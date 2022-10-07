#include "config.h"

!activating debug output
!#define KS_DEBUG

!-------------------------------------------------------------------------
! Module for multi-timestep
!-------------------------------------------------------------------------
module timestep
  use grid
  implicit none
  private
  integer,save,private :: STEP_MODE          ! predictor or corrector?
  integer,save,private :: CurrentLevel       ! 現在のレベル
  integer,save,private :: CurrentIndex       ! CurrentGridId = Wlist(CurrentIndex)
  ! CurrentLevel の時間
  real(kind=DBL_KIND),save :: Ctime    ! time
  real(kind=DBL_KIND),save :: Cdtime   ! 時間進み幅
  integer(kind=LLONG_KIND),save :: Cstep    ! step number
  integer(kind=LLONG_KIND),save :: Cdstep   ! delta step number
  ! 親レベルの時間
  real(kind=DBL_KIND),save :: Ptime    ! time
  real(kind=DBL_KIND),save :: Pdtime   ! 時間進み幅
  integer(kind=LLONG_KIND),save :: Pstep    ! step number
  integer(kind=LLONG_KIND),save :: Pdstep   ! delta step number
  ! 流束
  real(kind=DBL_KIND),save,dimension(ARRAYSIZE_IJKMGH,MX:MZ),target :: F
  real(kind=DBL_KIND),save,dimension(:,:,:,:,:),pointer :: Fp
  ! 作業変数のリスト
  integer,save,dimension(Gidmin:Gidmax) :: Wlist, Wnlist, Ulist
  integer,save :: ListMax

  !KS ADDED
#if MODEL_ART > 0
  real(kind=DBL_KIND),save,dimension(ARRAYSIZE_IJKGH) :: gam
#endif !MODEL_ART    
  
  public ::  step_all_grid
#if defined(SINGLE_STEP) && defined(WITH_SELFGRAVITY)
  public :: source_g_all_level
#endif
contains
  !-----------------------------------------------------------------------
  ! initialize timestep module
  !-----------------------------------------------------------------------
  subroutine timestep_init
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
       Wlist(n) = alloc_block()
       Wnlist(n) = alloc_block()
       Ulist(n) = GidList(n,  CurrentLevel )
    enddo
    ListMax = GidListMax( CurrentLevel )
  end subroutine timestep_init
  !-----------------------------------------------------------------------
  ! finalize module timestep
  !-----------------------------------------------------------------------
  subroutine timestep_finalize
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
       call dealloc_block(Wlist(n))
       call dealloc_block(Wnlist(n))
    enddo
  end subroutine timestep_finalize
  !-----------------------------------------------------------------------
  ! get w pointer
  !-----------------------------------------------------------------------
  function get_wp( n ) result( w )
    integer,intent(IN) :: n
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: w
    w => get_Up(Wlist(n))
  end function get_wp
  !-----------------------------------------------------------------------
  ! get wn pointer
  !-----------------------------------------------------------------------
  function get_wnp( n ) result( wn )
    integer,intent(IN) :: n
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: wn
    wn => get_Up(Wnlist(n))
  end function get_wnp
  !-----------------------------------------------------------------------
  ! ひとつのレベルを時間推進
  !-----------------------------------------------------------------------
  subroutine step_all_grid( level )
    use reflux
    use rescue
    use eos !KS ADDED
    integer,intent(IN) :: level
    integer :: id
    CurrentLevel = level

#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank) print *, "(step_all_grid, KS DEBUG) ",level     !KS DEBUG
#endif !KS_DEBUG
    call timestep_init
    call update_timestep
    ! predictor stage
    STEP_MODE = PREDICTOR
    call boundary_cond
    
    !たまにrefineでできた直後のセルを救済しなきゃいけないときがある
    call rescueLev( CurrentLevel ) ! KS ADDED

    do CurrentIndex = Gidmin, ListMax !各グリッド
       !--------------- KS DEBUG -------------!
       globdbg_mygid = Ulist(CurrentIndex)
       !--------------- KS DEBUG -------------!

#if MODEL_ART > 0
       call get_gamma(Ulist(CurrentIndex),gam)        !get cell-dependent gamma (KS ADDED)
#endif !MODEL_ART    
       
       call backup_u_2order
       call convert_u2w         ! u -> w
    enddo
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank)  print *, "(KS DEBUG A) "     !KS DEBUG
#endif !KS_DEBUG
    do CurrentIndex = Gidmin, ListMax !各グリッド
       !--------------- KS DEBUG -------------!
       globdbg_mygid = Ulist(CurrentIndex)
       !--------------- KS DEBUG -------------!

#if MODEL_ART > 0
       call get_gamma(Ulist(CurrentIndex),gam)         !get cell-dependent gamma (KS ADDED)
#endif !MODEL_ART    

       call get_flux            ! F
#ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid)  print *, "(KS DEBUG B) "     !KS DEBUG
#endif !KS_DEBUG
       call w_update(Cdtime/2)  ! w -(F)-> wn
#ifdef WITH_SELFGRAVITY
       call source_g(Cdtime/2)    ! wn
#endif
#ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid)  print *, "(KS DEBUG C) "     !KS DEBUG
#endif !KS_DEBUG
       call convert_w2u         ! wn -> u
    enddo

    ! t+dt/2 で壊れたセルを救済
    call rescueLev( CurrentLevel )

    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1

    ! corrector stage
    STEP_MODE = CORRECTOR
    call boundary_cond

    ! boundary_condで壊れる（ゴースト）セルもあるっぽい
    call rescueLev( CurrentLevel ) ! KS ADDED

    do CurrentIndex = Gidmin, ListMax !各グリッド
       call backup_u_1order
    enddo
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank)  print *, "(KS DEBUG D) "     !KS DEBUG
#endif !KS_DEBUG
    do CurrentIndex = Gidmin, ListMax !各グリッド
       !--------------- KS DEBUG -------------!
       globdbg_mygid = Ulist(CurrentIndex)
       !--------------- KS DEBUG -------------!

#if MODEL_ART > 0
       call get_gamma(Ulist(CurrentIndex),gam) !get cell-dependent gamma (KS ADDED)
#endif !MODEL_ART    

       call get_flux            ! F
#ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid)  print *, "(KS DEBUG E) "     !KS DEBUG
#endif !KS_DEBUG
       call w_update(Cdtime)    ! w -(F)-> wn
#if defined(WITH_SELFGRAVITY) && !defined(SINGLE_STEP)
       call source_g(Cdtime)      ! wn
#endif
#ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid)  print *, "(KS DEBUG F) "     !KS DEBUG
#endif !KS_DEBUG
       call convert_w2u         ! wn -> u
       Fp => F
       call save_flux( Ulist(CurrentIndex), Fp )
    enddo
    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1
    STEP_MODE = COMPLETE
    call boundary_cond          ! rescue bounary
    call rescueLev( CurrentLevel )
    call timestep_finalize
  end subroutine step_all_grid
  !-----------------------------------------------------------------------
  ! update timestep module parameter
  !-----------------------------------------------------------------------
  subroutine update_timestep
    Ctime  = Time(  CurrentLevel )
    Cdtime = Dtime( CurrentLevel )
    Cstep  = Step(  CurrentLevel )
    Cdstep = Dstep( CurrentLevel )
    if ( CurrentLevel > 0 ) then
       Ptime  = Time(  CurrentLevel-1 )
       Pdtime = Dtime( CurrentLevel-1 )
       Pstep  = Step(  CurrentLevel-1 )
       Pdstep = Dstep( CurrentLevel-1 )
    endif
    Time( CurrentLevel ) = Ctime + Cdtime
    Step( CurrentLevel ) = Cstep + Cdstep
  end subroutine update_timestep
  !-------------------------------------------------------------------------
  ! bakup to second order
  !-------------------------------------------------------------------------
  subroutine backup_u_2order
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, u2
    u => get_up( Ulist(CurrentIndex ) )
    u2 => get_u2orderp( Ulist(CurrentIndex) )
    u2 = u
    U2_StepNumber(CurrentLevel)          = U_StepNumber(CurrentLevel)
    U2_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_2order
  !-------------------------------------------------------------------------
  ! bakup to second order
  !-------------------------------------------------------------------------
  subroutine backup_u_1order
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, u1
    u => get_up( Ulist(CurrentIndex) )
    u1 => get_u1orderp( Ulist(CurrentIndex) )
    u1 = u
    U1_StepNumber(CurrentLevel)          = U_StepNumber(CurrentLevel)
    U1_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_1order
  !-----------------------------------------------------------------------
  ! convert u2w
  !-----------------------------------------------------------------------
  subroutine convert_u2w
    use eos
    use grid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: w
    u => get_up( Ulist(CurrentIndex) )
    w => get_wp( CurrentIndex )
    ! call u2w( u, w, get_dv(CurrentLevel) )
    call u2w_withgam( u, w, get_dv(CurrentLevel), gam ) !KS MODIFIED
  end subroutine convert_u2w
  !-----------------------------------------------------------------------
  ! convert u2w
  !-----------------------------------------------------------------------
  subroutine convert_w2u
    use eos
    use grid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: wn
    u => get_up( Ulist(CurrentIndex) )
    wn => get_wnp( CurrentIndex )
    ! call w2u( wn, u, get_dv(CurrentLevel) )
    call w2u_withgam( wn, u, get_dv(CurrentLevel), gam ) !KS MODIFIED
  end subroutine convert_w2u
  !-----------------------------------------------------------------------
  ! boundary condition
  !-----------------------------------------------------------------------
  subroutine boundary_cond
    use grid_boundary
    use boundary
    integer :: n
    call boundary_grid( CurrentLevel, STEP_MODE )
    do n = Gidmin, ListMax !各グリッド
       call boundary_u( Ulist(n), STEP_MODE)
    enddo
  end subroutine boundary_cond
  !-----------------------------------------------------------------------
  ! get numerical flux in 3dim components
  !-----------------------------------------------------------------------
  ! minmod limiter
#define MINMOD(x, y) (max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x)))
#define FLMT(x, y) MINMOD(x, y)
  subroutine get_flux
    use eos
    use util, only : util_arroffset
#ifdef FLUX_ROEM2
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH,MX:MZ) :: ul, ur
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH,MX:MZ) :: pratio
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: f1d, ul1d, ur1d
#else !FLUX_ROEM2
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: f1d, ul, ur
#endif !FLUX_ROEM2
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: n, io,jo,ko,i2,j2,k2,i,j,k,m

#ifdef RECONSTRUCTION_MUSCL3
    real(kind=DBL_KIND),parameter :: ETA = 1.d0/3.d0 ! 1/3 for 3rd order accuracy in space
    real(kind=DBL_KIND),parameter :: BW = (3.d0-ETA)/(1.d0-ETA)
!!$    real(kind=DBL_KIND),parameter :: ETA = -1.d0
!!$    real(kind=DBL_KIND),parameter :: BW = 1.d0
    real(kind=DBL_KIND) :: dva, dvb
#endif !RECONSTRUCTION_MUSCL3
    ! -------------------------
    ! definition of UL and UR
    ! -------------------------
#ifdef FLUX_ROEM2
#define UL_ ul(i,j,k,m,n)
#define UR_ ur(i,j,k,m,n)
#else !FLUX_ROEM2
#define UL_ ul(i,j,k,m)
#define UR_ ur(i,j,k,m)
#endif !FLUX_ROEM2

    u => get_up( Ulist(CurrentIndex) )

#ifdef MHD_ROE
    do n = MX, MZ
       call flux(u, f1d, n)
       F(:,:,:,:,n) = f1d(:,:,:,:)
    enddo
#else !MHD_ROE
#ifdef FLUX_ROEM2
    ul = Undefd
    ur = Undefd
#endif !FLUX_ROEM2
    do n = MX, MZ
       call util_arroffset(n,io,jo,ko)
       i2 = io*2
       j2 = jo*2
       k2 = ko*2
#ifndef FLUX_ROEM2
       ul = Undefd
       ur = Undefd
#endif !FLUX_ROEM2
       do m = Mmin, Mmax
          do k = Kmin-ko, Kmax
             do j = Jmin-jo, Jmax
                do i = Imin-io, Imax
#if defined(RECONSTRUCTION_NONE)
                   UL_ = u(i,j,k,m)
                   UR_ = u(i+io,j+jo,k+ko,m)
#elif defined(RECONSTRUCTION_MUSCL3)
                   dva = u(i+io,j+jo,k+ko,m) - u(i,j,k,m)
                   dvb = u(i,j,k,m) - u(i-io,j-jo,k-ko,m)
                   UL_ = u(i,j,k,m) &
                        + (1.d0 - ETA)/4.d0 * MINMOD(dvb, BW*dva) &
                        + (1.d0 + ETA)/4.d0 * MINMOD(dva, BW*dvb)
                   dva = u(i+i2,j+j2,k+k2,m) - u(i+io,j+jo,k+ko,m)
                   dvb = u(i+io,j+jo,k+ko,m) - u(i,j,k,m)
                   UR_ = u(i+io,j+jo,k+ko,m) &
                        - (1.d0 - ETA)/4.d0 * MINMOD(dva, BW*dvb) &
                        - (1.d0 + ETA)/4.d0 * MINMOD(dvb, BW*dva)
#else !RECONSTRUCTION_MUSCL2 (default)
                   UL_ = u(i,j,k,m) &
                        + (FLMT(u(i+io,j+jo,k+ko,m)-u(i,j,k,m), u(i,j,k,m)-u(i-io,j-jo,k-ko,m)))*0.5d0
                   UR_ = u(i+io,j+jo,k+ko,m) &
                        - (FLMT(u(i+io,j+jo,k+ko,m)-u(i,j,k,m), u(i+i2,j+j2,k+k2,m)-u(i+io,j+jo,k+ko,m)))*0.5d0
#endif

                enddo
             enddo
          enddo
       enddo
#ifdef FLUX_ROEM2
       !------------- KS DEBUG ---------------!
  #ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
          i=globdbg_i+lbound(u,1)
          j=globdbg_j+lbound(u,2)
          k=globdbg_k+lbound(u,1)
          print '(A,6I8,/,1P6E15.7)', "(get_flux, KS DEBUG)" , &
               globdbg_myrank, globdbg_mygid, i-lbound(u,1), j-lbound(u,2), k-lbound(u,3),n,&
               ul(i,j,k,MRHO,n),ur(i,j,k,MRHO,n),u(i-io,j-jo,k-ko,MRHO),&
               u(i,j,k,MRHO),u(i+io,j+jo,k+ko,MRHO),u(i+i2,j+j2,k+k2,MRHO)
          print '(1P6E15.7)', &
               ul(i,j,k,MP,n),ur(i,j,k,MP,n),u(i-io,j-jo,k-ko,MP),&
               u(i,j,k,MP),u(i+io,j+jo,k+ko,MP),u(i+i2,j+j2,k+k2,MP)
          print '(1P27E12.4,/,1P27E12.4,/,1P27E12.4)', &
               u(i-1:i+1,j-1:j+1,k-1:k+1,MP),&
               ul(i-1:i+1,j-1:j+1,k-1:k+1,MP,n),&
               ur(i-1:i+1,j-1:j+1,k-1:k+1,MP,n)
       end if
  #endif !KS_DEBUG
       !------------- KS DEBUG ---------------!
    enddo
    pratio = min(ul(:,:,:,MP,:)/ur(:,:,:,MP,:),ur(:,:,:,MP,:)/ul(:,:,:,MP,:))
    do n = MX, MZ
#endif !FLUX_ROEM2
       mcycle = cyclecomp( n )
#ifdef FLUX_ROEM2
       ul1d = ul(:,:,:,mcycle,n)
       ur1d = ur(:,:,:,mcycle,n)
!KS MODIFIED
  #if MODEL_ART > 0
       call flux(ul1d, ur1d, pratio, f1d, n, gam) !gammaを渡す
  #else !MODEL_ART
       call flux(ul1d, ur1d, pratio, f1d, n)
  #endif !MODEL_ART


#else !FLUX_ROEM2p
       ul = ul(:,:,:,mcycle)
       ur = ur(:,:,:,mcycle)

  #if MODEL_ART > 0
       call flux(ul, ur, f1d, gam)
  #else !MODEL_ART
       call flux(ul, ur, f1d)
  #endif !MODEL_ART


#endif !FLUX_ROEM2
       F(:,:,:,mcycle,n) = f1d(:,:,:,:)
    enddo
#endif !MHD_ROE
#ifdef Emulate_1Dim
    F(:,:,:,:,MY) = 0.d0
#endif !Emulate_1Dim
#ifdef EMULATE_2DIM
    F(:,:,:,:,MZ) = 0.d0
#endif !EMULATE_2DIM
  end subroutine get_flux
#undef FLMT
#undef UL_
#undef UR_
  !-----------------------------------------------------------------------
  subroutine w_update( dt )
    use grid
    use util, only : util_arroffset
#ifdef MHD
    use eos, only : source_b
#endif
#ifdef EXTERNALFORCE
    use externalForce
#endif
    use mpilib ! KS DEBUG
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(MX:MZ) :: ds
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: w
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: wn
    integer :: is,ie,js,je,ks,ke, i,j,k,m,n, io,jo,ko
    w => get_wp( CurrentIndex )
    wn => get_wnp( CurrentIndex )
    ds = get_ds( CurrentLevel )
    wn = w

    do n=MX,MZ
       call util_arroffset(n,io,jo,ko)
       !------------- KS DEBUG ---------------!
#ifdef KS_DEBUG
       if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
          i=globdbg_i+lbound(w,1)
          j=globdbg_j+lbound(w,2)
          k=globdbg_k+lbound(w,1)
          m=MRHO
          print '(A,6I8,/,1P6E15.7)', "(w_update, KS DEBUG)",&
               globdbg_myrank, globdbg_mygid, i-lbound(w,1), j-lbound(w,2), k-lbound(w,3),n,&
               wn(i,j,k,m),wn(i,j,k,m)-dt*ds(n)*(F(i,j,k,m,n)-F(i-io,j-jo,k-ko,m,n)),&
               dt, ds(n), F(i,j,k,m,n),F(i-io,j-jo,k-ko,m,n)
          m=MP
          print '(1P6E15.7)',wn(i,j,k,m),wn(i,j,k,m)-dt*ds(n)*(F(i,j,k,m,n)-F(i-io,j-jo,k-ko,m,n)),&
               dt, ds(n), F(i,j,k,m,n),F(i-io,j-jo,k-ko,m,n)
       end if
#endif !KS_DEBUG
       !------------- KS DEBUG ---------------!

       do m = Mmin, Mmax
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   wn(i,j,k,m) = wn(i,j,k,m) - dt*ds(n)*(F(i,j,k,m,n)-F(i-io,j-jo,k-ko,m,n))
                enddo
             enddo
          enddo
       enddo
    enddo
#ifdef MHD
    call source_b(F, wn, dt, Ulist(CurrentIndex))
#endif
#ifdef EXTERNALFORCE
    Fp => F
    call source_externalForce(wn, dt, Ulist(CurrentIndex), Fp)
#endif
  end subroutine w_update
  !-----------------------------------------------------------------------
  ! add source term to wn
  !-----------------------------------------------------------------------
#ifdef WITH_SELFGRAVITY
  subroutine source_g( dt )
    use grid
!!$#ifdef SINKPARTICLE
!!$    use sinkParticle, only : sp_gravity
!!$#endif !SINKPARTICLE
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: wn
    real(kind=DBL_KIND) :: dv, dtdvrho
    integer :: i,j,k,m,n
!!$#ifdef SINKPARTICLE
!!$    call sp_gravity(Ulist(CurrentIndex), activate=1)
!!$#endif !SINKPARTICLE
    u => get_up( Ulist(CurrentIndex) )
    wn => get_wnp( CurrentIndex )
    dv = get_dv( CurrentLevel )
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             dtdvrho = dt*dv*u(i,j,k,MRHO)
             wn(i,j,k,MVX) = wn(i,j,k,MVX) + u(i,j,k,MGX)*dtdvrho
             wn(i,j,k,MVY) = wn(i,j,k,MVY) + u(i,j,k,MGY)*dtdvrho
             wn(i,j,k,MVZ) = wn(i,j,k,MVZ) + u(i,j,k,MGZ)*dtdvrho
#ifndef OMIT_WORK_OF_GRAVITY
#ifdef MP
             wn(i,j,k,MP)  = wn(i,j,k,MP) + dtdvrho*( &
                  u(i,j,k,MVX)*u(i,j,k,MGX) &
                  +u(i,j,k,MVY)*u(i,j,k,MGY) &
                  +u(i,j,k,MVZ)*u(i,j,k,MGZ))
#endif !MP
#endif !OMIT_WORK_OF_GRAVITY


          enddo
       enddo
    enddo
!!$#ifdef SINKPARTICLE
!!$    call sp_gravity(Ulist(CurrentIndex), activate=-1)
!!$#endif !SINKPARTICLE
  end subroutine source_g
#endif !WITH_SELFGRAVITY
#if defined(WITH_SELFGRAVITY) && defined(SINGLE_STEP)
  ! -----------------------------------------------------------------
  ! 以下は、SINGLE_STEPが定義されている場合に用いるコード。
  ! レベルごとに共通のタイムステップを刻む場合。
  ! corrector step で時間２次精度自己重力を評価する。
  !-----------------------------------------------------------------------
  subroutine source_g_all_level
    use eos
    use boundary
!!$#ifdef SINKPARTICLE
!!$    use sinkParticle, only : sp_gravity
!!$#endif !SINKPARTICLE
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, u1
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: wn
    integer :: level, n, gid, ncrd, i, j, k, io, jo, ko
    real(kind=DBL_KIND) :: dtdvrho, dv, dt
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          call boundary_g(gid) ! baundary condition for gravity vector
!!$#ifdef SINKPARTICLE
!!$          call sp_gravity(gid, activate=1)
!!$#endif !SINKPARTICLE
          ! ----------------------
          ! source term
          ! ----------------------
          u1 => get_u1orderp( gid ) ! u at predictor stage
          u => get_Up( gid )
          dv = get_dv( level )
          call u2w( u, wn, dv )
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   dtdvrho = Dtime(level)*dv*u1(i,j,k,MRHO)
                   wn(i,j,k,MVX) = wn(i,j,k,MVX) + u(i,j,k,MGX)*dtdvrho
                   wn(i,j,k,MVY) = wn(i,j,k,MVY) + u(i,j,k,MGY)*dtdvrho
                   wn(i,j,k,MVZ) = wn(i,j,k,MVZ) + u(i,j,k,MGZ)*dtdvrho
#ifndef OMIT_WORK_OF_GRAVITY
#ifdef MP
                   wn(i,j,k,MP)  = wn(i,j,k,MP) + dtdvrho*( &
                         u1(i,j,k,MVX)*u(i,j,k,MGX) &
                        +u1(i,j,k,MVY)*u(i,j,k,MGY) &
                        +u1(i,j,k,MVZ)*u(i,j,k,MGZ))
#endif !MP
#endif !OMIT_WORK_OF_GRAVITY


                enddo
             enddo
          enddo
!!$#ifdef SINKPARTICLE
!!$          call sp_gravity(gid, activate=-1)
!!$#endif !SINKPARTICLE
          call w2u( wn, u, dv )
       enddo
    enddo
  end subroutine source_g_all_level
#endif  !defined(WITH_SELFGRAVITY) && defined(SINGLE_STEP)
  !-----------------------------------------------------------------------
end module timestep
