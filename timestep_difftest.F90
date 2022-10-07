#include "config.h"
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
  public ::  step_all_grid
#ifdef SINGLE_STEP
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
    integer,intent(IN) :: level
    integer :: id
    CurrentLevel = level
    call timestep_init
    call update_timestep
    ! predictor stage
    STEP_MODE = PREDICTOR
    call boundary_cond
    do CurrentIndex = Gidmin, ListMax !各グリッド
       call backup_u_2order
       call convert_u2w         ! u -> w
    enddo
!!$    do CurrentIndex = Gidmin, ListMax !各グリッド
!!$       call get_flux            ! F
!!$       call w_update(Cdtime/2)  ! w -(F)-> wn
!!$#ifdef WITH_SELFGRAVITY
!!$       call source_g(Cdtime/2)    ! wn
!!$#endif
!!$       call convert_w2u         ! wn -> u
!!$       call rescueLev( CurrentLevel )
!!$    enddo
    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1

    ! corrector stage
    STEP_MODE = CORRECTOR
!!$    call boundary_cond
    do CurrentIndex = Gidmin, ListMax !各グリッド
       call backup_u_1order
    enddo
!!$    do CurrentIndex = Gidmin, ListMax !各グリッド
!!$       call get_flux            ! F
!!$       call w_update(Cdtime)    ! w -(F)-> wn
!!$#if defined(WITH_SELFGRAVITY) && !defined(SINGLE_STEP)
!!$       call source_g(Cdtime)      ! wn
!!$#endif
!!$       call convert_w2u         ! wn -> u
!!$       Fp => F
!!$       call save_flux( Ulist(CurrentIndex), Fp )
!!$       call rescueLev( CurrentLevel )
!!$    enddo
    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1
    STEP_MODE = COMPLETE
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
    call u2w( u, w, get_dv(CurrentLevel) )
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
    call w2u( wn, u, get_dv(CurrentLevel) )
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
  ! update w by numerical flux
  !-----------------------------------------------------------------------
  ! get numerical flux in 3dim components
  !-----------------------------------------------------------------------
#define FLMT(x,y) max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x))
  subroutine get_flux
    use eos
    use util, only : util_arroffset
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: f1d, ul, ur
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: n, io,jo,ko,i2,j2,k2,i,j,k,m
    u => get_up( Ulist(CurrentIndex) )


#ifdef MHD
    do n = MX, MZ
       call flux(u, f1d, n)
       F(:,:,:,:,n) = f1d(:,:,:,:)
    enddo
#else !MHD
    do n = MX, MZ
       call util_arroffset(n,io,jo,ko)
       i2 = io*2
       j2 = jo*2
       k2 = ko*2
       ul = Undefd
       ur = Undefd
       do m = Mmin, Mmax
          do k = Kmin-ko, Kmax
             do j = Jmin-jo, Jmax
                do i = Imin-io, Imax
                   ul(i,j,k,m) = u(i,j,k,m) &
                        + (FLMT(u(i+io,j+jo,k+ko,m)-u(i,j,k,m), u(i,j,k,m)-u(i-io,j-jo,k-ko,m)))/2
                   ur(i,j,k,m) = u(i+io,j+jo,k+ko,m) &
                        - (FLMT(u(i+io,j+jo,k+ko,m)-u(i,j,k,m), u(i+i2,j+j2,k+k2,m)-u(i+io,j+jo,k+ko,m)))/2
                enddo
             enddo
          enddo
       enddo
       mcycle = cyclecomp( n )
       ul = ul(:,:,:,mcycle)
       ur = ur(:,:,:,mcycle)
       call flux(ul, ur, f1d)
       F(:,:,:,mcycle,n) = f1d(:,:,:,:)
    enddo
#endif !MHD
#ifdef Emulate_1Dim
    F(:,:,:,:,MY) = 0.d0
#endif !Emulate_1Dim
#ifdef EMULATE_2DIM
    F(:,:,:,:,MZ) = 0.d0
#endif !EMULATE_2DIM
  end subroutine get_flux
#undef FLMT
  !-----------------------------------------------------------------------
  subroutine w_update( dt )
    use grid
    use util, only : util_arroffset
#ifdef MHD
    use eos, only : source_b
#endif
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
  end subroutine w_update
  !-----------------------------------------------------------------------
  ! add source term to wn
  !-----------------------------------------------------------------------
#ifdef WITH_SELFGRAVITY
  subroutine source_g( dt )
    use grid
#ifdef SINKPARTICLE
    use sinkParticle, only : sp_gravity
#endif !SINKPARTICLE
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: wn
    real(kind=DBL_KIND) :: dv, dtdvrho
    integer :: i,j,k,m,n
#ifdef SINKPARTICLE
    call sp_gravity(Ulist(CurrentIndex), activate=1)
#endif !SINKPARTICLE
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
#ifdef MP
             wn(i,j,k,MP)  = wn(i,j,k,MP) + dtdvrho*( &
                  u(i,j,k,MVX)*u(i,j,k,MGX) &
                  +u(i,j,k,MVY)*u(i,j,k,MGY) &
                  +u(i,j,k,MVZ)*u(i,j,k,MGZ))
#endif !MP
          enddo
       enddo
    enddo
#ifdef SINKPARTICLE
    call sp_gravity(Ulist(CurrentIndex), activate=-1)
#endif !SINKPARTICLE
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
#ifdef SINKPARTICLE
    use sinkParticle, only : sp_gravity
#endif !SINKPARTICLE
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, u1
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: wn
    integer :: level, n, gid, ncrd, i, j, k, io, jo, ko
    real(kind=DBL_KIND) :: dtdvrho, dv
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          call boundary_g(gid) ! baundary condition for gravity vector
#ifdef SINKPARTICLE
          call sp_gravity(gid, activate=1)
#endif !SINKPARTICLE
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
#ifdef MP
                   wn(i,j,k,MP)  = wn(i,j,k,MP) + dtdvrho*( &
                         u1(i,j,k,MVX)*u(i,j,k,MGX) &
                        +u1(i,j,k,MVY)*u(i,j,k,MGY) &
                        +u1(i,j,k,MVZ)*u(i,j,k,MGZ))
#endif !MP
                enddo
             enddo
          enddo
#ifdef SINKPARTICLE
          call sp_gravity(gid, activate=-1)
#endif !SINKPARTICLE
          call w2u( wn, u, dv )
       enddo
    enddo
  end subroutine source_g_all_level
#endif  !defined(WITH_SELFGRAVITY) && defined(SINGLE_STEP)
  !-----------------------------------------------------------------------
end module timestep
