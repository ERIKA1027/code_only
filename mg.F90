#include "config.h"
#include "debug_fmg.h"
! Show error
! #define VERBOSE
!-------------------------------------------------------------------------
! Module for uniform Multigrid
!-------------------------------------------------------------------------
module mg
  use mg_data
  use fmg_boundary_phys
  implicit none
  real(kind=DBL_KIND),save :: Resmaxg
  private
  public :: mg_lin, mg_cycle, mg_boundary_fill0, mg_fill0, mg_nonlin, mg_nonlin_vcycle, mg_init, mg_finalize
contains

#ifdef FMG_POISSON
#include "mg_poisson.F90"
#endif !FMG_POISSON

#ifdef FMG_OHMIC_DISSIPATION
#include "mg_od.F90"
#endif !FMG_OHMIC_DISSIPATION

#ifdef FMG_DIFFUSION
#include "mg_diffusion.F90"
#endif !FMG_DIFFUSION

#ifdef FMG_AMBIPOLAR_DIFFUSION
#include "mg_ad.F90"
#endif !FMG_AMBIPOLAR_DIFFUSION

  !-------------------------------------------------------------------------
  ! IN .... PSI/U (initial guess), SRC/RHO (source)
  ! OUT ... PSI/U
  ! Given the initial guess of ju and jrho, this routine returns ju.
  ! This routine uses SRC and PSI templarily.
  !
  ! Parameter: nmax, errormax
  !  vmg の下請け業務としては、nmg = 1 程度で問題ない。
  !-------------------------------------------------------------------------
  subroutine mg_cycle(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: ju, jrho
    real(kind=DBL_KIND),parameter :: errormax = 1.e-2
!!$    real(kind=DBL_KIND),parameter :: errormax = 1.e-15
    integer,parameter :: nmg = 10
    real(kind=DBL_KIND) :: rhomax
    integer ::  n, mglev
    mglev = MG_LevelMin
    ! ----------------------------------
    ! 入力データを IPSI, ISRC に統一する
    ! ----------------------------------
    if (ju /= IPSI) then
       call mg_alloc_arr(mglev, IPSI)
       call mg_copy(mglev, IPSI,  ju)   ! backup u -> psi
    end if
    if (jrho /= ISRC) then
       call mg_alloc_arr(mglev, ISRC)
       call mg_copy(mglev, ISRC, jrho)   ! backup rho -> src
    end if
    call mg_boundary_fill0(mglev, IPSI) ! 孤立境界に必要。
    call mg_boundary_u(mglev, IPSI)     ! 周期境界等で上書き
    ! ----------------
    ! Multigrid cycle
    ! ----------------
    call mg_alloc_arr(mglev, IRHO)
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       ! nothing to do
    elseif ( &
         FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call mg_alloc_arr(mglev, IETA)
       call mg_boundary_u(mglev, IETA)
    else
       print *, '*** mg_cycle: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    do n = 1, nmg
       call mg_resid(mglev, IRHO, IPSI, ISRC)
       if (n == 1) rhomax = mg_get_absmax(mglev, IRHO)
       call mg_lin                 ! solve IU by IRHO
       call mg_add(mglev, IPSI, IU)
#ifdef VERBOSE
       print *, 'MG error', Resmaxg / rhomax
#endif !VERBOSE
       if ( Resmaxg / rhomax < errormax ) exit
    end do
    if (ju /= IPSI) call mg_copy(mglev, ju, IPSI)       ! restore U
  end subroutine mg_cycle
  !-------------------------------------------------------------------------
  ! solve Poisson eq. of u, given rho (source).
  ! This subroutine is given by `Mutigrid Methods for Boundary Value
  ! Problems. I. W. H.Press & S.A.Teukolsky, Computer in Physics, 5
  !  514-519(1991)
  !
  ! PARAMETER:
  !   NG = the number of grids hierarchy.
  !   Npre = the number of pre-smoothing iteration.
  !   Npre = the number of post-smoothing iteration.
  !-------------------------------------------------------------------------
  subroutine mg_lin
    use string
    use io_util
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
!!$    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    INTEGER :: Ncycle, Npre, Npost
    integer :: mglev, jcycle, jpre, jpost, vlevel
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       Ncycle=2
       Npre=2
       Npost=2
    else
       Ncycle=2
       Npre=2
       Npost=2
    endif
    ! --------------------------
    ! 最細レベルの係数の境界条件
    ! --------------------------
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
       ! nothing to do
    elseif ( &
         FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or.&
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    else
       print *, '*** mg_lin: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    ! --------------------------
    ! rho を各レベルへremap
    ! --------------------------
    do mglev = MG_LevelMin+1, MG_LevelMax
       call mg_alloc_arr(mglev, IRHO)
       call mg_rstrct(mglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
          ! nothing to do
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       else
          print *, '*** mg_lin: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
#ifdef DEBUG_VMG_OUTPUT_RES
    call mg_alloc_arr(MG_LevelMin, IDBG)
#endif !DEBUG_VMG_OUTPUT_RES
    ! ----------------
    ! 最粗レベルで解く
    ! ----------------
    mglev = MG_LevelMax
    call mg_alloc_arr(mglev, IU)
    call mg_alloc_arr(mglev, IRHS)
    call mg_fill0(mglev, IU)    ! initial guess of IU
    call mg_slvsml(mglev, IU, IRHO)
    ! ----------------------
    ! 粗レベルから細レベルへ
    ! ----------------------
    do mglev = MG_LevelMax-1, MG_LevelMin, -1
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRES)
       call mg_interp(mglev, IU, IU, cubic=.TRUE.)
       call mg_copy(mglev, IRHS, IRHO)

       do jcycle = 1, Ncycle     ! Ncycle
          do vlevel = mglev, MG_LevelMax-1 ! 下り
             do jpre=1,Npre
                call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
             enddo
             call mg_resid(vlevel, IRES, IU, IRHS)
             call mg_rstrct(vlevel+1, IRHS, IRES)
             call mg_fill0(vlevel+1, IU)
          enddo
          vlevel = MG_LevelMax
#ifdef VERBOSE
          call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
          call mg_slvsml(vlevel, IU, IRHS)
#ifdef VERBOSE
          call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
          do vlevel = MG_LevelMax-1, mglev, -1 ! 上り
             call mg_addint(vlevel, IU, IU, IRES)
             do jpost = 1, Npost
                call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
             enddo
          enddo
#ifdef VERBOSE
          call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' at level = '// trim(num2char(mglev)) )
#endif !VERBOSE
       enddo
    enddo

!!$    print *, 'resmax', Resmaxg
  end subroutine mg_lin
  !-------------------------------------------------------------------------
  ! IN .... PSI/U (initial guess), SRC/RHO (source)
  ! OUT ... PSI/U
  ! Given the initial guess of ju and jrho, this routine returns ju.
  !
  ! solve nonlinear PDE for u, given rho (source) by FMG-FAS.
  ! This subroutine is given by `Mutigrid Methods for Boundary Value
  ! Problems. I. W. H.Press & S.A.Teukolsky, Computer in Physics, 5
  !  526-529(1991)
  !
  ! PARAMETER:
  !   NG = the number of grids hierarchy.
  !   Npre = the number of pre-smoothing iteration.
  !   Npre = the number of post-smoothing iteration.
  !-------------------------------------------------------------------------
  subroutine mg_nonlin(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
!!$    INTEGER,parameter :: Ncycle=10, Npre=10, Npost=10
    integer :: mglev, jcycle, jpre, jpost, vlevel, icode
!!$    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d-20
    real(kind=DBL_KIND) :: trerr
    ! ----------------------------------
    ! 入力データを IU, IRHO に統一する
    ! ----------------------------------
    mglev = MG_LevelMin
    if (ju /= IU) then
       call mg_alloc_arr(mglev, IU)
       call mg_copy(mglev, IU,  ju)
    end if
    if (jrho /= IRHO) then
       call mg_alloc_arr(mglev, IRHO)
       call mg_copy(mglev, IRHO, jrho)
    end if
    ! --------------------------
    ! 最細レベルの係数の境界条件
    ! --------------------------
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call mg_boundary_u(MG_LevelMin, IDOD)
       call mg_boundary_u(MG_LevelMin, IDHE)
       call mg_boundary_u(MG_LevelMin, IDAD)
    else
       print *, '*** mg_nonlin: this type of PDE is not supported', FMG_PDE_TYPE
    endif
#ifdef DEBUG_VMG_OUTPUT_RES
    call mg_alloc_arr(MG_LevelMin, IDBG)
#endif !DEBUG_VMG_OUTPUT_RES
    ! ---------------------------
    ! rho, 係数 を各レベルへremap
    ! ---------------------------
    do mglev = MG_LevelMin+1, MG_LevelMax
       call mg_alloc_arr(mglev, IRHO)
       call mg_rstrct(mglev, IRHO, IRHO)
       call mg_alloc_arr(mglev, IU) !初期値もrestrict
       call mg_rstrct(mglev, IU, IU)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
          ! nothing to do
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
          do icode = IDOD, IDAD
             call mg_alloc_arr(mglev, icode)
             call mg_rstrct(mglev, icode, icode)
             call mg_boundary_u(mglev, icode)
          end do
       else
          print *, '*** mg_nonlin: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
    ! ----------------
    ! 最粗レベルで解く
    ! ----------------
    mglev = MG_LevelMax
    call mg_slvsml(mglev, IU, IRHO)
    ! ----------------------
    ! 粗レベルから細レベルへ
    ! ----------------------
    call mg_alloc_arr(mglev, IRHS)
    call mg_alloc_arr(mglev, IRES)
    do mglev = MG_LevelMax-1, MG_LevelMin, -1
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRES)
       call mg_alloc_arr(mglev+1, IRUF)
       call mg_interp(mglev, IU, IU, cubic=.TRUE.)
       call mg_copy(mglev, IRHS, IRHO)
       do jcycle = 1, Ncycle     ! Ncycle
          do vlevel = mglev, MG_LevelMax-1 ! 下り
             do jpre=1,Npre
                call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
             enddo
             call mg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) ! Uc = R Uf
             call mg_copy(vlevel+1, IRUF, IU)      ! backup Uc (=R Uf) for GCG
             call mg_rstrct(vlevel+1, IRHS, IRHS)  ! RHSc = R RHSf
             call mg_tau(vlevel+1, IRES, IU)       ! tau-correction (RES is tempolary used for tau).
             call mg_add(vlevel+1, IRHS, IRES)     ! RHS = RHS + tau
             if (vlevel == mglev) trerr = ALPHA*mg_get_absmax(vlevel+1, IRES)
          enddo
          vlevel = MG_LevelMax
#ifdef VERBOSE
          call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
          call mg_slvsml(vlevel, IU, IRHS)
#ifdef VERBOSE
          call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
          do vlevel = MG_LevelMax-1, mglev, -1 ! 上り
             call mg_sub(vlevel+1, IRES, IU, IRUF) ! RESc = Uc - R Uf
             call mg_interp(vlevel, IRES, IRES, cubic=.TRUE.)    ! RESf = P (RESc) = P (Uc - R Uf)
             call mg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
             do jpost = 1, Npost
                call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(mglev)))
#endif !VERBOSE
             enddo
          enddo
#ifdef VERBOSE
          call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
               // ' at level = '// trim(num2char(mglev)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
          ! error estimate: exit jcycle loop
          if (Resmaxg < trerr) exit
          if (trerr <= tiny(trerr)) exit
       enddo
    enddo
! #ifdef VERBOSE
    call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
         // ' at level = 0' // ', jcycle = '//trim(num2char(jcycle)) // ' (final)')
! #endif !VERBOSE
    if (ju /= IU) call mg_copy(MG_LevelMin, ju, IU)       ! restore U
  end subroutine mg_nonlin
  !-------------------------------------------------------------------------
  ! IN .... PSI/U (initial guess), SRC/RHO (source)
  ! OUT ... PSI/U
  ! Given the initial guess of ju and jrho, this routine returns ju.
  !
  ! solve nonlinear PDE for u, given rho (source) by FAS-v-cycle.
  ! This subroutine is given by `Mutigrid Methods for Boundary Value
  ! Problems. I. W. H.Press & S.A.Teukolsky, Computer in Physics, 5
  !  526-529(1991)
  !
  ! PARAMETER:
  !   NG = the number of grids hierarchy.
  !   Npre = the number of pre-smoothing iteration.
  !   Npre = the number of post-smoothing iteration.
  !-------------------------------------------------------------------------
  subroutine mg_nonlin_vcycle(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    INTEGER,parameter :: Ncycle=1, Npre=2, Npost=2
!!$    INTEGER,parameter :: Ncycle=10, Npre=10, Npost=10
    integer :: mglev, jcycle, jpre, jpost, vlevel, icode
!!$    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d-20
    real(kind=DBL_KIND) :: trerr
    ! ----------------------------------
    ! 入力データを IU, IRHO に統一する
    ! ----------------------------------
    mglev = MG_LevelMin
    if (ju /= IU) then
       call mg_alloc_arr(mglev, IU)
       call mg_copy(mglev, IU,  ju)
    end if
    if (jrho /= IRHO) then
       call mg_alloc_arr(mglev, IRHO)
       call mg_copy(mglev, IRHO, jrho)
    end if
    ! --------------------------
    ! 最細レベルの係数の境界条件
    ! --------------------------
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call mg_boundary_u(MG_LevelMin, IDOD)
       call mg_boundary_u(MG_LevelMin, IDHE)
       call mg_boundary_u(MG_LevelMin, IDAD)
    else
       print *, '*** mg_nonlin_vcycle: this type of PDE is not supported', FMG_PDE_TYPE
    endif
#ifdef DEBUG_VMG_OUTPUT_RES
    call mg_alloc_arr(MG_LevelMin, IDBG)
#endif !DEBUG_VMG_OUTPUT_RES
    ! ---------------------------
    ! 係数 を各レベルへremap
    ! ---------------------------
    do mglev = MG_LevelMin+1, MG_LevelMax
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
          ! nothing to do
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
          do icode = IDOD, IDAD
             call mg_alloc_arr(mglev, icode)
             call mg_rstrct(mglev, icode, icode)
             call mg_boundary_u(mglev, icode)
          end do
       else
          print *, '*** mg_nonlin_vcycle: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
    ! --------------------
    ! allocation of arrays
    ! --------------------
    do mglev = MG_LevelMin, MG_LevelMax
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRUF)
       call mg_alloc_arr(mglev, IRES)
    end do
    ! ------------
    ! v-cycle
    ! ------------
    call mg_copy(MG_LevelMin, IRHS, IRHO) ! set source to rhs
    do jcycle = 1, Ncycle     ! Ncycle
       do vlevel = MG_LevelMin, MG_LevelMax-1 ! 下り
          do jpre=1,Npre
             call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
             call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
                  ' at level = '// trim(num2char(vlevel)))
#endif !VERBOSE
          enddo
          call mg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) ! Uc = R Uf
          call mg_copy(vlevel+1, IRUF, IU)      ! backup Uc (=R Uf) for GCG
          call mg_rstrct(vlevel+1, IRHS, IRHS)  ! RHSc = R RHSf
          call mg_tau(vlevel+1, IRES, IU)       ! tau-correction (RES is tempolary used for tau).
          call mg_add(vlevel+1, IRHS, IRES)     ! RHS = RHS + tau
          if (vlevel == MG_LevelMin) trerr = ALPHA*mg_get_absmax(vlevel+1, IRES)
       enddo
       vlevel = MG_LevelMax
#ifdef VERBOSE
       call print_msg( 'MG error in pre-smooth  = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)))
#endif !VERBOSE
       call mg_slvsml(vlevel, IU, IRHS)
#ifdef VERBOSE
       call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
               ' at level = '// trim(num2char(vlevel)))
#endif !VERBOSE
       do vlevel = MG_LevelMax-1, MG_LevelMin, -1 ! 上り
          call mg_sub(vlevel+1, IRES, IU, IRUF) ! RESc = Uc - R Uf
          call mg_interp(vlevel, IRES, IRES, cubic=.TRUE.)    ! RESf = P (RESc) = P (Uc - R Uf)
          call mg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
          do jpost = 1, Npost
             call mg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
             call print_msg( 'MG error in post-smooth = ' // trim(num2char(Resmaxg)) // &
                  ' at level = '// trim(num2char(vlevel)))
#endif !VERBOSE
          enddo
       enddo
#ifdef VERBOSE
       call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
            // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
       ! error estimate: exit jcycle loop
       if (Resmaxg < trerr) exit
       if (trerr <= tiny(trerr)) exit
    enddo
! #ifdef VERBOSE
    call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
         // ' at level = 0' // ', jcycle = '//trim(num2char(jcycle)) // ' (final)')
! #endif !VERBOSE
    if (ju /= IU) call mg_copy(MG_LevelMin, ju, IU)       ! restore U
  end subroutine mg_nonlin_vcycle
  !-------------------------------------------------------------------------
  ! initialize
  !-------------------------------------------------------------------------
  subroutine mg_init(fmglev)
    integer,intent(IN) :: fmglev
    call mg_data_init(fmglev)
  end subroutine mg_init
  !-------------------------------------------------------------------------
  ! finalize
  !-------------------------------------------------------------------------
  subroutine mg_finalize
    use mg_interpol_cubic
    call mg_data_finalize
    call mg_interp_cubic_finalize
  end subroutine mg_finalize
  !-------------------------------------------------------------------------
  ! solve in the coarest level
  !-------------------------------------------------------------------------
  subroutine mg_slvsml(mglev, ju, jrhs)
    integer,intent(IN) :: mglev, ju, jrhs
    integer :: n
    do n = 1, 2
       call mg_relax(mglev, ju, jrhs)
    end do
  end subroutine mg_slvsml
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine mg_resid(mglev, jres, ju, jrhs)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=DBL_KIND) :: h2i

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call mg_poisson_resid(mglev, jres, ju, jrhs)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call mg_diff_resid(mglev, jres, ju, jrhs)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call mg_od_resid(mglev, jres, ju, jrhs)
#endif !FMG_OHMIC_DISSIPATION
    else
       print *, '*** mg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine mg_resid
  !-------------------------------------------------------------------------
  ! Restrction
  ! mglevel = mg level of coarse grid
  ! icode    = identifier of data set
  !-------------------------------------------------------------------------
  subroutine mg_rstrct(mglev, iuc, iuf, cubic)
    use restriction
    use fmg_data, only : FMG_PDE_LINEAR
    integer,intent(IN) :: mglev, iuc, iuf
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: uf, uc
    logical,intent(IN),optional :: cubic
    logical :: bool_cubic
    integer :: ics, jcs, kcs, ice, jce, kce, lf, lc
    bool_cubic = .FALSE.             ! default
    if (present(cubic)) bool_cubic = cubic
!!$    if (.not. FMG_PDE_LINEAR) bool_cubic = .TRUE. ! 非線型ではcubicを適用

    lf = mglev-1
    lc = mglev
    call mg_get_gridsize(lc, ics,jcs,kcs,ice,jce,kce)
    uf => mg_get_arrp(lf,iuf)
    uc => mg_get_arrp(lc,iuc)
    if (bool_cubic) call mg_boundary_u(lf, iuf)
    call rstrct(uc, uf, &
            ics, jcs, kcs, ice, jce, kce, &
            ics, jcs, kcs, cubic)
  end subroutine mg_rstrct
  !-------------------------------------------------------------------------
  ! Interpolation
  ! mglevel = mg level of original data (fine level)
  ! icode    = identifier of data set
  !-------------------------------------------------------------------------
  subroutine mg_interp(mglev, iuf, iuc, cubic)
    use mg_interpol_cubic
    integer,intent(IN) :: mglev, iuf, iuc
    logical,intent(IN),optional :: cubic
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: uf, uc
    logical :: bool_cubic
    integer :: lf, lc
    lf = mglev
    lc = mglev+1
    bool_cubic = .FALSE.             ! default
    if (present(cubic)) bool_cubic = cubic

    call mg_boundary_extrap(lc, iuc)
    call mg_boundary_u(lc,iuc)

    if (bool_cubic) then
       call mg_interp_cubic(lf, iuf, iuc)
    else
       call mg_interp_linear(lf, iuf, iuc)
    end if

    call mg_boundary_fill0(lf, iuf)
    call mg_boundary_fill0(lc, iuc)
    call mg_boundary_u(lf,iuf)
    call mg_boundary_u(lc,iuc)
  end subroutine mg_interp
  !-------------------------------------------------------------------------
  subroutine mg_interp_linear(mglev, iuf, iuc)
    use interpolation
    integer,intent(IN) :: mglev, iuf, iuc
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, &
         if, jf, kf, &
         ic, jc, kc, &
         ic1, jc1, kc1, &
         ifs, jfs, kfs, ife, jfe, kfe, m
    lf = mglev
    lc = mglev+1
    call mg_get_gridsize(lf, ifs,jfs,kfs,ife,jfe,kfe)
    uf => mg_get_arrp(lf,iuf)
    uc => mg_get_arrp(lc,iuc)
    call interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, ifs, jfs, kfs)
  end subroutine mg_interp_linear
  ! ----------------------------------------------------------------
  ! fill zero
  ! ----------------------------------------------------------------
  subroutine mg_fill0(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    real(kind=DBL_KIND),parameter :: zero = 0.d0
    u => mg_get_arrp(mglev, icode)
    u = zero
  end subroutine mg_fill0
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine mg_relax(mglev, ju, jrhs)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, ju, jrhs
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call mg_poisson_relax(mglev, ju, jrhs)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call mg_diff_relax(mglev, ju, jrhs)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call mg_od_relax(mglev, ju, jrhs)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call mg_ad_relax(mglev, ju, jrhs)
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** mg_relax: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine mg_relax
  !-------------------------------------------------------------------------
  ! copy arr
  !-------------------------------------------------------------------------
  subroutine mg_copy(mglev, iout, iin)
    integer,intent(IN) :: mglev, iout, iin
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ain, aout
    aout => mg_get_arrp(mglev,iout)
    ain  => mg_get_arrp(mglev,iin)
    aout = ain
  end subroutine mg_copy
  !-------------------------------------------------------------------------
  ! substract
  !-------------------------------------------------------------------------
  subroutine mg_sub(mglev, iout, ia, ib)
    integer,intent(IN) :: mglev, iout, ia, ib
    integer :: amrlev, gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: out, a, b
    out  => mg_get_arrp(mglev,iout)
    a    => mg_get_arrp(mglev,ia)
    b    => mg_get_arrp(mglev,ib)
    out = a - b
  end subroutine mg_sub
  !-------------------------------------------------------------------------
  ! add
  !-------------------------------------------------------------------------
  subroutine mg_add(mglev, ia, ida)
    integer,intent(IN) :: mglev, ia, ida
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: a, da
    a  => mg_get_arrp(mglev,ia)
    da => mg_get_arrp(mglev,ida)
    a = a + da
  end subroutine mg_add
  !-------------------------------------------------------------------------
  ! add int
  !-------------------------------------------------------------------------
  subroutine mg_addint(mglev, juf, juc, jres)
    integer,intent(IN) :: mglev, juf, juc, jres
    call mg_interp(mglev, jres, juc)
    call mg_add(mglev, juf, jres)
  end subroutine mg_addint
  !-------------------------------------------------------------------------
  ! fill zero in ghost cell for all the grids.
  !-------------------------------------------------------------------------
  subroutine mg_boundary_fill0(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    real(kind=DBL_KIND) :: zero = 0.d0
    integer :: is, ie, js, je, ks, ke

    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    u => mg_get_arrp(mglev,icode)
    if (TouchBoundary(Left, MX)) u(lbound(u,1):is-1,:,:,:) = zero
    if (TouchBoundary(Right,MX)) u(ie+1:ubound(u,1),:,:,:) = zero
    if (TouchBoundary(Left, MY)) u(:,lbound(u,2):js-1,:,:) = zero
    if (TouchBoundary(Right,MY)) u(:,je+1:ubound(u,2),:,:) = zero
    if (TouchBoundary(Left, MZ)) u(:,:,lbound(u,3):ks-1,:) = zero
    if (TouchBoundary(Right,MZ)) u(:,:,ke+1:ubound(u,3),:) = zero

  end subroutine mg_boundary_fill0
  !-------------------------------------------------------------------------
  ! fill zero in ghost cell for all the grids.
  !-------------------------------------------------------------------------
  subroutine mg_boundary_extrap(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    real(kind=DBL_KIND) :: zero = 0.d0
    integer :: is, ie, js, je, ks, ke
    integer :: isg, ieg, jsg, jeg, ksg, keg

    call mg_get_gridsize  (mglev, is,js,ks,ie,je,ke)
    call mg_get_gridsizeGh(mglev, isg,jsg,ksg,ieg,jeg,keg)
    u => mg_get_arrp(mglev,icode)
    if (TouchBoundary(Left, MX)) u(isg:is-1,:,:,:) = -u(is:is-1+Ngh:-1,:,:,:)
    if (TouchBoundary(Right,MX)) u(ie+1:ieg,:,:,:) = -u(ie+1-Ngh:ie:-1,:,:,:)
    if (TouchBoundary(Left, MY)) u(:,jsg:js-1,:,:) = -u(:,js:js-1+Ngh:-1,:,:)
    if (TouchBoundary(Right,MY)) u(:,je+1:jeg,:,:) = -u(:,je+1-Ngh:je:-1,:,:)
    if (TouchBoundary(Left, MZ)) u(:,:,ksg:ks-1,:) = -u(:,:,ks:ks-1+Ngh:-1,:)
    if (TouchBoundary(Right,MZ)) u(:,:,ke+1:keg,:) = -u(:,:,ke+1-Ngh:ke:-1,:)

  end subroutine mg_boundary_extrap
!!$  !-------------------------------------------------------------------------
!!$  ! fill zero in ghost cell for all the grids.
!!$  !-------------------------------------------------------------------------
!!$  subroutine mg_boundary_extrap(mglev, icode)
!!$    integer,intent(IN) :: mglev, icode
!!$    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
!!$    real(kind=DBL_KIND) :: zero = 0.d0
!!$    integer :: is, ie, js, je, ks, ke
!!$
!!$    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
!!$    u => mg_get_arrp(mglev,icode)
!!$    if (TouchBoundary(Left, MX)) u(is-1,:,:,:) = -u(is,:,:,:)
!!$    if (TouchBoundary(Right,MX)) u(ie+1,:,:,:) = -u(ie,:,:,:)
!!$    if (TouchBoundary(Left, MY)) u(:,js-1,:,:) = -u(:,js,:,:)
!!$    if (TouchBoundary(Right,MY)) u(:,je+1,:,:) = -u(:,je,:,:)
!!$    if (TouchBoundary(Left, MZ)) u(:,:,ks-1,:) = -u(:,:,ks,:)
!!$    if (TouchBoundary(Right,MZ)) u(:,:,ke+1,:) = -u(:,:,ke,:)
!!$
!!$  end subroutine mg_boundary_extrap
  !-------------------------------------------------------------------------
  ! check minmax val
  !-------------------------------------------------------------------------
  function mg_get_absmax(mglev, icode) result(umax)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND) :: umax
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    umax = maxval(abs(u(SZ)))
#undef SZ
  end function mg_get_absmax
  !-------------------------------------------------------------------------
  subroutine mg_minmax(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    print *, 'minmax', minval(u(SZ)), maxval(u(SZ))
#undef SZ
  end subroutine mg_minmax
  !-------------------------------------------------------------------------
  subroutine mg_minmaxall(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    u => mg_get_arrp(mglev,icode)
    print *, 'minmaxall', minval(u), maxval(u)
  end subroutine mg_minmaxall
  !-------------------------------------------------------------------------
  subroutine mg_minmaxbnd(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: umin, umax
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    umin = minval( u(lbound(u,1):is-1,:,:,:) )
    umin = min(umin, minval(u(ie+1:ubound(u,1),:,:,:)))
    umin = min(umin, minval(u(:,lbound(u,2):js-1,:,:)))
    umin = min(umin, minval(u(:,je+1:ubound(u,2),:,:)))
    umin = min(umin, minval(u(:,:,lbound(u,3):ks-1,:)))
    umin = min(umin, minval(u(:,:,ke+1:ubound(u,3),:)))

    umax = maxval( u(lbound(u,1):is-1,:,:,:) )
    umax = max(umax, maxval(u(ie+1:ubound(u,1),:,:,:)))
    umax = max(umax, maxval(u(:,lbound(u,2):js-1,:,:)))
    umax = max(umax, maxval(u(:,je+1:ubound(u,2),:,:)))
    umax = max(umax, maxval(u(:,:,lbound(u,3):ks-1,:)))
    umax = max(umax, maxval(u(:,:,ke+1:ubound(u,3),:)))


    print *, 'minmaxbnd', umin, umax
  end subroutine mg_minmaxbnd
  !-------------------------------------------------------------------------
  ! tau correction for nonlinear multigrid
  !-------------------------------------------------------------------------
  subroutine mg_tau(mglev, jtau, ju)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call mg_od_tau(mglev, jtau, ju)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call mg_ad_tau(mglev, jtau, ju)
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** mg_tau: this type of PDE is not supported'
    endif
  end subroutine mg_tau
end module mg
