#include "config.h"
#include "debug_fmg.h"
! Show error log
!#define VERBOSE
!
! Following CPP macros are available.  Set them in config.h
!   FMG_POISSON .............. Poisson equation for selfgravity
!   FMG_OHMIC_DISSIPATION .... Ohmic dissipation for resistive MHD
!   FMG_AMBIPOLAR_DIFFUSION .. Ambipolar diffusion
!   FMG_DIFFUSION ............ Diffusion equation
!-------------------------------------------------------------------------
! Module for FMG cycle
!-------------------------------------------------------------------------
module fmg
  use fmg_data
  implicit none
  private
  real(kind=DBL_KIND),save :: Resh2maxg
  real(kind=DBL_KIND),save,allocatable :: Resmaxg(:), Resmaxl(:), Trerr(:) ! Termination criterion
  integer,save :: N_fmg_cycle

#ifdef FMG_POISSON
  public :: fmg_poisson
#endif !FMG_POISSON

#ifdef FMG_OHMIC_DISSIPATION
  public :: fmg_ohmic_dissipation
#endif !FMG_OHMIC_DISSIPATION

#ifdef FMG_AMBIPOLAR_DIFFUSION
  public :: fmg_ambipolar_diffusion
#endif !FMG_AMBIPOLAR_DIFFUSION

#ifdef FMG_DIFFUSION
  public :: fmg_diffusion
#endif !FMG_DIFFUSION
  

contains

#ifdef FMG_POISSON
#include "fmg_poisson.F90"
#endif !FMG_POISSON

#ifdef FMG_OHMIC_DISSIPATION
#include "fmg_od.F90"
#endif !FMG_OHMIC_DISSIPATION

#ifdef FMG_AMBIPOLAR_DIFFUSION
#include "fmg_ad.F90"
#endif !FMG_AMBIPOLAR_DIFFUSION

#ifdef FMG_DIFFUSION
#include "fmg_diffusion.F90"
#endif !FMG_DIFFUSION

#include "fmg_data_if.F90"
! #include "fmg_test.F90"
  !-------------------------------------------------------------------------
  ! Frontend of multigrid iteration
  !-------------------------------------------------------------------------
#ifdef FMG_POISSON
  subroutine fmg_poisson
    FMG_PDE_TYPE = FMG_PDE_TYPE_POISSON_EQUATION
    call fmg_data_init              ! initialize fmg_data
    call fmg_cycle
    ! -------------------------
    ! test for multigrid engine
    ! -------------------------
!!$    call fmg_test
!!$    call vmg_test
!!$    call mg_test
!!$    call vmg_test_nonlinear
  end subroutine fmg_poisson
#endif !FMG_POISSON
#ifdef FMG_DIFFUSION
  !-------------------------------------------------------------------------
  ! Frontend of multigrid iteration
  !-------------------------------------------------------------------------
  subroutine fmg_diffusion
    FMG_PDE_TYPE = FMG_PDE_TYPE_POISSON_EQUATION
    call fmg_data_init              ! initialize fmg_data
    call fmg_cycle
  end subroutine fmg_diffusion
#endif !FMG_DIFFUSION
#ifdef FMG_OHMIC_DISSIPATION
  !-------------------------------------------------------------------------
  ! Frontend of multigrid iteration
  !-------------------------------------------------------------------------
  subroutine fmg_ohmic_dissipation
    FMG_PDE_TYPE = FMG_PDE_TYPE_OHMIC_DISSIPATION
    call fmg_data_init(MX,MZ)        ! initialize fmg_data
    call fmg_setVector(IETA, .FALSE.)
    N_fmg_cycle = 0
    do
       call fmg_cycle
       if (N_fmg_cycle >= Dt_refineRatio) exit
    enddo
    ! -------------------------
    ! test for multigrid engine
    ! -------------------------
!!$    FMG_PDE_LINEAR = .FALSE.
!!$    call fmg_nonlin
!!$    call fmg_test
!!$    call vmg_test
!!$    call mg_test
!!$    call vmg_test_nonlinear
!!$    call mg_test_nonlinear
  end subroutine fmg_ohmic_dissipation
#endif !FMG_OHMIC_DISSIPATION
#ifdef FMG_AMBIPOLAR_DIFFUSION
  !-------------------------------------------------------------------------
  ! Frontend of multigrid iteration
  !-------------------------------------------------------------------------
  subroutine fmg_ambipolar_diffusion
    FMG_PDE_TYPE = FMG_PDE_TYPE_AMBIP_DIFFUSION
    call fmg_data_init(MX,MZ, NghostCell=2)        ! initialize fmg_data
    call fmg_setVector(IDOD, .FALSE.)
    call fmg_setVector(IDHE, .FALSE.)
    call fmg_setVector(IDAD, .FALSE.)
    FMG_PDE_LINEAR = .FALSE.
    call fmg_nonlin
!!$    call fmg_nonlin_vcycle
    ! -------------------------
    ! test for multigrid engine
    ! -------------------------
!!$    call vmg_test_nonlinear
!!$    call mg_test_nonlinear
  end subroutine fmg_ambipolar_diffusion
#endif !FMG_AMBIPOLAR_DIFFUSION
  !-------------------------------------------------------------------------
  ! initial guess IU と density IRHO が与えられると、IU を更新する。
  !-------------------------------------------------------------------------
  subroutine fmg_cycle
    use mpilib
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_converge
    use fmg_ghostcell
    use vmg
    use mg
    use io_util
    use io
    use string
!!$    real(kind=DBL_KIND),parameter :: errormax = 1.d-7
    real(kind=DBL_KIND),parameter :: errormax = 1.d-3
!!$    real(kind=DBL_KIND),parameter :: errormax = 1.d-15
!!$    integer,parameter :: nfmgmax = 10
    integer :: nfmgmax
    real(kind=DBL_KIND) :: resh2max, ratio_resh2max, ratio_resh2max_prev
    integer :: fmglev, n
    logical :: bool_converge
    logical :: isNotFinite

    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data           ! copy u and rho to LoL

    call vmg_init(FMG_LevelMax)     ! call it every FMG call (depends on AMR level)
    call mg_init(FMG_LevelMax)      ! call it onece at the run
    fmglev = FMG_LevelMin
    call fmg_alloc_arr(fmglev, IPSI)
    call fmg_alloc_arr(fmglev, ISRC)
    call fmg_alloc_f(fmglev)
    call fmg_copy(fmglev, IPSI,  IU)    ! backup u -> psi
    call fmg_copy(fmglev, ISRC, IRHO)   ! backup rho -> src
    call fmg_converge_c2p(fmglev, ISRC)
    call fmg_converge_c2p(fmglev, IPSI)
    call fmg_boundary_physical_grid(IPSI, ISRC)
    call fmg_boundary_u(fmglev,IPSI)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then ! estimate initial error
       resh2max = fmg_get_h2absmax(FMG_LevelMin, ISRC)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       resh2max = fmg_get_absmax(FMG_LevelMin, ISRC)
    else
       resh2max = fmg_get_h2absmax(FMG_LevelMin, ISRC)
    endif
    ratio_resh2max_prev = HUGE(ratio_resh2max_prev)
    nfmgmax = 10
    if (abs(resh2max) <= TINY(resh2max)) nfmgmax = 0
    bool_converge = .true.
    do n = 1, nfmgmax
       call fmg_resid(fmglev, IRHO, IPSI, ISRC)
       call fmg_lin                 ! solve IU by IRHO
       call fmg_add(fmglev, IPSI, IU)

! #ifdef VERBOSE
       ratio_resh2max = Resh2maxg/resh2max
       call print_msg( 'fmg error = ' // trim(num2char(ratio_resh2max)) // ' (PDE_TYPE = ' // trim(num2char(FMG_PDE_TYPE)) // ')' )
! #endif
       if ( ratio_resh2max < errormax ) exit
       if (n == nfmgmax) bool_converge = .false.
       if (ratio_resh2max > ratio_resh2max_prev * 0.5d0 ) bool_converge = .false.
       if ( .not. bool_converge ) exit
       ratio_resh2max_prev = ratio_resh2max
    end do
    if ((.not. bool_converge) .and. (FMG_PDE_TYPE /= FMG_PDE_TYPE_POISSON_EQUATION)) then
       Dt_refineRatio = Dt_refineRatio * 2
       N_fmg_cycle = N_fmg_cycle * 2
       call fmg_set_dtime(fmg_get_dtime()/Dt_refineRatio)
       call print_msg( '** fmg can not converge residual.' )
       call print_msg( '** Dtime = ' // trim(num2char(fmg_get_dtime())) // ' (Dt_refiementRatio = ' // trim(num2char(Dt_refineRatio)) // ')' )
       call mg_finalize
       call vmg_finalize
       call fmg_finalize
       return
    endif
    N_fmg_cycle = N_fmg_cycle + 1
    call fmg_copy(fmglev, IU, IPSI, noskip=.true.) ! restore U (物理境界もコピーする)
    call fmg_converge_c2p(fmglev,IU)
    call fmg_boundary_u(fmglev, IU)
#ifdef FMG_POISSON
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call fmg_psi2g(IU)
       call fmg_restore_u((/MPSI/), IU)
    endif
#endif !FMG_POISSON
#ifdef FMG_DIFFUSION
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       call fmg_restore_u((/MPSI/), IU)
    endif
#endif !FMG_DIFFUSION
#ifdef FMG_OHMIC_DISSIPATION
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
    endif
#endif !FMG_OHMIC_DISSIPATION
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
!!$    call dumpdata
!!$    HALT
  end subroutine fmg_cycle
  !-------------------------------------------------------------------------
  ! initialize
  !-------------------------------------------------------------------------
  subroutine fmg_init
    allocate( Resmaxg(AMR_LevelMin:AMR_LevelMax) )
    allocate( Resmaxl(AMR_LevelMin:AMR_LevelMax) )
    allocate(   Trerr(AMR_LevelMin:AMR_LevelMax) )
  end subroutine fmg_init
  !-------------------------------------------------------------------------
  ! finalize
  !-------------------------------------------------------------------------
  subroutine fmg_finalize
    use fmg_interpol_cubic
    use fmg_converge
    use fmg_ghostcell
    use fmg_reflux
    deallocate(Resmaxg)
    deallocate(Resmaxl)
    deallocate(Trerr)
    call fmg_interp_cubic_finalize
    call fmg_converge_finalize
    call fmg_ghostcell_finalize
    call fmg_reflux_finalize
    call fmg_dealloc_LoL
  end subroutine fmg_finalize
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
  subroutine fmg_lin
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use string
!!$    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
!!$    INTEGER,parameter :: Ncycle=2, Npre=4, Npost=4
    INTEGER :: Ncycle, Npre, Npost
    integer :: fmglev, jcycle, jpre, jpost, vlevel
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
    ! rho を各レベルへremap
    ! --------------------------
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA)
       end if
    enddo
    ! ----------------
    ! 最粗レベルで解く
    ! ----------------
    fmglev = FMG_LevelMax
    call fmg_alloc_arr(fmglev, IU)
    call fmg_alloc_arr(fmglev, IRHS)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_f(fmglev)
    call fmg_fill0(fmglev, IU)  ! initial guess of IU
    call fmg_slvsml(fmglev, IU, IRHO)
#ifdef DEBUG_FMG_OUTPUT_RES
    call fmg_alloc_arr(FMG_LevelMin, IDBG)
#endif !DEBUG_FMG_OUTPUT_RES
    ! ----------------------
    ! 粗レベルから細レベルへ
    ! ----------------------
    do fmglev = FMG_LevelMax-1, FMG_LevelMin, -1
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_f(fmglev)
       call fmg_interp(fmglev, IU, IU, cubic=.TRUE.) ! FMG interpolation
       call fmg_copy(fmglev, IRHS, IRHO)
       do jcycle = 1, Ncycle     ! Ncycle
          do vlevel = fmglev, FMG_LevelMax-1 ! 下り
             do jpre=1,Npre
                call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
#ifdef VERBOSE
                call print_msg( 'error in pre-smooth = ' // trim(num2char(Resh2maxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(fmglev)))
#endif !VERBOSE
             enddo
             call fmg_resid(vlevel, IRES, IU, IRHS, boundary_fill0=.true.)
             call fmg_rstrct(vlevel+1, IRHS, IRES)
             call fmg_fill0(vlevel+1, IU)
             FmgLevel_fill0 = vlevel+1
          enddo
          vlevel = FMG_LevelMax
!!$          do jpre=1,Npre
!!$             call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
!!$#ifdef VERBOSE
!!$             call print_msg( 'error in pre-smooth = ' // trim(num2char(Resh2maxg)) // &
!!$                  ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(fmglev)))
!!$#endif !VERBOSE
!!$          end do
          call fmg_slvsml(vlevel, IU, IRHS)
!!$          do jpost = 1, Npost
!!$             call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
!!$#ifdef VERBOSE
!!$             call print_msg( 'error in post-smooth = ' // trim(num2char(Resh2maxg)) // &
!!$                  ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(fmglev)))
!!$#endif !VERBOSE
!!$          end do
          do vlevel = FMG_LevelMax-1, fmglev, -1 ! 上り
             call fmg_addint(vlevel, IU, IU, IRES)
             do jpost = 1, Npost
                call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
#ifdef VERBOSE
                call print_msg( 'error in post-smooth = ' // trim(num2char(Resh2maxg)) // &
                     ' at level = '// trim(num2char(vlevel)) // ' ' // trim(num2char(fmglev)))
#endif !VERBOSE
             enddo
          enddo
       enddo
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
#ifdef DEBUG_FMG_OUTPUT_RES
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_FMG_OUTPUT_RES
  end subroutine fmg_lin
  !-------------------------------------------------------------------------
  ! Multigrid for nonlinear PDE.
  ! FAS FMG-cycle is adopted
  ! 非線型MGでは残差は使わないので、RESを一時配列として利用する。
  !-------------------------------------------------------------------------
  subroutine fmg_nonlin
    use vmg
    use mg
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use io
    use string
    use mpilib
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    integer :: jcycle, jpre, jpost, vlevel, fmglev
    integer :: icode
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0
!!$    real(kind=DBL_KIND),parameter :: ALPHA = 1.d-5
    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data           ! copy u and rho to LoL
    call vmg_init(FMG_LevelMax)     ! call it every FMG call (depends on AMR level)
    call mg_init(FMG_LevelMax)      ! call it onece at the run
    ! -----------------------------
    ! rho と eta を各レベルへremap
    ! -----------------------------
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       call fmg_alloc_arr(fmglev, IU) !初期値もrestrict
       call fmg_rstrct(fmglev, IU, IU)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA, cubic=.TRUE.)
       end if
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
          do icode = IDOD, IDAD
             call fmg_alloc_arr(fmglev, icode)
             call fmg_rstrct(fmglev, icode, icode)
             call fmg_converge_c2p(fmglev, icode)
             call fmg_boundary_u(fmglev, icode)
             call fmg_ghostcell_fix(fmglev, icode, cubic=.TRUE.)
          end do
       end if
    enddo
#ifdef DEBUG_FMG_OUTPUT_RES
    call fmg_alloc_arr(FMG_LevelMin, IDBG)
#endif !DEBUG_FMG_OUTPUT_RES
#ifdef DEBUG_FMG_OUTPUT_TAU
    call fmg_alloc_arr(FMG_LevelMin, IDBG)
    call fmg_alloc_arr(FMG_LevelMin+1, IDBG)
#endif !DEBUG_FMG_OUTPUT_TAU
    ! ----------------
    ! 最粗レベルで解く
    ! ----------------
    fmglev = FMG_LevelMax
    call fmg_alloc_arr(fmglev, IU)
    call fmg_alloc_arr(fmglev, IRHS)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_f(fmglev)
    call fmg_slvsml(fmglev, IU, IRHO, fmg=.TRUE.)
!!$    do jpre=1,Npost
!!$       call fmg_relax(fmglev, IU, IRHO)
!!$#ifdef VERBOSE
!!$       call fmg_print_prepostlog('post', fmglev, 0)
!!$#endif !VERBOSE
!!$    enddo
    ! ----------------------
    ! 粗レベルから細レベルへ
    ! ----------------------
    do fmglev = FMG_LevelMax-1, FMG_LevelMin, -1
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_arr(fmglev+1, IRUF)
       call fmg_alloc_f(fmglev)
       call fmg_interp(fmglev, IU, IU, cubic=.TRUE.) ! FMG interpolation
       call fmg_copy(fmglev, IRHS, IRHO)
       do jcycle = 1, Ncycle     ! Ncycle
          do vlevel = fmglev, FMG_LevelMax-1 ! 下り
             do jpre=1,Npre
                call fmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call fmg_print_prepostlog(' pre', vlevel, jcycle)
#endif !VERBOSE
             enddo
             call fmg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) ! Uc = R Uf
             call fmg_copy(vlevel+1, IRUF, IU)      ! backup Uc (=R Uf) for GCG
             call fmg_rstrct(vlevel+1, IRHS, IRHS)  ! RHSc = R RHSf
             call fmg_tau(vlevel+1, IRES, IU)       ! tau-correction (RES is tempolary used for tau)
             call fmg_add(vlevel+1, IRHS, IRES)     ! RHS = RHS + tau
             if (vlevel == fmglev) then             ! estimate truncation error
                call fmg_max_forAMRLevel(Trerr, vlevel+1, IRES, absolute=.TRUE., skip=.TRUE.)
                Trerr = ALPHA*Trerr
             end if
          enddo
          vlevel = FMG_LevelMax
!!$          do jpre=1,Npre
!!$             call fmg_relax(vlevel, IU, IRHS)
!!$#ifdef VERBOSE
!!$             call fmg_print_prepostlog(' pre', vlevel, jcycle)
!!$#endif !VERBOSE
!!$          end do
          call fmg_slvsml(vlevel, IU, IRHS)
!!$          do jpost = 1, Npost
!!$             call fmg_relax(vlevel, IU, IRHS)
!!$#ifdef VERBOSE
!!$             call fmg_print_prepostlog('post', vlevel, jcycle)
!!$#endif !VERBOSE
!!$          end do
          do vlevel = FMG_LevelMax-1, fmglev, -1 ! 上り
             call fmg_sub(vlevel+1, IRES, IU, IRUF) ! RESc = Uc - R Uf
             call fmg_interp(vlevel, IRES, IRES, cubic=.TRUE.) ! RESf = P (RESc) = P (Uc - R Uf)
             call fmg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
             do jpost = 1, Npost
                call fmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call fmg_print_prepostlog('post', vlevel, jcycle)
#endif !VERBOSE
             enddo
          enddo
          call fmg_print_errlog(jcycle)
          if (all(Resmaxg < Trerr)) exit
       enddo
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
    ! --------------------
    ! データをAMR-gridに戻す
    ! --------------------
#if defined(FMG_OHMIC_DISSIPATION) || defined(FMG_AMBIPOLAR_DIFFUSION)
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
    endif
#endif !FMG_OHMIC_DISSIPATION or FMG_AMBIPOLAR_DIFFUSION
#ifdef DEBUG_FMG_OUTPUT_RES
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_FMG_OUTPUT_RES
#ifdef DEBUG_FMG_OUTPUT_TAU
    call fmg_interp(FMG_LevelMin, IDBG, IDBG, cubic=.TRUE.)
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_FMG_OUTPUT_TAU
    ! ---------------------
    ! finalize of multigrid
    ! ---------------------
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine fmg_nonlin
  !-------------------------------------------------------------------------
  ! Multigrid for nonlinear PDE.
  ! FAS V-cycle is adopted
  ! 非線型MGでは残差は使わないので、RESを一時配列として利用する。
  !-------------------------------------------------------------------------
  subroutine fmg_nonlin_vcycle
    use vmg
    use mg
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use io
    use string
    use mpilib
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    integer :: jcycle, jpre, jpost, vlevel, fmglev
    integer :: icode
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0
!!$    real(kind=DBL_KIND),parameter :: ALPHA = 1.d-5
    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data           ! copy u and rho to LoL
    call vmg_init(FMG_LevelMax)     ! call it every FMG call (depends on AMR level)
    call mg_init(FMG_LevelMax)      ! call it onece at the run
    ! -----------------------------
    ! rho と eta を各レベルへremap
    ! -----------------------------
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA, tricubic=.TRUE.)
       end if
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
          do icode = IDOD, IDAD
             call fmg_alloc_arr(fmglev, icode)
             call fmg_rstrct(fmglev, icode, icode)
             call fmg_converge_c2p(fmglev, icode)
             call fmg_boundary_u(fmglev, icode)
             call fmg_ghostcell_fix(fmglev, icode, tricubic=.TRUE.)
          end do
       end if
    enddo
    ! --------------------
    ! allocation of arrays
    ! --------------------
    do fmglev = FMG_LevelMin, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRUF)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_f(fmglev)
    end do
#ifdef DEBUG_FMG_OUTPUT_RES
    call fmg_alloc_arr(FMG_LevelMin, IDBG)
#endif !DEBUG_FMG_OUTPUT_RES
#ifdef DEBUG_FMG_OUTPUT_TAU
    call fmg_alloc_arr(FMG_LevelMin, IDBG)
    call fmg_alloc_arr(FMG_LevelMin+1, IDBG)
#endif !DEBUG_FMG_OUTPUT_TAU
    ! ----------------------
    ! V cycle
    ! ----------------------
    call fmg_copy(FMG_LevelMin, IRHS, IRHO) ! set source to rhs
    do jcycle = 1, Ncycle     ! Ncycle
       do vlevel = FMG_LevelMin, FMG_LevelMax-1 ! 下り
          do jpre=1,Npre
             call fmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
          call fmg_print_prepostlog(' pre', vlevel, jcycle)
#endif !VERBOSE
          enddo
          call fmg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) ! Uc = R Uf
          call fmg_copy(vlevel+1, IRUF, IU)      ! backup Uc (=R Uf) for GCG
          call fmg_rstrct(vlevel+1, IRHS, IRHS)  ! RHSc = R RHSf
          call fmg_tau(vlevel+1, IRES, IU)       ! tau-correction (RES is tempolary used for tau)
          call fmg_add(vlevel+1, IRHS, IRES)     ! RHS = RHS + tau
          if (vlevel == FMG_LevelMin) then       ! estimate truncation error
             call fmg_max_forAMRLevel(Trerr, vlevel+1, IRES, absolute=.TRUE., skip=.TRUE.)
             Trerr = ALPHA*Trerr
          end if
       enddo
       vlevel = FMG_LevelMax
!!$       do jpre=1,Npre
!!$          call fmg_relax(vlevel, IU, IRHS)
!!$#ifdef VERBOSE
!!$          call fmg_print_prepostlog(' pre', vlevel, jcycle)
!!$#endif !VERBOSE
!!$       enddo
       call fmg_slvsml(vlevel, IU, IRHS)
!!$       do jpost = 1, Npost
!!$          call fmg_relax(vlevel, IU, IRHS)
!!$#ifdef VERBOSE
!!$          call fmg_print_prepostlog('post', vlevel, jcycle)
!!$#endif !VERBOSE
!!$       enddo
       do vlevel = FMG_LevelMax-1, FMG_LevelMin, -1    ! 上り
          call fmg_sub(vlevel+1, IRES, IU, IRUF) ! RESc = Uc - R Uf
          call fmg_interp(vlevel, IRES, IRES, cubic=.TRUE.) ! RESf = P (RESc) = P (Uc - R Uf)
          call fmg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
          do jpost = 1, Npost
             call fmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
             call fmg_print_prepostlog('post', vlevel, jcycle)
#endif !VERBOSE
          enddo
       enddo
#ifdef VERBOSE
       call fmg_print_errlog(jcycle)
#endif !VERBOSE
       ! termination criterion
       if (all(Resmaxg < Trerr)) exit
       if (jcycle == Ncycle) &
            call print_msg('***** fmg_nonlin was not converged')
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
    ! --------------------
    ! データをAMR-gridに戻す
    ! --------------------
#if defined(FMG_OHMIC_DISSIPATION) || defined(FMG_AMBIPOLAR_DIFFUSION)
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
    endif
#endif !FMG_OHMIC_DISSIPATION or FMG_AMBIPOLAR_DIFFUSION
#ifdef DEBUG_FMG_OUTPUT_RES
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_FMG_OUTPUT_RES
#ifdef DEBUG_FMG_OUTPUT_TAU
    call fmg_interp(FMG_LevelMin, IDBG, IDBG, cubic=.TRUE.)
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_FMG_OUTPUT_TAU
    ! ---------------------
    ! finalize of multigrid
    ! ---------------------
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine fmg_nonlin_vcycle
  !-------------------------------------------------------------------------
  ! solve in the coarest level
  !-------------------------------------------------------------------------
  subroutine fmg_slvsml(fmglev, ju, jrhs, fmg)
    use vmg
    integer,intent(IN) :: fmglev, ju, jrhs
    logical,optional :: fmg
    logical :: bool_fmg
    integer :: n
    bool_fmg = .FALSE. ; if (present(fmg)) bool_fmg = fmg
    if (FMG_PDE_LINEAR) then
       call vmg_fas(ju, jrhs)
    else
       if (bool_fmg) then
          call vmg_fas_fmg(ju, jrhs)
       else
          call vmg_fas(ju, jrhs)
       end if
    end if

!!$    do n = 1, 100
!!$       call fmg_relax(fmglev, ju, jrhs)
!!$    end do
  end subroutine fmg_slvsml
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine fmg_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call fmg_poisson_resid(fmglev, jres, ju, jrhs, boundary_fill0)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call fmg_diff_resid(fmglev, jres, ju, jrhs, boundary_fill0)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call fmg_od_resid(fmglev, jres, ju, jrhs, boundary_fill0)
#endif !FMG_OHMIC_DISSIPATION
    else
       print *, '*** fmg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_resid
  !-------------------------------------------------------------------------
  ! Restrction
  ! fmglevel = fmg level of coarse grid
  ! icode    = identifier of data set
  !-------------------------------------------------------------------------
  subroutine fmg_rstrct(fmglev, iuc, iuf, cubic)
    use fmg_ghostcell
    use restriction
    integer,intent(IN) :: fmglev, iuc, iuf
    logical,intent(IN),optional :: cubic
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, gid, la, if, jf, kf, ic, jc, kc, m, &
         ics, jcs, kcs, ice, jce, kce
    real(kind=DBL_KIND) :: rdv
    logical :: bool_cubic
    bool_cubic = .FALSE.
    if (present(cubic)) bool_cubic = cubic
    lf = fmglev-1
    lc = fmglev
    call fmg_get_gridsize(lc, ics, jcs, kcs, ice, jce, kce)
    if (bool_cubic) call fmg_ghostcell_fix(lf,iuf,tricubic=cubic)
    do la = AMR_LevelMin, AMR_LevelMax
       rdv = fmg_get_dv(la,lf)/fmg_get_dv(la,lc)
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)
          if ( fmg_skip_grid(gid, la, fmglev) ) cycle
          uf => fmg_get_arrp(la,lf,gid,iuf)
          uc => fmg_get_arrp(la,lc,gid,iuc)
          call rstrct( uc, uf, &
               ics, jcs, kcs, ice, jce, kce, &
               ics, jcs, kcs, cubic)
       end do
    end do
!!$    call fmg_boundary_fill0(fmglev, iuc)
  end subroutine fmg_rstrct
  !-------------------------------------------------------------------------
  ! Interpolation
  ! fmglev  = fmg level of finer level
  ! iuf     = icode of fine array (OUT)
  ! iuc     = icode of coarse array (IN)
  ! cubic   = set .true. for cubic interpolation
  !-------------------------------------------------------------------------
  subroutine fmg_interp(fmglev, iuf, iuc, cubic)
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_interpol_cubic
    use io_util
    integer,intent(IN) :: fmglev, iuf, iuc
    logical,intent(IN),optional :: cubic
    logical :: bool_cubic
    integer :: lf, lc
    integer,parameter :: MINBLOCKSIZE = 4

    lf = fmglev
    lc = fmglev+1

    bool_cubic = .FALSE.             ! default
    if (present(cubic)) bool_cubic = cubic

    call fmg_converge_c2p(lc, iuc, cubic) !子と隣接する親の袖
    call fmg_boundary_u(lc, iuc)
    call fmg_ghostcell_fix(lc, iuc, tricubic=bool_cubic)
    call fmg_boundary_extrap(lc, iuc)
    call fmg_boundary_u(lc,iuc)

    if (bool_cubic) then
       call fmg_interp_cubic(fmglev, iuf, iuc)
    else
       call fmg_interp_linear(fmglev, iuf, iuc)
    end if

    call fmg_boundary_fill0(lc, iuc)
    call fmg_boundary_fill0(lf, iuf)
    call fmg_boundary_u(lc,iuc)
    call fmg_boundary_u(lf,iuf)

  end subroutine fmg_interp
  ! ----------------------------------------------------------------
  ! Trilinear interpolation
  ! ----------------------------------------------------------------
  subroutine fmg_interp_linear(fmglev, iuf, iuc)
    use interpolation
    integer,intent(IN) :: fmglev, iuf, iuc
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, gid, la, &
         ifs, jfs, kfs, ife, jfe, kfe
    lf = fmglev
    lc = fmglev+1
    call fmg_get_gridsize(lf, ifs, jfs, kfs, ife, jfe, kfe)
    do la = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)
          if ( fmg_skip_grid(gid, la, lf) ) cycle
          uf => fmg_get_arrp(la,lf,gid,iuf)
          uc => fmg_get_arrp(la,lc,gid,iuc)
          call interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, ifs, jfs, kfs)
       end do
    end do
  end subroutine fmg_interp_linear
  ! ----------------------------------------------------------------
  ! fill zero
  ! ----------------------------------------------------------------
  subroutine fmg_fill0(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    real(kind=DBL_KIND),parameter :: zero = 0.d0
    integer :: amrlev, gid
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          u => fmg_get_arrp(amrlev, fmglev, gid, icode)
          u(:,:,:,:) = zero
       end do
    end do
  end subroutine fmg_fill0
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine fmg_relax(fmglev, ju, jrhs, boundary_fill0)
    use mpilib
    integer,intent(IN) :: fmglev, ju, jrhs
    logical,optional :: boundary_fill0
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call fmg_poisson_relax(fmglev, ju, jrhs, Resh2maxg, boundary_fill0)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call fmg_diff_relax(fmglev, ju, jrhs, Resh2maxg, boundary_fill0)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call fmg_od_relax(fmglev, ju, jrhs, Resh2maxg, boundary_fill0)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call fmg_ad_relax(fmglev, ju, jrhs, Resh2maxg) ! 非線型PDEでは boundary_fill0は不要
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** fmg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_relax
  !-------------------------------------------------------------------------
  ! copy arr
  !-------------------------------------------------------------------------
  subroutine fmg_copy(fmglev, iout, iin, noskip)
    integer,intent(IN) :: fmglev, iout, iin
    logical,optional :: noskip
    logical :: bool_skip
    integer :: amrlev, gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ain, aout
    bool_skip = .TRUE.; if (present(noskip)) bool_skip = .not. noskip
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          aout => fmg_get_arrp(amrlev,fmglev,gid,iout)
          ain  => fmg_get_arrp(amrlev,fmglev,gid,iin)
          aout = ain
       end do
    end do
  end subroutine fmg_copy
  !-------------------------------------------------------------------------
  ! substract
  !-------------------------------------------------------------------------
  subroutine fmg_sub(fmglev, iout, ia, ib)
    integer,intent(IN) :: fmglev, iout, ia, ib
    integer :: amrlev, gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: out, a, b
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          out  => fmg_get_arrp(amrlev,fmglev,gid,iout)
          a    => fmg_get_arrp(amrlev,fmglev,gid,ia)
          b    => fmg_get_arrp(amrlev,fmglev,gid,ib)
          out = a - b
       end do
    end do
  end subroutine fmg_sub
  !-------------------------------------------------------------------------
  ! add
  !-------------------------------------------------------------------------
  subroutine fmg_add(fmglev, ia, ida)
    integer,intent(IN) :: fmglev, ia, ida
    integer :: amrlev, gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: a, da
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          a  => fmg_get_arrp(amrlev,fmglev,gid,ia)
          da => fmg_get_arrp(amrlev,fmglev,gid,ida)
          a = a + da
       end do
    end do
  end subroutine fmg_add
  !-------------------------------------------------------------------------
  ! add int
  !-------------------------------------------------------------------------
  subroutine fmg_addint(fmglev, juf, juc, jres)
    integer,intent(IN) :: fmglev, juf, juc, jres
    call fmg_interp(fmglev, jres, juc)
    call fmg_add(fmglev, juf, jres)
  end subroutine fmg_addint
  !-------------------------------------------------------------------------
  ! tau correction for nonlinear multigrid
  !-------------------------------------------------------------------------
  subroutine fmg_tau(fmglev, jtau, ju)
    integer,intent(IN) :: fmglev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call fmg_od_tau(fmglev, jtau, ju)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call fmg_ad_tau(fmglev, jtau, ju)
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** fmg_tau: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_tau
  !-------------------------------------------------------------------------
  ! shift cell for debug
  !-------------------------------------------------------------------------
  subroutine fmg_shift(fmglev, icodenew, icode)
    integer,intent(IN) :: fmglev, icodenew, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, unew
    integer :: i,j,k,is,js,ks,ie,je,ke,if,jf,kf,amrlev,gid, m
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          u => fmg_get_arrp(amrlev,fmglev,gid,icode)
          unew => fmg_get_arrp(amrlev,fmglev,gid,icodenew)
          do m = lbound(u,4), ubound(u,4)
             do k = ks, ke
                do j = js, je
                   do i = is, ie
                      if = i-1
                      jf = j-1
                      kf = k-1
                      unew(i,j,k,m) = u(if,jf,kf,m)
                   end do
                enddo
             enddo
          enddo
       end do
    end do
  end subroutine fmg_shift
  !-------------------------------------------------------------------------
  ! fill zero in ghost cell for all the grids.
  !-------------------------------------------------------------------------
  subroutine fmg_fill0_ghcell(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, unew
    integer :: amrlev,gid
!!$    real(kind=DBL_KIND) :: foo = 1.d30
    real(kind=DBL_KIND) :: foo = 0.d0

    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          u => fmg_get_arrp(amrlev,fmglev,gid,icode)
          u(lbound(u,1),:,:,:) = foo
          u(ubound(u,1),:,:,:) = foo
          u(:,lbound(u,2),:,:) = foo
          u(:,ubound(u,2),:,:) = foo
          u(:,:,lbound(u,3),:) = foo
          u(:,:,ubound(u,3),:) = foo
       end do
    end do
  end subroutine fmg_fill0_ghcell
  !-------------------------------------------------------------------------
  ! L∞ノルムを得る。 h**2の重みあり。
  !-------------------------------------------------------------------------
  function fmg_get_h2absmax(fmglev, icode) result(normmax)
    use mpilib
    integer,intent(IN) :: icode
    real(kind=DBL_KIND) :: normmax
    real(kind=DBL_KIND) :: h2, normmaxl
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arr
    integer :: amrlev, gid, fmglev
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    normmaxl = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       h2 = fmg_get_h(amrlev, fmglev) **2
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          arr => fmg_get_arrp(amrlev,fmglev,gid,icode)
          normmaxl =  max(normmaxl, maxval(abs(arr(SZ)))*h2)
       end do
    end do
    call mpi_allreduce(normmaxl, normmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
#undef SZ
!!$    print *, normmax
  end function fmg_get_h2absmax
  !-------------------------------------------------------------------------
  ! L∞ノルムを得る。 h**2の重みなし
  !-------------------------------------------------------------------------
  function fmg_get_absmax(fmglev, icode) result(normmax)
    use mpilib
    integer,intent(IN) :: fmglev, icode
    real(kind=DBL_KIND) :: normmax
    real(kind=DBL_KIND) :: normmaxl
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arr
    integer :: amrlev, gid
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    normmaxl = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          arr => fmg_get_arrp(amrlev,fmglev,gid,icode)
          normmaxl =  max(normmaxl, maxval(abs(arr(SZ))))
       end do
    end do
    call mpi_allreduce(normmaxl, normmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
#undef SZ
  end function fmg_get_absmax
  !-------------------------------------------------------------------------
  ! get maximum for each AMR level
  !-------------------------------------------------------------------------
  subroutine fmg_max_forAMRLevel(arrmax, fmglev, icode, absolute, skip)
    use mpilib
    real(kind=DBL_KIND),dimension(AMR_LevelMin:),intent(OUT) :: arrmax
    integer,intent(IN) :: fmglev, icode
    logical,intent(IN),optional :: absolute, skip
    real(kind=DBL_KIND),dimension(ARRAYSIZE1(arrmax)) :: arrmaxlocal
    integer :: amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call fmg_max(arrmaxlocal(amrlev), amrlev, fmglev, icode, absolute, skip, mpireduce=.FALSE.)
    end do
    call mpi_allreduce(arrmaxlocal, arrmax, size(arrmax), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
!!$    print *, '*** truncation error is estimated'
  end subroutine fmg_max_forAMRLevel
  !-------------------------------------------------------------------------
  ! print error log
  !-------------------------------------------------------------------------
  subroutine fmg_print_errlog(jcycle)
    use io_util
    use string
    integer,intent(IN) :: jcycle
    integer :: amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call print_msg( 'error in FMG = ' // trim(num2char(Resmaxg(amrlev))) // ' / ' // trim(num2char(Trerr(amrlev))) &
            // ' at amrlevel = '//trim(num2char(amrlev)) &
            // ', jcycle = '//trim(num2char(jcycle)) &
            )
    end do
  end subroutine fmg_print_errlog
  !-------------------------------------------------------------------------
  ! print error log for pre- and post-smoothing
  !-------------------------------------------------------------------------
  subroutine fmg_print_prepostlog(prepost, vlevel, jcycle)
    use io_util
    use string
    character(len=*),intent(IN) :: prepost
    integer,intent(IN) :: vlevel, jcycle
    integer,dimension(1) :: maxlev
    maxlev = maxloc(Resmaxg)-1+lbound(Resmaxg,1)
    call print_msg( 'FMG error in '//trim(prepost)//'-smooth = '// trim(num2char(maxval(Resmaxg))) &
         // ' max at amrlev = '//trim(num2char(maxlev(1))) &
         // ', vlevel = '//trim(num2char(vlevel)) &
         // ', jcycle = '//trim(num2char(jcycle)) &
         )
  end subroutine fmg_print_prepostlog
end module fmg
