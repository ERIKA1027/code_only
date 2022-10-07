!  test suit for multigrid
#include "config.h"
#include "debug_fmg.h"
#define VERBOSE
  !-------------------------------------------------------------------------
  ! Test for module "mg" for linear PDE.
  !-------------------------------------------------------------------------
  subroutine mg_test
    use mpilib
    use fmg_converge
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_ghostcell
    use mg_data
    use mg
    use io
    integer :: fmglev,mglev, n
    fmglev = 0
    mglev = 0
    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data          ! copy u and rho to LoL
    call fmg_alloc_f(fmglev)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_arr(fmglev, IRHS)

    ! 多重極展開などの境界条件を設定する場合,残差を計算する。
    ! mg_cycleではこのような境界条件を仮定していないため。
    call fmg_fill0(fmglev, IU)
    call fmg_boundary_physical_grid(IU, IRHO)
    call fmg_boundary_u(fmglev,IU)
    call fmg_resid(fmglev, IRES, IU, IRHO)
    call fmg_fill0(fmglev, IU)
    call fmg_converge_c2p(fmglev, IRES)
    call fmg_converge_c2p(fmglev, IRHO)
    call fmg_copy(fmglev, IRHO, IRES)

    call mg_init(fmglev)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_dataPrepare((/IRHO, IU/), (/IRHO, IU/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       call mg_dataPrepare((/IRHO, IU, IETA/), (/IRHO, IU, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call mg_dataPrepare((/IRHO, IU, IETA/), (/IRHO, IU, IETA/))
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    if ( get_myrank() == PRIMARY_RANK) then
!!$       call mg_fill0(mglev, IU)
!!$       call mg_lin
       call mg_fill0(mglev, IU)
       call mg_cycle(IU, IRHO)
    end if
#ifdef IU
    call mg_dataRestore((/IU/), (/IU/))
#endif !IU
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef MPI
       call fmg_restore_u((/MPSI/), IU)
#endif !MPI
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif

    call mg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine mg_test
  !-------------------------------------------------------------------------
  ! Test for module "mg" for nonlinear PDE.
  !-------------------------------------------------------------------------
  subroutine mg_test_nonlinear
    use mpilib
    use mg_data
    use fmg_converge
    use vmg
    use mg
    use io
    integer :: fmglev, n
    fmglev = 0
    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data           ! copy u and rho to LoL
    call fmg_converge_c2p(fmglev, IRHO) !子を持つブロックではrhoが計算されないため。
    call mg_init(fmglev)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_dataPrepare((/IRHO, IU/), (/IRHO, IU/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       call mg_dataPrepare((/IRHO, IU, IETA/), (/IRHO, IU, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call mg_dataPrepare((/IRHO, IU, IETA/), (/IRHO, IU, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call mg_dataPrepare((/IRHO, IU, IDOD, IDHE, IDAD/), (/IRHO, IU, IDOD, IDHE, IDAD/))
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    if (get_myrank() == PRIMARY_RANK) then
       call mg_nonlin(IU, IRHO)
!!$       call mg_nonlin_vcycle(IU, IRHO)
    end if

    call mg_dataRestore((/IU/), (/IU/))

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef MPSI
       call fmg_restore_u((/MPSI/), IU)
#endif !MPSI
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif

!!$    ! debug 最細グリッドの残差を出力する
!!$    call fmg_alloc_arr(fmglev, IETA)
!!$    call mg_dataRestore((/IETA/), (/IETA/))
!!$    call fmg_restore_u((/MVX, MVY, MVZ/), IETA)

    call mg_finalize            ! clear all buffer
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine mg_test_nonlinear
  !-------------------------------------------------------------------------
  ! Test for vmg_fas for *linear* PDE.
  !-------------------------------------------------------------------------
  subroutine vmg_test
    use mpilib
    use fmg_converge
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_ghostcell
    use vmg
    use mg
    use io
    integer :: fmglev, n
    fmglev = 0

    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data          ! copy u and rho to LoL
    call fmg_alloc_f(fmglev)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_arr(fmglev, IRHS)

    ! vmg test

    call fmg_fill0(fmglev, IU)  ! initial guess is zeo
    call fmg_boundary_physical_grid(IU, IRHO)
    call fmg_boundary_u(fmglev,IU)
    call fmg_resid(fmglev, IRES, IU, IRHO)
    call fmg_fill0(fmglev, IU)
    call vmg_init(fmglev)
    call mg_init(fmglev)
    call vmg_fas(IU, IRES)

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef MPSI
       call fmg_restore_u((/MPSI/), IU)
#endif !MPSI
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine vmg_test
  !-------------------------------------------------------------------------
  ! Test for vmg_fas for *nonlinear* PDE.
  !-------------------------------------------------------------------------
  subroutine vmg_test_nonlinear
    use mpilib
    use vmg
    use mg
    use io
    integer :: fmglev, n
    fmglev = 0

    call fmg_init
    call fmg_alloc_LoL          ! make LoL data structure
    call fmg_prepare_data       ! copy u and rho to LoL
    call fmg_alloc_f(fmglev)

    ! vmg test
    call vmg_init(fmglev)
    call mg_init(fmglev)
!!$    call vmg_fas(IU, IRHO)
    call vmg_fas_fmg(IU, IRHO)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef MPSI
       call fmg_restore_u((/MPSI/), IU)
#endif !MPSI
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif
#if defined(DEBUG_VMG_OUTPUT_RES) || defined(DEBUG_VMG_OUTPUT_TAU)
    if (fmg_isVector(IDBG)) then
       call fmg_restore_u((/MVX, MVY, MVZ/), IDBG)
    else
       call fmg_restore_u((/MVX/), IDBG)
    end if
#endif !DEBUG_VMG_OUTPUT_RES
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine vmg_test_nonlinear
  !-------------------------------------------------------------------------
  ! Test for fmg_cycle for linear PDE.
  ! For nonlinear PDE, call fmg_nonlin directly
  !-------------------------------------------------------------------------
  subroutine fmg_test
    use mpilib
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_converge
    use vmg
    use mg
    use io
    real(kind=DBL_KIND),parameter :: errormax = 1.d-5
    integer,parameter :: nfmgmax = 10
    real(kind=DBL_KIND) :: resh2max
    integer :: fmglev, n
    call fmg_init
    call fmg_alloc_LoL              ! make LoL data structure
    call fmg_prepare_data          ! copy u and rho to LoL
    call vmg_init(FMG_LevelMax)     ! call it every FMG call (depends on AMR level)
    call mg_init(FMG_LevelMax)      ! call it onece at the run
    fmglev = FMG_LevelMin
    call fmg_alloc_arr(fmglev, IPSI)
    call fmg_alloc_arr(fmglev, ISRC)
    call fmg_alloc_f(fmglev)

    call fmg_fill0(fmglev, IU)        ! initial guess
    call fmg_copy(fmglev, IPSI,  IU)    ! backup u -> psi
    call fmg_copy(fmglev, ISRC, IRHO)   ! backup rho -> src
    call fmg_converge_c2p(fmglev, ISRC)
    call fmg_converge_c2p(fmglev, IPSI)
    call fmg_boundary_physical_grid(IPSI, ISRC)
    call fmg_boundary_u(fmglev,IPSI)

    do n = 1, nfmgmax
       call fmg_resid(fmglev, IRHO, IPSI, ISRC)
       if (n == 1) resh2max = fmg_get_h2absmax(fmglev, IRHO)
       call fmg_lin                 ! solve IU by IRHO
       call fmg_add(fmglev, IPSI, IU)
#ifdef VERBOSE
       if (myrank == 0) print *, 'fmg error', Resh2maxg/resh2max, '*'
#endif !VERBOSE
       if ( Resh2maxg/resh2max < errormax ) exit
    end do

    call fmg_copy(fmglev, IU, IPSI)       ! restore U
!!$    call fmg_copy(fmglev, IU, IRHO)       ! debug moniter residual
    call fmg_converge_c2p(fmglev,IU)

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef MPSI
       call fmg_restore_u((/MPSI/), IU)
#endif !MPSI
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#if defined(MBX) && defined(MBY) && defined(MBZ)
       call fmg_restore_u((/MBX, MBY, MBZ/), IU)
#endif ! MBX MBY MBZ
    else
       print *, '*** mg_test: this type of PDE is not supported', FMG_PDE_TYPE
    endif

    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    HALT
  end subroutine fmg_test
