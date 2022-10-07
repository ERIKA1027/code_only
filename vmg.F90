#include "config.h"
#include "debug_fmg.h"
! Show error log
! #define VERBOSE
! ポテンシャルの境界値を外挿してから、子グリッドへ内装するか？
! U を解く場合に、SOLVE_U を定義する。
! 残差を解く場合には、未定義とする。境界値はゼロなので。
!-------------------------------------------------------------------------
! Module for V cycle of FAS
!-------------------------------------------------------------------------
module vmg
  use mpilib
  use fmg_data
  use vmg_interpol
  use fmg_converge
  implicit none
  private
  integer,save :: FMG_Level
  integer,save :: Imin, Jmin, Kmin, Imax, Jmax, Kmax
  integer,save :: Imingh, Jmingh, Kmingh, Imaxgh, Jmaxgh, Kmaxgh
  integer,save :: Icmin, Jcmin, Kcmin, Icmax, Jcmax, Kcmax
  real(kind=DBL_KIND),save,allocatable :: Resmaxg(:), Trerr(:)
  logical,save :: CubicInterp
  logical,save :: UseSerialMg
  public :: vmg_fas, vmg_fas_fmg, &
       vmg_init, vmg_finalize, vmg_rstrct, vmg_interp, vmg_relax
contains

#ifdef FMG_POISSON
#include "vmg_poisson.F90"
#endif !FMG_POISSON

#ifdef FMG_OHMIC_DISSIPATION
#include "vmg_od.F90"
#endif !FMG_OHMIC_DISSIPATION

#ifdef FMG_DIFFUSION
#include "vmg_diffusion.F90"
#endif !FMG_DIFFUSION

#ifdef FMG_AMBIPOLAR_DIFFUSION
#include "vmg_ad.F90"
#endif !FMG_AMBIPOLAR_DIFFUSION

  !-------------------------------------------------------------------------
  ! initialize (call every FMG. VMG_init depends on AMR level)
  !-------------------------------------------------------------------------
  subroutine vmg_init(fmglev)
    use util
    integer,intent(IN) :: fmglev
    FMG_Level = fmglev
    Imin = GridSize(FMG_Level)%Imin
    Jmin = GridSize(FMG_Level)%Jmin
    kmin = GridSize(FMG_Level)%Kmin
    Imax = GridSize(FMG_Level)%Imax
    Jmax = GridSize(FMG_Level)%Jmax
    Kmax = GridSize(FMG_Level)%Kmax

    Imingh = fmg_get_imingh(FMG_Level)
    Jmingh = fmg_get_jmingh(FMG_Level)
    Kmingh = fmg_get_kmingh(FMG_Level)
    Imaxgh = fmg_get_imaxgh(FMG_Level)
    Jmaxgh = fmg_get_jmaxgh(FMG_Level)
    Kmaxgh = fmg_get_kmaxgh(FMG_Level)

    Icmin = Imin
    Jcmin = Jmin
    Kcmin = Kmin
    Icmax = IJKC(Imax,Imin)
    Jcmax = IJKC(Jmax,Jmin)
    Kcmax = IJKC(Kmax,Kmin)

    if (Ngh >= 2) then
       CubicInterp = .TRUE.
    else
       CubicInterp = .FALSE.
    end if

    allocate( Resmaxg(AMR_LevelMin:AMR_LevelMax) )
    allocate( Trerr(AMR_LevelMin+1:AMR_LevelMax) )

    call vmg_interp_init(FMG_Level)
    call fmg_converge_init(FMG_Level)

    if ( util_isPowerOfTow(NGI_BASE) .and. &
         util_isPowerOfTow(NGJ_BASE) .and. &
         util_isPowerOfTow(NGK_BASE) ) then
       UseSerialMg = .TRUE.
    else
       UseSerialMg = .FALSE.
    end if

  end subroutine vmg_init
  !-------------------------------------------------------------------------
  ! finalize
  !-------------------------------------------------------------------------
  subroutine vmg_finalize
    deallocate(Resmaxg)
    deallocate(Trerr)
    call vmg_interp_finalize
    call fmg_converge_finalize
  end subroutine vmg_finalize
  !-------------------------------------------------------------------------
  ! solve Poisson eq. of u, given u (initial guess) and rho (source)
  ! This subroutine is given by `Mutigrid Methods for Boundary Value
  ! Problems. I. W. H.Press & S.A.Teukolsky, Computer in Physics, 5
  !  514-519(1991)
  !
  ! PARAMETER:
  !   NG = the number of grids hierarchy.
  !   Npre = the number of pre-smoothing iteration.
  !   Npre = the number of post-smoothing iteration.
  !-------------------------------------------------------------------------
! vmg.F90 の IRUF と競合するので、IRUFの代りにIPSIを使う。FASではIPSIは未使用なので。
#define IRUF IPSI
  subroutine vmg_fas(ju, jrho)
    use fmg_boundary
    use fmg_boundary_phys
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    integer :: vmglev, jcycle, jpre, jpost, vlevel
!!$    INTEGER,parameter :: Ncycle=1, Npre=2, Npost=2
    INTEGER :: Ncycle, Npre, Npost
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       Ncycle=1
       Npre=2
       Npost=2
    else
       Ncycle=1
       Npre=2
       Npost=2
    endif

#ifdef DEBUG_VMG_OUTPUT_RES
    call fmg_alloc_arr(FMG_Level, IDBG) ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES

    call fmg_alloc_arr(FMG_Level, IRHO)
    call fmg_alloc_arr(FMG_Level, IRHS)
    call fmg_alloc_arr(FMG_Level, IRES)
    call fmg_alloc_arr(FMG_Level, IU)
    call fmg_alloc_arr(FMG_Level, IRUF)
    call fmg_alloc_f(FMG_Level)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_alloc_arr(FMG_Level, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call fmg_alloc_arr(FMG_Level, IDOD)
       call fmg_alloc_arr(FMG_Level, IDHE)
       call fmg_alloc_arr(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (jrho /= IRHO) then
       do vmglev = AMR_LevelMin, AMR_LevelMax
          call vmg_copy(vmglev, IRHO, jrho)
       enddo
    endif
    if (ju /= IU) then
       print *, '**** error in mg_dataIF'
       call mpi_finalize(ierr)
       stop
    end if
    call fmg_converge_c2p(FMG_Level, IU)
    call fmg_converge_c2p(FMG_Level, IRHO)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_converge_c2p(FMG_Level, IETA)
       call fmg_boundary_u(FMG_Level, IETA)
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_converge_c2p(FMG_Level, IDOD)
       call fmg_boundary_u(FMG_Level, IDOD)
       call fmg_converge_c2p(FMG_Level, IDHE)
       call fmg_boundary_u(FMG_Level, IDHE)
       call fmg_converge_c2p(FMG_Level, IDAD)
       call fmg_boundary_u(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    call fmg_boundary_fill0(FMG_Level, IU)
    call fmg_boundary_u(FMG_Level, IU)

    ! ----------------------
    ! V cycle
    ! ----------------------
    do vmglev = AMR_LevelMin, AMR_LevelMax ! set source to rhs
       call vmg_copy(vmglev, IRHS, IRHO)
    enddo
    do jcycle = 1, Ncycle     ! Ncycle
       do vlevel = AMR_LevelMax, AMR_LevelMin+1, -1 ! 下り
          call vmg_gridboundary(vlevel, IU)
          do jpre = 1, Npre
             call vmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
             call print_msg( 'VMG error in pre-smooth =  ' // trim(num2char(Resmaxg(vlevel))) // &
                  ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
          enddo

          call vmg_rstrct(vlevel-1, IU, cubic=.TRUE.) ! restrict U to parent grid
          call vmg_copy(vlevel-1, IRUF, IU) ! backup restricted U for CGC
          call vmg_rstrct(vlevel-1, IRHS)        ! restrict RHS to parent grid
          call vmg_tau(vlevel-1, IRES, IU)       ! tau-correction (RES is tempolary used for tau)
          call vmg_add_onov(vlevel-1, IRHS, IRES)  ! RHS = RHS + tau (for overlaped blocks)
          Trerr(vlevel) = ALPHA* vmg_get_max(vlevel-1, IRES, absolute=.TRUE., xskip=.TRUE.) ! estimate truncation error
       enddo
#ifdef ZERO_AVERAGE_DENSITY
       if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) &
            call vmg_zeroAverage(AMR_LevelMin, IRHS)
#endif !ZERO_AVERAGE_DENSITY

       call vmg_slvsml(IU, IRHS)

       do vlevel = AMR_LevelMin+1, AMR_LevelMax  ! 上り
          call vmg_sub(vlevel-1, IRES, IU, IRUF) ! RESc = Uc - R Uf
          call vmg_interp(vlevel, IRES, IRES)    ! RESf = P (RESc) = P (Uc - R Uf)
          call vmg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
          call vmg_gridboundary(vlevel, IU)
          do jpost = 1, Npost
             call vmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
             call print_msg( 'VMG error in post-smooth = ' // trim(num2char(Resmaxg(vlevel))) // &
                  ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
          enddo
#ifdef VERBOSE
          call print_msg( 'error in VMG = ' // trim(num2char(Resmaxg(vlevel))) // ' / ' // trim(num2char(Trerr(vlevel))) &
               // ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
       enddo
       ! termination criterion
       if (all(Resmaxg(AMR_LevelMin+1:AMR_LevelMax) < Trerr(AMR_LevelMin+1:AMR_LevelMax))) exit
!!$       if( vmg_get_error(IRHS) < errormax ) return
    enddo
#undef IRUF
  end subroutine vmg_fas
  !-------------------------------------------------------------------------
  ! solve Poisson eq. of u, given u (initial guess) and rho (source)
  ! This subroutine is given by `Mutigrid Methods for Boundary Value
  ! Problems. I. W. H.Press & S.A.Teukolsky, Computer in Physics, 5
  !  514-519(1991)
  !
  ! PARAMETER:
  !   NG = the number of grids hierarchy.
  !   Npre = the number of pre-smoothing iteration.
  !   Npre = the number of post-smoothing iteration.
  !-------------------------------------------------------------------------
! vmg.F90 の IRUF と競合するので、IRUFの代りにIPSIを使う。FASではIPSIは未使用なので。
#define IRUF IPSI
  subroutine vmg_fas_fmg(ju, jrho)
    use fmg_boundary
    use fmg_boundary_phys
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    integer :: vmglev, jcycle, jpre, jpost, vlevel
    ! original
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    real(kind=DBL_KIND),parameter :: ALPHA = 1.d0/3.d0
#if defined(DEBUG_VMG_OUTPUT_RES) || defined(DEBUG_VMG_OUTPUT_TAU)
    call fmg_alloc_arr(FMG_Level, IDBG) ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
    call fmg_alloc_arr(FMG_Level, IRHO)
    call fmg_alloc_arr(FMG_Level, IRHS)
    call fmg_alloc_arr(FMG_Level, IRES)
    call fmg_alloc_arr(FMG_Level, IU)
    call fmg_alloc_arr(FMG_Level, IRUF)
    call fmg_alloc_f(FMG_Level)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_alloc_arr(FMG_Level, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call fmg_alloc_arr(FMG_Level, IDOD)
       call fmg_alloc_arr(FMG_Level, IDHE)
       call fmg_alloc_arr(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (jrho /= IRHO) then
       do vmglev = AMR_LevelMin, AMR_LevelMax
          call vmg_copy(vmglev, IRHO, jrho)
       enddo
    endif
    if (ju /= IU) then
       print *, '**** error in mg_dataIF'
       call mpi_finalize(ierr)
       stop
    end if
    call fmg_converge_c2p(FMG_Level, IU)
    call fmg_converge_c2p(FMG_Level, IRHO)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       ! nothing to do
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_converge_c2p(FMG_Level, IETA)
       call fmg_boundary_u(FMG_Level, IETA)
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_converge_c2p(FMG_Level, IDOD)
       call fmg_boundary_u(FMG_Level, IDOD)
       call fmg_converge_c2p(FMG_Level, IDHE)
       call fmg_boundary_u(FMG_Level, IDHE)
       call fmg_converge_c2p(FMG_Level, IDAD)
       call fmg_boundary_u(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    call fmg_boundary_fill0(FMG_Level, IU)
    call fmg_boundary_u(FMG_Level, IU)

    ! ----------------------
    ! set source to rhs
    ! ----------------------
    do vmglev = AMR_LevelMin, AMR_LevelMax
       call vmg_copy(vmglev, IRHS, IRHO)
    enddo
    ! ----------------
    ! 最粗レベルで解く
    ! ----------------
    call vmg_slvsml(IU, IRHO,fmg=.TRUE.)
    ! ----------------------
    ! 粗レベルから細レベルへ
    ! ----------------------
    do vmglev = AMR_LevelMin+1, AMR_LevelMax
       call vmg_interp(vmglev, IU, IU)
       do jcycle = 1, Ncycle     ! Ncycle
          do vlevel = vmglev, AMR_LevelMin+1, -1 ! 下り
             call vmg_gridboundary(vlevel, IU)
             do jpre = 1, Npre
                call vmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'VMG error in pre-smooth =  ' // trim(num2char(Resmaxg(vlevel))) // &
                     ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
             enddo

             call vmg_rstrct(vlevel-1, IU, cubic=.TRUE.) ! restrict U to parent grid
             call vmg_copy(vlevel-1, IRUF, IU) ! backup restricted U for CGC
             call vmg_rstrct(vlevel-1, IRHS)        ! restrict RHS to parent grid
             call vmg_tau(vlevel-1, IRES, IU)       ! tau-correction (RES is tempolary used for tau)
             call vmg_add_onov(vlevel-1, IRHS, IRES)  ! RHS = RHS + tau (for overlaped blocks)
             if (vlevel == vmglev) &
                  Trerr(vlevel) = ALPHA* vmg_get_max(vlevel-1, IRES, absolute=.TRUE., xskip=.TRUE.) ! estimate truncation error
          enddo
#ifdef ZERO_AVERAGE_DENSITY
          if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) &
               call vmg_zeroAverage(AMR_LevelMin, IRHS)
#endif !ZERO_AVERAGE_DENSITY

          call vmg_slvsml(IU, IRHS)

          do vlevel = AMR_LevelMin+1, vmglev  ! 上り
             call vmg_sub(vlevel-1, IRES, IU, IRUF) ! RESc = Uc - R Uf
             call vmg_interp(vlevel, IRES, IRES)    ! RESf = P (RESc) = P (Uc - R Uf)
             call vmg_add(vlevel, IU, IRES)         ! CGC: U = U + RESf = U + P (Uc - R Uf)
             call vmg_gridboundary(vlevel, IU)
             do jpost = 1, Npost
                call vmg_relax(vlevel, IU, IRHS)
#ifdef VERBOSE
                call print_msg( 'VMG error in post-smooth = ' // trim(num2char(Resmaxg(vlevel))) // &
                  ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
#endif !VERBOSE
             enddo
! #ifdef VERBOSE
             call print_msg( 'error in VMG = ' // trim(num2char(Resmaxg(vlevel))) // ' / ' // trim(num2char(Trerr(vlevel))) &
                  // ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
! #endif !VERBOSE
          enddo
          ! termination criterion
          if (all(Resmaxg(AMR_LevelMin+1:AMR_LevelMax) < Trerr(AMR_LevelMin+1:AMR_LevelMax))) exit
!!$       if( vmg_get_error(IRHS) < errormax ) return
       enddo
    end do
#undef IRUF
  end subroutine vmg_fas_fmg
  !-------------------------------------------------------------------------
  ! solve in the coarest level
  !-------------------------------------------------------------------------
  subroutine vmg_slvsml(ju, jrhs, fmg)
    use mg_data
    use mg
    use fmg_boundary_phys
    integer,intent(IN) :: ju, jrhs
    logical,optional :: fmg
    logical :: bool_fmg
    integer :: n
    integer,parameter :: N_VMGRELAX = 10

    if (.not. UseSerialMg) then
       do n = 1, N_VMGRELAX
          call vmg_relax( AMR_LevelMin, ju, jrhs)
       end do
       return
    endif

    bool_fmg = .FALSE. ; if (present(fmg)) bool_fmg = fmg

    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_dataPrepare((/ISRC, IPSI/), (/jrhs, ju/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       call mg_dataPrepare((/ISRC, IPSI, IETA/), (/jrhs, ju, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call mg_dataPrepare((/ISRC, IPSI, IETA/), (/jrhs, ju, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call mg_dataPrepare((/ISRC, IPSI, IDOD, IDHE, IDAD/), (/jrhs, ju, IDOD, IDHE, IDAD/))
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    if (get_myrank() == PRIMARY_RANK) then
       if (FMG_PDE_LINEAR) then
          call mg_cycle(IPSI, ISRC)
       else
          if (bool_fmg) then
             call mg_nonlin(IPSI, ISRC)
          else
             call mg_nonlin_vcycle(IPSI, ISRC)
          end if
       endif
    endif
    call mg_dataRestore((/ju/), (/IPSI/))
#if defined(DEBUG_VMG_OUTPUT_RES)
    call mg_dataRestore((/IDBG/), (/IDBG/))
#endif !DEBUG_VMG_OUTPUT_RES
  end subroutine vmg_slvsml
  !-------------------------------------------------------------------------
  ! tau correction (tau = L R uf - R L uf)
  !-------------------------------------------------------------------------
  subroutine vmg_tau(amrlev, jtau, ju)
    integer,intent(IN) :: amrlev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call vmg_poisson_tau(amrlev, jtau, ju)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call vmg_diff_tau(amrlev, jtau, ju)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call vmg_od_tau(amrlev, jtau, ju)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call vmg_ad_tau(amrlev, jtau, ju)
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** vmg_tau: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (.not. CubicInterp) call vmg_tau_crop(amrlev, jtau) ! 境界を３次補間すると不要
  end subroutine vmg_tau
  !-------------------------------------------------------------------------
  ! crop tau into interrior of grid-level.
  ! Fill zero into interface point. See 362 of Multigrid text book Trottengberg et al.
  !-------------------------------------------------------------------------
  subroutine vmg_tau_crop(amrlev, jtau)
    use mpilib
    integer,intent(IN) :: amrlev, jtau
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau
    integer :: ndir, lr, cgid, lc, lp, crank
    real(kind=DBL_KIND),parameter :: zero = 0.D0
    type(t_ablock),pointer :: cblock ! child block
    myrank = get_myrank()
    lp = amrlev                 ! parent AMR level
    lc = lp + 1                 ! child AMR level
    do crank = 0, (NPE)-1
       do cgid = lbound(Ancestry(lc,crank)%Block,1), ubound(Ancestry(lc,crank)%Block,1) !子供から親を探す
          cblock => Ancestry(lc,crank)%Block(cgid)

          if ( cblock%NeighborParentLevel(Left, MX) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(Imin,:,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Left, MY) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,Jmin,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Left, MZ) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,:,Kmin,:) = zero
          end if

          if ( cblock%NeighborParentLevel(Right, MX) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(Imax,:,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Right, MY) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,Jmax,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Right, MZ) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,:,Kmax,:) = zero
          end if

       end do
    end do
  end subroutine vmg_tau_crop
  !-------------------------------------------------------------------------
  ! substract iout = ia - ib on blocks overlaped with child blocks
  !-------------------------------------------------------------------------
  subroutine vmg_sub_onov(amrlev, iout, ia, ib)
    integer,intent(IN) :: amrlev, iout, ia, ib
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: out, a, b
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       out  => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       a    => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       b    => fmg_get_arrp(amrlev,FMG_Level,gid,ib)
       out = a - b
    end do
  end subroutine vmg_sub_onov
  !-------------------------------------------------------------------------
  ! substract iout = ia - ib on blocks overlaped with child blocks
  !-------------------------------------------------------------------------
  subroutine vmg_sub(amrlev, iout, ia, ib)
    integer,intent(IN) :: amrlev, iout, ia, ib
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: out, a, b
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       out  => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       a    => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       b    => fmg_get_arrp(amrlev,FMG_Level,gid,ib)
       out = a - b
    end do
  end subroutine vmg_sub
  !-------------------------------------------------------------------------
  ! add jtau to jrhs on blocks overlaped with child blocks
  !-------------------------------------------------------------------------
  subroutine vmg_add_onov(amrlev, ia, ida)
    integer,intent(IN) :: amrlev, ia, ida
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: a, da
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       a  => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       da => fmg_get_arrp(amrlev,FMG_Level,gid,ida)
       a = a + da
    end do
  end subroutine vmg_add_onov
  !-------------------------------------------------------------------------
  ! add
  !-------------------------------------------------------------------------
  subroutine vmg_add(amrlev, ia, ida)
    integer,intent(IN) :: amrlev, ia, ida
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: a, da
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       a  => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       da => fmg_get_arrp(amrlev,FMG_Level,gid,ida)
       a = a + da
    end do
  end subroutine vmg_add
  !-------------------------------------------------------------------------
  ! Restrction
  ! amrlev = vmg level of coarse grid
  ! icode  = identifier of data set
  !-------------------------------------------------------------------------
  subroutine vmg_rstrct(amrlev, icode, cubic) ! debug EXP
    use fmg_ghostcell
    integer,intent(IN) :: amrlev, icode
    logical,intent(IN),optional :: cubic
    integer :: ndir
    if (present(cubic)) then
       if (cubic)  then         !cubicの場合には,細かいグリッドの境界条件を整える
          call fmg_ghfix_parentlev_init(FMG_Level)
          call fmg_ghfix_samelev_init(FMG_Level)
          do ndir = MX, MZ
             call fmg_ghfix_parentlev(amrlev+1, FMG_Level, ndir, icode, tricubic=cubic)
             call fmg_ghfix_samelev(amrlev+1, FMG_Level, ndir,icode)
          end do
       end if
    end if
    call fmg_converge_c2p_lev(amrlev,FMG_Level,icode, cubic)
  end subroutine vmg_rstrct
  !-------------------------------------------------------------------------
  ! Interpolation
  ! amrlev = vmg level of original data (fine level)
  ! juf    = output data (fine)
  ! juc    = input data (coarse)
  !-------------------------------------------------------------------------
  subroutine vmg_interp(amrlev, juf, juc)
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, juf, juc
    call vmg_ghostcell(amrlev-1, juc)
    call vmg_boundary_extrap(amrlev-1, juc) !１次的に袖を外挿する。
    call vmg_boundary_u(amrlev-1, FMG_Level, juc) ! 折り返しなどの境界条件
    call vmg_interp_p2c_lev(amrlev, juf, juc, cubic=CubicInterp)
    call vmg_boundary_fill0(amrlev-1, juc) !外装した袖を元に戻す。
    call vmg_boundary_u(amrlev-1, FMG_Level, juc)      ! 折り返しなどの境界条件
  end subroutine vmg_interp
  ! ----------------------------------------------------------------
  ! interp and CGC (Coarse Grid Correction) in FAS
  ! U_new = U_guess + P ( U_H - R U_h ) = U_guess + P U_H - P RUF
  ! ----------------------------------------------------------------
#define SZ Imin:Imax,Jmin:Jmax,Kmin:Kmax,:
  subroutine vmg_cgc(amrlev, ju, jruf, jtmp)
    integer,intent(IN) :: amrlev, ju, jruf, jtmp
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: du, uf
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tmp, ruf, uc
!!$    ! 親グリッドで
!!$    do gid = fmg_get_gidmin(amrlev-1), fmg_get_gidmax(amrlev-1)
!!$       if (Ancestry(amrlev-1,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
!!$       tmp  => fmg_get_arrp(amrlev-1, FMG_Level, gid, jtmp)
!!$       uc   => fmg_get_arrp(amrlev-1, FMG_Level, gid, ju)
!!$       ruf  => fmg_get_arrp(amrlev-1, FMG_Level, gid, jruf)
!!$       tmp = uc - ruf
!!$    enddo
!!$    call vmg_interp(amrlev, jtmp, jtmp) ! tmpの境界は決まっていない
!!$    ! 子グリッドで
!!$    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
!!$       du  => fmg_get_arrp(amrlev, FMG_Level, gid, jtmp)
!!$       uf  => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
!!$       uf = uf + du
!!$    enddo
    call vmg_interp(amrlev, jtmp, ju)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       du  => fmg_get_arrp(amrlev, FMG_Level, gid, jtmp)
       uf  => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
       uf(SZ) = uf(SZ) + du(SZ) ! correction (+ P U_H)
    enddo
    call vmg_interp(amrlev, jtmp, jruf)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       du  => fmg_get_arrp(amrlev, FMG_Level, gid, jtmp)
       uf  => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
       uf(SZ) = uf(SZ) - du(SZ) ! correction (- P RUF)
    enddo
  end subroutine vmg_cgc
#undef SZ
  ! ----------------------------------------------------------------
  ! fill0
  ! ----------------------------------------------------------------
  subroutine vmg_fill0(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       u  => fmg_get_arrp(amrlev, FMG_Level, gid, icode)
       u = 0.d0
    enddo
  end subroutine vmg_fill0
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine vmg_relax(amrlev, ju, jrhs)
    use mpilib
    integer,intent(IN) :: amrlev, ju, jrhs
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
#ifdef FMG_POISSON
       call vmg_poisson_relax(amrlev, ju, jrhs)
#endif !FMG_POISSON
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
#ifdef FMG_DIFFUSION
       call vmg_diff_relax(amrlev, ju, jrhs)
#endif !FMG_DIFFUSION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
#ifdef FMG_OHMIC_DISSIPATION
       call vmg_od_relax(amrlev, ju, jrhs)
#endif !FMG_OHMIC_DISSIPATION
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
#ifdef FMG_AMBIPOLAR_DIFFUSION
       call vmg_ad_relax(amrlev, ju, jrhs)
#endif !FMG_AMBIPOLAR_DIFFUSION
    else
       print *, '*** vmg_relax: this type of PDE is not supported', FMG_PDE_TYPE
    end if
  end subroutine vmg_relax
  !-------------------------------------------------------------------------
  ! copy arr
  !-------------------------------------------------------------------------
  subroutine vmg_copy(amrlev, iout, iin)
    integer,intent(IN) :: amrlev, iout, iin
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ain, aout
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       aout => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       ain  => fmg_get_arrp(amrlev,FMG_Level,gid,iin)
       aout = ain
    end do
  end subroutine vmg_copy
  !-------------------------------------------------------------------------
  ! copy arr
  !-------------------------------------------------------------------------
  subroutine vmg_copy_onov(amrlev, iout, iin)
    integer,intent(IN) :: amrlev, iout, iin
    integer :: gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ain, aout
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       aout => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       ain  => fmg_get_arrp(amrlev,FMG_Level,gid,iin)
       aout = ain
    end do
  end subroutine vmg_copy_onov
  !-------------------------------------------------------------------------
  ! fix the ghost cell from parent grid level
  !-------------------------------------------------------------------------
  subroutine vmg_gridboundary(amrlev, icode)
    use fmg_ghostcell
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, icode
    integer :: ndir
!!$    call vmg_ghostcell(amrlev-1, icode) ! additional6/2
    call vmg_boundary_u(amrlev-1, FMG_Level, icode)
    call fmg_ghfix_parentlev_init(FMG_Level)
!!$    call fmg_ghfix_samelev_init(FMG_Level) ! debug EXP
    do ndir = MX, MZ
!!$       call fmg_ghfix_samelev(amrlev-1, FMG_Level, ndir,icode) ! debug EXP
       call fmg_ghfix_parentlev(amrlev, FMG_Level, ndir, icode, cubic=CubicInterp)
    enddo
    call vmg_boundary_fill0(amrlev, icode)
    call vmg_boundary_u(amrlev, FMG_Level, icode)
  end subroutine vmg_gridboundary
  !-------------------------------------------------------------------------
  ! fix the ghost cell between the same grid levels
  !-------------------------------------------------------------------------
  subroutine vmg_ghostcell(amrlev, icode)
    use fmg_ghostcell
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, icode
    integer :: ndir
    call vmg_boundary_u(amrlev, FMG_Level, icode)
    call fmg_ghfix_samelev_init(FMG_Level)
!!$    call fmg_ghfix_parentlev_init(FMG_Level) ! debug EXP
    do ndir = MX, MZ
!!$       call fmg_ghfix_parentlev(amrlev, FMG_Level, ndir, icode, cubic=CubicInterp) ! debug EXP
       call fmg_ghfix_samelev(amrlev, FMG_Level, ndir,icode)
    enddo
    call vmg_boundary_fill0(amrlev, icode)
    call vmg_boundary_u(amrlev, FMG_Level, icode)
  end subroutine vmg_ghostcell
  !-------------------------------------------------------------------------
  ! set zero boundary condition for U, RES, ....
  !-------------------------------------------------------------------------
  subroutine vmg_boundary_fill0(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    real(kind=DBL_KIND),parameter :: zero = 0.d0
    integer :: id, fmglev
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: up

    fmglev = FMG_Level

    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MX) ) then
          up(Imingh:Imin-1,:,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MX) ) then
          up(Imax+1:Imaxgh,:,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MY) ) then
          up(:,Jmingh:Jmin-1,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MY) ) then
          up(:,Jmax+1:Jmaxgh,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MZ) ) then
          up(:,:,Kmingh:Kmin-1,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MZ) ) then
          up(:,:,Kmax+1:Kmaxgh,:) = zero
       endif
    enddo
  end subroutine vmg_boundary_fill0
  !-------------------------------------------------------------------------
  ! extrapolate to the grid boundaries
  !-------------------------------------------------------------------------
  subroutine vmg_boundary_extrap(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: id, fmglev
    integer :: i, j, k
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: up

    fmglev = FMG_Level

    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MX) ) then
          do i = Imin-1, Imingh, -1
             up(i,:,:,:) = -up(i+1,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MX) ) then
          do i = Imax+1, Imaxgh
             up(i,:,:,:) = -up(i-1,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MY) ) then
          do j = Jmin-1, Jmingh, -1
             up(:,j,:,:) = -up(:,j+1,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MY) ) then
          do j = Jmax+1, Jmaxgh
             up(:,j,:,:) = -up(:,j-1,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MZ) ) then
          do k = Kmin-1, Kmingh, -1
             up(:,:,k,:) = -up(:,:,k+1,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MZ) ) then
          do k = Kmax+1, Kmaxgh
             up(:,:,k,:) = -up(:,:,k-1,:)
          end do
       end if
    enddo
  end subroutine vmg_boundary_extrap
  !-------------------------------------------------------------------------
  ! extrapolate to the grid boundaries
  !-------------------------------------------------------------------------
  subroutine vmg_boundary_extrap_BAK(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: id, fmglev
    integer :: i, j, k
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: up

    fmglev = FMG_Level

    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MX) ) then
          do i = Imin-1, Imingh, -1
             up(i,:,:,:) = 2*up(i+1,:,:,:) - up(i+2,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MX) ) then
          do i = Imax+1, Imaxgh
             up(i,:,:,:) = 2*up(i-1,:,:,:) - up(i-2,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MY) ) then
          do j = Jmin-1, Jmingh, -1
             up(:,j,:,:) = 2*up(:,j+1,:,:) - up(:,j+2,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MY) ) then
          do j = Jmax+1, Jmaxgh
             up(:,j,:,:) = 2*up(:,j-1,:,:) - up(:,j-2,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,MZ) ) then
          do k = Kmin-1, Kmingh, -1
             up(:,:,k,:) = 2*up(:,:,k+1,:) - up(:,:,k+2,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,MZ) ) then
          do k = Kmax+1, Kmaxgh
             up(:,:,k,:) = 2*up(:,:,k-1,:) - up(:,:,k-2,:)
          end do
       end if
    enddo
  end subroutine vmg_boundary_extrap_BAK
  !-------------------------------------------------------------------------
  ! get maximum
  !-------------------------------------------------------------------------
  function vmg_get_max(amrlev, icode, absolute, skip, xskip, mpireduce) result(maxarr)
    integer,intent(IN) :: amrlev, icode
    real(kind=DBL_KIND) :: maxarr
    logical,intent(IN),optional :: absolute, skip, xskip, mpireduce
    call fmg_max(maxarr, amrlev, FMG_Level, icode, absolute, skip, xskip, mpireduce)
  end function vmg_get_max
  !-------------------------------------------------------------------------
  ! get sum of array
  !-------------------------------------------------------------------------
  subroutine vmg_sum(sumarr, amrlev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, icode
    real(kind=DBL_KIND),intent(OUT) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_sum(sumarr, amrlev, FMG_Level, icode, absolute, skip, mpireduce)
  end subroutine vmg_sum
  !-------------------------------------------------------------------------
  ! get sum of array
  !-------------------------------------------------------------------------
  subroutine vmg_ave(avearr, amrlev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, icode
    real(kind=DBL_KIND),intent(OUT) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_ave(avearr, amrlev, FMG_Level, icode, absolute, skip, mpireduce)
  end subroutine vmg_ave
  !-------------------------------------------------------------------------
  ! Array is normalized so that the average equals zero.
  ! 周期境界条件でポアソン方程式を解く場合、tauコレクションによって、密度の平均値が 0 ではなくなる。
  ! これを強制するために、ベースの密度場をゼロとする。
  !-------------------------------------------------------------------------
#define SZ is:ie,js:je,ks:ke,Mmin:Mmin
  subroutine vmg_zeroAverage(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    real(kind=DBL_KIND) :: avearr
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: arr
    integer :: fmglev, gid, is,js,ks,ie,je,ke
    fmglev = FMG_Level
    call vmg_ave(avearr, amrlev, icode, skip=.FALSE.)
!!$    print *, 'average density', avearr
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       arr => fmg_get_arrp(amrlev, fmglev, gid, icode)
       arr(SZ) = arr(SZ) - avearr
    enddo
    ! check
!!$    call vmg_ave(avearr, amrlev, icode, skip=.FALSE.)
!!$    print *, 'average density', avearr
  end subroutine vmg_zeroAverage
#undef SZ
  !-------------------------------------------------------------------------
  ! Poisson 方程式の誤差 max| res | / max| Δrho | の最大値を得る
  !-------------------------------------------------------------------------
  function vmg_get_error(jrhs) result(error)
    use mpilib
    integer,intent(IN) :: jrhs
    real(kind=DBL_KIND) :: error
    real(kind=DBL_KIND) :: rhoave, drho, drhog
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: rho
    real(kind=DBL_KIND),parameter :: eps = 1.D-7
    integer :: amrlev, gid
#define SZ Imin:Imax,Jmin:Jmax,Kmin:Kmax,:
    error = 0.d0
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call vmg_ave(rhoave, amrlev, jrhs) ! get rhoave
       ! get SUM( |drho| )
       drho = 0.d0
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          rho => fmg_get_arrp(amrlev,FMG_Level,gid,jrhs)
          drho = max(drho, maxval(ABS(rho(SZ) - rhoave)))
       end do
       call mpi_allreduce(drho, drhog, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       error = max(error, Resmaxg(amrlev)/(drhog*(1.d0+eps)))
    end do
#ifdef VERBOSE
    if (myrank == PRIMARY_RANK) print *, 'vmg error',error
#endif !VERBOSE

#undef SZ
  end function vmg_get_error
end module vmg
