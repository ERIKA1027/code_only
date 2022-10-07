#include "config.h"
!-------------------------------------------------------------------------
!
! Module for controlling output data
!
!-------------------------------------------------------------------------
module outputdata
  implicit none
  private
  public :: output_data
contains
  ! -----------------------------------------------------------------
  ! output procedures 
  ! -----------------------------------------------------------------
  subroutine output_data
    use io
#ifndef MP
    use eos, only : Rhocr
#endif !MP
    use analysis, only : RhoMax
#ifdef SINKPARTICLE
    use sinkParticle
#endif !SINKPARTICLE
    use writeSnap
    use uniformgrid, only : uniformgrid_write
    use grid, only : Lmin, LevelMax
    use modelParameter, only : MP_Boxsize ! KS ADDED
    use mpilib ! KS ADDED

    integer :: level
    real(kind=DBL_KIND) :: hw

    if (.not. bool_output() ) return

    ! output whole computational domain
    ! if (get_myrank() == PRIMARY_RANK) then
    !    print '(/,A,I0)', "writeSnap_whole: level = ", -1
    ! end if
    ! call writeSnap_whole

    ! output central box
    call output_centralbox

    ! output regions around sink particles 
    ! call output_sinkParticles ! (comment out for rad_test)

  end subroutine output_data

  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming.
  ! -----------------------------------------------------------------
  function bool_output() result(bool)
    use grid
    use eos
    use analysis, only : RhoMax
    use modelParameter, only : MP_Dstep
    logical :: bool
#ifndef SINGLE_STEP
    integer,parameter :: ioskip = 100           ! output interval
    integer(kind=LLONG_KIND),save :: ioStep = 0
#else
    real(kind=DBL_KIND),parameter :: logrho_skip = 0.2d0 ! interval for output in unit of log(rho_max)
    integer(kind=LLONG_KIND),parameter :: step_skip = 500 ! interval for output in unit of step number
    integer,parameter :: level_io = 10                    ! grid level where step_skip is applied.
    real(kind=DBL_KIND),save :: logrhoio
    real(kind=DBL_KIND),save :: rhomax_prev = Huge(rhomax_prev)
    logical,save :: bool_output_initialized = .false.
#endif

    ! -------------
    ! for debug 
    ! -------------
    bool = .false.
    !MP_Dstep毎にデータを出力
    if (level_sync() == Lmin .and. mod(Step(Lmin),MP_Dstep) == 0 ) then !KS MODIFIED
       bool = .true.
    endif
    return
  end function bool_output

  ! -----------------------------------------------------------------
  ! Output central box.
  ! From level = Lmin-1 to LevelMax
  ! File name is cb step . level . d
  ! -----------------------------------------------------------------
  subroutine output_centralbox()
    use grid, only : Lmin, LevelMax
    use parameter
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use mpilib

    
    character(len=2),parameter :: prefix = 'cb'
    real(kind=DBL_KIND) :: halfwidth
    !integer,parameter :: NIug=32 ! data points in each direction (KS MODIFIED, ~ 2 MB/file)
    integer,parameter :: NIug=32 ! KS DEBUG
    integer,parameter :: l_range=20 ! num of levels to be output
    integer :: level
    

    do level = Lmin-1,LevelMax
       if (level < LevelMax - l_range + 1) cycle
       halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction       
       if (halfwidth > MP_boxsize) cycle
       if (get_myrank() == PRIMARY_RANK) then
          print '(/,A,I0,A,(1P1E9.2),A)', "output_centralbox: level = ", &
               level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
       end if
       call uniformgrid_write(-halfwidth,-halfwidth,-halfwidth,halfwidth,halfwidth,halfwidth, level, interpolate=.false.,prefix=prefix)
    enddo
  end subroutine output_centralbox

  ! -----------------------------------------------------------------
  ! Output regions around sink particles.
  ! From level = Lmin-1 to LevelMax
  ! File name is sp. pid . step . level . d
  ! -----------------------------------------------------------------
  subroutine output_sinkParticles()
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use unit, only : Unit_msun
    use sinkParticle
    use mpilib
    use string, only : CHARLEN, num2char, concat

    integer,parameter :: NIug=32, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 2 MB/file)
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=DBL_KIND) :: halfwidth
    integer :: level, np, n, dummy
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    character(len=CHARLEN) :: prefix
    integer,dimension(:),allocatable :: pid
    real(kind=DBL_KIND),dimension(:),allocatable :: pmass
    real(kind=DBL_KIND),dimension(:,:),allocatable :: pr

    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )

    ! pr ... location of sinkParticle
    ! np ... number of sinkParticle
    np = sp_getNparticle()

    if (np > 0) then     ! if sink particles exist
       allocate(pid(np), pmass(np),  pr(MX:MZ, np))
       call sp_sinkdata2array(dummy, pmass, pr=pr, pid=pid) ! from sink particle
    end if

    do n = 1, np
       !massが小さい場合はスキップ (ひとまず1 M_sunを境界にする)
       if (pmass(n)*Unit_msun < 1d0) then
          if (get_myrank() == PRIMARY_RANK) then
             print '(/,A,I0,A,1P1E9.2,A)', "pid = ", pid(n), ", pmass = ",pmass(n)*Unit_msun, &
             ", skip sinkParticle with mass < 1 Msun"
          end if
          cycle
       end if

       do level = LevelMax-5,LevelMax !最大レベルから6レベル分をプロット
          if (level <= -2) cycle 
          halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction in level -1
          if (get_myrank() == PRIMARY_RANK) then
             print '(/,A,I0,A,1P1E9.2,A,I0,A,(1P1E9.2),A)', "output_sinkParticles: pid = ", pid(n), " (pmass = ", &
             pmass(n)*Unit_msun ,"), level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
          end if

          xmin = coordPhys(MX)
          ymin = coordPhys(MY)
          zmin = coordPhys(MZ)
          xmax = coordPhys(MZ+1+MX)
          ymax = coordPhys(MZ+1+MY)
          zmax = coordPhys(MZ+1+MZ)
          ! 領域の AND をとる
          xmin = max(xmin, pr(MX,n)-halfwidth)
          ymin = max(ymin, pr(MY,n)-halfwidth)
          zmin = max(zmin, pr(MZ,n)-halfwidth)
          xmax = min(xmax, pr(MX,n)+halfwidth)
          ymax = min(ymax, pr(MY,n)+halfwidth)
          zmax = min(zmax, pr(MZ,n)+halfwidth)

          ! speicfy prefix
          prefix = 'sp.'
          prefix = concat(concat(prefix, num2char(pid(n))),'.')

          call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)

       end do ! level
    end do ! n
  end subroutine output_sinkParticles


   end module outputdata
