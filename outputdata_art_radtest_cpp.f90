module outputdata
  implicit none
  private
  public :: output_data
contains
  subroutine output_data
    use io
    use analysis, only : RhoMax
    use sinkParticle
    use writeSnap
    use uniformgrid, only : uniformgrid_write
    use grid, only : Lmin, LevelMax
    use modelParameter, only : MP_Boxsize 
    use mpilib 
    integer :: level
    real(kind=8) :: hw
    if (.not. bool_output() ) return
    call output_centralbox
  end subroutine output_data
  function bool_output() result(bool)
    use grid
    use eos
    use analysis, only : RhoMax
    use modelParameter, only : MP_Dstep
    logical :: bool
    real(kind=8),parameter :: logrho_skip = 0.2d0 
    integer(kind=8),parameter :: step_skip = 500 
    integer,parameter :: level_io = 10 
    real(kind=8),save :: logrhoio
    real(kind=8),save :: rhomax_prev = Huge(rhomax_prev)
    logical,save :: bool_output_initialized = .false.
    bool = .false.
    if (level_sync() == Lmin .and. mod(Step(Lmin),MP_Dstep) == 0 ) then 
       bool = .true.
    endif
    return
  end function bool_output
  subroutine output_centralbox()
    use grid, only : Lmin, LevelMax
    use parameter
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use mpilib
    character(len=2),parameter :: prefix = 'cb'
    real(kind=8) :: halfwidth
    integer,parameter :: NIug=32 
    integer,parameter :: l_range=20 
    integer :: level
    do level = Lmin-1,LevelMax
       if (level < LevelMax - l_range + 1) cycle
       halfwidth = MP_Boxsize / 2.**(level-Lmin) / (8*8/dble(NIug)) 
       if (halfwidth > MP_boxsize) cycle
       if (get_myrank() == 0) then
          print '(/,A,I0,A,(1P1E9.2),A)', "output_centralbox: level = ", &
               level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
       end if
       call uniformgrid_write(-halfwidth,-halfwidth,-halfwidth,halfwidth,halfwidth,halfwidth, level, interpolate=.false.,prefix=pre&
&fix)
    enddo
  end subroutine output_centralbox
  subroutine output_sinkParticles()
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use unit, only : Unit_msun
    use sinkParticle
    use mpilib
    use string, only : CHARLEN, num2char, concat
    integer,parameter :: NIug=32, NJug=NIug, NKug=NIug 
    real(kind=8) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=8) :: halfwidth
    integer :: level, np, n, dummy
    real(kind=8) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    character(len=CHARLEN) :: prefix
    integer,dimension(:),allocatable :: pid
    real(kind=8),dimension(:),allocatable :: pmass
    real(kind=8),dimension(:,:),allocatable :: pr
    call ob_computationBoxOfCoordPhys( coordPhys )
    np = sp_getNparticle()
    if (np > 0) then 
       allocate(pid(np), pmass(np), pr(0:2, np))
       call sp_sinkdata2array(dummy, pmass, pr=pr, pid=pid) 
    end if
    do n = 1, np
       if (pmass(n)*Unit_msun < 1d0) then
          if (get_myrank() == 0) then
             print '(/,A,I0,A,1P1E9.2,A)', "pid = ", pid(n), ", pmass = ",pmass(n)*Unit_msun, &
             ", skip sinkParticle with mass < 1 Msun"
          end if
          cycle
       end if
       do level = LevelMax-5,LevelMax 
          if (level <= -2) cycle
          halfwidth = MP_Boxsize / 2.**(level-Lmin) / (8*8/dble(NIug)) 
          if (get_myrank() == 0) then
             print '(/,A,I0,A,1P1E9.2,A,I0,A,(1P1E9.2),A)', "output_sinkParticles: pid = ", pid(n), " (pmass = ", &
             pmass(n)*Unit_msun ,"), level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
          end if
          xmin = coordPhys(0)
          ymin = coordPhys(1)
          zmin = coordPhys(2)
          xmax = coordPhys(2 +1+0)
          ymax = coordPhys(2 +1+1)
          zmax = coordPhys(2 +1+2)
          xmin = max(xmin, pr(0,n)-halfwidth)
          ymin = max(ymin, pr(1,n)-halfwidth)
          zmin = max(zmin, pr(2,n)-halfwidth)
          xmax = min(xmax, pr(0,n)+halfwidth)
          ymax = min(ymax, pr(1,n)+halfwidth)
          zmax = min(zmax, pr(2,n)+halfwidth)
          prefix = 'sp.'
          prefix = concat(concat(prefix, num2char(pid(n))),'.')
          call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)
       end do 
    end do 
  end subroutine output_sinkParticles
   end module outputdata
