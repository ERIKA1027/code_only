#include "config.h"
#include "planet.h"

!-------------------------------------------------------------------------
! Module for MHD parameters at the locations of the planets
! Calculation of planets' position is done only in the primary rank.
!-------------------------------------------------------------------------
module planet_insitu
!!$  use planet, only : t_planet
  implicit none
  private
!!$  type t_planet_insitu
!!$     type(t_planet),pointer :: planetp => null()
!!$     real(kind=DBL_KIND) :: x, y, z     ! position
!!$     real(kind=DBL_KIND) :: u(0:NM-1)   ! state vector
!!$  end type t_planet_insitu
!!$  logical :: Initialized
  integer :: CurrentLevel

  public :: planet_insitu_all
contains
  !-------------------------------------------------------------------------
  ! 
  !-------------------------------------------------------------------------
  subroutine planet_insitu_all(level)
    use dates, only : dates_timetoJd2000
    use planet, only : LIST_OF_PLANETS, planet_set_time
    use grid, only : Time
    integer,intent(IN) :: level ! current grid level for solving MHD
    CurrentLevel = level        ! module variable
    call planet_set_time(dates_timetoJd2000(Time(level)))
    call insitu_log(Mercury)
    call insitu_log(Venus)
    call insitu_log(Earth)
    call insitu_log(Mars)
    call insitu_log(Jupiter)
    call insitu_log(Saturn)
    call insitu_log(Uranus)
    call insitu_log(Neptune)
    call insitu_log(Pluto)
  end subroutine planet_insitu_all
  !-------------------------------------------------------------------------
  ! In situ observation of solar wind at a given planet.
  ! Results are written into log files. 
  !-------------------------------------------------------------------------
  subroutine insitu_log(plnt)
    use mpilib
    use string, only : CHARLEN
    use unit, only : Unit_au
    use grid, only : Mmin, Mmax, Step, Time, Undefi
    use overBlockCoordinates, only : t_obPointPhys, t_obRectPhys, ob_assignCoordPhysToPointPhys, ob_computationBoxOfRectPhys, ob_PointPhysWithinRectPhys
    use ob_interp, only : ob_interpolatedU
    use planet, only : t_planet, planet_get_pos, planet_set_pos
    type(t_planet),intent(INOUT) :: plnt
    real(kind=DBL_KIND) :: x_au, y_au, z_au
    real(kind=DBL_KIND),dimension(MX:MZ) :: coords
    integer,dimension(Mmin:Mmax) :: mlist
    real(kind=DBL_KIND),dimension(Mmin:Mmax) :: u
    integer,parameter :: LUN = 11
    integer :: level
    type(t_obPointPhys) :: pos
    type(t_obRectPhys) :: compbox
    if (get_myrank() == PRIMARY_RANK) then
       call planet_set_pos(plnt)
       call planet_get_pos(plnt, x_au, y_au, z_au) ! in au
       coords = (/x_au, y_au, z_au/) / Unit_au ! in simulation unit
    endif
    call mpi_bcast(coords, size(coords), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call ob_assignCoordPhysToPointPhys(coords, pos)
    call ob_computationBoxOfRectPhys(compbox)
    if (.not. ob_PointPhysWithinRectPhys(pos, compbox)) return
    mlist = (/MRHO, MVX, MVY, MVZ, MBX, MBY, MBZ, MP, MDB/)
    call ob_interpolatedU(pos, mlist, u, levelp=level)
    if (level /= CurrentLevel) return
    if (level == Undefi) return
    if (get_myrank() == PRIMARY_RANK) then
       open(LUN, file=get_filename(plnt), position='append')
       write(LUN, '(I12, I4, 1PE12.5, 3(1PE12.5), 9(1PE12.5))') Step(level), level, Time(level), coords, u
       close(LUN)
    endif
  end subroutine insitu_log
  !-------------------------------------------------------------------------
  ! get file name of log fle for each planet.
  !-------------------------------------------------------------------------
  function get_filename(plnt) result(filename)
    use string, only : CHARLEN
    use planet, only : t_planet
    use io_util, only : read_env
    type(t_planet),intent(IN) :: plnt
    character(len=CHARLEN) :: filename, dir
    call read_env('DIR', dir)
    filename = trim(adjustl(dir)) // 'insitu_' // trim(adjustl( plnt%name ))
  end function get_filename
  !-------------------------------------------------------------------------
  ! initialize module
  !-------------------------------------------------------------------------
!!$  subroutine planet_insitu_init
!!$    use planet, only : LIST_OF_PLANETS, t_planet
!!$    implicit none
!!$    integer :: n
!!$    if (Initialized) return
!!$    Initialized = .TRUE.
!!$  end subroutine planet_insitu_init

end module planet_insitu
