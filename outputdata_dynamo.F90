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
    use writeSnap
!!$    use uniformgrid, only : uniformgrid_write

    if (.not. bool_output() ) return

    call writeSnap_whole

  end subroutine output_data
  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming.
  ! -----------------------------------------------------------------
  function bool_output() result(bool)
    use grid
    use eos
    logical :: bool
    integer,parameter :: ioskip = 100           ! output interval

    bool = .false.

    if ( level_sync() /= Lmin ) return

    if (maxval(Step) == 0) bool = .true. ! initial codition

    if ( Step(Lmin) > 100 .and. Step(Lmin) < 200 .and. mod(Step(Lmin),10) == 0 ) then
       bool = .true.
    endif
    return

    if ( mod(Step(Lmin),ioskip) == 0 ) then
       bool = .true.
    endif
    return

  end function bool_output
end module outputdata
