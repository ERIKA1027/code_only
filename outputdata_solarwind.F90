#include "config.h"
! #define BENCHMARK
#define DEBUG
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
    use grid, only : Lmin, LevelMax
    use io
    use writeSnap
    use uniformgrid, only : uniformgrid_write, uniformgrid_write
    use modelParameter, only : MP_rOuterBoundary
    integer :: lev
    real(kind=DBL_KIND) :: wlen
#ifdef BENCHMARK
    return
#endif !BENCHMARK

    if (bool_output()) then     ! daily output
       do lev = Lmin, LevelMax
          wlen = MP_rOuterBoundary / 2**(lev-Lmin)
          call uniformgrid_write(-wlen,-wlen,-wlen,wlen,wlen,wlen,lev, interpolate=.true.)
       enddo
    endif

    if (bool_output_monthly()) then ! monthly output
       call dumpdata
    endif

!!$    call writeSnap_whole

  end subroutine output_data
  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming.
  ! -----------------------------------------------------------------
  function bool_output() result(bool)
    use modelParameter, only : MP_CarringtonRotDay
    use unit, only : Unit_day
    use grid
    use eos
    use analysis
    logical :: bool
    real(kind=DBL_KIND) :: ioTimeSkip
    integer,save :: nioPrev = -1 ! counter for prevous IO
    integer,save :: nio = 0      ! counter for IO

#if defined(DEBUG)
    ! -------------
    ! for debug
    ! -------------
    bool = .false.
    if (level_sync() == Lmin ) then
       nio = nio + 1
       if ( mod(nio, 5) == 0 ) bool = .true.
    endif
    if (maxval(Step) == 0) bool = .true. ! initial codition
    return
#else

!!$    ioTimeSkip = MP_CarringtonRotDay / Unit_day / 12 ! 12 frame/revolution
    ioTimeSkip = 1.d0 / Unit_day ! every day

    bool = .FALSE.
    if (level_sync() /= Lmin ) return

    nio = int(Time(LevelMax)/ioTimeSkip)
    if ( nioPrev < 0 ) nioPrev = nio
    if ( nio >= nioPrev  + 1  ) then
       nioPrev = nio
       bool = .TRUE.
    endif
#endif
    return
  end function bool_output
  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming. (monthly)
  ! -----------------------------------------------------------------
  function bool_output_monthly() result(bool)
    use dates, only : dates_time2date
    use grid
    logical :: bool
    integer :: year, month, day, hour, mnt, sec
    integer,save :: year_done=0, month_done=0
    logical,save :: initialized = .FALSE.
    
    bool = .FALSE.
    if (maxval(Step) == 0) return ! initial codition
    if (.not. initialized) then  
       call dates_time2date(Time(Lmin), year, month, day, hour, mnt, sec)
       year_done  = year
       month_done = month
       initialized = .TRUE.
       return
    endif

    if (level_sync() /= Lmin ) return

    call dates_time2date(Time(Lmin), year, month, day, hour, mnt, sec)
    if ( year == year_done .and. month == month_done ) return ! already outputed
    year_done  = year
    month_done = month
    bool = .TRUE.

  end function bool_output_monthly
end module outputdata
