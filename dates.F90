#include "config.h"
! #define USE_MAIN_PROGRAM
!-------------------------------------------------------------------------
! A module for dates
!-------------------------------------------------------------------------
module dates
  implicit none
  private
  ! 
  integer,parameter :: JDN_15821015 = 2299161
  ! Julian date number of 1582/10/15.  It is a date when changing a calender 
  ! from the Julian calendaer to the Gregorian calender.

  real(kind=DBL_KIND),parameter :: JD2000_OFFSET = 2451545.d0
  ! Offset of JD2000 from JD, corresponding to Julian day of 2000/01/01
  ! Franz & Harper 2002 (DOI:10.1016/S0032-0633(01)00119-2)

  real(kind=DBL_KIND),parameter :: ORIGIN_OF_SIM_TIME_IN_JD2K = 0.d0
  ! The origin of simlation time in unit of JD2000.

  integer,parameter :: CHARLEN_YYYYMMDD = 8
  ! length of characters for YYYYMMDD

  logical,save :: Initialized = .FALSE.

  public :: CHARLEN_YYYYMMDD, ORIGIN_OF_SIM_TIME_IN_JD2K

  public :: &
       dates_jd2000toDate, &
       dates_datetoJd2000, &
       dates_jdtoJd2000, &
       dates_jd2000toJd, &
       dates_datetoJd, &
       dates_jdtoDate, &
       dates_timetoJd2000, &
       dates_jd2000toTime, &
       dates_time2yyyymmdd, &
       dates_yyyymmdd2time, &
       dates_time2date, &
       dates_date2time
contains
  !-------------------------------------------------------------------------
  ! This routine converts JD200 to variables of dates (year, month, day,
  ! hour, minute, second).
  ! Julian day 2000 (JD2000) is JD whose base is 1200UT January 1st, 2000.
  !   JD2000 = JD - 2451545.0
  !-------------------------------------------------------------------------
  subroutine dates_jd2000toDate(jd2000, year, month, day, hour, mnt, sec)
    real(kind=DBL_KIND), intent(in) :: jd2000
    integer, intent(out) :: year, month, day, hour, mnt, sec
    real(kind=DBL_KIND) :: jd
    jd = dates_jd2000toJd(jd2000)
    call dates_jdtoDate(jd, year, month, day, hour, mnt, sec)
  end subroutine dates_jd2000toDate
  !-------------------------------------------------------------------------
  subroutine dates_datetoJd2000(year, month, day, hour, mnt, sec, jd2000)
    integer, intent(in) :: year, month, day, hour, mnt, sec
    real(kind=DBL_KIND), intent(out) :: jd2000
    real(kind=DBL_KIND) :: jd
    call dates_datetoJd(year, month, day, hour, mnt, sec, jd)
    jd2000 = dates_jdtoJd2000(jd)
  end subroutine dates_datetoJd2000
  !-------------------------------------------------------------------------
  ! Convertion between JD and JD2000
  !-------------------------------------------------------------------------
  function dates_jdtoJd2000(jd) result(jd2000)
    real(kind=DBL_KIND),intent(IN) :: jd
    real(kind=DBL_KIND) :: jd2000
    jd2000 = jd - JD2000_OFFSET
  end function dates_jdtoJd2000
  !-------------------------------------------------------------------------
  function dates_jd2000toJd(jd2000) result(jd)
    real(kind=DBL_KIND),intent(IN) :: jd2000
    real(kind=DBL_KIND) :: jd
    jd = jd2000 + JD2000_OFFSET
  end function dates_jd2000toJd
  !-------------------------------------------------------------------------
  ! Convert date to JD
  ! https://en.wikipedia.org/wiki/Julian_day#Calculation
  !-------------------------------------------------------------------------
  subroutine dates_datetoJd(year, month, day, hour, mnt, sec, jd)
    integer,intent(IN) :: year, month, day, hour, mnt, sec
    real(kind=DBL_KIND),intent(OUT) :: jd
    integer :: a, y, m, jdn
    a = int((14-month)/12)
    y = year + 4800 - a
    m = month + 12*a - 3
    ! Gregorian calender
    jdn = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045
    ! Julian calendaer
    if (jdn < JDN_15821015) jdn = day + (153*m+2)/5 + 365*y + y/4 - 32083
    jd = jdn + (hour - 12.d0)/24.d0 + mnt/1440.d0 + sec/86400.d0
  end subroutine dates_datetoJd
  !-------------------------------------------------------------------------
  ! Convert JD to date
  ! https://en.wikipedia.org/wiki/Julian_day#Calculation
  !-------------------------------------------------------------------------
  subroutine dates_jdtoDate(jd, year, month, day, hour, mnt, sec)
    real(kind=DBL_KIND),intent(IN) :: jd
    integer,intent(OUT) :: year, month, day, hour, mnt, sec
    integer :: f, e, g, h, jdn
    real(kind=DBL_KIND) :: resid
    jdn = int(jd + 0.5d0)
    if (jdn < JDN_15821015) then ! Julian calender
       f = jdn + 1401
    else                        ! Gregorian calender
       f = jdn + 1401 + (((4 * jdn + 274277)/146097)*3)/4-38
    endif
    e = 4 * f + 3
    g = mod(e, 1461) / 4
    h = 5*g + 2
    day = (mod(h, 153))/5 + 1
    month = mod(h/153 + 2, 12) + 1
    year = e/1461 - 4716 + (12 + 2 - month)/12
    resid = jd + 0.5d0 - jdn    ! residual within a day
    hour = int(resid*24)
    mnt = int((resid*24 - hour)*60)
    sec = int(((resid*24 - hour)*60 - mnt)*60 + 0.5d0) ! 0.5d0 is a margin
  end subroutine dates_jdtoDate
  !-------------------------------------------------------------------------
  ! Conversion between simulation time and JD2000
  !-------------------------------------------------------------------------
  function dates_timetoJd2000(time) result(jd2000)
    use unit, only : Unit_day
    real(kind=DBL_KIND),intent(IN) :: time
    real(kind=DBL_KIND) :: jd2000
    call dates_init
    jd2000 = time * Unit_day + ORIGIN_OF_SIM_TIME_IN_JD2K
  end function dates_timetoJd2000
  !-------------------------------------------------------------------------
  function dates_jd2000toTime(jd2000) result(time)
    use unit, only : Unit_day
    real(kind=DBL_KIND),intent(IN) :: jd2000
    real(kind=DBL_KIND) :: time
    call dates_init
    time = (jd2000 - ORIGIN_OF_SIM_TIME_IN_JD2K)/ Unit_day
  end function dates_jd2000toTime
  !-------------------------------------------------------------------------
  ! Conversion between simulation time and a string of YYYYMMDD
  !-------------------------------------------------------------------------
  function dates_time2yyyymmdd(time) result(yyyymmdd)
    real(kind=DBL_KIND),intent(IN) :: time
    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd    
    integer :: year, month, day, hour, mnt, sec
    call dates_jd2000toDate(dates_timetoJd2000(time), year, month, day, hour, mnt, sec)
    write(yyyymmdd, '(I4.4,2I2.2)') year, month, day
  end function dates_time2yyyymmdd
  !-------------------------------------------------------------------------
  function dates_yyyymmdd2time(yyyymmdd) result(time)
    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
    real(kind=DBL_KIND) :: time
    integer :: year, month, day, hour, mnt, sec
    real(kind=DBL_KIND) :: jd2000
    read(yyyymmdd(1:4),*) year
    read(yyyymmdd(5:6),*) month
    read(yyyymmdd(7:8),*) day
    hour = 0
    mnt = 0
    sec = 0
!!$    print *,  '**',year, month, day, hour, mnt, sec
    call dates_datetoJd2000(year, month, day, hour, mnt, sec, jd2000)
    time = dates_jd2000toTime(jd2000)
  end function dates_yyyymmdd2time
  !-------------------------------------------------------------------------
  ! Conversion between simulation time and a date
  !-------------------------------------------------------------------------
  subroutine dates_time2date(time, year, month, day, hour, mnt, sec)
    real(kind=DBL_KIND), intent(IN) :: time
    integer, intent(OUT) :: year, month, day, hour, mnt, sec
    call dates_jd2000toDate(dates_timetoJd2000(time), year, month, day, hour, mnt, sec)
  end subroutine dates_time2date
  !-------------------------------------------------------------------------
  subroutine dates_date2time(year, month, day, hour, mnt, sec, time)
    integer, intent(IN) :: year, month, day, hour, mnt, sec
    real(kind=DBL_KIND), intent(OUT) :: time
    real(kind=DBL_KIND) :: jd2000
    call dates_datetoJd2000(year, month, day, hour, mnt, sec, jd2000)
    time = dates_jd2000toTime(jd2000)
  end subroutine dates_date2time
  !-------------------------------------------------------------------------
  ! Initialized dependent routine.
  !-------------------------------------------------------------------------
  subroutine dates_init
    use parameter
    use unit
    if (Initialized) return
    Initialized = .TRUE.
    call parameter_init
    call unit_init
  end subroutine dates_init
end module dates

#ifdef USE_MAIN_PROGRAM
!-------------------------------------------------------------------------
! A main program for test the module.
! Utilty for JD is available here:
!  http://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/date2jd.cgi
!-------------------------------------------------------------------------
program main
  use dates
  implicit none
  real(kind=DBL_KIND) :: jd2000
  real(kind=DBL_KIND) :: jd
  integer :: year, month, day, hour, mnt, sec, length
  character(len=CHARLEN_YYYYMMDD) :: yyyymmdd
  character(:), allocatable :: arg

!!$  jd = 2457877.34309d0          !2017/05/03	20:14:03
!!$  jd = 2448015.34309d0          !1990/05/03	20:14:03
!!$  jd = 2451545.00000d0          !2000/01/01	12:00:00	
  jd = 2488128.50000d0         !2100/03/01	0:00:00
  jd2000 = dates_jdtoJd2000(jd)
  call dates_jd2000toDate(jd2000, year, month, day, hour, mnt, sec)
  print *, jd, year, month, day, hour, mnt, sec


!!$  call dates_datetoJd(year, month, day, hour, mnt, sec, jd)
  call dates_datetoJd(1000, 10, 15, 13, 34, 9, jd)
  print *, jd
  call dates_jdtoDate(jd, year, month, day, hour, mnt, sec)
  print *, year, month, day, hour, mnt, sec

  yyyymmdd = dates_time2yyyymmdd(jd2000*24.)
  print *, yyyymmdd
  jd2000 = dates_yyyymmdd2time(yyyymmdd) / 24.d0
  print *, dates_jd2000toJd(jd2000)

  call get_command_argument(0, length=length)
  allocate(character(length) :: arg)
  call get_command_argument(0, arg)
  print *, arg

end program main
#endif !USE_MAIN_PROGRAM

