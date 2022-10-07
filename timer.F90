
#include "config.h"
!-------------------------------------------------------------------------
! Module for timer.
! Estimate elapse times.
! Example
! call tm_init  ! Initialize timer module at the beginning of program
!   ...
! call tm_start(TIMER_SECT_IO)     ! estimate exec time of code A
!  code A
! call tm_end  !TIMER_SECT_IO
!  ...
! call tm_start(TIMER_SECT_PKSV)   ! estimate exec time of codes B + D
!  code B
! call tm_start(TIMER_SECT_PKWAIT) ! estimate exec time of codes C
!  code C
! call tm_end ! TIMER_SECT_PKWAIT
!  code D
! call tm_end ! TIMER_SECT_PKSV
!  ...
! call tm_finalize   ! Finalize timer module at the end of program
!-------------------------------------------------------------------------
module timer
  implicit none
  private
  ! Constant parameters for section
  integer,parameter :: TIMER_SECT_MAIN = 0
  integer,parameter :: TIMER_SECT_IO   = 1
  integer,parameter :: TIMER_SECT_HD   = 2
  integer,parameter :: TIMER_SECT_GRAV = 3
  integer,parameter :: TIMER_SECT_FMG  = 4
  integer,parameter :: TIMER_SECT_VMG  = 5
  integer,parameter :: TIMER_SECT_MG   = 6
  integer,parameter :: TIMER_SECT_FMGIF= 7
  integer,parameter :: TIMER_SECT_FMGLIN= 8
  integer,parameter :: TIMER_SECT_MGIF  = 9
  integer,parameter :: TIMER_SECT_FMGCONV  = 10
  integer,parameter :: TIMER_SECT_FMGCONVF  = 11
  integer,parameter :: TIMER_SECT_FMGCONVI  = 12
  integer,parameter :: TIMER_SECT_FMGCONVL  = 13
  integer,parameter :: TIMER_SECT_PKSV      = 14
  integer,parameter :: TIMER_SECT_PKPUSH    = 15
  integer,parameter :: TIMER_SECT_PKPOP     = 16
  integer,parameter :: TIMER_SECT_PKLEN     = 17
  integer,parameter :: TIMER_SECT_PKTEST    = 18
  integer,parameter :: TIMER_SECT_PKWAIT    = 19
  ! Constant parameters for XCLOCK
  integer,parameter :: TIMER_GETTIME_STR = 1
  integer,parameter :: TIMER_GETTIME_INT = 2
  integer,parameter :: TIMER_CPUTIME_START = 3
  integer,parameter :: TIMER_CPUTIME_SEC = 4
  integer,parameter :: TIMER_CPUTIME_LAP_SEC = 5
  integer,parameter :: TIMER_CPUTIME_REM_SEC = 6
  integer,parameter :: TIMER_ETIME_START = 7
  integer,parameter :: TIMER_ETIME_SEC = 8
  ! stack for level of nest
  integer,parameter :: NSTACK = 100
  integer,save,dimension(0:NSTACK) :: NestStack
  integer,save :: NestLev = 0   ! stack pointer
  ! Time for each section
  integer,parameter :: KIND_TIME = 8
  integer,parameter :: NSECTION = 100
  real(kind=KIND_TIME),save,dimension(0:NSECTION) :: Time
  real(kind=KIND_TIME),save :: TimeStamp
  public :: TIMER_SECT_MAIN, TIMER_SECT_IO, TIMER_SECT_HD, &
       TIMER_SECT_GRAV, TIMER_SECT_FMG, TIMER_SECT_VMG, TIMER_SECT_MG, &
       TIMER_SECT_FMGIF, TIMER_SECT_FMGLIN, TIMER_SECT_MGIF, &
       TIMER_SECT_FMGCONV, TIMER_SECT_FMGCONVF, TIMER_SECT_FMGCONVI, TIMER_SECT_FMGCONVL, &
       TIMER_SECT_PKSV, TIMER_SECT_PKPUSH, TIMER_SECT_PKPOP, TIMER_SECT_PKLEN, &
       TIMER_SECT_PKTEST, TIMER_SECT_PKWAIT, &
       tm_init, tm_finalize, tm_start, tm_end, tm_showTime
contains
  !-------------------------------------------------------------------------
  ! Initialize timer
  !-------------------------------------------------------------------------
  subroutine tm_init
    real(kind=KIND_TIME) :: dummy
    NestLev = TIMER_SECT_MAIN
    NestStack(:) = 0
    TimeStamp = 0.0
    call xclock(dummy, TIMER_ETIME_START)
    call xclock(TimeStamp, TIMER_ETIME_SEC)
  end subroutine tm_init
  !-------------------------------------------------------------------------
  ! Finalize timer
  !-------------------------------------------------------------------------
  subroutine tm_finalize
    real(kind=KIND_TIME) :: t
    call xclock(t, TIMER_ETIME_SEC)
    Time(TIMER_SECT_MAIN) = Time(TIMER_SECT_MAIN) + t - TimeStamp
    if (get_sectionNow() /= 0) print *, '*** error in tm_finalize. section is', get_sectionNow()
    call tm_showTime
  end subroutine tm_finalize
  !-------------------------------------------------------------------------
  ! start timer for a given section
  !-------------------------------------------------------------------------
  subroutine tm_start(section)
    integer,intent(IN) :: section
    real(kind=KIND_TIME) :: t
    integer :: sectionParent, sectionChild
    sectionParent = get_sectionNow()
    sectionChild  = section
    call xclock(t, TIMER_ETIME_SEC)
    Time(sectionParent) = Time(sectionParent) + t - TimeStamp
    TimeStamp = t
    call pushSection(sectionChild)
  end subroutine tm_start
  !-------------------------------------------------------------------------
  ! end timer for a given section
  !-------------------------------------------------------------------------
  subroutine tm_end
    real(kind=KIND_TIME) :: t
    integer :: sectionChild
    call xclock(t, TIMER_ETIME_SEC)
    call popSection(sectionChild)
    Time(sectionChild) = Time(sectionChild) + t - TimeStamp
    TimeStamp = t
  end subroutine tm_end
  !-------------------------------------------------------------------------
  ! push section to stack
  !-------------------------------------------------------------------------
  subroutine pushSection(section)
    integer,intent(IN) :: section
    NestLev = NestLev + 1
    if (NestLev > NSTACK) print *, '*** error in timer::pushSection. NestStack is overflow.',NestLev
    NestStack(NestLev) =  section
  end subroutine pushSection
  !-------------------------------------------------------------------------
  ! pop from stack
  !-------------------------------------------------------------------------
  subroutine popSection(section)
    integer,intent(OUT) :: section
    section = NestStack(NestLev)
    NestLev = NestLev - 1
    if (NestLev < 0) print *, '*** error in timer::popSection. NestStack is underflow.',NestLev
  end subroutine popSection
  !-------------------------------------------------------------------------
  ! pop from stack
  !-------------------------------------------------------------------------
  function get_sectionNow() result(section)
    integer :: section
    section = NestStack(NestLev)
  end function get_sectionNow
  !-------------------------------------------------------------------------
  ! print times for all sections
  !-------------------------------------------------------------------------
  subroutine tm_showTime
    use mpilib
    integer :: n
    if (get_myrank() == PRIMARY_RANK) then
       do n = 0, NSECTION
          if (Time(n) /= 0.0) then
             print *, 'section', n, Time(n)
          endif
       end do
       print *, 'NestLev', NestLev
    end if
  end subroutine tm_showTime
  !-------------------------------------------------------------------------
  ! emulate xclock
  !-------------------------------------------------------------------------
  subroutine xclock(p, q)
    real(kind=KIND_TIME),intent(OUT) :: p
    integer,intent(IN) :: q
    if (q == TIMER_ETIME_START) then
       p = 0
    else
       p = p + 1
    endif
  end subroutine xclock
end module timer
