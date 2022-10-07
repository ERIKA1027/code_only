
#include "config.h"
!-------------------------------------------------------------------------
! Module for general utility
!-------------------------------------------------------------------------
module util
  implicit none
contains
  !-----------------------------------------------------------------------
  ! arrey offset
  ! when ncrd = MX, io = 1
  ! when ncrd = MY, jo = 1
  ! when ncrd = MZ, jo = 1
  !-----------------------------------------------------------------------
  subroutine util_arroffset(ncrd,io,jo,ko)
    integer,intent(IN) :: ncrd
    integer,intent(OUT) :: io,jo,ko
    io = 0
    jo = 0
    ko = 0
    select case(ncrd)
    case ( MX )                 ! x-direction
       io = 1
    case ( MY )                 ! y-direction
       jo = 1
    case ( MZ )                 ! z-direction
       ko = 1
    case default
       write(*,*) '*** bad ncrd in subroutine arroffset'
       stop
    end select
  end subroutine util_arroffset
  !-------------------------------------------------------------------------
  ! bool が 真である index を見つける
  !-------------------------------------------------------------------------
  function util_whereInt1D( bool , last) result( index )
    logical,dimension(0:),intent(IN) :: bool
    logical,intent(IN),optional :: last
    integer,dimension( ARRAYSIZE1( bool ) ) :: unit
    integer :: index, loc(1), n
    ! defalut is last=.FALSE.
    if ( present( last ) ) then
       do n = lbound(bool,1), ubound(bool,1)
          unit(n) = n
       enddo
    else
       unit(:) = 1
    endif
    if ( any( bool ) ) then
       loc = maxloc(unit, bool)
       index = loc(1) -1 + lbound(bool, 1)
    else
       index = lbound(bool, 1) - 1
    endif
  end function util_whereInt1D
  !-------------------------------------------------------------------------
  ! n が2のべき乗の場合 真を返す。それ以外は偽を返す。
  !-------------------------------------------------------------------------
  function util_isPowerOfTow(n) result(bool)
    integer,intent(IN) :: n
    logical :: bool
    bool = .FALSE.
    if (n /= 0 .and. IAND(n, n-1) == 0) bool = .TRUE.
  end function util_isPowerOfTow
end module util
!-------------------------------------------------------------------------
subroutine util_print_hwm
  integer,parameter :: lun=92
  character(len=80) :: line
  integer(kind=8) :: hwm,rss
  open(lun,file='/proc/self/status')
  hwm=0
  rss=0
  do while(.true.)
     read(lun,'(a)',end=99) line
     if (line(1:6).eq.'VmHWM:') read(line(8:80),*) hwm
     if (line(1:6).eq.'VmRSS:') read(line(8:80),*) rss
  enddo
   99 close(lun)
  print *,'HWM = ',hwm,'kB =',hwm*1024,'bytes'
  print *,'RSS = ',rss,'kB =',rss*1024,'bytes'
end subroutine util_print_hwm
