#include "config.h"
!-------------------------------------------------------------------------
! module for parameter
!-------------------------------------------------------------------------
module parameter
  implicit none
  private
  real(kind=DBL_KIND),save :: Pi, Pi2, Pi4, Pi8, Pii, Pi2i, Pi4i, Pi8i
  public :: parameter_init
  public :: Pi, Pi2, Pi4, Pi8, Pii, Pi2i, Pi4i, Pi8i
contains
  subroutine parameter_init
    implicit none
    Pi = 4 * atan(1.D0)
    Pi2 = 2 * Pi
    Pi4 = 2 * Pi2
    Pi8 = 2 * Pi4
    Pii = 1 / Pi
    Pi2i = 1 / Pi2
    Pi4i = 1 / Pi4
    Pi8i = 1 / Pi8

!!$    MP_Gconst = Pi4i
!!$    call read_env('T_LAST', MP_T_Last)
!!$  T_LAST = 0.22360680d0         ! ½ªÎ»»þ¹ï
!!$  T_LAST = 6.d0            ! ½ªÎ»»þ¹ï
!!$  T_LAST = 1.d0            ! ½ªÎ»»þ¹ï
  end subroutine parameter_init
end module parameter
