module parameter
  implicit none
  private
  real(kind=8),save :: Pi, Pi2, Pi4, Pi8, Pii, Pi2i, Pi4i, Pi8i
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
  end subroutine parameter_init
end module parameter
