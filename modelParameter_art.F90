#include "config.h"
!-------------------------------------------------------------------------
! module for parameter for a sample code
!-------------------------------------------------------------------------
module modelParameter
  implicit none
  private
#define PUBLIC_VARIABLES MP_Gconst, MP_T_Last, MP_CloudDensity, MP_CloudPressure, MP_CloudCs
  real(kind=DBL_KIND),save :: PUBLIC_VARIABLES
  public :: PUBLIC_VARIABLES
  public :: modelParameter_init
contains
  subroutine modelParameter_init
    use mpilib
    use parameter, only :  Pi4i
    use io_util, only : read_env
    use unit
    real(kind=DBL_KIND) :: cs

    MP_Gconst = Pi4i
    call read_env('T_LAST', MP_T_Last)

    ! Background cloud
    MP_CloudDensity = ModelParam_rho / Unit_rho ! cloud density
    MP_CloudCs = 1.d0                            ! sound speed (= Unit_v)
    MP_CloudPressure = MP_CloudCs**2 * MP_CloudDensity / ModelParam_gamma ! pressure


    if (get_myrank() == PRIMARY_RANK) then
       PRINTV(MP_CloudDensity)
       PRINTV(MP_CloudCs)
       PRINTV(MP_CloudPressure)
    endif
  end subroutine modelParameter_init
end module modelParameter
