!=============================================================
!                       KINZOKU.F90
!             subroutine for metal effects         
!   
!
!=============================================================
#include "config.h"

module kinzoku
  use unit

  implicit none
  private

  !metallicity========
   real(kind=DBL_KIND),parameter :: Metallicity = 1.e0          ! Metallicity
  !===================

  !parameter for raytracing
  real(kind=DBL_KIND),parameter :: tau_max_dust = 1d2           ! maximum optical depth for raytracing
  real(kind=DBL_KIND),parameter :: frac_C_solar = 0.927d-4      ! fraction of C 
  real(kind=DBL_KIND),parameter :: frac_O_solar = 3.568d-4      ! fraction of O
  real(kind=DBL_KIND),parameter :: CONST_AC     = 12.011        ! Atomic weight of Carbon 
  real(kind=DBL_KIND),parameter :: CONST_AO     = 15.999        ! Atomic weight of Oxygen 

  real(kind=DBL_KIND),parameter :: frac_AHe     = 0.0972        ! fraction of He
  real(kind=DBL_KIND),parameter :: CONST_AHe    = 4.004         ! Atomic weight of He

  !constant for calculations
  real(kind=DBL_KIND),save :: frac_C, frac_O, yA_tot

  !Parameters for dust grains
  real(kind=DBL_KIND),parameter :: Qd_uv        = 1.d0          ! effective parameter for UV photons
  real(kind=DBL_KIND),parameter :: a_dust       = 1.d-5         ! typical dust grain radius [cm]
  real(kind=DBL_KIND),parameter :: rho_dust     = 3.e0          ! dust grain density [g cm^-3]
  real(kind=DBL_KIND),parameter :: f_dust       = 0.01          ! dust to gas mass ratio
    
  real(kind=DBL_KIND),save :: sig_duv, NcnH_max 


  public :: init_kinzoku
  public :: Metallicity, CONST_AC, CONST_AO, frac_AHe, CONST_AHe 
  public :: frac_C, frac_O, yA_tot
  public :: sig_duv, NcnH_max

contains
   
  !-----------------------------------------
  ! subroutine for initialize kinzoku 
  !-----------------------------------------

  subroutine init_kinzoku
    implicit none

    !fracion of each metal
    frac_C      = frac_C_solar * Metallicity
    frac_O      = frac_O_solar * Metallicity

    !total 
    yA_tot      = 1.d0 + frac_AHe * CONST_AHe + frac_C * CONST_AC + frac_O * CONST_AO 

    !cross section of dust grains
    sig_duv     = 3.d0*cgs_mp*Qd_uv*f_dust*Metallicity/(4.d0*a_dust*rho_dust) * yA_tot ![cm^2] 


    !max
    NcnH_max    = tau_max_dust / sig_duv ! maximum colomun density for raytracing 



  end subroutine init_kinzoku



end module kinzoku

