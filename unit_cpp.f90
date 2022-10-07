module unit
  implicit none
  private
  real(kind=8),parameter :: ModelParam_rho = 1.D-18 
  real(kind=8),parameter :: ModelParam_temp = 10000.d0 
  real(kind=8),parameter :: ModelParam_gamma = 1.4d0 
  real(kind=8),parameter :: MP_hnu_IR = 0.01d0*1.602176463158D-12 
  real(kind=8),save :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc &
    , Unit_yr, Unit_msun, Unit_gauss, Unit_ugauss
  real(kind=8),save :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc
  real(kind=8),save :: cgs_au, cgs_msun, cgs_rsun, cgs_sigma, cgs_lsun, cgs_ev, cgs_hpk, cgs_asb 
  real(kind=8),save :: Unit_acc, Unit_invmu, Unit_rhomu, Unit_lrhomu, Unit_l3, Unit_erg, cgs_yHe
  public :: unit_init
  public :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_gauss, Unit&
&_ugauss
  public :: Unit_msun,Unit_acc, Unit_invmu,Unit_rhomu, Unit_lrhomu, Unit_l3
  public :: ModelParam_rho, ModelParam_temp, ModelParam_gamma
  public :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun, cgs_sigma, &
&cgs_lsun, cgs_ev 
  public :: cgs_hpk, cgs_asb, Unit_erg
  public :: cgs_yHe, MP_hnu_IR
contains
  subroutine unit_init
    use mpilib
    use parameter
    implicit none
    cgs_c = 2.99792458D10 
    cgs_amu = 1.660538921D-24 
    cgs_mp = 1.672621777D-24 
    cgs_me = 9.10938291D-28 
    cgs_mh = cgs_mp + cgs_me 
    cgs_mu = 2.3 
    cgs_yr = 3.1556926D7 
    cgs_gc = 6.67384D-8 
    cgs_kb = 1.3806488D-16 
    cgs_au = 1.495978707D13 
    cgs_pc = cgs_au * 180.d0 * 3600.d0 / Pi 
    cgs_msun = 1.9884D33 
    cgs_rsun = 6.960D10 
    cgs_sigma = 5.67051D-5 
    cgs_asb =cgs_sigma*4.d0/cgs_c 
    cgs_lsun = 3.839D33 
    cgs_ev = 1.602176463158D-12 
    cgs_hpk = 6.62606876e-27 
    cgs_yHe = 9.7222222d-2 
    Unit_rho = ModelParam_rho 
    Unit_v = 1.d5
    Unit_l = 9.36d-2 * cgs_pc
    Unit_t = Unit_l / Unit_v
    Unit_l3 = Unit_l**3.d0 
    Unit_m = Unit_rho * Unit_l**3 
    Unit_e = Unit_m / Unit_l / Unit_t**2 
    Unit_b = sqrt(Unit_e) 
    Unit_n = Unit_rho/(cgs_mu*cgs_amu) 
    Unit_kms = Unit_v/1.D5 
    Unit_au = Unit_l/cgs_au 
    Unit_pc = Unit_l/cgs_pc 
    Unit_yr = Unit_t/cgs_yr 
    Unit_msun = Unit_m/cgs_msun 
    Unit_acc = Unit_l / Unit_t**2.d0 
    Unit_erg = Unit_m * Unit_l**2 / Unit_t**2 
    Unit_gauss = Unit_b 
    Unit_ugauss = Unit_b/1.d-6 
    Unit_invmu = 1.d0/((1.d0 + 4.d0*cgs_yHe)*cgs_mh)
    Unit_rhomu = Unit_rho / ((1.d0 + 4.d0*cgs_yHe)*cgs_mh)
    Unit_lrhomu = Unit_l*Unit_rho / ((1.d0 + 4.d0*cgs_yHe)*cgs_mh)
    myrank = get_myrank()
    if (myrank==0) then
       print *, "Unit_rho=", Unit_rho,"g/cm^3"
       print *, "Unit_v=", Unit_v,"cm/s"
       print *, "Unit_t=", Unit_t/cgs_yr,"yr"
       print *, "Unit_l=", Unit_l/cgs_au,"au"
       print *, "Unit_invmu", Unit_invmu
    end if
  end subroutine unit_init
end module unit
