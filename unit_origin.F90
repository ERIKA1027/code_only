#include "config.h"
!-------------------------------------------------------------------------
! module for physical units
!-------------------------------------------------------------------------
module unit
  implicit none
  private
  ! ================
  ! model parameters
  ! ================
  real(kind=DBL_KIND),parameter :: ModelParam_rho = 1.D-18 ! central density (in g cm^{-3}) for turbulent molecular core with sink
  real(kind=DBL_KIND),parameter :: ModelParam_temp = 10000.d0  ! temperature (in k)
  real(kind=DBL_KIND),parameter :: ModelParam_gamma = 1.4d0  ! polytropic index
  real(kind=DBL_KIND),parameter :: MP_hnu_IR  = 0.01d0*1.602176463158D-12 ![erg] energy of invidial IR photons
!!$  real(kind=DBL_KIND),parameter :: ModelParam_gamma = 1.d0  ! polytropic index
  real(kind=DBL_KIND),save :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc &
    , Unit_yr, Unit_msun, Unit_gauss, Unit_ugauss
  real(kind=DBL_KIND),save :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc
  real(kind=DBL_KIND),save :: cgs_au, cgs_msun, cgs_rsun, cgs_sigma, cgs_lsun, cgs_ev, cgs_hpk, cgs_asb !KS MODIFIED
  real(kind=DBL_KIND),save :: Unit_acc, Unit_invmu, Unit_rhomu, Unit_lrhomu, Unit_l3, Unit_erg, cgs_yHe
  public :: unit_init
  public :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_gauss, Unit_ugauss
  public :: Unit_msun,Unit_acc, Unit_invmu,Unit_rhomu, Unit_lrhomu, Unit_l3
  public :: ModelParam_rho, ModelParam_temp, ModelParam_gamma
  public :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun, cgs_sigma, cgs_lsun, cgs_ev  !KS MODIFIED
  public :: cgs_hpk, cgs_asb, Unit_erg
  public :: cgs_yHe, MP_hnu_IR
contains
  !-------------------------------------------------------------------------
  ! module for physical units
  !-------------------------------------------------------------------------
  subroutine unit_init
    use mpilib
    use parameter
    implicit none
!    real(kind=DBL_KIND) :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun  !KS MODIFIED

    ! ======================
    ! 物理定数の設定 in cgs
    ! ======================
    cgs_c = 2.99792458D10         ! 光速 (in cm/s)
    cgs_amu = 1.660538921D-24     ! 原子質量単位 (in g)
    cgs_mp = 1.672621777D-24      ! 陽子質量 (in g)
    cgs_me = 9.10938291D-28       ! 電子質量 (in g)
    cgs_mh =  cgs_mp + cgs_me     ! 水素原子質量 (in g) 近似的に正しい。
    cgs_mu = 2.3                  ! 平均分子量
    cgs_yr = 3.1556926D7          ! 1年 (in s)
    cgs_gc = 6.67384D-8           ! 重力定数 (in dyn cm^2 g^-2)
    cgs_kb = 1.3806488D-16        ! ボルツマン定数 (in erg deg^{-1})
    cgs_au = 1.495978707D13       ! AU (in cm)
    cgs_pc = cgs_au * 180.d0 * 3600.d0 / Pi  ! pc (in cm)
    cgs_msun = 1.9884D33          ! 太陽質量 (in g)
    cgs_rsun = 6.960D10           ! 太陽半径 (in cm)
    cgs_sigma = 5.67051D-5        ! シュテファン=ボルツマン定数 (in erg cm^-2 K^-4 s^-1) (KS ADDED)
    cgs_asb =cgs_sigma*4.d0/cgs_c ! radiation constant
    cgs_lsun = 3.839D33           ! 太陽光度 (in erg s^-1) (KS ADDED)
    cgs_ev = 1.602176463158D-12   ! eV (in erg) (KS ADDED)
    cgs_hpk = 6.62606876e-27      ! plank constant
    cgs_yHe  = 9.7222222d-2       ! Helium abundance  

    ! =================================
    ! Units of scales
    !    physical unit / non-dimension
    ! =================================
    ! Example:
    !  x ....... length in non-dimension
    !  x_au .... length in physical unit, AU
    !    x   = x_au / Unit_au
    !    x_au = x * Unit_au

    Unit_rho    = ModelParam_rho     ! density [g/cm^3]
!    Unit_v      = sqrt(ModelParam_gamma * cgs_kb * ModelParam_temp / (cgs_mu *cgs_amu)) ! velocity [cm/s]
!    Unit_t      = 1.d0 / sqrt( Pi4 * cgs_gc * Unit_rho) ! time [s]
!    Unit_l      = Unit_v * Unit_t                       ! length [cm]
 
    !----------EO_added----------!
    Unit_v = 1.d5 
    ! Unit_l = 2.d0 * cgs_gc * 1.d4 * cgs_msun / (2.d1*Unit_v* 2.d1*Unit_v)
    Unit_l = 1.4d5 * cgs_au
    Unit_t = Unit_l / Unit_v
    !----------EO_added----------!

    Unit_l3     = Unit_l**3.d0                          ! volume [cm^3]
    Unit_m      = Unit_rho * Unit_l**3                  ! mass [g]
    Unit_e      = Unit_m / Unit_l / Unit_t**2 ! energy density [erg/cm^3]
    Unit_b      = sqrt(Unit_e)                ! magnetic field [Gauss]
    Unit_n      = Unit_rho/(cgs_mu*cgs_amu)    ! number density [cm^-3]
    Unit_kms    = Unit_v/1.D5                 ! velocity [km/s]
    Unit_au     = Unit_l/cgs_au                   ! lenght [AU]
    Unit_pc     = Unit_l/cgs_pc                   ! length [pc]
    Unit_yr     = Unit_t/cgs_yr                   ! time  [yr]
    Unit_msun   = Unit_m/cgs_msun                ! mass [M_sun]
    Unit_acc    = Unit_l / Unit_t**2.d0           ! acc  [cm s^{-2}]
    Unit_erg    = Unit_m * Unit_l**2 / Unit_t**2 ! [erg]
    Unit_gauss  = Unit_b                    ! magnetic field [Gauss]
    Unit_ugauss = Unit_b/1.d-6             ! magnetic field [micro Gauss]

    Unit_invmu  = 1.d0/((1.d0 + 4.d0*cgs_yHe)*cgs_mh)
    Unit_rhomu  = Unit_rho / ((1.d0 + 4.d0*cgs_yHe)*cgs_mh)
    Unit_lrhomu = Unit_l*Unit_rho / ((1.d0 + 4.d0*cgs_yHe)*cgs_mh)


    ! KS ADDED
    myrank = get_myrank()
    if (myrank==PRIMARY_RANK) then
       print *, "Unit_rho=", Unit_rho,"g/cm^3"
       print *, "Unit_v=", Unit_v,"cm/s"
       print *, "Unit_t=", Unit_t/cgs_yr,"yr"
       print *, "Unit_l=", Unit_l/cgs_au,"au"
       print *, "Unit_invmu", Unit_invmu
    end if
    

  end subroutine unit_init
end module unit
