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
  real(kind=DBL_KIND),parameter :: ModelParam_rho = 1.44867883D-19 ! central density (in g cm^{-3}) for a filament with 0.1 pc full width
  real(kind=DBL_KIND),parameter :: ModelParam_temp = 10.d0  ! temperature (in k)
  real(kind=DBL_KIND),parameter :: ModelParam_gamma = 1.d0  ! polytropic index
  real(kind=DBL_KIND),save :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_msun
  public :: unit_init
  public :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_msun
  public :: ModelParam_rho, ModelParam_temp, ModelParam_gamma
contains
  !-------------------------------------------------------------------------
  ! module for physical units
  !-------------------------------------------------------------------------
  subroutine unit_init
    use mpilib
    use parameter
    implicit none
    real(kind=DBL_KIND) :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun

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

    ! =================================
    ! Units of scales
    !    physical unit / non-dimension
    ! =================================
    ! Example:
    !  x ....... length in non-dimension
    !  x_au .... length in physical unit, AU
    !    x   = x_au / Unit_au
    !    x_au = x * Unit_au

    Unit_rho = ModelParam_rho     ! density [g/cm^3]
    Unit_v   = sqrt(ModelParam_gamma * cgs_kb * ModelParam_temp / (cgs_mu *cgs_amu)) ! velocity [cm/s]
    Unit_t   = 1.d0 / sqrt( Pi4 * cgs_gc * Unit_rho) ! time [s]
    Unit_l   = Unit_v * Unit_t                       ! length [cm]
    Unit_m   = Unit_rho * Unit_l**3                  ! mass [g]
    Unit_e   = Unit_m / Unit_l / Unit_t**2 ! energy density [erg/cm^3]
    Unit_b   = sqrt(Unit_e)                ! magnetic field [Gauss]
    Unit_n   = Unit_rho/(cgs_mu*cgs_amu)    ! number density [cm^-3]
    Unit_kms = Unit_v/1.D5                 ! velocity [km/s]
    Unit_au  = Unit_l/cgs_au                   ! lenght [AU]
    Unit_pc  = Unit_l/cgs_pc                   ! length [pc]
    Unit_yr  = Unit_t/cgs_yr                   ! time  [yr]
    Unit_msun = Unit_m/cgs_msun                ! mass [M_sun]

  end subroutine unit_init
end module unit
