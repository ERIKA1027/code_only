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
  real(kind=DBL_KIND),save :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_msun, Unit_rsun, Unit_sec, Unit_hour, Unit_day, Unit_gauss, Unit_ugauss, Unit_tesla, Unit_utesla, Unit_ntesla
  real(kind=DBL_KIND),save :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun, cgs_hour, cgs_day

  public :: unit_init
  public :: Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, Unit_b, Unit_n, Unit_kms, Unit_au, Unit_pc, Unit_yr, Unit_msun, Unit_rsun, Unit_sec, Unit_hour, Unit_day, Unit_gauss, Unit_ugauss, Unit_tesla, Unit_utesla, Unit_ntesla
  public :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun, cgs_hour, cgs_day

contains
  !-------------------------------------------------------------------------
  ! module for physical units
  !-------------------------------------------------------------------------
  subroutine unit_init
    use parameter
    implicit none
!!$    real(kind=DBL_KIND) :: cgs_c, cgs_amu, cgs_mp, cgs_me, cgs_mh, cgs_mu, cgs_yr, cgs_gc, cgs_kb, cgs_pc, cgs_au, cgs_msun, cgs_rsun, cgs_hour, cgs_day

    ! ======================
    ! 物理定数の設定 in cgs
    ! ======================
    cgs_c = 2.99792458D10         ! 光速 (in cm/s)
    cgs_amu = 1.660538921D-24     ! 原子質量単位 (in g)
    cgs_mp = 1.672621777D-24      ! 陽子質量 (in g)
    cgs_me = 9.10938291D-28       ! 電子質量 (in g)
    cgs_mh = cgs_mp + cgs_me      ! 水素原子質量 (in g) 近似的に正しい。
    cgs_mu = 0.5d0                ! 平均分子量
    cgs_hour = 3600.d0            ! 1時間 (in s)
    cgs_day = 24.d0 * cgs_hour    ! 1日 (in s)
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

    Unit_n   = 1.d0       ! number density of proton [cm^-3]
    Unit_rho = cgs_mh * Unit_n     ! density [g/cm^3]
    Unit_l   = cgs_rsun   ! length [cm]
    Unit_t   = cgs_hour   ! time [s]
    Unit_v   = Unit_l/Unit_t ! velocity [cm/s]
    Unit_m   = Unit_rho * Unit_l**3                  ! mass [g]
    Unit_e   = Unit_m / Unit_l / Unit_t**2 ! energy density [erg/cm^3]
    Unit_b   = sqrt(Unit_e)                ! magnetic field [Gauss]
    Unit_kms = Unit_v/1.D5                 ! velocity [km/s]
    Unit_rsun = Unit_l/cgs_rsun            ! length [R_sun]
    Unit_au  = Unit_l/cgs_au               ! lenght [AU]
    Unit_pc  = Unit_l/cgs_pc               ! length [pc]
    Unit_sec = Unit_t                      ! time [sec]
    Unit_hour = Unit_t/cgs_hour            ! time [hour]
    Unit_day = Unit_t/cgs_day              ! time [day]
    Unit_yr  = Unit_t/cgs_yr               ! time  [yr]
    Unit_msun = Unit_m/cgs_msun            ! mass [M_sun]
    Unit_gauss = Unit_b                    ! magnetic field [Gauss]
    Unit_ugauss = Unit_b/1.d-6             ! magnetic field [micro Gauss]
    Unit_tesla = Unit_b / 1.D4             ! magnetic field [Tesla]
    Unit_utesla = Unit_tesla/1.D-6         ! magnetic field [micro Tesla]
    Unit_ntesla = Unit_tesla/1.D-9         ! magnetic field [nano Tesla]

  end subroutine unit_init
end module unit
