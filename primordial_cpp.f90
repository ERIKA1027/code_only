module primordial
  use unit, only: cgs_amu
  use kinzoku
  use kinzoku2
    use modelParameter, only: MP_frac_C_solar, MP_frac_O_solar, MP_mu, MP_frac_COsum, MP_AC, MP_AO
  implicit none
  private
  real(kind=8), parameter :: boltz = 1.380662d-16 
  real(kind=8), parameter :: evolt = 1.6022d-12 
  real(kind=8), parameter :: yHe = 9.7222222d-2 
  real(kind=8), parameter :: zred = 0.d0 
  real(kind=8), parameter :: MP_Tcmb = 2.73 
  real(kind=8), parameter :: CONST_G = 6.6726e-8
  real(kind=8), parameter :: CONST_amu= 1.66053886e-24
  real(kind=8), parameter :: CONST_kB = 1.3806505e-16
  real(kind=8), parameter :: cgs_sigma = 5.67051D-5
  type chem_vals
    real(kind=8) :: nH 
    real(kind=8) :: Tg 
    real(kind=8) :: Td 
    real(kind=8),dimension(0:6 -1) :: ychem 
    real(kind=8) :: yco 
    real(kind=8) :: rHpi 
    real(kind=8) :: rH2pd 
    real(kind=8) :: rHmpd 
    real(kind=8) :: heat 
    real(kind=8) :: rgfuv 
    real(kind=8) :: rdph 
    real(kind=8) :: rcopd 
    real(kind=8) :: rOII 
    real(kind=8) :: xlmbdj 
    real(kind=8) :: dt 
    real(kind=8) :: fd 
    real(kind=8) :: EradIR 
    real(kind=8) :: chi_d 
    real(kind=8) :: metal 
    real(kind=8) :: xNcH 
    real(kind=8), dimension(0:2) :: dvdr 
  end type chem_vals
  type rad_others
    real(kind=8) :: xeuv 
    real(kind=8) :: xfuv 
    real(kind=8) :: alpha_euv 
    real(kind=8) :: heat_euv 
    real(kind=8) :: hhm 
    real(kind=8) :: lumeuv 
    real(kind=8) :: lumfuv 
    real(kind=8) :: sig_euv 
    real(kind=8) :: sig_fuv 
    real(kind=8) :: rOII 
  end type rad_others
  integer,save :: dbg_flg_prm = 0 
  public :: chemreact, HUVAHHH, yHe, zred, c_H2, c_H2_2, CoolSolverExplicit, CoolSolverImplicit &
    , adjust_abundance, chemreact_adptv
  public :: dbg_flg_prm, H2cool, H2cool_HM, H2cool_LTEfit, H2cool_LowD_G15, H2cool_Omukai
  public :: chem_vals, rad_others
  public :: get_xmu
  public :: ProstFit2
  external DGESV
contains
  subroutine react_rat_highspeed2(xk, xnH, y, r_f_tot)
    implicit none
    real(kind=8),intent(IN) :: xk(0:25 -1),xnH,y(0:6 -1)
    real(kind=8),intent(OUT) :: r_f_tot(0:6 -1)
    integer :: isp, ire
    real (kind=8) :: yhmn, rate, rrate
    r_f_tot(:) = 0.d0
    rate = xk(0)*y(0)*y(2)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(2) = r_f_tot(2) + rate
    r_f_tot(3) = r_f_tot(3) + rate
    rate = xk(1)*y(3)*y(2)*xnH
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(0) = r_f_tot(0) + rate
    rate = xk(2)*y(4)*y(0)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(1) = r_f_tot(1) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(3)*y(1)*y(3)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(5) = r_f_tot(5) + rate
    r_f_tot(0) = r_f_tot(0) + rate
    rate = xk(4)*y(1)*y(2)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(5)*y(1)*y(0)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(6)*y(0)*y(0)*y(0)*xnH*xnH
    r_f_tot(0) = r_f_tot(0) - 2.d0*rate
    r_f_tot(1) = r_f_tot(1) + rate
    rate = xk(7)*y(0)*y(0)*y(1)*xnH*xnH
    r_f_tot(0) = r_f_tot(0) - 2.d0*rate
    r_f_tot(1) = r_f_tot(1) + rate
    rate = xk(8)*y(1)*y(1)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(9)*y(0)*y(2)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(4) = r_f_tot(4) + rate
    rate = xk(10)*y(0)
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(11)*y(1)
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(3) = r_f_tot(3) + 2.d0*rate
    r_f_tot(2) = r_f_tot(2) + 2.d0*rate
    rate = xk(12)*y(1)
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(13)*y(0)*y(0)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(14)*y(4)*y(2)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(2) = r_f_tot(2) + rate
    r_f_tot(0) = r_f_tot(0) + rate
    rate = xk(15)*y(4)*y(3)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(5) = r_f_tot(5) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(16)*y(4)*y(3)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(17)*y(4)
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(0) = r_f_tot(0) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(18)*y(0)*y(3)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(5) = r_f_tot(5) + rate
    rate = xk(19)*y(5)*y(0)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(1) = r_f_tot(1) + rate
    r_f_tot(3) = r_f_tot(3) + rate
    rate = xk(20)*y(5)*y(2)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(0) = r_f_tot(0) + 2.d0*rate
    rate = xk(21)*y(5)*y(4)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(1) = r_f_tot(1) + rate
    r_f_tot(0) = r_f_tot(0) + rate
    rate = xk(22)*y(0)*xnH
    r_f_tot(0) = r_f_tot(0) - 2.d0*rate
    r_f_tot(1) = r_f_tot(1) + rate
    rate = xk(23)*y(0)
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) + rate
    r_f_tot(2) = r_f_tot(2) + rate
    rate = xk(24)*y(1)
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(5) = r_f_tot(5) + rate
    r_f_tot(2) = r_f_tot(2) + rate
  end subroutine react_rat_highspeed2
  subroutine react_drdy(xk, xnH, y, dr_fdy)
    real(kind=8),intent(IN) :: xk(0:25 -1),xnH,y(0:6 -1)
    real(kind=8),dimension(0:6 -1,0:6 -1),intent(OUT) :: dr_fdy
    real(kind=8) :: xnH_2, ddr_A, ddr_B, ddr_C
    dr_fdy(:,:) = 0.d0
    xnH_2 = xnH*xnH
    ddr_A = xk(0)*y(2)*xnH
    ddr_B = xk(0)*y(0)*xnH
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_A
    dr_fdy(2,0) = dr_fdy(2,0) + ddr_A
    dr_fdy(3,0) = dr_fdy(3,0) + ddr_A
    dr_fdy(0,2) = dr_fdy(0,2) - ddr_B
    dr_fdy(2,2) = dr_fdy(2,2) + ddr_B
    dr_fdy(3,2) = dr_fdy(3,2) + ddr_B
    ddr_A = xk(1)*y(2)*xnH
    ddr_B = xk(1)*y(3)*xnH
    dr_fdy(3,3) = dr_fdy(3,3) - ddr_A
    dr_fdy(2,3) = dr_fdy(2,3) - ddr_A
    dr_fdy(0,3) = dr_fdy(0,3) + ddr_A
    dr_fdy(3,2) = dr_fdy(3,2) - ddr_B
    dr_fdy(2,2) = dr_fdy(2,2) - ddr_B
    dr_fdy(0,2) = dr_fdy(0,2) + ddr_B
    ddr_A = xk(2)*y(0)*xnH
    ddr_B = xk(2)*y(4)*xnH
    dr_fdy(4,4) = dr_fdy(4,4) - ddr_A
    dr_fdy(0,4) = dr_fdy(0,4) - ddr_A
    dr_fdy(1,4) = dr_fdy(1,4) + ddr_A
    dr_fdy(2,4) = dr_fdy(2,4) + ddr_A
    dr_fdy(4,0) = dr_fdy(4,0) - ddr_B
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_B
    dr_fdy(1,0) = dr_fdy(1,0) + ddr_B
    dr_fdy(2,0) = dr_fdy(2,0) + ddr_B
    ddr_A = xk(3)*y(3)*xnH
    ddr_B = xk(3)*y(1)*xnH
    dr_fdy(1,1) = dr_fdy(1,1) - ddr_A
    dr_fdy(3,1) = dr_fdy(3,1) - ddr_A
    dr_fdy(5,1) = dr_fdy(5,1) + ddr_A
    dr_fdy(0,1) = dr_fdy(0,1) + ddr_A
    dr_fdy(1,3) = dr_fdy(1,3) - ddr_B
    dr_fdy(3,3) = dr_fdy(3,3) - ddr_B
    dr_fdy(5,3) = dr_fdy(5,3) + ddr_B
    dr_fdy(0,3) = dr_fdy(0,3) + ddr_B
    ddr_A = xk(4)*y(2)*xnH
    ddr_B = xk(4)*y(1)*xnH
    dr_fdy(1,1) = dr_fdy(1,1) - ddr_A
    dr_fdy(0,1) = dr_fdy(0,1) + 2.d0*ddr_A
    dr_fdy(1,2) = dr_fdy(1,2) - ddr_B
    dr_fdy(0,2) = dr_fdy(0,2) + 2.d0*ddr_B
    ddr_A = xk(5)*y(0)*xnH
    ddr_B = xk(5)*y(1)*xnH
    dr_fdy(1,1) = dr_fdy(1,1) - ddr_A
    dr_fdy(0,1) = dr_fdy(0,1) + 2.d0*ddr_A
    dr_fdy(1,0) = dr_fdy(1,0) - ddr_B
    dr_fdy(0,0) = dr_fdy(0,0) + 2.d0*ddr_B
    ddr_A = 3.d0*xk(6)*y(0)*y(0)*xnH_2
    dr_fdy(0,0) = dr_fdy(0,0) - 2.d0*ddr_A
    dr_fdy(1,0) = dr_fdy(1,0) + ddr_A
    ddr_A = 2.d0*xk(7)*y(1)*y(0)*xnH_2
    ddr_B = xk(7)*y(0)*y(0)*xnH_2
    dr_fdy(0,0) = dr_fdy(0,0) - 2.d0*ddr_A
    dr_fdy(1,0) = dr_fdy(1,0) + ddr_A
    dr_fdy(0,1) = dr_fdy(0,1) - 2.d0*ddr_B
    dr_fdy(1,1) = dr_fdy(1,1) + ddr_B
    ddr_A = 2.d0*xk(8)*y(1)*xnH
    dr_fdy(1,1) = dr_fdy(1,1) - ddr_A
    dr_fdy(0,1) = dr_fdy(0,1) + 2.d0*ddr_A
    ddr_A = xk(9)*y(2)*xnH
    ddr_B = xk(9)*y(0)*xnH
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_A
    dr_fdy(2,0) = dr_fdy(2,0) - ddr_A
    dr_fdy(4,0) = dr_fdy(4,0) + ddr_A
    dr_fdy(0,2) = dr_fdy(0,2) - ddr_B
    dr_fdy(2,2) = dr_fdy(2,2) - ddr_B
    dr_fdy(4,2) = dr_fdy(4,2) + ddr_B
    dr_fdy(0,0) = dr_fdy(0,0) - xk(10)
    dr_fdy(3,0) = dr_fdy(3,0) + xk(10)
    dr_fdy(2,0) = dr_fdy(2,0) + xk(10)
    dr_fdy(1,1) = dr_fdy(1,1) - xk(11)
    dr_fdy(3,1) = dr_fdy(3,1) + 2.d0*xk(11)
    dr_fdy(2,1) = dr_fdy(2,1) + 2.d0*xk(11)
    dr_fdy(1,1) = dr_fdy(1,1) - xk(12)
    dr_fdy(0,1) = dr_fdy(0,1) + 2.d0*xk(12)
    ddr_A = 2.d0*xk(13)*y(0)*xnH
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_A
    dr_fdy(3,0) = dr_fdy(3,0) + ddr_A
    dr_fdy(2,0) = dr_fdy(2,0) + ddr_A
    ddr_A = xk(14)*y(2)*xnH
    ddr_B = xk(14)*y(4)*xnH
    dr_fdy(4,4) = dr_fdy(4,4) - ddr_A
    dr_fdy(2,4) = dr_fdy(2,4) + ddr_A
    dr_fdy(0,4) = dr_fdy(0,4) + ddr_A
    dr_fdy(4,2) = dr_fdy(4,2) - ddr_B
    dr_fdy(2,2) = dr_fdy(2,2) + ddr_B
    dr_fdy(0,2) = dr_fdy(0,2) + ddr_B
    ddr_A = xk(15)*y(3)*xnH
    ddr_B = xk(15)*y(4)*xnH
    dr_fdy(4,4) = dr_fdy(4,4) - ddr_A
    dr_fdy(3,4) = dr_fdy(3,4) - ddr_A
    dr_fdy(5,4) = dr_fdy(5,4) + ddr_A
    dr_fdy(2,4) = dr_fdy(2,4) + ddr_A
    dr_fdy(4,3) = dr_fdy(4,3) - ddr_B
    dr_fdy(3,3) = dr_fdy(3,3) - ddr_B
    dr_fdy(5,3) = dr_fdy(5,3) + ddr_B
    dr_fdy(2,3) = dr_fdy(2,3) + ddr_B
    ddr_A = xk(16)*y(3)*xnH
    ddr_B = xk(16)*y(4)*xnH
    dr_fdy(4,4) = dr_fdy(4,4) - ddr_A
    dr_fdy(3,4) = dr_fdy(3,4) - ddr_A
    dr_fdy(0,4) = dr_fdy(0,4) + 2.d0*ddr_A
    dr_fdy(4,3) = dr_fdy(4,3) - ddr_B
    dr_fdy(3,3) = dr_fdy(3,3) - ddr_B
    dr_fdy(0,3) = dr_fdy(0,3) + 2.d0*ddr_B
    dr_fdy(4,4) = dr_fdy(4,4) - xk(17)
    dr_fdy(0,4) = dr_fdy(0,4) + xk(17)
    dr_fdy(2,4) = dr_fdy(2,4) + xk(17)
    ddr_A = xk(18)*y(3)*xnH
    ddr_B = xk(18)*y(0)*xnH
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_A
    dr_fdy(3,0) = dr_fdy(3,0) - ddr_A
    dr_fdy(5,0) = dr_fdy(5,0) + ddr_A
    dr_fdy(0,3) = dr_fdy(0,3) - ddr_B
    dr_fdy(3,3) = dr_fdy(3,3) - ddr_B
    dr_fdy(5,3) = dr_fdy(5,3) + ddr_B
    ddr_A = xk(19)*y(0)*xnH
    ddr_B = xk(19)*y(5)*xnH
    dr_fdy(5,5) = dr_fdy(5,5) - ddr_A
    dr_fdy(0,5) = dr_fdy(0,5) - ddr_A
    dr_fdy(1,5) = dr_fdy(1,5) + ddr_A
    dr_fdy(3,5) = dr_fdy(3,5) + ddr_A
    dr_fdy(5,0) = dr_fdy(5,0) - ddr_B
    dr_fdy(0,0) = dr_fdy(0,0) - ddr_B
    dr_fdy(1,0) = dr_fdy(1,0) + ddr_B
    dr_fdy(3,0) = dr_fdy(3,0) + ddr_B
    ddr_A = xk(20)*y(2)*xnH
    ddr_B = xk(20)*y(5)*xnH
    dr_fdy(5,5) = dr_fdy(5,5) - ddr_A
    dr_fdy(2,5) = dr_fdy(2,5) - ddr_A
    dr_fdy(0,5) = dr_fdy(0,5) + 2.d0*ddr_A
    dr_fdy(5,2) = dr_fdy(5,2) - ddr_B
    dr_fdy(2,2) = dr_fdy(2,2) - ddr_B
    dr_fdy(0,2) = dr_fdy(0,2) + 2.d0*ddr_B
    ddr_A = xk(21)*y(4)*xnH
    ddr_B = xk(21)*y(5)*xnH
    dr_fdy(5,5) = dr_fdy(5,5) - ddr_A
    dr_fdy(4,5) = dr_fdy(4,5) - ddr_A
    dr_fdy(1,5) = dr_fdy(1,5) + ddr_A
    dr_fdy(0,5) = dr_fdy(0,5) + ddr_A
    dr_fdy(5,4) = dr_fdy(5,4) - ddr_B
    dr_fdy(4,4) = dr_fdy(4,4) - ddr_B
    dr_fdy(1,4) = dr_fdy(1,4) + ddr_B
    dr_fdy(0,4) = dr_fdy(0,4) + ddr_B
    dr_fdy(0,0) = dr_fdy(0,0) - 2.d0*xk(22)*xnH
    dr_fdy(1,0) = dr_fdy(1,0) + xk(22)*xnH
    dr_fdy(0,0) = dr_fdy(0,0) - xk(23)
    dr_fdy(3,0) = dr_fdy(3,0) + xk(23)
    dr_fdy(2,0) = dr_fdy(2,0) + xk(23)
    dr_fdy(1,1) = dr_fdy(1,1) - xk(24)
    dr_fdy(5,1) = dr_fdy(5,1) + xk(24)
    dr_fdy(2,1) = dr_fdy(2,1) + xk(24)
  end subroutine react_drdy
  subroutine adjust_abundance(y, yco, metal)
    implicit none
    real(kind=8),intent(IN) :: metal
    real(kind=8),intent(INOUT) :: y(0:6 -1), yco
    real(kind=8) :: ycII
    real(kind=8) :: yHtot, yCtot, yOtot, yDtot, chrgtot, nchrg, pchrg
    real(kind=8), parameter :: min_value = 1d-20 
    real(kind=8) :: frac_C, frac_O
    frac_C = metal*MP_frac_C_solar
    frac_O = metal*MP_frac_O_solar
    y(0) = MAX(MIN(y(0), 1.d0), min_value)
    y(1) = MAX(MIN(y(1), 0.5d0), min_value)
    y(3) = MAX(MIN(y(3), 1.d0), min_value)
    y(4) = MAX(MIN(y(4), 1.d0), min_value)
    y(5) = MAX(MIN(y(5), 1.d0), min_value)
    y(2) = MAX(min_value, y(2))
    yco = MAX(MIN(yco, frac_C), min_value)
    ycII= MAX(frac_C-yco, min_value)
    yHtot = y(0)+2.d0*y(1)+y(3)+y(4)+2.d0*y(5)
    if(abs(yHtot-1) > 1.d-10) then
        if(yHtot < 1.d0) then
            y(0) = y(0)+1.d0-yHtot
        else
           y(0)=y(0)/yHtot
           y(1)=y(1)/yHtot
           y(3)=y(3)/yHtot
           y(4)=y(4)/yHtot
           y(5)=y(5)/yHtot
        endif
    endif
    pchrg = ycII +y(3)+y(5)
    nchrg = y(2)+y(4)
    chrgtot = 2 * (pchrg-nchrg) / (nchrg+pchrg)
    if(abs(chrgtot) > 1.d-10) then
        y(2) = MAX(min_value, y(2)+pchrg-nchrg)
    endif
  end subroutine adjust_abundance
  subroutine setting_ymax_ymin(y_max, y_min, metal)
      implicit none
      real(kind=8),dimension(0:6 -1) :: y_max, y_min
      real(kind=8), parameter :: min_value = 1d-20 
      real(kind=8) :: frac_C, frac_O, metal
      frac_C = metal*MP_frac_C_solar
      frac_O = metal*MP_frac_O_solar
      y_min(:) = min_value
      y_max(1) = 0.5d0
      y_max(2) = 3.d0 
      y_max(0) = 1.d0
      y_max(3) = 1.d0
      y_max(4) = 1.d0
      y_max(5) = 1.d0
  end subroutine setting_ymax_ymin
  subroutine chemreact(xch,tchem,xk_in)
    type(chem_vals) :: xch
    real(kind=8) :: y(0:6 -1)
    real(kind=8),intent(OUT) :: tchem
    real(kind=8),dimension(0:25 -1),intent(IN),optional :: xk_in 
    real(kind=8),dimension(0:6 -1) :: y_init,y_tmp,dy,ddy
    real(kind=8),dimension(0:25 -1) :: xk
    real(kind=8),dimension(0:6 -1) :: r_f,r_f_fw,r_f_bw
    real(kind=8),dimension(0:6 -1,0:6 -1) :: dr_fdy,A
    real(kind=8),dimension(0:6 -1,0:6 -1) :: dr_fdy_dbg 
    real(kind=8),dimension(0:6 -1) :: y_max, y_min
    integer :: isp,jsp,ia,info,itr,iguess
    real(kind=8) :: delta_y,dr_f,err,err_max,r_f_big,tch
    real(kind=8),parameter :: eps=1.d-5 
    real(kind=8),parameter :: eps_y=1.d-10 
    real(kind=8),parameter :: eps_conv = 1.d-5 
    integer, parameter :: itrmax = 100 
    integer, parameter :: itr_change1 = 20
    integer, parameter :: itr_change2 = 40
    integer,parameter :: inc=1, n=6
    integer :: ipiv(0:6 -1)
    logical :: nan_flg
    integer, dimension(0:6 -1) :: indx
    real(kind=8) :: d
    y(:) = xch%ychem(:)
    y_init(:) = y(:)
    if (present(xk_in)) then
       xk(:) = xk_in(:)
    else
       call react_coef(xch, xk)
    end if
    call setting_ymax_ymin(y_max, y_min, xch%metal)
    do iguess=0,2
       if (iguess==0) then
          y(:) = y_init(:) 
       else if (iguess==1) then
          y(:) = 1d-12 
          y(0) = 1d0
          call adjust_abundance(y &
            , xch%yco &
            , xch%metal)
       else if (iguess==2) then
          y(:) = 1d-12 
          y(2) = 1d0
          y(3) = 1d0
          call adjust_abundance(y &
            , xch%yco &
            ,xch%metal)
       end if
       dy(:) = y(:) - y_init(:)
       nan_flg = .false.
       do itr=0,itrmax-1
          call react_rat_highspeed2(xk,xch%nH,y,r_f)
          if (dbg_flg_prm == 1) then
             print '(/,A,I0,A,(1P10E15.7))', "itr: ",itr, ", xnH T_K rHpi rH2pd rHmpd dt: ",&
                  xch%nH,xch%Tg,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
             print '(A,(1P10E15.7))', "y_init: ",y_init(:)
             print '(A,(1P10E15.7))', "y:",y(:)
             print '(A,(1P10E15.7))', "r_f:",r_f(:)
          end if
          call react_drdy(xk, xch%nH, y, dr_fdy)
          do jsp=0,6 -1
            do isp=0,6 -1
                   A(isp,jsp) = -xch%dt*dr_fdy(isp,jsp)
            enddo
          enddo
          do jsp=0,6 -1
            A(jsp,jsp) = A(jsp,jsp) + 1.d0
          enddo
          do isp=0,6 -1
             ddy(isp) = r_f(isp)*xch%dt-dy(isp)
          enddo
          call ludcmp(A, indx, d)
          call lubksb(A, indx, ddy)
          if(itr < itr_change1) then
            do isp=0,6 -1
               dy(isp) = dy(isp) + ddy(isp)
               y(isp) = y(isp) + ddy(isp) 
            enddo
          else
            if(itr < itr_change2) then
              do isp=0,6 -1
                 dy(isp) = dy(isp) + ddy(isp)*0.5d0
                 y(isp) = y(isp) + ddy(isp)*0.5d0 
              enddo
            else
              do isp=0,6 -1
                if(ddy(isp) < 0.d0) then
                  y_max(isp) = y(isp)
                  ddy(isp) = max((y_max(isp)+y_min(isp))*0.5d0-y(isp), ddy(isp))
                  dy(isp) = dy(isp) + ddy(isp)
                  y(isp) = y(isp) + ddy(isp) 
                else
                  y_min(isp) = y(isp)
                  ddy(isp) = min((y_max(isp)+y_min(isp))*0.5d0-y(isp), ddy(isp))
                  dy(isp) = dy(isp) + ddy(isp)
                  y(isp) = y(isp) + ddy(isp) 
                endif
              enddo
            endif
          endif
          if (dbg_flg_prm == 1) then
             print '(A,(1P10E15.7))', "(Before_adjust) y:",y(:) 
          end if
          if(itr < itr_change2) then
            call adjust_abundance(y &
            , xch%yco &
            , xch%metal)
          endif
          dy(:) = y(:) - y_init(:) 
          if (dbg_flg_prm == 1) then
             print '(A,(1P10E15.7))', "(After_adjust) y:",y(:) 
          end if
          do isp=0,6 -1
             if (y(isp)/=y(isp) .or. y(isp)*0/=0 .or. y(isp)==Huge(1d0)) then 
                nan_flg = .true.
                print *, "(chemreact) NaN/Inf/Huge found, move to next initial guess"
                exit
             end if
          end do
          if (nan_flg) exit 
          err_max = 0.d0
          do isp=0,6 -1
             if (y(isp) == 0d0) then
                err_max=Huge(1d0) 
                exit
             end if
             if(y(isp) > 1d-20 .or. y_init(isp) > 1d-20) then 
                err = dabs(ddy(isp)/y(isp))
                err_max = max(err,err_max)
             end if
         enddo
         if (dbg_flg_prm == 1) then
             print *,"itr, (error_y, error_max) y:", itr,err_max
             do isp=0,6 -1
              print *, isp, y(isp), ddy(isp)/y(isp), y_max(isp), y_min(isp)
             enddo
         end if
          if(err_max.lt.eps_conv) then
             call adjust_abundance(y &
               , xch%yco &
               , xch%metal)
             tchem = HUGE(1d0)
             do isp=0,6 -1
                if (y(isp) /= y_init(isp)) then
                   if ((y(isp)+y_init(isp))/2. < 1d-5) then 
                      cycle
                   endif
                   tch = dabs((y(isp) + y_init(isp)) &
                        /(2.d0*(y(isp) - y_init(isp))))*xch%dt
                   tchem = min(tch,tchem)
                end if
             enddo
             xch%ychem(:) = y(:)
             return
          endif
          if(itr.eq.itrmax-1 .and. dbg_flg_prm == 1) then
             if (iguess==0) then
                print *, "chemreact: iteration not converged (initial guess of previous abundance)"
                print *, "try initial guess of fully dissociated gas"
             else if (iguess==1) then
                print *, "chemreact: iteration not converged (initial guess of fully dissociated gas)"
                print *, "try initial guess of fully ionized gas"
             else if (iguess==2) then
                print *, "chemreact: iteration not converged (initial guess of fully ionized gas)"
                print *, "try another timestep update"
             end if
          end if
       end do
    end do 
    tchem = -1d0 
  end subroutine chemreact
  subroutine chemreact_adptv(xch,tchem,force_substep,xk_in)
    type(chem_vals) :: xch
    type(chem_vals) :: xch_tmp
    real(kind=8),intent(OUT) :: tchem
    logical,intent(IN),optional :: force_substep
    real(kind=8),dimension(0:25 -1),intent(IN),optional :: xk_in 
    real(kind=8) :: t_sub
    integer :: i,n
    integer, parameter :: imax = 10 
    integer, parameter :: nmax = 20 
    logical :: bool_force_substep
    real(kind=8) :: y_o(0:6 -1), y_dbg(0:20,0:6 -1), dt_sub_n(0:20)
    integer :: nn, i_itr_n(0:20)
    xch_tmp = xch
    y_o(:) = xch%ychem(:)
    y_dbg(:,:) = -1d10
    dt_sub_n(:) = -1d10
    i_itr_n(:) = -100
    bool_force_substep = .false.
    if ( present( force_substep) ) then
       bool_force_substep = force_substep
    endif
    t_sub=0d0 
    do n=0,nmax
       do i=0,imax
          xch_tmp%ychem(:) = xch%ychem(:)
          xch_tmp%dt = (xch%dt-t_sub)*10.d0**(-i) 
          if (n==0 .and. i==0 .and. bool_force_substep) then
             xch_tmp%dt = 0.5d0 * xch_tmp%dt 
          end if
          if (present(xk_in)) then
             call chemreact(xch_tmp,tchem,xk_in=xk_in) 
          else
             call chemreact(xch_tmp,tchem)
          end if
          i_itr_n(n) = i
          if (tchem>0d0) exit 
          if (i==imax) then
             print *, "(chemreact_adptv), i == imax ", i, imax
             print '(2(A,1P1E15.7))', "(chemreact_adptv), chemreact not converged w/ dt_sub=", xch_tmp%dt, ", dt=", xch%dt
             print '(/,A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt_sub: ",&
                  xch%nH,xch%Tg,xch%rHpi,xch%rH2pd,xch%rHmpd,xch_tmp%dt
             print '(A,(1P10E15.7))', "y_o: ",y_o(:)
             print *, "stopping..."
             stop
          end if
       end do 
       xch%ychem(:) = xch_tmp%ychem(:)
       t_sub = t_sub + xch_tmp%dt
       y_dbg(n,:) = xch_tmp%ychem(:)
       dt_sub_n(n) = xch_tmp%dt
       if ((xch%dt-t_sub) < 1d-3*xch%dt) exit
       if (n==nmax) then
          print *, "(chemreact_adptv), n == nmax ", n, nmax
          print '(3(A,1P1E15.7))', "(chemreact_adptv), max itr. reaches before reaching dt w/ t_sub=",&
               t_sub, ", dt=", xch%dt, ", dt_sub=", xch_tmp%dt
          print '(/,A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt_sub: ",&
               xch%nH,xch%Tg,xch%rHpi,xch%rH2pd,xch%rHmpd, xch_tmp%dt
          print '(A,(1P10E15.7))', "y_o: ",y_o(:)
          do nn=0, nmax
             print *, "n, i_itr_n, dt_sub_n, y_n: ",nn,i_itr_n(nn),dt_sub_n(nn), y_dbg(nn,:)
          end do
          print *, "stopping..."
          stop
       end if
    end do 
  end subroutine chemreact_adptv
  subroutine react_coef(xch, xk)
    type(chem_vals) :: xch
    real(kind=8),intent(OUT) :: xk(0:25 -1)
    integer :: i
    real(kind=8) :: T_eV, xlnT_eV, log_T, log_T2, log_T3,&
       xk_L, xk_H, rlgT4, rn_cr, a, xkcid, xkdt, T300
    real(kind=8) :: f_a, r_a, inv_T, sqrt_T, exp_Av
    T_eV = 8.61735d-5*xch%Tg

    xlnT_eV = dlog(T_eV)

    log_T = dlog10(xch%Tg)

    log_T2 = log_T*log_T

    log_T3 = log_T2*log_T

    inv_T = 1.d0/xch%Tg
    sqrt_T = dsqrt(xch%Tg)
    T300 = xch%Tg/300.d0
    xk(0) = dexp(-32.71396786 + (13.536556 &
         + (-5.73932875 + (1.56315498 + (-0.2877056 &
         + (3.48255977d-2 + (-2.63197617d-3 &
         + (1.11954395d-4 - 2.03914985d-6*xlnT_eV)*xlnT_eV) &
         *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
      xk(1) = 1.269d-13 * (315614.d0*inv_T)**1.503 *(1+(604625.d0*inv_T)**0.470)**(-1.923) 
    xk(2) = 1.35d-9* (xch%Tg**0.098493 + 0.32852*xch%Tg**0.5561 + 2.771d-7*xch%Tg**2.1826) &
         /(1.0 + 6.191d-3*xch%Tg**1.0461 + 8.9712d-11*xch%Tg**3.0424 + 3.2576d-14*xch%Tg**3.7741)
    if(xch%Tg <= 1.e4) then
       xk(3) = 3.d-10*dexp(-2.1050d4*inv_T)

    else
       xk(3) = 1.5d-10*dexp(-1.4000d4*inv_T)

    end if
    xk(4) = 4.4d-10*(xch%Tg**3.5d-1)*dexp(-1.02000d5*inv_T)
    xkcid = xkcidM96(xch%nH, xch%Tg, log_T, log_T2, log_T3)
    xkdt = xkdtM96 (xch%nH, xch%Tg, log_T, log_T2, log_T3)
    xk(5) = xkcid + xkdt
    xk(6) = 6.d-32*xch%Tg**(-0.25) + 2.d-31*xch%Tg**(-0.5)
    xk(7) = xk(6)*0.125d0
    xk_L = 1.18d-10*dexp(-6.95d+4*inv_T)
    xk_H = 8.125d-8/sqrt_T*dexp(-5.2e+4*inv_T)*(1.0 - dexp(-6.d3*inv_T))
    rlgT4 = log_T - 4.d0 
    rn_cr = 1.d1**(4.845d0 - 1.3d0*rlgT4 + 1.62d0*rlgT4**2)
    a = 1.d0/(1.d0 + xch%nH/rn_cr)
    if(a.eq.1.d0) then
       xk(8) = xk_L
    elseif(a.eq.0.d0) then
       xk(8) = xk_H
    else
       xk(8) = xk_H**(1.d0-a)*xk_L**a
    endif
    xk(9) = 1.4d-18*(xch%Tg**9.28d-1)*dexp(-xch%Tg/1.62d4)
    xk(10) = xch%rHpi
    xk(11) = xch%rHpi
    xk(12) = xch%rH2pd
    xk(13) = 1.7d-4*xk(0)
    xk(14) = dexp(- 18.01849334d0 + (2.3608522d0 &
         + (-0.28274430d0+(1.62331664d-2+(-3.36501203d-2 &
         + (1.17832978d-2+(-1.65619470d-3+(1.06827520d-4 &
         - 2.63128581d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV) &
         *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
    if(xch%Tg.le.8.d3) then
       xk(15) = 6.9d-9/(xch%Tg**3.5d-1)
    else
       xk(15) = 9.6d-7/(xch%Tg**9.d-1)
    endif
    xk(16) = 6.3d-8 + 5.7d-6/sqrt_T - 9.2d-11*sqrt_T + 4.4d-13*xch%Tg
    xk(17) = xch%rHmpd
    xk(18) = 1.d1**(-19.38d0 - 1.523d0*log_T &
         + 1.118d0*log_T**2 - 0.1269d0*log_T**3)
    xk(19) = 6.4d-10

    xk(20) = 2.d-7/sqrt_T

    xk(21) = 39.8371685741d-7/sqrt_T 
    f_a=1.d0/(1.d0+dexp(7.5d2*(1.d0/75.d0-1.d0/xch%Td)))
    r_a=0.3464101615d-17*sqrt_T*f_a/(1.d0+4.0d-2*dsqrt(xch%Tg+xch%Td)+2.0d-3*xch%Tg+ &
      8.0d-6*xch%Tg*xch%Tg)*xch%fd 
    xk(22) = r_a
   exp_Av = dexp(-xch%xNcH/(4.3d25)) 
   xk(23) = 1.5d0*2.d-16*exp_Av
      xk(24) = 2.d0*xk(23)
  end subroutine react_coef
  function xkcidM96(xnH, T, logT, logT2, logT3)
    real(kind=8) :: xkcidM96
    real(kind=8),intent(IN) :: xnH, T, logT, logT2, logT3
    real(kind=8),dimension(21) :: a
    real(kind=8) :: log_g_h1, log_g_h2, log_g_l1, log_g_l2,&
         log_n_c1, n_c1, log_n_c2, n_c2, p, log_gamma_cd
    data a /-1.784239e2, -6.842243e1, 4.320243e1, -4.633167, 6.970086e1,&
         4.087038e4, -2.370570e4, 1.288953e2, -5.391334e1, 5.315517,&
         -1.973427e1, 1.678095e4, -2.578611e4, 1.482123e1, -4.890915,& 
         4.749030e-1, -1.338283e2, -1.164408, 8.227443e-1, 5.864073e-1, -2.056313/
  log_g_h1 = a(1) + a(2)*logT + a(3)*logT2 + a(4)*logT3 + a(5)*log10(1+a(6)/T)
  log_g_h2 = a(7)/T
  log_g_l1 = a(8) + a(9)*logT + a(10)*logT2 + a(11)*log10(1+a(12)/T)
  log_g_l2 = a(13)/T
  log_n_c1 = a(14) + a(15)*logT + a(16)*logT2 + a(17)/T
  n_c1 = 10.0**log_n_c1
  log_n_c2 = a(18) + log_n_c1
  n_c2 = 10.0**log_n_c2
  p = a(19) + a(20)*exp(-T/1850.0) + a(21)*exp(-T/440.0)
  log_gamma_cd = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+(xnH/n_c1)**p) &
                 + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+(xnH/n_c2)**p)
  xkcidM96 = 10.0**log_gamma_cd
  end function xkcidM96
  function xkdtM96(xnH, T, logT, logT2, logT3)
    real(kind=8) :: xkdtM96
    real(kind=8),intent(IN) :: xnH, T, logT, logT2, logT3
    real(kind=8),dimension(21) :: a
    real(kind=8) :: log_g_h1, log_g_h2, log_g_l1, log_g_l2,&
         log_n_c1, n_c1, log_n_c2, n_c2, p, log_gamma_dt
    data a /-1.427664e2, 4.270741e1, -2.027365, -2.582097e-1,&
         2.136094e1, 2.753531e4, -2.146779e4, 6.034928e1, -2.743096e1,&
         2.676150, -1.128215e1, 1.425455e4, -2.312520e4, 9.305564, &
         -2.464009, 1.985955e-1, 7.430600e2, -1.174242, 7.502286e-1,&
         2.358848e-1, 2.937507/
  log_g_h1 = a(1) + a(2)*logT + a(3)*logT2 + a(4)*logT3 + a(5)*log10(1+a(6)/T)
  log_g_h2 = a(7)/T
  log_g_l1 = a(8) + a(9)*logT + a(10)*logT2 + a(11)*log10(1+a(12)/T)
  log_g_l2 = a(13)/T
  log_n_c1 = a(14) + a(15)*logT + a(16)*logT2 + a(17)/T
  if(log_n_c1 > 99.d0) then
    n_c1 = 1.d99
  else
    n_c1 = 10.0**log_n_c1
  endif
  log_n_c2 = a(18) + log_n_c1
  if(log_n_c2 > 99.d0) then
    n_c2 = 1.d99
  else
    n_c2 = 10.0**log_n_c2
  endif
  p = a(19) + a(20)*exp(-T/1850.0) + a(21)*exp(-T/440.0)
  log_gamma_dt = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+(xnH/n_c1)**p) &
                 + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+(xnH/n_c2)**p)
  xkdtM96 = 10.0**log_gamma_dt
  end function xkdtM96
  subroutine tot_cool(xch,xLmbd_tot,xk_in)
    type(chem_vals) :: xch
    real(kind=8),intent(OUT) :: xLmbd_tot 
    real(kind=8),dimension(0:25 -1),intent(IN),optional :: xk_in 
    real(kind=8) :: xkd, tau_cnt, esc_cont
    real(kind=8) :: xk(0:25 -1), chi_d
    real(kind=8) :: xLmbd_chem,xLmbd_line, dust_cool, ph_heat, ne
    real(kind=8) :: inside_exp
    logical :: isNotFinite
    real(kind=8) :: cr_heat, qcr_H, qcr_H2, lognH
    if (present(xk_in)) then
       xk(:) = xk_in(:)
    else
       call react_coef(xch, xk)
    end if
    call chem_cool(xch, xk, xLmbd_chem)
    call line_cool(xch, xLmbd_line)
    xLmbd_tot = xLmbd_chem + xLmbd_line
      call dtemp_radtr(xch%nH, xch%Tg, xch%Td, xch%EradIR, chi_d &
        ,dph_in=xch%rdph)
    dust_cool = 5.83e-8*xch%nH*(MP_mu*cgs_amu)*xch%nH*sqrt(xch%Tg*1.d-3) &
        *(1.d0-0.8d0*exp(-75.d0/xch%Tg))*(xch%Tg-xch%Td)*xch%fd 
    xLmbd_tot = xLmbd_tot + dust_cool
    ne = xch%ychem(2) * xch%nH
    call PhotoelectricHeating(ne, xch%nH, xch%Tg, xch%rgfuv, ph_heat, xch%fd)
    xLmbd_tot = xLmbd_tot - ph_heat 
    lognH = log10(xch%nH)
    qcr_H = (6.5d0+26.4d0*sqrt(xch%ychem(2)/(xch%ychem(2)+0.07d0)))*evolt 
    if(lognH < 2.d0) then
      qcr_H2 = 10.d0*evolt
    else if(lognH < 4.d0) then
      qcr_H2 = (10.d0+3.d0*(lognH-2.d0)*0.5d0)*evolt
    else if(lognH < 7.d0) then
      qcr_H2 = (13.d0+4.d0*(lognH-4.d0)/3.d0) *evolt
    else if(lognH < 10.d0) then
      qcr_H2 = (17.d0+(lognH - 7.d0)/3.d0) *evolt
    else
      qcr_H2 = 18.d0*evolt
    endif
    cr_heat = (xk(23)*qcr_H*xch%ychem(0)+xk(24)*qcr_H2*xch%ychem(1))*xch%nH 
    xLmbd_tot = xLmbd_tot - cr_heat
    if (dbg_flg_prm == 1) then
       print '(A,(1P7E15.7))', "xLmbd_chem, xLmbd_line, xLmbd_tot, dust_cool, ph_heat: " &
         ,xch%Tg, xLmbd_chem, xLmbd_line, xLmbd_tot, dust_cool, ph_heat, xch%rgfuv 
       print '(A,/)', "#-----------------------------------------------------------#"
    end if
  end subroutine tot_cool
    subroutine chem_cool(xch, xk, xLmbd_chem)
    type(chem_vals) :: xch
    real(kind=8),intent(IN) :: xk(0:25 -1)
    real(kind=8),intent(OUT) :: xLmbd_chem 
    real(kind=8) :: xn_cr,crit,rtHm,rt3b,rtdis,rtHci,rtHra,rtHpi,xL_H2,xL_ci,xL_rec,xL_HM,xG_pi
    real(kind=8) :: heat_H2disso, rtH2d
    real(kind=8) :: log_T,log_T2,log_T3,log_T4,log_T5
    real(kind=8),dimension(0:6 -1) :: r_f 
    real(kind=8) :: xL_H2d
    log_T = log10(xch%Tg)
    log_T2 = log_T*log_T
    log_T3 = log_T2*log_T
    log_T4 = log_T3*log_T
    log_T5 = log_T4*log_T
    xn_cr = 1.d6/dsqrt(xch%Tg) &
         /(1.6d0*xch%ychem(0)*dexp(-(4.d2/xch%Tg)**2) + 1.4d0*xch%ychem(1)*dexp(-1.2d4/(xch%Tg+1.2d3)))
    crit = 1.d0/(1.d0 + xn_cr/xch%nH)
    rtHm = xk(2)*xch%ychem(4)*xch%ychem(0) 
    rt3b = (xk(6)*(xch%ychem(0)**3) + xk(7)*(xch%ychem(0)**2)*xch%ychem(1))*xch%nH 
    rtdis = xk(5)*xch%ychem(1)*xch%ychem(0) + xk(8)*(xch%ychem(1)**2) 
    xL_H2 = - (3.73*rtHm*crit + 4.48*(rt3b*crit-rtdis)) * xch%nH*xch%nH*evolt 
    rtH2d = xk(22)*xch%ychem(0)
    xL_H2d = -(0.2d0+4.2d0*crit)*rtH2d*xch%nH*xch%nH*evolt
    rtHci = xk(0)*xch%ychem(0)*xch%ychem(2)
    xL_ci = rtHci*13.6*evolt*xch%nH*xch%nH
    xL_rec = 1d1 ** (-26.02 + 9.187d-1*log_T - 3.733d-1*log_T2&
         + 1.174d-1*log_T3 - 1.707d-2*log_T4 + 8.119d-4*log_T5)*xch%ychem(3)*xch%ychem(2)*xch%nH*xch%nH 
    rtHra = xk(9)*xch%ychem(0)*xch%ychem(2)
    xL_HM = boltz*xch%Tg*rtHra*xch%nH*xch%nH 
    rtHpi = xk(10)*xch%ychem(0)+xk(11)*xch%ychem(1)
    xG_pi = rtHpi*xch%heat*xch%nH
    rtH2d = xk(12)*xch%ychem(1) 
    heat_H2disso = 0.4d0*evolt*rtH2d*xch%nH + 9.d0*rtH2d*(2.2d0*evolt*crit)*xch%nH 
    xLmbd_chem = xL_H2 + xL_ci + xL_rec + xL_HM - xG_pi - heat_H2disso
    xLmbd_chem = xLmbd_chem + xL_H2d
    end subroutine chem_cool
  subroutine line_cool(xch, xLmbd_line)
    type(chem_vals) :: xch
    real(kind=8),intent(OUT) :: xLmbd_line 
    real(kind=8) :: xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, y_H2, y_H, y_e, y_Hp, T_5, tau_cnt
    real(kind=8) :: gff1, xLmbd_ff, xNc_H2, cooling_CO, heating_CO, xLmbd_CO
    real(kind=8) :: y_CII, nCII, ne, nHn, cool_rate_CII, xLmbd_CII
    real(kind=8) :: y_OI, nOI, cool_rate_OI, xLmbd_OI, y_OII, nOII,&
       cool_rate_OII, xLmbd_OII, y_OIII, nOIII, cool_rate_OIII, xLmbd_OIII
    real(kind=8),dimension(0:2) :: b_cont, dvdr
    real(kind=8) :: rho, t_ff, radius, v_bulk, dv_d, kappa
    real(kind=8) :: nCO, xNc_CO
    integer :: ii
    real(kind=8), parameter :: min_value = 1.d-50
    real(kind=8) :: frac_C, frac_O
    real(kind=8) :: dvdr_av, l_sbl, l_sbl_H2, l_sbl_CI, l_sbl_CO, v_th
    real(kind=8), save :: esc_OI(1:3)
    real(kind=8) :: xNc_OI, l_sbl_OI
    integer :: ifirst = 0
    real(kind=8), save :: esc_CII
    real(kind=8) :: l_sbl_CII, xNc_CII
    real(kind=8) :: xLmbd_CI, tau_sbl, cool_rate_CI
    real(kind=8) :: xNc_HD, xLmbd_HD, cool_rate_HD, l_sbl_HD
    real(kind=8) :: xNc_OH, l_sbl_OH, xLmbd_OH
    real(kind=8) ::xNc_H2O, xLmbd_H2O, l_sbl_H2O
    real(kind=8) :: pi
    data pi/3.14159265358979d0/
    logical :: isNotFinite
    if (ifirst == 0) then
      esc_OI(:) = 1.d0
      esc_CII = 1.d0
      ifirst = 1
    endif
    frac_C = xch%metal*MP_frac_C_solar
    frac_O = xch%metal*MP_frac_O_solar
    y_H = xch%ychem(0) 
    y_H2 = xch%ychem(1) 
    y_e = xch%ychem(2) 
    y_Hp = xch%ychem(3) 
    dvdr_av = (xch%dvdr(0)+xch%dvdr(1)+xch%dvdr(2))/3.d0+1.d-30
    v_th = sqrt(2.d0*CONST_kB*xch%Tg/CONST_amu)
    l_sbl = v_th / dvdr_av 
    l_sbl_H2 = l_sbl/sqrt(2.d0)
    xNc_H2 = y_H2*xch%nH*min(0.5*xch%xlmbdj,l_sbl_H2)
    ne = y_e*xch%nH
    nHn = y_H*xch%nH
    rho = (MP_mu*cgs_amu)*xch%nH 
    t_ff = sqrt(3.0*pi/(32.0*CONST_G*rho)) 
    radius = 0.5d0*xch%xlmbdj
    v_bulk = radius/3.0/t_ff
    y_CII = ( frac_C - xch%yco)
    y_CII = MAX(MIN(frac_C,y_CII),min_value)
    nCII = y_CII*xch%nH
    nCO = xch%yco * xch%nH
    l_sbl_CO = l_sbl/sqrt(MP_AC+MP_AO)
    xNc_CO = xch%yco*xch%nH*min(radius, l_sbl_CO)
    y_OI = frac_O*(1.d0 - y_Hp - xch%yco) 
    y_OI = MAX(MIN(frac_O,y_OI),min_value)
    nOI = y_OI*xch%nH
    call OII_OII_ratio(ne, y_Hp, xch%Tg, xch%rOII, y_OII, y_OIII, xch%metal)
    nOII = y_OII*xch%nH
    nOIII = y_OIII*xch%nH
    call find_dop(xch%Td, kappa) 
    tau_cnt = (kappa*xch%fd)*rho*radius 
    do ii = 0, 2
        b_cont(ii) = exp(-tau_cnt)
    enddo
    call H2cool(xch%nH,xch%Tg,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)
    T_5 = 1.d-5*xch%Tg
    xLmbd_Lya = y_e*y_H*7.50d-19/(1.d0 + dsqrt(T_5)) &
         *dexp(-1.18348d5/xch%Tg) *xch%nH*xch%nH
    xLmbd_cpt = 5.65d-36*((1.d0 + zred)**4)*(xch%Tg - MP_Tcmb)*xch%nH*y_e
    xLmbd_Hep = 5.54d-17/(1.d0 + sqrt(T_5))/(xch%Tg**0.397) &
         *exp(-4.73638d+5/xch%Tg)*yHe*(y_e*xch%nH)**2
    if (xch%Tg<3.2d5) then
       gff1 = 0.79464 + 0.1243 * dlog10(xch%Tg)
    else
       gff1 = 2.13164 - 0.1240 * dlog10(xch%Tg)
    endif
    xLmbd_ff = 1.426d-27 * sqrt(xch%Tg) * gff1* y_Hp * y_e * xch%nH * xch%nH
    xLmbd_CI = 0.d0
    l_sbl_CII = min(l_sbl/sqrt(MP_AC), radius)
    tau_sbl = (kappa*xch%fd) *rho *l_sbl_CII
    xNc_CII = nCII*l_sbl_CII
    call CIIcool(xch%nH, xch%Tg, MP_Tcmb ,xNc_CII,y_H2,y_H,y_e,esc_CII,tau_sbl,cool_rate_CII)
    xLmbd_CII = nCII * cool_rate_CII
    call COcool(xch%nH , xch%Tg ,y_H2,xNc_CO,tau_cnt,cooling_CO)
    call COcool(xch%nH ,MP_Tcmb ,y_H2,xNc_CO,tau_cnt,heating_CO)
    xLmbd_CO = nCO*xch%nH*(cooling_co - heating_co) 
    l_sbl_OI = min(l_sbl/sqrt(MP_AO), radius)
    xNc_OI = nOI*l_sbl_OI
    tau_sbl = (kappa*xch%fd)*rho*l_sbl_OI
    call OIcool(xch%nH, xch%Tg, MP_Tcmb, xNc_OI, y_H2, y_H, y_e, esc_OI, tau_sbl, cool_rate_OI)
    xLmbd_OI = nOI * cool_rate_OI 
    dv_d = sqrt(2.0*CONST_kB*xch%Tg/(MP_AO)/CONST_amu) 
    do ii = 0, 2
        dvdr(ii) = max(dv_d/radius, xch%dvdr(ii))
    end do
    call oxygenII_cooling(nOII, ne, xch%Tg, MP_Tcmb, dvdr, b_cont, cool_rate_OII)
    xLmbd_OII = nOII* cool_rate_OII 
    call oxygenIII_cooling(nOIII, ne, nHn, xch%Tg, MP_Tcmb, dvdr, b_cont, cool_rate_OIII)
    xLmbd_OIII = nOIII * cool_rate_OIII
    xLmbd_HD = 0.d0
    xLmbd_OH = 0.d0
    xLmbd_H2O = 0.d0
    xLmbd_line = xLmbd_H2 + xLmbd_Lya + xLmbd_cpt + xLmbd_Hep + xLmbd_ff &
               + xLmbd_CI + xLmbd_CII + xLmbd_OI + xLmbd_CO &
               + xLmbd_OII + xLmbd_OIII+ xLmbd_HD + xLmbd_OH + xLmbd_H2O
    if (dbg_flg_prm == 1) then
print '(A,(1P5E15.7))', "xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, xLmbd_ff: ",&
     xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, xLmbd_ff 
    end if
  end subroutine line_cool
  subroutine H2cool(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)
    real(kind=8),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt
    real(kind=8),intent(OUT) :: xLmbd_H2
    real(kind=8) :: xLdH2_n0,xLdH2_LTE,y_He,T_rad, xLdH2_dbg,xLdH2_dbg2
    logical :: isNotFinite
    y_He = yHe 
    T_rad = 0d0 
    if(xnH.lt.1.d+8) then
       call H2cool_LowD_G15(xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He,xLdH2_n0)
       call H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau_cnt,T_rad,xLdH2_LTE)
       xLmbd_H2 = xLdH2_LTE/(1.d0+xLdH2_LTE/xLdH2_n0)
    else
       call H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau_cnt,T_rad,xLdH2_LTE)
       xLmbd_H2 = xLdH2_LTE
    end if
  end subroutine H2cool
  subroutine H2cool_HM(xnH,T_K,y_H,y_H2,xLmbd_H2)
    real(kind=8),intent(IN) :: xnH,T_K,y_H,y_H2
    real(kind=8),intent(OUT) :: xLmbd_H2
    real(kind=8), parameter :: E_r_20 = 512.0d0*boltz
    real(kind=8), parameter :: E_r_31 = 853.3d0*boltz
    real(kind=8), parameter :: E_v_10 = 5860.0d0*boltz
    real(kind=8), parameter :: E_v_20 = 11720.0d0*boltz
    real(kind=8) :: lgT, sqT, T3,&
         kd_a, kd_m, L_r_H_LTE, L_r_H2_LTE, L_v_H_LTE, L_v_H2_LTE,&
         lg_L_r_H_n0, L_r_H_n0, L_r_H2_n0, L_v_H_n0, L_v_H2_n0,&
         kd_cr_r_H, kd_cr_v_H, kd_cr_r_H2, kd_cr_v_H2, L_vr_H,&
         L_vr_H2, Lda_rv,&
         Tiv, xkTiv, T3iv, kdaiv, kdmiv,&
         L_r_LTE, L_v_LTE, gamma_r_H2_2, gamma_r_H2_3, gamma_v_10_H,&
         gamma_v_20_H, gamma_v_10_H2
    if(y_H2.gt.0.d0 .and. T_K.gt.3.d+1) then
      kd_a = xnH*y_H
      kd_m = xnH*y_H2
      kdaiv = 1.d0/kd_a
      kdmiv = 1.d0/kd_m
      Tiv = 1.d0/T_K
      xkTiv = Tiv/boltz
      lgT = log10(T_K)
      sqT = sqrt(T_K)
      T3 = 1.d-3*T_K
      T3iv = 1.d+3*Tiv
      L_r_LTE = ( 9.5d-22*(T3**3.76)/(1.d0 + 0.12d0*(T3**2.1)))*exp(-(0.13*T3iv)**3)&
           + 3.0d-24*exp(-0.51d0*T3iv)
      L_v_LTE = ( 6.7d-19*exp(-5.86d0*T3iv) + 1.6d-18*exp(-1.17d1*T3iv))
      L_r_H_LTE = L_r_LTE*kdaiv
      L_r_H2_LTE = L_r_LTE*kdmiv
      L_v_H_LTE = L_v_LTE*kdaiv
      L_v_H2_LTE = L_v_LTE*kdmiv
      lg_L_r_H_n0 = - 1.03d2 + 9.759d+1*lgT - 4.805d+1*(lgT*lgT)& 
           + 1.08d+1*(lgT*lgT*lgT) - 0.9032d0*(lgT**4)
      L_r_H_n0 = 10.0**lg_L_r_H_n0
      gamma_r_H2_2 = ( 3.3d-12 + 6.6d-12*T3 )*7.0071412d-1
      gamma_r_H2_3 = ( 3.3d-12 + 6.6d-12*T3 )*1.0041880d+0
      L_r_H2_n0 = 0.25d0*( 5.0d0*gamma_r_H2_2 *exp(-E_r_20*xkTiv)*E_r_20 )&
           + 0.75d0*( 2.333d0*gamma_r_H2_3 *exp(-E_r_31*xkTiv)*E_r_31 )
      gamma_v_10_H = 1.0d-12*sqT*exp(-1.0d+3*Tiv)
      gamma_v_20_H = 1.6d-12*sqT*exp(-(4.0d+2*Tiv)**2)
      gamma_v_10_H2 = 1.4d-12*sqT*exp(-1.2d+4/(T_K + 1.2d+3))
      L_v_H_n0 = gamma_v_10_H*exp(-E_v_10*xkTiv)*E_v_10 + gamma_v_20_H*exp(-E_v_20*xkTiv)*E_v_20
      L_v_H2_n0 = gamma_v_10_H2*exp(-E_v_10*xkTiv)*E_v_10
      kd_cr_r_H = kd_a*L_r_H_LTE/L_r_H_n0
      kd_cr_v_H = kd_a*L_v_H_LTE/L_v_H_n0
      kd_cr_r_H2 = kd_m*L_r_H2_LTE/L_r_H2_n0
      kd_cr_v_H2 = kd_m*L_v_H2_LTE/L_v_H2_n0
      L_vr_H = L_r_H_LTE/( 1.d0 + kd_cr_r_H*kdaiv ) + L_v_H_LTE/( 1.d0 + kd_cr_v_H*kdaiv )
      L_vr_H2 = L_r_H2_LTE/( 1.d0 + kd_cr_r_H2*kdmiv ) + L_v_H2_LTE/( 1.d0 + kd_cr_v_H2*kdmiv )
      Lda_rv = kd_m*( y_H*L_vr_H + y_H2*L_vr_H2 )
      xLmbd_H2 = 0.55d0*Lda_rv*xnH 
      else
        xLmbd_H2 = 0.d0
      end if
      return
  end subroutine H2cool_HM
  function Q_bg(T_nu)
    real(kind=8) :: Q_bg
    real(kind=8),intent(IN) :: T_nu
    real(kind=8) :: T_rad, x, Q_bg_CMB
    T_rad = MP_Tcmb
    x = T_nu/T_rad
    if(x.gt.1.d2) then
       Q_bg_CMB = 0.d0
    else
       Q_bg_CMB = 1.d0/(dexp(x)-1.d0)
    endif
    Q_bg = Q_bg_CMB
    return
  end function Q_bg
  function beta_esc(tau_L,tau_C)
    real(kind=8) :: beta_esc
    real(kind=8) :: tau_L, tau_C
    if(tau_L.lt.0.d0) then
       beta_esc=1.d0
    elseif(tau_L.lt.1.d-5) then
       beta_esc = dexp(-tau_C)
    else
       beta_esc = dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
    endif
    return
  end function beta_esc
  function c_H2(T_K)
    real(kind=8),intent(IN) :: T_K
    real(kind=8) :: c_H2, dT
    integer :: ii
    real(kind=8) :: logTk
    real(kind=8),dimension(50) :: ca, xlTa
    data xlTa/ 1.00000d0,1.08163d0,1.16327d0,1.24490d0,1.32653d0,1.40816d0 &
    ,1.48980d0,1.57143d0,1.65306d0,1.73469d0,1.81633d0,1.89796d0,1.97959d0 &
    ,2.06122d0,2.14286d0,2.22449d0,2.30612d0,2.38776d0,2.46939d0,2.55102d0 &
    ,2.63265d0,2.71429d0,2.79592d0,2.87755d0,2.95918d0,3.04082d0,3.12245d0 &
    ,3.20408d0,3.28571d0,3.36735d0,3.44898d0,3.53061d0,3.61224d0,3.69388d0 &
    ,3.77551d0,3.85714d0,3.93878d0,4.02041d0,4.10204d0,4.18367d0,4.26531d0 &
    ,4.34694d0,4.42857d0,4.51020d0,4.59184d0,4.67347d0,4.75510d0,4.83673d0 &
    ,4.91837d0,5.00000d0 /
    data ca/ 1.50000d0,1.50000d0,1.50000d0,1.50000d0,1.50000d0,1.50000d0 &
    ,1.50000d0,1.50000d0,1.50482d0,1.51494d0,1.53443d0,1.59157d0,1.68857d0 &
    ,1.81833d0,1.97800d0,2.14694d0,2.27310d0,2.37931d0,2.44969d0,2.48382d0 &
    ,2.49866d0,2.50000d0,2.50559d0,2.52365d0,2.55735d0,2.62265d0,2.71562d0 &
    ,2.83053d0,2.95857d0,3.07868d0,3.17516d0,3.26061d0,3.32957d0,3.38085d0 &
    ,3.41796d0,3.44314d0,3.46176d0,3.47367d0,3.48120d0,3.48937d0,3.49000d0 &
    ,3.49518d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0 &
    ,3.50000d0,3.50000d0 /
    logTk = dlog10(T_K)
    ii = int((logTk-1.d0)*49.d0/4.d0) + 1
    if(ii < 1) then
      c_H2 = ca(1)
    else if(ii>49) then
      c_H2 = ca(50)
    else
      c_H2=(ca(ii+1)-ca(ii))/(xlTa(ii+1)-xlTa(ii))*(logTk-xlTa(ii))+ca(ii)
    endif
  end function c_H2
  function c_H2_2(T_K)
    real(kind=8),intent(IN) :: T_K
    real(kind=8) :: c_H2_2
    real(kind=8),dimension(30) :: xlTa, ca
    real(kind=8) :: xlT
      data xlta / 0.150d1, 0.160d1, 0.171d1, 0.181d1, 0.191d1,&
                  0.202d1, 0.212d1, 0.222d1, 0.233d1, 0.243d1,&
                  0.253d1, 0.264d1, 0.274d1, 0.284d1, 0.295d1,&
                  0.305d1, 0.316d1, 0.326d1, 0.336d1, 0.347d1,&
                  0.357d1, 0.367d1, 0.378d1, 0.388d1, 0.398d1,&
                  0.409d1, 0.419d1, 0.429d1, 0.440d1, 0.450d1 /
      data ca / 0.150d1, 0.150d1, 0.151d1, 0.153d1, 0.160d1,&
                0.174d1, 0.193d1, 0.214d1, 0.231d1, 0.243d1,&
                0.248d1, 0.250d1, 0.250d1, 0.251d1, 0.255d1,&
                0.263d1, 0.276d1, 0.292d1, 0.307d1, 0.320d1,&
                0.330d1, 0.337d1, 0.342d1, 0.345d1, 0.347d1,&
                0.348d1, 0.349d1, 0.349d1, 0.350d1, 0.350d1 /
      xlT = dlog10(T_K)
      call linear(xlTa,ca,30,xlT,c_H2_2)
      return
    contains
    subroutine linear(xa,ya,m,x,y)
      integer,intent(IN) :: m
      real(kind=8),dimension(m),intent(IN) :: xa,ya
      real(kind=8),intent(IN) :: x
      real(kind=8),intent(OUT) :: y
      integer :: ms, i
      real(kind=8) :: t, y1, y2
      ms = 0
      do i=1,m
         if(x-xa(i).le.0.d0 .or. i.eq.m) then
            ms = i
            exit
         endif
      end do
      if(ms.eq.0) then
        print *, 'zero ms in linear (primordial)'
        stop
      end if
      if(ms.eq.1) ms=2
      y1 = ya(ms-1)
      y2 = ya(ms)
      t = (x - xa(ms-1))/(xa(ms) - xa(ms-1))
      y = (1.d0-t)*y1 + t*y2
      return
    end subroutine linear
  end function c_H2_2
  subroutine H2cool_Omukai(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)
    real(kind=8),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt
    real(kind=8),intent(OUT) :: xLmbd_H2
    real(kind=8) :: A(0:2,0:2,0:22,0:22),ET(0:2,0:22),p(0:22),f(0:2,0:22),f_v(0:2)
    real(kind=8) :: xk_B,h_P,pi,xm_p,c_light
    data xk_B/1.380662d-16/,h_P/6.626176d-27/,pi/3.14159265358979d0/,xm_p/1.67d-24/,c_light/2.99792458d10/
    real(kind=8) :: gamma_h_j2down_omukai(0:18)
    real(kind=8) g_0, g_1, g_2, A_10, A_20, A_21,&
         gamma_10H, gamma_20H, gamma_21H, gamma_10H2, gamma_20H2, gamma_21H2,&
         gamma_10e, gamma_20e, gamma_21e, gamma_10Hp, gamma_20Hp, gamma_21Hp,&
         xn_a, xn_m,&
         C_10, C_20, C_21, DT_10, DT_20, DT_21, C_01, C_02, C_12,&
         Q_10, Q_20, Q_21,&
         R_10, R_20, R_21, R_01, R_02, R_12,&
         f_0, f_1, f_2, f_para, f_ortho, z_para, z_ortho, &
         gamma_H, gamma_H2, gamma_e, gamma_Hp, xj, xj2, xj3, C_ul, DT, Q_ul, g_u, g_l, r,&
         v_th, DE, xnu_Hz, xnu_Hz3, xNc_u, xNc_l, tau_ul, esc, x, xmeV, xn_e, xn_Hp
    integer :: i,ii,iv,ivf,ivi,j,jf,ji
    data (((A(i,ii,j,j+2),ii=0,i-1),i=1,2),j=0,20) &
         /8.54d-7 ,3.47d-7 ,1.29d-6 &
         ,4.23d-7 ,1.61d-7 ,6.40d-7 &
         ,2.90d-7 ,1.03d-7 ,4.41d-7 &
         ,2.09d-7 ,6.98d-8 ,3.18d-7 &
         ,1.50d-7 ,4.72d-8 ,2.28d-7 &
         ,1.06d-7 ,3.15d-8 ,1.62d-7 &
         ,7.38d-8 ,2.05d-8 ,1.12d-7 &
         ,4.98d-8 ,1.31d-8 ,7.50d-8 &
         ,3.27d-8 ,8.07d-9 ,4.88d-8 &
         ,2.07d-8 ,4.83d-9 ,3.06d-8 &
         ,1.27d-8 ,2.79d-9 ,1.85d-8 &
         ,7.44d-9 ,1.55d-9 ,1.07d-8 &
         ,4.16d-9 ,8.27d-10,5.84d-9 &
         ,2.19d-9 ,4.21d-10,3.00d-9 &
         ,1.08d-9 ,2.04d-10,1.43d-9 &
         ,4.87d-10,9.45d-11,6.17d-10 &
         ,1.96d-10,4.23d-11,2.33d-10 &
         ,6.71d-11,1.91d-11,7.32d-11 &
         ,1.82d-11,9.63d-12,1.73d-11 &
         ,3.38d-12,6.26d-12,2.46d-12 &
         ,2.96d-13,5.78d-12,1.16d-13/
    data (((A(i,ii,j,j),ii=0,i-1),i=1,2),j=1,20) &
         /4.29d-7 ,1.94d-7 ,6.37d-7 &
         ,3.03d-7 ,1.38d-7 ,4.50d-7 &
         ,2.78d-7 ,1.29d-7 ,4.12d-7 &
         ,2.65d-7 ,1.25d-7 ,3.91d-7 &
         ,2.55d-7 ,1.23d-7 ,3.74d-7 &
         ,2.45d-7 ,1.21d-7 ,3.58d-7 &
         ,2.34d-7 ,1.20d-7 ,3.40d-7 &
         ,2.23d-7 ,1.18d-7 ,3.22d-7 &
         ,2.12d-7 ,1.17d-7 ,3.03d-7 &
         ,1.99d-7 ,1.15d-7 ,2.84d-7 &
         ,1.87d-7 ,1.13d-7 ,2.63d-7 &
         ,1.74d-7 ,1.11d-7 ,2.42d-7 &
         ,1.61d-7 ,1.08d-7 ,2.21d-7 &
         ,1.47d-7 ,1.05d-7 ,2.01d-7 &
         ,1.34d-7 ,1.02d-7 ,1.80d-7 &
         ,1.21d-7 ,9.88d-8 ,1.60d-7 &
         ,1.08d-7 ,9.50d-8 ,1.41d-7 &
         ,9.61d-8 ,9.07d-8 ,1.22d-7 &
         ,8.43d-8 ,8.62d-8 ,1.05d-7 &
         ,7.32d-8 ,8.14d-8 ,8.85d-8/
    data (((A(i,ii,j,j-2),ii=0,i),i=0,2),j=2,20) &
         /2.94d-11,2.53d-7 ,2.79d-11 &
         ,1.27d-7 ,3.68d-7 ,2.56d-11 &
         ,4.76d-10,3.47d-7,4.50d-10 &
         ,1.90d-7 ,4.98d-7 ,4.12d-10 &
         ,2.76d-9 ,3.98d-7 ,2.59d-9 &
         ,2.38d-7 ,5.60d-7 ,2.37d-9 &
         ,9.84d-9 ,4.21d-7 ,9.21d-9 &
         ,2.77d-7 ,5.77d-7 ,8.37d-9 &
         ,2.64d-8 ,4.19d-7 ,2.46d-8 &
         ,3.07d-7 ,5.57d-7 ,2.22d-8 &
         ,5.88d-8 ,3.96d-7 ,5.44d-8 &
         ,3.28d-7 ,5.05d-7 ,4.88d-8 &
         ,1.14d-7 ,3.54d-7 ,1.05d-7 &
         ,3.39d-7 ,4.30d-7 ,9.33d-8 &
         ,2.00d-7 ,2.98d-7 ,1.82d-7 &
         ,3.40d-7 ,3.38d-7 ,1.61d-7 &
         ,3.24d-7 ,2.34d-7 ,2.92d-7 &
         ,3.30d-7 ,2.41d-7 ,2.55d-7 &
         ,4.90d-7 ,1.68d-7 ,4.38d-7 &
         ,3.12d-7 ,1.49d-7 ,3.79d-7 &
         ,7.03d-7 ,1.05d-7 ,6.21d-7 &
         ,2.85d-7 ,7.20d-8 ,5.32d-7 &
         ,9.64d-7 ,5.30d-8 ,8.42d-7 &
         ,2.53d-7 ,1.96d-8 ,7.13d-7 &
         ,1.27d-6 ,1.65d-8 ,1.10d-6 &
         ,2.16d-7 ,1.49d-11,9.18d-7 &
         ,1.62d-6 ,4.26d-10,1.38d-6 &
         ,1.78d-7 ,1.91d-8 ,1.14d-6 &
         ,2.00d-6 ,8.38d-9 ,1.69d-6 &
         ,1.39d-7 ,8.07d-8 ,1.37d-6 &
         ,2.41d-6 ,4.27d-8 ,2.00d-6 &
         ,1.01d-7 ,1.86d-7 ,1.61d-6 &
         ,2.83d-6 ,1.04d-7 ,2.32d-6 &
         ,6.80d-8 ,3.35d-7 ,1.84d-6 &
         ,3.26d-6 ,1.93d-7 ,2.64d-6 &
         ,3.98d-8 ,5.24d-7 ,2.05d-6 &
         ,3.68d-6 ,3.08d-7 ,2.93d-6 &
         ,1.84d-8 ,7.49d-7 ,2.23d-6/
    g_0=1.d0
    g_1=1.d0
    g_2=1.d0
    A_10=8.3d-7
    A_20=4.1d-7
    A_21=1.1d-6
    gamma_10H=1.0d-12*dsqrt(T_K)*dexp(-1.d3/T_K)
    gamma_20H=1.6d-12*dsqrt(T_K)*dexp(-(4.d2/T_K)**2)
    gamma_21H=4.5d-12*dsqrt(T_K)*dexp(-(5.d2/T_K)**2)
    gamma_10H2=1.4d-12*dsqrt(T_K)*dexp(-1.81d4/(T_K+1.2d3))
    gamma_20H2=0.d0
    gamma_21H2=gamma_10H2
    do iv=0,2
       do j=0,22
          ET(iv,j)=E_H2_BFM(iv,J)
       end do
    end do
    DT_10=ET(1,0)-ET(0,0)
    DT_20=ET(2,0)-ET(0,0)
    DT_21=ET(2,0)-ET(1,0)
    gamma_10e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_10/T_K)
    gamma_20e=2.5d-12*(T_K**0.5d0)/(1.d0+0.5d0*DT_20/T_K)
    gamma_21e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_21/T_K)
    gamma_10Hp=1.4d-4*T_K**(-1.344d0) &
         *(1.d0+4.005d-9*T_K**2.066d0) &
         *dexp((DT_10-9589d0)/T_K)
    gamma_20Hp=4.585d-5*T_K**(-1.291d0) &
         *(1.d0+1.378d-8*T_K**1.903d0) &
         *dexp((DT_20-1.933d4)/T_K)
    gamma_21Hp=5.593d-3*T_K**(-1.770d0) &
         *(1.d0+3.505d-10*T_K**2.406d0) &
         *dexp((DT_21-1.460d4)/T_K)
    xn_a=y_H*xnH
    xn_m=y_H2*xnH
    xn_e=y_e*xnH
    xn_Hp=y_Hp*xnH
    C_10=gamma_10H*xn_a+gamma_10H2*xn_m &
         +gamma_10e*xn_e+gamma_10Hp*xn_Hp
    C_20=gamma_20H*xn_a+gamma_20H2*xn_m &
         +gamma_20e*xn_e+gamma_20Hp*xn_Hp
    C_21=gamma_21H*xn_a+gamma_21H2*xn_m &
         +gamma_21e*xn_e+gamma_21Hp*xn_Hp
    C_01=C_10*dexp(-DT_10/T_K)
    C_02=C_20*dexp(-DT_20/T_K)
    C_12=C_21*dexp(-DT_21/T_K)
    Q_10=Q_bg(DT_10)
    Q_20=Q_bg(DT_20)
    Q_21=Q_bg(DT_21)
    R_10=A_10*(1.d0+Q_10)+C_10
    R_20=A_20*(1.d0+Q_20)+C_20
    R_21=A_21*(1.d0+Q_21)+C_21
    R_01=(g_1/g_0)*A_10*Q_10+C_01
    R_02=(g_2/g_0)*A_20*Q_20+C_02
    R_12=(g_2/g_1)*A_21*Q_21+C_12
    f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) ) &
         /( (R_01+R_02+R_20)*(R_10+R_12+R_21) &
         -(R_01-R_21)*(R_10-R_20) )
    f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
    f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)
    f_v(0)=f_0
    f_v(1)=f_1
    f_v(2)=f_2
    f_para=0.25d0
    f_ortho=0.75d0
    do iv=0,2
       z_para=0.d0
       z_ortho=0.d0
       p(0)=1.d0
       p(1)=1.d0
       do j=0,18
          if(mod(j,2).eq.0) then
             z_para=z_para+p(j)
          else
             z_ortho=z_ortho+p(j)
          endif
          if(j == 0) then
             gamma_H=2.93d-14+1.21d-15*T_K &
                  +2.16d-19*T_K**2+1.32d-21*T_K**3
          elseif(j == 1) then
             gamma_H=8.34d-14+5.97d-16*T_K &
                  +7.76d-19*T_K**2+1.72d-21*T_K**3
          elseif(j == 2) then
             gamma_H=7.54d-14+2.65d-16*T_K &
                  +2.14d-19*T_K**2+2.65d-21*T_K**3
          elseif(j == 3) then
             gamma_H=2.95d-14+1.21d-16*T_K &
                  +5.53d-19*T_K**2+2.51d-21*T_K**3
          else
             gamma_H=4.6d-12*(2.d0*dble(j+2)-3.d0)*dsqrt(T_K) &
                  *dsqrt(1.d0+2.d0*85.25d0*(2.d0*dble(j+2)-1.d0)/T_K) &
                  *dexp(-(1.d1*85.25d0*(2.d0*dble(j+2)-1.d0)) &
                  /(T_K+85.25d0*dble(j+2)*dble(j+3)) &
                  -0.1187d0*(4.d0*dble(j+2)-2.d0))
          endif
          gamma_H_j2down_omukai(j)=gamma_H
          gamma_H2=(3.3d-12+6.6d-15*T_K) &
               *0.276d0*dble((j+2)**2) &
               *dexp(-(dble(j+2)/3.18d0)**1.7d0)
          xmeV=11.604d0
          if(j == 0) then
             gamma_Hp=2.889d-10*dexp(0.171d0*xmeV/T_K)
          elseif(j == 1) then
             gamma_Hp=8.202d-10*dexp(0.250d0*xmeV/T_K)
          elseif(j == 2) then
             gamma_Hp=5.700d-10*dexp(0.159d0*xmeV/T_K)
          elseif(j == 3) then
             gamma_Hp=8.883d-10*dexp(0.223d0*xmeV/T_K)
          elseif(j == 4) then
             gamma_Hp=5.209d-10*dexp(0.124d0*xmeV/T_K)
          elseif(j == 5) then
             gamma_Hp=7.977d-10*dexp(0.177d0*xmeV/T_K)
          elseif(j == 6) then
             gamma_Hp=4.489d-10*dexp(0.101d0*xmeV/T_K)
          else
             gamma_Hp=5.915d-10*dexp(0.094d0*xmeV/T_K)
          endif
          xj=dble(j)
          x=(ET(iv,j+2)-ET(iv,j))/T_K
          gamma_e=1.d-10*(xj+2.d0)*(xj+1.d0)/(2.d0*xj+5.d0)*(1.d0+1.5d0/x)
          C_ul=gamma_H*xn_a+gamma_H2*xn_m &
               +gamma_Hp*xn_Hp+gamma_e*xn_e
          DT=ET(iv,j+2)-ET(iv,j)
          Q_ul=Q_bg(DT)
          g_u=2.d0*dble(j)+5.d0
          g_l=2.d0*dble(j)+1.d0
          r=(g_u/g_l)*(A(iv,iv,j+2,j)*Q_ul+C_ul*dexp(-DT/T_K)) &
               /(A(iv,iv,j+2,j)*(1.d0+Q_ul)+C_ul)
          p(j+2)=p(j)*r
       enddo
       do j=0,20
          if(mod(j,2).eq.0) then
             f(iv,j)=(p(j)/z_para)*f_para*f_v(iv)
          else
             f(iv,j)=(p(j)/z_ortho)*f_ortho*f_v(iv)
          endif
       enddo
    enddo
    v_th=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))
    xLmbd_H2=0.d0
    do ji=0,18
       do ivi=1,2
          do ivf=0,ivi-1
             jf=ji+2
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l) &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do
    do ji=1,20
       do ivi=1,2
          do ivf=0,ivi-1
             jf=ji
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l) &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do
    do ji=2,20
       do ivi=0,2
          do ivf=0,ivi
             jf=ji-2
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l) &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do
    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 
    return
  end subroutine H2cool_Omukai
  function E_H2_BFM(iv,J)
    integer,intent(IN) :: iv,J
    real(kind=8) :: E_H2_BFM
    real(kind=8) :: vv,rrot,Ev,Ev0,Ev1,Ev2,Ev3,Ev4,Bv,Bv0,Bv1,Bv2,Bv3,Bv4,&
         Dv,Dv0,Dv1,Dv2,Dv3,Dv4,Fv,Fv0,Fv1,Fv2,Fv3,Fv4,Gv,Gv0,Gv1,Gv2,Gv3,Gv4,&
         Hv,Hv0,Hv1,Hv2,Hv3,Hv4,Ov,Ov0,Ov1,Ov2,Ov3,Ov4,Pv,Pv0,Pv1,Pv2,Pv3,Pv4
    Ev0=0.38496d0
    Ev1=-0.04609d0
    Ev2=0.00178d0
    Ev3=-7.d-5
    Ev4=2.9511d-6
    Bv0=54.438d0
    Bv1=4.6063d0
    Bv2=-2.0050d0
    Bv3=0.19260d0
    Bv4=-6.2953d-3
    Dv0=-0.72593d0
    Dv1=5.9990d0
    Dv2=-1.5187d0
    Dv3=0.12721d0
    Dv4=-3.3391d-3
    Fv0=-12.662d0
    Fv1=20.047d0
    Fv2=-4.7873d0
    Fv3=0.36900d0
    Fv4=-9.003d-3
    Gv0=-24.006d0
    Gv1=33.989d0
    Gv2=-7.6841d0
    Gv3=0.52413d0
    Gv4=-1.1297d-2
    Hv0=-22.384d0
    Hv1=31.150d0
    Hv2=-6.6139d0
    Hv3=0.39226d0
    Hv4=-9.496d-3
    Ov0=-10.541d0
    Ov1=14.746d0
    Ov2=-2.9476d0
    Ov3=0.16016d0
    Ov4=-6.5005d-3
    Pv0=-2.0021d0
    Pv1=2.8400d0
    Pv2=-0.54654d0
    Pv3=0.031636d0
    Pv4=-2.4398d-3
    vv=dble(iv)+0.5d0
    Ev=Ev0+Ev1*vv+Ev2*vv**2+Ev3*vv**3+Ev4*vv**4
    Bv=Bv0+Bv1*vv+Bv2*vv**2+Bv3*vv**3+Bv4*vv**4
    Dv=Dv0+Dv1*vv+Dv2*vv**2+Dv3*vv**3+Dv4*vv**4
    Fv=Fv0+Fv1*vv+Fv2*vv**2+Fv3*vv**3+Fv4*vv**4
    Gv=Gv0+Gv1*vv+Gv2*vv**2+Gv3*vv**3+Gv4*vv**4
    Hv=Hv0+Hv1*vv+Hv2*vv**2+Hv3*vv**3+Hv4*vv**4
    Ov=Ov0+Ov1*vv+Ov2*vv**2+Ov3*vv**3+Ov4*vv**4
    Pv=Pv0+Pv1*vv+Pv2*vv**2+Pv3*vv**3+Pv4*vv**4
    rrot=dble(J*(J+1))
    E_H2_BFM=1.43879d0* &
         (-Ev*1.d5+Bv*rrot-Dv*1.d-2*rrot**2+Fv*1.d-5*rrot**3 &
         -Gv*1.d-8*rrot**4+Hv*1.d-11*rrot**5-Ov*1.d-14*rrot**6 &
         +Pv*1.d-17*rrot**7)
    return
  end function E_H2_BFM
  subroutine H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau,T_rad,xLmbd_H2)
    real(kind=8),intent(IN) :: xnH,T_K,y_H2,xNc_H2,tau,T_rad
    real(kind=8),intent(OUT) :: xLmbd_H2
    real(kind=8) :: fesc,xlgT3,a,xLmbd,ar,xLmbdr
    logical :: isNotFinite
    fesc=0.d0
    if(T_K > 1.d2) then
       xlgT3=dlog10(T_K*1.d-3)
       a=-20.584225d0+5.0194035d0*xlgT3-1.5738805d0*xlgT3**2 &
            -4.7155769d0*xlgT3**3+2.4714161d0*xlgT3**4 &
            +5.4710750d0*xlgT3**5-3.9467356d0*xlgT3**6 &
            -2.2148338d0*xlgT3**7+1.8161874d0*xlgT3**8
       xLmbd=10.d0**a
    else
       xLmbd=0.d0
    endif
    if(T_rad > 1.d2) then
       xlgT3=dlog10(T_rad*1.d-3)
       ar=-20.584225d0+5.0194035d0*xlgT3-1.5738805d0*xlgT3**2 &
            -4.7155769d0*xlgT3**3+2.4714161d0*xlgT3**4 &
            +5.4710750d0*xlgT3**5-3.9467356d0*xlgT3**6 &
            -2.2148338d0*xlgT3**7+1.8161874d0*xlgT3**8
       xLmbdr=10.d0**ar
    else
       xLmbdr=0.d0
    endif
    xLmbd_H2=(xLmbd-xLmbdr)/xnH
    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 
    call f_fit(xNc_H2,T_K,fesc)
    xLmbd_H2=xLmbd_H2*fesc*dexp(-tau)
  end subroutine H2cool_LTEfit
  subroutine f_fit(xNc_H2,T_K,fesc)
    real(kind=8),intent(IN) :: xNc_H2,T_K
    real(kind=8),intent(OUT) :: fesc
    real(kind=8) :: a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5,alpha,beta,xlnf
    if(T_K < 1.d3) then
       a0= 0.978382d0
       a1= 0.000399572d0
       a2= -1.73108d-006
       a3= 1.15363d-009
       a4= 8.24607d-013
       a5= -7.65975d-016
    elseif(T_K <4.d3) then
       a0=0.827472d0
       a1= 0.000107697d0
       a2= -8.25123d-008
       a3= 1.92812d-012
       a4= 5.7192d-015
       a5= -7.869d-019
    else
       a0= 1.11569d0
       a1= -0.000329302d0
       a2= 1.01846d-007
       a3= -1.46666d-011
       a4= 1.00764d-015
       a5= -2.68873d-020
    endif
    alpha=a0+a1*T_K+a2*T_K**2+a3*T_K**3+a4*T_K**4+a5*T_K**5
    b0 = 24.0561d0
    b1 = 0.00110043d0
    b2 = -2.87224d-007
    b3 = 6.11525d-011
    b4 = -6.55034d-015
    b5 = 2.54997d-019
    beta=b0+b1*T_K+b2*T_K**2+b3*T_K**3+b4*T_K**4+b5*T_K**5
    if(beta < 200.d0) then
      xlnf=-alpha*dlog10(1.d0+xNc_H2/10.d0**beta)
      fesc=10.d0**xlnf
    else
      fesc = 1.d0
    endif
  end subroutine f_fit
  subroutine H2cool_LowD_G15(xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He,xLmbd_H2)
    real(kind=8),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He
    real(kind=8),intent(OUT) :: xLmbd_H2
    real(kind=8) :: x_o,x_p,a0,a1,a2,a3,a4,a5,a6,a7,a8, &
         T3,xlgT3,xlgT3_2,xlgT3_3,xlgT3_4,xlgT3_5,xlgT3_6,xlgT3_7,xlgT3_8, &
         xLd_pH2H,xlgLd_pH2H,xLd_oH2H,xlgLd_oH2H,xLd_H2H,xlgLd_pH2pH2,xLd_pH2pH2, &
         xlgLd_pH2oH2,xLd_pH2oH2,xlgLd_oH2pH2,xLd_oH2pH2,&
         xlgLd_oH2oH2,xLd_oH2oH2,xLd_H2H2,&
         xlgLd_pH2He,xLd_pH2He,xlgLd_oH2He,xLd_oH2He,xLd_H2He, &
         xlgLd_H2Hp,xLd_H2Hp,xlgLd_H2e,xLd_H2e
    logical :: isNotFinite
    x_o=0.75d0
    x_p=0.25d0
    if(T_K < 1000.d0) then
       a0=-24.330855d0
       a1=4.4404496d0
       a2=-4.0460989d0
       a3=-1.1390725d0
       a4=9.8094223d0
       a5=8.6273872d0
       a6=0.0d0
       a7=0.0d0
    else
       a0=-24.329086d0
       a1=4.6105087d0
       a2=-3.9505350d0
       a3=12.363818d0
       a4=-32.403165d0
       a5=48.853562d0
       a6=-38.542008d0
       a7=12.066770d0
    endif
    T3=T_K*1.d-3
    xlgT3=dlog10(T3)
    xlgT3_2=xlgT3**2
    xlgT3_3=xlgT3_2*xlgT3
    xlgT3_4=xlgT3_3*xlgT3
    xlgT3_5=xlgT3_4*xlgT3
    xlgT3_6=xlgT3_5*xlgT3
    xlgT3_7=xlgT3_6*xlgT3
    xlgT3_8=xlgT3_7*xlgT3
    if(T_K < 100.d0) then
       xLd_oH2H=5.09d-27*(T3**0.5d0)*dexp(-852.5d0/T_K)
    else
       xlgLd_oH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7
       xLd_oH2H=10.d0**xlgLd_oH2H
    endif
    if(T_K < 1000.d0) then
       a0=-24.216387d0
       a1=3.3237480d0
       a2=-11.642384d0
       a3=-35.553366d0
       a4=-35.105689d0
       a5=-10.922078d0
       a6=0.0d0
       a7=0.0d0
    else
       a0=-24.216387d0
       a1=4.2046488d0
       a2=-1.3155285d0
       a3=-1.6552763d0
       a4=4.1780102d0
       a5=-0.56949697d0
       a6=-3.3824407d0
       a7=1.0904027d0
    endif
    if(T_K < 100.d0) then
       xLd_pH2H=8.16d-26*(T3**0.5d0)*dexp(-509.85d0/T_K)
    else
       xlgLd_pH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7
       xLd_pH2H=10.d0**xlgLd_pH2H
    endif
    xLd_H2H=(x_o/(x_o+x_p))*xLd_oH2H +(x_p/(x_o+x_p))*xLd_pH2H
    a0=-23.889798d0
    a1=1.8550774d0
    a2=-0.55593388d0
    a3=0.28429361d0
    a4=-0.20581113d0
    a5=0.13112378d0
    xlgLd_pH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2pH2=10.d0**xlgLd_pH2pH2
    a0=-23.748534d0
    a1=1.76676480d0
    a2=-0.58634325d0
    a3=0.31074159d0
    a4=-0.17455629d0
    a5=0.18530758d0
    xlgLd_pH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2oH2=10.d0**xlgLd_pH2oH2
    a0=-24.126177d0
    a1=2.3258217d0
    a2=-1.0082491d0
    a3=0.54823768d0
    a4=-0.33679759d0
    a5=0.20771406d0
    xlgLd_oH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2pH2=10.d0**xlgLd_oH2pH2
    a0=-24.020047d0
    a1=2.2687566d0
    a2=-1.0200304d0
    a3=0.83561432d0
    a4=-0.40772247d0
    a5=0.096025713d0
    xlgLd_oH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2oH2=10.d0**xlgLd_oH2oH2
    xLd_H2H2=x_p**2*xLd_pH2pH2+x_p*x_o*xLd_pH2oH2 &
         +x_o*x_p*xLd_oH2pH2+x_o**2*xLd_oH2oH2
    a0=-23.489029d0
    a1=1.8210825d0
    a2=-0.59110559d0
    a3=0.42280623d0
    a4=-0.30171138d0
    a5=0.12872839d0
    xlgLd_pH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2He=10.d0**xlgLd_pH2He
    a0=-23.7749d0
    a1=2.40654d0
    a2=-1.23449d0
    a3=0.739874d0
    a4=-0.258940d0
    a5=0.120573d0
    xlgLd_oH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2He=10.d0**xlgLd_oH2He
    xLd_H2He=(x_o/(x_o+x_p))*xLd_oH2He &
         +(x_p/(x_o+x_p))*xLd_pH2He
    if(10.d0<T_K) then
       a0=-22.089523d0
       a1=1.5714711d0
       a2=0.015391166d0
       a3=-0.23619985d0
       a4=-0.51002221d0
       a5=0.32168730d0
       xlgLd_H2Hp=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5
       xLd_H2Hp=10.d0**xlgLd_H2Hp
    else
       xLd_H2Hp=0.d0
    endif
    if(100.d0<T_K .and. T_K<500.d0) then
       a0=-21.928796d0
       a1=16.815730d0
       a2=96.743155d0
       a3=343.19180d0
       a4=734.71651d0
       a5=983.67576d0
       a6=801.81247d0
       a7=364.14446d0
       a8=70.609154d0
       xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 &
            +a7*xlgT3_7+a8*xlgT3_8
       xLd_H2e=10.d0**xlgLd_H2e
    else if(500.d0<T_K) then
       a0=-22.921189d0
       a1=1.6802758d0
       a2=0.93310622d0
       a3=4.0406627d0
       a4=-4.7274036d0
       a5=-8.8077017d0
       a6=8.9167183d0
       a7=6.4380698d0
       a8=-6.3701156d0
       xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 &
            +a7*xlgT3_7+a8*xlgT3_8
       xLd_H2e=10.d0**xlgLd_H2e
    else
       xLd_H2e=0.d0
    endif
    xLmbd_H2=xLd_H2H*y_H+xLd_H2H2*y_H2+xLd_H2He*y_He &
         +xLd_H2Hp*y_Hp+xLd_H2e*y_e
    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 
  end subroutine H2cool_LowD_G15
  subroutine CoolSolverExplicit(xch,t_chemcool,t_cool)
    type(chem_vals) :: xch
    real(kind=8),intent(INOUT) :: t_chemcool, t_cool
    real(kind=8) :: xmu, gamma, c_v, t_chem, xLmbd_tot, en, pr, cs, T_new, rho, radius, chi_d
    real(kind=8),dimension(0:25 -1) :: xk
    xmu = get_xmu(xch%ychem)
    pr = xch%nH*MP_mu/xmu * boltz * xch%Tg 
    gamma = 1.d0+(1.d0+4.d0*yHe) & 
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))
    en = pr / (gamma-1.) 
    rho = (MP_mu*cgs_amu)*xch%nH 
    radius = 0.5d0*xch%xlmbdj
    call dtemp_radtr(xch%nH, xch%Tg, xch%Td, xch%EradIR, xch%chi_d &
    , dph_in=xch%rdph)
    call react_coef(xch, xk)
    call chemreact_adptv(xch, t_chem, xk_in=xk)
    call update_yco(xch%nH,xch%Tg,xch%ychem,xch%yco,xch%rcopd,xch%metal,xch%dt)
    call tot_cool(xch,xLmbd_tot,xk_in=xk)
    en = en - xLmbd_tot * xch%dt
    xmu = get_xmu(xch%ychem)
    gamma = 1.d0+(1.d0+4.d0*yHe) & 
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))
    pr = (gamma-1.) * en
    T_new = pr / (xch%nH*MP_mu/xmu * boltz) 
    xch%Tg = T_new
    t_cool = abs(en/xLmbd_tot) 
    t_chemcool = min(t_cool,t_chem) 
  end subroutine CoolSolverExplicit
 subroutine CoolSolverImplicit(xch)
    type(chem_vals) :: xch
    real(kind=8),parameter :: err_eps=1d-2 
    real(kind=8),parameter :: T_min=2.d0, T_max=1d5 
    integer,parameter :: k_nr_max = 30 
    integer,parameter :: k_srch_max = 100 
    integer,parameter :: k_bs_max = 100 
    integer :: k_nr, k_bs, k_srch
    real(kind=8) :: G, G0, G1,G_o, err_G, deno, T_o, T_n, T_n0, T_n1, dT_K, err_Tn, y_o(0:6 -1)
    logical :: force_substep 
    real(kind=8) :: gfuv, dph, drcopd
    real(kind=8) :: yco_o
    real(kind=8),dimension(0:25 -1) :: xk
    force_substep = .False.
    T_o = xch%Tg
    y_o(:) = xch%ychem(:)
    yco_o = xch%yco
    gfuv = xch%rgfuv
    dph = xch%rdph
    drcopd = xch%rcopd
    call UpdateTY_with_Tn(xch, T_o)
    G_o = T_o - xch%Tg
    if (dbg_flg_prm == 1) then
       print '(A,(1P3E15.7),/)', "(CoolSolverImplicit, beginning) T_o, T_K, G_o",&
            T_o, xch%Tg, G_o
    end if
    if (abs(G_o/T_o) > 0.5) then
       if (dbg_flg_prm == 1) print '(A,(1P1E15.7),A)', "(CoolSolverImplicit) abs(G_o/T_o) = ", &
            abs(G_o/T_o), " > 0.5    =>    force_substep"
       force_substep = .True.
       xch%Tg = T_o
       xch%ychem(:) = y_o(:)
       xch%yco = yco_o
       call UpdateTY_with_Tn(xch,T_o,force_substep=force_substep)
       G_o = T_o - xch%Tg
    end if
    T_n0 = T_o
    dT_K = T_o*1.d-3
    T_n = T_o + dT_K
    G0 = G_o
    G = G_o
    do k_nr=0, k_nr_max -1
       xch%Tg = T_o
       xch%ychem(:) = y_o(:)
       xch%yco = yco_o
       call UpdateTY_with_Tn(xch,T_n,force_substep=force_substep)
       G = T_n - xch%Tg 
       if (abs(G/T_o) > 1. .and. k_nr >= 2) then
          if (dbg_flg_prm == 1) &
               print '(A,/,I0,(1P4E15.7),A)',"(CoolSolverImplicit, NR) abs(G/T_o) >1. with k_nr >= 2",&
               k_nr, T_o, T_n, xch%Tg, G, "    =>     move to bisection method"
          exit
       end if
       deno = G - G0
       if(abs(deno) < 1.d-50) then
         deno = 1.d-50
         G = G0 - 1.d-50
       endif
       dT_K = G*(T_n0 - T_n)/deno
       T_n0 = T_n
       G0 = G
       T_n = T_n + dT_K 
       if(T_n > T_max) T_n = T_max
       if(T_n < T_min) T_n = T_min
       if (T_n == T_n0) then
          if (dbg_flg_prm == 1) &
               print '(A,/,I0,(1P4E15.7),A)',"(CoolSolverImplicit, NR) T_n == T_n0 occurs during NR iteration",&
               k_nr, T_o, T_n, xch%Tg, G, "    =>     move to bisection method"
          exit
       end if
       err_Tn = abs((T_n - T_n0)/T_n) 
       err_G = abs(G/T_o) 
       if (dbg_flg_prm == 1) then
          print '(A,I0,(1P4E15.7))', "(CoolSolverImplicit, NR) k_nr, T_o, T_K, T_n0, T_n: ",&
               k_nr, T_o, xch%Tg, T_n0, T_n
          print '(A,(1P3E15.7),/)', "G, err_Tn, err_G: ",G, err_Tn, err_G
       end if
       if (err_Tn < err_eps .and. err_G < err_eps) then
          return
       end if
    end do
    if (G_o>0.) then
       T_n0 = T_o

       do k_srch = 0, k_srch_max-1
          T_n1 = T_n0 
          T_n0 = T_n0*0.9
          xch%ychem(:) = y_o(:)
          xch%yco = yco_o
          xch%Tg = T_o
          call UpdateTY_with_Tn(xch,T_n0,force_substep=force_substep)
          G = T_n0 - xch%Tg
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec_lb) k_srch, T_o, T_n0, T_n1, G: ",k_srch, T_o, T_n0, T_n1, G
       end if
          if (G < 0.) then
             exit 
          end if
          if (T_n0 < T_min .or. k_srch == k_srch_max-1) then
             print *, "CLODE_IMP_bisec: lower_bound_not_found", T_o, T_n0, xch%Td, T_min
             print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
             print '(A,(1P10E15.7))', "y_init: ",y_o(:)
             print *, 'gfuv, dph, rcopd, rOII', gfuv, dph, xch%rcopd, xch%rOII
             print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
             print *, 'fd, EradIR, chi_d', xch%fd, xch%EradIR, xch%chi_d
             print *, "stopping..."
             if(k_srch == k_srch_max-1) then 
                stop 
             end if
          end if
       end do
    else
       T_n1 = T_o
       do k_srch = 0, k_srch_max-1
          T_n0 = T_n1 
          T_n1 = T_n1*1.1
          xch%ychem(:) = y_o(:)
          xch%yco = yco_o
          xch%Tg = T_o
          call UpdateTY_with_Tn(xch,T_n1,force_substep=force_substep)
          G = T_n1 - xch%Tg
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec_ub) k_srch, T_o, T_n0, T_n1, G: ",k_srch, T_o, T_n0, T_n1, G
       end if
          if (G > 0.) then
             exit 
          end if
          if (T_n1 > T_max .or. k_srch == k_srch_max-1) then
             print *, "CLODE_IMP_bisec: upper_bound_not_found", T_o, T_n1, T_max
             print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
             print '(A,(1P10E15.7))', "y_init: ",y_o(:)
             print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
             print *, "stopping..."
             stop 
          end if
       end do
    end if
    do k_bs = 0, k_bs_max-1
       T_n = 0.5*(T_n1 + T_n0)
       xch%ychem(:) = y_o(:)
       xch%yco = yco_o
       xch%Tg = T_o
       call UpdateTY_with_Tn(xch,T_n,force_substep=force_substep)
       G = T_n - xch%Tg
       if (G < 0.) then
          T_n0 = T_n 
       else
          T_n1 = T_n 
       end if
       err_Tn = abs((T_n1-T_n0)/T_n0) 
       err_G = abs(G/T_o) 
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec) k_bs, T_o, T_K, err_Tn, err_G: ", k_bs, T_o, xch%Tg, err_Tn, err&
&_G
          print '(A,(1P5E20.12))', " T_n0, T_n1, T_n, T_K, G: ", T_n0, T_n1, T_n, xch%Tg, G
       end if
       if (err_Tn < err_eps) then
          xch%Tg = T_n
          return
       end if
       if (k_bs == k_bs_max-1) then
          print '(A,(1P3E15.7))', "CLODE_IMP_bisec: hitting the iteration limit",T_o, T_n0, T_n1
          print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
          print '(A,(1P10E15.7))', "y_init: ",y_o(:)
          print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
          print *, "stopping..."
          stop 
       end if
    end do
  end subroutine CoolSolverImplicit
  subroutine UpdateTY_with_Tn(xch, T_n, force_substep)
    type(chem_vals) :: xch
    real(kind=8),intent(IN) :: T_n
    logical,intent(IN),optional :: force_substep
    real(kind=8),dimension(0:25 -1) :: xk 
    real(kind=8) :: t_chem,xLmbd_tot,xmu,gamma,T_new,pr,en,t_cool, t_chemcool
    real(kind=8) :: rho, radius, Tdust, Tg
    logical :: bool_force_substep
    bool_force_substep = .false.
    if ( present( force_substep) ) then
       bool_force_substep = force_substep
    endif
    Tg = xch%Tg
    rho = (MP_mu*cgs_amu)*xch%nH 
    radius = 0.5d0*xch%xlmbdj
    call dtemp_radtr(xch%nH, Tg, xch%Td, xch%EradIR, xch%chi_d &
    , dph_in=xch%rdph)
    xch%Tg = T_n
    call react_coef(xch, xk)
    xmu = get_xmu(xch%ychem)
    pr = xch%nH*MP_mu/xmu * boltz * Tg 
    gamma = 1.d0+(1.d0+4.d0*yHe) & 
          /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(Tg)*xch%ychem(1)))
    en = pr / (gamma-1.) 
    if (dbg_flg_prm == 1) then
       print '(A,(1P7E15.7))', "(BEFORE chemreact) T_K,T_n,en,pr,xmu,gamma,xnH: ",Tg,T_n,en,pr,xmu,gamma,xch%nH 
    end if
    if(xch%nH > 1e11 .or. bool_force_substep) then 
       call chemreact_adptv(xch, t_chem, force_substep=.True., xk_in=xk) 
    else
       call chemreact_adptv(xch, t_chem, xk_in=xk)
    end if
    call update_yco(xch%nH,T_n,xch%ychem,xch%yco,xch%rcopd,xch%metal,xch%dt)
    call tot_cool(xch,xLmbd_tot,xk_in=xk)
    en = en - xLmbd_tot * xch%dt
    xmu = get_xmu(xch%ychem)
    gamma = 1.d0+(1.d0+4.d0*yHe) & 
          /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(Tg)*xch%ychem(1)))
    pr = (gamma-1.) * en
    T_new = pr / (xch%nH*MP_mu/xmu * boltz) 
    xch%Tg = T_new
    if (dbg_flg_prm == 1) then
       print '(A,(1P7E15.7))', "(AFTER chemreact) T_K,T_n,en,pr,xmu,gamma,xnH: ",xch%Tg,T_n,en,pr,xmu,gamma,xch%nH 
    end if
    t_cool = abs(en/xLmbd_tot) 
    t_chemcool = min(t_cool,t_chem) 
  end subroutine UpdateTY_with_Tn
  subroutine HUVAHHH(T, huv, alpha, heat, hve, hhm)
    real(kind=8),intent(IN) :: T
    real(kind=8),intent(OUT) :: huv, alpha, heat, hve, hhm
    real(kind=8),dimension(54) :: tbhuv,tbalph,tbheat,tbhved,tbhmds
    real(kind=8) :: tl1,tl2,dtl,tlg
    integer ::nt,it,it1
    data tbhuv /-5.21206d+01,-4.52041d+01,-3.90139d+01,-3.34752d+01,&
         -2.85209d+01,-2.40910d+01,-2.01313d+01,-1.65935d+01,-1.34341d+01,&
         -1.06141d+01,-8.09864d+00,-5.85626d+00,-3.85882d+00,-2.08108d+00,&
         -5.00370d-01, 9.03649d-01, 2.14923d+00, 3.25274d+00, 4.22889d+00,&
         5.09087d+00, 5.85053d+00, 6.51851d+00, 7.10437d+00, 7.61669d+00,&
         8.06318d+00, 8.45079d+00, 8.78575d+00, 9.07370d+00, 9.31969d+00,&
         9.52832d+00, 9.70371d+00, 9.84963d+00, 9.96948d+00, 1.00663d+01,&
         1.01430d+01, 1.02020d+01, 1.02455d+01, 1.02757d+01, 1.02943d+01,&
         1.03028d+01, 1.03028d+01, 1.02954d+01, 1.02817d+01, 1.02627d+01,&
         1.02390d+01, 1.02116d+01, 1.01808d+01, 1.01473d+01, 1.01115d+01,&
         1.00737d+01, 1.00344d+01, 9.99361d+00, 9.95174d+00, 9.90893d+00/
    data tbalph /-1.72075d+01,-1.72084d+01,-1.72093d+01,-1.72103d+01,&
         -1.72115d+01,-1.72128d+01,-1.72142d+01,-1.72158d+01,-1.72176d+01,&
         -1.72196d+01,-1.72218d+01,-1.72243d+01,-1.72271d+01,-1.72301d+01,&
         -1.72335d+01,-1.72373d+01,-1.72415d+01,-1.72462d+01,-1.72515d+01,&
         -1.72573d+01,-1.72637d+01,-1.72709d+01,-1.72789d+01,-1.72877d+01,&
         -1.72976d+01,-1.73085d+01,-1.73206d+01,-1.73339d+01,-1.73487d+01,&
         -1.73651d+01,-1.73831d+01,-1.74030d+01,-1.74248d+01,-1.74486d+01,&
         -1.74747d+01,-1.75030d+01,-1.75336d+01,-1.75667d+01,-1.76022d+01,&
         -1.76403d+01,-1.76808d+01,-1.77239d+01,-1.77695d+01,-1.78176d+01,&
         -1.78681d+01,-1.79210d+01,-1.79762d+01,-1.80337d+01,-1.80935d+01,&
         -1.81554d+01,-1.82193d+01,-1.82852d+01,-1.83530d+01,-1.84226d+01/
    data tbheat /-1.28617d+01,-1.28148d+01,-1.27678d+01,-1.27209d+01,&
         -1.26740d+01,-1.26272d+01,-1.25803d+01,-1.25336d+01,-1.24868d+01,&
         -1.24401d+01,-1.23935d+01,-1.23469d+01,-1.23004d+01,-1.22540d+01,&
         -1.22076d+01,-1.21613d+01,-1.21151d+01,-1.20690d+01,-1.20230d+01,&
         -1.19772d+01,-1.19314d+01,-1.18858d+01,-1.18404d+01,-1.17951d+01,&
         -1.17500d+01,-1.17051d+01,-1.16604d+01,-1.16160d+01,-1.15718d+01,&
         -1.15278d+01,-1.14842d+01,-1.14410d+01,-1.13981d+01,-1.13558d+01,&
         -1.13140d+01,-1.12729d+01,-1.12325d+01,-1.11928d+01,-1.11539d+01,&
         -1.11159d+01,-1.10788d+01,-1.10426d+01,-1.10074d+01,-1.09731d+01,&
         -1.09397d+01,-1.09073d+01,-1.08758d+01,-1.08452d+01,-1.08155d+01,&
         -1.07867d+01,-1.07587d+01,-1.07315d+01,-1.07052d+01,-1.06796d+01/
    data tbhved /-4.00350E+01,-3.43797E+01,-2.93209E+01,-2.47971E+01,&
         -2.07532E+01,-1.71399E+01,-1.39128E+01,-1.10321E+01,-8.46216E+00,&
         -6.17092E+00,-4.12969E+00,-2.31271E+00,-6.96923E-01, 7.38349E-01,&
         2.01159E+00, 3.13933E+00, 4.13631E+00, 5.01571E+00, 5.78929E+00,&
         6.46755E+00, 7.05991E+00, 7.57479E+00, 8.01978E+00, 8.40170E+00,&
         8.72672E+00, 9.00041E+00, 9.22781E+00, 9.41352E+00, 9.56170E+00,&
         9.67618E+00, 9.76043E+00, 9.81768E+00, 9.85086E+00, 9.86267E+00,&
         9.85558E+00, 9.83184E+00, 9.79345E+00, 9.74225E+00, 9.67985E+00,&
         9.60767E+00, 9.52697E+00, 9.43886E+00, 9.34428E+00, 9.24407E+00,&
         9.13894E+00, 9.02954E+00, 8.91638E+00, 8.79995E+00, 8.68064E+00,&
         8.55880E+00, 8.43475E+00, 8.30875E+00, 8.18102E+00, 8.05178E+00/
    data tbhmds / 1.1613E+01, 1.2113E+01, 1.2571E+01, 1.2992E+01,&
        1.3380E+01, 1.3739E+01, 1.4070E+01, 1.4376E+01, 1.4662E+01,&
        1.4927E+01, 1.5175E+01, 1.5407E+01, 1.5624E+01, 1.5829E+01,&
        1.6023E+01, 1.6205E+01, 1.6379E+01, 1.6544E+01, 1.6701E+01,&
        1.6851E+01, 1.6995E+01, 1.7134E+01, 1.7267E+01, 1.7395E+01,&
        1.7519E+01, 1.7639E+01, 1.7756E+01, 1.7869E+01, 1.7979E+01,&
        1.8087E+01, 1.8191E+01, 1.8293E+01, 1.8391E+01, 1.8487E+01,&
        1.8579E+01, 1.8668E+01, 1.8754E+01, 1.8837E+01, 1.8917E+01,&
        1.8994E+01, 1.9068E+01, 1.9139E+01, 1.9208E+01, 1.9275E+01,&
        1.9340E+01, 1.9404E+01, 1.9465E+01, 1.9525E+01, 1.9584E+01,&
        1.9642E+01, 1.9698E+01, 1.9754E+01, 1.9808E+01, 1.9862E+01/
    data tl1,tl2,dtl,nt/ 3.0000d+00, 5.5000d+00, 4.716981132d-02,54/
    tlg = (min(max(tl1,log10(t+1.0e-37)),tl2)-tl1)/dtl
    it = min(int(tlg)+1,nt)
    it1 = min(it+1,nt)
    tlg = tlg-int(tlg)
    huv = 10**(tbhuv(it)*(1.d0-tlg) + tbhuv(it1)*tlg)
    alpha = 10**(tbalph(it)*(1.d0-tlg) + tbalph(it1)*tlg)
    heat = 10**(tbheat(it)*(1.d0-tlg) + tbheat(it1)*tlg)
    hve = 10**(tbhved(it)*(1.d0-tlg) + tbhved(it1)*tlg)
    hhm = 10**(tbhmds(it)*(1.d0-tlg) + tbhmds(it1)*tlg)
  end subroutine HUVAHHH
  subroutine ProstFit2(Mass, Mdot, tage, Rad, Lum, Trad, rda)
    use mpilib
    use string, only : concat
    use io_util, only : read_env
    real(kind=8),intent(IN) :: Mass 
    real(kind=8),intent(IN) :: Mdot 
    real(kind=8),intent(IN) :: tage 
    real(kind=8),intent(Out) :: Rad 
    real(kind=8),intent(Out) :: Lum 
    real(kind=8),intent(Out) :: Trad 
    type(rad_others), intent(Out), optional :: rda
    integer,parameter :: nmd=11 
    integer,parameter :: npt=250 
    character(100) :: fname(0:nmd-1) = (/"1e-6", "1e-5", "1e-4", "3e-4", "1e-3", "3e-3" &
      , "6e-3", "1e-2", "3e-2", "6e-2", "1e-1"/)
    character(100) :: path2fitdat = "./ProstFit/"//trim("Z1")//"/"
    character(len=128) :: filename, dummy, dir, file
    integer,parameter :: FH = 11, FH2 = 12
    character(len=100) :: ffn
    integer :: err
    real(kind=8),save :: xm(0:npt-1) 
    real(kind=8),save :: xr(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xls(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xtrad(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xxeuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xxfuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xaleuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xheat(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xhhm(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xlumeuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xlumfuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xsigeuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xsigfuv(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xrOII(0:npt-1,0:nmd-1) 
    real(kind=8),parameter :: xmd(0:nmd-1) = (/1.d-6, 1.d-5, 1.d-4, 3.d-4, 1.d-3, 3.d-3, 6.d-3,&
         1.d-2, 3.d-2, 6.d-2, 1.d-1/) 
    integer :: i, j, id, ix, iy, n
    real(kind=8) :: x_in, x_in0, y_in, x1, x2, y1, y2, z_y1, z_y2, hi_R, hi_L
    real(kind=8) :: tap1, tap2, tap3, tap4, mstar, tlife, tmin, tm_1, tm_2, tp_1, tp_2
    real(kind=8), dimension(:,:), allocatable :: buf
    integer,save :: ifirst=0
    real(kind=8),save :: mass_max, dlm_log
    type stellar_evol
      integer :: nmodel
      real(kind=8), dimension(:), allocatable :: Age 
      real(kind=8), dimension(:), allocatable :: Rstar 
      real(kind=8), dimension(:), allocatable :: Lumi 
      real(kind=8), dimension(:), allocatable :: Teff 
      real(kind=8), dimension(:), allocatable :: xeuv 
      real(kind=8), dimension(:), allocatable :: xfuv 
      real(kind=8), dimension(:), allocatable :: alpha_euv 
      real(kind=8), dimension(:), allocatable :: heat_euv 
      real(kind=8), dimension(:), allocatable :: lumeuv 
      real(kind=8), dimension(:), allocatable :: lumfuv 
      real(kind=8), dimension(:), allocatable :: sigd_euv 
      real(kind=8), dimension(:), allocatable :: sigd_fuv 
      real(kind=8), dimension(:), allocatable :: rOII 
    end type stellar_evol
    type (stellar_evol), allocatable, dimension(:), save :: StarTable
    real(kind=8), allocatable, dimension(:), save :: stmass
    real(kind=8) :: tap_m1, tap_m2, tap_p1, tap_p2, ym, yp
    integer, save :: nstarmdl
    integer :: index_m, index_p, nmodel, nmodel_m, nmodel_p, index_tm, index_tp
    logical :: check_acc
    if (ifirst == 0) then
       allocate(buf(0:npt-1, 14))
       do id = 0, nmd-1
          ffn= trim(path2fitdat)//trim(fname(id))
          if(get_myrank() == 0) then
            print *, "PROST_FIT: reading ", ffn
            open(FH, file=ffn,status='old',iostat=err)
            if (err > 0) then
              print '(A,/,A)', "PROST_FIT: file not found","stopping..."
              stop
            end if
            do i = 0, npt-1
              read(FH,*,iostat=err) buf(i, 1:14)
              if (err > 0) then
                print *, "PROST_FIT: **WARNING** data size might be inconsistent"
                exit 
              end if
            end do
            close(FH)
          endif
          call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          xm(0:npt-1) = buf(0:npt-1, 1)
          xr(0:npt-1,id) = buf(0:npt-1, 2)
          xls(0:npt-1,id) = buf(0:npt-1, 3)
          xtrad(0:npt-1,id) = buf(0:npt-1, 4)
          xxeuv(0:npt-1,id) = buf(0:npt-1, 5)
          xxfuv(0:npt-1,id) = buf(0:npt-1, 6)
          xaleuv(0:npt-1,id) = buf(0:npt-1, 7)
          xheat(0:npt-1,id) = buf(0:npt-1, 8)
          xhhm(0:npt-1,id) = buf(0:npt-1, 9)
          xlumeuv(0:npt-1,id) = buf(0:npt-1,10)
          xlumfuv(0:npt-1,id) = buf(0:npt-1,11)
          xsigeuv(0:npt-1,id) = buf(0:npt-1,12)
          xsigfuv(0:npt-1,id) = buf(0:npt-1,13)
          xrOII(0:npt-1,id) = buf(0:npt-1,14)
       end do
       deallocate(buf)
       mass_max = 10.d0**xm(npt-1)
       if(get_myrank() == 0) then
          call read_env('DIR', dir)
          filename = trim("ProstFit/tsevolv_0.0.dat")
          file=concat(dir,filename)
          open(FH2, file=file)
          read(FH2, fmt=*) nstarmdl 
       endif
       call mpi_bcast(nstarmdl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       allocate(StarTable(nstarmdl), stmass(nstarmdl))
       do i=1, nstarmdl
          if(get_myrank() == 0) then
            read(FH2, fmt=*) dummy
            read(FH2, fmt=*) mstar, nmodel
          endif
          call mpi_bcast(mstar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          call mpi_bcast(nmodel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          stmass(i) = mstar 
          StarTable(i)%nmodel = nmodel
          allocate(StarTable(i)%Age(nmodel) &
                  ,StarTable(i)%Rstar(nmodel) &
                  ,StarTable(i)%Lumi(nmodel) &
                  ,StarTable(i)%Teff(nmodel) &
                  ,StarTable(i)%xeuv(nmodel) &
                  ,StarTable(i)%xfuv(nmodel) &
                  ,StarTable(i)%alpha_euv(nmodel) &
                  ,StarTable(i)%heat_euv(nmodel) &
                  ,StarTable(i)%lumeuv(nmodel) &
                  ,StarTable(i)%lumfuv(nmodel) &
                  ,StarTable(i)%sigd_euv(nmodel) &
                  ,StarTable(i)%sigd_fuv(nmodel) &
                  ,StarTable(i)%rOII(nmodel))
          allocate(buf(13, nmodel))
          if(get_myrank() == 0) then
            do n=1, nmodel
              read(FH2, fmt=*) buf(1:13, n)
            enddo
          endif
          call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          do n=1, nmodel
            StarTable(i)%Age(n) = buf(1 , n) 
            StarTable(i)%Rstar(n) = buf(2 , n) 
            StarTable(i)%Lumi(n) = buf(3 , n) 
            StarTable(i)%Teff(n) = buf(4 , n) 
            StarTable(i)%xeuv(n) = buf(5 , n) 
            StarTable(i)%xfuv(n) = buf(6 , n) 
            StarTable(i)%alpha_euv(n) = buf(7 , n) 
            StarTable(i)%heat_euv(n) = buf(8 , n) 
            StarTable(i)%lumeuv(n) = buf(9 , n) 
            StarTable(i)%lumfuv(n) = buf(10, n) 
            StarTable(i)%sigd_euv(n) = buf(11, n) 
            StarTable(i)%sigd_fuv(n) = buf(12, n) 
            StarTable(i)%rOII(n) = buf(13, n) 
          enddo
          deallocate(buf)
       enddo
       if(get_myrank() == 0) then
         close(FH2)
       endif
       dlm_log = (log10(stmass(nstarmdl))-log10(stmass(1)))/dble(nstarmdl-1)
       ifirst = 1
    end if
    check_acc = .False. 
    if(Mass < stmass(1)) then
      index_m = 1
      check_acc = .True.
    elseif (Mass > stmass(nstarmdl)) then
      index_m = nstarmdl-1
      check_acc = .True.
    else
      index_m = int((log10(Mass)-log10(stmass(1)))/dlm_log)+1
      if (Mdot > xmd(0)) then 
        check_acc = .True.
      endif
    endif
    index_p = index_m+1
    nmodel_m = StarTable(index_m)%nmodel
    nmodel_p = StarTable(index_p)%nmodel
    tlife = min(StarTable(index_m)%Age(nmodel_m), StarTable(index_p)%Age(nmodel_p))
    tmin = max(StarTable(index_m)%Age(1), StarTable(index_p)%Age(1))
    if(tage >= tlife) then
      Rad = 0.d0
      Lum = 0.d0
      Trad = 0.d0
      rda%xeuv = 0.d0
      rda%xfuv = 0.d0
      rda%alpha_euv = 0.d0
      rda%heat_euv = 0.d0
      rda%hhm = 0.d0
      rda%lumeuv = 0.d0
      rda%lumfuv = 0.d0
      rda%sig_euv = 0.d0
      rda%sig_fuv = 0.d0
      rda%rOII = 0.d0
      return
    endif
    if(tage < tmin) check_acc = .True.
    if(check_acc) then
      x_in0 = log10(Mass)
      x_in = min(max(x_in0, xm(0)), xm(npt-1))
      y_in = min(max(Mdot,xmd(0)), xmd(nmd-1))
      ix = int((x_in+2.0)/0.02d0)
      ix = max(min(ix, npt-2),0)
      iy = nmd-2
      do j = 0, nmd-2
        if (y_in <= xmd(j+1)) then
          iy = j
          exit
        endif
      enddo
      if(x_in0 > xm(npt-1)) then
        hi_L = (Mass/mass_max)
        hi_R = hi_L**0.5d0
      else
        hi_L = 1.d0
        hi_R = 1.d0
      endif
      x1 = xm(ix)
      x2 = xm(ix+1)
      y1 = xmd(iy)
      y2 = xmd(iy+1)
      tap1 = (x2-x_in)/(x2-x1)
      tap2 = (x_in-x1)/(x2-x1)
      tap3 = (y2-y_in)/(y2-y1)
      tap4 = (y_in-y1)/(y2-y1)
      z_y1 = tap1*xr(ix,iy) + tap2*xr(ix+1,iy) 
 z_y2 = tap1*xr(ix,iy+1)+ tap2*xr(ix+1,iy+1) 
 Rad = tap3*z_y1 + tap4*z_y2
      Rad = Rad * hi_R
      z_y1 = tap1*xls(ix,iy) + tap2*xls(ix+1,iy) 
 z_y2 = tap1*xls(ix,iy+1)+ tap2*xls(ix+1,iy+1) 
 Lum = tap3*z_y1 + tap4*z_y2
      Lum = Lum * hi_L
      z_y1 = tap1*xtrad(ix,iy) + tap2*xtrad(ix+1,iy) 
 z_y2 = tap1*xtrad(ix,iy+1)+ tap2*xtrad(ix+1,iy+1) 
 Trad = tap3*z_y1 + tap4*z_y2
      if(present(rda)) then
        z_y1 = tap1*xxeuv(ix,iy) + tap2*xxeuv(ix+1,iy) 
 z_y2 = tap1*xxeuv(ix,iy+1)+ tap2*xxeuv(ix+1,iy+1) 
 rda%xeuv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xxfuv(ix,iy) + tap2*xxfuv(ix+1,iy) 
 z_y2 = tap1*xxfuv(ix,iy+1)+ tap2*xxfuv(ix+1,iy+1) 
 rda%xfuv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xaleuv(ix,iy) + tap2*xaleuv(ix+1,iy) 
 z_y2 = tap1*xaleuv(ix,iy+1)+ tap2*xaleuv(ix+1,iy+1) 
 rda%alpha_euv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xheat(ix,iy) + tap2*xheat(ix+1,iy) 
 z_y2 = tap1*xheat(ix,iy+1)+ tap2*xheat(ix+1,iy+1) 
 rda%heat_euv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xhhm(ix,iy) + tap2*xhhm(ix+1,iy) 
 z_y2 = tap1*xhhm(ix,iy+1)+ tap2*xhhm(ix+1,iy+1) 
 rda%hhm = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xlumeuv(ix,iy) + tap2*xlumeuv(ix+1,iy) 
 z_y2 = tap1*xlumeuv(ix,iy+1)+ tap2*xlumeuv(ix+1,iy+1) 
 rda%lumeuv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xlumfuv(ix,iy) + tap2*xlumfuv(ix+1,iy) 
 z_y2 = tap1*xlumfuv(ix,iy+1)+ tap2*xlumfuv(ix+1,iy+1) 
 rda%lumfuv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xsigeuv(ix,iy) + tap2*xsigeuv(ix+1,iy) 
 z_y2 = tap1*xsigeuv(ix,iy+1)+ tap2*xsigeuv(ix+1,iy+1) 
 rda%sig_euv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xsigfuv(ix,iy) + tap2*xsigfuv(ix+1,iy) 
 z_y2 = tap1*xsigfuv(ix,iy+1)+ tap2*xsigfuv(ix+1,iy+1) 
 rda%sig_fuv = tap3*z_y1 + tap4*z_y2
        z_y1 = tap1*xrOII(ix,iy) + tap2*xrOII(ix+1,iy) 
 z_y2 = tap1*xrOII(ix,iy+1)+ tap2*xrOII(ix+1,iy+1) 
 rda%rOII = tap3*z_y1 + tap4*z_y2
      endif
      return
    else
      index_tm = nmodel_m-1
      do n = 1, nmodel_m-1
        if(tage < StarTable(index_m)%Age(n+1)) then
          index_tm = n
          exit
        endif
      enddo
      index_tp = nmodel_p-1
      do n = 1, nmodel_p-1
        if(tage < StarTable(index_p)%Age(n+1)) then
          index_tp = n
          exit
        endif
      enddo
      x1 = log10(stmass(index_m))
      x2 = log10(stmass(index_p))
      x_in = log10(Mass)
      tm_1 = StarTable(index_m)%Age(index_tm)
      tm_2 = StarTable(index_m)%Age(index_tm+1)
      tp_1 = StarTable(index_p)%Age(index_tp)
      tp_2 = StarTable(index_p)%Age(index_tp+1)
      tap_m1 = (tm_2-tage)/(tm_2-tm_1)
      tap_m2 = (tage-tm_1)/(tm_2-tm_1)
      tap_p1 = (tp_2-tage)/(tp_2-tp_1)
      tap_p2 = (tage-tp_1)/(tp_2-tp_1)
      tap1 = (x2-x_in)/(x2-x1)
      tap2 = (x_in-x1)/(x2-x1)
        ym = StarTable(index_m)%Rstar(index_tm)*tap_m1+StarTable(index_m)%Rstar(index_tm+1)*tap_m2
 yp = StarTable(index_p)%Rstar(index_tp)*tap_p1+StarTable(index_p)%Rstar(index_tp+1)*tap_p2
 Rad = ym*tap1 + yp*tap2
        ym = StarTable(index_m)%Lumi(index_tm)*tap_m1+StarTable(index_m)%Lumi(index_tm+1)*tap_m2
 yp = StarTable(index_p)%Lumi(index_tp)*tap_p1+StarTable(index_p)%Lumi(index_tp+1)*tap_p2
 Lum = ym*tap1 + yp*tap2
        ym = StarTable(index_m)%Teff(index_tm)*tap_m1+StarTable(index_m)%Teff(index_tm+1)*tap_m2
 yp = StarTable(index_p)%Teff(index_tp)*tap_p1+StarTable(index_p)%Teff(index_tp+1)*tap_p2
 Trad = ym*tap1 + yp*tap2
        if(present(rda)) then
          ym = StarTable(index_m)%xeuv(index_tm)*tap_m1+StarTable(index_m)%xeuv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%xeuv(index_tp)*tap_p1+StarTable(index_p)%xeuv(index_tp+1)*tap_p2
 rda%xeuv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%xfuv(index_tm)*tap_m1+StarTable(index_m)%xfuv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%xfuv(index_tp)*tap_p1+StarTable(index_p)%xfuv(index_tp+1)*tap_p2
 rda%xfuv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%alpha_euv(index_tm)*tap_m1+StarTable(index_m)%alpha_euv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%alpha_euv(index_tp)*tap_p1+StarTable(index_p)%alpha_euv(index_tp+1)*tap_p2
 rda%alpha_euv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%heat_euv(index_tm)*tap_m1+StarTable(index_m)%heat_euv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%heat_euv(index_tp)*tap_p1+StarTable(index_p)%heat_euv(index_tp+1)*tap_p2
 rda%heat_euv = ym*tap1 + yp*tap2
          rda%hhm = 0.d0
          ym = StarTable(index_m)%lumeuv(index_tm)*tap_m1+StarTable(index_m)%lumeuv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%lumeuv(index_tp)*tap_p1+StarTable(index_p)%lumeuv(index_tp+1)*tap_p2
 rda%lumeuv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%lumfuv(index_tm)*tap_m1+StarTable(index_m)%lumfuv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%lumfuv(index_tp)*tap_p1+StarTable(index_p)%lumfuv(index_tp+1)*tap_p2
 rda%lumfuv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%sigd_euv(index_tm)*tap_m1+StarTable(index_m)%sigd_euv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%sigd_euv(index_tp)*tap_p1+StarTable(index_p)%sigd_euv(index_tp+1)*tap_p2
 rda%sig_euv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%sigd_fuv(index_tm)*tap_m1+StarTable(index_m)%sigd_fuv(index_tm+1)*tap_m2
 yp = StarTable(index_p)%sigd_fuv(index_tp)*tap_p1+StarTable(index_p)%sigd_fuv(index_tp+1)*tap_p2
 rda%sig_fuv = ym*tap1 + yp*tap2
          ym = StarTable(index_m)%rOII(index_tm)*tap_m1+StarTable(index_m)%rOII(index_tm+1)*tap_m2
 yp = StarTable(index_p)%rOII(index_tp)*tap_p1+StarTable(index_p)%rOII(index_tp+1)*tap_p2
 rda%rOII = ym*tap1 + yp*tap2
        endif
    endif
  end subroutine ProstFit2
  subroutine ludcmp(A, indx, d)
    real(kind=8), dimension(0:6 -1, 0:6 -1),intent(INOUT):: A
    integer, dimension(0:6 -1),intent(OUT) :: indx
    real(kind=8), intent(OUT) :: d
    integer :: i, j, kk, imax_1 = 0
    real(kind=8) :: big,dum,summ, tempora
    real(kind=8), dimension(0:6 -1) :: vv
    imax_1 = 0
    tempora= 0.d0
    vv = 0.d0
    d = 1.d0
    do i=0, 6 -1
      big = 0.d0
      do j=0, 6 -1
        tempora = dabs(A(i,j))
        if (tempora > big) big = tempora
      enddo
      if(big == 0.d0) then
        print *, "Singular matrix in rourine ludcmp"
        print *, "stopping ...."
        stop
      endif
      vv(i) = 1.d0/big
    enddo
    do j = 0, 6 -1
      do i=0, j-1
        summ = A(i,j)
        do kk=0, i-1
          summ = summ - A(i,kk)*A(kk,j)
        enddo
        A(i,j) = summ
      enddo
      big = 0.d0
      do i=j, 6 -1
        summ = A(i,j)
        do kk=0, j-1
          summ = summ - A(i, kk)*A(kk,j)
        enddo
        A(i,j) = summ
        dum = vv(i)*dabs(summ)
        if(dum >= big) then
          big = dum
          imax_1 = i
        endif
      enddo
      if (j .ne. imax_1) then
        do kk = 0, 6 -1
          dum = A(imax_1, kk)
          A(imax_1, kk) = A(j,kk)
          A(j,kk) = dum
        enddo
        d = -d
        vv(imax_1) = vv(j)
      endif
      indx(j) = imax_1
      if( j .ne. 6 -1) then
        dum = 1.d0 / (A(j,j))
        do i = j+1, 6 -1
          A(i, j) = A(i, j)*dum
        enddo
      endif
    enddo
  end subroutine ludcmp
  subroutine lubksb(A, indx, b)
    real(kind=8), dimension(0:6 -1, 0:6 -1),intent(IN):: A
    integer, dimension(0:6 -1),intent(IN) :: indx
    real(kind=8), dimension(0:6 -1), intent(INOUT) :: b
    integer :: i, ii, ip, j
    real(kind=8) :: summ
    ii = 0
    do i=0, 6 -1
      ip = indx(i)
      summ = b(ip)
      b(ip)= b(i)
      if (ii .ne. 0) then
        do j=ii-1, i-1
          summ = summ - a(i,j)*b(j)
        enddo
      else if (summ .ne. 0.d0) then
        ii = i+1
      endif
      b(i) = summ
    enddo
    do i = 6 -1, 0, -1
      summ = b(i)
      do j=i+1, 6 -1
        summ = summ - a(i,j)*b(j)
      enddo
      b(i) = summ/a(i,i)
    enddo
  end subroutine lubksb
  function get_xmu(ychem)
    real(kind=8), dimension(0:6 -1), intent(IN) :: ychem
    real(kind=8) :: get_xmu, abund_tot, weight_tot
    integer :: ichem
    abund_tot = 0.d0
    weight_tot = 0.d0
    abund_tot = ychem(0)+ ychem(1)+ ychem(2)+ ychem(3)+yHe
    get_xmu = MP_mu / abund_tot
  end function get_xmu
end module primordial
