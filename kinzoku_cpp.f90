module kinzoku
  use unit
  use parameter
  use mpilib
  use modelParameter
  implicit none
  private
  real(kind=8),parameter :: tau_max_dust = 1d2 
  real(kind=8),parameter :: frac_AHe = 0.0972 
  real(kind=8),parameter :: CONST_AHe = 4.004 
  real(kind=8),parameter :: MP_FISRF = 1.6d-3 
  real(kind=8),parameter :: Qd_uv = 1.d0 
  real(kind=8),parameter :: Qd_ir = 1.d-2 
  real(kind=8),parameter :: a_dust = 1.d-5 
  real(kind=8),parameter :: rho_dust = 3.d0 
  real(kind=8),parameter :: fdust_solar = 0.01d0 
  real(kind=8),save :: NcnH_max
  real(kind=8),parameter :: g1_CII = 2.0
  real(kind=8),parameter :: g2_CII = 4.0
  real(kind=8),parameter :: delE_21_CII = 1.259d-14 
  real(kind=8),save :: nu_21_CII, delT_21_CII, A21_CII, delE_21_kB_CII
  real(kind=8),save :: delT_21_OI, delT_31_OI, delT_41_OI, delT_51_OI
  real(kind=8),save :: delT_32_OI, delT_42_OI, delT_52_OI, delT_43_OI, delT_53_OI, delT_54_OI
  real(kind=8),save :: gamma41H_OI, gamma51H_OI
  real(kind=8),save :: gamma21e_OI, gamma31e_OI, gamma41e_OI, gamma51e_OI
  real(kind=8),save :: gamma42H_OI, gamma52H_OI
  real(kind=8),save :: gamma32e_OI, gamma42e_OI, gamma52e_OI
  real(kind=8),save :: gamma43H_OI, gamma53H_OI
  real(kind=8),save :: gamma43e_OI, gamma53e_OI
  real(kind=8),save :: gamma54H_OI, gamma54e_OI
  real(kind=8),save,dimension(1:5) :: g_OI, inv_g_OI
  real(kind=8),save,dimension(1:5,1:5) :: A_OI, delE_OI, nu_OI
  real(kind=8),save,dimension(0:99999) :: data_dp, data_temp, data_ros
  real(kind=8),save :: excite_const_OII, g0_OII, g1_OII, g2_OII, DT10_OII, DT20_OII, DT21_OII &
    ,DE10_OII, DE20_OII, DE21_OII, nu10_OII, nu20_OII, nu21_OII, A10_OII, A20_OII, A21_OII
  real(kind=8),save,dimension(1:5,1:5) :: A_OIII, delE_OIII, nu_OIII
  real(kind=8),save,dimension(1:5) :: g_OIII, inv_g_OIII
  real(kind=8),save :: delT_21_OIII, delT_31_OIII, delT_41_OIII, delT_51_OIII, delT_32_OIII &
    , delT_42_OIII, delT_52_OIII, delT_43_OIII, delT_53_OIII, delT_54_OIII
  real(kind=8),parameter :: MP_Tcmb = 2.73
  public :: init_kinzoku
  public :: frac_AHe, CONST_AHe, MP_FISRF
  public :: NcnH_max
  public :: oxygenI_cooling, carbonII_cooling, oxygenII_cooling, oxygenIII_cooling
  public :: PhotoelectricHeating, depsilon, OII_OII_ratio
  public :: update_yco, COcool, COdissociation_rate
  public :: find_dop, fdust_solar, find_dros, dtemp_radtr
contains
  subroutine init_kinzoku
    use io_util, only : readenv, print_msg, read_env
    use unit
    implicit none
    integer i
    real(kind=8) :: Td
    NcnH_max = tau_max_dust / (1.d-21) 
    NcnH_max = max(NcnH_max, tau_max_dust / (5.34d-22*MP_Metallicity*2.5)) 
    if(get_myrank() == 0) then
        print *, 'NcnH_max            =', NcnH_max
    endif
    nu_21_CII = delE_21_CII/cgs_hpk 
    delT_21_CII = delE_21_CII/cgs_kb 
    A21_CII = 2.4e-6
    delE_21_kB_CII = delE_21_CII / cgs_kb
    g_OI(1) = 5.0 
    g_OI(2) = 3.0 
    g_OI(3) = 1.0 
    g_OI(4) = 5.0 
    g_OI(5) = 1.0 
    do i=1, 5
        inv_g_OI(i) = 1.e0/g_OI(i)
    end do
    gamma41H_OI = 1.0e-12
    gamma42H_OI = 1.0e-12
    gamma43H_OI = 1.0e-12
    gamma51H_OI = 1.0e-12
    gamma52H_OI = 1.0e-12
    gamma53H_OI = 1.0e-12
    gamma54H_OI = 0.0
    gamma21e_OI = 1.4e-8
    gamma31e_OI = 1.4e-8
    gamma32e_OI = 5.0e-9
    gamma41e_OI = 1.0e-10
    gamma42e_OI = 1.0e-10
    gamma43e_OI = 1.0e-10
    gamma51e_OI = 1.0e-10
    gamma52e_OI = 1.0e-10
    gamma53e_OI = 1.0e-10
    gamma54e_OI = 0.0
    A_OI(2,1) = 8.9e-5
    A_OI(3,1) = 1.0e-10
    A_OI(3,2) = 1.7e-5
    A_OI(4,1) = 6.3e-3
    A_OI(4,2) = 2.1e-3
    A_OI(4,3) = 7.3e-7
    A_OI(5,1) = 2.9e-4
    A_OI(5,2) = 7.3e-2
    A_OI(5,3) = 0.0 
    A_OI(5,4) = 1.2
    delE_OI(2,1) = 3.144e-14 
    delE_OI(3,2) = 1.365e-14 
    delE_OI(4,3) = 3.14e-12 
    delE_OI(5,4) = 3.56e-12 
    delE_OI(3,1) = delE_OI(3,2) + delE_OI(2,1)
    delE_OI(4,1) = delE_OI(4,3) + delE_OI(3,2) + delE_OI(2,1)
    delE_OI(5,1) = delE_OI(5,4) + delE_OI(4,3) + delE_OI(3,2) + delE_OI(2,1)
    delE_OI(4,2) = delE_OI(4,3) + delE_OI(3,2)
    delE_OI(5,2) = delE_OI(5,4) + delE_OI(4,3) + delE_OI(3,2)
    delE_OI(5,3) = delE_OI(5,4) + delE_OI(4,3)
    delT_21_OI = delE_OI(2,1) / cgs_kb
    delT_31_OI = delE_OI(3,1) / cgs_kb
    delT_32_OI = delE_OI(3,2) / cgs_kb
    delT_41_OI = delE_OI(4,1) / cgs_kb
    delT_42_OI = delE_OI(4,2) / cgs_kb
    delT_43_OI = delE_OI(4,3) / cgs_kb
    delT_51_OI = delE_OI(5,1) / cgs_kb
    delT_52_OI = delE_OI(5,2) / cgs_kb
    delT_53_OI = delE_OI(5,3) / cgs_kb
    delT_54_OI = delE_OI(5,4) / cgs_kb
    nu_OI(2,1) = delE_OI(2,1) / cgs_hpk
    nu_OI(3,1) = delE_OI(3,1) / cgs_hpk
    nu_OI(3,2) = delE_OI(3,2) / cgs_hpk
    nu_OI(4,1) = delE_OI(4,1) / cgs_hpk
    nu_OI(4,2) = delE_OI(4,2) / cgs_hpk
    nu_OI(4,3) = delE_OI(4,3) / cgs_hpk
    nu_OI(5,1) = delE_OI(5,1) / cgs_hpk
    nu_OI(5,2) = delE_OI(5,2) / cgs_hpk
    nu_OI(5,3) = delE_OI(5,3) / cgs_hpk
    nu_OI(5,4) = delE_OI(5,4) / cgs_hpk
    excite_const_OII = cgs_hpk**2.e0 / (2.e0*Pi*cgs_me)**1.5d0 / (cgs_kb)**0.5
    g0_OII = 4.d0
    g1_OII = 6.d0
    g2_OII = 4.d0
    DT10_OII = 3.9d4
    DT20_OII = 3.9d4
    DT21_OII = 30.d0
    DE10_OII = DT10_OII*cgs_kb
    DE20_OII = DT20_OII*cgs_kb
    DE21_OII = DT21_OII*cgs_kb
    nu10_OII = DE10_OII / cgs_hpk
    nu20_OII = DE20_OII / cgs_hpk
    nu21_OII = DE21_OII / cgs_hpk
    A10_OII = 5.1d-5
    A20_OII = 1.7d-4
    A21_OII = 1.3d-7
    g_OIII(1) = 1.d0 
    g_OIII(2) = 3.d0 
    g_OIII(3) = 5.d0 
    g_OIII(4) = 5.d0 
    g_OIII(5) = 1.d0 
    do i=1, 5
        inv_g_OIII(i) = 1.d0 / g_OIII(i)
    enddo
    A_OIII(5,4) = 1.8
    A_OIII(5,3) = 7.8d-4
    A_OIII(5,2) = 2.2d-1
    A_OIII(5,1) = 0.d0
    A_OIII(4,3) = 2.d-2
    A_OIII(4,2) = 6.7d-3
    A_OIII(4,1) = 2.7d-6
    A_OIII(3,2) = 9.8d-5
    A_OIII(3,1) = 3.d-11
    A_OIII(2,1) = 2.6d-5
    delE_OIII(2,1) = 2.247177d-14 
    delE_OIII(3,2) = 3.834950d-14
    delE_OIII(4,3) = 3.967530d-12
    delE_OIII(5,4) = 4.55290d-12
    delE_OIII(3,1) = delE_OIII(3,2) + delE_OIII(2,1)
    delE_OIII(4,1) = delE_OIII(4,3) + delE_OIII(3,2) + delE_OIII(2,1)
    delE_OIII(5,1) = delE_OIII(5,4) + delE_OIII(4,3) + delE_OIII(3,2) + delE_OIII(2,1)
    delE_OIII(4,2) = delE_OIII(4,3) + delE_OIII(3,2)
    delE_OIII(5,2) = delE_OIII(5,4) + delE_OIII(4,3) + delE_OIII(3,2)
    delE_OIII(5,3) = delE_OIII(5,4) + delE_OIII(4,3)
    delT_21_OIII = delE_OIII(2,1) / cgs_kb
    delT_31_OIII = delE_OIII(3,1) / cgs_kb
    delT_32_OIII = delE_OIII(3,2) / cgs_kb
    delT_41_OIII = delE_OIII(4,1) / cgs_kb
    delT_42_OIII = delE_OIII(4,2) / cgs_kb
    delT_43_OIII = delE_OIII(4,3) / cgs_kb
    delT_51_OIII = delE_OIII(5,1) / cgs_kb
    delT_52_OIII = delE_OIII(5,2) / cgs_kb
    delT_53_OIII = delE_OIII(5,3) / cgs_kb
    delT_54_OIII = delE_OIII(5,4) / cgs_kb
    nu_OIII(2,1) = delE_OIII(2,1) / cgs_hpk
    nu_OIII(3,1) = delE_OIII(3,1) / cgs_hpk
    nu_OIII(3,2) = delE_OIII(3,2) / cgs_hpk
    nu_OIII(4,1) = delE_OIII(4,1) / cgs_hpk
    nu_OIII(4,2) = delE_OIII(4,2) / cgs_hpk
    nu_OIII(4,3) = delE_OIII(4,3) / cgs_hpk
    nu_OIII(5,1) = delE_OIII(5,1) / cgs_hpk
    nu_OIII(5,2) = delE_OIII(5,2) / cgs_hpk
    nu_OIII(5,3) = delE_OIII(5,3) / cgs_hpk
    nu_OIII(5,4) = delE_OIII(5,4) / cgs_hpk
    call read_plankop
  end subroutine init_kinzoku
  subroutine carbonII_cooling(nCII, ne, nH, T_gas, T_rad, dvdr, b_cont, cool_rate_CII)
    implicit none
    real(kind=8),intent(IN) :: nCII, ne, nH, T_gas, T_rad
    real(kind=8), dimension(0:2), intent(IN) :: dvdr, b_cont
    real(kind=8),intent(OUT) :: cool_rate_CII
    real(kind=8) :: gamma12H, gamma12e, gamma21H, gamma21e
    real(kind=8) :: C21,C12, R21, R12, f1,f2, Q21
    real(kind=8) :: T100, exp_21, coef, beta21
    integer :: idirec
    real(kind=8), dimension(0:2) :: tau21, esc21
    Q21 = compute_Q_bg(T_gas, T_rad)
    gamma21H = 8.0e-10 * (T_gas*1e-2)**(0.07)
    gamma21e = 2.8e-7 * (T_gas*1e-2)**(-0.5)
    exp_21 = dexp(-delE_21_kB_CII/(T_gas))
    gamma12H = (g2_CII/g1_CII) * gamma21H * exp_21
    gamma12e = (g2_CII/g1_CII) * gamma21e * exp_21
    C21 = gamma21H * nH + gamma12e * ne
    C12 = g2_CII/g1_CII * C21 * exp_21
    R21 = A21_CII*(1.0+Q21) + C21
    R12 = g2_CII/g1_CII * A21_CII*Q21 + C12
    f1 = R21/(R21+R12)
    f2 = R12/(R21+R12)
    coef = (A21_CII/8.0/Pi)*(cgs_c/nu_21_CII)**3*(f1*g2_CII/g1_CII-f2)*nCII
    do idirec=0, 2
        tau21(idirec)=coef/dvdr(idirec)
        if(tau21(idirec)<300.0) then
            esc21(idirec) = compute_b_esc(tau21(idirec),b_cont(idirec))
        else
            esc21(idirec) = 0.0
        end if
        if(esc21(idirec)>1.0) then
            esc21(idirec) = 1.0
        end if
    end do
    beta21 = (esc21(0)+esc21(1)+esc21(2))/3.0
    cool_rate_CII = A21_CII * delE_21_CII * beta21*(f2-Q21*(g2_CII/g1_CII*f1-f2))
    cool_rate_CII = max(0.0, cool_rate_CII)
  contains
    function compute_Q_bg(Tnu, Trad)
        implicit none
        real(kind=8) :: compute_Q_bg
        real(kind=8),intent(IN) :: Tnu, Trad
        real(kind=8) :: xx, Q_bg
        xx = Tnu/Trad
        if(xx>100.e0) then
            compute_Q_bg = 0.e0
        else
            compute_Q_bg = 1.0/(dexp(xx)-1.0)
        endif
    end function compute_Q_bg
    function compute_b_esc(tau_line, b_cont)
        implicit none
        real(kind=8),intent(IN) :: tau_line, b_cont
        real(kind=8) :: compute_b_esc
            if(tau_line<0.0) then
                compute_b_esc = 1e0
            else if(tau_line<1.0e-5) then
                compute_b_esc = b_cont
            else
                compute_b_esc = (1.0 - dexp(-tau_line)) / tau_line * b_cont
            endif
    end function compute_b_esc
  end subroutine carbonII_cooling
  subroutine oxygenI_cooling(nOI, ne, nH, T_gas, T_rad, dvdr, b_cont, cool_rate_OI)
    implicit none
    real(kind=8),intent(IN) :: nOI, ne, nH, T_gas, T_rad
    real(kind=8), dimension(0:2), intent(IN) :: dvdr, b_cont
    real(kind=8),intent(OUT) :: cool_rate_OI
    real(kind=8) :: R21, R31, R41, R51, R32, R42, R52, R43, R53, R54
    real(kind=8) :: R12, R13, R14, R15, R23, R24, R25, R34, R35, R45
    real(kind=8), dimension(1:5, 1:5) :: Q
    real(kind=8), dimension(0:4, 0:4) :: Coeff
    real(kind=8), dimension(0:4, 0:0) :: b
    real(kind=8), dimension(1:5) :: f
    real(kind=8) :: gamma21H_OI,gamma31H_OI,gamma32H_OI, inv_T
    real(kind=8) :: exp_31,exp_21,exp_41,exp_51,exp_32,exp_42
    real(kind=8) :: exp_52,exp_43,exp_53,exp_54
    real(kind=8) :: gamma12H_OI, gamma13H_OI, gamma14H_OI, gamma15H_OI
    real(kind=8) :: gamma12e_OI, gamma13e_OI, gamma14e_OI, gamma15e_OI
    real(kind=8) :: gamma23H_OI, gamma24H_OI, gamma25H_OI
    real(kind=8) :: gamma23e_OI, gamma24e_OI, gamma25e_OI
    real(kind=8) :: gamma34H_OI, gamma35H_OI
    real(kind=8) :: gamma34e_OI, gamma35e_OI
    real(kind=8) :: gamma45H_OI
    real(kind=8) :: gamma45e_OI
    integer :: idirec,upp,low
    real(kind=8) :: coef, beta
    real(kind=8), dimension(0:2) :: tau, esc
    inv_T = 1.e0 / T_gas
    gamma21H_OI = 9.2e-11 * (T_gas*1.e-2)**(0.67)
    gamma31H_OI = 4.3e-11 * (T_gas*1.e-2)**(0.8)
    gamma32H_OI = 1.1e-10 * (T_gas*1.e-2)**(0.44)
    Q(2,1) = compute_Q_bg(delT_21_OI, T_rad)
    Q(3,1) = compute_Q_bg(delT_31_OI, T_rad)
    Q(3,2) = compute_Q_bg(delT_32_OI, T_rad)
    Q(4,1) = compute_Q_bg(delT_41_OI, T_rad)
    Q(4,2) = compute_Q_bg(delT_42_OI, T_rad)
    Q(4,3) = compute_Q_bg(delT_43_OI, T_rad)
    Q(5,1) = compute_Q_bg(delT_51_OI, T_rad)
    Q(5,2) = compute_Q_bg(delT_52_OI, T_rad)
    Q(5,3) = compute_Q_bg(delT_53_OI, T_rad)
    Q(5,4) = compute_Q_bg(delT_54_OI, T_rad)
    exp_21 = dexp(-delT_21_OI*inv_T)
    exp_31 = dexp(-delT_31_OI*inv_T)
    exp_41 = dexp(-delT_41_OI*inv_T)
    exp_51 = dexp(-delT_51_OI*inv_T)
    exp_32 = dexp(-delT_32_OI*inv_T)
    exp_42 = dexp(-delT_42_OI*inv_T)
    exp_52 = dexp(-delT_52_OI*inv_T)
    exp_43 = dexp(-delT_43_OI*inv_T)
    exp_53 = dexp(-delT_53_OI*inv_T)
    exp_54 = dexp(-delT_54_OI*inv_T)
    gamma12H_OI = (g_OI(2)*inv_g_OI(1)) * gamma21H_OI * exp_21
    gamma13H_OI = (g_OI(3)*inv_g_OI(1)) * gamma31H_OI * exp_31
    gamma14H_OI = (g_OI(4)*inv_g_OI(1)) * gamma41H_OI * exp_41
    gamma15H_OI = (g_OI(5)*inv_g_OI(1)) * gamma51H_OI * exp_51
    gamma23H_OI = (g_OI(3)*inv_g_OI(2)) * gamma32H_OI * exp_32
    gamma24H_OI = (g_OI(4)*inv_g_OI(2)) * gamma42H_OI * exp_42
    gamma25H_OI = (g_OI(5)*inv_g_OI(2)) * gamma52H_OI * exp_52
    gamma34H_OI = (g_OI(4)*inv_g_OI(3)) * gamma43H_OI * exp_43
    gamma35H_OI = (g_OI(5)*inv_g_OI(3)) * gamma53H_OI * exp_53
    gamma45H_OI = (g_OI(5)*inv_g_OI(4)) * gamma54H_OI * exp_54
    gamma12e_OI = (g_OI(2)*inv_g_OI(1)) * gamma21e_OI * exp_21
    gamma13e_OI = (g_OI(3)*inv_g_OI(1)) * gamma31e_OI * exp_31
    gamma14e_OI = (g_OI(4)*inv_g_OI(1)) * gamma41e_OI * exp_41
    gamma15e_OI = (g_OI(5)*inv_g_OI(1)) * gamma51e_OI * exp_51
    gamma23e_OI = (g_OI(3)*inv_g_OI(2)) * gamma32e_OI * exp_32
    gamma24e_OI = (g_OI(4)*inv_g_OI(2)) * gamma42e_OI * exp_42
    gamma25e_OI = (g_OI(5)*inv_g_OI(2)) * gamma52e_OI * exp_52
    gamma34e_OI = (g_OI(4)*inv_g_OI(3)) * gamma43e_OI * exp_43
    gamma35e_OI = (g_OI(5)*inv_g_OI(3)) * gamma53e_OI * exp_53
    gamma45e_OI = (g_OI(5)*inv_g_OI(4)) * gamma54e_OI * exp_54
    R21 = A_OI(2,1)*(1.0+Q(2,1)) + gamma21H_OI * nH + gamma21e_OI * ne
    R31 = A_OI(3,1)*(1.0+Q(3,1)) + gamma31H_OI * nH + gamma31e_OI * ne
    R41 = A_OI(4,1)*(1.0+Q(4,1)) + gamma41H_OI * nH + gamma41e_OI * ne
    R51 = A_OI(5,1)*(1.0+Q(5,1)) + gamma51H_OI * nH + gamma51e_OI * ne
    R32 = A_OI(3,2)*(1.0+Q(3,2)) + gamma32H_OI * nH + gamma32e_OI * ne
    R42 = A_OI(4,2)*(1.0+Q(4,2)) + gamma42H_OI * nH + gamma42e_OI * ne
    R52 = A_OI(5,2)*(1.0+Q(5,2)) + gamma52H_OI * nH + gamma52e_OI * ne
    R43 = A_OI(4,3)*(1.0+Q(4,3)) + gamma43H_OI * nH + gamma43e_OI * ne
    R53 = A_OI(5,3)*(1.0+Q(5,3)) + gamma53H_OI * nH + gamma53e_OI * ne
    R54 = A_OI(5,4)*(1.0+Q(5,4)) + gamma54H_OI * nH + gamma54e_OI * ne
    R12 = (g_OI(2)*inv_g_OI(1))*A_OI(2,1)*Q(2,1) + gamma12H_OI * nH + gamma12e_OI * ne
    R13 = (g_OI(3)*inv_g_OI(1))*A_OI(3,1)*Q(3,1) + gamma13H_OI * nH + gamma13e_OI * ne
    R14 = (g_OI(4)*inv_g_OI(1))*A_OI(4,1)*Q(4,1) + gamma14H_OI * nH + gamma14e_OI * ne
    R15 = (g_OI(5)*inv_g_OI(1))*A_OI(5,1)*Q(5,1) + gamma15H_OI * nH + gamma15e_OI * ne
    R23 = (g_OI(3)*inv_g_OI(2))*A_OI(3,2)*Q(3,2) + gamma23H_OI * nH + gamma23e_OI * ne
    R24 = (g_OI(4)*inv_g_OI(2))*A_OI(4,2)*Q(4,2) + gamma24H_OI * nH + gamma24e_OI * ne
    R25 = (g_OI(5)*inv_g_OI(2))*A_OI(5,2)*Q(5,2) + gamma25H_OI * nH + gamma25e_OI * ne
    R34 = (g_OI(4)*inv_g_OI(3))*A_OI(4,3)*Q(4,3) + gamma34H_OI * nH + gamma34e_OI * ne
    R35 = (g_OI(5)*inv_g_OI(3))*A_OI(5,3)*Q(5,3) + gamma35H_OI * nH + gamma35e_OI * ne
    R45 = (g_OI(5)*inv_g_OI(4))*A_OI(5,4)*Q(5,4) + gamma45H_OI * nH + gamma45e_OI * ne
    Coeff(0,0) = 1.0
    Coeff(0,1) = 1.0
    Coeff(0,2) = 1.0
    Coeff(0,3) = 1.0
    Coeff(0,4) = 1.0
    Coeff(1,0) = R12
    Coeff(1,1) = -(R21 + R23 + R24 + R25)
    Coeff(1,2) = R32
    Coeff(1,3) = R42
    Coeff(1,4) = R52
    Coeff(2,0) = R13
    Coeff(2,1) = R23
    Coeff(2,2) = -(R31 + R32 + R34 + R35)
    Coeff(2,3) = R43
    Coeff(2,4) = R53
    Coeff(3,0) = R14
    Coeff(3,1) = R24
    Coeff(3,2) = R34
    Coeff(3,3) = -(R41 + R42 + R43 + R45)
    Coeff(3,4) = R54
    Coeff(4,0) = R15
    Coeff(4,1) = R25
    Coeff(4,2) = R35
    Coeff(4,3) = R45
    Coeff(4,4) = -(R51 + R52 + R53 + R54)
    b(0,0) = 1.0
    b(1,0) = 0
    b(2,0) = 0
    b(3,0) = 0
    b(4,0) = 0
    call gaussj5(Coeff, 5, b, 1)
    f(1) = b(0,0)
    f(2) = b(1,0)
    f(3) = b(2,0)
    f(4) = b(3,0)
    f(5) = b(4,0)
    cool_rate_OI = 0.0
    do upp=1, 5
        do low=1, upp-1
            coef = (A_OI(upp,low)/8.0/Pi)*(cgs_c/nu_OI(upp,low))**(3.0)*(g_OI(upp)/g_OI(low)*f(low)-f(upp))*nOI
            do idirec=0, 2
                tau(idirec)=coef/dvdr(idirec)
                if(tau(idirec)<300.0) then
                    esc(idirec) = compute_b_esc(tau(idirec),b_cont(idirec))
                else
                    esc(idirec) = 0.0
                end if
                if(esc(idirec)>1.0) then
                    esc(idirec) = 1.0
                end if
            end do
            beta = (esc(0)+esc(1)+esc(2))/3.0
            cool_rate_OI = cool_rate_OI + A_OI(upp,low) * delE_OI(upp,low) * beta*(f(upp)-Q(upp,low)*(g_OI(upp)/g_OI(low)*f(low)-f(&
&upp)))
        end do
    end do
    cool_rate_OI = max(cool_rate_OI, 0.0)
  contains
    function compute_Q_bg(Tnu, Trad)
        implicit none
        real(kind=8) :: compute_Q_bg
        real(kind=8),intent(IN) :: Tnu, Trad
        real(kind=8) :: xx, Q_bg
        xx = Tnu/Trad
        if(xx>100.e0) then
            compute_Q_bg = 0.e0
        else
            compute_Q_bg = 1.0/(dexp(xx)-1.0)
        endif
    end function compute_Q_bg
    function compute_b_esc(tau_line, b_cont)
        implicit none
        real(kind=8),intent(IN) :: tau_line, b_cont
        real(kind=8) :: compute_b_esc
            if(tau_line<0.0) then
                compute_b_esc = 1e0
            else if(tau_line<1.0e-5) then
                compute_b_esc = b_cont
            else
                compute_b_esc = (1.0 - dexp(-tau_line)) / tau_line * b_cont
            endif
    end function compute_b_esc
    subroutine gaussj5(a, n, b, m)
    implicit none
    real(kind=8), dimension(0:4,0:4), intent(INOUT) :: a
    real(kind=8), dimension(0:4,0:0), intent(INOUT) :: b
    real(kind=8) :: xyz
    integer, intent(IN) :: n, m
    integer :: i,icol,irow,j,k,l,ll
    real(kind=8) :: big,dum,pivinv,temp
    integer, dimension(0:4) :: indxc, indxr, ipiv
    ipiv = 0.d0
    do i = 0, n-1
        big=0.0
        do j=0, n-1
            if (ipiv(j) .ne. 1) then
                do k=0, n-1
                    if (ipiv(k) .eq. 0) then
                        if (abs(a(j,k)) >= big) then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    endif
                enddo
            endif
        enddo
        ipiv(icol) = ipiv(icol) + 1
        if (irow .ne. icol) then
            do l=0, n-1
                temp = a(irow,l)
 a(irow,l) = a(icol,l)
 a(icol,l) = temp
            enddo
            do l=0, m-1
                temp = b(irow,l)
 b(irow,l) = b(icol,l)
 b(icol,l) = temp
            enddo
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol) .eq. 0.d0) then
            print *, "Error gaussj5 Singular martix"
            stop
        endif
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do l = 0, n-1
            a(icol,l) = a(icol,l) * pivinv
        enddo
        do l = 0, m-1
            b(icol,l) = b(icol,l) * pivinv
        enddo
        do ll = 0, n-1
            if (ll .ne. icol) then
                dum = a(ll,icol)
                a(ll,icol)=0.d0
                do l = 0, n-1
                    a(ll,l) = a(ll,l) - a(icol,l)*dum
                enddo
                do l = 0, m-1
                    b(ll,l) = b(ll,l) - b(icol,l)*dum
                enddo
            endif
        enddo
    enddo
    do l = n - 1, 0, -1
        if (indxr(l) .ne. indxc(l)) then
            do k=0,n-1
                temp = a(k,indxr(l))
 a(k,indxr(l)) = a(k,indxc(l))
 a(k,indxc(l)) = temp
            enddo
        endif
    enddo
    end subroutine gaussj5
    end subroutine oxygenI_cooling
    subroutine oxygenII_cooling(nOII, ne, T_gas, T_rad, dvdr, b_cont, cool_rate_OII)
        implicit none
        real(kind=8), intent(IN) :: nOII, ne, T_gas, T_rad
        real(kind=8), dimension(0:2), intent(IN) :: dvdr, b_cont
        real(kind=8), intent(OUT) :: cool_rate_OII
        real(kind=8) :: T4, T_sq, lnT4, T4_sq
        real(kind=8) :: Q10, Q20, Q21, gamma10e, gamma20e, gamma21e
        real(kind=8) :: C10, C20, C21, C01, C02, C12, R10, R20, R21, R01, R02, R12 &
          , f0, f1, f2
        integer idirec
        real(kind=8) :: coef, beta10,beta20,beta21
        real(kind=8), dimension(0:2) :: tau10, tau20, tau21, esc10, esc20, esc21
        T4 = T_gas * 1.d-4
        T_sq = (T_gas)**0.5d0
        lnT4 = dlog(T4)
        T4_sq= (T4)**(-0.5d0)
        Q10 = compute_Q_bg(DT10_OII, T_rad)
        Q20 = compute_Q_bg(DT20_OII, T_rad)
        Q21 = compute_Q_bg(DT21_OII, T_rad)
        gamma10e = 1.3d-8 * T4_sq
        gamma20e = 1.3d-8 * T4_sq
        gamma21e = 2.5d-8 * T4_sq
        C10 = ne*gamma10e
        C20 = ne*gamma20e
        C21 = ne*gamma21e
        C01 = g1_OII/g0_OII * C10 * dexp(-DT10_OII/T_gas)
        C02 = g2_OII/g0_OII * C20 * dexp(-DT20_OII/T_gas)
        C12 = g2_OII/g1_OII * C21 * dexp(-DT21_OII/T_gas)
        R10 = A10_OII*(1.0+Q10) + C10
        R20 = A20_OII*(1.0+Q20) + C20
        R21 = A21_OII*(1.0+Q21) + C21
        R01 = g1_OII/g0_OII*A10_OII*Q10 + C01
        R02 = g2_OII/g0_OII*A20_OII*Q20 + C02
        R12 = g2_OII/g1_OII*A21_OII*Q21 + C12
        f0 = ( R21*(R10-R20) + R20*(R10+R12+R21) )/ ( (R01+R02+R20)*(R10+R12+R21) - (R01-R21)*(R10-R20) )
        f1 = (f0*(R01-R21)+R21)/(R10+R12+R21)
        f2 = (f0*R02+f1*R12)/(R21+R20)
        cool_rate_OII = 0.d0
        coef = (A10_OII/8.0/Pi)*(cgs_c/nu10_OII)**3.d0*(f0*g1_OII/g0_OII-f1)*nOII
        do idirec=0, 2
            tau10(idirec)=coef/dvdr(idirec)
            if(tau10(idirec)<300.0) then
                esc10(idirec) = b_esc(tau10(idirec), b_cont(idirec))
            else
                esc10(idirec) = 0.d0
            endif
            if(esc10(idirec) > 1.d0) esc10(idirec) = 1.d0
        enddo
        beta10 = (esc10(0)+esc10(1)+esc10(2))/3.d0
        cool_rate_OII = cool_rate_OII + A10_OII * DE10_OII * beta10*(f1-Q10*(g1_OII/g0_OII*f0-f1))
        coef = (A20_OII/8.0/Pi)*(cgs_c/nu20_OII)**3.d0*(f0*g2_OII/g0_OII-f2)*nOII
        do idirec=0, 2
            tau20(idirec)=coef/dvdr(idirec)
            if(tau20(idirec)<300.0) then
                esc20(idirec) = b_esc(tau20(idirec), b_cont(idirec))
            else
                esc20(idirec) = 0.d0
            endif
            if(esc20(idirec) > 1.d0) esc20(idirec) = 1.d0
        enddo
        beta20 = (esc20(0)+esc20(1)+esc20(2))/3.d0
        cool_rate_OII = cool_rate_OII + A20_OII * DE20_OII * beta20*(f2-Q20*(g2_OII/g0_OII*f0-f2))
        coef = (A21_OII/8.0/Pi)*(cgs_c/nu21_OII)**3.d0*(f1*g2_OII/g1_OII-f2)*nOII
        do idirec=0, 2
            tau21(idirec)=coef/dvdr(idirec)
            if(tau21(idirec)<300.0) then
                esc21(idirec) = b_esc(tau21(idirec), b_cont(idirec))
            else
                esc21(idirec) = 0.d0
            endif
            if(esc21(idirec) > 1.d0) esc21(idirec) = 1.d0
        enddo
        beta21 = (esc21(0)+esc21(1)+esc21(2))/3.d0
        cool_rate_OII = cool_rate_OII + A21_OII * DE21_OII * beta21*(f2-Q21*(g2_OII/g1_OII*f1-f2))
  contains
    function compute_Q_bg(Tnu, Trad)
        implicit none
        real(kind=8) :: compute_Q_bg
        real(kind=8),intent(IN) :: Tnu, Trad
        real(kind=8) :: xx, Q_bg
        xx = Tnu/Trad
        if(xx>100.e0) then
            compute_Q_bg = 0.e0
        else
            compute_Q_bg = 1.0/(dexp(xx)-1.0)
        endif
    end function compute_Q_bg
    function b_esc(tau_line, b_cont)
        implicit none
        real(kind=8),intent(IN) :: tau_line, b_cont
        real(kind=8) :: b_esc
        if(tau_line < 0.d0) then
            b_esc = 1.d0
        else if (tau_line < 1.d-5) then
            b_esc = b_cont
        else
            b_esc = (1.0 - dexp(-tau_line)) / tau_line * b_cont
        endif
    end function b_esc
    end subroutine oxygenII_cooling
    subroutine oxygenIII_cooling(nOIII, ne, nH, T_gas, T_rad, dvdr, b_cont, cool_rate_OIII)
        implicit none
        real(kind=8), intent(IN) :: nOIII, ne, T_gas, T_rad, nH
        real(kind=8), dimension(0:2), intent(IN) :: dvdr, b_cont
        real(kind=8), intent(OUT) :: cool_rate_OIII
        real(kind=8) :: gamma21e, gamma31e, gamma41e, gamma51e, gamma32e, gamma42e, gamma52e &
          , gamma43e, gamma53e, gamma54e, gamma12H, gamma13H, gamma14H, gamma15H, gamma12e, gamma13e &
          , gamma14e, gamma15e, gamma23H, gamma24H, gamma25H, gamma23e, gamma24e, gamma25e, gamma34H &
          , gamma35H, gamma34e, gamma35e, gamma45H, gamma45e
        real(kind=8) :: omega21, omega31, omega41, omega51, omega32, omega42, omega52 &
          , omega43, omega53, omega54, gamma_ce
        real(kind=8) :: R21, R31, R41, R51, R32, R42, R52, R43, R53, R54, &
          R12, R13, R14, R15, R23, R24, R25, R34, R35, R45
        real(kind=8), dimension(1:5,1:5) :: Q
        real(kind=8), dimension(0:4,0:4) :: Coeff
        real(kind=8), dimension(0:4, 0:0) :: b
        real(kind=8), dimension(1:5) :: f
        integer :: i
        real(kind=8) :: exp_31,exp_21,exp_41,exp_51,exp_32,exp_42,exp_52,exp_43,exp_53,exp_54
        real(kind=8) :: T4,lnT4,tap,inv_T,beta
        integer :: idirec,upp,low
        real(kind=8) :: coef
        real(kind=8), dimension(0:2) :: tau,esc
        real(kind=8) :: gamma21H, gamma31H, gamma32H, gamma41H, gamma42H, gamma43H, gamma51H, gamma52H, gamma53H, gamma54H
        gamma21H = 0.d0
        gamma31H = 0.d0
        gamma32H = 0.d0
        gamma41H = 0.d0
        gamma42H = 0.d0
        gamma43H = 0.d0
        gamma51H = 0.d0
        gamma52H = 0.d0
        gamma53H = 0.d0
        gamma54H = 0.d0
        T4 = T_gas*1.d-4
        lnT4 = dlog(T4)
        inv_T = 1.d0 / T_gas
        gamma_ce = 8.629d-8 / dsqrt(T4)
        omega54 = 0.523 * T4 ** (0.210-0.099*lnT4)
        tap = T4 ** (0.118+0.057*lnT4)
        omega53 = 0.1605*tap
        omega52 = 0.0963*tap
        omega51 = 0.0321*tap
        tap = T4 ** (0.120+0.031*lnT4)
        omega43 = 1.215*tap
        omega42 = 0.729*tap
        omega41 = 0.243*tap
        omega32 = 1.23 * T4 ** (0.053+0.007*lnT4)
        omega31 = 0.257* T4 ** (0.081+0.017*lnT4)
        omega21 = 0.522* T4 ** (0.033-0.009*lnT4)
        gamma54e = gamma_ce*omega54*inv_g_OIII(5)
        gamma53e = gamma_ce*omega53*inv_g_OIII(5)
        gamma52e = gamma_ce*omega52*inv_g_OIII(5)
        gamma51e = gamma_ce*omega51*inv_g_OIII(5)
        gamma43e = gamma_ce*omega43*inv_g_OIII(4)
        gamma42e = gamma_ce*omega42*inv_g_OIII(4)
        gamma41e = gamma_ce*omega41*inv_g_OIII(4)
        gamma32e = gamma_ce*omega32*inv_g_OIII(3)
        gamma31e = gamma_ce*omega31*inv_g_OIII(3)
        gamma21e = gamma_ce*omega32*inv_g_OIII(2)
        Q(2,1) = compute_Q_bg(delT_21_OIII, T_rad)
        Q(3,1) = compute_Q_bg(delT_31_OIII, T_rad)
        Q(3,2) = compute_Q_bg(delT_32_OIII, T_rad)
        Q(4,1) = compute_Q_bg(delT_41_OIII, T_rad)
        Q(4,2) = compute_Q_bg(delT_42_OIII, T_rad)
        Q(4,3) = compute_Q_bg(delT_43_OIII, T_rad)
        Q(5,1) = compute_Q_bg(delT_51_OIII, T_rad)
        Q(5,2) = compute_Q_bg(delT_52_OIII, T_rad)
        Q(5,3) = compute_Q_bg(delT_53_OIII, T_rad)
        Q(5,4) = compute_Q_bg(delT_54_OIII, T_rad)
        exp_21 = dexp(-delT_21_OIII*inv_T)
        exp_31 = dexp(-delT_31_OIII*inv_T)
        exp_41 = dexp(-delT_41_OIII*inv_T)
        exp_51 = dexp(-delT_51_OIII*inv_T)
        exp_32 = dexp(-delT_32_OIII*inv_T)
        exp_42 = dexp(-delT_42_OIII*inv_T)
        exp_52 = dexp(-delT_52_OIII*inv_T)
        exp_43 = dexp(-delT_43_OIII*inv_T)
        exp_53 = dexp(-delT_53_OIII*inv_T)
        exp_54 = dexp(-delT_54_OIII*inv_T)
        gamma12H = (g_OIII(2)*inv_g_OIII(1)) * gamma21H * exp_21
        gamma13H = (g_OIII(3)*inv_g_OIII(1)) * gamma31H * exp_31
        gamma14H = (g_OIII(4)*inv_g_OIII(1)) * gamma41H * exp_41
        gamma15H = (g_OIII(5)*inv_g_OIII(1)) * gamma51H * exp_51
        gamma23H = (g_OIII(3)*inv_g_OIII(2)) * gamma32H * exp_32
        gamma24H = (g_OIII(4)*inv_g_OIII(2)) * gamma42H * exp_42
        gamma25H = (g_OIII(5)*inv_g_OIII(2)) * gamma52H * exp_52
        gamma34H = (g_OIII(4)*inv_g_OIII(3)) * gamma43H * exp_43
        gamma35H = (g_OIII(5)*inv_g_OIII(3)) * gamma53H * exp_53
        gamma45H = (g_OIII(5)*inv_g_OIII(4)) * gamma54H * exp_54
        gamma12e = (g_OIII(2)*inv_g_OIII(1)) * gamma21e * exp_21
        gamma13e = (g_OIII(3)*inv_g_OIII(1)) * gamma31e * exp_31
        gamma14e = (g_OIII(4)*inv_g_OIII(1)) * gamma41e * exp_41
        gamma15e = (g_OIII(5)*inv_g_OIII(1)) * gamma51e * exp_51
        gamma23e = (g_OIII(3)*inv_g_OIII(2)) * gamma32e * exp_32
        gamma24e = (g_OIII(4)*inv_g_OIII(2)) * gamma42e * exp_42
        gamma25e = (g_OIII(5)*inv_g_OIII(2)) * gamma52e * exp_52
        gamma34e = (g_OIII(4)*inv_g_OIII(3)) * gamma43e * exp_43
        gamma35e = (g_OIII(5)*inv_g_OIII(3)) * gamma53e * exp_53
        gamma45e = (g_OIII(5)*inv_g_OIII(4)) * gamma54e * exp_54
        R21 = A_OIII(2,1)*(1.0+Q(2,1)) + gamma21H * nH + gamma21e * ne
        R31 = A_OIII(3,1)*(1.0+Q(3,1)) + gamma31H * nH + gamma31e * ne
        R41 = A_OIII(4,1)*(1.0+Q(4,1)) + gamma41H * nH + gamma41e * ne
        R51 = A_OIII(5,1)*(1.0+Q(5,1)) + gamma51H * nH + gamma51e * ne
        R32 = A_OIII(3,2)*(1.0+Q(3,2)) + gamma32H * nH + gamma32e * ne
        R42 = A_OIII(4,2)*(1.0+Q(4,2)) + gamma42H * nH + gamma42e * ne
        R52 = A_OIII(5,2)*(1.0+Q(5,2)) + gamma52H * nH + gamma52e * ne
        R43 = A_OIII(4,3)*(1.0+Q(4,3)) + gamma43H * nH + gamma43e * ne
        R53 = A_OIII(5,3)*(1.0+Q(5,3)) + gamma53H * nH + gamma53e * ne
        R54 = A_OIII(5,4)*(1.0+Q(5,4)) + gamma54H * nH + gamma54e * ne
        R12 = (g_OIII(2)*inv_g_OIII(1))*A_OIII(2,1)*Q(2,1) + gamma12H * nH + gamma12e * ne
        R13 = (g_OIII(3)*inv_g_OIII(1))*A_OIII(3,1)*Q(3,1) + gamma13H * nH + gamma13e * ne
        R14 = (g_OIII(4)*inv_g_OIII(1))*A_OIII(4,1)*Q(4,1) + gamma14H * nH + gamma14e * ne
        R15 = (g_OIII(5)*inv_g_OIII(1))*A_OIII(5,1)*Q(5,1) + gamma15H * nH + gamma15e * ne
        R23 = (g_OIII(3)*inv_g_OIII(2))*A_OIII(3,2)*Q(3,2) + gamma23H * nH + gamma23e * ne
        R24 = (g_OIII(4)*inv_g_OIII(2))*A_OIII(4,2)*Q(4,2) + gamma24H * nH + gamma24e * ne
        R25 = (g_OIII(5)*inv_g_OIII(2))*A_OIII(5,2)*Q(5,2) + gamma25H * nH + gamma25e * ne
        R34 = (g_OIII(4)*inv_g_OIII(3))*A_OIII(4,3)*Q(4,3) + gamma34H * nH + gamma34e * ne
        R35 = (g_OIII(5)*inv_g_OIII(3))*A_OIII(5,3)*Q(5,3) + gamma35H * nH + gamma35e * ne
        R45 = (g_OIII(5)*inv_g_OIII(4))*A_OIII(5,4)*Q(5,4) + gamma45H * nH + gamma45e * ne
        Coeff(0,0) = 1.0
        Coeff(0,1) = 1.0
        Coeff(0,2) = 1.0
        Coeff(0,3) = 1.0
        Coeff(0,4) = 1.0
        Coeff(1,0) = R12
        Coeff(1,1) = -(R21 + R23 + R24 + R25)
        Coeff(1,2) = R32
        Coeff(1,3) = R42
        Coeff(1,4) = R52
        Coeff(2,0) = R13
        Coeff(2,1) = R23
        Coeff(2,2) = -(R31 + R32 + R34 + R35)
        Coeff(2,3) = R43
        Coeff(2,4) = R53
        Coeff(3,0) = R14
        Coeff(3,1) = R24
        Coeff(3,2) = R34
        Coeff(3,3) = -(R41 + R42 + R43 + R45)
        Coeff(3,4) = R54
        Coeff(4,0) = R15
        Coeff(4,1) = R25
        Coeff(4,2) = R35
        Coeff(4,3) = R45
        Coeff(4,4) = -(R51 + R52 + R53 + R54)
        b(0,0) = 1.0
        b(1,0) = 0.0
        b(2,0) = 0.0
        b(3,0) = 0.0
        b(4,0) = 0.0
        call gaussj5(Coeff, 5, b, 1)
        f(1) = b(0,0)
        f(2) = b(1,0)
        f(3) = b(2,0)
        f(4) = b(3,0)
        f(5) = b(4,0)
        cool_rate_OIII = 0.0
        do upp=2, 5
            do low=1, upp-1
                coef = (A_OIII(upp,low)/8.0/Pi)*(cgs_c/nu_OIII(upp,low))**3.d0*(g_OIII(upp)/g_OIII(low)*f(low)-f(upp))*nOIII
                do idirec = 0, 2
                    tau(idirec) = coef / dvdr(idirec)
                    if(tau(idirec)<300.d0) then
                        esc(idirec) = compute_b_esc(tau(idirec),b_cont(idirec))
                    else
                        esc(idirec) = 0.e0
                    endif
                    if (tau(idirec)>1.d0) esc(idirec) = 1.e0
                    beta = (esc(0) + esc(1) + esc(2) )/3.d0
                enddo
                 cool_rate_OIII = cool_rate_OIII + A_OIII(upp,low) * delE_OIII(upp,low) * beta*(f(upp)-Q(upp,low) &
                        *(g_OIII(upp)/g_OIII(low)*f(low)-f(upp)))
            enddo
        enddo
        cool_rate_OIII = max(cool_rate_OIII, 0.d0)
  contains
    function compute_Q_bg(Tnu, Trad)
        implicit none
        real(kind=8) :: compute_Q_bg
        real(kind=8),intent(IN) :: Tnu, Trad
        real(kind=8) :: xx, Q_bg
        xx = Tnu/Trad
        if(xx>100.e0) then
            compute_Q_bg = 0.e0
        else
            compute_Q_bg = 1.0/(dexp(xx)-1.0)
        endif
    end function compute_Q_bg
    function compute_b_esc(tau_line, b_cont)
        implicit none
        real(kind=8),intent(IN) :: tau_line, b_cont
        real(kind=8) :: compute_b_esc
            if(tau_line<0.0) then
                compute_b_esc = 1e0
            else if(tau_line<1.0e-5) then
                compute_b_esc = b_cont
            else
                compute_b_esc = (1.0 - dexp(-tau_line)) / tau_line * b_cont
            endif
    end function compute_b_esc
    subroutine gaussj5(a, n, b, m)
    implicit none
    real(kind=8), dimension(0:4,0:4), intent(INOUT) :: a
    real(kind=8), dimension(0:4,0:0), intent(INOUT) :: b
    real(kind=8) :: xyz
    integer, intent(IN) :: n, m
    integer :: i,icol,irow,j,k,l,ll
    real(kind=8) :: big,dum,pivinv,temp
    integer, dimension(0:4) :: indxc, indxr, ipiv
    ipiv = 0.d0
    do i = 0, n-1
        big=0.0
        do j=0, n-1
            if (ipiv(j) .ne. 1) then
                do k=0, n-1
                    if (ipiv(k) .eq. 0) then
                        if (abs(a(j,k)) >= big) then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    endif
                enddo
            endif
        enddo
        ipiv(icol) = ipiv(icol) + 1
        if (irow .ne. icol) then
            do l=0, n-1
                temp = a(irow,l)
 a(irow,l) = a(icol,l)
 a(icol,l) = temp
            enddo
            do l=0, m-1
                temp = b(irow,l)
 b(irow,l) = b(icol,l)
 b(icol,l) = temp
            enddo
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol) .eq. 0.d0) then
            print *, "Error gaussj5 Singular martix"
            stop
        endif
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do l = 0, n-1
            a(icol,l) = a(icol,l) * pivinv
        enddo
        do l = 0, m-1
            b(icol,l) = b(icol,l) * pivinv
        enddo
        do ll = 0, n-1
            if (ll .ne. icol) then
                dum = a(ll,icol)
                a(ll,icol)=0.d0
                do l = 0, n-1
                    a(ll,l) = a(ll,l) - a(icol,l)*dum
                enddo
                do l = 0, m-1
                    b(ll,l) = b(ll,l) - b(icol,l)*dum
                enddo
            endif
        enddo
    enddo
    do l = n - 1, 0, -1
        if (indxr(l) .ne. indxc(l)) then
            do k=0,n-1
                temp = a(k,indxr(l))
 a(k,indxr(l)) = a(k,indxc(l))
 a(k,indxc(l)) = temp
            enddo
        endif
    enddo
    end subroutine gaussj5
    end subroutine oxygenIII_cooling
    subroutine dtemp_radtr(N_H, Tg, Td, Erad, chi_d, dph_in)
        use modelParameter, only : MP_Crd, MP_Ctil
        implicit none
        real(kind=8),intent(IN) :: N_H, Tg, Erad
        real(kind=8),intent(IN),optional :: dph_in
        real(kind=8),intent(INOUT) :: Td
        real(kind=8),intent(OUT) :: chi_d
        real(kind=8), parameter :: Rd = 5.83e-8 
        real(kind=8), parameter :: h = 1.d-5
        integer, parameter :: iteration_max = 20
        real(kind=8) :: col_g
        real(kind=8) :: gas_col, rad_heat, xd_rad, kapb_rad
        integer :: notconverge, iteration
        real(kind=8) :: derror, err_f, Trad
        real(kind=8) :: kappa, kapb_d, f, f_pl, Td_pl , max_v
        real(kind=8) :: dif, fu, fl, fm, Td_l, Td_m, Td_u, dTul, dT_lm, T_ini
        logical :: isNotFinite
        col_g = Rd * N_H * dsqrt(Tg*1.d-3)*(1.e0-0.8*dexp(-75.d0/Tg))
        Trad = (MP_Crd*Erad /cgs_asb)**(0.25d0)
        call find_dop(Trad, xd_rad)
        rad_heat = xd_rad*MP_Ctil*Erad 
        if (present(dph_in)) then
          rad_heat = rad_heat + dph_in 
        endif
        Td = Trad
        if(Td < MP_Tcmb) then
          Td = Tg
        endif
        T_ini = Td
        notconverge = 0
        iteration = 0
        derror = 1.d0
        err_f = 1.d0
        do while( derror > 1.d-5 )
            iteration = iteration + 1
            call find_dop(Td, kappa)
            kapb_d = kappa*cgs_asb*cgs_c*Td**4.d0
            f = col_g * (Tg - Td) + rad_heat - kapb_d
            Td_pl = Td*(1.d0+h)
            call find_dop(Td_pl, kappa)
            kapb_d = kappa*cgs_asb*cgs_c*Td_pl**4.d0
            f_pl = col_g * (Tg - Td_pl) + rad_heat - kapb_d
            dif = (f_pl - f) / (Td*h)
            derror = abs(-f/dif) / Td
            max_v = max(rad_heat, kapb_d)
            err_f = abs(f / max_v)
            if(derror < 0.5d0) then
                Td = Td - f / dif
            else
                Td = Td- f / dif * 0.5d0
            endif
            if(Td < 0.d0) then
                notconverge = 2
                exit
            end if
            if(iteration > iteration_max) then
                notconverge = 1
                exit
            end if
        end do 
        if(notconverge .ne. 0) then
            Td = T_ini
            call find_dop(Td, kappa)
            kapb_d = kappa*cgs_asb*cgs_c*Td**4.d0
            f = col_g * (Tg - Td) + rad_heat - kapb_d
            fu = f
            fl = f
            Td_l = Td
            Td_u = Td
            iteration = 0
            dTul = 1.d-1
            if(f >= 0.d0) then
                do while (fl >= 0.d0)
                    iteration = iteration + 1
                    Td_l = Td * (1.d0 + dTul)
                    call find_dop(Td_l, kappa)
                    kapb_d = kappa*cgs_asb*cgs_c*Td_l**4.d0
                    fl = col_g * (Tg - Td_l) + rad_heat - kapb_d
                    dTul = dTul * 5.d0
                    if(iteration .eq. 30) then
                        Td_l = 1.d5
                        exit
                    endif
                end do
            else
                do while (fu < 0.d0)
                    iteration = iteration + 1
                    Td_u = Td / (1.d0 + dTul)
                    call find_dop(Td_u, kappa)
                    kapb_d = kappa*cgs_asb*cgs_c*Td_u**4.d0
                    fu = col_g * (Tg - Td_u) + rad_heat - kapb_d
                    dTul = dTul * 5.d0
                    if(iteration .eq. 30) then
                        Td_u = 0.1d0
                        exit
                    endif
                end do
            end if
            derror = 1.d0
            Td_m = 0.5d0 * (Td_u + Td_l)
            iteration = 0
            dT_lm = 1.d0
            do while(dT_lm > 1.d-5)
                iteration = iteration + 1
                call find_dop(Td_m, kappa)
                kapb_d = kappa*cgs_asb*cgs_c*Td_m**4.d0
                fm = col_g * (Tg - Td_m) + rad_heat - kapb_d
                if(fm >= 0.d0) then
                    Td_u = Td_m
                else
                    Td_l = Td_m
                endif
                Td_m = 0.5d0 * (Td_u + Td_l)
                Td = Td_m
                max_v = max(rad_heat, kapb_d)
                derror = abs(fm) / max_v
                dT_lm = abs((Td_u - Td_l)/Td_m)
                if(iteration == 30) then
                    notconverge = 101
                    Td = Td_m
                    exit
                end if
            end do
        end if 
        if(Td < 0.0)then
            print *, 'Td < 0 K appear and rescue Td'
            print *, 'duste_temp', Td, notconverge, col_g * (Tg - Td) ,rad_heat , kapb_d
            Td = 1.d0
        endif
        Td = max(Td,MP_Tcmb)
        call find_dop(Td, kappa)
        chi_d = 4.0*kappa*cgs_asb*cgs_c*Td**3.d0/col_g
    end subroutine dtemp_radtr
    subroutine find_dop(Td, xkd)
        implicit none
        real(kind=8),intent(IN) :: Td
        real(kind=8),intent(OUT):: xkd
        integer ii
        if(Td >= 1.d5) then
            xkd = data_dp(99999)
            return
        elseif(Td <= 0.1d0) then
            ii = 0
        else
            ii = int(1.d5 * (dlog10(Td)+1.d0 )/ 6.d0)
            ii = min(99998, ii)
        endif
        xkd = data_dp(ii) + (data_dp(ii+1)-data_dp(ii))/(data_temp(ii+1)-data_temp(ii))*(Td-data_temp(ii))
    end subroutine find_dop
    subroutine find_dros(Td, xkd)
        implicit none
        real(kind=8),intent(IN) :: Td
        real(kind=8),intent(OUT):: xkd
        integer ii
        if(Td >= 1.d5) then
            xkd = data_ros(99999)
            return
        elseif(Td <= 0.1d0) then
            ii = 0
        else
            ii = int(1.d5 * (dlog10(Td)+1.d0 )/ 6.d0)
            ii = min(99998, ii)
        endif
        xkd = data_ros(ii) + (data_ros(ii+1)-data_ros(ii))/(data_temp(ii+1)-data_temp(ii))*(Td-data_temp(ii))
    end subroutine find_dros
    subroutine read_plankop
        implicit none
        integer,parameter :: FH = 13
        integer,parameter :: FH2 = 15
        character(100) :: path2dust = "./Dust_op/plank_op.dat"
        character(100) :: path2dust2 = "./Dust_op/ross_op.dat"
        character(len=100) :: ffn
        integer, parameter :: inum = 100000
        integer :: i, err, ierr1, ierr2
        real(kind=8) :: xx, xop
        ffn= trim(path2dust)
        open(FH, file=ffn,status='old',iostat=err)
        if (err > 0) then
            if(get_myrank() == 0) print '(A,/,A)', "dustop: file not found","stopping..."
            stop
        end if
        do i = 0, inum-1
            read(FH,*,iostat=err) xx, xop
            if (err > 0) then
                print *, "dustop: **WARNING** data size might be inconsistent"
                exit
            endif
            data_temp(i) = xx
            data_dp(i) = xop*fdust_solar
        enddo
        close(FH)
        ffn= trim(path2dust2)
        open(FH2, file=ffn,status='old',iostat=err)
        if (err > 0) then
            if(get_myrank() == 0) print '(A,/,A)', "dustop: file not found","stopping..."
            stop
        end if
        do i = 0, inum-1
            read(FH2,*,iostat=err) xx, xop
            if (err > 0) then
                print *, "dustop: **WARNING** data size might be inconsistent"
                exit
            endif
            data_ros(i) = xop*fdust_solar
        enddo
        close(FH2)
    end subroutine read_plankop
    subroutine PhotoelectricHeating(ne0, nH, Tg, Gfuv, gamma_pe, fd)
        implicit none
        real(kind=8),intent(IN) :: ne0 
        real(kind=8),intent(IN) :: nH 
        real(kind=8),intent(IN) :: Tg 
        real(kind=8),intent(IN) :: Gfuv, fd 
        real(kind=8),intent(OUT):: gamma_pe 
        real(kind=8) :: ne, epsi, ggtg, heat, cool, beta, Tg_fit
        ne = ne0 
        Tg_fit = Tg
        if(ne > 0.d0) then
            ggtg = Gfuv*dsqrt(Tg_fit)/ne
            epsi = 4.87d-2/(1.0 + 4.d-3*(ggtg)**0.73) + 3.65d-2*(Tg_fit*1.d-4)**0.7/(1.0 + 2.d-4*ggtg)
            heat = 1.d-24 * epsi * Gfuv * nH * fd
            beta = 0.735d0/Tg_fit**0.068d0
            cool = 3.49d-30*Tg_fit**0.944d0*ggtg**beta*ne*nH *fd
            gamma_pe = heat - cool
        else
            gamma_pe = 0.d0
        endif
    end subroutine PhotoelectricHeating
    subroutine depsilon(Td, epsid)
        implicit none
        real(kind=8),intent(IN) :: Td 
        real(kind=8),intent(OUT) :: epsid 
        real(kind=8), parameter :: Tevap = 1200.e0
        real(kind=8) :: xx
        xx = (Td - Tevap)/50.d0
        epsid = 1.d0
    end subroutine depsilon
    subroutine OII_OII_ratio(ne, y_HII, T_gas, kpiOII, y_OII, y_OIII, mmetal)
        implicit none
        real(kind=8),intent(IN) :: ne, y_HII, T_gas, kpiOII, mmetal
        real(kind=8),intent(OUT) :: y_OII, y_OIII
        real(kind=8) :: y_Oion, y_Oion2, y_Oion3, alpha_RR, alpha_DR, tRR, tDR, yy_OII,yy_OIII
        real(kind=8) :: tDR2, x, Nealpha, xi_col, U_col,beta
        real(kind=8), save :: Z_OIII, a_RR, b_RR, c_RR, d_RR, a_DR, b_DR, c_DR, d_DR, f_DR &
          , A_col, P_col, X_col, K_col, delta_E
        integer,save :: ifirst = 0
        if(ifirst == 0) then
            Z_OIII = 2.d0
            a_RR = 4.092
            b_RR = -0.6413d0
            c_RR = 0.d0
            d_RR = 1.0d0
            a_DR = -0.0036d0
            b_DR = 0.7519d0
            c_DR = 1.5252d0
            d_DR = -0.0838d0
            f_DR = 0.2769d0
            A_col = 0.139d-7
            delta_E = 35.1e0 * cgs_ev
            P_col = 1.d0
            X_col = 0.212
            K_col = 0.22
            ifirst = 1
        endif
        y_Oion = MP_frac_O_solar*mmetal*y_HII
        tRR = T_gas * 1.d-4 / (Z_OIII*Z_OIII)
        a_RR = 1.d-13 * Z_OIII * a_RR * tRR**b_RR / (1.d0 + c_RR * tRR**d_RR)
        tDR = 1.d-4 * T_gas
        tDR2 = tDR*tDR
        a_DR = 1.d-12 * (a_DR/tDR + b_DR + c_DR*tDR + d_DR*tDR2) * tDR**(-1.5d0) * dexp(-f_DR/tDR)
        U_col = delta_E / (T_gas*cgs_kb)
        xi_col = A_col * (1.d0 + P_col*dsqrt(U_col))/(X_col+U_col) * U_col**K_col*dexp(-U_col)
        Nealpha = ne * (a_RR + a_DR) 
        yy_OII = Nealpha / (kpiOII + ne* xi_col + Nealpha) * y_Oion
        yy_OIII = y_Oion - yy_OII
        yy_OII = max(yy_OII ,0.e0)
        yy_OII = min(yy_OII ,y_Oion)
        yy_OIII = max(yy_OIII,0.e0)
        yy_OIII = min(yy_OIII,y_Oion)
        y_OII = yy_OII
        y_OIII = yy_OIII
    end subroutine OII_OII_ratio
    subroutine update_yco(xnH, T_K, y, yco, rcopd, mmetal, dt)
        implicit none
        real(kind=8),intent(IN) :: xnH, T_K, y(0:6 -1), rcopd, mmetal, dt
        real(kind=8),intent(INOUT) :: yco
        real(kind=8) :: frac_O, frac_C
        real(kind=8),parameter :: k0 = 5.d-16 
        real(kind=8),parameter :: k1 = 5.d-10 
        real(kind=8),parameter :: min_value = 1.d-50
        real(kind=8) :: nH2, nOI, gam_chx, rf, knOI, rd, yo, yco_o
        frac_O = MP_frac_O_solar*mmetal
        frac_C = MP_frac_C_solar*mmetal
        yco_o = yco
        nH2 = y(1)*xnH 
        yo = (1.d0 - y(3)) *frac_O - yco_o
        yo = MAX(MIN(frac_O,yo),min_value)
        nOI = yo*xnH 
        rd = rcopd 
        gam_chx = 5.d0 * rd 
        knOI = k1 * nOI
        rf = k0 * nH2 * knOI / (knOI + gam_chx) 
        yco = ( yco_o + rf*dt*frac_C ) / (1.d0 + rf*dt + rd*dt)
        yco = MAX(MIN(frac_C,yco),min_value)
    end subroutine update_yco
    subroutine COcool(xnH,T_K,y_H2,xNc_CO,tau_cnt,xLd_CO)
      implicit none
      real(kind=8),intent(IN) :: xnH,T_K,y_H2,tau_cnt
      real(kind=8),intent(OUT):: xLd_CO
      real(kind=8) :: xNc_CO
      real(kind=8) :: xlTa(1:11)=(/1.000d0, 1.301d0, 1.477d0 ,1.699d0, 1.903d0 &
        , 2.000d0, 2.477d0, 2.778d0,3.000d0, 3.176d0, 3.301d0/)
      real(kind=8) :: xlNa(1:11)=(/14.0d0, 14.5d0, 15.0d0, 15.5d0, 16.0d0 &
        ,16.5d0, 17.0d0, 17.5d0, 18.0d0, 18.5d0, 19.0d0/)
      real(kind=8) :: aL0a(1:11)=(/0.2477d2, 0.2438d2, 0.2421d2, 0.2403d2, 0.2389d2, 0.2382d2 &
        ,0.2342d2, 0.2313d2, 0.2291d2, 0.2263d2, 0.2228d2/)
      real(kind=8) :: aLLTEa(1:11,1:11) = reshape((/0.2108d2, 0.2035d2 &
        , 0.1994d2, 0.1945d2, 0.1901d2, 0.1880d2,&
          0.1781d2, 0.1723d2, 0.1686d2, 0.1666d2, 0.1655d2, &
          0.2109d2, 0.2035d2, 0.1995d2, 0.1945d2, 0.1901d2, 0.1880d2, &
          0.1781d2, 0.1723d2, 0.1686d2, 0.1666d2, 0.1655d2, &
          0.2111d2, 0.2037d2, 0.1996d2, 0.1946d2, 0.1901d2, 0.1880d2, &
          0.1781d2, 0.1723d2, 0.1686d2, 0.1666d2, 0.1655d2, &
          0.2118d2, 0.2040d2, 0.1998d2, 0.1947d2, 0.1902d2, 0.1881d2, &
          0.1782d2, 0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, &
          0.2137d2, 0.2051d2, 0.2005d2, 0.1952d2, 0.1905d2, 0.1883d2, &
          0.1782d2, 0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, &
          0.2167d2, 0.2073d2, 0.2023d2, 0.1964d2, 0.1913d2, 0.1890d2, &
          0.1785d2, 0.1725d2, 0.1688d2, 0.1667d2, 0.1656d2, &
          0.2204d2, 0.2105d2, 0.2052d2, 0.1987d2, 0.1932d2, 0.1906d2, &
          0.1792d2, 0.1728d2, 0.1690d2, 0.1669d2, 0.1658d2, &
          0.2244d2, 0.2142d2, 0.2086d2, 0.2019d2, 0.1960d2, 0.1933d2, &
          0.1808d2, 0.1738d2, 0.1697d2, 0.1675d2, 0.1663d2, &
          0.2287d2, 0.2182d2, 0.2124d2, 0.2055d2, 0.1995d2, 0.1966d2, &
          0.1834d2, 0.1759d2, 0.1715d2, 0.1691d2, 0.1678d2, &
          0.2330d2, 0.2223d2, 0.2165d2, 0.2094d2, 0.2032d2, 0.2003d2, &
          0.1867d2, 0.1789d2, 0.1748d2, 0.1726d2, 0.1712d2, &
          0.2376d2, 0.2266d2, 0.2206d2, 0.2135d2, 0.2071d2, 0.2042d2, &
          0.1903d2, 0.1826d2, 0.1793d2, 0.1774d2, 0.1761d2/),(/11,11/))
       real(kind=8) :: xlnha(1:11,1:11) =reshape((/0.329d1, &
          0.349d1, 0.367d1, 0.397d1, 0.430d1, 0.446d1, &
          0.517d1, 0.547d1, 0.553d1, 0.530d1, 0.470d1, &
          0.327d1, 0.348d1, 0.366d1, 0.396d1, 0.430d1, 0.445d1, &
          0.516d1, 0.547d1, 0.553d1, 0.530d1, 0.470d1, &
          0.322d1, 0.345d1, 0.364d1, 0.394d1, 0.429d1, 0.445d1, &
          0.516d1, 0.547d1, 0.553d1, 0.530d1, 0.470d1, &
          0.307d1, 0.334d1, 0.356d1, 0.389d1, 0.426d1, 0.442d1, &
          0.515d1, 0.546d1, 0.552d1, 0.530d1, 0.470d1, &
          0.272d1, 0.309d1, 0.335d1, 0.374d1, 0.416d1, 0.434d1, &
          0.513d1, 0.545d1, 0.551d1, 0.529d1, 0.468d1, &
          0.224d1, 0.265d1, 0.295d1, 0.342d1, 0.392d1, 0.414d1, &
          0.506d1, 0.541d1, 0.548d1, 0.526d1, 0.464d1, &
          0.174d1, 0.215d1, 0.247d1, 0.295d1, 0.349d1, 0.374d1, &
          0.486d1, 0.530d1, 0.539d1, 0.517d1, 0.453d1, &
          0.124d1, 0.165d1, 0.197d1, 0.245d1, 0.300d1, 0.325d1, &
          0.447d1, 0.502d1, 0.516d1, 0.494d1, 0.427d1, &
          0.742d0, 0.115d1, 0.147d1, 0.195d1, 0.250d1, 0.275d1, &
          0.398d1, 0.457d1, 0.473d1, 0.452d1, 0.384d1, &
          0.242d0, 0.652d0, 0.966d0, 0.145d1, 0.200d1, 0.225d1, &
          0.348d1, 0.407d1, 0.424d1, 0.403d1, 0.335d1, &
         -0.258d0, 0.152d0, 0.466d0, 0.954d0, 0.150d1, 0.175d1, &
          0.298d1, 0.357d1, 0.374d1, 0.353d1, 0.285d1/),(/11,11/))
      real(kind=8):: alphaa(1:11,1:11) =reshape((/0.439d0, &
          0.409d0, 0.392d0, 0.370d0, 0.361d0, 0.357d0, &
          0.385d0, 0.437d0, 0.428d0, 0.354d0, 0.322d0, &
          0.436d0, 0.407d0, 0.391d0, 0.368d0, 0.359d0, 0.356d0, &
          0.385d0, 0.437d0, 0.427d0, 0.354d0, 0.322d0, &
          0.428d0, 0.401d0, 0.385d0, 0.364d0, 0.356d0, 0.352d0, &
          0.383d0, 0.436d0, 0.427d0, 0.352d0, 0.320d0, &
          0.416d0, 0.388d0, 0.373d0, 0.353d0, 0.347d0, 0.345d0, &
          0.380d0, 0.434d0, 0.425d0, 0.349d0, 0.316d0, &
          0.416d0, 0.378d0, 0.360d0, 0.338d0, 0.332d0, 0.330d0, &
          0.371d0, 0.429d0, 0.421d0, 0.341d0, 0.307d0, &
          0.450d0, 0.396d0, 0.367d0, 0.334d0, 0.322d0, 0.317d0, &
          0.355d0, 0.419d0, 0.414d0, 0.329d0, 0.292d0, &
          0.492d0, 0.435d0, 0.403d0, 0.362d0, 0.339d0, 0.329d0, &
          0.343d0, 0.406d0, 0.401d0, 0.317d0, 0.276d0, &
          0.529d0, 0.473d0, 0.441d0, 0.404d0, 0.381d0, 0.370d0, &
          0.362d0, 0.410d0, 0.392d0, 0.316d0, 0.272d0, &
          0.555d0, 0.503d0, 0.473d0, 0.440d0, 0.423d0, 0.414d0, &
          0.418d0, 0.446d0, 0.404d0, 0.335d0, 0.289d0, &
          0.582d0, 0.528d0, 0.499d0, 0.469d0, 0.457d0, 0.451d0, &
          0.470d0, 0.487d0, 0.432d0, 0.364d0, 0.310d0, &
          0.596d0, 0.546d0, 0.519d0, 0.492d0, 0.483d0, 0.479d0, &
          0.510d0, 0.516d0, 0.448d0, 0.372d0, 0.313d0/),(/11,11/))
      real(kind=8) :: xn_c, xlT, cs, xlN, aL0, aLLTE, xlnh, alpha
      real(kind=8) :: xL0inv,xLLTEinv,xn_h, xLinv,xL
      xn_c=y_H2*xnH
      xlT=dlog10(T_K)
      cs=1.d-5*dsqrt(2.d0*cgs_kb*T_K/(28.d0*cgs_mp))
      xNc_CO=xNc_CO+1.d-10
      xlN=dlog10(xNc_CO/cs)
      call linear2(xlTa,aL0a,11,xlT,aL0)
      call bilinear2(xlTa,xlNa,aLLTEa,11,11,xlT,xlN,aLLTE)
      call bilinear2(xlTa,xlNa,xlnha,11,11,xlT,xlN,xlnh)
      call bilinear2(xlTa,xlNa,alphaa,11,11,xlT,xlN,alpha)
      xL0inv=10.d0**aL0
      xLLTEinv=10.d0**aLLTE
      xn_h=10.d0**xlnh
      xLinv=xL0inv+xn_c*xLLTEinv &
        +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
      xL=1.d0/xLinv
      xLd_CO=xL*dexp(-tau_cnt)
contains
      subroutine linear2(xa,ya,m,x,y)
        implicit none
        integer :: m, ms,i
        real(kind=8), dimension(m) :: xa, ya
        real(kind=8) :: x, y, y1, y2, t
        do 11 i=1,m
            if(x-xa(i).le.0.d0) then
                ms=i
                go to 12
            endif
 11 continue
        ms=m
 12 continue
        if(ms.eq.1) ms=2
        y1=ya(ms-1)
        y2=ya(ms)
        t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
        y=(1.d0-t)*y1+t*y2
        return
      end subroutine linear2
      subroutine bilinear2(x1a,x2a,ya,m,n,x1,x2,y)
        implicit none
        integer :: m, n, i, ms, ns
        real(kind=8) :: x1, x2, y
        real(kind=8), dimension(m) :: x1a
        real(kind=8), dimension(n) :: x2a
        real(kind=8), dimension(m,n) :: ya
        real(kind=8) :: y1, y2, y3, y4, t, u
        do 11 i=1,m
            if(x1-x1a(i).le.0.d0) then
                ms=i
                go to 12
            endif
 11 continue
        ms=m
 12 continue
        do 13 i=1,n
            if(x2-x2a(i).le.0.d0) then
                ns=i
                go to 14
            endif
 13 continue
        ns=n
 14 continue
        if(ms.eq.1) ms=2
        if(ns.eq.1) ns=2
        y1=ya(ms-1,ns-1)
        y2=ya(ms,ns-1)
        y3=ya(ms,ns)
        y4=ya(ms-1,ns)
        t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
        u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
        y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
      return
      end subroutine bilinear2
      end subroutine COcool
      subroutine COdissociation_rate(NcCO, NcH2, Av, Gfuv2, rcod)
        real(kind=8),intent(IN) :: NcCO, NcH2, Av, Gfuv2
        real(kind=8),intent(OUT):: rcod
        real(kind=8), parameter :: pdis = 1.03d-10 
        real(kind=8) :: theta1, theta2, theta3, logNco, logNH2, logtheta, logAv
        if(NcCO > 1.d12) then
            logNco = dlog10(NcCO)
            logtheta = (((((-1.62507133d-4*logNco+1.43588538d-2)*logNco &
                        -5.21418820d-1)*logNco+9.95253514d0)*logNco-1.05308265d2)*logNco &
                        +5.85860800d2)*logNco-1.33950326d3
            theta1 = 10.d0**logtheta
        else
            theta1 = 1.d0
        endif
        if(NcH2 > 1.e13) then
            logNH2 = dlog10(NcH2)
            logtheta = -2.09263955d-18*dexp(1.89035939*logNH2)
            theta2 = 10.d0**logtheta
        else
            theta2 = 1.d0
        endif
        if(Av < 40.d0) then
            theta3 = dexp(-2.5*Av)
        else
            theta3 = 0.d0
        endif
        rcod = Gfuv2 * pdis * theta1 * theta2 * theta3 
      end subroutine COdissociation_rate
end module kinzoku
