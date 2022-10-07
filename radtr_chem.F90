#include "config.h"
#include "chemistry_label.h"
#define DEBUG_RADTR   NO

!--------------------------------------------------------------------------------
! 
!  Module for chemistry with radtr
!
!--------------------------------------------------------------------------------
module radtr_chem
  use unit
  use parameter, only : Pi
  use modelParameter
  use grid
  use mpilib
  use sinkParticle, only : sp_getSinkRadius, sp_getLevel
#if MODEL_ART > 0
  use radiationSource
  use primordial,only: CoolSolverExplicit, CoolSolverImplicit, adjust_abundance &
    , yHe, c_H2, dbg_flg_prm, chem_vals, chem_vals, get_xmu
#endif !MODEL_ART    


#ifdef METAL_TRANSFER
  use kinzoku, only : dtemp_radtr, find_dop, find_dros, fdust_solar
#else
  use kinzoku, only : dtemp_radtr, find_dop, find_dros, fdust_solar
#endif

  implicit none
  private

  ! parameter for photoionization
  real(kind=DBL_KIND) :: SinkRadius
  type(t_rs_info),save,pointer :: rs_info    !info for radiation source

#ifdef M1CLOSER_EUV_TRANSFER
  type rad_mean
    real(kind=DBL_KIND) :: alpha_EUV  ! crossection of EUV for HI [cm^2]
    real(kind=DBL_KIND) :: heat_EUV   ! heating rate of EUV per one ionizing event [erg]
    real(kind=DBL_KIND) :: sig_EUV    ! dust cross section for EUV photons [cm^2]
    real(kind=DBL_KIND) :: sig_FUV    ! dust cross section for FUV photons [cm^2]
    real(kind=DBL_KIND) :: rOII       ! ionization rate of OII to HI
    real(kind=DBL_KIND) :: erg_EUV    ! total energy of per one EUV photons [erg]
    real(kind=DBL_KIND) :: erg_FUV    ! total energy of per one FUV photons [erg]
    integer :: num_rad                ! number of radiation source
  end type rad_mean
#endif


#ifdef M1CLOSER_EUV_TRANSFER
  public :: ch_radtrchem, rad_mean
#endif
#ifdef M1CLOSER_IR_TRANSFER
  public :: radtr_IRdust
#endif

contains

  !-----------------------------------------
  !   chemistry  更新
  !----------------------------------------

#ifdef M1CLOSER_EUV_TRANSFER
  subroutine ch_radtrchem(dt_code, dt_hyd, rsm)

    use modelParameter,only : MP_Tmin,MP_Tmax
    use chemistry, only : ch_CellChemCool

    real(kind=DBL_KIND),intent(IN) :: dt_code,dt_hyd
    type(rad_mean), intent(IN) :: rsm
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    !----------EO_added----------!
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: vx, vy, vz 
    !----------EO_added----------!
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p
#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: yco
#endif

    real(kind=DBL_KIND) :: l_jeans, l_sobolev, dop_b5, cf_FUV  
    real(kind=DBL_KIND),dimension(MX:MZ) :: dvdr
    real(kind=DBL_KIND), dimension(6)    :: dvdr_s, NH2_sbv, lsb_s, NCO_sbv, lsb_rho_s

  #ifdef M1CLOSER_EUV_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad_FUV, Frx_FUV, Fry_FUV, Frz_FUV
    real(kind=DBL_KIND) :: lfuv, fsH2_1, fsH2_2, fsCO_1, fsCO_2, Gfuv_0
    real(kind=DBL_KIND),parameter :: sig_CO = 1.594880d-18 ![cm^2]
    real(kind=DBL_KIND),parameter :: sigd_FUV = 2.23212d-21

    ! Drain & Bertldi
    real(kind=DBL_KIND),parameter :: Efuv_0  = 5.29d-14    ! erg em^{-3}
    real(kind=DBL_KIND),parameter :: xi_d0   = 4.27d-11  ! s^{-1} photodissociation rate per Efuv_0

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad_FUV_CO, Frx_FUV_CO, Fry_FUV_CO, Frz_FUV_CO
      real(kind=DBL_KIND) :: Nph_fuv_co
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad_FUV_DUST, Frx_FUV_DUST, Fry_FUV_DUST, Frz_FUV_DUST
      real(kind=DBL_KIND) :: Nph_fuv_dust, cf_FUV
    #endif


  #endif




  #ifdef M1CLOSER_IR_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad_IR, Frx_IR, Fry_IR, Frz_IR 
  #endif
  #ifdef METAL
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Tdust
  #endif
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: mu, khmpd
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem0
    real(kind=DBL_KIND) :: xmu,  cs, T_K0, yco0, sig_euv, sig_fuv
    real(kind=DBL_KIND) :: Nph_euv, Nph_fuv, Eph_fuv, dtau, dtau_exp, rho_cgs, xNcH2_1, xNcCO_1, xNcH2_2, xNcCO_2
    real(kind=DBL_KIND) :: xNcH_1, xNcH_2, v_th, dvdr_av, xNcH_H2, xNcH_CO
    integer :: level, n, gid, i, j, k, ic
  #ifdef METAL
    #ifdef DUST_NOTCONSTANT
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
    #endif
  #endif
    real(kind=DBL_KIND) :: Trad, rho_dust, xk_dust, xk_rad, xk_ros, tap1, tap2, d_EradIR, d_NradIR

  #ifndef NO_RADFORCE
    #ifdef EXTERNALFORCE
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    #endif
  #endif !NO_RADFORCE
  #ifdef EXTERNALFORCE
    real(kind=DBL_KIND) :: cf_EUV
  #endif
  #ifdef DUST_NOTCONSTANT
    real(kind=DBL_KIND) :: Tevap, t_dyn  
  #endif
  #ifdef METAL_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
  #endif
    integer :: sp_Level

    logical :: isNotFinite 
    integer :: count_output=0, max_count=100 !KS DEBUG
    logical :: naninf_flg
    logical :: sourceCell_flag
    type(chem_vals) :: xch

    real(kind=DBL_KIND), dimension(Imin-1:Imax+1,Jmin-1:Jmax+1,Kmin-1:Kmax+1) :: yh2_old, yco_old

    ! ----------------------------------------------------------
    type chemsp 
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(NCEHM_MIN:NCEHM_MAX)
    integer :: ichem
    ! ----------------------------------------------------------


    !------------------------------
    ! get information of RS
    !------------------------------
    call rs_GetSourceInfo(rs_info)
    SinkRadius = sp_getSinkRadius()

    ! -----------------------------
    ! get splevel
    ! -----------------------------
    sp_Level = sp_getLevel()

    ! timestep in cgs
    xch%dt = dt_code*Unit_t

    naninf_flg = .False.

    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        
        ! 子グリッドを持つ場合はchemistryは解かなくてよい (あとでconvergeを呼ぶ)
        if (has_child_grid(gid))   cycle

        !座標
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)

        vx => get_Ucomp(MVX, gid)
        vy => get_Ucomp(MVY, gid)
        vz => get_Ucomp(MVZ, gid)

                  !変数の取り出し
        rho => get_Ucomp(MRHO,gid)
        p => get_Ucomp(MP,gid)


        ! get chemistry pointer
        do ichem = NCEHM_MIN, NCEHM_MAX
          chemary3(ichem)%y => get_Ucomp(ichem,gid)
        enddo

  #ifdef CHEM_MODEL_HF2020
        yco  => get_Ucomp(MCO, gid)
  #endif


  #ifdef METAL
        Tdust => get_Ucomp(MTD,gid)
  #endif

  #ifndef NO_RADFORCE
    #ifdef EXTERNALFORCE
        fx_hpi => get_Ucomp(MXPI,gid)
        fy_hpi => get_Ucomp(MYPI,gid)
        fz_hpi => get_Ucomp(MZPI,gid)
    #endif
  #endif !NO_RADFORCE


  #ifdef M1CLOSER_EUV_TRANSFER
        Nrad => get_Ucomp(MER, gid)
        Frx  => get_Ucomp(MFRX, gid)
        Fry  => get_Ucomp(MFRY, gid)
        Frz  => get_Ucomp(MFRZ, gid)

  #endif

  #ifdef M1CLOSER_FUV_TRANSFER
        Nrad_FUV => get_Ucomp(MEF, gid)
        Frx_FUV  => get_Ucomp(MFRFX, gid)
        Fry_FUV  => get_Ucomp(MFRFY, gid)
        Frz_FUV  => get_Ucomp(MFRFZ, gid)

    #ifdef  M1CLOSER_SEPARATE_FUV_TRANS 
        Nrad_FUV_CO => get_Ucomp(MECOF, gid)
        Frx_FUV_CO  => get_Ucomp(MFCORFX, gid)
        Fry_FUV_CO  => get_Ucomp(MFCORFY, gid)
        Frz_FUV_CO  => get_Ucomp(MFCORFZ, gid)

        Nrad_FUV_DUST => get_Ucomp(MEDUSTF, gid)
        Frx_FUV_DUST  => get_Ucomp(MFDUSTRFX, gid)
        Fry_FUV_DUST  => get_Ucomp(MFDUSTRFY, gid)
        Frz_FUV_DUST  => get_Ucomp(MFDUSTRFZ, gid)
    #endif
  #endif

  #ifdef M1CLOSER_IR_TRANSFER
        Nrad_IR => get_Ucomp(MEIR, gid)
        Frx_IR  => get_Ucomp(MFIRX, gid)
        Fry_IR  => get_Ucomp(MFIRY, gid)
        Frz_IR  => get_Ucomp(MFIRZ, gid)
  #endif

  #ifdef defined(METAL) && defined(DUST_NOTCONSTANT)
        rhod => get_Ucomp(MDRHO,gid)
  #endif

  #ifdef METAL_TRANSFER
         mmetal => get_Ucomp(MMET,gid)
  #endif

        call GetKHmpdThin(x,y,z,khmpd) !get khmpd


        ! -------------------------------
        !        save old values 
        ! -------------------------------
        do k=Kmin-1, Kmax+1
          do j=Jmin-1, Jmax+1
            do  i=Imin-1, Imax+1
              yh2_old(i,j,k) = chemary3(MH2)%y(i,j,k)
  #ifdef CHEM_MODEL_HF2020
              yco_old(i,j,k) = yco(i,j,k)
  #else
              yco_old(i,j,k) = chemary3(MCO)%y(i,j,k)
  #endif
            enddo
          enddo
        enddo

        do i=Imin, Imax
          do j=Jmin, Jmax
            do k=Kmin, Kmax


        !----------EO_added----------!
       ! if(Nrad(i,j,k) /= 0.d0) then
       !   print*, "*********Nrad is not zero begin ch_radtrchem********"
       ! end if

        !----------EO_added----------!


              ! -----------------------
              !    set metallicity
              ! -----------------------
  #ifdef METAL_TRANSFER
              xch%metal = mmetal(i,j,k)
  #else
              xch%metal = MP_Metallicity
  #endif

  #ifdef METAL
              xch%Td = Tdust(i,j,k)
  #endif
              

              ! --------------------------
              !   化学組成変数の詰め替え
              ! --------------------------
              do ichem = 0, NCHEM-1
                 xch%ychem(ichem) = chemary3(ichem+NCEHM_MIN)%y(i,j,k)
              enddo 
  #ifdef CHEM_MODEL_HF2020
              xch%yco = yco(i,j,k)
  #endif


              !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
              call adjust_abundance(xch%ychem &
#ifdef CHEM_MODEL_HF2020
                ,xch%yco &
#endif
                ,xch%metal)

              !derived variables
              rho_cgs= rho(i,j,k)*Unit_rho
              xch%nH = rho_cgs*Unit_invmu                !水素原子核の数密度
              xmu    = get_xmu(xch%ychem) ! 平均分子量
              xch%Tg = p(i,j,k)*Unit_e*cgs_amu*xmu /rho_cgs/cgs_kb      !温度 [K]


              !derived variables
              if (xch%Tg < MP_Tmin*0.999d0 .or. xch%Tg > MP_Tmax*1.001d0) then !表示を抑制するため係数をかける
              
                !output warning
                if (count_output < max_count) then
                   print '(A, 6(1PE12.4))', &
                        "*** before chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                        xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                   count_output = count_output+1
                   !suppress warning
                   if (count_output == max_count) then
                      print *, 'count_output reaches max_count: no more warning will be shown'
                      count_output = count_output+1
                   end if
                end if
              end if
              xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)

              ! ------------------------------
              !        check sink cell
              ! ------------------------------
              if (level .eq. sp_Level) then
                call check_is_sourceCell(x(i), y(j), z(k), sourceCell_flag)
              else
                sourceCell_flag = .false.
              endif

              !---------------------- dust grain abundance ----------------------!
  #if defined(METAL) && defined(DUST_NOTCONSTANT)
              rhod(i,j,k)= min( rhod(i,j,k), fdust_solar*MP_Metallicity)
              xch%fd     = rhod(i,j,k) / fdust_solar ! dust abundance including sublimation & metallicity    
              xch%fd     = max(xch%fd, 0.d0)
  #elif defined(METAL) && defined(METAL_TRANSFER)
              xch%fd     = mmetal(i,j,k)
  #elif defined(METAL)
              xch%fd     = MP_Metallicity
  #else
              xch%fd     = 0.d0
  #endif
              !------------------------------------------------------------------!


              ! ------------------------------------------------
              !       H2 & CO photodisocciation rates
              ! ------------------------------------------------
              cs      = sqrt(cgs_kb*xch%Tg/(xmu*cgs_amu))       !等温音速 [cm/s]
              l_jeans = cs*sqrt(Pi/(cgs_gc*rho_cgs))            !ジーンズ長さ [cm]

              ! sobolev length each direction 
              dvdr_s(1)    = dabs(vx(i+1,j,k)-vx(i,j,k))/(CellWidth(MX,level)*Unit_t)+1.d-30 ! [s^-1]
              dvdr_s(2)    = dabs(vx(i,j,k)-vx(i-1,j,k))/(CellWidth(MX,level)*Unit_t)+1.d-30 ! [s^-1]
              dvdr_s(3)    = dabs(vy(i,j+1,k)-vy(i,j,k))/(CellWidth(MY,level)*Unit_t)+1.d-30 ! [s^-1]
              dvdr_s(4)    = dabs(vy(i,j,k)-vy(i,j-1,k))/(CellWidth(MY,level)*Unit_t)+1.d-30 ! [s^-1]
              dvdr_s(5)    = dabs(vz(i,j,k+1)-vz(i,j,k))/(CellWidth(MZ,level)*Unit_t)+1.d-30 ! [s^-1]
              dvdr_s(6)    = dabs(vz(i,j,k)-vz(i,j,k-1))/(CellWidth(MZ,level)*Unit_t)+1.d-30 ! [s^-1]

              ! sobolev length for H2 
              v_th         = sqrt(cgs_kb*xch%Tg/cgs_mh) 
              lsb_s(:)     = v_th/dvdr_s(:)

              NH2_sbv(1)   = min(yh2_old(i+1,j,k), yh2_old(i,j,k))*lsb_s(1)
              NH2_sbv(2)   = min(yh2_old(i,j,k), yh2_old(i-1,j,k))*lsb_s(2)
              NH2_sbv(3)   = min(yh2_old(i,j+1,k), yh2_old(i,j,k))*lsb_s(3)
              NH2_sbv(4)   = min(yh2_old(i,j,k), yh2_old(i,j-1,k))*lsb_s(4)
              NH2_sbv(5)   = min(yh2_old(i,j,k+1), yh2_old(i,j,k))*lsb_s(5)
              NH2_sbv(6)   = min(yh2_old(i,j,k), yh2_old(i,j,k-1))*lsb_s(6)

              ! sobolev length for CO
              v_th         = sqrt(2.d0*cgs_kb*xch%Tg/((MP_AC+MP_AO)*cgs_amu))
              lsb_s(:)     = v_th/dvdr_s(:)

              NCO_sbv(1)   = min(yco_old(i+1,j,k), yco_old(i,j,k))*lsb_s(1)
              NCO_sbv(2)   = min(yco_old(i,j,k), yco_old(i-1,j,k))*lsb_s(2)
              NCO_sbv(3)   = min(yco_old(i,j+1,k), yco_old(i,j,k))*lsb_s(3)
              NCO_sbv(4)   = min(yco_old(i,j,k), yco_old(i,j-1,k))*lsb_s(4)
              NCO_sbv(5)   = min(yco_old(i,j,k+1), yco_old(i,j,k))*lsb_s(5)
              NCO_sbv(6)   = min(yco_old(i,j,k), yco_old(i,j,k-1))*lsb_s(6)

              ! sobolev length
              xch%dvdr(MX) = (dvdr_s(1)+dvdr_s(2))*0.5d0
              xch%dvdr(MY) = (dvdr_s(3)+dvdr_s(4))*0.5d0  
              xch%dvdr(MZ) = (dvdr_s(5)+dvdr_s(6))*0.5d0 

              xch%xlmbdj = l_jeans                ! jeans長のうち短い方を採用
              xch%xNcH   = xch%nH*0.5d0*l_jeans   ! column density of hydrogen


  #ifdef M1CLOSER_FUV_TRANSFER

              ! ---------------------
              !    column density 
              ! ---------------------
    #if FUV_COLUMNDENS_OPTION == 1     
              lfuv    = xch%xlmbdj   ! jeans length         
              xNcH_1  = xch%nH*lfuv                        ! 水素原子の柱密度 [cm^-2] at step i
              xNcH2_1 = xch%ychem(X_H2)*xNcH_1  !水素分子の柱密度の概算 [cm^-2]
      #ifdef CHEM_MODEL_HF2020
              xNcCO_1 = xch%yco*xNcH_1          !CO分子の柱密度 [cm^-2]
      #else
              xNcCO_1 = xch%ychem(X_CO)*xNcH_1          !CO分子の柱密度 [cm^-2]
      #endif

              v_th      = sqrt(cgs_kb*xch%Tg/cgs_mh)  ! [cm s^{-1}]
              dop_b5    = v_th*1.d-5                  ! Doppler broadening parameter

              call selfshield_H2_CO(xch%nH, xch%Tg, xNcH_1, xNcH_1, xNcH2_1, xNcCO_1, xch%fd, dop_b5, fsH2_1, fsCO_1)

    #elif FUV_COLUMNDENS_OPTION == 2


              ! -------------------
              !  H2 column density
              ! -------------------
              dvdr_av = (xch%dvdr(MX)+xch%dvdr(MY)+xch%dvdr(MZ))/3.d0 + 1.d-30  ! lower limit [s^{-1}]
              
              v_th      = sqrt(cgs_kb*xch%Tg/cgs_mh)  ! [cm s^{-1}]
              dop_b5    = v_th*1.d-5                  ! Doppler broadening parameter
              lfuv      = min(v_th/dvdr_av,xch%xlmbdj)  ! minimum of jeans length & sobolev               ! [ cm ]
              xNcH_H2   = xch%nH*lfuv 
              xNcH2_1   = min(minval(NH2_sbv), xch%xlmbdj*xch%ychem(X_H2)*xch%nH)

              v_th    = sqrt(2.d0*cgs_kb*xch%Tg/((MP_AC+MP_AO)*cgs_amu))
              lfuv    = min(v_th/dvdr_av,xch%xlmbdj)                ! [ cm ]
              xNcH_CO = xch%nH*lfuv
      #ifdef CHEM_MODEL_HF2020
              xNcCO_1   = min(minval(NCO_sbv), xch%xlmbdj*xch%yco*xch%nH)
      #else
              xNcCO_1   = min(minval(NCO_sbv), xch%xlmbdj*xch%ychem(X_CO)*xch%nH)
      #endif

              call selfshield_H2_CO(xch%nH, xch%Tg, xNcH_H2, xNcH_CO, xNcH2_1, xNcCO_1, xch%fd, dop_b5, fsH2_1, fsCO_1)

    #else
              print *, "this option for FUV_COLUMNDENS_OPTION has not defined yet"
              stop
    #endif
      
              ! shielding rate
              if(isNotFinite(fsCO_1)) then
                fsCO_1 = 0.d0
              endif

              ! photodissociation rate
              Nph_fuv   = Nrad_FUV(i,j,k)/Unit_l3*MP_PHON !number density of photon in cgs [cm^{-3}]
              Nph_fuv   = max(Nph_fuv, 0.d0)

              Eph_fuv   = Nph_fuv*rsm%erg_FUV*MP_Crd ! energy density of FUV [erg cm^{-3}] 
                                                     ! 高速を制限しているため、流速が滞り光子密度が上がっている分を補正する
              Gfuv_0    = Eph_fuv/Efuv_0             ! Habig's unit

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
              Nph_fuv_co = Nrad_FUV_CO(i,j,k)/Unit_l3*MP_PHON !number density of photon in cgs [cm^{-3}]
              Nph_fuv_co = max(Nph_fuv_co, 0.d0)

              xch%rH2pd = xi_d0*Gfuv_0               ! [s^{-1}]
              xch%rcopd = sig_CO*MP_Ctil*Nph_fuv_co
    #else
              xch%rH2pd = xi_d0 *Gfuv_0         *fsH2_1 ! [s^{-1}] local self shielding fractor
              xch%rcopd = sig_CO*MP_Ctil*Nph_fuv*fsCO_1
    #endif


    #ifndef METAL 
      xch%rcopd = 0.d0
    #endif

  #else
                xch%rH2pd = 0.d0
                xch%rcopd = 0.d0  
  #endif
              !---------------------------------------------------------------!


              ! ----------------------- HI photoionization rate  ------------------!
  #if defined(NO_IONIZATION)
              Nph_euv  = 0.d0
              xch%heat = 0.d0
              xch%rHpi = 0.d0
  #elif defined(M1CLOSER_EUV_TRANSFER)
              Nph_euv  = Nrad(i,j,k)/Unit_l3*MP_PHON    !number density of photon in cgs [cm^{-3}]
              Nph_euv  = max(Nph_euv, 0.d0)
              xch%rHpi = rsm%alpha_EUV*MP_Ctil*Nph_euv  !ionization rate per H atom [s^{-1}]
              xch%heat = rsm%heat_EUV                   ! heating deposit per ionization[erg]
  #else
              Nph_euv  = 0.d0
              xch%rHpi = 0.d0
              xch%heat = 0.d0
  #endif
              ! -----------------------------------------------------------------!
  
              
  #ifdef METAL                      
              ! -------------------------------------!
              !       photoelectric heating 
              !  primordial.F90: PHOTOELECTRIC_HEATING == YES is needed for including PEH
              ! -------------------------------------!

    #if defined(M1CLOSER_FUV_TRANSFER) && defined(M1CLOSER_SEPARATE_FUV_TRANS)
              Nph_fuv_dust = Nrad_FUV_DUST(i,j,k)/Unit_l3*MP_PHON
              xch%rgfuv    = Nph_fuv_dust*rsm%erg_FUV*MP_Crd/Efuv_0 ! G0
    #elif defined(M1CLOSER_FUV_TRANSFER) 
              xch%rgfuv    = Gfuv_0 
    #else
              xch%rgfuv   = 0.d0  !TODO 
    #endif


              ! --------------------------------------!
              !    dust heating rate from EUV & FUV
              ! ---------------------------------------

              ! dust heating rate
    #ifndef SET_NODUST_ATTENUATION
              ! set sigma -------------------
              sig_euv = xch%fd * rsm%sig_EUV  ! 融解 & Metallicity 込み
              sig_fuv = xch%fd * rsm%sig_FUV
              ! -----------------------------

      #if defined(M1CLOSER_FUV_TRANSFER) && defined(M1CLOSER_SEPARATE_FUV_TRANS)
              xch%rdph=(rsm%erg_EUV*sig_euv*Nph_euv + rsm%erg_FUV*sig_fuv*Nph_fuv_dust)*MP_Ctil*Unit_invmu ! ダスト別に計算
      #elif defined(M1CLOSER_FUV_TRANSFER) 
              xch%rdph=(rsm%erg_EUV*sig_euv*Nph_euv + rsm%erg_FUV*sig_fuv*Nph_fuv)*MP_Ctil*Unit_invmu ! FUV一色
      #else
              xch%rdph=rsm%erg_EUV*sig_euv*Nph_euv*MP_Ctil*Unit_invmu 
      #endif

    #else
              ! set sigma -------------------
              sig_euv = 0.d0  
              sig_fuv = 0.d0
              xch%rdph=0.d0
              ! -----------------------------
    #endif

              ! -------------------------------------------
              !         ionization rate of OII
              ! -------------------------------------------
              xch%rOII  = rsm%rOII 

  #else
              xch%rgfuv   = 0.d0
              xch%rdph    = 0.d0
              xch%rOII    = 0.d0
  #endif

              ! set -------------------
              xch%rHmpd = khmpd(i,j,k) ! H- photo-dissociation rate
              ! -----------------------

              ! IR photon density ------
  #ifdef M1CLOSER_IR_TRANSFER
              xch%EradIR = Nrad_IR(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR ! [erg cm^-3]
              xch%EradIR = max(xch%EradIR, 0.d0)
  #else
              xch%EradIR = 0.d0
  #endif
              ! ------------------------


  #ifdef IGNORING_SINK_TRANS_M1
              ! ------------ inside the sink cells ----------------!
              if(sourceCell_flag) then
                Nph_euv   = 0.d0
                xch%rHpi  = 0.d0
                xch%heat  = 0.d0
                xch%rH2pd = 0.d0    
                xch%rcopd = 0.d0
              endif
              ! ---------------------------------------------------!
  #endif






              !-------------        check Nan/Inf     -------------!
              ychem0(:) = xch%ychem(:)
  #ifdef CHEM_MODEL_HF2020
              yco0 = xch%yco
  #endif
              T_K0 = xch%Tg

#if DEBUG_RADTR == YES
              if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                 naninf_flg = .True.
              end if
              do ic=0,NCHEM-1
                 if (isNotFinite(xch%ychem(ic))) then
                    naninf_flg = .True.
                 end if
              end do
  #ifdef CHEM_MODEL_HF2020
              if (isNotFinite(xch%yco)) then
                 naninf_flg = .True.
              endif
  #endif
              if (isNotFinite(xch%rHpi)   .or. &
                  isNotFinite(xch%rH2pd)  .or. &
                  isNotFinite(xch%rHmpd)  .or. &
                  isNotFinite(xch%heat)   .or. &
                  isNotFinite(xch%rgfuv)  .or. &
                  isNotFinite(xch%rdph)   .or. &
                  isNotFinite(xch%rcopd)  .or. &
                  isNotFinite(xch%rOII)   .or. &
                  isNotFinite(xch%xlmbdj) .or. &
                  isNotFinite(xch%fd)     .or. &
                  isNotFinite(xch%EradIR) .or. &
                  isNotFinite(xch%metal)        &
                  ) then
                  naninf_flg = .True.
              endif

              if (naninf_flg) then
                 print '(A, (1P6E15.7))', "(chemistry) NaN/Inf found before chem update: ", xch%ychem(:)
                 print '((1P8E15.7))', xch%nH, xch%Tg, xch%Td, xch%yco, xch%rHpi, xch%rH2pd, xch%rHmpd, xch%heat
                 print '((1P9E15.7))', xch%rgfuv, xch%rdph, xch%rcopd, xch%rOII, xch%xlmbdj, xch%fd, xch%EradIR, xch%metal
                 print *, "stopping..."
                 stop
              end if
#endif      

              !-----------    セルの時間発展 (KS NOTE: heat_hpi/khpi = energy deposit per ionization) -------!
              call ch_CellChemCool(xch)
              !--------------------------------------------------------------------------------------------!

              !-------------        check Nan/Inf     -------------!
              if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                 naninf_flg = .True.
              end if
  #ifdef CHEM_MODEL_HF2020
              if (isNotFinite(xch%yco)) then
                 naninf_flg = .True.
              endif
  #endif
              do ic=0,NCHEM-1
                 if (isNotFinite(xch%ychem(ic))) then
                    naninf_flg = .True.
                 end if
              end do
              if(isNotFinite(xch%Td))then
                 naninf_flg = .True.  
              endif
              if (naninf_flg) then
                 print *, "(chemistry) NaN/Inf found after chem update: ", xch%ychem(:)
  #ifdef CHEM_MODEL_HF2020
                 print *, "yco:", xch%yco
  #endif
                 print *,  xch%nH,xch%Tg,xch%rHpi,xch%heat,x(i),y(j),z(k), xch%dt
                 print *,  ychem0(:)
                 print *,  xch%rdph, xch%rgfuv
                 print *, "Td", xch%Td
                 print *, "stopping..."
                 stop
              end if


              !impose T floor (after chem_update)
              if (xch%Tg < MP_Tmin .or. xch%Tg > MP_Tmax) then
                 !output warning
                 if (count_output < max_count) then
                    print '(A, 6(1PE12.4))', &
                         "*** after chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                         xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                    count_output = count_output+1
                    !suppress warning
                    if (count_output == max_count) then
                       print *, 'count_output reaches max_count: no more warning will be shown'
                       count_output = count_output+1
                    end if
                 end if
              end if
              xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)


              ! ----------------------------------
              !     packing updated chemistry 
              ! ----------------------------------
              do ichem = 0, NCHEM-1
                 chemary3(ichem+NCEHM_MIN)%y(i,j,k) = xch%ychem(ichem)
              enddo 
  #ifdef CHEM_MODEL_HF2020
              yco(i,j,k) = xch%yco
  #endif

              xmu      = get_xmu(xch%ychem) ! 平均分子量
              p(i,j,k) = (cgs_kb*xch%Tg)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e      !温度 [K]


  #ifdef METAL
              Tdust(i,j,k) = xch%Td
  #endif


              !************************************************************
              !  
              !                 Radiation Transport Part 
              !
              !************************************************************

              ! ====================================
              !          EUV components
              ! ====================================

              if(.not. sourceCell_flag) then

  #ifdef M1CLOSER_EUV_TRANSFER

                ! ----------------------------
                !     radiation pressure 
                ! ----------------------------

                ! ----------------------------
                !         absorption 
                ! ----------------------------
    #ifdef NO_IONIZATION
                dtau = 0.d0
    #else
                dtau = (xch%ychem(X_HI)+xch%ychem(X_H2))*xch%nH*rsm%alpha_EUV*MP_Ctil*xch%dt
    #endif
                ! dust absorption ---------------------------------
    #ifndef SET_NODUST_ATTENUATION
                dtau = dtau + xch%nH*sig_euv*MP_Ctil*xch%dt  ! [cm^-3 cm^2 cm s^-1 s ] = [noD] 
    #endif
                dtau_exp = 1.d0/(1.d0+dtau)

    #ifdef RECOM_ION_PHOTO
                Nrad(i,j,k) = (Nrad(i,j,k)+Nph_recom(xch%Tg, xch%nH, xch%ychem(X_HII), xch%ychem(X_EL), xch%dt))*dtau_exp
    #else
                Nrad(i,j,k) = Nrad(i,j,k)*dtau_exp  ! on the spot approximation
    #endif
                Frx(i,j,k)  = Frx(i,j,k)*dtau_exp
                Fry(i,j,k)  = Fry(i,j,k)*dtau_exp
                Frz(i,j,k)  = Frz(i,j,k)*dtau_exp

    !----------EO_added----------!
   ! if(Frx(i,j,k)/= 0.d0 .or. Fry(i,j,k)/=0.d0 .or. Frz(i,j,k)/=0.d0) then
   ! print*, "**********Frxyz(i,j,k) is not zero**********"
   ! end if

   ! if(Nrad(i,j,k) /= 0.d0) then
   ! print*, "**********Nrad(i,j,k) is not zero**********"
   ! end if 
 
   ! if(dtau_exp == 0.d0) then
   ! print*, "**********dtau_exp is zero**********"
   ! end if 
    !----------EO_added----------!


    #ifdef EXTERNALFORCE
                ! parameter for EUV force 
      #ifdef NO_IONIZATION
                cf_EUV = MP_PHON*rsm%erg_EUV*(sig_euv)*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hyd) ! erg cm^2 g^-1 cm^-1 s => cm^3 s^-1
      !----------EO_added----------!
     ! if(cf_EUV == 0.d0) then
     ! print*, "**********cf_EUV is zero *************"
     ! end if
      !----------EO_added----------!

      #else
                cf_EUV = MP_PHON*rsm%erg_EUV*(rsm%alpha_EUV*(xch%ychem(X_HI)+xch%ychem(X_H2))+sig_euv)*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hyd) ! erg cm^2 g^-1 cm^-1 s => cm^3 s^-1
      #endif
                fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx(i,j,k)*cf_EUV ![noD]
                fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry(i,j,k)*cf_EUV ![noD]
                fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz(i,j,k)*cf_EUV ![noD]

      !----------EO_added_willremove----------!
     ! if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_M1CLOSER_EUV_TRANSFER is not zero.******"
     ! end if

     ! if (fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
      !   print*, "******fxyz_hpi(i,j,k)_M1CLOSER_EUV_TRANSFER is zero*******"
      !end if

     ! if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
      !   print*, "******fxyz_hpi(i,j,k)_M1CLOSER_EUV_TRANSFER is not finite******"
     ! end if
      !----------EO_added----------!


    #endif


    #if DEBUG_RADTR == YES
                if (isNotFinite(Nrad(i,j,k)) .or. isNotFinite(Frx(i,j,k)) .or. isNotFinite(Fry(i,j,k)) &
                  .or. isNotFinite(Frz(i,j,k)) ) then
                  print *, "Nan appear at EUV", Nrad(i,j,k), Frx(i,j,k), Fry(i,j,k), Frz(i,j,k)
                  print *, Nph_recom(xch%Tg, xch%nH, xch%ychem(X_HII), xch%ychem(X_EL), xch%dt), dtau_exp
                  stop
                end if
    #endif
  

  #endif !M1CLOSER_EUV_TRANSFER

              ! ====================================
              !          FUV components
              ! ====================================

  #if defined(M1CLOSER_FUV_TRANSFER) && defined(M1CLOSER_SEPARATE_FUV_TRANS)
                !---------------------- H2 LINE TRAPPING -----------------------!

                cs  = sqrt(cgs_kb*xch%Tg/(xmu*cgs_amu))       !等温音速 [cm/s]
                xch%xlmbdj = cs*sqrt(Pi/(cgs_gc*rho_cgs))    !ジーンズ長さ [cm]
                lfuv    = min(xch%xlmbdj, r_FUV)             ! ジーンズ長か光源までの平均距離で短い方を採用する
                xNcH_1  = xch%nH*lfuv                        ! 水素原子の柱密度 [cm^-2] at step i
                xNcH_2  = xch%nH*(lfuv+MP_Ctil*xch%dt) ! 水素原子の柱密度 [cm^-2] at step i+1

                xNcH2_1 = xch%ychem(X_H2)*xNcH_1  !水素分子の柱密度の概算 [cm^-2]
                xNcH2_2 = xch%ychem(X_H2)*xNcH_2  !水素分子の柱密度の概算 [cm^-2] i+1 step

                call selfshield_H2(xNcH2_1, fsH2_1)
                call selfshield_H2(xNcH2_2, fsH2_2)

                dtau = xch%nH*sigd_FUV*MP_Ctil*xch%dt  ! extinction

                if(fsH2_1 > 0.d0) then
                  dtau_exp = fsH2_2/fsH2_1/(1.d0+dtau)
                else
                  dtau_exp = 0.d0
                endif

                Nrad_FUV(i,j,k) = Nrad_FUV(i,j,k)*dtau_exp
                Frx_FUV(i,j,k)  = Frx_FUV(i,j,k) *dtau_exp
                Fry_FUV(i,j,k)  = Fry_FUV(i,j,k) *dtau_exp 
                Frz_FUV(i,j,k)  = Frz_FUV(i,j,k) *dtau_exp
                !---------------------------------------------------------------!

    #if DEBUG_RADTR == YES
                if (isNotFinite(Nrad_FUV(i,j,k)) .or. isNotFinite(Frx_FUV(i,j,k)) .or. isNotFinite(Fry_FUV(i,j,k)) &
                  .or. isNotFinite(Frz_FUV(i,j,k)) ) then
                  print *, "Nan appear at FUV", Nrad_FUV(i,j,k), Frx_FUV(i,j,k), Fry_FUV(i,j,k), Frz_FUV(i,j,k)
                  print *, dtau_exp
                  stop
                end if
    #endif

                !-------------------- CO LINE TRAPPING -------------------------!
                xNcCO_1 = xch%yco*xNcH_1  !CO分子の柱密度の概算 [cm^-2]
                xNcCO_2 = xch%yco*xNcH_2  !CO分子の柱密度の概算 [cm^-2] i+1 step

                call selfshield_CO(xNcH2_1, xNcCO_1, fsCO_1)
                call selfshield_CO(xNcH2_2, xNcCO_2, fsCO_2)

                dtau = dtau*0.9282296651d0

                if(fsCO_1 > 0.d0) then
                  dtau_exp = fsCO_2/fsCO_1/(1.d0+dtau)
                else
                  dtau_exp = 0.d0
                endif

                Nrad_FUV_CO(i,j,k) = Nrad_FUV_CO(i,j,k)*dtau_exp
                Frx_FUV_CO(i,j,k)  = Frx_FUV_CO(i,j,k) *dtau_exp
                Fry_FUV_CO(i,j,k)  = Fry_FUV_CO(i,j,k) *dtau_exp 
                Frz_FUV_CO(i,j,k)  = Frz_FUV_CO(i,j,k) *dtau_exp
                !---------------------------------------------------------------!

    #if DEBUG_RADTR == YES
                if (isNotFinite(Nrad_FUV_CO(i,j,k)) .or. isNotFinite(Frx_FUV_CO(i,j,k)) .or. isNotFinite(Fry_FUV_CO(i,j,k)) &
                  .or. isNotFinite(Frz_FUV_CO(i,j,k)) ) then
                  print *, "Nan appear at FUV_CO", Nrad_FUV_CO(i,j,k), Frx_FUV_CO(i,j,k), Fry_FUV_CO(i,j,k), Frz_FUV_CO(i,j,k)
                  print *, dtau_exp
                  stop
                end if
    #endif

                ! ----------------------------
                !     radiation pressure 
                ! ----------------------------



    #ifndef SET_NODUST_ATTENUATION
                dtau = xch%nH*sig_fuv*MP_Ctil*xch%dt  ! [cm^-3 cm^2 cm s^-1 s ] = [noD] 
    #else
                dtau = 0.d0
    #endif
                dtau_exp = 1.d0/(1.d0+dtau)

                Nrad_FUV_DUST(i,j,k) = Nrad_FUV_DUST(i,j,k)*dtau_exp
                Frx_FUV_DUST(i,j,k)  = Frx_FUV_DUST(i,j,k) *dtau_exp
                Fry_FUV_DUST(i,j,k)  = Fry_FUV_DUST(i,j,k) *dtau_exp 
                Frz_FUV_DUST(i,j,k)  = Frz_FUV_DUST(i,j,k) *dtau_exp
                ! --------------------------------------------------------------!

    #ifdef EXTERNALFORCE
                ! parameter for EUV force 
                cf_FUV = MP_PHON*rsm%erg_FUV*sig_fuv*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hyd) ! erg cm^2 g^-1 cm^-1 s => cm^3 s^-1
                fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx_FUV_DUST(i,j,k)*cf_FUV ![noD]
                fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry_FUV_DUST(i,j,k)*cf_FUV ![noD]
                fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz_FUV_DUST(i,j,k)*cf_FUV ![noD]
    #endif

      !----------EO_added----------!
     ! if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_(M1CLOSER_FUV_TRANSFER,M1CLOSER_SEPARATE_FUV_TRANS) is not zero******"
     ! end if

     ! if(fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
      !   print*, "******fxyz_hpi(i,j,k)_(M1CLOSER_FUV_TRANSFER,M1CLOSER_SEPARATE_FUV_TRANS) is zero*******"
     ! end if


!      if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
!         print*, "******fxyz_hpi(i,j,k)_(M1CLOSER_FUV_TRANSFER,M1CLOSER_SEPARATE_FUV_TRANS) is not finite******"
!      end if
      !----------EO_added----------!



  #elif defined(M1CLOSER_FUV_TRANSFER) 


                ! ----------------------------
                !     radiation pressure 
                !  単純にダストの減光を考慮すればよい
                ! ----------------------------

                !---------------------------------------------------------------!
                dtau     = xch%nH*sigd_FUV*MP_Ctil*xch%dt  ! extinction
                dtau_exp = 1.d0/(1.d0+dtau)
                Nrad_FUV(i,j,k) = Nrad_FUV(i,j,k)*dtau_exp
                Frx_FUV(i,j,k)  = Frx_FUV(i,j,k) *dtau_exp
                Fry_FUV(i,j,k)  = Fry_FUV(i,j,k) *dtau_exp 
                Frz_FUV(i,j,k)  = Frz_FUV(i,j,k) *dtau_exp
                !---------------------------------------------------------------!

    #ifdef EXTERNALFORCE
                ! parameter for EUV force 
                cf_FUV = MP_PHON*rsm%erg_FUV*sig_fuv*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hyd) ! erg cm^2 g^-1 cm^-1 s => cm^3 s^-1
                fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx_FUV(i,j,k)*cf_FUV ![noD]
                fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry_FUV(i,j,k)*cf_FUV ![noD]
                fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz_FUV(i,j,k)*cf_FUV ![noD]
    #endif

      !----------EO_added----------!
!      if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
!         print*, "******fxyz_hpi(i,j,k)_M1CLOSER_FUV_TRANSFER is not zero******"
!      end if

     ! if(fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_M1CLOSER_FUV_TRANSFER is zero*******"
     ! end if


!      if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
!         print*, "******fxyz_hpi(i,j,k)_M1CLOSER_FUV_TRANSFER is not finite******"
!      end if
      !----------EO_added----------!


    #if DEBUG_RADTR == YES
                if (isNotFinite(Nrad_FUV(i,j,k)) .or. isNotFinite(Frx_FUV(i,j,k)) .or. isNotFinite(Fry_FUV(i,j,k)) &
                  .or. isNotFinite(Frz_FUV(i,j,k)) ) then
                  print *, "Nan appear at FUV", Nrad_FUV(i,j,k), Frx_FUV(i,j,k), Fry_FUV(i,j,k), Frz_FUV(i,j,k)
                  print *, dtau_exp
                  stop
                end if
    #endif

  #endif
              endif!if(.not. sourceCell_flag)


              ! ====================================
              !          IR components
              ! ====================================

  #ifdef M1CLOSER_IR_TRANSFER

              ! ----------------------------
              !    dust abundance & opacity
              ! ----------------------------
    #ifdef METAL
      #ifdef DUST_NOTCONSTANT
              rho_dust = rhod(i,j,k)
      #else

        #ifdef METAL_TRANSFER
              rho_dust = fdust_solar*mmetal(i,j,k)
        #else
              rho_dust = fdust_solar*MP_Metallicity
        #endif
      #endif
    #else 
              rho_dust = 0.d0
    #endif  

              Trad = (MP_Crd*xch%EradIR/cgs_asb)**(0.25d0)    ! temperature of radiation fieild of IR photon

              !read opacity
              call pdop_radtr(xch%Td  , rho_dust, xk_dust)
              call pdop_radtr(Trad    , rho_dust, xk_rad )
              call rdop_radtr(Trad    , rho_dust, xk_ros )


              ! ----------------------------
              !     radiation pressure 
              ! ----------------------------
              ! ----------------------------
              !         absorption 
              ! ----------------------------

              ! set Delta rad
              tap1 = rho_cgs*xch%dt*(xk_dust*cgs_c*cgs_asb*xch%Td**4.d0 - xk_rad*MP_Ctil*xch%EradIR) ! [ erg cm^{-3}]
              tap2 = 1.d0 + rho_cgs*xch%dt*xk_rad*MP_Ctil/(1.d0+xch%chi_d) ! [NoD]

              d_EradIR = tap1/tap2 ! [erg cm^{-3}]
              d_NradIR = d_EradIR*Unit_l3/(MP_PHON*MP_hnu_IR) ! [noD]

              Nrad_IR(i,j,k) = Nrad_IR(i,j,k) + d_NradIR      ! update NradIR

              dtau = xk_ros*rho_cgs*MP_Ctil*xch%dt ![noD]
              Frx_IR(i,j,k)  = Frx_IR(i,j,k)/(1.d0+dtau)
              Fry_IR(i,j,k)  = Fry_IR(i,j,k)/(1.d0+dtau)
              Frz_IR(i,j,k)  = Frz_IR(i,j,k)/(1.d0+dtau)

    #ifdef EXTERNALFORCE
              tap1 = MP_PHON*MP_hnu_IR*xk_ros/cgs_c/(Unit_v*Unit_l**2)*dt_code/dt_hyd
              fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx_IR(i,j,k)*tap1 ![noD]
              fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry_IR(i,j,k)*tap1 ![noD]
              fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz_IR(i,j,k)*tap1 ![noD]
    #endif

      !----------EO_added----------!
     ! if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_M1CLOSER_IR_TRANSFER is not zero******"
     ! end if

     ! if(fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_M1CLOSER_IR_TRANSFER is zero*******"
     ! end if


     ! if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
     !    print*, "******fxyz_hpi(i,j,k)_M1CLOSER_IR_TRANSFER is not finite******"
     ! end if
      !----------EO_added----------!


      #if DEBUG_RADTR == YES
                if (isNotFinite(Nrad_IR(i,j,k)) .or. isNotFinite(Frx_IR(i,j,k)) .or. isNotFinite(Fry_IR(i,j,k)) &
                  .or. isNotFinite(Frz_IR(i,j,k)) ) then
                  print *, "Nan appear at IR", Nrad_IR(i,j,k), Frx_IR(i,j,k), Fry_IR(i,j,k), Frz_IR(i,j,k)
                  print *, tap1, tap2, d_EradIR, d_NradIR
                  stop
                end if
      #endif
  #endif

              ! ------------------------------------
              !          dust sublimation
              ! ------------------------------------

  #ifdef DUST_NOTCONSTANT
              ! Tevap を超えたらダスト密度をゆっくり減らす(改善の余地あり)
              Tevap = 2.d3*rho_cgs**(0.0195d0)
              if(xch%Td > Tevap)then

                if(rhod(i,j,k) < 1.d-50)then
                  rhod(i,j,k) = 0.d0
                else
                  t_dyn = 2.228175d9/Unit_t
                  tap1 = dexp(-1.d1*dt_code/t_dyn)
                  rhod(i,j,k) = rhod(i,j,k)*tap1
                endif
              endif
              if(rhod(i,j,k) < 0.d0) rhod(i,j,k) = 0.d0
  #endif
    
  #ifdef EXTERNALFORCE
              if(rsm%num_rad == 0) then
                fx_hpi(i,j,k) = 0.d0
                fy_hpi(i,j,k) = 0.d0
                fz_hpi(i,j,k) = 0.d0
              endif

              if(isNotFinite(fx_hpi(i,j,k))) then
                fx_hpi(i,j,k) = 0.d0
              endif
              if(isNotFinite(fy_hpi(i,j,k))) then
                fy_hpi(i,j,k) = 0.d0
              endif
              if(isNotFinite(fz_hpi(i,j,k))) then
                fz_hpi(i,j,k) = 0.d0
              endif

  #endif

       !----------EO_added----------!
     ! if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_(end_ch_radtrchem) is not zero******"
     ! end if

     ! if(fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_(end_ch_radtrchem) is zero*******"
     ! end if


!      if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
!         print*, "******fxyz_hpi(i,j,k)_(end_ch_radtrchem) is not finite******"
!      end if
      !----------EO_added----------!

           enddo
          enddo
        enddo


      enddo
    enddo

  end subroutine ch_radtrchem
#endif ! M1CLOSER_EUV_TRANSFER


  subroutine check_is_sourceCell(xi, yj, zk, sourceCell_flag)
    real(kind=DBL_KIND),intent(IN) :: xi, yj, zk
    logical, intent(OUT) :: sourceCell_flag
    integer :: nsource_glob, isrc
    real(kind=DBL_KIND),dimension(MX:MZ) :: spos
    real(kind=DBL_KIND) :: r2


    nsource_glob = rs_info%nsource

    sourceCell_flag = .false.

    do isrc = 0, nsource_glob-1
      spos = rs_info%spos(:,isrc)
      r2   = (xi-spos(MX))**2 + (yj-spos(MY))**2 + (zk-spos(MZ))**2
      
      if (r2 < SinkRadius**2) then
        sourceCell_flag = .true.
      endif
    enddo


  end subroutine check_is_sourceCell

#ifdef M1CLOSER_IR_TRANSFER
  !--------------------------------------------
  ! calculate dust temperature and update Erad
  !--------------------------------------------
  subroutine radtr_IRdust(dt_code, dt_hyd)

    real(kind=DBL_KIND), intent(IN) :: dt_code, dt_hyd
    integer :: level, gid
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, Nrad, Frx, Fry, Frz
#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: yco
#endif
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
    real(kind=DBL_KIND) :: yco_l, xnH, xmu, T_K, Td, Erad, dt, rho_cgs, t_dyn
    real(kind=DBL_KIND) :: chi_d, Trad, xk_dust, xk_rad, xk_ros, d_Erad, tap1, tap2, d_Nrad, tau
#ifdef METAL
    real(kind=DBL_KIND) :: rho_dust
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Tdust
  #ifdef DUST_NOTCONSTANT
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
  #endif
#endif
#ifdef EXTERNALFORCE
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
#endif
    real(kind=DBL_KIND) :: Tevap
#ifdef RADTR_DIRECT
    real(kind=DBL_KIND),dimension(:,:,:), pointer :: kdph
    real(kind=DBL_KIND) :: dph
#endif
#ifdef METAL_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
#endif
    integer :: i, j, k, n
    logical :: isNotFinite 

    ! ----------------------------------------------------------
    type chemsp 
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(NCEHM_MIN:NCEHM_MAX)
    integer :: ichem
    ! ----------------------------------------------------------

    ! timestep in cgs
    dt = dt_code * Unit_t

    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        
        ! 子グリッドを持つ場合はchemistryは解かなくてよい (あとでconvergeを呼ぶ)
        if (has_child_grid(gid))   cycle

                  !変数の取り出し
        rho => get_Ucomp(MRHO,gid)
        p => get_Ucomp(MP,gid)

        ! get chemistry pointer
        do ichem = NCEHM_MIN, NCEHM_MAX
          chemary3(ichem)%y => get_Ucomp(ichem,gid)
        enddo
#ifdef CHEM_MODEL_HF2020
        yco  => get_Ucomp(MCO, gid)
#endif

        Nrad => get_Ucomp(MEIR, gid)
        Frx  => get_Ucomp(MFIRX, gid)
        Fry  => get_Ucomp(MFIRY, gid)
        Frz  => get_Ucomp(MFIRZ, gid)
        Tdust => get_Ucomp(MTD,gid)
#ifdef DUST_NOTCONSTANT
        rhod => get_Ucomp(MDRHO,gid)
#endif

#ifdef EXTERNALFORCE
        fx_hpi => get_Ucomp(MXPI,gid)
        fy_hpi => get_Ucomp(MYPI,gid)
        fz_hpi => get_Ucomp(MZPI,gid)
#endif

#ifdef RADTR_DIRECT
        kdph  => get_Ucomp(MDPH,gid)
#endif

#ifdef METAL_TRANSFER
        mmetal => get_Ucomp(MMET,gid)
#endif

        do i=Imin, Imax
          do j=Jmin, Jmax
            do k=Kmin, Kmax

              ! --------------------------
              !   化学組成変数の詰め替え
              ! --------------------------
              do ichem = 0, NCHEM-1
                 ychem(ichem) = chemary3(ichem+NCEHM_MIN)%y(i,j,k)
              enddo 
#ifdef CHEM_MODEL_HF2020
              yco_l = yco(i,j,k)
#endif

              !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
              call adjust_abundance(ychem &
#ifdef CHEM_MODEL_HF2020
                , yco_l &
#endif
#ifdef METAL_TRANSFER
                , mmetal(i,j,k))
#else
                , MP_Metallicity)
#endif

              !derived variables
              rho_cgs = rho(i,j,k)*Unit_rho
              xnH     = rho_cgs/(MP_mu*cgs_amu)                !水素原子核の数密度
              xmu     = get_xmu(ychem) ! 平均分子量
              T_K     = p(i,j,k)*Unit_e*cgs_amu*xmu/rho_cgs/cgs_kb      !温度 [K]

              ! set max and minimum
              T_K = min(max(T_K,MP_Tmin),MP_Tmax)

              ! radiation energy
              Erad = Nrad(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR ! [erg cm^-3]
              Erad = max(Erad, 0.d0)

              ! 直接光の影響
#ifdef RADTR_DIRECT
              dph = kdph(i,j,k)/rho_cgs ! [erg s^-1 cm^-3 * g^-1 cm^3] = [erg g^-1 s^-1]
#endif

              ! ==================================
              ! get dust temperature at each bin 
              ! ==================================
              Td = Tdust(i,j,k)
#ifdef RADTR_DIRECT
              call dtemp_radtr(xnH, T_K, Td, Erad, chi_d, dph)
#else
              call dtemp_radtr(xnH, T_K, Td, Erad, chi_d)
#endif

              ! input dust temperature
              Tdust(i,j,k) = Td
      
              ! temperature of radiation field
              Trad = (MP_Crd*Erad /cgs_asb)**(0.25d0) 

#ifdef DUST_NOTCONSTANT
              rhod(i,j,k)= min( rhod(i,j,k), fdust_solar*MP_Metallicity)
              rho_dust = rhod(i,j,k)
#else
  #ifdef METAL_TRANSFER
              rho_dust = fdust_solar*mmetal(i,j,k)
  #else
              rho_dust = fdust_solar*MP_Metallicity
  #endif
#endif

              ! read opacity
              call pdop_radtr(Td  , rho_dust, xk_dust)
              call pdop_radtr(Trad, rho_dust, xk_rad )
              call rdop_radtr(Trad, rho_dust, xk_ros )


#ifdef EXTERNALFORCE
              tap1 = MP_PHON*MP_hnu_IR*xk_ros/cgs_c/(Unit_v*Unit_l**2)*dt_code/dt_hyd ! erg * cm^2 g^-1 /(cm s^-1) = cm^3 s^-1 = (cm s^-1) cm^2
              fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx(i,j,k)*tap1 ![noD]
              fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry(i,j,k)*tap1 ![noD]
              fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz(i,j,k)*tap1 ![noD]
#endif

      !----------EO_added----------!
     ! if (fx_hpi(i,j,k)/=0.d0 .or. fy_hpi(i,j,k)/=0.d0 .or. fz_hpi(i,j,k)/=0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_(radtr_IRdust) is not zero******"
     ! end if

     ! if(fx_hpi(i,j,k)==0.d0 .or. fy_hpi(i,j,k)==0.d0 .or. fz_hpi(i,j,k)==0.d0) then
     !    print*, "******fxyz_hpi(i,j,k)_(radtr_IRdust) is zero*******"
     ! end if


     ! if (isNotFinite(fx_hpi(i,j,k)) .or. isNotFinite(fy_hpi(i,j,k)) .or. isNotFinite(fz_hpi(i,j,k))) then
     !    print*, "******fxyz_hpi(i,j,k)_(radtr_IRdust) is not finite******"
     ! end if
      !----------EO_added----------!


              ! set Delta rad
              tap1 = rho_cgs*dt*(xk_dust*cgs_c*cgs_asb*Td**4.d0 - xk_rad*MP_Ctil*Erad) ! [ erg cm^{-3}]
              tap2 = 1.d0 + rho_cgs*dt*xk_rad*MP_Ctil/(1.d0+chi_d) ! [NoD]

              d_Erad = tap1/tap2 ! [erg cm^{-3}]
              d_Nrad = d_Erad*Unit_l3/(MP_PHON*MP_hnu_IR) ! [noD]

              ! update
              Nrad(i,j,k) = Nrad(i,j,k) + d_Nrad

              tau = xk_ros*rho_cgs*MP_Ctil*dt ![noD]
              Frx(i,j,k)  = Frx(i,j,k)/(1.d0+tau)
              Fry(i,j,k)  = Fry(i,j,k)/(1.d0+tau)
              Frz(i,j,k)  = Frz(i,j,k)/(1.d0+tau)

#if DEBUG_RADTR == YES
              if (isNotFinite(Nrad(i,j,k)) .or. isNotFinite(Frx(i,j,k)) .or. isNotFinite(Fry(i,j,k)) &
                .or. isNotFinite(Frz(i,j,k)) ) then
                print *, "Nan appear at IR", Nrad(i,j,k), Frx(i,j,k), Fry(i,j,k), Frz(i,j,k)
                print *, tap1, tap2, d_Erad, d_Nrad
                stop
              end if
#endif


#ifdef DUST_NOTCONSTANT
              ! Tevap を超えたらダスト密度をゆっくり減らす(改善の余地あり)
              Tevap = 2.d3*rho_cgs**(0.0195d0)
              if(Td > Tevap)then

                if(rhod(i,j,k) < 1.d-50)then
                  rhod(i,j,k) = 0.d0
                else
                  t_dyn = 2.228175d9/Unit_t 
                  tap1 = dexp(-1.d1*dt_code/t_dyn)
                  rhod(i,j,k) = rhod(i,j,k)*tap1
                endif
              endif
#endif

            enddo ! k=Kmin, Kmax
          enddo ! j=Jmin, Jmax
        enddo ! i=Imin, Imax


      enddo ! n = Gidmin, GidListMax( level )
    enddo ! level = Lmin, Lmax


  
  end subroutine radtr_IRdust
#endif ! M1CLOSER_IR_TRANSFER




  ! --------------------------
  !   plank dust opacity
  ! -------------------------- 
  subroutine pdop_radtr(Td, rhod, xkd)

    real(kind=DBL_KIND),intent(IN)  :: Td, rhod 
    real(kind=DBL_KIND),intent(OUT) :: xkd
    real(kind=DBL_KIND) :: f_dust, xop

    f_dust = rhod/(fdust_solar) 
    call find_dop(Td, xop)

    xkd = f_dust*xop

  end subroutine pdop_radtr

  ! --------------------------
  !   rosseland dust opacity
  ! -------------------------- 
  subroutine rdop_radtr(Td, rhod, xkd)

    real(kind=DBL_KIND),intent(IN)  :: Td, rhod 
    real(kind=DBL_KIND),intent(OUT) :: xkd
    real(kind=DBL_KIND) :: f_dust, xop

    f_dust = rhod/(fdust_solar) 
    call find_dros(Td, xop)

    xkd = f_dust*xop

  end subroutine rdop_radtr


  !-------------------------------------------------------------------------
  ! Get Hm photo-detachemet rate [s^-1] assuming optically thin NIR flux
  !-------------------------------------------------------------------------
  subroutine GetKHmpdThin(x,y,z,khmpd)
    real(kind=DBL_KIND),dimension(:),pointer,intent(IN) :: x, y, z
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK),intent(OUT) :: khmpd

    integer :: isrc,i,j,k
    real(kind=DBL_KIND) :: r2,r_star

    khmpd(:,:,:)=0.d0

    return !今は考慮しない 

    do isrc=0, rs_info%nsource-1
       r_star = sqrt(rs_info%lum(isrc)/(4.*Pi * cgs_sigma*rs_info%Trad(isrc)**4)) !radius
       do i = Imin, Imax
          do j = Jmin, Jmax
             do k = Kmin, Kmax
                r2 = (x(i) - rs_info%spos(MX,isrc))**2 + &
                     (y(j) - rs_info%spos(MY,isrc))**2 + &
                     (z(k) - rs_info%spos(MZ,isrc))**2
                khmpd(i,j,k)=khmpd(i,j,k)+rs_info%hhm(isrc)*r_star**2/(r2*Unit_l**2) 
             end do
          end do
       end do
    end do

  end subroutine GetKHmpdThin

  !---------------------------------------------------------------------
  ! number of recombination photons
  !--------------------------------------------------------------------

  function Nph_recom(Tg, xnH, yHII, ye, dt)
    
    real(kind=DBL_KIND) :: Nph_recom
    real(kind=DBL_KIND) :: Tg, xnH, yHII, ye, dt
    real(kind=DBL_KIND) :: alpha_A, alpha_B, inv_T, nHII, ne

    inv_T = 1.d0/Tg

    alpha_A = 1.269d-13 * (315614.d0*inv_T)**1.503 *(1+(604625.d0*inv_T)**0.470)**(-1.923)  !Case A [cm^3 s^-1]
    alpha_B = 2.753d-14 * (315614.d0*inv_T)**1.5   *(1+(115188.d0*inv_T)**0.407)**(-2.242)  !Case B

    nHII = xnH * yHII
    ne   = xnH * ye

    Nph_recom = (alpha_A - alpha_B)*nHII*ne*dt/MP_PHON *Unit_l**3.d0 ! [cm^-3] => [nD]

    return

  end function Nph_recom

  !-------------------------------------------------------------
  !             self-shielding rate for H2 & CO
  !-------------------------------------------------------------

  subroutine selfshield_H2_CO(nH, Tg, NcH_H2, NcH_CO, NcH2, NcCO, fd, dop_b5, fH2, fCO)

    real(kind=DBL_KIND), intent(IN) :: nH, Tg, NcH_H2, NcH_CO, NcH2, NcCO, fd, dop_b5
    real(kind=DBL_KIND), intent(OUT):: fH2, fCO
    real(kind=DBL_KIND),parameter :: NcH2_cr = 1d14, NcH2_max = 1d22
    real(kind=DBL_KIND) :: x,y,z ! temp vars
    real(kind=DBL_KIND) :: logNco, logtheta, theta1, logNH2, theta2, theta3, Av_H2, Av_CO
    logical :: isNotFinite 
#if SELF_SHIELDING_H2_FORMULA == 2
    real(kind=DBL_KIND) :: A1, A2, alpha, x1_sqrt, log10T
#endif

    ! extinction
    Av_H2 = 5.34d-22*NcH_H2*fd
    Av_CO = 5.34d-22*NcH_CO*fd

    ! extinction is too much ---
    if(Av_H2 > 20.d0) then
      fH2 = 0.d0
      fCO = 0.d0
      return
    endif
    ! --------------------------


    ! shelf_shielding rate of H2 -------------
#if SELF_SHIELDING_H2_FORMULA == 1
    if(NcH2 > NcH2_max) then
      fH2 = 0.d0
    else if (NcH2 < NcH2_cr) then
      fH2 = 1.d0
    else
      x = sqrt(NcH2_cr/NcH2)
      y = sqrt(x)
      fH2 = x*y
    endif
#elif SELF_SHIELDING_H2_FORMULA == 2
    ! formula of Wolcott-Green and Haiman 2019
    if(NcH2 > 1.d24 ) then
      fH2 = 0.d0
    else if (NcH2 < 1.d13) then
      fH2 = 1.d0
    else
      
      log10T  =  log10(Tg)
      A1      =  0.8711d0*log10T-1.928d0
      A2      = -0.9639d0*log10T+3.892d0
      alpha   =  A1*exp(-0.2856d0*log10(nH))+A2
      x       =  NcH2*2.d0*1.d-15
      x1_sqrt =  sqrt(1.d0+x)

      fH2 = 0.965d0/(1.d0+x/dop_b5)**alpha + 0.035d0/x1_sqrt*exp(-8.5d-4*x1_sqrt)
    endif
#else
    print *, "this option of SELF_SHIELDING_H2_FORMULA has not been defined yet"
    stop
#endif
    ! ----------------------------------------

#ifdef M1CLOSER_SEPARATE_FUV_TRANS
    fH2 = fH2*dexp(-4.18d0*Av_H2) 
#endif

    ! shelf_shielding rate of CO -------------
    if(NcCO > 1.d12) then
        logNco = dlog10(NcCO)
        logtheta = (((((-1.62507133d-4*logNco+1.43588538d-2)*logNco &
                    -5.21418820d-1)*logNco+9.95253514d0)*logNco-1.05308265d2)*logNco &
                    +5.85860800d2)*logNco-1.33950326d3
        theta1 = 10.d0**logtheta
    else
        theta1 = 1.d0
    endif

    ! Theta2
    if(NcH2 > 1.e13) then
        logNH2 = dlog10(NcH2)
        logtheta = -2.09263955d-18*dexp(1.89035939*logNH2)
        theta2 = 10.d0**logtheta
    else
        theta2 = 1.d0
    endif

#ifdef M1CLOSER_SEPARATE_FUV_TRANS
    theta3 = dexp(-3.88*Av_CO)
    fCO = theta1*theta2*theta3
#else
    fCO = theta1*theta2
#endif


    ! ----------------------------------------

 
  end subroutine selfshield_H2_CO



  subroutine selfshield_H2(NcH2, fH2)

    real(kind=DBL_KIND), intent(IN) :: NcH2
    real(kind=DBL_KIND), intent(OUT):: fH2
    real(kind=DBL_KIND),parameter :: NcH2_cr = 1d14, NcH2_max = 1d22
    real(kind=DBL_KIND) :: x,y,z ! temp vars

    ! shelf_shielding rate of H2 -------------
    if(NcH2 > NcH2_max) then
      fH2 = 0.d0
    else if (NcH2 < NcH2_cr) then
      fH2 = 1.d0
    else
      x = sqrt(NcH2_cr/NcH2)
      y = sqrt(x)
      fH2 = x*y
    endif
    ! ----------------------------------------

  end subroutine selfshield_H2

  subroutine selfshield_CO(NcH2, NcCO, fCO)

    real(kind=DBL_KIND), intent(IN) :: NcH2, NcCO
    real(kind=DBL_KIND), intent(OUT):: fCO 
    real(kind=DBL_KIND) :: logNH2,logNco, logtheta, theta1, theta2
    logical :: isNotFinite

    ! shelf_shielding rate of CO -------------
    if(NcCO > 1.d12) then
        logNco = dlog10(NcCO)
        logtheta = (((((-1.62507133d-4*logNco+1.43588538d-2)*logNco &
                    -5.21418820d-1)*logNco+9.95253514d0)*logNco-1.05308265d2)*logNco &
                    +5.85860800d2)*logNco-1.33950326d3
        theta1 = 10.d0**logtheta
    else
        theta1 = 1.d0
    endif

    ! Theta2
    if(NcH2 > 1.e13) then
        logNH2 = dlog10(NcH2)
        logtheta = -2.09263955d-18*dexp(1.89035939*logNH2)
        theta2 = 10.d0**logtheta
    else
        theta2 = 1.d0
    endif
    fCO = theta1*theta2

    if(isNotFinite(fCO)) then
      fCO = 0.d0
    endif

  end subroutine selfshield_CO 

end module radtr_chem
