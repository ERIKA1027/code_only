module radtr_chem
  use unit
  use parameter, only : Pi
  use modelParameter
  use grid
  use mpilib
  use sinkParticle, only : sp_getSinkRadius, sp_getLevel
  use radiationSource
  use primordial,only: CoolSolverExplicit, CoolSolverImplicit, adjust_abundance &
    , yHe, c_H2, dbg_flg_prm, chem_vals, chem_vals, get_xmu
  use kinzoku, only : dtemp_radtr, find_dop, find_dros, fdust_solar
  implicit none
  private
  real(kind=8) :: SinkRadius
  type(t_rs_info),save,pointer :: rs_info 
  type rad_mean
    real(kind=8) :: alpha_EUV 
    real(kind=8) :: heat_EUV 
    real(kind=8) :: sig_EUV 
    real(kind=8) :: sig_FUV 
    real(kind=8) :: rOII 
    real(kind=8) :: erg_EUV 
    real(kind=8) :: erg_FUV 
    integer :: num_rad 
  end type rad_mean
  public :: ch_radtrchem, rad_mean
  public :: radtr_IRdust
contains
  subroutine ch_radtrchem(dt_code, dt_hyd, rsm)
    use modelParameter,only : MP_Tmin,MP_Tmax
    use chemistry, only : ch_CellChemCool
    real(kind=8),intent(IN) :: dt_code,dt_hyd
    type(rad_mean), intent(IN) :: rsm
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:),pointer :: vx, vy, vz
    real(kind=8),dimension(:,:,:),pointer :: rho, p
    real(kind=8),dimension(:,:,:),pointer :: yco
    real(kind=8) :: l_jeans, l_sobolev, dop_b5, cf_FUV
    real(kind=8),dimension(0:2) :: dvdr
    real(kind=8), dimension(6) :: dvdr_s, NH2_sbv, lsb_s, NCO_sbv, lsb_rho_s
    real(kind=8),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
    real(kind=8),dimension(:,:,:),pointer :: Nrad_FUV, Frx_FUV, Fry_FUV, Frz_FUV
    real(kind=8) :: lfuv, fsH2_1, fsH2_2, fsCO_1, fsCO_2, Gfuv_0
    real(kind=8),parameter :: sig_CO = 1.594880d-18 
    real(kind=8),parameter :: sigd_FUV = 2.23212d-21
    real(kind=8),parameter :: Efuv_0 = 5.29d-14 
    real(kind=8),parameter :: xi_d0 = 4.27d-11 
    real(kind=8),dimension(:,:,:),pointer :: Nrad_IR, Frx_IR, Fry_IR, Frz_IR
    real(kind=8),dimension(:,:,:),pointer :: Tdust
    real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax) :: mu, khmpd
    real(kind=8),dimension(0:6 -1) :: ychem0
    real(kind=8) :: xmu, cs, T_K0, yco0, sig_euv, sig_fuv
    real(kind=8) :: Nph_euv, Nph_fuv, Eph_fuv, dtau, dtau_exp, rho_cgs, xNcH2_1, xNcCO_1, xNcH2_2, xNcCO_2
    real(kind=8) :: xNcH_1, xNcH_2, v_th, dvdr_av, xNcH_H2, xNcH_CO
    integer :: level, n, gid, i, j, k, ic
    real(kind=8) :: Trad, rho_dust, xk_dust, xk_rad, xk_ros, tap1, tap2, d_EradIR, d_NradIR
    real(kind=8),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    real(kind=8) :: cf_EUV
    integer :: sp_Level
    logical :: isNotFinite
    integer :: count_output=0, max_count=100 
    logical :: naninf_flg
    logical :: sourceCell_flag
    type(chem_vals) :: xch
    real(kind=8), dimension(Imin-1:Imax+1,Jmin-1:Jmax+1,Kmin-1:Kmax+1) :: yh2_old, yco_old
    type chemsp
      real(kind=8),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(5:10)
    integer :: ichem
    call rs_GetSourceInfo(rs_info)
    SinkRadius = sp_getSinkRadius()
    sp_Level = sp_getLevel()
    xch%dt = dt_code*Unit_t
    naninf_flg = .False.
    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        if (has_child_grid(gid)) cycle
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)
        vx => get_Ucomp(1, gid)
        vy => get_Ucomp(2, gid)
        vz => get_Ucomp(3, gid)
        rho => get_Ucomp(0,gid)
        p => get_Ucomp(4,gid)
        do ichem = 5, 10
          chemary3(ichem)%y => get_Ucomp(ichem,gid)
        enddo
        yco => get_Ucomp(11, gid)
        Tdust => get_Ucomp(21,gid)
        fx_hpi => get_Ucomp(14,gid)
        fy_hpi => get_Ucomp(15,gid)
        fz_hpi => get_Ucomp(16,gid)
        Nrad => get_Ucomp(22, gid)
        Frx => get_Ucomp(23, gid)
        Fry => get_Ucomp(24, gid)
        Frz => get_Ucomp(25, gid)
        Nrad_FUV => get_Ucomp(26, gid)
        Frx_FUV => get_Ucomp(27, gid)
        Fry_FUV => get_Ucomp(28, gid)
        Frz_FUV => get_Ucomp(29, gid)
        Nrad_IR => get_Ucomp(30, gid)
        Frx_IR => get_Ucomp(31, gid)
        Fry_IR => get_Ucomp(32, gid)
        Frz_IR => get_Ucomp(33, gid)
        call GetKHmpdThin(x,y,z,khmpd) 
        do k=Kmin-1, Kmax+1
          do j=Jmin-1, Jmax+1
            do i=Imin-1, Imax+1
              yh2_old(i,j,k) = chemary3(6)%y(i,j,k)
              yco_old(i,j,k) = yco(i,j,k)
            enddo
          enddo
        enddo
        do i=Imin, Imax
          do j=Jmin, Jmax
            do k=Kmin, Kmax
              xch%metal = MP_Metallicity
              xch%Td = Tdust(i,j,k)
              do ichem = 0, 6 -1
                 xch%ychem(ichem) = chemary3(ichem+5)%y(i,j,k)
              enddo
              xch%yco = yco(i,j,k)
              call adjust_abundance(xch%ychem &
                ,xch%yco &
                ,xch%metal)
              rho_cgs= rho(i,j,k)*Unit_rho
              xch%nH = rho_cgs*Unit_invmu 
              xmu = get_xmu(xch%ychem) 
              xch%Tg = p(i,j,k)*Unit_e*cgs_amu*xmu /rho_cgs/cgs_kb 
              if (xch%Tg < MP_Tmin*0.999d0 .or. xch%Tg > MP_Tmax*1.001d0) then 
                if (count_output < max_count) then
                   print '(A, 6(1PE12.4))', &
                        "*** before chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                        xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                   count_output = count_output+1
                   if (count_output == max_count) then
                      print *, 'count_output reaches max_count: no more warning will be shown'
                      count_output = count_output+1
                   end if
                end if
              end if
              xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)
              if (level .eq. sp_Level) then
                call check_is_sourceCell(x(i), y(j), z(k), sourceCell_flag)
              else
                sourceCell_flag = .false.
              endif
              xch%fd = MP_Metallicity
              cs = sqrt(cgs_kb*xch%Tg/(xmu*cgs_amu)) 
              l_jeans = cs*sqrt(Pi/(cgs_gc*rho_cgs)) 
              dvdr_s(1) = dabs(vx(i+1,j,k)-vx(i,j,k))/(CellWidth(0,level)*Unit_t)+1.d-30 
              dvdr_s(2) = dabs(vx(i,j,k)-vx(i-1,j,k))/(CellWidth(0,level)*Unit_t)+1.d-30 
              dvdr_s(3) = dabs(vy(i,j+1,k)-vy(i,j,k))/(CellWidth(1,level)*Unit_t)+1.d-30 
              dvdr_s(4) = dabs(vy(i,j,k)-vy(i,j-1,k))/(CellWidth(1,level)*Unit_t)+1.d-30 
              dvdr_s(5) = dabs(vz(i,j,k+1)-vz(i,j,k))/(CellWidth(2,level)*Unit_t)+1.d-30 
              dvdr_s(6) = dabs(vz(i,j,k)-vz(i,j,k-1))/(CellWidth(2,level)*Unit_t)+1.d-30 
              v_th = sqrt(cgs_kb*xch%Tg/cgs_mh)
              lsb_s(:) = v_th/dvdr_s(:)
              NH2_sbv(1) = min(yh2_old(i+1,j,k), yh2_old(i,j,k))*lsb_s(1)
              NH2_sbv(2) = min(yh2_old(i,j,k), yh2_old(i-1,j,k))*lsb_s(2)
              NH2_sbv(3) = min(yh2_old(i,j+1,k), yh2_old(i,j,k))*lsb_s(3)
              NH2_sbv(4) = min(yh2_old(i,j,k), yh2_old(i,j-1,k))*lsb_s(4)
              NH2_sbv(5) = min(yh2_old(i,j,k+1), yh2_old(i,j,k))*lsb_s(5)
              NH2_sbv(6) = min(yh2_old(i,j,k), yh2_old(i,j,k-1))*lsb_s(6)
              v_th = sqrt(2.d0*cgs_kb*xch%Tg/((MP_AC+MP_AO)*cgs_amu))
              lsb_s(:) = v_th/dvdr_s(:)
              NCO_sbv(1) = min(yco_old(i+1,j,k), yco_old(i,j,k))*lsb_s(1)
              NCO_sbv(2) = min(yco_old(i,j,k), yco_old(i-1,j,k))*lsb_s(2)
              NCO_sbv(3) = min(yco_old(i,j+1,k), yco_old(i,j,k))*lsb_s(3)
              NCO_sbv(4) = min(yco_old(i,j,k), yco_old(i,j-1,k))*lsb_s(4)
              NCO_sbv(5) = min(yco_old(i,j,k+1), yco_old(i,j,k))*lsb_s(5)
              NCO_sbv(6) = min(yco_old(i,j,k), yco_old(i,j,k-1))*lsb_s(6)
              xch%dvdr(0) = (dvdr_s(1)+dvdr_s(2))*0.5d0
              xch%dvdr(1) = (dvdr_s(3)+dvdr_s(4))*0.5d0
              xch%dvdr(2) = (dvdr_s(5)+dvdr_s(6))*0.5d0
              xch%xlmbdj = l_jeans 
              xch%xNcH = xch%nH*0.5d0*l_jeans 
              dvdr_av = (xch%dvdr(0)+xch%dvdr(1)+xch%dvdr(2))/3.d0 + 1.d-30 
              v_th = sqrt(cgs_kb*xch%Tg/cgs_mh) 
              dop_b5 = v_th*1.d-5 
              lfuv = min(v_th/dvdr_av,xch%xlmbdj) 
              xNcH_H2 = xch%nH*lfuv
              xNcH2_1 = min(minval(NH2_sbv), xch%xlmbdj*xch%ychem(1)*xch%nH)
              v_th = sqrt(2.d0*cgs_kb*xch%Tg/((MP_AC+MP_AO)*cgs_amu))
              lfuv = min(v_th/dvdr_av,xch%xlmbdj) 
              xNcH_CO = xch%nH*lfuv
              xNcCO_1 = min(minval(NCO_sbv), xch%xlmbdj*xch%yco*xch%nH)
              call selfshield_H2_CO(xch%nH, xch%Tg, xNcH_H2, xNcH_CO, xNcH2_1, xNcCO_1, xch%fd, dop_b5, fsH2_1, fsCO_1)
              if(isNotFinite(fsCO_1)) then
                fsCO_1 = 0.d0
              endif
              Nph_fuv = Nrad_FUV(i,j,k)/Unit_l3*MP_PHON 
              Nph_fuv = max(Nph_fuv, 0.d0)
              Eph_fuv = Nph_fuv*rsm%erg_FUV*MP_Crd 
              Gfuv_0 = Eph_fuv/Efuv_0 
              xch%rH2pd = xi_d0 *Gfuv_0 *fsH2_1 
              xch%rcopd = sig_CO*MP_Ctil*Nph_fuv*fsCO_1
              Nph_euv = Nrad(i,j,k)/Unit_l3*MP_PHON 
              Nph_euv = max(Nph_euv, 0.d0)
              xch%rHpi = rsm%alpha_EUV*MP_Ctil*Nph_euv 
              xch%heat = rsm%heat_EUV 
              xch%rgfuv = Gfuv_0
              sig_euv = xch%fd * rsm%sig_EUV 
              sig_fuv = xch%fd * rsm%sig_FUV
              if(xch%Td<=MP_Tsubl) then
                sig_euv = xch%fd * rsm%sig_EUV 
                sig_fuv = xch%fd * rsm%sig_FUV
              else
              sig_euv = 0.d0 
              sig_fuv = 0.d0
              endif
              xch%rdph=(rsm%erg_EUV*rsm%sig_EUV*Nph_euv + rsm%erg_FUV*rsm%sig_FUV*Nph_fuv)*MP_Ctil*Unit_invmu 
              xch%rOII = rsm%rOII
              xch%rHmpd = khmpd(i,j,k) 
              xch%EradIR = Nrad_IR(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR 
              xch%EradIR = max(xch%EradIR, 0.d0)
              ychem0(:) = xch%ychem(:)
              yco0 = xch%yco
              T_K0 = xch%Tg
              call dtemp_radtr(xch%nH, xch%Tg, xch%Td, xch%EradIR, xch%chi_d ,dph_in=xch%rdph)
              if(xch%Td>MP_Tsubl) then
                xch%fd=0.d0
              endif
              call ch_CellChemCool(xch)
              if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                 naninf_flg = .True.
              end if
              if (isNotFinite(xch%yco)) then
                 naninf_flg = .True.
              endif
              do ic=0,6 -1
                 if (isNotFinite(xch%ychem(ic))) then
                    naninf_flg = .True.
                 end if
              end do
              if(isNotFinite(xch%Td))then
                 naninf_flg = .True.
              endif
              if (naninf_flg) then
                 print *, "(chemistry) NaN/Inf found after chem update: ", xch%ychem(:)
                 print *, "yco:", xch%yco
                 print *, xch%nH,xch%Tg,xch%rHpi,xch%heat,x(i),y(j),z(k), xch%dt
                 print *, ychem0(:)
                 print *, xch%rdph, xch%rgfuv
                 print *, "Td", xch%Td
                 print *, "stopping..."
                 stop
              end if
              if (xch%Tg < MP_Tmin .or. xch%Tg > MP_Tmax) then
                 if (count_output < max_count) then
                    print '(A, 6(1PE12.4))', &
                         "*** after chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                         xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                    count_output = count_output+1
                    if (count_output == max_count) then
                       print *, 'count_output reaches max_count: no more warning will be shown'
                       count_output = count_output+1
                    end if
                 end if
              end if
              xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)
              do ichem = 0, 6 -1
                 chemary3(ichem+5)%y(i,j,k) = xch%ychem(ichem)
              enddo
              yco(i,j,k) = xch%yco
              xmu = get_xmu(xch%ychem) 
              p(i,j,k) = (cgs_kb*xch%Tg)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
              Tdust(i,j,k) = xch%Td
              if(.not. sourceCell_flag) then
                dtau = (xch%ychem(0)+xch%ychem(1))*xch%nH*rsm%alpha_EUV*MP_Ctil*xch%dt
                if(xch%Td<=MP_Tsubl) then
                  dtau = dtau + xch%nH*sig_euv*MP_Ctil*xch%dt 
                else
                  dtau = dtau
                endif
                dtau_exp = 1.d0/(1.d0+dtau)
                Nrad(i,j,k) = (Nrad(i,j,k)+Nph_recom(xch%Tg, xch%nH, xch%ychem(3), xch%ychem(2), xch%dt))*dtau_exp
                Frx(i,j,k) = Frx(i,j,k)*dtau_exp
                Fry(i,j,k) = Fry(i,j,k)*dtau_exp
                Frz(i,j,k) = Frz(i,j,k)*dtau_exp
      if(xch%Td<=MP_Tsubl) then
        cf_EUV = MP_PHON*rsm%erg_EUV*(rsm%alpha_EUV*(xch%ychem(0)+xch%ychem(1))+sig_euv)*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_co&
&de/dt_hyd) 
      else
        cf_EUV = MP_PHON*rsm%erg_EUV*(rsm%alpha_EUV*(xch%ychem(0)+xch%ychem(1)))*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hy&
&d) 
      endif
                fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx(i,j,k)*cf_EUV 
                fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry(i,j,k)*cf_EUV 
                fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz(i,j,k)*cf_EUV 
                if(xch%Td<=MP_Tsubl) then
                  dtau = xch%nH*sigd_FUV*MP_Ctil*xch%dt 
                else
                  dtau = 0.d0
                endif
                dtau_exp = 1.d0/(1.d0+dtau)
                Nrad_FUV(i,j,k) = Nrad_FUV(i,j,k)*dtau_exp
                Frx_FUV(i,j,k) = Frx_FUV(i,j,k) *dtau_exp
                Fry_FUV(i,j,k) = Fry_FUV(i,j,k) *dtau_exp
                Frz_FUV(i,j,k) = Frz_FUV(i,j,k) *dtau_exp
                if(xch%Td<=MP_Tsubl) then
                  cf_FUV = MP_PHON*rsm%erg_FUV*sig_fuv*Unit_invmu/cgs_c/(Unit_v*Unit_l**2)*(dt_code/dt_hyd) 
                else
                  cf_FUV = 0.d0
                endif
                fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx_FUV(i,j,k)*cf_FUV 
                fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry_FUV(i,j,k)*cf_FUV 
                fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz_FUV(i,j,k)*cf_FUV 
              endif
              rho_dust = fdust_solar*MP_Metallicity
              Trad = (MP_Crd*xch%EradIR/cgs_asb)**(0.25d0) 
              call pdop_radtr(xch%Td , rho_dust, xk_dust)
              call pdop_radtr(Trad , rho_dust, xk_rad )
              call rdop_radtr(Trad , rho_dust, xk_ros )
                if(xch%Td<=MP_Tsubl) then
                  tap1 = rho_cgs*xch%dt*(xk_dust*cgs_c*cgs_asb*xch%Td**4.d0 - xk_rad*MP_Ctil*xch%EradIR) 
                else
                  tap1 = 0.d0
                endif
              tap2 = 1.d0 + rho_cgs*xch%dt*xk_rad*MP_Ctil/(1.d0+xch%chi_d) 
              d_EradIR = tap1/tap2 
              d_NradIR = d_EradIR*Unit_l3/(MP_PHON*MP_hnu_IR) 
              Nrad_IR(i,j,k) = Nrad_IR(i,j,k) + d_NradIR 
              if(xch%Td<=MP_Tsubl) then
                dtau = xk_ros*rho_cgs*MP_Ctil*xch%dt 
              else
                dtau = 0.d0
              endif
              Frx_IR(i,j,k) = Frx_IR(i,j,k)/(1.d0+dtau)
              Fry_IR(i,j,k) = Fry_IR(i,j,k)/(1.d0+dtau)
              Frz_IR(i,j,k) = Frz_IR(i,j,k)/(1.d0+dtau)
              tap1 = MP_PHON*MP_hnu_IR*xk_ros/cgs_c/(Unit_v*Unit_l**2)*dt_code/dt_hyd
                if(xch%Td<=MP_Tsubl) then
                  tap1 = MP_PHON*MP_hnu_IR*xk_ros/cgs_c/(Unit_v*Unit_l**2)*dt_code/dt_hyd
                else
                  tap1 = 0.d0
                endif
              fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx_IR(i,j,k)*tap1 
              fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry_IR(i,j,k)*tap1 
              fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz_IR(i,j,k)*tap1 
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
           enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine ch_radtrchem
  subroutine check_is_sourceCell(xi, yj, zk, sourceCell_flag)
    real(kind=8),intent(IN) :: xi, yj, zk
    logical, intent(OUT) :: sourceCell_flag
    integer :: nsource_glob, isrc
    real(kind=8),dimension(0:2) :: spos
    real(kind=8) :: r2
    nsource_glob = rs_info%nsource
    sourceCell_flag = .false.
    do isrc = 0, nsource_glob-1
      spos = rs_info%spos(:,isrc)
      r2 = (xi-spos(0))**2 + (yj-spos(1))**2 + (zk-spos(2))**2
      if (r2 < SinkRadius**2) then
        sourceCell_flag = .true.
      endif
    enddo
  end subroutine check_is_sourceCell
  subroutine radtr_IRdust(dt_code, dt_hyd)
    real(kind=8), intent(IN) :: dt_code, dt_hyd
    integer :: level, gid
    real(kind=8),dimension(:,:,:),pointer :: rho, p, Nrad, Frx, Fry, Frz
    real(kind=8),dimension(:,:,:),pointer :: yco
    real(kind=8),dimension(0:6 -1) :: ychem
    real(kind=8) :: yco_l, xnH, xmu, T_K, Td, Erad, dt, rho_cgs, t_dyn
    real(kind=8) :: chi_d, Trad, xk_dust, xk_rad, xk_ros, d_Erad, tap1, tap2, d_Nrad, tau
    real(kind=8) :: rho_dust
    real(kind=8),dimension(:,:,:),pointer :: Tdust
    real(kind=8),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    real(kind=8) :: Tevap
    integer :: i, j, k, n
    logical :: isNotFinite
    type chemsp
      real(kind=8),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(5:10)
    integer :: ichem
    dt = dt_code * Unit_t
    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        if (has_child_grid(gid)) cycle
        rho => get_Ucomp(0,gid)
        p => get_Ucomp(4,gid)
        do ichem = 5, 10
          chemary3(ichem)%y => get_Ucomp(ichem,gid)
        enddo
        yco => get_Ucomp(11, gid)
        Nrad => get_Ucomp(30, gid)
        Frx => get_Ucomp(31, gid)
        Fry => get_Ucomp(32, gid)
        Frz => get_Ucomp(33, gid)
        Tdust => get_Ucomp(21,gid)
        fx_hpi => get_Ucomp(14,gid)
        fy_hpi => get_Ucomp(15,gid)
        fz_hpi => get_Ucomp(16,gid)
        do i=Imin, Imax
          do j=Jmin, Jmax
            do k=Kmin, Kmax
              do ichem = 0, 6 -1
                 ychem(ichem) = chemary3(ichem+5)%y(i,j,k)
              enddo
              yco_l = yco(i,j,k)
              call adjust_abundance(ychem &
                , yco_l &
                , MP_Metallicity)
              rho_cgs = rho(i,j,k)*Unit_rho
              xnH = rho_cgs/(MP_mu*cgs_amu) 
              xmu = get_xmu(ychem) 
              T_K = p(i,j,k)*Unit_e*cgs_amu*xmu/rho_cgs/cgs_kb 
              T_K = min(max(T_K,MP_Tmin),MP_Tmax)
              Erad = Nrad(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR 
              Erad = max(Erad, 0.d0)
              Td = Tdust(i,j,k)
              call dtemp_radtr(xnH, T_K, Td, Erad, chi_d)
              Tdust(i,j,k) = Td
              Trad = (MP_Crd*Erad /cgs_asb)**(0.25d0)
              rho_dust = fdust_solar*MP_Metallicity
              call pdop_radtr(Td , rho_dust, xk_dust)
              call pdop_radtr(Trad, rho_dust, xk_rad )
              call rdop_radtr(Trad, rho_dust, xk_ros )
              tap1 = MP_PHON*MP_hnu_IR*xk_ros/cgs_c/(Unit_v*Unit_l**2)*dt_code/dt_hyd 
              fx_hpi(i,j,k)=fx_hpi(i,j,k)+Frx(i,j,k)*tap1 
              fy_hpi(i,j,k)=fy_hpi(i,j,k)+Fry(i,j,k)*tap1 
              fz_hpi(i,j,k)=fz_hpi(i,j,k)+Frz(i,j,k)*tap1 
              tap1 = rho_cgs*dt*(xk_dust*cgs_c*cgs_asb*Td**4.d0 - xk_rad*MP_Ctil*Erad) 
              tap2 = 1.d0 + rho_cgs*dt*xk_rad*MP_Ctil/(1.d0+chi_d) 
              d_Erad = tap1/tap2 
              d_Nrad = d_Erad*Unit_l3/(MP_PHON*MP_hnu_IR) 
              Nrad(i,j,k) = Nrad(i,j,k) + d_Nrad
              tau = xk_ros*rho_cgs*MP_Ctil*dt 
              Frx(i,j,k) = Frx(i,j,k)/(1.d0+tau)
              Fry(i,j,k) = Fry(i,j,k)/(1.d0+tau)
              Frz(i,j,k) = Frz(i,j,k)/(1.d0+tau)
            enddo 
          enddo 
        enddo 
      enddo 
    enddo 
  end subroutine radtr_IRdust
  subroutine pdop_radtr(Td, rhod, xkd)
    real(kind=8),intent(IN) :: Td, rhod
    real(kind=8),intent(OUT) :: xkd
    real(kind=8) :: f_dust, xop
    f_dust = rhod/(fdust_solar)
    call find_dop(Td, xop)
    xkd = f_dust*xop
  end subroutine pdop_radtr
  subroutine rdop_radtr(Td, rhod, xkd)
    real(kind=8),intent(IN) :: Td, rhod
    real(kind=8),intent(OUT) :: xkd
    real(kind=8) :: f_dust, xop
    f_dust = rhod/(fdust_solar)
    call find_dros(Td, xop)
    xkd = f_dust*xop
  end subroutine rdop_radtr
  subroutine GetKHmpdThin(x,y,z,khmpd)
    real(kind=8),dimension(:),pointer,intent(IN) :: x, y, z
    real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax),intent(OUT) :: khmpd
    integer :: isrc,i,j,k
    real(kind=8) :: r2,r_star
    khmpd(:,:,:)=0.d0
    return 
    do isrc=0, rs_info%nsource-1
       r_star = sqrt(rs_info%lum(isrc)/(4.*Pi * cgs_sigma*rs_info%Trad(isrc)**4)) 
       do i = Imin, Imax
          do j = Jmin, Jmax
             do k = Kmin, Kmax
                r2 = (x(i) - rs_info%spos(0,isrc))**2 + &
                     (y(j) - rs_info%spos(1,isrc))**2 + &
                     (z(k) - rs_info%spos(2,isrc))**2
                khmpd(i,j,k)=khmpd(i,j,k)+rs_info%hhm(isrc)*r_star**2/(r2*Unit_l**2)
             end do
          end do
       end do
    end do
  end subroutine GetKHmpdThin
  function Nph_recom(Tg, xnH, yHII, ye, dt)
    real(kind=8) :: Nph_recom
    real(kind=8) :: Tg, xnH, yHII, ye, dt
    real(kind=8) :: alpha_A, alpha_B, inv_T, nHII, ne
    inv_T = 1.d0/Tg
    alpha_A = 1.269d-13 * (315614.d0*inv_T)**1.503 *(1+(604625.d0*inv_T)**0.470)**(-1.923) 
    alpha_B = 2.753d-14 * (315614.d0*inv_T)**1.5 *(1+(115188.d0*inv_T)**0.407)**(-2.242) 
    nHII = xnH * yHII
    ne = xnH * ye
    Nph_recom = (alpha_A - alpha_B)*nHII*ne*dt/MP_PHON *Unit_l**3.d0 
    return
  end function Nph_recom
  subroutine selfshield_H2_CO(nH, Tg, NcH_H2, NcH_CO, NcH2, NcCO, fd, dop_b5, fH2, fCO)
    real(kind=8), intent(IN) :: nH, Tg, NcH_H2, NcH_CO, NcH2, NcCO, fd, dop_b5
    real(kind=8), intent(OUT):: fH2, fCO
    real(kind=8),parameter :: NcH2_cr = 1d14, NcH2_max = 1d22
    real(kind=8) :: x,y,z 
    real(kind=8) :: logNco, logtheta, theta1, logNH2, theta2, theta3, Av_H2, Av_CO
    logical :: isNotFinite
    real(kind=8) :: A1, A2, alpha, x1_sqrt, log10T
    Av_H2 = 5.34d-22*NcH_H2*fd
    Av_CO = 5.34d-22*NcH_CO*fd
    if(Av_H2 > 20.d0) then
      fH2 = 0.d0
      fCO = 0.d0
      return
    endif
    if(NcH2 > 1.d24 ) then
      fH2 = 0.d0
    else if (NcH2 < 1.d13) then
      fH2 = 1.d0
    else
      log10T = log10(Tg)
      A1 = 0.8711d0*log10T-1.928d0
      A2 = -0.9639d0*log10T+3.892d0
      alpha = A1*exp(-0.2856d0*log10(nH))+A2
      x = NcH2*2.d0*1.d-15
      x1_sqrt = sqrt(1.d0+x)
      fH2 = 0.965d0/(1.d0+x/dop_b5)**alpha + 0.035d0/x1_sqrt*exp(-8.5d-4*x1_sqrt)
    endif
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
    fCO = theta1*theta2
  end subroutine selfshield_H2_CO
  subroutine selfshield_H2(NcH2, fH2)
    real(kind=8), intent(IN) :: NcH2
    real(kind=8), intent(OUT):: fH2
    real(kind=8),parameter :: NcH2_cr = 1d14, NcH2_max = 1d22
    real(kind=8) :: x,y,z 
    if(NcH2 > NcH2_max) then
      fH2 = 0.d0
    else if (NcH2 < NcH2_cr) then
      fH2 = 1.d0
    else
      x = sqrt(NcH2_cr/NcH2)
      y = sqrt(x)
      fH2 = x*y
    endif
  end subroutine selfshield_H2
  subroutine selfshield_CO(NcH2, NcCO, fCO)
    real(kind=8), intent(IN) :: NcH2, NcCO
    real(kind=8), intent(OUT):: fCO
    real(kind=8) :: logNH2,logNco, logtheta, theta1, theta2
    logical :: isNotFinite
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
    fCO = theta1*theta2
    if(isNotFinite(fCO)) then
      fCO = 0.d0
    endif
  end subroutine selfshield_CO
end module radtr_chem
