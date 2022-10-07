module modelParameter
  implicit none
  private
  real(kind=8),save :: MP_ElapseLimit ,MP_T_LAST, MP_N0, MP_R0, MP_Tau0, MP_T0, MP_spNcr, MP_spRadius_lamJ,MP_SubDiskNcr
  integer,save :: MP_Dstep, MP_Lmax0, MP_JeansConst, MP_spRadius_cell, MP_GasType, MP_SubDiskAbs,MP_SubAbsCell,MP_SubStrmAbs
  real(kind=8),save :: MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_spCs,MP_Boxsize, MP_Gconst, MP_Mstar, MP_Lstar, MP_Mdot, MP_Bondi_rad&
&ius, Bondi_radius,sound_speed, yhni, yh2i, yhpi, yhmi, yh2pi, yeli, ycoi, MP_BHLII_radius,MP_BHLI_radius, MP_T20sink, MP_dg, MP_lu&
&m_edd, MP_lum_edd_dg, MP_lum_star, MP_kappa_pl, MP_Tsubl, MP_sinksubl
  integer,save :: MP_CONNECTION_RUN
  real(kind=8) :: MP_Ctil, MP_Ctil_nD, MP_PHON, MP_Crd
  real(kind=8),parameter :: yHe = 9.7222222d-2
  real(kind=8),parameter :: MP_pi = 3.14159265359d0
  real(kind=8),parameter :: MP_kappa_es = 4.0d-1
  real(kind=8),save :: MP_Metallicity
  real(kind=8),save :: MP_frac_C_solar, MP_frac_O_solar, MP_AC, MP_AO &
    ,MP_frac_C, MP_frac_O, MP_fracC_AC, MP_fracO_AO, MP_mu, MP_frac_COsum
  public :: MP_ElapseLimit ,MP_T_LAST, MP_N0, MP_R0, MP_Tau0, MP_T0, MP_spNcr, MP_spRadius_lamJ,MP_SubDiskNcr
  public :: MP_Dstep, MP_Lmax0, MP_JeansConst, MP_spRadius_cell, MP_GasType, MP_SubDiskAbs,MP_SubAbsCell,MP_SubStrmAbs
  public :: MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_spCs, MP_Boxsize, MP_Gconst, MP_Mstar, MP_Lstar, MP_Mdot, MP_Bondi_radius, Bondi&
&_radius, sound_speed, MP_BHLII_radius, MP_BHLI_radius, MP_T20sink, MP_kappa_pl, MP_Tsubl
  public :: modelParameter_init
  public :: MP_CONNECTION_RUN
  public :: MP_Ctil, MP_Ctil_nD, MP_PHON, MP_Crd
  public :: MP_Metallicity
  public :: MP_frac_C_solar, MP_frac_O_solar, MP_frac_C, MP_frac_O, MP_AC, MP_AO &
    , MP_fracC_AC, MP_fracO_AO, MP_mu, MP_frac_COsum
contains
  subroutine modelParameter_init
    use mpilib
    use io_util, only : readenv, print_msg, read_env
    use parameter, only : Pi, Pi4, Pi4i
    use unit
    real(kind=8) :: h0, hmax, csp, jlength, sp_Rhocr
    call print_msg('start initializing modelParameters (prolbem: radtest)')
    if ( &
         readenv('ELAPSELIMIT', MP_ElapseLimit) .and. &
         readenv('Dstep', MP_Dstep) .and. &
         readenv('T_LAST', MP_T_LAST) .and. &
         readenv('N0', MP_N0) .and. &
         readenv('T0', MP_T0) .and. &
         readenv('Lmax0', MP_Lmax0) .and. &
         readenv('sp_radius_in_cell', MP_spRadius_cell) .and. &
         readenv('sp_Ncr' , MP_spNcr) .and. &
         readenv('Mstar', MP_Mstar) .and. &
         readenv('Mdot', MP_Mdot) .and. &
         readenv('Boxsize', MP_Boxsize) .and. &
         readenv('GasType', MP_GasType) .and. &
         readenv('Tsubl', MP_Tsubl) .and. &
         readenv('Lstar', MP_Lstar) &
         ) then
       if (get_myrank() == 0) then
          print *, 'ELAPSELIMIT       =', MP_ElapseLimit, '[h]' 
          print *, 'Dstep             =', MP_Dstep 
          print *, 'T_LAST            =', MP_T_LAST, '[yr]' 
          print *, 'N0                =', MP_N0, '[cm^-3]' 
          print *, 'T0                =', MP_T0, '[K]' 
          print *, 'Lmax0             =', MP_Lmax0 
          print *, 'sp_radius_in_cell =', MP_spRadius_cell 
          print *, 'sp_Ncr            =', MP_spNcr, '[cm^-3]'
          print *, 'Mstar             =', MP_Mstar, '[M_sun]' 
          print *, 'Mdot              =', MP_Mdot, '[M_sun/yr]' 
          print *, 'Boxsize           =', MP_Boxsize, '[BHL]' 
          print *, 'GasType           =', MP_GasType 
          print *, 'Tsubl             =', MP_Tsubl, '[K]'
          print *, 'Lstar             =', MP_Lstar
       endif
    else
       print *, '****', 'error in modelParameter_init.'
       stop
    end if
    call flush(6)
    MP_Tmax = 1D5 
    MP_Tmin = 1D1 
    MP_Nmin = 1D-2 
    MP_Vmax = 2D3 
    MP_spCs = 1.9D5 
    MP_JeansConst=-1 
    MP_spRadius_lamJ = 0.5d0
   MP_Bondi_radius = 1.4d5 * cgs_au / Unit_l
   MP_Mstar = MP_Mstar*cgs_msun / Unit_m
   MP_BHLII_radius = 9.36d-2 * cgs_pc / Unit_l
   MP_BHLI_radius = 2.1d-1 * cgs_pc / Unit_l
   MP_Boxsize = MP_Boxsize * MP_BHLI_radius
      if (readenv('Metallicity', MP_Metallicity)) then
         if(get_myrank() == 0) then
           print *, 'Metallicity       =', MP_Metallicity 
         endif
       else
         if(get_myrank() == 0) print *, 'metallicity is not defined'
       endif
   if(MP_Tsubl<1.3d3) then
     MP_kappa_pl = 5.880285d2*1.d-2 
   else
     MP_kappa_pl = 6.535156d2*1.d-2 
   end if
    MP_Gconst = cgs_gc*Unit_rho*Unit_t**2 
    if (readenv('SubDiskAbs', MP_SubDiskAbs)) then
       if (get_myrank() == 0) &
            print *, 'SubDiskAbs       =', MP_SubDiskAbs
    else
       MP_SubDiskAbs = 1 
       if (get_myrank() == 0) &
            print *, 'SubDiskAbs       =', MP_SubDiskAbs, ' (Default)'
    end if
    if (MP_SubDiskAbs >= 1) then
       if (MP_SubDiskAbs == 1) then
          MP_SubDiskNcr = 1d-1 * MP_spNcr 
       else if (MP_SubDiskAbs == 2) then
          MP_SubDiskNcr = 1d-2 * MP_spNcr 
       end if
       if (get_myrank() == 0) &
            print *, 'SubDiskNcr       =', MP_SubDiskNcr, '[cm^-3] (Default)'
    end if
    if (readenv('SubStrmAbs', MP_SubStrmAbs)) then
       if (get_myrank() == 0) &
            print *, 'SubStrmAbs       =', MP_SubStrmAbs
    else
       MP_SubStrmAbs = 1 
       if (get_myrank() == 0) &
            print *, 'SubStrmAbs       =', MP_SubStrmAbs, ' (Default)'
    end if
    if (MP_SubDiskAbs >= 1 .or. MP_SubStrmAbs >= 1) then
       MP_SubAbsCell = 2
       if (get_myrank() == 0) &
            print *, 'SubAbsCell       =', MP_SubAbsCell, ' (Default)'
    end if
    MP_T_LAST = MP_T_LAST / Unit_yr
    if (readenv('CONNECTION_RUN', MP_CONNECTION_RUN)) then
       if (get_myrank() == 0) &
            print *, 'CONNECTION_RUN       =', MP_CONNECTION_RUN
    else
       MP_CONNECTION_RUN = 0 
       if (get_myrank() == 0) &
            print *, 'CONNECTION_RUN       =', MP_CONNECTION_RUN, ' (Default)'
    end if
     MP_Crd = 6.d-3 
    MP_Ctil = MP_Crd*cgs_c
    MP_Ctil_nD = MP_Ctil/Unit_v
    MP_PHON = 1.d49
    MP_frac_C_solar = 0.927d-4
    MP_frac_O_solar = 3.568d-4
    MP_frac_C = MP_frac_C_solar * MP_Metallicity
    MP_frac_O = MP_frac_O_solar * MP_Metallicity
    MP_frac_COsum = MP_frac_C + MP_frac_O
    MP_AC = 12.011d0
    MP_AO = 15.999d0
    MP_fracC_AC = MP_frac_C * MP_AC
    MP_fracO_AO = MP_frac_O * MP_AO
    MP_mu = 1.008d0 + 4.004d0*yHe + MP_fracC_AC + MP_fracO_AO
    call print_msg('end initializing modelParameters')
    call flush(6)
  end subroutine modelParameter_init
end module modelParameter
