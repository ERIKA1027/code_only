#include "config.h"
!-------------------------------------------------------------------------
! module for set-up of parameters in radiation test
!-------------------------------------------------------------------------
module modelParameter
  implicit none
  private
  real(kind=DBL_KIND),save :: MP_ElapseLimit ,MP_T_LAST, MP_N0, MP_R0, MP_Tau0, MP_T0, MP_spNcr, MP_spRadius_lamJ,MP_SubDiskNcr
  integer,save :: MP_Dstep, MP_Lmax0, MP_JeansConst, MP_spRadius_cell, MP_GasType, MP_SubDiskAbs,MP_SubAbsCell,MP_SubStrmAbs
  real(kind=DBL_KIND),save :: MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_spCs,MP_Boxsize, MP_Gconst, MP_Mstar, MP_Lstar, MP_Mdot, MP_Bondi_radius, Bondi_radius,sound_speed, yhni, yh2i, yhpi, yhmi, yh2pi, yeli, ycoi, MP_BHLII_radius,MP_BHLI_radius, MP_T20sink, MP_dg, MP_lum_edd, MP_lum_edd_dg, MP_lum_star, MP_kappa_pl, MP_Tsubl, MP_sinksubl
  integer,save :: MP_CONNECTION_RUN
  real(kind=DBL_KIND) :: MP_Ctil, MP_Ctil_nD, MP_PHON, MP_Crd
  real(kind=DBL_KIND),parameter :: yHe      = 9.7222222d-2 
  real(kind=DBL_KIND),parameter :: MP_pi    = 3.14159265359d0
  real(kind=DBL_KIND),parameter :: MP_kappa_es    = 4.0d-1


  !-------------------------------- chemistry --------------------------------!
#ifndef METAL_TRANSFER
  real(kind=DBL_KIND),save :: MP_Metallicity
#endif
  real(kind=DBL_KIND),save :: MP_frac_C_solar, MP_frac_O_solar, MP_AC, MP_AO &
    ,MP_frac_C, MP_frac_O, MP_fracC_AC, MP_fracO_AO, MP_mu, MP_frac_COsum
  !---------------------------------------------------------------------------!

  public :: MP_ElapseLimit ,MP_T_LAST, MP_N0, MP_R0, MP_Tau0, MP_T0, MP_spNcr, MP_spRadius_lamJ,MP_SubDiskNcr
  public :: MP_Dstep, MP_Lmax0, MP_JeansConst, MP_spRadius_cell, MP_GasType, MP_SubDiskAbs,MP_SubAbsCell,MP_SubStrmAbs
  
  public :: MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_spCs, MP_Boxsize, MP_Gconst, MP_Mstar, MP_Lstar, MP_Mdot, MP_Bondi_radius, Bondi_radius, sound_speed, MP_BHLII_radius, MP_BHLI_radius, MP_T20sink, MP_kappa_pl, MP_Tsubl
  public :: modelParameter_init
  public :: MP_CONNECTION_RUN
  public :: MP_Ctil, MP_Ctil_nD, MP_PHON, MP_Crd

  !-------------------------------- chemistry --------------------------------!
#ifndef METAL_TRANSFER
  public :: MP_Metallicity
#endif
  public :: MP_frac_C_solar, MP_frac_O_solar, MP_frac_C, MP_frac_O, MP_AC, MP_AO &
    , MP_fracC_AC, MP_fracO_AO, MP_mu, MP_frac_COsum
  !---------------------------------------------------------------------------!
contains
  subroutine modelParameter_init
    use mpilib
    use io_util, only : readenv, print_msg, read_env
    use parameter, only :  Pi, Pi4, Pi4i
    use unit
    !use primordial

    real(kind=DBL_KIND) :: h0, hmax, csp, jlength, sp_Rhocr
!    real(kind=DBL_KIND) :: yhni, yh2i, yeli, yhpi, yhmi, yh2pi, ycoi
!    real(kind=DBL_KIND) :: yhni2, yh2i2, yeli2, yhpi2, yhmi2, yh2pi2, rho02, cs02, radius
 
    call print_msg('start initializing modelParameters (prolbem: radtest)')

    ! ------------------------------
    ! read from enviroment variables
    ! ------------------------------
    if ( &
         readenv('ELAPSELIMIT', MP_ElapseLimit) .and. &
         readenv('Dstep', MP_Dstep) .and. &
         readenv('T_LAST', MP_T_LAST) .and. &
         readenv('N0', MP_N0) .and. &
        ! readenv('R0', MP_R0) .and. &
        ! readenv('Tau0', MP_Tau0) .and. &
         readenv('T0', MP_T0) .and. &
         readenv('Lmax0', MP_Lmax0) .and. &         
         readenv('sp_radius_in_cell', MP_spRadius_cell) .and.  &
         readenv('sp_Ncr' , MP_spNcr) .and. &
         readenv('Mstar', MP_Mstar) .and. &
         readenv('Mdot', MP_Mdot)  .and. &
         readenv('Boxsize', MP_Boxsize) .and. &
         readenv('GasType', MP_GasType) .and. &
         readenv('Tsubl', MP_Tsubl) .and. &
!         readenv('sinksubl', MP_sinksubl) .and. &
         readenv('Lstar', MP_Lstar) &
         ) then
       if (get_myrank() == PRIMARY_RANK) then
          print *, 'ELAPSELIMIT       =', MP_ElapseLimit, '[h]' ! maximum computational time in hour
          print *, 'Dstep             =', MP_Dstep              ! data output interval
          print *, 'T_LAST            =', MP_T_LAST, '[yr]'     ! maximum simulation time in yr
          print *, 'N0                =', MP_N0, '[cm^-3]'      ! central density of cloud in cm^-3
         ! print *, 'R0                =', MP_R0,  '[g/cm3]'
         ! print *, 'Tau0              =', MP_Tau0
          print *, 'T0                =', MP_T0, '[K]'          ! temperature of cloud in K
          print *, 'Lmax0             =', MP_Lmax0              ! initial maximum nest level
          print *, 'sp_radius_in_cell =', MP_spRadius_cell      ! radius of sink particles wrt finest cell size
          print *, 'sp_Ncr            =', MP_spNcr, '[cm^-3]' 
          print *, 'Mstar             =', MP_Mstar, '[M_sun]'   ! mass of protostar
          print *, 'Mdot              =', MP_Mdot, '[M_sun/yr]' ! mdot of protostar
!          print *, 'Boxsize           =', MP_Boxsize, '[au]'    ! half of side length of computational box in au
          print *, 'Boxsize           =', MP_Boxsize, '[BHL]'    ! half of side length of computational box in bhl
          print *, 'GasType           =', MP_GasType            ! type of ambient gas: 0 -> HI, 1 -> H2, 2 -> H2&disk-like
          print *, 'Tsubl             =', MP_Tsubl, '[K]'
!          print *, 'sinksubl          =', MP_sinksubl
          print *, 'Lstar             =', MP_Lstar
       endif
    else
       print *, '****', 'error in modelParameter_init.'
       stop
    end if
    call flush(6)


    ! ------------------------------
    ! parameters determined here
    ! ------------------------------
    MP_Tmax = 1D5     ! 計算中に課す最大温度 [K]
    MP_Tmin = 1D1     ! 計算中に課す最小温度 [K]
    MP_Nmin = 1D-2    ! 計算中に課す最小数密度 [cm-3]
    MP_Vmax = 2D3     ! 計算中に課す最大速度 [km/s]
    MP_spCs = 1.9D5   ! ジーンズ長からシンク半径を求める際に仮定する音速 [cm/s] (sound speed for gas at nH = 1e12 cm^-3 (yH2=0.5, T=1e3 K) is ~ 1.9 km/s)

    !dummy
    MP_JeansConst=-1             !Jeans長を何波長でという条件でrefineしない
    MP_spRadius_lamJ = 0.5d0
  
    !-------------------------------------------
    ! Bondi_raius for Boxsize 
    !------------------------------------------- 

   !----------EO_added----------!
   ! Bondi_radius=2.d0*cgs_gc*MP_Mstar*cgs_msun / ((2.d6)*(2.d6)) 
   ! MP_Bondi_radius=Bondi_radius / Unit_l
   !MP_Boxsize = MP_Boxsize*MP_Bondi_radius

   MP_Bondi_radius = 1.4d5 * cgs_au / Unit_l
 !  MP_Boxsize      = MP_Boxsize * MP_Bondi_radius
  
   !Setting Sink cell!
 !  h0 = 2d0*MP_Boxsize / ((NGI_BASE) * (NI))     !ベースグリッドのセルサイズ
 !  hmax = h0/(2d0**MP_Lmax0)                     !最大レベルでのセルサイズ
 !  MP_spRadius_cell = int((MP_Bondi_radius/3.d1) / hmax) 

  ! print *, 'sp_radius_in_cell =', MP_spRadius_cell      ! radius of sink particles wrt finest cell size 
   !----------EO_added----------!

   MP_Mstar = MP_Mstar*cgs_msun / Unit_m  
 
   !----------EO_added----------!
   MP_BHLII_radius = 9.36d-2 * cgs_pc / Unit_l
   MP_BHLI_radius  = 2.1d-1 * cgs_pc / Unit_l
   MP_Boxsize      = MP_Boxsize * MP_BHLI_radius

   !dust_sublimation
!    MP_dustsubl = (2.0**(0.5)) * 8.4d-3 * cgs_pc / Unit_l
!    MP_dustsubl = 5.52d-3 * cgs_pc / Unit_l
   #ifndef METAL_TRANSFER
      if (readenv('Metallicity', MP_Metallicity)) then
         if(get_myrank() == PRIMARY_RANK) then
           print *, 'Metallicity       =', MP_Metallicity    ! Metallicity
         endif
       else
         if(get_myrank() == PRIMARY_RANK) print *, 'metallicity is not defined'
       endif
   #endif


  ! MP_dg          = 1.0 /(1.0+710.0*MP_Metallicity)
  ! MP_lum_edd     = 4.d0 * MP_pi * cgs_c * cgs_gc * MP_Mstar * Unit_m / MP_kappa_es ![cgs]
  ! MP_lum_edd_dg     = MP_dg * MP_lum_edd ![cgs]
  ! MP_lum_star       = MP_Lstar * MP_lum_edd_dg ![cgs]

   if(MP_Tsubl<1.3d3) then
    ! MP_kappa_pl = 5.935908d2*1.d-2  ![cgs], solormetal, 1200K
     MP_kappa_pl = 5.880285d2*1.d-2  ![cgs], solormetal, 1000K
   else
     MP_kappa_pl = 6.535156d2*1.d-2  ![cgs], solormetal, 1500K
   end if

   !昇華半径は金属度依存性なし(kappa_plと2.8d2で打ち消し合う)。ただし、光度には金属度依存性入っている。
!   MP_dustsubl = (MP_lum_star * 2.8d2 / (16.0 * MP_pi * MP_kappa_pl * cgs_sigma *((MP_Tsubl)**(4.d0)) ))**(0.5)
!   MP_dustsubl = MP_dustsubl / Unit_l 


  ! MP_T20sink = 2.24d-2 * cgs_pc / Unit_l

   !Setting Sink cell!
!   h0 = 2d0*MP_Boxsize / ((NGI_BASE) * (NI))     !ベースグリッドのセルサイズ
!   hmax = h0/(2d0**MP_Lmax0)                     !最大レベルでのセルサイズ
!   MP_spRadius_cell = int((MP_sinksubl*MP_dustsubl) / hmax) 
!   MP_spRadius_cell = 10  
   !----------EO_added----------!



    !sound_speed
!    MP_csiso = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
!       /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))    ! isothermal sound speed

    !Bondi radius
!    MP_Bondi_radius=2.d0*MP_Gconst*MP_Mstar*cgs_msun / (MP_csiso**2 + (1.d1*Unit_v)**2) 
    !MP_Mstar
!    MP_Mstar = MP_Mstar*cgs_msun / Unit_m  
   

    ! gconst
    MP_Gconst = cgs_gc*Unit_rho*Unit_t**2 ! (=Pi4i)

    !sink粒子がちょうどLmax0に作成されるように MP_spNcr をfine tune
!    h0 = 2d0*MP_Boxsize / ((NGI_BASE) * (NI))     !ベースグリッドのセルサイズ
!    hmax = h0/(2d0**MP_Lmax0)                     !最大レベルでのセルサイズ
!    csp = MP_spCs/Unit_v                          !sink半径決定時に仮定される音速
!    jlength= MP_spRadius_cell * hmax / MP_spRadius_lamJ !ジーンズ長がこの値ならLmax0がシンク粒子のレベルになる
!    sp_Rhocr = Pi / MP_Gconst * (csp/jlength)**2
!    MP_spNcr = sp_Rhocr / (cgs_mh * (1.+4.*cgs_yHe) / Unit_rho)

    !-----------------------------------------------------
    ! whether subgrid model for disk model of absoprtion is used
    ! 0 -> no subgrid disk absoprtion,
    ! 1 -> activate subgrid disk with critical density 1/10 of sink density (default),
    ! 2 -> activate subgrid disk with critical density 1/100 of sink density
    !-----------------------------------------------------
    if (readenv('SubDiskAbs', MP_SubDiskAbs)) then 
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubDiskAbs       =', MP_SubDiskAbs
    else
       MP_SubDiskAbs = 1 ! Default value
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubDiskAbs       =', MP_SubDiskAbs, ' (Default)'
    end if
       
    ! set critical density of subgrid disk
    if (MP_SubDiskAbs >= 1) then
       if (MP_SubDiskAbs == 1) then
          MP_SubDiskNcr = 1d-1 * MP_spNcr ! ひとまず、sink粒子密度の1/10をサブグリッド円盤吸収の閾値とする
       else if (MP_SubDiskAbs == 2) then
          MP_SubDiskNcr = 1d-2 * MP_spNcr ! sink粒子密度の1/100をサブグリッド円盤吸収の閾値とする
       end if
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubDiskNcr       =', MP_SubDiskNcr, '[cm^-3] (Default)'
    end if

    !-----------------------------------------------------
    ! switch for subgrid model of EUV absoprtion based on Stromgren radius
    ! 0 -> de-activate subgrid EUV absoprtion,
    ! 1 -> activate subgrid EUV absorption
    !-----------------------------------------------------
    if (readenv('SubStrmAbs', MP_SubStrmAbs)) then 
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubStrmAbs       =', MP_SubStrmAbs
    else
       MP_SubStrmAbs = 1 ! Default value
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubStrmAbs       =', MP_SubStrmAbs, ' (Default)'
    end if
    
    ! set number of buffer cells for subgrid absorption
    if (MP_SubDiskAbs >= 1 .or. MP_SubStrmAbs >= 1) then
       MP_SubAbsCell = 2
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'SubAbsCell       =', MP_SubAbsCell, ' (Default)'
    end if

    ! --------------------------------------------
    ! convert from physical unit to code unit
    ! --------------------------------------------
    MP_T_LAST = MP_T_LAST / Unit_yr


    !------------------------------------------------------------------------------
    ! flag to activate the mode to connect two different resultion runs
    !    (mdot_disk is assumed to be constant to avoid artificial enhancement)
    !  0 -> de-activate (default), >1 -> # of steps for transitional period
    !------------------------------------------------------------------------------
    if (readenv('CONNECTION_RUN', MP_CONNECTION_RUN)) then 
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'CONNECTION_RUN       =', MP_CONNECTION_RUN
    else
       MP_CONNECTION_RUN = 0 ! Default value
       if (get_myrank() == PRIMARY_RANK) &
            print *, 'CONNECTION_RUN       =', MP_CONNECTION_RUN, ' (Default)'
    end if

    !-------------------------------------------
    ! reduction rate of speed light in M1-closer
    !-------------------------------------------
!    MP_Crd  = 3.d-4 ! reduction rate of light speed
     MP_Crd  = 6.d-3 ! reduction rate of light speed
  !   MP_Crd  = 9.d-3
 !    MP_Crd  = 3.d-3

    !-----------EO_added-----------!
!    MP_Crd  = 6.d-3 ! reduction rate of light speed
    !-----------EO_added-----------!


    MP_Ctil = MP_Crd*cgs_c
    MP_Ctil_nD = MP_Ctil/Unit_v
    MP_PHON = 1.d49

    !-------------------------------------------
    ! parameter for chemistry 
    !-------------------------------------------
!#ifndef METAL_TRANSFER
!    if (readenv('Metallicity', MP_Metallicity)) then
!      if(get_myrank() == PRIMARY_RANK) then
!        print *, 'Metallicity       =', MP_Metallicity    ! Metallicity
!      endif
!    else
!      if(get_myrank() == PRIMARY_RANK) print *, 'metallicity is not defined'
!    endif
!#endif

    ! abundance
    MP_frac_C_solar = 0.927d-4
    MP_frac_O_solar = 3.568d-4

    MP_frac_C = MP_frac_C_solar * MP_Metallicity
    MP_frac_O = MP_frac_O_solar * MP_Metallicity 
    MP_frac_COsum = MP_frac_C + MP_frac_O

    ! weight
    MP_AC = 12.011d0
    MP_AO = 15.999d0

    ! O & C
    MP_fracC_AC = MP_frac_C * MP_AC
    MP_fracO_AO = MP_frac_O * MP_AO

    ! mean molecular weight (/amu)
    MP_mu = 1.008d0 + 4.004d0*yHe + MP_fracC_AC + MP_fracO_AO

    !-------------------------------------------
    ! Bondi_raius for Boxsize 
    !------------------------------------------- 
!     if (MP_GasType == 0) then ! HI gas ¤Î¾ì¹ç
!        yh2i=1.d-8
!        yhpi=1.d-8
!        yhmi=1.d-20
!        yh2pi=1.d-20
!        yhni=1.d0 - (yhpi+yhmi) - 2*(yh2i+yh2pi)
!        yeli=yhpi - yhmi + yh2pi
!        ycoi = 0.927d-4
!     else if (MP_GasType == 1 .or. MP_GasType == 2) then ! H2 gas ¤Î¾ì¹ç
!             yhni=1.d-8
!             yhpi=1.d-8
!             yhmi=1.d-20
!             yh2pi=1.d-20
!             yh2i= (1.d0 - (yhni+yhpi+yhmi) - 2*yh2pi)/2.d0
!             yeli=yhpi - yhmi + yh2pi
!             ycoi = 0.927d-4 
!           end if


!           sound_speed = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
!                /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))
     

!    Bondi_radius=2.d0*cgs_gc*MP_Mstar*cgs_msun / (sound_speed**2 + (5.d6)**2) 
!    Bondi_radius=2.d0*cgs_gc*MP_Mstar*cgs_msun / ((5.d6)*(5.d6)) 
! 
!    MP_Bondi_radius=Bondi_radius / Unit_l
 
!    MP_Boxsize = MP_Boxsize*MP_Bondi_radius

       


    call print_msg('end initializing modelParameters')
    call flush(6)

  end subroutine modelParameter_init
end module modelParameter
