#include "config.h"

!#define CREATE_RADSOURCE_BY_HAND

!control debug messages
!#define RS_DEBUG

!-----------------------------------------------------------------------
! module for radiation sources
!-----------------------------------------------------------------------
module radiationSource
  use grid
  use unit
  use parameter
  use mpilib
  use overBlockCoordinates  
#ifdef METAL 
  use kinzoku
#endif

  implicit none
  private

  !parameter
  integer,parameter :: MAX_RADIATION_SOUCE = 10000 !maximum number of radiation source in each rank

  ! type for information about radition sources
  type t_rs_info
     integer :: nsource                                                     ! number of sources
     integer,dimension(0:MAX_RADIATION_SOUCE-1) :: sid                      ! ID of sources (use pid for sink particles; ToDo)
     real(kind=DBL_KIND),dimension(MX:MZ,0:MAX_RADIATION_SOUCE-1) :: spos   ! positions of sources
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: lum          ! luminosities of sources
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: Trad         ! radiaiton temperature of sources
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: x_euv !x_euv*lum is emissivity of stellar EUV photons
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: x_fuv !x_fuv*lum is emissivity of stellar FUV photons
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: alpha_euv !mean cross section for EUV photons (cm^2)
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: heat_euv  !mean energy transferred to gas per one photoionizaiton event (erg)
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: hhm !nu-averaged H- dissociation rate 
#ifdef METAL 
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: lumeuv   !ratio of FUV LUMINOSITY
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: lumfuv   !ratio of FUV LUMINOSITY
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: sig_euv  !dust cross section for EUV
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: sig_fuv  !dust cross section for FUV
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: rOII     !ration of OII photoionization rate to that of H ionization
#endif
#ifdef SUB_GRID_MODEL_DIRECTLIGHT
     real(kind=DBL_KIND),dimension(0:MAX_RADIATION_SOUCE-1) :: mskr
#endif

  end type t_rs_info

  type(t_rs_info), target :: rs_info ! information about radition sources

  public :: t_rs_info, rs_GetSourceInfo
#ifdef STOCHASTIC_STELLAR_MODEL
  public :: stochastic_radiation_source 
#endif
#ifndef STOCHASTIC_STELLAR_MODEL
  public :: radiation_source
#endif
contains

!-----------------------------------------------------------------------
! subroutine to calculate radiation sources
!-----------------------------------------------------------------------
#ifndef STOCHASTIC_STELLAR_MODEL
subroutine radiation_source  
  use modelParameter
  use sinkParticle
#ifndef RADSOURCE_SC
  use primordial,only : ProstFit2, rad_others
#endif
 
  !for logfile
  use string, only : concat, num2char, CHARLEN
  use io_util, only : read_env
  
  ! real(kind=DBL_KIND),parameter :: mass_cr = 1d-2 !in M_sun
  !--------------------------- KS DEBUG --------------------------!
#ifndef MINMASS_RS
  real(kind=DBL_KIND),parameter :: mass_cr = 5d0 !in M_sun
#else
  real(kind=DBL_KIND),parameter :: mass_cr = MINMASS_RS !in M_sun
#endif

  !--------------------------- KS DEBUG --------------------------!  

  integer :: nsource, sid,n,np, nparticle,sp_Level,ns
  real(kind=DBL_KIND),dimension(MX:MZ) :: pos
  real(kind=DBL_KIND) :: lum,radius,Trad,mass,dmass
  integer,dimension(:),allocatable :: pid
  real(kind=DBL_KIND),dimension(:),allocatable :: pmass,pmdot,pmdot_disk,pt_prev,pdm_disk
#ifndef RADSOURCE_SC
  real(kind=DBL_KIND),dimension(:),allocatable :: ptcrt
  real(kind=DBL_KIND) :: tage
#endif
  real(kind=DBL_KIND),dimension(:,:),allocatable :: pr,pJ_disk
  real :: mdot
  integer :: n_rs !radiation sourceのindexはsink粒子のindexと別に用意する必要あり (KS ADDED)

  !for log file
  integer,parameter :: LUN = 11
  character(len=CHARLEN) :: file, dir
  character(len=CHARLEN),parameter :: FILENAME_LOG='logRadiationSource'
  integer,parameter :: log_skip = 10


  !----------EO_added----------!
  ! real(kind=DBL_KIND),parameter :: eta  
   real(kind=DBL_KIND),parameter :: pi  = 3.14159265   
   real(kind=DBL_KIND),parameter :: sigma_tom  = 6.65d-25  
   real(kind=DBL_KIND) :: eta_edd, lum_edd
  !----------EO_added----------!


#ifdef RADSOURCE_SC
  real(kind=DBL_KIND) :: x_euv,x_fuv,alpha_euv,heat_euv,hhm, lumeuv, lumfuv,sig_euv, sig_fuv, rOII
#else
  type(rad_others) :: rda
#endif

  !for initial call
  integer,save :: ifirst = 0

#ifdef NO_RADIATION
  rs_info%nsource = 0
  if(get_myrank() == PRIMARY_RANK) &
     print '(A)', "NO_RADIATION option is chosen  -->   rs_info%nsource = 0"
  return
#endif !NO_RADIATION

  !get sink particle info
  nparticle = sp_getNparticle()

!#ifdef CREATE_RADSOURCE_BY_HAND
  !--------------    create sink particle by hand (KS DEBUG)    -----------------!
  if (ifirst == 0) then
     if(nparticle == 0) then
        sp_Level = sp_getLevel()
        pos = (/0.d0, 0.d0, 0.d0/)
        mass = MP_Mstar                         !mass of particle in code unit
        call sp_newParticle(mass, pos, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), &
         mdot_disk=MP_Mdot/(Unit_msun/Unit_yr) , J_disk=(/0.d0, 0.d0, 1.d0/))        ! give disk mdot/J
        nparticle = sp_getNparticle()
     end if
     ifirst = 1
  end if
  !------------------------------------------------------------------------------!
!#endif

  if (nparticle > 0) then     ! sink particles exist
     ! allocate(pid(nparticle), pmass(nparticle), pmdot(nparticle), pr(MX:MZ, nparticle), pmdot_disk(nparticle), pJ_disk(MX:MZ, nparticle))
     ! call sp_sinkdata2array(np, pmass, pmdot=pmdot, pr=pr, pid=pid, pmdot_disk=pmdot_disk, pJ_disk=pJ_disk)

#ifdef RADSOURCE_SC
     allocate(pid(nparticle), pmass(nparticle), pr(MX:MZ, nparticle), pmdot_disk(nparticle), pJ_disk(MX:MZ, nparticle),&
      pt_prev(nparticle),  pdm_disk(nparticle))
     call sp_sinkdata2array(np, pmass, pr=pr, pid=pid, pmdot_disk=pmdot_disk, pJ_disk=pJ_disk, pt_prev=pt_prev, pdm_disk=pdm_disk)
#else
     allocate(pid(nparticle), pmass(nparticle), pr(MX:MZ, nparticle), pmdot_disk(nparticle), pJ_disk(MX:MZ, nparticle),&
      pt_prev(nparticle),  pdm_disk(nparticle), ptcrt(nparticle))
     call sp_sinkdata2array(np, pmass, pr=pr, pid=pid, pmdot_disk=pmdot_disk, pJ_disk=pJ_disk, pt_prev=pt_prev &
       , pdm_disk=pdm_disk, ptcrt=ptcrt)
#endif


  else
     return ! return if sink does not exist
  end if


  if(get_myrank() == PRIMARY_RANK) &
       print '(/,A)', "(RadSrc) Begin evaluating radiation sources -----------------------------------"


  !create radiation sources from sink particles
  rs_info%nsource = 0
  n_rs = 0 !radiation sourceのindexは別に定義する必要あり
  do n = 0, nparticle-1
     ns = n + 1 !sink index starts from 1 (otherwise index starts from 0)

     !------------------------!
     !        get mdot        !        
     !------------------------!
     !基本的にはdiskの降着率をmdotとする
     mdot = pmdot_disk(ns)
     !ただし、「平均降着率がまだ評価されていない」 or 「前に求めた平均降着率が非常に小さかった」場合は
     if (mdot*Unit_msun/Unit_yr < 1d-20) then 
        mdot =  pdm_disk(ns)/(Time(Lmin) - pt_prev(ns)) ! t_prevから現在までの暫定的な平均値を利用
        if(get_myrank() == PRIMARY_RANK) &
             print '(A,I6,A,3((1P1E12.4),A))', "mdot_disk is zero -> use tentative average: pid = ",pid(ns), &
             ", dt = ", (Time(Lmin)-pt_prev(ns))*Unit_yr," yr, t_prev = ",pt_prev(ns)*Unit_yr," yr, dm_disk = ", pdm_disk(ns)*Unit_msun, " M_sun"
     end if


     !assume only particles with mass > 1e-2 M_sun emit radiation
     if (pmass(ns)*Unit_msun < mass_cr) then
        if(get_myrank() == PRIMARY_RANK) &
             print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,1((1P1E12.4),A))', "pid =",pid(ns),", mass=",pmass(ns)*Unit_msun, &
             " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,ns),&
             " (in code unit)","   ==> skipping because of its small mass (<", mass_cr ," M_sun)"
        cycle 
     end if
     
     !------- WARNING WARNING  WARNING WARNING (KS DEBUG) WARNING WARNING WARNING WARNING -------!
     ! if(get_myrank() == PRIMARY_RANK) &
     !      print '(A)', "*** KS WARNING *** M and Mdot are artifically set to 200 M_sun and 1e-3 M_sun/yr, respectively"
     ! pmass(ns)= 2d2/Unit_msun
     ! pmdot_disk(ns)= 1d-3/(Unit_msun/Unit_yr)
     ! mdot= 1d-3/(Unit_msun/Unit_yr)
     !------- WARNING WARNING  WARNING WARNING (KS DEBUG) WARNING WARNING WARNING WARNING -------!     

#ifdef RADSOURCE_SC


     call rs_star_cluster(pmass(ns)*Unit_msun, lum, Trad, x_euv, alpha_euv, heat_euv, x_fuv, hhm, &
            lumeuv, lumfuv, sig_euv, sig_fuv, rOII)

     if(get_myrank() == PRIMARY_RANK) &
          print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,4((1P1E12.4),A))',  "pid =",pid(ns),", mass=",pmass(ns)*Unit_msun, &
             " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,ns),&
             " (in code unit)","   ==>   L=", lum/cgs_lsun, " L_sun, T_rad=",Trad," K, Ndot_ion=", lum*x_euv, " s^-1, Ndot_LW=", lum*x_fuv, " s^-1"
          
     ! ----------------------- write radiation source data to logfile ------------------- !
     if(get_myrank() == PRIMARY_RANK) then
        if (mod(Step(Lmin), log_skip) == 0) then
           call read_env('DIR', dir)
           file = concat(dir,FILENAME_LOG)
           open(LUN, file=file, position='APPEND')       
           write(LUN, '(I14, 1P1E17.9, I5, 1P7E17.9)') Step(Lmin), Time(Lmin)*Unit_yr, &
                pid(ns), pmass(ns)*Unit_msun, mdot*Unit_msun/Unit_yr, &
                lum/cgs_lsun, 0.e0, Trad, lum*x_euv, lum*x_fuv
           call flush(LUN)
           close(LUN)
        end if
     end if
     ! ---------------------------------------------------------------------------------- !

     !number of radiation sources
     rs_info%nsource = rs_info%nsource + 1
     if (rs_info%nsource == MAX_RADIATION_SOUCE) then
        print '(A,/,A)','(RadSrc) number of radiation sources exceeds the maximum', 'stopping...'
        stop
     end if

     !set rs_info
     rs_info%sid(n_rs)        = pid(ns)
     rs_info%spos(:,n_rs)     = pr(:,ns)
     rs_info%lum(n_rs)        = lum         ! luminosity [erg s^-1]
     rs_info%Trad(n_rs)       = Trad        ! effective temperature [K]
     rs_info%x_euv(n_rs)      = x_euv       ! the emissivity of euv photons
     rs_info%x_fuv(n_rs)      = x_fuv       !                   fuv
     rs_info%alpha_euv(n_rs)  = alpha_euv   ! mean cross section for EUV photons [cm^-2]
     rs_info%heat_euv(n_rs)   = heat_euv    ! mean heating rate [erg]
     rs_info%hhm(n_rs)        = hhm         ! H- dissociation rate at stellar surface [s^-1]
  #ifdef METAL
     rs_info%lumeuv(n_rs)     = lumeuv      ! ratio of euv photons
     rs_info%lumfuv(n_rs)     = lumfuv      ! ratio of fuv photons
     rs_info%sig_euv(n_rs)    = sig_euv     ! dust cross section for euv photons [cm^-2]
     rs_info%sig_fuv(n_rs)    = sig_fuv     ! dust cross section for fuv photons [cm^-2]
     rs_info%rOII(n_rs)       = rOII        ! fraction of OII photoionization rate to that of H
  #endif



     
     !increment index for next radiaiton source
     n_rs = n_rs + 1
    
#else

    ! ===============================
    !         Editting here 
    ! ===============================


     !(mass, mdot) -> (lum, Trad) 
     !call ProstFit(pmass(ns)*Unit_msun, mdot*Unit_msun/Unit_yr, radius, lum)     

     !----------EO_added----------!
    ! eta_edd = 100.0d0 * 5.33d-4  !when dusty-gas
    ! eta_edd = 0.6              !when only gas(no dust)
    ! lum_edd = 4.d0 * pi * cgs_c * cgs_gc * pmass(ns) * Unit_m  * cgs_mp / sigma_tom     
    ! eta_edd = 1.d0 / (1.d0 + 7.1d1)
    ! lum_edd = 4.d0 * pi * cgs_c * cgs_gc * 1.d4 * cgs_msun * cgs_mp / sigma_tom     
    ! lum     = eta_edd * lum_edd    
     !----------EO_added----------!


  !----------EO_removed----------!
  #ifdef REDUCE_SMRATE
     call ProstFit2(pmass(ns)*Unit_msun*REDUCE_SMRATE,mdot*Unit_msun/Unit_yr,(Time(Lmin)-ptcrt(ns))*Unit_yr,radius,lum,Trad,rda)

  #else
     call ProstFit2(pmass(ns)*Unit_msun              ,mdot*Unit_msun/Unit_yr,(Time(Lmin)-ptcrt(ns))*Unit_yr,radius,lum,Trad,rda)
  #endif

    ! lum = lum * cgs_lsun                                           !luminosity in [erg s^-1]
  !----------EO_removed----------!

      !----------EO_added----------!
    ! eta_edd = 100.0d0 * 5.33d-4  !when dusty-gas
    ! eta_edd = 0.6              !when only gas(no dust)
    ! lum_edd = 4.d0 * pi * cgs_c * cgs_gc * pmass(ns) * Unit_m  * cgs_mp / sigma_tom     

     eta_edd = 0.03
    ! eta_edd = 1.d0 / (1.d0 + 7.1d1)
     lum_edd = 4.d0 * pi * cgs_c * cgs_gc * 1.d4 * cgs_msun * cgs_mp / sigma_tom     
     lum     = eta_edd * lum_edd    
     !----------EO_added----------!

  


     if(get_myrank() == PRIMARY_RANK) &
          print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,4((1P1E12.4),A))',  "pid =",pid(ns),", mass=",pmass(ns)*Unit_msun, &
             " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,ns),&
             " (in code unit)","   ==>   L=", lum/cgs_lsun, " L_sun, T_rad=",Trad," K, Ndot_ion=" &
             , lum*rda%xeuv, " s^-1, Ndot_LW=", lum*rda%xfuv, " s^-1"
          
     ! ----------------------- write radiation source data to logfile ------------------- !
     if(get_myrank() == PRIMARY_RANK) then
        if (mod(Step(Lmin), log_skip) == 0) then
           call read_env('DIR', dir)
           file = concat(dir,FILENAME_LOG)
           open(LUN, file=file, position='APPEND')       
           write(LUN, '(I14, 1P1E17.9, I5, 1P9E17.9)') Step(Lmin), Time(Lmin)*Unit_yr, &
                pid(ns), pmass(ns)*Unit_msun, mdot*Unit_msun/Unit_yr, &
                lum/cgs_lsun, radius, Trad, lum*rda%xeuv, lum*rda%xfuv, rda%lumeuv, rda%lumfuv
           call flush(LUN)
           close(LUN)
        end if
     end if
     ! ---------------------------------------------------------------------------------- !

     !number of radiation sources
     rs_info%nsource = rs_info%nsource + 1
     if (rs_info%nsource == MAX_RADIATION_SOUCE) then
        print '(A,/,A)','(RadSrc) number of radiation sources exceeds the maximum', 'stopping...'
        stop
     end if

     !set rs_info
     !----------EO_removed----------!
     rs_info%sid(n_rs)        = pid(ns)
     rs_info%spos(:,n_rs)     = pr(:,ns)
     rs_info%lum(n_rs)        = lum         ! luminosity [erg s^-1]
     rs_info%Trad(n_rs)       = Trad        ! effective temperature [K]
     rs_info%x_euv(n_rs)      = rda%xeuv       ! the emissivity of euv photons
     rs_info%x_fuv(n_rs)      = rda%xfuv       !                   fuv
     rs_info%alpha_euv(n_rs)  = rda%alpha_euv   ! mean cross section for EUV photons [cm^-2]
     rs_info%heat_euv(n_rs)   = rda%heat_euv    ! mean heating rate [erg]
     rs_info%hhm(n_rs)        = rda%hhm         ! H- dissociation rate at stellar surface [s^-1]
  #ifdef METAL
     rs_info%lumeuv(n_rs)     = rda%lumeuv      ! ratio of euv photons
     rs_info%lumfuv(n_rs)     = rda%lumfuv      ! ratio of fuv photons
     rs_info%sig_euv(n_rs)    = rda%sig_euv     ! dust cross section for euv photons [cm^-2]
     rs_info%sig_fuv(n_rs)    = rda%sig_fuv     ! dust cross section for fuv photons [cm^-2]
     rs_info%rOII(n_rs)       = rda%rOII        ! fraction of OII photoionization rate to that of H
  #endif
     !----------EO_removed----------!


     !----------EO_added----------!
!     rs_info%sid(n_rs)        = pid(ns)
!     rs_info%spos(:,n_rs)     = pr(:,ns)
!     rs_info%lum(n_rs)        = lum         ! luminosity [erg s^-1]
!     rs_info%Trad(n_rs)       = 1.d0        ! effective temperature [K]
!     rs_info%x_euv(n_rs)      = 1.d0 / (1.36d1 * cgs_ev)       ! the emissivity of euv photons
!     rs_info%x_fuv(n_rs)      = 0.d0       !                   fuv
!     rs_info%alpha_euv(n_rs)  = 6.3d0 * 1.d-18   ! mean cross section for EUV photons [cm^-2]
!     rs_info%alpha_euv(n_rs)  = 1.d-21   ! mean cross section for EUV photons [cm^-2]


!     rs_info%heat_euv(n_rs)   = 0.d0    ! mean heating rate [erg]
!     rs_info%hhm(n_rs)        = 0.d0         ! H- dissociation rate at stellar surface [s^-1]
!  #ifdef METAL
!     rs_info%lumeuv(n_rs)     = 1.d0      ! ratio of euv photons
!     rs_info%lumfuv(n_rs)     = 0.d0      ! ratio of fuv photons
!     rs_info%sig_euv(n_rs)    = 1.d-21     ! dust cross section for euv photons [cm^-2]
!     rs_info%sig_fuv(n_rs)    = 1.d-21     ! dust cross section for fuv photons [cm^-2]
!     rs_info%rOII(n_rs)       = 0.d0        ! fraction of OII photoionization rate to that of H
!  #endif
     !----------EO_added-----------!     

     !increment index for next radiaiton source
     n_rs = n_rs + 1

#endif

  end do
  if(get_myrank() == PRIMARY_RANK) &
       print '(A)', "(RadSrc) End evaluating radiation sources -----------------------------------"  

  ! deallocate ----------
  deallocate(pid, pmass, pr, pmdot_disk, pJ_disk, pt_prev, pdm_disk)
  ! ---------------------

  end subroutine radiation_source
#endif !STOCHASTIC_STELLAR_MODEL

!-----------------------------------------------------------------------
! subroutine to calculate radiation sources in stochastic
!-----------------------------------------------------------------------
#ifdef STOCHASTIC_STELLAR_MODEL
subroutine stochastic_radiation_source  
  use modelParameter
  use sinkParticle
  use stochastic_star
 
  !for logfile
  use string, only : concat, num2char, CHARLEN
  use io_util, only : read_env
  
  !--------------------------- KS DEBUG --------------------------!  
  integer :: np, nparticle,sp_Level,ns,i, rank
  real(kind=DBL_KIND),dimension(MX:MZ) :: pos
  real(kind=DBL_KIND) :: lum,radius,Trad,mass
  integer,dimension(:),allocatable :: pid
  real(kind=DBL_KIND),dimension(:),allocatable :: pmass,pmdot_disk,pt_prev,pdm_disk
  real(kind=DBL_KIND),dimension(:,:),allocatable :: pr,pJ_disk
  real(kind=DBL_KIND) :: mdot
  integer :: n_rs !radiation sourceのindexはsink粒子のindexと別に用意する必要あり (KS ADDED)
  integer :: nplc, np_local, lnp, id_sink, my_rank, il, il_int, npy
  integer :: recvcounts(0:NPE-1), dispals(0:NPE-1)
  integer :: recvcounts_int(0:NPE-1), dispals_int(0:NPE-1)
  integer, parameter :: ninf = 11
  real(kind=DBL_KIND) :: stime
  real(kind=DBL_KIND), allocatable, dimension(:) :: sbuf, rbuf
  integer, allocatable, dimension(:) :: sbuf_int, rbuf_int
  type(st_rs_info) :: rsinf_stc
  type(st_rs_info), allocatable, dimension(:) :: rsinf_stc_all
  integer :: pid_max, tsink, local_pid, np_tot
  integer, allocatable, dimension(:) :: pid2ns

  !for log file
  integer,parameter :: LUN = 11
  character(len=CHARLEN) :: file, dir
  character(len=CHARLEN),parameter :: FILENAME_LOG='logRadiationSource'
  integer,parameter :: log_skip = 10
  real(kind=DBL_KIND) :: x_euv,x_fuv,alpha_euv,heat_euv, lumeuv, lumfuv,sig_euv, sig_fuv, rOII

  !for initial call
  integer,save :: ifirst = 0

  ! get my_rank ----------
  my_rank = get_myrank()
  ! ----------------------

#ifdef NO_RADIATION
  rs_info%nsource = 0
  if(my_rank == PRIMARY_RANK) &
     print '(A)', "NO_RADIATION option is chosen  -->   rs_info%nsource = 0"
  return
#endif !NO_RADIATION

  !get sink particle info
  nparticle = sp_getNparticle()


!#ifdef CREATE_RADSOURCE_BY_HAND
  !--------------    create sink particle by hand (KS DEBUG)    -----------------!
  if (ifirst == 0) then
     if(nparticle == 0) then
        sp_Level = sp_getLevel()
        pos = (/0.d0, 0.d0, 0.d0/)
        mass = MP_Mstar*cgs_msun/Unit_m                         !mass of particle in code unit
        call sp_newParticle(mass, pos, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), &
         mdot_disk=MP_Mdot/(Unit_msun/Unit_yr) , J_disk=(/0.d0, 0.d0, 1.d0/))        ! give disk mdot/J
        nparticle = sp_getNparticle()

        if(my_rank == PRIMARY_RANK) then
          call set_radsc_stchast(my_rank, 1, 1, MP_Mstar, 0.d0, -1.d6, rsinf_stc)  
        endif
     end if
     ifirst = 1
  end if
  !------------------------------------------------------------------------------!
!#endif

  if (nparticle > 0) then     ! sink particles exist
     ! allocate(pid(nparticle), pmass(nparticle), pmdot(nparticle), pr(MX:MZ, nparticle), pmdot_disk(nparticle), pJ_disk(MX:MZ, nparticle))
     allocate(pid(nparticle), pmass(nparticle), pr(MX:MZ, nparticle), pmdot_disk(nparticle), pJ_disk(MX:MZ, nparticle),&
      pt_prev(nparticle),  pdm_disk(nparticle))
     call sp_sinkdata2array(np_tot, pmass, pr=pr, pid=pid, pmdot_disk=pmdot_disk, pJ_disk=pJ_disk, pt_prev=pt_prev, pdm_disk=pdm_disk)

     ! ---------------------------------------------
     pid_max = 0
     do np = 1, np_tot 
        pid_max = max(pid(np), pid_max)
     enddo
     allocate(pid2ns(0:pid_max))
     
     pid2ns(:) = -1
  
     do np = 1, np_tot
        pid2ns(pid(np)) = np
     enddo

     tsink = pid_max+1
     ! ---------------------------------------------

  else
     return ! return if sink does not exist
  end if


  if(my_rank == PRIMARY_RANK) &
       print '(/,A)', "(RadSrc) Begin evaluating radiation sources -----------------------------------"


  ! set labels for read stellar data
  nplc     = int(tsink/NPE)
  npy      = tsink-nplc*NPE
  np_local = nplc
  if(my_rank < npy) then
    np_local = np_local + 1
  endif

  ! ------------------------------------------
  do i=0, NPE-1
    recvcounts(i)     = nplc*ninf
    recvcounts_int(i) = nplc
  enddo

  do i=0, npy-1
      recvcounts(i)     = recvcounts(i) + ninf
      recvcounts_int(i) = recvcounts_int(i) + 1
  enddo
  dispals(0)     = 0
  dispals_int(0) = 0
  do rank = 1, NPE-1
    dispals(rank)       = dispals(rank-1)     + recvcounts(rank-1)
    dispals_int(rank)   = dispals_int(rank-1) + recvcounts_int(rank-1) 
  enddo  
  ! ------------------------------------------

  ! allocate 
  allocate(sbuf(np_local*ninf), rbuf(tsink*ninf))
  allocate(sbuf_int(np_local), rbuf_int(tsink))


  ! -----------------------------!
  !  stochastic stellar evolv
  ! -----------------------------!
  do np =0, np_local-1

    lnp       = np + 1
    local_pid = my_rank + NPE*(np)
    id_sink   = pid2ns(local_pid) 

    if ( id_sink > 0) then

      mdot    = pmdot_disk(id_sink)*Unit_msun/Unit_yr ! [Msun/yr] 
      stime   = Time(Lmin)*Unit_yr ! [yr] elapsed time in simulations

      call set_radsc_stchast(my_rank, local_pid, lnp, pmass(id_sink)*Unit_msun , mdot, stime, rsinf_stc)  

      il = np*ninf
      sbuf(il+1)      = rsinf_stc%lum 
      sbuf(il+2)      = rsinf_stc%Teff
      sbuf(il+3)      = rsinf_stc%xeuv
      sbuf(il+4)      = rsinf_stc%xfuv
      sbuf(il+5)      = rsinf_stc%alpha_euv
      sbuf(il+6)      = rsinf_stc%heat_euv
      sbuf(il+7)      = rsinf_stc%lumeuv
      sbuf(il+8)      = rsinf_stc%lumfuv
      sbuf(il+9)      = rsinf_stc%sigd_euv
      sbuf(il+10)     = rsinf_stc%sigd_fuv
      sbuf(il+11)     = rsinf_stc%rOII
      sbuf_int(np+1)  = rsinf_stc%nstar
    else
      il = np*ninf
      sbuf(il+1)      = 0.d0 
      sbuf(il+2)      = 0.d0 
      sbuf(il+3)      = 0.d0 
      sbuf(il+4)      = 0.d0 
      sbuf(il+5)      = 0.d0 
      sbuf(il+6)      = 0.d0 
      sbuf(il+7)      = 0.d0 
      sbuf(il+8)      = 0.d0 
      sbuf(il+9)      = 0.d0 
      sbuf(il+10)     = 0.d0 
      sbuf(il+11)     = 0.d0 
      sbuf_int(np+1)  = 0
    endif

  enddo ! np =0, np_local-1  

  ! -----------------------------!
  !       MPI transport          !
  ! -----------------------------!
  call mpi_allgatherv( sbuf, np_local*ninf, MPI_DOUBLE_PRECISION, rbuf, recvcounts, dispals, MPI_DOUBLE_PRECISION &
    , MPI_COMM_WORLD, ierr)

  call mpi_allgatherv( sbuf_int, np_local, MPI_INTEGER, rbuf_int, recvcounts_int, dispals_int, MPI_INTEGER &
    , MPI_COMM_WORLD, ierr)



  ! deallocate -------------
  deallocate(sbuf, sbuf_int)
  allocate(rsinf_stc_all(np_tot))
  ! ------------------------

  il     = 0
  il_int = 1
  do i=0, NPE-1
    do np=0, nplc-1

      local_pid= i + NPE*(np)
      id_sink = pid2ns(local_pid) 

      if ( id_sink > 0) then
        rsinf_stc_all(id_sink)%lum       = rbuf(il+1)  
        rsinf_stc_all(id_sink)%Teff      = rbuf(il+2) 
        rsinf_stc_all(id_sink)%xeuv      = rbuf(il+3) 
        rsinf_stc_all(id_sink)%xfuv      = rbuf(il+4) 
        rsinf_stc_all(id_sink)%alpha_euv = rbuf(il+5) 
        rsinf_stc_all(id_sink)%heat_euv  = rbuf(il+6) 
        rsinf_stc_all(id_sink)%lumeuv    = rbuf(il+7) 
        rsinf_stc_all(id_sink)%lumfuv    = rbuf(il+8) 
        rsinf_stc_all(id_sink)%sigd_euv  = rbuf(il+9) 
        rsinf_stc_all(id_sink)%sigd_fuv  = rbuf(il+10)
        rsinf_stc_all(id_sink)%rOII      = rbuf(il+11)
        rsinf_stc_all(id_sink)%nstar     = rbuf_int(il_int)
      endif

      ! ----------------
      il     = il + ninf
      il_int = il_int + 1
      ! ----------------

    enddo

    if(i < npy) then

      local_pid= i + NPE*(nplc)
      id_sink = pid2ns(local_pid) 

      if ( id_sink > 0) then
        rsinf_stc_all(id_sink)%lum       = rbuf(il+1)  
        rsinf_stc_all(id_sink)%Teff      = rbuf(il+2) 
        rsinf_stc_all(id_sink)%xeuv      = rbuf(il+3) 
        rsinf_stc_all(id_sink)%xfuv      = rbuf(il+4) 
        rsinf_stc_all(id_sink)%alpha_euv = rbuf(il+5) 
        rsinf_stc_all(id_sink)%heat_euv  = rbuf(il+6) 
        rsinf_stc_all(id_sink)%lumeuv    = rbuf(il+7) 
        rsinf_stc_all(id_sink)%lumfuv    = rbuf(il+8) 
        rsinf_stc_all(id_sink)%sigd_euv  = rbuf(il+9) 
        rsinf_stc_all(id_sink)%sigd_fuv  = rbuf(il+10)
        rsinf_stc_all(id_sink)%rOII      = rbuf(il+11)
        rsinf_stc_all(id_sink)%nstar     = rbuf_int(il_int)
      endif

      ! ----------------
      il     = il + ninf
      il_int = il_int + 1
      ! ----------------
    endif

  enddo


  ! -----------------------------!
  !     radiation  source        !
  ! -----------------------------!
  n_rs = 0
  rs_info%nsource = 0
  do id_sink =1, np_tot

    if(rsinf_stc_all(id_sink)%nstar > 1) then ! 光源1つ以上ある場合だけ記述

        rs_info%nsource = rs_info%nsource + 1

        ! check maximum !
        if (rs_info%nsource > MAX_RADIATION_SOUCE) then
           print '(A,/,A)','(RadSrc) number of radiation sources exceeds the maximum', 'stopping...'
           stop
        end if

        !set rs_info
        rs_info%sid(n_rs)        = pid(id_sink)
        rs_info%spos(:,n_rs)     = pr(:,id_sink)
        rs_info%lum(n_rs)        = rsinf_stc_all(id_sink)%lum       ! luminosity [erg s^-1]
        rs_info%Trad(n_rs)       = rsinf_stc_all(id_sink)%Teff      ! effective temperature [K]
        rs_info%x_euv(n_rs)      = rsinf_stc_all(id_sink)%xeuv      ! the emissivity of euv photons
        rs_info%x_fuv(n_rs)      = rsinf_stc_all(id_sink)%xfuv      !                   fuv
        rs_info%alpha_euv(n_rs)  = rsinf_stc_all(id_sink)%alpha_euv ! mean cross section for EUV photons [cm^-2]
        rs_info%heat_euv(n_rs)   = rsinf_stc_all(id_sink)%heat_euv  ! mean heating rate [erg]
        rs_info%hhm(n_rs)        = 0.d0             ! H- dissociation rate at stellar surface [s^-1]
#ifdef METAL
        rs_info%lumeuv(n_rs)     = rsinf_stc_all(id_sink)%lumeuv    ! ratio of euv photons
        rs_info%lumfuv(n_rs)     = rsinf_stc_all(id_sink)%lumfuv    ! ratio of fuv photons
        rs_info%sig_euv(n_rs)    = rsinf_stc_all(id_sink)%sigd_euv  ! dust cross section for euv photons [cm^-2]
        rs_info%sig_fuv(n_rs)    = rsinf_stc_all(id_sink)%sigd_fuv  ! dust cross section for fuv photons [cm^-2]
        rs_info%rOII(n_rs)       = rsinf_stc_all(id_sink)%rOII      ! fraction of OII photoionization rate to that of H
#endif
        n_rs = n_rs + 1

    endif
  enddo


  ! ---------------------------------
  !           logoutput 
  ! ---------------------------------
  if(my_rank == PRIMARY_RANK) then

    n_rs = 0
    do id_sink =1, np_tot

      if(rsinf_stc_all(id_sink)%nstar < 2) then
        print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,1(I6,A))', "pid =",pid(id_sink),", mass=",pmass(id_sink)*Unit_msun, &
        " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,id_sink),&
        " (in code unit)","   ==> skipping because of stellar number is less than (< 2", rsinf_stc_all(id_sink)%nstar ," ko)"
        cycle
      endif

      mdot = pmdot_disk(id_sink)
      lum  = rs_info%lum(n_rs) 
      Trad = rs_info%Trad(n_rs)  
      x_euv= rs_info%x_euv(n_rs)
      x_fuv= rs_info%x_fuv(n_rs) 
      print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,4((1P1E12.4),A))',  "pid =",pid(id_sink),", mass=",pmass(id_sink)*Unit_msun, &
         " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,id_sink),&
         " (in code unit)","   ==>   L=", lum/cgs_lsun, " L_sun, T_rad=",Trad," K, Ndot_ion=", lum*x_euv, " s^-1, Ndot_LW=", lum*x_fuv, " s^-1"

      ! -----------------------------------------------------------------------------
      if (mod(Step(Lmin), log_skip) == 0) then
        call read_env('DIR', dir)
        file = concat(dir,FILENAME_LOG)
        open(LUN, file=file, position='APPEND')       
        write(LUN, '(I14, 1P1E17.9, I5, 1P12E17.9)') Step(Lmin), Time(Lmin)*Unit_yr, &
             pid(id_sink), pmass(id_sink)*Unit_msun, mdot*Unit_msun/Unit_yr, &
             lum/cgs_lsun, 0.e0, Trad, lum*x_euv, lum*x_fuv, &
             rs_info%alpha_euv(n_rs), rs_info%heat_euv(n_rs), rs_info%sig_euv(n_rs), rs_info%sig_fuv(n_rs), rs_info%rOII(n_rs)
        call flush(LUN)
        close(LUN)
      endif
      ! -----------------------------------------------------------------------------

      n_rs = n_rs + 1
    enddo

    print '(A)', "(RadSrc) End evaluating radiation sources -----------------------------------"  

  endif


  ! deallocate ----------
  deallocate(pid2ns)
  deallocate(rsinf_stc_all)
  deallocate(rbuf, rbuf_int)
  deallocate(pid, pmass, pr, pmdot_disk, pJ_disk, pt_prev, pdm_disk)
  ! ---------------------

  end subroutine stochastic_radiation_source  
#endif !STOCHASTIC_STELLAR_MODEL



  !-------------------------------------------------------------------------
  ! get radiation source info
  ! OUTPUT:
  !  p_rs_info ...... information about all radiation sources (pointer)
  !-------------------------------------------------------------------------
  subroutine rs_GetSourceInfo(p_rs_info)
    type(t_rs_info),pointer,intent(OUT) :: p_rs_info
    p_rs_info => rs_info
  end subroutine rs_GetSourceInfo
  !-------------------------------------------------------------------------
  ! search isrc for given sid
  !-------------------------------------------------------------------------
  subroutine rs_sid2isrc(sid,isrc)
    integer,intent(IN)::sid
    integer,intent(OUT)::isrc
    integer :: nsource,i
    nsource = rs_info%nsource
    do i = 0, nsource-1
       if (rs_info%sid(i) == sid) then
          isrc = i
          exit
       end if
    end do
  end subroutine rs_sid2isrc

  !-------------------------------------------------------------------------
  ! get R_* and L_* from M_* and Mdot using a fitting of Prost data
  !-------------------------------------------------------------------------
  subroutine ProstFit(Mass, Mdot, Rad, Lum)
    real(kind=DBL_KIND),intent(IN) :: Mass !mass in [M_sun]
    real(kind=DBL_KIND),intent(IN) :: Mdot !mdot in [M_sun/yr]
    real(kind=DBL_KIND),intent(Out) :: Rad !radius in [R_sun]
    real(kind=DBL_KIND),intent(Out) :: Lum !luminosity in [L_sun]

    !data file
    integer,parameter :: nmd=11  !number of mdot in data
    integer,parameter :: npt=250 !number of mass in data
    character(100) :: fname(0:nmd-1) = (/"1em5fit.dat", "1em4fit.dat", "3em4fit.dat", "1em3fit.dat", &
         "3em3fit.dat", "6em3fit.dat", "1em2fit.dat", "3em2fit.dat", &
         "6em2fit.dat", "1em1fit.dat", "3em1fit.dat"/)    
    character(100) :: path2fitdat = "./ProstFit/"
    integer,parameter :: FH = 11
    character(len=100) :: ffn
    integer :: err

    !data container
    real(kind=DBL_KIND),save :: xm(0:npt-1)                !mass in [M_sun]
    real(kind=DBL_KIND),save :: xr(0:npt-1,0:nmd-1)        !radius in [R_sun]
    real(kind=DBL_KIND),save :: xls(0:npt-1,0:nmd-1)       !luminosity in [L_sun]
    real(kind=DBL_KIND),parameter ::  xmd(0:nmd-1) = (/1.d-5, 1.d-4, 3.d-4, 1.d-3, 3.d-3, 6.d-3,& 
         1.d-2, 3.d-2, 6.d-2, 1.d-1, 3.d-1/)          !mdot in [M_sun/yr]

    !other variables
    real(kind=DBL_KIND) :: xxm
    integer,save :: ifirst=0
    integer :: i, id

    !initialization (read data)
    if (ifirst == 0) then
       do id = 0, nmd-1
          ffn= trim(path2fitdat)//trim(fname(id))
          if(get_myrank() == PRIMARY_RANK) print *, "PROST_FIT: reading ", ffn
          open(FH, file=ffn,status='old',iostat=err)
          if (err > 0) then
             if(get_myrank() == PRIMARY_RANK) print '(A,/,A)', "PROST_FIT: file not found","stopping..."
             stop
          end if
          do i = 0, npt-1
             read(FH,*,iostat=err) xxm, xr(i,id), xls(i,id)
             if (id == 0) xm(i) = 1d1**xxm;
             if (err > 0) then
                print *, "PROST_FIT: **WARNING** data size might be inconsistent"
                exit !ファイルの最後に達したら終了
             end if
          end do
          close(FH)
       end do
       ifirst = 1
    end if

    !output warning
    if(get_myrank() == PRIMARY_RANK) then
       if (Mass < xm(0)) &
            print '(2(A, (1P1E12.4)))', "** ProstFIT: Wargning **  M=",Mass, "is set to Mmin=",xm(0)
       if (Mass > xm(npt-1)) &
            print '(2(A, (1P1E12.4)))', "** ProstFIT: Wargning **  M=",Mass, "is set to Mmax=",xm(npt-1)
       if (Mdot < xmd(0)) &
            print '(2(A, (1P1E12.4)))', "** ProstFIT: Wargning **  Mdot=",Mdot, "is set to Mdotmin=",xmd(0)
       if (Mdot > xmd(nmd-1)) &
         print '(2(A, (1P1E12.4)))', "** ProstFIT: Wargning **  M=",Mdot, "is set to Mdotmax=",xmd(nmd-1)
    end if

    call lin2D(xm,npt,xmd,nmd,xr,Mass,Mdot,Rad)
    call lin2D(xm,npt,xmd,nmd,xls,Mass,Mdot,Lum)    
    !   print *, "KS DEBUG: Mass,Mdot,Rad,Lum = ",Mass,Mdot,Rad,Lum
  end subroutine ProstFit

!  !-------------------------------------------------------------------------
!  ! get R_* and L_* from M_* and Mdot using a fitting of Prost data
!  !-------------------------------------------------------------------------
!#ifndef RADSOURCE_SC
!  subroutine ProstFit2(Mass, Mdot, Rad, Lum, Trad, xeuv, xfuv, alpha_euv, heat_euv, hhm, &
!      lumeuv, lumfuv, sig_euv, sig_fuv, rOII)
!    real(kind=DBL_KIND),intent(IN) :: Mass        ! mass in [M_sun]
!    real(kind=DBL_KIND),intent(IN) :: Mdot        ! mdot in [M_sun/yr]
!    real(kind=DBL_KIND),intent(Out) :: Rad        ! radius in [R_sun]
!    real(kind=DBL_KIND),intent(Out) :: Lum        ! luminosity in [L_sun]
!    real(kind=DBL_KIND),intent(Out) :: Trad       ! Effective temeprature [K]
!    real(kind=DBL_KIND),intent(Out) :: xeuv       ! the emissivity of stellar EUV photons xeuv * lum 
!                                                  !            = emissivity of EUV photons [s^-1]
!    real(kind=DBL_KIND),intent(Out) :: xfuv       ! the emissivity of stellar FUV photons
!    real(kind=DBL_KIND),intent(Out) :: alpha_euv  ! mean cross section for EUV [cm^-2]
!    real(kind=DBL_KIND),intent(Out) :: heat_euv   ! mean heating rate per one ionization [erg per one ionization]
!    real(kind=DBL_KIND),intent(Out) :: hhm        ! nu averaged H^- cross section at stellar surface [s^-1] 
!                                                  ! you need to multiply (r/R*) in code
!    real(kind=DBL_KIND),intent(Out) :: lumeuv     ! the raio of EUV luminosity
!    real(kind=DBL_KIND),intent(Out) :: lumfuv     ! the raio of EUV luminosity
!    real(kind=DBL_KIND),intent(Out) :: sig_euv    ! dust cross section for EUV at Z=Zsun [cm^2]
!    real(kind=DBL_KIND),intent(Out) :: sig_fuv    ! dust cross section for FUV at Z=Zsun [cm^2]
!    real(kind=DBL_KIND),intent(Out) :: rOII       ! ratio of photodisociation rate of OII and HI
!
!    !data file
!    integer,parameter :: nmd=11  !number of mdot in data
!    integer,parameter :: npt=250 !number of mass in data
!    character(100) :: fname(0:nmd-1) = (/"0e00", "1e-5", "1e-4", "3e-4", "1e-3", "3e-3" &
!      , "6e-3", "1e-2", "3e-2", "6e-2", "1e-1"/)
!    character(100) :: path2fitdat = "./ProstFit/"//trim(FOL_RADSOURCE)//"/"
!    integer,parameter :: FH = 11
!    character(len=100) :: ffn 
!    integer :: err 
!
!        !data container
!    real(kind=DBL_KIND),save :: xm(0:npt-1)                !mass in [M_sun]
!    real(kind=DBL_KIND),save :: xr(0:npt-1,0:nmd-1)        !radius in [R_sun]
!    real(kind=DBL_KIND),save :: xls(0:npt-1,0:nmd-1)       !luminosity in [L_sun]
!    real(kind=DBL_KIND),save :: xtrad(0:npt-1,0:nmd-1)     !Trad [K]
!    real(kind=DBL_KIND),save :: xxeuv(0:npt-1,0:nmd-1)     !xeuv
!    real(kind=DBL_KIND),save :: xxfuv(0:npt-1,0:nmd-1)     !xfuv
!    real(kind=DBL_KIND),save :: xaleuv(0:npt-1,0:nmd-1)    !alpha_euv
!    real(kind=DBL_KIND),save :: xheat(0:npt-1,0:nmd-1)     !heat_euv
!    real(kind=DBL_KIND),save :: xhhm(0:npt-1,0:nmd-1)      !hhm
!    real(kind=DBL_KIND),save :: xlumeuv(0:npt-1,0:nmd-1)   !lumeuv
!    real(kind=DBL_KIND),save :: xlumfuv(0:npt-1,0:nmd-1)   !lumfuv
!    real(kind=DBL_KIND),save :: xsigeuv(0:npt-1,0:nmd-1)   !sigeuv
!    real(kind=DBL_KIND),save :: xsigfuv(0:npt-1,0:nmd-1)   !sigfuv
!    real(kind=DBL_KIND),save :: xrOII(0:npt-1,0:nmd-1)     !rOII
!    real(kind=DBL_KIND),parameter ::  xmd(0:nmd-1) = (/0.d0, 1.d-5, 1.d-4, 3.d-4, 1.d-3, 3.d-3, 6.d-3,&
!         1.d-2, 3.d-2, 6.d-2, 1.d-1/)          !mdot in [M_sun/yr]
!
!    !other variables
!    integer,save :: ifirst=0
!    integer :: i, j, id, ix, iy
!    real(kind=DBL_KIND) :: x_in, y_in, x1, x2, y1, y2, z_y1, z_y2
!    real(kind=DBL_KIND) :: tap1, tap2, tap3, tap4
!
!    !initialization (read data)
!    if (ifirst == 0) then
!       do id = 0, nmd-1
!          ffn= trim(path2fitdat)//trim(fname(id))
!          if(get_myrank() == PRIMARY_RANK) print *, "PROST_FIT: reading ", ffn
!          open(FH, file=ffn,status='old',iostat=err)
!          if (err > 0) then
!             print '(A,/,A)', "PROST_FIT: file not found","stopping..."
!             stop
!          end if
!          do i = 0, npt-1
!             read(FH,*,iostat=err) xm(i), xr(i,id), xls(i,id), xtrad(i,id), xxeuv(i,id), xxfuv(i,id), xaleuv(i,id)&
!               , xheat(i,id), xhhm(i,id), xlumeuv(i,id), xlumfuv(i,id), xsigeuv(i,id), xsigfuv(i,id), xrOII(i,id)
!             if (err > 0) then
!                print *, "PROST_FIT: **WARNING** data size might be inconsistent"
!                exit !真真真真真真真
!             end if
!          end do
!          close(FH)
!       end do
!       ifirst = 1
!    end if
!
!
!    ! linear interpolation
!    x_in    = dlog10(Mass)
!    x_in    = min(max(x_in, xm(0)),  xm(npt-1))
!    y_in    = min(max(Mdot,xmd(0)), xmd(nmd-1))
!
!    ! get index of m & mdot
!    ix  = int((x_in+2.0)/0.02d0)
!    ix  = max(min(ix, npt-2),0)
!    iy = nmd-2
!    do j = 0, nmd-2
!      if (y_in <= xmd(j+1)) then
!        iy = j
!        exit
!      endif
!    enddo
!
!
!    !interpolation
!    x1 = xm(ix)
!    x2 = xm(ix+1)
!    y1 = xmd(iy)
!    y2 = xmd(iy+1)
!
!    tap1 = (x2-x_in)/(x2-x1)
!    tap2 = (x_in-x1)/(x2-x1)
!    tap3 = (y2-y_in)/(y2-y1)
!    tap4 = (y_in-y1)/(y2-y1)
!
!#define INTERPOLATION_LIN(zval, zarr) \
!    z_y1 = tap1*zarr(ix,iy)  + tap2*zarr(ix+1,iy)    ;\
!    z_y2 = tap1*zarr(ix,iy+1)+ tap2*zarr(ix+1,iy+1)  ;\
!    zval = tap3*z_y1 + tap4*z_y2
!
!    INTERPOLATION_LIN( Rad , xr   )
!    INTERPOLATION_LIN( Lum , xls  )
!    INTERPOLATION_LIN( Trad, xtrad)
!    INTERPOLATION_LIN( xeuv, xxeuv)
!    INTERPOLATION_LIN( xfuv, xxfuv)
!    INTERPOLATION_LIN( alpha_euv, xaleuv)
!    INTERPOLATION_LIN( heat_euv, xheat)
!    INTERPOLATION_LIN( hhm, xhhm)
!    INTERPOLATION_LIN( lumeuv, xlumeuv)
!    INTERPOLATION_LIN( lumfuv, xlumfuv)
!    INTERPOLATION_LIN( sig_euv, xsigeuv)
!    INTERPOLATION_LIN( sig_fuv, xsigfuv)
!    INTERPOLATION_LIN( rOII, xrOII)
!
!    !print *, x_in, x1, x2, ix, y_in, iy, xeuv, xxeuv(ix, iy)!, x_in -x1, x2 - x_in
!      !(x2-x_in)/(x2-x1)*xxeuv(ix,iy)+(x_in-x1)/(x2-x1)*xxeuv(ix+1,iy), &
!      !(x2-x_in)/(x2-x1)*xxeuv(ix,iy+1)+(x_in-x1)/(x2-x1)*xxeuv(ix+1,iy+1)
!
!#undef INTERPOLATION_LIN
!
!  end subroutine ProstFit2
!#endif

  !--------------------------------------------------------
  !        get source properties of star closter particels
  !--------------------------------------------------------
  subroutine rs_star_cluster(Mass, lumi, Trad, x_euv, alpha, heat_euv, x_fuv, hhm, &
      lumeuv, lumfuv, sig_euv, sig_fuv, rOII)
    real(kind=DBL_KIND),intent(IN) :: Mass  !mass in [M_sun]
    real(kind=DBL_KIND),intent(OUT):: lumi  !in cgs
    real(kind=DBL_KIND),intent(OUT):: Trad  !in K
    real(kind=DBL_KIND),intent(OUT):: x_euv
    real(kind=DBL_KIND),intent(OUT):: alpha ! mean cross section
    real(kind=DBL_KIND),intent(OUT):: heat_euv ! heating rate per one ionization event
    real(kind=DBL_KIND),intent(OUT):: x_fuv
    real(kind=DBL_KIND),intent(OUT):: hhm
    real(kind=DBL_KIND),intent(OUT):: lumeuv, lumfuv, sig_euv, sig_fuv, rOII


    character(100) :: path2fitdat = "./ProstFit/rs_chabrier.dat"
    character(100) :: path2fitdat2 = "./ProstFit/rs_chabrier2.dat"
    character(100) :: path2fitdat3 = "./ProstFit/rs_chabrier3.dat"

    integer,parameter :: FH = 11
    integer,save :: ifirst=0
    character(len=100) :: ffn
    integer :: err

    real(kind=DBL_KIND),save :: L_rec, S_rec, Teff_rec, huv_rec, alpha_rec, heat_rec, hfuv_rec
    real(kind=DBL_KIND),save :: lumeuv_rec, lumfuv_rec, sig_euv_rec, sig_fuv_rec, rOII_rec

    if(ifirst == 0)then

#if FOL_RADSOURCE == 0
        ffn= trim(path2fitdat)
        if(get_myrank() == PRIMARY_RANK) print *, 'use stellar table at 1Zsun'
#elif FOL_RADSOURCE == 1
        ffn= trim(path2fitdat3)
        if(get_myrank() == PRIMARY_RANK) print *, 'use stellar table at 0.1Zsun'
#elif FOL_RADSOURCE == 2
        ffn= trim(path2fitdat2)
        if(get_myrank() == PRIMARY_RANK) print *, 'use stellar table at 0.01Zsun'
#else
        print *, 'rs_chabrier is not defined with FOL_RADSOURCE=', FOL_RADSOURCE
        stop
#endif

        if(get_myrank() == PRIMARY_RANK) print *, "reading: ", ffn
        open(FH, file=ffn,status='old',iostat=err)
        if (err > 0) then
           if(get_myrank() == PRIMARY_RANK) print '(A,/,A)', "RS_STAR_CLUSTER: file not found","stopping..."
           stop
        end if
        read(FH,*,iostat=err) L_rec, S_rec, Teff_rec, huv_rec, alpha_rec, heat_rec, hfuv_rec
        read(FH,*,iostat=err) lumeuv_rec, lumfuv_rec, sig_euv_rec, sig_fuv_rec, rOII_rec
        close(FH)

        !if(get_myrank() == PRIMARY_RANK) then
    
        !  print *, "-------------- rs values ---------------"
        !  print *, "L:", L_rec, "[ Lsun g^-1 ]"
        !  print *, "Teff:", Teff_rec, "[ K ]"
        !  print *, "x_euv", huv_rec
        !  print *, "x_fuv", hfuv_rec 
        !  print *, "heat_euv"

        !  print *, "read rs..", L_rec, S_rec, Teff_rec, huv_rec &
        !  , alpha_rec, heat_rec, hfuv_rec
        !endif

        ! L_rec [Lsun/Msun], S_rec [ 1/s/Msun]

        ifirst = 1
    endif

    
    lumi      = L_rec*Mass*cgs_lsun ! luminosity in cgs [cgs]
    Trad      = Teff_rec            ! effective temperature 
    x_euv     = huv_rec             ! emissivity of uv photons  [1.e0/erg]
    alpha     = alpha_rec           ! mean cross section [cm^-2]
    heat_euv  = heat_rec            ! heating rate  [erg]
    x_fuv     = hfuv_rec            ! fuv photons  [1.0/erg]
    hhm       = 0.e0                ! we neglect ionization of H- here
    lumeuv    = lumeuv_rec          ! euv luminosity rate
    lumfuv    = lumfuv_rec          ! fuv luminosity rate
    sig_euv   = sig_euv_rec         ! cross section for euv [cm^-2]
    sig_fuv   = sig_fuv_rec         ! cross section for fuv [cm^-2]
    rOII      = rOII_rec



  end subroutine rs_star_cluster

    !------------------------------------------------
    ! linear interpolation for 2D array    
    !------------------------------------------------
    subroutine lin2D(x_arr,nx,y_arr,ny,z_arr,x,y,z)
      integer,intent(IN)  :: nx, ny
      real(kind=DBL_KIND),intent(IN) :: x_arr(0:nx-1),y_arr(0:ny-1),z_arr(0:nx-1,0:ny-1),x,y
      real(kind=DBL_KIND),intent(Out) :: z

    real(kind=DBL_KIND) :: x1,x2,y1,y2,z_y1,z_y2,x_in,y_in
    integer :: i,ix,j,iy

    !no extrapolation
    x_in = min(max(x,x_arr(0)),x_arr(nx-1))
    y_in = min(max(y,y_arr(0)),y_arr(ny-1))    


    !find x bin
    ix = nx-2 !ix = nx-2 when x > x_arr(nx-1)
    do i = 0, nx-2
       if (x_in < x_arr(i+1)) then
          ix = i !x_arr(ix) < x < x_arr(ix+1)
          exit
       end if
    end do

    !find y bin
    iy = ny-2 !iy = ny-2 when y > y_arr(nx-1)
    do j = 0, ny-2
       if (y_in < y_arr(j+1)) then
          iy = j !y_arr(iy) < y < y_arr(iy+1)
          exit
       end if
    end do

    !interpolation
    x1 = x_arr(ix)
    x2 = x_arr(ix+1)
    y1 = y_arr(iy)
    y2 = y_arr(iy+1)

    !interp. in x-direction
    z_y1 = (x2-x_in)/(x2-x1) * z_arr(ix,iy) + (x_in-x1)/(x2-x1) * z_arr(ix+1,iy)
    z_y2 = (x2-x_in)/(x2-x1) * z_arr(ix,iy+1) + (x_in-x1)/(x2-x1) * z_arr(ix+1,iy+1)

    !interp. in y-direction
    z = (y2-y_in)/(y2-y1) * z_y1 + (y_in-y1)/(y2-y1) * z_y2

  end subroutine lin2D

    !------------------------------------------------
    ! linear interpolation/extraporation for 2D array    
    !------------------------------------------------
    subroutine lin2D_ex(x_arr,nx,y_arr,ny,z_arr,x,y,z)
      integer,intent(IN)  :: nx, ny
      real(kind=DBL_KIND),intent(IN) :: x_arr(0:nx-1),y_arr(0:ny-1),z_arr(0:nx-1,0:ny-1),x,y
      real(kind=DBL_KIND),intent(Out) :: z

    real(kind=DBL_KIND) :: x1,x2,y1,y2,z_y1,z_y2
    integer :: i,ix,j,iy

    !find x bin
    ix = nx-2 !ix = nx-2 when x > x_arr(nx-1)
    do i = 0, nx-2
       if (x < x_arr(i+1)) then
          ix = i !x_arr(ix) < x < x_arr(ix+1)
          exit
       end if
    end do

    !find y bin
    iy = ny-2 !iy = ny-2 when y > y_arr(nx-1)
    do j = 0, ny-2
       if (y < y_arr(j+1)) then
          iy = j !y_arr(iy) < y < y_arr(iy+1)
          exit
       end if
    end do

    !interpolation
    x1 = x_arr(ix)
    x2 = x_arr(ix+1)
    y1 = y_arr(iy)
    y2 = y_arr(iy+1)

    !interp. in x-direction
    z_y1 = (x2-x)/(x2-x1) * z_arr(ix,iy) + (x-x1)/(x2-x1) * z_arr(ix+1,iy)
    z_y2 = (x2-x)/(x2-x1) * z_arr(ix,iy+1) + (x-x1)/(x2-x1) * z_arr(ix+1,iy+1)

    !interp. in y-direction
    z = (y2-y)/(y2-y1) * z_y1 + (y-y1)/(y2-y1) * z_y2

  end subroutine lin2D_ex



end module radiationSource
