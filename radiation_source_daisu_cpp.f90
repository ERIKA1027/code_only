module radiationSource
  use grid
  use unit
  use parameter
  use mpilib
  use overBlockCoordinates
  use kinzoku
  implicit none
  private
  integer,parameter :: MAX_RADIATION_SOUCE = 10000 
  type t_rs_info
     integer :: nsource 
     integer,dimension(0:MAX_RADIATION_SOUCE-1) :: sid 
     real(kind=8),dimension(0:2,0:MAX_RADIATION_SOUCE-1) :: spos 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: lum 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: Trad 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: x_euv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: x_fuv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: alpha_euv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: heat_euv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: hhm 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: lumeuv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: lumfuv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: sig_euv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: sig_fuv 
     real(kind=8),dimension(0:MAX_RADIATION_SOUCE-1) :: rOII 
  end type t_rs_info
  type(t_rs_info), target :: rs_info 
  public :: t_rs_info, rs_GetSourceInfo
  public :: radiation_source
   real(kind=8) :: MP_dustsubl
   public :: MP_dustsubl
contains
subroutine radiation_source
  use modelParameter
  use sinkParticle
  use primordial,only : ProstFit2, rad_others
  use string, only : concat, num2char, CHARLEN
  use io_util, only : read_env
  real(kind=8),parameter :: mass_cr = 0.1d0 
  integer :: nsource, sid,n,np, nparticle,sp_Level,ns
  real(kind=8),dimension(0:2) :: pos
  real(kind=8) :: lum,radius,Trad,mass,dmass
  integer,dimension(:),allocatable :: pid
  real(kind=8),dimension(:),allocatable :: pmass,pmdot,pmdot_disk,pt_prev,pdm_disk
  real(kind=8),dimension(:),allocatable :: ptcrt
  real(kind=8) :: tage
  real(kind=8),dimension(:,:),allocatable :: pr,pJ_disk
  integer :: n_rs 
  integer,parameter :: LUN = 11
  character(len=CHARLEN) :: file, dir
  character(len=CHARLEN),parameter :: FILENAME_LOG='logRadiationSource'
  integer,parameter :: log_skip = 10
  type(rad_others) :: rda
  integer,save :: ifirst = 0
   real(kind=8),parameter :: eta = 1.d-1
   real(kind=8),parameter :: pi = 3.14159265
   real(kind=8),parameter :: sigma_tom = 6.65d-25
   real(kind=8) :: p, nu0
   real(kind=8) :: nu_euv_min, nu_euv_max, nu_fuv_min, nu_fuv_max, nu_max, nu_min
   real(kind=8) :: S_euv, S_fuv, sigma_mean, heat_cal1, heat_cal2, lumi_euv, lumi_fuv
   real(kind=8) :: eta_edd, lum_edd, lum0, x_dg, lum_edd_dg
   real(kind=8) :: acrate_bhl, lambda
   real(kind=8) :: mdot_norm,mdot
   real(kind=8) :: time_HL, sim_time
  nparticle = sp_getNparticle()
  if (ifirst == 0) then
     if(nparticle == 0) then
        sp_Level = sp_getLevel()
        pos = (/0.d0, 0.d0, 0.d0/)
        mass = MP_Mstar 
        call sp_newParticle(mass, pos, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), &
         mdot_disk=MP_Mdot/(Unit_msun/Unit_yr) , J_disk=(/0.d0, 0.d0, 1.d0/)) 
        nparticle = sp_getNparticle()
     end if
     ifirst = 1
  end if
  if (nparticle > 0) then 
     allocate(pid(nparticle), pmass(nparticle), pr(0:2, nparticle), pmdot(nparticle), pJ_disk(0:2, nparticle),&
      pt_prev(nparticle), pdm_disk(nparticle), ptcrt(nparticle))
     call sp_sinkdata2array(np, pmass, pr=pr, pid=pid, pmdot=pmdot, pJ_disk=pJ_disk, pt_prev=pt_prev &
       , pdm_disk=pdm_disk, ptcrt=ptcrt)
  else
     return 
  end if
  if(get_myrank() == 0) &
       print '(/,A)', "(RadSrc) Begin evaluating radiation sources -----------------------------------"
  rs_info%nsource = 0
  n_rs = 0 
  do n = 0, nparticle-1
     ns = n + 1 
     mdot = pmdot(ns)
     mass = pmass(ns)
     lum_edd = 4.d0 * pi * cgs_c * cgs_gc * 1.d4 * cgs_msun * cgs_mp / sigma_tom
     mdot_norm = mdot*(Unit_m/Unit_t) / (lum_edd/(eta*cgs_c*cgs_c))
     if (mdot_norm>2.0) then
        lum = 2.0 * lum_edd * (1.0+log(mdot_norm/2.0))
     else
        lum = lum_edd * mdot_norm
     end if
     MP_dustsubl = (lum* 2.8d2 / (16.0 * pi * MP_kappa_pl * cgs_sigma *((MP_Tsubl)**(4.d0)) ))**(0.5)
     MP_dustsubl = MP_dustsubl / Unit_l
     if(get_myrank() == 0) &
          print '(A,I6,A,2((1P1E12.4),A), 1P3E10.2,A,/,A,4((1P1E12.4),A))', "pid =",pid(ns),", mass=",pmass(ns)*Unit_msun, &
             " M_sun, mdot_disk=",mdot*Unit_msun/Unit_yr," M_sun yr^-1, J_disk=",pJ_disk(:,ns),&
             " (in code unit)","   ==>   L=", lum/cgs_lsun, " L_sun, T_rad=",Trad," K, Ndot_ion=" &
             , lum*rda%xeuv, " s^-1, Ndot_LW=", lum*rda%xfuv, " s^-1"
     if(get_myrank() == 0) then
        if (mod(Step(Lmin), log_skip) == 0) then
           call read_env('DIR', dir)
           file = concat(dir,FILENAME_LOG)
           open(LUN, file=file, position='APPEND')
           write(LUN, '(I14, 1P1E17.9, I5, 1P9E17.9)') Step(Lmin), Time(Lmin)*Unit_yr, &
                pid(ns), pmass(ns)*Unit_msun, mdot*Unit_msun/Unit_yr, &
                lum/cgs_lsun, radius, Trad, lum*rda%xeuv, lum*rda%xfuv, rda%lumeuv, rda%lumfuv, MP_dustsubl*Unit_pc
           call flush(LUN)
           close(LUN)
        end if
     end if
     rs_info%nsource = rs_info%nsource + 1
     if (rs_info%nsource == MAX_RADIATION_SOUCE) then
        print '(A,/,A)','(RadSrc) number of radiation sources exceeds the maximum', 'stopping...'
        stop
     end if
     p = - 1.5d0 
     nu_min = 6.0d0 * cgs_ev / cgs_hpk
     nu_max = 1.d3 * cgs_ev / cgs_hpk
     nu0 = nu_min
     lum0 = lum * nu0**p * (p + 1.d0) / (nu_max**(p + 1.d0) - nu_min**(p + 1.d0))
     nu_euv_min = 1.36d1 * cgs_ev / cgs_hpk
     nu_euv_max = 1.d3 * cgs_ev / cgs_hpk
     S_euv = lum0 *(nu_euv_max**p - nu_euv_min**p) / (cgs_hpk * p * nu0**p)
     rda%xeuv = S_euv / lum
     nu_fuv_min = 6.0d0 * cgs_ev / cgs_hpk
     nu_fuv_max = 1.36d1 * cgs_ev / cgs_hpk
     S_fuv = lum0 *(nu_fuv_max**p - nu_fuv_min**p) / (cgs_hpk * p * nu0**p)
     rda%xfuv = S_fuv / lum
     sigma_mean = 6.3d0 * 1.d-18 * lum0 * nu_euv_min**3.d0 *(nu_euv_max**(p-3.d0) - nu_euv_min**(p - 3.d0)) / (cgs_hpk * (nu0**p) *&
& (p - 3.d0))
     rda%alpha_euv = sigma_mean / S_euv
     heat_cal1 = lum0 * (nu_euv_max**(p + 1.d0) - nu_euv_min**(p + 1.d0)) / ((p + 1.d0) * nu0**p)
     heat_cal2 = lum0 * nu_euv_min * (nu_euv_max**p - nu_euv_min**p) / (p * nu0**p)
     rda%heat_euv = (heat_cal1 - heat_cal2) / S_euv
     lumi_euv = lum0 * (nu_euv_max**(p + 1.d0) - nu_euv_min**(p + 1.d0)) / ((p + 1.d0) * nu0**p)
     rda%lumeuv = lumi_euv / lum
     lumi_fuv = lum0 * (nu_fuv_max**(p + 1.d0) - nu_fuv_min**(p + 1.d0)) / ((p + 1.d0) * nu0**p)
     rda%lumfuv = lumi_fuv / lum
     rda%hhm = 0.d0
     rda%sig_euv = 2.8d2 * cgs_mh
     rda%sig_fuv = 2.8d2 * cgs_mh
     rda%rOII = 0.d0
     rs_info%sid(n_rs) = pid(ns)
     rs_info%spos(:,n_rs) = pr(:,ns)
     rs_info%lum(n_rs) = lum 
     rs_info%x_euv(n_rs) = S_euv / lum 
     rs_info%x_fuv(n_rs) = S_fuv / lum 
     rs_info%alpha_euv(n_rs) = sigma_mean / S_euv 
     rs_info%heat_euv(n_rs) = (heat_cal1 - heat_cal2) / S_euv 
     rs_info%hhm(n_rs) = 0.d0 
     rs_info%lumeuv(n_rs) = lumi_euv / lum 
     rs_info%lumfuv(n_rs) = lumi_fuv / lum 
     rs_info%sig_euv(n_rs) = 2.8d2 * cgs_mh 
     rs_info%sig_fuv(n_rs) = 2.8d2 * cgs_mh 
     rs_info%rOII(n_rs) = 0.d0 
     n_rs = n_rs + 1
  end do
  if(get_myrank() == 0) &
       print '(A)', "(RadSrc) End evaluating radiation sources -----------------------------------"
  deallocate(pid, pmass, pr, pmdot, pJ_disk, pt_prev, pdm_disk)
  end subroutine radiation_source
  subroutine rs_GetSourceInfo(p_rs_info)
    type(t_rs_info),pointer,intent(OUT) :: p_rs_info
    p_rs_info => rs_info
  end subroutine rs_GetSourceInfo
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
  subroutine ProstFit(Mass, Mdot, Rad, Lum)
    real(kind=8),intent(IN) :: Mass 
    real(kind=8),intent(IN) :: Mdot 
    real(kind=8),intent(Out) :: Rad 
    real(kind=8),intent(Out) :: Lum 
    integer,parameter :: nmd=11 
    integer,parameter :: npt=250 
    character(100) :: fname(0:nmd-1) = (/"1em5fit.dat", "1em4fit.dat", "3em4fit.dat", "1em3fit.dat", &
         "3em3fit.dat", "6em3fit.dat", "1em2fit.dat", "3em2fit.dat", &
         "6em2fit.dat", "1em1fit.dat", "3em1fit.dat"/)
    character(100) :: path2fitdat = "./ProstFit/"
    integer,parameter :: FH = 11
    character(len=100) :: ffn
    integer :: err
    real(kind=8),save :: xm(0:npt-1) 
    real(kind=8),save :: xr(0:npt-1,0:nmd-1) 
    real(kind=8),save :: xls(0:npt-1,0:nmd-1) 
    real(kind=8),parameter :: xmd(0:nmd-1) = (/1.d-5, 1.d-4, 3.d-4, 1.d-3, 3.d-3, 6.d-3,&
         1.d-2, 3.d-2, 6.d-2, 1.d-1, 3.d-1/) 
    real(kind=8) :: xxm
    integer,save :: ifirst=0
    integer :: i, id
    if (ifirst == 0) then
       do id = 0, nmd-1
          ffn= trim(path2fitdat)//trim(fname(id))
          if(get_myrank() == 0) print *, "PROST_FIT: reading ", ffn
          open(FH, file=ffn,status='old',iostat=err)
          if (err > 0) then
             if(get_myrank() == 0) print '(A,/,A)', "PROST_FIT: file not found","stopping..."
             stop
          end if
          do i = 0, npt-1
             read(FH,*,iostat=err) xxm, xr(i,id), xls(i,id)
             if (id == 0) xm(i) = 1d1**xxm

             if (err > 0) then
                print *, "PROST_FIT: **WARNING** data size might be inconsistent"
                exit 
             end if
          end do
          close(FH)
       end do
       ifirst = 1
    end if
    if(get_myrank() == 0) then
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
  end subroutine ProstFit
  subroutine rs_star_cluster(Mass, lumi, Trad, x_euv, alpha, heat_euv, x_fuv, hhm, &
      lumeuv, lumfuv, sig_euv, sig_fuv, rOII)
    real(kind=8),intent(IN) :: Mass 
    real(kind=8),intent(OUT):: lumi 
    real(kind=8),intent(OUT):: Trad 
    real(kind=8),intent(OUT):: x_euv
    real(kind=8),intent(OUT):: alpha 
    real(kind=8),intent(OUT):: heat_euv 
    real(kind=8),intent(OUT):: x_fuv
    real(kind=8),intent(OUT):: hhm
    real(kind=8),intent(OUT):: lumeuv, lumfuv, sig_euv, sig_fuv, rOII
    character(100) :: path2fitdat = "./ProstFit/rs_chabrier.dat"
    character(100) :: path2fitdat2 = "./ProstFit/rs_chabrier2.dat"
    character(100) :: path2fitdat3 = "./ProstFit/rs_chabrier3.dat"
    integer,parameter :: FH = 11
    integer,save :: ifirst=0
    character(len=100) :: ffn
    integer :: err
    real(kind=8),save :: L_rec, S_rec, Teff_rec, huv_rec, alpha_rec, heat_rec, hfuv_rec
    real(kind=8),save :: lumeuv_rec, lumfuv_rec, sig_euv_rec, sig_fuv_rec, rOII_rec
    if(ifirst == 0)then
        print *, 'rs_chabrier is not defined with FOL_RADSOURCE=', "Z1"
        stop
        if(get_myrank() == 0) print *, "reading: ", ffn
        open(FH, file=ffn,status='old',iostat=err)
        if (err > 0) then
           if(get_myrank() == 0) print '(A,/,A)', "RS_STAR_CLUSTER: file not found","stopping..."
           stop
        end if
        read(FH,*,iostat=err) L_rec, S_rec, Teff_rec, huv_rec, alpha_rec, heat_rec, hfuv_rec
        read(FH,*,iostat=err) lumeuv_rec, lumfuv_rec, sig_euv_rec, sig_fuv_rec, rOII_rec
        close(FH)
        ifirst = 1
    endif
    lumi = L_rec*Mass*cgs_lsun 
    Trad = Teff_rec 
    x_euv = huv_rec 
    alpha = alpha_rec 
    heat_euv = heat_rec 
    x_fuv = hfuv_rec 
    hhm = 0.e0 
    lumeuv = lumeuv_rec 
    lumfuv = lumfuv_rec 
    sig_euv = sig_euv_rec 
    sig_fuv = sig_fuv_rec 
    rOII = rOII_rec
  end subroutine rs_star_cluster
    subroutine lin2D(x_arr,nx,y_arr,ny,z_arr,x,y,z)
      integer,intent(IN) :: nx, ny
      real(kind=8),intent(IN) :: x_arr(0:nx-1),y_arr(0:ny-1),z_arr(0:nx-1,0:ny-1),x,y
      real(kind=8),intent(Out) :: z
    real(kind=8) :: x1,x2,y1,y2,z_y1,z_y2,x_in,y_in
    integer :: i,ix,j,iy
    x_in = min(max(x,x_arr(0)),x_arr(nx-1))
    y_in = min(max(y,y_arr(0)),y_arr(ny-1))
    ix = nx-2 
    do i = 0, nx-2
       if (x_in < x_arr(i+1)) then
          ix = i 
          exit
       end if
    end do
    iy = ny-2 
    do j = 0, ny-2
       if (y_in < y_arr(j+1)) then
          iy = j 
          exit
       end if
    end do
    x1 = x_arr(ix)
    x2 = x_arr(ix+1)
    y1 = y_arr(iy)
    y2 = y_arr(iy+1)
    z_y1 = (x2-x_in)/(x2-x1) * z_arr(ix,iy) + (x_in-x1)/(x2-x1) * z_arr(ix+1,iy)
    z_y2 = (x2-x_in)/(x2-x1) * z_arr(ix,iy+1) + (x_in-x1)/(x2-x1) * z_arr(ix+1,iy+1)
    z = (y2-y_in)/(y2-y1) * z_y1 + (y_in-y1)/(y2-y1) * z_y2
  end subroutine lin2D
    subroutine lin2D_ex(x_arr,nx,y_arr,ny,z_arr,x,y,z)
      integer,intent(IN) :: nx, ny
      real(kind=8),intent(IN) :: x_arr(0:nx-1),y_arr(0:ny-1),z_arr(0:nx-1,0:ny-1),x,y
      real(kind=8),intent(Out) :: z
    real(kind=8) :: x1,x2,y1,y2,z_y1,z_y2
    integer :: i,ix,j,iy
    ix = nx-2 
    do i = 0, nx-2
       if (x < x_arr(i+1)) then
          ix = i 
          exit
       end if
    end do
    iy = ny-2 
    do j = 0, ny-2
       if (y < y_arr(j+1)) then
          iy = j 
          exit
       end if
    end do
    x1 = x_arr(ix)
    x2 = x_arr(ix+1)
    y1 = y_arr(iy)
    y2 = y_arr(iy+1)
    z_y1 = (x2-x)/(x2-x1) * z_arr(ix,iy) + (x-x1)/(x2-x1) * z_arr(ix+1,iy)
    z_y2 = (x2-x)/(x2-x1) * z_arr(ix,iy+1) + (x-x1)/(x2-x1) * z_arr(ix+1,iy+1)
    z = (y2-y)/(y2-y1) * z_y1 + (y-y1)/(y2-y1) * z_y2
  end subroutine lin2D_ex
end module radiationSource
