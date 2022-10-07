
! HLLD flux
! reference: Miyoshi and Kusano, JCP, 208, (2005), 315-344
!
! coding rule
! S_L ........ sl
! S_R ........ sr
! S_L^* ...... sal
! S_R^* ...... sar
! S_M ........ sm
! rho_L ...... rhol
! rho_L^* .... rhoal
! rho_L^** ... rhoml

#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state.
!-------------------------------------------------------------------------
module eos

#ifndef RADTR_M1closer
  use modelParameter, only : MP_Tmin, MP_Tmax, MP_mu, MP_frac_COsum   !KS ADDED
#else
  use modelParameter, only : MP_Tmin, MP_Tmax, MP_Ctil_nD, MP_mu, MP_frac_COsum ! HF ADDED
#endif

  implicit none
#ifdef MODEL_SOLARWIND
  real(kind=DBL_KIND),parameter :: Gamma = 1.46d0 ! for solar wind
#else !MODEL_SOLARWIND
  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0
!!$  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0
!!$  real(kind=DBL_KIND),parameter :: Gamma = 1.01d0
#endif !MODEL_SOLARWIND
  real(kind=DBL_KIND),parameter :: Cr = 0.18d0 ! Cr := cp^2/ch (for divB clean)
  real(kind=DBL_KIND),save :: Ch
contains
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  function cyclecomp(ncrd, invert) result( mcycle )
    use grid , only : Mmin, Mmax
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(Mmin:Mmax) :: mcycle
    integer,dimension(MX:MZ) :: mcycleV = (/ MVX, MVY, MVZ /)
    integer,dimension(MX:MZ) :: mcycleB = (/ MBX, MBY, MBZ /)
    integer :: m
    do m = Mmin, Mmax
       mcycle(m) = m
    enddo
    if ( present( invert ) ) then
       mcycle(MVX:MVX+size(mcycleV)-1) = cshift( mcycleV, -ncrd)
       mcycle(MBX:MBX+size(mcycleB)-1) = cshift( mcycleB, -ncrd)
    else
       mcycle(MVX:MVX+size(mcycleV)-1) = cshift( mcycleV, ncrd )
       mcycle(MBX:MBX+size(mcycleB)-1) = cshift( mcycleB, ncrd )
    endif
  end function cyclecomp
  !-----------------------------------------------------------------------
  ! get numerical flux in one dimension (Miyoshi & Kusano, 2005, JCP 208, 315-344
  !-----------------------------------------------------------------------

  ! SWD(A, B) = 1.d0  if A >= B. otherwise 0.d0
#define SWD(A, B)  (0.5d0 + 0.5d0*sign(1.d0, (A)-(B)))
  ! SWR(A, B) = 1.d0 if A * B < 0, otherwise 0.d0
#define SWR(A, B)  (0.5d0 - 0.5d0*sign(1.d0, (A))*sign(1.d0, (B)))

#if MODEL_ART > 0 || MODEL_ART < 0
 ! give a gamma 
  subroutine flux( ql, qr, f, gam)
#else !MODEL_ART
  subroutine flux( ql, qr, f)
#endif !MODEL_ART
    use grid, only : Imin,Jmin,Kmin,Mmin, Imax,Jmax,Kmax,Mmax,Imingh, Jmingh, Kmingh, &
                     CellWidth, LevelMax, Lmin, Dtime, &
         globdbg_myrank, globdbg_mygid,globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
   ! use util, only : util_arroffset            
    use parameter, only : Pi4i, Pi4
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: qr  ! (IN)
  !  real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,MX:) :: pratio  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: f   ! (OUT)
    real(kind=DBL_KIND) :: sw_deg_l, sw_deg_r, swl, swr, swal, swar, swml, swmr, swt
    real(kind=DBL_KIND),parameter :: eps = 1.D-6
    integer :: i,j,k, is,js,ks,ie,je,ke,io,jo,ko
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,bxl,byl,bzl,pbl,ptl,el, &
         rhoal,ual,val,wal,pal,byal,bzal,pbal,ptal,eal, &
         rhoml,uml,vml,wml,pml,byml,bzml,pbml,ptml,eml, &
         cfl,csl2,cal2,caxl2,sl,sal,sql,psil,rhoil, &
         denom_l, &
                                !
         rhor,ur,vr,wr,pr,bxr,byr,bzr,pbr,ptr,er, &
         rhoar,uar,var,war,par,byar,bzar,pbar,ptar,ear, &
         rhomr,umr,vmr,wmr,pmr,bymr,bzmr,pbmr,ptmr,emr, &
         cfr,csr2,car2,caxr2,sr,sar,sqr,psir,rhoir, &
         denom_r, &
                                !
         gm1i,sqrtpi4i,sqrtpi4,bxm,psim,sm,pta,sqlri,vm,wm,bym,bzm,gzm,signbxm,bxm2, &
         rho, u, v, w, bx, by, bz, p, e, pt, &
         gam1, cl, ub
    real(kind=DBL_KIND) :: sql1, sqr1, sqa1, sqa
    real(kind=DBL_KIND), dimension(NCEHM_MIN:NCEHM_MAX) :: rho_cheml, rho_chemr, rhoacheml, rhoachemr, rhomcheml, rhomchemr &
                                                           , rhochem
#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND) :: rhocol, rhocor, rhoacol, rhoacor, rhomcol, rhomcor, rhoco
#endif


   
#if MODEL_ART > 0
    real(kind=DBL_KIND),dimension(:,:,:),intent(IN) :: gam
#endif !MODEL_ART
 
#ifdef SINGLE_STEP
    h = CellWidth( :, LevelMax )
#else  !SINGLE_STEP
    h = CellWidth( :, Lmin )
#endif !SINGLE_STEP
    Ch = (CFL) * minval(h) / Dtime(Lmin) /3.d0

    is = max( lbound(ql,1), lbound(qr,1), lbound(f,1) )
    js = max( lbound(ql,2), lbound(qr,2), lbound(f,2) )
    ks = max( lbound(ql,3), lbound(qr,3), lbound(f,3) )
    ie = min( ubound(ql,1), ubound(qr,1), ubound(f,1) )
    je = min( ubound(ql,2), ubound(qr,2), ubound(f,2) )
    ke = min( ubound(ql,3), ubound(qr,3), ubound(f,3) )

    !gm1i = 1.d0/(Gamma - 1)
    sqrtpi4i = sqrt(Pi4i)
    sqrtpi4 = sqrt(Pi4)
    do k = ks, ke
       do j = js, je
          do i = is, ie
#if MODEL_ART > 0 || MODEL_ART < 0
         gam1 = gam(i,j,k) - 1
         gm1i = 1.d0/(gam(i,j,k) - 1)
#else !MODEL_ART
         gam1 = Gamma - 1
         gm1i = 1.d0/(Gamma - 1)
#endif !MODEL_ART
             ! -------------------------
             ! add-on for div B cleaning
             ! -------------------------
             bxl = ql(i,j,k,MBX)*sqrtpi4i
             bxr = qr(i,j,k,MBX)*sqrtpi4i
             psil = ql(i,j,k,MDB)*sqrtpi4i
             psir = qr(i,j,k,MDB)*sqrtpi4i
             bxm  = 0.5d0*(bxl+bxr - (psir-psil)/Ch)
             psim = 0.5d0*(psil+psir - (bxr-bxl)*Ch)
             f(i,j,k,MBX) = psim         * sqrtpi4
             f(i,j,k,MDB) = bxm*Ch**2    * sqrtpi4
             bxm2 = bxm**2
             ! --------------------
             ! U_L; left variables
             ! --------------------
             rhol = abs(ql(i,j,k,MRHO)) ! 真真真真真真真真真?
#if MODEL_ART > 0
             rho_cheml(NCEHM_MIN:NCEHM_MAX) = abs(rhol*ql(i,j,k,NCEHM_MIN:NCEHM_MAX))
  #ifdef CHEM_MODEL_HF2020
             rhocol = abs(rhol*ql(i,j,k,MCO)) 
  #endif

#endif !MODEL_ART

             ul = ql(i,j,k,MVX)
             vl = ql(i,j,k,MVY)
             wl = ql(i,j,k,MVZ)
             byl = ql(i,j,k,MBY)*sqrtpi4i
             bzl = ql(i,j,k,MBZ)*sqrtpi4i
             pl = abs(ql(i,j,k,MP))
             pbl = (bxm2 + byl**2 + bzl**2)*0.5d0
             ptl = pl + pbl
             el = rhol * (ul**2 + vl**2 + wl**2) * 0.5d0 + pbl + pl*gm1i
             ! ---------------------
             ! U_R; right variables
             ! ---------------------
             rhor = abs(qr(i,j,k,MRHO)) ! 真真真真真真真真真?

#if MODEL_ART > 0
             rho_chemr(NCEHM_MIN:NCEHM_MAX) = abs(rhor*qr(i,j,k,NCEHM_MIN:NCEHM_MAX))
  #ifdef CHEM_MODEL_HF2020
             rhocor = abs(rhor*qr(i,j,k,MCO))
  #endif
#endif !MODEL_ART

             ur = qr(i,j,k,MVX)
             vr = qr(i,j,k,MVY)
             wr = qr(i,j,k,MVZ)
             byr = qr(i,j,k,MBY)*sqrtpi4i
             bzr = qr(i,j,k,MBZ)*sqrtpi4i
             pr = abs(qr(i,j,k,MP))
             pbr = (bxm2 + byr**2 + bzr**2)*0.5d0
             ptr = pr + pbr
             er = rhor * (ur**2 + vr**2 + wr**2) * 0.5d0 + pbr + pr*gm1i
             ! -----------
             ! wave speed
             ! -----------
             rhoil = 1.d0/rhol
             rhoir = 1.d0/rhor
#if MODEL_ART > 0 || MODEL_ART < 0
             csl2 = gam(i,j,k) * pl * rhoil ! sound speed
             csr2 = gam(i,j,k) * pr * rhoir ! sound speed
#else !MODEL_ART
             csl2 = Gamma * pl * rhoil ! sound speed
             csr2 = Gamma * pr * rhoir ! sound speed
#endif !MODEL_ART
             cal2 = pbl * rhoil * 2.d0 ! alfven speed
             car2 = pbr * rhoir * 2.d0 ! alfven speed
             caxl2 = bxm2 * rhoil    ! alfven speed in x-direction
             caxr2 = bxm2 * rhoir    ! alfven speed in x-direction
             cfl = sqrt( 0.5d0*( csl2 + cal2 + sqrt( (csl2 + cal2)**2 - 4.d0*csl2*caxl2) ) ) ! fast wave
             cfr = sqrt( 0.5d0*( csr2 + car2 + sqrt( (csr2 + car2)**2 - 4.d0*csr2*caxr2) ) ) ! fast wave
             sl = min(ul-cfl, ur-cfr)
             sr = max(ul+cfl, ur+cfr)
             sm = ((sr - ur)*rhor*ur - (sl - ul)*rhol*ul - ptr + ptl)/((sr - ur)*rhor - (sl - ul)*rhol)


             rhoal = rhol * (sl - ul)/(sl - sm)
             rhoar = rhor * (sr - ur)/(sr - sm)
#if MODEL_ART > 0
             ! left
             rhoacheml(NCEHM_MIN:NCEHM_MAX) = rho_cheml(NCEHM_MIN:NCEHM_MAX) * (sl - ul)/(sl - sm)
             ! right
             rhoachemr(NCEHM_MIN:NCEHM_MAX) = rho_chemr(NCEHM_MIN:NCEHM_MAX) * (sr - ur)/(sr - sm) 

  #ifdef CHEM_MODEL_HF2020
             rhoacol  = rhocol * (sl - ul)/(sl - sm)        
             rhoacor  = rhocor * (sr - ur)/(sr - sm)
  #endif
#endif !MODEL_ART

             sql = sqrt(rhoal)
             sqr = sqrt(rhoar)
             sal = sm - abs(bxm) / sql
             sar = sm + abs(bxm) / sqr
             ! ---------------------------
             ! switching flux
             ! ---------------------------
             swl  = SWD(sl, 0.d0)
             swal = SWR(sl, sal)
             swml = SWR(sal, sm)
             swmr = SWR(sm, sar)
             swar = SWR(sar, sr)
             swr  = SWD(0.d0, sr)
             swt = 1.d0/(swl+swr+swal+swar+swml+swmr)
             ! ---------------
             ! U*_L and U*_R
             ! ---------------
             pta = ((sr-ur)*rhor*ptl - (sl-ul)*rhol*ptr + rhol*rhor*(sr-ur)*(sl-ul)*(ur-ul))/((sr-ur)*rhor-(sl-ul)*rhol)

             denom_l = rhol*(sl-ul)*(sl-sm)-bxm2
             sw_deg_l = SWD(eps, abs(denom_l)) ! switch for degenerate
             denom_l = (1.d0-sw_deg_l) / (denom_l + sw_deg_l)
             ual = sm
             val = vl - byl*bxm*(sm-ul) * denom_l
             wal = wl - bzl*bxm*(sm-ul) * denom_l
             byal = byl*(rhol*(sl-ul)**2-bxm2) * denom_l
             bzal = bzl*(rhol*(sl-ul)**2-bxm2) * denom_l
             ptal = pta
             eal = ((sl-ul)*el-ptl*ul+ptal*sm+bxm*(ul*bxm+vl*byl+wl*bzl-ual*bxm-val*byal-wal*bzal))/(sl-sm)

             denom_r = rhor*(sr-ur)*(sr-sm)-bxm2
             sw_deg_r = SWD(eps, abs(denom_r))
             denom_r = (1.d0-sw_deg_r) / (denom_r + sw_deg_r)
             uar = sm
             var = vr - bxm*byr*(sm-ur) * denom_r
             war = wr - bxm*bzr*(sm-ur) * denom_r
             byar = byr*(rhor*(sr-ur)**2-bxm2) * denom_r
             bzar = bzr*(rhor*(sr-ur)**2-bxm2) * denom_r
             ptar = pta
             ear = ((sr-ur)*er-ptr*ur+ptar*sm+bxm*(ur*bxm+vr*byr+wr*bzr-uar*bxm-var*byar-war*bzar))/(sr-sm)
             ! --------------
             ! U**_L, U**_R
             ! --------------
             sqlri = 1.d0/(sql+sqr)
             signbxm = sign(1.d0,bxm)
             vm = (sql*val+sqr*var+(byar-byal)*signbxm)*sqlri
             wm = (sql*wal+sqr*war+(bzar-bzal)*signbxm)*sqlri
             bym = (sql*byar+sqr*byal+sql*sqr*(var-val)*signbxm)*sqlri
             bzm = (sql*bzar+sqr*bzal+sql*sqr*(war-wal)*signbxm)*sqlri

             rhoml = rhoal

#if MODEL_ART > 0
             rhomcheml(NCEHM_MIN:NCEHM_MAX) = rhoacheml(NCEHM_MIN:NCEHM_MAX) 
  #ifdef CHEM_MODEL_HF2020
             rhomcol = rhoacol
  #endif
#endif !MODEL_ART

             ptml = ptal
             uml = sm
             vml = vm
             wml = wm
             byml = bym
             bzml = bzm
             eml = eal - sql*(ual*bxm+val*byal+wal*bzal - uml*bxm-vml*byml-wml*bzml)*signbxm

             rhomr = rhoar

#if MODEL_ART > 0
             rhomchemr(NCEHM_MIN:NCEHM_MAX) = rhoachemr(NCEHM_MIN:NCEHM_MAX)
  #ifdef CHEM_MODEL_HF2020
             rhomcor = rhoacor
  #endif
#endif !MODEL_ART


             ptmr = ptar
             umr = sm
             vmr = vm
             wmr = wm
             bymr = bym
             bzmr = bzm
             emr = ear + sqr*(uar*bxm+var*byar+war*bzar - umr*bxm-vmr*bymr-wmr*bzmr)*signbxm
             ! -------------
             ! switching U
             ! -------------
             rho = (rhol*swl + rhor*swr + rhoal*swal + rhoar*swar + rhoml*swml + rhomr*swmr)*swt

#if MODEL_ART > 0
             rhochem(NCEHM_MIN:NCEHM_MAX)=(rho_cheml(NCEHM_MIN:NCEHM_MAX)*swl+rho_chemr(NCEHM_MIN:NCEHM_MAX)*swr &
                      + rhoacheml(NCEHM_MIN:NCEHM_MAX)*swal+rhoachemr(NCEHM_MIN:NCEHM_MAX)*swar                  &
                      + rhomcheml(NCEHM_MIN:NCEHM_MAX)*swml+rhomchemr(NCEHM_MIN:NCEHM_MAX)*swmr)*swt
  #ifdef CHEM_MODEL_HF2020
             rhoco=(rhocol*swl+rhocor*swr+rhoacol*swal+rhoacor*swar+rhomcol*swml+rhomcor*swmr)*swt
  #endif
#endif !MODEL_ART


             u = (ul*swl + ur*swr + ual*swal + uar*swar + uml*swml + umr*swmr)*swt
             v = (vl*swl + vr*swr + val*swal + var*swar + vml*swml + vmr*swmr)*swt
             w = (wl*swl + wr*swr + wal*swal + war*swar + wml*swml + wmr*swmr)*swt
             bx = bxm
             by = (byl*swl + byr*swr + byal*swal + byar*swar + byml*swml + bymr*swmr)*swt
             bz = (bzl*swl + bzr*swr + bzal*swal + bzar*swar + bzml*swml + bzmr*swmr)*swt
             pt = (ptl*swl + ptr*swr + ptal*swal + ptar*swar + ptml*swml + ptmr*swmr)*swt
             e = (el*swl + er*swr + eal*swal + ear*swar + eml*swml + emr*swmr)*swt
             ! --------
             ! flux
             ! --------
             f(i,j,k,MRHO) = rho*u
#if MODEL_ART > 0
             f(i,j,k,NCEHM_MIN:NCEHM_MAX) = rhochem(NCEHM_MIN:NCEHM_MAX)*u 
  #ifdef CHEM_MODEL_HF2020
             f(i,j,k,MCO)  = rhoco*u
  #endif
#endif !MODEL_ART
             f(i,j,k, MVX) = rho*u**2 + pt - bx**2
             f(i,j,k, MVY) = rho*v*u - bx*by
             f(i,j,k, MVZ) = rho*w*u - bx*bz
             f(i,j,k, MBY) = (by*u - bx*v)*sqrtpi4
             f(i,j,k, MBZ) = (bz*u - bx*w)*sqrtpi4
             f(i,j,k,  MP) = (e + pt)*u - bx*(u*bx+v*by+w*bz)
             ! Up-wind flux of chemical abundances (ART CHEM)             
!#if MODEL_ART > 0

             !sql1 = sqrt(rhol)
             !sqr1 = sqrt(rhor)
             !sqa1 = sql1 + sqr1
             !ub = (sql1 * ul + sqr1 * ur) / sqa1

           !  sqa = sql + sqr              !ske try
           !  ub = (sql * ul + sqr * ur) / sqa
!prevoius method (I took from Sugimura-san code)
           !  f(i,j,k,MHN) = 0.5d0 * ub * (ql(i,j,k,MHN)*rhol + qr(i,j,k,MHN)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MHN)*rhor - ql(i,j,k,MHN)*rhol)
           !  f(i,j,k,MHP) = 0.5d0 * ub * (ql(i,j,k,MHP)*rhol + qr(i,j,k,MHP)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MHP)*rhor - ql(i,j,k,MHP)*rhol)
           !  f(i,j,k,MH2) = 0.5d0 * ub * (ql(i,j,k,MH2)*rhol + qr(i,j,k,MH2)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MH2)*rhor - ql(i,j,k,MH2)*rhol)
           !  f(i,j,k,MEL) = 0.5d0 * ub * (ql(i,j,k,MEL)*rhol + qr(i,j,k,MEL)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MEL)*rhor - ql(i,j,k,MEL)*rhol)
           !  f(i,j,k,MHM) = 0.5d0 * ub * (ql(i,j,k,MHM)*rhol + qr(i,j,k,MHM)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MHM)*rhor - ql(i,j,k,MHM)*rhol)
           !  f(i,j,k,MH2P) = 0.5d0 * ub * (ql(i,j,k,MH2P)*rhol + qr(i,j,k,MH2P)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MH2P)*rhor - ql(i,j,k,MH2P)*rhol)
!#ifdef CHEM_D
           !  f(i,j,k,MDN) = 0.5d0 * ub * (ql(i,j,k,MDN)*rhol + qr(i,j,k,MDN)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MDN)*rhor - ql(i,j,k,MDN)*rhol)
           !  f(i,j,k,MDP) = 0.5d0 * ub * (ql(i,j,k,MDP)*rhol + qr(i,j,k,MDP)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MDP)*rhor - ql(i,j,k,MDP)*rhol)
           !  f(i,j,k,MHDN) = 0.5d0 * ub * (ql(i,j,k,MHDN)*rhol + qr(i,j,k,MHDN)*rhor) &
           !       - 0.5d0 * abs(ub) * (qr(i,j,k,MHDN)*rhor - ql(i,j,k,MHDN)*rhol)
!#endif !CHEM_D
!#endif !MODEL_ART
          enddo
       enddo
    enddo
  end subroutine flux

  !-----------------------------------------------------------------------
  ! source term due to div B
  !-----------------------------------------------------------------------
  subroutine source_b(f, w, dt, gid)
    use grid, only : get_dv, get_ds, get_level, get_up, Imingh, Jmingh, Kmingh, Mmin, &
         Imin, Imax, Jmin, Jmax, Kmin, Kmax, CellWidth, Lmin
    use parameter
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:,MX:) :: f  ! (IN)
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: w  ! (OUT)
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=DBL_KIND) :: dv
    real(kind=DBL_KIND),dimension(MX:MZ) :: ds
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: i,j,k, level
    real(kind=DBL_KIND) :: divbdv, dtchcr, cr_len
    level = get_level(gid)
    u => get_up(gid)
    ds = get_ds(level)
    dv = get_dv(level)
    cr_len = Cr * maxval(CellWidth(:,Lmin) * (/(NI)*(NGI_BASE), (NJ)*(NGJ_BASE), (NK)*(NGK_BASE)/))
    dtchcr = dt*Ch/cr_len
    ! f(i,j,k,MBX,MX) = Psi @ x-surface
    ! f(i,j,k,MBY,MY) = Psi @ y-surface
    ! f(i,j,k,MBZ,MZ) = Psi @ z-surface
    ! f(i,j,k,MDB,MX) = ch^2 Bxm

    do k=Kmin,Kmax
       do j=Jmin,Jmax
          do i=Imin,Imax
             divbdv = &
                  ((f(i,j,k,MDB,MX) - f(i-1,j,k,MDB,MX))*ds(MX) &
                  +(f(i,j,k,MDB,MY) - f(i,j-1,k,MDB,MY))*ds(MY) &
                  +(f(i,j,k,MDB,MZ) - f(i,j,k-1,MDB,MZ))*ds(MZ)) &
                  /(4.d0*pi*Ch**2)
             w(i,j,k,MVX)=w(i,j,k,MVX)-dt*divbdv*u(i,j,k,MBX)
             w(i,j,k,MVY)=w(i,j,k,MVY)-dt*divbdv*u(i,j,k,MBY)
             w(i,j,k,MVZ)=w(i,j,k,MVZ)-dt*divbdv*u(i,j,k,MBZ)
             w(i,j,k,MDB)=w(i,j,k,MDB)*exp(-dtchcr)
             w(i,j,k,MP)=w(i,j,k,MP) -dt &
                  *(u(i,j,k,MBX)*(f(i,j,k,MBX,MX)-f(i-1,j,k,MBX,MX))*ds(MX) &
                  + u(i,j,k,MBY)*(f(i,j,k,MBY,MY)-f(i,j-1,k,MBY,MY))*ds(MY) &
                  + u(i,j,k,MBZ)*(f(i,j,k,MBZ,MZ)-f(i,j,k-1,MBZ,MZ))*ds(MZ) ) &
                  /(4.d0*pi)
          enddo
       enddo
    enddo
#ifdef EMULATE_2DIM
    w(:,:,:,MVZ)=0.d0
    w(:,:,:,MBZ)=0.d0
#endif
  end subroutine source_b
  !-----------------------------------------------------------------------
  ! find dt according to CFL condtion
  !-----------------------------------------------------------------------
  function get_dt_by_cflcond( id ) result( dt )
    use parameter
    use grid
    integer,intent(IN) :: id
    real(kind=DBL_KIND) :: dt
    integer :: i,j,k,level
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH) :: ca
#if MODEL_ART > 0
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH) :: gam
#endif !MODEL_ART
#ifdef RADTR_M1closer
    real(kind=DBL_KIND):: dt_radtr 
#endif
    level = get_level(id)
    h = CellWidth( :, level )
    rho => get_Ucomp(MRHO,id)
    vx => get_Ucomp(MVX,id)
    vy => get_Ucomp(MVY,id)
    vz => get_Ucomp(MVZ,id)
    bx => get_Ucomp(MBX,id)
    by => get_Ucomp(MBY,id)
    bz => get_Ucomp(MBZ,id)
    p => get_Ucomp(MP,id)
#ifdef WITH_SELFGRAVITY
    gx => get_Ucomp(MGX,id)
    gy => get_Ucomp(MGY,id)
    gz => get_Ucomp(MGZ,id)
#endif !WITH_SELFGRAVITY
    x => get_Xp(id)
    y => get_Yp(id)
    z => get_Zp(id)
#if MODEL_ART > 0
    call get_gamma(id,gam)
    ca= sqrt( abs( ( gam*p +(bx**2+by**2+bz**2)*pi4i )/rho ) ) ! fast wave
#else !MODEL_ART
    ca= sqrt( abs( ( Gamma*p +(bx**2+by**2+bz**2)*pi4i )/rho ) ) ! fast wave
#endif !MODEL_ART
#ifdef WITH_SELFGRAVITY
#define DTI_(V_,G_,H_) sqrt((V_)**2+2.d0*(G_)*(H_)) / (H_)
    dt = (CFL) / maxval( &
         DTI_( abs(vx)+ca, abs(gx), h(MX) ) + &
         DTI_( abs(vy)+ca, abs(gy), h(MY) ) + &
         DTI_( abs(vz)+ca, abs(gz), h(MZ) ), mask=GridMask )
#else   !WITH_SELFGRAVITY
    dt = (CFL) / maxval( &
         (abs(vx)+ca)/h(MX) + &
         (abs(vy)+ca)/h(MY) + &
         (abs(vz)+ca)/h(MZ), mask=GridMask )
#endif  !WITH_SELFGRAVITY

#ifdef RADTR_M1closer
    dt_radtr = minval(h)*CFLfac_radtr/ (3.d0*MP_Ctil_nD)*MAX_STEP_RADTR_TO_HYDRO
    dt = min(dt, dt_radtr)
#endif

  end function get_dt_by_cflcond

  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable w
  !
  !-----------------------------------------------------------------------
  subroutine u2w(u,w,dv)
    use parameter
    use grid, only : Mmin
    use unit
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART    
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv
    integer :: ichem
    
#if MODEL_ART > 0
!    real(kind=DBL_KIND),dimension(ARRAYSIZE3(u)) :: gam
    real(kind=DBL_KIND) :: T_K, xmu, gam
    integer :: i,j,k
#endif        

    w(:,:,:,MRHO)=u(:,:,:,MRHO)*dv
    w(:,:,:,MVX)=u(:,:,:,MVX)*w(:,:,:,MRHO)
    w(:,:,:,MVY)=u(:,:,:,MVY)*w(:,:,:,MRHO)
    w(:,:,:,MVZ)=u(:,:,:,MVZ)*w(:,:,:,MRHO)
    w(:,:,:,MBX)= u(:,:,:,MBX)*dv
    w(:,:,:,MBY)= u(:,:,:,MBY)*dv
    w(:,:,:,MBZ)= u(:,:,:,MBZ)*dv
    w(:,:,:,MDB) = u(:,:,:,MDB)*dv


#ifdef DM_POTENTIAL
    w(:,:,:,MDMRHO)=u(:,:,:,MDMRHO)*dv
#endif

    !ART CHEM
#if MODEL_ART > 0
    do ichem = NCEHM_MIN, NCEHM_MAX
      w(:,:,:,ichem)=u(:,:,:,ichem)*w(:,:,:,MRHO)
    enddo
  #ifdef CHEM_MODEL_HF2020
    w(:,:,:,MCO)=u(:,:,:,MCO)*w(:,:,:,MRHO)
  #endif

    w(:,:,:,MKPI)=u(:,:,:,MKPI) * dv
    w(:,:,:,MHPI)=u(:,:,:,MHPI) * dv
    w(:,:,:,MXPI)=u(:,:,:,MXPI)*w(:,:,:,MRHO)
    w(:,:,:,MYPI)=u(:,:,:,MYPI)*w(:,:,:,MRHO)
    w(:,:,:,MZPI)=u(:,:,:,MZPI)*w(:,:,:,MRHO) ! original 真真真真? HF
    !w(:,:,:,MXPI)=u(:,:,:,MXPI) * dv
    !w(:,:,:,MYPI)=u(:,:,:,MYPI) * dv
    !w(:,:,:,MZPI)=u(:,:,:,MZPI) * dv

  #ifdef METAL
    w(:,:,:,MTD )=u(:,:,:,MTD ) * dv
  #endif
  #ifdef DUST_NOTCONSTANT
    w(:,:,:,MDRHO)=u(:,:,:,MDRHO)*w(:,:,:,MRHO)
  #endif
  #ifdef METAL_TRANSFER
    w(:,:,:,MMET)=u(:,:,:,MMET)*w(:,:,:,MRHO)
  #endif

#endif !MODEL_ART

#ifdef RADTR_DIRECT
  #if MODEL_ART == 2
    w(:,:,:,MKPD)=u(:,:,:,MKPD) * dv
    #ifdef METAL
      w(:,:,:,MGFUV)=u(:,:,:,MGFUV) * dv
      w(:,:,:,MDPH) =u(:,:,:,MDPH) * dv
      w(:,:,:,MDPCO) =u(:,:,:,MDPCO) * dv
      w(:,:,:,MKPOII) =u(:,:,:,MKPOII) * dv
    #endif
  #endif !MODEL_ART
#endif

#ifdef MPSI
    w(:,:,:,MPSI) = u(:,:,:,MPSI) * dv
    w(:,:,:,MGX) = u(:,:,:,MGX) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGY) = u(:,:,:,MGY) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGZ) = u(:,:,:,MGZ) * u(:,:,:,MRHO) * dv
#endif !MPSI

#if MODEL_ART > 0
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !+MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
             T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb      !真 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax)                            !c_H2(T_K)真真真真真真真真真真真真?
             gam = 1.d0+(1.d0+4.d0*yHe) &                            !真? (Hm, H2p真?)
                  /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
             w(i,j,k,MP) =dv*(u(i,j,k,MP)/(gam-1.d0) &
                  +u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2)/2.d0 &
                  +(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*pi8i)
                  
             ! if (i==5 .and. j==5 .and. k==5)  then
             !    print *, "KS DEBUG1, gam(5,5,5)", gam(i,j,k),u(i,j,k,MHN),u(i,j,k,MH2),u(i,j,k,MEL),u(i,j,k,MHP),u(i,j,k,MRHO) ,xmu, T_K
             ! end if
          end do
       end do
    end do
    ! w(:,:,:,MP) =dv*(u(:,:,:,MP)/(gam(:,:,:)-1.d0) &
    !      +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2)
#else !MODEL_ART
    w(:,:,:,MP) =dv*(u(:,:,:,MP)/(Gamma-1.d0) &
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0 &
         +(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i)

#endif !MODEL_ART


#ifdef M1CLOSER_EUV_TRANSFER ! EUV photons
    w(:,:,:,MER)  = u(:,:,:,MER) *dv
    w(:,:,:,MFRX) = u(:,:,:,MFRX)*dv
    w(:,:,:,MFRY) = u(:,:,:,MFRY)*dv
    w(:,:,:,MFRZ) = u(:,:,:,MFRZ)*dv
#endif

#ifdef M1CLOSER_FUV_TRANSFER !FUV photons
    w(:,:,:,MEF)   = u(:,:,:,MEF)  *dv
    w(:,:,:,MFRFX) = u(:,:,:,MFRFX)*dv
    w(:,:,:,MFRFY) = u(:,:,:,MFRFY)*dv
    w(:,:,:,MFRFZ) = u(:,:,:,MFRFZ)*dv

  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    w(:,:,:,MECOF)     = u(:,:,:,MECOF)  *dv
    w(:,:,:,MFCORFX)   = u(:,:,:,MFCORFX)*dv
    w(:,:,:,MFCORFY)   = u(:,:,:,MFCORFY)*dv
    w(:,:,:,MFCORFZ)   = u(:,:,:,MFCORFZ)*dv
    w(:,:,:,MEDUSTF)   = u(:,:,:,MEDUSTF)  *dv
    w(:,:,:,MFDUSTRFX) = u(:,:,:,MFDUSTRFX)*dv
    w(:,:,:,MFDUSTRFY) = u(:,:,:,MFDUSTRFY)*dv
    w(:,:,:,MFDUSTRFZ) = u(:,:,:,MFDUSTRFZ)*dv
  #endif

#endif

#ifdef M1CLOSER_IR_TRANSFER ! IR photons
    w(:,:,:,MEIR)  = u(:,:,:,MEIR) *dv
    w(:,:,:,MFIRX) = u(:,:,:,MFIRX)*dv
    w(:,:,:,MFIRY) = u(:,:,:,MFIRY)*dv
    w(:,:,:,MFIRZ) = u(:,:,:,MFIRZ)*dv
#endif


  end subroutine u2w

  !-----------------------------------------------------------------------
  ! MVX, MVY, MVZ, MPのみ変換する
  !-----------------------------------------------------------------------
  subroutine u2w_4(u,w,dv)
#define M_MRHO  0
#define M_MVX   1
#define M_MVY   2
#define M_MVZ   3
#define M_MP    4

    ! use grid, only : Mmin
    use grid, only : Mmin
    use unit
    use parameter
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART    
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv
    real(kind=DBL_KIND) :: tap

#if MODEL_ART > 0
!    real(kind=DBL_KIND),dimension(ARRAYSIZE3(u)) :: gam
    real(kind=DBL_KIND) :: T_K, xmu, gam
    integer :: i,j,k

#endif        

    tap = MP_mu

    w(:,:,:,M_MRHO)=u(:,:,:,MRHO)*dv
    w(:,:,:,M_MVX) =u(:,:,:,MVX) *w(:,:,:,M_MRHO)
    w(:,:,:,M_MVY) =u(:,:,:,MVY) *w(:,:,:,M_MRHO)
    w(:,:,:,M_MVZ) =u(:,:,:,MVZ) *w(:,:,:,M_MRHO)

#if MODEL_ART > 0
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) !+MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
             T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb      !真 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax)                            !c_H2(T_K)真真真真真真真真真真真真?
             gam = 1.d0+tap &                            !真? (Hm, H2p真?)
                  /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
             w(i,j,k,M_MP) =dv*(u(i,j,k,MP)/(gam-1.d0) &
                  +u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2)/2.d0 &
                  +(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*pi8i)
          end do
       end do
    end do
#else !MODEL_ART
    w(:,:,:,M_MP) =dv*(u(:,:,:,MP)/(Gamma-1.d0) &
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0 &
         +(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i)
#endif !MODEL_ART

  end subroutine u2w_4

  ! convert conservative variables w to primitive variable u
  !-----------------------------------------------------------------------
  subroutine w2u(w,u,dv)
    use parameter
    use grid, only : Mmin
    use unit
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

#if MODEL_ART > 0
!    real(kind=DBL_KIND),dimension(ARRAYSIZE3(u)) :: gam
    real(kind=DBL_KIND) :: T_K, xmu, gam, p_gam_1, gam0 
    ! real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
    integer :: i,j,k
    integer :: nske = 0.d0
    integer :: itr, n_itr=30
    real(kind=DBL_KIND) :: err_gam, max_err_gam=1d-10 !gamma真?1e-6真真真?
    ! integer :: count_output=0, max_count=100 !KS DEBUG    
    real(kind=DBL_KIND) :: f1, f2, gam1,dif
    integer :: notConverge
#endif
    integer :: ichem

    u(:,:,:,MRHO) =w(:,:,:,MRHO)/dv
    u(:,:,:,MVX)  =w(:,:,:,MVX) /w(:,:,:,MRHO)
    u(:,:,:,MVY)  =w(:,:,:,MVY) /w(:,:,:,MRHO)
    u(:,:,:,MVZ)  =w(:,:,:,MVZ) /w(:,:,:,MRHO)
    u(:,:,:,MBX)  =w(:,:,:,MBX)/dv
    u(:,:,:,MBY)  =w(:,:,:,MBY)/dv
    u(:,:,:,MBZ)  =w(:,:,:,MBZ)/dv
    u(:,:,:,MDB)  = w(:,:,:,MDB)/dv

#ifdef DM_POTENTIAL
    u(:,:,:,MDMRHO)= w(:,:,:,MDMRHO)/dv 
#endif
    !ART CHEM
#if MODEL_ART > 0
    do ichem = NCEHM_MIN, NCEHM_MAX
      u(:,:,:,ichem)=w(:,:,:,ichem)/w(:,:,:,MRHO)
    enddo

  #ifdef CHEM_MODEL_HF2020
    u(:,:,:,MCO)=w(:,:,:,MCO)/w(:,:,:,MRHO)
  #endif

    u(:,:,:,MKPI)=w(:,:,:,MKPI) / dv
    u(:,:,:,MHPI)=w(:,:,:,MHPI) / dv
    u(:,:,:,MXPI)=w(:,:,:,MXPI) / w(:,:,:,MRHO)
    u(:,:,:,MYPI)=w(:,:,:,MYPI) / w(:,:,:,MRHO)
    u(:,:,:,MZPI)=w(:,:,:,MZPI) / w(:,:,:,MRHO)

  #ifdef METAL
    u(:,:,:,MTD )=w(:,:,:,MTD ) / dv
  #endif

  #ifdef DUST_NOTCONSTANT
    u(:,:,:,MDRHO)=w(:,:,:,MDRHO)/w(:,:,:,MRHO)
  #endif
  #ifdef METAL_TRANSFER
    u(:,:,:,MMET)=w(:,:,:,MMET)/w(:,:,:,MRHO)
  #endif

#endif !MODEL_ART
#ifdef RADTR_DIRECT
  #if MODEL_ART == 2
    u(:,:,:,MKPD)=w(:,:,:,MKPD) / dv
    #ifdef METAL
      u(:,:,:,MGFUV)= w(:,:,:,MGFUV) / dv
      u(:,:,:,MDPH) = w(:,:,:,MDPH)  / dv
      u(:,:,:,MDPCO) = w(:,:,:,MDPCO)  / dv
      u(:,:,:,MKPOII) = w(:,:,:,MKPOII)  / dv
    #endif
  #endif !MODEL_ART
#endif
#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
    u(:,:,:,MGX) = w(:,:,:,MGX) / w(:,:,:,MRHO)
    u(:,:,:,MGY) = w(:,:,:,MGY) / w(:,:,:,MRHO)
    u(:,:,:,MGZ) = w(:,:,:,MGZ) / w(:,:,:,MRHO)
#endif !MPSI

#if MODEL_ART > 0
    !まずgammaを求める

    ! initialize ----------
    gam = 5d0/3d0                             ! initial guess of gamma
    ! ---------------------

    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             p_gam_1 = w(i,j,k,MP)/dv &
                  - 0.5d0*u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2) & 
                   -(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*pi8i 
             xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) ! +MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
             u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

             !refinement真真真真真真真真真?w2u真真真真真真真真真?
             !真真真真真真? => gam真?iteration真真?
             if ( .not. (u(i,j,k,MRHO) > 0d0 .and. u(i,j,k,MP) > 0d0 .and. &
                  u(i,j,k,MHN)  > 0d0 .and. u(i,j,k,MH2) > 0d0 .and. &
                  u(i,j,k,MEL)  > 0d0 .and. u(i,j,k,MHP) > 0d0  .and. &
                  u(i,j,k,MHN)  < 2d0 .and. u(i,j,k,MH2) < 1d0 .and. &
                  u(i,j,k,MEL)  < 2d0 .and. u(i,j,k,MHP) < 2d0 )) cycle
             
             ! --------------------------------------------- Newton-Raphson ------------------------------------------------------
             ! ---------------
             notConverge = 1
             ! ---------------

             do  itr = 0, n_itr-1 
                gam0 = gam
                T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb     !温度 [K]
                gam = 1.d0+(1.d0+4.d0*yHe)/(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
                f1 = gam - gam0

                err_gam = abs(f1)/max(gam,gam0)

                if (err_gam < max_err_gam) then
                  notConverge = 0
                  exit
                endif

                gam = gam0*(1.d0+1d-5)
                gam1= gam
                u(i,j,k,MP) = p_gam_1*(gam-1)
                T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb
                gam = 1.d0+(1.d0+4.d0*yHe)/(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
                f2  = gam - gam1
                dif = (f2-f1)/(gam0*1d-5)
                gam = gam0 - f1/dif
                u(i,j,k,MP) = p_gam_1*(gam-1)

                if (itr == n_itr-1) then
                   !print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam0, err_gam
                   exit
                end if
             enddo
             ! ------------------------------------------------------------------------------------------------------------------

             ! ------------------------------------- 地道な方法 --------------------------------------!
             if(notConverge == 1) then
                gam = 5d0/3d0                             ! initial guess of gamma
                xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) !+MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
                u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

                do itr = 0, n_itr-1   !gammaが収束するまでiteration
                   T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb     !温度 [K]
                   !T_K = min(max(T_K,MP_Tmin),MP_Tmax)                                      !c_H2(T_K)がおかしくならないよう上下の温度フロアを課しておく
                   gam0 = gam
                   gam = 1.d0+(1.d0+4.d0*yHe) &                                             !比熱比 (Hm, H2pは無視)
                        /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2_2(T_K)*u(i,j,k,MH2)))
                   err_gam = abs(gam-gam0)/max(gam,gam0)
                   u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

                   ! if (i==5 .and. j==5 .and. k==5)  then
                   !    print *, "KS DEBUG, gam(5,5,5)", itr, gam0,gam,err_gam, u(i,j,k,MHN),u(i,j,k,MH2),u(i,j,k,MHP) ,xmu, T_K
                   ! end if

                   if (err_gam < max_err_gam) exit ! gammaが収束していたらiterationから抜ける
                   ! if (err_gam /= err_gam) exit ! gammaがNaNなら意味ないので抜ける (すでに非物理的な値の場合はチェック済み)
                   if (itr == n_itr-1) then
                      print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam0, err_gam

                      ! initialize ----------
                      gam = 5d0/3d0                             ! initial guess of gamma
                      ! ---------------------

                      p_gam_1 = (cgs_kb*MP_Tmin)*(u(i,j,k,MRHO)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                      u(i,j,k,MP) = p_gam_1*(gam-1)  

                   end if
                end do
             endif
             ! ------------------------------------------------------------------------------------------------------------------


          end do
       end do
    end do


    
#else !MODEL_ART
    u(:,:,:,MP) = ( w(:,:,:,MP)/dv -0.5d0*u(:,:,:,MRHO) &
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) &
         -(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i &
         )*(Gamma-1.0d0)
#endif



#ifdef M1CLOSER_EUV_TRANSFER
    u(:,:,:,MER)  =  w(:,:,:,MER)  / dv 
    u(:,:,:,MFRX) =  w(:,:,:,MFRX) / dv
    u(:,:,:,MFRY) =  w(:,:,:,MFRY) / dv
    u(:,:,:,MFRZ) =  w(:,:,:,MFRZ) / dv
#endif

#ifdef M1CLOSER_FUV_TRANSFER
    u(:,:,:,MEF)   = w(:,:,:,MEF)  / dv
    u(:,:,:,MFRFX) = w(:,:,:,MFRFX)/ dv
    u(:,:,:,MFRFY) = w(:,:,:,MFRFY)/ dv
    u(:,:,:,MFRFZ) = w(:,:,:,MFRFZ)/ dv

  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    u(:,:,:,MECOF)     = w(:,:,:,MECOF)  / dv
    u(:,:,:,MFCORFX)   = w(:,:,:,MFCORFX)/ dv
    u(:,:,:,MFCORFY)   = w(:,:,:,MFCORFY)/ dv
    u(:,:,:,MFCORFZ)   = w(:,:,:,MFCORFZ)/ dv
    u(:,:,:,MEDUSTF)   = w(:,:,:,MEDUSTF)  / dv
    u(:,:,:,MFDUSTRFX) = w(:,:,:,MFDUSTRFX)/ dv
    u(:,:,:,MFDUSTRFY) = w(:,:,:,MFDUSTRFY)/ dv
    u(:,:,:,MFDUSTRFZ) = w(:,:,:,MFDUSTRFZ)/ dv
  #endif

#endif

#ifdef M1CLOSER_IR_TRANSFER
    u(:,:,:,MEIR)  =  w(:,:,:,MEIR)  / dv 
    u(:,:,:,MFIRX) =  w(:,:,:,MFIRX) / dv
    u(:,:,:,MFIRY) =  w(:,:,:,MFIRY) / dv
    u(:,:,:,MFIRZ) =  w(:,:,:,MFIRZ) / dv
#endif
  end subroutine w2u





  subroutine w2u_4(w,u,dv)
    ! use grid, only : Mmin
    use parameter
    use grid, only : Mmin
    use unit
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

#if MODEL_ART > 0
!    real(kind=DBL_KIND),dimension(ARRAYSIZE3(u)) :: gam
    real(kind=DBL_KIND) :: T_K, xmu, gam, p_gam_1, gam0 
    ! real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
    integer :: i,j,k
    integer :: itr, n_itr=30
    real(kind=DBL_KIND) :: err_gam, max_err_gam=1d-10 !gammaにして1e-6の誤差まで許す
    ! integer :: count_output=0, max_count=100 !KS DEBUG    
    real(kind=DBL_KIND) :: f1, f2, gam1,dif
    integer :: notConverge
#endif
    real(kind=DBL_KIND) :: tap

    tap = MP_mu


#ifdef OUTFLOW_ON
    u(:,:,:,MRHO) =w(:,:,:,M_MRHO)/dv
#endif

    u(:,:,:,MVX) =w(:,:,:,M_MVX)/w(:,:,:,M_MRHO)
    u(:,:,:,MVY) =w(:,:,:,M_MVY)/w(:,:,:,M_MRHO)
    u(:,:,:,MVZ) =w(:,:,:,M_MVZ)/w(:,:,:,M_MRHO)

#if MODEL_ART > 0
    !まずgammaを求める

    ! initialize ----------
    gam = 5d0/3d0                             ! initial guess of gamma
    ! ---------------------

    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             p_gam_1 = w(i,j,k,M_MP)/dv &
                  - 0.5d0*u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2) & 
                   -(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*pi8i 
             xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) ! +MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
             u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

             !refinement真真真真真真真真真?w2u真真真真真真真真真?
             !真真真真真真? => gam真?iteration真真?
             if ( .not. (u(i,j,k,MRHO) > 0d0 .and. u(i,j,k,MP) > 0d0 .and. &
                  u(i,j,k,MHN)  > 0d0 .and. u(i,j,k,MH2) > 0d0 .and. &
                  u(i,j,k,MEL)  > 0d0 .and. u(i,j,k,MHP) > 0d0  .and. &
                  u(i,j,k,MHN)  < 2d0 .and. u(i,j,k,MH2) < 1d0 .and. &
                  u(i,j,k,MEL)  < 2d0 .and. u(i,j,k,MHP) < 2d0 )) cycle
             
             ! --------------------------------------------- Newton-Raphson ------------------------------------------------------
             ! ---------------
             notConverge = 1
             ! ---------------

             do  itr = 0, n_itr-1 
                gam0 = gam
                T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb     !温度 [K]
                gam = 1.d0+tap/(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
                f1 = gam - gam0

                err_gam = abs(f1)/max(gam,gam0)

                if (err_gam < max_err_gam) then
                  notConverge = 0
                  exit
                endif

                gam = gam0*(1.d0+1d-5)
                gam1= gam
                u(i,j,k,MP) = p_gam_1*(gam-1)
                T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb
                gam = 1.d0+tap/(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
                f2  = gam - gam1
                dif = (f2-f1)/(gam0*1d-5)
                gam = gam0 - f1/dif
                u(i,j,k,MP) = p_gam_1*(gam-1)

                if (itr == n_itr-1) then
                   !print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam0, err_gam
                   exit
                end if
             enddo
             ! ------------------------------------------------------------------------------------------------------------------

             ! ------------------------------------- 地道な方法 --------------------------------------!
             if(notConverge == 1) then
                gam = 5d0/3d0                             ! initial guess of gamma
                xmu = MP_mu/ &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) ! +MP_frac_COsum-u(i,j,k,MCO))  !Hm, H2p 真真真真真
                u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

                do itr = 0, n_itr-1   !gammaが収束するまでiteration
                   T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb     !温度 [K]
                   !T_K = min(max(T_K,MP_Tmin),MP_Tmax)                                      !c_H2(T_K)がおかしくならないよう上下の温度フロアを課しておく
                   gam0 = gam
                   gam = 1.d0+tap &                                             !比熱比 (Hm, H2pは無視)
                        /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2_2(T_K)*u(i,j,k,MH2)))
                   err_gam = abs(gam-gam0)/max(gam,gam0)
                   u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

                   ! if (i==5 .and. j==5 .and. k==5)  then
                   !    print *, "KS DEBUG, gam(5,5,5)", itr, gam0,gam,err_gam, u(i,j,k,MHN),u(i,j,k,MH2),u(i,j,k,MHP) ,xmu, T_K
                   ! end if

                   if (err_gam < max_err_gam) exit ! gammaが収束していたらiterationから抜ける
                   ! if (err_gam /= err_gam) exit ! gammaがNaNなら意味ないので抜ける (すでに非物理的な値の場合はチェック済み)
                   if (itr == n_itr-1) then
                      print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam0, err_gam

                      ! initialize ----------
                      gam = 5d0/3d0                             ! initial guess of gamma
                      ! ---------------------
                      p_gam_1 = (cgs_kb*MP_Tmin)*(u(i,j,k,MRHO)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                      u(i,j,k,MP) = p_gam_1*(gam-1)  

                   end if
                end do
             endif
             ! ------------------------------------------------------------------------------------------------------------------


          end do
       end do
    end do


    
#else !MODEL_ART
    u(:,:,:,MP) = ( w(:,:,:,M_MP)/dv -0.5d0*u(:,:,:,MRHO) &
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) &
         -(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i &
         )*(Gamma-1.0d0)
#endif

  end subroutine w2u_4
  
#undef M_MRHO  
#undef M_MVX   
#undef M_MVY   
#undef M_MVZ   
#undef M_MP    







#if MODEL_ART > 0
  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable w
  ! gamma真真真真真真真?refinement真真真  
  !-----------------------------------------------------------------------
  subroutine u2w_withgam(u,w,dv,gam)
    use parameter
    ! use grid, only : Mmin
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
    use unit
    use primordial
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
    real(kind=DBL_KIND),dimension(:,:,:),intent(IN) :: gam
    real(kind=DBL_KIND),intent(IN) :: dv

    real(kind=DBL_KIND) :: T_K, xmu
    integer :: i,j,k, ichem

    w(:,:,:,MRHO)=u(:,:,:,MRHO)*dv
    w(:,:,:,MVX)=u(:,:,:,MVX)*w(:,:,:,MRHO)
    w(:,:,:,MVY)=u(:,:,:,MVY)*w(:,:,:,MRHO)
    w(:,:,:,MVZ)=u(:,:,:,MVZ)*w(:,:,:,MRHO)
    w(:,:,:,MBX)= u(:,:,:,MBX)*dv
    w(:,:,:,MBY)= u(:,:,:,MBY)*dv
    w(:,:,:,MBZ)= u(:,:,:,MBZ)*dv
    w(:,:,:,MDB) = u(:,:,:,MDB)*dv

#ifdef DM_POTENTIAL
    w(:,:,:,MDMRHO)=u(:,:,:,MDMRHO)*dv 
#endif

    !ART CHEM
    do ichem = NCEHM_MIN, NCEHM_MAX
      w(:,:,:,ichem)=u(:,:,:,ichem)*w(:,:,:,MRHO)
    enddo

  #ifdef CHEM_MODEL_HF2020
    w(:,:,:,MCO)=u(:,:,:,MCO)*w(:,:,:,MRHO)
  #endif

    w(:,:,:,MKPI)=u(:,:,:,MKPI) * dv
    w(:,:,:,MHPI)=u(:,:,:,MHPI) * dv
    w(:,:,:,MXPI)=u(:,:,:,MXPI) * w(:,:,:,MRHO)
    w(:,:,:,MYPI)=u(:,:,:,MYPI) * w(:,:,:,MRHO)
    w(:,:,:,MZPI)=u(:,:,:,MZPI) * w(:,:,:,MRHO)

#ifdef METAL
    w(:,:,:,MTD )=u(:,:,:,MTD ) * dv
#endif

#ifdef DUST_NOTCONSTANT
    w(:,:,:,MDRHO)=u(:,:,:,MDRHO)*w(:,:,:,MRHO)
#endif
#ifdef METAL_TRANSFER
    w(:,:,:,MMET)=u(:,:,:,MMET)*w(:,:,:,MRHO)
#endif



#ifdef RADTR_DIRECT
  #if MODEL_ART == 2
    w(:,:,:,MKPD)=u(:,:,:,MKPD) * dv
    #ifdef METAL
      w(:,:,:,MGFUV) = u(:,:,:,MGFUV) * dv
      w(:,:,:,MDPH)  = u(:,:,:,MDPH)  * dv
      w(:,:,:,MDPCO)  = u(:,:,:,MDPCO)  * dv
      w(:,:,:,MKPOII)  = u(:,:,:,MKPOII)  * dv
    #endif
  #endif !MODEL_ART
#endif
#ifdef MPSI
    w(:,:,:,MPSI) = u(:,:,:,MPSI) * dv
    w(:,:,:,MGX) = u(:,:,:,MGX) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGY) = u(:,:,:,MGY) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGZ) = u(:,:,:,MGZ) * u(:,:,:,MRHO) * dv
#endif !MPSI

    w(:,:,:,MP) =dv*(u(:,:,:,MP)/(gam(:,:,:)-1.d0) &
                  +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0 &
                  +(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i)
                 
#ifdef M1CLOSER_EUV_TRANSFER
    w(:,:,:,MER)  = u(:,:,:,MER) *dv
    w(:,:,:,MFRX) = u(:,:,:,MFRX)*dv
    w(:,:,:,MFRY) = u(:,:,:,MFRY)*dv
    w(:,:,:,MFRZ) = u(:,:,:,MFRZ)*dv
#endif

#ifdef M1CLOSER_FUV_TRANSFER !FUV photons
    w(:,:,:,MEF)   = u(:,:,:,MEF)  *dv
    w(:,:,:,MFRFX) = u(:,:,:,MFRFX)*dv
    w(:,:,:,MFRFY) = u(:,:,:,MFRFY)*dv
    w(:,:,:,MFRFZ) = u(:,:,:,MFRFZ)*dv

  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    w(:,:,:,MECOF)     = u(:,:,:,MECOF)  *dv
    w(:,:,:,MFCORFX)   = u(:,:,:,MFCORFX)*dv
    w(:,:,:,MFCORFY)   = u(:,:,:,MFCORFY)*dv
    w(:,:,:,MFCORFZ)   = u(:,:,:,MFCORFZ)*dv
    w(:,:,:,MEDUSTF)   = u(:,:,:,MEDUSTF)  *dv
    w(:,:,:,MFDUSTRFX) = u(:,:,:,MFDUSTRFX)*dv
    w(:,:,:,MFDUSTRFY) = u(:,:,:,MFDUSTRFY)*dv
    w(:,:,:,MFDUSTRFZ) = u(:,:,:,MFDUSTRFZ)*dv
  #endif

#endif

#ifdef M1CLOSER_IR_TRANSFER
    w(:,:,:,MEIR)  = u(:,:,:,MEIR) *dv
    w(:,:,:,MFIRX) = u(:,:,:,MFIRX)*dv
    w(:,:,:,MFIRY) = u(:,:,:,MFIRY)*dv
    w(:,:,:,MFIRZ) = u(:,:,:,MFIRZ)*dv
#endif
    
  end subroutine u2w_withgam
 
  !-----------------------------------------------------------------------
  ! convert conservative variables w to primitive variable u
  ! gamma真真真真真真真真真真真真真真?
  !-----------------------------------------------------------------------
  subroutine w2u_withgam(w,u,dv,gam)
    use parameter
    ! use grid, only : Mmin
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
    use unit

    use primordial
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(OUT)
    real(kind=DBL_KIND),dimension(:,:,:),intent(IN) :: gam
    real(kind=DBL_KIND),intent(IN) :: dv

    real(kind=DBL_KIND) :: T_K, xmu
    integer :: i,j,k,ichem    

    u(:,:,:,MRHO) =w(:,:,:,MRHO)/dv
    u(:,:,:,MVX) =w(:,:,:,MVX)/w(:,:,:,MRHO)
    u(:,:,:,MVY) =w(:,:,:,MVY)/w(:,:,:,MRHO)
    u(:,:,:,MVZ) =w(:,:,:,MVZ)/w(:,:,:,MRHO)
    u(:,:,:,MBX) =w(:,:,:,MBX)/dv
    u(:,:,:,MBY) =w(:,:,:,MBY)/dv
    u(:,:,:,MBZ) =w(:,:,:,MBZ)/dv
    u(:,:,:,MDB) = w(:,:,:,MDB)/dv

#ifdef DM_POTENTIAL
    u(:,:,:,MDMRHO) = w(:,:,:,MDMRHO) / dv 
#endif

    do ichem = NCEHM_MIN, NCEHM_MAX
      u(:,:,:,ichem)=w(:,:,:,ichem)/w(:,:,:,MRHO)
    enddo

#ifdef CHEM_MODEL_HF2020
    u(:,:,:,MCO)=w(:,:,:,MCO)/w(:,:,:,MRHO)
#endif

    u(:,:,:,MKPI)=w(:,:,:,MKPI) / dv
    u(:,:,:,MHPI)=w(:,:,:,MHPI) / dv
    u(:,:,:,MXPI)=w(:,:,:,MXPI) / w(:,:,:,MRHO)
    u(:,:,:,MYPI)=w(:,:,:,MYPI) / w(:,:,:,MRHO)
    u(:,:,:,MZPI)=w(:,:,:,MZPI) / w(:,:,:,MRHO)

#ifdef METAL
    u(:,:,:,MTD )=w(:,:,:,MTD ) / dv
#endif

#ifdef DUST_NOTCONSTANT
    u(:,:,:,MDRHO)=w(:,:,:,MDRHO)/w(:,:,:,MRHO)
#endif
#ifdef METAL_TRANSFER
    u(:,:,:,MMET)=w(:,:,:,MMET)/w(:,:,:,MRHO)
#endif


#ifdef RADTR_DIRECT
  #if MODEL_ART == 2
    u(:,:,:,MKPD)=w(:,:,:,MKPD) / dv
    #ifdef METAL
      u(:,:,:,MGFUV) = w(:,:,:,MGFUV)/  dv
      u(:,:,:,MDPH)  = w(:,:,:,MDPH) /  dv
      u(:,:,:,MDPCO)  = w(:,:,:,MDPCO) /  dv
      u(:,:,:,MKPOII)  = w(:,:,:,MKPOII) /  dv
    #endif
  #endif !MODEL_ART
#endif

#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
    u(:,:,:,MGX) = w(:,:,:,MGX) / w(:,:,:,MRHO)
    u(:,:,:,MGY) = w(:,:,:,MGY) / w(:,:,:,MRHO)
    u(:,:,:,MGZ) = w(:,:,:,MGZ) / w(:,:,:,MRHO)
#endif !MPSI

    u(:,:,:,MP) = ( w(:,:,:,MP)/dv -0.5d0*u(:,:,:,MRHO) &
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)  &
         -(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i &
         )*(gam(:,:,:)-1.d0)

#ifdef M1CLOSER_EUV_TRANSFER
    u(:,:,:,MER)  =  w(:,:,:,MER) / dv 
    u(:,:,:,MFRX) =  w(:,:,:,MFRX)/ dv 
    u(:,:,:,MFRY) =  w(:,:,:,MFRY)/ dv
    u(:,:,:,MFRZ) =  w(:,:,:,MFRZ)/ dv
#endif

#ifdef M1CLOSER_FUV_TRANSFER
    u(:,:,:,MEF)   = w(:,:,:,MEF)  / dv
    u(:,:,:,MFRFX) = w(:,:,:,MFRFX)/ dv
    u(:,:,:,MFRFY) = w(:,:,:,MFRFY)/ dv
    u(:,:,:,MFRFZ) = w(:,:,:,MFRFZ)/ dv

  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    u(:,:,:,MECOF)     = w(:,:,:,MECOF)  / dv
    u(:,:,:,MFCORFX)   = w(:,:,:,MFCORFX)/ dv
    u(:,:,:,MFCORFY)   = w(:,:,:,MFCORFY)/ dv
    u(:,:,:,MFCORFZ)   = w(:,:,:,MFCORFZ)/ dv
    u(:,:,:,MEDUSTF)   = w(:,:,:,MEDUSTF)  / dv
    u(:,:,:,MFDUSTRFX) = w(:,:,:,MFDUSTRFX)/ dv
    u(:,:,:,MFDUSTRFY) = w(:,:,:,MFDUSTRFY)/ dv
    u(:,:,:,MFDUSTRFZ) = w(:,:,:,MFDUSTRFZ)/ dv
  #endif

#endif

#ifdef M1CLOSER_IR_TRANSFER
    u(:,:,:,MEIR)  =  w(:,:,:,MEIR) / dv 
    u(:,:,:,MFIRX) =  w(:,:,:,MFIRX)/ dv 
    u(:,:,:,MFIRY) =  w(:,:,:,MFIRY)/ dv
    u(:,:,:,MFIRZ) =  w(:,:,:,MFIRZ)/ dv
#endif

!---------skeDEBUG------------------------------------------------------    
 !   do i = lbound(u,1),ubound(u,1)
 !      do j = lbound(u,2),ubound(u,2)
 !         do k = lbound(u,3),ubound(u,3)
               ! if(u(i,j,k,MP)<0.d0) then
                !   print*,'----skeDEBUG---w2u_withgam------'
                !   print'(A,(1P4E15.7))',' i,j,k   : ',i,j,k
                !   print'(A,(1P4E15.7))',' n , p   : ',u(i,j,k,MRHO)*Unit_rho/(xmu*cgs_mh),u(i,j,k,MP)
                !   print'(A,(1P4E15.7))','   B     : ',dsqrt(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*dsqrt(Unit_e)
                !   print'(A,(1P4E15.7))','Et,Ek,Em : ',u(i,j,k,MP)/(gam-1)  &
  !                                                    ,u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2)*0.5d0  &
  !                                                    ,(u(i,j,k,MBX)**2+u(i,j,k,MBY)**2+u(i,j,k,MBZ)**2)*pi8i
                !  print*, '---skeENDBUG-------------------' 
                !endif
   !       enddo
   !    enddo
   ! enddo
!--------------skeENDBUG---------------------------------------------
  end subroutine w2u_withgam
#endif !MODEL_ART

  !-----------------------------------------------------------------------
  ! change u to w in one array
  !-----------------------------------------------------------------------
  subroutine conv_u2w(uw, dv)
    real(kind=DBL_KIND),dimension(:,:,:,:) :: uw !(INOUT)
    real(kind=DBL_KIND),intent(IN) :: dv
    real(kind=DBL_KIND),dimension(ARRAYSIZE4(uw)) :: swap
    call u2w(uw,swap,dv)
    uw = swap
  end subroutine conv_u2w
  !-----------------------------------------------------------------------
  ! change w to u in one array
  !-----------------------------------------------------------------------
  subroutine conv_w2u(wu, dv)
    real(kind=DBL_KIND),dimension(:,:,:,:) :: wu !(INOUT)
    real(kind=DBL_KIND),intent(IN) :: dv
    real(kind=DBL_KIND),dimension(ARRAYSIZE4(wu)) :: swap
    call w2u(wu,swap,dv)
    wu = swap
  end subroutine conv_w2u
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
#if MODEL_ART > 0
  subroutine get_gamma(gid, gam)
    use parameter
    use grid, only : Imingh,Jmingh,Kmingh,Imaxgh,Jmaxgh,Kmaxgh, get_Ucomp
    use unit
    use primordial
    
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH),intent(OUT) :: gam

    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp

   ! real(kind=DBL_KIND),dimension(:,:,:),pointer :: bx, by, bz
    real(kind=DBL_KIND) :: xmu, T_K
    integer :: i,j,k

    !chemical abundance
    rho => get_Ucomp(MRHO,gid)
    p => get_Ucomp(MP,gid)
    yhn => get_Ucomp(MHN,gid)
    yh2 => get_Ucomp(MH2,gid)
    yel => get_Ucomp(MEL,gid)
    yhp => get_Ucomp(MHP,gid)

    do i = Imingh, Imaxgh
       do j = Jmingh, Jmaxgh
          do k = Kmingh, Kmaxgh
             xmu = MP_mu /(yhn(i,j,k)+yh2(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe) !+MP_frac_COsum-yco(i,j,k))  !Hm, H2p 真真真真真
             T_K = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb      !真 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax) !真真真真真?             
             gam(i,j,k) = 1.d0+(1.d0+4.d0*yHe) &                            !真? (Hm, H2p真?)
                  /(xmu*(1.5d0*(yhn(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe) + c_H2(T_K)*yh2(i,j,k)))
             ! if (i==5 .and. j==5 .and. k==5)  then
             !    print *, "KS DEBUG3, gid, gam(5,5,5)", gid, gam(i,j,k),yhn(i,j,k),yh2(i,j,k),yel(i,j,k),yhp(i,j,k) ,xmu, T_K
             ! end if
          end do
       end do
    end do
  end subroutine get_gamma
#endif !MODEL_ART  
end module eos
subroutine eos_init
#if MODEL_ART == 0
  use eos, only : Gamma
  use unit, only: ModelParam_gamma
  real(kind=8) :: csd
  Gamma = ModelParam_gamma
#endif !MODEL_ART
end subroutine eos_init

