#include "config.h"

!activating debug output
!#define KS_DEBUG

!-------------------------------------------------------------------------
! module depending on equation of state.
!-------------------------------------------------------------------------
module eos
 
#ifndef RADTR_M1closer
  use modelParameter, only : MP_Tmin, MP_Tmax, MP_mu, MP_frac_COsum ! KS ADDED
#else
  use modelParameter, only : MP_Tmin, MP_Tmax, MP_Ctil_nD, MP_mu, MP_frac_COsum ! HF ADDED
#endif
  implicit none
  ! KS MODIFIED
  ! gamma is given in a consistent way with chemical abundances

  ! Gamma is given by eos_init()
#if MODEL_ART == 0
  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0
  ! real(kind=DBL_KIND),save :: Gamma
#endif

!!$  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0
!!$  real(kind=DBL_KIND),parameter :: Gamma = 1.01d0
contains
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  function cyclecomp(ncrd, invert) result( mcycle )
    use grid , only : Mmin, Mmax
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(Mmin:Mmax) :: mcycle
    integer,dimension(MX:MZ) :: mcycle3 = (/ MVX, MVY, MVZ /)
    integer :: m
    do m = Mmin, Mmax
       mcycle(m) = m
    enddo
    if ( present( invert ) ) then
       mcycle(MVX:MVX+size(mcycle3)-1) = cshift( mcycle3, -ncrd)
    else
       mcycle(MVX:MVX+size(mcycle3)-1) = cshift( mcycle3, ncrd )
    endif
  end function cyclecomp
  !-----------------------------------------------------------------------
  ! get numerical flux in one dimension, RoeM2 (Kim et al. 2003, JCP, 185, 342)
  !-----------------------------------------------------------------------
  ! macro for entropy condition
#define ELMOD(EL,ELBAR,EL1, EL2) \
  eps=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-eps) ;\
  x2=1.d0-x1 ;\
  EL= x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(eps+x1)+eps)

#if MODEL_ART > 0
  !gammaを渡す
  subroutine flux( ql, qr, pratio, f, ncrd, gam )
#else !MODEL_ART
  subroutine flux( ql, qr, pratio, f, ncrd )
#endif !MODEL_ART    
    use grid, only : Imin,Jmin,Kmin,Mmin, Imax,Jmax,Kmax,Mmax, Imingh,Jmingh,Kmingh, Imaxgh,Jmaxgh,Kmaxgh, &
         globdbg_myrank, globdbg_mygid,globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
    use util, only : util_arroffset
    integer,intent(IN) :: ncrd
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,MX:) :: pratio  ! (IN)
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: f   ! (OUT)

    integer :: i,j,k, is,js,ks,ie,je,ke,io,jo,ko
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,hl,el,cl, &
         rhor,ur,vr,wr,pr,hr,er,cr, &
         drho, drhou, drhov, drhow, drhoh, dp, dh, du, dv, dw, &
         sql,sqr,sqa,rhob,ub,vb,wb,hb,qb2,cb2,cb,ub2,vb2,wb2, &
         gm1,mn, b1, b2, b1b2, db12, ff, gg, rbdq, hh, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5
    real(kind=DBL_KIND) :: el4c, a7c, a8c, a9c, w1c
    real(kind=DBL_KIND) :: eps, x1, x2
    integer :: ichem


#if MODEL_ART > 0
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:),intent(IN) :: gam
#endif !MODEL_ART

    f(:,:,:,:) = 0.d0

    call util_arroffset(ncrd,io,jo,ko)
    ks = Kmin-ko
    ke = Kmax
    js = Jmin-jo
    je = Jmax
    is = Imin-io
    ie = Imax
#ifdef EMULATE_2DIM
    ks = Kmin
    ke = Kmin
#endif !EMULATE_2DIM

    do k = ks, ke
       do j = js, je
          do i = is, ie
#if MODEL_ART > 0
             gm1 = gam(i,j,k) - 1
#else !MODEL_ART
             gm1 = Gamma - 1
#endif !MODEL_ART
             ! --------------
             ! left variables
             ! --------------
             rhol = abs(ql(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
             ul = ql(i,j,k,MVX)
             vl = ql(i,j,k,MVY)
             wl = ql(i,j,k,MVZ)
             pl = abs(ql(i,j,k,MP))
             el = rhol * (ul**2 + vl**2 + wl**2) / 2 + pl / gm1
             hl = ( el + pl )/rhol
#if MODEL_ART > 0
             cl = sqrt( gam(i,j,k) * pl / rhol )
#else !MODEL_ART
             cl = sqrt( Gamma * pl / rhol )
#endif !MODEL_ART
             ! --------------
             ! right variables
             ! --------------
             rhor = abs(qr(i,j,k,MRHO))
             ur = qr(i,j,k,MVX)
             vr = qr(i,j,k,MVY)
             wr = qr(i,j,k,MVZ)
             pr = abs(qr(i,j,k,MP))
             er = rhor * (ur**2 + vr**2 + wr**2) / 2 + pr / gm1
             hr = ( er + pr )/rhor
#if MODEL_ART > 0
             cr = sqrt( gam(i,j,k) * pr / rhor )
#else !MODEL_ART
             cr = sqrt( Gamma * pr / rhor )
#endif!MODEL_ART
             ! -------
             ! delta Q
             ! -------
             drho = rhor - rhol
             drhou = rhor * ur - rhol * ul
             drhov = rhor * vr - rhol * vl
             drhow = rhor * wr - rhol * wl
             drhoh = rhor * hr - rhol * hl
             du = ur - ul
             dv = vr - vl
             dw = wr - wl
             dp = pr - pl
             dh = hr - hl
             ! ----------------------
             ! intermidiate variables
             ! ----------------------
             sql = sqrt(rhol)
             sqr = sqrt(rhor)
             sqa = sql + sqr
             rhob = sql * sqr
             ub = (sql * ul + sqr * ur) / sqa
             vb = (sql * vl + sqr * vr) / sqa
             wb = (sql * wl + sqr * wr) / sqa
             hb = (sql * hl + sqr * hr) / sqa
             ub2 = ub**2
             vb2 = vb**2
             wb2 = wb**2
             qb2 = ub2 + vb2 + wb2
             cb2 = gm1*(hb - qb2/2)
             cb  =  sqrt( cb2 )
             mn = abs(ub/cb)         !Mach number

             !
             b1 = max(0.d0, ub+cb, ur+cb)
             b2 = min(0.d0, ub-cb, ul-cb)
             b1b2 = b1 * b2
             db12 = b1 - b2

             if (ncrd == MX) then
                hh = 1.d0 - min(pratio(i,j,k,MX), &
                     pratio(i  ,j-1,k  ,MY), pratio(i+1,j-1,k  ,MY), &
                     pratio(i  ,j  ,k  ,MY), pratio(i+1,j  ,k  ,MY), &
                     pratio(i  ,j  ,k-1,MZ), pratio(i+1,j  ,k-1,MZ), &
                     pratio(i  ,j  ,k  ,MZ), pratio(i+1,j  ,k  ,MZ))
             else if (ncrd == MY) then
                hh = 1.d0 - min(pratio(i,j,k,MY), &
                     pratio(i  ,j  ,k-1,MZ), pratio(i  ,j+1,k-1,MZ), &
                     pratio(i  ,j  ,k  ,MZ), pratio(i  ,j+1,k  ,MZ), &
                     pratio(i-1,j  ,k  ,MX), pratio(i-1,j+1,k  ,MX), &
                     pratio(i  ,j  ,k  ,MX), pratio(i  ,j+1,k  ,MX))
             else
                hh = 1.d0 - min(pratio(i,j,k,MZ), &
                     pratio(i-1,j  ,k  ,MX), pratio(i-1,j  ,k+1,MX), &
                     pratio(i  ,j  ,k  ,MX), pratio(i  ,j  ,k+1,MX), &
                     pratio(i  ,j-1,k  ,MY), pratio(i  ,j-1,k+1,MY), &
                     pratio(i  ,j  ,k  ,MY), pratio(i  ,j  ,k+1,MY))
             endif
             if (qb2 == 0.d0 ) then
                ff = 1.d0
             else
                ff = mn**hh
             endif
             if (mn == 0.d0) then
                gg = 1.d0
             else
                gg = mn**(1.d0-pratio(i,j,k,ncrd))
             endif
             rbdq = drho - ff * dp / cb2
             ! ---------
             ! FL and FR
             ! ---------
             fl1 = rhol * ul
             fl2 = fl1 * ul + pl
             fl3 = fl1 * vl
             fl4 = fl1 * wl
             fl5 = ul * ( el + pl )

             fr1 = rhor * ur
             fr2 = fr1 * ur + pr
             fr3 = fr1 * vr
             fr4 = fr1 * wr
             fr5 = ur * ( er + pr )
             ! -----------------
             ! intermidiate flux
             ! -----------------
             f(i,j,k,MRHO) = (b1*fl1 - b2*fr1)/db12 &
                  + b1b2/db12 * drho &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq
             f(i,j,k,MVX)  = (b1*fl2 - b2*fr2)/db12 &
                  + b1b2/db12 * drhou &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq * ub
             f(i,j,k,MVY)  = (b1*fl3 - b2*fr3)/db12 &
                  + b1b2/db12 * drhov &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * vb + rhob * dv)
             f(i,j,k,MVZ)  = (b1*fl4 - b2*fr4)/db12 &
                  + b1b2/db12 * drhow &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * wb + rhob * dw)
             f(i,j,k,MP)  = (b1*fl5 - b2*fr5)/db12 &
                  + b1b2/db12 * drhoh &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * hb + rhob * dh)

             ! Up-wind flux of chemical abundances (ART CHEM)             
#if MODEL_ART > 0
             do ichem = NCEHM_MIN, NCEHM_MAX
                f(i,j,k,ichem) = 0.5d0 * ub * (ql(i,j,k,ichem)*rhol + qr(i,j,k,ichem)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,ichem)*rhor - ql(i,j,k,ichem)*rhol)
             enddo
  #ifdef CHEM_MODEL_HF2020
             f(i,j,k,MCO) = 0.5d0 * ub * (ql(i,j,k,MCO)*rhol + qr(i,j,k,MCO)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,MCO)*rhor - ql(i,j,k,MCO)*rhol)
  #endif
#endif !MODEL_ART


#ifdef DUST_NOTCONSTANT
             f(i,j,k,MDRHO) = 0.5d0 * ub * (ql(i,j,k,MDRHO)*rhol + qr(i,j,k,MDRHO)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,MDRHO)*rhor - ql(i,j,k,MDRHO)*rhol)
#endif
#ifdef METAL_TRANSFER
             f(i,j,k,MMET) = 0.5d0 * ub * (ql(i,j,k,MMET)*rhol + qr(i,j,k,MMET)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,MMET)*rhor - ql(i,j,k,MMET)*rhol)
#endif

             !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
             if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid .and. &
                  i==globdbg_i+lbound(f,1) .and. j==globdbg_j+lbound(f,2) .and.k==globdbg_k+lbound(f,3)) then
                print *, "(KS DEBUG) f: ",i-lbound(f,1),j-lbound(f,2),k-lbound(f,3),&
                f(i,j,k,MRHO), b1,fl1,b2,fr1,db12, &
                   b1b2,db12, drho, &
                   gg, b1b2,db12, 1.d0, mn, rbdq
                print *, "(KS DEBUG) rbdq: ",&
                     rbdq,drho,ff,dp,cb2
                print *, "(KS DEBUG) ff: ",ff, mn,hh
             end if
#endif !KS_DEBUG
             !------------------- KS DEBUG ---------------------!
          enddo
       enddo
    enddo
#ifdef EMULATE_2DIM
    do k = Kmin-1, Kmax
       if (k /= Kmin) f(:,:,k,:) = f(:,:,Kmin,:)
    enddo
#endif !EMULATE_2DIM
    !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
       i=globdbg_i+lbound(f,1)
       j=globdbg_j+lbound(f,2)
       k=globdbg_k+lbound(f,3)
       print '(A,6I8,/,1P6E15.7)', "(flux, KS DEBUG)",&
            globdbg_myrank, globdbg_mygid, i-lbound(f,1), j-lbound(f,2), k-lbound(f,3),ncrd,&
            ql(i,j,k,MRHO),qr(i,j,k,MRHO),f(i,j,k,MRHO),f(i-1,j,k,MRHO),f(i,j-1,k,MRHO), f(i,j,k-1,MRHO)
       print '(1P6E15.7)',ql(i,j,k,MRHO),qr(i,j,k,MRHO),ql(i,j,k,MP),qr(i,j,k,MP),ql(i,j,k,MVX),qr(i,j,k,MVX)
       print '(1P7E15.7)',ql(i,j,k,MP),qr(i,j,k,MP),gam(i,j,k), f(i,j,k,MRHO),f(i-1,j,k,MRHO),f(i,j-1,k,MRHO), f(i,j,k-1,MRHO)
    end if
#endif !KS_DEBUG
    !------------------- KS DEBUG ---------------------!       

  end subroutine flux
  !-----------------------------------------------------------------------
  ! find dt according to CFL condtion
  !-----------------------------------------------------------------------
  function get_dt_by_cflcond( id ) result( dt )
    use grid
    integer,intent(IN) :: id
    real(kind=DBL_KIND) :: dt
    integer :: i,j,k,level
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, p, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH) :: cs
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
    call get_gamma(id, gam)
    cs = sqrt(abs(gam*p/rho))
#else !MODEL_ART
    cs = sqrt(abs(Gamma*p/rho))
#endif !MODEL_ART
#ifdef WITH_SELFGRAVITY
#define DTI_(V_,G_,H_) sqrt((V_)**2+2.d0*(G_)*(H_)) / (H_)
    dt = (CFL) / maxval( &
         DTI_( abs(vx)+cs, abs(gx), h(MX) ) + &
         DTI_( abs(vy)+cs, abs(gy), h(MY) ) + &
         DTI_( abs(vz)+cs, abs(gz), h(MZ) ), mask=GridMask )
#else   !WITH_SELFGRAVITY
    dt = (CFL) / maxval( &
         (abs(vx)+cs)/h(MX) + &
         (abs(vy)+cs)/h(MY) + &
         (abs(vz)+cs)/h(MZ), mask=GridMask )
#endif  !WITH_SELFGRAVITY

#ifdef RADTR_M1closer
    dt_radtr = minval(h)*CFLfac_radtr/ (3.d0*MP_Ctil_nD)*MAX_STEP_RADTR_TO_HYDRO
    dt = min(dt, dt_radtr)
#endif

  end function get_dt_by_cflcond

  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable w
  ! gammaが与えられていない場合（例えばrefinementの際に利用）  
  !-----------------------------------------------------------------------
  subroutine u2w(u,w,dv)
    ! use grid, only : Mmin
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
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
    w(:,:,:,MZPI)=u(:,:,:,MZPI)*w(:,:,:,MRHO)
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
             xmu = MP_mu / &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
             T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb      !温度 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax)                            !c_H2(T_K)がおかしくならないよう上下の温度フロアを課しておく
             gam = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
                  /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
             w(i,j,k,MP) =dv*(u(i,j,k,MP)/(gam-1.d0) &
                  +u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2)/2.d0)
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
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0)
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


    
    !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
       i=globdbg_i+lbound(w,1)
       j=globdbg_j+lbound(w,2)
       k=globdbg_k+lbound(w,3)
       print '(A,5I8,/,1P6E15.7,/,1P6E15.7)', "(u2w, KS DEBUG)",&
            globdbg_myrank, globdbg_mygid, i-lbound(w,1), j-lbound(w,2), k-lbound(w,3),&
            u(i,j,k,MRHO),u(i,j,k,MP),u(i,j,k,MVX),u(i,j,k,MVY),u(i,j,k,MVZ),u(i,j,k,MH2), &
            w(i,j,k,MRHO),w(i,j,k,MP),w(i,j,k,MVX),w(i,j,k,MVY),w(i,j,k,MVZ),w(i,j,k,MH2)
    end if
#endif !KS DEBUG
    !------------------- KS DEBUG ---------------------!       
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
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
    use unit
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
             xmu = MP_mu / &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
             T_K = u(i,j,k,MP)*Unit_e*cgs_amu*xmu /(u(i,j,k,MRHO)*Unit_rho)/cgs_kb      !温度 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax)                            !c_H2(T_K)がおかしくならないよう上下の温度フロアを課しておく
             gam = 1.d0+tap &                            !比熱比 (Hm, H2pは無視)
                  /(xmu*(1.5d0*(u(i,j,k,MHN)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe) + c_H2(T_K)*u(i,j,k,MH2)))
             w(i,j,k,M_MP) =dv*(u(i,j,k,MP)/(gam-1.d0) &
                  +u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2)/2.d0)
          end do
       end do
    end do
#else !MODEL_ART
    w(:,:,:,M_MP) =dv*(u(:,:,:,MP)/(Gamma-1.d0) &
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0)
#endif !MODEL_ART

  end subroutine u2w_4

  !-----------------------------------------------------------------------
  ! convert conservative variables w to primitive variable u
  ! gammaが与えられていない場合（例えばrefinementの際に利用）
  ! ときどきunphysicalな値が渡される場合があるみたいなので（ゴーストセルとかで値を使わない？？）、あまり固執しすぎないように対応
  !-----------------------------------------------------------------------
  subroutine w2u(w,u,dv)
    ! use grid, only : Mmin
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
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
    integer :: ichem

    u(:,:,:,MRHO) =w(:,:,:,MRHO)/dv
    u(:,:,:,MVX)  =w(:,:,:,MVX) /w(:,:,:,MRHO)
    u(:,:,:,MVY)  =w(:,:,:,MVY) /w(:,:,:,MRHO)
    u(:,:,:,MVZ)  =w(:,:,:,MVZ) /w(:,:,:,MRHO)

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
                  - 0.5d0*u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2) !(gam-1)をかけてpになる量
             xmu = MP_mu / &
                  (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
             u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

             !refinementのときなどに、後で使わないセルについてw2uが呼ばれて変な値が入っていることがある
             !明らかに変な値を見つけたら => gam探査のiterationをスキップ
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
                xmu = MP_mu / &
                    (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
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
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) )*(Gamma-1.0d0)
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
    

    !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
       i=globdbg_i+lbound(w,1)
       j=globdbg_j+lbound(w,2)
       k=globdbg_k+lbound(w,3)
       print '(A,5I8,/,1P6E15.7,/,1P6E15.7)', "(w2u, KS DEBUG)" , &
            globdbg_myrank, globdbg_mygid, i-lbound(w,1), j-lbound(w,2), k-lbound(w,3),&
            w(i,j,k,MRHO),w(i,j,k,MP),w(i,j,k,MVX),w(i,j,k,MVY),w(i,j,k,MVZ),w(i,j,k,MH2), &
            u(i,j,k,MRHO),u(i,j,k,MP),u(i,j,k,MVX),u(i,j,k,MVY),u(i,j,k,MVZ),u(i,j,k,MH2)
    end if
#endif !KS_DEBUG
    !------------------- KS DEBUG ---------------------!       
  end subroutine w2u





  subroutine w2u_4(w,u,dv)
    ! use grid, only : Mmin
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k ! KS DEBUG
    use unit
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
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
                  - 0.5d0*u(i,j,k,MRHO)*(u(i,j,k,MVX)**2+u(i,j,k,MVY)**2+u(i,j,k,MVZ)**2) !(gam-1)をかけてpになる量
             xmu = MP_mu / &
                    (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
             u(i,j,k,MP) = p_gam_1*(gam-1)               ! pressure

             !refinementのときなどに、後で使わないセルについてw2uが呼ばれて変な値が入っていることがある
             !明らかに変な値を見つけたら => gam探査のiterationをスキップ
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
                xmu = MP_mu / &
                    (u(i,j,k,MHN)+u(i,j,k,MH2)+u(i,j,k,MEL)+u(i,j,k,MHP)+yHe)  !Hm, H2p は量が少ないので無視
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
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) )*(Gamma-1.0d0)
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
  ! gammaが与えられていない場合（例えばrefinementの際に利用）  
  !-----------------------------------------------------------------------
  subroutine u2w_withgam(u,w,dv,gam)
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
                  +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0)

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
    
    !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
       i=globdbg_i+lbound(w,1)
       j=globdbg_j+lbound(w,2)
       k=globdbg_k+lbound(w,3)
       print '(A,5I8,/,1P6E15.7,/,1P6E15.7)', "(u2w, KS DEBUG)",&
            globdbg_myrank, globdbg_mygid, i-lbound(w,1), j-lbound(w,2), k-lbound(w,3),&
            u(i,j,k,MRHO),u(i,j,k,MP),u(i,j,k,MVX),u(i,j,k,MVY),u(i,j,k,MVZ),u(i,j,k,MH2), &
            w(i,j,k,MRHO),w(i,j,k,MP),w(i,j,k,MVX),w(i,j,k,MVY),w(i,j,k,MVZ),w(i,j,k,MH2)
    end if
#endif !KS DEBUG
    !------------------- KS DEBUG ---------------------!       
  end subroutine u2w_withgam

  !-----------------------------------------------------------------------
  ! convert conservative variables w to primitive variable u
  ! gammaが与えられている場合（例えば流体のアップデートの際に利用）
  !-----------------------------------------------------------------------
  subroutine w2u_withgam(w,u,dv,gam)
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
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) )*(gam(:,:,:)-1.0d0)

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

    !------------------- KS DEBUG ---------------------!
#ifdef KS_DEBUG
    if (globdbg_myrank==globdbg_rank .and. globdbg_mygid==globdbg_gid) then
       i=globdbg_i+lbound(w,1)
       j=globdbg_j+lbound(w,2)
       k=globdbg_k+lbound(w,3)
       print '(A,5I8,/,1P6E15.7,/,1P6E15.7)', "(w2u, KS DEBUG)" , &
            globdbg_myrank, globdbg_mygid, i-lbound(w,1), j-lbound(w,2), k-lbound(w,3),&
            w(i,j,k,MRHO),w(i,j,k,MP),w(i,j,k,MVX),w(i,j,k,MVY),w(i,j,k,MVZ),w(i,j,k,MH2), &
            u(i,j,k,MRHO),u(i,j,k,MP),u(i,j,k,MVX),u(i,j,k,MVY),u(i,j,k,MVZ),u(i,j,k,MH2)
    end if
#endif !KS_DEBUG
    !------------------- KS DEBUG ---------------------!       
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
  !get gamma consistent with chemical abundances
  !-----------------------------------------------------------------------
#if MODEL_ART > 0
  subroutine get_gamma(gid, gam)
    use grid, only : Imingh,Jmingh,Kmingh,Imaxgh,Jmaxgh,Kmaxgh, get_Ucomp
    use unit
    use primordial
    
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKGH),intent(OUT) :: gam

    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp
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
             xmu = MP_mu /(yhn(i,j,k)+yh2(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe)  !Hm, H2p は量が少ないので無視
             T_K = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb      !温度 [K]
             T_K = min(max(T_K,MP_Tmin),MP_Tmax) !上下の温度フロアを課す             
             gam(i,j,k) = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
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
