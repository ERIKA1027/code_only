
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

  subroutine flux( ql, qr, f)
    use grid, only : Mmin, CellWidth, LevelMax, Lmin, Dtime
    use parameter, only : Pi4i, Pi4
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: f   ! (OUT)
    real(kind=DBL_KIND) :: sw_deg_l, sw_deg_r, swl, swr, swal, swar, swml, swmr, swt
    real(kind=DBL_KIND),parameter :: eps = 1.D-6
    integer :: i,j,k, is,js,ks,ie,je,ke
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
         rho, u, v, w, bx, by, bz, p, e, pt

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

    gm1i = 1.d0/(Gamma - 1)
    sqrtpi4i = sqrt(Pi4i)
    sqrtpi4 = sqrt(Pi4)
    do k = ks, ke
       do j = js, je
          do i = is, ie
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
             rhol = abs(ql(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
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
             rhor = abs(qr(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
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
             csl2 = Gamma * pl * rhoil ! sound speed
             csr2 = Gamma * pr * rhoir ! sound speed
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
             ptml = ptal
             uml = sm
             vml = vm
             wml = wm
             byml = bym
             bzml = bzm
             eml = eal - sql*(ual*bxm+val*byal+wal*bzal - uml*bxm-vml*byml-wml*bzml)*signbxm

             rhomr = rhoar
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
             f(i,j,k, MVX) = rho*u**2 + pt - bx**2
             f(i,j,k, MVY) = rho*v*u - bx*by
             f(i,j,k, MVZ) = rho*w*u - bx*bz
             f(i,j,k, MBY) = (by*u - bx*v)*sqrtpi4
             f(i,j,k, MBZ) = (bz*u - bx*w)*sqrtpi4
             f(i,j,k,  MP) = (e + pt)*u - bx*(u*bx+v*by+w*bz)
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
    ca= sqrt( abs( ( Gamma*p +(bx**2+by**2+bz**2)*pi4i )/rho ) ) ! fast wave
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
  end function get_dt_by_cflcond
  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable w
  !-----------------------------------------------------------------------
  subroutine u2w(u,w,dv)
    use parameter
    use grid, only : Mmin
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

    w(:,:,:,MRHO)=u(:,:,:,MRHO)*dv
    w(:,:,:,MVX)=u(:,:,:,MVX)*w(:,:,:,MRHO)
    w(:,:,:,MVY)=u(:,:,:,MVY)*w(:,:,:,MRHO)
    w(:,:,:,MVZ)=u(:,:,:,MVZ)*w(:,:,:,MRHO)
    w(:,:,:,MBX)= u(:,:,:,MBX)*dv
    w(:,:,:,MBY)= u(:,:,:,MBY)*dv
    w(:,:,:,MBZ)= u(:,:,:,MBZ)*dv
    w(:,:,:,MP) =dv*(u(:,:,:,MP)/(Gamma-1.d0) &
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2.d0 &
         +(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i)
    w(:,:,:,MDB) = u(:,:,:,MDB)*dv
#ifdef MPSI
    w(:,:,:,MPSI) = u(:,:,:,MPSI) * dv
    w(:,:,:,MGX) = u(:,:,:,MGX) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGY) = u(:,:,:,MGY) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGZ) = u(:,:,:,MGZ) * u(:,:,:,MRHO) * dv
#endif !MPSI
  end subroutine u2w
  !-----------------------------------------------------------------------
  ! convert conservative variables w to primitive variable u
  !-----------------------------------------------------------------------
  subroutine w2u(w,u,dv)
    use parameter
    use grid, only : Mmin
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

    u(:,:,:,MRHO) =w(:,:,:,MRHO)/dv
    u(:,:,:,MVX) =w(:,:,:,MVX)/w(:,:,:,MRHO)
    u(:,:,:,MVY) =w(:,:,:,MVY)/w(:,:,:,MRHO)
    u(:,:,:,MVZ) =w(:,:,:,MVZ)/w(:,:,:,MRHO)
    u(:,:,:,MBX) =w(:,:,:,MBX)/dv
    u(:,:,:,MBY) =w(:,:,:,MBY)/dv
    u(:,:,:,MBZ) =w(:,:,:,MBZ)/dv
    u(:,:,:,MP) = ( w(:,:,:,MP)/dv -0.5d0*u(:,:,:,MRHO) &
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) &
         -(u(:,:,:,MBX)**2+u(:,:,:,MBY)**2+u(:,:,:,MBZ)**2)*pi8i &
         )*(Gamma-1.0d0)
    u(:,:,:,MDB) = w(:,:,:,MDB)/dv
#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
    u(:,:,:,MGX) = w(:,:,:,MGX) / (u(:,:,:,MRHO) * dv)
    u(:,:,:,MGY) = w(:,:,:,MGY) / (u(:,:,:,MRHO) * dv)
    u(:,:,:,MGZ) = w(:,:,:,MGZ) / (u(:,:,:,MRHO) * dv)
#endif !MPSI
  end subroutine w2u
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
end module eos
subroutine eos_init
  use eos
  real(kind=8) :: csd

end subroutine eos_init
