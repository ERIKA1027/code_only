#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state.
!-------------------------------------------------------------------------
module eos
  implicit none
  ! Gamma is given by eos_init()
  real(kind=DBL_KIND),save :: Gamma

!!$  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0
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

  subroutine flux( ql, qr, pratio, f, ncrd )
    use grid, only : Imin,Jmin,Kmin,Mmin, Imax,Jmax,Kmax,Mmax, Imingh,Jmingh,Kmingh
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
             gm1 = Gamma - 1
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
             cl = sqrt( Gamma * pl / rhol )
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
             cr = sqrt( Gamma * pr / rhor )
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

             
             f(i,j,k,MHN) = 0.5d0 * ub * (ql(i,j,k,MHN)*rhol + qr(i,j,k,MHN)*rhor) - 0.5d0 * abs(ub) * (qr(i,j,k,MHN)*rhor - ql(i,j,k,MHN)*rhol) ! Up-wind
             f(i,j,k,MHP) = 0.5d0 * ub * (ql(i,j,k,MHP)*rhol + qr(i,j,k,MHP)*rhor) - 0.5d0 * abs(ub) * (qr(i,j,k,MHP)*rhor - ql(i,j,k,MHP)*rhol) ! Up-wind             
          enddo
       enddo
    enddo
#ifdef EMULATE_2DIM
    do k = Kmin-1, Kmax
       if (k /= Kmin) f(:,:,k,:) = f(:,:,Kmin,:)
    enddo
#endif !EMULATE_2DIM
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
    cs = sqrt(abs(Gamma*p/rho))
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
  end function get_dt_by_cflcond
  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable w
  !-----------------------------------------------------------------------
  subroutine u2w(u,w,dv)
    use grid, only : Mmin
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

    w(:,:,:,MRHO)=u(:,:,:,MRHO)*dv
    w(:,:,:,MVX)=u(:,:,:,MVX)*w(:,:,:,MRHO)
    w(:,:,:,MVY)=u(:,:,:,MVY)*w(:,:,:,MRHO)
    w(:,:,:,MVZ)=u(:,:,:,MVZ)*w(:,:,:,MRHO)
    w(:,:,:,MP) =dv*(u(:,:,:,MP)/(Gamma-1.d0) &
         +u(:,:,:,MRHO)*(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2)/2)
    w(:,:,:,MHN)=u(:,:,:,MHN)*w(:,:,:,MRHO) !
    w(:,:,:,MHP)=u(:,:,:,MHP)*w(:,:,:,MRHO) !
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
    use grid, only : Mmin
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: w !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: u !(OUT)
    real(kind=DBL_KIND),intent(IN) :: dv

    u(:,:,:,MRHO) =w(:,:,:,MRHO)/dv
    u(:,:,:,MVX) =w(:,:,:,MVX)/w(:,:,:,MRHO)
    u(:,:,:,MVY) =w(:,:,:,MVY)/w(:,:,:,MRHO)
    u(:,:,:,MVZ) =w(:,:,:,MVZ)/w(:,:,:,MRHO)
    u(:,:,:,MP) = ( w(:,:,:,MP)/dv -0.5d0*u(:,:,:,MRHO) &
         *(u(:,:,:,MVX)**2+u(:,:,:,MVY)**2+u(:,:,:,MVZ)**2) )*(Gamma-1.0d0)
    u(:,:,:,MHN)=w(:,:,:,MHN)/w(:,:,:,MRHO) !
    u(:,:,:,MHP)=w(:,:,:,MHP)/w(:,:,:,MRHO) !
#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
    u(:,:,:,MGX) = w(:,:,:,MGX) / w(:,:,:,MRHO)
    u(:,:,:,MGY) = w(:,:,:,MGY) / w(:,:,:,MRHO)
    u(:,:,:,MGZ) = w(:,:,:,MGZ) / w(:,:,:,MRHO)
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
  use eos, only : Gamma
  use unit, only: ModelParam_gamma
  real(kind=8) :: csd
  Gamma = ModelParam_gamma
end subroutine eos_init
