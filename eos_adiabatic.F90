
#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state.
!-------------------------------------------------------------------------
module eos
  implicit none
!!$  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0
  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0
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
  ! get numerical flux in one dimension (Roe's original method) Roe 1981, JCP, 43, 357
  !-----------------------------------------------------------------------
  ! macro for entropy condition
#define ELMOD(EL,ELBAR,EL1, EL2) \
  eps=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-eps) ;\
  x2=1.d0-x1 ;\
  EL= x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(eps+x1)+eps)

  subroutine flux( ql, qr, f)
    use grid, only : Mmin, Imingh,Jmingh,Kmingh
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: f   ! (OUT)
    integer :: i,j,k, is,js,ks,ie,je,ke
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,hl,el,cl, &
         rhor,ur,vr,wr,pr,hr,er,cr, &
         dq1,dq2,dq3,dq4,dq5, &
         sql,sqr,sqa,rhob,ub,vb,wb,hb,qb2,cb2,cb,cb24i,ub2,vb2,wb2, &
         gm1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10, &
         w1,w2,w3,w4,w5, &
         el1,el2,el3,el4,el5, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5 
    real(kind=DBL_KIND) :: eps, x1, x2
    is = max( lbound(ql,1), lbound(qr,1), lbound(f,1) )
    js = max( lbound(ql,2), lbound(qr,2), lbound(f,2) )
    ks = max( lbound(ql,3), lbound(qr,3), lbound(f,3) )
    ie = min( ubound(ql,1), ubound(qr,1), ubound(f,1) )
    je = min( ubound(ql,2), ubound(qr,2), ubound(f,2) )
    ke = min( ubound(ql,3), ubound(qr,3), ubound(f,3) )

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
             dq1 = rhor - rhol
             dq2 = rhor * ur - rhol * ul
             dq3 = rhor * vr - rhol * vl
             dq4 = rhor * wr - rhol * wl
             dq5 = er - el
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
             cb24i = 1/(cb2 * 4)
             ! ----------------------
             ! w = R |Lambda| R^(-1)
             ! ----------------------
             ELMOD(el1, ub-cb, ul-cl, ur-cr)
             ELMOD(el2, ub,    ul,    ur)
             ELMOD(el3, ub,    ul,    ur)
             ELMOD(el4, ub,    ul,    ur)
             ELMOD(el5, ub+cb, ul+cl, ur+cr)

             a1 = vb*dq3
             a2 = wb*dq4
             a3 = dq1*qb2
             a4 = (2*dq5 + a3 - 2*(dq2*ub + a1 + a2))
             a5 = cb*(el1 - el5)
             a6 = (el1 + el5)
             a7 = (a6 - 2*el4)
             a8 = gm1*a7*a4
             a9 = a8 - a5*(dq2 - dq1*ub)*2
             a10 = el4*qb2

             w1 = (4*cb2*dq1*el4 + a9)*cb24i

             w2 = (2*cb2*(dq2*a6 - dq1*a7*ub) & 
                  - a5*( gm1*(a3 + 2*(dq5 - a1 - a2)) &
                  + 2*( (2 - Gamma)*dq2*ub - dq1*ub2)) &
                  + a8*ub)*cb24i

             w3 = (a9*vb + 4*cb2*(dq3*el3 + dq1*(-el3 + el4)*vb))*cb24i

             w4 = (a9*wb + 4*cb2*(dq4*el2 + dq1*(-el2 + el4)*wb))*cb24i

             w5 = (gm1*(hb*a6 - a10)*a4 &
                  + 2*cb2*(dq2*a6*ub + 2*(el3*a1 + el2*a2) &
                  + dq1*(a10 - a6*ub2 - 2*(el3*vb2 + el2*wb2))) &
                  - a5*(2*dq2*(hb - gm1*ub2) &
                  + ub*(2*dq5*gm1 - dq1*(2*hb - qb2*gm1) &
                  - 2*gm1*(a1 + a2))))*cb24i

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
             f(i,j,k,MRHO) = (fl1 + fr1 - w1)/2
             f(i,j,k,MVX)  = (fl2 + fr2 - w2)/2
             f(i,j,k,MVY)  = (fl3 + fr3 - w3)/2
             f(i,j,k,MVZ)  = (fl4 + fr4 - w4)/2
             f(i,j,k,MP)   = (fl5 + fr5 - w5)/2

          enddo
       enddo
    enddo
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
  ! return sound speed (c) when density (rho) is given
  !-----------------------------------------------------------------------
  subroutine get_cs(c, rho, p)
    real(kind=DBL_KIND),intent(OUT) :: c
    real(kind=DBL_KIND),intent(IN) :: rho, p
    real(kind=DBL_KIND) :: flag
    c = sqrt(abs(Gamma*p/rho))
  end subroutine get_cs
!!$  !-----------------------------------------------------------------------
!!$  ! return pressure (p) when density (rho) is given
!!$  !-----------------------------------------------------------------------
!!$  subroutine get_p(p, rho)
!!$    real(kind=DBL_KIND),intent(OUT) :: p
!!$    real(kind=DBL_KIND),intent(IN) :: rho
!!$    real(kind=DBL_KIND) :: flag
!!$    GETP(p, rho)
!!$  end subroutine get_p
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
end module eos
subroutine eos_init
  use eos
  real(kind=8) :: csd

end subroutine eos_init
