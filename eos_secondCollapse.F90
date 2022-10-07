#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state. (Barotropic equation of state)
!-------------------------------------------------------------------------
module eos
  implicit none
  real(kind=DBL_KIND),parameter :: Cs = 1.0d0    ! isothermal sound speed
  real(kind=DBL_KIND),parameter :: Gamma1 = 1.4d0 ! specific heat for polytrope
  real(kind=DBL_KIND),parameter :: Gamma2 = 1.1d0 ! specific heat for polytrope
  real(kind=DBL_KIND),parameter :: Rhocr1 = 1.D-13/1.D-16  ! critical density betweeen isothermal and polytrope of Gamma1
  real(kind=DBL_KIND),parameter :: Rhocr2 = 1.D-9/1.D-16  ! critical density betweeen polytrope gas of Gamma1 and Gamma2
!!$  real(kind=DBL_KIND),parameter :: Rhocr1 = 2.D-13/1.D-19  ! critical density betweeen isothermal and polytrope of Gamma1
!!$  real(kind=DBL_KIND),parameter :: Rhocr2 = 1.D-9/1.D-19  ! critical density betweeen polytrope gas of Gamma1 and Gamma2
  real(kind=DBL_KIND),save :: Kappa1, Kappa2
  real(kind=DBL_KIND),parameter :: Rhocr = Rhocr1 ! backward comatibility
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
  ! get numerical flux in one dimension
  !-----------------------------------------------------------------------
  ! macro for entropy condition
  ! ----------------------------
#define ELMOD(EL,ELBAR,EL1, EL2) \
  eps=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-eps) ;\
  x2=1.d0-x1 ;\
  EL= x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(eps+x1)+eps)
  subroutine flux( ql, qr, f)
    use grid, only : Mmin
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,Mmin:) :: f   ! (OUT)
    integer :: i,j,k, is,js,ks,ie,je,ke
    real(kind=DBL_KIND) :: &
         rhobar, csbar, ubar, vbar, wbar, sr0, sr1, sri, csi, &
         rhol,vxl,vyl,vzl,pl,csl, &
         rhor,vxr,vyr,vzr,pr,csr, &
         el1,el2,el3,el4, &
         w1,w2,w3,w4, &
         fl1,fl2,fl3,fl4, &
         fr1,fr2,fr3,fr4
    real(kind=DBL_KIND) :: eps, x1, x2
    real(kind=DBL_KIND) :: flag0, flag1, flag2
    real(kind=DBL_KIND) :: rhomax, rhomin, rho_bar, delrho
    is = max( lbound(ql,1), lbound(qr,1), lbound(f,1) )
    js = max( lbound(ql,2), lbound(qr,2), lbound(f,2) )
    ks = max( lbound(ql,3), lbound(qr,3), lbound(f,3) )
    ie = min( ubound(ql,1), ubound(qr,1), ubound(f,1) )
    je = min( ubound(ql,2), ubound(qr,2), ubound(f,2) )
    ke = min( ubound(ql,3), ubound(qr,3), ubound(f,3) )

    do k = ks, ke
       do j = js, je
          do i = is, ie
             ! left variable
             rhol = abs(ql(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
             vxl  = ql(i,j,k,MVX)
             vyl  = ql(i,j,k,MVY)
             vzl  = ql(i,j,k,MVZ)
             ! right variable
             rhor = abs(qr(i,j,k,MRHO))
             vxr  = qr(i,j,k,MVX)
             vyr  = qr(i,j,k,MVY)
             vzr  = qr(i,j,k,MVZ)
             ! get cs and P
             GETCSBAR(csbar, rhol, rhor)
             GETCS(csl, rhol)
             GETCS(csr, rhor)
             GETP(pl, rhol)
             GETP(pr, rhor)
             ! this routine returns the numerical fluxes (x-direction)
             ! of 1st order acuracy, given left ans right physical value.
             ! computaiotn of fluxes (left)
             fl1=rhol*vxl
             fl2=fl1*vxl+pl
             fl3=fl1*vyl
             fl4=fl1*vzl
             ! computaiotn of fluxes (right)
             fr1=rhor*vxr
             fr2=fr1*vxr+pr
             fr3=fr1*vyr
             fr4=fr1*vzr
             ! bared variables
             csi = 1.d0/csbar
             sr0=sqrt(rhol)
             sr1=sqrt(rhor)
             sri=1.0d0/(sr0+sr1)
             rhobar=sr0*sr1
             ubar=(sr0*vxl+sr1*vxr)*sri
             vbar=(sr0*vyl+sr1*vyr)*sri
             wbar=(sr0*vzl+sr1*vzr)*sri
             ! entropy condition & cal. eigen value of marix a ; eli
             ELMOD(el1,ubar+csbar,vxl+csl,vxr+csr)
             ELMOD(el2,ubar-csbar,vxl-csl,vxr-csr)
             ELMOD(el3,ubar,      vxl,    vxr    )
             ELMOD(el4,ubar,      vxl,    vxr    )
             ! amplitude; w = sigma wi*ri; ri are eigenvector
             w1=( vxr-vxl+csbar*(rhor-rhol)/rhobar)*rhobar*5.0d-1*csi
             w2=(-vxr+vxl+csbar*(rhor-rhol)/rhobar)*rhobar*5.0d-1*csi
             w3=vyr-vyl
             w4=vzr-vzl
             ! fulx spliting; computation of f(i+1/2,j)
             f(i,j,k,MRHO) = (fl1+fr1-el1*w1-el2*w2)*0.5d0
             f(i,j,k,MVX)  = (fl2+fr2-el1*w1*(ubar+csbar)-el2*w2*(ubar-csbar))*0.5d0
             f(i,j,k,MVY)  = (fl3+fr3-el1*w1*vbar-el2*w2*vbar-el3*w3*rhobar)*0.5d0
             f(i,j,k,MVZ)  = (fl4+fr4-el1*w1*wbar-el2*w2*wbar-el4*w4*rhobar)*0.5d0
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
    real(kind=DBL_KIND) :: flag0, flag1, flag2, csp
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, p, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    level = get_level(id)
    h = CellWidth( :, level )
    rho => get_Ucomp(MRHO,id)
    vx => get_Ucomp(MVX,id)
    vy => get_Ucomp(MVY,id)
    vz => get_Ucomp(MVZ,id)
    gx => get_Ucomp(MGX,id)
    gy => get_Ucomp(MGY,id)
    gz => get_Ucomp(MGZ,id)
    x => get_Xp(id)
    y => get_Yp(id)
    z => get_Zp(id)
    dt = HUGE(dt)
#define DTI_(V_,G_,H_) (H_)/sqrt((V_)**2+2.d0*(G_)*(H_))
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             GETCS(csp, rho(i,j,k))
             dt = min(dt, &
                  (CFL) /( &
                  DTI_( abs(vx(i,j,k))+csp, abs(gx(i,j,k)), h(MX) ) + &
                  DTI_( abs(vy(i,j,k))+csp, abs(gy(i,j,k)), h(MY) ) + &
                  DTI_( abs(vz(i,j,k))+csp, abs(gz(i,j,k)), h(MZ) ) ))
          enddo
       enddo
    enddo
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
#ifdef MPSI
    w(:,:,:,MPSI) = u(:,:,:,MPSI) * dv
    w(:,:,:,MGX) = u(:,:,:,MGX) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGY) = u(:,:,:,MGY) * u(:,:,:,MRHO) * dv
    w(:,:,:,MGZ) = u(:,:,:,MGZ) * u(:,:,:,MRHO) * dv
#endif
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
#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
    u(:,:,:,MGX) = w(:,:,:,MGX) / w(:,:,:,MRHO)
    u(:,:,:,MGY) = w(:,:,:,MGY) / w(:,:,:,MRHO)
    u(:,:,:,MGZ) = w(:,:,:,MGZ) / w(:,:,:,MRHO)
#endif
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
  subroutine get_cs(c, rho)
    real(kind=DBL_KIND),intent(OUT) :: c
    real(kind=DBL_KIND),intent(IN) :: rho
    real(kind=DBL_KIND) :: flag0, flag1, flag2
    GETCS(c, rho)
  end subroutine get_cs
  !-----------------------------------------------------------------------
  ! return pressure (p) when density (rho) is given
  !-----------------------------------------------------------------------
  subroutine get_p(p, rho)
    real(kind=DBL_KIND),intent(OUT) :: p
    real(kind=DBL_KIND),intent(IN) :: rho
    real(kind=DBL_KIND) :: flag0, flag1, flag2
    GETP(p, rho)
  end subroutine get_p
end module eos
subroutine eos_init
  use eos
  Kappa1 = Cs**2 * Rhocr1 ** (1.d0-Gamma1) ! kappa of polytrope
  Kappa2 = Kappa1 * Rhocr2 ** (Gamma1 - Gamma2)
end subroutine eos_init
