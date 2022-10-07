#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state. (Barotropic equation of state)
!-------------------------------------------------------------------------
module eos
  implicit none
  real(kind=DBL_KIND),parameter :: Cs = 1.0d0    ! isothermal sound speed
  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0 ! specific heat for polytrope
  real(kind=DBL_KIND),parameter :: Rhocr = 2.D-13/1.D-19  ! critical density betweeen isothermal and polytrope
  real(kind=DBL_KIND),save :: Kappa
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
  ! ------------
  ! macro for cs
  ! ------------
#define GETCS(CS_,RHO_) \
  flag = ( 1.d0 + sign(1.d0, (RHO_)-Rhocr) ) /2.d0 ;\
  CS_ = (1.d0-flag) * Cs + flag * sqrt( Gamma * Kappa * (RHO_) ** (Gamma -1.d0) )
  ! ------------------------------------------------------------
  ! macro for csbar
  ! flagi for isothermal (rhomin < rhomax < Rhocr)
  ! flagp for polytrope  (Rhocr < rhomin < rhomax)
  ! flagb for isothermal and polytrope (rhomin < Rhocr < rhomax)
  ! ------------------------------------------------------------
#define GETCSBAR(CSBAR_,RHOL_,RHOR_) \
  rhomax=max((RHOL_),(RHOR_)) ;\
  rhomin=min((RHOL_),(RHOR_)) ;\
  rho_bar=sqrt((RHOL_)*(RHOR_)) ;\
  delrho=max(5.0d-1*dlog(rhomax/rhomin),1.0d-5) ;\
  flagi = 0 ;\
  flagp = 0 ;\
  flagb = 0 ;\
  flagi = (1.d0+sign(1.d0, Rhocr-rhomax))/2 ;\
  flagp = (1.d0+sign(1.d0, rhomin-Rhocr))/2 ;\
  flagb = (1.d0-sign(1.d0, (rhomax-Rhocr)*(rhomin-Rhocr)))/2 ;\
  CSBAR_ = \
      (flagi * Cs \
      +flagp * sqrt(Kappa*(rho_bar**(Gamma-1.0d0))*sinh(Gamma*delrho)/sinh(delrho)) \
      +flagb * sqrt(abs((Kappa*(rhomax**Gamma)-(Cs**2)*rhomin))/((rhomax-rhomin)+1.0d-20*rho_bar)) \
      )/(flagi + flagp + flagb)
  ! --------------------
  ! macro for Pressure
  ! --------------------
#define GETP(P_,RHO_) \
  flag = ( 1.d0 + sign(1.d0, (RHO_)-Rhocr) ) /2 ;\
  P_ = (1-flag) * Cs**2 * (RHO_) + flag * Kappa*(RHO_)**Gamma

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
    real(kind=DBL_KIND) :: flag
    real(kind=DBL_KIND) :: rhomax, rhomin, rho_bar, flagi, flagp, flagb, delrho
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
    real(kind=DBL_KIND) :: flag, csp
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, p
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    level = get_level(id)
    h = CellWidth( :, level )
    rho => get_Ucomp(MRHO,id)
    vx => get_Ucomp(MVX,id)
    vy => get_Ucomp(MVY,id)
    vz => get_Ucomp(MVZ,id)
    x => get_Xp(id)
    y => get_Yp(id)
    z => get_Zp(id)
    dt = HUGE(dt)
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             GETCS(csp, rho(i,j,k))
             dt = (CFL)*min( &
                  h(MX)/(abs(vx(i,j,k))+csp), &
                  h(MY)/(abs(vy(i,j,k))+csp), &
                  h(MZ)/(abs(vz(i,j,k))+csp) )
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
    w(:,:,:,MVX)=u(:,:,:,MVX)
    w(:,:,:,MVY)=u(:,:,:,MVY)
    w(:,:,:,MVZ)=u(:,:,:,MVZ)
!!$    w(:,:,:,MVX)=u(:,:,:,MVX)*w(:,:,:,MRHO)
!!$    w(:,:,:,MVY)=u(:,:,:,MVY)*w(:,:,:,MRHO)
!!$    w(:,:,:,MVZ)=u(:,:,:,MVZ)*w(:,:,:,MRHO)
#ifdef MPSI
    w(:,:,:,MPSI) = u(:,:,:,MPSI) * dv
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
    u(:,:,:,MVX) =w(:,:,:,MVX)
    u(:,:,:,MVY) =w(:,:,:,MVY)
    u(:,:,:,MVZ) =w(:,:,:,MVZ)
!!$    u(:,:,:,MVX) =w(:,:,:,MVX)/w(:,:,:,MRHO)
!!$    u(:,:,:,MVY) =w(:,:,:,MVY)/w(:,:,:,MRHO)
!!$    u(:,:,:,MVZ) =w(:,:,:,MVZ)/w(:,:,:,MRHO)
#ifdef MPSI
    u(:,:,:,MPSI) = w(:,:,:,MPSI) / dv
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
  !
  !-----------------------------------------------------------------------
end module eos
subroutine eos_init
  use eos
  Kappa = Cs**2 * Rhocr ** (1.d0-Gamma) ! kappa of polytrope
end subroutine eos_init
