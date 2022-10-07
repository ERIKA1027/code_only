#include "config.h"
#include "barotropic.h"
!-------------------------------------------------------------------------
! module depending on equation of state (and system equations).
! Adiabatic MHD
!  Numerical flux is based on Hanawa & Fukuda (19...)
!-------------------------------------------------------------------------
module eos
  implicit none
  real(kind=DBL_KIND),parameter :: Cs = 1.0d0    ! isothermal sound speed
  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0 ! specific heat for polytrope
!  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0 ! specific heat for polytrope
!!$  real(kind=DBL_KIND),parameter :: Rhocr = 2.D-13/1.D-19  ! critical density betweeen isothermal and polytrope
!!$  real(kind=DBL_KIND),parameter :: Rhocr = 1.D-13/1.D-16  ! critical density betweeen isothermal and polytrope
  real(kind=DBL_KIND),parameter :: Rhocr = 1.D10
  real(kind=DBL_KIND),save :: Kappa
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
  ! get numerical flux in one dimension
  !-----------------------------------------------------------------------
  ! macro for entropy condition
  ! ----------------------------
#define ELMOD(EL,ELBAR,EL1, EL2) \
  epsel=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-epsel) ;\
  x2=1.d0-x1 ;\
  EL= -(x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(epsel+x1)+epsel))
#define FLMT(x,y) max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x))
!!$#define FLMT(x,y) (x)*0.d0+(y)*0.d0
  subroutine flux(u, f1d, ncrd)
    use parameter
    use grid, only : Imingh, Jmingh, Kmingh, Imaxgh, Jmaxgh, Kmaxgh, &
         Imin, Jmin, Kmin, Mmin, Imax, Jmax, Kmax, Mmax, &
         Lmin, CellWidth, Dtime, BlockMem, LevelMax
    use util, only : util_arroffset
    integer,intent(IN) :: ncrd
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: u    ! (IN)
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: f1d  ! (OUT)
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: q, f
    integer :: i,j,k, is,js,ks,ie,je,ke,io,jo,ko,io2,jo2,ko2
    real(kind=DBL_KIND),parameter :: eps = 1.0d-5
    real(kind=DBL_KIND),parameter :: eps2=2.0d-5
    real(kind=DBL_KIND) :: &
         alphf, alphs, astar2, betay, betaz, bxave, bxbar, byave,  &
         bybar, bybz, bzave, bzbar, ca, ca2, cfa, cfast, cfast2, &
         cfca, cfcs, cs2x, csca, cslow, cslow2, csx, dby, dbz, detp,&
         detq, drho, dvx, dvy, dvz, ea1, ea2, ea3, ea4, ea5, ea6, &
         ea7, eb1, eb2, eb3, eb4, eb5, eb6, eb7, eeps, el1, el2, el3, &
         el4, el5, el6, el7, fl1, fl2, fl3, fl4, fl5, fl6, &
         fr1, fr2, fr3, fr4, fr5, fr6, p11, p12, p21, p22, pl, pr, &
         q11, q12, q21, q22, re1, re2, re3, re4, re5, re6, rhoave, &
         rhoavi, rhobar, rma1, rma2, rma3, rma4, rma5, rma6, rmf1, &
         rmf2, rmf3, rmf4, rmf5, rmf6, rms1, rms2, rms3, rms4, rms5, &
         rms6, rpa1, rpa2, rpa3, rpa4, rpa5, rpa6, rpf1, rpf2, rpf3, &
         rpf4, rpf5, rpf6, rps1, rps2, rps3, rps4, rps5, rps6, s1, &
         s2, s3, sgnbx, sgr, sgr2, sp, sp2, sr0, sr1, sri, t1, t3, &
         vxave, vxbar, vyave, vybar, vzave, vzbar, w1, w2, w3, w4, &
         w5, w6, w7, wl1, wl2, wl3, wl4, wl5, wl6, wl7, wm1, wm2, &
         wm3, wm4, wm5, wm6, wm7, wr1, wr2, wr3, wr4, wr5, wr6, wr7, &
         rhol,vxl,vyl,vzl,bxl,byl,bzl,psil,cl, &
         rhor,vxr,vyr,vzr,bxr,byr,bzr,psir,cr, &
         rholl,vxll,vyll,vzll,bxll,byll,bzll,psill,cll, &
         rhorr,vxrr,vyrr,vzrr,bxrr,byrr,bzrr,psirr,crr, &
         csbar, cs2, &
         bxm, psim, bxlm, bxrm, psilm, psirm
    real(kind=DBL_KIND) :: flag
    real(kind=DBL_KIND) :: rhomax, rhomin, rho_bar, flagi, flagp, flagb, delrho
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND) :: epsel, x1, x2
    integer,dimension(Mmin:Mmax) :: mcycle
    integer :: m

    mcycle = cyclecomp( ncrd )
    q(:,:,:,:) = u(:,:,:,mcycle)
    f(:,:,:,:) = 0.d0           ! initialize

#ifdef SINGLE_STEP
    h = CellWidth( :, LevelMax )
#else  !SINGLE_STEP
    h = CellWidth( :, Lmin )
#endif !SINGLE_STEP
    Ch = (CFL) * minval(h) / Dtime(Lmin) / 3.d0

    call util_arroffset(ncrd,io,jo,ko)
    io2 = io*2
    jo2 = jo*2
    ko2 = ko*2
    is = Imin-io
    js = Jmin-jo
    ks = Kmin-ko
    ie = Imax
    je = Jmax
    ke = Kmax

    do k = ks, ke
       do j = js, je
          do i = is, ie
             ! -------------------------
             ! left and right variables
             ! -------------------------
             rholl = q(i-io,j-jo,k-ko,MRHO)
             vxll =  q(i-io,j-jo,k-ko,MVX)
             vyll =  q(i-io,j-jo,k-ko,MVY)
             vzll =  q(i-io,j-jo,k-ko,MVZ)
             bxll =  q(i-io,j-jo,k-ko,MBX)
             byll =  q(i-io,j-jo,k-ko,MBY)
             bzll =  q(i-io,j-jo,k-ko,MBZ)
             psill = q(i-io,j-jo,k-ko,MDB)
             rhol =  q(i,j,k,MRHO)
             vxl =   q(i,j,k,MVX)
             vyl =   q(i,j,k,MVY)
             vzl =   q(i,j,k,MVZ)
             bxl =   q(i,j,k,MBX)
             byl =   q(i,j,k,MBY)
             bzl =   q(i,j,k,MBZ)
             psil =  q(i,j,k,MDB)
             rhor =  q(i+io,j+jo,k+ko,MRHO)
             vxr =   q(i+io,j+jo,k+ko,MVX)
             vyr =   q(i+io,j+jo,k+ko,MVY)
             vzr =   q(i+io,j+jo,k+ko,MVZ)
             bxr =   q(i+io,j+jo,k+ko,MBX)
             byr =   q(i+io,j+jo,k+ko,MBY)
             bzr =   q(i+io,j+jo,k+ko,MBZ)
             psir =  q(i+io,j+jo,k+ko,MDB)
             rhorr = q(i+io2,j+jo2,k+ko2,MRHO)
             vxrr =  q(i+io2,j+jo2,k+ko2,MVX)
             vyrr =  q(i+io2,j+jo2,k+ko2,MVY)
             vzrr =  q(i+io2,j+jo2,k+ko2,MVZ)
             bxrr =  q(i+io2,j+jo2,k+ko2,MBX)
             byrr =  q(i+io2,j+jo2,k+ko2,MBY)
             bzrr =  q(i+io2,j+jo2,k+ko2,MBZ)
             psirr = q(i+io2,j+jo2,k+ko2,MDB)

             GETCSBAR(csbar,rhol,rhor)
             cs2 = csbar**2
             GETP(pl, rhol)
             GETP(pr, rhor)
             !-----barred variables
             sr0=sqrt(rhol)
             sr1=sqrt(rhor)
             sri=1.0d0/(sr0+sr1)
             rhobar=sr0*sr1
             vxbar=(sr0*vxl+sr1*vxr)*sri
             vybar=(sr0*vyl+sr1*vyr)*sri
             vzbar=(sr0*vzl+sr1*vzr)*sri
             bxbar=(sr0*bxr+sr1*bxl)*sri
             bybar=(sr0*byr+sr1*byl)*sri
             bzbar=(sr0*bzr+sr1*bzl)*sri
             !-----arithmetric average value
             rhoave=0.5*(rhol+rhor)
             rhoavi=1.0d0/rhoave
             vxave=0.5*(vxl+vxr)
             vyave=0.5*(vyl+vyr)
             vzave=0.5*(vzl+vzr)
             bxave=0.5*(bxl+bxr)
             byave=(byl+byr)/2.0d0
             bzave=(bzl+bzr)/2.0d0

             ! ---------------------
             ! add-on for div B free
             ! ---------------------
             bxlm =bxl +0.5*FLMT(bxr -bxl, bxl-bxll)
             bxrm =bxr -0.5*FLMT(bxr -bxl, bxrr-bxr)
             psilm=psil+0.5*FLMT(psir-psil,psil-psill)
             psirm=psir-0.5*FLMT(psir-psil,psirr-psir)
             bxm  = 0.5d0*(bxlm+bxrm - (psirm-psilm)/Ch)
             psim = 0.5d0*(psilm+psirm - (bxrm-bxlm)*Ch)
             f(i,j,k,MBX) = psim
             f(i,j,k,MDB) = bxm*Ch**2
             bxbar=bxm
             bxave=bxbar
             !-----characteristic speed
             cs2x=cs2+((byl-byr)**2+(bzl-bzr)**2)*sri**2*pi8i
             astar2=cs2x+(bxbar**2+bybar**2+bzbar**2)*pi4i/rhobar
             ca2=bxbar**2/(pi4*rhobar)
             bybz=(bybar**2+bzbar**2)*pi4i/rhobar
             csca=cs2x-ca2
             cfcs=sqrt(csca**2+bybz*(2.0d0*(cs2x+ca2)+bybz))
             cfast2=0.5d0*(astar2+cfcs)
             cslow2=cs2x*ca2/cfast2
             cfast=sqrt(cfast2)
             cslow=sqrt(cslow2)
             ca=sqrt(ca2)
             csx=sqrt(cs2x)
             !-----flux_l
             fl1=rhol*vxl
             fl2=fl1*vxl+pl+(-bxbar**2+byl**2+bzl**2)*pi8i
             fl3=fl1*vyl-bxbar*byl*pi4i
             fl4=fl1*vzl-bxbar*bzl*pi4i
             fl5=vxl*byl-vyl*bxbar
             fl6=vxl*bzl-vzl*bxbar
             !-----flux_r
             fr1=rhor*vxr
             fr2=fr1*vxr+pr+(-bxbar**2+byr**2+bzr**2)*pi8i
             fr3=fr1*vyr-bxbar*byr*pi4i
             fr4=fr1*vzr-bxbar*bzr*pi4i
             fr5=vxr*byr-vyr*bxbar
             fr6=vxr*bzr-vzr*bxbar
             !-----for singular points
             sgr=bybar**2+bzbar**2-eps
             sp=0.5d0+sign(0.5,sgr)
             betay=sp*bybar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp)) &
                  +sqrt(0.5d0)*(1.0d0-sp)
             betaz=sp*bzbar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp)) &
                  +sqrt(0.5d0)*(1.0d0-sp)
             sgr2=cfcs-eps2
             sp2=0.5d0+sign(0.5,sgr2)
             cfca=(bybz*(2.0*(ca2+cs2x)+bybz))/(abs(csca)+cfcs) +max(csca,0.0)
             cfa=(bybz*(2.0*(ca2+cs2x)+bybz))/(abs(csca)+cfcs) -min(csca,0.0)
             alphf=sp2*sqrt(cfca/(cfcs+1.0d0-sp2)) +1.0d0-sp2
             alphs=sp2*sqrt((cfa)/(cfcs+1.0d0-sp2))
             sgnbx=sign(1.0,bxbar)
             !-----eigen values
             ELMOD(el1,vxbar+cfast, vxl+cfast, vxr+cfast)
             ELMOD(el7,vxbar-cfast, vxl-cfast, vxr-cfast)
             ELMOD(el3,vxbar+cslow, vxl+cslow, vxr+cslow)
             ELMOD(el5,vxbar-cslow, vxl-cslow, vxr-cslow)
             ELMOD(el2,vxbar+ca,    vxl+ca,    vxr+ca)
             ELMOD(el6,vxbar-ca,    vxl-ca,    vxr-ca)
             ELMOD(el4,vxbar,       vxl,       vxr)
             !-----calc amplitude
             p11=alphs*cfast*sqrt(pi4/rhobar)
             p12=-alphf*cs2x/cfast*sqrt(pi4/rhobar)
             p21=alphf
             p22=alphs
             !
             q11=alphf*cfast
             q12=alphs*cslow
             q21=-alphs*ca*sgnbx
             q22=alphf*csx*sgnbx
             !
             detp=p11*p22-p12*p21
             detq=q11*q22-q12*q21
             !-----midium amplitude
             drho=rhor-rhol
             dvx=rhobar*rhoavi*((rhor*vxr-rhol*vxl)-vxave*drho)
             dvy=rhobar*rhoavi*((rhor*vyr-rhol*vyl)-vyave*drho)
             dvz=rhobar*rhoavi*((rhor*vzr-rhol*vzl)-vzave*drho)
             dby=byr-byl
             dbz=bzr-bzl
             !
             t1=betay*dby+betaz*dbz
             t3=betaz*dby-betay*dbz
             !
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             !
             wm1=0.5d0*((p22*t1-p12*drho)/detp+(q22*s1-q12*s2)/detq)
             wm7=0.5d0*((p22*t1-p12*drho)/detp-(q22*s1-q12*s2)/detq)
             wm3=0.5d0*((-p21*t1+p11*drho)/detp+(-q21*s1+q11*s2)/detq)
             wm5=0.5d0*((-p21*t1+p11*drho)/detp-(-q21*s1+q11*s2)/detq)
             wm2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wm6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wm4=drho-alphf*(wm1+wm7)-alphs*(wm3+wm5)
             !-----left amplitude
             drho=rhol-rholl
             dvx=rhobar*rhoavi*((rhol*vxl-rholl*vxll)-vxave*drho)
             dvy=rhobar*rhoavi*((rhol*vyl-rholl*vyll)-vyave*drho)
             dvz=rhobar*rhoavi*((rhol*vzl-rholl*vzll)-vzave*drho)
             dby=byl-byll
             dbz=bzl-bzll
             !
             t1=betay*dby+betaz*dbz
             t3=betaz*dby-betay*dbz
             !
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             !
             wl1=0.5d0*((p22*t1-p12*drho)/detp+(q22*s1-q12*s2)/detq)
             wl7=0.5d0*((p22*t1-p12*drho)/detp-(q22*s1-q12*s2)/detq)
             wl3=0.5d0*((-p21*t1+p11*drho)/detp+(-q21*s1+q11*s2)/detq)
             wl5=0.5d0*((-p21*t1+p11*drho)/detp-(-q21*s1+q11*s2)/detq)
             wl2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wl6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wl4=drho-alphf*(wl1+wl7)-alphs*(wl3+wl5)
             !-----right amplitude
             drho=rhorr-rhor
             dvx=rhobar*rhoavi*((rhorr*vxrr-rhor*vxr)-vxave*drho)
             dvy=rhobar*rhoavi*((rhorr*vyrr-rhor*vyr)-vyave*drho)
             dvz=rhobar*rhoavi*((rhorr*vzrr-rhor*vzr)-vzave*drho)
             dby=byrr-byr
             dbz=bzrr-bzr
             !
             t1=betay*dby+betaz*dbz
             t3=betaz*dby-betay*dbz
             !
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             !
             wr1=0.5d0*((p22*t1-p12*drho)/detp+(q22*s1-q12*s2)/detq)
             wr7=0.5d0*((p22*t1-p12*drho)/detp-(q22*s1-q12*s2)/detq)
             wr3=0.5d0*((-p21*t1+p11*drho)/detp+(-q21*s1+q11*s2)/detq)
             wr5=0.5d0*((-p21*t1+p11*drho)/detp-(-q21*s1+q11*s2)/detq)
             wr2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wr6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wr4=drho-alphf*(wr1+wr7)-alphs*(wr3+wr5)
             !------muscl in amplitude
             !------upwind calculation
             ea1=0.5d0-sign(0.5,vxbar+cfast)
             ea7=0.5d0-sign(0.5,vxbar-cfast)
             ea3=0.5d0-sign(0.5,vxbar+cslow)
             ea5=0.5d0-sign(0.5,vxbar-cslow)
             ea2=0.5d0-sign(0.5,vxbar+ca)
             ea6=0.5d0-sign(0.5,vxbar-ca)
             ea4=0.5d0-sign(0.5,vxbar)
             eb1=0.5d0-sign(0.5,-(vxbar+cfast))
             eb7=0.5d0-sign(0.5,-(vxbar-cfast))
             eb3=0.5d0-sign(0.5,-(vxbar+cslow))
             eb5=0.5d0-sign(0.5,-(vxbar-cslow))
             eb2=0.5d0-sign(0.5,-(vxbar+ca))
             eb6=0.5d0-sign(0.5,-(vxbar-ca))
             eb4=0.5d0-sign(0.5,-vxbar)
             !-------muscl in amplitude
             w1=wm1-FLMT(wm1,wl1)*eb1-FLMT(wm1,wr1)*ea1
             w2=wm2-FLMT(wm2,wl2)*eb2-FLMT(wm2,wr2)*ea2
             w3=wm3-FLMT(wm3,wl3)*eb3-FLMT(wm3,wr3)*ea3
             w4=wm4-FLMT(wm4,wl4)*eb4-FLMT(wm4,wr4)*ea4
             w5=wm5-FLMT(wm5,wl5)*eb5-FLMT(wm5,wr5)*ea5
             w6=wm6-FLMT(wm6,wl6)*eb6-FLMT(wm6,wr6)*ea6
             w7=wm7-FLMT(wm7,wl7)*eb7-FLMT(wm7,wr7)*ea7
             !-----components of the eigen vectors
             rpf1=alphf
             rpf2=alphf*(vxbar+cfast)
             rpf3=alphf*vybar-alphs*betay*ca*sgnbx
             rpf4=alphf*vzbar-alphs*betaz*ca*sgnbx
             rpf5=alphs*betay*cfast*sqrt(pi4/rhobar)
             rpf6=alphs*betaz*cfast*sqrt(pi4/rhobar)
             !
             rmf1=alphf
             rmf2=alphf*(vxbar-cfast)
             rmf3=alphf*vybar+alphs*betay*ca*sgnbx
             rmf4=alphf*vzbar+alphs*betaz*ca*sgnbx
             rmf5=rpf5
             rmf6=rpf6
             !
             rps1=alphs
             rps2=alphs*(vxbar+cslow)
             rps3=alphs*vybar+csx*sgnbx*alphf*betay
             rps4=alphs*vzbar+csx*sgnbx*alphf*betaz
             rps5=-sqrt(pi4/rhobar)*cs2x*alphf*betay/cfast
             rps6=-sqrt(pi4/rhobar)*cs2x*alphf*betaz/cfast
             !
             rms1=alphs
             rms2=alphs*(vxbar-cslow)
             rms3=alphs*vybar-csx*sgnbx*alphf*betay
             rms4=alphs*vzbar-csx*sgnbx*alphf*betaz
             rms5=rps5
             rms6=rps6
             !
             rpa1=0.0d0
             rpa2=0.0d0
             rpa3=-sgnbx*betaz
             rpa4=sgnbx*betay
             rpa5=sqrt(pi4/rhobar)*betaz
             rpa6=-sqrt(pi4/rhobar)*betay
             !
             rma1=0.0d0
             rma2=0.0d0
             rma3=-rpa3
             rma4=-rpa4
             rma5=rpa5
             rma6=rpa6
             !
             re1=1.0d0
             re2=vxbar
             re3=vybar
             re4=vzbar
             re5=0.0d0
             re6=0.0d0
             !
             !-----computation of f(i+1/2,j)
             f(i,j,k,MRHO)=(fl1+fr1 &
                  +el1*w1*rpf1 &
                  +el7*w7*rmf1 &
                  +el3*w3*rps1 &
                  +el5*w5*rms1 &
                  +el4*w4*re1)*5.0d-1
             f(i,j,k,MVX)=(fl2+fr2 &
                  +el1*w1*rpf2 &
                  +el7*w7*rmf2 &
                  +el3*w3*rps2 &
                  +el5*w5*rms2 &
                  +el4*w4*re2)*0.5d0
             f(i,j,k,MVY)=(fl3+fr3 &
                  +el1*w1*rpf3 &
                  +el7*w7*rmf3 &
                  +el3*w3*rps3 &
                  +el5*w5*rms3 &
                  +el2*w2*rpa3 &
                  +el6*w6*rma3 &
                  +el4*w4*re3)*0.5d0
             f(i,j,k,MVZ)=(fl4+fr4 &
                  +el1*w1*rpf4 &
                  +el7*w7*rmf4 &
                  +el3*w3*rps4 &
                  +el5*w5*rms4 &
                  +el2*w2*rpa4 &
                  +el6*w6*rma4 &
                  +el4*w4*re4)*0.5d0
             f(i,j,k,MBY)=(fl5+fr5 &
                  +el1*w1*rpf5 &
                  +el7*w7*rmf5 &
                  +el3*w3*rps5 &
                  +el5*w5*rms5 &
                  +el2*w2*rpa5 &
                  +el6*w6*rma5)*0.5d0
             f(i,j,k,MBZ)=(fl6+fr6 &
                  +el1*w1*rpf6 &
                  +el7*w7*rmf6 &
                  +el3*w3*rps6 &
                  +el5*w5*rms6 &
                  +el2*w2*rpa6 &
                  +el6*w6*rma6)*0.5d0
          enddo
       enddo
    enddo
    f1d(:,:,:,mcycle) = f(:,:,:,:)
  end subroutine flux
  !-----------------------------------------------------------------------
  ! source term due to div B
  !-----------------------------------------------------------------------
  subroutine source_b(f, w, dt, gid)
    use parameter
    use grid, only : get_dv, get_ds, get_level, get_up, Imingh, Jmingh, Kmingh, Mmin, &
         Imin, Imax, Jmin, Jmax, Kmin, Kmax, CellWidth, Lmin
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
          enddo
       enddo
    enddo
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
    real(kind=DBL_KIND) :: flag, csp, ca
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    level = get_level(id)
    h = CellWidth( :, level )
    rho => get_Ucomp(MRHO,id)
    vx => get_Ucomp(MVX,id)
    vy => get_Ucomp(MVY,id)
    vz => get_Ucomp(MVZ,id)
    bx => get_Ucomp(MBX,id)
    by => get_Ucomp(MBY,id)
    bz => get_Ucomp(MBZ,id)
#ifdef WITH_SELFGRAVITY
    gx => get_Ucomp(MGX,id)
    gy => get_Ucomp(MGY,id)
    gz => get_Ucomp(MGZ,id)
#endif !WITH_SELFGRAVITY
    x => get_Xp(id)
    y => get_Yp(id)
    z => get_Zp(id)
    dt = HUGE(dt)
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             GETCS(csp, rho(i,j,k))
             ca = sqrt( csp**2 + (bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)*pi4i/rho(i,j,k) )
#ifdef WITH_SELFGRAVITY
#define DTI_(V_,G_,H_) sqrt((V_)**2+2.d0*(G_)*(H_))/(H_)
             dt = min(dt, &
                  (CFL) / ( &
                  DTI_( abs(vx(i,j,k))+ca, abs(gx(i,j,k)), h(MX) ) + &
                  DTI_( abs(vy(i,j,k))+ca, abs(gy(i,j,k)), h(MY) ) + &
                  DTI_( abs(vz(i,j,k))+ca, abs(gz(i,j,k)), h(MZ) ) ))
#else   !WITH_SELFGRAVITY
             dt = min(dt, &
                  (CFL) / ( &
                  (abs(vx(i,j,k))+ca)/h(MX) + &
                  (abs(vy(i,j,k))+ca)/h(MY) + &
                  (abs(vz(i,j,k))+ca)/h(MZ) ) )
#endif  !WITH_SELFGRAVITY
          enddo
       enddo
    enddo
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
    u(:,:,:,MDB) = w(:,:,:,MDB)/dv
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
  subroutine get_cs(c, rho)
    real(kind=DBL_KIND),intent(OUT) :: c
    real(kind=DBL_KIND),intent(IN) :: rho
    real(kind=DBL_KIND) :: flag
    GETCS(c, rho)
  end subroutine get_cs
  !-----------------------------------------------------------------------
  ! return pressure (p) when density (rho) is given
  !-----------------------------------------------------------------------
  subroutine get_p(p, rho)
    real(kind=DBL_KIND),intent(OUT) :: p
    real(kind=DBL_KIND),intent(IN) :: rho
    real(kind=DBL_KIND) :: flag
    GETP(p, rho)
  end subroutine get_p
end module eos
subroutine eos_init
  use eos
  Kappa = Cs**2 * Rhocr ** (1.d0-Gamma) ! kappa of polytrope
end subroutine eos_init
