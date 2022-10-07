#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state (and system equations).
! Adiabatic MHD
!  Numerical flux is based on Hanawa & Fukuda (19...)
!-------------------------------------------------------------------------
module eos
  implicit none
  real(kind=DBL_KIND),parameter :: Gamma = 5.d0/3.d0
!!$  real(kind=DBL_KIND),parameter :: Gamma = 1.4d0
!!$  real(kind=DBL_KIND),parameter :: Gamma = 2.d0
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
#define ELMOD(EL,ELBAR,EL1, EL2) \
  epsel=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-epsel) ;\
  x2=1.d0-x1 ;\
  EL= -(x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(epsel+x1)+epsel))

#define FLMT(x,y) max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x))
! #define FLMT(x,y) (x)*0.d0+(y)*0.d0
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
!!$    real(kind=DBL_KIND),parameter :: eps = 1.0d-6
    real(kind=DBL_KIND),parameter :: eps = 1.0d-10
    real(kind=DBL_KIND),parameter :: eps2=eps*2.d0
    real(kind=DBL_KIND) :: &
         rhol,vxl,vyl,vzl,pl,hl,el,bxl,byl,bzl,psil,cl, &
         rhor,vxr,vyr,vzr,pr,hr,er,bxr,byr,bzr,psir,cr, &
         rholl,vxll,vyll,vzll,pll,hll,ell,bxll,byll,bzll,psill,cll, &
         rhorr,vxrr,vyrr,vzrr,prr,hrr,err,bxrr,byrr,bzrr,psirr,crr, &
         gm1, gm2, gm21, &
         alphf, alphs, astar2, betay, betaz, bxave, bxbar, byave,     &
         bybar, bybz, bzave, bzbar, ca, ca2, cfa, cfast, cfast2,      &
         cfca, cfcs, cs2x, csca, cslow, cslow2, csx, dby, dbz, detp,  &
         detq, drho, dvx, dvy, dvz, ea1, ea2, ea3, ea4, ea5, ea6,     &
         ea7, eb1, eb2, eb3, eb4, eb5, eb6, eb7, eeps, el1, el2, el3, &
         el4, el5, el6, el7, fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8,  &
         fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8, p11, p12, p21, p22,  &
         q11, q12, q21, q22, re1, re2, re3, re4, re5, re6, rhoave,    &
         rhoavi, rhobar, rma1, rma2, rma3, rma4, rma5, rma6, rmf1,    &
         rmf2, rmf3, rmf4, rmf5, rmf6, rms1, rms2, rms3, rms4, rms5,  &
         rms6, rpa1, rpa2, rpa3, rpa4, rpa5, rpa6, rpf1, rpf2, rpf3,  &
         rpf4, rpf5, rpf6, rps1, rps2, rps3, rps4, rps5, rps6, s1,    &
         s2, s3, sgnbx, sgr, sgr2, sp, sp2, sr0, sr1, sri, t1, t3,    &
         vxave, vxbar, vyave, vybar, vzave, vzbar, w1, w2, w3, w4,    &
         w5, w6, w7, wl1, wl2, wl3, wl4, wl5, wl6, wl7, wm1, wm2,     &
         wm3, wm4, wm5, wm6, wm7, wr1, wr2, wr3, wr4, wr5, wr6, wr7,  &
         cs2,csi,                                                     &
         hbar,vx2ave,vy2ave,vz2ave,enl,enll,enr,enrr,                 &
         cbr2,cs,rpf7,rpf8,rmf7,rmf8,rps7,rps8,rms7,rms8,rpa7,        &
         rpa8,rma7,rma8,re7,re8,dbx,delb2,dp,t2, bxm, psim, bxlm, bxrm, psilm, psirm
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND) :: epsel, x1, x2
    integer,dimension(Mmin:Mmax) :: mcycle

    mcycle = cyclecomp( ncrd )
    q(:,:,:,:) = u(:,:,:,mcycle)

#ifdef SINGLE_STEP
    h = CellWidth( :, LevelMax )
#else  !SINGLE_STEP
    h = CellWidth( :, Lmin )
#endif !SINGLE_STEP
    Ch = (CFL) * minval(h) / Dtime(Lmin) / 3.d0

    gm1 = Gamma - 1.d0
    gm2 = Gamma - 2.d0
    gm21 = gm2/gm1

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
#ifdef EMULATE_2DIM
    ks = Kmin
    ke = Kmin
#endif !EMULATE_2DIM


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
             pll =   q(i-io,j-jo,k-ko,MP)
             psill = q(i-io,j-jo,k-ko,MDB)
             rhol =  q(i,j,k,MRHO)
             vxl =   q(i,j,k,MVX)
             vyl =   q(i,j,k,MVY)
             vzl =   q(i,j,k,MVZ)
             bxl =   q(i,j,k,MBX)
             byl =   q(i,j,k,MBY)
             bzl =   q(i,j,k,MBZ)
             pl =    q(i,j,k,MP)
             psil =  q(i,j,k,MDB)
             rhor =  q(i+io,j+jo,k+ko,MRHO)
             vxr =   q(i+io,j+jo,k+ko,MVX)
             vyr =   q(i+io,j+jo,k+ko,MVY)
             vzr =   q(i+io,j+jo,k+ko,MVZ)
             bxr =   q(i+io,j+jo,k+ko,MBX)
             byr =   q(i+io,j+jo,k+ko,MBY)
             bzr =   q(i+io,j+jo,k+ko,MBZ)
             pr =    q(i+io,j+jo,k+ko,MP)
             psir =  q(i+io,j+jo,k+ko,MDB)
             rhorr = q(i+io2,j+jo2,k+ko2,MRHO)
             vxrr =  q(i+io2,j+jo2,k+ko2,MVX)
             vyrr =  q(i+io2,j+jo2,k+ko2,MVY)
             vzrr =  q(i+io2,j+jo2,k+ko2,MVZ)
             bxrr =  q(i+io2,j+jo2,k+ko2,MBX)
             byrr =  q(i+io2,j+jo2,k+ko2,MBY)
             bzrr =  q(i+io2,j+jo2,k+ko2,MBZ)
             prr =   q(i+io2,j+jo2,k+ko2,MP)
             psirr = q(i+io2,j+jo2,k+ko2,MDB)
             ! ------------------------------------------
             ! add-on for div B free (Dedner et al. 2002)
             ! ------------------------------------------
             bxlm =bxl +0.5*FLMT(bxr -bxl, bxl-bxll)
             bxrm =bxr -0.5*FLMT(bxr -bxl, bxrr-bxr)
             psilm=psil+0.5*FLMT(psir-psil,psil-psill)
             psirm=psir-0.5*FLMT(psir-psil,psirr-psir)
             bxm  = 0.5d0*(bxlm+bxrm - (psirm-psilm)/Ch)
             psim = 0.5d0*(psilm+psirm - (bxrm-bxlm)*Ch)
             f(i,j,k,MBX) = psim
             f(i,j,k,MDB) = bxm*Ch**2
             bxbar=bxm          ! value at interface
             bxave=bxbar

             ! -------------
             ! Roe variable
             ! -------------
             sr0=sqrt(rhol)
             sr1=sqrt(rhor)
             sri=1.0d0/(sr0+sr1)
             rhobar=sr0*sr1
             vxbar=(sr0*vxl+sr1*vxr)*sri
             vybar=(sr0*vyl+sr1*vyr)*sri
             vzbar=(sr0*vzl+sr1*vzr)*sri
!!$             bxbar=(sr0*bxr+sr1*bxl)*sri
             bybar=(sr0*byr+sr1*byl)*sri
             bzbar=(sr0*bzr+sr1*bzl)*sri
             hl=0.5d0*(vxl**2+vyl**2+vzl**2)+gamma*pl/(gm1*rhol) &
                  +(bxbar**2+byl**2+bzl**2)/(pi4*rhol)
             hr=0.5d0*(vxr**2+vyr**2+vzr**2)+gamma*pr/(gm1*rhor) &
                  +(bxbar**2+byr**2+bzr**2)/(pi4*rhor)
             hbar=(sr0*hl+sr1*hr)*sri
             rhoave=0.5d0*(rhol+rhor)
             rhoavi=1.0d0/rhoave
             vxave=0.5d0*(vxl+vxr)
             vyave=0.5d0*(vyl+vyr)
             vzave=0.5d0*(vzl+vzr)
             vx2ave=0.5d0*(vxl**2+vxr**2)
             vy2ave=0.5d0*(vyl**2+vyr**2)
             vz2ave=0.5d0*(vzl**2+vzr**2)
!!$             bxave=0.5d0*(bxl+bxr)
             byave=0.5d0*(byl+byr)
             bzave=0.5d0*(bzl+bzr)
             enl=0.5d0*rhol*(vxl**2+vyl**2+vzl**2)+pl/gm1 &
                  +pi8i*(bxbar**2+byl**2+bzl**2)
             enll=0.5d0*rholl*(vxll**2+vyll**2+vzll**2)+pll/gm1 &
                  +pi8i*(bxbar**2+byll**2+bzll**2)
             enr=0.5d0*rhor*(vxr**2+vyr**2+vzr**2)+pr/gm1 &
                  +pi8i*(bxbar**2+byr**2+bzr**2)
             enrr=0.5d0*rhorr*(vxrr**2+vyrr**2+vzrr**2)+prr/gm1 &
                  +pi8i*(bxbar**2+byrr**2+bzrr**2)

             !-----characteristic speed
             delb2=gm21*((byr-byl)**2+(bzr-bzl)**2)*sri**2*pi8i
             cs2=gm1*(hbar-0.5d0*(vxbar**2+vybar**2+vzbar**2) &
                  -delb2-(bxbar**2+bybar**2+bzbar**2)*pi4i/rhobar)
             astar2=gm1*(hbar-0.5d0*(vxbar**2+vybar**2+vzbar**2)-delb2) &
                  -gm2*(bxbar**2+bybar**2+bzbar**2)*pi4i/rhobar
             ca2=bxbar**2/(pi4*rhobar)
             !      cfast2=0.5d0*(astar2+sqrt(astar2**2-4.0d0*cs2*ca2))
             cbr2=(bybar**2+bzbar**2)*pi4i/rhobar
             cfast2=0.5d0*(astar2+sqrt(cbr2*(astar2+cs2+ca2)+(cs2-ca2)**2))
             cslow2=cs2*ca2/cfast2
             cfast=sqrt(cfast2)
             cslow=sqrt(cslow2)
             ca=sqrt(ca2)
             cs=sqrt(cs2)


             !----- for singular points
             sgr=bybar**2+bzbar**2-eps ! 2007/1/18
!!$             sgr=bybar**2+bzbar**2 -eps*bxbar
             sp=0.5d0+sign(0.5d0,sgr)
             betay=sp*bybar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp))+sqrt(0.5d0)*(1.0d0-sp)
             betaz=sp*bzbar*sqrt(1.0d0/(bybar**2+bzbar**2+1.0d0-sp))+sqrt(0.5d0)*(1.0d0-sp)
             sgr2=(bybar**2+bzbar**2)/(pi4*rhobar)+abs(ca2-cs2)-eps2 !2007/1/18
!!$             sgr2=(bybar**2+bzbar**2)/(pi4*rhobar)+abs(ca2-cs2)-eps2*bxbar
             sp2=0.5d0+sign(0.5d0,sgr2)
             cfca=max(0.0d0,cfast2-ca2)
             cfcs=max(0.0d0,cfast2-cslow2)
             cfa=max(0.0d0,cfast2-cs2)
             alphf=sp2*sqrt(cfca/(cfcs+1.0d0-sp2))+1.0d0-sp2
             alphs=sp2*sqrt((cfa)/(cfcs+1.0d0-sp2))
             sgnbx=sign(1.0d0,bxbar)
             !----- eigen value & entropy condition
             ELMOD(el1,vxbar+cfast, vxl+cfast, vxr+cfast)
             ELMOD(el7,vxbar-cfast, vxl-cfast, vxr-cfast)
             ELMOD(el3,vxbar+cslow, vxl+cslow, vxr+cslow)
             ELMOD(el5,vxbar-cslow, vxl-cslow, vxr-cslow)
             ELMOD(el2,vxbar+ca,    vxl+ca,    vxr+ca)
             ELMOD(el6,vxbar-ca,    vxl-ca,    vxr-ca)
             ELMOD(el4,vxbar,       vxl,       vxr)
             !----- components of the eigen vectors
             !----- fast mode +(w1)
             rpf1=alphf
             rpf2=alphf*(vxbar+cfast)
             rpf3=alphf*vybar-alphs*betay*ca*sgnbx
             rpf4=alphf*vzbar-alphs*betaz*ca*sgnbx
             rpf6=alphs*betay*cfast*sqrt(pi4/rhobar)
             rpf7=alphs*betaz*cfast*sqrt(pi4/rhobar)
             rpf8=alphf*(0.5d0*(vxbar**2+vybar**2+vzbar**2)+delb2+cfast*vxbar &
                  +cfast2/gm1 &
                  +(cfast2-cs2)*gm21) &
                  -alphs*ca*(betay*vybar+betaz*vzbar)*sgnbx
             !----- fast mode -(w7)
             rmf1=alphf
             rmf2=alphf*(vxbar-cfast)
             rmf3=alphf*vybar+alphs*betay*ca*sgnbx
             rmf4=alphf*vzbar+alphs*betaz*ca*sgnbx
             rmf6=rpf6
             rmf7=rpf7
             rmf8=alphf*(0.5d0*(vxbar**2+vybar**2+vzbar**2) &
                  +delb2-cfast*vxbar+cfast2/gm1 &
                  +(cfast2-cs2)*gm21) &
                  +alphs*ca*(betay*vybar+betaz*vzbar)*sgnbx
             !---- slow mode +(w3)
             rps1=alphs
             rps2=alphs*(vxbar+cslow)
             rps3=alphs*vybar+cs*sgnbx*alphf*betay
             rps4=alphs*vzbar+cs*sgnbx*alphf*betaz
             rps6=-sqrt(pi4/rhobar)*cs2*alphf*betay/cfast
             rps7=-sqrt(pi4/rhobar)*cs2*alphf*betaz/cfast
             rps8=alphs*(0.5d0*(vxbar**2+vybar**2+vzbar**2) &
                  +delb2+cslow*vxbar+cslow2/gm1 &
                  +(cslow2-cs2)*gm21) &
                  +alphf*cs*(betay*vybar+betaz*vzbar)*sgnbx
             !----- slow mode -(w5)
             rms1=alphs
             rms2=alphs*(vxbar-cslow)
             rms3=alphs*vybar-cs*sgnbx*alphf*betay
             rms4=alphs*vzbar-cs*sgnbx*alphf*betaz
             rms6=rps6
             rms7=rps7
             rms8=alphs*(0.5d0*(vxbar**2+vybar**2+vzbar**2) &
                  +delb2-cslow*vxbar+cslow2/gm1 &
                  +(cslow2-cs2)*gm21) &
                  -alphf*cs*(betay*vybar+betaz*vzbar)*sgnbx
             !----- alfven mode +(w2)
             rpa1=0.0d0
             rpa2=0.0d0
             rpa3=-sgnbx*betaz
             rpa4=sgnbx*betay
             rpa6=sqrt(pi4/rhobar)*betaz
             rpa7=-sqrt(pi4/rhobar)*betay
             rpa8=-(betaz*vybar-betay*vzbar)*sgnbx
             !----- alfven mode -(w6)
             rma1=0.0d0
             rma2=0.0d0
             rma3=-rpa3
             rma4=-rpa4
             rma6=rpa6
             rma7=rpa7
             rma8=-rpa8
             !----- entropy wave (w4)
             re1=1.0d0
             re2=vxbar
             re3=vybar
             re4=vzbar
             re6=0.0d0
             re7=0.0d0
             re8=0.5d0*(vxbar**2+vybar**2+vzbar**2)+delb2
             !----- calc amplitude
             p11=alphs*cfast*sqrt(pi4/rhobar)
             p12=-alphf*cs2/cfast*sqrt(pi4/rhobar)
             p21=alphf*(cfast2-cs2*gm21)
             p22=alphs*(cslow2-cs2*gm21)
             q11=alphf*cfast
             q12=alphs*cslow
             q21=-alphs*ca*sgnbx
             q22=alphf*cs*sgnbx
             detp=p11*p22-p12*p21
             detq=q11*q22-q12*q21
             !----- medium amplitude
             drho=rhor-rhol
             dvx=rhobar*rhoavi*((rhor*vxr-rhol*vxl)-vxave*(rhor-rhol))
             dvy=rhobar*rhoavi*((rhor*vyr-rhol*vyl)-vyave*(rhor-rhol))
             dvz=rhobar*rhoavi*((rhor*vzr-rhol*vzl)-vzave*(rhor-rhol))
             dbx=0.0d0
             dby=byr-byl
             dbz=bzr-bzl
             dp=gm1*(enr-enl-pi4i*(bxbar*dbx+byave*dby+bzave*dbz) &
                  -vxave*(rhor*vxr-rhol*vxl)-vyave*(rhor*vyr-rhol*vyl) &
                  -vzave*(rhor*vzr-rhol*vzl) &
                  +(rhor-rhol)*(vxave**2+vyave**2+vzave**2 &
                  -0.5d0*(vx2ave+vy2ave+vz2ave)))
             t1=betay*dby+betaz*dbz
             t2=(dp+(byave*dby+bzave*dbz)*pi4i &
                  +gm2*(bybar*dby+bzbar*dbz)*pi4i)/gm1
             t3=betaz*dby-betay*dbz
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             wm1=0.5d0*((p22*t1-p12*t2)/detp+(q22*s1-q12*s2)/detq)
             wm7=0.5d0*((p22*t1-p12*t2)/detp-(q22*s1-q12*s2)/detq)
             wm3=0.5d0*((-p21*t1+p11*t2)/detp+(-q21*s1+q11*s2)/detq)
             wm5=0.5d0*((-p21*t1+p11*t2)/detp-(-q21*s1+q11*s2)/detq)
             wm2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wm6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wm4=drho-alphf*(wm1+wm7)-alphs*(wm3+wm5)
             !----- left amplitude
             drho=rhol-rholl
             dvx=rhobar*rhoavi*((rhol*vxl-rholl*vxll)-vxave*(rhol-rholl))
             dvy=rhobar*rhoavi*((rhol*vyl-rholl*vyll)-vyave*(rhol-rholl))
             dvz=rhobar*rhoavi*((rhol*vzl-rholl*vzll)-vzave*(rhol-rholl))
             dbx=0.0d0
             dby=byl-byll
             dbz=bzl-bzll
             dp=gm1*(enl-enll-pi4i*(bxbar*dbx+byave*dby+bzave*dbz) &
                  -vxave*(rhol*vxl-rholl*vxll)-vyave*(rhol*vyl-rholl*vyll) &
                  -vzave*(rhol*vzl-rholl*vzll) &
                  +(rhol-rholl)*(vxave**2+vyave**2+vzave**2 &
                  -0.5d0*(vx2ave+vy2ave+vz2ave)))
             t1=betay*dby+betaz*dbz
             t2=(dp+(byave*dby+bzave*dbz)*pi4i &
                  +gm2*(bybar*dby+bzbar*dbz)*pi4i)/gm1
             t3=betaz*dby-betay*dbz
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             wl1=0.5d0*((p22*t1-p12*t2)/detp+(q22*s1-q12*s2)/detq)
             wl7=0.5d0*((p22*t1-p12*t2)/detp-(q22*s1-q12*s2)/detq)
             wl3=0.5d0*((-p21*t1+p11*t2)/detp+(-q21*s1+q11*s2)/detq)
             wl5=0.5d0*((-p21*t1+p11*t2)/detp-(-q21*s1+q11*s2)/detq)
             wl2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wl6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wl4=drho-alphf*(wl1+wl7)-alphs*(wl3+wl5)
             !----- right amplitude
             drho=rhorr-rhor
             dvx=rhobar*rhoavi*((rhorr*vxrr-rhor*vxr)-vxave*(rhorr-rhor))
             dvy=rhobar*rhoavi*((rhorr*vyrr-rhor*vyr)-vyave*(rhorr-rhor))
             dvz=rhobar*rhoavi*((rhorr*vzrr-rhor*vzr)-vzave*(rhorr-rhor))
             dbx=0.0d0
             dby=byrr-byr
             dbz=bzrr-bzr
             dp=gm1*(enrr-enr-pi4i*(bxbar*dbx+byave*dby+bzave*dbz) &
                  -vxave*(rhorr*vxrr-rhor*vxr)-vyave*(rhorr*vyrr-rhor*vyr) &
                  -vzave*(rhorr*vzrr-rhor*vzr) &
                  +(rhorr-rhor)*(vxave**2+vyave**2+vzave**2 &
                  -0.5d0*(vx2ave+vy2ave+vz2ave)))
             t1=betay*dby+betaz*dbz
             t2=(dp+(byave*dby+bzave*dbz)*pi4i &
                  +gm2*(bybar*dby+bzbar*dbz)*pi4i)/gm1
             t3=betaz*dby-betay*dbz
             s1=dvx
             s2=betay*dvy+betaz*dvz
             s3=betaz*dvy-betay*dvz
             wr1=0.5d0*((p22*t1-p12*t2)/detp+(q22*s1-q12*s2)/detq)
             wr7=0.5d0*((p22*t1-p12*t2)/detp-(q22*s1-q12*s2)/detq)
             wr3=0.5d0*((-p21*t1+p11*t2)/detp+(-q21*s1+q11*s2)/detq)
             wr5=0.5d0*((-p21*t1+p11*t2)/detp-(-q21*s1+q11*s2)/detq)
             wr2=0.5d0*(sqrt(rhobar*pi4i)*t3-sgnbx*s3)
             wr6=0.5d0*(sqrt(rhobar*pi4i)*t3+sgnbx*s3)
             wr4=drho-alphf*(wr1+wr7)-alphs*(wr3+wr5)
             !----- muscl in amplitude
             !----- upwind calculation
             ea1=0.5d0-sign(0.5d0,vxbar+cfast)
             ea7=0.5d0-sign(0.5d0,vxbar-cfast)
             ea3=0.5d0-sign(0.5d0,vxbar+cslow)
             ea5=0.5d0-sign(0.5d0,vxbar-cslow)
             ea2=0.5d0-sign(0.5d0,vxbar+ca)
             ea6=0.5d0-sign(0.5d0,vxbar-ca)
             ea4=0.5d0-sign(0.5d0,vxbar)
             eb1=0.5d0-sign(0.5d0,-(vxbar+cfast))
             eb7=0.5d0-sign(0.5d0,-(vxbar-cfast))
             eb3=0.5d0-sign(0.5d0,-(vxbar+cslow))
             eb5=0.5d0-sign(0.5d0,-(vxbar-cslow))
             eb2=0.5d0-sign(0.5d0,-(vxbar+ca))
             eb6=0.5d0-sign(0.5d0,-(vxbar-ca))
             eb4=0.5d0-sign(0.5d0,-(vxbar))
             !----- muscl in amplitude
             w1=wm1-FLMT(wm1,wl1)*eb1-FLMT(wm1,wr1)*ea1
             w2=wm2-FLMT(wm2,wl2)*eb2-FLMT(wm2,wr2)*ea2
             w3=wm3-FLMT(wm3,wl3)*eb3-FLMT(wm3,wr3)*ea3
             w4=wm4-FLMT(wm4,wl4)*eb4-FLMT(wm4,wr4)*ea4
             w5=wm5-FLMT(wm5,wl5)*eb5-FLMT(wm5,wr5)*ea5
             w6=wm6-FLMT(wm6,wl6)*eb6-FLMT(wm6,wr6)*ea6
             w7=wm7-FLMT(wm7,wl7)*eb7-FLMT(wm7,wr7)*ea7
             !----- flux_l
             fl1=rhol*vxl
             fl2=fl1*vxl+pl+(-bxbar**2+byl**2+bzl**2)*pi8i
             fl3=fl1*vyl-bxbar*byl*pi4i
             fl4=fl1*vzl-bxbar*bzl*pi4i
             fl5=0.0d0
             fl6=vxl*byl-vyl*bxbar
             fl7=vxl*bzl-vzl*bxbar
             fl8=fl1*hl-bxbar*(bxbar*vxl+byl*vyl+bzl*vzl)*pi4i
             !----- flux_r
             fr1=rhor*vxr
             fr2=fr1*vxr+pr+(-bxbar**2+byr**2+bzr**2)*pi8i
             fr3=fr1*vyr-bxbar*byr*pi4i
             fr4=fr1*vzr-bxbar*bzr*pi4i
             fr5=0.0d0
             fr6=vxr*byr-vyr*bxbar
             fr7=vxr*bzr-vzr*bxbar
             fr8=fr1*hr-bxbar*(bxbar*vxr+byr*vyr+bzr*vzr)*pi4i
             !----- computation of f(i+1/2,j)
             f(i,j,k,MRHO) =0.5d0*(fl1+fr1 &
                  +el1*w1*rpf1+el7*w7*rmf1+el3*w3*rps1+el5*w5*rms1 &
                  +el4*w4*re1)
             f(i,j,k,MVX) =0.5d0*(fl2+fr2 &
                  +el1*w1*rpf2+el7*w7*rmf2+el3*w3*rps2+el5*w5*rms2 &
                  +el4*w4*re2)
             f(i,j,k,MVY) =0.5d0*(fl3+fr3 &
                  +el1*w1*rpf3+el7*w7*rmf3+el3*w3*rps3+el5*w5*rms3 &
                  +el2*w2*rpa3+el6*w6*rma3+el4*w4*re3)
             f(i,j,k,MVZ) =0.5d0*(fl4+fr4 &
                  +el1*w1*rpf4+el7*w7*rmf4+el3*w3*rps4+el5*w5*rms4 &
                  +el2*w2*rpa4+el6*w6*rma4+el4*w4*re4)
             f(i,j,k,MBY) =0.5d0*(fl6+fr6 &
                  +el1*w1*rpf6+el7*w7*rmf6+el3*w3*rps6+el5*w5*rms6 &
                  +el2*w2*rpa6+el6*w6*rma6)
             f(i,j,k,MBZ) =0.5d0*(fl7+fr7 &
                  +el1*w1*rpf7+el7*w7*rmf7+el3*w3*rps7+el5*w5*rms7 &
                  +el2*w2*rpa7+el6*w6*rma7)
             f(i,j,k,MP) =0.5d0*(fl8+fr8 &
                  +el1*w1*rpf8+el7*w7*rmf8+el3*w3*rps8+el5*w5*rms8 &
                  +el2*w2*rpa8+el6*w6*rma8+el4*w4*re8)
          enddo
       enddo
    enddo
#ifdef EMULATE_2DIM
    do k = Kmin-1, Kmax
       if (k /= Kmin) f(:,:,k,:) = f(:,:,Kmin,:)
    enddo
#endif !EMULATE_2DIM
    f1d(:,:,:,mcycle) = f(:,:,:,:)
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
