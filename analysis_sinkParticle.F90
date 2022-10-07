#include "config.h"
#include "barotropic.h"
!-------------------------------------------------------------------------
! Module for analysis for using sink particles.
! This routine should be called in at the end of subroutine
! step_all_level of multi_timestep.F90
! -------------------------------------------------------------------------
module analysis
  use mpilib
  use grid
  implicit none
  private
  integer,parameter :: IntervalStep = 10
  real(kind=DBL_KIND),save :: RhoMax, RhoMin
  real(kind=DBL_KIND),save,dimension(Lmin:Lmax) :: RhoMaxLev, RhoMinLev
  !
  public :: RhoMax, RhoMin
  public :: analysis_keyparam
contains
  !-----------------------------------------------------------------------
  ! Front end of analysis
  !-----------------------------------------------------------------------
  subroutine analysis_keyparam
    use eos, only : Rhocr
    use unit
    use sinkParticle, only : sp_getRhocr
    real(kind=DBL_KIND) :: rhoCritical
    real(kind=DBL_KIND),parameter :: LargeRadius = 1.D10
    call an_rhomax              ! find RhoMax and RhoMin
    if (.not. bool_analysis()) return ! skip

!     rhoCritical = Rhocr
!     if (an_skip(rhoCritical)) call an_sumvalues(rhoCritical, 'logfileRhocr')
!     rhoCritical = 1.d0
!     if (an_skip(rhoCritical)) call an_sumvalues(rhoCritical, 'logfile1d0')
!     rhoCritical = sp_getRhocr()*0.1d0
!     if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile01rhosp')

    rhoCritical = 0.d0          ! whole
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile0')

    rhoCritical = 1.d3          ! 1.d3
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile1d3')

    rhoCritical = 1.d5          ! 1.d5 (Rhocr)
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile1d5')

    rhoCritical = 1.d6          ! 1.d6
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile1d6')

    rhoCritical = 1.d8          ! 1.d8
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile1d8')

    rhoCritical = RhoMax*0.01d0 ! 1.d8
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile001rhomax')

    rhoCritical = RhoMax*0.1d0  ! 1.d9
    if (an_skip(rhoCritical)) call an_sinkAround(rhoCritical, LargeRadius, 'logfile01rhomax')

  end subroutine analysis_keyparam
  !-----------------------------------------------------------------------
  ! return ture when analysis shoudl be executed.
  !-----------------------------------------------------------------------
  function bool_analysis() result(bool)
    use grid, only : level_sync, Step, Lmin
    logical :: bool
    bool = .FALSE.
    if ( level_sync() == Lmin .and. mod(Step(Lmin), IntervalStep) == 0 ) bool = .TRUE.
  end function bool_analysis
  !-----------------------------------------------------------------------
  ! find rhomax
  !-----------------------------------------------------------------------
  subroutine an_rhomax
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    real(kind=DBL_KIND) :: rhmx, rhmn
    integer :: n, level, gid
    myrank = get_myrank()
    ! find rhomax
    do level = Lmin, LevelMax
       rhmx = TINY(rhmx)
       rhmn = HUGE(rhmn)
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( MRHO, gid )
          rhmx = max(rhmx, MAXVAL(rho, mask=GridMask))
          rhmn = min(rhmn, MINVAL(rho, mask=GridMask))
       enddo
       call mpi_allreduce(rhmx, RhoMaxLev(level), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(rhmn, RhoMinLev(level), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    enddo
    RhoMax = MAXVAL(RhoMaxLev(Lmin:LevelMax))
    RhoMin = MINVAL(RhoMinLev(Lmin:LevelMax))

  end subroutine an_rhomax
  !-----------------------------------------------------------------------
  ! check critical density
  !-----------------------------------------------------------------------
  function an_skip(rhoCritical) result(bool)
    real(kind=DBL_KIND),intent(IN) :: rhoCritical
    logical :: bool
    bool = .false.
    if ( rhoCritical > RhoMax ) return
    if ( level_sync() > level_rhocr(rhoCritical) ) return
    bool = .true.
  end function an_skip
  ! -----------------------------------------------------------------
  ! rhoCritical を内包するグリッド
  ! -----------------------------------------------------------------
  function level_rhocr(rhoCritical) result(levelcr)
    real(kind=DBL_KIND),intent(IN) :: rhoCritical
    integer :: levelcr, level
    do level = Lmin, LevelMax
       if ( rhoCritical <= RhoMaxLev(level)) then
          levelcr = level
          exit
       end if
    end do
  end function level_rhocr
  !-----------------------------------------------------------------------
  ! get key parameters
  ! rhomax, volume, mass, I, Ekin, Ethermal, Emag, Egrav
  ! bin : all, rho > rhomax/2, rho > rhomax/10, rho > rhocr
  !-----------------------------------------------------------------------
  subroutine an_sumvalues(rhoCritical, logfile)
    use io_util, only : read_env
    use string, only : CHARLEN, concat
    use eos, only : Rhocr, Cs, Kappa, Gamma
    use parameter
    real(kind=DBL_KIND),intent(IN) :: rhoCritical
    character(len=*),intent(IN) :: logfile
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    logical,dimension(ARRAYSIZE3(GridMask)) :: mask
    real(kind=DBL_KIND),dimension(7) :: bufs, bufr
    real(kind=DBL_KIND) :: vol, mass, eKin, eMag, eGrav, inertia, eTh, dv, p
    integer :: n, level, gid, i, j, k, istat
    integer,parameter ::  FH=11
    character(len=CHARLEN) :: dir, fn

    real(kind=DBL_KIND) :: flag

    myrank = get_myrank()
    vol = 0.d0
    mass = 0.d0
    eKin = 0.d0
    eMag = 0.d0
    eGrav = 0.d0
    inertia = 0.d0
    eTh = 0.d0
    do level = Lmin, LevelMax
       dv = get_dv( level )
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( MRHO, gid )
          vx => get_Ucomp( MVX, gid )
          vy => get_Ucomp( MVY, gid )
          vz => get_Ucomp( MVZ, gid )
#ifdef MHD
          bx => get_Ucomp( MBX, gid )
          by => get_Ucomp( MBY, gid )
          bz => get_Ucomp( MBZ, gid )
#endif !MHD
          gx => get_Ucomp( MGX, gid )
          gy => get_Ucomp( MGY, gid )
          gz => get_Ucomp( MGZ, gid )
          x => get_Xp( gid )
          y => get_Yp( gid )
          z => get_Zp( gid )
          mask = GridMask .and. ( rho >= rhoCritical )
          if ( count(mask) == 0 ) cycle
          vol  = vol  + COUNT( mask ) * dv
          mass = mass + SUM(rho, mask=mask) * dv
          eKin = eKin + SUM(rho*(vx**2+vy**2+vz**2), mask=mask) * dv
#ifdef MHD
          eMag = eMag + SUM(bx**2+by**2+bz**2, mask=mask) * dv
#endif !MHD
          eGrav = eGrav + SUM(gx**2+gy**2+gz**2, mask=mask) * dv
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   if ( mask(i,j,k) ) then
                      GETP(p, rho(i,j,k))
                      eTh = eTh + p * dv
                      inertia = inertia + rho(i,j,k) * (x(i)**2+y(j)**2+z(k)**2) * dv
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    bufs = (/ vol, mass, eKin, eMag, eGrav, eTh, inertia /)
    call mpi_allreduce(bufs, bufr, size(bufs), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    vol = bufr(1)
    mass = bufr(2)
    eKin = bufr(3) * 0.5d0
    eMag = bufr(4) * pi8i
    eGrav = bufr(5) * 0.5d0
    eTh = bufr(6)
    inertia = bufr(7)
    if ( get_myrank() == PRIMARY_RANK ) then
       call read_env('DIR', dir)
       fn = concat(dir, logfile)
       open(FH, file=fn, position='APPEND')
       write(FH,'(I12, 9(1PE12.5))') Step(Lmin), Time(Lmin), RhoMax, vol, mass, eKin, eMag, eGrav, eTh, inertia
       call flush(FH)
       close(FH)
    endif
  end subroutine an_sumvalues
  !-----------------------------------------------------------------------
  ! estimate parameters for every sink particles
  !-----------------------------------------------------------------------
  subroutine an_sinkAround(rhoCritical, radius, logfile)
    use sinkParticle
    real(kind=DBL_KIND),intent(IN) :: rhoCritical, radius
    character(len=*),intent(IN) :: logfile
    integer :: nparticle, np, n
    integer,dimension(:),allocatable :: pid
    real(kind=DBL_KIND),dimension(:),allocatable :: pmass
    real(kind=DBL_KIND),dimension(:,:),allocatable :: pr, pv, ps
    nparticle = sp_getNparticle()
    if (nparticle > 0) then     ! sink particles exist
       allocate(pid(nparticle), pmass(nparticle), pr(MX:MZ, nparticle), pv(MX:MZ, nparticle), ps(MX:MZ, nparticle))
       call sp_sinkdata2array(np, pmass, pr=pr, pv=pv, ps=ps, pid=pid)
       if (np /= nparticle) &   ! error check
            print *, '*** error in an_sinkAround: np is not equal to nparticle', np, nparticle
    else                        ! no sink particle
       nparticle = 1
       allocate(pid(nparticle), pmass(nparticle), pr(MX:MZ, nparticle), pv(MX:MZ, nparticle), ps(MX:MZ, nparticle))
       pmass = 0.d0             ! fill dummy values
       pr = 0.d0
       pv = 0.d0
       ps = 0.d0
       pid = -1                 ! mark for no sink particle
    endif
    do n = lbound(pid,1), ubound(pid,1)
       call an_sinkAround_eachParticle(rhoCritical, radius, logfile, pmass(n), pr(:,n), pv(:,n), ps(:,n), pid(n))
    enddo
    deallocate(pid, pmass, pr, pv, ps)
  end subroutine an_sinkAround
  !-----------------------------------------------------------------------
  ! estimate parameters for each sink particle
  !-----------------------------------------------------------------------
  subroutine an_sinkAround_eachParticle(rhoCritical, radius, logfile, pmass, pr, pv, ps, pid)
    use io_util, only : read_env
    use string, only : CHARLEN, concat
    use eos, only : Rhocr, Cs, Kappa, Gamma
    use parameter
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    real(kind=DBL_KIND),intent(IN) :: rhoCritical, radius
    character(len=*),intent(IN) :: logfile
    real(kind=DBL_KIND),intent(IN) :: pmass, pr(MX:MZ), pv(MX:MZ), ps(MX:MZ)
    integer,intent(IN) :: pid
    real(kind=DBL_KIND) :: dv, vol, mass, kt, kg, ki, buf(6), eGrav, eTh, eMag, skt
    real(kind=DBL_KIND),dimension(MX:MZ) :: rg, vg, jt, js, jg, bave, ksi, ksg, sjt, sksi, sksg, sbsi, sbsg
    real(kind=DBL_KIND),dimension(MX:MZ,MX:MZ) :: it, sit
    integer :: level, n, gid, i, j, k
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, gx, gy, gz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(ARRAYSIZE3(GridMask)) :: xx, yy, zz, dvx, dvy, dvz
    real(kind=DBL_KIND),dimension(lbound(GridMask,1):ubound(GridMask,1)) :: dx, xn
    real(kind=DBL_KIND),dimension(lbound(GridMask,2):ubound(GridMask,2)) :: dy, yn
    real(kind=DBL_KIND),dimension(lbound(GridMask,3):ubound(GridMask,3)) :: dz, zn
    logical,dimension(ARRAYSIZE3(GridMask)) :: mask
    integer,parameter ::  FH=11
    character(len=CHARLEN) :: dir, fn
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    real(kind=DBL_KIND) :: xlenCompBox, ylenCompBox, zlenCompBox, xshift, yshift, zshift
    real(kind=DBL_KIND) :: rs, rc, sint, cost, sinp, cosp, vr, vt,vp, eps, br, bt, bp
    logical :: bool_touch, bbuf
    integer :: ibool_touch, istat
    real(kind=DBL_KIND) :: buf33(3,3), buf6(6), buf3(3), buf1

    real(kind=DBL_KIND) :: flag, p

    ! size of computation box
    call ob_computationBoxOfCoordPhys( coordPhys )
    xlenCompBox = coordPhys(MZ+1+MX) - coordPhys(MX)
    ylenCompBox = coordPhys(MZ+1+MY) - coordPhys(MY)
    zlenCompBox = coordPhys(MZ+1+MZ) - coordPhys(MZ)
    xshift = coordPhys(MX) + pr(MX)
    yshift = coordPhys(MY) + pr(MY)
    zshift = coordPhys(MZ) + pr(MZ)

    myrank = get_myrank()
    vol = 0.d0
    mass = 0.d0
    eMag = 0.d0
    kt = 0.d0
    bave = 0.d0
    rg = 0.d0
    vg = 0.d0
    jt = 0.d0
    it = 0.d0
    ksi = 0.d0
    ksg = 0.d0
    eGrav = 0.d0
    eTh = 0.d0
    bool_touch = .false.
    do level = Lmin, LevelMax
       dv = get_dv( level )
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( MRHO, gid )
          vx => get_Ucomp( MVX, gid )
          vy => get_Ucomp( MVY, gid )
          vz => get_Ucomp( MVZ, gid )
#ifdef MHD
          bx => get_Ucomp( MBX, gid )
          by => get_Ucomp( MBY, gid )
          bz => get_Ucomp( MBZ, gid )
#endif !MHD
          gx => get_Ucomp( MGX, gid )
          gy => get_Ucomp( MGY, gid )
          gz => get_Ucomp( MGZ, gid )
          x => get_Xp( gid )
          y => get_Yp( gid )
          z => get_Zp( gid )
          ! shift coordinates due to periodic boundary condition.
          xn = modulo(x - xshift, xlenCompBox) + xshift
          yn = modulo(y - yshift, ylenCompBox) + yshift
          zn = modulo(z - zshift, zlenCompBox) + zshift
          ! spread from 1D to 3D arrary
          xx = spread(spread(xn, 2, size(y)), 3, size(z))
          yy = spread(spread(yn, 1, size(x)), 3, size(z))
          zz = spread(spread(zn, 1, size(x)), 2, size(y))

          mask = GridMask .and. ( rho >= rhoCritical ) .and. ((xx-pr(MX))**2 + (yy-pr(MY))**2 + (zz-pr(MZ))**2 < radius**2)
          if ( count(mask) == 0 ) cycle

          bool_touch = bool_touch .or. bool_maskoverflow()

          vol  = vol  + COUNT( mask ) * dv
          mass = mass + SUM(rho, mask=mask) * dv
          kt = kt + SUM(rho*(vx**2+vy**2+vz**2), mask=mask) * dv
#ifdef MHD
          eMag = eMag + SUM(bx**2+by**2+bz**2, mask=mask) * dv
          bave(MX) = bave(MX) + SUM(bx, mask=mask)*dv
          bave(MY) = bave(MY) + SUM(by, mask=mask)*dv
          bave(MZ) = bave(MZ) + SUM(bz, mask=mask)*dv
#endif !MHD

          rg(MX) = rg(MX) + SUM(rho*xx, mask=mask)*dv
          rg(MY) = rg(MY) + SUM(rho*yy, mask=mask)*dv
          rg(MZ) = rg(MZ) + SUM(rho*zz, mask=mask)*dv

          vg(MX) = vg(MX) + SUM(rho*vx, mask=mask)*dv
          vg(MY) = vg(MY) + SUM(rho*vy, mask=mask)*dv
          vg(MZ) = vg(MZ) + SUM(rho*vz, mask=mask)*dv

          jt(MX) = jt(MX) + SUM(rho*(yy*vz-zz*vy), mask=mask)*dv
          jt(MY) = jt(MY) + SUM(rho*(zz*vx-xx*vz), mask=mask)*dv
          jt(MZ) = jt(MZ) + SUM(rho*(xx*vy-yy*vx), mask=mask)*dv

          eGrav = eGrav + SUM(gx**2+gy**2+gz**2, mask=mask) * dv

          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   if ( mask(i,j,k) ) then
                      GETP(p, rho(i,j,k))
                      eTh = eTh + p / (Gamma-1.d0) * dv
                   endif
                enddo
             enddo
          enddo
       end do
    end do
    buf = (/vol, mass, eMag, kt, eGrav, eTh/)
    call mpi_allreduce(buf, buf6, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = buf6
    vol = buf(1); mass = buf(2); eMag = buf(3); kt = buf(4); eGrav = buf(5); eTh = buf(6)
    call mpi_allreduce(bave, buf3, size(bave), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    bave = buf3
    call mpi_allreduce(rg, buf3, size(rg), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    rg = buf3
    call mpi_allreduce(vg, buf3, size(vg), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    vg = buf3
    call mpi_allreduce(jt, buf3, size(jt), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    jt = buf3
    call mpi_allreduce(bool_touch, bbuf, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    bool_touch = bbuf
    eMag = eMag/pi8             ! B^2 / 8 pi
    bave = bave/vol             ! bar(B)
    rg = rg/mass                ! wieght of position (baricenter)
    vg = vg/mass                ! velocity of baricenter
    kt = kt/mass                ! < v^2 > total kinetic energy
    kg = sum(vg**2)             ! vg^2    kinetic energy of baricenter
    ki = kt - kg                ! < (v - vg)^2 > internal turbulent energy
    jt = jt/mass                ! total specific angular momentum
    jg(MX) = rg(MY)*vg(MZ) - rg(MZ)*vg(MY) ! specific angular momentum of baricenter
    jg(MY) = rg(MZ)*vg(MX) - rg(MX)*vg(MZ)
    jg(MZ) = rg(MX)*vg(MY) - rg(MY)*vg(MX)
    js = jt - jg                ! spin specific angular momentum
    eGrav = eGrav * 0.5d0
    if (bool_touch) then
       ibool_touch = 1
    else
       ibool_touch = 0
    endif

    do level = Lmin, LevelMax
       dv = get_dv( level )
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( MRHO, gid )
          vx => get_Ucomp( MVX, gid )
          vy => get_Ucomp( MVY, gid )
          vz => get_Ucomp( MVZ, gid )
          x => get_Xp( gid )
          y => get_Yp( gid )
          z => get_Zp( gid )
          ! shift coordinates due to periodic boundary condition.
          xn = modulo(x - xshift, xlenCompBox) + xshift
          yn = modulo(y - yshift, ylenCompBox) + yshift
          zn = modulo(z - zshift, zlenCompBox) + zshift
          ! deviation from baricenter
          dx = xn - rg(MX)
          dy = yn - rg(MY)
          dz = zn - rg(MZ)
          ! spread from 1D to 3D arrary
          xx = spread(spread(dx, 2, size(y)), 3, size(z))
          yy = spread(spread(dy, 1, size(x)), 3, size(z))
          zz = spread(spread(dz, 1, size(x)), 2, size(y))

          mask = GridMask .and. ( rho >= rhoCritical ) .and. (xx**2 + yy**2 + zz**2 < radius**2)
          if ( count(mask) == 0 ) cycle

          it(MX,MX) = it(MX,MX) + SUM(rho*xx**2, mask=mask)*dv
          it(MY,MY) = it(MY,MY) + SUM(rho*yy**2, mask=mask)*dv
          it(MZ,MZ) = it(MZ,MZ) + SUM(rho*zz**2, mask=mask)*dv
          it(MX,MY) = it(MX,MY) + SUM(rho*xx*yy, mask=mask)*dv
          it(MX,MZ) = it(MX,MZ) + SUM(rho*xx*zz, mask=mask)*dv
          it(MY,MZ) = it(MY,MZ) + SUM(rho*yy*zz, mask=mask)*dv

          ! cartesian to spherical coocrdinates
#define MR_ MX
#define MT_ MY
#define MP_ MZ
          eps = minval(CellWidth(:,level))*1.D-6
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   if ( mask(i,j,k) ) then
                      rs = sqrt(dx(i)**2+dy(j)**2+dz(k)**2) ! spherical radius
                      rc = sqrt(dx(i)**2+dy(j)**2)          ! cylindrical radius
                      sint=rc/(rs+eps)
                      cost=dz(k) /(rs+eps)
                      sinp=dy(j) /(rc+eps)
                      cosp=dx(i) /(rc+eps)
                      vr= (vx(i,j,k)-vg(MX))*sint*cosp+(vy(i,j,k)-vg(MY))*sint*sinp+(vz(i,j,k)-vg(MZ))*cost
                      vt= (vx(i,j,k)-vg(MX))*cost*cosp+(vy(i,j,k)-vg(MY))*cost*sinp-(vz(i,j,k)-vg(MZ))*sint
                      vp=-(vx(i,j,k)-vg(MX))*sinp+(vy(i,j,k)-vg(MY))*cosp
                      ksi(MR_) =  ksi(MR_) + vr**2 * rho(i,j,k)*dv ! <vs^2>
                      ksi(MT_) =  ksi(MT_) + vt**2 * rho(i,j,k)*dv
                      ksi(MP_) =  ksi(MP_) + vp**2 * rho(i,j,k)*dv
                      ksg(MR_) =  ksg(MR_) + vr * rho(i,j,k)*dv ! <vs>^2
                      ksg(MT_) =  ksg(MT_) + vt * rho(i,j,k)*dv
                      ksg(MP_) =  ksg(MP_) + vp * rho(i,j,k)*dv
                   endif
                enddo
             enddo
          enddo
       end do
    end do

    call mpi_allreduce(it, buf33,  size(it), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    it = buf33
    it(MY,MX) = it(MX,MY)
    it(MZ,MX) = it(MX,MZ)
    it(MZ,MY) = it(MY,MZ)
    it = it/mass                ! moment of inertia

    call mpi_allreduce(ksi, buf3, size(ksi), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ksi = buf3
    call mpi_allreduce(ksg, buf3, size(ksg), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ksg = buf3
    ksi = ksi / mass            ! <vs^2>
    ksg = (ksg / mass)**2       ! <vs>^2

    ! -----------------------------
    ! for around sink particles
    ! -----------------------------
    sjt = 0.d0
    skt = 0.d0
    sit = 0.d0
    sksi = 0.d0
    sksg = 0.d0
    sbsi = 0.d0
    sbsg = 0.d0
    if (pid /= -1) then
       ! integral with respective to sink particle
       do level = Lmin, LevelMax
          dv = get_dv( level )
          do n = Gidmin, GidListMax( level )
             gid = GidList(n, level)
             if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
             rho => get_Ucomp( MRHO, gid )
             vx => get_Ucomp( MVX, gid )
             vy => get_Ucomp( MVY, gid )
             vz => get_Ucomp( MVZ, gid )
             x => get_Xp( gid )
             y => get_Yp( gid )
             z => get_Zp( gid )
#ifdef MHD
             bx => get_Ucomp( MBX, gid )
             by => get_Ucomp( MBY, gid )
             bz => get_Ucomp( MBZ, gid )
#endif !MHD
             ! shift coordinates due to periodic boundary condition.
             xn = modulo(x - xshift, xlenCompBox) + xshift
             yn = modulo(y - yshift, ylenCompBox) + yshift
             zn = modulo(z - zshift, zlenCompBox) + zshift
             ! spread from 1D to 3D arrary
             xx = spread(spread(xn, 2, size(y)), 3, size(z))
             yy = spread(spread(yn, 1, size(x)), 3, size(z))
             zz = spread(spread(zn, 1, size(x)), 2, size(y))

             ! from sink particle
             xx = xx - pr(MX)
             yy = yy - pr(MY)
             zz = zz - pr(MZ)

             mask = GridMask .and. ( rho >= rhoCritical ) .and. (xx**2 + yy**2 + zz**2 < radius**2)
             if ( count(mask) == 0 ) cycle

             ! relative velocity
             dvx = vx - pv(MX)
             dvy = vy - pv(MY)
             dvz = vz - pv(MZ)
             ! angular momentum
             sjt(MX) = sjt(MX) + SUM(rho*(yy*dvz-zz*dvy), mask=mask)*dv
             sjt(MY) = sjt(MY) + SUM(rho*(zz*dvx-xx*dvz), mask=mask)*dv
             sjt(MZ) = sjt(MZ) + SUM(rho*(xx*dvy-yy*dvx), mask=mask)*dv
             ! <(v-pv)^2>
             skt = skt + SUM(rho*(dvx**2+dvy**2+dvz**2), mask=mask) * dv
             ! moment of intertia
             sit(MX,MX) = sit(MX,MX) + SUM(rho*xx**2, mask=mask)*dv
             sit(MY,MY) = sit(MY,MY) + SUM(rho*yy**2, mask=mask)*dv
             sit(MZ,MZ) = sit(MZ,MZ) + SUM(rho*zz**2, mask=mask)*dv
             sit(MX,MY) = sit(MX,MY) + SUM(rho*xx*yy, mask=mask)*dv
             sit(MX,MZ) = sit(MX,MZ) + SUM(rho*xx*zz, mask=mask)*dv
             sit(MY,MZ) = sit(MY,MZ) + SUM(rho*yy*zz, mask=mask)*dv
             ! cartesian to spherical coocrdinates
#define MR_ MX
#define MT_ MY
#define MP_ MZ
             eps = minval(CellWidth(:,level))*1.D-6
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i = Imin, Imax
                      if ( mask(i,j,k) ) then
                         rs = sqrt(xx(i,j,k)**2+yy(i,j,k)**2+zz(i,j,k)**2) ! spherical radius
                         rc = sqrt(xx(i,j,k)**2+yy(i,j,k)**2)          ! cylindrical radius
                         sint=rc/(rs+eps)
                         cost=zz(i,j,k) /(rs+eps)
                         sinp=yy(i,j,k) /(rc+eps)
                         cosp=xx(i,j,k) /(rc+eps)
                         vr= dvx(i,j,k)*sint*cosp+dvy(i,j,k)*sint*sinp+dvz(i,j,k)*cost
                         vt= dvx(i,j,k)*cost*cosp+dvy(i,j,k)*cost*sinp-dvz(i,j,k)*sint
                         vp=-dvx(i,j,k)*sinp+dvy(i,j,k)*cosp
                         sksi(MR_) =  sksi(MR_) + vr**2 * rho(i,j,k)*dv ! <vs^2>
                         sksi(MT_) =  sksi(MT_) + vt**2 * rho(i,j,k)*dv
                         sksi(MP_) =  sksi(MP_) + vp**2 * rho(i,j,k)*dv
                         sksg(MR_) =  sksg(MR_) + vr * rho(i,j,k)*dv ! <vs>
                         sksg(MT_) =  sksg(MT_) + vt * rho(i,j,k)*dv
                         sksg(MP_) =  sksg(MP_) + vp * rho(i,j,k)*dv
#ifdef MHD
                         br= bx(i,j,k)*sint*cosp+by(i,j,k)*sint*sinp+bz(i,j,k)*cost
                         bt= bx(i,j,k)*cost*cosp+by(i,j,k)*cost*sinp-bz(i,j,k)*sint
                         bp=-bx(i,j,k)*sinp+by(i,j,k)*cosp
                         sbsi(MR_) =  sbsi(MR_) + br**2 * dv ! <bs^2>
                         sbsi(MT_) =  sbsi(MT_) + bt**2 * dv
                         sbsi(MP_) =  sbsi(MP_) + bp**2 * dv
                         sbsg(MR_) =  sbsg(MR_) + br * dv ! <bs>
                         sbsg(MT_) =  sbsg(MT_) + bt * dv
                         sbsg(MP_) =  sbsg(MP_) + bp * dv
#endif !MHD
                      endif
                   enddo
                enddo
             enddo
          end do
       end do
       ! all reduce
       call mpi_allreduce(sjt, buf3, size(sjt), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sjt = buf3
       call mpi_allreduce(skt, buf1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       skt = buf1
       call mpi_allreduce(sit, buf33, size(sit), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sit = buf33
       call mpi_allreduce(sksi, buf3, size(sksi), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sksi = buf3
       call mpi_allreduce(sksg, buf3, size(sksg), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sksg = buf3
       call mpi_allreduce(sbsi, buf3, size(sbsi), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sbsi = buf3
       call mpi_allreduce(sbsg, buf3, size(sbsg), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       sbsg = buf3

       sjt = sjt/mass              ! angular momentum
       skt = skt/mass              ! <(v-pv)^2>
       sit(MY,MX) = sit(MX,MY)
       sit(MZ,MX) = sit(MX,MZ)
       sit(MZ,MY) = sit(MY,MZ)
       sit = sit/mass                ! moment of inertia
       sksi = sksi / mass            ! <vr^2>, <vp^2>, <vt^2>
       sksg = (sksg / mass)**2       ! <vr>^2, <vp>^2, <vt>^2
       sbsi = sbsi / vol             ! <Br^2>, <Bp^2>, <Bt^2>
       sbsg = (sbsg / vol)**2        ! <Br>^2, <Bp>^2, <Bt>^2
    endif

    ! -----------------------------------
    ! output data
    ! -----------------------------------
    if ( get_myrank() == PRIMARY_RANK ) then
       call read_env('DIR', dir)
       fn = concat(dir, logfile)
       open(FH, file=fn, position='APPEND')
       write(FH,'(I12, 1PE12.5, I3, 77(1PE12.5), I3)') Step(Lmin), Time(Lmin), pid, pmass, pr(:), pv(:), ps(:), RhoMax, vol, mass, eMag, &
            bave(:), rg(:), vg(:), kt, kg, ki, jt(:), jg(:), js(:), it(:,:), ksi(:), ksg(:), eGrav, eTh, &
            sjt(:), skt, sit(:,:), sksi(:), sksg(:), sbsi(:), sbsg(:), &
            ibool_touch
       call flush(FH)
       close(FH)
    endif
  contains
    ! -------------------------------------------------------------------------
    ! check mask for overflowing region
    ! -------------------------------------------------------------------------
    function bool_maskoverflow() result(bool)
      logical :: bool
      real(kind=DBL_KIND) :: xlmax, ylmax, zlmax, xlmin, ylmin, zlmin
      bool = .false.
      xlmin = minval(xx, mask=mask)
      xlmax = maxval(xx, mask=mask)
      ylmin = minval(yy, mask=mask)
      ylmax = maxval(yy, mask=mask)
      zlmin = minval(zz, mask=mask)
      zlmax = maxval(zz, mask=mask)
      if ( xlmax - pr(MX) > xlenCompBox/2 ) bool = .true.
      if ( pr(MX) - xlmax > xlenCompBox/2 ) bool = .true.
      if ( ylmax - pr(MY) > ylenCompBox/2 ) bool = .true.
      if ( pr(MY) - ylmax > ylenCompBox/2 ) bool = .true.
      if ( zlmax - pr(MZ) > zlenCompBox/2 ) bool = .true.
      if ( pr(MZ) - zlmax > zlenCompBox/2 ) bool = .true.
    end function bool_maskoverflow
  end subroutine an_sinkAround_eachParticle

end module analysis
