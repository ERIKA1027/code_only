#include "config.h"
! ----------------------------------------------------------------------------
! Module for potential field source surface.
! ----------------------------------------------------------------------------
module pfss
  implicit none
  private
  integer,parameter :: ThetaStart = 5 ! colatitude
  integer,parameter :: ThetaEnd = 175
  integer,parameter :: ThetaSkip = 5
  integer,parameter :: PhiStart = 0 ! longitude
  integer,parameter :: PhiEnd = 360
  integer,parameter :: PhiSkip = 5
  integer,parameter :: NTheta = (ThetaEnd-ThetaStart)/ThetaSkip+1
  integer,parameter :: NPhi = (PhiEnd-PhiStart)/PhiSkip+1
  real(kind=DBL_KIND),dimension(NTheta) :: ThetaSS ! theta at source surface
  real(kind=DBL_KIND),dimension(NPhi) :: PhiSS     ! phi at source surface
  ! Superradial expansion factor, Fs
  ! Radial magnetic field on the source surface, BrSS
  real(kind=DBL_KIND),dimension(Nphi, NTheta) :: Fs, BrSS
  real(kind=DBL_KIND),dimension(Nphi, NTheta) :: Vrsw, Vpsw, Brsw, Bpsw, Rhosw, Psw
  !
  logical,save :: Bool_synop_chart_provide = .FALSE.
  public :: pfss_synop_chart_provide, pfss_boundary_gid, pfss_init_cond_gid
contains
  ! ----------------------------------------------------------------------------
  ! Innner boundary condition according to solar wind model
  ! ----------------------------------------------------------------------------
  subroutine pfss_boundary_gid(gid, t)
    use modelParameter, only : MP_rInnerBoundary, MP_drInnerBoundary, MP_CarringtonRotDay
    use unit, only : Unit_day, Unit_n, Unit_e, cgs_mu, cgs_kb
    use parameter
    use grid
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),intent(IN) :: t
    integer :: i, j, k, ipha, iphb, jtha, jthb
    real(kind=DBL_KIND) :: vr, vp, vt, br, bp, bt, r, theta, phi, dtheta, dphi, dpha, dphb, dtha, dthb, s2d, sint, cost, sinp, cosp, omega, nh0, rho0, p0, temp0
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, db

    if (.not. Bool_synop_chart_provide) then
       print *, '*** error Bool_synop_chart_provide is', Bool_synop_chart_provide
       stop
    end if
    
    if (.not. pfss_touches_innerboundary(gid)) return
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    rho => get_Ucomp(MRHO, gid)
    vx => get_Ucomp(MVX, gid)
    vy => get_Ucomp(MVY, gid)
    vz => get_Ucomp(MVZ, gid)
    bx => get_Ucomp(MBX, gid)
    by => get_Ucomp(MBY, gid)
    bz => get_Ucomp(MBZ, gid)
    p  => get_Ucomp(MP,  gid)
    db => get_Ucomp(MDB, gid)
    dtheta = ThetaSkip * Pi/180.d0
    dphi = PhiSkip * Pi/180.d0
    omega = Pi2 / (MP_CarringtonRotDay/Unit_day)
    nh0 = Unit_n                ! default hydrogen number density
    temp0 = 1.e6                ! default temperature in K
    rho0 = nh0 / Unit_n
    p0 = nh0/cgs_mu * cgs_kb * temp0 / Unit_e
    do k = Kmingh, Kmaxgh
       do j = Jmingh, Jmaxgh
          do i = Imingh, Imaxgh
             r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
             if (r <= MP_rInnerBoundary .and. r >= MP_rInnerBoundary - MP_drInnerBoundary) then
                theta = atan2(sqrt(x(i)**2 + y(j)**2), z(k))            ! [0, pi]
                phi = atan2(y(j), x(i))

                sint = sin(theta)
                cost = cos(theta)
                sinp = sin(phi)
                cosp = cos(phi)

                ! shift according to rotation
                phi = phi - t * omega
                phi = modulo(phi, Pi2) ! [0, 2pi]

                ipha = int((phi - PhiSS(1))/dphi) + 1
                ipha = min(max(ipha,1), NPhi-1)
                iphb = ipha + 1
                jtha = int((theta - ThetaSS(1))/dtheta) + 1
                jtha = min(max(jtha, 1), NTheta-1)
                jthb = jtha + 1

                dpha = phi - PhiSS(ipha)
                dphb = PhiSS(iphb) - phi
                dtha = theta - ThetaSS(jtha)
                dthb = ThetaSS(jthb) - theta
                s2d = dphi * dtheta

                vr = ( &
                     Vrsw(iphb,jthb)*dpha*dtha + Vrsw(ipha,jthb)*dphb*dtha + &
                     Vrsw(iphb,jtha)*dpha*dthb + Vrsw(ipha,jtha)*dphb*dthb ) / s2d
                vp = ( &
                     Vpsw(iphb,jthb)*dpha*dtha + Vpsw(ipha,jthb)*dphb*dtha + &
                     Vpsw(iphb,jtha)*dpha*dthb + Vpsw(ipha,jtha)*dphb*dthb ) / s2d
                vt = 0.d0
                vx(i,j,k) = vr*sint*cosp + vt*cost*cosp - vp*sinp
                vy(i,j,k) = vr*sint*sinp + vt*cost*sinp + vp*cosp
                vz(i,j,k) = vr*cost - vt*sint


                br = ( &
                     Brsw(iphb,jthb)*dpha*dtha + Brsw(ipha,jthb)*dphb*dtha + &
                     Brsw(iphb,jtha)*dpha*dthb + Brsw(ipha,jtha)*dphb*dthb ) / s2d
                bp = ( &
                     Bpsw(iphb,jthb)*dpha*dtha + Bpsw(ipha,jthb)*dphb*dtha + &
                     Bpsw(iphb,jtha)*dpha*dthb + Bpsw(ipha,jtha)*dphb*dthb ) / s2d
                bt = 0.d0
                bx(i,j,k) = br*sint*cosp + bt*cost*cosp - bp*sinp
                by(i,j,k) = br*sint*sinp + bt*cost*sinp + bp*cosp
                bz(i,j,k) = br*cost - bt*sint

                rho(i,j,k) = ( &
                     Rhosw(iphb,jthb)*dpha*dtha + Rhosw(ipha,jthb)*dphb*dtha + &
                     Rhosw(iphb,jtha)*dpha*dthb + Rhosw(ipha,jtha)*dphb*dthb ) / s2d
                p(i,j,k) = ( &
                     Psw(iphb,jthb)*dpha*dtha + Psw(ipha,jthb)*dphb*dtha + &
                     Psw(iphb,jtha)*dpha*dthb + Psw(ipha,jtha)*dphb*dthb ) / s2d
                db(i,j,k) = 0.d0
             else if ( r < MP_rInnerBoundary - MP_drInnerBoundary ) then
                rho(i,j,k) = rho0
                vx(i,j,k) = 0.d0
                vy(i,j,k) = 0.d0
                vz(i,j,k) = 0.d0
                bx(i,j,k) = 0.d0
                by(i,j,k) = 0.d0
                bz(i,j,k) = 0.d0
                p(i,j,k) = p0
                db(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
    
  end subroutine pfss_boundary_gid
  ! ----------------------------------------------------------------------------
  ! provide boundary condition by PFSS
  ! ----------------------------------------------------------------------------
  function pfss_touches_innerboundary(gid) result(bool)
    use grid
    use modelParameter, only : MP_rInnerBoundary
    integer,intent(IN) :: gid
    logical :: bool
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: xedge, yedge, zedge, dist
    integer :: level
    bool = .FALSE.              ! default value
    if ( gid == Undefi ) return
    level = get_level(gid)
    if ( level /= LevelMax ) return ! solar wind is only for finest grid.
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    h = CellWidth(:,level)
    ! inner edges of the block
    xedge = min(abs(x(Imin)-h(MX)*0.5d0), abs(x(Imax)+h(MX)*0.5d0))
    yedge = min(abs(y(Jmin)-h(MY)*0.5d0), abs(y(Jmax)+h(MY)*0.5d0))
    zedge = min(abs(z(Kmin)-h(MZ)*0.5d0), abs(z(Kmax)+h(MZ)*0.5d0))
    dist = sqrt(xedge**2 + yedge**2 + zedge**2)
    if (dist > MP_rInnerBoundary) return ! block is outside the inner boundary
    bool = .TRUE.
  end function pfss_touches_innerboundary
  ! ----------------------------------------------------------------------------
  ! load synoptic chart for all vars
  ! ----------------------------------------------------------------------------
  subroutine pfss_synop_chart_provide(file)
    use mpilib
    character(len=*),intent(IN) :: file
    if (get_myrank() == PRIMARY_RANK) then
       call pfss_synop_chart_read(file)
    end if
    call pfss_synop_chart_bcast
    call pfss_synop_chart_make_sw
    Bool_synop_chart_provide = .TRUE.
  end subroutine pfss_synop_chart_provide
  ! ----------------------------------------------------------------------------
  ! provide solar wind variables
  ! ----------------------------------------------------------------------------
  subroutine pfss_synop_chart_make_sw
    use parameter
    use modelParameter
    use unit
    integer :: i, j
    real(kind=DBL_KIND) :: vrsw_kms, rsun50, ncc, tsw, omega
    rsun50 = 50.d0 * MP_SolarRadius
    omega = Pi2 / (MP_CarringtonRotDay/Unit_day)
    do j = 1, NTheta
       do i = 1, Nphi
          ! Shiota et al. (2014), eq. (9), or Arge & Pizzo (2000), eq. (4)
          vrsw_kms = 267.5d0 + 410.d0/Fs(i,j)**0.4
          Vrsw(i,j) = vrsw_kms / Unit_kms
!!$          Vpsw(i,j) = 2 * MP_rInnerBoundary * omega
          Vpsw(i,j) = 0.d0

          ! Kataoka et al. (2009), eq. (3)
!!$          Brsw(i,j) = 10.d0* sign(1.d0, BrSS(i,j))*sqrt(abs(BrSS(i,j))/Unit_utesla) / (MP_rInnerBoundary/MP_SourceSurface)**2
          ! Simple formula  Br ~ r^(-2)
          Brsw(i,j) = BrSS(i,j)/Unit_utesla / (MP_rInnerBoundary/MP_SourceSurface)**2

          ! Kataoka et al. (2009), eq. (6)
          Bpsw(i,j) = -Brsw(i,j) * MP_rInnerBoundary * omega * sin(ThetaSS(j))/Vrsw(i,j)
!!$          Bpsw(i,j) = 0.d0

          ! Shiota et al. (2014), eq. (11), or Hayashi et al. (2003), eq. (10)
          ncc = (rsun50/MP_rInnerBoundary)**2 * (62.98d0 + 866.4d0 * (vrsw_kms/100 - 1.549)**(-3.402))
          Rhosw(i,j) = ncc / Unit_n

          ! Shiota et al. (2014), eq. (12), or Hayashi et al. (2003), eq. (11)
          tsw = (rsun50/MP_rInnerBoundary)**(2*(MP_Gamma-1)) * (-0.455d0 + 0.1943 * vrsw_kms/100)*1.D6 ! in K
          Psw(i,j) = ncc/cgs_mu * cgs_kb * tsw / Unit_e
       end do
    end do
  end subroutine pfss_synop_chart_make_sw
  ! ----------------------------------------------------------------------------
  ! read synop chart for all ranks
  ! ----------------------------------------------------------------------------
  subroutine pfss_synop_chart_bcast
    use mpilib
    call mpi_bcast(ThetaSS, size(ThetaSS), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(PhiSS, size(PhiSS), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(Fs, size(Fs), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(BrSS, size(BrSS), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
  end subroutine pfss_synop_chart_bcast
  ! ----------------------------------------------------------------------------
  ! read data from file
  ! ----------------------------------------------------------------------------
  subroutine pfss_synop_chart_read(file)
    character(len=*),intent(IN) :: file
    integer,parameter :: UNIT = 11
    integer :: nphi_r, ntheta_r
    print *, 'pfss file reading... ', trim(file)
    open(UNIT, file=file, form='unformatted')
    read(UNIT) nphi_r, ntheta_r
    if (nphi_r /= NPhi) print *, '*** nphi_r /= NPhi'
    if (ntheta_r /= NTheta) print *, '*** ntheta_r /= NTheta'
    read(UNIT) PhiSS
    read(UNIT) ThetaSS
    read(UNIT) Fs
    read(UNIT) Brss
    call flush(UNIT)
    close(UNIT)
  end subroutine pfss_synop_chart_read
  ! ----------------------------------------------------------------------------
  ! Initial condition
  ! ----------------------------------------------------------------------------
  subroutine pfss_init_cond_gid(gid)
    use modelParameter, only : MP_rInnerBoundary, MP_drInnerBoundary, MP_CarringtonRotDay, MP_Gamma
    use unit, only : Unit_day, Unit_n, Unit_e, cgs_mu, cgs_kb, Unit_kms, Unit_au
    use parameter
    use grid
    integer,intent(IN) :: gid
    real(kind=DBL_KIND) :: t
    real(kind=DBL_KIND),parameter :: VR_INIT_ASYMPTOTE = 1000.d0 ! in km/s
    real(kind=DBL_KIND),parameter :: R_VR_INIT_ASYMPTOTE = 200.d0 ! in au
    integer :: i, j, k, ipha, iphb, jtha, jthb
    real(kind=DBL_KIND) :: vr, vp, vt, br, bp, bt, r, theta, phi, dtheta, dphi, dpha, dphb, dtha, dthb, s2d, sint, cost, sinp, cosp, omega
    real(kind=DBL_KIND) :: rho0, p0, vr0, br0, vra, ra
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, db

    if (.not. Bool_synop_chart_provide) then
       print *, '*** error Bool_synop_chart_provide is', Bool_synop_chart_provide
       stop
    end if


    t = 0.d0                    ! initial condition
    vra = VR_INIT_ASYMPTOTE / Unit_kms
    ra = R_VR_INIT_ASYMPTOTE / Unit_au
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    rho => get_Ucomp(MRHO, gid)
    vx => get_Ucomp(MVX, gid)
    vy => get_Ucomp(MVY, gid)
    vz => get_Ucomp(MVZ, gid)
    bx => get_Ucomp(MBX, gid)
    by => get_Ucomp(MBY, gid)
    bz => get_Ucomp(MBZ, gid)
    p  => get_Ucomp(MP,  gid)
    db => get_Ucomp(MDB, gid)
    dtheta = ThetaSkip * Pi/180.d0
    dphi = PhiSkip * Pi/180.d0
    omega = Pi2 / (MP_CarringtonRotDay/Unit_day)
    do k = Kmingh, Kmaxgh
       do j = Jmingh, Jmaxgh
          do i = Imingh, Imaxgh
             r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

             theta = atan2(sqrt(x(i)**2 + y(j)**2), z(k))            ! [0, pi]
             phi = atan2(y(j), x(i))

             sint = sin(theta)
             cost = cos(theta)
             sinp = sin(phi)
             cosp = cos(phi)

             ! shift according to rotation
             phi = phi - t * omega
             phi = modulo(phi, Pi2) ! [0, 2pi]

             ipha = int((phi - PhiSS(1))/dphi) + 1
             ipha = min(max(ipha,1), NPhi-1)
             iphb = ipha + 1
             jtha = int((theta - ThetaSS(1))/dtheta) + 1
             jtha = min(max(jtha, 1), NTheta-1)
             jthb = jtha + 1

             dpha = phi - PhiSS(ipha)
             dphb = PhiSS(iphb) - phi
             dtha = theta - ThetaSS(jtha)
             dthb = ThetaSS(jthb) - theta
             s2d = dphi * dtheta

             vr0 = ( &
                  Vrsw(iphb,jthb)*dpha*dtha + Vrsw(ipha,jthb)*dphb*dtha + &
                  Vrsw(iphb,jtha)*dpha*dthb + Vrsw(ipha,jtha)*dphb*dthb ) / s2d
             vr = ( vr0 * (ra - r) + vra * (r - MP_rInnerBoundary))/(ra - MP_rInnerBoundary)
             vp = 0.d0
             vt = 0.d0

             vx(i,j,k) = vr*sint*cosp + vt*cost*cosp - vp*sinp
             vy(i,j,k) = vr*sint*sinp + vt*cost*sinp + vp*cosp
             vz(i,j,k) = vr*cost - vt*sint


             br0 = ( &
                  Brsw(iphb,jthb)*dpha*dtha + Brsw(ipha,jthb)*dphb*dtha + &
                  Brsw(iphb,jtha)*dpha*dthb + Brsw(ipha,jtha)*dphb*dthb ) / s2d
             br = br0 * (MP_rInnerBoundary/r)**2
             bp = 0.d0
             bt = 0.d0

             bx(i,j,k) = br*sint*cosp + bt*cost*cosp - bp*sinp
             by(i,j,k) = br*sint*sinp + bt*cost*sinp + bp*cosp
             bz(i,j,k) = br*cost - bt*sint

             rho0 = ( &
                  Rhosw(iphb,jthb)*dpha*dtha + Rhosw(ipha,jthb)*dphb*dtha + &
                  Rhosw(iphb,jtha)*dpha*dthb + Rhosw(ipha,jtha)*dphb*dthb ) / s2d
             rho(i,j,k) = rho0 * MP_rInnerBoundary**2 * vr0 / (r**2 * vr)
             p0 = ( &
                  Psw(iphb,jthb)*dpha*dtha + Psw(ipha,jthb)*dphb*dtha + &
                  Psw(iphb,jtha)*dpha*dthb + Psw(ipha,jtha)*dphb*dthb ) / s2d
             p(i,j,k) = p0 * (rho(i,j,k)/rho0) ! **MP_Gamma
             db(i,j,k) = 0.d0
          end do
       end do
    end do

  end subroutine pfss_init_cond_gid
end module pfss
