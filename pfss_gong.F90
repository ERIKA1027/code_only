
#include "config.h"
! check_update_keystages
!    if time in [stage_a, stage_b)
! 
! if (bool_update) then
!    stage A (yyyymmdd) を用意     
!         if (yyyymmdd+1 exist?) then
!              copy from stage B
!         else
!              read new file (初回だけ実行される)
!    stage B (yyyymmdd+1) を用意
!         read new file
!
! linear interpolation in time
!     interplation is performed in the unit of simulation time
! linear interpolation in space
!     rotation of sun
! 
! ----------------------------------------------------------------------------
! Module for potential field source surface.
! ----------------------------------------------------------------------------
module pfss
  use string, only : CHARLEN
  use grid, only : Undefi
  implicit none
  private
  integer,parameter :: ThetaStart = 0 ! colatitude (degree)
  integer,parameter :: ThetaEnd = 180
  integer,parameter :: ThetaSkip = 1
  integer,parameter :: PhiStart = 0 ! longitude (degree)
  integer,parameter :: PhiEnd = 360
  integer,parameter :: PhiSkip = 1
  integer,parameter :: NTheta = (ThetaEnd-ThetaStart)/ThetaSkip+1
  integer,parameter :: NPhi = (PhiEnd-PhiStart)/PhiSkip+1
  real(kind=DBL_KIND),dimension(NTheta) :: ThetaSS ! theta at source surface (radian)
  real(kind=DBL_KIND),dimension(NPhi) :: PhiSS     ! phi at source surface (radian)
  !
  logical,save :: Bool_Initialized = .FALSE.
  integer,parameter :: CHARLEN_YYYYMMDD = 8
  character(len=CHARLEN),save :: PfssDir
  real(kind=DBL_KIND),save :: Phi_Offset ! offset of longitude for synoptic map for Carrington longitude =0 when t=0
  real(kind=DBL_KIND),save :: Time_Offset ! offset of time of simulation in terms of JD2000
  real(kind=DBL_KIND),parameter :: BR_ENHANCE = 1.d1 ! enhancement factor of Br

!!$  integer(kind=LLONG_KIND),save :: StartDateInUnixTime
!!$  character(len=CHARLEN_YYYYMMDD) :: CARRINGTON_ROTATION_START_DATE = '18531109'
!!$  real(kind=DBL_KIND),save :: CarringtonRotationStartTime
  character(len=CHARLEN_YYYYMMDD) :: DATE_UNIDEF = '        '
  integer(kind=LLONG_KIND),save :: StepForSw = Undefi ! step of Sw
!!$  character(len=CHARLEN) :: UNIX_CMD_DATE = '/bin/date'
  
  ! Type for Soar wind
  type t_synop_chart
     character(len=CHARLEN_YYYYMMDD) :: date
     real(Kind=DBL_KIND) :: time
     real(kind=DBL_KIND),dimension(Nphi, NTheta) :: vr, vp, br, bp, rho, p
  end type t_synop_chart
  ! Sw_a, Sw_b ... solar wind at the key stages
  ! Sw ........... temporlally inerplated solar wind
  type(t_synop_chart),save :: Sw_a, Sw_b, Sw


  public :: pfss_boundary_gid, pfss_init_cond_gid
  public :: pfss_test           ! debug
contains
  ! ----------------------------------------------------------------------------
  ! Innner boundary condition according to solar wind model
  ! ----------------------------------------------------------------------------
  subroutine pfss_boundary_gid(gid, t)
    use grid, only : Step, LevelMax, get_level
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),intent(IN) :: t
    call pfss_init
    
    if ( get_level(gid) /= LevelMax ) return ! solar wind is only for finest grid.
    if (Step(LevelMax) /= StepForSw) then ! if this call is the first time for this time level
       call provide_sw_by_time(t)
       StepForSw = Step(LevelMax)
    end if
    
    call sw_interp_in_space(gid, t)
  end subroutine pfss_boundary_gid
  ! ----------------------------------------------------------------------------
  ! Returns .TRUE. if grid touches the inner boundary.
  ! ----------------------------------------------------------------------------
  function touches_innerboundary(gid) result(bool)
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
  end function touches_innerboundary
  ! ----------------------------------------------------------------------------
  ! Check wheather key stages should be updated.
  ! ----------------------------------------------------------------------------
  function check_update_keystages(t) result(bool)
    real(kind=DBL_KIND),intent(IN) :: t
    logical,save :: bool_initial_stage = .TRUE.
    logical :: bool
    if ( bool_initial_stage ) then
       bool = .TRUE.
       bool_initial_stage = .FALSE.
       return
    end if
    if ( Sw_a%time <= t .and. t < Sw_b%time ) then ! time is between key strages?
       bool = .FALSE.
    else
       bool = .TRUE.
    end if
!!$    print *, bool, t, Sw_a%time , Sw_b%time 
  end function check_update_keystages
  ! ----------------------------------------------------------------------------
  ! Provide solar wind, Sw, in accordance with time level.  The key stages are
  ! also updated if necessary.
  ! ----------------------------------------------------------------------------
  subroutine provide_sw_by_time(t)
    use dates, only : dates_time2yyyymmdd
    real(kind=DBL_KIND),intent(IN) :: t
    integer(kind=LLONG_KIND),parameter ::  oneDay = 1
    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd_a, yyyymmdd_b
    if (check_update_keystages(t)) then
       yyyymmdd_a = dates_time2yyyymmdd(t)
       yyyymmdd_b = add_yyyymmdd(yyyymmdd_a, oneDay)

       ! provide stage a
       if (Sw_b%date == DATE_UNIDEF) then
          call keystage_provide(yyyymmdd_a, Sw_a)
       else
          Sw_a = Sw_b              ! copy
       end if

       ! provide stage b
       call keystage_provide(yyyymmdd_b, Sw_b)
    end if

    ! interpolation in time of Sw_a and Sw_b to make Sw
    call sw_interp_in_time(t)

  end subroutine provide_sw_by_time
  ! ----------------------------------------------------------------------------
  ! Provide all the variable of synoptic chart for a key stage.
  ! ----------------------------------------------------------------------------
  subroutine keystage_provide(yyyymmdd, sw)
    use mpilib
    use dates, only : dates_yyyymmdd2time
    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
    character(len=CHARLEN) :: file
    type(t_synop_chart),intent(OUT) :: sw
    real(kind=DBL_KIND),dimension(Nphi, NTheta) :: brss_gauss, vrss_kms
    if (get_myrank() == PRIMARY_RANK) then
       file = get_filename_by_yyyymmdd(yyyymmdd)
       call keystage_read(file, brss_gauss, vrss_kms)
    end if
    call keystage_bcast(brss_gauss, vrss_kms)
    call keystage_make_sw(brss_gauss, vrss_kms, sw)
    sw%date = yyyymmdd
    sw%time = dates_yyyymmdd2time(yyyymmdd)
  end subroutine keystage_provide
  ! ----------------------------------------------------------------------------
  ! Provide solar wind variables on the synoptic chart.
  ! ----------------------------------------------------------------------------
  subroutine keystage_make_sw(brss_gauss, vrss_kms, sw)
    use parameter
    use modelParameter
    use unit
    real(kind=DBL_KIND),dimension(Nphi, NTheta),intent(IN) :: brss_gauss, vrss_kms
    type(t_synop_chart),intent(OUT) :: sw
    integer :: i, j
    real(kind=DBL_KIND) :: rsun50, ncc, tsw, omega
    rsun50 = 50.d0 * MP_SolarRadius
    omega = Pi2 / (MP_CarringtonRotDay/Unit_day)
    do j = 1, NTheta
       do i = 1, Nphi
          sw%vr(i,j) = vrss_kms(i,j) / Unit_kms
!!$          Vpsw(i,j) = 2 * MP_rInnerBoundary * omega
          sw%vp(i,j) = 0.d0

          ! Kataoka et al. (2009), eq. (3)
!!$          Brsw(i,j) = 10.d0* sign(1.d0, BrSS(i,j))*sqrt(abs(BrSS(i,j))/Unit_utesla) / (MP_rInnerBoundary/MP_SourceSurface)**2
          ! Simple formula  Br ~ r^(-2)
          sw%br(i,j) = brss_gauss(i,j) / Unit_gauss / (MP_rInnerBoundary/MP_SourceSurface)**2
          sw%br(i,j) = sw%br(i,j) * BR_ENHANCE

          ! Kataoka et al. (2009), eq. (6)
          sw%bp(i,j) = -sw%br(i,j) * MP_rInnerBoundary * omega * sin(ThetaSS(j))/sw%vr(i,j)

          ! Shiota et al. (2014), eq. (11), or Hayashi et al. (2003), eq. (10)
          ncc = (rsun50/MP_rInnerBoundary)**2 * (62.98d0 + 866.4d0 * (vrss_kms(i,j)/100 - 1.549)**(-3.402))
          sw%rho(i,j) = ncc / Unit_n

          ! Shiota et al. (2014), eq. (12), or Hayashi et al. (2003), eq. (11)
          tsw = (rsun50/MP_rInnerBoundary)**(2*(MP_Gamma-1)) * (-0.455d0 + 0.1943 * vrss_kms(i,j)/100)*1.D6 ! in K
          sw%p(i,j) = ncc/cgs_mu * cgs_kb * tsw / Unit_e
       end do
    end do
  end subroutine keystage_make_sw
  ! ----------------------------------------------------------------------------
  ! read synop chart for all ranks
  ! ----------------------------------------------------------------------------
  subroutine keystage_bcast(br_ss, vr_ss)
    use mpilib
    real(kind=DBL_KIND),dimension(Nphi, NTheta),intent(INOUT) :: br_ss, vr_ss
    call mpi_bcast(br_ss, size(br_ss), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vr_ss, size(vr_ss), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
  end subroutine keystage_bcast
  ! ----------------------------------------------------------------------------
  ! Read synoptic chart from file
  ! br_ss .... Br on source surface in unit of Gauss
  ! vr_ss .... Vr on source surface in unit of km/s
  ! ----------------------------------------------------------------------------
  subroutine keystage_read(file, br_ss, vr_ss)
    character(len=*),intent(IN) :: file
    integer,parameter :: UNIT = 11
    real(kind=DBL_KIND),dimension(Nphi, NTheta),intent(OUT) :: br_ss, vr_ss
    real(kind=DBL_KIND),dimension(Nphi, NTheta) :: ff_ss, br_ps, bb_ps
    print *, 'pfss file reading... ', trim(file)
    open(UNIT, file=file, form='unformatted',convert='big_endian')
    read(UNIT) br_ss, vr_ss, ff_ss, br_ps, bb_ps
    call flush(UNIT)
    close(UNIT)
  end subroutine keystage_read
  ! ----------------------------------------------------------------------------
  ! temporal interpolation of solar winds at the key stages, Sw_a, and Sw_b, to make Sw
  ! ----------------------------------------------------------------------------
  subroutine sw_interp_in_time(t)
    real(kind=DBL_KIND),intent(IN) :: t
    real(kind=DBL_KIND) :: a, b, dt
    dt = Sw_b%time - Sw_a%time
    a = (Sw_b%time - t)/dt
    b = (t - Sw_a%time)/dt
    Sw%time = t
    Sw%vr  = a * Sw_a%vr  + b * Sw_b%vr
    Sw%vp  = a * Sw_a%vp  + b * Sw_b%vp
    Sw%br  = a * Sw_a%br  + b * Sw_b%br
    Sw%bp  = a * Sw_a%bp  + b * Sw_b%bp
    Sw%rho = a * Sw_a%rho + b * Sw_b%rho
    Sw%p   = a * Sw_a%p   + b * Sw_b%p
  end subroutine sw_interp_in_time
  ! ----------------------------------------------------------------------------
  ! Spatail interpolation of solar wind to remap to the inner boundary
  ! ----------------------------------------------------------------------------
  subroutine sw_interp_in_space(gid, t)
    use modelParameter, only : MP_rInnerBoundary, MP_drInnerBoundary, MP_CarringtonRotDay
    use unit, only : Unit_day, Unit_n, Unit_e, cgs_mu, cgs_kb
    use parameter
    use grid
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),intent(IN) :: t
    integer :: i, j, k, ipha, iphb, jtha, jthb
    real(kind=DBL_KIND) :: vr, vp, vt, br, bp, bt, r, theta, phi, dtheta, dphi, dpha, dphb, dtha, dthb, s2d, sint, cost, sinp, cosp, omega, nh0, rho0, p0, temp0, phi0
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, db
    
    if (.not. touches_innerboundary(gid)) return
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
    phi0 = - Phi_Offset - (t - Time_Offset) * omega !  offset of phi and shift according to rotation
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
                phi = phi + phi0
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
                     Sw%vr(iphb,jthb)*dpha*dtha + Sw%vr(ipha,jthb)*dphb*dtha + &
                     Sw%vr(iphb,jtha)*dpha*dthb + Sw%vr(ipha,jtha)*dphb*dthb ) / s2d
                vp = ( &
                     Sw%vp(iphb,jthb)*dpha*dtha + Sw%vp(ipha,jthb)*dphb*dtha + &
                     Sw%vp(iphb,jtha)*dpha*dthb + Sw%vp(ipha,jtha)*dphb*dthb ) / s2d
                vt = 0.d0
                vx(i,j,k) = vr*sint*cosp + vt*cost*cosp - vp*sinp
                vy(i,j,k) = vr*sint*sinp + vt*cost*sinp + vp*cosp
                vz(i,j,k) = vr*cost - vt*sint

                br = ( &
                     Sw%br(iphb,jthb)*dpha*dtha + Sw%br(ipha,jthb)*dphb*dtha + &
                     Sw%br(iphb,jtha)*dpha*dthb + Sw%br(ipha,jtha)*dphb*dthb ) / s2d
                bp = ( &
                     Sw%bp(iphb,jthb)*dpha*dtha + Sw%bp(ipha,jthb)*dphb*dtha + &
                     Sw%bp(iphb,jtha)*dpha*dthb + Sw%bp(ipha,jtha)*dphb*dthb ) / s2d
                bt = 0.d0
                bx(i,j,k) = br*sint*cosp + bt*cost*cosp - bp*sinp
                by(i,j,k) = br*sint*sinp + bt*cost*sinp + bp*cosp
                bz(i,j,k) = br*cost - bt*sint

                rho(i,j,k) = ( &
                     Sw%rho(iphb,jthb)*dpha*dtha + Sw%rho(ipha,jthb)*dphb*dtha + &
                     Sw%rho(iphb,jtha)*dpha*dthb + Sw%rho(ipha,jtha)*dphb*dthb ) / s2d
                p(i,j,k) = ( &
                     Sw%p(iphb,jthb)*dpha*dtha + Sw%p(ipha,jthb)*dphb*dtha + &
                     Sw%p(iphb,jtha)*dpha*dthb + Sw%p(ipha,jtha)*dphb*dthb ) / s2d
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
    
  end subroutine sw_interp_in_space
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
    real(kind=DBL_KIND) :: rho0, p0, vr0, br0, vra, ra, phi0
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, db

    call pfss_init

    t = Time(LevelMax)                    ! initial condition

    if (Step(LevelMax) /= StepForSw) then ! if this call is the first time for this time level
       call provide_sw_by_time(t)
       StepForSw = Step(LevelMax)
    end if

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
    phi0 = - Phi_Offset - (t - Time_Offset) * omega !  offset of phi and shift according to rotation
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
             phi = phi + phi0
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
                  Sw%vr(iphb,jthb)*dpha*dtha + Sw%vr(ipha,jthb)*dphb*dtha + &
                  Sw%vr(iphb,jtha)*dpha*dthb + Sw%vr(ipha,jtha)*dphb*dthb ) / s2d
             vr = ( vr0 * (ra - r) + vra * (r - MP_rInnerBoundary))/(ra - MP_rInnerBoundary)
             vp = 0.d0
             vt = 0.d0

             vx(i,j,k) = vr*sint*cosp + vt*cost*cosp - vp*sinp
             vy(i,j,k) = vr*sint*sinp + vt*cost*sinp + vp*cosp
             vz(i,j,k) = vr*cost - vt*sint


             br0 = ( &
                  Sw%br(iphb,jthb)*dpha*dtha + Sw%br(ipha,jthb)*dphb*dtha + &
                  Sw%br(iphb,jtha)*dpha*dthb + Sw%br(ipha,jtha)*dphb*dthb ) / s2d
             br = br0 * (MP_rInnerBoundary/r)**2
             bp = 0.d0
             bt = 0.d0

             bx(i,j,k) = br*sint*cosp + bt*cost*cosp - bp*sinp
             by(i,j,k) = br*sint*sinp + bt*cost*sinp + bp*cosp
             bz(i,j,k) = br*cost - bt*sint

             rho0 = ( &
                  Sw%rho(iphb,jthb)*dpha*dtha + Sw%rho(ipha,jthb)*dphb*dtha + &
                  Sw%rho(iphb,jtha)*dpha*dthb + Sw%rho(ipha,jtha)*dphb*dthb ) / s2d
             rho(i,j,k) = rho0 * MP_rInnerBoundary**2 * vr0 / (r**2 * vr)
             p0 = ( &
                  Sw%p(iphb,jthb)*dpha*dtha + Sw%p(ipha,jthb)*dphb*dtha + &
                  Sw%p(iphb,jtha)*dpha*dthb + Sw%p(ipha,jtha)*dphb*dthb ) / s2d
             p(i,j,k) = p0 * (rho(i,j,k)/rho0) ! **MP_Gamma
             db(i,j,k) = 0.d0
          end do
       end do
    end do

!!$    rho = 1.d0                   ! debug
!!$    vx = 0.d0
!!$    vy = 0.d0
!!$    vz = 0.d0
!!$    bx = 0.d0
!!$    by = 0.d0
!!$    bz = 0.d0
!!$    p = 0.d0
  end subroutine pfss_init_cond_gid
!!$  ! ----------------------------------------------------------------------------
!!$  ! conversion among three types of time
!!$  !   time .... A simulation time in unit of hour. The clock starts at the initial condition.
!!$  !             A offset between time and unixtime is StartDateInUnixTime (in unit of second).
!!$  ! 
!!$  !   unixtime .... A unixtime in unit of second. It begins 1970/01/01 00:00:00z. 
!!$  ! 
!!$  !   yyyymmdd .... Characters showing a date. The unix date command is used for conversion
!!$  !                 between yyyymmdd to the other types of time.
!!$  ! ----------------------------------------------------------------------------
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert a stirng yyyymmdd to an unitxtime in sec.
!!$  ! ----------------------------------------------------------------------------
!!$  function yyyymmdd2unixtime(yyyymmdd) result(unixtime)
!!$    use string, only : CHARLEN
!!$    use systemcall, only : systemcall_command
!!$    use io_util, only : read_env
!!$    use mpilib
!!$    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
!!$    integer(kind=LLONG_KIND) :: unixtime
!!$    character(len=*),parameter :: FDATE = 'date_stdout.txt'
!!$    character(len=CHARLEN) :: dir, tmpfile
!!$    integer,parameter :: UNIT = 11
!!$    logical :: tsafe
!!$    if (get_myrank() == PRIMARY_RANK) then
!!$       call read_env('DIR', dir)
!!$       tmpfile = trim(dir)//FDATE
!!$       call systemcall_command( trim(UNIX_CMD_DATE)//" -u -d "//yyyymmdd//"  +'%s' > "//trim(tmpfile))
!!$       open(UNIT, file=tmpfile)
!!$       read(UNIT, *) unixtime
!!$       call flush(UNIT)
!!$       close(unit=UNIT, status='DELETE')
!!$    endif
!!$    call mpi_bcast( unixtime, 1, MPI_LLONG_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
!!$  end function yyyymmdd2unixtime
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert an unix time to a string yyyymmdd
!!$  ! ----------------------------------------------------------------------------
!!$  function unixtime2yyyymmdd(unixtime) result(yyyymmdd)
!!$    use string, only : CHARLEN, num2char
!!$    use systemcall, only : systemcall_command
!!$    use io_util, only : read_env
!!$    use mpilib
!!$    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd
!!$    integer(kind=LLONG_KIND),intent(IN) :: unixtime
!!$    character(len=*),parameter :: FDATE = 'date_stdout.txt'
!!$    character(len=CHARLEN) :: dir, tmpfile
!!$    integer,parameter :: UNIT = 11
!!$    logical :: tsafe
!!$    if (get_myrank() == PRIMARY_RANK) then 
!!$       call read_env('DIR', dir)
!!$       tmpfile = trim(dir)//FDATE
!!$       call systemcall_command( trim(UNIX_CMD_DATE)//" -u -d @"//trim(num2char(unixtime))//"  +'%Y%m%d' > "//trim(tmpfile))
!!$       open(UNIT, file=tmpfile)
!!$       read(UNIT, *) yyyymmdd
!!$       call flush(UNIT)
!!$       close(unit=UNIT, status='DELETE')
!!$    endif
!!$    call mpi_bcast( yyyymmdd, len(yyyymmdd), MPI_CHARACTER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
!!$  end function unixtime2yyyymmdd
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert time in the simulation unit to a date in the form of yyyymmdd
!!$  ! The resul is rounded within the unit of a day.
!!$  ! ----------------------------------------------------------------------------
!!$  function time2yyyymmdd(time) result(yyyymmdd)
!!$    use unit, only : Unit_sec
!!$    real(kind=DBL_KIND),intent(IN) :: time
!!$    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd
!!$    integer(kind=LLONG_KIND) :: unixtime
!!$    unixtime = time2unixtime(time)
!!$    yyyymmdd = unixtime2yyyymmdd(unixtime)
!!$  end function time2yyyymmdd
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert a string yyyymmdd to simulation time
!!$  ! ----------------------------------------------------------------------------
!!$  function yyyymmdd2time(yyyymmdd) result(time)
!!$    use unit, only : Unit_sec
!!$    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
!!$    real(kind=DBL_KIND) :: time
!!$    integer(kind=LLONG_KIND) :: unixtime
!!$    unixtime = yyyymmdd2unixtime(yyyymmdd)
!!$    time = unixtime2time(unixtime)
!!$  end function yyyymmdd2time
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert time (hours from simulation start) to unixtime (sec)
!!$  ! ----------------------------------------------------------------------------
!!$  function time2unixtime(time) result(unixtime)
!!$    use unit, only : Unit_sec
!!$    real(kind=DBL_KIND),intent(IN) :: time
!!$    integer(kind=LLONG_KIND) :: time_in_sec, unixtime
!!$    time_in_sec = time * Unit_sec ! cast from double precision to long integer
!!$    unixtime = time_in_sec + StartDateInUnixTime
!!$  end function time2unixtime
!!$  ! ----------------------------------------------------------------------------
!!$  ! Convert unixtime to time (hours from simulation start)
!!$  ! ----------------------------------------------------------------------------
!!$  function unixtime2time(unixtime) result(time)
!!$    use unit, only : Unit_sec
!!$    integer(kind=LLONG_KIND),intent(IN) :: unixtime
!!$    real(kind=DBL_KIND) :: time
!!$    time = dble(unixtime - StartDateInUnixTime) / Unit_sec
!!$  end function unixtime2time
  ! ----------------------------------------------------------------------------
  ! add days to yyyymmdd
  ! ----------------------------------------------------------------------------
  function add_yyyymmdd(yyyymmdd, days) result(yyyymmdd_out)
    use unit, only : Unit_sec, Unit_day
    use dates, only : dates_time2yyyymmdd, dates_yyyymmdd2time
    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
    integer(kind=LLONG_KIND),intent(IN) :: days
    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd_out
    integer,parameter :: margin_sec = 2 ! for leap second
    real(kind=DBL_KIND) :: margin
    margin = margin_sec / Unit_sec ! to non-dimensional time (simulation time)
    yyyymmdd_out = dates_time2yyyymmdd( dates_yyyymmdd2time(yyyymmdd) + days/Unit_day + margin)
  end function add_yyyymmdd
  ! ----------------------------------------------------------------------------
  ! get filename (path) of a synoptic chart given by a date of yyyymmdd.
  ! ----------------------------------------------------------------------------
  function get_filename_by_yyyymmdd(yyyymmdd) result(fn)
    character(len=CHARLEN_YYYYMMDD),intent(IN) :: yyyymmdd
    character(len=CHARLEN) :: fn
    character(len=*),parameter :: PREFIX = 'PFN_OPEN.GONG.'
    character(len=*),parameter :: SUFFIX = '.dat'
    fn = trim(PfssDir)//PREFIX//yyyymmdd(3:CHARLEN_YYYYMMDD)//SUFFIX
  end function get_filename_by_yyyymmdd
  ! ----------------------------------------------------------------------------
  ! Provide spherical coordinates for a synoptic map
  ! ----------------------------------------------------------------------------
  subroutine spherical_coordinates_provide
    use parameter, only : Pi
    integer :: i, j
    do i = lbound(ThetaSS, 1), ubound(ThetaSS, 1)
       ThetaSS(i) = (i - 1)* ThetaSkip + ThetaStart
    end do
    do j = lbound(PhiSS, 1), ubound(PhiSS, 1)
       PhiSS(j) = (j - 1)* PhiSkip + PhiStart
    end do
    ! conversion unit from degree to radian
    ThetaSS = ThetaSS * Pi/180.d0
    PhiSS = PhiSS * Pi/180.d0
  end subroutine spherical_coordinates_provide
  ! ----------------------------------------------------------------------------
  ! initialize this module.
  ! ----------------------------------------------------------------------------
  subroutine pfss_init()
    use parameter, only : Pi, Pi2
    use unit, only : Unit_kms, Unit_day
    use modelParameter, only : MP_rInnerBoundary, MP_SourceSurface, MP_CarringtonRotDay, MP_T_Start, MP_T_Last
    use dates, only : ORIGIN_OF_SIM_TIME_IN_JD2K, dates_yyyymmdd2time
    use io_util, only : read_env, print_msg
    real(kind=DBL_KIND),parameter :: VSW_KMS = 4.d2 ! typical solar wind speed in km/s
    if (Bool_Initialized) return
    Bool_Initialized = .TRUE.
    call print_msg('pfss initialized')
    call read_env('PFSS_DIR', PfssDir)
    
    Phi_Offset = ( 180D0 - 153.1044D0 + 2.35D0 ) /180D0 *Pi  &
         - ( MP_rInnerBoundary - MP_SourceSurface ) * Unit_kms / VSW_KMS * Pi2 / ( MP_CarringtonRotDay / Unit_day )
    ! 
    ! Longitude_offset is the separation angle between HGI longitude=0 and Carrington longitude =0 when t=0
    ! It is ported from SUSANOO, bdc_commn_real_daily_rdc.f90 : L113

    Time_Offset = ORIGIN_OF_SIM_TIME_IN_JD2K / Unit_day
    ! Offset of time from JD2000=0

    ! initialize Sw_a and Sw_b
    Sw_a%date = DATE_UNIDEF
    Sw_b%date = DATE_UNIDEF
    call spherical_coordinates_provide
  end subroutine pfss_init
  ! ----------------------------------------------------------------------------
  ! initialize this module.
  ! ----------------------------------------------------------------------------
  subroutine pfss_test
    use mpilib
    character(len=CHARLEN_YYYYMMDD) :: yyyymmdd = '20110228'
!!$    character(len=CHARLEN) :: file
!!$    real(kind=DBL_KIND),dimension(Nphi, NTheta) :: br_ss, vr_ss
    call pfss_init
!!$    call keystage_provide(yyyymmdd, Sw_a)
!!$    call keystage_provide(add_yyyymmdd(yyyymmdd, 1), Sw_b)
!!$    print *, check_update_keystages(1.d0)
    call provide_sw_by_time(0.d0)

  end subroutine pfss_test
  
end module pfss
