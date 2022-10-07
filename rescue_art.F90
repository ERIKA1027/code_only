!-------------------------------------------------------------------------
! Module for rescue arrays
!-------------------------------------------------------------------------
#include "config.h"
#include "chemistry_label.h"
!-------------------------------------------------------------------------
! Use of rescue
! Define if rescue module is used (recommended).
! 
#define USE_RESCUE
! 
!-------------------------------------------------------------------------
! Level of warning
! VERBOSE ... Output warnings the most frequently.  When this macro is
! not defined, the warnings are shown only for the real grid (the grid
! having no chlid grid).
! SILENT .... Suppress all the warnings.
! 
#define VERBOSE
! #define SILENT
!
!-------------------------------------------------------------------------
! Turn on the sub-step method, restricting CFL condition and taking
! smaller time-step.  This macro exists also in multi_timestep.F90.
! Both the two macro should be define for use of this sub-step method.
! 
! #define USE_RESCUE_RESTRICTCFL
! 
!-------------------------------------------------------------------------
module rescue
  use grid, only: Lmin, Lmax
  implicit none
  private
!!$  real(kind=DBL_KIND),parameter :: V_UPPER_LIMIT = 1.d1 ! for solar wind
!!$  real(kind=DBL_KIND),parameter :: V_UPPER_LIMIT = 2.d2 ! for NL=7 (binary accretion)
  real(kind=DBL_KIND),parameter :: V_UPPER_LIMIT = 1.d2 ! for NL=5

#ifdef USE_RESCUE_RESTRICTCFL
  logical,dimension(Lmin:Lmax),save :: Bool_SubstepOn = .FALSE. ! sw for substep on/off
  logical,dimension(Lmin:Lmax),save :: Bool_SubstepTurnOn = .FALSE. ! true if sw is turned on
  integer,dimension(Lmin:Lmax),save :: N_SubstepCounter = 0     ! counter for substeps
#endif ! USE_RESCUE_RESTRICTCFL

  public :: rescueLev, rescueAllLev, rescue_rhopsi_for_poisson ! KS MODIFIED

#ifdef USE_RESCUE_RESTRICTCFL
  public :: rescue_restrictCFL
#endif ! USE_RESCUE_RESTRICTCFL

contains
  !-----------------------------------------------------------------------
  ! rescue all the component of U in a given grid level
  !-----------------------------------------------------------------------
  subroutine rescueLev(level)
    integer :: level
#ifndef USE_RESCUE
    return
#endif
#ifdef EMULATE_2DIM
    call recover_2dim( level )
#endif !EMULATE_2DIM
    call rescueNanInf( level )
    !call setDiffusion( level )
    call rescueByFloor ( level ) ! KS MODIFIED
#ifdef USE_RESCUE_RESTRICTCFL
    call allreduce_bool_substepon(level)
#endif ! USE_RESCUE_RESTRICTCFL

!!$    call setUpperLimitVelocity( level )

  end subroutine rescueLev
  !-----------------------------------------------------------------------
  ! rescue all the component of U in all the grid level
  !-----------------------------------------------------------------------
  subroutine rescueAllLev
    use grid, only : Lmin, LevelMax
    integer :: level
#ifndef USE_RESCUE
    return
#endif
    do level = Lmin, LevelMax
       call rescueLev( level )
    enddo
  end subroutine rescueAllLev
  !-------------------------------------------------------------------------
  ! rescue data if values have non-finite data
  !-------------------------------------------------------------------------
  subroutine rescueNanInf( level )
    use grid, only : Gidmin, GidListMax, GidList, Mmin, Mmax, get_Ucomp, has_child_grid, get_Xp, get_Yp, get_Zp, get_Up
    use mpilib ! KS ADDED
    integer,intent(IN) :: level
    integer :: n, gid, m, ncell, i, j, k, is, js, ks, ie, je, ke, ii, jj, kk, nw
    logical :: allIsFinite, isFinite ! function in naninf
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: a
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: ave
    logical :: bool_realgrid, bool_write
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u ! KS ADDED

    logical :: isNotFinite ! function in naninf (KS DEBUG)

    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level) ! gid for U
       bool_realgrid = .not. has_child_grid(gid)

       ! control message
#if defined(VERBOSE)
       bool_write = bool_realgrid
#elif defined(SILENT)
       bool_write = .FALSE.
#else
       bool_write = .TRUE.
#endif
       do m = Mmin, Mmax
          a => get_Ucomp(m, gid)
          if ( allIsFinite( a , size(a) ) ) cycle
          ks = lbound(a,3)
          ke = ubound(a,3)
          js = lbound(a,2)
          je = ubound(a,2)
          is = lbound(a,1)
          ie = ubound(a,1)
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)
          do k = ks, ke
             do j = js, je
                do i = is, ie
                   if ( isFinite(a(i,j,k)) ) cycle
                !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!
                   u => get_Up(gid)
                   print '(A, 5I5, (1P11E15.7))', "(rescueNanInf) KS DEBUG", &
                        get_myrank(),gid,i-lbound(u,1),j-lbound(u,2),k-lbound(u,3), &
                        u(i,j,k,MP), u(i,j,k,MRHO), u(i,j,k,MVX), u(i,j,k,MVY), u(i,j,k,MVZ),&
                        x(i), y(j), z(k),u(i,j,k,MXPI), u(i,j,k,MYPI),u(i,j,k,MXPI)
                   if (isNotFinite(u(i,j,k,MP)) .or.  isNotFinite(u(i,j,k,MRHO)) .or. isNotFinite(u(i,j,k,MVX)) &
                        .or. isNotFinite(u(i,j,k,MVY)) .or.  isNotFinite(u(i,j,k,MVZ))) then
                      print '(A,/,A)', '(rescue KS DEBUG) NaN/Inf found','stopping'
                      stop
                   end if
                !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!

#ifdef DM_POTENTIAL
                  if(isNotFinite(u(i,j,k,MDMRHO))) then
                    print *, "DM: NaN/Inf found", get_myrank(),gid,i-lbound(u,1),j-lbound(u,2),k-lbound(u,3)
                    stop
                  endif
#endif


                   
                   do nw = 1, max(ie-is+1, je-js+1, ke-ks+1) ! size of window
                      ave = 0.d0
                      ncell = 0
                      do kk = max(k-nw, ks), min(k+nw, ke)
                         do jj = max(j-nw, js), min(j+nw, je)
                            do ii = max(i-nw, is), min(i+nw, ie)
                               if ( isFinite(a(ii,jj,kk)) ) then
                                  ave = ave + a(ii,jj,kk)
                                  ncell = ncell + 1
                               endif
                            enddo
                         enddo
                      enddo
                      if ( ncell > 0 ) exit
                   enddo
                   if (ncell == 0) then
                      write(*, *) "*** Error in rescue3D, ncell == 0"
                      stop
                   endif
                   a(i,j,k) = ave/ncell
                   if (bool_write) &
                        write(*,'("*** rescueNanInf (realg, level, gid, x, y, z, m, nw) ", L, I3, I4, 3(1PE14.6), I3, I3)') &
                        bool_realgrid, level, gid, x(i), y(j), z(k), m, nw
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine rescueNanInf
  !-------------------------------------------------------------------------
  ! set uppper limit of velocities
  !-------------------------------------------------------------------------
  subroutine setUpperLimitVelocity( level )
    use grid, only : Gidmin, GidListMax, GidList, get_Ucomp, has_child_grid
    integer,intent(IN) :: level
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: v
    real(kind=DBL_KIND) :: vmax
    integer :: n, gid, m
    real(kind=DBL_KIND) :: xloc, yloc, zloc
    logical :: bool_realgrid, bool_write
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       bool_realgrid = .not. has_child_grid(gid)

       ! control message
#if defined(VERBOSE)
       bool_write = bool_realgrid
#elif defined(SILENT)
       bool_write = .FALSE.
#else
       bool_write = .TRUE.
#endif

       do m = MVX, MVZ
          v => get_Ucomp(m, gid)
          vmax = maxval(abs(v))
          if ( vmax >= V_UPPER_LIMIT ) then
             call gridloc(gid, xloc, yloc, zloc)
             if (bool_write) then
                write(*, '("*** rescue: V limitted (realg, level, gid, x, y, z) ", L, I3, I4, 3(1PE14.6))') &
                     bool_realgrid, level, gid, xloc, yloc, zloc
             endif
             where( v > V_UPPER_LIMIT  ) v =  V_UPPER_LIMIT
             where( v < -V_UPPER_LIMIT ) v = -V_UPPER_LIMIT
          end if
       end do
    end do
  end subroutine setUpperLimitVelocity
  !-------------------------------------------------------------------------
  ! cure upper negative density and pressure, and limit velocity by diffusion
  !-------------------------------------------------------------------------
  subroutine setDiffusion( level )
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax, Mmin, Mmax, get_Ucomp, get_Up, Gidmin, GidListMax, GidList, get_Xp, get_Yp, get_Zp, has_child_grid, Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh
#ifdef MODEL_SOLARWIND
    use modelParameter, only: MP_rInnerBoundary
#endif !MODEL_SOLARWIND
    !----------------- KS ADDED ------------------!
    use modelParameter, only : MP_Tmin, MP_Tmax !KS ADDED
    use unit
#if MODEL_ART > 0
    use primordial, only : yHe
#endif    
    !---------------------------------------------!
    integer,intent(IN) :: level
    integer,dimension(3) :: mlistvel = (/MVX, MVY, MVZ/)
#ifdef MP
    integer,dimension(1) :: mlistpos_r = (/MRHO/)
    integer,dimension(1) :: mlistpos_p = (/MP/)
#else !MP
    integer,dimension(1) :: mlistpos_r = (/MRHO/)
#endif !MP
#ifdef MODEL_SOLARWIND
    integer,parameter :: NPHASE_MIN = 2
#else !MODEL_SOLARWIND
    integer,parameter :: NPHASE_MIN = 1
#endif !MODEL_SOLARWIND
    integer,parameter :: NPHASE_MAX = 3
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJKMGH) :: w
    real(kind=DBL_KIND),dimension(Mmin:Mmax) :: wave, uave
    integer :: n, gid, i, j, k, ii, jj, kk, is, ie, js, je, ks, ke, nw, m, nphase, np, ncell, ntrouble, ntrouble_r, ntrouble_p, ntrouble_v
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    logical :: bool_realgrid, bool_write, bool_allok

    logical :: isNotFinite ! function in naninf (KS DEBUG)

    ! real(kind=DBL_KIND) :: xmu, T_K !KS ADDED
!    integer :: count_output=0, max_count=100 !KS ADDED
    ! integer :: count_output=0, max_count=10000 !KS DEBUG

    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       u => get_Up(gid)

       if (&
#ifdef MP
            all(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistpos_p(:)) > 0.d0) .and. & !! check the values are all positive for p
#endif !MP
            all(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistpos_r(:)) > 0.d0) .and. & ! check the values are all positive for rho
            all(abs(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistvel(:))) < V_UPPER_LIMIT)) cycle  ! check the velocities are under limitation.
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       bool_realgrid = .not. has_child_grid(gid)

#ifdef USE_RESCUE_RESTRICTCFL
       if (bool_realgrid) then
          Bool_SubstepTurnOn(level) = .TRUE.
       endif
#endif ! USE_RESCUE_RESTRICTCFL

       ! control message
#if defined(VERBOSE)
       bool_write = bool_realgrid
#elif defined(SILENT)
       bool_write = .FALSE.
#else
       bool_write = .TRUE.
#endif
       ! convert to conservative variables
       w(:,:,:,:) = u(:,:,:,:)
       w(:,:,:,MVX) = u(:,:,:,MVX)*u(:,:,:,MRHO)
       w(:,:,:,MVY) = u(:,:,:,MVY)*u(:,:,:,MRHO)
       w(:,:,:,MVZ) = u(:,:,:,MVZ)*u(:,:,:,MRHO)

       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
#ifdef MODEL_SOLARWIND
                if (x(i)**2 + y(j)**2 + z(k)**2 < MP_rInnerBoundary**2) cycle  !skip if boundary region
#endif !MODEL_SOLARWIND
                ntrouble_p = 0    ! type of trouble
                ntrouble_r = 0    ! type of trouble
                ntrouble_v = 0    ! type of trouble               
#ifdef MP
                if (.not. all(u(i,j,k,mlistpos_p(:)) > 0.d0)) ntrouble_p = 1
#endif !MP
                if (.not. all(u(i,j,k,mlistpos_r(:)) > 0.d0)) ntrouble_r = 1
                if (.not. all(abs(u(i,j,k,mlistvel(:))) < V_UPPER_LIMIT)) ntrouble_v = 1
                ntrouble = ntrouble_p + ntrouble_r * 2 + ntrouble_v * 4
                if (ntrouble == 0) cycle
                !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!
                print '(A, 4I5,(1P8E15.7))', "(rescue) KS DEBUG", &
                     gid, i-lbound(u,1),j-lbound(u,2),k-lbound(u,3),&
                     u(i,j,k,MRHO), u(i,j,k,MP), u(i,j,k,MVX), u(i,j,k,MVY), u(i,j,k,MVZ),&
                     x(i), y(j), z(k)
                if (isNotFinite(u(i,j,k,MP)) .or.  isNotFinite(u(i,j,k,MRHO)) .or. isNotFinite(u(i,j,k,MVX)) &
                     .or. isNotFinite(u(i,j,k,MVY)) .or.  isNotFinite(u(i,j,k,MVZ))) then
                   print '(A,/,A)', '(rescue KS DEBUG) NaN/Inf found','stopping'
                   stop
                end if
                !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!                
                do nw = 1, max(Imax-Imin+1, Jmax-Jmin+1, Kmax-Kmin+1) ! expand wind size
                   do np = NPHASE_MIN, NPHASE_MAX ! phase counter
                      nphase = np
                      is = max(i-nw, Imin)
                      ie = min(i+nw, Imax)
                      js = max(j-nw, Jmin)
                      je = min(j+nw, Jmax)
                      ks = max(k-nw, Kmin)
                      ke = min(k+nw, Kmax)

                      if (nphase == 1) then ! simple average around the cell (i,j,k)
                         do m = Mmin, Mmax
                            wave(m) = sum( w(is:ie, js:je, ks:ke, m) )
                         end do
                         wave(:) = wave(:) / ((ie-is+1)*(je-js+1)*(ke-ks+1))
                         uave = wave
                         uave(MVX:MVZ) = wave(MVX:MVZ)/wave(MRHO)
                      else if (nphase == 2) then ! same as phase 1 but excluding the central point from the average 
                         do m = Mmin, Mmax
                            wave(m) = sum( w(is:ie, js:je, ks:ke, m) )
                         end do
                         wave(:) = wave(:) - w(i,j,k,:)
                         wave(:) = wave(:) / ((ie-is+1)*(je-js+1)*(ke-ks+1)-1)
                         uave = wave
                         uave(MVX:MVZ) = wave(MVX:MVZ)/wave(MRHO)
                      else if (nphase == 3) then ! average only for cells with regular values. use u directly
                         uave = 0.d0
                         ncell = 0
                         do kk = ks, ke
                            do jj = js, je
                               do ii = is, ie
                                  if (ii == i .and. jj == j .and. kk == k) cycle
                                  if (&
#ifdef MP
                                       all(u(ii,jj,kk,mlistpos_p(:)) > 0.d0) .and. &
#endif !MP
                                       all(u(ii,jj,kk,mlistpos_r(:)) > 0.d0) .and. &
                                       all(abs(u(ii,jj,kk,mlistvel(:))) < V_UPPER_LIMIT) &
                                       ) then
                                     uave(:) = uave(:) + u(ii,jj,kk,:)
                                     ncell = ncell + 1
                                  end if
                               end do
                            end do
                         end do
                         if (ncell > 0) then
!!$                            print *, 'ncell', ncell
                            uave(:) = uave(:)/ncell
                         endif
                      end if
                      bool_allok = ( &
#ifdef MP
                           all(uave(mlistpos_p(:)) > 0.d0) .and. &
#endif !MP
                           all(uave(mlistpos_r(:)) > 0.d0) .and. &
                           all(abs(uave(mlistvel(:))) < V_UPPER_LIMIT))
                      if ( bool_allok ) exit
                   end do       ! try next phase

                   if ( bool_allok ) then
                      if (nphase == 1) then
                         do m = Mmin, Mmax ! fill all the cells which are used in average
                            u(is:ie,js:je,ks:ke,m)  = uave(m)
                         end do
                      else      ! cure only a cell
#ifdef MP
                         if (ntrouble_p == 1) u(i,j,k,mlistpos_p) = uave(mlistpos_p) ! cure only p
#endif !MP
                         if (ntrouble_r == 1 .or. ntrouble_v == 1) u(i,j,k,:) = uave(:) ! cure all component
                      endif
                      exit
                   endif

                end do          ! try larger window

                if ( .not. bool_allok  ) then ! phase 4
                   nphase = nphase + 1
#ifdef MP
                   if (.not. all(uave(mlistpos_p(:)) > 0.d0)) then
                      u(i,j,k,mlistpos_p) = abs(w(i,j,k,mlistpos_p))
                   endif
#endif !MP
                   if (.not. all(uave(mlistpos_r(:)) > 0.d0)) then
                      u(i,j,k,mlistpos_r) = abs(w(i,j,k,mlistpos_r))
                      u(i,j,k,MVX:MVZ) = w(i,j,k,MVX:MVZ)/u(i,j,k,MRHO)
                   endif
                   if (.not. all(abs(u(i,j,k,mlistvel(:))) < V_UPPER_LIMIT) ) then
                      u(i,j,k,mlistvel) = w(i,j,k,mlistvel)/sqrt(sum(w(i,j,k,mlistvel)**2)) * V_UPPER_LIMIT
                   endif
                   nw = 1
                endif
                if (bool_write) &
                     write(*,'("*** phase",I2," rescue diff (realg, level, gid, x, y, z, nw, ntb) ", L, I3, I4, 3(1PE14.6), I3, I3)') &
                     nphase, bool_realgrid, level, gid, x(i), y(j), z(k), nw, ntrouble
             end do
          end do
       end do
    end do

    ! bool and counter should have same value in all the node

  end subroutine setDiffusion
  !-------------------------------------------------------------------------
  ! get location of grid
  !-------------------------------------------------------------------------
  subroutine gridloc(gid, xloc, yloc, zloc)
    use grid
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),intent(OUT) :: xloc, yloc, zloc
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    xloc = (x(Imin)+x(Imax))*0.5d0
    yloc = (y(Jmin)+y(Jmax))*0.5d0
    zloc = (z(Kmin)+z(Kmax))*0.5d0
  end subroutine gridloc

#ifdef EMULATE_2DIM
  !-------------------------------------------------------------------------
  ! Emulation for 2-dimension
  !-------------------------------------------------------------------------
  subroutine recover_2dim(level)
    use grid, only : Gidmin, GidListMax, GidList, get_Up, Kmin, Kmaxgh
    integer,intent(IN) :: level
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: gid, n, k
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level) ! gid for U
       u => get_Up(gid)
       do k=Kmin+1, Kmaxgh
          u(:,:,k,:) = u(:,:,Kmin,:)
       end do
    end do
  end subroutine recover_2dim
#endif !EMULATE_2DIM

#ifdef USE_RESCUE_RESTRICTCFL

  !-------------------------------------------------------------------------
  ! all reduce bools for substepon
  !-------------------------------------------------------------------------
  subroutine allreduce_bool_substepon(level)
    use mpilib
    use grid, only: Lmin, Lmax
    integer,intent(IN) :: level                  ! the grid level
    call MPI_Allreduce(MPI_IN_PLACE, Bool_SubstepTurnOn(level), 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    if (Bool_SubstepTurnOn(level)) then
       Bool_SubstepOn(level) = .TRUE.
       N_SubstepCounter(level) = 0
       if ( get_myrank() == PRIMARY_RANK ) then
          write(*,'("*** Switch turn ON for restrictCFL ", I3, I3)') level, N_SubstepCounter(level)
       end if
    end if
    Bool_SubstepTurnOn(level) = .FALSE.
  end subroutine allreduce_bool_substepon
  !-------------------------------------------------------------------------
  ! restrict cfl condition
  !-------------------------------------------------------------------------
  subroutine rescue_restrictCFL(dtlocal, level)
    use mpilib
    real(kind=DBL_KIND),intent(INOUT) :: dtlocal ! dt on the grid level
    integer,intent(IN) :: level                  ! the grid level
    integer,parameter :: NSUBSTEP = 10 ! number of substep
    if (.not. Bool_SubstepOn(level)) return
    dtlocal = dtlocal / NSUBSTEP
    N_SubstepCounter(level) = N_SubstepCounter(level) + 1

    ! Last substep?
    if (N_SubstepCounter(level) >= NSUBSTEP) then
       N_SubstepCounter(level) = 0
       Bool_SubstepOn(level) = .FALSE.
    endif

    ! messages
    if ( get_myrank() == PRIMARY_RANK ) then
       write(*,'("*** Switch ON for restrictCFL ", I3, I3)') level, N_SubstepCounter(level)
    end if
  end subroutine rescue_restrictCFL
  
#endif ! USE_RESCUE_RESTRICTCFL

  !----------------------------------------------------------------------------------------!
  ! rescue U1 (KS ADEED)                                                                   !
  !                                                                                        !
  ! same as original rescue subroutines but for U1 (quantities at t_(n+1/2) = tn + dt/2)   !
  !----------------------------------------------------------------------------------------!
  !
  !-----------------------------------------------------------------------
  ! rescue all the component of U1 in all the grid level
  !-----------------------------------------------------------------------
  subroutine rescue_rhopsi_for_poisson
    use grid, only : Lmin, LevelMax
    integer :: level
    do level = Lmin, LevelMax
       call check_rhopsi_for_poisson( level )
    enddo
  end subroutine rescue_rhopsi_for_poisson

  !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!
  subroutine check_rhopsi_for_poisson( level )
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax, get_Ucomp, get_U1comp, get_U2comp, Gidmin, GidListMax, GidList, Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh,&
         get_Xp, get_Yp, get_Zp
    use mpilib
    use unit
    integer,intent(IN) :: level
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho,rho1,rho2,psi
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    integer :: n, gid, i, j, k
    logical :: isNotFinite ! function in naninf (KS DEBUG)

    ! if ( get_myrank() == PRIMARY_RANK ) &
    !      print *, "KS DEBUG: check_rhopsi_for_poisson, level=",level
    ! if ( get_myrank() == PRIMARY_RANK ) then
    !    print *, "KS DEBUG ** WARNING ** change psi == 0 for initial guess of poisson eq."
    !    print *, "KS DEBUG ** WARNING ** change rho1 = rho0*r_au**-2 for poisson eq."
    ! end if
    
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       ! rho => get_Ucomp(MRHO,gid)
       rho1 => get_U1comp(MRHO,gid)
       ! rho2 => get_U2comp(MRHO,gid)
       psi => get_Ucomp(MPSI,gid)
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)       
       do k = Kmingh, Kmaxgh
          do j = Jmingh, Jmaxgh
             do i = Imingh, Imaxgh
                !------ WARNING ------- WARNING ------- WARNING ------- WARNING -------!
                ! psi(i,j,k) = 0.d0 * psi(i,j,k)
                ! rho1(i,j,k) = max(1d12*cgs_mh/Unit_rho / ((x(i)**2+y(j)**2+z(k)**2)*Unit_au**2), 1d1*cgs_mh/Unit_rho)
                !------ WARNING ------- WARNING ------- WARNING ------- WARNING -------!
                !check rho of U1
                if (.not. (rho1(i,j,k)*Unit_rho/cgs_mh > 1d-2 .and. rho1(i,j,k)*Unit_rho/cgs_mh < 1d16)) then !密度が変になっていそうなとき
                   print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "** WARNING ** (check U1_rho) ", &
                        rho1(i,j,k)*Unit_rho/cgs_mh, i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                end if
                !check psi of U
                if (.not. (psi(i,j,k) > -1d2 .and. psi(i,j,k) < 1d0)) then          !重力場がおかしくなっていそうなとき
                   print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "** WARNING ** (check U_psi) ", &
                        psi(i,j,k), i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                end if
                ! !check rho of U
                ! if (.not. (rho(i,j,k)*Unit_rho/cgs_mh > 1d0 .and. rho(i,j,k)*Unit_rho/cgs_mh < 1d13)) then !密度が変になっていそうなとき
                !    print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "(rescue U_rho) ", &
                !         rho(i,j,k)*Unit_rho/cgs_mh, i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                ! end if
                ! !check rho of U2
                ! if (.not. (rho2(i,j,k)*Unit_rho/cgs_mh > 1d0 .and. rho2(i,j,k)*Unit_rho/cgs_mh < 1d13)) then !密度が変になっていそうなとき
                !    print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "(rescue U_rho) ", &
                !         rho2(i,j,k)*Unit_rho/cgs_mh, i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                ! end if

             end do
          end do
       end do
    end do
  end subroutine check_rhopsi_for_poisson
  !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!

  !-------------------------------------------------------------------------
  ! rescue data by imposing floor
  !-------------------------------------------------------------------------
  subroutine rescueByFloor( level )
    use grid, only : Gidmin, GidListMax, GidList, get_Ucomp, has_child_grid, get_Xp, get_Yp, get_Zp, Imin, Imax, Jmin, Jmax, Kmin, Kmax,&
         Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh
#ifdef METAL_TRANSFER
    use modelParameter, only : MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_mu, MP_frac_COsum
#else
    use modelParameter, only : MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_metallicity, MP_mu, MP_frac_COsum
#endif
    use unit
    use primordial
    use mpilib
    integer,intent(IN) :: level
    integer :: n, gid, i, j, k, is, ie, js, je, ks, ke
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    logical :: bool_realgrid, bool_gh
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, vx, vy, vz
#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: yco
#endif
    real(kind=DBL_KIND) :: T_K, xmu, v2, xnH, yco_l
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
#ifdef METAL_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
#endif
    integer :: count_output_NL=0, count_output_TL=0, count_output_TU=0, count_output_VU=0, max_count=100 !KS DEBUG
    ! ----------------------------------------------------------
    type chemsp 
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(NCEHM_MIN:NCEHM_MAX)
    integer :: ichem
    ! ----------------------------------------------------------
    

    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level) ! gid for U
       bool_realgrid = .not. has_child_grid(gid)
       !座標
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)

       !変数の取り出し
       rho => get_Ucomp(MRHO,gid)
       p => get_Ucomp(MP,gid)

       ! get chemistry pointer
       do ichem = NCEHM_MIN, NCEHM_MAX
         chemary3(ichem)%y => get_Ucomp(ichem,gid)
       enddo

#ifdef CHEM_MODEL_HF2020
       yco  => get_Ucomp(MCO, gid)
#endif

       vx => get_Ucomp(MVX,gid)
       vy => get_Ucomp(MVY,gid)
       vz => get_Ucomp(MVZ,gid)

#ifdef METAL_TRANSFER
       mmetal => get_Ucomp(MMET,gid)
#endif

       ! ks = lbound(rho,3)
       ! ke = ubound(rho,3)
       ! js = lbound(rho,2)
       ! je = ubound(rho,2)
       ! is = lbound(rho,1)
       ! ie = ubound(rho,1)


       ! do i = Imin, Imax
       !    do j = Jmin, Jmax
       !       do k = Kmin, Kmax
       ! do i = is,ie
       !    do j = js,je
       !       do k = ks,ke
       do i = Imingh, Imaxgh
          do j = Jmingh, Jmaxgh
             do k = Kmingh, Kmaxgh
                
                bool_gh = .not. (i<Imin .or. i>Imax .or. j<Jmin .or. j>Jmax .or. k<Kmin .or. k>Kmax)
                      
                !-------------------------------
                ! impose conservation law in chemistry
                !-------------------------------
                ! --------------------------
                !   化学組成変数の詰め替え
                ! --------------------------
                do ichem = 0, NCHEM-1
                   ychem(ichem) = chemary3(ichem+NCEHM_MIN)%y(i,j,k)
                enddo 
#ifdef CHEM_MODEL_HF2020
                yco_l = yco(i,j,k)
#endif

                !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
                call adjust_abundance(ychem &
#ifdef CHEM_MODEL_HF2020
                , yco_l &
#endif
#ifdef METAL_TRANSFER
                , mmetal(i,j,k))
#else
                , MP_Metallicity)
#endif

                ! ----------------------------------
                !     packing updated chemistry 
                ! ----------------------------------
                do ichem = 0, NCHEM-1
                   chemary3(ichem+NCEHM_MIN)%y(i,j,k) = ychem(ichem)
                enddo 
#ifdef CHEM_MODEL_HF2020
                yco(i,j,k) = yco_l
#endif

                !------------------------------
                ! metallicity floor
                !------------------------------
#ifdef METAL_TRANSFER
                mmetal(i,j,k) = MAX(mmetal(i,j,k), 1.d-99)
#endif

                !-------------------------------
                ! impose density floor
                !-------------------------------
                xnH = rho(i,j,k)*Unit_rho/(MP_mu*cgs_amu)                !水素原子核の数密度
                if (.not. xnH >= MP_Nmin) then
                   !output warning
                   if (count_output_NL <= max_count) then
                      print '(A, 5I5, 2L, I3, 5(1PE14.6))', "*** rescueByFloor (N lower floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, xnH, MP_Nmin, x(i),y(j),z(k)
                      if (count_output_NL == max_count) &
                           print *, 'count_output_NL reaches max_count -> no more warning shown'
                      count_output_NL = count_output_NL+1
                   end if
                   !impose floor
                   xnH = MP_Nmin
                   rho(i,j,k) = xnH * (MP_mu*cgs_amu) / Unit_rho

                   !密度をfixして圧力そのままだと温度がめちゃくちゃになる -> Fixしてしまう
                   xmu = get_xmu(ychem) ! 平均分子量  
                   T_K = sqrt(MP_Tmin*MP_Tmax)                                         ! 幾何平均
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                end if


                !-------------------------------
                ! impose temperature floor
                !-------------------------------
                !温度の取得
                xmu = get_xmu(ychem) ! 平均分子量  
                T_K = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb      !温度 [K]

                if (.not. T_K >= MP_Tmin) then
                   !output warning
                   if (count_output_TL <= max_count) then
                      print '(A, 5I5, 2L, I3, 8(1PE14.6))', "*** rescueByFloor (T lower floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, p(i,j,k), rho(i,j,k), xmu, T_K, MP_Tmin, x(i),y(j),z(k)
                      if (count_output_TL == max_count) &
                           print *, 'count_output_TL reaches max_count -> no more warning shown'
                      count_output_TL = count_output_TL+1
                   end if
                   !impose floor
                   T_K = MP_Tmin
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                end if

                if (.not. T_K <= MP_Tmax) then
                   !output warning
                   if (count_output_TU <= max_count) then
                      print '(A, 5I5, 2L, I3, 8(1PE14.6))', "*** rescueByFloor (T upper floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, p(i,j,k), rho(i,j,k), xmu, T_K, MP_Tmax, x(i),y(j),z(k)
                      if (count_output_TU == max_count) &
                           print *, 'count_output_TU reaches max_count -> no more warning shown'
                      count_output_TU = count_output_TU+1
                   end if
                   !impose floor
                   T_K = MP_Tmax
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                end if

                !-------------------------------
                ! impose velocity floor
                !-------------------------------
                v2 = (vx(i,j,k)*vx(i,j,k) + vy(i,j,k)*vy(i,j,k) + vz(i,j,k)*vz(i,j,k))*Unit_kms*Unit_kms

                if (.not. v2 <= MP_Vmax*MP_Vmax) then
                   !output warning
                   if (count_output_VU <= max_count) then
                      print '(A, 5I5, 2L, I3, 7(1PE14.6))', "*** rescueByFloor (V upper floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, vx(i,j,k)*Unit_kms,vy(i,j,k)*Unit_kms,vz(i,j,k)*Unit_kms, MP_Vmax, x(i),y(j),z(k)
                      if (count_output_VU == max_count) &
                           print *, 'count_output_VU reaches max_count -> no more warning shown'
                      count_output_VU = count_output_VU+1
                   end if
                   !impose floor                   
                   vx(i,j,k) = MP_Vmax/sqrt(v2) * vx(i,j,k)
                   vy(i,j,k) = MP_Vmax/sqrt(v2) * vy(i,j,k)
                   vz(i,j,k) = MP_Vmax/sqrt(v2) * vz(i,j,k)
                end if


                
             end do
          end do
       end do
    enddo
  end subroutine rescueByFloor

end module rescue
