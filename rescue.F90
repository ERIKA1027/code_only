!-------------------------------------------------------------------------
! Module for rescue arrays
!-------------------------------------------------------------------------
#include "config.h"
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

  public :: rescueLev, rescueAllLev
  public :: V_UPPER_LIMIT

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
    call setDiffusion( level )
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
    use grid, only : Gidmin, GidListMax, GidList, Mmin, Mmax, get_Ucomp, has_child_grid, get_Xp, get_Yp, get_Zp
    integer,intent(IN) :: level
    integer :: n, gid, m, ncell, i, j, k, is, js, ks, ie, je, ke, ii, jj, kk, nw
    logical :: allIsFinite, isFinite ! function in naninf
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: a
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: ave
    logical :: bool_realgrid, bool_write

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
                        write(*,'("*** rescueNanInf (realg, level, gid, x, y, z, m, nw) ", L, I3, I4, 3(1PE12.4), I3, I3)') &
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
                write(*, '("*** rescue: V limitted (realg, level, gid, x, y, z) ", L, I3, I4, 3(1PE12.4))') &
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
                print *, "(rescue) KS DEBUG", u(i,j,k,MP), u(i,j,k,MRHO), u(i,j,k,MVX), u(i,j,k,MVY), u(i,j,k,MZ),&
                     x(i), y(j), z(k)
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
                     write(*,'("*** phase",I2," rescue diff (realg, level, gid, x, y, z, nw, ntb) ", L, I3, I4, 3(1PE12.4), I3, I3)') &
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
  
end module rescue
