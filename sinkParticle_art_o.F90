#include "config.h"
#define FNSPP_  "sp_param.txt"
!-------------------------------------------------------------------------
! module for sink particle
! This module supports only a single-timestep mode.
!-------------------------------------------------------------------------
! Particles are test particles.
#define ACCRETION_ON
#define ACCRETION_RHOCR
! #define ACCRETION_BONDI
#define CREATION_ON
#define ADVECTION_ON
#define MERGER_ON
! #define TESTPARTICLE

#ifdef TESTPARTICLE
#define ADVECTION_ON
#undef ACCRETION_ON
#undef ACCRETION_RHOCR
#undef ACCRETION_BONDI
#undef MERGER_ON
#endif !TESTPARTICLE

! Particle can be created in boundary planes.
! otherwise particle creation is forbidden in the boundary plane.
#define PARTICLECREATION_IN_BOUNDARY

! Gas does not have gravity
! #define FLUID_GRAVITY_OFF

! Linear momentum is conserved between gas and particle.
! The gravity of gas acting on particles are evaluated as back reaction of gravity from particle to gas.
!#define MOMENTUM_CONSERVATION

! density foor in Bondi accretion
#define RHO_FLOOR 1.d0

module sinkParticle
  use overBlockCoordinates
  use grid
  use string, only : CHARLEN
  use unit
  implicit none
  private
  type t_spParticle
     integer :: pid                  ! id of particle
     real(kind=DBL_KIND) :: mass     ! mass of particle
     real(kind=DBL_KIND) :: dmass    ! delta mass (due to accretion)
     real(kind=DBL_KIND) :: r(MX:MZ) ! location of particle
     real(kind=DBL_KIND) :: dr(MX:MZ) ! delta r (due to gravity)
     real(kind=DBL_KIND) :: v(MX:MZ) ! velocity of particle
     real(kind=DBL_KIND) :: dv(MX:MZ) ! delta velocity (due to gravity)
     real(kind=DBL_KIND) :: dp(MX:MZ) ! delta momentum (due to accretion)
     real(kind=DBL_KIND) :: ds(MX:MZ) ! delta specific spin angular momentum (due to accreion)
     real(kind=DBL_KIND) :: s(MX:MZ)  ! specific spin angular momentum
     !----- KS ADDED ----!
     real(kind=DBL_KIND) :: t_prev         ! previous time (starting point of averaging for disk accretion)
     real(kind=DBL_KIND) :: dm_disk        ! accreted mass after t_prev
     real(kind=DBL_KIND) :: dJ_disk(MX:MZ) ! accreted ang. mom. after t_prev
     real(kind=DBL_KIND) :: mdot_disk      ! acc. rate through disk
     real(kind=DBL_KIND) :: J_disk(MX:MZ)  ! ang. mom. accreted onto disk     
     !----- KS ADDED ----!     
     type(t_spParticle),pointer :: next => null()
  end type t_spParticle
  ! Type for link list of position of new particles
  type t_spPos
     real(DBL_KIND),dimension(MX:MZ) :: r
     type(t_spPos),pointer :: next => null()
  end type t_spPos
  integer,save :: Nparticle = 0
  integer,save :: PidMax = -1
  integer,save :: sp_Level = Undefi
  logical,save :: Initialized = .false.
  type(t_spParticle),save,pointer :: Particle => null() ! link list
  real(kind=DBL_KIND),save :: sp_Rhocr, sp_RhocrCreate ! KS MODIFIED
#ifdef MP
  real(kind=DBL_KIND),save :: sp_Cs ! KS MODIFIED
!  real(kind=DBL_KIND),parameter :: sp_Cs = 1.d0 ! sound speed of sink particle creation. used for estimation of sink radius through Jeans length.  
#endif !MP
  real(kind=DBL_KIND),save :: sp_SinkRadius ! radius of sink particle, defined at sp_init
  real(kind=DBL_KIND),save :: sp_SofteningRadius ! softening radius, defined at sp_init
  real(kind=DBL_KIND),save :: sp_MaskRadius      ! masking of particle creation, defined at sp_init
  real(kind=DBL_KIND),parameter :: sp_f_rAcc = 4.d0
  real(kind=DBL_KIND),save :: sp_DtimeCFL
  real(kind=DBL_KIND),parameter :: sp_CGS = 1.d-1
  real(kind=DBL_KIND),parameter :: sp_CVS = 5.d-1 ! limit of timestep
  character(len=CHARLEN),parameter :: FILENAME_LOG='logSinkparticle'
  character(len=CHARLEN),parameter :: FILENAME='sinkparticle.d'
#ifdef ACCRETION_BONDI
  ! for bondi accretion
  integer,save :: sp_BONDI_NMAX
  real(KIND=DBL_KIND),save,dimension(:),allocatable :: sp_BONDI_RHO
  real(kind=DBL_KIND),save :: sp_BONDI_DX
#endif !ACCRETION_BONDI

  ! accretion time scale = t_visc at outer boundary of disk (KS ADDED, cf. Hosokawa+16)
  real(kind=DBL_KIND),parameter :: dt_acc = 3d2 ! in yr

  !
  public :: sp_update, sp_write, sp_read, sp_writeLog, sp_restrictCFL, sp_refineCond, sp_getLevel, &
       sp_sinkdata2array, sp_getNparticle, sp_getRhocr, sp_newParticle, sp_gravityOfParticle, &
       sp_advector, sp_getSinkRadius, sp_is_inside_sink, sp_refineCond_KS !KS MODIFIED
contains
  !-------------------------------------------------------------------------
  ! Update sink particles.
  ! Exported subroutine.
  !-------------------------------------------------------------------------
  subroutine sp_update
    use parameter

    call sp_init

#ifdef CREATION_ON
    call sp_create
#endif !CREATION_ON
    if (Nparticle == 0) return
    if (LevelMax /= sp_Level) then
       print *, '*** waring : LevelMax /= sp_Level', LevelMax, sp_Level
       call flush(6)
    endif
#ifdef ACCRETION_ON
    call sp_accretion_init
#ifdef ACCRETION_BONDI
    call sp_accretionBondi( Dtime(sp_Level) )
#endif !ACCRETION_BONDI
#ifdef ACCRETION_RHOCR
   call sp_accretion
#endif !ACCRETION_RHOCR
    call sp_accretion_comm
    call sp_accretionUpdate
#endif !ACCRETION_ON
#ifdef ADVECTION_ON
    call sp_advector( Dtime(sp_Level) ) ! also make CFL condition for restricted timestep of HD.
#endif !ADVECTION_ON
#ifdef MERGER_ON
    call sp_merge
#endif !MERGER_ON

    !----- KS ADDED ----!
    call sp_update_subdisk
    !----- KS ADDED ----!

    call sp_boundary_eject
    call sp_writeLog
!!$    call sp_testIO                 ! debug

  end subroutine sp_update
#ifdef CREATION_ON
  !-------------------------------------------------------------------------
  ! Create Particle.
  ! This routine determs *only* the locations of new particles.
  !-------------------------------------------------------------------------
  subroutine sp_create
    use mpilib
    type(t_spPos),pointer :: list_newPos
    if (LevelMax < sp_Level) return
    call sp_getListNewPos( list_newPos )
#ifndef PARTICLECREATION_IN_BOUNDARY
    call sp_createConditionBoundary( list_newPos )
#endif !PARTICLECREATION_IN_BOUNDARY
    call sp_createConditionEnergies( list_newPos )
    call sp_createGather( list_newPos )
    if (get_myrank() == PRIMARY_RANK .and. associated(list_newPos)) print *, 'new particles created'
    call sp_makeParticle( list_newPos )
    call sp_deallocate_list_newPos( list_newPos )
  end subroutine sp_create
  !-------------------------------------------------------------------------
  ! Get list of new particles  (MPI global)
  ! OUTPUT:
  !   list_newPos = a link list for positions of new particles.
  !-------------------------------------------------------------------------
  subroutine sp_getListNewPos(list_newPos)
    type(t_spPos),pointer :: list_newPos       ! MPI global
    type(t_spPos),pointer :: list_newPos_Local ! MPI local
    ! make list locally
    call sp_getListNewPos_Local(list_newPos_Local)
    ! allreduce a list over the every nodes to make list_newPos
    call sp_listNewPos_Allreduce(list_newPos, list_newPos_Local)
    call sp_deallocate_list_newPos(list_newPos_Local)
  end subroutine sp_getListNewPos
  !-------------------------------------------------------------------------
  ! Get list of new particles (MPI local)
  ! OUTPUT:
  !   list_newPos_Local = link list of newPos made in each node individually.
  !-------------------------------------------------------------------------
  subroutine sp_getListNewPos_Local( list_newPos )
    use mpilib
    use grid_boundary
    type(t_spPos),pointer :: list_newPos
    integer :: level, gid, n
    myrank = get_myrank()
    ! make serial list
    nullify(list_newPos)
    level = sp_Level
    call boundary_grid( level, COMPLETE )
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       call sp_findNewParticle_local(list_newPos, gid)
    end do
  end subroutine sp_getListNewPos_Local
  !-------------------------------------------------------------------------
  ! Evaluate condition of particle creation and returns new position of particle (MPI local).
  ! Returned value:
  !   true if a particle is created. elsewhere false.
  ! INPUT:
  !   gid = grid id
  ! OUTPUT:
  !   list_newPos = a list of positions of particles (type: newPos)
  ! Conditions:
  !   1. gravitational potential is local minimum.
  !   2. density is higher than sp_Rhocr.
  !   3. divergence velocity is compressed
  !-------------------------------------------------------------------------
  subroutine sp_findNewParticle_local(list_newPos, gid)
    use mpilib
    type(t_spPos),pointer :: list_newPos
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: psi, rho, vx, vy, vz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    logical,dimension(ARRAYSIZE3(GridMask)) :: maskpos
    type(t_spParticle),pointer :: ptr
    type(t_spPos),pointer :: newPos, ptrn
    real(kind=DBL_KIND),dimension(MX:MZ) :: r
    integer :: i, j, k
    rho => get_Ucomp(MRHO, gid)
    if (maxval(rho, mask=GridMask) < sp_Rhocr) return
    if ( get_level(gid) < sp_Level ) return
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    psi => get_Ucomp(MPSI, gid)
    vx => get_Ucomp(MVX, gid)
    vy => get_Ucomp(MVY, gid)
    vz => get_Ucomp(MVZ, gid)
    ! ------------------
    ! initialize maskpos
    ! ------------------
    maskpos = GridMask
    ptr => Particle
    do
       if (.not. associated(ptr)) exit

       if ( &              ! out of region of gid
            ptr%r(MX) < x(Imin) - sp_MaskRadius .or. ptr%r(MX) > x(Imax) + sp_MaskRadius .or. &
            ptr%r(MY) < y(Jmin) - sp_MaskRadius .or. ptr%r(MY) > y(Jmax) + sp_MaskRadius .or. &
            ptr%r(MZ) < z(Kmin) - sp_MaskRadius .or. ptr%r(MZ) > z(Kmax) + sp_MaskRadius) then
          ptr => ptr%next
          cycle
       end if

       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                r = (/x(i), y(j), z(k)/)
                if (sum((r - ptr%r)**2) <= sp_MaskRadius**2) maskpos(i,j,k) = .false.
             end do
          end do
       end do
       ptr => ptr%next
    end do
    ! -------------------
    ! find local minimum
    ! -------------------
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             if (.not. maskpos(i,j,k)) cycle
             if ( rho(i,j,k) < sp_RhocrCreate ) cycle ! exceed critical density?
             if ( psi(i,j,k) > minval(psi(i-1:i+1,j-1:j+1,k-1:k+1))) cycle ! local minimum?
             if ( .not. sp_createConditionDivergenceV() ) cycle
             ! ----------------
             ! assign to newPos
             ! ----------------
             allocate( newPos )
             newPos%r = (/x(i), y(j), z(k)/)
             ptrn => list_newPos
             list_newPos => newPos
             list_newPos%next => ptrn
             maskpos(i,j,k) = .false.
             print *, 'allocate new particle', gid; call flush(6)
          end do
       end do
    end do
    contains
      !-------------------------------------------------------------------------
      ! check divergence of velocity tensor for sp creation
      !-------------------------------------------------------------------------
      function sp_createConditionDivergenceV() result(bool)
        logical :: bool
        real(kind=DBL_KIND) :: vm(MX:MZ, MX:MZ), eigenvalues(MX:MZ)
        bool = .false.
        ! ordinal divergence v
        if (vx(i+1,j,k)-vx(i-1,j,k)+vy(i,j+1,k)-vy(i,j-1,k)+vz(i,j,k+1)-vz(i,j,k-1) > 0) then
           print *, '* particle creation is rejected by divergence v'
           return
        end if
        ! principal axis of divergence v
        vm(MX,MX) = vx(i+1,j,k)-vx(i-1,j,k)
        vm(MX,MY) = (vy(i+1,j,k)-vy(i-1,j,k)+vx(i,j+1,k)-vx(i,j-1,k))*0.5d0
        vm(MX,MZ) = (vz(i+1,j,k)-vz(i-1,j,k)+vx(i,j,k+1)-vx(i,j,k-1))*0.5d0
        vm(MY,MY) = vy(i,j+1,k)-vy(i,j-1,k)
        vm(MY,MZ) = (vz(i,j+1,k)-vz(i,j-1,k)+vy(i,j,k+1)-vy(i,j,k-1))*0.5d0
        vm(MZ,MZ) = vz(i,j,k+1)-vz(i,j,k-1)
        call DSYEVC3(vm, eigenvalues) ! get eigenvalues
        if (eigenvalues(MX) > 0 .or. eigenvalues(MY) > 0 .or. eigenvalues(MZ) > 0) then
           print *, '* particle creation is rejected by eigenvalues of v-tensor', eigenvalues
           return
        endif
        bool = .true.
      end function sp_createConditionDivergenceV
  end subroutine sp_findNewParticle_local
  !-------------------------------------------------------------------------
  ! create condition for Jeans condition
  !-------------------------------------------------------------------------
  subroutine sp_createConditionEnergies(list_newPos)
    use eos
    use parameter, only : Pi8i
    use modelParameter, only : MP_Gconst
    use mpilib
    type(t_spPos),pointer :: list_newPos
    !
    type(t_spPos),pointer :: ptr, prev, dell
    integer :: level, gid, rank, i,j,k, ig, jg, kg
    integer,dimension(MX:MZ) :: ijkgL, ijkgR
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz
#ifdef MP
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: p
#else !MP
    real(kind=DBL_KIND) :: p
#endif !MP
    real(kind=DBL_KIND) :: dv, Eth, Egrav, Ekin, Emag, Etot, mass, vxg, vyg, vzg, Ekinr
    real(kind=DBL_KIND),dimension(6) :: buf, bufr
    !
    real(kind=DBL_KIND),parameter :: sp_Gamma = 1.1d0
#ifndef MP
    ! Estimate of kappa for barotropic gas
    real(kind=DBL_KIND) :: sp_Kappa
    sp_Kappa = Kappa * sp_Rhocr ** (Gamma - sp_Gamma)
#endif !MP
    level = sp_Level
    dv = get_dv(level)

    ptr => list_newPos
!     prev => null()
    nullify(prev)
    do
       if (.not. associated(ptr)) exit
       ! evaluate blocks
       posL(:) = ptr%r(:) - sp_SinkRadius
       posR(:) = ptr%r(:) + sp_SinkRadius
#ifdef PARTICLECREATION_IN_BOUNDARY
       call sp_restrict_within_boundary(posL, posR)
#endif !PARTICLECREATION_IN_BOUNDARY
       Eth = 0.d0
       Egrav = 0.d0
       Ekin = 0.d0
       Emag = 0.d0
       mass = 0.d0
       vxg = 0.d0
       vyg = 0.d0
       vzg = 0.d0

       call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
       call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
       do kg = ijkgL(MZ), ijkgR(MZ)
          do jg = ijkgL(MY), ijkgR(MY)
             do ig = ijkgL(MX), ijkgR(MX)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(MRHO, gid)
                vx => get_Ucomp(MVX, gid)
                vy => get_Ucomp(MVY, gid)
                vz => get_Ucomp(MVZ, gid)
#ifdef MP
                p => get_Ucomp(MP, gid)
#endif !MP
#ifdef MHD
                bx => get_Ucomp(MBX, gid)
                by => get_Ucomp(MBY, gid)
                bz => get_Ucomp(MBZ, gid)
#endif !MHD
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         if ( (x(i)-ptr%r(MX))**2 + (y(j)-ptr%r(MY))**2 + (z(k)-ptr%r(MZ))**2 <= sp_SinkRadius**2 ) then
!!$                            Eth = Eth + rho(i,j,k) * dv
!!$                            call get_p(p, rho(i,j,k))
#ifdef MP
                            Eth = Eth + p(i,j,k)*dv ! Eth = (3/2) p dV. for a factor (3/2) see below
#else !MP
                            p = sp_Kappa * rho(i,j,k) ** sp_Gamma
                            Eth = Eth + p * dv
#endif !MP
#ifdef MHD
                            Emag = Emag + (bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2) * Pi8i *dv
#endif !MHD
                            mass = mass + rho(i,j,k)*dv
                            vxg = vxg + rho(i,j,k)*vx(i,j,k)*dv
                            vyg = vyg + rho(i,j,k)*vy(i,j,k)*dv
                            vzg = vzg + rho(i,j,k)*vz(i,j,k)*dv
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
       buf(:) = (/Eth, Emag, mass, vxg, vyg, vzg/)
!        call mpi_allreduce( MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce( buf, bufr, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       buf = bufr
       !Eth = buf(1) * 1.5d0     ! (3/2) p dV
       Eth = buf(1)      ! sinkができずに周囲にモノが溜まり過ぎるように見えたので、条件を変えてやってみる (KS TODO)
       Emag = buf(2)
       mass = buf(3)
       vxg = buf(4)
       vyg = buf(5)
       vzg = buf(6)
       Egrav = - 3 * MP_Gconst *mass**2 / (5.d0 * sp_SinkRadius)
       vxg = vxg/mass
       vyg = vyg/mass
       vzg = vzg/mass
       do kg = ijkgL(MZ), ijkgR(MZ)
          do jg = ijkgL(MY), ijkgR(MY)
             do ig = ijkgL(MX), ijkgR(MX)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(MRHO, gid)
                vx => get_Ucomp(MVX, gid)
                vy => get_Ucomp(MVY, gid)
                vz => get_Ucomp(MVZ, gid)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         if ( (x(i)-ptr%r(MX))**2 + (y(j)-ptr%r(MY))**2 + (z(k)-ptr%r(MZ))**2 <= sp_SinkRadius**2 ) then
                            Ekin = Ekin + 0.5d0*rho(i,j,k)*( &
                                 (vx(i,j,k)-vxg)**2 + (vy(i,j,k)-vyg)**2 + (vz(i,j,k)-vzg)**2)*dv
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
!        call mpi_allreduce( MPI_IN_PLACE, Ekin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce( Ekin, Ekinr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       Ekin = Ekinr

!--- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING ---!
       Etot = Eth + Ekin + Emag + Egrav ! inward velocityの成分を取り除いてもいいかも (KS TODO)
       !Etot = Eth + Emag + Egrav      !KS MODIFEID (TENTATIVE?)
!--- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING ---!       
       if (Etot > 0) then
          if (get_myrank() == PRIMARY_RANK) print *, '* particle creation is rejected by energy criterion ', Etot, Eth, Ekin, Emag, Egrav
          if ( associated(ptr, list_newPos) ) then ! delete first element
             list_newPos => ptr%next
          else                  ! delete second or later elements
             prev%next => ptr%next
          end if
          dell => ptr
          ptr => ptr%next
          deallocate( dell )
          cycle
       end if
       prev => ptr
       ptr => ptr%next
    end do
  end subroutine sp_createConditionEnergies
  !-------------------------------------------------------------------------
  ! create condition for boundaries.
  ! particle creation over boundary surface is rejected.
  !
  ! INPUT :
  !   list_newPos = a link list for positions of new particles.
  !-------------------------------------------------------------------------
  subroutine sp_createConditionBoundary( list_newPos )
    type(t_spPos),pointer :: list_newPos
    type(t_spPos),pointer :: ptr, prev, dell
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    type(t_obRectPhys) :: compDomain, particleDomain
    type(t_obPointPhys) :: pL, pR
    call ob_computationBoxOfRectPhys( compDomain )
    ptr => list_newPos
    nullify(prev)
    do
       if (.not. associated(ptr)) exit
       posL(:) = ptr%r(:) - sp_SinkRadius
       posR(:) = ptr%r(:) + sp_SinkRadius
       call ob_assignCoordPhysToPointPhys(posL, pL)
       call ob_assignCoordPhysToPointPhys(posR, pR)
       call ob_assignPointPhysToRectPhys(pL, particleDomain, 'L')
       call ob_assignPointPhysToRectPhys(pR, particleDomain, 'R')
       if ( ob_rectPhysComp( particleDomain, compDomain ) /= OB_RECT_INNER ) then
          if (get_myrank() == PRIMARY_RANK) print *, '* particle creation is rejected by boundary.'
          if ( associated(ptr, list_newPos) ) then ! delete first element
             list_newPos => ptr%next
          else                  ! delete second or later elements
             prev%next => ptr%next
          end if
          dell => ptr
          ptr => ptr%next
          deallocate( dell )
          cycle
       end if
       prev => ptr
       ptr => ptr%next
    end do
  end subroutine sp_createConditionBoundary
  !-------------------------------------------------------------------------
  ! deallocate list_newPos
  ! INPUT :
  !   list_newPos = link list of t_newPos
  !-------------------------------------------------------------------------
  subroutine sp_deallocate_list_newPos( list_newPos )
    type(t_spPos),pointer :: list_newPos, next, ptr
    ptr => list_newPos
    do
       if (.not. associated( ptr ) ) exit
       next => ptr%next
       deallocate( ptr )
       ptr => next
    end do
    nullify( list_newPos )
  end subroutine sp_deallocate_list_newPos
  !-------------------------------------------------------------------------
  ! get list of new particles (MPI local)
  ! INPUT:
  !   list_newPos_Local = link list of newPos made in each node individually.
  ! OUTPU:
  !   list_newPos = link list of rallreduced list_newPos_Local
  !-------------------------------------------------------------------------
  subroutine sp_listNewPos_Allreduce(list_newPos, list_newPos_Local)
    use mpilib
    type(t_spPos),pointer :: list_newPos, list_newPos_Local
    integer :: n, rank, n_list_newPos, n_list_newPos_Local
    real(DBL_KIND),dimension(:,:),allocatable :: buf
    type(t_spPos),pointer :: elem, ptr, ptrg
    ! count n_list_newPos_Local
    n_list_newPos_Local = 0
    ptr => list_newPos_Local
    do
       if (.not. associated(ptr)) exit
       n_list_newPos_Local = n_list_newPos_Local + 1
       ptr => ptr%next
    end do
    ! allreduce
    myrank = get_myrank()
    nullify(list_newPos)
    do rank = 0, NPE-1
       ! bcast n_list_newPos_Local to n_list
       if (rank == myrank) n_list_newPos = n_list_newPos_Local
       call mpi_bcast(n_list_newPos, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (n_list_newPos == 0) cycle
       ! push to buffer
       allocate(buf(MX:MZ, n_list_newPos))
       if (rank == myrank) then
          ptr => list_newPos_Local
          do n = 1, n_list_newPos
             buf(:,n) = ptr%r
             ptr => ptr%next
          end do
       end if
       call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
       ! push to list_newPos
       do n = 1, n_list_newPos
          allocate( elem )
          elem%r(:) = buf(:,n)
          ptrg => list_newPos
          list_newPos => elem
          list_newPos%next => ptrg
       end do
       deallocate( buf )
    end do
  end subroutine sp_listNewPos_Allreduce
#endif !CREATION_ON
  !-------------------------------------------------------------------------
  ! gather particles which locate close each other
  !-------------------------------------------------------------------------
  subroutine sp_createGather(list_newPos)
    use mpilib
    type(t_spPos),pointer :: list_newPos, ptra, ptrb, prev
    type t_spw              ! list of weght
       integer :: w         ! weight
       type(t_spw),pointer :: next => null()
    end type t_spw
    type(t_spw),pointer :: weight, list_weight, ptrwa, ptrwb, prevw
    integer :: pida, pidb       ! particle id
    integer,dimension(2) :: pair
    real(kind=DBL_KIND) :: dmin2, dist2
    if (.not. associated(list_newPos)) return
    ! -----------------
    ! make list_weight
    ! -----------------
    nullify(list_weight)
    ptra => list_newPos
    do
       if (.not. associated( ptra ) ) exit
       allocate(weight)
       weight%w = 1
       if (associated(list_weight)) then
          ptrwa%next => weight
          ptrwa => ptrwa%next
       else
          list_weight => weight        !1st element
          ptrwa => list_weight
       endif
       ptra => ptra%next
    end do
    do
       ! ---------------------
       ! search a closest pair
       ! ---------------------
       pair(:) = Undefi
       dmin2 = huge(dmin2)
       ptra => list_newPos
       pida = 0
       do                          ! loop of ptra
          if (.not. associated( ptra ) ) exit

          ptrb => list_newPos
          pidb = 0
          do                       ! loop of ptrb
             if (.not. associated( ptrb ) ) exit
             if (pidb >= pida) exit ! compare with pidb (< pida)
             ! if (get_myrank() == PRIMARY_RANK) print *, 'pair', pida, pidb
             dist2 = sum((ptra%r - ptrb%r)**2)
             if (dist2 < dmin2) then
                dmin2 = dist2
                pair(:) = (/ min(pida, pidb), max(pida, pidb) /)
             endif
             pidb = pidb + 1
             ptrb => ptrb%next
          end do

          pida = pida + 1
          ptra => ptra%next

       end do
       ! ----------------------
       ! return from subroutine
       ! ----------------------
       if (dmin2 > sp_MaskRadius**2) then
          call sp_createGather_dealloc_listweight(list_weight)
          return
       end if
       ! if (get_myrank() == PRIMARY_RANK) print *, 'merge pair, dmin2', pair, dmin2
       ! ----------------------------------------------------------
       ! marge pair(1) and pair(2) into pair(1), and remove pair(2)
       ! ----------------------------------------------------------
       nullify(prev)
       nullify(prevw)
       ptra => list_newPos
       pida = 0
       ptrwa => list_weight
       do                          ! loop of ptra
          if (.not. associated( ptra ) ) exit
          if (pida == pair(1)) then
             ptrb => list_newPos
             pidb = 0
             ptrwb => list_weight
             do                       ! loop of ptrb
                if (.not. associated( ptrb ) ) exit
                if (pidb == pair(2)) then
                   ! if (get_myrank() == PRIMARY_RANK) print *, 'pair ave', pida, pidb, ptrwa%w, ptrwb%w
                   ptra%r = (ptra%r * ptrwa%w + ptrb%r * ptrwb%w)/(ptrwa%w + ptrwb%w) ! average with weight
                   ptrwa%w = ptrwa%w + ptrwb%w
                   if (get_myrank() == PRIMARY_RANK) print "('merge particles: pos, weight ',3E12.5, I2)", ptra%r, ptrwa%w
                   ! delete ptrb (ptrb is not 1st element. ptra can be 1st element.)
                   prev%next => ptrb%next
                   prevw%next => ptrwb%next
                   deallocate( ptrb )
                   deallocate( ptrwb )
                   exit
                end if
                prev => ptrb
                pidb = pidb + 1
                ptrb => ptrb%next
                prevw => ptrwb
                ptrwb => ptrwb%next
             end do
          end if
          pida = pida + 1
          ptra => ptra%next
          ptrwa => ptrwa%next
       end do
    end do
    call sp_createGather_dealloc_listweight(list_weight)
    contains
      !-------------------------------------------------------------------------
      ! deallocate list of weight
      !-------------------------------------------------------------------------
      subroutine sp_createGather_dealloc_listweight(list_weight)
        type(t_spw),pointer :: list_weight, next, ptr
        ptr => list_weight
        do
           if (.not. associated(ptr)) exit
           next => ptr%next
           deallocate( ptr )
           ptr => next
        end do
      end subroutine sp_createGather_dealloc_listweight

  end subroutine sp_createGather
  !-------------------------------------------------------------------------
  ! make particle because of list_newPos
  ! INPUT:
  !   list_newPos = link list of t_newPos (MPI global)
  !-------------------------------------------------------------------------
  subroutine sp_makeParticle(list_newPos)
    type(t_spPos),pointer :: list_newPos, ptr
    ptr => list_newPos
    do
       if (.not. associated( ptr ) ) exit
       !call sp_newParticle( 0.d0, ptr%r, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/) )
       call sp_newParticle( 0.d0, ptr%r, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), t_prev=Time(sp_Level) ) !KS MODIFIED
       ptr => ptr%next
    end do
  end subroutine sp_makeParticle
#ifdef ADVECTION_ON
  !-------------------------------------------------------------------------
  ! advector of sink particles
  ! dtGlobal = timestep to be advect.
  !-------------------------------------------------------------------------
  subroutine sp_advector(dtGlobal)
    use mpilib
    real(kind=DBL_KIND),intent(IN) :: dtGlobal
    real(kind=DBL_KIND) :: dtLocal, t
    real(kind=DBL_KIND),dimension(:,:),allocatable :: gravFluid
    logical :: bool_last, bool_subcycle_on
    t = 0.d0                    ! local closk for subcycle: 0 <= t <= dtGlobal
    bool_last = .false.
    bool_subcycle_on = .false.
    allocate(gravFluid(MX:MZ, Nparticle))
    do
       call sp_timestepOfSubcycle(dtLocal, gravFluid, bool_subcycle_on)
       if ( dtLocal < dtGlobal ) then
          if (.not. bool_subcycle_on) then ! Initialize for subcycling.
             bool_subcycle_on = .true.
             call sp_gravityF2P_all(gravFluid, dtGlobal) ! dt is dtGlobal. 
             ! Gravity of gas, gravFluid, is obtained here and it is used through the subcycles.
             ! Particles affects gas during dt = dtGlobal.
             ! When subcycling, here after, particles never affect gas.
             if (get_myrank() == PRIMARY_RANK) print *, 'Sink particle: subcycle ON'
          endif
       endif
       if ( t + dtLocal >= dtGlobal ) then
          dtLocal = dtGlobal - t
          bool_last = .true.
       end if
       t = t + dtLocal
       call sp_advectorSubcycle(dtLocal, gravFluid, bool_subcycle_on)
       call sp_advectionUpdate
       call sp_boundary_eject
       if (get_myrank() == PRIMARY_RANK) print '("* advector:",1P3E12.5)', dtLocal, dtGlobal, t
       if (bool_last) exit
    end do
    deallocate(gravFluid)
  end subroutine sp_advector
  !-------------------------------------------------------------------------
  ! get fluid gravity for all particles
  !-------------------------------------------------------------------------
  subroutine sp_gravityF2P_all(gravFluid, dt)
    real(kind=DBL_KIND),dimension(MX:,:),intent(OUT) :: gravFluid
    real(kind=DBL_KIND),intent(IN) :: dt
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:MZ) :: gravF
    integer :: np
    np = 0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       np = np + 1
       call sp_gravityFluid2particle(ptr, gravF, dt) 
       gravFluid(:,np) = gravF
       ptr => ptr%next
    end do
  end subroutine sp_gravityF2P_all
  !-------------------------------------------------------------------------
  ! advect sink particle
  ! dt = delta t
  !-------------------------------------------------------------------------
  subroutine sp_advectorSubcycle(dt, gravFluid, bool_subcycle_on)
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(MX:,:),intent(IN) :: gravFluid
    logical,intent(IN) :: bool_subcycle_on
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:MZ) :: gravP, gravF
    real(kind=DBL_KIND),dimension(:,:),allocatable :: grav
    integer,dimension(MX:MZ) :: mlist
    integer :: np
    myrank = get_myrank()
    mlist = (/MGX, MGY, MGZ/)
    sp_DtimeCFL = HUGE(1.d0)
    allocate(grav(MX:MZ, Nparticle))
    ! -------------------
    ! predictor step (dr)
    ! -------------------
    np = 0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       np = np + 1
       call sp_gravityP2P(ptr, gravP) ! gravity of particle
       if ( bool_subcycle_on ) then
          gravF = gravFluid(:,np)
       else
          call sp_gravityFluid2particle(ptr, gravF, dt*0.5d0)
       endif
       grav(:, np) = gravP + gravF
       ptr%dr(:) = ptr%v(:)*dt + grav(:, np) *dt**2 * 0.5d0
       ! CFL condition
       sp_DtimeCFL = min( sp_DtimeCFL, &
            minval( CellWidth(:,sp_Level) /abs(ptr%v(:) + grav(:, np)*dt)) &
            )
       ptr => ptr%next
    end do
    ! -------------------
    ! update r temporally
    ! -------------------
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ptr%r = ptr%r + ptr%dr
       ptr => ptr%next
    end do
    ! -------------------
    ! corrector step (dv)
    ! -------------------
    np = 0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       np = np + 1
       call sp_gravityP2P(ptr, gravP) ! gravity of particle
       if ( bool_subcycle_on ) then
          gravF = gravFluid(:,np)
       else
          call sp_gravityFluid2particle(ptr, gravF, dt*0.5d0)
       endif
       grav(:, np) = (grav(:, np) + gravP + gravF)*0.5d0
       ptr%dv(:) = grav(:, np) *dt
       ! CFL condition
       sp_DtimeCFL = min( sp_DtimeCFL, &
            minval( CellWidth(:,sp_Level) /abs(ptr%v(:) + ptr%dv(:))) &
            )
       ptr => ptr%next
    end do
    ! -------------------
    ! restore r
    ! -------------------
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ptr%r = ptr%r - ptr%dr
       ptr => ptr%next
    end do

    deallocate( grav )
  end subroutine sp_advectorSubcycle
  !-------------------------------------------------------------------------
  ! Get time step of subcycle for advector of particles.
  ! OUTPUT:
  !    dtsub ... timestep of subcycle
  ! INPUT:
  !    gravFluid ... gravity of gas, which should be given if bool_subcycle_on is true.
  !    bool_subcycle_on (optional) 
  !-------------------------------------------------------------------------
  subroutine sp_timestepOfSubcycle(dtsub, gravFluid, bool_subcycle_on)
    real(kind=DBL_KIND),intent(OUT) :: dtsub
    real(kind=DBL_KIND),dimension(MX:,:),intent(IN) :: gravFluid
    logical,intent(IN) :: bool_subcycle_on
    type(t_spParticle),pointer :: ptri, ptrj
    real(kind=DBL_KIND) :: dr2min, grav, h
    integer,dimension(MX:MZ) :: mlist
    real(kind=DBL_KIND),dimension(MX:MZ) :: gravP, gravF
    integer :: np
    np = 0
    mlist = (/MGX, MGY, MGZ/)
    dtsub = HUGE(dtsub)
    ptri => Particle
    do                          ! for ptri
       if (.not. associated(ptri)) exit
       np = np + 1
       dr2min = HUGE(dr2min)
       ptrj => Particle
       do
          if (.not. associated(ptrj)) exit
          if (ptri%pid == ptrj%pid) then
             ptrj => ptrj%next
             cycle
          end if
          dr2min = min(dr2min, sum((ptri%r - ptrj%r)**2))
          ptrj => ptrj%next
       end do
       call sp_gravityP2P(ptri, gravP)
       if (bool_subcycle_on) then
          gravF = gravFluid(:,np)
       else
          call sp_gravityOfFluid(ptri, gravF)
       endif
       grav =  sqrt(sum((gravP + gravF)**2))
       h = minval(CellWidth(:, sp_Level))
       dtsub = min(dtsub, sqrt(min(sqrt(dr2min), h)/grav) )
       ptri => ptri%next
    end do
    dtsub = sp_CGS * dtsub
  end subroutine sp_timestepOfSubcycle
#endif !ADVECTION_ON
  !-------------------------------------------------------------------------
  ! gravity from particles to a particle
  !-------------------------------------------------------------------------
  subroutine sp_gravityP2P(ptri, grav)
    type(t_spParticle),pointer :: ptri
    real(kind=DBL_KIND),dimension(MX:MZ),intent(OUT) :: grav
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND) :: gx, gy, gz
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr
    ptr => Particle
    grav = 0
    do
       if (.not. associated(ptr)) exit
       if (ptr%pid == ptri%pid) then
          ptr => ptr%next
          cycle
       end if
       dr = ptri%r - ptr%r
       call sp_gravityOfParticle( dr, ptr%mass, gx, gy, gz, no_softening=.true.)
       grav(MX) = grav(MX) + gx
       grav(MY) = grav(MY) + gy
       grav(MZ) = grav(MZ) + gz
       ptr => ptr%next
    end do
  end subroutine sp_gravityP2P
#ifdef ACCRETION_ON
  !-------------------------------------------------------------------------
  ! Initialize accretion
  !-------------------------------------------------------------------------
  subroutine sp_accretion_init
    type(t_spParticle),pointer :: ptr
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ptr%dmass = 0.d0
       ptr%dp(:) = 0.d0
       ptr%ds(:) = 0.d0
       ptr => ptr%next
    end do
  end subroutine sp_accretion_init
  !-------------------------------------------------------------------------
  ! all reduce accretion
  !-------------------------------------------------------------------------
  subroutine sp_accretion_comm
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(:,:),allocatable :: buf, bufr
    integer :: n
    allocate(buf(7, Nparticle), bufr(7, Nparticle))
    ptr => Particle
    do n = 1, Nparticle
       if (.not. associated(ptr)) then
          print *, '*** error in sp_accretion: ptr is no associated.'
       end if
       buf(1, n) = ptr%dmass
       buf(2:4, n) = ptr%dp
       buf(5:7, n) = ptr%ds
       ptr => ptr%next
    end do
    call mpi_allreduce( buf, bufr , size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    ptr => Particle
    do n = 1, Nparticle
       if (.not. associated(ptr)) then
          print *, '*** error in sp_accretion: ptr is no associated.'
       end if
       ptr%dmass = buf(1, n)
       ptr%dp = buf(2:4, n)
       ptr%ds = buf(5:7, n)
!!$       if (get_myrank() == PRIMARY_RANK) print *, ptr%pid, 'dmass =', ptr%dmass
       ptr => ptr%next
    end do
    deallocate(buf, bufr)
  end subroutine sp_accretion_comm
#endif !ACCRETION_ON
#ifdef ACCRETION_RHOCR
  !-------------------------------------------------------------------------
  ! accretion to sink particle
  !-------------------------------------------------------------------------
  subroutine sp_accretion
    use mpilib
    type(t_spParticle),pointer :: ptr
    integer :: level, gid, rank, i,j,k, ig, jg, kg
    integer,dimension(MX:MZ) :: ijkgL, ijkgR
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, pr !KS ADDED
    real(kind=DBL_KIND) :: dv, rhoacc
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr, du
    ptr => Particle
    level = LevelMax
    do
       if (.not. associated(ptr)) exit
       ! evaluate blocks
       posL(:) = ptr%r(:) - sp_SinkRadius
       posR(:) = ptr%r(:) + sp_SinkRadius
#ifdef PARTICLECREATION_IN_BOUNDARY
       call sp_restrict_within_boundary(posL, posR)
#endif !PARTICLECREATION_IN_BOUNDARY
       call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
       call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
       do kg = ijkgL(MZ), ijkgR(MZ)
          do jg = ijkgL(MY), ijkgR(MY)
             do ig = ijkgL(MX), ijkgR(MX)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(MRHO, gid)
                pr =>  get_Ucomp(MP, gid) !KS ADDED
                vx => get_Ucomp(MVX,gid)
                vy => get_Ucomp(MVY,gid)
                vz => get_Ucomp(MVZ,gid)
                dv = get_dv(level)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         if ( (x(i)-ptr%r(MX))**2 + (y(j)-ptr%r(MY))**2 + (z(k)-ptr%r(MZ))**2 < sp_SinkRadius**2 ) then
                            rhoacc = max(rho(i,j,k) - sp_Rhocr, 0.d0) ! mass of accretion
                            ! particle gets mass and momentum
                            ptr%dmass = ptr%dmass +  rhoacc*dv
                            ptr%dp(MX) = ptr%dp(MX) + vx(i,j,k)*rhoacc*dv
                            ptr%dp(MY) = ptr%dp(MY) + vy(i,j,k)*rhoacc*dv
                            ptr%dp(MZ) = ptr%dp(MZ) + vz(i,j,k)*rhoacc*dv
                            dr = (/x(i), y(j), z(k)/) - ptr%r
                            du = (/vx(i,j,k), vy(i,j,k), vz(i,j,k)/) - ptr%v
                            ptr%ds(MX) = ptr%ds(MX) + (dr(MY)*du(MZ)-dr(MZ)*du(MY))*rhoacc*dv
                            ptr%ds(MY) = ptr%ds(MY) + (dr(MZ)*du(MX)-dr(MX)*du(MZ))*rhoacc*dv
                            ptr%ds(MZ) = ptr%ds(MZ) + (dr(MX)*du(MY)-dr(MY)*du(MX))*rhoacc*dv
                            ! mass and momentum of gas is reduced by accretion
                            ! note that velocity of gas keeps remains by accretion
                            !--------------------------- KS ADDED ------------------------!
                            pr(i,j,k) = pr(i,j,k) * (rho(i,j,k)-rhoacc)/rho(i,j,k) !密度と同じ比率で減少させる
                            !-------------------------------------------------------------!
                            rho(i,j,k) = rho(i,j,k) - rhoacc
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do

!!$       call mpi_allreduce( ptr%dmass, dmass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       if (get_myrank() == PRIMARY_RANK) then
!!$          print *, 'dm (normal) = ', dmass
!!$          call flush(6)
!!$       endif

       ptr => ptr%next
    end do
    ! all reduce of change due to accretion and update particle variables
  end subroutine sp_accretion
#endif !ACCRETION_RHOCR
#ifdef ACCRETION_BONDI
  !-------------------------------------------------------------------------
  ! emulating accretion of Bondi 1952
  !-------------------------------------------------------------------------
  subroutine sp_accretionBondi( dt )
    use mpilib
    use parameter, only : Pi4
    use modelParameter, only : MP_Gconst
    use eos
    use ob_interp
    real(kind=DBL_KIND),intent(IN) :: dt
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND) :: rBH, rAcc, rKernel
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: vx, vy, vz, rho
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: csp, venv2, lambda, rhobar, rhoenv, mdot, dw, wtotal, r2, dv, rhoacc, hmin, mpibuf
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    integer,dimension(MX:MZ) :: ijkgL, ijkgR
    integer :: rank, gid, level, ig, jg, kg, i, j, k
    integer,dimension(MRHO:MVZ) :: mlist
#ifdef MP
    real(kind=DBL_KIND),dimension(MRHO:MVZ+1) :: u0
    real(kind=DBL_KIND) :: p0
#else !MP
    real(kind=DBL_KIND),dimension(MRHO:MVZ) :: u0
#endif !MP
    real(kind=DBL_KIND) :: rho0
    real(kind=DBL_KIND),dimension(MVX:MVZ) :: v0
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr, du
    lambda = exp(1.5d0)*0.25d0
    level = LevelMax
    rAcc = sp_f_rAcc * minval(CellWidth(:,level)) ! sink radius
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       if (ptr%mass < TINY(ptr%mass)) then
          ptr => ptr%next
          cycle
       endif

       ! --------------------------------
       ! rBH, rKernel, venv2, csp (given from point)
       ! --------------------------------
#ifdef MP
       mlist = (/MRHO, MVX, MVY, MVZ, MP/)
       call ob_interpolatedUbyCoordPhys(ptr%r, mlist, u0)
       rho0 = u0(MRHO)
       v0 = u0(MVX:MVZ)
       p0 = u0(MP)
       csp = sqrt(Gamma * p0/rho0)
#else !MP
       mlist = (/MRHO, MVX, MVY, MVZ/)
       call ob_interpolatedUbyCoordPhys(ptr%r, mlist, u0)
       rho0 = u0(MRHO)
       v0 = u0(MVX:MVZ)
       call get_cs(csp, rho0)
#endif !MP
       venv2 = sum( (ptr%v - v0)**2 )
       rBH = MP_Gconst * ptr%mass / (venv2 + csp**2)
       rKernel = min(max(rBH, minval(CellWidth(:,level))/4.d0), rAcc/2)
       ! -------------------------------------
       ! rhobar (weighted average of density)
       ! wtotal (weight of kernel)
       ! -------------------------------------
       rhobar = 0.d0
       wtotal = 0.d0
       posL(:) = ptr%r(:) - rAcc
       posR(:) = ptr%r(:) + rAcc
#ifdef PARTICLECREATION_IN_BOUNDARY
       call sp_restrict_within_boundary(posL, posR)
#endif !PARTICLECREATION_IN_BOUNDARY
       call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
       call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
       do kg = ijkgL(MZ), ijkgR(MZ)
          do jg = ijkgL(MY), ijkgR(MY)
             do ig = ijkgL(MX), ijkgR(MX)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(MRHO, gid)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         r2 = (x(i)-ptr%r(MX))**2 + (y(j)-ptr%r(MY))**2 + (z(k)-ptr%r(MZ))**2
                         if ( r2 < rAcc**2 ) then
                            dw =  exp( -r2/rKernel**2 )
                            rhobar =  rhobar + rho(i,j,k) * dw
                            wtotal = wtotal + dw
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
       call mpi_allreduce(wtotal, mpibuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       wtotal = mpibuf
       call mpi_allreduce(rhobar, mpibuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       rhobar = mpibuf
       rhobar = rhobar / wtotal
       ! ----------------------------------------------
       ! rhoenv (rho_infinity), mdot (accretion rate)
       ! ---------------------------------------------
       rhoenv = rhobar / get_alpha(1.2d0 * minval(CellWidth(:,level))/rBH)
       mdot = Pi4 * rhoenv * rBH**2 * sqrt(lambda**2 * csp**2 + venv2)

       ! debug
!!$       if (get_myrank() == PRIMARY_RANK) then
!!$          print *, 'mdot, venv2, rhobar alpha rhoenv = ', mdot, venv2, rhobar, rhobar/rhoenv, rhoenv
!!$          call flush(6)
!!$       endif
       if ( mdot < 0) then
          HALT
       endif

       ! ------------
       ! accretion
       ! ------------
       hmin = minval(CellWidth(:,level))
       do kg = ijkgL(MZ), ijkgR(MZ)
          do jg = ijkgL(MY), ijkgR(MY)
             do ig = ijkgL(MX), ijkgR(MX)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(MRHO, gid)
                vx => get_Ucomp(MVX, gid)
                vy => get_Ucomp(MVY, gid)
                vz => get_Ucomp(MVZ, gid)
                dv = get_dv(level)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         r2 = (x(i)-ptr%r(MX))**2 + (y(j)-ptr%r(MY))**2 + (z(k)-ptr%r(MZ))**2
                         if ( r2 < rAcc**2 ) then
                            dw =  exp( -r2/rKernel**2 ) /wtotal
                            rhoacc = mdot * dt * dw / dv
                            ! reduce rhoacc due to rotation
                            call reduce_accretion_by_rotation(rhoacc, ptr, vx(i,j,k), vy(i,j,k), vz(i,j,k), x(i), y(j), z(k), hmin, rKernel)
                            rhoacc = min(rhoacc, rho(i,j,k)*0.25d0)
#ifdef RHO_FLOOR
                            rhoacc = min(rhoacc, max(rho(i,j,k)-RHO_FLOOR, 0.d0))
#endif !RHO_FLOOR
                            ! particle gets mass and momentum
                            ptr%dmass = ptr%dmass +  rhoacc*dv
                            ptr%dp(MX) = ptr%dp(MX) + vx(i,j,k)*rhoacc*dv
                            ptr%dp(MY) = ptr%dp(MY) + vy(i,j,k)*rhoacc*dv
                            ptr%dp(MZ) = ptr%dp(MZ) + vz(i,j,k)*rhoacc*dv
                            dr = (/x(i), y(j), z(k)/) - ptr%r
                            du = (/vx(i,j,k), vy(i,j,k), vz(i,j,k)/) - ptr%v
                            ptr%ds(MX) = ptr%ds(MX) + (dr(MY)*du(MZ)-dr(MZ)*du(MY))*rhoacc*dv
                            ptr%ds(MY) = ptr%ds(MY) + (dr(MZ)*du(MX)-dr(MX)*du(MZ))*rhoacc*dv
                            ptr%ds(MZ) = ptr%ds(MZ) + (dr(MX)*du(MY)-dr(MY)*du(MX))*rhoacc*dv
                            ! mass and momentum of gas is reduced by accretion
                            ! note that velocity of gas keeps remains by accretion
                            rho(i,j,k) = rho(i,j,k) - rhoacc
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do

!!$       ! temporally allreduce for outpu dm
!!$       call mpi_allreduce( ptr%dmass, dmass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       if (get_myrank() == PRIMARY_RANK) then
!!$          print *, 'dm (Bondi) =  ', dmass, mdot, rBH, rKernel
!!$          call flush(6)
!!$       endif

       ptr => ptr%next

    end do
  contains
    !-------------------------------------------------------------------------
    ! get alpha
    !-------------------------------------------------------------------------
    function get_alpha(xp) result(alpha)
      real(kind=DBL_KIND),intent(IN) :: xp
      real(kind=DBL_KIND) :: alpha, x0, x1
      integer :: n0, n1
      n0 = min(max(int(xp/sp_BONDI_DX),1), sp_BONDI_NMAX-1)
      n1 = n0 + 1
      x0 = n0*sp_BONDI_DX
      x1 = n1*sp_BONDI_DX
      alpha = ( (x1 - xp)*sp_BONDI_RHO(n0) + (xp - x0)*sp_BONDI_RHO(n1) ) / (x1 - x0)
      alpha = max(alpha, 1.d0)  ! for a very large xp
    end function get_alpha
    !-------------------------------------------------------------------------
    ! reduce accretion by rotation
    !-------------------------------------------------------------------------
    subroutine reduce_accretion_by_rotation(rhoacc, ptr, vx, vy, vz, x, y, z, hmin, rKernel)
      use modelParameter, only : MP_Gconst
      real(kind=DBL_KIND),intent(INOUT) :: rhoacc
      real(kind=DBL_KIND),intent(IN) :: vx, vy, vz, x, y, z, hmin, rKernel
      type(t_spParticle),pointer :: ptr ! (IN) particle
      real(kind=DBL_KIND) :: esp2, jsp2, ekin, rmin, xc, yc, zc, gx, gy, gz, psi, hsc, sinkSize, dx, dy, dz, gm, gm2
      real(kind=DBL_KIND),dimension(MX:MZ) :: dr
      integer,parameter :: N_SubCell = 8
      integer :: i, j, k, n
      sinkSize = hmin * 0.25d0
      ! check : particle is located inside the host cell and its bordering cells?
      if ( (x - ptr%r(MX))**2 + (y - ptr%r(MY))**2 + (z - ptr%r(MZ))**2 < (hmin*1.5d0)**2 ) then
         if ( rKernel >= sinkSize ) rhoacc = 0.d0 ! no-accretion
         return
      endif

      ekin = 0.5d0*( (vx-ptr%v(MX))**2 + (vy-ptr%v(MY))**2 + (vz-ptr%v(MZ))**2 )
      gm = MP_Gconst * ptr%mass
      gm2 = gm**2
      hsc = hmin / N_SubCell    ! cell width of subcell
      n = 0
      do k = 0, N_SubCell-1
         do j = 0, N_SubCell-1
            do i = 0, N_SubCell-1
               ! position of subcell
               xc = x + ( i - (N_SubCell-1)*0.5d0 ) * hsc
               yc = y + ( j - (N_SubCell-1)*0.5d0 ) * hsc
               zc = z + ( k - (N_SubCell-1)*0.5d0 ) * hsc
               dx = xc - ptr%r(MX)
               dy = yc - ptr%r(MY)
               dz = zc - ptr%r(MZ)
               dr = (/xc - ptr%r(MX), yc - ptr%r(MY), zc - ptr%r(MZ)/)
               call sp_gravityOfParticle(dr, ptr%mass, gx, gy, gz, psi) ! gravitational postential
               esp2 = (ekin + psi)*2.d0 ! total energy of subcell * 2
               jsp2 = (dy*vz-dz*vy)**2 + (dz*vx-dx*vz)**2 + (dx*vy-dy*vx)**2 !specific angular momentum^2
               if ( esp2 >= 0.d0 ) then
                  rmin = HUGE(rmin)
               else
                  rmin = -gm / esp2 * ( 1.d0 - sqrt( 1.d0 + jsp2 * esp2 /gm2) )
               end if
               if ( rmin > sinkSize ) n = n + 1
            end do
         end do
      end do
      rhoacc = rhoacc * ( 1.d0 - dble(n)/N_SubCell**3)
    end subroutine reduce_accretion_by_rotation
  end subroutine sp_accretionBondi
#endif !ACCRETION_BONDI
  !-------------------------------------------------------------------------
  ! assign new particle
  !-------------------------------------------------------------------------
  !----------------------------- KS MODIFIED --------------------------!    
  !subroutine sp_newParticle(mass, r, v, s, pid, dmass, dr, dv, dp, ds)
  subroutine sp_newParticle(mass, r, v, s, pid, dmass, dr, dv, dp, ds, &
    t_prev, dm_disk, dJ_disk, mdot_disk, J_disk)
    !----------------------------- KS MODIFIED --------------------------!    
    real(kind=DBL_KIND),intent(IN) :: mass
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: r, v, s
    real(kind=DBL_KIND),intent(IN),optional :: dmass
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN),optional :: dr, dv, dp, ds
    integer,intent(IN),optional :: pid
    !----------------------------- KS ADDED --------------------------!
    real(kind=DBL_KIND),intent(IN),optional :: t_prev, dm_disk, mdot_disk
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN),optional :: dJ_disk, J_disk
    !----------------------------- KS ADDED --------------------------!        
    type(t_spParticle),pointer :: newParticle
    ! prepare a node of link list
    allocate(newParticle)
    newParticle%r = r
    newParticle%v = v
    newParticle%s = s
    newParticle%mass = mass
    if (present(dmass)) then
       newParticle%dmass = dmass
    else
       newParticle%dmass = 0.d0
    end if
    if (present(dp)) then
       newParticle%dp =  dp
    else
       newParticle%dp =  0.d0
    end if
    if (present(dr)) then
       newParticle%dr = dr
    else
       newParticle%dr = 0.d0
    end if
    if (present(dv)) then
       newParticle%dv = dv
    else
       newParticle%dv = 0.d0
    end if
    if (present(ds)) then
       newParticle%ds = ds
    else
       newParticle%ds = 0.d0
    end if
    !----------------------------- KS ADDED --------------------------!
    if (present(t_prev)) then
       newParticle%t_prev = t_prev
    else
       newParticle%t_prev = 0.d0
    end if
    if (present(dm_disk)) then
       newParticle%dm_disk = dm_disk
    else
       newParticle%dm_disk = 0.d0
    end if
    if (present(dJ_disk)) then
       newParticle%dJ_disk = dJ_disk
    else
       newParticle%dJ_disk = 0.d0
    end if    
    if (present(mdot_disk)) then
       newParticle%mdot_disk = mdot_disk
    else
       newParticle%mdot_disk = 0.d0
    end if
    if (present(J_disk)) then
       newParticle%J_disk = J_disk
    else
       newParticle%J_disk = 0.d0
    end if
    !----------------------------- KS ADDED --------------------------!    
    if (present(pid)) then
       if (bool_checkPidDuplicated (pid)) &
            print *, '*** error in sp_newParticle: pid is duplicated.', pid
       newParticle%pid = pid
    else
       newParticle%pid = PidMax + 1
    end if
    PidMax = max(PidMax, newParticle%pid)
    Nparticle = Nparticle + 1
    ! assign newParticle to a link list, Particle
    call sp_appendParticleList(Particle, newParticle)
    !-------------------------------------------------------------------------
    ! returns true if pid is duplicated.
    !-------------------------------------------------------------------------
  contains
    function bool_checkPidDuplicated(pid) result(bool)
      integer,intent(IN) :: pid
      logical :: bool
      type(t_spParticle),pointer :: ptr
      bool = .false.
      ptr => Particle
      do
         if (.not. associated(ptr)) exit
         if ( pid  == ptr%pid ) then
            bool = .true.
            exit
         end if
         ptr => ptr%next
      end do
    end function bool_checkPidDuplicated
  end subroutine sp_newParticle
  !-------------------------------------------------------------------------
  ! delete particle of pid
  !-------------------------------------------------------------------------
  subroutine sp_deleteParticle(pid)
    integer,intent(IN) :: pid
    type(t_spParticle),pointer :: ptr, ptrdel, next
    if (.not. associated(Particle)) then
       print *, '** error in sp_deleteParticle: No particle is assigned', pid
       call flush(6)
       return
    end if

    Nparticle = Nparticle - 1
    ptr => Particle
    next => ptr%next
    if (Particle%pid == pid) then
       deallocate( Particle )
       Particle => next
       return
    end if
    do
       if (.not. associated(ptr%next)) then
          print *, '*** No particle is deleted. Any pid is not match.'
          call flush(6)
          return
       endif
       if (ptr%next%pid == pid) then
          ptrdel => ptr%next
          ptr%next => ptr%next%next
          deallocate( ptrdel )
          return
       end if
       ptr => ptr%next
    end do
  end subroutine sp_deleteParticle
  !-------------------------------------------------------------------------
  ! delete all particles
  !-------------------------------------------------------------------------
  subroutine sp_deleteAllParticles
    type(t_spParticle),pointer :: ptr, next
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       next => ptr%next
       call sp_deleteParticle(ptr%pid)
       ptr => next
    end do
  end subroutine sp_deleteAllParticles
  !-------------------------------------------------------------------------
  ! append a to b. result is a = (a, b)
  ! INOUT: a (list or element)
  ! IN :   b   (list or element)
  ! if both a and b is null pointers, this routine returns error.
  !-------------------------------------------------------------------------
  subroutine sp_appendParticleList(a, b)
    type(t_spParticle),pointer :: a, b
    type(t_spParticle),pointer :: ptr
    if ((.not. associated(a)) .and. (.not. associated(b))) then
       print *, '*** error in sp_appendParticleList: Both pointers are not associated.'
       return
    end if
    if (.not. associated(a)) then
       a => b
       return
    end if
    if (.not. associated(b)) then
       b => a
       return
    end if
    ptr => a
    do
       if (.not. associated(ptr%next)) then
          ptr%next => b
          exit
       end if
       ptr => ptr%next
    end do
  end subroutine sp_appendParticleList
  !-------------------------------------------------------------------------
  ! boundary condition for particles
  !-------------------------------------------------------------------------
  subroutine sp_boundary_eject
    type(t_spParticle),pointer :: ptr, next
    type(t_obRectPhys) :: compDomain, particleDomain
    type(t_obPointPhys) :: pL, pR
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    call ob_computationBoxOfRectPhys( compDomain )
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       posL(:) = ptr%r(:) - sp_SinkRadius
       posR(:) = ptr%r(:) + sp_SinkRadius
       call ob_assignCoordPhysToPointPhys(posL, pL)
       call ob_assignCoordPhysToPointPhys(posR, pR)
       call ob_assignPointPhysToRectPhys(pL, particleDomain, 'L')
       call ob_assignPointPhysToRectPhys(pR, particleDomain, 'R')
       next => ptr%next
       if ( ob_rectPhysComp( particleDomain, compDomain ) /= OB_RECT_INNER ) then
          call sp_boundary_eject_write(ptr)
          call sp_deleteParticle(ptr%pid)
       end if
       ptr => next
    end do
  contains
    subroutine sp_boundary_eject_write(ptr)
      use mpilib
      use string, only : concat
      use io_util, only : read_env
      type(t_spParticle),pointer :: ptr
      character(len=CHARLEN),parameter :: FILENAME_EJECT='logSinkparticle_ejected'
      character(len=CHARLEN) :: file, dir
      integer,parameter :: LUN = 11
      if (get_myrank() /= PRIMARY_RANK) return
      if (.not. associated(Particle)) return
      call read_env('DIR', dir)
      file = concat(dir,FILENAME_EJECT)
      open(LUN, file=file, position='APPEND')
      write(LUN, '(I14, E17.9, I7, 10(1PE17.9))') Step(sp_Level), Time(sp_Level), ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%s
      call flush(LUN)
      close(LUN)
    end subroutine sp_boundary_eject_write
  end subroutine sp_boundary_eject
  !-------------------------------------------------------------------------
  ! merge particles within the sink radius.
  !-------------------------------------------------------------------------
  subroutine sp_merge
    use parameter
    use eos

    use modelParameter, only : MP_CONNECTION_RUN ! KS ADDED
    
    type(t_spParticle),pointer :: ptri, ptrj
    real(kind=DBL_KIND) :: distanceMerge
    integer :: piddel
    !----------------------------- KS ADDED --------------------------------!
    real(kind=DBL_KIND) :: rcm(MX:MZ), vcm(MX:MZ) ! cneter of mass
    real(kind=DBL_KIND) :: dri(MX:MZ), drj(MX:MZ), dui(MX:MZ), duj(MX:MZ)
    !----------------------------- KS ADDED (END) --------------------------!

    distanceMerge = sp_SinkRadius * 2.d0
    ptri => Particle
    do
       if (.not. associated(ptri)) exit
       ptrj => Particle
       do
          if (.not. associated(ptrj)) exit
          if ( ptri%pid == ptrj%pid ) then
             ptrj => ptrj%next
             cycle
          endif
          if ( sum((ptri%r - ptrj%r)**2) > distanceMerge**2 ) then
             ptrj => ptrj%next
             cycle
          end if
          if (ptri%mass < TINY(ptri%mass) .and. ptrj%mass < TINY(ptrj%mass)) then
             ptrj => ptrj%next
             cycle
          end if
          if (get_myrank() == PRIMARY_RANK) then
             print *, '*** particles are merged ', ptri%pid, ptrj%pid
             call flush(6)
          end if
          !----------------------------- KS ADDED --------------------------!
          if (MP_CONNECTION_RUN == 0) then
             ! center of mass position&velocity
             rcm = (ptri%r*ptri%mass + ptrj%r*ptrj%mass)/(ptri%mass + ptrj%mass)
             vcm = (ptri%v*ptri%mass + ptrj%v*ptrj%mass)/(ptri%mass + ptrj%mass)

             ! position&velocity with respect to the center of mass
             dri(:) = ptri%r(:) - rcm(:)
             drj(:) = ptrj%r(:) - rcm(:)
             dui(:) = ptri%v(:) - vcm(:)
             duj(:) = ptrj%v(:) - vcm(:)

             !軽い方のシンクが円盤経由で降着とモデル化
             ptri%dm_disk = ptri%dm_disk + ptrj%mass

             !合体時の重心から見たシンク粒子の角運動量をdJ_diskに加える -> 実質的に円盤面と合体の軌道面を一致させる
             ptri%dJ_disk(MX) = ptri%dJ_disk(MX) + &
                  (dri(MY)*dui(MZ)-dri(MZ)*dui(MY))*ptri%mass + (drj(MY)*duj(MZ)-drj(MZ)*duj(MY))*ptrj%mass
             ptri%dJ_disk(MY) = ptri%dJ_disk(MY) + &
                  (dri(MZ)*dui(MX)-dri(MX)*dui(MZ))*ptri%mass + (drj(MZ)*duj(MX)-drj(MX)*duj(MZ))*ptrj%mass
             ptri%dJ_disk(MZ) = ptri%dJ_disk(MZ) + &
                  (dri(MX)*dui(MY)-dri(MY)*dui(MX))*ptri%mass + (drj(MX)*duj(MY)-drj(MY)*duj(MX))*ptrj%mass
          end if

          !現状の円盤の降着率/角運動量は二つの円盤の和を利用
          ptri%mdot_disk = ptri%mdot_disk + ptrj%mdot_disk
          ptri%J_disk = ptri%J_disk + ptrj%J_disk
          !----------------------------- KS ADDED --------------------------!

          ptri%r = (ptri%r*ptri%mass + ptrj%r*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%v = (ptri%v*ptri%mass + ptrj%v*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%s = (ptri%s*ptri%mass + ptrj%s*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dr = (ptri%dr*ptri%mass + ptrj%dr*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dv = (ptri%dv*ptri%mass + ptrj%dv*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dp = ptri%dp + ptrj%dp
          ptri%ds = ptri%ds + ptrj%ds
          ptri%mass = ptri%mass + ptrj%mass
          ptri%dmass = ptri%dmass + ptrj%dmass
          piddel = ptrj%pid
          ptrj => ptrj%next
          call sp_deleteParticle(piddel)
       end do
       ptri => ptri%next
    end do
  end subroutine sp_merge
  !-------------------------------------------------------------------------
  ! write sink particle data to binary file
  !-------------------------------------------------------------------------
  subroutine sp_write
    use systemcall
    use mpilib
    use string, only : concat, num2char
    use io_util, only : read_env
    integer,parameter :: LUN = 11
    character(len=CHARLEN) :: file, dir, filen
    character(len=2),parameter :: PREFIX = 'sp', SUFFIX = '.d'
    type(t_spParticle),pointer :: ptr
    logical :: exist
    call sp_init
    if (Nparticle == 0) return
    if (get_myrank() /= PRIMARY_RANK) return
    call read_env('DIR', dir)
    file = concat(dir,FILENAME)                   ! sinkparticle.d
    filen = trim(dir)//PREFIX//trim(num2char(Step(LevelMax)))//SUFFIX ! sp1234.d
    open(LUN, file=filen, form='unformatted')
    write(LUN) Nparticle, PidMax
    write(LUN) sp_DtimeCFL
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ! write(LUN) ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds
            !----------------------------- KS MODIFIED --------------------------!
       write(LUN) ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds, &
            ptr%t_prev, ptr%dm_disk, ptr%dJ_disk, ptr%mdot_disk, ptr%J_disk 
            !----------------------------- KS MODIFIED --------------------------!
       ptr => ptr%next
    end do
    call flush(LUN)
    close(LUN)
    inquire(file=filen, exist=exist)
    if (.not. exist) print *, '*** ERROR in sp_write'
    inquire(file=file, exist=exist)
    if (exist) call systemcall_unlink(file)
    call systemcall_symlink(filen, file)
  end subroutine sp_write
  !-------------------------------------------------------------------------
  ! write sink particle data to binary file
  !-------------------------------------------------------------------------
  subroutine sp_read
    use mpilib
    use string, only : concat
    use io_util, only : read_env
    integer,parameter :: LUN = 11
    character(len=CHARLEN) :: file, dir
    integer :: n, nptcle, pid
    real(kind=DBL_KIND),dimension(:,:),allocatable :: buf
    integer,dimension(:),allocatable :: bufi
    real(kind=DBL_KIND),dimension(MX:MZ) :: r, v, dr, dv, dp, s, ds
    real(kind=DBL_KIND) :: mass, dmass
    logical :: exist
    !----------------------------- KS ADDED --------------------------------!
    real(kind=DBL_KIND) :: t_prev, dm_disk, mdot_disk
    real(kind=DBL_KIND),dimension(MX:MZ) :: dJ_disk, J_disk
    !----------------------------- KS ADDED (END) --------------------------!
    
    call sp_init
    if (get_myrank() == PRIMARY_RANK ) then
       call read_env('DIR', dir)
       file = concat(dir,FILENAME)
       inquire(file=file, exist=exist)
    end if
    call mpi_bcast(exist, 1, MPI_LOGICAL, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    if (.not. exist) return
    if (get_myrank() == PRIMARY_RANK ) then
       open(LUN, file=file, form='unformatted')
       read(LUN) nptcle, PidMax
       read(LUN) sp_DtimeCFL
    end if
    call mpi_bcast(nptcle, 1, MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(PidMax, 1, MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sp_DtimeCFL, 1, MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    if (nptcle == 0) return
    allocate( bufi(nptcle) )
    !----------------------------- KS MODIFIED --------------------------!    
    !allocate( buf(1+1+3*7, nptcle ) )
    allocate( buf(1+1+3*7+1*3+3*2, nptcle ) )
    !----------------------------- KS MODIFIED --------------------------!
    if (get_myrank() == PRIMARY_RANK) then
       do n = 1, nptcle
          read(LUN) bufi(n), buf(:, n)
       end do
       close(LUN)
    end if
    call mpi_bcast(bufi, size(bufi), MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(buf,  size(buf),  MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    do n = 1, nptcle
       pid = bufi(n)
       mass = buf(1,n)
       r = buf(2:2+2,n)
       v = buf(5:5+2,n)
       dmass = buf(8,n)
       dr = buf(9:9+2,n)
       dv = buf(12:12+2,n)
       dp = buf(15:15+2,n)
       s  = buf(18:18+2,n)
       ds = buf(21:21+2,n)
       !----------------------------- KS ADDED --------------------------!
       t_prev = buf(24,n)
       dm_disk = buf(25,n)
       dJ_disk = buf(26:26+2,n)
       mdot_disk = buf(29,n)
       J_disk = buf(30:30+2,n)
       !----------------------------- KS ADDED --------------------------!
       !----------------------------- KS MODIFIED --------------------------!    
       ! call sp_newParticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds )
       call sp_newParticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds, &
            t_prev=t_prev, dm_disk=dm_disk, dJ_disk=dJ_disk, mdot_disk=mdot_disk, J_disk=J_disk)
       !----------------------------- KS MODIFIED --------------------------!
    end do
    if (nptcle /= Nparticle) print *, '*** error in sp_read: Nparticle is not consistent', Nparticle, nptcle
    deallocate( bufi, buf )
  end subroutine sp_read
  !-------------------------------------------------------------------------
  ! write sink particle data to logfile
  !-------------------------------------------------------------------------
  subroutine sp_writeLog
    use mpilib
    ! use string, only : concat
    use string, only : concat, num2char  
    use io_util, only : read_env
    integer,parameter :: LUN = 11

    character(len=CHARLEN) :: file, dir
    character(len=CHARLEN) :: file_each ! KS ADDED
    character(len=CHARLEN),parameter :: prefix_log='spLog' !KS ADDED
    character(len=CHARLEN),parameter :: prefix_log2='logsp' !KS ADDED
    integer,parameter :: LUN2 = 12      !KS ADDED
    
    type(t_spParticle),pointer :: ptr
    integer,parameter :: skip = 10

    if (mod(Step(sp_Level), skip) /= 0) return

    call sp_init
    if (Nparticle == 0) return
    if (get_myrank() /= PRIMARY_RANK) return
    if (.not. associated(Particle)) return
    call read_env('DIR', dir)
    file = concat(dir,FILENAME_LOG)

    !------------   KS ADDED (to include parameter info) ------------!
 !    if (access(file) /= 0) then
 !       open(LUN, file=file, status='new')
 !       write(LUN, '(A,/, 2(1PE17.9), I7, 23(1PE17.9))') '#Header:  sp_SinkRadius =   2.938469073116814E-004
 ! sp_Level =           13
 ! sp_SofteningRadius =   2.938469073116814E-004
 ! sp_MaskRadius =   5.876938146233628E-004'
 !    else
 !       open(LUN, file=file, position='APPEND')
 !    end if
    !---------------------------------------------------------------!

    open(LUN, file=file, position='APPEND')       
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       write(LUN, '(I14, 2(1PE17.9), I7, 23(1PE17.9))') Step(sp_Level), Time(sp_Level), Dtime(sp_Level), &
            ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds

       !----------------------------- KS ADDED --------------------------!
       file_each = concat(dir,prefix_log)                   ! /dir/spLog
       file_each = concat(file_each,num2char(ptr%pid))          ! /dir/spLog10
       file_each = concat(file_each,'.dat')                 ! /dir/spLog10.dat
       open(LUN2, file=file_each, position='APPEND')

       !mdot_diskがゼロで無ければ
       if (ptr%mdot_disk > tiny(ptr%mdot_disk)) then
          write(LUN2, '(I10, 14(1PE15.7))') Step(sp_Level), Time(sp_Level)*Unit_yr, &
               ptr%mass*Unit_msun, ptr%dmass/Dtime(sp_Level)*Unit_msun/Unit_yr,&
               ptr%r*Unit_au, ptr%v*Unit_kms, ptr%t_prev*Unit_yr, &
               ptr%mdot_disk*Unit_msun/Unit_yr, ptr%J_disk*Unit_msun*Unit_kms*Unit_au !1*5 + 3*3 = 14 floats

       !mdot_diskがゼロのとき
       else
          !プロットの際に軸を決めるときに困るので、mdot_diskがゼロでJ_diskもゼロベクトルになる時はシンク粒子のスピンを渡す
          write(LUN2, '(I10, 14(1PE15.7))') Step(sp_Level), Time(sp_Level)*Unit_yr, &
               ptr%mass*Unit_msun, ptr%dmass/Dtime(sp_Level)*Unit_msun/Unit_yr,&
               ptr%r*Unit_au, ptr%v*Unit_kms, ptr%t_prev*Unit_yr, &
               ptr%mdot_disk*Unit_msun/Unit_yr, ptr%s !1*5 + 3*3 = 14 floats
       end if
       call flush(LUN2)
       close(LUN2)
       !----------------------------- KS ADDED (END) --------------------------!

       !----------------------------- KS ADDED2 --------------------------!
       !spLogをやめてlogspに移行する予定
       file_each = concat(dir,prefix_log2)                   ! /dir/logsp
       file_each = concat(file_each,num2char(ptr%pid))          ! /dir/logsp10
       file_each = concat(file_each,'.dat')                 ! /dir/logsp10.dat
       open(LUN2, file=file_each, position='APPEND')

       write(LUN2, '(I10, 17(1PE15.7))') Step(sp_Level), Time(sp_Level)*Unit_yr, &
               ptr%mass*Unit_msun, ptr%dmass/Dtime(sp_Level)*Unit_msun/Unit_yr,&
               ptr%r*Unit_au, ptr%v*Unit_kms, ptr%s, &
               ptr%t_prev*Unit_yr, ptr%mdot_disk*Unit_msun/Unit_yr, &
               ptr%J_disk*Unit_msun*Unit_kms*Unit_au !1*5 + 3*4 = 17 floats
       call flush(LUN2)
       close(LUN2)
       !----------------------------- KS ADDED2 (END) --------------------------!       
       ptr => ptr%next
    end do
    call flush(LUN)
    close(LUN)
  end subroutine sp_writeLog
  !-------------------------------------------------------------------------
  ! gravity of particle. This subroutine should be inline-expanded.
  !-------------------------------------------------------------------------
  subroutine sp_gravityOfParticle(dr, mass, gx, gy, gz, psi, no_softening)
    use modelParameter, only : MP_Gconst
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: dr
    real(kind=DBL_KIND),intent(IN) :: mass
    real(kind=DBL_KIND),intent(OUT) :: gx, gy, gz
    real(kind=DBL_KIND),intent(OUT),optional :: psi
    logical,optional :: no_softening
    real(kind=DBL_KIND) :: flagIn, flagOut, gabs, gm, sradii, sradii2, sradii3, drabs
    logical :: no_soft
    if (present(no_softening)) then
       no_soft = no_softening
    else
       no_soft = .false.
    endif

    sradii = sp_SofteningRadius
    sradii2 = sradii**2
    sradii3 = sradii**3
    gm = MP_Gconst * mass
    drabs = sqrt(sum(dr**2))
    if (no_soft) then
       flagOut = 1
    else
       flagOut = min(int( drabs /sradii), 1) ! outside sink radius
    endif
    flagIn  = 1 - flagOut            ! inside sink radius
    gabs = - gm * ( flagIn / sradii3 + flagOut / (drabs**3 + flagIn) )
    gx = gabs * dr(MX)
    gy = gabs * dr(MY)
    gz = gabs * dr(MZ)
    if (present(psi)) then
       psi = - gm * ( &
            flagIn * (1.5d0 - 0.5d0*drabs**2/sradii2) / sradii + &
            flagOut / (drabs + flagIn) &
         )
    endif
  end subroutine sp_gravityOfParticle
  !-------------------------------------------------------------------------
  ! Interaction of gravity between fluid and a particle.
  ! Gas velocity is updated bacause of gravity of a particle.
  ! Gravity of gas acting on the particle is obtained.
  ! 
  ! INPUT:
  !    ptr = a pointer of particle
  !    dt  = time step
  ! OUTPUT:
  !    gravF = gravity from fluid to a particle
  !-------------------------------------------------------------------------
  subroutine sp_gravityFluid2particle(ptr, gravF, dt)
    use mpilib
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:MZ),intent(OUT) :: gravF
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(MX:MZ) :: gravBlock
    integer :: level, n, gid
    gravF = 0.d0
#ifdef FLUID_GRAVITY_OFF
    return
#else !FLUID_GRAVITY_OFF
    myrank = get_myrank()
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          call sp_gravityFluid2particleForBlck(ptr, gravBlock, dt, gid)
#ifdef MOMENTUM_CONSERVATION
          ! gravity of gas is evaluated as back reaction of gravity of particle acting on gas.
          if (ChildGid(Left, Left, Left, gid, myrank) == Undefi) then
             gravF = gravF + gravBlock
          end if
#endif !MOMENTUM_CONSERVATION
       end do
    end do
#ifdef MOMENTUM_CONSERVATION
    call mpi_allreduce(MPI_IN_PLACE, gravF, size(gravF), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else !MOMENTUM_CONSERVATION
    call sp_gravityOfFluid(ptr, gravF) ! get gravity of gas acting on ptr.
#endif !MOMENTUM_CONSERVATION

#endif !FLUID_GRAVITY_OFF
  end subroutine sp_gravityFluid2particle
#ifndef FLUID_GRAVITY_OFF
  !-------------------------------------------------------------------------
  ! Interaction of gravity between a block and a particle.
  ! INPUT:
  !    ptr = a pointer of particle
  !    dt  = time step
  !    gid = block id
  ! OUTPUT:
  !    grav = gravity from fluid of given block to a particle
  !-------------------------------------------------------------------------
  subroutine sp_gravityFluid2particleForBlck(ptr, grav, dt, gid)
    use eos, only : w2u, u2w
    use mpilib ! KS DEBUG
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:MZ),intent(OUT) :: grav
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, w
    integer :: level, i, j, k, ic, jc, kc
    real(kind=DBL_KIND) :: dv, dvsc, xc, yc, zc, dtdvrho, gxp, gyp, gzp
    real(kind=DBL_KIND),dimension(MX:MZ) :: h, hsc, bl, br, pl, pr, cl, cr, dr
    integer,parameter :: N_SubCell = 8

    logical :: isNotFinite ! defined in naninf.F90  (KS DEBUG)


    level = get_level(gid)
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    u => get_Up(gid)
    dv = get_dv(level)
    allocate( w(ARRAYSIZE4(u)) )
    call u2w(u, w, dv)

    h(:) = CellWidth(:,level)
    bl = (/x(Imin), y(Jmin), z(Kmin)/) - h*0.5d0 ! left side of block
    br = (/x(Imax), y(Jmax), z(Kmax)/) + h*0.5d0 ! right side of block
    hsc(:) = CellWidth(:, level) / N_SubCell    ! cell width of subcell
    dvsc = dv/N_SubCell**3
    grav = 0.d0

    pl = ptr%r - sp_SofteningRadius
    pr = ptr%r + sp_SofteningRadius
    if ( max(pl(MX),bl(MX)) < min(pr(MX),br(MX)) .and. &
         max(pl(MY),bl(MY)) < min(pr(MY),br(MY)) .and. &
         max(pl(MZ),bl(MZ)) < min(pr(MZ),br(MZ)) ) then ! block is overlaped with particle radius
       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                cl = (/x(i), y(j), z(k)/) - h*0.5d0
                cr = (/x(i), y(j), z(k)/) + h*0.5d0
                if ( max(pl(MX),cl(MX)) < min(pr(MX),cr(MX)) .and. &
                     max(pl(MY),cl(MY)) < min(pr(MY),cr(MY)) .and. &
                     max(pl(MZ),cl(MZ)) < min(pr(MZ),cr(MZ)) ) then ! cell is overlaped with particle radius
                   dtdvrho = dt * dvsc * u(i,j,k,MRHO)
                   do kc = 0, N_SubCell-1
                      do jc = 0, N_SubCell-1
                         do ic = 0, N_SubCell-1
                            xc = x(i) + ( ic - (N_SubCell-1)*0.5d0 ) * hsc(MX)
                            yc = y(j) + ( jc - (N_SubCell-1)*0.5d0 ) * hsc(MY)
                            zc = z(k) + ( kc - (N_SubCell-1)*0.5d0 ) * hsc(MZ)
                            dr = (/ xc-ptr%r(MX), yc-ptr%r(MY), zc-ptr%r(MZ) /)
                            call sp_gravityOfParticle(dr, ptr%mass, gxp, gyp, gzp)
#ifndef TESTPARTICLE
                            w(i,j,k,MVX) = w(i,j,k,MVX) + gxp*dtdvrho
                            w(i,j,k,MVY) = w(i,j,k,MVY) + gyp*dtdvrho
                            w(i,j,k,MVZ) = w(i,j,k,MVZ) + gzp*dtdvrho
#ifdef MP
                            w(i,j,k,MP)  = w(i,j,k,MP) + dtdvrho*( &
                                  u(i,j,k,MVX)*gxp &
                                 +u(i,j,k,MVY)*gyp &
                                 +u(i,j,k,MVZ)*gzp )
#endif !MP
#endif !TESTPARTICLE
                            grav = grav - (/gxp, gyp, gzp/) * dvsc * u(i,j,k,MRHO)
                         end do
                      end do
                   end do
                else         ! cell is outside particle radius
                   dtdvrho = dt * dv * u(i,j,k,MRHO)
                   dr = (/ x(i)-ptr%r(MX), y(j)-ptr%r(MY), z(k)-ptr%r(MZ) /)
                   call sp_gravityOfParticle(dr, ptr%mass, gxp, gyp, gzp)
#ifndef TESTPARTICLE
                   w(i,j,k,MVX) = w(i,j,k,MVX) + gxp*dtdvrho
                   w(i,j,k,MVY) = w(i,j,k,MVY) + gyp*dtdvrho
                   w(i,j,k,MVZ) = w(i,j,k,MVZ) + gzp*dtdvrho
#ifdef MP
                   w(i,j,k,MP)  = w(i,j,k,MP) + dtdvrho*( &
                         u(i,j,k,MVX)*gxp &
                        +u(i,j,k,MVY)*gyp &
                        +u(i,j,k,MVZ)*gzp )
#endif !MP
#endif !TESTPARTICLE
                   grav = grav - (/gxp, gyp, gzp/) * dv * u(i,j,k,MRHO)
                end if
             end do
          end do
       end do
    else                     ! block is outside particle radius
       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                dtdvrho = dt * dv * u(i,j,k,MRHO)
                dr = (/ x(i)-ptr%r(MX), y(j)-ptr%r(MY), z(k)-ptr%r(MZ) /)
                call sp_gravityOfParticle(dr, ptr%mass, gxp, gyp, gzp)
#ifndef TESTPARTICLE
                   w(i,j,k,MVX) = w(i,j,k,MVX) + gxp*dtdvrho
                   w(i,j,k,MVY) = w(i,j,k,MVY) + gyp*dtdvrho
                   w(i,j,k,MVZ) = w(i,j,k,MVZ) + gzp*dtdvrho
#ifdef MP
                   w(i,j,k,MP)  = w(i,j,k,MP) + dtdvrho*( &
                         u(i,j,k,MVX)*gxp &
                        +u(i,j,k,MVY)*gyp &
                        +u(i,j,k,MVZ)*gzp )
#endif !MP
#endif !TESTPARTICLE
                   grav = grav - (/gxp, gyp, gzp/) * dv * u(i,j,k,MRHO)
             end do
          end do
       end do
    end if
    call w2u(w, u, dv)

    !NaN/Inf check for pressure (KS DEBUG)
    do i = Imin,Imax
       do j = Jmin,Jmax
          do k = Kmin,Kmax
             if (isNotFinite(u(i,j,k,MP))) then
                print '(A,5I6, 1P3E15.7)', "(spgf2pfb) NaN/Inf found in u(i,j,k,MP)",&
                     get_myrank(),gid,i,j,k,x(i),y(j),z(k)
                print '(A,1P21E15.7)','w: ',w(i,j,k,:)
                print '(A,1P21E15.7)','u: ',u(i,j,k,:)
                print '(A)','stopping...'
                stop
             end if
          end do
       end do
    end do
    
    deallocate( w )
    if (ptr%mass > TINY(ptr%mass)) then
       grav = grav / ptr%mass
    else
       grav = 0.d0
    end if
  end subroutine sp_gravityFluid2particleForBlck
#endif !FLUID_GRAVITY_OFF
  !-------------------------------------------------------------------------
  ! get gravity of fluid
  !-------------------------------------------------------------------------
  subroutine sp_gravityOfFluid(ptr, grav)
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:MZ),intent(OUT) :: grav
    real(kind=DBL_KIND) :: vol
    real(kind=DBL_KIND) :: vol_spsph ! KS ADDED
    
#ifndef FLUID_GRAVITY_OFF
    call sp_sumUinSphere(ptr%r, sp_SofteningRadius, (/MGX, MGY, MGZ/), grav, vol)
    grav = grav/vol
    !------------------------- KS ADDED -----------------------------!
    ! check whether vol is consistent with sink volume
    vol_spsph = 4.*3.14/3.*sp_SofteningRadius**3 ! volume of sink particle sphere 
    if (abs((vol-vol_spsph)/vol_spsph) > 0.1) then
       if (get_myrank() == PRIMARY_RANK) &
            print '(A,1P10E15.7)', "(sp_gravityOfFluid) ** WARNING ** sp volume is not consistent: ", &
            ptr%r, vol, vol_spsph
    end if
    !----------------------------------------------------------------!
    
#else !FLUID_GRAVITY_OFF
    grav = 0.d0
#endif !FLUID_GRAVITY_OFF
  end subroutine sp_gravityOfFluid
  !-------------------------------------------------------------------------
  ! get sum value of u within radius
  !-------------------------------------------------------------------------
  subroutine sp_sumUinSphere(point, radius, mlist, usum, vol)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: point
    real(kind=DBL_KIND),intent(IN) :: radius
    integer,intent(IN),dimension(:) :: mlist
    real(kind=DBL_KIND),dimension(size(mlist)),intent(OUT) :: usum
    real(kind=DBL_KIND),intent(OUT) :: vol
    integer,parameter :: N_SubCell = 8
    integer :: level, gid, rank, i,j,k, m, ic, jc, kc, ig, jg, kg
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    integer,dimension(MX:MZ) :: ijkgL, ijkgR
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: u
    real(kind=DBL_KIND) :: dv, hsc(MX:MZ), xc, yc, zc
    real(kind=DBL_KIND),dimension(size(mlist)+1) :: buf, bufr
    level = sp_Level
    dv = get_dv(level) / N_SubCell**3
    hsc(:) = CellWidth(:, level) / N_SubCell    ! cell width of subcell
    usum(:) = 0.d0
    vol = 0.d0
    ! evaluate blocks
    posL(:) = point(:) - radius
    posR(:) = point(:) + radius
#ifdef PARTICLECREATION_IN_BOUNDARY
    call sp_restrict_within_boundary(posL, posR)
#endif !PARTICLECREATION_IN_BOUNDARY
    call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
    call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
             ! !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!
             ! if (get_myrank() == PRIMARY_RANK) then
             !    print '(A, 6I6, 1P4E15.7)', "(sp_sumu) KS DEBUG", ig,jg,kg,level,gid,rank, point(:), radius
             ! end if
             ! !---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---- KS DEBUG ---!       
             if (gid == Undefi) cycle
             if (rank /= get_myrank() ) cycle
             x => get_Xp(gid)
             y => get_Yp(gid)
             z => get_Zp(gid)
             do m = lbound(mlist, 1), ubound(mlist, 1)
                u => get_Ucomp(mlist(m), gid)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         ! loop for subcell
                         do kc = 0, N_SubCell-1
                            do jc = 0, N_SubCell-1
                               do ic = 0, N_SubCell-1
                                  ! position of subcell
                                  xc = x(i) + ( ic - (N_SubCell-1)*0.5d0 ) * hsc(MX)
                                  yc = y(j) + ( jc - (N_SubCell-1)*0.5d0 ) * hsc(MY)
                                  zc = z(k) + ( kc - (N_SubCell-1)*0.5d0 ) * hsc(MZ)
                                  if ( (xc-point(MX))**2 + (yc-point(MY))**2 + (zc-point(MZ))**2 <= radius**2 ) then
                                     usum(m) = usum(m) + u(i,j,k)*dv
                                     if (m == 1) vol = vol + dv
                                  end if
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    buf(1) = vol
    buf(2:size(mlist)+1) = usum(:)
    call mpi_allreduce( buf, bufr, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    vol = buf(1)
    usum(:) = buf(2:size(mlist)+1)
  end subroutine sp_sumUinSphere
  !-------------------------------------------------------------------------
  ! update accretion
  !-------------------------------------------------------------------------
  subroutine sp_accretionUpdate
    use modelParameter, only : MP_CONNECTION_RUN ! KS ADDED

    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND) :: mass_new
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ! accretion is updated (ptr%dp, ptr%dmass, ptr%ds)
       mass_new = ptr%mass + ptr%dmass
       ptr%v = (ptr%v * ptr%mass + ptr%dp) / mass_new
       ptr%s = (ptr%s * ptr%mass + ptr%ds) / mass_new
       ptr%mass = mass_new

       !----- KS ADDED ----!
       if (MP_CONNECTION_RUN == 0) then       !通常の場合
          ptr%dm_disk = ptr%dm_disk + ptr%dmass
          ptr%dJ_disk = ptr%dJ_disk + ptr%ds
       else                                   !connection runの場合 (人工的に降着率が上がるのを防止)
          ptr%dm_disk = ptr%dm_disk + ptr%mdot_disk * Dtime(sp_Level) !円盤の降着率を仮定
          ptr%dJ_disk = ptr%dJ_disk + ptr%ds !円盤の角運動量は目安としてしか使っていないので、あまり気にしない
       end if
          
       !----- KS ADDED ----!

!!$       ! for L1551NE
!!$       ! dr ... mass
!!$       ! dp ... momentum
!!$       ! ds ... spin angular momentum
!!$       ptr%dr(MX) = ptr%dr(MX) + ptr%dmass
!!$       ptr%dv = ptr%dv + ptr%dp
!!$       ptr%s  = ptr%s  + ptr%ds

       ptr => ptr%next
    end do
  end subroutine sp_accretionUpdate
  !-------------------------------------------------------------------------
  ! update advection
  !-------------------------------------------------------------------------
  subroutine sp_advectionUpdate
    type(t_spParticle),pointer :: ptr
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ! advection is updated (ptr%dr, ptr%dv)
       ptr%r = ptr%r + ptr%dr
       ptr%v = ptr%v + ptr%dv
       ptr => ptr%next
    end do
  end subroutine sp_advectionUpdate
  !-------------------------------------------------------------------------
  ! restrict Dtime in Hyddrodynamics
  !-------------------------------------------------------------------------
  subroutine sp_restrictCFL(dt)
    real(kind=DBL_KIND),intent(INOUT) :: dt
    if (Nparticle == 0) return
#ifdef ADVECTION_ON
    dt = min(dt, sp_CVS * sp_DtimeCFL)
#endif !ADVECTION_ON
  end subroutine sp_restrictCFL
  !-------------------------------------------------------------------------
  ! test for IO
  !-------------------------------------------------------------------------
  subroutine sp_testIO
    use mpilib
    type(t_spParticle),pointer :: ptr
    if (.not. associated(Particle)) return
    call sp_write
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** a', ptr%pid, ptr%mass ; call flush(6)
       ptr => ptr%next
    end do
    print *, '** a', Nparticle, PidMax ; call flush(6)
    call sp_deleteAllParticles
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** b', ptr%pid, ptr%mass ; call flush(6)
       ptr => ptr%next
    end do
    print *, '** b', Nparticle, PidMax ; call flush(6)
    call sp_read
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** c', ptr%pid, ptr%mass ; call flush(6)
       ptr => ptr%next
    end do
    print *, '** c', Nparticle, PidMax ; call flush(6)
  end subroutine sp_testIO
#ifdef ACCRETION_BONDI
  !-------------------------------------------------------------------------
  ! read Bobdi data (sp_BONDI_NMAX, sp_BONDI_RHO, sp_BONDI_DX)
  !-------------------------------------------------------------------------
  subroutine sp_readBondi
    use mpilib
    character(len=*),parameter :: BONDI_FILE  = '../Bondi/bondi.dat'
    integer,parameter :: Lun = 11
    real(kind=DBL_KIND) :: x_dummy, y_dummy, x0
    integer :: n
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'BONDI_FILE = '//BONDI_FILE
       open(LUN, file=BONDI_FILE)
       read(LUN, '(I8)') sp_BONDI_NMAX
       print *, 'sp_BONDI_NMAX =', sp_BONDI_NMAX
    endif
    call mpi_bcast(sp_BONDI_NMAX, 1, MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    allocate(sp_BONDI_RHO(sp_BONDI_NMAX))
    if (get_myrank() == PRIMARY_RANK) then
       read(LUN, '(3(1PE12.5))') x0, y_dummy, sp_BONDI_RHO(1)
       read(LUN, '(3(1PE12.5))') x_dummy, y_dummy, sp_BONDI_RHO(2)
       sp_BONDI_DX = x_dummy - x0
       do n = 3, sp_BONDI_NMAX
          read(LUN, '(3(1PE12.5))') x_dummy, y_dummy, sp_BONDI_RHO(n)
       end do
       close(LUN)
       print *, 'sp_BONDI_DX = ', sp_BONDI_DX
    endif
    call mpi_bcast(sp_BONDI_RHO, size(sp_BONDI_RHO), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sp_BONDI_DX, 1, MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
  end subroutine sp_readBondi
#endif !ACCRETION_BONDI
  !-------------------------------------------------------------------------
  ! refinement condition due to sink particles
  ! INPUT:
  !   gid = grid id for evaluated grid.
  ! INOUT:
  !   bool = .true. if this grid should be refinement.
  !-------------------------------------------------------------------------
  subroutine sp_refineCond( gid, bool )
    integer,intent(IN) :: gid
    logical,intent(INOUT) :: bool
    integer :: thisLevel
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: margin
    real(kind=DBL_KIND),dimension(MX:MZ) :: hh
    call sp_init
    thisLevel = get_level(gid)
    if (thisLevel + 1 > sp_Level) then ! over refinement
       bool = .false.
       return
    end if
    if ( bool ) return                      ! if true return

    margin = sp_SinkRadius
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    hh = CellWidth(:,get_level(gid)) * 0.5d0 ! half cell width
    ptr => Particle
    do
       if (.not. associated(ptr)) exit

       if ( &
            max(ptr%r(MX) - margin, x(Imin)-hh(MX)) <= min(ptr%r(MX) + margin, x(Imax)+hh(MX)) .and. &
            max(ptr%r(MY) - margin, y(Jmin)-hh(MY)) <= min(ptr%r(MY) + margin, y(Jmax)+hh(MY)) .and. &
            max(ptr%r(MZ) - margin, z(Kmin)-hh(MZ)) <= min(ptr%r(MZ) + margin, z(Kmax)+hh(MZ))) then
          bool = .true.
          return
       end if
       ptr => ptr%next
    end do
  end subroutine sp_refineCond
  !-------------------------------------------------------------------------
  ! refinement condition due to sink particles (MODIFIED BY KS)
  ! INPUT:
  !   gid = grid id for evaluated grid.
  ! INOUT:
  !   bool = .true. if this grid should be refinement.
  !-------------------------------------------------------------------------
  subroutine sp_refineCond_KS( gid, bool )
    integer,intent(IN) :: gid
    logical,intent(INOUT) :: bool
    integer :: thisLevel
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: margin
    real(kind=DBL_KIND),dimension(MX:MZ) :: hh
    call sp_init
    thisLevel = get_level(gid)

    if (thisLevel + 1 > sp_Level) then ! over refinement
       bool = .false.
       return
    end if

    if ( bool ) return                      ! if true return

    !sink粒子を含むグリッドは孫が必ず存在するので、それによってrefinementが要請されている
    !不良グリッドの除外で想定が外れる可能性があるので、下のレベルにもsink粒子を持つことによる細分化を波及させる
    ! if ( thisLevel + 1 /= sp_Level ) return ! Keep sp_Level to be LevelMax

    !marginを復活させる（一時期margin無しで計算していた）
    margin = sp_SinkRadius
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    hh = CellWidth(:,get_level(gid))*0.5d0
    !箱の内外の判定の際にセル境界とセル中心のズレを考慮する
    
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       if ( &
            max(ptr%r(MX) - margin, x(Imin)-hh(MX)) <= min(ptr%r(MX) + margin, x(Imax)+hh(MX)) .and. &
            max(ptr%r(MY) - margin, y(Jmin)-hh(MY)) <= min(ptr%r(MY) + margin, y(Jmax)+hh(MY)) .and. &
            max(ptr%r(MZ) - margin, z(Kmin)-hh(MZ)) <= min(ptr%r(MZ) + margin, z(Kmax)+hh(MZ))) then
          bool = .true.
          return
       end if
       ! if ( &
       !      x(Imin)-0.5*h(MX) <= ptr%r(MX) .and. ptr%r(MX) <= x(Imax)+0.5*h(MX) .and. &
       !      y(Jmin)-0.5*h(MY) <= ptr%r(MY) .and. ptr%r(MY) <= y(Jmax)+0.5*h(MY) .and. &
       !      z(Kmin)-0.5*h(MY) <= ptr%r(MZ) .and. ptr%r(MZ) <= z(Kmax)+0.5*h(MZ)) then
       !    bool = .true.
       !    return
       ! end if
       ptr => ptr%next
    end do
  end subroutine sp_refineCond_KS
  ! ----------------------------------------------------------------------------
  ! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
  ! analytical algorithm.
  ! Only the diagonal and upper triangular parts of A are accessed. The access
  ! is read-only.
  ! ----------------------------------------------------------------------------
  ! Parameters:
  !   A: The symmetric input matrix
  !   W: Storage buffer for eigenvalues
  ! ----------------------------------------------------------------------------
  SUBROUTINE DSYEVC3(A, W)
    !     .. Arguments ..
    real(kind=DBL_KIND),intent(IN) :: A(3,3)
    real(kind=DBL_KIND),intent(OUT) :: W(3)

    !     .. Parameters ..
    real(kind=DBL_KIND),parameter :: SQRT3 = 1.73205080756887729352744634151D0

    !     .. Local Variables ..
    DOUBLE PRECISION M, C1, C0
    DOUBLE PRECISION DE, DD, EE, FF
    DOUBLE PRECISION P, SQRTP, Q, C, S, PHI

    !     Determine coefficients of characteristic poynomial. We write
    !           | A   D   F  |
    !      A =  | D*  B   E  |
    !           | F*  E*  C  |
    DE    = A(1,2) * A(2,3)
    DD    = A(1,2)**2
    EE    = A(2,3)**2
    FF    = A(1,3)**2
    M     = A(1,1) + A(2,2) + A(3,3)
    C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
    C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0D0 * A(1,3)*DE
    P     = M**2 - 3.0D0 * C1
    Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
    SQRTP = SQRT(ABS(P))
    PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1)  + C0 * (Q + (27.0D0/4.0D0)*C0) )
    PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)
    C     = SQRTP * COS(PHI)
    S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)
    W(2) = (1.0D0/3.0D0) * (M - C)
    W(3) = W(2) + S
    W(1) = W(2) + C
    W(2) = W(2) - S
  END SUBROUTINE DSYEVC3
  !-------------------------------------------------------------------------
  ! get sp_Level
  !-------------------------------------------------------------------------
  function sp_getLevel() result(level)
    integer :: level
    call sp_init
    level = sp_Level
  end function sp_getLevel
  !-------------------------------------------------------------------------
  !  get sp_Rhocr
  !-------------------------------------------------------------------------
  function sp_getRhocr() result(rhocr)
    real(kind=DBL_KIND) :: rhocr
    rhocr = sp_Rhocr
  end function sp_getRhocr
  !-------------------------------------------------------------------------
  ! get number of particles
  !-------------------------------------------------------------------------
  function sp_getNparticle() result(number)
    integer :: number
    number = Nparticle
  end function sp_getNparticle
  !-------------------------------------------------------------------------
  ! get sp_SinkRadius (KS ADDED)
  !-------------------------------------------------------------------------
  function sp_getSinkRadius() result(SinkRadius)
    real(kind=DBL_KIND) :: SinkRadius
    SinkRadius = sp_SinkRadius
  end function sp_getSinkRadius
  !-------------------------------------------------------------------------
  ! pack sink particle data into array
  !  Modified to inlude pmdot (KS MDOIFIED)
  !
  ! OUTPUT:
  !  np ... number of sink particle
  !  pmass ... a vector having masses of sink paricles
  !  pmdot ... a vector having mdot of sink paricles (optional)
  !  pr ...... a vector having positions of sink particles (optional)
  !  pv ...... a vector having velocities of sink particles (optional)
  !  ps ...... a vector having spin of sink particles (optional)
  !  pid ..... a vector having particle ID number (optional)
  !-------------------------------------------------------------------------
  !----------------------------- KS MODIFIED --------------------------!    
  ! subroutine sp_sinkdata2array(np, pmass, pmdot, pr, pv, ps, pid)
  subroutine sp_sinkdata2array(np, pmass, pmdot, pr, pv, ps, pid, pmdot_disk, pJ_disk, pmass_prev, pt_prev, pdm_disk)  
    !----------------------------- KS MODIFIED --------------------------!        
    integer,intent(OUT) :: np
    real(kind=DBL_KIND),dimension(:),intent(OUT) :: pmass
    real(kind=DBL_KIND),dimension(MX:MZ,size(pmass)),intent(OUT),optional :: pr, pv, ps
    integer,dimension(size(pmass)),intent(OUT),optional :: pid
    !----------------------------- KS ADDED --------------------------!
    real(kind=DBL_KIND),dimension(:),intent(OUT),optional :: pmdot ! obsolete
    real(kind=DBL_KIND),dimension(size(pmass)),intent(OUT),optional :: pmdot_disk, pmass_prev, pt_prev, pdm_disk
    real(kind=DBL_KIND),dimension(MX:MZ,size(pmass)),intent(OUT),optional :: pJ_disk
    !----------------------------- KS ADDED --------------------------!    
    type(t_spParticle),pointer :: ptr
    np = 0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       np = np + 1
       pmass(np) = ptr%mass

       if (present(pr)) pr(:, np) = ptr%r(:)
       if (present(pv)) pv(:, np) = ptr%v(:)
       if (present(ps)) ps(:, np) = ptr%s(:)
       if (present(pid)) pid(np) = ptr%pid

       !----------------------------- KS ADDED --------------------------!
       if (present(pmdot)) then
          if (Dtime(sp_Level) == 0.d0) then
             print *, "sp_sinkdata2array: Dtime(sp_Level) is zero, stopping..."
             stop
          end if
          pmdot(np) = ptr%dmass/Dtime(sp_Level)
       end if
       if (present(pmdot_disk)) pmdot_disk(np) = ptr%mdot_disk
       if (present(pJ_disk)) pJ_disk(:, np) = ptr%J_disk(:)
       if (present(pmass_prev)) pmass_prev(np) = ptr%mass - ptr%dm_disk ! diskに滞まっているmassを差っ引くとt=t_prevでのmassが得られる
       if (present(pt_prev)) pt_prev(np) = ptr%t_prev
       if (present(pdm_disk)) pdm_disk(np) = ptr%dm_disk
       !----------------------------- KS ADDED --------------------------!       

       ptr => ptr%next
    end do
  end subroutine sp_sinkdata2array
  !-------------------------------------------------------------------------
  ! check whether given position is inside a sink cell (KS ADDED)
  ! INPUT:
  !   pos ... given position
  ! OUTPUT:
  !   true if inside a sink cell
  !-------------------------------------------------------------------------
  function sp_is_inside_sink(pos) result(bool)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: pos
    logical :: bool
    
    integer :: np
    real(kind=DBL_KIND),dimension(MX:MZ) :: pr
    type(t_spParticle),pointer :: ptr
    np = 0
    ptr => Particle
    bool = .False.
    
    !loop over sink particles
    do
       if (.not. associated(ptr)) exit
       if ( (pos(MX)-ptr%r(MX))**2 + (pos(MY)-ptr%r(MY))**2 + (pos(MZ)-ptr%r(MZ))**2 < sp_SinkRadius**2 ) &
            bool = .True.
       ptr => ptr%next
    end do
  end function sp_is_inside_sink
#ifdef PARTICLECREATION_IN_BOUNDARY
  !-------------------------------------------------------------------------
  ! restrict sink region within computational domain
  !-------------------------------------------------------------------------
  subroutine sp_restrict_within_boundary(posL, posR)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(INOUT) :: posL, posR
    type(t_obRectPhys) :: compDomain, particleDomain
    type(t_obPointPhys) :: pL, pR
    integer :: n
    call ob_computationBoxOfRectPhys( compDomain )
    call ob_extractPointPhysFromRectPhys(pL, compDomain, 'L')
    call ob_extractPointPhysFromRectPhys(pR, compDomain, 'R')
    posL = max(posL, pL%p + CellWidth(:, sp_Level)*0.5d0)
    posR = min(posR, pR%p - CellWidth(:, sp_Level)*0.5d0)
  end subroutine sp_restrict_within_boundary
#endif !PARTICLECREATION_IN_BOUNDARY
  !-------------------------------------------------------------------------
  ! initialize this module
  !-------------------------------------------------------------------------
  subroutine sp_init
    use mpilib
    use parameter, only : Pi
    use modelParameter, only : MP_spNcr,MP_spCs,MP_spRadius_cell,MP_spRadius_lamJ, MP_Gconst
#if MODEL_ART > 0
    use primordial,only: yhe    !KS MODIFIED
#endif !MODEL_ART    
    use eos
    real(kind=DBL_KIND) :: csp, jlength, hmax
    integer :: level
    if (Initialized) return
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'initialize sink particle'
       call flush(6)
    end if

    ! ------------------------------
    ! set sp_Rhocr and sp_RhocrCreate
    ! ------------------------------
#if MODEL_ART > 0
    sp_Rhocr       = MP_spNcr * cgs_mh * (1.+4.*yHe) / Unit_rho  ! KS MODIFIED
#else !MODEL_ART    
    sp_Rhocr       = MP_spNcr * cgs_mh / Unit_rho  ! KS MODIFIED
#endif !MODEL_ART    
    sp_RhocrCreate = sp_Rhocr    ! KS MODIFIED

    ! --------------------
    ! define sp_SinkRadius
    ! --------------------
    ! sound speed for gas at sp_Ncr
    csp = MP_spCs/Unit_v    ! KS MODIFIED
    jlength = csp * sqrt(Pi / (MP_Gconst * sp_Rhocr))
    do level = Lmin, Lmax
       hmax = maxval(CellWidth(:,Lmin))/2.d0**(level-Lmin)
       !if ( jlength >= JeansConst * hmax ) then ! jeans condition

       !------------- KS DEBUG -------------!
       ! if (get_myrank() == PRIMARY_RANK) then
       !    print *, "KS DEBUG", level, hmax, jlength * MP_spRadius_lamJ, MP_spRadius_cell * hmax * 0.99
       ! end if
       !------------- KS DEBUG (END) ---------!
       !       if ( jlength * MP_spRadius_lamJ >= MP_spRadius_cell * hmax * 0.99) then !KS MODIFIED (0.99は実数の大小判定における安全係数)
       if ( jlength * MP_spRadius_lamJ * sqrt(2.d0) >= MP_spRadius_cell * hmax ) then !KS MODIFIED (2**-0.5 < sp_SinkRadius/ (jlength * MP_spRadius_lamJ) < 2**0.5)
          sp_SinkRadius = MP_spRadius_cell * hmax
          sp_Level = level
          if (get_myrank() == PRIMARY_RANK) then
             print '(A,I4,1P3E12.4)', "level, hmax, jlength, jlength/hmax: ", level, hmax, jlength, jlength/hmax
          end if
          exit
       end if
    end do

    !KS ADDED
    if (sp_Level == Undefi) then !もしLmaxが小さい等でsp_Levelが見つからなかったらstop
       if (get_myrank() == PRIMARY_RANK) &
            print '(/,A,/,A,/)', "cannot determine sp_Level.", "stopping..."
       stop
    endif
!!$    sp_Level = Lmax
!!$    hmax = maxval(CellWidth(:,Lmax))
!!$    sp_SinkRadius = sp_f_rAcc * hmax
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'sp_SinkRadius = ', sp_SinkRadius
       print *, 'sp_SinkRadius[pc]= ', sp_SinkRadius*Unit_pc
       print *, 'sp_SinkRadius[au]= ', sp_SinkRadius*Unit_au
       print *, 'sp_Level = ', sp_Level
       call flush(6)
    end if

    ! ---------------------
    ! defne sp_SofteningRadius
    ! ---------------------
    sp_SofteningRadius = sp_SinkRadius
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'sp_SofteningRadius = ', sp_SofteningRadius
       call flush(6)
    end if

    ! ---------------------
    ! define sp_MaskRadius
    ! ---------------------
    sp_MaskRadius = sp_SinkRadius * 2
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'sp_MaskRadius = ', sp_MaskRadius
       call flush(6)
    end if

#ifdef ACCRETION_BONDI
    ! ---------------
    ! read Bondi file
    ! ---------------
    call sp_readBondi
#endif !ACCRETION_BONDI

    call write_spparam

    Initialized = .true.
contains

  ! -----------------------------------------------------------------
  ! write parameter to a file
  ! -----------------------------------------------------------------
  subroutine write_spparam
    use string, only : concat, CHARLEN
    use io_util, only : readenv, wchar
    
    character(len=CHARLEN) :: fn, dir
    if (.not. readenv('DIR', dir) ) stop
    if (get_myrank() == PRIMARY_RANK) then
       fn = concat(dir, FNSPP_)
       call wchar(6,'write parameter file = '//fn)
       open(1, file=fn)

       write(1,*) sp_SinkRadius
       write(1,*) sp_Level
       write(1,*) sp_SofteningRadius
       write(1,*) sp_MaskRadius
       close(1)

    end if
  end subroutine write_spparam
    
end subroutine sp_init

! -----------------------------------------------------------------
! update subgrid disk
! -----------------------------------------------------------------
subroutine sp_update_subdisk
  type(t_spParticle),pointer :: ptr
  real(kind=DBL_KIND) :: t

  t = Time(sp_Level) !現在時刻
  ptr => Particle
  do
     if (.not. associated(ptr)) exit

     !前回のt_prevからdt_acc以上経過していたらdisk情報を更新
     if (t-ptr%t_prev > dt_acc/Unit_yr) then
        
        !disk情報を更新
        ptr%mdot_disk = ptr%dm_disk/(t-ptr%t_prev) ! 平均のmdot
        ptr%J_disk = ptr%dJ_disk                   ! Jの合計

        !次の計測に向けてリセット
        ptr%t_prev = t                             ! t_prevを更新
        ptr%dm_disk = 0d0                          ! dmをリセット
        ptr%dJ_disk = 0d0                          ! dJをリセット
     end if
     ptr => ptr%next
  end do
  
end subroutine sp_update_subdisk
  

end module sinkParticle
