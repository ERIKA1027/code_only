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
!#define CREATION_ON
!#define ADVECTION_ON
!#define MERGER_ON
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

! user parallel cal at P2P
#define PARALLEL_SP_GRAVITYP2P

#define NOT_REFINEMENT_SP

#ifdef NOT_REFINEMENT_SP
  #define MOMENTUM_CONSERVATION
#endif

! separate timestep for subcycle
!#define SEPARATE_TIME_STEP 

! use maximum level to estimate CFL condition
#define USE_SPLEVEL_SP_CFL

! softning in P2P ?
#define SOFTNING_P2P_INCLUDE YES 



! density foor in Bondi accretion
#define RHO_FLOOR 1.d0


! outflow threshold (Msun) 
#define MINIMUM_MASS_SPOUTFLOW 0.1d0


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
     !----- HF ADDED ----!
     real(kind=DBL_KIND) :: t_crt          ! create time of sink cell
     integer :: lev                        ! maximum level of grid on sink cell
     real(kind=DBL_KIND) :: h(MX:MZ)       ! cell size 
     !----- HF ADDED ----!
     type(t_spParticle),pointer :: next => null()
  end type t_spParticle
  ! Type for link list of position of new particles
  type t_spPos
     real(DBL_KIND),dimension(MX:MZ) :: r
     type(t_spPos),pointer :: next => null()
  end type t_spPos
 
  !-----------EO_added----------!
  real(kind=DBL_KIND),save :: mass_disk=0.0
  real(kind=DBL_KIND),save :: mass_subdisk=0.0
  !-----------EO_added----------!

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


#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
  real(kind=DBL_KIND), dimension(MX:MZ), save :: hh_sp
#endif

  !
  public :: sp_update, sp_write, sp_read, sp_writeLog, sp_restrictCFL, sp_refineCond, sp_getLevel, &
       sp_sinkdata2array, sp_getNparticle, sp_getRhocr, sp_newParticle, sp_gravityOfParticle, &
       sp_getSinkRadius, sp_is_inside_sink, sp_refineCond_KS !KS MODIFIED

#ifdef ADVECTION_ON
  public :: sp_advector
#endif
contains
  !-------------------------------------------------------------------------
  ! Update sink particles.
  ! Exported subroutine.
  !-------------------------------------------------------------------------
  subroutine sp_update
    use parameter
    use modelParameter

    integer, save :: ifirst = 0
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, v
    real(kind=DBL_KIND) :: a, vr, mass

    !--------------    create sink particle by hand (KS DEBUG)    -----------------!
    !if (ifirst == 0) then 
    !    a    = 2.d0*MP_Rcloud / Unit_pc
    !    mass = 10.d0*cgs_msun / Unit_m
    !    vr   = sqrt(cgs_gc*MP_Mcl*cgs_msun/(2.d0*MP_Rcloud*cgs_pc)) / (Unit_l/Unit_t)

    !    pos  = 0.d0
    !    pos(MX) = a
    !    v    = (/0.d0, vr, 0.d0/)
    !    call sp_newParticle(mass, pos, v, (/0.d0, 0.d0, 0.d0/), mdot_disk=0.d0, J_disk=(/0.d0, 0.d0, 0.d0/) )

    !    Nparticle = 1
    !    ifirst    = 1
    !endif
    ! -------------------------------------------------------------------------------


    call sp_init

#ifdef CREATION_ON
    call sp_create
#endif !CREATION_ON
    if (Nparticle == 0) return

    ! ----------------
    call get_plev
    ! ----------------


#ifndef NOT_REFINEMENT_SP
    if (LevelMax /= sp_Level) then
       print *, '*** waring : LevelMax /= sp_Level', LevelMax, sp_Level
       call flush(6)
    endif
#endif

#ifdef ACCRETION_ON
    call sp_accretion_init
  #ifdef ACCRETION_BONDI
    #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      call sp_accretionBondi( Dtime(LevelMax) )
    #else
      call sp_accretionBondi( Dtime(sp_Level) )
    #endif
  #endif !ACCRETION_BONDI
  #ifdef ACCRETION_RHOCR
   call sp_accretion
  #endif !ACCRETION_RHOCR
    call sp_accretion_comm
    call sp_accretionUpdate
#endif !ACCRETION_ON


#ifdef ADVECTION_ON
  #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
    call sp_advector( Dtime(LevelMax) ) ! also make CFL condition for restricted timestep of HD.
  #else
    call sp_advector( Dtime(sp_Level) ) ! also make CFL condition for restricted timestep of HD.
  #endif
#endif !ADVECTION_ON
#ifdef MERGER_ON
    call sp_merge
#endif !MERGER_ON

    !----- KS ADDED ----!
    call sp_update_subdisk
    !----- KS ADDED ----!

#ifdef OUTFLOW_ON
  #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
    call sp_outflow( Dtime(LevelMax) ) ! also make CFL condition for restricted timestep of HD.
  #else
    call sp_outflow( Dtime(sp_Level) ) ! also make CFL condition for restricted timestep of HD.
  #endif
#endif

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

#ifdef OUTPUT_SP_CRTEPOCH
    call output_spcrt(list_newPos)
#endif

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
       call sp_newParticle( 0.d0, ptr%r, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), t_prev=Time(sp_Level) &
         , t_crt=Time(sp_Level) ) !KS MODIFIED
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

       !if ( dtLocal < dtGlobal ) then
       !   if (.not. bool_subcycle_on) then ! Initialize for subcycling.
       !      bool_subcycle_on = .true.
       !      call sp_gravityF2P_all(gravFluid, dtGlobal) ! dt is dtGlobal. 
       !      ! Gravity of gas, gravFluid, is obtained here and it is used through the subcycles.
       !      ! Particles affects gas during dt = dtGlobal.
       !      ! When subcycling, here after, particles never affect gas.
       !      if (get_myrank() == PRIMARY_RANK) print *, 'Sink particle: subcycle ON'
       !   endif
       !endif

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
    integer :: np

    !allocate(gravF(MX:MZ, Nparticle))

    if(Nparticle > 0) then
      call sp_gravityFluid2particle(gravFluid, dt) 
    endif

    !gravFluid(:,:) = gravF(:,:)

    !deallocate(gravF)

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
    real(kind=DBL_KIND),dimension(:,:),allocatable :: grav, gravP, gravF
#ifdef DM_NFW_PROFILE
    real(kind=DBL_KIND),dimension(:,:),allocatable :: gravDM
#endif
    integer,dimension(MX:MZ) :: mlist
    integer :: np
    myrank = get_myrank()
    mlist = (/MGX, MGY, MGZ/)
    sp_DtimeCFL = HUGE(1.d0)
    allocate(grav(MX:MZ, Nparticle), gravP(MX:MZ, Nparticle), gravF(MX:MZ, Nparticle))
#ifdef DM_NFW_PROFILE
    allocate(gravDM(MX:MZ, Nparticle))
#endif

    ! -------------------
    ! predictor step (dr)
    ! -------------------
    call sp_gravityP2P_chi(gravP) ! gravity of particle

    if ( bool_subcycle_on ) then
      gravF = gravFluid
    else
      call sp_gravityFluid2particle(gravF, dt*0.5d0)
    endif

#ifdef DM_NFW_PROFILE
    call sp_gravityOfDM(gravDM, dt*0.5d0) ! get gravity of DM
#endif

#ifdef DM_NFW_PROFILE
    grav = gravP + gravF + gravDM
#else
    grav = gravP + gravF
#endif

    ! --------------------
    ! update r temporally
    ! --------------------
    np = 0
    ptr => Particle
    do 
      if (.not. associated(ptr)) exit
      np = np + 1

      ptr%dr(:) = ptr%v(:)*dt + grav(:, np) *dt**2 * 0.5d0

      ! -------------------
      ! CFL condition
      ! -------------------

#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      sp_DtimeCFL = min( sp_DtimeCFL, &
           minval( ptr%h(:) / abs(ptr%v(:) + grav(:, np)*dt)) &
           )
#else
      sp_DtimeCFL = min( sp_DtimeCFL, &
           minval( CellWidth(:,sp_Level) /abs(ptr%v(:) + grav(:, np)*dt)) &
           )
#endif
      ptr%r = ptr%r + ptr%dr
      call get_plev_kobetu(ptr)


      ptr => ptr%next
    enddo


    ! -------------------
    ! corrector step (dv)
    ! -------------------
    call sp_gravityP2P_chi(gravP) ! gravity of particle

    if ( bool_subcycle_on ) then
      gravF = gravFluid
    else
      call sp_gravityFluid2particle(gravF, dt*0.5d0)
    endif

#ifdef DM_NFW_PROFILE
    call sp_gravityOfDM(gravDM, dt*0.5d0) ! get gravity of DM
#endif

#ifdef DM_NFW_PROFILE
    grav = (grav + gravP + gravF + gravDM)*0.5d0
#else
    grav = (grav + gravP + gravF)*0.5d0
#endif

    np = 0
    ptr => Particle
    do 
      if (.not. associated(ptr)) exit
      np = np + 1

      ptr%dv(:) = grav(:, np) *dt

      ! -------------------
      ! CFL condition
      ! -------------------

#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      sp_DtimeCFL = min( sp_DtimeCFL, &
           minval( ptr%h(:) /abs(ptr%v(:) + ptr%dv(:))) &
           )
#else
      sp_DtimeCFL = min( sp_DtimeCFL, &
           minval( CellWidth(:,sp_Level) /abs(ptr%v(:) + ptr%dv(:))) &
           )
#endif

       ptr%r = ptr%r - ptr%dr
       call get_plev_kobetu(ptr)

      ptr => ptr%next
    enddo

    deallocate( grav, gravP, gravF )
#ifdef DM_NFW_PROFILE
    deallocate(gravDM)
#endif
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
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
       h = minval(ptri%h(:))
#else
       h = minval(CellWidth(:, sp_Level))
#endif

#ifdef SEPARATE_TIME_STEP 
       dr2min = max(dr2min, sp_SofteningRadius**2.d0)
#endif
       !dtsub = min(dtsub, sqrt(min(sqrt(dr2min), h)/grav) )
       dtsub = min(dtsub, sqrt(h/grav) ) ! causion !!!! HF modified

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
#if SOFTNING_P2P_INCLUDE == YES 
       call sp_gravityOfParticle( dr, ptr%mass, gx, gy, gz)
#else
       call sp_gravityOfParticle( dr, ptr%mass, gx, gy, gz, no_softening=.true.)
#endif
       grav(MX) = grav(MX) + gx
       grav(MY) = grav(MY) + gy
       grav(MZ) = grav(MZ) + gz
       ptr => ptr%next
    end do
  end subroutine sp_gravityP2P

  !-------------------------------------------------------------------------
  ! gravity from particles to a particle chi.ver
  !-------------------------------------------------------------------------

#ifdef PARALLEL_SP_GRAVITYP2P
  subroutine sp_gravityP2P_chi(grav)
    type(t_spParticle),pointer :: ptri
    real(kind=DBL_KIND),dimension(MX:, :),intent(OUT) :: grav
    real(kind=DBL_KIND),dimension(:,:),allocatable :: grav_r
    type(t_spParticle),pointer :: ptr1, ptr2
    real(kind=DBL_KIND) :: gx, gy, gz
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr
    integer :: np, next_np

    allocate(grav_r(lbound(grav,1):ubound(grav,1), lbound(grav,2):ubound(grav,2)))

    grav   = 0.d0
    grav_r = 0.d0
    np   = 0
    ptr1 => Particle
    next_np = get_myrank() + 1

    do
      if (.not. associated(ptr1)) exit
      np = np + 1

      if(np == next_np) then

        ptr2 => Particle
        do
           if (.not. associated(ptr2)) exit
           if (ptr2%pid == ptr1%pid) then
              ptr2 => ptr2%next
              cycle
           end if
           dr = ptr1%r - ptr2%r
#if SOFTNING_P2P_INCLUDE == YES 
           call sp_gravityOfParticle( dr, ptr2%mass, gx, gy, gz)
#else
           call sp_gravityOfParticle( dr, ptr2%mass, gx, gy, gz, no_softening=.true.)
#endif
           grav_r(MX,np) = grav_r(MX,np) + gx
           grav_r(MY,np) = grav_r(MY,np) + gy
           grav_r(MZ,np) = grav_r(MZ,np) + gz
           ptr2 => ptr2%next
        end do

        ! next_np
        next_np = next_np + NPE 
      endif

      ptr1 => ptr1%next
    enddo

    call mpi_allreduce( grav_r, grav, size(grav_r), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    deallocate(grav_r)
  end subroutine sp_gravityP2P_chi

#else

  subroutine sp_gravityP2P_chi(grav)
    type(t_spParticle),pointer :: ptri
    real(kind=DBL_KIND),dimension(MX:, :),intent(OUT) :: grav
    type(t_spParticle),pointer :: ptr1, ptr2
    real(kind=DBL_KIND) :: gx, gy, gz
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr
    integer :: np

    grav = 0
    np   = 0
    ptr1 => Particle

    do
      if (.not. associated(ptr1)) exit
      np = np + 1

      ptr2 => Particle
      do
         if (.not. associated(ptr2)) exit
         if (ptr2%pid == ptr1%pid) then
            ptr2 => ptr2%next
            cycle
         end if
         dr = ptr1%r - ptr2%r
         call sp_gravityOfParticle( dr, ptr2%mass, gx, gy, gz, no_softening=.true.)
         grav(MX,np) = grav(MX,np) + gx
         grav(MY,np) = grav(MY,np) + gy
         grav(MZ,np) = grav(MZ,np) + gz
         ptr2 => ptr2%next
      end do

      ptr1 => ptr1%next
    enddo
  end subroutine sp_gravityP2P_chi

#endif


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

    !----------EO_added----------!
    real(kind=DBL_KIND) :: r_diskout, v_diskr, t_vis, real_t, t_free, v_free
    real(kind=DBL_KIND) :: dmass_sink, dmass_bh, dmass_disk, dotm_sink
    logical :: isNotFinite
    !----------EO_added----------!


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

       !----------EO_added----------!
      ! r_diskout = 4.615991693d-4 * cgs_pc
      ! v_diskr   = 0.01 * sqrt(cgs_gc*cgs_msun*1.d4 / r_diskout) / Unit_v
      ! r_diskout = r_diskout/Unit_l
      ! t_vis     = r_diskout/v_diskr
      ! v_free    = sqrt(2.0*cgs_gc*1.d4*cgs_msun/(sp_SinkRadius*Unit_l))/Unit_v
      ! t_free    = (sp_SinkRadius-r_diskout)/v_free 

      ! #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      ! real_t = Time(LevelMax) !現在時刻
      ! #else
      ! real_t = Time(sp_Level) !現在時刻
      ! #endif




!       dmass_sink = ptr%dmass

!       #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
!       mass_subdisk = (mass_subdisk + dmass_sink) / (1.0+Dtime(LevelMax)/t_free)
!       dmass_disk   = mass_subdisk / t_free * Dtime(LevelMax)
!       mass_disk    = (mass_disk + dmass_disk) / (1.0 + Dtime(LevelMax)/t_vis)
!       dmass_bh   = mass_disk/t_vis*Dtime(LevelMax)
!       #else
!       mass_subdisk = (mass_subdisk + dmass_sink) / (1.0+Dtime(sp_Level)/t_free)
!       dmass_disk   = mass_subdisk / t_free * Dtime(sp_Level)
!       mass_disk    = (mass_disk + dmass_disk) / (1.0 + Dtime(sp_Level)/t_vis)
!       dmass_bh   = mass_disk/t_vis*Dtime(sp_Level)
!       #endif

!       ptr%dmass  = max(dmass_bh,0.0)
!       buf(1, n)  = ptr%dmass


       !dmass_sink = ptr%dmass
      ! #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      ! dmass_bh   = mass_disk / t_vis*Dtime(LevelMax) 
      ! mass_disk  = mass_disk + dmass_sink - dmass_bh
      ! #else
      ! dmass_bh   = mass_disk / t_vis*Dtime(sp_Level) 
      ! mass_disk  = mass_disk + dmass_sink - dmass_bh
      ! #endif

       !ptr%dmass  = max(dmass_bh,0.0)
       !buf(1, n)  = ptr%dmass




       !----------EO_added----------!

       !----------EO_added----------!
       ! print*, "dmass_bh=", dmass_bh, t_vis, mass_disk, Dtime(LevelMax)
       !----------EO_added----------!


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
#ifdef OUTFLOW_ON
    real(kind=DBL_KIND) :: acc_rate
#endif
    ptr => Particle

#ifdef NOT_REFINEMENT_SP
    if (LevelMax < sp_Level) return
    level = sp_Level
#else
    level = LevelMax
#endif

    do
       ! -------------------------------
       if (.not. associated(ptr)) exit
       ! -------------------------------
#ifdef OUTFLOW_ON
       if(ptr%mass*Unit_msun > MINIMUM_MASS_SPOUTFLOW) then
          acc_rate = 1.d0/(1.d0+OUTFLOW_FW)
       else
          acc_rate = 1.d0
       endif
#endif

       if (ptr%lev == sp_Level) then
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

#ifdef OUTFLOW_ON
                            ! particle gets mass and momentum
                            ptr%dmass = ptr%dmass +  rhoacc*dv*acc_rate
                            ptr%dp(MX) = ptr%dp(MX) + vx(i,j,k)*rhoacc*dv  ! outflowのtotal momentumは0
                            ptr%dp(MY) = ptr%dp(MY) + vy(i,j,k)*rhoacc*dv
                            ptr%dp(MZ) = ptr%dp(MZ) + vz(i,j,k)*rhoacc*dv
                            dr = (/x(i), y(j), z(k)/) - ptr%r
                            du = (/vx(i,j,k), vy(i,j,k), vz(i,j,k)/) - ptr%v
                            ptr%ds(MX) = ptr%ds(MX) + (dr(MY)*du(MZ)-dr(MZ)*du(MY))*rhoacc*dv ! angular momentumも持ち出さないとする
                            ptr%ds(MY) = ptr%ds(MY) + (dr(MZ)*du(MX)-dr(MX)*du(MZ))*rhoacc*dv
                            ptr%ds(MZ) = ptr%ds(MZ) + (dr(MX)*du(MY)-dr(MY)*du(MX))*rhoacc*dv
#else
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


                          
#endif




                            ! mass and momentum of gas is reduced by accretion
                            ! note that velocity of gas keeps remains by accretion
                            !--------------------------- KS ADDED ------------------------!
                            pr(i,j,k) = pr(i,j,k) * (rho(i,j,k)-rhoacc)/rho(i,j,k) !密度と同じ比率で減少させる
                            !-------------------------------------------------------------!
                            rho(i,j,k) = rho(i,j,k) - rhoacc

                            !----------EO_added----------!
                          !  vx(i,j,k) = 0.d0
                          !  vy(i,j,k) = 0.d0
                          !  vz(i,j,k) = 0.d0

                           !---sink most outer grid inflow??---!
                         !  if (kg==ijkgL(MZ) .or. kg ==ijkgR(MZ) .or. jg==ijkgL(MY) .or. jg==ijkgR(MY) .or. ig==ijkgL(MX) .or. ig==ijkgR(MX)) then
                         !     if (x(i) > 0.d0) then
                         !        if (vx(i,j,k) > 0.d0) then
                         !           vx(i,j,k) = -vx(i,j,k)
                         !        end if
                         !     else
                         !        if (vx(i,j,k) < 0.d0) then
                         !           vx(i,j,k) = -vx(i,j,k)
                         !        end if
                         !     end if
 
                          !    if (y(j) > 0.d0) then
                          !       if (vy(i,j,k) > 0.d0) then
                          !          vy(i,j,k) = -vy(i,j,k)
                          !       end if
                          !    else
                          !       if (vy(i,j,k) < 0.d0) then
                          !          vy(i,j,k) = -vy(i,j,k)
                          !       end if
                          !    end if
 
                          !    if (z(k) > 0.d0) then
                          !       if (vz(i,j,k) > 0.d0) then
                          !          vz(i,j,k) = -vz(i,j,k)
                          !       end if
                          !    else
                          !       if (vz(i,j,k) < 0.d0) then
                          !          vz(i,j,k) = -vz(i,j,k)
                          !       end if
                          !    end if
           
                          !else
                          !    vx(i,j,k) = 0.d0
                          !    vy(i,j,k) = 0.d0
                          !    vz(i,j,k) = 0.d0
                          !end if
                           !----------EO_added----------!

                         end if
                      end do  ! i = Imin, Imax
                   end do ! j = Jmin, Jmax
                end do ! k = Kmin, Kmax
              end do ! ig = ijkgL(MX), ijkgR(MX)
            end do ! jg = ijkgL(MY), ijkgR(MY)
         end do ! kg = ijkgL(MZ), ijkgR(MZ)
       endif ! if (ptr%lev == sp_Level) then

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
#ifdef NOT_REFINEMENT_SP
    if (LevelMax < sp_Level) return
    level = sp_Level
#else
    level = LevelMax
#endif
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
    t_prev, dm_disk, dJ_disk, mdot_disk, J_disk, t_crt)
    !----------------------------- KS MODIFIED --------------------------!    
    real(kind=DBL_KIND),intent(IN) :: mass
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: r, v, s
    real(kind=DBL_KIND),intent(IN),optional :: dmass
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN),optional :: dr, dv, dp, ds
    integer,intent(IN),optional :: pid
    !----------------------------- KS ADDED --------------------------!
    real(kind=DBL_KIND),intent(IN),optional :: t_prev, dm_disk, mdot_disk, t_crt
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
    if (present(t_crt)) then
       newParticle%t_crt = t_crt
    else
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
       newParticle%t_crt = Time(LevelMax)
#else
       newParticle%t_crt = Time(sp_Level)
#endif
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
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
      write(LUN, '(I14, E17.9, I7, 10(1PE17.9))') Step(LevelMax), Time(LevelMax), ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%s
#else
      write(LUN, '(I14, E17.9, I7, 10(1PE17.9))') Step(sp_Level), Time(sp_Level), ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%s
#endif
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

          ! 形成時間は古い方
          ptri%t_crt = min(ptri%t_crt, ptrj%t_crt)

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
            ptr%t_prev, ptr%dm_disk, ptr%dJ_disk, ptr%mdot_disk, ptr%J_disk, ptr%t_crt
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
    
    !----------------------------- HF ADDED --------------------------------!
    real(kind=DBL_KIND) :: t_crt
    !-----------------------------------------------------------------------!
    
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
    !allocate( buf(1+1+3*7+1*3+3*2, nptcle ) ) 
    allocate( buf(1+1+3*7+1*3+3*2+1, nptcle ) ) !HF modified
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

       !----------------------------- HF ADDED --------------------------!
       t_crt  = buf(33,n)
       !----------------------------- HF ADDED --------------------------!

       !----------------------------- KS MODIFIED --------------------------!    
       ! call sp_newParticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds )
       !call sp_newparticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds, &
       !     t_prev=t_prev, dm_disk=dm_disk, dj_disk=dj_disk, mdot_disk=mdot_disk, j_disk=j_disk)
       call sp_newparticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds, &
            t_prev=t_prev, dm_disk=dm_disk, dj_disk=dj_disk, mdot_disk=mdot_disk,            &
            j_disk=j_disk, t_crt=t_crt)
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
    integer :: slevel

#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
    slevel = LevelMax
#else
    slevel = sp_Level
#endif

    if (mod(Step(slevel), skip) /= 0) return

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
       write(LUN, '(I14, 2(1PE17.9), I7, 24(1PE23.15))') Step(slevel), Time(slevel), Dtime(slevel), &
            ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds,Time(slevel)-ptr%t_crt

       !----------------------------- KS ADDED --------------------------!
       !file_each = concat(dir,prefix_log)                   ! /dir/spLog
       !file_each = concat(file_each,num2char(ptr%pid))          ! /dir/spLog10
       !file_each = concat(file_each,'.dat')                 ! /dir/spLog10.dat
       !open(LUN2, file=file_each, position='APPEND')

       !!mdot_diskがゼロで無ければ
       !if (ptr%mdot_disk > tiny(ptr%mdot_disk)) then
       !   write(LUN2, '(I10, 14(1PE15.7))') Step(slevel), Time(slevel)*Unit_yr, &
       !        ptr%mass*Unit_msun, ptr%dmass/Dtime(slevel)*Unit_msun/Unit_yr,&
       !        ptr%r*Unit_au, ptr%v*Unit_kms, ptr%t_prev*Unit_yr, &
       !        ptr%mdot_disk*Unit_msun/Unit_yr, ptr%J_disk*Unit_msun*Unit_kms*Unit_au !1*5 + 3*3 = 14 floats

       !!mdot_diskがゼロのとき
       !else
       !   !プロットの際に軸を決めるときに困るので、mdot_diskがゼロでJ_diskもゼロベクトルになる時はシンク粒子のスピンを渡す
       !   write(LUN2, '(I10, 14(1PE15.7))') Step(slevel), Time(slevel)*Unit_yr, &
       !        ptr%mass*Unit_msun, ptr%dmass/Dtime(slevel)*Unit_msun/Unit_yr,&
       !        ptr%r*Unit_au, ptr%v*Unit_kms, ptr%t_prev*Unit_yr, &
       !        ptr%mdot_disk*Unit_msun/Unit_yr, ptr%s !1*5 + 3*3 = 14 floats
       !end if
       !call flush(LUN2)
       !close(LUN2)
       !----------------------------- KS ADDED (END) --------------------------!

       !----------------------------- KS ADDED2 --------------------------!
       !spLogをやめてlogspに移行する予定
       file_each = concat(dir,prefix_log2)                   ! /dir/logsp
       file_each = concat(file_each,num2char(ptr%pid))          ! /dir/logsp10
       file_each = concat(file_each,'.dat')                 ! /dir/logsp10.dat
       open(LUN2, file=file_each, position='APPEND')

       write(LUN2, '(I10, 18(1PE15.7))') Step(slevel), Time(slevel)*Unit_yr, &
               ptr%mass*Unit_msun, ptr%dmass/Dtime(slevel)*Unit_msun/Unit_yr,&
               ptr%r*Unit_au, ptr%v*Unit_kms, ptr%s, &
               ptr%t_prev*Unit_yr, ptr%mdot_disk*Unit_msun/Unit_yr, &
               ptr%J_disk*Unit_msun*Unit_kms*Unit_au,(Time(slevel)-ptr%t_crt)*Unit_yr  !1*5 + 3*4+1 = 18 floats
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
  ! gravity of particl to fluide
  !-------------------------------------------------------------------------
  subroutine sp_gravityOfParticle2fluid(dr, mass, gx, gy, gz, soft_length, psi, no_softening)
    use modelParameter, only : MP_Gconst
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: dr
    real(kind=DBL_KIND),intent(IN) :: mass, soft_length
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

    sradii = soft_length
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
  end subroutine sp_gravityOfParticle2fluid

  !-------------------------------------------------------------------------
  ! gravity of particle. This subroutine should be inline-expanded.
  !-------------------------------------------------------------------------
  subroutine sp_gravityOfParticle_out_chi(dr, mass, gx, gy, gz)
    use modelParameter, only : MP_Gconst
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: dr
    real(kind=DBL_KIND),intent(IN) :: mass
    real(kind=DBL_KIND),intent(OUT) :: gx, gy, gz
    real(kind=DBL_KIND) :: gabs, gm, drabs2

    gm = MP_Gconst * mass
    drabs2 = sum(dr**2)
    gabs = -gm/(drabs2*sqrt(drabs2))
    
    gx = gabs * dr(MX)
    gy = gabs * dr(MY)
    gz = gabs * dr(MZ)

  end subroutine sp_gravityOfParticle_out_chi
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
  subroutine sp_gravityFluid2particle(gravF, dt)
    use mpilib
    real(kind=DBL_KIND),dimension(MX:,:),intent(OUT) :: gravF
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(:,:), allocatable :: gravBlock
    integer :: level, n, gid

    allocate(gravBlock(MX:MZ,Nparticle))

    gravF = 0.d0

#ifdef FLUID_GRAVITY_OFF
    return
#else !FLUID_GRAVITY_OFF

    myrank = get_myrank()
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )

          gid = GidList(n, level)
          call sp_gravityFluid2particleForBlck(gravBlock, dt, gid)


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
    call sp_gravityOfFluid_chi(gravF) ! get gravity of gas acting on ptr.
#endif !MOMENTUM_CONSERVATION

#endif !FLUID_GRAVITY_OFF


    deallocate(gravBlock)

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
#define M_MRHO  0
#define M_MVX   1
#define M_MVY   2
#define M_MVZ   3
#define M_MP    4

  subroutine sp_gravityFluid2particleForBlck(grav, dt, gid)
    use eos, only : w2u_4, u2w_4
    use mpilib ! KS DEBUG
    use modelParameter, only : MP_Gconst
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:, :),intent(OUT) :: grav
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, w
    integer :: level, i, j, k, ic, jc, kc, ilev, ilgid, rank, levsp
    real(kind=DBL_KIND) :: dv, dvsc, xc, yc, zc, dtdvrho, gxp, gyp, gzp, nsubm1, dtdv, dtdvsc
    real(kind=DBL_KIND),dimension(MX:MZ) :: h, hh, hsc, bl, br, pl, pr, cl, cr, dr
    real(kind=DBL_KIND) :: gm, drabs2, gabs, dryz_2, drz_2, inv_dt
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    real(kind=DBL_KIND) :: soft_length
    integer,parameter :: N_SubCell = 8
    integer :: np

    logical :: isNotFinite ! defined in naninf.F90  (KS DEBUG)


    level = get_level(gid)
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    u => get_Up(gid)
    dv = get_dv(level)
    allocate( w(ARRAYSIZE4_4(u)) )
    call u2w_4(u, w, dv)

    h(:) = CellWidth(:,level)
    hh = 0.5d0*h
    bl = (/x(Imin), y(Jmin), z(Kmin)/) - hh ! left side of block
    br = (/x(Imax), y(Jmax), z(Kmax)/) + hh ! right side of block
    hsc(:) = CellWidth(:, level) / N_SubCell    ! cell width of subcell
    dvsc = dv/N_SubCell**3
    grav = 0.d0

    dtdv   = dt*dv
    dtdvsc = dt*dvsc 
    inv_dt = 1.d0/dt

    nsubm1 = (N_SubCell-1)*0.5d0

    ! -------------
    ! step for np
    ! -------------
    np = 0
    ptr => Particle
    
    do
      if (.not. associated(ptr)) exit
      np = np + 1

      ! mass -----------------
      gm = MP_Gconst* ptr%mass
      ! ----------------------

      ! softning length -
      soft_length = sp_SofteningRadius*2.d0**(sp_Level-ptr%lev)

      pl = ptr%r - soft_length
      pr = ptr%r + soft_length
      if ( max(pl(MX),bl(MX)) < min(pr(MX),br(MX)) .and. &
           max(pl(MY),bl(MY)) < min(pr(MY),br(MY)) .and. &
           max(pl(MZ),bl(MZ)) < min(pr(MZ),br(MZ)) ) then ! block is overlaped with particle radius
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  cl = (/x(i), y(j), z(k)/) - hh
                  cr = (/x(i), y(j), z(k)/) + hh
                  if ( max(pl(MX),cl(MX)) < min(pr(MX),cr(MX)) .and. &
                       max(pl(MY),cl(MY)) < min(pr(MY),cr(MY)) .and. &
                       max(pl(MZ),cl(MZ)) < min(pr(MZ),cr(MZ)) ) then ! cell is overlaped with particle radius
                     dtdvrho = dtdvsc*u(i,j,k,MRHO)
                     do kc = 0, N_SubCell-1
                        do jc = 0, N_SubCell-1
                           do ic = 0, N_SubCell-1
                              xc = x(i) + ( ic - nsubm1 ) * hsc(MX)
                              yc = y(j) + ( jc - nsubm1 ) * hsc(MY)
                              zc = z(k) + ( kc - nsubm1 ) * hsc(MZ)
                              dr = (/ xc-ptr%r(MX), yc-ptr%r(MY), zc-ptr%r(MZ) /)
                              call sp_gravityOfParticle2fluid(dr, ptr%mass, gxp, gyp, gzp, soft_length)
  #ifndef TESTPARTICLE
                              w(i,j,k,M_MVX) = w(i,j,k,M_MVX) + gxp*dtdvrho
                              w(i,j,k,M_MVY) = w(i,j,k,M_MVY) + gyp*dtdvrho
                              w(i,j,k,M_MVZ) = w(i,j,k,M_MVZ) + gzp*dtdvrho
  #ifdef MP
                              w(i,j,k,M_MP)  = w(i,j,k,M_MP) + dtdvrho*( &
                                    u(i,j,k,MVX)*gxp &
                                   +u(i,j,k,MVY)*gyp &
                                   +u(i,j,k,MVZ)*gzp )
  #endif !MP
  #endif !TESTPARTICLE
                              grav(:,np) = grav(:,np) - (/gxp, gyp, gzp/) * dvsc * u(i,j,k,MRHO)
                           end do
                        end do
                     end do
                  else         ! cell is outside particle radius
                     dtdvrho = dtdv*u(i,j,k,MRHO)
                     dr = (/ x(i)-ptr%r(MX), y(j)-ptr%r(MY), z(k)-ptr%r(MZ) /)
                     call sp_gravityOfParticle_out_chi(dr, ptr%mass, gxp, gyp, gzp)
  #ifndef TESTPARTICLE
                     w(i,j,k,M_MVX) = w(i,j,k,M_MVX) + gxp*dtdvrho
                     w(i,j,k,M_MVY) = w(i,j,k,M_MVY) + gyp*dtdvrho
                     w(i,j,k,M_MVZ) = w(i,j,k,M_MVZ) + gzp*dtdvrho
    #ifdef MP
                     w(i,j,k,M_MP)  = w(i,j,k,M_MP) + dtdvrho*( &
                           u(i,j,k,MVX)*gxp &
                          +u(i,j,k,MVY)*gyp &
                          +u(i,j,k,MVZ)*gzp )
    #endif !MP
  #endif !TESTPARTICLE
                     grav(:,np)=grav(:,np)-(/gxp, gyp, gzp/)*dv*u(i,j,k,MRHO)
                  end if
               end do
            end do
         end do
      else                     ! block is outside particle radius
         do k = Kmin, Kmax 
            do j = Jmin, Jmax 
               do i = Imin, Imax 
                  dtdvrho = dtdv*u(i,j,k,MRHO)
                  dr = (/ x(i)-ptr%r(MX), y(j)-ptr%r(MY), z(k)-ptr%r(MZ) /)
                  call sp_gravityOfParticle_out_chi(dr, ptr%mass, gxp, gyp, gzp) 
  #ifndef TESTPARTICLE
                     w(i,j,k,M_MVX) = w(i,j,k,M_MVX) + gxp*dtdvrho
                     w(i,j,k,M_MVY) = w(i,j,k,M_MVY) + gyp*dtdvrho
                     w(i,j,k,M_MVZ) = w(i,j,k,M_MVZ) + gzp*dtdvrho
  #ifdef MP
                     w(i,j,k,M_MP)  = w(i,j,k,M_MP) + dtdvrho*( &
                           u(i,j,k,MVX)*gxp &
                          +u(i,j,k,MVY)*gyp &
                          +u(i,j,k,MVZ)*gzp )
  #endif !MP
  #endif !TESTPARTICLE
                     grav(:,np) = grav(:,np) - (/gxp, gyp, gzp/) * dv * u(i,j,k,MRHO)
               end do
            end do
         end do
      end if

      if (ptr%mass > TINY(ptr%mass)) then
        grav(:,np) = grav(:,np) / ptr%mass
      else
        grav(:,np) = 0.d0
      end if

      ptr => ptr%next
  
    enddo ! do


    call w2u_4(w, u, dv)

    !NaN/Inf check for pressure (KS DEBUG)
    !do i = Imin,Imax
    !   do j = Jmin,Jmax
    !      do k = Kmin,Kmax
    !         if (isNotFinite(u(i,j,k,MP))) then
    !            print '(A,5I6, 1P3E15.7)', "(spgf2pfb) NaN/Inf found in u(i,j,k,MP)",&
    !                 get_myrank(),gid,i,j,k,x(i),y(j),z(k)
    !            !print '(A,1P21E15.7)','w: ',w(i,j,k,:)
    !            print '(A,1P21E15.7)','u: ',u(i,j,k,:)
    !            print '(A)','stopping...'
    !            stop
    !         end if
    !      end do
    !   end do
    !end do
    
    deallocate( w )
  end subroutine sp_gravityFluid2particleForBlck

#undef M_MRHO  
#undef M_MVX   
#undef M_MVY   
#undef M_MVZ   
#undef M_MP    



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
    call sp_sumUinSphere_chi(ptr%r, sp_SofteningRadius, (/MGX, MGY, MGZ/), grav, vol) 
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
  ! get gravity of fluid
  !-------------------------------------------------------------------------
  subroutine sp_gravityOfFluid_chi(grav)
    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND),dimension(MX:,:),intent(OUT) :: grav
    real(kind=DBL_KIND),dimension(MX:MZ) :: grav_ptr
    real(kind=DBL_KIND) :: vol
    real(kind=DBL_KIND) :: vol_spsph ! KS ADDED
    integer :: np

#ifndef FLUID_GRAVITY_OFF
  
    !-----------
    ! loop sp
    !----------
    np  = 0
    ptr => Particle

    do
      if (.not. associated(ptr)) exit
      np = np + 1

      call sp_sumUinSphere_chi(ptr%r, sp_SofteningRadius, (/MGX, MGY, MGZ/), grav_ptr, vol)
      grav(:,np) = grav_ptr(:)/vol
      !------------------------- KS ADDED -----------------------------!
      ! check whether vol is consistent with sink volume
      vol_spsph = 4.*3.14/3.*sp_SofteningRadius**3 ! volume of sink particle sphere 
      if (abs((vol-vol_spsph)/vol_spsph) > 0.1) then
         if (get_myrank() == PRIMARY_RANK) &
              print '(A,1P10E15.7)', "(sp_gravityOfFluid) ** WARNING ** sp volume is not consistent: ", &
              ptr%r, vol, vol_spsph
      end if
      !----------------------------------------------------------------!
      ptr => ptr%next
    enddo
    
#else !FLUID_GRAVITY_OFF
    grav = 0.d0
#endif !FLUID_GRAVITY_OFF
  end subroutine sp_gravityOfFluid_chi
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
  ! get sum value of u within radius chi
  !-------------------------------------------------------------------------

  subroutine sp_sumUinSphere_chi(point, radius, mlist, usum, vol)
    use modelParameter, only : MP_spRadius_cell
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: point
    real(kind=DBL_KIND),intent(IN) :: radius
    integer,intent(IN),dimension(:) :: mlist
    real(kind=DBL_KIND),dimension(size(mlist)),intent(OUT) :: usum
    real(kind=DBL_KIND),intent(OUT) :: vol
    integer,parameter :: N_SubCell = 8
    integer :: level, gid, rank, i,j,k, m, ic, jc, kc, ig, jg, kg
    integer :: N_spsub
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR, posSP, posbt, bl, br, psb
    integer,dimension(MX:MZ) :: ijkg, ijkgL, ijkgR
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: u
    real(kind=DBL_KIND) :: dv, hsc(MX:MZ), xc, yc, zc
    real(kind=DBL_KIND) :: h(MX:MZ), hspace, radius2
    real(kind=DBL_KIND),dimension(size(mlist)+1) :: buf, bufr
    real(kind=DBL_KIND) :: ti, uj, vk, tim, ujm, vkm, ulin
    integer :: srank, slevel, sgid, splevel_max

    radius2 = radius*radius

    hsc(:) = CellWidth(:, Lmin) / (N_SubCell*2**(sp_Level-Lmin))    ! cell width of subcell
    dv = hsc(MX)*hsc(MY)*hsc(MZ) 

    usum(:) = 0.d0
    vol = 0.d0
    ! evaluate blocks
    posL(:) = point(:) - radius
    posR(:) = point(:) + radius
#ifdef PARTICLECREATION_IN_BOUNDARY
    call sp_restrict_within_boundary(posL, posR)
#endif !PARTICLECREATION_IN_BOUNDARY

    N_spsub = 2*MP_spRadius_cell*N_Subcell

    posSP(MX)=hsc(MX)*dble(nint(point(MX)/hsc(MX)))
    posSP(MY)=hsc(MY)*dble(nint(point(MY)/hsc(MY)))
    posSP(MZ)=hsc(MZ)*dble(nint(point(MZ)/hsc(MZ)))

    ! get center cell
    posbt(MX)=hsc(MX)*(dble(nint(point(MX)/hsc(MX)))-dble(N_spsub)*0.5d0-0.5d0)
    posbt(MY)=hsc(MY)*(dble(nint(point(MY)/hsc(MY)))-dble(N_spsub)*0.5d0-0.5d0)
    posbt(MZ)=hsc(MZ)*(dble(nint(point(MZ)/hsc(MZ)))-dble(N_spsub)*0.5d0-0.5d0)

    !print *, "x coordinate of sink", hsc(MX)*dble(nint(point(MX)/hsc(MX))), point(MX), hsc(MX)
    !print *, "check radius spacing", 2.d0*radius/dble(N_spsub), hsc(MX)
    !print *, "positions of bt", posbt(MX), hsc(MX)*dble(nint(point(MX)/hsc(MX)))-radius-0.5d0*hsc(MX)
    !print *, "LevelMax, Lmin", LevelMax, Lmin


    
    do kc = 1, N_spsub
      zc = posbt(MZ)+kc*hsc(MZ)
      do jc = 1, N_spsub
        yc = posbt(MY)+jc*hsc(MY)
        do ic = 1, N_spsub    
          xc = posbt(MX)+ic*hsc(MX)

          if ( (xc-posSP(MX))**2 + (yc-posSP(MY))**2 + (zc-posSP(MZ))**2 <= radius2 ) then

             psb =  (/xc, yc, zc/)

             ! ----------
             sgid  = Undefi
             srank = Undefi
             splevel_max = 0
             ! ----------

             do level = LevelMax, Lmin, -1 ! sinkが存在する点に該当する最深のgridを探査
                 
                 call ob_getIjkgridFromCoordPhys(ijkg, level, psb)
                 call get_gid_from_ijkgrid(ijkg(MX),ijkg(MY),ijkg(MZ),level,gid,rank)
                 
                 if (gid /= Undefi) then
                   sgid   = gid
                   slevel = level
                   srank  = rank
                   splevel_max = max(splevel_max, slevel)
                   exit ! 探査終わり  
                 endif
             enddo
             
             ! ------------------------------------------------------
             !if (srank == get_myrank()) then
             !    x => get_Xp(sgid)
             !    y => get_Yp(sgid)
             !    z => get_Zp(sgid)
             !    h = CellWidth(:,slevel)

             !    ! point of cell
             !    i = int((psb(MX)-x(Imin))/h(MX) + 0.5d0)+Imin
             !    j = int((psb(MY)-y(Jmin))/h(MY) + 0.5d0)+Jmin
             !    k = int((psb(MZ)-z(Kmin))/h(MZ) + 0.5d0)+Kmin
             !    do m = lbound(mlist, 1), ubound(mlist, 1)
             !        u => get_Ucomp(mlist(m), sgid)
             !        usum(m) = usum(m) + u(i,j,k)*dv
             !        if (m == 1) vol = vol + dv
             !    enddo
             !endif
             ! ------------------------------------------------------

             ! ------------------------------------------------------
             if (srank == get_myrank()) then
                 x => get_Xp(sgid)
                 y => get_Yp(sgid)
                 z => get_Zp(sgid)
                 h = CellWidth(:,slevel)


                 ! point of cell
                 i = int((psb(MX)-x(Imingh))/h(MX))+Imingh
                 j = int((psb(MY)-y(Jmingh))/h(MY))+Jmingh
                 k = int((psb(MZ)-z(Kmingh))/h(MZ))+Kmingh

                 ti = (psb(MX) - x(i))/(x(i+1) - x(i))
                 uj = (psb(MY) - y(j))/(y(j+1) - y(j))
                 vk = (psb(MZ) - z(k))/(z(k+1) - z(k))

                 tim = 1.d0 - ti
                 ujm = 1.d0 - uj
                 vkm = 1.d0 - vk

                 !print *, "ti, uj, vk", ti, uj, vk, i, j, k

                 do m = lbound(mlist, 1), ubound(mlist, 1)
                     u => get_Ucomp(mlist(m), sgid)

                     ulin = tim*ujm*vkm*u(i,j,k)    + ti *ujm*vkm*u(i+1,j,k)   &
                          + tim*uj *vkm*u(i,j+1,k)  + ti *uj *vkm*u(i+1,j+1,k) &
                          + tim*ujm*vk *u(i,j,k+1)  + ti *ujm*vk *u(i+1,j,k+1) &
                          + tim*uj *vk *u(i,j+1,k+1)+ ti *uj *vk *u(i+1,j+1,k+1)

                     usum(m) = usum(m) + ulin*dv
                     if (m == 1) vol = vol + dv
                 enddo
             endif
             ! ------------------------------------------------------

          endif

        enddo
      enddo
    enddo


    buf(1) = vol
    buf(2:size(mlist)+1) = usum(:)
    call mpi_allreduce( buf, bufr, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    vol = buf(1)
    usum(:) = buf(2:size(mlist)+1)


    !if(get_myrank() == 0 ) print *, "sinkpos:", splevel_max, point(MX), point(MY), point(MZ) &
    !  , (sqrt(sum(point**2))-2.d0*MP_Rcloud/Unit_pc)/(2.d0*MP_Rcloud/Unit_pc)
  end subroutine sp_sumUinSphere_chi













  !-------------------------------------------------------------------------
  ! update accretion
  !-------------------------------------------------------------------------
  subroutine sp_accretionUpdate
    use modelParameter, only : MP_CONNECTION_RUN ! KS ADDED

    !----------EO_added----------!
    real(kind=DBL_KIND) :: r_diskout, v_diskr, t_vis, real_t, t_free, v_free
    real(kind=DBL_KIND) :: dmass_sink, dmass_bh, dmass_disk, dotm_sink
    logical :: isNotFinite
    !----------EO_added----------!



    type(t_spParticle),pointer :: ptr
    real(kind=DBL_KIND) :: mass_new
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ! accretion is updated (ptr%dp, ptr%dmass, ptr%ds)


       !----------EO_added----------!
      ! r_diskout = 4.3d-5 * cgs_pc
       r_diskout = 4.615991693d-4 * cgs_pc
       v_diskr   = 0.01 * sqrt(cgs_gc*cgs_msun*1.d4 / r_diskout) / Unit_v
       r_diskout = r_diskout/Unit_l
       t_vis     = r_diskout/v_diskr
       v_free    = sqrt(2.0*cgs_gc*1.d4*cgs_msun/(sp_SinkRadius*Unit_l))/Unit_v
      ! v_free    = 0.01*2.d6 / Unit_v
       t_free    = (sp_SinkRadius-r_diskout)/v_free 

       #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
       real_t = Time(LevelMax) !現在時刻
       #else
       real_t = Time(sp_Level) !現在時刻
       #endif

      ! if(real_t<=t_vis) then
      !   ptr%dmass  = 1.d-8
      ! else
       dmass_sink = ptr%dmass
       #if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
       dmass_bh   = mass_disk / t_vis*Dtime(LevelMax) 
       mass_disk  = mass_disk + dmass_sink - dmass_bh
       #else
       dmass_bh   = mass_disk / t_vis*Dtime(sp_Level) 
       mass_disk  = mass_disk + dmass_sink - dmass_bh
       #endif

       ptr%dmass  = max(dmass_bh,0.0)
      ! end if
       !----------EO_added----------!
 



       mass_new = ptr%mass + ptr%dmass
       ptr%v = (ptr%v * ptr%mass + ptr%dp) / mass_new
       ptr%s = (ptr%s * ptr%mass + ptr%ds) / mass_new
       ptr%mass = mass_new

       !----- KS ADDED ----!
       if (MP_CONNECTION_RUN == 0) then       !通常の場合
          ptr%dm_disk = ptr%dm_disk + ptr%dmass
          ptr%dJ_disk = ptr%dJ_disk + ptr%ds
       else                                   !connection runの場合 (人工的に降着率が上がるのを防止)
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
          ptr%dm_disk = ptr%dm_disk + ptr%mdot_disk * Dtime(LevelMax) !円盤の降着率を仮定
#else
          ptr%dm_disk = ptr%dm_disk + ptr%mdot_disk * Dtime(sp_Level) !円盤の降着率を仮定
#endif
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
       call get_plev_kobetu(ptr)
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

#ifdef NOT_REFINEMENT_SP
    return
#endif

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

#ifdef NOT_REFINEMENT_SP
    return
#endif

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
#ifdef OUTFLOW_ON
    margin = 2.d0*sp_SinkRadius ! outflowを吹かせる場合は2倍のmarginが必要
#else
    margin = sp_SinkRadius      ! outflowを吹かせる場合は2倍のmarginが必要
#endif

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
  subroutine sp_sinkdata2array(np, pmass, pmdot, pr, pv, ps, pid, pmdot_disk, pJ_disk, pmass_prev  &
      , pt_prev, pdm_disk, ptcrt)  
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

    !----------------------------- HF ADDED --------------------------!
    real(kind=DBL_KIND),dimension(size(pmass)),intent(OUT),optional :: ptcrt
    !-----------------------------------------------------------------!

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
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
          if (Dtime(LevelMax) == 0.d0) then
             print *, "sp_sinkdata2array: Dtime(LevelMax) is zero, stopping..."
             stop
          end if
          pmdot(np) = ptr%dmass/Dtime(LevelMax)
#else
          if (Dtime(sp_Level) == 0.d0) then
             print *, "sp_sinkdata2array: Dtime(sp_Level) is zero, stopping..."
             stop
          end if
          pmdot(np) = ptr%dmass/Dtime(sp_Level)
#endif
       end if
       if (present(pmdot_disk)) pmdot_disk(np) = ptr%mdot_disk
       if (present(pJ_disk)) pJ_disk(:, np) = ptr%J_disk(:)
       if (present(pmass_prev)) pmass_prev(np) = ptr%mass - ptr%dm_disk ! diskに滞まっているmassを差っ引くとt=t_prevでのmassが得られる
       if (present(pt_prev)) pt_prev(np) = ptr%t_prev
       if (present(pdm_disk)) pdm_disk(np) = ptr%dm_disk
       if (present(ptcrt)) ptcrt(np) = ptr%t_crt
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
    use modelParameter, only : MP_spNcr,MP_spCs,MP_spRadius_cell,MP_spRadius_lamJ &
      , MP_Gconst, MP_mu, MP_Bondi_radius, MP_Lmax0, MP_Boxsize
    use eos
    real(kind=DBL_KIND) :: csp, jlength, hmax, mjeans, h0
    integer :: level
    if (Initialized) return
    if (get_myrank() == PRIMARY_RANK) then
       print *, 'initialize sink particle'
       call flush(6)
    end if

    !----------EO_added----------!
    !mass_disk = 0.0
    !----------EO_added----------!


    ! ------------------------------
    ! set sp_Rhocr and sp_RhocrCreate
    ! ------------------------------
#if MODEL_ART > 0
    sp_Rhocr       = MP_spNcr * cgs_amu * MP_mu / Unit_rho  ! KS MODIFIED
#else !MODEL_ART    
    sp_Rhocr       = MP_spNcr * cgs_mh / Unit_rho  ! KS MODIFIED
#endif !MODEL_ART    
    sp_RhocrCreate = sp_Rhocr    ! KS MODIFIED

    ! --------------------
    ! define sp_SinkRadius
    ! --------------------
    ! sound speed for gas at sp_Ncr

    !----------EO_removed----------!
   ! csp = MP_spCs/Unit_v    ! KS MODIFIED
   ! jlength = csp * sqrt(Pi / (MP_Gconst * sp_Rhocr))
   ! do level = Lmin, Lmax
   !    hmax = maxval(CellWidth(:,Lmin))/2.d0**(level-Lmin)
   !    if ( jlength >= JeansConst * hmax ) then ! jeans condition
    !----------EO_removed----------!

       !------------- KS DEBUG -------------!
       ! if (get_myrank() == PRIMARY_RANK) then
       !    print *, "KS DEBUG", level, hmax, jlength * MP_spRadius_lamJ, MP_spRadius_cell * hmax * 0.99
       ! end if
       !------------- KS DEBUG (END) ---------!
       !       if ( jlength * MP_spRadius_lamJ >= MP_spRadius_cell * hmax * 0.99) then !KS MODIFIED (0.99は実数の大小判定における安全係数)
!       if ( jlength * MP_spRadius_lamJ * sqrt(2.d0) >= MP_spRadius_cell * hmax ) then !KS MODIFIED (2**-0.5 < sp_SinkRadius/ (jlength * MP_spRadius_lamJ) < 2**0.5)

      ! if ( 1.d-1*MP_Bondi_radius >= MP_spRadius_cell * hmax ) then
      ! if ( MP_Bondi_radius / 3.d1 >= MP_spRadius_cell * hmax ) then

      !    sp_SinkRadius = MP_spRadius_cell * hmax
      !    sp_Level = level
      !    mjeans   = 4.d0/3.d0*pi*(jlength*0.5d0)**3*sp_Rhocr*Unit_msun
      !    if (get_myrank() == PRIMARY_RANK) then
      !       print '(A,I4,1P4E12.4)', "level, hmax, jlength, jlength/hmax, Mjeans: ", level, hmax, jlength, jlength/hmax, mjeans
      !    end if
      !    exit
      ! end if
   ! end do
      !----------EO_removed----------!


      !-----------EO_added----------!
      level = MP_Lmax0
      hmax = maxval(CellWidth(:,Lmin))/2.d0**(level-Lmin)
     ! h0 = 2d0*MP_Boxsize / ((NGI_BASE) * (NI))     !ベースグリッドのセルサイズ
     ! hmax = h0/(2d0**MP_Lmax0)                     !最大レベルでのセルサイズ
      sp_SinkRadius = MP_spRadius_cell * hmax
      sp_Level = MP_Lmax0
      if (get_myrank() == PRIMARY_RANK) then
         print '(A,I4,1P4E12.4)', "level, hmax, jlength, jlength/hmax, Mjeans: ", level, hmax, jlength, jlength/hmax, mjeans
      end if 
      !----------EO_added-----------!

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
       print *, 'MP_spRadius_cell = ', MP_spRadius_cell
       print *, 'LevelMax = ', LevelMax
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

    ! -----------------------                                                                                                                                             
    !  cellwidth of levelmax                                                                                                                                              
    ! -----------------------                                                                                                                                             
#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)                                                                                                                    
    hh_sp(:) = CellWidth(:, Lmin) / (2.d0**(sp_Level-Lmin))    ! cell width of subcell                                                                                    
#endif     


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
#ifdef USE_TKH_TACC
  use primordial, only: ProstFit2
  real(kind=DBL_KIND) :: mass, mdot, radius, lum, Trad, tacc, tkh, dtaccc, tage
#endif
  type(t_spParticle),pointer :: ptr
  real(kind=DBL_KIND) :: t
  logical :: isNotFinite
  


#if defined(SINGLE_STEP) && defined(NOT_REFINEMENT_SP)
  t = Time(LevelMax) !現在時刻
#else
  t = Time(sp_Level) !現在時刻
#endif

  ptr => Particle
  do
     if (.not. associated(ptr)) exit

#ifdef USE_TKH_TACC
     mass = ptr%mass*Unit_msun
     mdot = ptr%mdot_disk*Unit_msun/Unit_yr
     tage = (Time(Lmin)-ptr%t_crt)*Unit_yr   ! [yr]

     call ProstFit2(mass, mdot, tage, radius, lum, Trad)



     tacc = mass/(mdot + 1.d-99) ! [yr]
     tkh  = cgs_gc*(mass*cgs_msun)**2.d0/(radius*cgs_rsun*lum*cgs_lsun)/cgs_yr ![yr]

     dtaccc = max(min(tacc, tkh)*0.1d0, 1.d0) ! 0.1倍の時間待つことにする
     dtaccc = min(dtaccc, 1.d4)               ! 1万年に一回はチェック        


     !前回のt_prevからdt_acc以上経過していたらdisk情報を更新
     if (t-ptr%t_prev > dtaccc/Unit_yr) then
        
        !disk情報を更新
        ptr%mdot_disk = ptr%dm_disk/(t-ptr%t_prev) ! 平均のmdot
        ptr%J_disk = ptr%dJ_disk                   ! Jの合計

        !次の計測に向けてリセット
        ptr%t_prev = t                             ! t_prevを更新
        ptr%dm_disk = 0d0                          ! dmをリセット
        ptr%dJ_disk = 0d0                          ! dJをリセット

        ! --------------------------------------------------------
        if (get_myrank() == PRIMARY_RANK) then
          if(isNotFinite(ptr%mdot_disk)) then
            print *, "NAN appear mdot:", ptr%dm_disk, t, ptr%t_prev
          endif
        endif
        ! --------------------------------------------------------

     end if

#else
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
#endif

     ptr => ptr%next
  end do
  
end subroutine sp_update_subdisk


#ifdef OUTFLOW_ON

  #define M_MRHO  0
  #define M_MVX   1
  #define M_MVY   2
  #define M_MVZ   3
  #define M_MP    4

  ! -----------------------------------------------------------------------------
  !                     Inject outflow momentum 
  ! -----------------------------------------------------------------------------
  subroutine sp_outflow(dtGlobal)
    use mpilib
    use modelParameter
    use primordial, only: ProstFit2
    use eos, only : w2u_4, u2w_4
    use parameter, only :  Pi
    real(kind=DBL_KIND), intent(IN) :: dtGlobal
    type(t_spParticle),pointer :: ptr
    integer :: kg, jg, ig, gid, rank, gg_ijk, N_sub
    integer :: i, j, k, cell_num
    integer :: ic, jc, kc
    integer,dimension(MX:MZ) :: ijkg
    integer :: sgid, srank, slevel, level, my_rank
    integer, dimension(:,:,:), allocatable :: lsrank, lslevel, lsgid
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, posbt, psb, h, r_i, jang, hsc
    real(kind=DBL_KIND) :: dv, dvsc, dist_outflw, dist_outflw2, abs_spin, rj_cross, theta, chi, xi, sqrt_theta2, del_theta
    real(kind=DBL_KIND) :: theta_p, theta_m, atan_p, atan_m, sum_eta, sum_local
    real(kind=DBL_KIND) :: wmdotdt, vkep, rstar
    real(kind=DBL_KIND) :: mass, mdot, radius, lum, Trad, dmom, c2_gam_w, gpt, gxp, gyp, gzp, dtdvdrho
    real(kind=DBL_KIND),dimension(:,:,:),allocatable :: eta
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND), parameter :: theta0 = 0.01d0
    real(kind=DBL_KIND), parameter :: T_w    = 1.d4       ! outflow temperature
    real(kind=DBL_KIND), parameter :: gam_w  = 5.d0/3.d0  
    real(kind=DBL_KIND), parameter :: xmu_w  = 6.622517d-1
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, w
    real(kind=DBL_KIND) :: xc, yc, zc
    real(kind=DBL_KIND) :: tage
    logical :: isNotFinite



    !-------------- pre loop --------------!
    hsc(:) = CellWidth(:, Lmin) / (2**(sp_Level-Lmin))    ! cell width at sp_Level
    dvsc   = hsc(MX)*hsc(MY)*hsc(MZ) 
    sqrt_theta2 = sqrt(theta0**2 + 1.d0)
    del_theta   = atan(1.d0/dble(2*MP_spRadius_cell))
    c2_gam_w = (cgs_kb*T_w/(xmu_w*cgs_mh))/Unit_v**2/(gam_w-1.d0) ![noD]
    my_rank  = get_myrank()
    !--------------------------------------!


    ! array for set weight
    cell_num = 2*MP_spRadius_cell ! number of cell
    N_sub = cell_num*2
    allocate(eta(N_sub, N_sub, N_sub))
    allocate(lsrank(N_sub, N_sub, N_sub), lslevel(N_sub, N_sub, N_sub), lsgid(N_sub, N_sub, N_sub))



    !--- pointer for sink ----!
    ptr => Particle
    !-------------------------!

  
    do 
      if (.not. associated(ptr)) exit

      ! ------------------------------
      if (ptr%lev .ne. sp_Level) then
        ptr => ptr%next         
        cycle
      endif
      ! ------------------------------

       if(ptr%mass*Unit_msun <= MINIMUM_MASS_SPOUTFLOW) then
        ptr => ptr%next         
        cycle
       endif


      ! outflow rate
      wmdotdt = OUTFLOW_FW*ptr%dmass 

      ! read table of stellar evolution
      mass = ptr%mass*Unit_msun
      mdot = ptr%mdot_disk*Unit_msun/Unit_yr
      tage = (Time(Lmin)-ptr%t_crt)*Unit_yr   ! [yr]
      call ProstFit2(mass, mdot, tage, radius, lum, Trad)
      rstar = radius*cgs_rsun                             ! convert [Rsun] => [cm]
      vkep  = sqrt(cgs_gc*ptr%mass*Unit_m/rstar)/Unit_v   ! convert [cm/s] => [noD]

      ! --------------------------------------------------------
      if (my_rank == PRIMARY_RANK) then
        if(isNotFinite(vkep)) then
          print *, "NAN appear at sp_outflowi:", mass, mdot, tage, rstar, vkep, ptr%mass
        endif
      endif
      ! --------------------------------------------------------

      !shift position to the central position
      pos(MX) = hsc(MX)*dble(nint(ptr%r(MX)/hsc(MX)))
      pos(MY) = hsc(MY)*dble(nint(ptr%r(MY)/hsc(MY)))
      pos(MZ) = hsc(MZ)*dble(nint(ptr%r(MZ)/hsc(MZ)))

      ! evaluate blocks
      posbt(:) = pos(:) - (dble(cell_num)+0.5d0)*hsc(:)

      ! spin vector --------------------------------------------------------
      abs_spin = dsqrt(ptr%s(MX)**2+ptr%s(MY)**2+ptr%s(MZ)**2)

      if(abs_spin <= 0.d0) then
        jang(MX) = 0.d0
        jang(MY) = 0.d0
        jang(MZ) = 1.d0
      else
        jang(:)  = ptr%s(:)/abs_spin
      endif
      ! --------------------------------------------------------------------

      ! reset
      eta = 0.d0

      ! loop
      do kc = 1, N_sub
         zc = posbt(MZ)+kc*hsc(MZ)
         do jc = 1, N_sub
            yc = posbt(MY)+jc*hsc(MY)
            do ic = 1, N_sub
               xc = posbt(MX)+ic*hsc(MX)

               ! distance to sink particle ----------------------------------------
               dist_outflw2 = (xc-pos(MX))**2 + (yc-pos(MY))**2 + (zc-pos(MZ))**2 
               dist_outflw  = dsqrt(dist_outflw2)
               ! ------------------------------------------------------------------

               if( dist_outflw < 2.d0*sp_SinkRadius .and. dist_outflw >= sp_SinkRadius ) then
             
                  psb =  (/xc, yc, zc/)

                  ! ----------
                  sgid  = Undefi
                  srank = Undefi
                  slevel= Undefi
                  ! ----------

                  ! -----------------------------------------------------------------------
                  do level = LevelMax, Lmin, -1 ! sinkが存在する点に該当する最深のgridを探査
               
                      call ob_getIjkgridFromCoordPhys(ijkg, level, psb)
                      call get_gid_from_ijkgrid(ijkg(MX),ijkg(MY),ijkg(MZ),level,gid,rank)
               
                      if (gid /= Undefi) then
                        sgid   = gid
                        slevel = level
                        srank  = rank
                        exit ! 探査終わり  
                      endif
                  enddo
                  ! ------------------------------------------------------------------------

                  lsgid(ic,jc,kc)   = sgid         
                  lsrank(ic,jc,kc)  = srank 
                  lslevel(ic,jc,kc) = slevel      

                  ! get eta --------------------------------------
                  if(srank == my_rank) then

                     ! vector --------------------------------
                     r_i(MX) = (xc - pos(MX))/dist_outflw
                     r_i(MY) = (yc - pos(MY))/dist_outflw
                     r_i(MZ) = (zc - pos(MZ))/dist_outflw
                     ! ---------------------------------------


                     ! obtain theta -----------------------------------------------
                     rj_cross = r_i(MX)*jang(MX)+r_i(MY)*jang(MY)+r_i(MZ)*jang(MZ)
                     rj_cross = max(min(rj_cross, 1.d0),-1.d0)
                     theta    = acos(rj_cross)
                     ! ------------------------------------------------------------
                    
                     !if(abs(sin(0.5*Pi - theta)) .ge. hsc(MX)/dist_outflw) then
                     if(abs(sin(0.5*Pi - theta)) .ge. 0.5d0) then
                       ! chi & xi -----------------------------------------
                       theta_p  = theta + 0.5d0*del_theta
                       theta_m  = theta - 0.5d0*del_theta
                       atan_p   = atan(sqrt_theta2*tan(theta_p)/theta0)
                       atan_m   = atan(sqrt_theta2*tan(theta_m)/theta0)
                       chi = 1.d0/dist_outflw2 
                       xi  = (atan_p-atan_m)/(theta0*sqrt_theta2)
                       ! --------------------------------------------------
                       eta(ic, jc, kc) = chi*xi*dvsc
                    endif


                  endif
               else
                  lsgid(ic,jc,kc)   = Undefi       
                  lsrank(ic,jc,kc)  = Undefi
                  lslevel(ic,jc,kc) = Undefi      
               endif

            enddo
         enddo
      enddo

      ! sum of kernel
      sum_local = sum(eta)
      call mpi_allreduce( sum_local, sum_eta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      eta = eta / sum_eta

      do kc = 1, N_sub
         do jc = 1, N_sub
            do ic = 1, N_sub

               if (lsrank(ic,jc,kc) /= my_rank ) cycle
               if (lsgid(ic,jc,kc) == Undefi) cycle

               ! ----------------------
               xc = posbt(MX)+ic*hsc(MX)
               yc = posbt(MY)+jc*hsc(MY)
               zc = posbt(MZ)+kc*hsc(MZ)
               psb =  (/xc, yc, zc/)
               ! ----------------------

               ! ----------------
               sgid = lsgid(ic,jc,kc) 
               x => get_Xp(sgid)
               y => get_Yp(sgid)
               z => get_Zp(sgid)
               u => get_Up(sgid)
               dv = get_dv(lslevel(ic,jc,kc))
               allocate( w(ARRAYSIZE4_4(u)) )
               call u2w_4(u, w, dv)
               ! ----------------

               ! point of cell
               h = CellWidth(:,lslevel(ic,jc,kc))
               i = int((psb(MX)-x(Imingh))/h(MX))+Imingh
               j = int((psb(MY)-y(Jmingh))/h(MY))+Jmingh
               k = int((psb(MZ)-z(Kmingh))/h(MZ))+Kmingh

               ! distance to sink particle ----------------------------------------
               dist_outflw2 = (xc-pos(MX))**2 + (yc-pos(MY))**2 + (zc-pos(MZ))**2 
               dist_outflw  = dsqrt(dist_outflw2)
               ! ------------------------------------------------------------------

               ! vector --------------------------------
               r_i(MX) = (xc - pos(MX))/dist_outflw
               r_i(MY) = (yc - pos(MY))/dist_outflw
               r_i(MZ) = (zc - pos(MZ))/dist_outflw
               ! ---------------------------------------

               ! --------------------
               ! mass injection rate
               ! --------------------
               w(i,j,k,M_MRHO) = w(i,j,k,M_MRHO) + wmdotdt*eta(ic, jc, kc)  ! [g in code unit]
               dtdvdrho = w(i,j,k,M_MRHO)*dtGlobal ! [ g s in code unit]
               
               ! --------------------
               ! momentum injection 
               ! --------------------
               dmom = OUTFLOW_FV*vkep*wmdotdt*eta(ic, jc, kc) ! [ g cm / s in code unit] 
               w(i,j,k,M_MVX) = w(i,j,k,M_MVX) + dmom*r_i(MX) 
               w(i,j,k,M_MVY) = w(i,j,k,M_MVY) + dmom*r_i(MY) 
               w(i,j,k,M_MVZ) = w(i,j,k,M_MVZ) + dmom*r_i(MZ) 

               gpt = dmom/dtdvdrho    ! [ cm / s^-2 in code unit] 
               gxp = gpt*r_i(MX) 
               gyp = gpt*r_i(MY) 
               gzp = gpt*r_i(MZ) 

               ! --------------------
               ! energy injection
               ! --------------------
               w(i,j,k,M_MP)  = w(i,j,k,M_MP) + wmdotdt*c2_gam_w*eta(ic, jc, kc) &
                                + dtdvdrho * (  u(i,j,k,MVX)*gxp  &
                                               +u(i,j,k,MVY)*gyp  &
                                               +u(i,j,k,MVZ)*gzp )
  
               call w2u_4(w, u, dv)

               deallocate(w)

            enddo
         enddo
      enddo

      ! --------------
      ptr => ptr%next         
      ! --------------


    enddo

    deallocate(eta)
    deallocate(lsrank, lslevel, lsgid)

  
  end subroutine sp_outflow

  #undef M_MRHO
  #undef M_MVX
  #undef M_MVY
  #undef M_MVZ
  #undef M_MP

#endif

subroutine get_plev

  real(kind=DBL_KIND),dimension(MX:MZ) :: pos
  type(t_spParticle),pointer :: ptr
  integer :: np, ilev, ilgid, rank, i, j, k

  ! --------------
  ! step for np
  ! --------------
  np  = 0
  ptr => Particle

  do 
    if (.not. associated(ptr)) exit
    np = np + 1

    ! --------------------------------------------------
    !  シンク粒子の存在するグリッドの最深レベルを探査
    ! --------------------------------------------------
    pos = ptr%r
    do ilev = LevelMax, Lmin, -1 ! 一番上のlevelから探査
      call ob_getBlockFromCoordPhys(pos, ilev, ilgid, rank, i, j, k)
      if (ilgid /= Undefi) then
        ptr%lev = ilev
        exit                       ! gidが見つかったら終わり
      endif

      if (ilgid == Undefi .and. ilev == Lmin) then    ! gidがLminでも見つからなかったとりあえずlevsp = Lminに
        ptr%lev = Lmin
        exit
      endif
    enddo

    ! --------------------------------------------------
    !     シンクのいる最深レベルのセルのセルサイズ
    ! --------------------------------------------------
#ifdef USE_SPLEVEL_SP_CFL
    ptr%h(:) = hh_sp(:)
#else
    ptr%h(:) = CellWidth(:, Lmin)/(2**(ptr%lev-Lmin))
#endif

    ! --------------
    ptr => ptr%next
    ! --------------

  enddo

end subroutine get_plev

subroutine get_plev_kobetu(ptr)

  real(kind=DBL_KIND),dimension(MX:MZ) :: pos
  type(t_spParticle),pointer :: ptr
  integer :: ilev, ilgid, rank, i, j, k


  ! --------------------------------------------------
  !  シンク粒子の存在するグリッドの最深レベルを探査
  ! --------------------------------------------------
  pos = ptr%r
  do ilev = LevelMax, Lmin, -1 ! 一番上のlevelから探査
    call ob_getBlockFromCoordPhys(pos, ilev, ilgid, rank, i, j, k)
    if (ilgid /= Undefi) then
      ptr%lev = ilev
      exit                       ! gidが見つかったら終わり
    endif

    if (ilgid == Undefi .and. ilev == Lmin) then    ! gidがLminでも見つからなかったとりあえずlevsp = Lminに
      ptr%lev = Lmin
      exit
    endif
  enddo

  ! --------------------------------------------------
  !     シンクのいる最深レベルのセルのセルサイズ
  ! --------------------------------------------------
#ifdef USE_SPLEVEL_SP_CFL
  ptr%h(:) = hh_sp(:)
#else
  ptr%h(:) = CellWidth(:, Lmin)/(2**(ptr%lev-Lmin))
#endif

end subroutine get_plev_kobetu




! ---------------------------------------------------------------------
!                         get gravity of DM
! ---------------------------------------------------------------------
#ifdef DM_NFW_PROFILE
subroutine sp_gravityOfDM(gravDM, dt) ! get gravity of DM

    use kinzoku, only: gravDMh

    real(kind=DBL_KIND),dimension(MX:,:),intent(OUT) :: gravDM
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND) :: gx, gy, gz
    type(t_spParticle),pointer :: ptr
    integer :: np

    ptr => Particle
    np = 0
    do
      if (.not. associated(ptr)) exit
      np = np + 1
      
      call gravDMh(ptr%r(MX),ptr%r(MY),ptr%r(MZ),gx,gy,gz)
      
      gravDM(MX,np) = gx
      gravDM(MY,np) = gy
      gravDM(MZ,np) = gz

      ptr => ptr%next
    enddo

end subroutine sp_gravityOfDM
#endif
  

  ! -------------------------------------------------
  !          output sink creation regions 
  !        From level = Lmin-1 to LevelMax
  !     File name is spcrt. pid . step . level . d
  ! -------------------------------------------------
  subroutine output_spcrt(list_newPos)
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use unit, only : Unit_msun
    use mpilib
    use string, only : CHARLEN, num2char, concat

    type(t_spPos),pointer :: list_newPos  ! positions of newparticles
    type(t_spPos),pointer :: ptr

    integer,parameter :: NIug=64, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 16 MB/file)    
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=DBL_KIND) :: halfwidth, output_width
    integer :: level, dummy
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    character(len=CHARLEN) :: prefix

    integer :: ntot, npid

#ifdef OUTPUT_SP_CRTBOXSIZE
    output_width = 2.d0*OUTPUT_SP_CRTBOXSIZE*cgs_pc/Unit_l
#else
    output_width = 2.d0*cgs_pc/Unit_l
#endif

    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )

    ! get number of new sink particle
    ! -----------------
    ptr  => list_newPos
    ntot = 0
    ! -----------------
    do 
      if (.not. associated(ptr)) exit
      ntot = ntot + 1
      ptr => ptr%next
    enddo

    ! -------------------
    if (ntot < 1) return
    ! -------------------

    ! id of fisrt new sink particle
    npid = Nparticle-ntot
       
    ptr  => list_newPos

    do 
      if (.not. associated(ptr)) exit

      ! output
      do level = Lmin-2,LevelMax       !level = -2まで許してみる

        ! ------------------------------------------------------------------
        halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))       
        if(halfwidth > output_width) cycle
        ! ------------------------------------------------------------------

        if (get_myrank() == PRIMARY_RANK) then
             print '(/,A,I0,A,I0,A,(1P1E9.2),A)', "outputnewparticle: pid = ", npid &
               , " level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
        endif

          xmin = coordPhys(MX)
          ymin = coordPhys(MY)
          zmin = coordPhys(MZ)
          xmax = coordPhys(MZ+1+MX)
          ymax = coordPhys(MZ+1+MY)
          zmax = coordPhys(MZ+1+MZ)
          ! 領域の AND をとる
          xmin = max(xmin, ptr%r(MX)-halfwidth)
          ymin = max(ymin, ptr%r(MY)-halfwidth)
          zmin = max(zmin, ptr%r(MZ)-halfwidth)
          xmax = min(xmax, ptr%r(MX)+halfwidth)
          ymax = min(ymax, ptr%r(MY)+halfwidth)
          zmax = min(zmax, ptr%r(MZ)+halfwidth)

          ! speicfy prefix
          prefix = 'crsp.'
          prefix = concat(concat(prefix, num2char(npid)),'.')

          call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)
      enddo


      npid = npid + 1
      ptr  => ptr%next
    enddo

  end subroutine output_spcrt




end module sinkParticle

