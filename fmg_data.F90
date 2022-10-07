#include "config.h"
!-------------------------------------------------------------------------
! Module for data of FMG cycle
!-------------------------------------------------------------------------
module fmg_data
  use grid, only : Left, Right, Gidmin, GidListMax, Undefi
  use mpilib
  implicit none

  ! Specify PDE type
  integer,parameter :: FMG_PDE_TYPE_POISSON_EQUATION = 0
  integer,parameter :: FMG_PDE_TYPE_DIFFUSION_EQUATION = 1
  integer,parameter :: FMG_PDE_TYPE_OHMIC_DISSIPATION = 2
  integer,parameter :: FMG_PDE_TYPE_AMBIP_DIFFUSION= 3
  integer,parameter :: FMG_PDE_NMAX = 3
  integer,save :: FMG_PDE_TYPE = FMG_PDE_TYPE_POISSON_EQUATION ! set default
  logical,save :: FMG_PDE_LINEAR = .TRUE.
#ifdef FMG_LAMBDA
  real(kind=DBL_KIND),save :: LAMBDA = FMG_LAMBDA
#endif !FMG_LAMBDA
  integer,save :: AMR_LevelMin
  integer,save :: AMR_LevelMax
  integer,save :: FMG_LevelMin
  integer,save :: FMG_LevelMax
  integer,save :: Mmin          ! size of vector components.
  integer,save :: Mmax
  integer,save :: Nmin          ! size of space
  integer,save :: Nmax
  ! ----------------
  ! data structures
  ! ----------------
  ! node local
  ! GridLevel(amr_level,fmg_level)%Block(grid_ID)%{U,Res,Rhs,Rho}(i,j,k,m)
  ! GridLevel(amr_level,fmg_level)%Block(grid_ID)%F(i,j,k,m,n)
  type t_block
     real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: Dbg=> null(), U => null(), Res => null(), Rhs => null(), Rho => null(), Src => null(), Psi => null(), Ruf => null(), Eta => null(), Dod => null(), Dhe => null(), Dad => null()
     real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: F => null()
  end type t_block
  type t_gridlevel
     type(t_block),pointer,dimension(:) :: Block => null()
     real(kind=DBL_KIND) :: h, ds, dv ! width, cross-section, and volume of a cell
  end type t_gridlevel
  type(t_gridlevel),save,pointer,dimension(:,:) :: GridLevel => null()

  ! node global
  ! Ancestry(amr_level,rank)%Block(grid_ID)%ParentGid
  ! Ancestry(amr_level,rank)%Block(grid_ID)%ChildGid(lr,lr,lr)
  ! Ancestry(amr_level,rank)%Block(grid_ID)%NeighborGid(lr,ndirection)
  ! Ancestry(amr_level,rank)%Block(grid_ID)%NeighborSameLevel(lr,ndirection)
  type t_ablock
     integer :: ParentGid=Undefi, CoarseGid=Undefi, ParentRank=MPI_PROC_NULL, CoarseRank=MPI_PROC_NULL
     integer,dimension(Left:Right,Left:Right,Left:Right) :: ChildGid=Undefi, ChildRank=MPI_PROC_NULL
     integer,dimension(Left:Right,MX:MZ) :: NeighborGid=Undefi, NeighborRank=MPI_PROC_NULL
     logical,dimension(Left:Right,MX:MZ) :: NeighborSameLevel=.FALSE. ! TRUE:same-level
     logical,dimension(Left:Right,MX:MZ) :: NeighborParentLevel=.FALSE. ! TRUE :parent-level
  end type t_ablock
  type t_ancestry
     type(t_ablock),pointer,dimension(:) :: Block => null()
  end type t_ancestry
  type(t_ancestry),save,pointer,dimension(:,:) :: Ancestry => null()
  ! ÎÙ¤¬Æ±¤¸¥ì¥Ù¥ë
  !   Ancestry()%Block()%NeighborSameLevel() == .TRUE.
  ! ÎÙ¤¬¿Æ¤Î¥ì¥Ù¥ë
  !   Ancestry()%Block()%NeighborParentLevel() == .TRUE.
  ! ¿Æ¤¬ÊªÍý¶­³¦
  !   Geom()%Block()%TouchBoundary() == .TRUE.
  !   ¤Þ¤¿¤Ï
  !  Ancestry()%Block()%NeighborSameLevel() == .FALSE.
  !  Ancestry()%Block()%NeighborParentLevel() == .FALSE.

  ! node local
  ! Geom(amr_level)%Block(grid_ID)%TouchBoundary()
  type t_gblock
     logical,dimension(Left:Right,MX:MZ) :: TouchBoundary=.FALSE.
  end type t_gblock
  type t_geom
     type(t_gblock),pointer,dimension(:) :: Block => null()
  end type t_geom
  type(t_geom),save,pointer,dimension(:) :: Geom => null()

  ! node local (common for each local)
  ! GridSize(fmg_level)%{Imin,Imax,Jmin,Jmax,Kmin,Kmax}
  type t_gridsize
     integer :: Imin, Jmin, Kmin, Imax, Jmax, Kmax
  end type t_gridsize
  type(t_gridsize),save,pointer,dimension(:) :: GridSize => null()

  integer,save :: Ngh   ! number of ghost cell
  integer,parameter :: IDBG=-1, IU=0, IRES=1, IRHS=2, IRHO=3, IFLX=4, ISRC=5, IPSI=6, IRUF=7, IETA=8, IDOD=9, IDHE=10, IDAD=11  ! code for data
  integer,parameter :: ICODEMIN = -1  ! min of icode
  integer,parameter :: ICODEMAX = 11 ! max of icode
  integer,save :: FmgLevel_fill0 = Undefi ! level filled by zero
  logical,save,dimension(ICODEMIN:ICODEMAX) :: FMG_bool_isVector
  integer,parameter :: MUX = MX, MUY = MY, MUZ = MZ ! component of U vector
  real(kind=DBL_KIND),save :: Dtime                      ! Delta time for implicit scheme
  integer,save :: Dt_refineRatio = 1 ! refinement ratio of Dtime
  ! generic interface
  interface fmg_arrp
     module procedure fmg_arrp_vect, fmg_arrp_scalar
  end interface fmg_arrp
  interface fmg_fp
     module procedure fmg_fp_vect, fmg_fp_scalar
  end interface fmg_fp
contains
  !-------------------------------------------------------------------------
  ! set size of fmg data
  ! (vcomponentMin, vcomponentMax) = Vector size. Default is (MX, MX), corresponding to scalar.
  ! (dimensionMin, dimensionMax) = Dimension size. Default is (MX, MZ), corresponding to 3-dimension.
  !-------------------------------------------------------------------------
  subroutine fmg_data_init(vcomponentMin, vcomponentMax, dimensionMin, dimensionMax, NghostCell)
    use grid, only: Lmin, LevelMax, Dtime
    integer,intent(IN),optional :: vcomponentMin, vcomponentMax, dimensionMin, dimensionMax, NghostCell
    FMG_LevelMax = log(dble(min((NI), (NJ), (NK)))) / log(2.d0) - 1 + 0.5
    FMG_LevelMin = 0
    AMR_LevelMin = Lmin
    AMR_LevelMax = LevelMax

    if (present(NghostCell)) then
       Ngh = NghostCell
    else
       Ngh = 1                  ! default
    end if

    if (present(vcomponentMin)) then
       Mmin = vcomponentMin
    else
       Mmin = MX
    endif
    if (present(vcomponentMax)) then
       Mmax = vcomponentMax
    else
       Mmax = MX
    endif
    if (present(dimensionMin)) then
       Nmin = dimensionMin
    else
       Nmin = MX
    endif
    if (present(dimensionMax)) then
       Nmax = dimensionMax
    else
       Nmax = MZ
    endif
    if ( Mmin == Mmax ) then    ! set default
       FMG_bool_isVector(:) = .FALSE.
    else
       FMG_bool_isVector(:) = .TRUE.
    endif

    ! set Dtime
!!$    Dt_refineRatio = 1
    call fmg_set_dtime(Dtime(AMR_LevelMax)/Dt_refineRatio)

  end subroutine fmg_data_init
  !-------------------------------------------------------------------------
  function fmg_get_arrp(amr_level, fmg_level, grid_id, icode) result(arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Dbg
    case (IU)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%U
    case (IRES)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Res
    case (IRHS)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Rhs
    case (IRHO)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Rho
    case (ISRC)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Src
    case (IPSI)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Psi
    case (IRUF)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Ruf
    case (IETA)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Eta
    case (IDOD)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Dod
    case (IDHE)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Dhe
    case (IDAD)
       arrp => GridLevel(amr_level,fmg_level)%Block(grid_id)%Dad
    case default
       print *, '*** This code is not supported', icode
       stop
    end select
  end function fmg_get_arrp
  !-------------------------------------------------------------------------
  subroutine fmg_arrp_vect(amr_level, fmg_level, grid_id, icode, arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    arrp => fmg_get_arrp(amr_level, fmg_level, grid_id, icode)
    if (.not. fmg_isVector(icode)) then
       print *, '*** error in fmg_arrp_vect: icode is scalar', icode
       stop
    endif
  end subroutine fmg_arrp_vect
  !-------------------------------------------------------------------------
  subroutine fmg_arrp_scalar(amr_level, fmg_level, grid_id, icode, arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: arrp
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrpv
    arrpv => fmg_get_arrp(amr_level, fmg_level, grid_id, icode)
    arrp => fmg_slice3d(arrpv(:,:,:,Mmin))
    if (fmg_isVector(icode)) then
       print *, '*** error in fmg_arrp_scalar: icode is vector', icode
       stop
    endif
  end subroutine fmg_arrp_scalar
  !-------------------------------------------------------------------------
  function fmg_slice3d(arr) result(ptra)
    real(kind=DBL_KIND),intent(IN),dimension(-Ngh:,-Ngh:,-Ngh:),target :: arr
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: ptra
    ptra => arr
  end function fmg_slice3d
  !-------------------------------------------------------------------------
  function fmg_get_fp(amr_level, fmg_level, grid_id) result(fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fp
    fp => GridLevel(amr_level,fmg_level)%Block(grid_id)%F
  end function fmg_get_fp
  !-------------------------------------------------------------------------
  subroutine fmg_fp_vect(amr_level, fmg_level, grid_id, fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fp
    fp => fmg_get_fp(amr_level, fmg_level, grid_id)
  end subroutine fmg_fp_vect
  !-------------------------------------------------------------------------
  subroutine fmg_fp_scalar(amr_level, fmg_level, grid_id, fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: fp
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fpv
    fpv => fmg_get_fp(amr_level, fmg_level, grid_id)
    fp => fmg_slice3d_f(fpv(:,:,:,:,Mmin))
  end subroutine fmg_fp_scalar
  !-------------------------------------------------------------------------
  function fmg_slice3d_f(arr) result(ptra)
    real(kind=DBL_KIND),intent(IN),dimension(-1:,-1:,-1:,Mmin:),target :: arr
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptra
    ptra => arr
  end function fmg_slice3d_f
  !-------------------------------------------------------------------------
  subroutine fmg_set_arrp(arrp, amr_level, fmg_level, grid_id, icode)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Dbg => arrp
    case (IU)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%U   => arrp
    case (IRES)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Res => arrp
    case (IRHS)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Rhs => arrp
    case (IRHO)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Rho => arrp
    case (ISRC)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Src => arrp
    case (IPSI)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Psi => arrp
    case (IRUF)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Ruf => arrp
    case (IETA)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Eta => arrp
    case (IDOD)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Dod => arrp
    case (IDHE)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Dhe => arrp
    case (IDAD)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Dad => arrp
    case default
       print *, '*** This code is not supported', icode
       stop
    end select
  end subroutine fmg_set_arrp
  !-------------------------------------------------------------------------
  ! make LoL
  !-------------------------------------------------------------------------
  subroutine fmg_alloc_LoL
    use grid, only : CellWidth
    integer :: n, m, rank, listmax
    ! make GridLevel()%Block(), node local
    allocate( GridLevel( AMR_LevelMin : AMR_LevelMax, FMG_LevelMin : FMG_LevelMax ) )
    do n = lbound(GridLevel,2), ubound(GridLevel,2)    ! for each FMG level
       do m = lbound(GridLevel,1), ubound(GridLevel,1) ! for each AMR level
          allocate(GridLevel(m,n)%Block(Gidmin:GidListMax(m)))
          GridLevel(m,n)%h  = CellWidth(MX, m) * 2**n
          GridLevel(m,n)%ds = GridLevel(m,n)%h**2
          GridLevel(m,n)%dv = GridLevel(m,n)%h**3
       enddo
    enddo
    ! make Ancestry()%Block(), node global
    allocate( Ancestry( AMR_LevelMin : AMR_LevelMax, 0:NPE-1 ) )
    myrank = get_myrank()
    do rank = lbound(Ancestry,2), ubound(Ancestry,2) ! for each rank
       do m = lbound(Ancestry,1), ubound(Ancestry,1) ! for each AMR level
          if ( rank == myrank ) listmax = GidListMax( m )
          call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
          allocate( Ancestry(m,rank)%Block(Gidmin:listmax) )
       enddo
    enddo
    ! make Geom, node local
    allocate( Geom( AMR_LevelMin : AMR_LevelMax ) )
    do m = lbound(Geom,1), ubound(Geom,1)
       allocate(Geom(m)%Block(Gidmin:GidListMax(m)))
    enddo

    ! make GridSize, node common
    allocate( GridSize( FMG_LevelMin : FMG_LevelMax ) )
    do n = lbound(GridSize,1), ubound(GridSize,1)
       GridSize(n)%Imin = 0
       GridSize(n)%Jmin = 0
       GridSize(n)%Kmin = 0
       GridSize(n)%Imax = (NI)/2**n -1
       GridSize(n)%Jmax = (NJ)/2**n -1
       GridSize(n)%Kmax = (NK)/2**n -1
    enddo
  end subroutine fmg_alloc_LoL
  !-------------------------------------------------------------------------
  ! delete LoL
  !-------------------------------------------------------------------------
#define FREE_POINTER_SAFE(ARRAY) \
             if (associated(ARRAY)) then ;\
                deallocate(ARRAY) ;\
                nullify(ARRAY) ;\
             endif

#define FREE_POINTER(ARRAY) \
             deallocate(ARRAY) ;\
             nullify(ARRAY)

  subroutine fmg_dealloc_LoL
    integer :: m, n, gid
    ! deallocate GridLevel
    do n = lbound(GridLevel,2), ubound(GridLevel,2)
       do m = lbound(GridLevel,1), ubound(GridLevel,1)
          do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Dbg)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%U)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Res)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Rhs)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Rho)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Src)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Psi)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Ruf)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Eta)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Dod)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Dhe)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%Dad)
             FREE_POINTER_SAFE(GridLevel(m,n)%Block(gid)%F)
          enddo
          FREE_POINTER(GridLevel(m,n)%Block)
       enddo
    enddo
    deallocate( GridLevel )
    ! deallocate Ancestry
    do n = lbound(Ancestry,2), ubound(Ancestry,2)
       do m = lbound(Ancestry,1), ubound(Ancestry,1)
          FREE_POINTER( Ancestry(m,n)%Block )
       enddo
    enddo
    FREE_POINTER( Ancestry )
    ! deallocate Geom
    do m = lbound(Geom,1), ubound(Geom,1)
       FREE_POINTER( Geom(m)%Block )
    enddo
    FREE_POINTER( Geom )
    ! deallocate GridSize
    FREE_POINTER( GridSize )
  end subroutine fmg_dealloc_LoL
  !-------------------------------------------------------------------------
  ! allocate array
  !-------------------------------------------------------------------------
  subroutine fmg_alloc_arr( fmglevel, icode )
    integer,intent(IN) :: fmglevel, icode
    integer :: imax, jmax, kmax, imin, jmin, kmin, n, m, gid
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    n = fmglevel
    imin = fmg_get_imingh( fmglevel )
    imax = fmg_get_imaxgh( fmglevel )
    jmin = fmg_get_jmingh( fmglevel )
    jmax = fmg_get_jmaxgh( fmglevel )
    kmin = fmg_get_kmingh( fmglevel )
    kmax = fmg_get_kmaxgh( fmglevel )
    do m = lbound(GridLevel,1), ubound(GridLevel,1) ! for each AMR level
       do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1) ! for each block

          arrp => fmg_get_arrp(m, fmglevel, gid, icode)
          if ( associated( arrp ) ) cycle
          if ( fmg_isVector(icode) ) then
             allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmax))
             arrp = 0.d0 !HF added allocate ¿¿¿Nan¿¿¿¿¿¿¿? 
          else
             allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmin))
             arrp = 0.d0 !HF added 
          end if
          call fmg_set_arrp(arrp, m, fmglevel, gid, icode)
       enddo
    enddo
  end subroutine fmg_alloc_arr
  !-------------------------------------------------------------------------
  ! allocate F for all block
  !-------------------------------------------------------------------------
  subroutine fmg_alloc_f( fmglevel )
    integer,intent(IN) :: fmglevel
    integer :: imax, jmax, kmax, imin, jmin, kmin, n, m, gid
    n = fmglevel
!!$    imin = GridSize(fmglevel)%Imin-1 ! no margin
!!$    imax = GridSize(fmglevel)%Imax
!!$    jmin = GridSize(fmglevel)%Jmin-1
!!$    jmax = GridSize(fmglevel)%Jmax
!!$    kmin = GridSize(fmglevel)%Kmin-1
!!$    kmax = GridSize(fmglevel)%Kmax
    imin = fmg_get_imingh( fmglevel )
    imax = fmg_get_imaxgh( fmglevel )
    jmin = fmg_get_jmingh( fmglevel )
    jmax = fmg_get_jmaxgh( fmglevel )
    kmin = fmg_get_kmingh( fmglevel )
    kmax = fmg_get_kmaxgh( fmglevel )
    do m = lbound(GridLevel,1), ubound(GridLevel,1) ! for each AMR level
       do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1) ! for each block
          if ( associated(GridLevel(m,n)%Block(gid)%F) ) cycle
          allocate(GridLevel(m,n)%Block(gid)%F(imin:imax,jmin:jmax,kmin:kmax,Nmin:Nmax,Mmin:Mmax))
       enddo
    enddo
  end subroutine fmg_alloc_f
  !-------------------------------------------------------------------------
  ! get minimum grid id for this rank (independent of FMGlevel)
  !-------------------------------------------------------------------------
  function fmg_get_gidmin(amrlevel) result(idmin)
    integer,intent(IN) :: amrlevel
    integer :: idmin
    idmin = lbound(GridLevel(amrlevel, FMG_LevelMin)%Block, 1)
  end function fmg_get_gidmin
  !-------------------------------------------------------------------------
  ! get maximum grid id for this rank (independent of FMGlevel)
  !-------------------------------------------------------------------------
  function fmg_get_gidmax(amrlevel) result(idmax)
    integer,intent(IN) :: amrlevel
    integer :: idmax
    idmax = ubound(GridLevel(amrlevel, FMG_LevelMin)%Block, 1)
  end function fmg_get_gidmax
  !-------------------------------------------------------------------------
  ! get minimum grid id for given amrlevel and rank
  !-------------------------------------------------------------------------
  function fmg_get_gidmin_rank(amrlevel, rank) result(idmin)
    integer,intent(IN) :: amrlevel, rank
    integer :: idmin
    idmin = lbound(Ancestry(amrlevel,rank)%Block, 1)
  end function fmg_get_gidmin_rank
  !-------------------------------------------------------------------------
  ! get maximum grid id for given amrlevel and rank
  !-------------------------------------------------------------------------
  function fmg_get_gidmax_rank(amrlevel, rank) result(idmax)
    integer,intent(IN) :: amrlevel, rank
    integer :: idmax
    idmax = ubound(Ancestry(amrlevel,rank)%Block, 1)
  end function fmg_get_gidmax_rank
  !-------------------------------------------------------------------------
  ! get array size without ghost cell
  !-------------------------------------------------------------------------
  function fmg_get_imin(fmglevel) result(imin)
    integer,intent(IN) :: fmglevel
    integer :: imin
    imin = GridSize(fmglevel)%Imin
  end function fmg_get_imin
  !-------------------------------------------------------------------------
  function fmg_get_jmin(fmglevel) result(jmin)
    integer,intent(IN) :: fmglevel
    integer :: jmin
    jmin = GridSize(fmglevel)%Jmin
  end function fmg_get_jmin
  !-------------------------------------------------------------------------
  function fmg_get_kmin(fmglevel) result(kmin)
    integer,intent(IN) :: fmglevel
    integer :: kmin
    kmin = GridSize(fmglevel)%Kmin
  end function fmg_get_kmin
  !-------------------------------------------------------------------------
  function fmg_get_imax(fmglevel) result(imax)
    integer,intent(IN) :: fmglevel
    integer :: imax
    imax = GridSize(fmglevel)%Imax
  end function fmg_get_imax
  !-------------------------------------------------------------------------
  function fmg_get_jmax(fmglevel) result(jmax)
    integer,intent(IN) :: fmglevel
    integer :: jmax
    jmax = GridSize(fmglevel)%Jmax
  end function fmg_get_jmax
  !-------------------------------------------------------------------------
  function fmg_get_kmax(fmglevel) result(kmax)
    integer,intent(IN) :: fmglevel
    integer :: kmax
    kmax = GridSize(fmglevel)%Kmax
  end function fmg_get_kmax
  !-------------------------------------------------------------------------
  ! get array size with ghost cell
  !-------------------------------------------------------------------------
  function fmg_get_imingh(fmglevel) result(imingh)
    integer,intent(IN) :: fmglevel
    integer :: imingh
    imingh = GridSize(fmglevel)%Imin-Ngh
  end function fmg_get_imingh
  !-------------------------------------------------------------------------
  function fmg_get_jmingh(fmglevel) result(jmingh)
    integer,intent(IN) :: fmglevel
    integer :: jmingh
    jmingh = GridSize(fmglevel)%Jmin-Ngh
  end function fmg_get_jmingh
  !-------------------------------------------------------------------------
  function fmg_get_kmingh(fmglevel) result(kmingh)
    integer,intent(IN) :: fmglevel
    integer :: kmingh
    kmingh = GridSize(fmglevel)%Kmin-Ngh
  end function fmg_get_kmingh
  !-------------------------------------------------------------------------
  function fmg_get_imaxgh(fmglevel) result(imaxgh)
    integer,intent(IN) :: fmglevel
    integer :: imaxgh
    imaxgh = GridSize(fmglevel)%Imax+Ngh
  end function fmg_get_imaxgh
  !-------------------------------------------------------------------------
  function fmg_get_jmaxgh(fmglevel) result(jmaxgh)
    integer,intent(IN) :: fmglevel
    integer :: jmaxgh
    jmaxgh = GridSize(fmglevel)%Jmax+Ngh
  end function fmg_get_jmaxgh
  !-------------------------------------------------------------------------
  function fmg_get_kmaxgh(fmglevel) result(kmaxgh)
    integer,intent(IN) :: fmglevel
    integer :: kmaxgh
    kmaxgh = GridSize(fmglevel)%Kmax+Ngh
  end function fmg_get_kmaxgh
  !-------------------------------------------------------------------------
  subroutine fmg_get_gridsize(fmglevel, imin,jmin,kmin,imax,jmax,kmax)
    integer,intent(IN) :: fmglevel
    integer,intent(OUT) :: imin,jmin,kmin,imax,jmax,kmax
    imin = GridSize(fmglevel)%Imin
    jmin = GridSize(fmglevel)%Jmin
    kmin = GridSize(fmglevel)%Kmin
    imax = GridSize(fmglevel)%Imax
    jmax = GridSize(fmglevel)%Jmax
    kmax = GridSize(fmglevel)%Kmax
  end subroutine fmg_get_gridsize
  !-------------------------------------------------------------------------
  !  get grid size with ghostcell
  !-------------------------------------------------------------------------
  subroutine fmg_get_gridsizeGh(fmglevel, imin,jmin,kmin,imax,jmax,kmax)
    integer,intent(IN) :: fmglevel
    integer,intent(OUT) :: imin,jmin,kmin,imax,jmax,kmax
    imin = fmg_get_imingh( fmglevel )
    imax = fmg_get_imaxgh( fmglevel )
    jmin = fmg_get_jmingh( fmglevel )
    jmax = fmg_get_jmaxgh( fmglevel )
    kmin = fmg_get_kmingh( fmglevel )
    kmax = fmg_get_kmaxgh( fmglevel )
  end subroutine fmg_get_gridsizeGh
  !-------------------------------------------------------------------------
  ! ¼«Ê¬¤Ï¿Æ¤ËÂÐ¤·¤Æ±¦¤«º¸¤«
  !-------------------------------------------------------------------------
  subroutine fmg_left_or_right(gid, rank, amr_level, lri, lrj, lrk)
    integer,intent(IN) :: gid, rank, amr_level
    integer,intent(OUT) :: lri, lrj, lrk
    integer :: i,j,k, pgid, prank
    lri = Undefi
    lrj = Undefi
    lrk = Undefi

    pgid  = Ancestry(amr_level,rank)%Block(gid)%ParentGid
    prank = Ancestry(amr_level,rank)%Block(gid)%ParentRank

    if ( pgid == Undefi .or. prank == MPI_PROC_NULL ) then
       print *, '*** invarid link list in left_or_right'
       stop
    endif

    do k = Left,Right
       do j = Left, Right
          do i = Left, Right
             if ( Ancestry(amr_level-1,prank)%Block(pgid)%ChildGid(i,j,k) == gid .and. &
                  Ancestry(amr_level-1,prank)%Block(pgid)%ChildRank(i,j,k) == rank ) then
                lri = i
                lrj = j
                lrk = k
             endif
          enddo
       enddo
    enddo
    if (lri < 0) then
       print *, '*** inconsistent link list in fmg_data'
       stop
    endif

  end subroutine fmg_left_or_right
  !-------------------------------------------------------------------------
  function fmg_get_h(amrlev,fmglev) result(h)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=DBL_KIND) :: h
    h = GridLevel(amrlev,fmglev)%h
  end function fmg_get_h
  !-------------------------------------------------------------------------
  function fmg_get_ds(amrlev,fmglev) result(ds)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=DBL_KIND) :: ds
    ds = GridLevel(amrlev,fmglev)%ds
  end function fmg_get_ds
  !-------------------------------------------------------------------------
  function fmg_get_dv(amrlev,fmglev) result(dv)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=DBL_KIND) :: dv
    dv = GridLevel(amrlev,fmglev)%dv
  end function fmg_get_dv
  !-------------------------------------------------------------------------
  ! skip grid having child grid
  !-------------------------------------------------------------------------
  function fmg_skip_grid(gid, amrlev, fmglev) result(bool)
    use mpilib
    integer,intent(IN) :: gid, amrlev, fmglev
    logical :: bool
!!$    bool = .false.
!!$    return
    myrank = get_myrank()
    bool = amrlev /= AMR_LevelMax &
         .and. &
         Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi
  end function fmg_skip_grid
  !-------------------------------------------------------------------------
  ! have child grid
  !-------------------------------------------------------------------------
  function fmg_have_child(gid, amrlev) result(bool)
    integer,intent(IN) :: gid, amrlev
    logical :: bool
    bool = Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi
  end function fmg_have_child
  !-------------------------------------------------------------------------
  ! get maximaum
  !-------------------------------------------------------------------------
  subroutine fmg_max(maxarr, amrlev, fmglev, icode, absolute, skip, xskip, mpireduce)
    real(kind=DBL_KIND),intent(OUT) :: maxarr
    integer,intent(IN) :: amrlev, fmglev, icode
    logical,intent(IN),optional :: absolute, skip, xskip, mpireduce
    logical :: bool_abs, bool_skip, bool_xskip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    real(kind=DBL_KIND) :: maxlocal
    bool_abs = .FALSE.;  if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.; if (present(skip)) bool_skip = skip
    bool_xskip = .FALSE.; if (present(xskip)) bool_xskip = xskip
    bool_mpi = .TRUE.;  if (present(mpireduce)) bool_mpi = mpireduce
    if (bool_abs) then
       maxlocal = 0.d0
    else
       maxlocal = -HUGE(maxlocal)
    end if
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       if (bool_xskip .and. .not. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       arrp => fmg_get_arrp(amrlev, fmglev, gid, icode)
       if (bool_abs) then
          maxlocal = max(maxlocal, MAXVAL(abs(arrp(SZ))))
       else
          maxlocal = max(maxlocal, MAXVAL(arrp(SZ)))
       end if
    end do
    if (bool_mpi) then
       call mpi_allreduce(maxlocal, maxarr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    else
       maxarr = maxlocal
    end if
#undef SZ
  end subroutine fmg_max
  !-------------------------------------------------------------------------
  function fmg_get_max(amrlev, fmglev, icode, absolute, skip, mpireduce) result(maxarr)
    integer,intent(IN) :: amrlev, fmglev, icode
    logical,intent(IN),optional :: absolute, skip, mpireduce
    real(kind=DBL_KIND) :: maxarr
    call fmg_max(maxarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_max
  !-------------------------------------------------------------------------
  ! get sum
  !-------------------------------------------------------------------------
  subroutine fmg_sum(sumarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=DBL_KIND),intent(OUT) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    logical :: bool_abs, bool_skip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke
#define SZ is:ie,js:je,ks:ke,:
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    real(kind=DBL_KIND) :: dv, sumlocal
    bool_abs = .FALSE.; if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.; if (present(skip)) bool_skip = skip
    bool_mpi = .TRUE.;  if (present(mpireduce)) bool_mpi = mpireduce
    dv = fmg_get_dv(amrlev,fmglev)
    sumlocal = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       arrp => fmg_get_arrp(amrlev, fmglev, gid, icode)
       if (bool_abs) then
          sumlocal = sumlocal + SUM(abs(arrp(SZ)))
       else
          sumlocal = sumlocal + SUM(arrp(SZ))
       end if
    end do
    if (bool_mpi) then
       call mpi_allreduce(sumlocal, sumarr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    else
       sumarr = sumlocal
    end if
    sumarr = sumarr * dv
#undef SZ
  end subroutine fmg_sum
  !-------------------------------------------------------------------------
  function fmg_get_sum(amrlev, fmglev, icode, absolute, skip, mpireduce) result(sumarr)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=DBL_KIND) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_sum(sumarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_sum
  !-------------------------------------------------------------------------
  ! get average
  !-------------------------------------------------------------------------
  subroutine fmg_ave(avearr, amrlev, fmglev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=DBL_KIND),intent(OUT) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    logical :: bool_abs, bool_skip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke, ngrid
#define SZ is:ie,js:je,ks:ke,:
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    real(kind=DBL_KIND) :: avelocal
    bool_abs = .FALSE.; if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.; if (present(skip)) bool_skip = skip
    bool_mpi = .TRUE.;  if (present(mpireduce)) bool_mpi = mpireduce
    avelocal = 0.d0
    ngrid = 0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       ngrid = ngrid + 1
       arrp => fmg_get_arrp(amrlev, fmglev, gid, icode)
       if (bool_abs) then
          avelocal = avelocal + SUM(abs(arrp(SZ)))
       else
          avelocal = avelocal + SUM(arrp(SZ))
       end if
    end do
    ! local average (³ºÅö¤¹¤ë¥°¥ê¥Ã¥É¤¬¤Ê¤±¤ì¤Ð¡¢¥¼¥í¤òÊÖ¤¹¡£)
    if ( ngrid == 0 ) then
       avelocal = 0.d0
    else
       avelocal = avelocal / ( (ie-is+1)*(je-js+1)*(je-js+1)*size(arrp,4) * ngrid )
    endif
    if (bool_mpi) then
       call mpi_allreduce(avelocal, avearr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       avearr = avearr / mpi_get_npe()
    else
       avearr = avelocal
    end if
#undef SZ
  end subroutine fmg_ave
  !-------------------------------------------------------------------------
  function fmg_get_ave(amrlev, fmglev, icode, absolute, skip, mpireduce) result(avearr)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=DBL_KIND) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_ave(avearr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_ave
  !-------------------------------------------------------------------------
  ! return icode arrary is scalar
  !-------------------------------------------------------------------------
  function fmg_isVector(icode) result(bool)
    integer,intent(IN) :: icode
    logical :: bool
    bool = FMG_bool_isVector(icode)
  end function fmg_isVector
  !-------------------------------------------------------------------------
  ! set if array is scalar
  !-------------------------------------------------------------------------
  subroutine fmg_setVector(icode, bool)
    integer,intent(IN) :: icode
    logical,intent(IN) :: bool
    FMG_bool_isVector(icode) = bool
  end subroutine fmg_setVector
  !-------------------------------------------------------------------------
  ! get number of componets
  !-------------------------------------------------------------------------
  function fmg_getNcomp( icode ) result( ncomp )
    integer,intent(IN) :: icode
    integer :: ncomp
    if ( fmg_isVector(icode) ) then
       ncomp = Mmax - Mmin + 1
    else
       ncomp = 1
    end if
  end function fmg_getNcomp
  !-------------------------------------------------------------------------
  ! get dtime
  !-------------------------------------------------------------------------
  function fmg_get_dtime() result( dt )
    real(kind=DBL_KIND) :: dt
    dt = Dtime
  end function fmg_get_dtime
  !-------------------------------------------------------------------------
  ! set dtime
  !-------------------------------------------------------------------------
  subroutine fmg_set_dtime(dt)
    real(kind=DBL_KIND),intent(IN) :: dt
    Dtime = dt
  end subroutine fmg_set_dtime

end module fmg_data
