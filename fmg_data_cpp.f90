module fmg_data
  use grid, only : Left, Right, Gidmin, GidListMax, Undefi
  use mpilib
  implicit none
  integer,parameter :: FMG_PDE_TYPE_POISSON_EQUATION = 0
  integer,parameter :: FMG_PDE_TYPE_DIFFUSION_EQUATION = 1
  integer,parameter :: FMG_PDE_TYPE_OHMIC_DISSIPATION = 2
  integer,parameter :: FMG_PDE_TYPE_AMBIP_DIFFUSION= 3
  integer,parameter :: FMG_PDE_NMAX = 3
  integer,save :: FMG_PDE_TYPE = FMG_PDE_TYPE_POISSON_EQUATION 
  logical,save :: FMG_PDE_LINEAR = .TRUE.
  integer,save :: AMR_LevelMin
  integer,save :: AMR_LevelMax
  integer,save :: FMG_LevelMin
  integer,save :: FMG_LevelMax
  integer,save :: Mmin 
  integer,save :: Mmax
  integer,save :: Nmin 
  integer,save :: Nmax
  type t_block
     real(kind=8),pointer,dimension(:,:,:,:) :: Dbg=> null(), U => null(), Res => null(), Rhs => null(), Rho => null(), Src => null&
&(), Psi => null(), Ruf => null(), Eta => null(), Dod => null(), Dhe => null(), Dad => null()
     real(kind=8),pointer,dimension(:,:,:,:,:) :: F => null()
  end type t_block
  type t_gridlevel
     type(t_block),pointer,dimension(:) :: Block => null()
     real(kind=8) :: h, ds, dv 
  end type t_gridlevel
  type(t_gridlevel),save,pointer,dimension(:,:) :: GridLevel => null()
  type t_ablock
     integer :: ParentGid=Undefi, CoarseGid=Undefi, ParentRank=MPI_PROC_NULL, CoarseRank=MPI_PROC_NULL
     integer,dimension(Left:Right,Left:Right,Left:Right) :: ChildGid=Undefi, ChildRank=MPI_PROC_NULL
     integer,dimension(Left:Right,0:2) :: NeighborGid=Undefi, NeighborRank=MPI_PROC_NULL
     logical,dimension(Left:Right,0:2) :: NeighborSameLevel=.FALSE. 
     logical,dimension(Left:Right,0:2) :: NeighborParentLevel=.FALSE. 
  end type t_ablock
  type t_ancestry
     type(t_ablock),pointer,dimension(:) :: Block => null()
  end type t_ancestry
  type(t_ancestry),save,pointer,dimension(:,:) :: Ancestry => null()
  type t_gblock
     logical,dimension(Left:Right,0:2) :: TouchBoundary=.FALSE.
  end type t_gblock
  type t_geom
     type(t_gblock),pointer,dimension(:) :: Block => null()
  end type t_geom
  type(t_geom),save,pointer,dimension(:) :: Geom => null()
  type t_gridsize
     integer :: Imin, Jmin, Kmin, Imax, Jmax, Kmax
  end type t_gridsize
  type(t_gridsize),save,pointer,dimension(:) :: GridSize => null()
  integer,save :: Ngh 
  integer,parameter :: IDBG=-1, IU=0, IRES=1, IRHS=2, IRHO=3, IFLX=4, ISRC=5, IPSI=6, IRUF=7, IETA=8, IDOD=9, IDHE=10, IDAD=11 
  integer,parameter :: ICODEMIN = -1 
  integer,parameter :: ICODEMAX = 11 
  integer,save :: FmgLevel_fill0 = Undefi 
  logical,save,dimension(ICODEMIN:ICODEMAX) :: FMG_bool_isVector
  integer,parameter :: MUX = 0, MUY = 1, MUZ = 2 
  real(kind=8),save :: Dtime 
  integer,save :: Dt_refineRatio = 1 
  interface fmg_arrp
     module procedure fmg_arrp_vect, fmg_arrp_scalar
  end interface fmg_arrp
  interface fmg_fp
     module procedure fmg_fp_vect, fmg_fp_scalar
  end interface fmg_fp
contains
  subroutine fmg_data_init(vcomponentMin, vcomponentMax, dimensionMin, dimensionMax, NghostCell)
    use grid, only: Lmin, LevelMax, Dtime
    integer,intent(IN),optional :: vcomponentMin, vcomponentMax, dimensionMin, dimensionMax, NghostCell
    FMG_LevelMax = log(dble(min((8), (8), (8)))) / log(2.d0) - 1 + 0.5
    FMG_LevelMin = 0
    AMR_LevelMin = Lmin
    AMR_LevelMax = LevelMax
    if (present(NghostCell)) then
       Ngh = NghostCell
    else
       Ngh = 1 
    end if
    if (present(vcomponentMin)) then
       Mmin = vcomponentMin
    else
       Mmin = 0
    endif
    if (present(vcomponentMax)) then
       Mmax = vcomponentMax
    else
       Mmax = 0
    endif
    if (present(dimensionMin)) then
       Nmin = dimensionMin
    else
       Nmin = 0
    endif
    if (present(dimensionMax)) then
       Nmax = dimensionMax
    else
       Nmax = 2
    endif
    if ( Mmin == Mmax ) then 
       FMG_bool_isVector(:) = .FALSE.
    else
       FMG_bool_isVector(:) = .TRUE.
    endif
    call fmg_set_dtime(Dtime(AMR_LevelMax)/Dt_refineRatio)
  end subroutine fmg_data_init
  function fmg_get_arrp(amr_level, fmg_level, grid_id, icode) result(arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
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
  subroutine fmg_arrp_vect(amr_level, fmg_level, grid_id, icode, arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    arrp => fmg_get_arrp(amr_level, fmg_level, grid_id, icode)
    if (.not. fmg_isVector(icode)) then
       print *, '*** error in fmg_arrp_vect: icode is scalar', icode
       stop
    endif
  end subroutine fmg_arrp_vect
  subroutine fmg_arrp_scalar(amr_level, fmg_level, grid_id, icode, arrp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=8),pointer,dimension(:,:,:) :: arrp
    real(kind=8),pointer,dimension(:,:,:,:) :: arrpv
    arrpv => fmg_get_arrp(amr_level, fmg_level, grid_id, icode)
    arrp => fmg_slice3d(arrpv(:,:,:,Mmin))
    if (fmg_isVector(icode)) then
       print *, '*** error in fmg_arrp_scalar: icode is vector', icode
       stop
    endif
  end subroutine fmg_arrp_scalar
  function fmg_slice3d(arr) result(ptra)
    real(kind=8),intent(IN),dimension(-Ngh:,-Ngh:,-Ngh:),target :: arr
    real(kind=8),dimension(:,:,:),pointer :: ptra
    ptra => arr
  end function fmg_slice3d
  function fmg_get_fp(amr_level, fmg_level, grid_id) result(fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fp
    fp => GridLevel(amr_level,fmg_level)%Block(grid_id)%F
  end function fmg_get_fp
  subroutine fmg_fp_vect(amr_level, fmg_level, grid_id, fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fp
    fp => fmg_get_fp(amr_level, fmg_level, grid_id)
  end subroutine fmg_fp_vect
  subroutine fmg_fp_scalar(amr_level, fmg_level, grid_id, fp)
    integer,intent(IN) :: amr_level, fmg_level, grid_id
    real(kind=8),pointer,dimension(:,:,:,:) :: fp
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fpv
    fpv => fmg_get_fp(amr_level, fmg_level, grid_id)
    fp => fmg_slice3d_f(fpv(:,:,:,:,Mmin))
  end subroutine fmg_fp_scalar
  function fmg_slice3d_f(arr) result(ptra)
    real(kind=8),intent(IN),dimension(-1:,-1:,-1:,Mmin:),target :: arr
    real(kind=8),dimension(:,:,:,:),pointer :: ptra
    ptra => arr
  end function fmg_slice3d_f
  subroutine fmg_set_arrp(arrp, amr_level, fmg_level, grid_id, icode)
    integer,intent(IN) :: amr_level, fmg_level, grid_id, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%Dbg => arrp
    case (IU)
       GridLevel(amr_level,fmg_level)%Block(grid_id)%U => arrp
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
  subroutine fmg_alloc_LoL
    use grid, only : CellWidth
    integer :: n, m, rank, listmax
    allocate( GridLevel( AMR_LevelMin : AMR_LevelMax, FMG_LevelMin : FMG_LevelMax ) )
    do n = lbound(GridLevel,2), ubound(GridLevel,2) 
       do m = lbound(GridLevel,1), ubound(GridLevel,1) 
          allocate(GridLevel(m,n)%Block(Gidmin:GidListMax(m)))
          GridLevel(m,n)%h = CellWidth(0, m) * 2**n
          GridLevel(m,n)%ds = GridLevel(m,n)%h**2
          GridLevel(m,n)%dv = GridLevel(m,n)%h**3
       enddo
    enddo
    allocate( Ancestry( AMR_LevelMin : AMR_LevelMax, 0:400 -1 ) )
    myrank = get_myrank()
    do rank = lbound(Ancestry,2), ubound(Ancestry,2) 
       do m = lbound(Ancestry,1), ubound(Ancestry,1) 
          if ( rank == myrank ) listmax = GidListMax( m )
          call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
          allocate( Ancestry(m,rank)%Block(Gidmin:listmax) )
       enddo
    enddo
    allocate( Geom( AMR_LevelMin : AMR_LevelMax ) )
    do m = lbound(Geom,1), ubound(Geom,1)
       allocate(Geom(m)%Block(Gidmin:GidListMax(m)))
    enddo
    allocate( GridSize( FMG_LevelMin : FMG_LevelMax ) )
    do n = lbound(GridSize,1), ubound(GridSize,1)
       GridSize(n)%Imin = 0
       GridSize(n)%Jmin = 0
       GridSize(n)%Kmin = 0
       GridSize(n)%Imax = (8)/2**n -1
       GridSize(n)%Jmax = (8)/2**n -1
       GridSize(n)%Kmax = (8)/2**n -1
    enddo
  end subroutine fmg_alloc_LoL
  subroutine fmg_dealloc_LoL
    integer :: m, n, gid
    do n = lbound(GridLevel,2), ubound(GridLevel,2)
       do m = lbound(GridLevel,1), ubound(GridLevel,1)
          do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1)
             if (associated(GridLevel(m,n)%Block(gid)%Dbg)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Dbg) 
 nullify(GridLevel(m,n)%Block(gid)%Dbg) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%U)) then 
 deallocate(GridLevel(m,n)%Block(gid)%U) 
 nullify(GridLevel(m,n)%Block(gid)%U) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Res)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Res) 
 nullify(GridLevel(m,n)%Block(gid)%Res) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Rhs)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Rhs) 
 nullify(GridLevel(m,n)%Block(gid)%Rhs) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Rho)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Rho) 
 nullify(GridLevel(m,n)%Block(gid)%Rho) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Src)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Src) 
 nullify(GridLevel(m,n)%Block(gid)%Src) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Psi)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Psi) 
 nullify(GridLevel(m,n)%Block(gid)%Psi) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Ruf)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Ruf) 
 nullify(GridLevel(m,n)%Block(gid)%Ruf) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Eta)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Eta) 
 nullify(GridLevel(m,n)%Block(gid)%Eta) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Dod)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Dod) 
 nullify(GridLevel(m,n)%Block(gid)%Dod) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Dhe)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Dhe) 
 nullify(GridLevel(m,n)%Block(gid)%Dhe) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%Dad)) then 
 deallocate(GridLevel(m,n)%Block(gid)%Dad) 
 nullify(GridLevel(m,n)%Block(gid)%Dad) 
 endif
             if (associated(GridLevel(m,n)%Block(gid)%F)) then 
 deallocate(GridLevel(m,n)%Block(gid)%F) 
 nullify(GridLevel(m,n)%Block(gid)%F) 
 endif
          enddo
          deallocate(GridLevel(m,n)%Block) 
 nullify(GridLevel(m,n)%Block)
       enddo
    enddo
    deallocate( GridLevel )
    do n = lbound(Ancestry,2), ubound(Ancestry,2)
       do m = lbound(Ancestry,1), ubound(Ancestry,1)
          deallocate(Ancestry(m,n)%Block) 
 nullify(Ancestry(m,n)%Block)
       enddo
    enddo
    deallocate(Ancestry) 
 nullify(Ancestry)
    do m = lbound(Geom,1), ubound(Geom,1)
       deallocate(Geom(m)%Block) 
 nullify(Geom(m)%Block)
    enddo
    deallocate(Geom) 
 nullify(Geom)
    deallocate(GridSize) 
 nullify(GridSize)
  end subroutine fmg_dealloc_LoL
  subroutine fmg_alloc_arr( fmglevel, icode )
    integer,intent(IN) :: fmglevel, icode
    integer :: imax, jmax, kmax, imin, jmin, kmin, n, m, gid
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    n = fmglevel
    imin = fmg_get_imingh( fmglevel )
    imax = fmg_get_imaxgh( fmglevel )
    jmin = fmg_get_jmingh( fmglevel )
    jmax = fmg_get_jmaxgh( fmglevel )
    kmin = fmg_get_kmingh( fmglevel )
    kmax = fmg_get_kmaxgh( fmglevel )
    do m = lbound(GridLevel,1), ubound(GridLevel,1) 
       do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1) 
          arrp => fmg_get_arrp(m, fmglevel, gid, icode)
          if ( associated( arrp ) ) cycle
          if ( fmg_isVector(icode) ) then
             allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmax))
             arrp = 0.d0 
          else
             allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmin))
             arrp = 0.d0 
          end if
          call fmg_set_arrp(arrp, m, fmglevel, gid, icode)
       enddo
    enddo
  end subroutine fmg_alloc_arr
  subroutine fmg_alloc_f( fmglevel )
    integer,intent(IN) :: fmglevel
    integer :: imax, jmax, kmax, imin, jmin, kmin, n, m, gid
    n = fmglevel
    imin = fmg_get_imingh( fmglevel )
    imax = fmg_get_imaxgh( fmglevel )
    jmin = fmg_get_jmingh( fmglevel )
    jmax = fmg_get_jmaxgh( fmglevel )
    kmin = fmg_get_kmingh( fmglevel )
    kmax = fmg_get_kmaxgh( fmglevel )
    do m = lbound(GridLevel,1), ubound(GridLevel,1) 
       do gid = lbound(GridLevel(m,n)%Block,1), ubound(GridLevel(m,n)%Block,1) 
          if ( associated(GridLevel(m,n)%Block(gid)%F) ) cycle
          allocate(GridLevel(m,n)%Block(gid)%F(imin:imax,jmin:jmax,kmin:kmax,Nmin:Nmax,Mmin:Mmax))
       enddo
    enddo
  end subroutine fmg_alloc_f
  function fmg_get_gidmin(amrlevel) result(idmin)
    integer,intent(IN) :: amrlevel
    integer :: idmin
    idmin = lbound(GridLevel(amrlevel, FMG_LevelMin)%Block, 1)
  end function fmg_get_gidmin
  function fmg_get_gidmax(amrlevel) result(idmax)
    integer,intent(IN) :: amrlevel
    integer :: idmax
    idmax = ubound(GridLevel(amrlevel, FMG_LevelMin)%Block, 1)
  end function fmg_get_gidmax
  function fmg_get_gidmin_rank(amrlevel, rank) result(idmin)
    integer,intent(IN) :: amrlevel, rank
    integer :: idmin
    idmin = lbound(Ancestry(amrlevel,rank)%Block, 1)
  end function fmg_get_gidmin_rank
  function fmg_get_gidmax_rank(amrlevel, rank) result(idmax)
    integer,intent(IN) :: amrlevel, rank
    integer :: idmax
    idmax = ubound(Ancestry(amrlevel,rank)%Block, 1)
  end function fmg_get_gidmax_rank
  function fmg_get_imin(fmglevel) result(imin)
    integer,intent(IN) :: fmglevel
    integer :: imin
    imin = GridSize(fmglevel)%Imin
  end function fmg_get_imin
  function fmg_get_jmin(fmglevel) result(jmin)
    integer,intent(IN) :: fmglevel
    integer :: jmin
    jmin = GridSize(fmglevel)%Jmin
  end function fmg_get_jmin
  function fmg_get_kmin(fmglevel) result(kmin)
    integer,intent(IN) :: fmglevel
    integer :: kmin
    kmin = GridSize(fmglevel)%Kmin
  end function fmg_get_kmin
  function fmg_get_imax(fmglevel) result(imax)
    integer,intent(IN) :: fmglevel
    integer :: imax
    imax = GridSize(fmglevel)%Imax
  end function fmg_get_imax
  function fmg_get_jmax(fmglevel) result(jmax)
    integer,intent(IN) :: fmglevel
    integer :: jmax
    jmax = GridSize(fmglevel)%Jmax
  end function fmg_get_jmax
  function fmg_get_kmax(fmglevel) result(kmax)
    integer,intent(IN) :: fmglevel
    integer :: kmax
    kmax = GridSize(fmglevel)%Kmax
  end function fmg_get_kmax
  function fmg_get_imingh(fmglevel) result(imingh)
    integer,intent(IN) :: fmglevel
    integer :: imingh
    imingh = GridSize(fmglevel)%Imin-Ngh
  end function fmg_get_imingh
  function fmg_get_jmingh(fmglevel) result(jmingh)
    integer,intent(IN) :: fmglevel
    integer :: jmingh
    jmingh = GridSize(fmglevel)%Jmin-Ngh
  end function fmg_get_jmingh
  function fmg_get_kmingh(fmglevel) result(kmingh)
    integer,intent(IN) :: fmglevel
    integer :: kmingh
    kmingh = GridSize(fmglevel)%Kmin-Ngh
  end function fmg_get_kmingh
  function fmg_get_imaxgh(fmglevel) result(imaxgh)
    integer,intent(IN) :: fmglevel
    integer :: imaxgh
    imaxgh = GridSize(fmglevel)%Imax+Ngh
  end function fmg_get_imaxgh
  function fmg_get_jmaxgh(fmglevel) result(jmaxgh)
    integer,intent(IN) :: fmglevel
    integer :: jmaxgh
    jmaxgh = GridSize(fmglevel)%Jmax+Ngh
  end function fmg_get_jmaxgh
  function fmg_get_kmaxgh(fmglevel) result(kmaxgh)
    integer,intent(IN) :: fmglevel
    integer :: kmaxgh
    kmaxgh = GridSize(fmglevel)%Kmax+Ngh
  end function fmg_get_kmaxgh
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
  subroutine fmg_left_or_right(gid, rank, amr_level, lri, lrj, lrk)
    integer,intent(IN) :: gid, rank, amr_level
    integer,intent(OUT) :: lri, lrj, lrk
    integer :: i,j,k, pgid, prank
    lri = Undefi
    lrj = Undefi
    lrk = Undefi
    pgid = Ancestry(amr_level,rank)%Block(gid)%ParentGid
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
  function fmg_get_h(amrlev,fmglev) result(h)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=8) :: h
    h = GridLevel(amrlev,fmglev)%h
  end function fmg_get_h
  function fmg_get_ds(amrlev,fmglev) result(ds)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=8) :: ds
    ds = GridLevel(amrlev,fmglev)%ds
  end function fmg_get_ds
  function fmg_get_dv(amrlev,fmglev) result(dv)
    integer,intent(IN) :: amrlev, fmglev
    real(kind=8) :: dv
    dv = GridLevel(amrlev,fmglev)%dv
  end function fmg_get_dv
  function fmg_skip_grid(gid, amrlev, fmglev) result(bool)
    use mpilib
    integer,intent(IN) :: gid, amrlev, fmglev
    logical :: bool
    myrank = get_myrank()
    bool = amrlev /= AMR_LevelMax &
         .and. &
         Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi
  end function fmg_skip_grid
  function fmg_have_child(gid, amrlev) result(bool)
    integer,intent(IN) :: gid, amrlev
    logical :: bool
    bool = Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi
  end function fmg_have_child
  subroutine fmg_max(maxarr, amrlev, fmglev, icode, absolute, skip, xskip, mpireduce)
    real(kind=8),intent(OUT) :: maxarr
    integer,intent(IN) :: amrlev, fmglev, icode
    logical,intent(IN),optional :: absolute, skip, xskip, mpireduce
    logical :: bool_abs, bool_skip, bool_xskip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    real(kind=8) :: maxlocal
    bool_abs = .FALSE.
 if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.
 if (present(skip)) bool_skip = skip
    bool_xskip = .FALSE.
 if (present(xskip)) bool_xskip = xskip
    bool_mpi = .TRUE.
 if (present(mpireduce)) bool_mpi = mpireduce
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
          maxlocal = max(maxlocal, MAXVAL(abs(arrp(is:ie,js:je,ks:ke,:))))
       else
          maxlocal = max(maxlocal, MAXVAL(arrp(is:ie,js:je,ks:ke,:)))
       end if
    end do
    if (bool_mpi) then
       call mpi_allreduce(maxlocal, maxarr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    else
       maxarr = maxlocal
    end if
  end subroutine fmg_max
  function fmg_get_max(amrlev, fmglev, icode, absolute, skip, mpireduce) result(maxarr)
    integer,intent(IN) :: amrlev, fmglev, icode
    logical,intent(IN),optional :: absolute, skip, mpireduce
    real(kind=8) :: maxarr
    call fmg_max(maxarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_max
  subroutine fmg_sum(sumarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=8),intent(OUT) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    logical :: bool_abs, bool_skip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    real(kind=8) :: dv, sumlocal
    bool_abs = .FALSE.
 if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.
 if (present(skip)) bool_skip = skip
    bool_mpi = .TRUE.
 if (present(mpireduce)) bool_mpi = mpireduce
    dv = fmg_get_dv(amrlev,fmglev)
    sumlocal = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       arrp => fmg_get_arrp(amrlev, fmglev, gid, icode)
       if (bool_abs) then
          sumlocal = sumlocal + SUM(abs(arrp(is:ie,js:je,ks:ke,:)))
       else
          sumlocal = sumlocal + SUM(arrp(is:ie,js:je,ks:ke,:))
       end if
    end do
    if (bool_mpi) then
       call mpi_allreduce(sumlocal, sumarr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    else
       sumarr = sumlocal
    end if
    sumarr = sumarr * dv
  end subroutine fmg_sum
  function fmg_get_sum(amrlev, fmglev, icode, absolute, skip, mpireduce) result(sumarr)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=8) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_sum(sumarr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_sum
  subroutine fmg_ave(avearr, amrlev, fmglev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=8),intent(OUT) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    logical :: bool_abs, bool_skip, bool_mpi
    integer :: gid
    integer :: is, ie, js, je, ks, ke, ngrid
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    real(kind=8) :: avelocal
    bool_abs = .FALSE.
 if (present(absolute)) bool_abs = absolute
    bool_skip = .FALSE.
 if (present(skip)) bool_skip = skip
    bool_mpi = .TRUE.
 if (present(mpireduce)) bool_mpi = mpireduce
    avelocal = 0.d0
    ngrid = 0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev)) cycle
       ngrid = ngrid + 1
       arrp => fmg_get_arrp(amrlev, fmglev, gid, icode)
       if (bool_abs) then
          avelocal = avelocal + SUM(abs(arrp(is:ie,js:je,ks:ke,:)))
       else
          avelocal = avelocal + SUM(arrp(is:ie,js:je,ks:ke,:))
       end if
    end do
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
  end subroutine fmg_ave
  function fmg_get_ave(amrlev, fmglev, icode, absolute, skip, mpireduce) result(avearr)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=8) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_ave(avearr, amrlev, fmglev, icode, absolute, skip, mpireduce)
  end function fmg_get_ave
  function fmg_isVector(icode) result(bool)
    integer,intent(IN) :: icode
    logical :: bool
    bool = FMG_bool_isVector(icode)
  end function fmg_isVector
  subroutine fmg_setVector(icode, bool)
    integer,intent(IN) :: icode
    logical,intent(IN) :: bool
    FMG_bool_isVector(icode) = bool
  end subroutine fmg_setVector
  function fmg_getNcomp( icode ) result( ncomp )
    integer,intent(IN) :: icode
    integer :: ncomp
    if ( fmg_isVector(icode) ) then
       ncomp = Mmax - Mmin + 1
    else
       ncomp = 1
    end if
  end function fmg_getNcomp
  function fmg_get_dtime() result( dt )
    real(kind=8) :: dt
    dt = Dtime
  end function fmg_get_dtime
  subroutine fmg_set_dtime(dt)
    real(kind=8),intent(IN) :: dt
    Dtime = dt
  end subroutine fmg_set_dtime
end module fmg_data
