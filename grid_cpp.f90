module grid
  use mpilib
  implicit none
  integer,parameter :: Undefi=-1 
  real(kind=8),parameter :: Undefd=1.2345678d0 
  real(kind=8),parameter :: Undefd0=0.d0 
  integer,parameter :: Left = 0, Right = 1
  integer,parameter :: COMPLETE = 0
  integer,parameter :: PREDICTOR = 1
  integer,parameter :: CORRECTOR = 2
  integer,parameter :: Imin=0, Jmin=0, Kmin=0, Mmin=0, Gidmin=0, Lmin = 0
  integer,parameter :: Imax=8 -1, Jmax=8 -1, Kmax=8 -1, Mmax=34 -1, Gidmax=5000 -1, Lmax = 16 -1
  integer,parameter :: Ngh = 2 
  integer,parameter :: Imingh=Imin-Ngh, Jmingh=Jmin-Ngh, Kmingh=Kmin-Ngh
  integer,parameter :: Imaxgh=Imax+Ngh, Jmaxgh=Jmax+Ngh, Kmaxgh=Kmax+Ngh
  type t_blockMem
     real(kind=8),dimension(:,:,:,:),pointer :: u => null()
     real(kind=8),dimension(:),pointer :: x => null()
     real(kind=8),dimension(:),pointer :: y => null()
     real(kind=8),dimension(:),pointer :: z => null()
  end type t_blockMem
  type(t_blockMem),save,dimension(:),allocatable :: BlockMem
  logical,save,dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: GridMask
  integer,save,dimension(Gidmin:Gidmax) :: Igrid=Undefi, Jgrid=Undefi, Kgrid=Undefi 
  logical,save,dimension(Gidmin:Gidmax) :: AllocMem=.FALSE.
  integer,save,dimension(Gidmin:Gidmax) :: Levels=Undefi
  integer,save,dimension(Gidmin:Gidmax) :: U1list=Undefi, U2list=Undefi
  integer,save,dimension(0:1,0:2,Gidmin:Gidmax,0:400 -1) :: NeighborGid=Undefi, NeighborRank=MPI_PROC_NULL
  integer,save,dimension(0:1,0:1,0:1,Gidmin:Gidmax,0:400 -1) :: ChildGid=Undefi, ChildRank=MPI_PROC_NULL
  integer,save,dimension(Gidmin:Gidmax,0:400 -1) :: ParentGid=Undefi, ParentRank=MPI_PROC_NULL
  integer,save,dimension(0:(8)-1,0:(8)-1,0:(8)-1) :: GidBase, RankBase
  integer,save :: LevelMax = Undefi 
  real(kind=8),save,dimension(Lmin:Lmax) :: Time 
  real(kind=8),save,dimension(Lmin:Lmax) :: Dtime 
  integer(kind=8),save,dimension(Lmin:Lmax) :: Step 
  integer(kind=8),save,dimension(Lmin:Lmax) :: Dstep 
  real(kind=8),save,dimension(0:2,Lmin:Lmax) :: CellWidth 
  integer,save,dimension(Gidmin:Gidmax,Lmin:Lmax) :: GidList = Undefi 
  integer,save,dimension(Lmin:Lmax) :: GidListMax = Undefi 
  integer,save,dimension(Gidmin:Gidmax,Lmin:Lmax,0:400 -1) :: GidListNode = Undefi 
  integer,save,dimension(Lmin:Lmax,0:400 -1) :: GidListNodeMax = Undefi 
  integer(kind=8),save,dimension(Lmin:Lmax) :: U_StepNumber = 0
  integer(kind=8),save,dimension(Lmin:Lmax) :: U_StepNumberGhostCell = 0
  integer(kind=8),save,dimension(Lmin:Lmax) :: U1_StepNumber = 0
  integer(kind=8),save,dimension(Lmin:Lmax) :: U1_StepNumberGhostCell = 0
  integer(kind=8),save,dimension(Lmin:Lmax) :: U2_StepNumber = 0
  integer(kind=8),save,dimension(Lmin:Lmax) :: U2_StepNumberGhostCell = 0
  integer(kind=8),save ::globdbg_flag = 0,globdbg_mygid = 0,globdbg_myrank = 0
  integer(kind=8),parameter:: globdbg_rank = 0, globdbg_gid = 6465, globdbg_i = 3, globdbg_j = 2, globdbg_k = 2
  logical,save :: globdbg_gridupdate
contains
  subroutine grid_init
    GridMask = .FALSE.
    GridMask(Imin:Imax,Jmin:Jmax,Kmin:Kmax) = .TRUE.
    allocate(BlockMem(Gidmin:Gidmax))
  end subroutine grid_init
  subroutine grid_finalize
    integer :: gid
    do gid=lbound(BlockMem,1), ubound(BlockMem,1)
       if (.not. AllocMem(gid)) cycle
       call dealloc_U(gid)
       call dealloc_U1order(gid)
       call dealloc_U2order(gid)
    enddo
    deallocate(BlockMem)
  end subroutine grid_finalize
  function get_level(gid) result( level )
    integer,intent(IN) :: gid
    integer :: level
    level = Levels(gid)
  end function get_level
  function get_Xp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%x
  end function get_Xp
  function get_Yp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%y
  end function get_Yp
  function get_Zp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%z
  end function get_Zp
  function get_Ucomp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=8),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(gid)%u(:,:,:,ncomp))
  end function get_Ucomp
  function slice3(arr) result(ptra)
    real(kind=8),intent(IN),dimension(Imingh:,Jmingh:,Kmingh:),target :: arr
    real(kind=8),dimension(:,:,:),pointer :: ptra
    ptra => arr
  end function slice3
  function get_Up(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(gid)%u
  end function get_Up
  function get_U1orderp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(U1list(gid))%u
  end function get_U1orderp
  function get_U1comp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=8),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(U1list(gid))%u(:,:,:,ncomp))
  end function get_U1comp
  function get_U2orderp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=8),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(U2list(gid))%u
  end function get_U2orderp
  function get_U2comp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=8),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(U2list(gid))%u(:,:,:,ncomp))
  end function get_U2comp
  function alloc_U(level) result(gid)
    integer,intent(IN) :: level
    integer :: gid
    gid = alloc_block()
    Levels( gid ) = level
  end function alloc_U
  subroutine alloc_U1order(gid)
    integer,intent(IN) :: gid
    U1list(gid) = alloc_block(bool_omit_coordinates=.TRUE.)
  end subroutine alloc_U1order
  subroutine alloc_U2order(gid)
    integer,intent(IN) :: gid
    U2list(gid) = alloc_block(bool_omit_coordinates=.TRUE.)
  end subroutine alloc_U2order
  function alloc_block(bool_omit_coordinates) result(gid)
    logical,intent(IN),optional :: bool_omit_coordinates
    integer :: gid, id, l
    gid = Undefi
    do id = lbound(AllocMem,1), ubound(AllocMem,1)
       if ( .not. AllocMem(id) ) then
          gid = id 
          call alloc_block_by_gid(gid, bool_omit_coordinates)
          return
       endif
    enddo
    print *, '*** gid reaches Gidmax.  Increase NGID in config.h.'
    do l = Lmin, LevelMax
       print *, get_myrank(), l, count( AllocMem )
    enddo
    stop
  end function alloc_block
  subroutine alloc_block_by_gid(gid, bool_omit_coordinates)
    integer,intent(IN) :: gid
    logical,intent(IN),optional :: bool_omit_coordinates
    logical :: bool_omit_coord
    integer :: err
    allocate(BlockMem(gid)%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax), stat=err)
    if (err /= 0) then
       print *, '*** Error in allocating BlockMem%u', get_myrank(), gid
       stop
    endif
    AllocMem(gid) = .TRUE.
    bool_omit_coord = .FALSE. 
    if (present(bool_omit_coordinates)) bool_omit_coord = bool_omit_coordinates
    if (bool_omit_coord) return
    allocate(BlockMem(gid)%x(Imingh:Imaxgh), stat=err)
    if (err /= 0) then
       print *, '*** Error in allocating BlockMem%x', get_myrank(), gid
       stop
    endif
    allocate(BlockMem(gid)%y(Jmingh:Jmaxgh), stat=err)
    if (err /= 0) then
       print *, '*** Error in allocating BlockMem%y', get_myrank(), gid
       stop
    endif
    allocate(BlockMem(gid)%z(Kmingh:Kmaxgh), stat=err)
    if (err /= 0) then
       print *, '*** Error in allocating BlockMem%z', get_myrank(), gid
       stop
    endif
  end subroutine alloc_block_by_gid
  subroutine dealloc_U(gid)
    integer,intent(IN) :: gid
    integer :: pgid, prank
    call dealloc_block( gid )
    myrank = get_myrank()
    pgid = ParentGid(gid,myrank)
    prank = ParentRank(gid,myrank)
    Igrid(gid)=Undefi
    Jgrid(gid)=Undefi
    Kgrid(gid)=Undefi
    Levels(gid)=Undefi
    NeighborGid(:,:,gid,myrank)=Undefi
    NeighborRank(:,:,gid,myrank)=MPI_PROC_NULL
    ChildGid(:,:,:,gid,myrank)=Undefi
    ChildRank(:,:,:,gid,myrank)=MPI_PROC_NULL
    ParentGid(gid,myrank)=Undefi
    ParentRank(gid,myrank)=MPI_PROC_NULL
  end subroutine dealloc_U
  subroutine dealloc_U1order(gid)
    integer,intent(IN) :: gid
    call dealloc_block( U1list(gid) )
    U1list(gid) = Undefi
  end subroutine dealloc_U1order
  subroutine dealloc_U2order(gid)
    integer,intent(IN) :: gid
    call dealloc_block( U2list(gid) )
    U2list(gid) = Undefi
  end subroutine dealloc_U2order
  subroutine dealloc_block(gid)
    integer,intent(IN) :: gid
    if (gid == Undefi) return
    if (associated(BlockMem(gid)%u)) deallocate(BlockMem(gid)%u) 
 nullify(BlockMem(gid)%u)
    if (associated(BlockMem(gid)%x)) deallocate(BlockMem(gid)%x) 
 nullify(BlockMem(gid)%x)
    if (associated(BlockMem(gid)%y)) deallocate(BlockMem(gid)%y) 
 nullify(BlockMem(gid)%y)
    if (associated(BlockMem(gid)%z)) deallocate(BlockMem(gid)%z) 
 nullify(BlockMem(gid)%z)
    AllocMem( gid ) = .FALSE.
  end subroutine dealloc_block
  subroutine update_gidlist(level)
    integer,intent(IN) :: level
    integer :: id, n
    GidList(:,level) = Undefi
    n = -1
    do id = lbound(Levels,1), ubound(Levels,1)
       if ( Levels(id) == level ) then
          n = n + 1
          GidList(n,level) = id
       endif
    enddo
    GidListMax( level ) = n
    call update_gidlist_node(level)
  end subroutine update_gidlist
  subroutine update_gidlist_node(level)
    integer,intent(IN) :: level
    integer :: rank, last
    integer,dimension(Gidmin:Gidmax) :: glist
    myrank = get_myrank()
    GidListNode(:,level,:) = Undefi
    do rank = 0, 400 -1
       if (rank == myrank) then
          last = GidListMax( level )
          glist(Gidmin:last) = GidList(Gidmin:last, level )
       endif
       call mpi_bcast(last, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       call mpi_bcast(glist, last-Gidmin+1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       GidListNodeMax(level, rank) = last
       GidListNode(Gidmin:last, level, rank) = glist(Gidmin:last)
    enddo
  end subroutine update_gidlist_node
  function grid_used_num() result(num)
    integer :: num, id
    num = 0
    do id = lbound(Levels,1), ubound(Levels,1)
       if ( Levels(id) /= Undefi ) then
          num = num + 1
       endif
    enddo
  end function grid_used_num
  function grid_used_num_global() result(num)
    integer :: num, num_node
    num_node = grid_used_num()
    call mpi_allreduce(num_node, num, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  end function grid_used_num_global
  function grid_used_num_level( level ) result(num)
    integer,intent(IN) :: level
    integer :: num
    num = GidListMax(level) - Gidmin + 1
  end function grid_used_num_level
  function grid_used_num_level_global( level ) result(num)
    integer,intent(IN) :: level
    integer :: num, num_node
    num_node = grid_used_num_level( level )
    call mpi_allreduce(num_node, num, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  end function grid_used_num_level_global
  function get_dv(level) result(dv)
    integer,intent(IN) :: level
    real(kind=8) :: dv
    dv = CellWidth(0,level) * CellWidth(1,level) * CellWidth(2,level)
  end function get_dv
  function get_ds(level) result(ds)
    integer,intent(IN) :: level
    real(kind=8) :: ds(0:2)
    ds = (/ &
         CellWidth(1,level) * CellWidth(2,level) , &
         CellWidth(2,level) * CellWidth(0,level) , &
         CellWidth(0,level) * CellWidth(1,level) /)
  end function get_ds
  subroutine left_or_right(gid, rank, lri, lrj, lrk)
    integer,intent(IN) :: gid, rank
    integer,intent(OUT) :: lri, lrj, lrk
    integer :: i,j,k, pgid, prank
    lri = Undefi
    lrj = Undefi
    lrk = Undefi
    pgid = ParentGid(gid,rank)
    prank = ParentRank(gid,rank)
    if ( pgid == Undefi .or. prank == MPI_PROC_NULL ) then
       print *, '*** invarid link list in left_or_right'
       stop
    endif
    do k = Left,Right
       do j = Left, Right
          do i = Left, Right
             if ( ChildGid(i,j,k,pgid,prank) == gid .and. &
                  ChildRank(i,j,k,pgid,prank) == rank ) then
                lri = i
                lrj = j
                lrk = k
             endif
          enddo
       enddo
    enddo
    if (lri < 0) then
       print *, '*** inconsistent link list'
       stop
    endif
  end subroutine left_or_right
  function encode_surf(lr, direc) result( code )
    integer,intent(IN) :: lr, direc
    integer :: code
    code = lr + direc * (2 -0)
  end function encode_surf
  function decode_direction(ncord) result (direc)
    integer,intent(IN) :: ncord
    integer :: direc
    direc = ncord / (2 -0)
  end function decode_direction
  function decode_LR(ncord) result (surf)
    integer,intent(IN) :: ncord
    integer :: surf
    surf = modulo(ncord , 2)
  end function decode_LR
  function swap_LR(surf) result(swap)
    integer,intent(IN) :: surf
    integer :: swap
    swap = modulo( surf+1, 2 )
  end function swap_LR
  subroutine convert_coord(gidfrom, gidto, if,jf,kf, it,jt,kt)
    integer,intent(IN) :: gidfrom, gidto, if,jf,kf
    integer,intent(OUT) :: it,jt,kt
    integer(kind=8) :: ifoffset,itoffset,jfoffset,jtoffset,kfoffset,ktoffset
    ifoffset = Igrid(gidfrom)*(Imax-Imin+1)
    itoffset = Igrid(gidto)*(Imax-Imin+1)
    jfoffset = Jgrid(gidfrom)*(Jmax-Jmin+1)
    jtoffset = Jgrid(gidto)*(Jmax-Jmin+1)
    kfoffset = Kgrid(gidfrom)*(Kmax-Kmin+1)
    ktoffset = Kgrid(gidto)*(Kmax-Kmin+1)
    if ( Levels(gidfrom) == Levels(gidto) ) then 
       it = if + ifoffset - itoffset
       jt = jf + jfoffset - jtoffset
       kt = kf + kfoffset - ktoffset
    elseif ( Levels(gidfrom) == Levels(gidto) + 1) then 
       it = ((if + ifoffset)-(0))/2 + mod(min((if + ifoffset)-(0),0),2) + (0) - itoffset
       jt = ((jf + jfoffset)-(0))/2 + mod(min((jf + jfoffset)-(0),0),2) + (0) - jtoffset
       kt = ((kf + kfoffset)-(0))/2 + mod(min((kf + kfoffset)-(0),0),2) + (0) - ktoffset
    elseif ( Levels(gidfrom) == Levels(gidto) - 1) then 
       it = ((if + ifoffset) - (0)) * 2 + (0) - itoffset
       jt = ((jf + jfoffset) - (0)) * 2 + (0) - jtoffset
       kt = ((kf + kfoffset) - (0)) * 2 + (0) - ktoffset
    else
       print *, 'It is not supported yet'
       stop
    endif
  end subroutine convert_coord
  subroutine convert_coord_cb(gidfrom, gidto, surf, if,jf,kf, it,jt,kt)
    integer,intent(IN) :: gidfrom, gidto, if,jf,kf, surf
    integer,intent(OUT) :: it,jt,kt
    integer(kind=8) :: ifoffset,itoffset,jfoffset,jtoffset,kfoffset,ktoffset
    call convert_coord(gidfrom, gidto, if,jf,kf, it,jt,kt)
    if ( Levels(gidfrom) == Levels(gidto) - 1) then 
       if ( surf == 0 ) it = it + 1
       if ( surf == 1 ) jt = jt + 1
       if ( surf == 2 ) kt = kt + 1
    endif
  end subroutine convert_coord_cb
  subroutine get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank, whichlevel)
    integer,intent(IN) :: ig, jg, kg, level
    integer,intent(IN),optional :: whichlevel
    integer,intent(OUT) :: gid, rank
    integer :: lri, lrj, lrk, npos, gidn, rankn, ig0, jg0, kg0, descend
    if ( present( whichlevel ) ) then
       descend = level - whichlevel
    else
       descend = 0
    endif
    if ( descend < 0 ) print *, '*** error in descend'
    ig0 = ishft(ig, -level)
    jg0 = ishft(jg, -level)
    kg0 = ishft(kg, -level)
    gid = GidBase(ig0, jg0, kg0)
    rank = RankBase(ig0, jg0, kg0)
    do npos = level-1, descend, -1
       if ( gid == Undefi ) exit
       if ( rank == MPI_PROC_NULL ) exit
       lri = ibits(ig,npos,1)
       lrj = ibits(jg,npos,1)
       lrk = ibits(kg,npos,1)
       gidn = ChildGid(lri, lrj, lrk, gid, rank)
       rankn = ChildRank(lri, lrj, lrk, gid, rank)
       gid = gidn
       rank = rankn
    enddo
  end subroutine get_gid_from_ijkgrid
  function has_child_grid(gid) result(bool)
    use mpilib
    integer,intent(IN) :: gid
    logical :: bool
    myrank = get_myrank()
    if (ChildGid(Left, Left, Left, gid, myrank) == Undefi) then
       bool = .FALSE.
    else
       bool = .TRUE.
    endif
  end function has_child_grid
  subroutine update_neighbor(gid, thisrank, ijkgrid, level)
    use mpilib
    integer,intent(IN) :: gid, thisrank, ijkgrid(0:2), level
    integer :: pgid, prank, pranklr, pgidlr, m
    integer,dimension(0:2) :: lra, lrb
    do m = 0, 2
       lra(m) = modulo(ijkgrid(m),2)
       lrb(m) = modulo(ijkgrid(m)+1,2)
    enddo
    pgid = ParentGid(gid, thisrank)
    prank = ParentRank(gid, thisrank)
    NeighborGid( lrb(0), 0, gid, thisrank ) = ChildGid( lrb(0), lra(1), lra(2), pgid, prank )
    NeighborGid( lrb(1), 1, gid, thisrank ) = ChildGid( lra(0), lrb(1), lra(2), pgid, prank )
    NeighborGid( lrb(2), 2, gid, thisrank ) = ChildGid( lra(0), lra(1), lrb(2), pgid, prank )
    NeighborRank( lrb(0), 0, gid, thisrank ) = ChildRank( lrb(0), lra(1), lra(2), pgid, prank )
    NeighborRank( lrb(1), 1, gid, thisrank ) = ChildRank( lra(0), lrb(1), lra(2), pgid, prank )
    NeighborRank( lrb(2), 2, gid, thisrank ) = ChildRank( lra(0), lra(1), lrb(2), pgid, prank )
    pgidlr = NeighborGid( lra(0), 0, pgid, prank )
    pranklr = NeighborRank( lra(0), 0, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(0), 0, gid, thisrank ) = ChildGid( lrb(0), lra(1), lra(2), pgidlr, pranklr )
       NeighborRank( lra(0), 0, gid, thisrank ) = ChildRank( lrb(0), lra(1), lra(2), pgidlr, pranklr )
    else
       NeighborGid( lra(0), 0, gid, thisrank ) = Undefi
       NeighborRank( lra(0), 0, gid, thisrank ) = MPI_PROC_NULL
    endif
    pgidlr = NeighborGid( lra(1), 1, pgid, prank )
    pranklr = NeighborRank( lra(1), 1, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(1), 1, gid, thisrank ) = ChildGid( lra(0), lrb(1), lra(2), pgidlr, pranklr )
       NeighborRank( lra(1), 1, gid, thisrank ) = ChildRank( lra(0), lrb(1), lra(2), pgidlr, pranklr )
    else
       NeighborGid( lra(1), 1, gid, thisrank ) = Undefi
       NeighborRank( lra(1), 1, gid, thisrank ) = MPI_PROC_NULL
    endif
    pgidlr = NeighborGid( lra(2), 2, pgid, prank )
    pranklr = NeighborRank( lra(2), 2, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(2), 2, gid, thisrank ) = ChildGid( lra(0), lra(1), lrb(2), pgidlr, pranklr )
       NeighborRank( lra(2), 2, gid, thisrank ) = ChildRank( lra(0), lra(1), lrb(2), pgidlr, pranklr )
    else
       NeighborGid( lra(2), 2, gid, thisrank ) = Undefi
       NeighborRank( lra(2), 2, gid, thisrank ) = MPI_PROC_NULL
    endif
  end subroutine update_neighbor
  subroutine update_neighbor_linklist(level)
    integer,intent(IN) :: level
    integer :: gid, n, rank, listmin, listmax
    integer,dimension(:),allocatable :: gidl
    integer,dimension(:,:),allocatable :: ijkgrid
    myrank = get_myrank()
    do rank = 0, 400 -1
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)
       allocate( gidl( listmin:listmax ) )
       allocate( ijkgrid( 0:2, listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level)
          do n = listmin, listmax
             ijkgrid(:,n) = (/ Igrid(gidl(n)), Jgrid(gidl(n)), Kgrid(gidl(n)) /)
          enddo
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       call mpi_bcast( ijkgrid, size(ijkgrid), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       do n = listmin, listmax
          call update_neighbor( gidl(n), rank, ijkgrid(:,n), level)
       enddo
       deallocate( gidl, ijkgrid )
    enddo
  end subroutine update_neighbor_linklist
  subroutine update_parent_linklist(level)
    integer,intent(IN) :: level
    integer :: n, cgid, pgid, crank, prank, rank, listmin, listmax, lri, lrj, lrk, ig,jg,kg
    integer,dimension(:),allocatable :: gidl
    integer,dimension(:,:),allocatable :: ijkgrid
    myrank = get_myrank()
    do rank = 0, 400 -1
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)
       allocate( gidl( listmin:listmax ) )
       allocate( ijkgrid( 0:2, listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level)
          do n = listmin, listmax
             ijkgrid(:,n) = (/ Igrid(gidl(n)), Jgrid(gidl(n)), Kgrid(gidl(n)) /)
          enddo
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       call mpi_bcast( ijkgrid, size(ijkgrid), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       do n = listmin, listmax
          cgid = gidl(n)
          crank = rank
          ig = ijkgrid(0,n)
          jg = ijkgrid(1,n)
          kg = ijkgrid(2,n)
          call get_gid_from_ijkgrid(ig,jg,kg,level, pgid,prank, whichlevel=level-1)
          if ( pgid == Undefi ) print *, '*** error pgid', pgid, prank, ig,jg,kg,level
          lri = modulo(ig,2)
          lrj = modulo(jg,2)
          lrk = modulo(kg,2)
          ChildGid(lri, lrj, lrk, pgid, prank) = cgid
          ChildRank(lri, lrj, lrk, pgid, prank) = crank
          ParentGid(cgid, crank) = pgid
          ParentRank(cgid, crank) = prank
       enddo
       deallocate( gidl, ijkgrid )
    enddo
  end subroutine update_parent_linklist
  subroutine clear_linklist(level)
    integer,intent(IN) :: level
    integer :: n, gid, cgid, pgid, crank, prank, rank, listmin, listmax, lri, lrj, lrk
    integer,dimension(:),allocatable :: gidl
    if ( level <= 0 ) return
    myrank = get_myrank()
    do rank = 0, 400 -1
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       do n = listmin, listmax
          gid = gidl(n)
          ParentGid(gid, rank) = Undefi
          ParentRank(gid, rank) = MPI_PROC_NULL
          ChildGid(:,:,:,gid,rank) = Undefi
          ChildRank(:,:,:,gid,rank) = MPI_PROC_NULL
          NeighborGid(:,:,gid,rank) = Undefi
          NeighborRank(:,:,gid,rank) = MPI_PROC_NULL
       enddo
       deallocate( gidl )
    enddo
    do rank = 0, 400 -1
       if ( rank == myrank ) &
            listmax = GidListMax( level-1 )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level-1)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       do n = listmin, listmax
          pgid = gidl(n)
          prank = rank
          ChildGid(:,:,:,pgid,prank) = Undefi
          ChildRank(:,:,:,pgid,prank) = MPI_PROC_NULL
       enddo
       deallocate( gidl )
    enddo
    if ( level + 1 > LevelMax ) return
    do rank = 0, 400 -1
       if ( rank == myrank ) &
            listmax = GidListMax( level+1 )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level+1)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       do n = listmin, listmax
          cgid = gidl(n)
          crank = rank
          ParentGid(cgid,crank) = Undefi
          ParentRank(cgid,crank) = MPI_PROC_NULL
       enddo
       deallocate( gidl )
    enddo
  end subroutine clear_linklist
  subroutine check_linklist(level)
    integer,intent(IN) :: level
    integer :: n, gid, pgid, prank, i,j,k, rank
    integer,dimension(Gidmin:Gidmax,0:400 -1) :: ptmp
    integer,dimension(0:1,0:1,0:1,Gidmin:Gidmax,0:400 -1) :: ctmp
    logical :: bool, boolt, boolr
    myrank = get_myrank()
    boolt = .TRUE.
    bool = .TRUE.
    do rank = 0, 400 -1
       ptmp = ParentGid
       call mpi_bcast(ptmp, size(ptmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ptmp-ParentGid)-minval(ptmp-ParentGid) /= 0) then
          print *, '** bad  consistency ParentGid, node=',rank
          bool = .false.
       end if
    enddo
    do rank = 0, 400 -1
       ptmp = ParentRank
       call mpi_bcast(ptmp, size(ptmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ptmp-ParentRank)-minval(ptmp-ParentRank) /= 0) then
          print *, '** bad  consistency RarentRank, node=',rank
          bool = .false.
       end if
    enddo
    do rank = 0, 400 -1
       ctmp = ChildGid
       call mpi_bcast(ctmp, size(ctmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ctmp-ChildGid)-minval(ctmp-ChildGid) /= 0) then
          print *, '** bad  consistency ChidlGid node=', rank
          bool = .false.
       end if
    enddo
    do rank = 0, 400 -1
       ctmp = ChildRank
       call mpi_bcast(ctmp, size(ctmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ctmp-ChildRank)-minval(ctmp-ChildRank) /= 0) then
          print *, '** bad  consistency ChildRank node =', rank
          bool = .false.
       end if
    enddo
    boolt = boolt .and. bool
    if ( level > Lmin ) then
       bool = .true.
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          if ( ParentGid(gid, myrank) == Undefi ) then
             print *, 'brake link ParentGid', gid, myrank
             bool = .false.
          end if
          if ( ParentRank(gid, myrank) == MPI_PROC_NULL ) then
             print *, 'brake link ParentRank', gid, myrank
             bool = .false.
          end if
       enddo
       if ( bool ) print *, myrank, level, 'OK: child -> parent'
       boolt = boolt .and. bool
    endif
    if ( level > Lmin .and. level <= LevelMax ) then
       bool = .true.
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          pgid = ParentGid(gid, myrank)
          prank = ParentRank(gid, myrank)
          k = modulo(Kgrid(gid),2)
          j = modulo(Jgrid(gid),2)
          i = modulo(Igrid(gid),2)
          if ( ChildGid(i,j,k,pgid,prank) /= gid ) then
             print *, 'brake link ChildGid', gid, myrank, ChildGid(i,j,k,pgid,prank)
             bool = .false.
          end if
          if ( ChildRank(i,j,k,pgid,prank) /= myrank ) then
             print *, 'brake link ChildRank', gid, myrank, ChildRank(i,j,k,pgid,prank)
             bool = .false.
          end if
       enddo
       if ( bool ) print *, myrank, level, 'OK: parent -> child'
       boolt = boolt .and. bool
    endif
    bool = .true.
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       call get_gid_from_ijkgrid(Igrid(gid),Jgrid(gid),Kgrid(gid),level,pgid,prank)
       if ( pgid /= gid .or. prank /= myrank ) then
          print *, 'brake link parent - child', gid, myrank, Igrid(gid),Jgrid(gid),Kgrid(gid)
          bool = .false.
       end if
    enddo
    if ( bool ) print *, myrank, level, 'OK: parent <-> child'
    boolt = boolt .and. bool
    call mpi_allreduce(boolt, boolr, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr )
    if ( .not. boolr ) then
       call mpi_finalize(ierr)
       print *, '*** stop'
       stop
    endif
  end subroutine check_linklist
  function level_sync() result(levsync)
    integer :: l, levsync, id
    integer(kind=8) :: steplmax
    id = 0
    steplmax = Step( LevelMax )
    levsync = LevelMax
    do l = LevelMax-1, Lmin, -1
       if ( Step( l ) == steplmax ) levsync = l
    enddo
  end function level_sync
  function get_current_level() result(level)
    integer :: lev, level
    level = LevelMax
    do lev = Lmin, LevelMax-1
       if (Step(lev) - Dstep(lev) == Step(lev+1)) then
          level = lev
          exit
       end if
    end do
  end function get_current_level
end module grid
