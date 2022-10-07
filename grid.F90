#include "config.h"
!-------------------------------------------------------------------------
! Module for grid configuration and operation
!-------------------------------------------------------------------------
module grid
  use mpilib
  implicit none
  !  private
  integer,parameter :: Undefi=-1 !未定義値（初期化時に代入）
  real(kind=DBL_KIND),parameter :: Undefd=1.2345678d0 !未定義値（初期化時に代入）
  real(kind=DBL_KIND),parameter :: Undefd0=0.d0 !未定義値（初期化時に代入）
  integer,parameter :: Left = 0, Right = 1
  ! STEP_MODE
  integer,parameter :: COMPLETE  = 0
  integer,parameter :: PREDICTOR = 1
  integer,parameter :: CORRECTOR = 2
  ! grid size
  integer,parameter :: Imin=0, Jmin=0, Kmin=0, Mmin=0, Gidmin=0, Lmin = 0
  integer,parameter :: Imax=NI-1, Jmax=NJ-1, Kmax=NK-1, Mmax=NM-1, Gidmax=NGID-1, Lmax = NL-1
  integer,parameter :: Ngh = N_GHOST_CELL  ! number of ghost cell
  ! grid size including ghost cells
  integer,parameter :: Imingh=Imin-Ngh, Jmingh=Jmin-Ngh, Kmingh=Kmin-Ngh
  integer,parameter :: Imaxgh=Imax+Ngh, Jmaxgh=Jmax+Ngh, Kmaxgh=Kmax+Ngh
  ! arrays for data, BlockMem(gid)%u(i,j,k,m), BlockMem(gid)%x(i), BlockMem(gid)%y(j), BlockMem(gid)%z(k)
  type t_blockMem
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u => null()
     real(kind=DBL_KIND),dimension(:),pointer :: x => null()
     real(kind=DBL_KIND),dimension(:),pointer :: y => null()
     real(kind=DBL_KIND),dimension(:),pointer :: z => null()
  end type t_blockMem
  type(t_blockMem),save,dimension(:),allocatable :: BlockMem
  ! gridmask
  logical,save,dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: GridMask
  ! ------------------------------
  !  property of grid (node local)
  ! ------------------------------
  integer,save,dimension(Gidmin:Gidmax) :: Igrid=Undefi, Jgrid=Undefi, Kgrid=Undefi ! coordinates of grid, key=gid
  ! The grid is free or used
  logical,save,dimension(Gidmin:Gidmax) :: AllocMem=.FALSE.
  ! list of level for U, key=gid
  integer,save,dimension(Gidmin:Gidmax) :: Levels=Undefi
  ! list of pointers indicating 1st order and 2nd order solutions, key=gid.
  integer,save,dimension(Gidmin:Gidmax) :: U1list=Undefi, U2list=Undefi
  ! -------------------------------
  !  property of grid (node global)
  ! -------------------------------
  integer,save,dimension(0:1,MX:MZ,Gidmin:Gidmax,0:NPE-1) :: NeighborGid=Undefi, NeighborRank=MPI_PROC_NULL
  ! NeighborGid( L/R, ijk-directions, grid id, rank )
  ! 斜めは i-方向転送 → j-方向転送 → k-方向転送 の順に行うと, 転送されるので、覚えておかなくて大丈夫
  integer,save,dimension(0:1,0:1,0:1,Gidmin:Gidmax,0:NPE-1) :: ChildGid=Undefi, ChildRank=MPI_PROC_NULL
  ! leftgid( i-direction, j-direction, k-direction, grid id, rank )
  integer,save,dimension(Gidmin:Gidmax,0:NPE-1) :: ParentGid=Undefi, ParentRank=MPI_PROC_NULL
  ! base grid
  integer,save,dimension(0:(NGI_BASE)-1,0:(NGJ_BASE)-1,0:(NGK_BASE)-1) :: GidBase, RankBase
  ! ----------------------
  ! property of grid level
  ! ----------------------
  integer,save :: LevelMax = Undefi          ! maximum grid level
  real(kind=DBL_KIND),save,dimension(Lmin:Lmax) :: Time        ! time
  real(kind=DBL_KIND),save,dimension(Lmin:Lmax) :: Dtime       ! delta t
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: Step        ! step number
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: Dstep       ! delta step number
  real(kind=DBL_KIND),save,dimension(MX:MZ,Lmin:Lmax) :: CellWidth  ! cell width
  integer,save,dimension(Gidmin:Gidmax,Lmin:Lmax) :: GidList = Undefi ! 自ノードの gid のリスト
  integer,save,dimension(Lmin:Lmax) :: GidListMax = Undefi !上記リストの最大要素番号
  integer,save,dimension(Gidmin:Gidmax,Lmin:Lmax,0:NPE-1) :: GidListNode = Undefi ! GidList の各ノード情報
  integer,save,dimension(Lmin:Lmax,0:NPE-1) :: GidListNodeMax = Undefi ! GidListMax の各ノード情報

  ! ---------------------------------------------------------------------
  ! block 内部と 袖のステップ数をカウントする変数(半ステップで1増加する)
  ! ---------------------------------------------------------------------
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U_StepNumber = 0
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U_StepNumberGhostCell = 0
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U1_StepNumber = 0
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U1_StepNumberGhostCell = 0
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U2_StepNumber = 0
  integer(kind=LLONG_KIND),save,dimension(Lmin:Lmax) :: U2_StepNumberGhostCell = 0
  ! level に含まれる gid のリスト。key1 = No, key2 = grid level。node local
  ! refinement ごとに更新されるべきもの。 levels(gid) から作られる。(左詰め)
  !  public :: Time, Dtime, Step, Dstep, CellWidth, GidList

  !---------- KS DEBUG --------------!
  integer(kind=LLONG_KIND),save ::globdbg_flag = 0,globdbg_mygid = 0,globdbg_myrank = 0
  !globdbg_i=i-lbound(u,1), etc.
  integer(kind=LLONG_KIND),parameter:: globdbg_rank = 0, globdbg_gid = 6465, globdbg_i = 3, globdbg_j = 2, globdbg_k = 2
  logical,save :: globdbg_gridupdate
  !---------- KS DEBUG --------------!

contains
  !-------------------------------------------------------------------------
  ! grid init
  !-------------------------------------------------------------------------
  subroutine grid_init
    GridMask = .FALSE.
    GridMask(ARRAYSIZE_IJK) = .TRUE.
    allocate(BlockMem(Gidmin:Gidmax))
  end subroutine grid_init
  !-------------------------------------------------------------------------
  ! grid finalize
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! get levels given by gid
  !-------------------------------------------------------------------------
  function get_level(gid) result( level )
    integer,intent(IN) :: gid
    integer :: level
    level = Levels(gid)
  end function get_level
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、X の pointer を返す。ex) x => get_Xp(gid)
  !-------------------------------------------------------------------------
  function get_Xp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%x
  end function get_Xp
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、Y の pointer を返す。ex) y => get_Yp(gid)
  !-------------------------------------------------------------------------
  function get_Yp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%y
  end function get_Yp
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、Z の pointer を返す。ex) y => get_Zp(gid)
  !-------------------------------------------------------------------------
  function get_Zp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:),pointer :: ptr
    ptr => BlockMem(gid)%z
  end function get_Zp
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、U の pointer を返す。ex) rho => get_Ucomp(MRHO,gid)
  !-------------------------------------------------------------------------
  function get_Ucomp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(gid)%u(:,:,:,ncomp))
  end function get_Ucomp
  !-------------------------------------------------------------------------
  ! 添字の下限値を指定して、引数配列を受け取る。ポインタを返す。
  !-------------------------------------------------------------------------
  function slice3(arr) result(ptra)
    real(kind=DBL_KIND),intent(IN),dimension(Imingh:,Jmingh:,Kmingh:),target :: arr
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: ptra
    ptra => arr
  end function slice3
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、U の pointer を返す。
  !-------------------------------------------------------------------------
  function get_Up(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(gid)%u
  end function get_Up
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、1st order solution の pointer を返す。
  !-------------------------------------------------------------------------
  function get_U1orderp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(U1list(gid))%u
  end function get_U1orderp
  !-------------------------------------------------------------------------
  function get_U1comp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(U1list(gid))%u(:,:,:,ncomp))
  end function get_U1comp
  !-------------------------------------------------------------------------
  ! gid が与えられたとき、2nd order solution の pointer を返す。
  !-------------------------------------------------------------------------
  function get_U2orderp(gid) result (ptr)
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptr
    ptr => BlockMem(U2list(gid))%u
  end function get_U2orderp
  !-------------------------------------------------------------------------
  function get_U2comp(ncomp, gid) result (ptr)
    integer,intent(IN) :: gid, ncomp
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: ptr
    ptr => slice3(BlockMem(U2list(gid))%u(:,:,:,ncomp))
  end function get_U2comp
  !-------------------------------------------------------------------------
  ! allocate memory for U
  !-------------------------------------------------------------------------
  function alloc_U(level) result(gid)
    integer,intent(IN) :: level
    integer :: gid
    gid = alloc_block()
    Levels( gid ) = level
  end function alloc_U
  !-------------------------------------------------------------------------
  ! allocate memory for 1st and 2nd order solutions
  !-------------------------------------------------------------------------
  subroutine alloc_U1order(gid)
    integer,intent(IN) :: gid
    U1list(gid) = alloc_block(bool_omit_coordinates=.TRUE.)
  end subroutine alloc_U1order
  !-------------------------------------------------------------------------
  subroutine alloc_U2order(gid)
    integer,intent(IN) :: gid
    U2list(gid) = alloc_block(bool_omit_coordinates=.TRUE.)
  end subroutine alloc_U2order
  !-------------------------------------------------------------------------
  ! メモリを確保し、gid を返す。
  ! bool_omit_coordinates が真のとき、座標の変数を割りつけるのを省略する。
  ! デフォルトは偽である（つまり座標も割りつける）。
  !-------------------------------------------------------------------------
  function alloc_block(bool_omit_coordinates) result(gid)
    logical,intent(IN),optional :: bool_omit_coordinates
    integer :: gid, id, l
    ! search free memory
    gid = Undefi
    do id = lbound(AllocMem,1), ubound(AllocMem,1)
       if ( .not. AllocMem(id) ) then
          gid = id              ! returned value
          call alloc_block_by_gid(gid, bool_omit_coordinates)
          return
       endif
    enddo

    ! if error
    print *, '*** gid reaches Gidmax.  Increase NGID in config.h.'
    do l = Lmin, LevelMax
       print *, get_myrank(), l, count( AllocMem )
    enddo
    stop
  end function alloc_block
  !-------------------------------------------------------------------------
  ! gid が与えられたら block のメモリを確保する。
  !-------------------------------------------------------------------------
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

    bool_omit_coord = .FALSE.   ! default
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
  !-------------------------------------------------------------------------
  ! Free and initialize grid id
  !-------------------------------------------------------------------------
  subroutine dealloc_U(gid)
    integer,intent(IN) :: gid
    integer :: pgid, prank
!!$    if (gid == Undefi) print *,'*** error U'
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
  !-------------------------------------------------------------------------
  ! Free U1 and U2 order
  !-------------------------------------------------------------------------
  subroutine dealloc_U1order(gid)
    integer,intent(IN) :: gid
    call dealloc_block( U1list(gid) )
    U1list(gid) = Undefi
  end subroutine dealloc_U1order
  !-------------------------------------------------------------------------
  subroutine dealloc_U2order(gid)
    integer,intent(IN) :: gid
    call dealloc_block( U2list(gid) )
    U2list(gid) = Undefi
  end subroutine dealloc_U2order
  !-------------------------------------------------------------------------
  ! Free block
  !-------------------------------------------------------------------------
  subroutine dealloc_block(gid)
    integer,intent(IN) :: gid
#define DEALLOC_SAFE(A) \
    if (associated(A)) deallocate(A) ;\
    nullify(A)
    if (gid == Undefi) return
    DEALLOC_SAFE(BlockMem(gid)%u)
    DEALLOC_SAFE(BlockMem(gid)%x)
    DEALLOC_SAFE(BlockMem(gid)%y)
    DEALLOC_SAFE(BlockMem(gid)%z)
    AllocMem( gid ) = .FALSE.
  end subroutine dealloc_block
  !-------------------------------------------------------------------------
  ! GidList を更新する
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! GidListNode を更新する
  !-------------------------------------------------------------------------
  subroutine update_gidlist_node(level)
    integer,intent(IN) :: level
    integer :: rank, last
    integer,dimension(Gidmin:Gidmax) :: glist
    myrank = get_myrank()
    GidListNode(:,level,:) = Undefi
    do rank = 0, NPE-1
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
  !-------------------------------------------------------------------------
  ! 合計のグリッド数
  !-------------------------------------------------------------------------
  function grid_used_num() result(num)
    integer :: num, id
    num = 0
    do id = lbound(Levels,1), ubound(Levels,1)
       if ( Levels(id) /=  Undefi ) then
          num = num + 1
       endif
    enddo
  end function grid_used_num
  !-------------------------------------------------------------------------
  ! 合計のグリッド数 (node global)
  !-------------------------------------------------------------------------
  function grid_used_num_global() result(num)
    integer :: num, num_node
    num_node = grid_used_num()
    call mpi_allreduce(num_node, num, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  end function grid_used_num_global
  !-------------------------------------------------------------------------
  ! level の合計のグリッド数 (node local)
  !-------------------------------------------------------------------------
  function grid_used_num_level( level ) result(num)
    integer,intent(IN) :: level
    integer :: num
    num = GidListMax(level) - Gidmin + 1
  end function grid_used_num_level
  !-------------------------------------------------------------------------
  ! level の合計のグリッド数 (node global)
  !-------------------------------------------------------------------------
  function grid_used_num_level_global( level ) result(num)
    integer,intent(IN) :: level
    integer :: num, num_node
    num_node = grid_used_num_level( level )
    call mpi_allreduce(num_node, num, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  end function grid_used_num_level_global
  !-------------------------------------------------------------------------
  ! dv
  !-------------------------------------------------------------------------
  function get_dv(level) result(dv)
    integer,intent(IN) :: level
    real(kind=DBL_KIND) :: dv
    dv = CellWidth(0,level) * CellWidth(1,level) * CellWidth(2,level)
  end function get_dv
  !-------------------------------------------------------------------------
  ! ds
  !-------------------------------------------------------------------------
  function get_ds(level) result(ds)
    integer,intent(IN) :: level
    real(kind=DBL_KIND) :: ds(MX:MZ)
    ds = (/ &
         CellWidth(1,level) * CellWidth(2,level) , &
         CellWidth(2,level) * CellWidth(0,level) , &
         CellWidth(0,level) * CellWidth(1,level) /)
!!$    ds(MZ) =0                   ! for uniform in z-direction
  end function get_ds
  !-------------------------------------------------------------------------
  ! 自分が親にたいして、右か左かを求める。
  !-------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------
  ! encode to surface code
  !-----------------------------------------------------------------------
  function encode_surf(lr, direc) result( code )
    integer,intent(IN) :: lr, direc
    integer :: code
    code = lr + direc * (MZ-MX)
  end function encode_surf
  ! 方向コードから direction を deocde
  function decode_direction(ncord) result (direc)
    integer,intent(IN) :: ncord
    integer :: direc
    direc = ncord / (MZ-MX)
  end function decode_direction
  ! 方向コードからLRを decode
  function decode_LR(ncord) result (surf)
    integer,intent(IN) :: ncord
    integer :: surf
    surf = modulo(ncord , 2)
  end function decode_LR
  ! swap L and R
  function swap_LR(surf) result(swap)
    integer,intent(IN) :: surf
    integer :: swap
    swap = modulo( surf+1, 2 )
  end function swap_LR
  !-----------------------------------------------------------------------
  ! convert coordinates for cell center
  !-----------------------------------------------------------------------
  subroutine convert_coord(gidfrom, gidto, if,jf,kf, it,jt,kt)
    integer,intent(IN) :: gidfrom, gidto, if,jf,kf
    integer,intent(OUT) :: it,jt,kt
    integer(kind=LLONG_KIND) :: ifoffset,itoffset,jfoffset,jtoffset,kfoffset,ktoffset

    ifoffset = Igrid(gidfrom)*(Imax-Imin+1)
    itoffset = Igrid(gidto)*(Imax-Imin+1)
    jfoffset = Jgrid(gidfrom)*(Jmax-Jmin+1)
    jtoffset = Jgrid(gidto)*(Jmax-Jmin+1)
    kfoffset = Kgrid(gidfrom)*(Kmax-Kmin+1)
    ktoffset = Kgrid(gidto)*(Kmax-Kmin+1)

    if ( Levels(gidfrom) == Levels(gidto) ) then ! within the same grid level
       it = if + ifoffset - itoffset
       jt = jf + jfoffset - jtoffset
       kt = kf + kfoffset - ktoffset
    elseif ( Levels(gidfrom) == Levels(gidto) + 1) then ! coarser
       it = IJKC( if + ifoffset ,0) - itoffset
       jt = IJKC( jf + jfoffset ,0) - jtoffset
       kt = IJKC( kf + kfoffset ,0) - ktoffset
    elseif ( Levels(gidfrom) == Levels(gidto) - 1) then ! finer
       it = IJKF( if + ifoffset ,0) - itoffset
       jt = IJKF( jf + jfoffset ,0) - jtoffset
       kt = IJKF( kf + kfoffset ,0) - ktoffset
    else
       print *, 'It is not supported yet'
       stop
    endif
  end subroutine convert_coord
  !-----------------------------------------------------------------------
  ! convert coordinates for cell boundary
  !-----------------------------------------------------------------------
  subroutine convert_coord_cb(gidfrom, gidto, surf, if,jf,kf, it,jt,kt)
    integer,intent(IN) :: gidfrom, gidto, if,jf,kf, surf
    integer,intent(OUT) :: it,jt,kt
    integer(kind=LLONG_KIND) :: ifoffset,itoffset,jfoffset,jtoffset,kfoffset,ktoffset

    call convert_coord(gidfrom, gidto, if,jf,kf, it,jt,kt)

    if ( Levels(gidfrom) == Levels(gidto) - 1) then ! finer
       if ( surf == MX ) it = it + 1
       if ( surf == MY ) jt = jt + 1
       if ( surf == MZ ) kt = kt + 1
    endif
  end subroutine convert_coord_cb
  !-----------------------------------------------------------------------
  ! get (gid, rank) given by (ig, jg, kg, level). Base grid からの血筋を追う。
  ! whichlevel (optional) is a level, at which (gid, rank) are defined.
  !-----------------------------------------------------------------------
  subroutine get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank, whichlevel)
    integer,intent(IN) :: ig, jg, kg, level
    integer,intent(IN),optional :: whichlevel
    integer,intent(OUT) :: gid, rank
    integer :: lri, lrj, lrk, npos, gidn, rankn, ig0, jg0, kg0, descend

    ! desend = どこまでさかのぼる？
    if ( present( whichlevel ) ) then
       descend = level - whichlevel
    else
       descend = 0
    endif
    if ( descend < 0 ) print *, '*** error in descend'
    ! get ijk in the base grid
    ig0 = ishft(ig, -level)
    jg0 = ishft(jg, -level)
    kg0 = ishft(kg, -level)
    gid = GidBase(ig0, jg0, kg0)
    rank = RankBase(ig0, jg0, kg0)
    ! 親から子へ血筋をたどる
    do npos = level-1, descend, -1
       if ( gid == Undefi ) exit
       if ( rank == MPI_PROC_NULL ) exit
!!$       print *, 'npos',npos
       lri = ibits(ig,npos,1)
       lrj = ibits(jg,npos,1)
       lrk = ibits(kg,npos,1)
       gidn = ChildGid(lri, lrj, lrk, gid, rank)
       rankn = ChildRank(lri, lrj, lrk, gid, rank)
       gid = gidn
       rank = rankn
    enddo
  end subroutine get_gid_from_ijkgrid
  !-------------------------------------------------------------------------
  ! Does a given grid have a child grid?
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! update NeighborGid and NeighborRank
  ! 親子関係(親のChildGid, 子のParentGid)と親の NeighborGid と NeighborRank は既知でなければならない
  !-------------------------------------------------------------------------
  subroutine update_neighbor(gid, thisrank, ijkgrid, level)
    use mpilib
    integer,intent(IN) :: gid, thisrank, ijkgrid(MX:MZ), level
    integer :: pgid, prank, pranklr, pgidlr, m
    integer,dimension(MX:MZ) :: lra, lrb

    do m = MX, MZ
       lra(m) = modulo(ijkgrid(m),2)
       lrb(m) = modulo(ijkgrid(m)+1,2)
    enddo

    pgid = ParentGid(gid, thisrank)
    prank = ParentRank(gid, thisrank)
    ! ChildGid(lra(0), lra(1), lra(2), pgid, prank) ! 自分自身
    ! 自分が R のとき、L 隣を知るには、親の L子
    NeighborGid( lrb(MX), MX, gid, thisrank ) = ChildGid( lrb(MX), lra(MY), lra(MZ), pgid, prank )
    NeighborGid( lrb(MY), MY, gid, thisrank ) = ChildGid( lra(MX), lrb(MY), lra(MZ), pgid, prank )
    NeighborGid( lrb(MZ), MZ, gid, thisrank ) = ChildGid( lra(MX), lra(MY), lrb(MZ), pgid, prank )

    NeighborRank( lrb(MX), MX, gid, thisrank ) = ChildRank( lrb(MX), lra(MY), lra(MZ), pgid, prank )
    NeighborRank( lrb(MY), MY, gid, thisrank ) = ChildRank( lra(MX), lrb(MY), lra(MZ), pgid, prank )
    NeighborRank( lrb(MZ), MZ, gid, thisrank ) = ChildRank( lra(MX), lra(MY), lrb(MZ), pgid, prank )

    ! 自分が R のとき、R 隣を知るには、親の R 隣の L子
    pgidlr = NeighborGid( lra(MX), MX, pgid, prank )
    pranklr = NeighborRank( lra(MX), MX, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(MX), MX, gid, thisrank ) = ChildGid( lrb(MX), lra(MY), lra(MZ), pgidlr, pranklr )
       NeighborRank( lra(MX), MX, gid, thisrank ) = ChildRank( lrb(MX), lra(MY), lra(MZ), pgidlr, pranklr )
    else
       NeighborGid( lra(MX), MX, gid, thisrank ) = Undefi
       NeighborRank( lra(MX), MX, gid, thisrank ) = MPI_PROC_NULL
    endif

    pgidlr = NeighborGid( lra(MY), MY, pgid, prank )
    pranklr = NeighborRank( lra(MY), MY, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(MY), MY, gid, thisrank ) = ChildGid( lra(MX), lrb(MY), lra(MZ), pgidlr, pranklr )
       NeighborRank( lra(MY), MY, gid, thisrank ) = ChildRank( lra(MX), lrb(MY), lra(MZ), pgidlr, pranklr )
    else
       NeighborGid( lra(MY), MY, gid, thisrank ) = Undefi
       NeighborRank( lra(MY), MY, gid, thisrank ) = MPI_PROC_NULL
    endif

    pgidlr = NeighborGid( lra(MZ), MZ, pgid, prank )
    pranklr = NeighborRank( lra(MZ), MZ, pgid, prank )
    if ( pgidlr /= Undefi ) then
       NeighborGid( lra(MZ), MZ, gid, thisrank ) = ChildGid( lra(MX), lra(MY), lrb(MZ), pgidlr, pranklr )
       NeighborRank( lra(MZ), MZ, gid, thisrank ) = ChildRank( lra(MX), lra(MY), lrb(MZ), pgidlr, pranklr )
    else
       NeighborGid( lra(MZ), MZ, gid, thisrank ) = Undefi
       NeighborRank( lra(MZ), MZ, gid, thisrank ) = MPI_PROC_NULL
    endif
  end subroutine update_neighbor
  !-------------------------------------------------------------------------
  ! あるレベルの全てのグリッドの neighbor list を update する。
  !-------------------------------------------------------------------------
  subroutine update_neighbor_linklist(level)
    integer,intent(IN) :: level
    integer :: gid, n, rank, listmin, listmax
    integer,dimension(:),allocatable :: gidl
    integer,dimension(:,:),allocatable :: ijkgrid
    myrank = get_myrank()
    do rank = 0, NPE-1
       ! rank から GidList の要素数 GidListMax( level ) を全 node に転送する。
       ! 結果は listmax
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)

       ! rank から GidList と IJKgrid を 全 node に転送する。
       ! 結果は gidl と ijkgrid
       allocate( gidl( listmin:listmax ) )
       allocate( ijkgrid( MX:MZ, listmin:listmax ) )
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
  !-----------------------------------------------------------------------
  ! 親子関係の link list を作り直す
  ! 親までのリンクリストが正しい必要がある。
  !-----------------------------------------------------------------------
  subroutine update_parent_linklist(level)
    integer,intent(IN) :: level
    integer :: n, cgid, pgid, crank, prank, rank, listmin, listmax, lri, lrj, lrk, ig,jg,kg
    integer,dimension(:),allocatable :: gidl
    integer,dimension(:,:),allocatable :: ijkgrid

    myrank = get_myrank()

    do rank = 0, NPE -1
       ! rank から GidList の要素数 GidListMax( level ) を全 node に転送する。
       ! 結果は listmax
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)

       ! rank から GidList と IJKgrid を 全 node に転送する。
       ! 結果は gidl と ijkgrid
       allocate( gidl( listmin:listmax ) )
       allocate( ijkgrid( MX:MZ, listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level)
          do n = listmin, listmax
             ijkgrid(:,n) = (/ Igrid(gidl(n)), Jgrid(gidl(n)), Kgrid(gidl(n)) /)
          enddo
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       call mpi_bcast( ijkgrid, size(ijkgrid), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       ! 順に, 親子の link を接続する
       do n = listmin, listmax
          cgid = gidl(n)
          crank = rank
          ig = ijkgrid(MX,n)
          jg = ijkgrid(MY,n)
          kg = ijkgrid(MZ,n)
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
  !-----------------------------------------------------------------------
  ! リンクリストをクリアする。
  ! 対象： (1) 自分のNeighbor*, (2) 自分から見た Parent* Child*,
  !        (3) 親から見た Child*, (4) 子から見た Parent*
  !        自分 = level, 子 = level + 1, 親 = level - 1
  !-----------------------------------------------------------------------
  subroutine clear_linklist(level)
    integer,intent(IN) :: level
    integer :: n, gid, cgid, pgid, crank, prank, rank, listmin, listmax, lri, lrj, lrk
    integer,dimension(:),allocatable :: gidl

    if ( level <= 0 ) return

    myrank = get_myrank()

    ! 自分のレベルからのリンク
    do rank = 0, NPE -1
       ! rank から GidList の要素数 GidListMax( level ) を全 node に転送する。
       ! 結果は listmax
       if ( rank == myrank ) &
            listmax = GidListMax( level )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)

       ! rank から GidList と IJKgrid を 全 node に転送する。
       ! 結果は gidl と ijkgrid
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       ! 順に, 親子の link をクリアする
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

    ! 親のレベルからのリンク
    do rank = 0, NPE -1
       ! rank から GidList の要素数 GidListMax( level-1 ) を全 node に転送する。
       ! 結果は listmax
       if ( rank == myrank ) &
            listmax = GidListMax( level-1 )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)

       ! rank から GidList と IJKgrid を 全 node に転送する。
       ! 結果は gidl と ijkgrid
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level-1)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       ! 順に, 親子の link をクリアする
       do n = listmin, listmax
          pgid = gidl(n)
          prank = rank
          ChildGid(:,:,:,pgid,prank) = Undefi
          ChildRank(:,:,:,pgid,prank) = MPI_PROC_NULL
       enddo
       deallocate( gidl )
    enddo

    ! 子のレベルからのリンク
    if ( level + 1 > LevelMax ) return
    do rank = 0, NPE -1
       ! rank から GidList の要素数 GidListMax( level-1 ) を全 node に転送する。
       ! 結果は listmax
       if ( rank == myrank ) &
            listmax = GidListMax( level+1 )
       call mpi_bcast(listmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       listmin = lbound(GidList,1)

       ! rank から GidList と IJKgrid を 全 node に転送する。
       ! 結果は gidl と ijkgrid
       allocate( gidl( listmin:listmax ) )
       if ( rank == myrank ) then
          gidl = GidList(listmin:listmax, level+1)
       endif
       call mpi_bcast( gidl, size(gidl), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       ! 順に, 親子の link をクリアする
       do n = listmin, listmax
          cgid = gidl(n)
          crank = rank
          ParentGid(cgid,crank) = Undefi
          ParentRank(cgid,crank) = MPI_PROC_NULL
       enddo
       deallocate( gidl )
    enddo

  end subroutine clear_linklist
  !-----------------------------------------------------------------------
  ! 親子link list をチェックする。
  !-----------------------------------------------------------------------
  subroutine check_linklist(level)
    integer,intent(IN) :: level
    integer :: n, gid, pgid, prank, i,j,k, rank
    integer,dimension(Gidmin:Gidmax,0:NPE-1) :: ptmp
    integer,dimension(0:1,0:1,0:1,Gidmin:Gidmax,0:NPE-1) :: ctmp
    logical :: bool, boolt, boolr

    myrank = get_myrank()
    boolt = .TRUE.

    ! check consisitency among nodes
    bool = .TRUE.
    do rank = 0, NPE-1
       ptmp = ParentGid
       call mpi_bcast(ptmp, size(ptmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ptmp-ParentGid)-minval(ptmp-ParentGid) /= 0) then
          print *, '** bad  consistency ParentGid, node=',rank
          bool = .false.
       end if
    enddo
    do rank = 0, NPE-1
       ptmp = ParentRank
       call mpi_bcast(ptmp, size(ptmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ptmp-ParentRank)-minval(ptmp-ParentRank) /= 0) then
          print *, '** bad  consistency RarentRank, node=',rank
          bool = .false.
       end if
    enddo
    do rank = 0, NPE-1
       ctmp = ChildGid
       call mpi_bcast(ctmp, size(ctmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ctmp-ChildGid)-minval(ctmp-ChildGid) /= 0) then
          print *, '** bad  consistency ChidlGid node=', rank
          bool = .false.
       end if
    enddo
    do rank = 0, NPE-1
       ctmp = ChildRank
       call mpi_bcast(ctmp, size(ctmp), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
       if (maxval(ctmp-ChildRank)-minval(ctmp-ChildRank) /= 0) then
          print *, '** bad  consistency ChildRank node =', rank
          bool = .false.
       end if
    enddo
    boolt = boolt .and. bool

    ! check link from child to parent
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
    ! check parent to child
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

    ! check coordinates
!!$    do n = Gidmin, GidListMax( level )
!!$       gid = GidList(n, level)
!!$       print *, Igrid(gid), Jgrid(gid), Kgrid(gid)
!!$    enddo

  end subroutine check_linklist

  ! -----------------------------------------------------------------
  ! どのグリッドまで同期?
  ! -----------------------------------------------------------------
  function level_sync() result(levsync)
    integer :: l, levsync, id
    integer(kind=LLONG_KIND) :: steplmax
    id = 0
    steplmax = Step( LevelMax )
    levsync = LevelMax
    do l = LevelMax-1, Lmin, -1
       if ( Step( l ) == steplmax ) levsync = l
    enddo
  end function level_sync
  ! -----------------------------------------------------------------
  ! get the level, in which time step proceeds
  ! -----------------------------------------------------------------
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
