
#include "config.h"
#include "packarr.h"
! ななめも refine するか？
#define NANAME
!!$#define WAVE_TEST
!-------------------------------------------------------------------------
! Module for refinement of grid
! 孫との隣接を禁止する。
!-------------------------------------------------------------------------
module refine
  use grid
  use mpilib
  implicit none
  private
  ! refinement bool Map
  ! 細分化（子グリッドの存在）が必要な親グリッドに .TRUE. のマークを付ける。
  logical,save,dimension(Gidmin:Gidmax,0:NPE-1) :: BoolMap
  ! list of grid id that should be refined. 虫食い左詰めリスト。 node global
  ! Finer   = サブグリッドを生成するべき grid id
  ! Fine    = サブグリッドを維持するべき grid id
  ! Coarser = サブグリッドを消すべき grid id
  ! Coarse  = サブグリッドがない grid id
  ! Newrank = refine するときの rank (新しいサブグリッドのrank)
  ! Ihash   = 逆引用index. Ihash(gid, rank) = リストの要素番号
  integer,save,dimension(Gidmin:Gidmax,0:NPE-1) :: Finer, Coarser, Fine, Coarse, Ihash
  integer,save,dimension(Left:Right,Left:Right,Left:Right,Gidmin:Gidmax,0:NPE-1) :: Newrank
  ! Finer の gid の要素で、転送に関する左詰めリスト. node global
  integer,save :: Levelf, Levelc     ! levels of subgrid and parent grid
  logical,save :: BoolRefine         ! refinement が必要か？
  integer,parameter :: Rgh = 1 ! numner of ghost cell for interpolation of refinement
  public :: refineLevel, refineAllLevel, refine_get_norder
contains

  !-------------------------------------------------------------------------
  ! refine all level
  !-------------------------------------------------------------------------
  subroutine refineALlLevel
    logical :: bool
    integer :: level
    do level = Lmin+1, Lmax
       call refineLevel(level, bool)
       if (.not. bool) exit
    enddo
  end subroutine refineALlLevel
  !-------------------------------------------------------------------------
  ! refine given level
  ! グリッドレベル level を作る。
  ! このレベルが必要なければ、bool = .false. を返す。
  ! それ以外は bool = .true. を返す。
  !-------------------------------------------------------------------------
  subroutine refineLevel(level, bool)
    use io_util, only : print_msg
    use string, only : num2char
    use grid, only : alloc_U1order, alloc_U2order, dealloc_U1order, dealloc_U2order
    use grid_boundary
    use boundary
#ifdef SINKPARTICLE
    use sinkParticle
#endif ! SINKPARTICLE
    use modelParameter, only : MP_CONNECTION_RUN  !KS ADDED

    integer,intent(IN) :: level
    logical,intent(OUT) :: bool ! refined?
    integer :: n

    !----------- KS DEBUG -----------!
    globdbg_gridupdate = .False.
    !----------- KS DEBUG -----------!

    Levelf = level
    Levelc = level - 1
    BoolRefine = .FALSE.
    bool = .FALSE.
    if ( Levelc < Lmin ) return
    if ( Levelf > Lmax ) return

    !----------------------------------------------------------
    ! level > sp_level may exist during transition period
    !----------------------------------------------------------
    if (MP_CONNECTION_RUN == 0) then !KS ADDED
#ifdef SINKPARTICLE
       if ( Levelf > sp_getLevel() ) return
#endif ! SINKPARTICLE
    end if ! KS ADDED

!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call print_msg( 'refined level = '//num2char(level) )
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)

    ! refinement の前に親の境界条件を整えておく
    call boundary_grid( Levelc , COMPLETE)
    do n = Gidmin, GidListMax( Levelc )
       call boundary_u( GidList(n, Levelc), COMPLETE )
    enddo

    call makeBoolRefine
    bool = BoolRefine

    !----------- KS DEBUG -----------!
    globdbg_gridupdate = bool_update_grid()
    !----------- KS DEBUG -----------!
    

    ! 実際に、グリッドを変更するか？
    if ( .not. bool_update_grid() ) return

    ! グリッドに付随するリストを壊す
    if ( Levelf <= LevelMax ) then
       do n = Gidmin, GidListMax( Levelf )
          call dealloc_U1order( GidList(n, Levelf) )
          call dealloc_U2order( GidList(n, Levelf) )
       enddo
    endif

    call do_coarser

    if ( bool ) then
       LevelMax = max(LevelMax,Levelf)
       call define_Newrank
       call do_fine
       call do_finer
    else
!!$       LevelMax = LevelMax -1
!!$       LevelMax = Levelc
       LevelMax = max(LevelMax-1, levelc)
    endif

    call update_gidlist(Levelf)
    call clear_linklist(Levelf)
    call update_parent_linklist(Levelf)
    if (Levelf + 1 <= LevelMax) then
       call update_parent_linklist(Levelf+1)
    endif
    call update_neighbor_linklist(Levelf)

    ! グリッドに付随するリストを作る
    if ( Levelf <= LevelMax ) then
       do n = Gidmin, GidListMax( Levelf )
          call alloc_U1order( GidList(n, Levelf) )
          call alloc_U2order( GidList(n, Levelf) )
       enddo
    endif

    ! check

!!$    myrank = get_myrank()
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    print *, 'Finer', count( Finer /= Undefi ), myrank
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    print *, 'Newrank', count( Newrank /= MPI_PROC_NULL ), myrank
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    print *, 'BufGid', count( BufGid /= Undefi ), LastB, myrank
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    print *, 'BoolTmp', count( BoolTmp ), myrank
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    print *, 'Levels',  count(Levels == Levelf), myrank
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)


!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$    call check_linklist(Levelf)
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr)

!!$       call mpi_finalize(ierr)
!!$       stop
!!$    call test
  end subroutine refineLevel
  !-------------------------------------------------------------------------
  ! Finer, Coarser, Fine, Coarse のリストを作る。(重い？)
  !-------------------------------------------------------------------------
  subroutine makeBoolRefine()
    use eos
    use refineCond
#ifdef SINKPARTICLE
    use sinkParticle
#endif ! SINKPARTICLE
    use unit ! KS DEBUG
    use modelParameter, only : MP_CONNECTION_RUN  !KS ADDED
    
    integer :: gid, n,m,lr, gidn, rankn, rank
    logical :: bufbool, boolself, boolgc
    logical,dimension(Gidmin:Gidmax,0:NPE-1) :: bufboolmap
    integer :: lr1, lr2, lr3, m1, m2, m3, &
         gid1, gid2, gid3, gid4, gid5, gid6, gid7, &
         rank1, rank2, rank3, rank4, rank5, rank6, rank7, &
         i, j, k, lbuf, ic, jc, kc, gidc, rankc

    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z    ! KS DEBUG
    
    
    myrank = get_myrank()
    ! ------------------
    ! BoolMap を定義する
    ! ------------------
    BoolMap(:,:) = .FALSE.

    do n = lbound(GidList,1), GidListMax( Levelc )
       boolself = .false.
       boolgc   = .false.
       gid = GidList(n, Levelc)

#ifdef Emulate_1Dim
       if (.not. Jgrid(gid)  == 0 ) cycle
#endif !Emulate_1Dim
#ifdef EMULATE_2DIM
       if (.not. Kgrid(gid)  == 0 ) cycle
#endif !EMULATE_2DIM

       ! 孫を調べる。ひとりでもいれば true
       do k = Left, Right
          do j = Left, Right
             do i = Left, Right
                gidc = ChildGid(i,j,k,gid,myrank)
                rankc = ChildRank(i,j,k,gid,myrank)
                if ( gidc == Undefi ) cycle
                if ( any( ChildGid(:,:,:, gidc, rankc) /= Undefi) )  &
                     boolgc = .true.
             enddo
          enddo
       enddo
       
       ! 自分を調べる。
       boolself = refineCond_eval(gid)
       ! ! KS DEBUG
       ! if (myrank == 232 .and. gid == 138) then
       !    print *, '(KS DEBUG) gc, jns: ', boolgc, boolself
       ! end if
       
#ifdef SINKPARTICLE
       ! call sp_refineCond(gid, boolself)
       call sp_refineCond_KS(gid, boolself) ! KS MODIFIED
#endif ! SINKPARTICLE
       ! ! KS DEBUG
       ! if (myrank == 232 .and. gid == 138) then
       !    print *, '(KS DEBUG) sp:', boolself
       ! end if

       ! このブロックは無関係なので次へ
       if (.not. (boolgc .or. boolself) ) cycle
       BoolMap(gid, myrank) = boolgc .or. boolself

       ! check  粗から細へ細分化で必要
!!$       if (.not. checkNeighbor(gid, myrank)) then
!!$          print *, '****ERROR in refinment'
!!$          cycle
!!$       endif

       ! 孫はマージンのために存在するのか？
#ifdef SINGLE_STEP
       !------ WARNING WARNING (KS MODIFIED) WARNING WARNING ------!
       ! if ( .not. ( boolself .or. have_allgrandchild(gid, myrank) ))  cycle
       !兄弟の孫のために子を持つ？ おじさんの存在を保証 (KS MODIFIED)
       if ( .not. ( boolself .or. have_allgrandchild(gid, myrank) ))  then
       ! 孫がいるサイドの兄弟に子を持つよう要請（マージンとは少し異なる要請）
          do k = Left, Right
             do j = Left, Right
                do i = Left, Right
                   gidc = ChildGid(i,j,k,gid,myrank)
                   rankc = ChildRank(i,j,k,gid,myrank)
                   if ( gidc == Undefi ) cycle
                   !孫がいるとき
                   if ( any( ChildGid(:,:,:, gidc, rankc) /= Undefi) )  then
                      !孫がいるサイドの兄弟をTrueにする
                      !x方向
                      gid1 = NeighborGid(i,MX,gid,myrank)
                      rank1 = NeighborRank(i,MX,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.
                      !y方向
                      gid1 = NeighborGid(j,MY,gid,myrank)
                      rank1 = NeighborRank(j,MY,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.                   
                      !z方向
                      gid1 = NeighborGid(k,MZ,gid,myrank)
                      rank1 = NeighborRank(k,MZ,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.                   
                   end if
                enddo
             enddo
          enddo
          cycle !これ以上のマージンは求めない
       end if
       !------ WARNING WARNING (KS MODIFIED) WARNING WARNING ------!       
#endif

#ifdef WAVE_TEST
       ! 最細グリッドはマージンなし。(衝撃波等の問題に有効?)
       if ( Lmax == Levelf ) cycle
#endif

#ifdef REFINE_NO_MARGIN
       cycle                    ! without margin. only for fixed grid
#endif

       do m1 = MX, MZ
          do lr1 = Left , Right     ! 隣を True にする
             gid1 = NeighborGid(lr1,m1,gid,myrank)
             rank1 = NeighborRank(lr1,m1,gid,myrank)
             if (gid1 == Undefi .or. rank1 == MPI_PROC_NULL ) cycle
!!$             ! 四方の隣接は兄弟でなければならない
!!$             if ( .not. checkNeighbor(gid1,rank1) ) cycle   ! 粗から細へ細分化で必要
             ! 登録する
             BoolMap( gid1, rank1 ) = .TRUE.
#ifdef NANAME
             do m2 = MX, MZ  ! 平面上の斜めも True にする
                if ( m2 == m1 ) cycle
                do lr2 = Left, Right
                   gid2 = NeighborGid(lr2,m2,gid1,rank1)
                   rank2 = NeighborRank(lr2,m2,gid1,rank1)
                   if (gid2 == Undefi .or. rank2 == MPI_PROC_NULL ) cycle
                   BoolMap( gid2, rank2 ) = .TRUE.
                   do m3 = MX, MZ  ! 立体の斜めも True にする
                      if ( m3 == m2 .or. m3 == m1 ) cycle
                      do lr3 = Left, Right
                         gid3 = NeighborGid(lr3,m3,gid2,rank2)
                         rank3 = NeighborRank(lr3,m3,gid2,rank2)
                         if (gid3 == Undefi .or. rank3 == MPI_PROC_NULL ) cycle
                         BoolMap( gid3, rank3 ) = .TRUE.
                      enddo
                   enddo
                enddo
             enddo
#endif !NANAME
          enddo
       enddo
    enddo


    call mpi_allreduce_BoolMap(MPI_LOR)

    ! 不良ブロックを除外する
    ! 2段飛びの不良グリッドはゴーストセルの値を正しく決められないので除外する必要があるらしい by 松本さん

    do n = lbound(GidList,1), GidListMax( Levelc )
       gid = GidList(n, Levelc)
       if (.not. checkNeighbor(gid, myrank)) then
          BoolMap(gid, myrank) = .FALSE.
       !-------------- KS MODIFIED (needs more check) --------------!
       ! if (.not. checkNeighbor_NANAME(gid, myrank)) then
       !    BoolMap(gid, myrank) = .FALSE.
       !-------------- KS MODIFIED (needs more check) --------------!          
          
!!$          print *, '*** Warning in refinement: bad block is detected.', gid, myrank,Levelc
          !------------------- KS DEBUG ------------------!
          ! x => get_Xp(gid)
          ! y => get_Yp(gid)
          ! z => get_Zp(gid)          
          ! print '(A,3I5, 1P6E15.7)', '(KS DEBUG) bad block in refinement: ', gid, myrank,Levelc,& !KS DEBUG
          !      x(Imin)*Unit_au, x(Imax)*Unit_au, y(Jmin)*Unit_au, y(Jmax)*Unit_au, z(Kmin)*Unit_au, z(Kmax)*Unit_au
          !------------------- KS DEBUG ------------------!          
       endif
    end do

    call mpi_allreduce_BoolMap(MPI_LAND) ! 今度はAND
    

#ifdef EMULATE_2DIM
    call Emulate_2Dim_BoolMap   ! Kgrid(gid) == 0 のグリッドを BoolMap = .false. とする
#endif !EMULATE_2DIM
    BoolRefine = any( BoolMap )
    call mpi_allreduce( BoolRefine, bufbool, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    BoolRefine = bufbool
    ! ------------------------
    ! Refiment List を定義する
    ! ------------------------
    Finer(:,:) = Undefi
    Coarser(:,:) = Undefi
    Fine(:,:) = Undefi
    Coarse(:,:) = Undefi
    Ihash(:,:) = Undefi
    Newrank(:,:,:,:,:) = MPI_PROC_NULL
    do rank = 0, NPE-1
       ! Lastelem(rank) -> GidListNodeMax(Levelc, rank)
       ! glist(n)       -> GidListNode(n, Levelc, rank)

       ! Finer, Fine, Coarse, or Coarser?
       ! node global process
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          gid = GidListNode(n, Levelc, rank)
          Ihash(gid,rank) = n   ! 逆引用インデックス
          if ( BoolMap(gid,rank) ) then
             if ( ChildGid(0,0,0,gid,rank) == Undefi ) then
                Finer(n,rank) = gid     !子がいないので、新たに refine
             else
                Fine(n,rank) = gid      !すでに子がいるので、現状維持
             endif
          else
             if ( ChildGid(0,0,0,gid,rank) == Undefi ) then
                Coarse(n,rank) = gid     !すでに子がいないので、現状維持
             else
                Coarser(n,rank) = gid    !子を殺す
             endif
          endif
       enddo
    enddo
    ! !KS DEBUG
    ! if (Levelc == 12 .and. myrank == 232) then
    !    print *, '(KS DEBUG) BM3, cgUndefi:', BoolMap( 138, 232 ), ChildGid(0,0,0,gid,rank)==Undefi
    ! end if

    !------------------- KS DEBUG --------------------------!
    if ( get_myrank() == PRIMARY_RANK ) then
       !if (mod(Step(Lmin),100) == 0) then
       if (mod(Step(Lmin),100) == 0 .or. MP_CONNECTION_RUN > 0) then          
          print '(A,4I6)', 'refine grids: Coarser, Coarse, Fine, Finer =>  ', &
               count( Coarser /= Undefi ), count( Coarse /= Undefi ), &
               count( Fine /= Undefi ), count( Finer /= Undefi )
       endif
    endif
    !------------------- KS DEBUG --------------------------!
  contains
    !-------------------------------------------------------------------------
    ! 四方全てに兄弟がいるかどうか、チェックする。
    !-------------------------------------------------------------------------
    function checkNeighbor(gid,rank) result(bool)
      integer,intent(IN) :: gid, rank
      logical :: bool
      integer :: lr, ndir, pgid, prank
      ! 境界に接するかどうかの判定：Igrid, Jgrid, Kgrid から決めても良いが、
      ! ここでは、親の隣が存在するかどうかで判定する。
      bool = .true.
      if (Levelc == Lmin) return
      pgid = ParentGid(gid,rank)
      prank = ParentRank(gid,rank)
      do ndir = MX, MZ
         do lr = Left, Right
            if ( NeighborGid(lr, ndir, pgid, prank) /= Undefi .and. & ! 境界に接しない
                 NeighborGid(lr, ndir, gid, rank) == Undefi ) &    ! 隣が未定義
                 bool = .false.
         enddo
      enddo
    end function checkNeighbor

    !-------------------------------------------------------------------------
    ! 斜めも含めて周囲全てに兄弟がいるかどうか、チェックする。 (KS ADDED)
    ! 境界に接しているかのチェック方法が正しいか自信なし
    ! *** あまり意味無さそうだったのでこの関数は使わない ***
    !-------------------------------------------------------------------------
    function checkNeighbor_NANAME(gid,rank) result(bool)
      integer,intent(IN) :: gid, rank
      logical :: bool
      integer :: pgid, prank, m1,m2,m3,lr1,lr2,lr3,gid1,gid2,gid3,rank1,rank2,rank3
      ! 境界に接するかどうかの判定：Igrid, Jgrid, Kgrid から決めても良いが、
      ! ここでは、親の隣が存在するかどうかで判定する。
      bool = .true.
      if (Levelc == Lmin) return
      pgid = ParentGid(gid,rank)
      prank = ParentRank(gid,rank)

      do m1 = MX, MZ   ! 隣を調べる
         do lr1 = Left , Right 
            ! if ( NeighborGid(lr1, m1, pgid, prank) == Undefi) cycle ! 境界に接していたらスキップ
            gid1 = NeighborGid(lr1,m1,gid,rank)
            rank1 = NeighborRank(lr1,m1,gid,rank)
            if (gid1 == Undefi )  bool = .false.    ! 未定義なら .false.

            do m2 = MX, MZ  ! 平面上の斜めを調べる
               if ( m2 == m1 ) cycle
               do lr2 = Left, Right
                  ! if ( NeighborGid(lr2, m2, pgid, prank) == Undefi) cycle ! 境界に接していたらスキップ
                  gid2 = NeighborGid(lr2,m2,gid1,rank1)
                  rank2 = NeighborRank(lr2,m2,gid1,rank1)
                  if (gid2 == Undefi )  bool = .false.    ! 未定義なら .false.

                   do m3 = MX, MZ  ! 立体の斜めを調べる
                      if ( m3 == m2 .or. m3 == m1 ) cycle
                      do lr3 = Left, Right
                         ! if ( NeighborGid(lr3, m3, pgid, prank) == Undefi) cycle ! 境界に接していたらスキップ
                         gid3 = NeighborGid(lr3,m3,gid2,rank2)
                         rank3 = NeighborRank(lr3,m3,gid2,rank2)
                         if (gid3 == Undefi )  bool = .false.    ! 未定義なら .false.

                      end do
                   end do
                end do
             end do
          end do
       end do
      
     end function checkNeighbor_NANAME
    
  end subroutine makeBoolRefine
  !-------------------------------------------------------------------------
  ! 実際にグリッドを更新する必要があるか？Finer, Coarserをもとに判断する。
  !-------------------------------------------------------------------------
  function bool_update_grid() result(bool)
    logical :: bool
    bool = ( any(Finer /= Undefi ) .or. any(Coarser /= Undefi ) )
  end function bool_update_grid
  !-------------------------------------------------------------------------
  ! ブロック(gid, rank)が覆う領域全てに孫が存在するか？
  ! 幾何学により、マージンではない本物の孫ブロックは、当該ブロック全ての領域を覆っている。
  !-------------------------------------------------------------------------
  function have_allgrandchild(gid, rank) result(bool)
    integer,intent(IN) :: gid, rank
    logical :: bool
    integer :: i, j, k, gidc, rankc
    bool = .true.
    do k = Left, Right ! for all child
       do j = Left, Right
          do i = Left, Right
             gidc = ChildGid(i,j,k,gid,rank)
             rankc = ChildRank(i,j,k,gid,rank)
             if ( gidc == Undefi ) then ! 子供すら不在
                bool = .false.
                return
             endif
             if ( any(ChildGid(:,:,:,gidc,rankc) == Undefi) ) then !孫不在
                bool = .false.
                return
             endif
          enddo
       enddo
    enddo
  end function have_allgrandchild
  !-------------------------------------------------------------------------
  ! allreduce for BoolMap
  ! This routine is equivalent to:
  ! mpi_allreduce( MPI_IN_PLACE, BoolMap, size(BoolMap), MPI_LOGICAL, operator_mpi, MPI_COMM_WORLD, ierr)
  !-------------------------------------------------------------------------
  subroutine mpi_allreduce_BoolMap(operator_mpi)
    integer,intent(IN) :: operator_mpi
    logical,dimension(:),allocatable :: buf
    integer :: bufsize, rank, pos, n, gid
    bufsize = sum( max((GidListNodeMax(Levelc, :)-Gidmin+1),0) )
    if (bufsize == 0) return
    allocate(buf(bufsize))
    ! pack BoolMap into buf
    pos = 1
    do rank = 0, NPE-1
       if ( GidListNodeMax(Levelc, rank) == Undefi ) cycle
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          gid = GidListNode(n, Levelc, rank)
          if (gid == Undefi) cycle
          buf(pos) = BoolMap(gid, rank)
          pos = pos + 1
       end do
    end do
    if (pos-1 /= bufsize) print *, '**** error in mpi_allreduce_BoolMap', pos-1, bufsize ! check
    call mpi_allreduce( MPI_IN_PLACE, buf, bufsize, MPI_LOGICAL, operator_mpi, MPI_COMM_WORLD, ierr)
    ! unpack buf into BoolMap
    pos = 1
    do rank = 0, NPE-1
       if ( GidListNodeMax(Levelc, rank) == Undefi ) cycle
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          gid = GidListNode(n, Levelc, rank)
          if (gid == Undefi) cycle
          BoolMap(gid, rank) = buf(pos)
          pos = pos + 1
       end do
    end do
    deallocate(buf)
  end subroutine mpi_allreduce_BoolMap
  !-------------------------------------------------------------------------
  ! 新しい rank リスト (Newrank) を定義する。
  !-------------------------------------------------------------------------
  subroutine define_Newrank
    integer :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer,parameter :: LARGE_LENGTH = 2**10
    ! 計算領域の範囲(子グリッド)
    call grid_range(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    ! define_Newrank を子グリッドの領域の広さに応じてスイッチする。
    ! 実験によると新方式(define_Newrank_largeArea)に固定しても十分高速であった。
    ! したがってこのスイッチは将来は不要かもしれない。
#ifdef SWITCH_DEFINE_NEWRANK_METHOD
    if (max(igmax-igmin+1, jgmax-jgmin+1, kgmax-kgmin+1) <= LARGE_LENGTH) then
       call define_Newrank_smallArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    else
       call define_Newrank_largeArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    endif
#else ! SWITCH_DEFINE_NEWRANK_METHOD
    call define_Newrank_largeArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
#endif ! SWITCH_DEFINE_NEWRANK_METHOD
  end subroutine define_Newrank
  !-------------------------------------------------------------------------
  ! 子グリッドの範囲(igmin, jgmin, kgmin, igmax, jgmax, kgmax)を求める
  !-------------------------------------------------------------------------
  subroutine grid_range(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(OUT) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: rank, n, gid
    integer,dimension(MX:MZ) :: buf, bufr ! buffer
    ! search grid range in the parent level
    rank = get_myrank()
    igmin = Huge(igmin)
    jgmin = Huge(jgmin)
    kgmin = Huge(kgmin)
    igmax = -Huge(igmax)
    jgmax = -Huge(jgmax)
    kgmax = -Huge(kgmax)
    do n = Gidmin, GidListNodeMax(Levelc, rank)
       if ( Fine(n, rank) /= Undefi ) then
          gid = Fine(n, rank)
       elseif ( Finer(n, rank) /= Undefi ) then
          gid = Finer(n, rank)
       else
          cycle
       endif
       igmin = min(igmin, Igrid(gid))
       jgmin = min(jgmin, Jgrid(gid))
       kgmin = min(kgmin, Kgrid(gid))
       igmax = max(igmax, Igrid(gid))
       jgmax = max(jgmax, Jgrid(gid))
       kgmax = max(kgmax, Kgrid(gid))
    enddo
    ! data transfer
    buf = (/igmin, jgmin, kgmin/)
    call mpi_allreduce(buf, bufr, size(buf), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr )
    igmin = bufr(MX);  jgmin = bufr(MY); kgmin = bufr(MZ)
    buf = (/igmax, jgmax, kgmax/)
    call mpi_allreduce(buf, bufr, size(buf), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
    igmax = bufr(MX);  jgmax = bufr(MY); kgmax = bufr(MZ)
!!$      print *, 'gid_range', igmin, jgmin, kgmin, igmax, jgmax, kgmax
    ! convert grid range in the parent level to that in the child level
    igmin = IJKF(igmin, 0)    ! left side
    jgmin = IJKF(jgmin, 0)
    kgmin = IJKF(kgmin, 0)
    igmax = IJKF(igmax, 0) + 1 ! right side
    jgmax = IJKF(jgmax, 0) + 1
    kgmax = IJKF(kgmax, 0) + 1
  end subroutine grid_range
  !-------------------------------------------------------------------------
  ! 各ノードが分担するグリッド数を求める。
  !-------------------------------------------------------------------------
  subroutine get_ngrid_node(ngrid_node)
    integer,intent(OUT) :: ngrid_node(0:NPE-1)
    integer :: ngrid_total, rank, n
    ! 子グリッド数の合計 (ngrid_total)
    ngrid_total = 0
    do rank = 0, NPE -1
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          if ( Fine(n, rank) /= Undefi .or. Finer(n, rank) /= Undefi ) &
               ngrid_total = ngrid_total + 1
       enddo
    enddo
    ngrid_total = ngrid_total * 8  ! 子供のグリッドは8個なので
    ! 各ノードの分担数(ngrid_node(:))
    ngrid_node(:) = int(ngrid_total/(NPE))
    do n = 1, mod( ngrid_total, NPE )
       ngrid_node(NPE-n) = ngrid_node(NPE-n) + 1
    enddo
  end subroutine get_ngrid_node
  !-------------------------------------------------------------------------
  ! get norder
  !-------------------------------------------------------------------------
  function refine_get_norder(igmin,jgmin,kgmin,igmax,jgmax,kgmax) result(norder)
    integer,intent(IN) :: igmin,jgmin,kgmin,igmax,jgmax,kgmax
    integer :: norder, n, width
    width = max(igmax-igmin+1, jgmax-jgmin+1, kgmax-kgmin+1 )
    n = 0
    do
       if ( 2**n >= width ) exit
       n = n + 1
    end do
    norder = n
  end function refine_get_norder
  !-------------------------------------------------------------------------
  ! 小さな領域のオーダリングを行う。(従来方式)
  !-------------------------------------------------------------------------
  subroutine define_Newrank_smallArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(IN) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: ranknew, gcount, norder
    integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
    integer :: ifc, jfc, kfc
    integer :: ngrid_node(0:NPE-1)
    ! 各ノードが分担するグリッド数
    call get_ngrid_node(ngrid_node)
    ! fill space
    ranknew = 0                   ! new rank id
    myrank = get_myrank()         ! myrank
    gcount = 0                    ! counter of grid for each rank
    ifc = igmin ; jfc = jgmin ; kfc = kgmin   ! point of filling curve
    call assign(ifc,jfc,kfc)
    norder = refine_get_norder(igmin,jgmin,kgmin,igmax,jgmax,kgmax)
    call fillcells(norder, R,U,L,D,B,F)
    if (maxval(NewRank) > NPE) then
       print *, 'error in define_Newrank (NewRank_max, NPE) =', maxval(NewRank), NPE
       stop
    endif
!!$    call checkNewrank
    !-------------------------------------------------------------------------
  contains
    !----------------------------------------------------------------------
    ! Peano-Hilbert filling curve
    !----------------------------------------------------------------------
    recursive subroutine fillcells(n, right, up, left, down, back, foward)
      integer,intent(IN) :: n, right, up, left, down, back, foward
      if (n == 0) return
      call fillcells(n-1,foward,right,back,left,down,up)
      call connect(up)
      call fillcells(n-1,up,foward,down,back,left,right)
      call connect(right)
      call fillcells(n-1,up,foward,down,back,left,right)
      call connect(down)
      call fillcells(n-1,left,down,right,up,back,foward)
      call connect(foward)
      call fillcells(n-1,left,down,right,up,back,foward)
      call connect(up)
      call fillcells(n-1,up,back,down,foward,right,left)
      call connect(left)
      call fillcells(n-1,up,back,down,foward,right,left)
      call connect(down)
      call fillcells(n-1,back,right,foward,left,up,down)
    end subroutine fillcells
    !-----------------------------------------------------------------------
    subroutine connect(dir)
      integer,intent(IN) :: dir
      select case (dir)
      case (R)
         ifc = ifc + 1
      case (L)
         ifc = ifc - 1
      case (F)
         jfc = jfc + 1
      case (B)
         jfc = jfc - 1
      case (U)
         kfc = kfc + 1
      case (D)
         kfc = kfc - 1
      end select
      call assign(ifc,jfc,kfc)
    end subroutine connect
    !-----------------------------------------------------------------------
    subroutine assign(ifc,jfc,kfc)
      integer,intent(IN) :: ifc,jfc,kfc
      integer :: gid, rank, i, j, k, ic, jc, kc
      if ( ifc < igmin .or. ifc > igmax .or. &
           jfc < jgmin .or. jfc > jgmax .or. &
           kfc < kgmin .or. kfc > kgmax ) return
      ! grid coordinates in the parent level
      ic = IJKC(ifc, 0)
      jc = IJKC(jfc, 0)
      kc = IJKC(kfc, 0)
      call get_gid_from_ijkgrid(ic,jc,kc,Levelc,gid,rank) ! get gid and rank from coordinates and level
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( Finer(Ihash(gid,rank),rank) == Undefi .and. &
           Fine(Ihash(gid,rank),rank) == Undefi ) return

      gcount = gcount + 1
      if ( gcount > ngrid_node(ranknew) ) then ! 次のランク
         ranknew = ranknew + 1
         gcount = 1
      endif
      i = modulo( ifc, 2 )      ! child cells are located in left or right on the parent grid
      j = modulo( jfc, 2 )
      k = modulo( kfc, 2 )
      Newrank(i,j,k,Ihash(gid,rank),rank) = ranknew
    end subroutine assign
    !-----------------------------------------------------------------------
    subroutine checkNewrank
      integer :: rank, n, gid, n_refined

      n_refined = (count(Finer /= Undefi) + count(Fine /=Undefi))*8
      if ( n_refined /= count(Newrank /=MPI_PROC_NULL) .or. &
           n_refined /= count(Ihash /= Undefi .and. BoolMap)*8 .or. &
           n_refined /= sum(ngrid_node) .or. &
           n_refined /= count(BoolMap)*8 ) then
         print *, 'error in checkNewrank'
         HALT
      end if

      if (count(Newrank /= MPI_PROC_NULL) /= sum(ngrid_node)) then
         print *, 'error: inconsisnten Newrank and ngrid_node',count(Newrank /= MPI_PROC_NULL),sum(ngrid_node)
         HALT
      endif

      do rank = 0, NPE -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            if ( Fine(n, rank) == Undefi .or. Finer(n, rank) == Undefi ) cycle
            if (any(Newrank(:,:,:,n,rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank'
               HALT
            endif
         end do
      end do

      do rank = 0, NPE -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            gid = Fine(n, rank)
            if (gid == Undefi) cycle
            if (.not. BoolMap(gid,rank)) cycle
            if (any(Newrank(:,:,:,Ihash(gid,rank),rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank - Fine', rank, gid
               HALT
            endif
         end do
      end do

      do rank = 0, NPE -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            gid = Finer(n, rank)
            if (gid == Undefi) cycle
            if (any(Newrank(:,:,:,Ihash(gid,rank),rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank - Finer', rank, gid
               HALT
            endif
         end do
      end do

    end subroutine checkNewrank
  end subroutine define_Newrank_smallArea
  !-------------------------------------------------------------------------
  ! get norder
  !-------------------------------------------------------------------------
  function refine_get_norder_L() result(norder)
    integer :: norder
    norder = refine_get_norder(0, 0, 0, (NGI_BASE)-1, (NGJ_BASE)-1, (NGK_BASE)-1)
    norder = norder + Levelf
  end function refine_get_norder_L
  !-------------------------------------------------------------------------
  ! 大きな領域のオーダリングを行う。(新方式。欠損領域のある hilbert 曲線を実装)
  !-------------------------------------------------------------------------
  subroutine define_Newrank_largeArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(IN) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: ranknew, gcount, norder
    integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
    integer,dimension(:),allocatable :: ifc, jfc, kfc
    integer :: ngrid_node(0:NPE-1)
    ! 各ノードが分担するグリッド数
    call get_ngrid_node(ngrid_node)
    ! fill space
    ranknew = 0                   ! new rank id
    myrank = get_myrank()         ! myrank
    gcount = 0                    ! counter of grid for each rank
    norder = refine_get_norder_L()
    allocate(ifc(norder+1), jfc(norder+1), kfc(norder+1))
    ifc(:) = 0 ; jfc(:) = 0 ; kfc(:) = 0   ! point of filling curve
    call fillcells_L(norder,R,U,L,D,B,F)
    deallocate(ifc, jfc, kfc)
    if (maxval(NewRank) > NPE) then
       print *, 'error in define_Newrank (NewRank_max, NPE) =', maxval(NewRank), NPE
       stop
    endif

  contains
    !----------------------------------------------------------------------
    ! Peano-Hilbert filling curve
    !----------------------------------------------------------------------
    recursive subroutine fillcells_L(n, right, up, left, down, back, foward)
      integer,intent(IN) :: n, right, up, left, down, back, foward
      if (n == 0) return

      call pendown(n, foward,right,back,left,down,up) ! get a starting point from paraent level

      ! skip 
      if (.not. boolBlockExist(n, ifc(n), jfc(n), kfc(n))) return

      if (n == 1) call assign_L(ifc(n),jfc(n),kfc(n))

      call fillcells_L(n-1,foward,right,back,left,down,up)
      call connect_L(n, up)

      call fillcells_L(n-1,up,foward,down,back,left,right)
      call connect_L(n, right)

      call fillcells_L(n-1,up,foward,down,back,left,right) 
      call connect_L(n, down) 

      call fillcells_L(n-1,left,down,right,up,back,foward) 
      call connect_L(n, foward) 

      call fillcells_L(n-1,left,down,right,up,back,foward) 
      call connect_L(n, up) 

      call fillcells_L(n-1,up,back,down,foward,right,left) 
      call connect_L(n, left) 

      call fillcells_L(n-1,up,back,down,foward,right,left) 
      call connect_L(n, down) 

      call fillcells_L(n-1,back,right,foward,left,up,down) 
    end subroutine fillcells_L
    !----------------------------------------------------------------------
    ! define the starting point from position on the paranet level.
    !----------------------------------------------------------------------
    subroutine pendown(n, foward,right,back,left,down,up)
      integer,intent(IN) :: n, foward,right,back,left,down,up
      ifc(n) = 2*ifc(n+1)
      jfc(n) = 2*jfc(n+1)
      kfc(n) = 2*kfc(n+1)
      if (up == L .or. right == L .or. foward == L) ifc(n) = ifc(n) + 1
      if (up == B .or. right == B .or. foward == B) jfc(n) = jfc(n) + 1
      if (up == D .or. right == D .or. foward == D) kfc(n) = kfc(n) + 1
    end subroutine pendown
    !----------------------------------------------------------------------
    ! move a pen inside an octant.
    !----------------------------------------------------------------------
    subroutine connect_L(n, dir)
      integer,intent(IN) :: n, dir
      select case (dir)
      case (R)
         ifc(n) = ifc(n) + 1
      case (L)
         ifc(n) = ifc(n) - 1
      case (F)
         jfc(n) = jfc(n) + 1
      case (B)
         jfc(n) = jfc(n) - 1
      case (U)
         kfc(n) = kfc(n) + 1
      case (D)
         kfc(n) = kfc(n) - 1
      end select
!!$      if (n == 1) print *, ifc(n),jfc(n),kfc(n)
      if (n == 1) call assign_L(ifc(n),jfc(n),kfc(n))
    end subroutine connect_L
    !-------------------------------------------------------------------------
    ! Check whether a point (if, jf, kf) in order n has a grid.
    !-------------------------------------------------------------------------
    function boolBlockExist(n, if, jf, kf) result(boolExist)
      integer,intent(IN) :: n, if, jf, kf
      integer :: ic, jc, kc, gid, rank, levelTest
      logical :: boolExist
      levelTest = Levelc-n+1    ! grid level to be test

      boolExist = .TRUE.
      if (Levelc == Lmin) return
      if (levelTest < Lmin) return

      boolExist = .FALSE.
      if ( n == 1 .and.  &
           (if < igmin .or. if > igmax .or. &
            jf < jgmin .or. jf > jgmax .or. &
            kf < kgmin .or. kf > kgmax )) return

      ! grid coordinates in the parent level
      ic = IJKC(if, 0)
      jc = IJKC(jf, 0)
      kc = IJKC(kf, 0)
      ! check whether (ic,jc,kc) is inside the base grid if levelTest equals Lmin
      if (levelTest == Lmin) then
         if ( ic < 0 .or. ic > (NGI_BASE)-1 .or. &
              jc < 0 .or. jc > (NGJ_BASE)-1 .or. &
              kc < 0 .or. kc > (NGK_BASE)-1 ) then
            return
         end if
      end if

      call get_gid_from_ijkgrid(ic,jc,kc,levelTest,gid,rank) ! get gid and rank from coordinates and level
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( n > 1 .and. ChildGid(Left, Left, Left, gid, rank) == Undefi ) return
      if ( n == 1 ) then
         if (Finer(Ihash(gid,rank),rank) == Undefi .and. &
              Fine(Ihash(gid,rank),rank) == Undefi ) return
      endif
      boolExist = .TRUE.
    end function boolBlockExist
    !-----------------------------------------------------------------------
    subroutine assign_L(if,jf,kf)
      integer,intent(IN) :: if,jf,kf
      integer :: gid, rank, i, j, k, ic, jc, kc
!!$      if (get_myrank() == PRIMARY_RANK) print *, if,jf,kf
!!$      return
      if ( if < igmin .or. if > igmax .or. &
           jf < jgmin .or. jf > jgmax .or. &
           kf < kgmin .or. kf > kgmax ) return
      ! grid coordinates in the parent level
      ic = IJKC(if, 0)
      jc = IJKC(jf, 0)
      kc = IJKC(kf, 0)
      call get_gid_from_ijkgrid(ic,jc,kc,Levelc,gid,rank) ! get gid and rank from coordinates and level
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( Finer(Ihash(gid,rank),rank) == Undefi .and. &
           Fine(Ihash(gid,rank),rank) == Undefi ) return

      gcount = gcount + 1
      if ( gcount > ngrid_node(ranknew) ) then ! 次のランク
         ranknew = ranknew + 1
         gcount = 1
      endif
      i = modulo( if, 2 )      ! child cells are located in left or right on the parent grid
      j = modulo( jf, 2 )
      k = modulo( kf, 2 )
      Newrank(i,j,k,Ihash(gid,rank),rank) = ranknew
    end subroutine assign_L
  end subroutine define_Newrank_largeArea
  !-------------------------------------------------------------------------
  ! make new fine grids (Finer)
  !-------------------------------------------------------------------------
  subroutine do_finer
    use io_util, only : print_msg
    use string, only : num2char
    use packarr
    integer,dimension(Left:Right) :: is, ie, js, je, ks, ke ! lower and higher bound of sub-block
    ! size of sub-block
    integer,parameter :: Iminr=Imin, Jminr=Jmin, Kminr=Kmin
    integer,parameter :: Imaxr=(Imax-Imin+1)/2+Imin-1, Jmaxr=(Jmax-Jmin+1)/2+Jmin-1, Kmaxr=(Kmax-Kmin+1)/2+Kmin-1
    integer,parameter :: Iminrgh=Iminr-Rgh, Jminrgh=Jminr-Rgh, Kminrgh=Kminr-Rgh
    integer,parameter :: Imaxrgh=Imaxr+Rgh, Jmaxrgh=Jmaxr+Rgh, Kmaxrgh=Kmaxr+Rgh
    ! buffer for sub-block
    real(kind=DBL_KIND),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh, Mmin:Mmax) :: ubuf
    real(kind=DBL_KIND),dimension(Iminrgh:Imaxrgh) :: xbuf
    real(kind=DBL_KIND),dimension(Jminrgh:Jmaxrgh) :: ybuf
    real(kind=DBL_KIND),dimension(Kminrgh:Kmaxrgh) :: zbuf
    integer :: igridbuf, jgridbuf, kgridbuf

    integer :: gid, n, gids, gidd, ngrid, prank, ranks, rankd, i, j, k
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z


    ngrid = count( Finer /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Finer)   ' // num2char(ngrid))
    myrank = get_myrank()

    ! -------------------
    ! parameter for level
    ! -------------------
    Time(levelf) = Time(Levelc)
    Step(levelf) = Step(Levelc)
    CellWidth(:,Levelf) = CellWidth(:,Levelc)/2

    ! -----------------------------------------
    ! 袖の値は未定義になるので、袖の時刻をずらしておく
    ! -----------------------------------------
    U_StepNumberGhostCell(Levelf) = U_StepNumber(Levelf) - 1

    ! ブロックの分割サイズ
    is = (/ Imin  -Rgh, Imaxr+1-Rgh /)
    ie = (/ Imaxr +Rgh, Imax   +Rgh /)
    js = (/ Jmin  -Rgh, Jmaxr+1-Rgh /)
    je = (/ Jmaxr +Rgh, Jmax   +Rgh /)
    ks = (/ Kmin  -Rgh, Kmaxr+1-Rgh /)
    ke = (/ Kmaxr +Rgh, Kmax   +Rgh /)

    call pkar_reset

    do prank = lbound(Finer,2), ubound(Finer,2)
       do n = lbound(Finer,1), GidListNodeMax(Levelc, prank)
          if ( Finer(n, prank) == Undefi ) cycle
          gids = Finer(n, prank)
          ranks = prank
          if (myrank == ranks) then
             u => get_Up(gids)
             x => get_Xp(gids)
             y => get_Yp(gids)
             z => get_Zp(gids)
          endif
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   rankd = Newrank(i,j,k,n, prank)
                   if (myrank == ranks) then
                      ubuf(:,:,:,:) = u(is(i):ie(i),js(j):je(j),ks(k):ke(k),:)
                      xbuf(:) = x(is(i):ie(i))
                      ybuf(:) = y(js(j):je(j))
                      zbuf(:) = z(ks(k):ke(k))
                   endif
                   PACK_SEND4(ubuf, myrank, ranks, rankd)
                   PACK_SEND1(xbuf, myrank, ranks, rankd)
                   PACK_SEND1(ybuf, myrank, ranks, rankd)
                   PACK_SEND1(zbuf, myrank, ranks, rankd)
                   PACK_SEND0(Igrid(gids), myrank, ranks, rankd)
                   PACK_SEND0(Jgrid(gids), myrank, ranks, rankd)
                   PACK_SEND0(Kgrid(gids), myrank, ranks, rankd)
                enddo
             enddo
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do prank = lbound(Finer,2), ubound(Finer,2)
       do n = lbound(Finer,1), GidListNodeMax(Levelc, prank)
          if ( Finer(n, prank) == Undefi ) cycle
          gids = Finer(n, prank)
          ranks = prank
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   rankd = Newrank(i,j,k,n, prank)
                   if ( myrank == rankd ) then ! 受信
                      UNPACK_RECV4(ubuf, myrank, ranks, rankd)
                      UNPACK_RECV1(xbuf, myrank, ranks, rankd)
                      UNPACK_RECV1(ybuf, myrank, ranks, rankd)
                      UNPACK_RECV1(zbuf, myrank, ranks, rankd)
                      UNPACK_RECV0(igridbuf, myrank, ranks, rankd)
                      UNPACK_RECV0(jgridbuf, myrank, ranks, rankd)
                      UNPACK_RECV0(kgridbuf, myrank, ranks, rankd)
                      call interpU( ubuf, xbuf, ybuf, zbuf, igridbuf, jgridbuf, kgridbuf, i, j, k )
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

  contains
    !-------------------------------------------------------------------------
    ! interp U by gid
    !-------------------------------------------------------------------------
    subroutine interpU( uc, xc, yc, zc, igc, jgc, kgc, lri, lrj, lrk )
      use eos
      real(kind=DBL_KIND),intent(IN),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh, Mmin:Mmax) ::  uc
      real(kind=DBL_KIND),intent(IN),dimension(Iminrgh:Imaxrgh) :: xc
      real(kind=DBL_KIND),intent(IN),dimension(Jminrgh:Jmaxrgh) :: yc
      real(kind=DBL_KIND),intent(IN),dimension(Kminrgh:Kmaxrgh) :: zc
      integer,intent(IN) :: igc, jgc, kgc, lri, lrj, lrk
      real(kind=DBL_KIND),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh,MX:MZ) :: grad
      real(kind=DBL_KIND) :: hf, hc, di, dj, dk, dvc, dvf, x0, y0, z0
      integer :: i,j,k,m, ic,jc,kc, gidf
      integer :: ig0, jg0, kg0                ! origin of grid coordinates

      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: uf
      real(kind=DBL_KIND),dimension(:),pointer :: xf, yf, zf
      real(kind=DBL_KIND) :: a, b, minmod
#ifdef MPSI
      integer :: ic1, jc1, kc1
#endif
      
#define MINMOD(a,b) ( sign(1.d0,(a))*max(0.d0,min(abs(a),sign(1.d0,(a))*(b))) )

      ! -----------------
      ! a new fine block
      ! -----------------
      gidf = alloc_U(Levelf)    ! allocate a new block in finer level
      uf => get_Up(gidf)
      ! --------------------------
      ! interp physical variables
      ! --------------------------
      hf = 1.d0                   ! fine cell width
      hc = 2* hf                  ! coarse cell width
      dvc = 1.d0                  ! dummy volume
      dvf = 1.d0                  ! dummy volume
      call conv_u2w( uc, dvc )
      ! get offset of fine and coarse grids
      ig0 = Imin
      jg0 = Jmin
      kg0 = Kmin

      do m = Mmin, Mmax
#ifdef MPSI
         if (m == MPSI .or. m == MGX .or. m == MGY .or. m == MGZ) cycle
#endif
         ! ------------
         ! 勾配を求める
         ! ------------
         do k = Kminr, Kmaxr
            do j = Jminr, Jmaxr
               do i = Iminr, Imaxr
                  grad(i,j,k,MX) = MINMOD(uc(i+1,j,k,m)-uc(i,j,k,m), uc(i,j,k,m)-uc(i-1,j,k,m))/hc
                  grad(i,j,k,MY) = MINMOD(uc(i,j+1,k,m)-uc(i,j,k,m), uc(i,j,k,m)-uc(i,j-1,k,m))/hc
                  grad(i,j,k,MZ) = MINMOD(uc(i,j,k+1,m)-uc(i,j,k,m), uc(i,j,k,m)-uc(i,j,k-1,m))/hc
                  ! for debug
!!$                grad(i,j,k,MX) = (uc(i+1,j,k,m)-uc(i-1,j,k,m))/hc/2
!!$                grad(i,j,k,MY) = (uc(i,j+1,k,m)-uc(i,j-1,k,m))/hc/2
!!$                grad(i,j,k,MZ) = (uc(i,j,k+1,m)-uc(i,j,k-1,m))/hc/2
               enddo
            enddo
         enddo
         ! This code assumes that odd and even indexes are left and right cells
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  ! indexes of coarse cell
                  ic = IJKC( i, ig0 )
                  jc = IJKC( j, jg0 )
                  kc = IJKC( k, kg0 )
                  ! deviation of fine grid from the coase cell
                  di = ( modulo(i,2) - 0.5d0 )*hf
                  dj = ( modulo(j,2) - 0.5d0 )*hf
                  dk = ( modulo(k,2) - 0.5d0 )*hf
                  uf(i,j,k,m) = uc(ic,jc,kc,m)+grad(ic,jc,kc,MX)*di+grad(ic,jc,kc,MY)*dj+grad(ic,jc,kc,MZ)*dk
               enddo
            end do
         end do
      end do
      ! -------------------------------------
      ! bi-linear for gravitational potential
      ! -------------------------------------
#define NEAR(n, nc)  ( (nc)+2*modulo((n),2)-1 )
#ifdef MPSI
      do m = Mmin, Mmax
         if (.not. (m == MPSI .or. m == MGX .or. m == MGY .or. m == MGZ) ) cycle
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  ! indexes of coarse cell
                  ic = IJKC( i, ig0 )
                  jc = IJKC( j, jg0 )
                  kc = IJKC( k, kg0 )
                  ! the next nearest cell
                  ic1 = NEAR( i, ic )
                  jc1 = NEAR( j, jc )
                  kc1 = NEAR( k, kc )
                  uf(i,j,k,m) = (27.d0 * uc(ic,jc,kc,m)  &
                       + 9.d0 * (uc(ic1,jc,kc,m) + uc(ic,jc1,kc,m) + uc(ic,jc,kc1,m)) &
                       + 3.d0 * (uc(ic,jc1,kc1,m) + uc(ic1,jc,kc1,m) + uc(ic1,jc1,kc,m)) &
                       + uc(ic1,jc1,kc1,m))/64.d0
               enddo
            enddo
         end do
      end do
#endif

    call conv_w2u( uf, dvf )

      ! -----------------------------
      ! interpolation of coordinates
      ! -----------------------------
      x0 = xc(Imin) - CellWidth(MX,Levelf)/2
      y0 = yc(Jmin) - CellWidth(MY,Levelf)/2
      z0 = zc(Kmin) - CellWidth(MZ,Levelf)/2
      xf => get_Xp( gidf )
      yf => get_Yp( gidf )
      zf => get_Zp( gidf )
      do i = Imingh, Imaxgh
         xf(i) = (i-Imin)*CellWidth(MX,Levelf) + x0
      enddo
      do j = Jmingh, Jmaxgh
         yf(j) = (j-Jmin)*CellWidth(MY,Levelf) + y0
      enddo
      do k = Kmingh, Kmaxgh
         zf(k) = (k-Kmin)*CellWidth(MZ,Levelf) + z0
      enddo

      ! ---------------------------
      ! define Igrid, Jgrid, Kgrid
      ! ---------------------------
      Igrid(gidf) = IJKF(igc,0) + lri ! position of chlid grid
      Jgrid(gidf) = IJKF(jgc,0) + lrj
      Kgrid(gidf) = IJKF(kgc,0) + lrk

    end subroutine interpU

  end subroutine do_finer
  !-------------------------------------------------------------------------
  ! free unsed grid (Coarser)
  !-------------------------------------------------------------------------
  subroutine do_coarser
    use io_util, only : print_msg
    use string, only : num2char
    use unit ! KS DEBUG
    integer :: pgid, cgid, crank, n, i, j, k, prank, ngrid, gid
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z    ! KS DEBUG
    !  node global process
    ngrid = count( Coarser /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Coarser) ' // num2char(ngrid) )
    myrank = get_myrank()
    do prank = lbound(Coarser,2), ubound(Coarser,2)
       do n = lbound(Coarser,1), GidListNodeMax(Levelc, prank)
          if ( Coarser(n, prank) == Undefi ) cycle
          pgid = Coarser(n, prank) ! parent gid

          !------------------- KS DEBUG --------------------------!
          ! if (myrank == prank) then
          !    x => get_Xp(pgid)
          !    y => get_Yp(pgid)
          !    z => get_Zp(pgid)          
          !    print '(A,2I6,1P6E15.7)', 'refine grids: Coarser prank, pgid =', prank, pgid,&
          !         x(Imin)*Unit_au, x(Imax)*Unit_au, y(Jmin)*Unit_au, y(Jmax)*Unit_au, &
          !        z(Kmin)*Unit_au, z(Kmax)*Unit_au
          ! end if
          
          ! if ((prank==318 .and. pgid == 106) .or. (prank==318 .and. pgid == 107) &
          !      .or. (prank==232 .and. pgid == 123) .or. (prank==232 .and. pgid == 138) &
          !      .or. (prank==275 .and. pgid == 110) .or. (prank==276 .and. pgid == 132)) then
          !    if (myrank == PRIMARY_RANK) &
          !    print '(A,2I6,A)', 'refine grids: prank, pgid =', prank, pgid,&
          !         ' => skip destruction'          
          !    cycle
          ! end if

          ! if (Levelc==11) then
          !    if (myrank == PRIMARY_RANK) &
          !     print '(A,2I6,A)', 'refine grids: prank, pgid =', prank, pgid, &
          !     ' => skip destruction'
          !    print *, prank, pgid, ChildGid(0,0,0,pgid,prank)==Undefi
          !    cycle
          ! end if
          !------------------- KS DEBUG --------------------------!

          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   cgid = ChildGid(i,j,k,pgid,prank)
                   crank = ChildRank(i,j,k,pgid,prank)
!!$                   ! 親から見た子へのリンクを消す
!!$                   ChildGid(i,j,k,pgid,prank) = Undefi
!!$                   ChildRank(i,j,k,pgid,prank) = MPI_PROC_NULL
!!$                   ! 子からみた親へのリンクを消す
!!$                   ParentGid(cgid,crank) = Undefi
!!$                   ParentRank(cgid,crank) = MPI_PROC_NULL
                   if ( myrank == crank ) call dealloc_U( cgid )
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine do_coarser
  !-------------------------------------------------------------------------
  ! Fine リストにしたがって,グリッドのノードを変更する。
  !-------------------------------------------------------------------------
  subroutine do_fine
    use io_util, only : print_msg
    use string, only : num2char
    use packarr
    integer :: pgid, prank, gids, gidd, ranks, rankd, n, i, j, k, ngrid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    !  node global process
    ngrid = count( Fine /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Fine)    ' // num2char( ngrid ) )


    myrank = get_myrank()
    call pkar_reset
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) ! parent gid
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   gids = ChildGid(i,j,k,pgid,prank)
                   ranks = ChildRank(i,j,k,pgid,prank)
                   rankd = Newrank(i,j,k,n,prank)
                   ! (gids,ranks) -> (new gid, rankd)
                   if ( ranks == rankd ) cycle
                   if ( ranks == myrank ) then ! 送信
                      u => get_Up(gids)
                      x => get_Xp(gids)
                      y => get_Yp(gids)
                      z => get_Zp(gids)
                   endif
                   PACK_SEND4_SZ(u, myrank, ranks, rankd, (Imaxgh-Imingh+1)*(Jmaxgh-Jmingh+1)*(Kmaxgh-Kmingh+1)*(Mmax-Mmin+1))
                   PACK_SEND1_SZ(x, myrank, ranks, rankd, (Imaxgh-Imingh+1))
                   PACK_SEND1_SZ(y, myrank, ranks, rankd, (Jmaxgh-Jmingh+1))
                   PACK_SEND1_SZ(z, myrank, ranks, rankd, (Kmaxgh-Kmingh+1))
                   PACK_SEND0(Igrid(gids), myrank, ranks, rankd)
                   PACK_SEND0(Jgrid(gids), myrank, ranks, rankd)
                   PACK_SEND0(Kgrid(gids), myrank, ranks, rankd)
                enddo
             enddo
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) ! parent gid
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   gids = ChildGid(i,j,k,pgid,prank)
                   ranks = ChildRank(i,j,k,pgid,prank)
                   rankd = Newrank(i,j,k,n,prank)
                   ! (gids,ranks) -> (new gid, rankd)
                   if ( ranks == rankd ) cycle
                   if ( rankd == myrank ) then ! 受信
                      gidd = alloc_U( Levelf )
                      u => get_Up(gidd)
                      x => get_Xp(gidd)
                      y => get_Yp(gidd)
                      z => get_Zp(gidd)
                      UNPACK_RECV4(u, myrank, ranks, rankd)
                      UNPACK_RECV1(x, myrank, ranks, rankd)
                      UNPACK_RECV1(y, myrank, ranks, rankd)
                      UNPACK_RECV1(z, myrank, ranks, rankd)
                      UNPACK_RECV0(Igrid(gidd), myrank, ranks, rankd)
                      UNPACK_RECV0(Jgrid(gidd), myrank, ranks, rankd)
                      UNPACK_RECV0(Kgrid(gidd), myrank, ranks, rankd)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    ! -----------------
    ! Free unused grid
    ! -----------------
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) ! parent gid
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   gids = ChildGid(i,j,k,pgid,prank)
                   ranks = ChildRank(i,j,k,pgid,prank)
                   rankd = Newrank(i,j,k,n,prank)
                   if ( ranks == rankd ) cycle
                   if ( myrank /= ranks ) cycle
                   call dealloc_U(gids)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine do_fine
  !-------------------------------------------------------------------------
  ! 2次元をエミュレーションするために、BoolMapを変更する (Kgrid = 0 のグリッドだけをtrue)
  !-------------------------------------------------------------------------
  subroutine Emulate_2Dim_BoolMap
    integer :: n, gid
    logical,dimension(Gidmin:Gidmax,0:NPE-1) :: bufboolmap
    myrank = get_myrank()
    do n = lbound(GidList,1), GidListMax( Levelc )
       gid = GidList(n, Levelc)
       if ( Kgrid(gid) == 0 ) cycle
       BoolMap( gid, myrank ) = .false.
    end do
    call mpi_allreduce_BoolMap(MPI_LAND)
  end subroutine Emulate_2Dim_BoolMap
end module refine

