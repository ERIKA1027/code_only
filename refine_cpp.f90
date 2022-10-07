module refine
  use grid
  use mpilib
  implicit none
  private
  logical,save,dimension(Gidmin:Gidmax,0:400 -1) :: BoolMap
  integer,save,dimension(Gidmin:Gidmax,0:400 -1) :: Finer, Coarser, Fine, Coarse, Ihash
  integer,save,dimension(Left:Right,Left:Right,Left:Right,Gidmin:Gidmax,0:400 -1) :: Newrank
  integer,save :: Levelf, Levelc 
  logical,save :: BoolRefine 
  integer,parameter :: Rgh = 1 
  public :: refineLevel, refineAllLevel, refine_get_norder
contains
  subroutine refineALlLevel
    logical :: bool
    integer :: level
    do level = Lmin+1, Lmax
       call refineLevel(level, bool)
       if (.not. bool) exit
    enddo
  end subroutine refineALlLevel
  subroutine refineLevel(level, bool)
    use io_util, only : print_msg
    use string, only : num2char
    use grid, only : alloc_U1order, alloc_U2order, dealloc_U1order, dealloc_U2order
    use grid_boundary
    use boundary
    use sinkParticle
    use modelParameter, only : MP_CONNECTION_RUN 
    integer,intent(IN) :: level
    logical,intent(OUT) :: bool 
    integer :: n
    globdbg_gridupdate = .False.
    Levelf = level
    Levelc = level - 1
    BoolRefine = .FALSE.
    bool = .FALSE.
    if ( Levelc < Lmin ) return
    if ( Levelf > Lmax ) return
    if (MP_CONNECTION_RUN == 0) then 
       if ( Levelf > sp_getLevel() ) return
    end if 
    call print_msg( 'refined level = '//num2char(level) )
    call boundary_grid( Levelc , COMPLETE)
    do n = Gidmin, GidListMax( Levelc )
       call boundary_u( GidList(n, Levelc), COMPLETE )
    enddo
    call makeBoolRefine
    bool = BoolRefine
    globdbg_gridupdate = bool_update_grid()
    if ( .not. bool_update_grid() ) return
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
       LevelMax = max(LevelMax-1, levelc)
    endif
    call update_gidlist(Levelf)
    call clear_linklist(Levelf)
    call update_parent_linklist(Levelf)
    if (Levelf + 1 <= LevelMax) then
       call update_parent_linklist(Levelf+1)
    endif
    call update_neighbor_linklist(Levelf)
    if ( Levelf <= LevelMax ) then
       do n = Gidmin, GidListMax( Levelf )
          call alloc_U1order( GidList(n, Levelf) )
          call alloc_U2order( GidList(n, Levelf) )
       enddo
    endif
  end subroutine refineLevel
  subroutine makeBoolRefine()
    use eos
    use refineCond
    use sinkParticle
    use unit 
    use modelParameter, only : MP_CONNECTION_RUN 
    integer :: gid, n,m,lr, gidn, rankn, rank
    logical :: bufbool, boolself, boolgc
    logical,dimension(Gidmin:Gidmax,0:400 -1) :: bufboolmap
    integer :: lr1, lr2, lr3, m1, m2, m3, &
         gid1, gid2, gid3, gid4, gid5, gid6, gid7, &
         rank1, rank2, rank3, rank4, rank5, rank6, rank7, &
         i, j, k, lbuf, ic, jc, kc, gidc, rankc
    real(kind=8),dimension(:),pointer :: x, y, z 
    myrank = get_myrank()
    BoolMap(:,:) = .FALSE.
    do n = lbound(GidList,1), GidListMax( Levelc )
       boolself = .false.
       boolgc = .false.
       gid = GidList(n, Levelc)
       do k = Left, Right
          do j = Left, Right
             do i = Left, Right
                gidc = ChildGid(i,j,k,gid,myrank)
                rankc = ChildRank(i,j,k,gid,myrank)
                if ( gidc == Undefi ) cycle
                if ( any( ChildGid(:,:,:, gidc, rankc) /= Undefi) ) &
                     boolgc = .true.
             enddo
          enddo
       enddo
       boolself = refineCond_eval(gid)
       call sp_refineCond_KS(gid, boolself) 
       if (.not. (boolgc .or. boolself) ) cycle
       BoolMap(gid, myrank) = boolgc .or. boolself
       if ( .not. ( boolself .or. have_allgrandchild(gid, myrank) )) then
          do k = Left, Right
             do j = Left, Right
                do i = Left, Right
                   gidc = ChildGid(i,j,k,gid,myrank)
                   rankc = ChildRank(i,j,k,gid,myrank)
                   if ( gidc == Undefi ) cycle
                   if ( any( ChildGid(:,:,:, gidc, rankc) /= Undefi) ) then
                      gid1 = NeighborGid(i,0,gid,myrank)
                      rank1 = NeighborRank(i,0,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.
                      gid1 = NeighborGid(j,1,gid,myrank)
                      rank1 = NeighborRank(j,1,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.
                      gid1 = NeighborGid(k,2,gid,myrank)
                      rank1 = NeighborRank(k,2,gid,myrank)
                      if (.not. (gid1 == Undefi .or. rank1 == MPI_PROC_NULL) ) &
                           BoolMap( gid1, rank1 ) = .TRUE.
                   end if
                enddo
             enddo
          enddo
          cycle 
       end if
       do m1 = 0, 2
          do lr1 = Left , Right 
             gid1 = NeighborGid(lr1,m1,gid,myrank)
             rank1 = NeighborRank(lr1,m1,gid,myrank)
             if (gid1 == Undefi .or. rank1 == MPI_PROC_NULL ) cycle
             BoolMap( gid1, rank1 ) = .TRUE.
             do m2 = 0, 2 
                if ( m2 == m1 ) cycle
                do lr2 = Left, Right
                   gid2 = NeighborGid(lr2,m2,gid1,rank1)
                   rank2 = NeighborRank(lr2,m2,gid1,rank1)
                   if (gid2 == Undefi .or. rank2 == MPI_PROC_NULL ) cycle
                   BoolMap( gid2, rank2 ) = .TRUE.
                   do m3 = 0, 2 
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
          enddo
       enddo
    enddo
    call mpi_allreduce_BoolMap(MPI_LOR)
    do n = lbound(GidList,1), GidListMax( Levelc )
       gid = GidList(n, Levelc)
       if (.not. checkNeighbor(gid, myrank)) then
          BoolMap(gid, myrank) = .FALSE.
       endif
    end do
    call mpi_allreduce_BoolMap(MPI_LAND) 
    BoolRefine = any( BoolMap )
    call mpi_allreduce( BoolRefine, bufbool, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    BoolRefine = bufbool
    Finer(:,:) = Undefi
    Coarser(:,:) = Undefi
    Fine(:,:) = Undefi
    Coarse(:,:) = Undefi
    Ihash(:,:) = Undefi
    Newrank(:,:,:,:,:) = MPI_PROC_NULL
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          gid = GidListNode(n, Levelc, rank)
          Ihash(gid,rank) = n 
          if ( BoolMap(gid,rank) ) then
             if ( ChildGid(0,0,0,gid,rank) == Undefi ) then
                Finer(n,rank) = gid 
             else
                Fine(n,rank) = gid 
             endif
          else
             if ( ChildGid(0,0,0,gid,rank) == Undefi ) then
                Coarse(n,rank) = gid 
             else
                Coarser(n,rank) = gid 
             endif
          endif
       enddo
    enddo
    if ( get_myrank() == 0 ) then
       if (mod(Step(Lmin),100) == 0 .or. MP_CONNECTION_RUN > 0) then
          print '(A,4I6)', 'refine grids: Coarser, Coarse, Fine, Finer =>  ', &
               count( Coarser /= Undefi ), count( Coarse /= Undefi ), &
               count( Fine /= Undefi ), count( Finer /= Undefi )
       endif
    endif
  contains
    function checkNeighbor(gid,rank) result(bool)
      integer,intent(IN) :: gid, rank
      logical :: bool
      integer :: lr, ndir, pgid, prank
      bool = .true.
      if (Levelc == Lmin) return
      pgid = ParentGid(gid,rank)
      prank = ParentRank(gid,rank)
      do ndir = 0, 2
         do lr = Left, Right
            if ( NeighborGid(lr, ndir, pgid, prank) /= Undefi .and. & 
                 NeighborGid(lr, ndir, gid, rank) == Undefi ) & 
                 bool = .false.
         enddo
      enddo
    end function checkNeighbor
    function checkNeighbor_NANAME(gid,rank) result(bool)
      integer,intent(IN) :: gid, rank
      logical :: bool
      integer :: pgid, prank, m1,m2,m3,lr1,lr2,lr3,gid1,gid2,gid3,rank1,rank2,rank3
      bool = .true.
      if (Levelc == Lmin) return
      pgid = ParentGid(gid,rank)
      prank = ParentRank(gid,rank)
      do m1 = 0, 2 
         do lr1 = Left , Right
            gid1 = NeighborGid(lr1,m1,gid,rank)
            rank1 = NeighborRank(lr1,m1,gid,rank)
            if (gid1 == Undefi ) bool = .false. 
            do m2 = 0, 2 
               if ( m2 == m1 ) cycle
               do lr2 = Left, Right
                  gid2 = NeighborGid(lr2,m2,gid1,rank1)
                  rank2 = NeighborRank(lr2,m2,gid1,rank1)
                  if (gid2 == Undefi ) bool = .false. 
                   do m3 = 0, 2 
                      if ( m3 == m2 .or. m3 == m1 ) cycle
                      do lr3 = Left, Right
                         gid3 = NeighborGid(lr3,m3,gid2,rank2)
                         rank3 = NeighborRank(lr3,m3,gid2,rank2)
                         if (gid3 == Undefi ) bool = .false. 
                      end do
                   end do
                end do
             end do
          end do
       end do
     end function checkNeighbor_NANAME
  end subroutine makeBoolRefine
  function bool_update_grid() result(bool)
    logical :: bool
    bool = ( any(Finer /= Undefi ) .or. any(Coarser /= Undefi ) )
  end function bool_update_grid
  function have_allgrandchild(gid, rank) result(bool)
    integer,intent(IN) :: gid, rank
    logical :: bool
    integer :: i, j, k, gidc, rankc
    bool = .true.
    do k = Left, Right 
       do j = Left, Right
          do i = Left, Right
             gidc = ChildGid(i,j,k,gid,rank)
             rankc = ChildRank(i,j,k,gid,rank)
             if ( gidc == Undefi ) then 
                bool = .false.
                return
             endif
             if ( any(ChildGid(:,:,:,gidc,rankc) == Undefi) ) then 
                bool = .false.
                return
             endif
          enddo
       enddo
    enddo
  end function have_allgrandchild
  subroutine mpi_allreduce_BoolMap(operator_mpi)
    integer,intent(IN) :: operator_mpi
    logical,dimension(:),allocatable :: buf
    integer :: bufsize, rank, pos, n, gid
    bufsize = sum( max((GidListNodeMax(Levelc, :)-Gidmin+1),0) )
    if (bufsize == 0) return
    allocate(buf(bufsize))
    pos = 1
    do rank = 0, 400 -1
       if ( GidListNodeMax(Levelc, rank) == Undefi ) cycle
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          gid = GidListNode(n, Levelc, rank)
          if (gid == Undefi) cycle
          buf(pos) = BoolMap(gid, rank)
          pos = pos + 1
       end do
    end do
    if (pos-1 /= bufsize) print *, '**** error in mpi_allreduce_BoolMap', pos-1, bufsize 
    call mpi_allreduce( MPI_IN_PLACE, buf, bufsize, MPI_LOGICAL, operator_mpi, MPI_COMM_WORLD, ierr)
    pos = 1
    do rank = 0, 400 -1
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
  subroutine define_Newrank
    integer :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer,parameter :: LARGE_LENGTH = 2**10
    call grid_range(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    call define_Newrank_largeArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
  end subroutine define_Newrank
  subroutine grid_range(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(OUT) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: rank, n, gid
    integer,dimension(0:2) :: buf, bufr 
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
    buf = (/igmin, jgmin, kgmin/)
    call mpi_allreduce(buf, bufr, size(buf), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr )
    igmin = bufr(0)
 jgmin = bufr(1)
 kgmin = bufr(2)
    buf = (/igmax, jgmax, kgmax/)
    call mpi_allreduce(buf, bufr, size(buf), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
    igmax = bufr(0)
 jgmax = bufr(1)
 kgmax = bufr(2)
    igmin = ((igmin) - (0)) * 2 + (0) 
    jgmin = ((jgmin) - (0)) * 2 + (0)
    kgmin = ((kgmin) - (0)) * 2 + (0)
    igmax = ((igmax) - (0)) * 2 + (0) + 1 
    jgmax = ((jgmax) - (0)) * 2 + (0) + 1
    kgmax = ((kgmax) - (0)) * 2 + (0) + 1
  end subroutine grid_range
  subroutine get_ngrid_node(ngrid_node)
    integer,intent(OUT) :: ngrid_node(0:400 -1)
    integer :: ngrid_total, rank, n
    ngrid_total = 0
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelc, rank)
          if ( Fine(n, rank) /= Undefi .or. Finer(n, rank) /= Undefi ) &
               ngrid_total = ngrid_total + 1
       enddo
    enddo
    ngrid_total = ngrid_total * 8 
    ngrid_node(:) = int(ngrid_total/(400))
    do n = 1, mod( ngrid_total, 400 )
       ngrid_node(400 -n) = ngrid_node(400 -n) + 1
    enddo
  end subroutine get_ngrid_node
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
  subroutine define_Newrank_smallArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(IN) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: ranknew, gcount, norder
    integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
    integer :: ifc, jfc, kfc
    integer :: ngrid_node(0:400 -1)
    call get_ngrid_node(ngrid_node)
    ranknew = 0 
    myrank = get_myrank() 
    gcount = 0 
    ifc = igmin 
 jfc = jgmin 
 kfc = kgmin 
    call assign(ifc,jfc,kfc)
    norder = refine_get_norder(igmin,jgmin,kgmin,igmax,jgmax,kgmax)
    call fillcells(norder, R,U,L,D,B,F)
    if (maxval(NewRank) > 400) then
       print *, 'error in define_Newrank (NewRank_max, NPE) =', maxval(NewRank), 400
       stop
    endif
  contains
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
    subroutine assign(ifc,jfc,kfc)
      integer,intent(IN) :: ifc,jfc,kfc
      integer :: gid, rank, i, j, k, ic, jc, kc
      if ( ifc < igmin .or. ifc > igmax .or. &
           jfc < jgmin .or. jfc > jgmax .or. &
           kfc < kgmin .or. kfc > kgmax ) return
      ic = ((ifc)-(0))/2 + mod(min((ifc)-(0),0),2) + (0)
      jc = ((jfc)-(0))/2 + mod(min((jfc)-(0),0),2) + (0)
      kc = ((kfc)-(0))/2 + mod(min((kfc)-(0),0),2) + (0)
      call get_gid_from_ijkgrid(ic,jc,kc,Levelc,gid,rank) 
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( Finer(Ihash(gid,rank),rank) == Undefi .and. &
           Fine(Ihash(gid,rank),rank) == Undefi ) return
      gcount = gcount + 1
      if ( gcount > ngrid_node(ranknew) ) then 
         ranknew = ranknew + 1
         gcount = 1
      endif
      i = modulo( ifc, 2 ) 
      j = modulo( jfc, 2 )
      k = modulo( kfc, 2 )
      Newrank(i,j,k,Ihash(gid,rank),rank) = ranknew
    end subroutine assign
    subroutine checkNewrank
      integer :: rank, n, gid, n_refined
      n_refined = (count(Finer /= Undefi) + count(Fine /=Undefi))*8
      if ( n_refined /= count(Newrank /=MPI_PROC_NULL) .or. &
           n_refined /= count(Ihash /= Undefi .and. BoolMap)*8 .or. &
           n_refined /= sum(ngrid_node) .or. &
           n_refined /= count(BoolMap)*8 ) then
         print *, 'error in checkNewrank'
         write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

      end if
      if (count(Newrank /= MPI_PROC_NULL) /= sum(ngrid_node)) then
         print *, 'error: inconsisnten Newrank and ngrid_node',count(Newrank /= MPI_PROC_NULL),sum(ngrid_node)
         write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

      endif
      do rank = 0, 400 -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            if ( Fine(n, rank) == Undefi .or. Finer(n, rank) == Undefi ) cycle
            if (any(Newrank(:,:,:,n,rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank'
               write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

            endif
         end do
      end do
      do rank = 0, 400 -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            gid = Fine(n, rank)
            if (gid == Undefi) cycle
            if (.not. BoolMap(gid,rank)) cycle
            if (any(Newrank(:,:,:,Ihash(gid,rank),rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank - Fine', rank, gid
               write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

            endif
         end do
      end do
      do rank = 0, 400 -1
         do n = Gidmin, GidListNodeMax(Levelc, rank)
            gid = Finer(n, rank)
            if (gid == Undefi) cycle
            if (any(Newrank(:,:,:,Ihash(gid,rank),rank) == MPI_PROC_NULL)) then
               print *, 'error in Newrank - Finer', rank, gid
               write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

            endif
         end do
      end do
    end subroutine checkNewrank
  end subroutine define_Newrank_smallArea
  function refine_get_norder_L() result(norder)
    integer :: norder
    norder = refine_get_norder(0, 0, 0, (8)-1, (8)-1, (8)-1)
    norder = norder + Levelf
  end function refine_get_norder_L
  subroutine define_Newrank_largeArea(igmin, jgmin, kgmin, igmax, jgmax, kgmax)
    integer,intent(IN) :: igmin, jgmin, kgmin, igmax, jgmax, kgmax
    integer :: ranknew, gcount, norder
    integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
    integer,dimension(:),allocatable :: ifc, jfc, kfc
    integer :: ngrid_node(0:400 -1)
    call get_ngrid_node(ngrid_node)
    ranknew = 0 
    myrank = get_myrank() 
    gcount = 0 
    norder = refine_get_norder_L()
    allocate(ifc(norder+1), jfc(norder+1), kfc(norder+1))
    ifc(:) = 0 
 jfc(:) = 0 
 kfc(:) = 0 
    call fillcells_L(norder,R,U,L,D,B,F)
    deallocate(ifc, jfc, kfc)
    if (maxval(NewRank) > 400) then
       print *, 'error in define_Newrank (NewRank_max, NPE) =', maxval(NewRank), 400
       stop
    endif
  contains
    recursive subroutine fillcells_L(n, right, up, left, down, back, foward)
      integer,intent(IN) :: n, right, up, left, down, back, foward
      if (n == 0) return
      call pendown(n, foward,right,back,left,down,up) 
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
    subroutine pendown(n, foward,right,back,left,down,up)
      integer,intent(IN) :: n, foward,right,back,left,down,up
      ifc(n) = 2*ifc(n+1)
      jfc(n) = 2*jfc(n+1)
      kfc(n) = 2*kfc(n+1)
      if (up == L .or. right == L .or. foward == L) ifc(n) = ifc(n) + 1
      if (up == B .or. right == B .or. foward == B) jfc(n) = jfc(n) + 1
      if (up == D .or. right == D .or. foward == D) kfc(n) = kfc(n) + 1
    end subroutine pendown
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
      if (n == 1) call assign_L(ifc(n),jfc(n),kfc(n))
    end subroutine connect_L
    function boolBlockExist(n, if, jf, kf) result(boolExist)
      integer,intent(IN) :: n, if, jf, kf
      integer :: ic, jc, kc, gid, rank, levelTest
      logical :: boolExist
      levelTest = Levelc-n+1 
      boolExist = .TRUE.
      if (Levelc == Lmin) return
      if (levelTest < Lmin) return
      boolExist = .FALSE.
      if ( n == 1 .and. &
           (if < igmin .or. if > igmax .or. &
            jf < jgmin .or. jf > jgmax .or. &
            kf < kgmin .or. kf > kgmax )) return
      ic = ((if)-(0))/2 + mod(min((if)-(0),0),2) + (0)
      jc = ((jf)-(0))/2 + mod(min((jf)-(0),0),2) + (0)
      kc = ((kf)-(0))/2 + mod(min((kf)-(0),0),2) + (0)
      if (levelTest == Lmin) then
         if ( ic < 0 .or. ic > (8)-1 .or. &
              jc < 0 .or. jc > (8)-1 .or. &
              kc < 0 .or. kc > (8)-1 ) then
            return
         end if
      end if
      call get_gid_from_ijkgrid(ic,jc,kc,levelTest,gid,rank) 
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( n > 1 .and. ChildGid(Left, Left, Left, gid, rank) == Undefi ) return
      if ( n == 1 ) then
         if (Finer(Ihash(gid,rank),rank) == Undefi .and. &
              Fine(Ihash(gid,rank),rank) == Undefi ) return
      endif
      boolExist = .TRUE.
    end function boolBlockExist
    subroutine assign_L(if,jf,kf)
      integer,intent(IN) :: if,jf,kf
      integer :: gid, rank, i, j, k, ic, jc, kc
      if ( if < igmin .or. if > igmax .or. &
           jf < jgmin .or. jf > jgmax .or. &
           kf < kgmin .or. kf > kgmax ) return
      ic = ((if)-(0))/2 + mod(min((if)-(0),0),2) + (0)
      jc = ((jf)-(0))/2 + mod(min((jf)-(0),0),2) + (0)
      kc = ((kf)-(0))/2 + mod(min((kf)-(0),0),2) + (0)
      call get_gid_from_ijkgrid(ic,jc,kc,Levelc,gid,rank) 
      if ( gid == Undefi .or. rank == MPI_PROC_NULL ) return
      if ( Finer(Ihash(gid,rank),rank) == Undefi .and. &
           Fine(Ihash(gid,rank),rank) == Undefi ) return
      gcount = gcount + 1
      if ( gcount > ngrid_node(ranknew) ) then 
         ranknew = ranknew + 1
         gcount = 1
      endif
      i = modulo( if, 2 ) 
      j = modulo( jf, 2 )
      k = modulo( kf, 2 )
      Newrank(i,j,k,Ihash(gid,rank),rank) = ranknew
    end subroutine assign_L
  end subroutine define_Newrank_largeArea
  subroutine do_finer
    use io_util, only : print_msg
    use string, only : num2char
    use packarr
    integer,dimension(Left:Right) :: is, ie, js, je, ks, ke 
    integer,parameter :: Iminr=Imin, Jminr=Jmin, Kminr=Kmin
    integer,parameter :: Imaxr=(Imax-Imin+1)/2+Imin-1, Jmaxr=(Jmax-Jmin+1)/2+Jmin-1, Kmaxr=(Kmax-Kmin+1)/2+Kmin-1
    integer,parameter :: Iminrgh=Iminr-Rgh, Jminrgh=Jminr-Rgh, Kminrgh=Kminr-Rgh
    integer,parameter :: Imaxrgh=Imaxr+Rgh, Jmaxrgh=Jmaxr+Rgh, Kmaxrgh=Kmaxr+Rgh
    real(kind=8),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh, Mmin:Mmax) :: ubuf
    real(kind=8),dimension(Iminrgh:Imaxrgh) :: xbuf
    real(kind=8),dimension(Jminrgh:Jmaxrgh) :: ybuf
    real(kind=8),dimension(Kminrgh:Kmaxrgh) :: zbuf
    integer :: igridbuf, jgridbuf, kgridbuf
    integer :: gid, n, gids, gidd, ngrid, prank, ranks, rankd, i, j, k
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(:),pointer :: x, y, z
    ngrid = count( Finer /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Finer)   ' // num2char(ngrid))
    myrank = get_myrank()
    Time(levelf) = Time(Levelc)
    Step(levelf) = Step(Levelc)
    CellWidth(:,Levelf) = CellWidth(:,Levelc)/2
    U_StepNumberGhostCell(Levelf) = U_StepNumber(Levelf) - 1
    is = (/ Imin -Rgh, Imaxr+1-Rgh /)
    ie = (/ Imaxr +Rgh, Imax +Rgh /)
    js = (/ Jmin -Rgh, Jmaxr+1-Rgh /)
    je = (/ Jmaxr +Rgh, Jmax +Rgh /)
    ks = (/ Kmin -Rgh, Kmaxr+1-Rgh /)
    ke = (/ Kmaxr +Rgh, Kmax +Rgh /)
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
                   if ((myrank) == (ranks)) call pkar_push(ubuf(lbound(ubuf,1),lbound(ubuf,2),lbound(ubuf,3),lbound(ubuf,4)), size(&
&ubuf), kind(ubuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(ubuf), kind(ubuf), ranks)
                   if ((myrank) == (ranks)) call pkar_push(xbuf(lbound(xbuf,1)), size(xbuf), kind(xbuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(xbuf), kind(xbuf), ranks)
                   if ((myrank) == (ranks)) call pkar_push(ybuf(lbound(ybuf,1)), size(ybuf), kind(ybuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(ybuf), kind(ybuf), ranks)
                   if ((myrank) == (ranks)) call pkar_push(zbuf(lbound(zbuf,1)), size(zbuf), kind(zbuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(zbuf), kind(zbuf), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Igrid(gids), 1, kind(Igrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Igrid(gids)), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Jgrid(gids), 1, kind(Jgrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Jgrid(gids)), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Kgrid(gids), 1, kind(Kgrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Kgrid(gids)), ranks)
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
                   if ( myrank == rankd ) then 
                      if ((myrank) == (rankd)) call pkar_pop(ubuf(lbound(ubuf,1),lbound(ubuf,2),lbound(ubuf,3),lbound(ubuf,4)), siz&
&e(ubuf), kind(ubuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(xbuf(lbound(xbuf,1)), size(xbuf), kind(xbuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(ybuf(lbound(ybuf,1)), size(ybuf), kind(ybuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(zbuf(lbound(zbuf,1)), size(zbuf), kind(zbuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(igridbuf, 1, kind(igridbuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(jgridbuf, 1, kind(jgridbuf), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(kgridbuf, 1, kind(kgridbuf), ranks)
                      call interpU( ubuf, xbuf, ybuf, zbuf, igridbuf, jgridbuf, kgridbuf, i, j, k )
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
  contains
    subroutine interpU( uc, xc, yc, zc, igc, jgc, kgc, lri, lrj, lrk )
      use eos
      real(kind=8),intent(IN),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh, Mmin:Mmax) :: uc
      real(kind=8),intent(IN),dimension(Iminrgh:Imaxrgh) :: xc
      real(kind=8),intent(IN),dimension(Jminrgh:Jmaxrgh) :: yc
      real(kind=8),intent(IN),dimension(Kminrgh:Kmaxrgh) :: zc
      integer,intent(IN) :: igc, jgc, kgc, lri, lrj, lrk
      real(kind=8),dimension(Iminrgh:Imaxrgh, Jminrgh:Jmaxrgh, Kminrgh:Kmaxrgh,0:2) :: grad
      real(kind=8) :: hf, hc, di, dj, dk, dvc, dvf, x0, y0, z0
      integer :: i,j,k,m, ic,jc,kc, gidf
      integer :: ig0, jg0, kg0 
      real(kind=8),dimension(:,:,:,:),pointer :: uf
      real(kind=8),dimension(:),pointer :: xf, yf, zf
      real(kind=8) :: a, b, minmod
      integer :: ic1, jc1, kc1
      gidf = alloc_U(Levelf) 
      uf => get_Up(gidf)
      hf = 1.d0 
      hc = 2* hf 
      dvc = 1.d0 
      dvf = 1.d0 
      call conv_u2w( uc, dvc )
      ig0 = Imin
      jg0 = Jmin
      kg0 = Kmin
      do m = Mmin, Mmax
         if (m == 17 .or. m == 18 .or. m == 19 .or. m == 20) cycle
         do k = Kminr, Kmaxr
            do j = Jminr, Jmaxr
               do i = Iminr, Imaxr
                  grad(i,j,k,0) = ( sign(1.d0,(uc(i+1,j,k,m)-uc(i,j,k,m)))*max(0.d0,min(abs(uc(i+1,j,k,m)-uc(i,j,k,m)),sign(1.d0,(u&
&c(i+1,j,k,m)-uc(i,j,k,m)))*(uc(i,j,k,m)-uc(i-1,j,k,m)))) )/hc
                  grad(i,j,k,1) = ( sign(1.d0,(uc(i,j+1,k,m)-uc(i,j,k,m)))*max(0.d0,min(abs(uc(i,j+1,k,m)-uc(i,j,k,m)),sign(1.d0,(u&
&c(i,j+1,k,m)-uc(i,j,k,m)))*(uc(i,j,k,m)-uc(i,j-1,k,m)))) )/hc
                  grad(i,j,k,2) = ( sign(1.d0,(uc(i,j,k+1,m)-uc(i,j,k,m)))*max(0.d0,min(abs(uc(i,j,k+1,m)-uc(i,j,k,m)),sign(1.d0,(u&
&c(i,j,k+1,m)-uc(i,j,k,m)))*(uc(i,j,k,m)-uc(i,j,k-1,m)))) )/hc
               enddo
            enddo
         enddo
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  ic = ((i)-(ig0))/2 + mod(min((i)-(ig0),0),2) + (ig0)
                  jc = ((j)-(jg0))/2 + mod(min((j)-(jg0),0),2) + (jg0)
                  kc = ((k)-(kg0))/2 + mod(min((k)-(kg0),0),2) + (kg0)
                  di = ( modulo(i,2) - 0.5d0 )*hf
                  dj = ( modulo(j,2) - 0.5d0 )*hf
                  dk = ( modulo(k,2) - 0.5d0 )*hf
                  uf(i,j,k,m) = uc(ic,jc,kc,m)+grad(ic,jc,kc,0)*di+grad(ic,jc,kc,1)*dj+grad(ic,jc,kc,2)*dk
               enddo
            end do
         end do
      end do
      do m = Mmin, Mmax
         if (.not. (m == 17 .or. m == 18 .or. m == 19 .or. m == 20) ) cycle
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  ic = ((i)-(ig0))/2 + mod(min((i)-(ig0),0),2) + (ig0)
                  jc = ((j)-(jg0))/2 + mod(min((j)-(jg0),0),2) + (jg0)
                  kc = ((k)-(kg0))/2 + mod(min((k)-(kg0),0),2) + (kg0)
                  ic1 = ( (ic)+2*modulo((i),2)-1 )
                  jc1 = ( (jc)+2*modulo((j),2)-1 )
                  kc1 = ( (kc)+2*modulo((k),2)-1 )
                  uf(i,j,k,m) = (27.d0 * uc(ic,jc,kc,m) &
                       + 9.d0 * (uc(ic1,jc,kc,m) + uc(ic,jc1,kc,m) + uc(ic,jc,kc1,m)) &
                       + 3.d0 * (uc(ic,jc1,kc1,m) + uc(ic1,jc,kc1,m) + uc(ic1,jc1,kc,m)) &
                       + uc(ic1,jc1,kc1,m))/64.d0
               enddo
            enddo
         end do
      end do
    call conv_w2u( uf, dvf )
      x0 = xc(Imin) - CellWidth(0,Levelf)/2
      y0 = yc(Jmin) - CellWidth(1,Levelf)/2
      z0 = zc(Kmin) - CellWidth(2,Levelf)/2
      xf => get_Xp( gidf )
      yf => get_Yp( gidf )
      zf => get_Zp( gidf )
      do i = Imingh, Imaxgh
         xf(i) = (i-Imin)*CellWidth(0,Levelf) + x0
      enddo
      do j = Jmingh, Jmaxgh
         yf(j) = (j-Jmin)*CellWidth(1,Levelf) + y0
      enddo
      do k = Kmingh, Kmaxgh
         zf(k) = (k-Kmin)*CellWidth(2,Levelf) + z0
      enddo
      Igrid(gidf) = ((igc) - (0)) * 2 + (0) + lri 
      Jgrid(gidf) = ((jgc) - (0)) * 2 + (0) + lrj
      Kgrid(gidf) = ((kgc) - (0)) * 2 + (0) + lrk
    end subroutine interpU
  end subroutine do_finer
  subroutine do_coarser
    use io_util, only : print_msg
    use string, only : num2char
    use unit 
    integer :: pgid, cgid, crank, n, i, j, k, prank, ngrid, gid
    real(kind=8),dimension(:),pointer :: x, y, z 
    ngrid = count( Coarser /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Coarser) ' // num2char(ngrid) )
    myrank = get_myrank()
    do prank = lbound(Coarser,2), ubound(Coarser,2)
       do n = lbound(Coarser,1), GidListNodeMax(Levelc, prank)
          if ( Coarser(n, prank) == Undefi ) cycle
          pgid = Coarser(n, prank) 
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   cgid = ChildGid(i,j,k,pgid,prank)
                   crank = ChildRank(i,j,k,pgid,prank)
                   if ( myrank == crank ) call dealloc_U( cgid )
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine do_coarser
  subroutine do_fine
    use io_util, only : print_msg
    use string, only : num2char
    use packarr
    integer :: pgid, prank, gids, gidd, ranks, rankd, n, i, j, k, ngrid
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(:),pointer :: x, y, z
    ngrid = count( Fine /= Undefi )
    if ( ngrid == 0 ) return
    call print_msg( 'refine grids (Fine)    ' // num2char( ngrid ) )
    myrank = get_myrank()
    call pkar_reset
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) 
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   gids = ChildGid(i,j,k,pgid,prank)
                   ranks = ChildRank(i,j,k,pgid,prank)
                   rankd = Newrank(i,j,k,n,prank)
                   if ( ranks == rankd ) cycle
                   if ( ranks == myrank ) then 
                      u => get_Up(gids)
                      x => get_Xp(gids)
                      y => get_Yp(gids)
                      z => get_Zp(gids)
                   endif
                   if ((myrank) == (ranks)) call pkar_push(u(lbound(u,1),lbound(u,2),lbound(u,3),lbound(u,4)), size(u), kind(u), ra&
&nkd) 
 if ((myrank) == (rankd)) call pkar_recvlen((Imaxgh-Imingh+1)*(Jmaxgh-Jmingh+1)*(Kmaxgh-Kmingh+1)*(Mmax-Mmin+1), kind(u), ranks)
                   if ((myrank) == (ranks)) call pkar_push(x(lbound(x,1)), size(x), kind(x), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen((Imaxgh-Imingh+1), kind(x), ranks)
                   if ((myrank) == (ranks)) call pkar_push(y(lbound(y,1)), size(y), kind(y), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen((Jmaxgh-Jmingh+1), kind(y), ranks)
                   if ((myrank) == (ranks)) call pkar_push(z(lbound(z,1)), size(z), kind(z), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen((Kmaxgh-Kmingh+1), kind(z), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Igrid(gids), 1, kind(Igrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Igrid(gids)), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Jgrid(gids), 1, kind(Jgrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Jgrid(gids)), ranks)
                   if ((myrank) == (ranks)) call pkar_push(Kgrid(gids), 1, kind(Kgrid(gids)), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(1, kind(Kgrid(gids)), ranks)
                enddo
             enddo
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) 
          do k=Left,Right
             do j=Left,Right
                do i=Left,Right
                   gids = ChildGid(i,j,k,pgid,prank)
                   ranks = ChildRank(i,j,k,pgid,prank)
                   rankd = Newrank(i,j,k,n,prank)
                   if ( ranks == rankd ) cycle
                   if ( rankd == myrank ) then 
                      gidd = alloc_U( Levelf )
                      u => get_Up(gidd)
                      x => get_Xp(gidd)
                      y => get_Yp(gidd)
                      z => get_Zp(gidd)
                      if ((myrank) == (rankd)) call pkar_pop(u(lbound(u,1),lbound(u,2),lbound(u,3),lbound(u,4)), size(u), kind(u), &
&ranks)
                      if ((myrank) == (rankd)) call pkar_pop(x(lbound(x,1)), size(x), kind(x), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(y(lbound(y,1)), size(y), kind(y), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(z(lbound(z,1)), size(z), kind(z), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(Igrid(gidd), 1, kind(Igrid(gidd)), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(Jgrid(gidd), 1, kind(Jgrid(gidd)), ranks)
                      if ((myrank) == (rankd)) call pkar_pop(Kgrid(gidd), 1, kind(Kgrid(gidd)), ranks)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    do prank = lbound(Fine,2), ubound(Fine,2)
       do n = lbound(Fine,1), GidListNodeMax(Levelc, prank)
          if ( Fine(n, prank) == Undefi ) cycle
          pgid = Fine(n, prank) 
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
  subroutine Emulate_2Dim_BoolMap
    integer :: n, gid
    logical,dimension(Gidmin:Gidmax,0:400 -1) :: bufboolmap
    myrank = get_myrank()
    do n = lbound(GidList,1), GidListMax( Levelc )
       gid = GidList(n, Levelc)
       if ( Kgrid(gid) == 0 ) cycle
       BoolMap( gid, myrank ) = .false.
    end do
    call mpi_allreduce_BoolMap(MPI_LAND)
  end subroutine Emulate_2Dim_BoolMap
end module refine
