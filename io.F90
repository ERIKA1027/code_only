#include "config.h"
#define TMP "tmp"
!
! For checking consistency of linklist when read data.
! #define CHECK_LINKLIST
!
! Use MPI_Barrier when read data from every node.
! #define IO_BLOCKING
!
!-------------------------------------------------------------------------
! Module for Input and Output
!-------------------------------------------------------------------------
module io
  use string, only : CHARLEN
  implicit none
  private
  integer,parameter :: IO_VERSION = 5
  integer,save :: Used_grid_num
#ifdef IO_BLOCKING
  logical,parameter :: BLOCKING_IO = .TRUE.
#else !IO_BLOCKING
  logical,parameter :: BLOCKING_IO = .FALSE.
#endif !IO_BLOCKING
  public :: dumpdata, restoredata, dumpslice, get_filename
contains
  !---------------------------------------------------------------------
  ! dump data PE順に書く
  !---------------------------------------------------------------------
  subroutine dumpdata
    use systemcall
    use mpilib
    use grid
    use io_util, only : print_msg
    use string, only : CHARLEN
    integer :: n
    character(len=CHARLEN) :: tmpfile
    logical :: exist
#ifdef FORBIT_WRITEDATA
    return
#endif

    call print_msg("write data")

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if ( get_myrank() == PRIMARY_RANK ) then
       tmpfile = get_tmpfile()
       inquire(file=tmpfile, exist=exist)
       if (exist) call systemcall_unlink(tmpfile)
    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call print_msg("(IO) KS DEBUG 1") ! KS DEBUG

    Used_grid_num = grid_used_num_global()
    do n = 0, mpi_get_npe() - 1
       if (BLOCKING_IO) call mpi_barrier(MPI_COMM_WORLD, ierr)
       if ( n == get_myrank() ) then
          call dumpdata_npe
       end if
    end do

    call print_msg("(IO) KS DEBUG 2") ! KS DEBUG

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if ( get_myrank() == PRIMARY_RANK ) then
       tmpfile = get_tmpfile()
       call print_msg("(IO) KS DEBUG "//"cat "// trim(tmpfile) // "* >" // trim(tmpfile)) ! KS DEBUG
       if (.not. BLOCKING_IO) call systemcall_command("cat "// trim(tmpfile) // "* >" // trim(tmpfile)) !gather tmpfile
       call print_msg("(IO) KS DEBUG "//"sh " // trim(tmpfile)) ! KS DEBUG
       call systemcall_command("sh " // trim(tmpfile))
    endif
    call print_msg("(IO) KS DEBUG 3") ! KS DEBUG

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (.not. BLOCKING_IO) call systemcall_unlink(get_tmpfile_npe())
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call print_msg("(IO) KS DEBUG 4") ! KS DEBUG

  end subroutine dumpdata
  !---------------------------------------------------------------------
  ! restore data PE順に読む
  !---------------------------------------------------------------------
  subroutine restoredata
    use mpilib
    use grid
    integer :: n, l
    GidListMax = Undefi
    GidList = Undefi
    NeighborGid = Undefi
    NeighborRank = MPI_PROC_NULL
    ChildGid = Undefi
    ChildRank = MPI_PROC_NULL
    ParentGid = Undefi
    ParentRank = MPI_PROC_NULL
    Levels = Undefi
    do n = 0, mpi_get_npe() - 1
       if (BLOCKING_IO) call mpi_barrier(MPI_COMM_WORLD, ierr)
       if ( n == get_myrank() ) then
          call restoredata_npe
       end if
    end do
    call restoreDatabase
    do l = Lmin, LevelMax
       call update_gidlist(l)
    end do
#ifdef CHECK_LINKLIST
    do l = Lmin, LevelMax
       call check_linklist(l)
    enddo
#endif
  end subroutine restoredata
  !---------------------------------------------------------------------
  ! restore database
  !---------------------------------------------------------------------
  subroutine restoreDatabase
    use grid
    use mpilib
    integer :: n, i, j, k, gid
    integer,dimension(ARRAYSIZE3(NeighborGid)) :: ngid
    integer,dimension(ARRAYSIZE3(NeighborRank)):: nrank
    integer,dimension(ARRAYSIZE4(ChildGid))    :: cgid
    integer,dimension(ARRAYSIZE4(ChildRank))   :: crank
    integer,dimension(ARRAYSIZE1(ParentGid))   :: pgid
    integer,dimension(ARRAYSIZE1(ParentRank))  :: prank
    integer,dimension(ARRAYSIZE3(GidBase))     :: gbase
    integer,dimension(ARRAYSIZE3(RankBase))    :: rbase

    myrank = get_myrank()

    ! restore GidBase and RankBase
    gbase(:,:,:) = Undefi
    rbase(:,:,:) = MPI_PROC_NULL
    do n = Gidmin, GidListMax(Lmin)
       gid = GidList(n, Lmin)
       i = Igrid(gid)
       j = Jgrid(gid)
       k = Kgrid(gid)
       gbase(i,j,k) = gid
       rbase(i,j,k) = myrank
    enddo
    call mpi_allreduce(gbase, GidBase, size(gbase), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
    call mpi_allreduce(rbase, RankBase, size(rbase), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )

    ! restore neighbor lists
    ngid = NeighborGid(:,:,:,myrank)
    call mpi_allgather(ngid , size(ngid), MPI_INTEGER, NeighborGid, size(ngid), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    nrank = NeighborRank(:,:,:,myrank)
    call mpi_allgather(nrank , size(nrank), MPI_INTEGER, NeighborRank, size(nrank), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    cgid = ChildGid(:,:,:,:,myrank)
    call mpi_allgather(cgid , size(cgid), MPI_INTEGER, ChildGid, size(cgid), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    crank = ChildRank(:,:,:,:,myrank)
    call mpi_allgather(crank , size(crank), MPI_INTEGER, ChildRank, size(crank), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    pgid = ParentGid(:,myrank)
    call mpi_allgather(pgid , size(pgid), MPI_INTEGER, ParentGid, size(pgid), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    prank = ParentRank(:,myrank)
    call mpi_allgather(prank , size(prank), MPI_INTEGER, ParentRank, size(prank), MPI_INTEGER, MPI_COMM_WORLD, ierr )

  end subroutine restoreDatabase
  !---------------------------------------------------------------------
  ! dump data for each PE
  !---------------------------------------------------------------------
#define FH 11
#define DUMP_PARAM(RW,FH) \
  RW(FH) Imin, Jmin, Kmin, Mmin, Gidmin, Lmin ;\
  RW(FH) Imax, Jmax, Kmax, Mmax, Gidmax, Lmax ;\
  RW(FH) Ngh ;\
  RW(FH) Imingh, Jmingh, Kmingh ;\
  RW(FH) Imaxgh, Jmaxgh, Kmaxgh

#define DUMP_HEADER(RW,FH) \
  RW(FH) LevelMax ;\
  RW(FH) GidListMax(Lmin:LevelMax) ;\
  RW(FH) GidList(Gidmin: maxval(GidListMax(Lmin:LevelMax)),Lmin:LevelMax) ;\
  RW(FH) Time(Lmin:LevelMax) ;\
  RW(FH) Dtime(Lmin:LevelMax) ;\
  RW(FH) Step(Lmin:LevelMax) ;\
  RW(FH) Dstep(Lmin:LevelMax) ;\
  RW(FH) CellWidth(:,Lmin:LevelMax)

#define DUMP_DATA(RW,FH) \
  RW(FH) gid ;\
  RW(FH) Levels(gid) ;\
  RW(FH) U1list(gid) ;\
  RW(FH) U2list(gid) ;\
  RW(FH) AllocMem(gid), AllocMem(U1list(gid)), AllocMem(U2list(gid)) ;\
  if (.not. associated(BlockMem(GID)%u)) allocate(BlockMem(GID)%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) ;\
  if (.not. associated(BlockMem(GID)%x)) allocate(BlockMem(GID)%x(Imingh:Imaxgh)) ;\
  if (.not. associated(BlockMem(GID)%y)) allocate(BlockMem(GID)%y(Jmingh:Jmaxgh)) ;\
  if (.not. associated(BlockMem(GID)%z)) allocate(BlockMem(GID)%z(Kmingh:Kmaxgh)) ;\
  if (.not. associated(BlockMem(U1list(gid))%u)) allocate(BlockMem(U1list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) ;\
  if (.not. associated(BlockMem(U2list(gid))%u)) allocate(BlockMem(U2list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) ;\
  RW(FH) BlockMem(gid)%u ;\
  RW(FH) BlockMem(U1list(gid))%u ;\
  RW(FH) BlockMem(U2list(gid))%u ;\
  RW(FH) BlockMem(gid)%x ;\
  RW(FH) BlockMem(gid)%y ;\
  RW(FH) BlockMem(gid)%z ;\
  RW(FH) Igrid(gid), Jgrid(gid), Kgrid(gid) ;\
  RW(FH) NeighborGid(:,:,gid,myrank), NeighborRank(:,:,gid,myrank) ;\
  RW(FH) ChildGid(:,:,:,gid,myrank), ChildRank(:,:,:,gid,myrank) ;\
  RW(FH) ParentGid(gid,myrank), ParentRank(gid, myrank)

#define DUMP_DATA_ALLOC_BLOCK(RW,GID) \
  if (RW == 'read')


  subroutine dumpdata_npe()
    use mpilib
    use grid
    use reflux, only : reflux_write
    use io_util, only : read_env, wchar
    character(len=CHARLEN) :: fn, dir, prefix, suffix, fns, fnt
    integer :: gid
    call read_env('DIR', dir)
    call read_env('PREFIX', prefix)
    call read_env('SUFFIX', suffix)
    fn = get_filename(dir, prefix, suffix) ! filename for IO
    open(FH, file=fn, form='unformatted')
    write(FH) IO_VERSION
    write(FH) Used_grid_num
    write(FH) grid_used_num()
    DUMP_PARAM(write,FH)
    DUMP_HEADER(write,FH)
    do gid = lbound(Levels,1), ubound(Levels,1)
       if ( Levels(gid) == Undefi ) cycle
       DUMP_DATA(write,FH)
       call reflux_write(FH, gid)
    enddo
    call flush(FH)
    close(FH)
    ! --------------------------
    ! make tmpfile for post dump
    ! --------------------------
    if (BLOCKING_IO) then
       fnt = get_tmpfile()
       open(FH, file=fnt, position='APPEND')
    else
       fnt = get_tmpfile_npe()
       open(FH, file=fnt)
    endif
    fns = get_dumpfile()
    call wchar(FH, "/bin/rm -f "//trim(fns))
    call wchar(FH, "ln -s "//trim(fn)//" "//trim(fns))
    call flush(FH)
    close(FH)
  end subroutine dumpdata_npe
  !---------------------------------------------------------------------
  ! restore data
  !---------------------------------------------------------------------
  subroutine restoredata_npe
    use mpilib
    use grid
    use reflux, only : reflux_read
    use string
    use io_util, only : wchar
    character(len=CHARLEN) :: fn
    integer :: io_v, i, j, k, m, gid, l, ngrid, eof

    fn = get_dumpfile()

    call wchar(6,'dumpfile = '//fn)

    open(FH, file=fn, form='unformatted')
    read(FH) io_v
    if (io_v /= IO_VERSION) then
       write(*,*)  '**** restoredata: IO_VERSION is not consistent'
       stop
    endif

    read(FH) ngrid
    read(FH) ngrid

    ! read min values
    read(FH) i, j, k, m, gid, l
    if ( i /= Imin .or. j /= Jmin .or. k /= Kmin .or. m /= Mmin ) then
       write(*,*)  '**** restoredata: IJKM MIN is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    if ( l < Lmin .or.  gid < Gidmin ) then
       write(*,*)  '**** restoredata: LMIN or NGID is not consistent'
       call mpi_finalize(ierr)
       stop
    endif

    ! read max values
    read(FH) i, j, k, m, gid, l
    if ( i /= Imax .or. j /= Jmax .or. k /= Kmax .or. m /= Mmax ) then
       write(*,*)  '**** restoredata: IJKM MAX is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
!!$    if ( l > Lmax .or.  gid > Gidmax ) then
!!$       write(*,*)  '**** restoredata: LMAX or NGID is not consistent'
!!$       call mpi_finalize(ierr)
!!$       stop
!!$    endif

    read(FH) i
    if ( i /= Ngh ) then
       write(*,*)  '**** restoredata: N_GHOST_CELL is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    !  this is dummy io
    read(FH) i,j,k
    read(FH) i,j,k


    DUMP_HEADER(read,FH)
#define FHEOF FH,iostat=eof
    do
       DUMP_DATA(read,FHEOF)
       call reflux_read(FH, gid, eof)

       if ( gid > Gidmax ) then
          print *, '*** error gid > Gidmax', gid, Gidmax
          stop
       endif
       if ( Levels(gid) > Lmax ) then
          print *, '*** error level > Lmax', gid, Gidmax
          stop
       endif

       if (eof /= 0) exit
    enddo
    call flush(FH)
    close(FH)

  end subroutine restoredata_npe
  !---------------------------------------------------------------------
  ! get file name
  !---------------------------------------------------------------------
  subroutine filename(fn,step,prefix,suffix)
    use string
    integer(kind=LLONG_KIND),intent(in) :: step
    character(len=*),intent(out) :: fn
    character(len=*),intent(in) :: prefix, suffix
    fn = num2char(step)
    fn = concat(prefix,fn)
    fn = concat(fn,suffix)
  end subroutine filename
  ! --------------------------------------------------------------------
  ! write slice data
  ! --------------------------------------------------------------------
  subroutine dumpslice(ncrd, xyz)
    use mpilib
    use grid
    use string
    use io_util, only : print_msg
    integer,intent(IN) :: ncrd
    real(kind=DBL_KIND),intent(IN) :: xyz
    integer,dimension(:),pointer :: iolist
    integer :: iolistmax
    integer :: n, ngridl, ngridg
#ifdef FORBIT_WRITEDATA
    return
#endif
    call print_msg("write slice data")
    ! 出力リストを作る
    call get_gid_list(ncrd, xyz, iolist, iolistmax)
    ngridl = iolistmax + 1
    call mpi_allreduce(ngridl, ngridg, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

    do n = 0, mpi_get_npe() - 1
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       if ( n == get_myrank() ) then
          call dumpslice_npe
       end if
    end do
    deallocate(iolist)
  contains
    ! ---------------------------------------------------
    ! iolist に登録されたデータを出力する
    ! ---------------------------------------------------
    subroutine dumpslice_npe
      use io_util, only : read_env
      character(len=CHARLEN) :: fn, dir, prefix, suffix
      character(len=1) :: dot = '.'
      integer :: gid, n
      call read_env('DIR', dir)
      call read_env('PREFIX', prefix)
      call read_env('SUFFIX', suffix)
      select case(ncrd)
      case (MX)
         prefix = concat(prefix, 'x')
      case (MY)
         prefix = concat(prefix, 'y')
      case (MZ)
         prefix = concat(prefix, 'z')
      end select
      fn = get_filename(dir, prefix, suffix)
      open(FH, file=fn, form='unformatted')
      write(FH) IO_VERSION
      write(FH) ngridg
      write(FH) ngridl
      DUMP_PARAM(write,FH)
      DUMP_HEADER(write,FH)
      do n = Gidmin, iolistmax
         gid = iolist(n)
         DUMP_DATA(write,FH)
      enddo
      call flush(FH)
      close(FH)
    end subroutine dumpslice_npe
    ! ---------------------------------------------------
    ! 出力するべきグリッドのリスト(iolist)を求める
    ! ---------------------------------------------------
    subroutine get_gid_list(ncrd, xyz, iolist, iolistmax)
      integer,intent(IN) :: ncrd
      real(kind=DBL_KIND),intent(IN) :: xyz
      integer,intent(OUT) :: iolistmax
      integer,dimension(:),pointer :: iolist
      integer :: gid
      allocate(iolist(Gidmin:grid_used_num()))
      iolistmax = Gidmin - 1
      do gid = lbound(Levels,1), ubound(Levels,1)
         if ( Levels(gid) == Undefi ) cycle
         if ( .not. bool_include(gid, ncrd, xyz)) cycle
         iolistmax = iolistmax + 1
         iolist(iolistmax) = gid
      end do
    end subroutine get_gid_list
    ! ---------------------------------------------------
    ! グリッド gid が xyz(ncrd) を含むかどうか？
    ! ---------------------------------------------------
    function bool_include(gid, ncrd, xyz) result(bool)
      integer,intent(IN) :: gid, ncrd
      real(kind=DBL_KIND),intent(IN) :: xyz
      logical :: bool
      real(kind=DBL_KIND),dimension(:),pointer :: coord
      real(kind=DBL_KIND) :: xyzmin, xyzmax
      integer :: lev
      bool = .false.
      lev = Levels(gid)
      if ( lev == Undefi ) return
      select case (ncrd)
      case (MX)
         coord => get_Xp(gid)
         xyzmin = coord(imin) - CellWidth(MX, lev)
         xyzmax = coord(imax) + CellWidth(MX, lev)
      case (MY)
         coord => get_Yp(gid)
         xyzmin = coord(jmin) - CellWidth(MY, lev)
         xyzmax = coord(jmax) + CellWidth(MY, lev)
      case (MZ)
         coord => get_Zp(gid)
         xyzmin = coord(kmin) - CellWidth(MZ, lev)
         xyzmax = coord(kmax) + CellWidth(MZ, lev)
      end select
      if ( xyz >= xyzmin .and. xyz <= xyzmax ) bool = .true.
    end function bool_include
  end subroutine dumpslice
  !-------------------------------------------------------------------------
  ! make file name for IO
  !-------------------------------------------------------------------------
  function get_filename(dir, prefix, suffix) result(fn)
    use string
    use grid, only : Step, LevelMax
    use mpilib
    character(len=*),intent(IN) :: dir, prefix, suffix
    character(len=CHARLEN) :: fn
    character(len=1) :: dot = '.'
    myrank = get_myrank()
    fn = dir                               ! /dir/
    fn = concat(fn,prefix)                 ! /dir/st
    fn = concat(fn,num2char(Step(LevelMax))) ! /dir/st12000
    fn = concat(fn,dot)                    ! /dir/st12000.
    fn = concat(fn,num2char(myrank))         ! /dir/st12000.0
    fn = concat(fn,dot)                    ! /dir/st12000.0.
    fn = concat(fn,suffix)                 ! /dir/st12000.0.d
  end function get_filename
  !-------------------------------------------------------------------------
  ! make file name for dumpfile
  !-------------------------------------------------------------------------
  function get_dumpfile() result(fn)
    use io_util, only : read_env
    use string
    use grid
    use mpilib
    character(len=CHARLEN) :: fn, dir, suffix, dfile
    character(len=1) :: dot = '.'
    myrank = get_myrank()
    call read_env('DIR', dir)
    call read_env('DUMPFILE', dfile)
    call read_env('SUFFIX', suffix)
    fn = dir                               ! /dir/
    fn = concat(fn,dfile)                  ! /dir/dump
    fn = concat(fn,dot)                    ! /dir/dump.
    fn = concat(fn,num2char(myrank))         ! /dir/dump.0
    fn = concat(fn,dot)                    ! /dir/dump.0.
    fn = concat(fn,suffix)                 ! /dir/st12000.0.d
  end function get_dumpfile
  !-------------------------------------------------------------------------
  ! make tmp file name for dump data
  !-------------------------------------------------------------------------
  function get_tmpfile() result(fn)
    use io_util, only : read_env
    use string
    use grid
    character(len=CHARLEN) :: fn, dir
    call read_env('DIR', dir)
    fn = concat(dir,TMP)
  end function get_tmpfile
  !-------------------------------------------------------------------------
  ! make tmp file name for each node
  !-------------------------------------------------------------------------
  function get_tmpfile_npe() result(fn)
    use string
    use mpilib
    character(len=CHARLEN) :: fn
    myrank = get_myrank()
    fn = concat(get_tmpfile(), num2char(myrank))
  end function get_tmpfile_npe
  !-------------------------------------------------------------------------
  ! get number of lines for given file
  !-------------------------------------------------------------------------
  function io_get_nline(fn) result(nline)
    use io_util, only : read_env
    use systemcall
    use string
    use mpilib
    character(len=CHARLEN) :: fn, dir, tmpfile
    integer :: nline
    integer,parameter :: UNIT=1
    logical :: exist
    if ( get_myrank() == PRIMARY_RANK ) then
       call read_env('DIR', dir)
       tmpfile = trim(dir) // 'temp_systemout'
       inquire(file=tmpfile, exist=exist)
       if (exist) call systemcall_unlink(tmpfile)

       call systemcall_command("out=(`wc -l " // trim(fn) // "`) && echo ${out[0]}>"// tmpfile  )
       open(UNIT, file=tmpfile)
       read(UNIT, *) nline
       call flush(UNIT)
       close(UNIT, status='DELETE')
    end if
    call mpi_bcast(nline, 1, MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
  end function io_get_nline
end module io
