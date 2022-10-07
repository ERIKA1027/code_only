module io
  use string, only : CHARLEN
  implicit none
  private
  integer,parameter :: IO_VERSION = 5
  integer,save :: Used_grid_num
  logical,parameter :: BLOCKING_IO = .FALSE.
  public :: dumpdata, restoredata, dumpslice, get_filename
contains
  subroutine dumpdata
    use systemcall
    use mpilib
    use grid
    use io_util, only : print_msg
    use string, only : CHARLEN
    integer :: n
    character(len=CHARLEN) :: tmpfile
    logical :: exist
    call print_msg("write data")
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if ( get_myrank() == 0 ) then
       tmpfile = get_tmpfile()
       inquire(file=tmpfile, exist=exist)
       if (exist) call systemcall_unlink(tmpfile)
    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call print_msg("(IO) KS DEBUG 1") 
    Used_grid_num = grid_used_num_global()
    do n = 0, mpi_get_npe() - 1
       if (BLOCKING_IO) call mpi_barrier(MPI_COMM_WORLD, ierr)
       if ( n == get_myrank() ) then
          call dumpdata_npe
       end if
    end do
    call print_msg("(IO) KS DEBUG 2") 
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if ( get_myrank() == 0 ) then
       tmpfile = get_tmpfile()
       call print_msg("(IO) KS DEBUG "//"cat "// trim(tmpfile) // "* >" // trim(tmpfile)) 
       if (.not. BLOCKING_IO) call systemcall_command("cat "// trim(tmpfile) // "* >" // trim(tmpfile)) 
       call print_msg("(IO) KS DEBUG "//"sh " // trim(tmpfile)) 
       call systemcall_command("sh " // trim(tmpfile))
    endif
    call print_msg("(IO) KS DEBUG 3") 
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (.not. BLOCKING_IO) call systemcall_unlink(get_tmpfile_npe())
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call print_msg("(IO) KS DEBUG 4") 
  end subroutine dumpdata
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
  end subroutine restoredata
  subroutine restoreDatabase
    use grid
    use mpilib
    integer :: n, i, j, k, gid
    integer,dimension(lbound(NeighborGid,1):ubound(NeighborGid,1),lbound(NeighborGid,2):ubound(NeighborGid,2),lbound(NeighborGid,3)&
&:ubound(NeighborGid,3)) :: ngid
    integer,dimension(lbound(NeighborRank,1):ubound(NeighborRank,1),lbound(NeighborRank,2):ubound(NeighborRank,2),lbound(NeighborRa&
&nk,3):ubound(NeighborRank,3)):: nrank
    integer,dimension(lbound(ChildGid,1):ubound(ChildGid,1),lbound(ChildGid,2):ubound(ChildGid,2),lbound(ChildGid,3):ubound(ChildGi&
&d,3),lbound(ChildGid,4):ubound(ChildGid,4)) :: cgid
    integer,dimension(lbound(ChildRank,1):ubound(ChildRank,1),lbound(ChildRank,2):ubound(ChildRank,2),lbound(ChildRank,3):ubound(Ch&
&ildRank,3),lbound(ChildRank,4):ubound(ChildRank,4)) :: crank
    integer,dimension(lbound(ParentGid,1):ubound(ParentGid,1)) :: pgid
    integer,dimension(lbound(ParentRank,1):ubound(ParentRank,1)) :: prank
    integer,dimension(lbound(GidBase,1):ubound(GidBase,1),lbound(GidBase,2):ubound(GidBase,2),lbound(GidBase,3):ubound(GidBase,3)) &
&:: gbase
    integer,dimension(lbound(RankBase,1):ubound(RankBase,1),lbound(RankBase,2):ubound(RankBase,2),lbound(RankBase,3):ubound(RankBas&
&e,3)) :: rbase
    myrank = get_myrank()
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
    fn = get_filename(dir, prefix, suffix) 
    open(11, file=fn, form='unformatted')
    write(11) IO_VERSION
    write(11) Used_grid_num
    write(11) grid_used_num()
    write(11) Imin, Jmin, Kmin, Mmin, Gidmin, Lmin 
 write(11) Imax, Jmax, Kmax, Mmax, Gidmax, Lmax 
 write(11) Ngh 
 write(11) Imingh, Jmingh, Kmingh 
 write(11) Imaxgh, Jmaxgh, Kmaxgh
    write(11) LevelMax 
 write(11) GidListMax(Lmin:LevelMax) 
 write(11) GidList(Gidmin: maxval(GidListMax(Lmin:LevelMax)),Lmin:LevelMax) 
 write(11) Time(Lmin:LevelMax) 
 write(11) Dtime(Lmin:LevelMax) 
 write(11) Step(Lmin:LevelMax) 
 write(11) Dstep(Lmin:LevelMax) 
 write(11) CellWidth(:,Lmin:LevelMax)
    do gid = lbound(Levels,1), ubound(Levels,1)
       if ( Levels(gid) == Undefi ) cycle
       write(11) gid 
 write(11) Levels(gid) 
 write(11) U1list(gid) 
 write(11) U2list(gid) 
 write(11) AllocMem(gid), AllocMem(U1list(gid)), AllocMem(U2list(gid)) 
 if (.not. associated(BlockMem(GID)%u)) allocate(BlockMem(GID)%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) 
 if (.not. associated(BlockMem(GID)%x)) allocate(BlockMem(GID)%x(Imingh:Imaxgh)) 
 if (.not. associated(BlockMem(GID)%y)) allocate(BlockMem(GID)%y(Jmingh:Jmaxgh)) 
 if (.not. associated(BlockMem(GID)%z)) allocate(BlockMem(GID)%z(Kmingh:Kmaxgh)) 
 if (.not. associated(BlockMem(U1list(gid))%u)) allocate(BlockMem(U1list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 if (.not. associated(BlockMem(U2list(gid))%u)) allocate(BlockMem(U2list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 write(11) BlockMem(gid)%u 
 write(11) BlockMem(U1list(gid))%u 
 write(11) BlockMem(U2list(gid))%u 
 write(11) BlockMem(gid)%x 
 write(11) BlockMem(gid)%y 
 write(11) BlockMem(gid)%z 
 write(11) Igrid(gid), Jgrid(gid), Kgrid(gid) 
 write(11) NeighborGid(:,:,gid,myrank), NeighborRank(:,:,gid,myrank) 
 write(11) ChildGid(:,:,:,gid,myrank), ChildRank(:,:,:,gid,myrank) 
 write(11) ParentGid(gid,myrank), ParentRank(gid, myrank)
       call reflux_write(11, gid)
    enddo
    call flush(11)
    close(11)
    if (BLOCKING_IO) then
       fnt = get_tmpfile()
       open(11, file=fnt, position='APPEND')
    else
       fnt = get_tmpfile_npe()
       open(11, file=fnt)
    endif
    fns = get_dumpfile()
    call wchar(11, "/bin/rm -f "//trim(fns))
    call wchar(11, "ln -s "//trim(fn)//" "//trim(fns))
    call flush(11)
    close(11)
  end subroutine dumpdata_npe
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
    open(11, file=fn, form='unformatted')
    read(11) io_v
    if (io_v /= IO_VERSION) then
       write(*,*) '**** restoredata: IO_VERSION is not consistent'
       stop
    endif
    read(11) ngrid
    read(11) ngrid
    read(11) i, j, k, m, gid, l
    if ( i /= Imin .or. j /= Jmin .or. k /= Kmin .or. m /= Mmin ) then
       write(*,*) '**** restoredata: IJKM MIN is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    if ( l < Lmin .or. gid < Gidmin ) then
       write(*,*) '**** restoredata: LMIN or NGID is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    read(11) i, j, k, m, gid, l
    if ( i /= Imax .or. j /= Jmax .or. k /= Kmax .or. m /= Mmax ) then
       write(*,*) '**** restoredata: IJKM MAX is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    read(11) i
    if ( i /= Ngh ) then
       write(*,*) '**** restoredata: N_GHOST_CELL is not consistent'
       call mpi_finalize(ierr)
       stop
    endif
    read(11) i,j,k
    read(11) i,j,k
    read(11) LevelMax 
 read(11) GidListMax(Lmin:LevelMax) 
 read(11) GidList(Gidmin: maxval(GidListMax(Lmin:LevelMax)),Lmin:LevelMax) 
 read(11) Time(Lmin:LevelMax) 
 read(11) Dtime(Lmin:LevelMax) 
 read(11) Step(Lmin:LevelMax) 
 read(11) Dstep(Lmin:LevelMax) 
 read(11) CellWidth(:,Lmin:LevelMax)
    do
       read(11,iostat=eof) gid 
 read(11,iostat=eof) Levels(gid) 
 read(11,iostat=eof) U1list(gid) 
 read(11,iostat=eof) U2list(gid) 
 read(11,iostat=eof) AllocMem(gid), AllocMem(U1list(gid)), AllocMem(U2list(gid)) 
 if (.not. associated(BlockMem(GID)%u)) allocate(BlockMem(GID)%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) 
 if (.not. associated(BlockMem(GID)%x)) allocate(BlockMem(GID)%x(Imingh:Imaxgh)) 
 if (.not. associated(BlockMem(GID)%y)) allocate(BlockMem(GID)%y(Jmingh:Jmaxgh)) 
 if (.not. associated(BlockMem(GID)%z)) allocate(BlockMem(GID)%z(Kmingh:Kmaxgh)) 
 if (.not. associated(BlockMem(U1list(gid))%u)) allocate(BlockMem(U1list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 if (.not. associated(BlockMem(U2list(gid))%u)) allocate(BlockMem(U2list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 read(11,iostat=eof) BlockMem(gid)%u 
 read(11,iostat=eof) BlockMem(U1list(gid))%u 
 read(11,iostat=eof) BlockMem(U2list(gid))%u 
 read(11,iostat=eof) BlockMem(gid)%x 
 read(11,iostat=eof) BlockMem(gid)%y 
 read(11,iostat=eof) BlockMem(gid)%z 
 read(11,iostat=eof) Igrid(gid), Jgrid(gid), Kgrid(gid) 
 read(11,iostat=eof) NeighborGid(:,:,gid,myrank), NeighborRank(:,:,gid,myrank) 
 read(11,iostat=eof) ChildGid(:,:,:,gid,myrank), ChildRank(:,:,:,gid,myrank) 
 read(11,iostat=eof) ParentGid(gid,myrank), ParentRank(gid, myrank)
       call reflux_read(11, gid, eof)
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
    call flush(11)
    close(11)
  end subroutine restoredata_npe
  subroutine filename(fn,step,prefix,suffix)
    use string
    integer(kind=8),intent(in) :: step
    character(len=*),intent(out) :: fn
    character(len=*),intent(in) :: prefix, suffix
    fn = num2char(step)
    fn = concat(prefix,fn)
    fn = concat(fn,suffix)
  end subroutine filename
  subroutine dumpslice(ncrd, xyz)
    use mpilib
    use grid
    use string
    use io_util, only : print_msg
    integer,intent(IN) :: ncrd
    real(kind=8),intent(IN) :: xyz
    integer,dimension(:),pointer :: iolist
    integer :: iolistmax
    integer :: n, ngridl, ngridg
    call print_msg("write slice data")
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
    subroutine dumpslice_npe
      use io_util, only : read_env
      character(len=CHARLEN) :: fn, dir, prefix, suffix
      character(len=1) :: dot = '.'
      integer :: gid, n
      call read_env('DIR', dir)
      call read_env('PREFIX', prefix)
      call read_env('SUFFIX', suffix)
      select case(ncrd)
      case (0)
         prefix = concat(prefix, 'x')
      case (1)
         prefix = concat(prefix, 'y')
      case (2)
         prefix = concat(prefix, 'z')
      end select
      fn = get_filename(dir, prefix, suffix)
      open(11, file=fn, form='unformatted')
      write(11) IO_VERSION
      write(11) ngridg
      write(11) ngridl
      write(11) Imin, Jmin, Kmin, Mmin, Gidmin, Lmin 
 write(11) Imax, Jmax, Kmax, Mmax, Gidmax, Lmax 
 write(11) Ngh 
 write(11) Imingh, Jmingh, Kmingh 
 write(11) Imaxgh, Jmaxgh, Kmaxgh
      write(11) LevelMax 
 write(11) GidListMax(Lmin:LevelMax) 
 write(11) GidList(Gidmin: maxval(GidListMax(Lmin:LevelMax)),Lmin:LevelMax) 
 write(11) Time(Lmin:LevelMax) 
 write(11) Dtime(Lmin:LevelMax) 
 write(11) Step(Lmin:LevelMax) 
 write(11) Dstep(Lmin:LevelMax) 
 write(11) CellWidth(:,Lmin:LevelMax)
      do n = Gidmin, iolistmax
         gid = iolist(n)
         write(11) gid 
 write(11) Levels(gid) 
 write(11) U1list(gid) 
 write(11) U2list(gid) 
 write(11) AllocMem(gid), AllocMem(U1list(gid)), AllocMem(U2list(gid)) 
 if (.not. associated(BlockMem(GID)%u)) allocate(BlockMem(GID)%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax)) 
 if (.not. associated(BlockMem(GID)%x)) allocate(BlockMem(GID)%x(Imingh:Imaxgh)) 
 if (.not. associated(BlockMem(GID)%y)) allocate(BlockMem(GID)%y(Jmingh:Jmaxgh)) 
 if (.not. associated(BlockMem(GID)%z)) allocate(BlockMem(GID)%z(Kmingh:Kmaxgh)) 
 if (.not. associated(BlockMem(U1list(gid))%u)) allocate(BlockMem(U1list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 if (.not. associated(BlockMem(U2list(gid))%u)) allocate(BlockMem(U2list(gid))%u(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mma&
&x)) 
 write(11) BlockMem(gid)%u 
 write(11) BlockMem(U1list(gid))%u 
 write(11) BlockMem(U2list(gid))%u 
 write(11) BlockMem(gid)%x 
 write(11) BlockMem(gid)%y 
 write(11) BlockMem(gid)%z 
 write(11) Igrid(gid), Jgrid(gid), Kgrid(gid) 
 write(11) NeighborGid(:,:,gid,myrank), NeighborRank(:,:,gid,myrank) 
 write(11) ChildGid(:,:,:,gid,myrank), ChildRank(:,:,:,gid,myrank) 
 write(11) ParentGid(gid,myrank), ParentRank(gid, myrank)
      enddo
      call flush(11)
      close(11)
    end subroutine dumpslice_npe
    subroutine get_gid_list(ncrd, xyz, iolist, iolistmax)
      integer,intent(IN) :: ncrd
      real(kind=8),intent(IN) :: xyz
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
    function bool_include(gid, ncrd, xyz) result(bool)
      integer,intent(IN) :: gid, ncrd
      real(kind=8),intent(IN) :: xyz
      logical :: bool
      real(kind=8),dimension(:),pointer :: coord
      real(kind=8) :: xyzmin, xyzmax
      integer :: lev
      bool = .false.
      lev = Levels(gid)
      if ( lev == Undefi ) return
      select case (ncrd)
      case (0)
         coord => get_Xp(gid)
         xyzmin = coord(imin) - CellWidth(0, lev)
         xyzmax = coord(imax) + CellWidth(0, lev)
      case (1)
         coord => get_Yp(gid)
         xyzmin = coord(jmin) - CellWidth(1, lev)
         xyzmax = coord(jmax) + CellWidth(1, lev)
      case (2)
         coord => get_Zp(gid)
         xyzmin = coord(kmin) - CellWidth(2, lev)
         xyzmax = coord(kmax) + CellWidth(2, lev)
      end select
      if ( xyz >= xyzmin .and. xyz <= xyzmax ) bool = .true.
    end function bool_include
  end subroutine dumpslice
  function get_filename(dir, prefix, suffix) result(fn)
    use string
    use grid, only : Step, LevelMax
    use mpilib
    character(len=*),intent(IN) :: dir, prefix, suffix
    character(len=CHARLEN) :: fn
    character(len=1) :: dot = '.'
    myrank = get_myrank()
    fn = dir 
    fn = concat(fn,prefix) 
    fn = concat(fn,num2char(Step(LevelMax))) 
    fn = concat(fn,dot) 
    fn = concat(fn,num2char(myrank)) 
    fn = concat(fn,dot) 
    fn = concat(fn,suffix) 
  end function get_filename
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
    fn = dir 
    fn = concat(fn,dfile) 
    fn = concat(fn,dot) 
    fn = concat(fn,num2char(myrank)) 
    fn = concat(fn,dot) 
    fn = concat(fn,suffix) 
  end function get_dumpfile
  function get_tmpfile() result(fn)
    use io_util, only : read_env
    use string
    use grid
    character(len=CHARLEN) :: fn, dir
    call read_env('DIR', dir)
    fn = concat(dir,"tmp")
  end function get_tmpfile
  function get_tmpfile_npe() result(fn)
    use string
    use mpilib
    character(len=CHARLEN) :: fn
    myrank = get_myrank()
    fn = concat(get_tmpfile(), num2char(myrank))
  end function get_tmpfile_npe
  function io_get_nline(fn) result(nline)
    use io_util, only : read_env
    use systemcall
    use string
    use mpilib
    character(len=CHARLEN) :: fn, dir, tmpfile
    integer :: nline
    integer,parameter :: UNIT=1
    logical :: exist
    if ( get_myrank() == 0 ) then
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
    call mpi_bcast(nline, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end function io_get_nline
end module io
