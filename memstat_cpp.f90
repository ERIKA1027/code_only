module memstat
  implicit none
  private
  public :: memstat_print
contains
  subroutine memstat_print
    use mpilib
    use io_util
    use string
    integer(kind=8),dimension(0:400 -1) :: hwmList, rssList
    integer(kind=8) :: hwm, rss, hwmmax, rssmax, hwmmin, rssmin
    integer :: hwmmaxrank, rssmaxrank, hwmminrank, rssminrank
    integer,dimension(1) :: loc
    call memstat_get_hwm_rss(hwm, rss)
    call mpi_gather( &
         hwm, 1, MPI_LLONG_INTEGER, &
         hwmList, 1, MPI_LLONG_INTEGER, &
         0, MPI_COMM_WORLD, ierr)
    call mpi_gather( &
         rss, 1, MPI_LLONG_INTEGER, &
         rssList, 1, MPI_LLONG_INTEGER, &
         0, MPI_COMM_WORLD, ierr)
    loc = maxloc(hwmList)
    hwmmaxrank = loc(1) - 1 + lbound(hwmList,1)
    hwmmax = hwmList(hwmmaxrank)
    loc = maxloc(rssList)
    rssmaxrank = loc(1) - 1 + lbound(rssList,1)
    rssmax = rssList(rssmaxrank)
    loc = minloc(hwmList)
    hwmminrank = loc(1) - 1 + lbound(hwmList,1)
    hwmmin = hwmList(hwmminrank)
    loc = minloc(rssList)
    rssminrank = loc(1) - 1 + lbound(rssList,1)
    rssmin = rssList(rssminrank)
    call print_msg( 'VmHWM = '// &
         trim(num2char(hwmmax)) // ' kB (rank = '// trim(num2char(hwmmaxrank)) // ') - ' // &
         trim(num2char(hwmmin)) // ' kB (rank = '// trim(num2char(hwmminrank)) // ')' )
    call print_msg( 'VmRSS = '// &
         trim(num2char(rssmax)) // ' kB (rank = '// trim(num2char(rssmaxrank)) // ') - ' // &
         trim(num2char(rssmin)) // ' kB (rank = '// trim(num2char(rssminrank)) // ')' )
  end subroutine memstat_print
  subroutine memstat_get_hwm_rss(hwm, rss)
    integer(kind=8),intent(OUT) :: hwm, rss
    integer,parameter :: LUN=92
    character(len=80) :: line
    open(LUN,file='/proc/self/status')
    hwm=0
    rss=0
    do while(.true.)
       read(LUN,'(a)',end=99) line
       if (line(1:6).eq.'VmHWM:') read(line(8:80),*) hwm
       if (line(1:6).eq.'VmRSS:') read(line(8:80),*) rss
    enddo
99 close(LUN)
  end subroutine memstat_get_hwm_rss
end module memstat
