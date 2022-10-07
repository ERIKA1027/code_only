module mpilib
  use mpi
  implicit none
  integer,save :: ierr,myrank,npe
  integer,save :: status(MPI_STATUS_SIZE)
  integer,save :: MPI_LLONG_INTEGER
contains
  subroutine mpilib_init
    call mpi_init(ierr)
    call mpi_type_contiguous(2, MPI_INTEGER, MPI_LLONG_INTEGER, ierr)
    call mpi_type_commit( MPI_LLONG_INTEGER, ierr )
  end subroutine mpilib_init
  subroutine mpi_param(npe,myrank)
    integer,intent(OUT) :: npe, myrank
    call mpi_comm_size(MPI_COMM_WORLD, npe, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  end subroutine mpi_param
  function get_myrank() result(myrank)
    integer :: myrank
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  end function get_myrank
  function mpi_get_npe() result(npe)
    integer :: npe
    call mpi_comm_size(MPI_COMM_WORLD, npe, ierr)
  end function mpi_get_npe
  subroutine mpi_halt(bool)
    logical,intent(IN) :: bool
    logical :: boolglobal
    call mpi_allreduce( bool, boolglobal, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    if (boolglobal) then
       call mpi_finalize(ierr)
       stop
    endif
  end subroutine mpi_halt
end module mpilib
