
#include "config.h"
!-------------------------------------------------------------------------
! MPI library for AMR code
!-------------------------------------------------------------------------
module mpilib
  use mpi
  implicit none
  integer,save :: ierr,myrank,npe
  integer,save :: status(MPI_STATUS_SIZE)
  integer,save :: MPI_LLONG_INTEGER
contains
  !-------------------------------------------------------------------------
  ! initialize and start MPI
  !-------------------------------------------------------------------------
  subroutine mpilib_init
    ! initialize MPI
    call mpi_init(ierr)
    ! define new type, LLONG_KIND
    call mpi_type_contiguous(2, MPI_INTEGER, MPI_LLONG_INTEGER, ierr)
    call mpi_type_commit( MPI_LLONG_INTEGER, ierr )
  end subroutine mpilib_init
  !-------------------------------------------------------------------------
  subroutine mpi_param(npe,myrank)
    integer,intent(OUT) :: npe, myrank
    call mpi_comm_size(MPI_COMM_WORLD, npe, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  end subroutine mpi_param
!-------------------------------------------------------------------------
  function get_myrank() result(myrank)
    integer :: myrank
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  end function get_myrank
!-------------------------------------------------------------------------
  function mpi_get_npe() result(npe)
    integer :: npe
    call mpi_comm_size(MPI_COMM_WORLD, npe, ierr)
  end function mpi_get_npe
!-------------------------------------------------------------------------
! 各プロセスのうちどれかが TRUE なら計算を止める
!-------------------------------------------------------------------------
  subroutine mpi_halt(bool)
    logical,intent(IN) :: bool
    logical :: boolglobal
    call mpi_allreduce( bool, boolglobal, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    if (boolglobal) then
       call mpi_finalize(ierr)
       stop
    endif
  end subroutine mpi_halt


!-------------------------------------------------------------------------
! 以下の試みはうまく行かない。 Fortan90 故に。。。
!-------------------------------------------------------------------------
!   include 'mpif.h'
!   integer :: ierr,myrank,npe
!   integer :: status(MPI_STATUS_SIZE)
!   call mpi_comm_size(MPI_COMM_WORLD, npe, ierr)
!   call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
!-------------------------------------------------------------------------
!!$  function mpilib_myrank() result(myrank)
!!$    include 'mpif.h'
!!$    integer :: ierr, myrank
!!$    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
!!$  end function mpilib_myrank
!!$!-------------------------------------------------------------------------
!!$  subroutine mpilib_allreduce_dbl_sum(arr)
!!$    include 'mpif.h'
!!$    real(kind=8),dimension(:) :: arr
!!$    real(kind=8),dimension(size(arr)) :: dest
!!$    integer :: ierr
!!$    call mpi_allreduce( arr, dest, size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
!!$    arr = dest
!!$  end subroutine mpilib_allreduce_dbl_sum
!!$!-------------------------------------------------------------------------
!!$  subroutine mpilib_allreduce_int_sum(arr)
!!$    include 'mpif.h'
!!$    integer,dimension(:) :: arr
!!$    integer,dimension(size(arr)) :: dest
!!$    integer :: ierr
!!$    call mpi_allreduce( arr, dest, size(arr), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
!!$    arr = dest
!!$  end subroutine mpilib_allreduce_int_sum
end module mpilib

