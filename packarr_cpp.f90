module packarr
  implicit none
  private
  integer,parameter :: BUFSIZE = 1024 
  integer,parameter :: KINDINT = kind(1) 
  integer,parameter :: IUNITY = 1 
  integer,parameter :: TAG = 1 
  type type_buflist
     integer,dimension(:),pointer :: pkg => null() 
     integer :: position = 1 
  end type type_buflist
  type(type_buflist),save,dimension(0:400 -1),target :: Bufs 
  type(type_buflist),save,dimension(0:400 -1),target :: Bufd 
  public :: KINDINT, pkar_reset, pkar_reset_all, pkar_sendrecv, pkar_push_int, pkar_pop_int, &
       pkar_recvlen_int
contains
  subroutine pkar_reset
    Bufs(:)%position = 1
    Bufd(:)%position = 1
  end subroutine pkar_reset
  subroutine pkar_reset_all
    integer :: rank
    Bufs(:)%position = 1
    Bufd(:)%position = 1
    do rank = 0, 400 -1
       if (associated(Bufs(rank)%pkg)) deallocate( Bufs(rank)%pkg )
       nullify( Bufs(rank)%pkg )
       if (associated(Bufd(rank)%pkg)) deallocate( Bufd(rank)%pkg )
       nullify( Bufd(rank)%pkg )
    end do
  end subroutine pkar_reset_all
  subroutine pkar_sendrecv()
    use mpilib
    integer,dimension(0:400 -1) :: reqs, reqd
    integer :: rank
    myrank = get_myrank()
    do rank = 0, 400 -1
       if (rank == myrank) cycle
       if (Bufs(rank)%position-1 /= 0) then
          call mpi_isend(Bufs(rank)%pkg(1), Bufs(rank)%position-1, MPI_INTEGER, &
               rank, TAG, MPI_COMM_WORLD, reqs(rank), ierr )
       endif
       if (Bufd(rank)%position-1 /= 0) then
          call pkar_test_buf(Bufd(rank), Bufd(rank)%position-1, init_position = 1)
          call mpi_irecv(Bufd(rank)%pkg(1), Bufd(rank)%position-1, MPI_INTEGER, &
               rank, TAG, MPI_COMM_WORLD, reqd(rank), ierr )
       endif
    enddo
    if (Bufs(myrank)%position-1 /= 0) then
       if (Bufd(myrank)%position /= Bufs(myrank)%position) then
          print*, '*** inconsistent data length', Bufd(myrank)%position, Bufs(myrank)%position, myrank
       endif
       call pkar_test_buf(Bufd(myrank), Bufd(myrank)%position-1, init_position = 1)
       Bufd(myrank)%pkg(1:Bufd(myrank)%position-1) = Bufs(myrank)%pkg(1:Bufs(myrank)%position-1)
    endif
    do rank = 0, 400 -1
       if (rank == myrank) cycle
       if (Bufs(rank)%position-1 /= 0) call mpi_wait(reqs(rank),status,ierr)
       if (Bufd(rank)%position-1 /= 0) call mpi_wait(reqd(rank),status,ierr)
    end do
    Bufs(:)%position = 1
    Bufd(:)%position = 1
  end subroutine pkar_sendrecv
  subroutine pkar_truncate_buf(buf, position)
    integer,intent(IN) :: position
    integer,dimension(:),pointer :: buf, bufnew
    integer :: sz
    sz = position -1
    allocate( bufnew(sz) )
    bufnew = buf(1:sz)
    deallocate( buf )
    buf => bufnew
  end subroutine pkar_truncate_buf
  subroutine pkar_test_buf(buft, sizea, init_position)
    type(type_buflist),intent(INOUT) :: buft
    integer,intent(IN) :: sizea
    integer,optional :: init_position
    integer,dimension(:),pointer :: buf, bufnew
    integer :: sza, szb, lb, ub, position
    if (present(init_position)) then
       position = init_position
    else
       position = buft%position
    endif
    if ( .not. associated( buft%pkg ) ) then
       allocate( buft%pkg(BUFSIZE) )
       if (.not. present(init_position)) buft%position = 1
    endif
    buf => buft%pkg
    sza = sizea
    szb = size(buf)
    lb = position
    ub = position + sza - 1
    do while (ub > szb)
       szb = szb + BUFSIZE
    end do
    if ( szb > size(buf) ) then
       allocate( bufnew(szb) )
       bufnew(lbound(buf,1):ubound(buf,1)) = buf(:)
       deallocate( buf )
       buft%pkg => bufnew
    end if
  end subroutine pkar_test_buf
  subroutine pkar_push_int(a, rankd)
    integer,dimension(:),intent(IN) :: a
    integer,intent(IN) :: rankd
    integer,dimension(:),pointer :: buf
    integer,pointer :: position
    integer :: lb, ub, sizea
    sizea = size(a)
    call pkar_test_buf(Bufs(rankd), sizea)
    position => Bufs(rankd)%position
    buf => Bufs(rankd)%pkg
    lb = position
    ub = position+sizea-1
    buf(lb:ub) = a(:)
    position = ub+1 
  end subroutine pkar_push_int
  subroutine pkar_pop_int(a, ranks)
    integer,dimension(:),intent(OUT) :: a
    integer,intent(IN) :: ranks
    integer,dimension(:),pointer :: buf
    integer,pointer :: position
    integer :: lb, ub
    position => Bufd(ranks)%position
    buf => Bufd(ranks)%pkg
    lb = position
    ub = position+size(a)-1
    if (ub > ubound(buf,1)) print *, 'error in pkar_pop_int', ub, ubound(buf,1)
    a(:) = buf(lb:ub)
    position = ub+1 
  end subroutine pkar_pop_int
  subroutine pkar_recvlen_int(len, ranks)
    integer,intent(IN) :: len, ranks
    Bufd(ranks)%position = Bufd(ranks)%position + len
  end subroutine pkar_recvlen_int
end module packarr
subroutine pkar_push(a, sizea, kinda, rankd)
  use packarr
  implicit none
  integer,intent(IN) :: sizea, kinda, rankd
  integer,intent(IN) :: a(sizea*kinda/KINDINT)
  call pkar_push_int(a, rankd)
end subroutine pkar_push
subroutine pkar_pop(a, sizea, kinda, ranks)
  use packarr
  implicit none
  integer,intent(IN) :: kinda, sizea, ranks
  integer,intent(OUT) :: a(sizea*kinda/KINDINT)
  call pkar_pop_int(a, ranks)
end subroutine pkar_pop
subroutine pkar_recvlen(len, kinda, ranks)
  use packarr
  implicit none
  integer,intent(IN) :: len, kinda, ranks
  call pkar_recvlen_int(len*kinda/KINDINT, ranks)
end subroutine pkar_recvlen
