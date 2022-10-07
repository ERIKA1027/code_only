
! Tomoaki Matsumoto <matsu@i.hosei.ac.jp>
#include "config.h"
! ---------------------------------------------------------------------
! Packing arrays to one-dimensional array for send and recv in MPI.
!     {array1, array2, ....}
!
! USAGE:
!   Define these CPP macros
!   #define PACK_SEND(A, DST) \
!       call pkar_push(A, size(A), kind(A), DST)
!   #define UNPACK_RECV(A, SRC) \
!       call pkar_pop(A, size(A), kind(A), SRC)
!
! EXAMPLE:
!    use packarr
!    #include "packarr.h"
!    integer :: position
!    [anytype] :: A(3,4,5), B(0:100, -4:10, 100), C(3,4,5)
!
!    call pkar_reset()                !! Reset modules before usage.
!    PACK_SEND(A, rankd1)             !! pack A into buffer. The internal position is slieded
!    PACK_SEND(B, rankd1)             !! pack B into buffer. The internal position is slieded
!    PACK_SEND(C, rankd2)             !! pack C into buffer. The internal position is slieded
!    call pkar_sendrecv()             !! transter packed arrays
!    UNPACK_RECV(A, ranks)            !! A is unpacked from buffer in ranks.
!    UNPACK_RECV(B, ranks)            !! B is unpacked from buffer in ranks.
!    UNPACK_RECV(C, ranks)            !! B is unpacked from buffer in ranks.
! ---------------------------------------------------------------------
module packarr
  implicit none
  private
  integer,parameter :: BUFSIZE = 1024  ! initial buffer size of pkg array
  integer,parameter :: KINDINT = kind(1) ! kind of integer
  integer,parameter :: IUNITY = 1        ! Unity
  integer,parameter :: TAG = 1           ! tag used in MPI
  ! buffer list
  type type_buflist
     integer,dimension(:),pointer :: pkg => null() ! buffer of array
     integer :: position = 1  ! position in pkg
  end type type_buflist
  type(type_buflist),save,dimension(0:NPE-1),target :: Bufs ! list of send buffer
  type(type_buflist),save,dimension(0:NPE-1),target :: Bufd ! list of recv buffer

  public :: KINDINT, pkar_reset, pkar_reset_all, pkar_sendrecv, pkar_push_int, pkar_pop_int, &
       pkar_recvlen_int
contains
  !---------------------------------------------------------------------
  ! reset this module
  !---------------------------------------------------------------------
  subroutine pkar_reset
    ! ---------------
    ! positionの初期化
    ! ---------------
    Bufs(:)%position = 1
    Bufd(:)%position = 1
  end subroutine pkar_reset
  !---------------------------------------------------------------------
  ! reset buffer and reset this module
  !---------------------------------------------------------------------
  subroutine pkar_reset_all
    integer :: rank
    Bufs(:)%position = 1
    Bufd(:)%position = 1
    do rank = 0, NPE-1
       if (associated(Bufs(rank)%pkg)) deallocate( Bufs(rank)%pkg )
       nullify( Bufs(rank)%pkg )
       if (associated(Bufd(rank)%pkg)) deallocate( Bufd(rank)%pkg )
       nullify( Bufd(rank)%pkg )
    end do
  end subroutine pkar_reset_all
  !---------------------------------------------------------------------
  ! send recv of buffers
  !---------------------------------------------------------------------
  subroutine pkar_sendrecv()
    use mpilib
    integer,dimension(0:NPE-1) :: reqs, reqd
    integer :: rank

    myrank = get_myrank()
    ! ---------------
    ! 本体を送受信する
    ! ---------------
    do rank = 0, NPE -1
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
    ! -----------------
    ! 同一ノード内はコピー
    ! -----------------
    if (Bufs(myrank)%position-1 /= 0) then
       if (Bufd(myrank)%position /= Bufs(myrank)%position) then
          print*, '*** inconsistent data length', Bufd(myrank)%position, Bufs(myrank)%position, myrank
       endif
       call pkar_test_buf(Bufd(myrank), Bufd(myrank)%position-1, init_position = 1)
       Bufd(myrank)%pkg(1:Bufd(myrank)%position-1) = Bufs(myrank)%pkg(1:Bufs(myrank)%position-1)
    endif
    ! -----------------
    ! 転送待ち
    ! -----------------
    do rank = 0, NPE -1
       if (rank == myrank) cycle
       if (Bufs(rank)%position-1 /= 0) call mpi_wait(reqs(rank),status,ierr)
       if (Bufd(rank)%position-1 /= 0) call mpi_wait(reqd(rank),status,ierr)
    end do
    ! --------------------------------
    ! positionの初期化 (この後のpopのため)
    ! --------------------------------
    Bufs(:)%position = 1
    Bufd(:)%position = 1
  end subroutine pkar_sendrecv
  ! ---------------------------------------------------------------------
  ! truncate buffer within size of position-1.
  ! *buf (inout) ...... buffer
  ! position .......... the last position.
  ! ---------------------------------------------------------------------
  subroutine pkar_truncate_buf(buf, position)
    integer,intent(IN) :: position
    integer,dimension(:),pointer :: buf, bufnew
    integer :: sz
    sz = position -1
    allocate( bufnew(sz) )
!VECT
    bufnew = buf(1:sz)
    deallocate( buf )
    buf => bufnew
  end subroutine pkar_truncate_buf
  ! ---------------------------------------------------------------------
  ! Test buffer and exted buffer size if neccesary.
  !  buft (inout) ....... buffer type
  !  sizea (in) ......... size of an array to store in unit of 4 BYTE (kind of int)
  !  init_position(IN, Optional) ... the initial position in buffer.
  !                       Buffer size is calculated by this argument,
  !                       but current position is not affected.
  ! ---------------------------------------------------------------------
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
    ! initial position
    if ( .not. associated( buft%pkg ) ) then
       allocate( buft%pkg(BUFSIZE) )
       if (.not. present(init_position)) buft%position = 1
    endif
    buf => buft%pkg
    sza = sizea
    szb = size(buf)
    ! lower and upper bound of position
    lb = position
    ub = position + sza - 1
    ! extend buffer size (buf)
    do while (ub > szb)
       szb = szb + BUFSIZE
    end do
    if ( szb >  size(buf) ) then
       allocate( bufnew(szb) )
!VECT
       bufnew(lbound(buf,1):ubound(buf,1)) = buf(:)
       deallocate( buf )
       buft%pkg => bufnew
    end if
  end subroutine pkar_test_buf
  ! ---------------------------------------------------------------------
  ! Push and pack a array of integer
  !  a (in) ............ array of integer type
  !  rankd(in) ......... rank for destination
  ! ---------------------------------------------------------------------
  subroutine pkar_push_int(a, rankd)
    integer,dimension(:),intent(IN) :: a
    integer,intent(IN) :: rankd
    integer,dimension(:),pointer :: buf
    integer,pointer :: position
    integer :: lb, ub, sizea
    ! test for buffer
    sizea = size(a)
    call pkar_test_buf(Bufs(rankd), sizea)
    ! push a to buf
    position => Bufs(rankd)%position
    buf => Bufs(rankd)%pkg
    lb = position
    ub = position+sizea-1
!VECT
    buf(lb:ub) = a(:)
    position = ub+1             ! shift to the next position
  end subroutine pkar_push_int
  ! ---------------------------------------------------------------------
  ! Pop and unpack array of integer
  !  a (out) ........... array of integer type
  !  ranks (in) ........ rank for destination
  ! ---------------------------------------------------------------------
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
!VECT
    a(:) = buf(lb:ub)
    position = ub+1             ! shift to the next position
  end subroutine pkar_pop_int
  ! ---------------------------------------------------------------------
  !  assgin receive length
  !  len (in) ...... length of array
  !  ranks (in) ....  rank of source
  ! ---------------------------------------------------------------------
  subroutine pkar_recvlen_int(len, ranks)
    integer,intent(IN) :: len, ranks
    Bufd(ranks)%position = Bufd(ranks)%position + len
  end subroutine pkar_recvlen_int
end module packarr
! ---------------------------------------------------------------------
! Belows are external subroutines for use of a consistent array
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! Pack and push an array of any type
! call pkar_push(a, sizea, buf, sizebuf, position, kinda)
!  a (in) ............ inputed 1-dimensional array of any type
!  sizea (in) ........ size of array in unit of element number
!  kinda (in) ........ kind of array a (size of an elements in unit of byte)
!  rankd (in) ........ rank for destination
! ---------------------------------------------------------------------
subroutine pkar_push(a, sizea, kinda, rankd)
  use packarr
  implicit none
  integer,intent(IN) :: sizea, kinda, rankd
  integer,intent(IN)  :: a(sizea*kinda/KINDINT)
  call pkar_push_int(a, rankd)
end subroutine pkar_push
! ---------------------------------------------------------------------
! unpack and pop array
! call pkar_pop(a, sizea, buf, sizebuf, position, kinda)
!  a (out) ........... outputed 1-dimensional array of any type
!  sizea (in) ........ size of array a (number of elements)
!  kinda (in) ........ kind of array a
!  ranks (in) ........ rank for source
! ---------------------------------------------------------------------
subroutine pkar_pop(a, sizea, kinda, ranks)
  use packarr
  implicit none
  integer,intent(IN) :: kinda, sizea, ranks
  integer,intent(OUT) :: a(sizea*kinda/KINDINT)
  call pkar_pop_int(a, ranks)
end subroutine pkar_pop
! ---------------------------------------------------------------------
! assign receive length
! len ..... length of array
! kinda ... kind of array
! ranks ... rank of source
! ---------------------------------------------------------------------
subroutine pkar_recvlen(len, kinda, ranks)
  use packarr
  implicit none
  integer,intent(IN) :: len, kinda, ranks
  call pkar_recvlen_int(len*kinda/KINDINT, ranks)
end subroutine pkar_recvlen
