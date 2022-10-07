#include "config.h"
!#include "config-gTV.h"
!#include "recl.h"
#define SZG 0:NGI-1,0:NGJ-1,0:NGK-1,MX:MZ
! #define USE_DIRECT_ACCESS
module initReadVelocity
  implicit none
  private
  ! mesh numbers in each direction
  integer,parameter :: NGI=(NI)*(NGI_BASE), NGJ=(NJ)*(NGJ_BASE), NGK=(NK)*(NGK_BASE)
  real(KIND=DBL_KIND),save,dimension(:,:,:,:),allocatable :: Vg
  public :: readVelocity, extractVelocity
contains
  !-------------------------------------------------------------------------
  ! read velocity field for all the node
  !-------------------------------------------------------------------------
  subroutine readVelocity(machNumber)
    use mpilib
    use io_util, only : wchar, readenv
    use string, only :  CHARLEN, concat
    real(KIND=DBL_KIND),intent(IN) :: machNumber
    character(len=CHARLEN) :: fn, dir
    !integer(KIND=LLONG_KIND) :: recl
    integer,parameter :: UNIT=11
    if ( get_myrank() /= PRIMARY_RANK ) return
    allocate( Vg(SZG) )
    !recl = int(size(Vg), LLONG_KIND)*int(DBL_KIND, LLONG_KIND)/int(RECL_UNIT, LLONG_KIND)
    if (.not. readenv('TURBVEL', fn)) then
       print *, '**** error in readVelocity: TURBVEL is not specified'
    end if
    if (.not. readenv('DIR', dir)) then
       print *, '**** error in readVelocity: DIR is not specified'
    end if
    fn = concat(dir,fn)
    call wchar(6,'read velocity file = '//fn)
!#ifdef USE_DIRECT_ACCESS
!!!$       print *, 'recl = ',recl
!    open(unit=UNIT, file=fn, form='unformatted', access='direct', recl=recl )
!    read(UNIT, rec=1) Vg
!    call flush(UNIT)
!    close(UNIT)
!#else !USE_DIRECT_ACCESS
    open(unit=UNIT, file=fn, form='unformatted', access='stream')
    read(UNIT) Vg
    call flush(UNIT)
    close(UNIT)
!#endif !USE_DIRECT_ACCESS
    Vg = Vg * machNumber
    call wchar(6,'... done')
  end subroutine readVelocity
  !-------------------------------------------------------------------------
  ! extract velocity
  !-------------------------------------------------------------------------
  subroutine extractVelocity(v, comp, gid, rank)
    use mpilib
    use grid, only : Igrid, Jgrid, Kgrid, Imin, Jmin, Kmin, Imax, Jmax, Kmax
    real(KIND=DBL_KIND),dimension(:,:,:),pointer :: v
    integer,intent(IN) :: comp, gid, rank
    integer :: is, js, ks, ie, je, ke, ijkbuf(MX:MZ)
    real(KIND=DBL_KIND),dimension(:,:,:),allocatable :: buf
    integer :: reqs, reqd
    integer,parameter :: TAG = 1           ! tag used in MPI

    myrank = get_myrank()
    if (myrank /= PRIMARY_RANK .and. myrank /= rank) then
       print *, '**** error in extractVelocity'
    end if


    if (myrank == rank) then
       is = Igrid(gid) * (NI)
       js = Jgrid(gid) * (NJ)
       ks = Kgrid(gid) * (NK)
       ijkbuf = (/is,js,ks/)
       call mpi_isend(ijkbuf, size(ijkbuf), MPI_INTEGER, PRIMARY_RANK, TAG, MPI_COMM_WORLD, reqs, ierr)
    end if
    if (myrank == PRIMARY_RANK ) then
       call mpi_irecv(ijkbuf, size(ijkbuf), MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, reqd, ierr)
    endif
    if (myrank == rank)         call mpi_wait(reqs,status,ierr)
    if (myrank == PRIMARY_RANK) call mpi_wait(reqd,status,ierr)

    if (myrank == PRIMARY_RANK ) then
       is = ijkbuf(MX)
       js = ijkbuf(MY)
       ks = ijkbuf(MZ)
       ie = is + (NI)-1
       je = js + (NJ)-1
       ke = ks + (NK)-1
    end if

    allocate(buf(NI, NJ, NK))
    if ( myrank == PRIMARY_RANK ) then
       buf = Vg(is:ie,js:je,ks:ke,comp)
       call mpi_isend(buf, size(buf), MPI_DOUBLE_PRECISION, rank, TAG, MPI_COMM_WORLD, reqs, ierr)
    end if
    if ( myrank == rank ) then
       call mpi_irecv(buf, size(buf), MPI_DOUBLE_PRECISION, PRIMARY_RANK, TAG, MPI_COMM_WORLD, reqd, ierr)
    end if
    if ( myrank == PRIMARY_RANK ) call mpi_wait(reqs,status,ierr)
    if ( myrank == rank )         call mpi_wait(reqd,status,ierr)

    if ( myrank == rank ) v(Imin:Imax,Jmin:Jmax,Kmin:Kmax) = buf

    deallocate(buf)

  end subroutine extractVelocity
end module initReadVelocity
