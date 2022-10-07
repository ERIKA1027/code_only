#include "config.h"
#include "packarr.h"
#include "overBlockCoordinates.h"
!-------------------------------------------------------------------------
! Modole: uniformpatch
! 
! Gather data from AMR-blocks to a uniform patch.
! Scatter data from a uniform patch to AMR-blocks.
! AMR-blocks ara in the same grid level.  
!-------------------------------------------------------------------------
module patchBlock
  use overBlockCoordinates
  use grid
  implicit none
  private
  type t_upatch
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
     real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
     integer,dimension(:),pointer :: comp
     real(kind=DBL_KIND),dimension(MX:MZ) :: cellwidth
     type(t_obRect) :: rect
  end type t_upatch
  public :: t_upatch, pb_gather, pb_scatter, pb_gather_byRectPhys, pb_gather_byCoordPhys, pb_isRectCoveredByBlocks
contains
  !-------------------------------------------------------------------------
  ! Gather data from AMR blocks to a uniform patch of rankd
  ! 
  ! INPUT:
  !   rect .... rectangle for patch region (type t_obRect)
  !   mlist ... list of components that the pactch stores
  !   rankd ... MPI rank where patch is created
  !   bool_checkRect ... check whether the rect is valid or not (MPI com intensive)
  ! OUTPUT:
  !   patch ... patch data, which is allocated in this subroutine.
  !             It should be deallocated somewhere.
  !-------------------------------------------------------------------------
  subroutine pb_gather(rect, patch, mlist, rankd, bool_checkRect)
    use mpilib
    use packarr
    type(t_obRect),intent(IN) :: rect
    type(t_upatch),pointer :: patch ! (OUT)
    integer,dimension(0:),intent(IN) :: mlist
    integer,intent(IN) :: rankd
    logical,intent(IN),optional :: bool_checkRect
    integer :: level, gid, ranks, ig, jg, kg, is, js, ks, ie, je, ke
    type(t_obPoint) :: posL, posR
    integer,dimension(MX:MZ) :: ijkgL, ijkgR, ijkgrid
    type(t_obRect) :: rectBlock, rectand
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: ubuf
    real(kind=DBL_KIND),dimension(:),allocatable :: xbuf, ybuf, zbuf
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    logical :: bool, bool_chk
    bool_chk = .FALSE.
    if (present(bool_checkRect)) bool_chk = bool_checkRect
    if (bool_chk) bool = pb_isRectCoveredByBlocks(rect)

    myrank = get_myrank()
    nullify(patch)                                     ! initialize
    call ob_extractLevelFromRect(level, rect) ! get level
    call ob_extractPointPairFromRect(posL, posR, rect) ! get posL, posR
    call ob_getIjkgridFromPoint(ijkgL, posL)           ! get ijkgL
    call ob_getIjkgridFromPoint(ijkgR, posR)           ! get ijkgR

    if (rankd == myrank) then
       call pb_alloc(rect, patch, size(mlist))
       patch%comp = mlist
       patch%cellwidth = CellWidth(:,level)
    endif
    call pkar_reset
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,ranks) !get gid, ranks
             if (ranks == myrank .or. rankd == myrank) then
                ijkgrid = (/ig, jg, kg/)
                call ob_getRectFromIjkgrid(rectBlock, ijkgrid, level)
                call ob_rectAnd( rect, rectBlock, rectand )
                is =  mod( rectand%left(MX),  Imax-Imin+1 ) ! index in block
                js =  mod( rectand%left(MY),  Jmax-Jmin+1 )
                ks =  mod( rectand%left(MZ),  Kmax-Kmin+1 )
                ie =  mod( rectand%right(MX), Imax-Imin+1 )
                je =  mod( rectand%right(MY), Jmax-Jmin+1 )
                ke =  mod( rectand%right(MZ), Kmax-Kmin+1 )
             endif
             if (ranks == myrank) then
                allocate(ubuf(is:ie,js:je,ks:ke,size(mlist)), xbuf(is:ie), ybuf(js:je), zbuf(ks:ke))
                u => get_Up(gid)
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                ubuf = u(is:ie,js:je,ks:ke,mlist)
                xbuf = x(is:ie)
                ybuf = y(js:je)
                zbuf = z(ks:ke)
             endif
             if (myrank == rankd .or. myrank == ranks) then
                PACK_SEND4_SZ(ubuf, myrank, ranks, rankd, (ie-is+1)*(je-js+1)*(ke-ks+1)*size(mlist))
                PACK_SEND1_SZ(xbuf, myrank, ranks, rankd, (ie-is+1))
                PACK_SEND1_SZ(ybuf, myrank, ranks, rankd, (je-js+1))
                PACK_SEND1_SZ(zbuf, myrank, ranks, rankd, (ke-ks+1))
             endif
             if (ranks == myrank) then
                deallocate(ubuf, xbuf, ybuf, zbuf)
             end if
          end do
       end do
    end do
    call pkar_sendrecv()
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,ranks) !get gid, ranks
             if (rankd == myrank) then
                ijkgrid = (/ig, jg, kg/)
                call ob_getRectFromIjkgrid(rectBlock, ijkgrid, level)
                call ob_rectAnd( rect, rectBlock, rectand )
                is = rectand%left(MX)  - posL%p(MX) ! index in patch
                js = rectand%left(MY)  - posL%p(MY)
                ks = rectand%left(MZ)  - posL%p(MZ)
                ie = rectand%right(MX) - posL%p(MX)
                je = rectand%right(MY) - posL%p(MY)
                ke = rectand%right(MZ) - posL%p(MZ)
                allocate(ubuf(is:ie,js:je,ks:ke,size(mlist)), xbuf(is:ie), ybuf(js:je), zbuf(ks:ke))
                UNPACK_RECV4(ubuf, myrank, ranks, rankd)
                UNPACK_RECV1(xbuf, myrank, ranks, rankd)
                UNPACK_RECV1(ybuf, myrank, ranks, rankd)
                UNPACK_RECV1(zbuf, myrank, ranks, rankd)
                patch%u(is:ie,js:je,ks:ke,:) = ubuf
                patch%x(is:ie) = xbuf
                patch%y(js:je) = ybuf
                patch%z(ks:ke) = zbuf
                deallocate(ubuf, xbuf, ybuf, zbuf)
             end if
          end do
       end do
    end do
  end subroutine pb_gather
  !-------------------------------------------------------------------------
  ! Wrapper of pb_gather
  !-------------------------------------------------------------------------
  subroutine pb_gather_byRectPhys(rectPhys, patch, mlist, rankd, level, bool_checkRect)
    type(t_obRectPhys),intent(IN) :: rectPhys
    type(t_upatch),pointer :: patch ! (OUT)
    integer,dimension(0:),intent(IN) :: mlist
    integer,intent(IN) :: rankd, level
    logical,intent(IN),optional :: bool_checkRect
    type(t_obRect) :: rect
    logical :: bool_chk
    call ob_RectPhys2RectOb( rectPhys, level, rect)
    bool_chk = .FALSE.
    if (present(bool_checkRect)) bool_chk = bool_checkRect
    call pb_gather(rect, patch, mlist, rankd, bool_chk)
  end subroutine pb_gather_byRectPhys
  !-------------------------------------------------------------------------
  ! Wrapper of pb_gather
  !-------------------------------------------------------------------------
  subroutine pb_gather_byCoordPhys(coordPhys, patch, mlist, rankd, level, bool_checkRect)
    real(kind=DBL_KIND),dimension(OB_COORDS_SZ),intent(IN) :: coordPhys
    type(t_upatch),pointer :: patch ! (OUT)
    integer,dimension(0:),intent(IN) :: mlist
    integer,intent(IN) :: rankd, level
    logical,intent(IN),optional :: bool_checkRect
    type(t_obRect) :: rect
    type(t_obRectPhys) :: rectPhys
    logical :: bool_chk
    call ob_assignCoordPhysToRectPhys(coordPhys, rectPhys)
    call ob_RectPhys2RectOb( rectPhys, level, rect)
    bool_chk = .FALSE.
    if (present(bool_checkRect)) bool_chk = bool_checkRect
    call pb_gather(rect, patch, mlist, rankd, bool_chk)
  end subroutine pb_gather_byCoordPhys
  !-------------------------------------------------------------------------
  ! Scatter data from a uniform patch in rank to AMR glocks
  !
  ! INPUT:
  !   patch ... patch data, which is allocated in this subroutine.
  !             It is deallocated in this subroutine.
  !   ranks ... MPI rank where patch exists
  ! OUTPUT:
  !   none
  !-------------------------------------------------------------------------
  subroutine pb_scatter(patch, ranks)
    use mpilib
    use packarr
    type(t_upatch),pointer :: patch ! (INOUT)
    integer,intent(IN) :: ranks
    integer :: level, gid, rankd, ig, jg, kg, is, js, ks, ie, je, ke, msize
    type(t_obPoint) :: posL, posR
    integer,dimension(MX:MZ) :: ijkgL, ijkgR, ijkgrid
    type(t_obRect) :: rectBlock, rectand, rect
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: ubuf
    real(kind=DBL_KIND),dimension(:),allocatable :: xbuf, ybuf, zbuf
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    integer,dimension(:),allocatable :: mlist
    myrank = get_myrank()
    if (ranks == myrank) then
       rect = patch%rect
       msize = size(patch%comp)
    end if
    call mpi_bcast(rect, 1, MPI_OB_RECT, ranks, MPI_COMM_WORLD, ierr)
    call mpi_bcast(msize, 1, MPI_INTEGER, ranks, MPI_COMM_WORLD, ierr)
    allocate(mlist(0:msize-1))
    if (ranks == myrank) then
       mlist = patch%comp
    end if
    call mpi_bcast(mlist, msize, MPI_INTEGER, ranks, MPI_COMM_WORLD, ierr)

    call ob_extractLevelFromRect(level, rect) ! get level
    call ob_extractPointPairFromRect(posL, posR, rect) ! get posL, posR
    call ob_getIjkgridFromPoint(ijkgL, posL)           ! get ijkgL
    call ob_getIjkgridFromPoint(ijkgR, posR)           ! get ijkgR

    call pkar_reset
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rankd) ! get gid, rankdd
             if (ranks == myrank .or. rankd == myrank) then
                ijkgrid = (/ig, jg, kg/)
                call ob_getRectFromIjkgrid(rectBlock, ijkgrid, level)
                call ob_rectAnd( rect, rectBlock, rectand )
                is = rectand%left(MX)  - posL%p(MX) ! index in patch
                js = rectand%left(MY)  - posL%p(MY)
                ks = rectand%left(MZ)  - posL%p(MZ)
                ie = rectand%right(MX) - posL%p(MX)
                je = rectand%right(MY) - posL%p(MY)
                ke = rectand%right(MZ) - posL%p(MZ)
             end if
             if (ranks == myrank) then
                allocate(ubuf(is:ie,js:je,ks:ke,size(mlist)), xbuf(is:ie), ybuf(js:je), zbuf(ks:ke))
                ubuf = patch%u(is:ie,js:je,ks:ke,:)
                xbuf = patch%x(is:ie)
                ybuf = patch%y(js:je)
                zbuf = patch%z(ks:ke)
             end if
             if (myrank == rankd .or. myrank == ranks) then             
                PACK_SEND4_SZ(ubuf, myrank, ranks, rankd, (ie-is+1)*(je-js+1)*(ke-ks+1)*size(mlist))
                PACK_SEND1_SZ(xbuf, myrank, ranks, rankd, (ie-is+1))
                PACK_SEND1_SZ(ybuf, myrank, ranks, rankd, (je-js+1))
                PACK_SEND1_SZ(zbuf, myrank, ranks, rankd, (ke-ks+1))
             endif
             if (ranks == myrank) then
                deallocate(ubuf, xbuf, ybuf, zbuf)
             end if
          end do
       end do
    end do
    call pkar_sendrecv()
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rankd) ! get gid, rankd
             if (rankd == myrank) then
                ijkgrid = (/ig, jg, kg/)
                call ob_getRectFromIjkgrid(rectBlock, ijkgrid, level)
                call ob_rectAnd( rect, rectBlock, rectand )
                is =  mod( rectand%left(MX),  Imax-Imin+1 ) ! index in block
                js =  mod( rectand%left(MY),  Jmax-Jmin+1 )
                ks =  mod( rectand%left(MZ),  Kmax-Kmin+1 )
                ie =  mod( rectand%right(MX), Imax-Imin+1 )
                je =  mod( rectand%right(MY), Jmax-Jmin+1 )
                ke =  mod( rectand%right(MZ), Kmax-Kmin+1 )
                allocate(ubuf(is:ie,js:je,ks:ke,size(mlist)), xbuf(is:ie), ybuf(js:je), zbuf(ks:ke))
                UNPACK_RECV4(ubuf, myrank, ranks, rankd)
                UNPACK_RECV1(xbuf, myrank, ranks, rankd)
                UNPACK_RECV1(ybuf, myrank, ranks, rankd)
                UNPACK_RECV1(zbuf, myrank, ranks, rankd)
                u => get_Up(gid)
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                u(is:ie,js:je,ks:ke,mlist) = ubuf 
                x(is:ie) = xbuf
                y(js:je) = ybuf
                z(ks:ke) = zbuf
                deallocate(ubuf, xbuf, ybuf, zbuf)
             end if
          end do
       end do
    end do
    deallocate(mlist)
    if (ranks == myrank) call pb_dealloc(patch)
  end subroutine pb_scatter
  !-------------------------------------------------------------------------
  ! check whether rect is covered by blocks or not
  !-------------------------------------------------------------------------
  function pb_isRectCoveredByBlocks(rect) result(bool)
    use mpilib
    use packarr
    type(t_obRect),intent(IN) :: rect
    logical :: bool             ! (OUT)
    type(t_obPoint) :: posL, posR
    integer,dimension(MX:MZ) :: ijkgL, ijkgR
    type(t_obRect) :: rectComp
    integer :: level, gid, ranks, ig, jg, kg
    logical :: bool_exist
    bool = .true.
    myrank = get_myrank()
    call ob_extractLevelFromRect(level, rect) ! get level
    call ob_extractPointPairFromRect(posL, posR, rect) ! get posL, posR
    call ob_getIjkgridFromPoint(ijkgL, posL)           ! get ijkgL
    call ob_getIjkgridFromPoint(ijkgR, posR)           ! get ijkgR
    ! check rect is within computational box
    call ob_computationBoxOfRect( rectComp, level )
    if (ob_rectsComp( rect,rectComp) /= OB_RECT_INNER .and. myrank == PRIMARY_RANK) then
       print *, '*** pb_isRectCoveredByBlocks: invalid rect', rect
    endif
    do kg = ijkgL(MZ), ijkgR(MZ)
       do jg = ijkgL(MY), ijkgR(MY)
          do ig = ijkgL(MX), ijkgR(MX)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,ranks) !get gid, ranks
             if (gid /= Undefi .and. ranks /= MPI_PROC_NULL) then
                bool_exist = .TRUE.
             else
                bool_exist = .FALSE.
             end if
             call mpi_allreduce(MPI_IN_PLACE, bool_exist, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
             if (.not. bool_exist .and. myrank == PRIMARY_RANK) then
                print *, '*** pb_isRectCoveredByBlocks: rect is not coveted by any blocks', rect, ig, jg, kg
             end if
             bool = bool .and. bool_exist
          end do
       end do
    end do
  end function pb_isRectCoveredByBlocks
  !-------------------------------------------------------------------------
  ! allocate patch
  !-------------------------------------------------------------------------
  subroutine pb_alloc(rect, patch, msize)
    type(t_obRect),intent(IN) :: rect ! rectangle
    type(t_upatch),pointer :: patch ! (INOUT)
    integer,intent(IN) :: msize     ! number of components
    integer,dimension(MX:MZ) :: sz

    allocate(patch)
    sz(:) = rect%right(:) - rect%left(:) + 1 ! cast from long int to int
    allocate(patch%u(0:sz(MX)-1, 0:sz(MY)-1, 0:sz(MZ)-1, 0:msize-1))
    allocate(patch%x(0:sz(MX)-1), patch%y(0:sz(MY)-1), patch%z(0:sz(MZ)-1))
    allocate(patch%comp(0:msize-1))
    patch%u = 0.d0
    patch%x = 0.d0
    patch%y = 0.d0
    patch%z = 0.d0
    patch%comp = -1
    patch%rect = rect
  end subroutine pb_alloc
  !-------------------------------------------------------------------------
  ! deallocate patch
  !-------------------------------------------------------------------------
#define DEALLOC(P)  deallocate(P); nullify(P)
  subroutine pb_dealloc(patch)
    type(t_upatch),pointer :: patch ! (INOUT)
    DEALLOC(patch%u)
    DEALLOC(patch%x)
    DEALLOC(patch%y)
    DEALLOC(patch%z)
    DEALLOC(patch%comp)
    DEALLOC(patch)
  end subroutine pb_dealloc

end module patchBlock
