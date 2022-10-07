#include "config.h"
#include "overBlockCoordinates.h"
!-------------------------------------------------------------------------
! Module for interpolation using overBlockCoordinates
!-------------------------------------------------------------------------
module ob_interp
  use overBlockCoordinates
  implicit none
  private
  public :: ob_interpolatedU, ob_interpolatedUbyCoordPhys
contains
  !-------------------------------------------------------------------------
  ! an wrapper of ob_interpolatedU.
  ! Input argument is a physical coordinates vector instead of physical point.
  !-------------------------------------------------------------------------
  subroutine ob_interpolatedUbyCoordPhys(coords, mlist, uave, levelp, gidp, rankp)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: coords
    integer,intent(IN),dimension(:) :: mlist
    real(kind=DBL_KIND),dimension(size(mlist)),intent(OUT) :: uave
    integer,intent(OUT),optional ::  levelp, gidp, rankp
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coords, pointPhys)
    call ob_interpolatedU(pointPhys, mlist, uave, levelp, gidp, rankp)
  end subroutine ob_interpolatedUbyCoordPhys
  !-------------------------------------------------------------------------
  ! get interpolated primitive variables at a given point
  ! INPUT:
  !   pointPhys = a point in physical coordinates.
  !   mlist     = list of components. i.e., mlist = (/MVX, MVY, MVZ/)
  ! OUTPUT:
  !   uave = averaged
  !   rankp = rank defined pointPhys (optional)
  !   gidp  = gid defined pointPhys  (optional)
  !   levelp = level defined pointPhys (optional)
  !-------------------------------------------------------------------------
  subroutine ob_interpolatedU(pointPhys, mlist, uave, levelp, gidp, rankp)
    use mpilib
    use grid_boundary
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN),dimension(:) :: mlist
    real(kind=DBL_KIND),dimension(size(mlist)),intent(OUT) :: uave
    integer,intent(OUT),optional ::  levelp, gidp, rankp
    integer:: ia, ja, ka, ib, jb, kb, m, gid, rank, level, plevel
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    myrank = get_myrank()
    do level =  LevelMax, Lmin, -1 ! from child to parent level
       call ob_getBlockFromPointPhys(pointPhys, level, gid, rank, ia,ja,ka)
       if ( ia == Imin .or. ia == Imax .or. &
            ja == Jmin .or. ja == Jmax .or. &
            ka == Kmin .or. ka == Kmax) &
            call boundary_grid( level, COMPLETE )
       ! (ia, ja, ka) がブロクの端の場合に袖にアクセスする*可能性*があるので、袖を fix する。
       ! 本当は ib が Imin-1 か Imax+1 の場合、
       ! または jb が Jmin-1 か Jmax+1 の場合、
       ! または kb が Kmin-1 か Kmax+1 の場合
       ! に boundary_grid を呼べば必要十分であるが、(ib, jb, kb) は rank ローカルな値のため、現状の仕様とした。
       if ( rank == myrank ) then
          ! next nearest index, in, jn, kn
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)
          ib = ia + int( sign(1.d0, pointPhys%p(MX) - x(ia) ) )
          jb = ja + int( sign(1.d0, pointPhys%p(MY) - y(ja) ) )
          kb = ka + int( sign(1.d0, pointPhys%p(MZ) - z(ka) ) )
          ! get gxave (gravity averaged in space) by linear interpolation in 3D [ia,ib][ja,jb][ka,kb]
          do m = lbound(mlist,1), ubound(mlist,1)
             uave(m) = trilinearInterp(ia, ja, ka, ib, jb, kb, &
                  pointPhys%p(MX), pointPhys%p(MY), pointPhys%p(MZ), mlist(m), gid)
          end do
       endif
       if (gid /= Undefi) then
          plevel = level
          exit ! exit if block is found
       endif
    end do
    if (gid == Undefi) then
       print *, '*** error in ob_interpolatedU.  No block is found.'
       call flush(6)
       uave = 0.d0
    end if
    call mpi_bcast(uave, size(uave), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    ! optional parameters
    if (present(gidp)) then
       gidp = gid
    end if
    if (present(levelp)) then
       levelp = plevel
    end if
    if (present(rankp)) rankp = rank
  contains
    !-------------------------------------------------------------------------
    ! tri-linear interpolation
    !-------------------------------------------------------------------------
    function trilinearInterp( ia, ja, ka, ib, jb, kb, x0, y0, z0, m, gid) result(u0)
      integer,intent(IN) :: ia, ja, ka, ib, jb, kb, m, gid
      real(kind=DBL_KIND),intent(IN) :: x0, y0, z0
      real(kind=DBL_KIND) :: u0, uaa, uba, uab, ubb, ua, ub
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: u
      u => get_Ucomp(m, gid)
      uaa = linearInterp(u(ia,ja,ka), u(ib,ja,ka), x(ia), x(ib), x0)
      uba = linearInterp(u(ia,jb,ka), u(ib,jb,ka), x(ia), x(ib), x0)
      uab = linearInterp(u(ia,ja,kb), u(ib,ja,kb), x(ia), x(ib), x0)
      ubb = linearInterp(u(ia,jb,kb), u(ib,jb,kb), x(ia), x(ib), x0)
      ua = linearInterp(uaa, uba, y(ja), y(jb), y0)
      ub = linearInterp(uab, ubb, y(ja), y(jb), y0)
      u0 = linearInterp(ua, ub, z(ka), z(kb), z0)
    end function trilinearInterp
    !-------------------------------------------------------------------------
    ! linear interpolation
    !-------------------------------------------------------------------------
    function linearInterp(ua, ub, xa, xb, x) result(u)
      real(kind=DBL_KIND),intent(IN) :: ua, ub, xa, xb, x
      real(kind=DBL_KIND) :: u
      u = (ub*(x - xa) + ua*(xb - x))/(xb-xa)
    end function linearInterp
  end subroutine ob_interpolatedU
  
end module ob_interp
