module ob_interp
  use overBlockCoordinates
  implicit none
  private
  public :: ob_interpolatedU, ob_interpolatedUbyCoordPhys
contains
  subroutine ob_interpolatedUbyCoordPhys(coords, mlist, uave, levelp, gidp, rankp)
    real(kind=8),dimension(0:2),intent(IN) :: coords
    integer,intent(IN),dimension(:) :: mlist
    real(kind=8),dimension(size(mlist)),intent(OUT) :: uave
    integer,intent(OUT),optional :: levelp, gidp, rankp
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coords, pointPhys)
    call ob_interpolatedU(pointPhys, mlist, uave, levelp, gidp, rankp)
  end subroutine ob_interpolatedUbyCoordPhys
  subroutine ob_interpolatedU(pointPhys, mlist, uave, levelp, gidp, rankp)
    use mpilib
    use grid_boundary
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN),dimension(:) :: mlist
    real(kind=8),dimension(size(mlist)),intent(OUT) :: uave
    integer,intent(OUT),optional :: levelp, gidp, rankp
    integer:: ia, ja, ka, ib, jb, kb, m, gid, rank, level, plevel
    real(kind=8),dimension(:),pointer :: x, y, z
    myrank = get_myrank()
    do level = LevelMax, Lmin, -1 
       call ob_getBlockFromPointPhys(pointPhys, level, gid, rank, ia,ja,ka)
       if ( ia == Imin .or. ia == Imax .or. &
            ja == Jmin .or. ja == Jmax .or. &
            ka == Kmin .or. ka == Kmax) &
            call boundary_grid( level, COMPLETE )
       if ( rank == myrank ) then
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)
          ib = ia + int( sign(1.d0, pointPhys%p(0) - x(ia) ) )
          jb = ja + int( sign(1.d0, pointPhys%p(1) - y(ja) ) )
          kb = ka + int( sign(1.d0, pointPhys%p(2) - z(ka) ) )
          do m = lbound(mlist,1), ubound(mlist,1)
             uave(m) = trilinearInterp(ia, ja, ka, ib, jb, kb, &
                  pointPhys%p(0), pointPhys%p(1), pointPhys%p(2), mlist(m), gid)
          end do
       endif
       if (gid /= Undefi) then
          plevel = level
          exit 
       endif
    end do
    if (gid == Undefi) then
       print *, '*** error in ob_interpolatedU.  No block is found.'
       call flush(6)
       uave = 0.d0
    end if
    call mpi_bcast(uave, size(uave), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    if (present(gidp)) then
       gidp = gid
    end if
    if (present(levelp)) then
       levelp = plevel
    end if
    if (present(rankp)) rankp = rank
  contains
    function trilinearInterp( ia, ja, ka, ib, jb, kb, x0, y0, z0, m, gid) result(u0)
      integer,intent(IN) :: ia, ja, ka, ib, jb, kb, m, gid
      real(kind=8),intent(IN) :: x0, y0, z0
      real(kind=8) :: u0, uaa, uba, uab, ubb, ua, ub
      real(kind=8),dimension(:,:,:),pointer :: u
      u => get_Ucomp(m, gid)
      uaa = linearInterp(u(ia,ja,ka), u(ib,ja,ka), x(ia), x(ib), x0)
      uba = linearInterp(u(ia,jb,ka), u(ib,jb,ka), x(ia), x(ib), x0)
      uab = linearInterp(u(ia,ja,kb), u(ib,ja,kb), x(ia), x(ib), x0)
      ubb = linearInterp(u(ia,jb,kb), u(ib,jb,kb), x(ia), x(ib), x0)
      ua = linearInterp(uaa, uba, y(ja), y(jb), y0)
      ub = linearInterp(uab, ubb, y(ja), y(jb), y0)
      u0 = linearInterp(ua, ub, z(ka), z(kb), z0)
    end function trilinearInterp
    function linearInterp(ua, ub, xa, xb, x) result(u)
      real(kind=8),intent(IN) :: ua, ub, xa, xb, x
      real(kind=8) :: u
      u = (ub*(x - xa) + ua*(xb - x))/(xb-xa)
    end function linearInterp
  end subroutine ob_interpolatedU
end module ob_interp
