! ----------------------------------------------------------------------

! module for soft-clustering
! ----------------------------------------------------------------------
#include "config.h"
module softClustering
  implicit none
  private
  public :: clst_make
contains
  ! ----------------------------------------------------------------------
  ! Clusterging points
  ! INPUT:
  !   rp(dim, npoint) = Points to be considered.
  !   mass(npoint)    = Masses of points.
  !   radiiCluster    = Radii of clusters.
  ! OUTPUT:
  !   c(npoint)        = Cluster numbers of points.
  !                      Cluster number begins with lbound(c,1).
  !   rg(dim,npoint)   = Baricenters of clusters.
  !   ncluster         = Number of clusters.
  ! ----------------------------------------------------------------------
  subroutine clst_make(rp, mass, radiiCluster, c, rg, ncluster)
    real(kind=DBL_KIND),dimension(:,:),intent(IN) :: rp
    real(kind=DBL_KIND),dimension(:),intent(IN) :: mass
    real(kind=DBL_KIND),intent(IN) :: radiiCluster
    integer,dimension(:),intent(OUT) :: c
    real(kind=DBL_KIND),dimension(:,:),intent(OUT) :: rg
    integer,intent(OUT) :: ncluster
    logical,dimension(size(rp,2),size(rp,2)) :: mask
    real(kind=DBL_KIND),dimension(size(rp,2),size(rp,2)) :: dr2
    real(kind=DBL_KIND),dimension(lbound(rp,1):ubound(rp,1)) :: rgij
    real(kind=DBL_KIND),dimension(size(rp,2)) :: dist2
    integer,dimension(2) :: pos
    logical :: bool_change
    integer :: i, j, n, inew, jnew, ci, cj, cmin, cmax, cn, co
    real(kind=DBL_KIND) :: dist2max

    ! check arrary sizes
    n = size(rp,2)
    if ( size(mass) /= n .or. size(c) /= n .or. size(rg,2) /= n) then
       print *, '*** errror in clst_make: array sizes are not consistent.'
       stop
    end if
    n = lbound(rp,2)
    if ( lbound(mass,1) /= n .or. lbound(c,1) /= n .or. lbound(rg,2) /= n) then
       print *, '*** errror in clst_make: array offsets are not consistent.'
       stop
    end if

    ! make mask initialize
    do j = lbound(mask,2), ubound(mask,2)
       do i = lbound(mask,1), ubound(mask,1)
          if ( i > j ) then
             mask(i,j) = .true.
          else
             mask(i,j) = .false.
          end if
       end do
    end do

    ! initial clustering
    do j = lbound(c,1), ubound(c,1)
       c(j) = j
    end do
    rg = rp

    bool_change = .true.
    do while (bool_change)
       if (.not. any(mask) ) exit

       ! update distance
       dr2 = 0
       do n = lbound(rg,1), ubound(rg,1) ! dimension
          do j = lbound(dr2,2), ubound(dr2,2)
             do i = lbound(dr2,1), ubound(dr2,1)
                if ( .not. mask(i,j) ) cycle
                dr2(i,j) = dr2(i,j) + (rg(n,i)-rg(n,j))**2
             end do
          end do
       end do
       pos = minloc(dr2, mask=mask)
       inew = pos(1)
       jnew = pos(2)
       ! baricenter of i - j
       rgij(:) = (rg(:,inew)*mass(inew) + rg(:,jnew)*mass(jnew))/(mass(inew) + mass(jnew))
       ci = c(inew)
       cj = c(jnew)
       ! maximum distance from baricenter
       dist2 = 0
       do n = lbound(rp,1), ubound(rp,1) !dimension
          dist2(:) = dist2(:) + (rp(n,:)-rgij(n))**2
       end do
       dist2max = maxval( dist2, mask = (c(:) == ci .or. c(:) == cj) )
       if (dist2max <= radiiCluster**2) then
!!$          print *, ci , ' -> ' , cj
          do j = lbound(c,1), ubound(c,1)
             if ( c(j) == ci .or. c(j) == cj ) then
                c(j) = cj
                rg(:,j) = rgij(:)
             end if
          end do
          ! update mask(inew, jnew)
          do j = lbound(mask,2), ubound(mask,2)
             do i = lbound(mask,1), ubound(mask,1)
                if (c(i) == cj .and.  c(j) == cj) then
                   mask(i,j) = .false.
                end if
             end do
          end do
          bool_change = .true.
       else
          bool_change = .false.
       end if
    end do

    ! packing cluster number
    cmin = minval(c)
    cmax = maxval(c)
    cn = lbound(c,1)
    do co = cmin, cmax
       where(c == co)
          c = cn
       end where
       if (any( c == cn )) cn = cn + 1
    end do
    ncluster = maxval(c)+lbound(c,1)-1

    ! returned values, rg
    do j = lbound(c,1), lbound(c,1) + ncluster - 1
       do i = lbound(c,1), ubound(c,1)
          if (c(i) == j) then
             rg(:,j) = rg(:,i)
             exit
          end if
       end do
    end do
  end subroutine clst_make
end module softClustering
