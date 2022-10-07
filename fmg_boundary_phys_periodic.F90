
#include "config.h"
!-------------------------------------------------------------------------
! Module for physical boundary condition of multigrid iteration.
! Periodic boundary condition in the x, y, z directions.
!-------------------------------------------------------------------------
module fmg_boundary_phys
  implicit none
  private
  public :: fmg_boundary_u, fmg_boundary_physical_grid, vmg_boundary_u, mg_boundary_u
contains
  ! ----------------------------------------------------------------
  ! base grids の boxsize (xmin, ymin, zmin, xmax, ymax, zmax) を求める。
  ! ----------------------------------------------------------------
  subroutine fmg_bndp_xyzminmax(xmin, ymin, zmin, xmax, ymax, zmax)
    use mpilib
    use grid, only : GidBase, RankBase, get_Xp, get_Yp, get_Zp, Imin, Jmin, Kmin, Imax, Jmax, Kmax
    real(kind=DBL_KIND),intent(OUT) :: xmin, ymin, zmin, xmax, ymax, zmax
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: buf(MX:MZ)
    integer :: gid, rank
    ! find xmin, ymin, zmin
    gid  =  GidBase(lbound(GidBase,1), lbound(GidBase,2), lbound(GidBase,3))
    rank = RankBase(lbound(RankBase,1), lbound(RankBase,2), lbound(RankBase,3))
    myrank = get_myrank()
    if (myrank == rank) then
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       buf(:) = (/x(Imin), y(Jmin), z(Kmin)/) ! min coordinates
    endif
    call mpi_bcast( buf, size(buf), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmin = buf(MX)
    ymin = buf(MY)
    zmin = buf(MZ)

    ! find xmax, ymax, zmax
    gid  =  GidBase(ubound(GidBase,1), ubound(GidBase,2), ubound(GidBase,3))
    rank = RankBase(ubound(RankBase,1), ubound(RankBase,2), ubound(RankBase,3))
    myrank = get_myrank()
    if (myrank == rank) then
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       buf(:) = (/x(Imax), y(Jmax), z(Kmax)/) ! max coordinates
    endif
    call mpi_bcast( buf, size(buf), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmax = buf(MX)
    ymax = buf(MY)
    zmax = buf(MZ)
  end subroutine fmg_bndp_xyzminmax
  ! ----------------------------------------------------------------
  ! 境界条件 (AMR FMG cycle)
  ! ----------------------------------------------------------------
  subroutine fmg_boundary_u(fmglev,ju)
    use fmg_data
    integer,intent(IN) :: fmglev, ju

  end subroutine fmg_boundary_u
  ! ----------------------------------------------------------------
  ! 境界条件 (AMR FAS V cycle)
  ! ----------------------------------------------------------------
  subroutine vmg_boundary_u(amrlev,fmglev,ju)
    use fmg_data
    integer,intent(IN) :: amrlev, fmglev, ju

  end subroutine vmg_boundary_u
  ! ----------------------------------------------------------------
  ! ポテンシャル残差の境界条件 (Base Grid FMG cycle)
  ! ----------------------------------------------------------------
  subroutine mg_boundary_u(mglev, ju)
    use mg_data
    integer,intent(IN) :: mglev, ju
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    integer :: imingh, imin, imaxgh, imax, i
    integer :: jmingh, jmin, jmaxgh, jmax, j
    integer :: kmingh, kmin, kmaxgh, kmax, k
    u => mg_get_arrp(mglev, ju)

    ! i-direction
    imingh =  mg_get_imingh(mglev)
    imin = GridSize(mglev)%Imin
    imaxgh =  mg_get_imaxgh(mglev)
    imax = GridSize(mglev)%Imax
    u(Imingh:Imin-1,:,:,:) = u(Imax+1-Ngh:Imax,:,:,:)
    u(Imax+1:Imaxgh,:,:,:) = u(Imin:Imin-1+Ngh,:,:,:)

    ! j-direction
    jmingh =  mg_get_jmingh(mglev)
    jmin = GridSize(mglev)%Jmin
    jmaxgh =  mg_get_jmaxgh(mglev)
    jmax = GridSize(mglev)%Jmax
    u(:,Jmingh:Jmin-1,:,:) = u(:,Jmax+1-Ngh:Jmax,:,:)
    u(:,Jmax+1:Jmaxgh,:,:) = u(:,Jmin:Jmin-1+Ngh,:,:)

    ! k-direction
    kmingh =  mg_get_kmingh(mglev)
    kmin = GridSize(mglev)%Kmin
    kmaxgh =  mg_get_kmaxgh(mglev)
    kmax = GridSize(mglev)%Kmax
    u(:,:,Kmingh:Kmin-1,:) = u(:,:,Kmax+1-Ngh:Kmax,:)
    u(:,:,Kmax+1:Kmaxgh,:) = u(:,:,Kmin:Kmin-1+Ngh,:)

  end subroutine mg_boundary_u
  ! ----------------------------------------------------------------
  ! 物理境界に接するグリッドに境界条件を与える
  ! 境界条件は最粗グリッドの 多重極展開により、求める。
  ! ----------------------------------------------------------------
  subroutine fmg_boundary_physical_grid(ju, jrho)
    integer,intent(IN) :: ju, jrho

  end subroutine fmg_boundary_physical_grid
end module fmg_boundary_phys
