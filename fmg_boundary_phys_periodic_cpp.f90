module fmg_boundary_phys
  implicit none
  private
  public :: fmg_boundary_u, fmg_boundary_physical_grid, vmg_boundary_u, mg_boundary_u
contains
  subroutine fmg_bndp_xyzminmax(xmin, ymin, zmin, xmax, ymax, zmax)
    use mpilib
    use grid, only : GidBase, RankBase, get_Xp, get_Yp, get_Zp, Imin, Jmin, Kmin, Imax, Jmax, Kmax
    real(kind=8),intent(OUT) :: xmin, ymin, zmin, xmax, ymax, zmax
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: buf(0:2)
    integer :: gid, rank
    gid = GidBase(lbound(GidBase,1), lbound(GidBase,2), lbound(GidBase,3))
    rank = RankBase(lbound(RankBase,1), lbound(RankBase,2), lbound(RankBase,3))
    myrank = get_myrank()
    if (myrank == rank) then
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       buf(:) = (/x(Imin), y(Jmin), z(Kmin)/) 
    endif
    call mpi_bcast( buf, size(buf), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmin = buf(0)
    ymin = buf(1)
    zmin = buf(2)
    gid = GidBase(ubound(GidBase,1), ubound(GidBase,2), ubound(GidBase,3))
    rank = RankBase(ubound(RankBase,1), ubound(RankBase,2), ubound(RankBase,3))
    myrank = get_myrank()
    if (myrank == rank) then
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       buf(:) = (/x(Imax), y(Jmax), z(Kmax)/) 
    endif
    call mpi_bcast( buf, size(buf), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmax = buf(0)
    ymax = buf(1)
    zmax = buf(2)
  end subroutine fmg_bndp_xyzminmax
  subroutine fmg_boundary_u(fmglev,ju)
    use fmg_data
    integer,intent(IN) :: fmglev, ju
  end subroutine fmg_boundary_u
  subroutine vmg_boundary_u(amrlev,fmglev,ju)
    use fmg_data
    integer,intent(IN) :: amrlev, fmglev, ju
  end subroutine vmg_boundary_u
  subroutine mg_boundary_u(mglev, ju)
    use mg_data
    integer,intent(IN) :: mglev, ju
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    integer :: imingh, imin, imaxgh, imax, i
    integer :: jmingh, jmin, jmaxgh, jmax, j
    integer :: kmingh, kmin, kmaxgh, kmax, k
    u => mg_get_arrp(mglev, ju)
    imingh = mg_get_imingh(mglev)
    imin = GridSize(mglev)%Imin
    imaxgh = mg_get_imaxgh(mglev)
    imax = GridSize(mglev)%Imax
    u(Imingh:Imin-1,:,:,:) = u(Imax+1-Ngh:Imax,:,:,:)
    u(Imax+1:Imaxgh,:,:,:) = u(Imin:Imin-1+Ngh,:,:,:)
    jmingh = mg_get_jmingh(mglev)
    jmin = GridSize(mglev)%Jmin
    jmaxgh = mg_get_jmaxgh(mglev)
    jmax = GridSize(mglev)%Jmax
    u(:,Jmingh:Jmin-1,:,:) = u(:,Jmax+1-Ngh:Jmax,:,:)
    u(:,Jmax+1:Jmaxgh,:,:) = u(:,Jmin:Jmin-1+Ngh,:,:)
    kmingh = mg_get_kmingh(mglev)
    kmin = GridSize(mglev)%Kmin
    kmaxgh = mg_get_kmaxgh(mglev)
    kmax = GridSize(mglev)%Kmax
    u(:,:,Kmingh:Kmin-1,:) = u(:,:,Kmax+1-Ngh:Kmax,:)
    u(:,:,Kmax+1:Kmaxgh,:) = u(:,:,Kmin:Kmin-1+Ngh,:)
  end subroutine mg_boundary_u
  subroutine fmg_boundary_physical_grid(ju, jrho)
    integer,intent(IN) :: ju, jrho
  end subroutine fmg_boundary_physical_grid
end module fmg_boundary_phys
