! #define DEBUG
#include "config.h"
!-------------------------------------------------------------------------
! Module for physical boundary condition of multigrid iteration.
! Periodic boundary condition in the z directions.
! Fix boundary on the cylindrical surface.
!-------------------------------------------------------------------------
module fmg_boundary_phys
  implicit none
  private
  public :: fmg_boundary_u, fmg_boundary_physical_grid, vmg_boundary_u, mg_boundary_u
contains
  ! ----------------------------------------------------------------
  ! base grids の boxsize (xmin, ymin, zmin, xmax, ymax, zmax) を求める。
  ! ----------------------------------------------------------------
  subroutine fmg_bndp_xyzminmax(xmin, ymin, zmin, xmax, ymax, zmax, xlen, ylen, zlen)
    use mpilib
    use grid, only : GidBase, RankBase, get_Xp, get_Yp, get_Zp, Imin, Jmin, Kmin, Imax, Jmax, Kmax, CellWidth, Lmin
    real(kind=DBL_KIND),intent(OUT) :: xmin, ymin, zmin, xmax, ymax, zmax, xlen, ylen, zlen
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

    xlen = xmax - xmin + CellWidth(MX, Lmin)
    ylen = ymax - ymin + CellWidth(MY, Lmin)
    zlen = zmax - zmin + CellWidth(MZ, Lmin)
  end subroutine fmg_bndp_xyzminmax
  ! ----------------------------------------------------------------
  ! 境界条件 (AMR FMG cycle)
  ! ----------------------------------------------------------------
  subroutine fmg_boundary_u(fmglev,ju)
    use fmg_data
    integer,intent(IN) :: fmglev, ju
    integer :: amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call vmg_boundary_u(amrlev, fmglev,ju)
    end do
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
    integer :: kmingh, kmin, kmaxgh, kmax, k
    u => mg_get_arrp(mglev, ju)
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
    use modelParameter, only : MP_rhomask
    use parameter, only : Pi2
    use fmg_data
    use grid, only : get_Xp, get_Yp, get_Zp, GidList, GidListMax, get_dv, &
         Imin, Jmin, Kmin, Imax, Jmax, Kmax, get_Ucomp
    integer,intent(IN) :: ju, jrho
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, u
    real(kind=DBL_KIND),dimension(1) :: abuf, abufd
    integer :: i,j,k, amrgid, amrlev, fmglev, gid, la
    real(kind=DBL_KIND) :: pi,ub,&
         a00,a10,a11r,a11i,a20,a21r,a21i,a22r,a22i,a30,a31r,a31i, &
         a32r,a32i,a33r,a33i, &
         r2,r3,r,rc,logrc, cost,sint,cosp,sinp,ri,r2i,r3i,r4i,rci,dv,rhodv, &
         xmin, ymin, zmin, xmax, ymax, zmax, xlen, ylen, zlen
    if (FMG_PDE_TYPE /= FMG_PDE_TYPE_POISSON_EQUATION) then
       print *, '**** error in fmg_boundary_phys_sphere. invalid FMG_PDE_TYPE', FMG_PDE_TYPE
    endif
    ! --------------
    ! define levels
    ! --------------
    fmglev = FMG_LevelMin       ! real grid only

    ! ----------------
    ! 密度によるマスク
    ! ----------------
    call fmg_bndp_xyzminmax(xmin, ymin, zmin, xmax, ymax, zmax, xlen, ylen, zlen)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = lbound(GidList,1), GidListMax( amrlev ) ! FMG gid
          amrgid = GidList(gid, amrlev)                 ! AMR gid
          call fmg_arrp(amrlev, fmglev, gid, jrho, rho)
          where (rho < MP_rhomask) rho = 0.d0
       enddo
    enddo
    ! --------
    ! 展開係数
    ! --------
    amrlev = AMR_LevelMin       ! coarsest grid only
    pi = 4.d0*atan(1.d0)
    dv = get_dv(amrlev)
    a00 = 0
    do gid = lbound(GidList,1), GidListMax( amrlev ) ! FMG gid
       amrgid = GidList(gid, amrlev)                 ! AMR gid
       call fmg_arrp(amrlev, fmglev, gid, jrho, rho)
       x => get_Xp(amrgid)
       y => get_Yp(amrgid)
       z => get_Zp(amrgid)
       do k=Kmin,Kmax
          do j=Jmin,Jmax
             do i=Imin,Imax
                rhodv = rho(i,j,k)*dv
                a00 = a00 +  rhodv
             enddo
          enddo
       enddo
    enddo
    abuf(1) =a00
    call mpi_allreduce(abuf, abufd, size(abuf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    a00  = abufd(1)
    a00 = a00 / zlen / Pi2 ! phi = 2GM ln(r) = (4 pi G M)/(2 pi) ln(r). Note that rho includes 4 pi G.
#ifdef DEBUG
    if ( myrank == 0 ) then
       PRINTV(a00)
       PRINTV(zlen)
       PRINTV(rad2)
    endif
#endif
    ! ---------------------
    ! m=±L の場合は2倍する
    ! ---------------------
#define U00 logrc*a00

    ! --------------------------
    ! グリッドに境界条件を与える
    ! --------------------------
#define UB_(I_,IVAL_,UBV_) \
    I_=IVAL_         ;\
    rc = sqrt(x(i)**2 + y(j)**2) ;\
    logrc = log(rc)      ;\
    UBV_ = U00
    do la = AMR_LevelMin, AMR_LevelMax  ! for all AMR level
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)  ! FMG gid
          amrgid = GidList(gid, la)                     ! AMR gid
          call fmg_arrp(la, fmglev, gid, ju, u)
          x => get_Xp(amrgid)
          y => get_Yp(amrgid)
          z => get_Zp(amrgid)
          if (Geom(la)%Block(gid)%TouchBoundary(Left,MX)) then
             do k=lbound(u,3),ubound(u,3)
                do j=lbound(u,2),ubound(u,2)
                   UB_(i,Imin-1,ub)
                   u(Imin-1,j,k) = ub
                enddo
             enddo
          endif
          if (Geom(la)%Block(gid)%TouchBoundary(Right,MX)) then
             do k=lbound(u,3),ubound(u,3)
                do j=lbound(u,2),ubound(u,2)
                   UB_(i,Imax+1,ub)
                   u(Imax+1,j,k) = ub
                enddo
             enddo
          endif

          if (Geom(la)%Block(gid)%TouchBoundary(Left,MY)) then
             do k=lbound(u,3),ubound(u,3)
                do i=lbound(u,1),ubound(u,1)
                   UB_(j,Jmin-1,ub)
                   u(i,Jmin-1,k) = ub
                enddo
             enddo
          endif
          if (Geom(la)%Block(gid)%TouchBoundary(Right,MY)) then
             do k=lbound(u,3),ubound(u,3)
                do i=lbound(u,1),ubound(u,1)
                   UB_(j,Jmax+1,ub)
                   u(i,Jmax+1,k) = ub
                enddo
             enddo
          endif

       end do
    enddo

  end subroutine fmg_boundary_physical_grid
end module fmg_boundary_phys
