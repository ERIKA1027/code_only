!
! Write Snapshot of data
!
#include "config.h"
module writeSnap
  use uniformgrid
  implicit none
  private
  public :: writeSnap_whole
#ifdef WITH_SELFGRAVITY
  public :: writeSnap_denseRegion
  public :: writeSnap_clusters
#endif !WITH_SELFGRAVITY
contains
  !---------------------------------------------------------------------
  ! write data of whole region of computational domain
  !---------------------------------------------------------------------
  subroutine writeSnap_whole
    use grid, only : Lmin
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    integer,parameter :: SampleRate = 2 ! sampling rate ( full sampling means 1)
    integer :: baseLevel
    real(KIND=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    call ob_computationBoxOfCoordPhys( coordPhys )
    xmin = coordPhys(MX)
    ymin = coordPhys(MY)
    zmin = coordPhys(MZ)
    xmax = coordPhys(MZ+1+MX)
    ymax = coordPhys(MZ+1+MY)
    zmax = coordPhys(MZ+1+MZ)
    baseLevel = Lmin - int( log10(dble(SampleRate))/log10(2.d0) + 0.5)
    if ( bool_skip( baseLevel ) ) return
    call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, baseLevel,interpolate=.true.)
  end subroutine writeSnap_whole
#ifdef WITH_SELFGRAVITY
  !---------------------------------------------------------------------
  ! write data of dense region
  ! INPUT:
  !   radius = radius of output region
  !   prefix = prefix of output filename (optional)
  !---------------------------------------------------------------------
  subroutine writeSnap_denseRegion( radius, prefix )
    use grid, only : LevelMax, Lmin, CellWidth, Undefi
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    real(kind=DBL_KIND),intent(IN) :: radius
    character(len=*),intent(IN),optional :: prefix
    integer,parameter :: NIug=64, NJug=NIug, NKug=NIug ! minimum resolution
    real(kind=DBL_KIND) :: dhReq                         ! mesh size required
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    integer :: baseLevel
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    !  check resolution
    dhReq = radius*2 / max(NIug, NJug, NKug)
    baseLevel = int(log(minval(CellWidth(:,Lmin))/dhReq)/log(2.d0)+0.5d0)
    if ( bool_skip( baseLevel ) ) return
    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )
    xmin = coordPhys(MX)
    ymin = coordPhys(MY)
    zmin = coordPhys(MZ)
    xmax = coordPhys(MZ+1+MX)
    ymax = coordPhys(MZ+1+MY)
    zmax = coordPhys(MZ+1+MZ)
    ! 最大密度の位置
    call rhomaxPosition(xp, yp, zp)
    ! 領域の AND をとる
    xmin = max(xmin, xp-radius)
    ymin = max(ymin, yp-radius)
    zmin = max(zmin, zp-radius)
    xmax = min(xmax, xp+radius)
    ymax = min(ymax, yp+radius)
    zmax = min(zmax, zp+radius)
    call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, baseLevel,interpolate=.true., prefix=prefix)
  contains
    !---------------------------------------------------------------------
    ! get a rectangle region which a given level covers
    !---------------------------------------------------------------------
    subroutine rhomaxPosition(xp, yp, zp)
      use mpilib
      use grid, only : get_Ucomp, get_Xp, get_Yp, get_Zp, Undefi, Gidmin, GidListMax, GidList, &
           Imin, Jmin, Kmin, Imax, Jmax, Kmax
      real(kind=DBL_KIND),intent(OUT) :: xp, yp, zp
      real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
      real(kind=DBL_KIND),dimension(0:NPE-1) :: rhomaxNode, rhomaxNoder
      real(kind=DBL_KIND),dimension(MX:MZ) ::  buf
      real(kind=DBL_KIND) :: rhomax, xmax, ymax, zmax
      integer :: n, gid, gidmax, maxl(MX:MZ), nodemax(1)
      ! find gid of maximum rho
      myrank = get_myrank()
      gidmax = Undefi
      rhomaxNode(:) = 0.d0
      do n = Gidmin, GidListMax(LevelMax)
         gid = GidList(n, LevelMax)
         rho => get_Ucomp(MRHO, gid)
         rhomax = maxval(rho(ARRAYSIZE_IJK))
         if (rhomax > rhomaxNode(myrank)) then
            rhomaxNode(myrank) = rhomax
            gidmax = gid
         endif
      enddo
      ! coordinates of rhomax for each node
      if ( gidmax /= Undefi ) then
         rho => get_Ucomp(MRHO, gidmax)
         x => get_Xp(gidmax)
         y => get_Yp(gidmax)
         z => get_Zp(gidmax)
         maxl = maxloc(rho(ARRAYSIZE_IJK))
         maxl = maxl - 1 + (/Imin, Jmin, Kmin/)
         xmax = x(maxl(MX))
         ymax = y(maxl(MY))
         zmax = z(maxl(MZ))
      endif
      ! find node of maximum rho
      call mpi_allreduce(rhomaxNode, rhomaxNoder, size(rhomaxNode), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      nodemax = maxloc(rhomaxNoder) - 1
      buf = (/ xmax, ymax, zmax /)
      call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, nodemax(1), MPI_COMM_WORLD, ierr)
      xp = buf(MX) ; yp = buf(MY) ; zp = buf(MZ)
    end subroutine rhomaxPosition
  end subroutine writeSnap_denseRegion
#endif !WITH_SELFGRAVITY
#ifdef WITH_SELFGRAVITY
  !---------------------------------------------------------------------
  ! write data of dense regions
  ! this routine can output **several** regions
  ! INPUT:
  !   radius = radius of output region
  !---------------------------------------------------------------------
  subroutine writeSnap_clusters ( radius, prefix )
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use string
    real(kind=DBL_KIND),intent(IN) :: radius
    character(len=*),intent(IN),optional :: prefix
    integer,parameter :: NREGIONS = 100                 ! maximum number of regions
    !integer,parameter :: NIug=128, NJug=NIug, NKug=NIug ! minimum resolution
    integer,parameter :: NIug=64, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 20 MB/file)
    !integer,parameter :: NIug=32, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 2 MB/file)
    real(kind=DBL_KIND) :: dhReq                         ! mesh size required
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    integer :: baseLevel, np, n
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    real(kind=DBL_KIND),dimension(MX:MZ,NREGIONS) :: pr
    character(len=CHARLEN) :: prf

    !if (Step(Lmin) <= 0) return
    if (Step(Lmin) < 0) return

    !  check resolution, and define baseLevel
    dhReq = radius*2 / max(NIug, NJug, NKug)
    baseLevel = floor(log(minval(CellWidth(:,Lmin))/dhReq)/log(2.d0)+0.5d0) !KS MODIFIED (自分としてはこっち)
    ! baseLevel = int(log(minval(CellWidth(:,Lmin))/dhReq)/log(2.d0)+0.5d0)    
    if ( bool_skip( baseLevel ) ) return
    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )
    ! pr ... location of cluster
    ! np ... number of clusters
    call mkPosition(pr, np)
    do n = 1, np
       xmin = coordPhys(MX)
       ymin = coordPhys(MY)
       zmin = coordPhys(MZ)
       xmax = coordPhys(MZ+1+MX)
       ymax = coordPhys(MZ+1+MY)
       zmax = coordPhys(MZ+1+MZ)
       ! 領域の AND をとる
       xmin = max(xmin, pr(MX,n)-radius)
       ymin = max(ymin, pr(MY,n)-radius)
       zmin = max(zmin, pr(MZ,n)-radius)
       xmax = min(xmax, pr(MX,n)+radius)
       ymax = min(ymax, pr(MY,n)+radius)
       zmax = min(zmax, pr(MZ,n)+radius)
       if (present( prefix )) then
          prf = prefix
       else
          prf = 'cl.'
       endif
       prf = concat(concat(prf, num2char(n)),'.')
       call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, baseLevel,interpolate=.true., prefix=prf)
    end do
  contains
    !---------------------------------------------------------------------
    ! get a positions where gravitational potential (PSI) has local minimum
    !---------------------------------------------------------------------
    subroutine mkPosition(pr, np)
      use mpilib
      use softClustering
#ifdef SINKPARTICLE
      use sinkparticle
#endif !SINKPARTICLE
      real(kind=DBL_KIND),dimension(MX:MZ,NREGIONS),intent(OUT) :: pr
      integer,intent(OUT) :: np
      !
      real(kind=DBL_KIND),dimension(MX:MZ,NREGIONS) :: rSink, rPsi
      real(kind=DBL_KIND),dimension(NREGIONS) :: mSink, mPsi
      integer :: nSink, nPsi
      !
      real(kind=DBL_KIND),dimension(:,:),allocatable :: rp, rg
      real(kind=DBL_KIND),dimension(:),allocatable :: mass
      integer,dimension(:),allocatable :: c
      integer :: ncluster
      integer :: n                        ! for debug
#ifdef SINKPARTICLE
      call sp_sinkdata2array(nSink, mSink, pr=rSink) ! from sink particle
#else
      nSink = 0
#endif !SINKPARTICLE
      
      call findPsiMin(rPsi, mPsi, nPsi)   ! from psi
      call prAllGather(rPsi, mPsi, nPsi)  ! allgather of rPsi, and mPsi
      ! merge list : (rSink, rPsi) => rSink
      if (nSink + nPsi > NREGIONS) then
         print *, '*** error in mkPosition: nSink + nPsi > NREGIONS', nSink, nPsi, NREGIONS
      end if
      if ( nPsi > 0 ) then
         rSink(:, nSink+1: nSink+nPsi) =rPsi(:,1:nPsi)
         mSink(   nSink+1: nSink+nPsi) =mPsi(1:nPsi)
         nSink = nSink + nPsi
      end if

      ! solve clusters
      allocate( rp(MX:MZ,nSink), mass(nSink), c(nSink), rg(MX:MZ,nSink) )
      rp(MX:MZ,1:nSink) = rSink(MX:MZ,1:nSink)
      mass(    1:nSink) = mSink(      1:nSink)
      call clst_make(rp, mass, radius, c, rg, ncluster)
      np = ncluster             ! returned value np
      pr(:,1:np) = rg(:,1:np)   ! returned value pr
      deallocate(rp, mass, c, rg)

    end subroutine mkPosition
    !---------------------------------------------------------------------
    ! find local minimum of gravitational potential for every node in node-local
    !---------------------------------------------------------------------
    subroutine findPsiMin(pr, pmass, np)
      use mpilib
      use parameter
      use grid, only : Undefi, get_Ucomp, get_Xp, get_Yp, get_Zp, Undefi, Gidmin, GidListMax, GidList, &
           Imin, Jmin, Kmin, Imax, Jmax, Kmax, ChildGid, Left, get_dv
      real(kind=DBL_KIND),dimension(MX:MZ,NREGIONS),intent(OUT) :: pr
      real(kind=DBL_KIND),dimension(NREGIONS),intent(OUT) :: pmass
      integer,intent(OUT) :: np
      real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: psi, rho
      real(kind=DBL_KIND) :: dv
      integer :: i, j, k, gid
      dv = get_dv(LevelMax)
      np = 0
      myrank = get_myrank()

      if (Step(Lmin) <= 0) return ! skip (KS MODIFIED)
      
      do n = Gidmin, GidListMax(LevelMax)
         gid = GidList(n, LevelMax)
         if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle ! gid is not finest at this point
         psi => get_Ucomp(MPSI, gid)
         rho => get_Ucomp(MRHO, gid)
         x => get_Xp(gid)
         y => get_Yp(gid)
         z => get_Zp(gid)
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  if ( psi(i,j,k) > minval(psi(i-1:i+1,j-1:j+1,k-1:k+1))) cycle ! local minimum?
                  np = np + 1
                  pr(:,np) = (/ x(i), y(j), z(k) /)
                  pmass(np) = rho(i,j,k) * dv ! add-hock
               end do
            end do
         end do
      end do
    end subroutine findPsiMin
    !---------------------------------------------------------------------
    ! allgather (pr, np) among every ranks
    !---------------------------------------------------------------------
    subroutine prAllGather(pr, pmass, np)
      use mpilib
      real(kind=DBL_KIND),dimension(MX:MZ,NREGIONS),intent(INOUT) :: pr
      real(kind=DBL_KIND),dimension(NREGIONS),intent(INOUT) :: pmass
      integer,intent(INOUT) :: np
      real(kind=DBL_KIND),dimension(MX:MZ+1,NREGIONS) :: rbuf, sbuf
      integer :: npoint, rank, nelem
      integer,dimension(0:NPE-1) :: recvcounts, displs
      call mpi_allreduce(np, npoint, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (npoint > NREGIONS) then
         print *, '*** error in prAllGather: npoint > NREGIONS', npoint, NREGIONS
      end if
      nelem = np * size(sbuf,1) ! number of elements in sbuf
      call mpi_allgather( &
           nelem,      1, MPI_INTEGER, &
           recvcounts, 1, MPI_INTEGER,  MPI_COMM_WORLD, ierr)
      displs(0) = 0
      do rank = 1, NPE-1
         displs(rank) = displs(rank-1) + recvcounts(rank-1)
      end do
      sbuf(MX:MZ,1:np) = pr(MX:MZ,1:np)
      sbuf(MZ+1,1:np) = pmass(1:np)
      call mpi_allgatherv( &
           sbuf, nelem, MPI_DOUBLE_PRECISION, &
           rbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      pr(MX:MZ,1:npoint)  = rbuf(MX:MZ,1:npoint)
      pmass(1:npoint) = rbuf(MZ+1,1:npoint)
      np = npoint
    end subroutine prAllGather

  end subroutine writeSnap_clusters
#endif !WITH_SELFGRAVITY
  !---------------------------------------------------------------------
  ! bool skip to write filename
  !---------------------------------------------------------------------
  function bool_skip( baseLevel ) result( bool )
    use grid, only : LevelMax
    integer,intent(IN) :: baseLevel
    logical :: bool
    integer,parameter :: DlevelMax = 3
    bool = .false.              ! defaut is no skip
    if ( baseLevel - LevelMax > DlevelMax ) bool = .true.
  end function bool_skip

end module writeSnap
