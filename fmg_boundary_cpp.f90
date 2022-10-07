module fmg_boundary
  use fmg_data
  implicit none
contains
  subroutine fmg_boundary_fill0(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    real(kind=8),parameter :: zero = 0.d0
    integer :: amrlev, id
    integer :: imin, jmin, kmin, imax, jmax, kmax
    integer :: imingh, jmingh, kmingh, imaxgh, jmaxgh, kmaxgh
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    imin = GridSize(fmglev)%Imin
    jmin = GridSize(fmglev)%Jmin
    kmin = GridSize(fmglev)%Kmin
    imax = GridSize(fmglev)%Imax
    jmax = GridSize(fmglev)%Jmax
    kmax = GridSize(fmglev)%Kmax
    imingh = fmg_get_imingh(fmglev)
    jmingh = fmg_get_jmingh(fmglev)
    kmingh = fmg_get_kmingh(fmglev)
    imaxgh = fmg_get_imaxgh(fmglev)
    jmaxgh = fmg_get_jmaxgh(fmglev)
    kmaxgh = fmg_get_kmaxgh(fmglev)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          up => fmg_get_arrp(amrlev, fmglev, id, icode)
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) &
               up(imingh:imin-1,:,:,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) &
               up(imax+1:imaxgh,:,:,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) &
               up(:,jmingh:jmin-1,:,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) &
               up(:,jmax+1:jmaxgh,:,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) &
               up(:,:,kmingh:kmin-1,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) &
               up(:,:,kmax+1:kmaxgh,:) = zero
       enddo
    enddo
  end subroutine fmg_boundary_fill0
  subroutine fmg_boundary_fill0_f(fmglev)
    integer,intent(IN) :: fmglev
    real(kind=8),parameter :: zero = 0.d0
    integer :: amrlev, id
    integer :: imin, jmin, kmin, imax, jmax, kmax
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fp
    call fmg_get_gridsize(fmglev, imin,jmin,kmin,imax,jmax,kmax)
    imin = imin - 1
    jmin = jmin - 1
    kmin = kmin - 1
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          fp => fmg_get_fp(amrlev, fmglev, id)
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) &
               fp(imin,:,:,0,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) &
               fp(imax,:,:,0,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) &
               fp(:,jmin,:,1,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) &
               fp(:,jmax,:,1,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) &
               fp(:,:,kmin,2,:) = zero
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) &
               fp(:,:,kmax,2,:) = zero
       enddo
    enddo
  end subroutine fmg_boundary_fill0_f
  subroutine fmg_boundary_extrap(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    integer :: id, amrlev
    integer :: i, j, k
    integer :: imin, jmin, kmin, imax, jmax, kmax
    integer :: imingh, jmingh, kmingh, imaxgh, jmaxgh, kmaxgh
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    imin = GridSize(fmglev)%Imin
    jmin = GridSize(fmglev)%Jmin
    kmin = GridSize(fmglev)%Kmin
    imax = GridSize(fmglev)%Imax
    jmax = GridSize(fmglev)%Jmax
    kmax = GridSize(fmglev)%Kmax
    imingh = fmg_get_imingh(fmglev)
    jmingh = fmg_get_jmingh(fmglev)
    kmingh = fmg_get_kmingh(fmglev)
    imaxgh = fmg_get_imaxgh(fmglev)
    jmaxgh = fmg_get_jmaxgh(fmglev)
    kmaxgh = fmg_get_kmaxgh(fmglev)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          up => fmg_get_arrp(amrlev, fmglev, id, icode)
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
             do i = imin-1, imingh, -1
                up(i,:,:,:) = -up(i+1,:,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
             do i = imax+1, imaxgh
                up(i,:,:,:) = -up(i-1,:,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
             do j = jmin-1, jmingh, -1
                up(:,j,:,:) = -up(:,j+1,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
             do j = jmax+1, jmaxgh
                up(:,j,:,:) = -up(:,j-1,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
             do k = kmin-1, kmingh, -1
                up(:,:,k,:) = -up(:,:,k+1,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
             do k = kmax+1, kmaxgh
                up(:,:,k,:) = -up(:,:,k-1,:)
             end do
          end if
       enddo
    end do
  end subroutine fmg_boundary_extrap
  subroutine fmg_boundary_extrap_BAK(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    integer :: id, amrlev
    integer :: i, j, k
    integer :: imin, jmin, kmin, imax, jmax, kmax
    integer :: imingh, jmingh, kmingh, imaxgh, jmaxgh, kmaxgh
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    imin = GridSize(fmglev)%Imin
    jmin = GridSize(fmglev)%Jmin
    kmin = GridSize(fmglev)%Kmin
    imax = GridSize(fmglev)%Imax
    jmax = GridSize(fmglev)%Jmax
    kmax = GridSize(fmglev)%Kmax
    imingh = fmg_get_imingh(fmglev)
    jmingh = fmg_get_jmingh(fmglev)
    kmingh = fmg_get_kmingh(fmglev)
    imaxgh = fmg_get_imaxgh(fmglev)
    jmaxgh = fmg_get_jmaxgh(fmglev)
    kmaxgh = fmg_get_kmaxgh(fmglev)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          up => fmg_get_arrp(amrlev, fmglev, id, icode)
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
             do i = imin-1, imingh, -1
                up(i,:,:,:) = 2*up(i+1,:,:,:) - up(i+2,:,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
             do i = imax+1, imaxgh
                up(i,:,:,:) = 2*up(i-1,:,:,:) - up(i-2,:,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
             do j = jmin-1, jmingh, -1
                up(:,j,:,:) = 2*up(:,j+1,:,:) - up(:,j+2,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
             do j = jmax+1, jmaxgh
                up(:,j,:,:) = 2*up(:,j-1,:,:) - up(:,j-2,:,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
             do k = kmin-1, kmingh, -1
                up(:,:,k,:) = 2*up(:,:,k+1,:) - up(:,:,k+2,:)
             end do
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
             do k = kmax+1, kmaxgh
                up(:,:,k,:) = 2*up(:,:,k-1,:) - up(:,:,k-2,:)
             end do
          end if
       enddo
    end do
  end subroutine fmg_boundary_extrap_BAK
  subroutine fmg_boundary_minmax(fmglev, icode)
    use mpilib
    integer,intent(IN) :: fmglev, icode
    real(kind=8) :: amin, amax, aming, amaxg
    integer :: amrlev, id
    integer :: imin, jmin, kmin, imax, jmax, kmax
    integer :: imingh, jmingh, kmingh, imaxgh, jmaxgh, kmaxgh
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    logical :: boolb, boolg
    imin = GridSize(fmglev)%Imin
    jmin = GridSize(fmglev)%Jmin
    kmin = GridSize(fmglev)%Kmin
    imax = GridSize(fmglev)%Imax
    jmax = GridSize(fmglev)%Jmax
    kmax = GridSize(fmglev)%Kmax
    imingh = fmg_get_imingh(fmglev)
    jmingh = fmg_get_jmingh(fmglev)
    kmingh = fmg_get_kmingh(fmglev)
    imaxgh = fmg_get_imaxgh(fmglev)
    jmaxgh = fmg_get_jmaxgh(fmglev)
    kmaxgh = fmg_get_kmaxgh(fmglev)
    amax = 0.d0
    amin = 0.d0
    boolb = .false.
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(id, amrlev, fmglev) ) cycle
          up => fmg_get_arrp(amrlev, fmglev, id, icode)
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
             amin = min(amin, minval(up(imingh:imin-1,:,:,:)))
             amax = max(amax, maxval(up(imingh:imin-1,:,:,:)))
             boolb = .true.
          endif
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
             amin = min(amin, minval(up(imax+1:imaxgh,:,:,:)))
             amax = max(amax, maxval(up(imax+1:imaxgh,:,:,:)))
             boolb = .true.
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
             amin = min(amin, minval(up(:,jmingh:jmin-1,:,:)))
             amax = max(amax, maxval(up(:,jmingh:jmin-1,:,:)))
             boolb = .true.
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
             amin = min(amin, minval(up(:,jmax+1:jmaxgh,:,:)))
             amax = max(amax, maxval(up(:,jmax+1:jmaxgh,:,:)))
             boolb = .true.
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
             amin = min(amin, minval(up(:,:,kmingh:kmin-1,:)))
             amax = max(amax, maxval(up(:,:,kmingh:kmin-1,:)))
             boolb = .true.
          end if
          if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
             amin = min(amin, minval(up(:,:,kmax+1:kmaxgh,:)))
             amax = max(amax, maxval(up(:,:,kmax+1:kmaxgh,:)))
             boolb = .true.
          end if
       enddo
    enddo
    call mpi_allreduce(boolb, boolg, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    if ( boolg ) then
       call mpi_allreduce(amin, aming, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
       call mpi_allreduce(amax, amaxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       print *, 'minmax boundary', amin,amax
    end if
    if (amin /= 0 .or. amax /= 0) then 
       write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

    endif
  end subroutine fmg_boundary_minmax
end module fmg_boundary
