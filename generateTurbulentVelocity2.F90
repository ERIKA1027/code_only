
!
! generate the initial turbulent velocity field
!
#include "config.h"
!#include "config-gTV.h"
!!$#undef SEED
!
program main
  use mpi
  use string, only : CHARLEN
  use io_util, only : readenv
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
#define SZG 0:NGI-1,0:NGJ-1,0:NGK-1
#define SZL 0:NGI-1,0:NGJ-1,0:NLK-1
#define FH 11
  ! mesh numbers in each direction
  integer,parameter :: NGI=(NI)*(NGI_BASE), NGJ=(NJ)*(NGJ_BASE), NGK=(NK)*(NGK_BASE)
  ! local mesh number in k-direction
  integer,parameter :: NLK=NGK
  ! vector potential (ax, ay, az)
!  complex(KIND=CMPLX_KIND),dimension(:,:,:),allocatable :: ax, ay, az, vx, vy, vz, work
  complex*16,dimension(:,:,:),allocatable :: ax, ay, az
  ! imaginary unit
  complex(KIND=CMPLX_KIND),parameter :: I_ = cmplx(0.d0, 1.d0)
  ! boxsize
  real(KIND=DBL_KIND) :: XMAX, YMAX, ZMAX
  ! k-space
  real(KIND=DBL_KIND),dimension(0:NGI-1) :: kx
  real(KIND=DBL_KIND),dimension(0:NGJ-1) :: ky
  real(KIND=DBL_KIND),dimension(0:NLK-1) :: kz

  real(KIND=DBL_KIND) :: r, phi
  ! sigma
  real(KIND=DBL_KIND) :: sigma
  real(KIND=DBL_KIND),parameter :: sigma0 = 1.d0
  ! power index (sigma_A^2 ~ k^(-n-2),  sigma_v^2 ~ k^(-n) )
  integer,parameter :: index=4

  real(KIND=DBL_KIND) :: PI
  integer :: i, j, k, kg, node, myrank, ierr
  real(KIND=DBL_KIND) :: randAx, randAy, randAz, randPhix, randPhiy, randPhiz
  real(KIND=DBL_KIND),dimension(:,:,:),allocatable :: vxr, vyr, vzr
  real(KIND=DBL_KIND) :: sigma_v_local, sigma_v, MP_Boxsize
  !
  character(len=CHARLEN) :: fn
  !
  type(C_PTR) :: plan, cvx, cvy, cvz
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: vx, vy, vz
  integer(C_INTPTR_T) ::  alloc_local, local_nz, local_k_start

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call fftw_mpi_init()

  !get Boxsize --------------------------------------------
  if(readenv('Boxsize', MP_Boxsize)) then
        print *, 'Boxsize  =', MP_Boxsize, '[au]'             ! (half) boxsize in au
  else
        print *, '***********', 'error in readind Boxsize'
        stop
  endif
  XMAX = MP_Boxsize
  YMAX = MP_Boxsize
  ZMAX = MP_Boxsize
  ! -------------------------------------------------------

  ! get size of array
  alloc_local = fftw_mpi_local_size_3d( NGK, NGJ, NGI, &
       MPI_COMM_WORLD, &
       local_nz, local_k_start)
  !if ( local_nz /= NLK) then
  !   print *, '*** error in MPI data layout', local_nz, NLK, NGK, NPE
  !   stop
  !endif
  ! allocate array
  cvx = fftw_alloc_complex(alloc_local)
  cvy = fftw_alloc_complex(alloc_local)
  cvz = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cvx, vx, [NGI, NGJ, NLK])
  call c_f_pointer(cvy, vy, [NGI, NGJ, NLK])
  call c_f_pointer(cvz, vz, [NGI, NGJ, NLK])

  ! create MPI plan
  plan = fftw_mpi_plan_dft_3d(NGK, NGJ, NGI, vx, vx, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_MEASURE)
  !
  allocate( ax(SZL), ay(SZL), az(SZL) )
  !
  ! define k-space
  !
  do i = 0, NGI/2
     kx(i) = dble(i)/XMAX
  enddo
  do i = NGI/2+1, NGI-1
     kx(i) = (dble(i)-NGI)/XMAX
  enddo

  do j = 0, NGJ/2
     ky(j) = dble(j)/YMAX
  enddo
  do j = NGJ/2+1, NGJ-1
     ky(j) = (dble(j)-NGJ)/YMAX
  end do

  do kg = 0, NGK/2              ! global K
     k = mod(kg,NLK)            ! local k
     node = kg/NLK
     if (myrank == node) kz(k) = dble(kg)/ZMAX
  end do
  do kg = NGK/2+1, NGK-1        ! global K
     k = mod(kg,NLK)            ! local k
     node = kg/NLK
     if (myrank == node) kz(k) = (dble(kg)-NGK)/ZMAX
  end do
  !
  ! vector potential (ax, ay, az) in k-space
  !
  call set_seed                 ! set common seed among nodes

  PI = 4.d0 * atan(1.d0)

  do kg = 0, NGK-1              ! global K
     k = mod(kg,NLK)            ! local k
     node = kg/NLK
     do j = 0, NGJ-1
        do i = 0, NGI-1
           ! produce random number even if this isnt my rank
           call random_number(randAx)
           call random_number(randAy)
           call random_number(randAz)
           call random_number(randPhix)
           call random_number(randPhiy)
           call random_number(randPhiz)
           if (myrank == node) then
              if ( kx(i) == 0 .and. ky(j) == 0 .and. kz(k) == 0 ) then
                 sigma = 0.d0
              else
                 sigma = sigma0/sqrt(kx(i)**2+ky(j)**2+kz(k)**2)**((index+2)/2.d0)
              endif
              r = sigma*sqrt(-2.d0*log(1.d0-randAx))
              phi = randPhix * 2.d0 * PI
              ax(i,j,k) = cmplx(r*cos(phi), r*sin(phi))
              r = sigma*sqrt(-2.d0*log(1.d0-randAy))
              phi = randPhiy * 2.d0 * PI
              ay(i,j,k) = cmplx(r*cos(phi), r*sin(phi))
              r = sigma*sqrt(-2.d0*log(1.d0-randAz))
              phi = randPhiz * 2.d0 * PI
              az(i,j,k) = cmplx(r*cos(phi), r*sin(phi))
           endif
        enddo
     enddo
  enddo
  !
  ! u(k) = i k x A(k)
  !
  do k = 0, NLK-1               ! local k
     do j = 0, NGJ-1
        do i = 0, NGI-1
           vx(i+1,j+1,k+1) = (ky(j) * az(i,j,k) - kz(k) * ay(i,j,k)) * I_
           vy(i+1,j+1,k+1) = (kz(k) * ax(i,j,k) - kx(i) * az(i,j,k)) * I_
           vz(i+1,j+1,k+1) = (kx(i) * ay(i,j,k) - ky(j) * ax(i,j,k)) * I_
        enddo
     enddo
  enddo
  !
  ! FFT u(k) -> v(r)
  !
  call fftw_mpi_execute_dft(plan, vx, vx)
  call fftw_mpi_execute_dft(plan, vy, vy)
  call fftw_mpi_execute_dft(plan, vz, vz)
!!$  vx = vx/(NGI*NGJ*NGK)
!!$  vy = vy/(NGI*NGJ*NGK)
!!$  vz = vz/(NGI*NGJ*NGK)
  !
  ! normalize
  ! v = v/sigma_v,  sigma_v^2 = <v^2>
  !
  allocate( vxr(SZL), vyr(SZL), vzr(SZL) )
  vxr = dble(vx)
  vyr = dble(vy)
  vzr = dble(vz)
  sigma_v_local = sum(vxr**2 + vyr**2 + vzr**2)
  call mpi_allreduce(sigma_v_local, sigma_v, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  sigma_v = sqrt(sigma_v/(NGI*NGJ*NGK))
  vxr = vxr/sigma_v
  vyr = vyr/sigma_v
  vzr = vzr/sigma_v
  call fftw_destroy_plan(plan)
  call fftw_free(cvx)
  call fftw_free(cvy)
  call fftw_free(cvz)
  !
  ! ouput in binary format
  !
  fn = get_filename()
  if (myrank == PRIMARY_RANK) then
     write(*,*) 'outputfile '// trim(fn)
  endif
  do node = 0, (NPE)-1
     call mpi_barrier(MPI_COMM_WORLD, ierr)
     if (myrank /= node ) cycle
     if ( node == 0 ) then
        open(unit=FH, file=fn, form='unformatted', access='stream')
     else
        open(unit=FH, file=fn, form='unformatted', access='stream', position='append')
     endif
     write(FH) vxr
     call flush(FH)
     close(FH)
  enddo
  do node = 0, (NPE)-1
     call mpi_barrier(MPI_COMM_WORLD, ierr)
     if (myrank /= node ) cycle
     open(unit=FH, file=fn, form='unformatted', access='stream', position='append')
     write(FH) vyr
     call flush(FH)
     close(FH)
  enddo
  do node = 0, (NPE)-1
     call mpi_barrier(MPI_COMM_WORLD, ierr)
     if (myrank /= node ) cycle
     open(unit=FH, file=fn, form='unformatted', access='stream', position='append')
     write(FH) vzr
     call flush(FH)
     close(FH)
  enddo

!!$  !
!!$  ! read  for test
!!$  !
!!$  allocate( vxr(SZL), vyr(SZL), vzr(SZL) )
!!$  sz = size(vxr)
!!$  do node = 0, (NNODE)-1
!!$     call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$     if (myrank /= node ) cycle
!!$     open(unit=FH, file=fn, form='unformatted', access='direct', recl=sz*(DBL_KIND) )
!!$     read(FH, rec=1+node) vxr
!!$     read(FH, rec=1+node+(NNODE)) vyr
!!$     read(FH, rec=1+node+(NNODE)*2) vzr
!!$     call flush(FH)
!!$     close(FH)
!!$  enddo
!!$  !
!!$  ! check sum
!!$  !
!!$  if (sum(vxr) == sum(dble(vx))) print *, 'vx OK', myrank
!!$  if (sum(vyr) == sum(dble(vy))) print *, 'vy OK', myrank
!!$  if (sum(vzr) == sum(dble(vz))) print *, 'vz OK', myrank
!!$  !
!!$  ! check all values
!!$  !
!!$  do k = 0, NLK-1               ! local k
!!$     do j = 0, NGJ-1
!!$        do i = 0, NGI-1
!!$           if  ( vxr(i,j,k) /= dble(vx(i,j,k)) ) print *, 'error in vx', i,j,k,myrank
!!$           if  ( vyr(i,j,k) /= dble(vy(i,j,k)) ) print *, 'error in vy', i,j,k,myrank
!!$           if  ( vzr(i,j,k) /= dble(vz(i,j,k)) ) print *, 'error in vz', i,j,k,myrank
!!$        enddo
!!$     enddo
!!$  enddo

  call mpi_finalize(ierr)
!-------------------------------------------------------------------------
contains
  !-------------------------------------------------------------------------
  ! set common seed among the nodes used
  !-------------------------------------------------------------------------
  subroutine set_seed
    integer,dimension(:),allocatable :: seed
    integer :: n
#ifndef SEED
    integer :: clock, i
#endif !SEED
    call random_seed( size=n )
    allocate( seed(n) )

    ! define seed at primary rank
    if ( myrank == PRIMARY_RANK ) then
#ifdef SEED
       seed(:) = SEED
#else ! SEED
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
#endif ! SEED
    end if

    call mpi_bcast(seed, size(seed), MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    call random_seed( put=seed )
    deallocate( seed )
  end subroutine set_seed
  !---------------------------------------------------------------------
  ! get file name
  !---------------------------------------------------------------------
  function get_filename() result(fn)
    use string
    character(len=CHARLEN) :: fn
    character(len=CHARLEN) :: dir
    character(len=2) :: suffix ='.d'
    character(len=CHARLEN) :: seed, size
    call getenv('DIR', dir)
    fn = concat(dir, FOUT)
#ifdef SEED
    seed = num2char(SEED)
    fn = concat(fn, 'Seed'//seed)
#endif !SEED
    size = num2char((NI)*(NGI_BASE))
    fn = concat(fn, 'Size'//size)
    fn = concat(fn, suffix)
  end function get_filename
end program main

! EOF
