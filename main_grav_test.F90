#include "config.h"
! AMRv3
! (c) 2004 Tomoaki Matsumoto <matsu@i.hosei.ac.jp>

program main
  use io
  use refine
!!$  use multi_timestep
  use fmg
  call init
  call restoredata              ! array should be allocated!!

  call cp_rho2rho1
  call fmg_multigrid

  call eval_error_g

!!$  call step_all_level
  call dumpslice(MZ, 0.d0)
  call dumpdata
  call finalize
contains
  !-------------------------------------------------------------------------
  ! Prepare rho
  !-------------------------------------------------------------------------
  subroutine cp_rho2rho1
    use grid
    use fg2cg
    integer :: level, n, gid
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, rho1, psi

!!$    do level = LevelMax, Lmin+1, -1
!!$       call fg2cg_u( level )
!!$    enddo
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          rho => get_Ucomp( MRHO, gid )
          rho1 => get_U1comp( MRHO, gid )
          rho1 = rho
          psi => get_Ucomp( MPSI, gid )
          psi = 0
       enddo
    enddo
  end subroutine cp_rho2rho1
  !-------------------------------------------------------------------------
  ! evaluate error of gravity
  !-------------------------------------------------------------------------
  subroutine eval_error_g
    use mpilib
    use grid
    use grid_boundary
    integer :: level, n, gid
    integer :: i,j,k, ncell
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: psi, gex, gey, gez
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: gx, gy, gz
!!$    real(kind=DBL_KIND) :: gx, gy, gz
    real(kind=DBL_KIND) :: dgx, dgy, dgz, ge_abs, dg
    real(kind=DBL_KIND),dimension(Lmin:Lmax) :: errormax, error_max, error2, error_2, error1, error_1
    myrank = get_myrank()
    errormax(:) = tiny(errormax)
    error2(:) = 0.d0
    error1(:) = 0.d0

    do level = Lmin, LevelMax
       call boundary_grid( level, 0 ) ! boundary fix for emulate of hydro
    enddo
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )

          gid = GidList(n, level) ! gid for U

          psi => get_Ucomp( MPSI, gid )
          psi = 0

          if ( ChildGid(Left, Left, Left, gid, myrank) /= Undefi ) cycle
!!$          psi => get_Ucomp( MPSI, gid )

          gex => get_Ucomp( MVX, gid )
          gey => get_Ucomp( MVY, gid )
          gez => get_Ucomp( MVZ, gid )

          gx => get_Ucomp( MGX, gid )
          gy => get_Ucomp( MGY, gid )
          gz => get_Ucomp( MGZ, gid )
          ncell = (Kmax-Kmin+1)*(Jmax-Jmin+1)*(Imax-Imin+1)*NGI_BASE*NGJ_BASE*NGK_BASE
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
!!$                   gx = -(psi(i+1,j,k)-psi(i-1,j,k))/CellWidth(MX, level)/2.d0
!!$                   gy = -(psi(i,j+1,k)-psi(i,j-1,k))/CellWidth(MY, level)/2.d0
!!$                   gz = -(psi(i,j,k+1)-psi(i,j,k-1))/CellWidth(MZ, level)/2.d0

                   ge_abs = sqrt(gex(i,j,k)**2+gey(i,j,k)**2+gez(i,j,k)**2)

                   dgx = gx(i,j,k) - gex(i,j,k)
                   dgy = gy(i,j,k) - gey(i,j,k)
                   dgz = gz(i,j,k) - gez(i,j,k)
!!$                   dgx = gx - gex(i,j,k)
!!$                   dgy = gy - gey(i,j,k)
!!$                   dgz = gz - gez(i,j,k)
                   dg = sqrt(dgx**2+dgy**2+dgz**2)

                   errormax(level) = max(errormax(level), dg/ge_abs)
                   error2(level) = error2(level) + (dg/ge_abs)**2 /ncell
                   error1(level) = error1(level) + dg/ge_abs /ncell

                   psi(i,j,k) = dg/ge_abs

!!$                   gex(i,j,k) = dg/ge_abs

                enddo
             enddo
          enddo
       enddo
       call mpi_allreduce(errormax(level), error_max(level), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       call mpi_allreduce(error2(level), error_2(level), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       call mpi_allreduce(error1(level), error_1(level), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    enddo

    do level = Lmin, LevelMax
       if (PRIMARY_RANK == myrank) &
            print *, 'errormax', level, error_max(level)
    enddo
    do level = Lmin, LevelMax
       if (PRIMARY_RANK == myrank) &
            print *, 'error1', level, sqrt(error_2(level))
    enddo
    do level = Lmin, LevelMax
       if (PRIMARY_RANK == myrank) &
            print *, 'error2', level, error_1(level)
    enddo

  end subroutine eval_error_g
end program main
!-------------------------------------------------------------------------
! Initialize program
!-------------------------------------------------------------------------
subroutine init
  use mpilib
  use grid
  call mpi_init(ierr)
  call hellomsg
  call grid_init
  call parameter_init
  call eos_init
end subroutine init
!-------------------------------------------------------------------------
subroutine hellomsg
  use systemcall
  use mpilib
  if (get_myrank() == PRIMARY_RANK) then
     call systemcall_command('date')
     print *, 'ZAKU-II (C) Tomoaki Matsumoto 2004,2005,2006'
  endif
  call mpi_barrier( MPI_COMM_WORLD, ierr )
end subroutine hellomsg
!-------------------------------------------------------------------------
! Initialize program
!-------------------------------------------------------------------------
subroutine finalize
  use mpilib
  call mpi_finalize(ierr)
end subroutine finalize
