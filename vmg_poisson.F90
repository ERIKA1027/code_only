#include "debug_fmg.h"
!-------------------------------------------------------------------------
! tau correction
! tau = L R uf - R L uf
!-------------------------------------------------------------------------
subroutine vmg_poisson_tau(amrlev, jtau, ju)
  use fmg_data
  integer,intent(IN) :: amrlev, jtau, ju
  real(kind=DBL_KIND),pointer,dimension(:,:,:) :: tau, u
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
  integer :: gid, i, j, k, amrlevf, m
  real(kind=DBL_KIND) :: hi
  amrlevf = amrlev+1 ! fine level
  if ( amrlevf > AMR_LevelMax ) return
  myrank = get_myrank()

  ! tau = - R L uf
  call vmg_poisson_flux(amrlevf, ju) ! in child (fine) grids.
  do gid = fmg_get_gidmin(amrlevf), fmg_get_gidmax(amrlevf)
     call fmg_arrp(amrlevf, FMG_Level, gid, jtau, tau)
     call fmg_arrp(amrlevf, FMG_Level, gid, ju,   u)
     call fmg_fp  (amrlevf, FMG_Level, gid, f)
     hi = 1.d0/fmg_get_h(amrlevf, FMG_Level)
     tau = 0.d0
     do k = Kmin, Kmax
        do j = Jmin, Jmax
           do i = Imin, Imax
              tau(i,j,k) = &
                   -(f(i,j,k,MX)-f(i-1,j,k,MX) &
                   + f(i,j,k,MY)-f(i,j-1,k,MY) &
                   + f(i,j,k,MZ)-f(i,j,k-1,MZ))*hi
           end do
        enddo
     enddo
  enddo
  call vmg_rstrct(amrlev, jtau)
!!$  print *, 'tau1', maxval(abs(tau(Imin:Imax,Jmin:Jmax,Kmin:Kmax))), amrlev

  ! tau = tau + L R uf = - R L uf + L R uf
  call vmg_poisson_flux(amrlev, ju) ! uf was already restricted (R uf)
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
     if (.not. fmg_have_child(gid, amrlev)) cycle
     call fmg_arrp(amrlev, FMG_Level, gid, jtau, tau)
     call fmg_arrp(amrlev, FMG_Level, gid, ju,   u)
     call fmg_fp  (amrlev, FMG_Level, gid, f)
     hi = 1.d0/fmg_get_h(amrlev, FMG_Level)
     do k = Kmin, Kmax
        do j = Jmin, Jmax
           do i = Imin, Imax
              tau(i,j,k) = tau(i,j,k) &
                   +(f(i,j,k,MX)-f(i-1,j,k,MX) &
                   + f(i,j,k,MY)-f(i,j-1,k,MY) &
                   + f(i,j,k,MZ)-f(i,j,k-1,MZ))*hi
           enddo
        enddo
     enddo
  end do
!!$  print *, 'tau2', maxval(abs(tau(Imin:Imax,Jmin:Jmax,Kmin:Kmax))), amrlev
end subroutine vmg_poisson_tau
! ----------------------------------------------------------------
! solve flux for one composit grid (for one VMG level)
! ----------------------------------------------------------------
subroutine vmg_poisson_flux(amrlev, ju)
  use fmg_boundary_phys
    integer,intent(IN) :: amrlev, ju
    real(kind=DBL_KIND) :: hi
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u
    integer :: gid, ndir, lr, i, j, k
    ! --------
    ! 境界条件
    ! --------
    call vmg_ghostcell(amrlev, ju)
    call vmg_boundary_u(amrlev, FMG_Level, ju)
    ! ----------
    ! フラックス
    ! ----------
    hi = 1.d0/fmg_get_h( amrlev, FMG_Level )
!OCL NOVREC
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       call fmg_fp(amrlev, FMG_Level, gid, f)
       call fmg_arrp(amrlev, FMG_Level, gid, ju, u)
       do k = Kmin-1, Kmax
          do j = Jmin-1, Jmax
!VECT
             do i = Imin-1, Imax
                f(i,j,k,MX)= (u(i+1,j,k)-u(i,j,k))*hi
                f(i,j,k,MY)= (u(i,j+1,k)-u(i,j,k))*hi
                f(i,j,k,MZ)= (u(i,j,k+1)-u(i,j,k))*hi
             end do
          end do
       end do
    end do
end subroutine vmg_poisson_flux
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine vmg_poisson_relax(amrlev, ju, jrhs)
    use mpilib
    integer,intent(IN) :: amrlev, ju, jrhs
    real(kind=DBL_KIND),parameter :: sixth = 1.d0/6.0d0
    real(kind=DBL_KIND) :: h, resh2, resmax, resmax_g
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, rhs
#ifdef DEBUG_VMG_OUTPUT_RES
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dbg ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
    integer :: gid, ndir, i, j, k, ipass
    myrank = get_myrank()
    resmax = 0.d0
    h = fmg_get_h( amrlev, FMG_Level )
    ! Red-black Gauss-Seidel iteration
!OCL NOVREC
    do ipass=1,2              ! 1-red, 2-black
       call vmg_poisson_flux(amrlev, ju)
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          call fmg_fp(amrlev, FMG_Level, gid, f)
          call fmg_arrp(amrlev, FMG_Level, gid, ju, u)
          call fmg_arrp(amrlev, FMG_Level, gid, jrhs, rhs)
#ifdef DEBUG_VMG_OUTPUT_RES
          call fmg_arrp(amrlev, FMG_Level, gid, IDBG, dbg) ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i= Imin + mod(j+k+ipass,2), Imax, 2
                   resh2 = &
                        (f(i,j,k,MX)-f(i-1,j,k,MX) &
                        +f(i,j,k,MY)-f(i,j-1,k,MY) &
                        +f(i,j,k,MZ)-f(i,j,k-1,MZ) &
                        -rhs(i,j,k) * h &
                        )* h
                   u(i,j,k) = u(i,j,k) + sixth * resh2
                   resmax = max(resmax, abs(resh2))
#ifdef DEBUG_VMG_OUTPUT_RES
                   dbg(i,j,k) = resh2/h**2 ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
                enddo
             enddo
          enddo
       enddo
    enddo
    resmax = resmax / h**2
    call mpi_allreduce(resmax, resmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    Resmaxg(amrlev) = resmax_g
  end subroutine vmg_poisson_relax
