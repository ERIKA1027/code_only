#include "config.h"
  !-------------------------------------------------------------------------
  ! prepare diffusion coefficient eta
  !-------------------------------------------------------------------------
  subroutine fmg_diff_prepare_eta(jeta)
    use fmg_data
    integer,intent(IN) :: jeta  ! IN rho, OUT eta
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: eta
    real(kind=DBL_KIND) :: dt
    integer :: fmglev, amrlev, gid
    fmglev = FMG_LevelMin
    dt = fmg_get_dtime()
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          call fmg_arrp(amrlev,fmglev,gid,jeta, eta)
          eta = sqrt(eta)       ! rho -> eta
!!$          eta = 1.d0
          eta = eta * dt * LAMBDA
       enddo
    end do
  end subroutine fmg_diff_prepare_eta
  !-------------------------------------------------------------------------
  ! prepare source term
  !-------------------------------------------------------------------------
  subroutine fmg_diff_prepare_source(jsrc, jeta)
    use fmg_data
    integer,intent(IN) :: jsrc  ! IN B, OUT S = -b -dt/2 eta ( grad B )
    integer,intent(IN) :: jeta  ! IN eta
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: src, eta
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: f
    real(kind=DBL_KIND) :: h, hi
    integer:: i, j, k, fmglev, amrlev, gid, is,js,ks,ie,je,ke
    fmglev = FMG_LevelMin
    call fmg_alloc_f(fmglev)
    call fmg_diff_flux(fmglev, jsrc) ! grad B
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          call fmg_arrp(amrlev,fmglev,gid,jsrc, src)
          call fmg_arrp(amrlev,fmglev,gid,jeta, eta)
          call fmg_fp(amrlev, fmglev, gid, f)
          h = fmg_get_h( amrlev, fmglev )
          hi = 1.d0/h
          do k = ks, ke
             do j = js, je
                do i= is, ie
                   src(i,j,k) = -src(i,j,k) - &
                        (f(i,j,k,MX)-f(i-1,j,k,MX) &
                        +f(i,j,k,MY)-f(i,j-1,k,MY) &
                        +f(i,j,k,MZ)-f(i,j,k-1,MZ)) * hi * (1.d0-LAMBDA)/LAMBDA
                end do
             end do
          end do
       end do
    end do
  end subroutine fmg_diff_prepare_source
  ! ----------------------------------------------------------------
  ! solve flux for one composit grid (for one FMG level)
  ! ----------------------------------------------------------------
  subroutine fmg_diff_flux(fmglev, ju, boundary_fill0)
    use fmg_data
    use fmg_boundary_phys, only : fmg_boundary_u
    use fmg_converge
    use fmg_reflux
    use fmg_ghostcell
    use fmg_boundary
    integer,intent(IN) :: fmglev, ju
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND) :: hi
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, eta
    real(kind=DBL_KIND) :: etaX, etaY, etaZ
    integer :: amrlev, gid, ndir, lr, i, j, k, imin,jmin,kmin,imax,jmax,kmax
    ! initial guess がゼロであるときは、境界の転送・設定は不要
    if (FmgLevel_fill0 /= fmglev) then
       ! ------------
       ! 袖の値の設定
       ! ------------
       call fmg_converge_c2p(fmglev,ju)
       call fmg_boundary_u(fmglev, ju)
       call fmg_ghostcell_fix(fmglev,ju)
       ! ------------------------------
       ! 外側の境界値にゼロを入れるか？
       ! ------------------------------
       if (present(boundary_fill0)) then
          bool_boundary_fill0 = boundary_fill0
       else
          bool_boundary_fill0 = .false.
       endif
       if (bool_boundary_fill0) call fmg_boundary_fill0(fmglev, ju)
       call fmg_boundary_u(fmglev, ju)
    endif
    FmgLevel_fill0 = Undefi
    ! ----------
    ! フラックス
    ! ----------
    myrank = get_myrank()
    call fmg_get_gridsize(fmglev, imin,jmin,kmin,imax,jmax,kmax)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          ! 最細レベルか子供がいなければ,以下に進む
          call fmg_fp(amrlev, fmglev, gid, f)
          call fmg_arrp(amrlev, fmglev, gid, ju, u)
          call fmg_arrp(amrlev,fmglev,gid,IETA, eta)
          do k = kmin-1, kmax
             do j = jmin-1, jmax
!VECT
                do i = imin-1, imax
                   etaX = (eta(i+1,j,k)+eta(i,j,k))*0.5d0
                   etaY = (eta(i,j+1,k)+eta(i,j,k))*0.5d0
                   etaZ = (eta(i,j,k+1)+eta(i,j,k))*0.5d0
                   f(i,j,k,MX)= (u(i+1,j,k)-u(i,j,k))*hi*etaX
                   f(i,j,k,MY)= (u(i,j+1,k)-u(i,j,k))*hi*etaY
                   f(i,j,k,MZ)= (u(i,j,k+1)-u(i,j,k))*hi*etaZ
                end do
             end do
          end do
       end do
    end do
    ! -------------------------------------
    ! フラックス保存 (子のフラックスを信用)
    ! -------------------------------------
    call fmg_fluxcorrection(fmglev)
  end subroutine fmg_diff_flux
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ju = a code of unknown variable
  ! jrhs = a code of right hand side
  ! ----------------------------------------------------------------
  subroutine fmg_diff_relax(fmglev, ju, jrhs, resh2maxg, boundary_fill0)
    use mpilib
    use fmg_data
    integer,intent(IN) :: fmglev, ju, jrhs
    real(kind=DBL_KIND),intent(OUT) :: resh2maxg
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND) :: h, hi, resh2, resh2max, alpha6
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, rhs, eta
    integer :: amrlev, gid, ndir, i, j, k, ipass, is,js,ks,ie,je,ke

    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif

    myrank = get_myrank()
    resh2max = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    ! Red-black Gauss-Seidel iteration
    do ipass=1,2              ! 1-red, 2-black
       call fmg_diff_flux(fmglev, ju, bool_boundary_fill0)
!!$       call fmg_boundary_minmax(fmglev, ju) ! debug
       do amrlev = AMR_LevelMin, AMR_LevelMax
          h = fmg_get_h( amrlev, fmglev )
          hi = 1.d0/h
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             ! 最細レベルか子供がいなければ,以下に進む
             call fmg_fp(amrlev, fmglev, gid, f)
             call fmg_arrp(amrlev, fmglev, gid, ju, u)
             call fmg_arrp(amrlev, fmglev, gid, jrhs, rhs)
             call fmg_arrp(amrlev, fmglev, gid, IETA, eta)
             do k = ks, ke
                do j = js, je
                   do i= is + mod(j+k+ipass,2), ie, 2
                      alpha6 = (eta(i+1,j,k) + eta(i-1,j,k) + eta(i,j+1,k) + eta(i,j-1,k) + eta(i,j,k+1) + eta(i,j,k-1) + eta(i,j,k) * 6.d0)*0.5d0*hi**2
                      resh2 = ( &
                           rhs(i,j,k) + u(i,j,k) - &
                           (f(i,j,k,MX)-f(i-1,j,k,MX) &
                           +f(i,j,k,MY)-f(i,j-1,k,MY) &
                           +f(i,j,k,MZ)-f(i,j,k-1,MZ))*hi )/(1.d0+ alpha6)
                      u(i,j,k) = u(i,j,k) - resh2
                      resh2max = max(resh2max,abs(resh2))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    call mpi_allreduce(resh2max, resh2maxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
!!$    print *, resh2maxg
  end subroutine fmg_diff_relax
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine fmg_diff_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    use fmg_data
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: res, u, rhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    integer :: amrlev, gid, i, j, k, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: hi
    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    call fmg_diff_flux(fmglev, ju, bool_boundary_fill0)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          call fmg_arrp(amrlev,fmglev,gid,jres, res)
          call fmg_arrp(amrlev,fmglev,gid,ju, u)
          call fmg_arrp(amrlev,fmglev,gid,jrhs, rhs)
          call fmg_fp(amrlev, fmglev, gid, f)
          hi = 1.d0/fmg_get_h( amrlev, fmglev )
          do k = ks, ke
             do j = js, je
                do i = is, ie
                   res(i,j,k) = &
                        rhs(i,j,k) + u(i,j,k) - &
                        (f(i,j,k,MX)-f(i-1,j,k,MX) &
                        +f(i,j,k,MY)-f(i,j-1,k,MY) &
                        +f(i,j,k,MZ)-f(i,j,k-1,MZ))*hi
                end do
             enddo
          enddo
       enddo
    enddo
  end subroutine fmg_diff_resid

