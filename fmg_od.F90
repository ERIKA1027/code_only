#include "config.h"
#include "barotropic.h"
  !-------------------------------------------------------------------------
  ! prepare diffusion coefficient eta
  !-------------------------------------------------------------------------
  subroutine fmg_od_prepare_eta_BAK(jeta)
    use unit
    use fmg_data
#ifndef MP
    use eos, only : Rhocr, Gamma, Cs, Kappa
#endif
    integer,intent(IN) :: jeta  ! IN rho, OUT eta
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: eta
    real(kind=DBL_KIND) :: xe, rho, nh2, temp, flag, dt
    integer :: fmglev, amrlev, gid
    fmglev = FMG_LevelMin
    dt = fmg_get_dtime()
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          eta => fmg_get_arrp(amrlev,fmglev,gid,jeta)
!!$          eta = 1.d0
          eta = eta * dt
       enddo
    end do
  end subroutine fmg_od_prepare_eta_BAK
  !-------------------------------------------------------------------------
  ! prepare diffusion coefficient eta
  !-------------------------------------------------------------------------
  subroutine fmg_od_prepare_eta(jeta)
    use unit
!    use grid, only : Dtime
    use fmg_data
#ifndef MP
    use eos, only : Rhocr, Gamma, Cs, Kappa
#endif !MP
    integer,intent(IN) :: jeta  ! IN rho, OUT eta
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: eta
    real(kind=DBL_KIND) :: xe, rho, nh2, temp, flag, dt
    real(kind=DBL_KIND) :: temp_init = ModelParam_temp ! initial temperature in Kelvin
    integer :: fmglev, amrlev, gid, i,j,k,is,js,ks,ie,je,ke
    ! upper limit for density
!    real(kind=DBL_KIND),parameter :: RHOLIMIT = 1.D-9/ModelParam_rho
    real(kind=DBL_KIND),parameter :: RHOLIMIT = HUGE(1.d0)
    fmglev = FMG_LevelMin
    dt = fmg_get_dtime()
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          eta => fmg_get_arrp(amrlev,fmglev,gid,jeta)
          do k = ks, ke
             do j = js, je
                do i= is, ie
                   rho = min(eta(i,j,k,Mmin), RHOLIMIT)
                   rho = abs(rho)     ! for sqrt(small value)
                   nh2 = rho * Unit_n      ! number dinsity in cgs
                   xe = 5.7d-4 / nh2              ! ionozation degree
#ifdef MP
                   temp = p/rho  ! adiabatic temperature
                   Ohmid dissipation for adiabatic EOS is not implemented.
#else !MP
                   GETTEMP( temp, rho )
#endif !MP
                   eta(i,j,k,Mmin) = 740.d0/xe*sqrt(temp/temp_init)*(1.d0-tanh(nh2/1.D15)) ! in cgs [cm^2 s^{-1}]
                   eta(i,j,k,Mmin) = abs(eta(i,j,k,Mmin)) * Unit_t / Unit_l**2 ! in non-dimensional
                end do
             end do
          end do
          eta = eta * dt
       enddo
    end do
  end subroutine fmg_od_prepare_eta
  !-------------------------------------------------------------------------
  ! prepare source term
  !-------------------------------------------------------------------------
  subroutine fmg_od_prepare_source(jsrc)
    use fmg_data
    integer,intent(IN) :: jsrc  ! IN B, OUT SRC
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: src
    real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: f
    real(kind=DBL_KIND) :: hi
    integer:: i, j, k, m, fmglev, amrlev, gid, is,js,ks,ie,je,ke
    fmglev = FMG_LevelMin
    call fmg_alloc_f(fmglev)
    call fmg_od_flux(fmglev, jsrc)
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          src => fmg_get_arrp(amrlev,fmglev,gid,jsrc)
          f => fmg_get_fp(amrlev, fmglev, gid)
          do m = Mmin, Mmax
             do k = ks, ke
                do j = js, je
                   do i= is, ie
                      src(i,j,k,m) = src(i,j,k,m) - &
                           (f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                           +f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                           +f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m)) * hi * (1.d0-LAMBDA)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine fmg_od_prepare_source
  ! ----------------------------------------------------------------
  ! solve flux for one composit grid (for one FMG level)
  ! ----------------------------------------------------------------
  subroutine fmg_od_flux(fmglev, ju, boundary_fill0)
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
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, eta
    real(kind=DBL_KIND) :: etahi
    integer :: amrlev, gid, i, j, k, n, is,js,ks,ie,je,ke
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
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          ! 最細レベルか子供がいなければ,以下に進む
          f => fmg_get_fp(amrlev, fmglev, gid)
          u => fmg_get_arrp(amrlev, fmglev, gid, ju)
          eta => fmg_get_arrp(amrlev,fmglev,gid,IETA)
          ! flux at i+1/2
          n = MX
          do k = ks, ke
             do j = js, je
!VECT
                do i = is-1, ie
                   etahi = (eta(i+1,j,k,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
                   f(i,j,k,n,MUX) = 0.d0
                   f(i,j,k,n,MUY) = (-u(i+1,j,k,MUY)+u(i,j,k,MUY) + (u(i,j+1,k,MUX)-u(i,j-1,k,MUX)+u(i+1,j+1,k,MUX)-u(i+1,j-1,k,MUX))*0.25d0)*etahi
                   f(i,j,k,n,MUZ) = (-u(i+1,j,k,MUZ)+u(i,j,k,MUZ) + (u(i,j,k+1,MUX)-u(i,j,k-1,MUX)+u(i+1,j,k+1,MUX)-u(i+1,j,k-1,MUX))*0.25d0)*etahi
                end do
             end do
          end do
          ! flux at j+1/2
          n = MY
          do k = ks, ke
             do j = js-1, je
!VECT
                do i = is, ie
                   etahi = (eta(i,j+1,k,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
                   f(i,j,k,n,MUX) = (-u(i,j+1,k,MUX)+u(i,j,k,MUX) + (u(i+1,j,k,MUY)-u(i-1,j,k,MUY)+u(i+1,j+1,k,MUY)-u(i-1,j+1,k,MUY))*0.25d0)*etahi
                   f(i,j,k,n,MUY) = 0.d0
                   f(i,j,k,n,MUZ) = (-u(i,j+1,k,MUZ)+u(i,j,k,MUZ) + (u(i,j,k+1,MUY)-u(i,j,k-1,MUY)+u(i,j+1,k+1,MUY)-u(i,j+1,k-1,MUY))*0.25d0)*etahi
                end do
             end do
          end do
          ! flux at k+1/2
          n = MZ
          do k = ks-1, ke
             do j = js, je
                do i = is, ie
!VECT
                   etahi = (eta(i,j,k+1,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
                   f(i,j,k,n,MUX) = (-u(i,j,k+1,MUX)+u(i,j,k,MUX) + (u(i+1,j,k,MUZ)-u(i-1,j,k,MUZ)+u(i+1,j,k+1,MUZ)-u(i-1,j,k+1,MUZ))*0.25d0)*etahi
                   f(i,j,k,n,MUY) = (-u(i,j,k+1,MUY)+u(i,j,k,MUY) + (u(i,j+1,k,MUZ)-u(i,j-1,k,MUZ)+u(i,j+1,k+1,MUZ)-u(i,j-1,k+1,MUZ))*0.25d0)*etahi
                   f(i,j,k,n,MUZ) = 0.d0
                end do
             end do
          end do
       end do
    end do
    ! -------------------------------------
    ! フラックス保存 (子のフラックスを信用)
    ! -------------------------------------
    call fmg_fluxcorrection(fmglev)
  end subroutine fmg_od_flux
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ju = a code of unknown variable
  ! jrhs = a code of right hand side
  ! ----------------------------------------------------------------
  subroutine fmg_od_relax(fmglev, ju, jrhs, resh2maxg, boundary_fill0)
    use mpilib
    use fmg_data
    integer,intent(IN) :: fmglev, ju, jrhs
    real(kind=DBL_KIND),intent(OUT) :: resh2maxg
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND) :: hi, resh2, resh2max, res
    real(kind=DBL_KIND) :: alpha, etaxl,etaxr,etayl,etayr,etazl,etazr
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, rhs, eta
    integer :: amrlev, gid, ndir, i, j, k, m, ipass, is,js,ks,ie,je,ke
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
       call fmg_od_flux(fmglev, ju, bool_boundary_fill0)
!!$       call fmg_boundary_minmax(fmglev, ju) ! debug
       do amrlev = AMR_LevelMin, AMR_LevelMax
          hi = 1.d0/fmg_get_h( amrlev, fmglev )
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             ! 最細レベルか子供がいなければ,以下に進む
             f   => fmg_get_fp(amrlev, fmglev, gid)
             u   => fmg_get_arrp(amrlev, fmglev, gid, ju)
             rhs => fmg_get_arrp(amrlev, fmglev, gid, jrhs)
             eta => fmg_get_arrp(amrlev, fmglev, gid, IETA)
             do m = Mmin, Mmax
                do k = ks, ke
                   do j = js, je
                      do i= is + mod(j+k+ipass,2), ie, 2
                         res = rhs(i,j,k,m) - u(i,j,k,m) &
                              -( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                              +  f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                              +  f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA
                         etaxl = (eta(i-1,j,k,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         etaxr = (eta(i+1,j,k,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         etayl = (eta(i,j-1,k,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         etayr = (eta(i,j+1,k,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         etazl = (eta(i,j,k-1,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         etazr = (eta(i,j,k+1,Mmin)+eta(i,j,k,Mmin))*0.5d0
                         if (m == MUX) then
                            alpha = etayl+etayr+etazl+etazr
                         elseif (m == MUY) then
                            alpha = etazl+etazr+etaxl+etaxr
                         elseif (m == MUZ) then
                            alpha = etaxl+etaxr+etayl+etayr
                         endif
                         alpha = alpha * LAMBDA * hi**2
                         resh2 = res/(alpha+1.d0)
                         u(i,j,k,m) = u(i,j,k,m) + resh2
                         resh2max = max(resh2max,abs(resh2))
                      enddo
                   end do
                enddo
             enddo
          enddo
       enddo
    enddo
    call mpi_allreduce(resh2max, resh2maxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
!!$    if (myrank == PRIMARY_RANK) print *, 'FMG sz', size(u, 1)-Ngh*2
  end subroutine fmg_od_relax
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine fmg_od_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    use fmg_data
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: res, u, rhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: amrlev, gid, i, j, k, m, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: hi
    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    call fmg_od_flux(fmglev, ju, bool_boundary_fill0)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          res => fmg_get_arrp(amrlev,fmglev,gid,jres)
          u   => fmg_get_arrp(amrlev,fmglev,gid,ju)
          rhs => fmg_get_arrp(amrlev,fmglev,gid,jrhs)
          f   => fmg_get_fp(amrlev, fmglev, gid)
          do m = Mmin, Mmax
             do k = ks, ke
                do j = js, je
!VECT
                   do i = is, ie
                      res(i,j,k,m) = rhs(i,j,k,m) - u(i,j,k,m) &
                           - (f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                           +  f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                           +  f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA
                   end do
                enddo
             enddo
          end do
       enddo
    enddo
  end subroutine fmg_od_resid
  !-------------------------------------------------------------------------
  ! tau correction for nonlinear multigrid
  !-------------------------------------------------------------------------
  subroutine fmg_od_tau(fmglev, jtau, ju)
    integer,intent(IN) :: fmglev, jtau, ju
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau, u
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: gid, i, j, k, is, ie, js, je, ks, ke, fmglevf, amrlev, m
    real(kind=DBL_KIND) :: hi

    fmglevf = fmglev-1 ! fine level

    if ( fmglev > FMG_LevelMax ) print *, '*** error: invarid fmglev in fmg_tau'
    if ( fmglevf < FMG_LevelMin ) print *, '*** error: invarid fmglevf in fmg_tau'
    myrank = get_myrank()

    ! tau = - R L uf
    call fmg_get_gridsize(fmglevf, is,js,ks,ie,je,ke)
    call fmg_od_flux(fmglevf, ju) ! in fine grids.
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglevf )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglevf) ) cycle
          tau => fmg_get_arrp(amrlev, fmglevf, gid, jtau)
          u   => fmg_get_arrp(amrlev, fmglevf, gid, ju)
          f   => fmg_get_fp  (amrlev, fmglevf, gid)
          do m = Mmin, Mmax
             do k = ks, ke
                do j = js, je
!VECT
                   do i = is, ie
                      tau(i,j,k,m) = -( &
                           ( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                           + f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                           + f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA &
                           + u(i,j,k,m))
                   end do
                enddo
             enddo
          enddo
       enddo
    end do
    call fmg_rstrct(fmglev, jtau, jtau)

    ! tau = tau + L R uf = -R L uf + L R uf
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    call fmg_od_flux(fmglev, ju) ! uf was already restricted (R uf)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          tau => fmg_get_arrp(amrlev, fmglev, gid, jtau)
          u   => fmg_get_arrp(amrlev, fmglev, gid, ju)
          f   => fmg_get_fp  (amrlev, fmglev, gid)
          do m = Mmin, Mmax
             do k = ks, ke
                do j = js, je
!VECT
                   do i = is, ie
                      tau(i,j,k,m) = tau(i,j,k,m) &
                           +( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                           +  f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                           +  f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA &
                           + u(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       end do
    end do
  end subroutine fmg_od_tau
