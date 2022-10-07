!-------------------------------------------------------------------------
! tau correction
! tau = L R uf - R L uf
!-------------------------------------------------------------------------
subroutine vmg_od_tau(amrlev, jtau, ju)
  use fmg_data
  integer,intent(IN) :: amrlev, jtau, ju
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau, u
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
  integer :: gid, i, j, k, amrlevf, m
  real(kind=DBL_KIND) :: hi
  amrlevf = amrlev+1 ! fine level
  if ( amrlevf > AMR_LevelMax ) return
  myrank = get_myrank()

  ! tau = - R L uf
  call vmg_od_flux(amrlevf, ju) ! in child (fine) grids.
  do gid = fmg_get_gidmin(amrlevf), fmg_get_gidmax(amrlevf)
     tau => fmg_get_arrp(amrlevf, FMG_Level, gid, jtau)
     u   => fmg_get_arrp(amrlevf, FMG_Level, gid, ju)
     f   => fmg_get_fp  (amrlevf, FMG_Level, gid)
     hi = 1.d0/fmg_get_h(amrlevf, FMG_Level)
     tau = 0.d0
     do m = Mmin, Mmax
        do k = Kmin, Kmax
           do j = Jmin, Jmax
              do i = Imin, Imax
                 tau(i,j,k,m) = -( &
                      ( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                      + f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                      +f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m))*hi*LAMBDA &
                      + u(i,j,k,m) )
              end do
           enddo
        enddo
     enddo
  enddo
  call vmg_rstrct(amrlev, jtau)
!!$  print *, 'tau1', maxval(abs(tau(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:))), amrlev

  ! tau = tau + L R uf = - R L uf + L R uf
  call vmg_od_flux(amrlev, ju) ! uf was already restricted (R uf)
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
     if (.not. fmg_have_child(gid, amrlev)) cycle
     tau => fmg_get_arrp(amrlev, FMG_Level, gid, jtau)
     u   => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
     f   => fmg_get_fp  (amrlev, FMG_Level, gid)
     hi = 1.d0/fmg_get_h(amrlev, FMG_Level)
     do m = Mmin, Mmax
        do k = Kmin, Kmax
           do j = Jmin, Jmax
              do i = Imin, Imax
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

!!$  print *, 'tau2', maxval(abs(tau(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:))), amrlev

end subroutine vmg_od_tau
! ----------------------------------------------------------------
! solve flux for one composit grid (for one VMG level)
! ----------------------------------------------------------------
subroutine vmg_od_flux(amrlev, ju)
  use fmg_boundary_phys
  integer,intent(IN) :: amrlev, ju
  real(kind=DBL_KIND) :: hi
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, eta
  integer :: gid, ndir, lr, i, j, k, n
  real(kind=DBL_KIND) :: etahi
  ! --------
  ! 境界条件
  ! --------
  call vmg_ghostcell(amrlev, ju)
  call vmg_boundary_u(amrlev, FMG_Level, ju)
  ! ----------
  ! フラックス
  ! ----------
  hi = 1.d0/fmg_get_h( amrlev, FMG_Level )
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
     f => fmg_get_fp(amrlev, FMG_Level, gid)
     u => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
     eta => fmg_get_arrp(amrlev,FMG_Level, gid, IETA)
     if ( .not. associated( f )) print *, '*** error'
     ! flux at i+1/2
     n = MX
     do k = Kmin, Kmax
        do j = Jmin, Jmax
           do i = Imin-1, Imax
              etahi = (eta(i+1,j,k,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
              f(i,j,k,n,MUX) = 0.d0
              f(i,j,k,n,MUY) = (-u(i+1,j,k,MUY)+u(i,j,k,MUY) + (u(i,j+1,k,MUX)-u(i,j-1,k,MUX)+u(i+1,j+1,k,MUX)-u(i+1,j-1,k,MUX))*0.25d0)*etahi
              f(i,j,k,n,MUZ) = (-u(i+1,j,k,MUZ)+u(i,j,k,MUZ) + (u(i,j,k+1,MUX)-u(i,j,k-1,MUX)+u(i+1,j,k+1,MUX)-u(i+1,j,k-1,MUX))*0.25d0)*etahi
           end do
        end do
     end do
     ! flux at j+1/2
     n = MY
     do k = Kmin, Kmax
        do j = Jmin-1, Jmax
           do i = Imin, Imax
              etahi = (eta(i,j+1,k,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
              f(i,j,k,n,MUX) = (-u(i,j+1,k,MUX)+u(i,j,k,MUX) + (u(i+1,j,k,MUY)-u(i-1,j,k,MUY)+u(i+1,j+1,k,MUY)-u(i-1,j+1,k,MUY))*0.25d0)*etahi
              f(i,j,k,n,MUY) = 0.d0
              f(i,j,k,n,MUZ) = (-u(i,j+1,k,MUZ)+u(i,j,k,MUZ) + (u(i,j,k+1,MUY)-u(i,j,k-1,MUY)+u(i,j+1,k+1,MUY)-u(i,j+1,k-1,MUY))*0.25d0)*etahi
           end do
        end do
     end do
     ! flux at k+1/2
     n = MZ
     do k = Kmin-1, Kmax
        do j = Jmin, Jmax
           do i = Imin, Imax
              etahi = (eta(i,j,k+1,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
              f(i,j,k,n,MUX) = (-u(i,j,k+1,MUX)+u(i,j,k,MUX) + (u(i+1,j,k,MUZ)-u(i-1,j,k,MUZ)+u(i+1,j,k+1,MUZ)-u(i-1,j,k+1,MUZ))*0.25d0)*etahi
              f(i,j,k,n,MUY) = (-u(i,j,k+1,MUY)+u(i,j,k,MUY) + (u(i,j+1,k,MUZ)-u(i,j-1,k,MUZ)+u(i,j+1,k+1,MUZ)-u(i,j-1,k+1,MUZ))*0.25d0)*etahi
              f(i,j,k,n,MUZ) = 0.d0
           end do
        end do
     end do
  end do
end subroutine vmg_od_flux
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine vmg_od_relax(amrlev, ju, jrhs)
    use mpilib
    integer,intent(IN) :: amrlev, ju, jrhs
    real(kind=DBL_KIND) :: res, resmax, resmax_g, hi
    real(kind=DBL_KIND) :: alpha, etaxl,etaxr,etayl,etayr,etazl,etazr
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, rhs, eta
    integer :: gid, ndir, i, j, k, ipass, m

    myrank = get_myrank()
    resmax = 0.d0
    hi = 1.d0/fmg_get_h( amrlev, FMG_Level )
    ! Red-black Gauss-Seidel iteration
    do ipass=1,2              ! 1-red, 2-black
       call vmg_od_flux(amrlev, ju)
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          f   => fmg_get_fp(amrlev, FMG_Level, gid)
          u   => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
          rhs => fmg_get_arrp(amrlev, FMG_Level, gid, jrhs)
          eta => fmg_get_arrp(amrlev, FMG_Level, gid, IETA)
          do m = Mmin, Mmax
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i= Imin + mod(j+k+ipass,2), Imax, 2
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
                      u(i,j,k,m) = u(i,j,k,m) + res/(alpha+1.d0)
                      resmax = max(resmax, abs(res))
                   end do
                enddo
             enddo
          enddo
       enddo
    enddo
    call mpi_allreduce(resmax, resmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    Resmaxg(amrlev) = resmax_g
  end subroutine vmg_od_relax
