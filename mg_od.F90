  !-------------------------------------------------------------------------
  ! get flux
  !-------------------------------------------------------------------------
  subroutine mg_od_flux(mglev, ju)
    use fmg_data, only : MUX, MUY, MUZ
    integer,intent(IN) :: mglev, ju
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, eta
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: i, j, k, is, ie, js, je, ks, ke, n
    real(kind=DBL_KIND) :: etahi, hi
    call mg_boundary_u(mglev,ju)
    call mg_alloc_f( mglev )
    f => mg_get_fp( mglev )
    u => mg_get_arrp(mglev, ju)
    eta => mg_get_arrp(mglev, IETA)
    hi = 1.d0/mg_get_h( mglev )
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    ! flux at i+1/2
    n = MX
    do k = ks, ke
       do j = js, je
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
             etahi = (eta(i,j,k+1,Mmin) + eta(i,j,k,Mmin))*0.5d0 * hi
             f(i,j,k,n,MUX) = (-u(i,j,k+1,MUX)+u(i,j,k,MUX) + (u(i+1,j,k,MUZ)-u(i-1,j,k,MUZ)+u(i+1,j,k+1,MUZ)-u(i-1,j,k+1,MUZ))*0.25d0)*etahi
             f(i,j,k,n,MUY) = (-u(i,j,k+1,MUY)+u(i,j,k,MUY) + (u(i,j+1,k,MUZ)-u(i,j-1,k,MUZ)+u(i,j+1,k+1,MUZ)-u(i,j-1,k+1,MUZ))*0.25d0)*etahi
             f(i,j,k,n,MUZ) = 0.d0
          end do
       end do
    end do
  end subroutine mg_od_flux
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine mg_od_resid(mglev, jres, ju, jrhs)
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: res, u, rhs, eta
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: i, j, k, is, ie, js, je, ks, ke, m
    real(kind=DBL_KIND) :: hi
    real(kind=DBL_KIND) :: fxl, fxr, fyl, fyr, fzl, fzr
    real(kind=DBL_KIND) :: etaxl, etaxr, etayl, etayr, etazl, etazr
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    call mg_od_flux( mglev, ju )
    res => mg_get_arrp(mglev,jres)
    u   => mg_get_arrp(mglev,ju)
    rhs => mg_get_arrp(mglev,jrhs)
    eta => mg_get_arrp(mglev,IETA)
    hi = 1.d0/mg_get_h( mglev )
    f   => mg_get_fp(mglev)
    do m = Mmin, Mmax
       do k = ks, ke
          do j = js, je
             do i = is, ie
                res(i,j,k,m) = rhs(i,j,k,m) - u(i,j,k,m) &
                     -( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                     +  f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                     +  f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA
             end do
          enddo
       enddo
    enddo
  end subroutine mg_od_resid
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine mg_od_relax(mglev, ju, jrhs)
    use fmg_data, only : MUX, MUY, MUZ
    integer,intent(IN) :: mglev, ju, jrhs
    real(kind=DBL_KIND) :: res, resmax
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, rhs, eta
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: i, j, k, ipass, is, ie, js, je, ks, ke, m
    real(kind=DBL_KIND) :: hi, hi2, alpha
    real(kind=DBL_KIND) :: fxl, fxr, fyl, fyr, fzl, fzr
    real(kind=DBL_KIND) :: etaxl, etaxr, etayl, etayr, etazl, etazr
    hi = 1.d0/mg_get_h( mglev )
    hi2 = hi**2
    call mg_alloc_f( mglev )
    u   => mg_get_arrp(mglev, ju)
    rhs => mg_get_arrp(mglev, jrhs)
    eta => mg_get_arrp(mglev,IETA)
    f   => mg_get_fp(mglev)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    ! Red-black Gauss-Seidel iteration
    resmax = 0.d0
    do ipass=1,2              ! 1-red, 2-black
       call mg_od_flux( mglev, ju )
       do m = Mmin, Mmax
          do k = ks, ke
             do j = js, je
                do i = is + mod(j+k+ipass,2), ie, 2
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
    Resmaxg = resmax
!!$#define SZ is:ie,js:je,ks:ke,:
!!$    print *, 'resmax', mglev, resmax, maxval(rhs(SZ))
!!$#undef SZ
!!$    print *, 'Resmaxg', Resmaxg
!!$    print *, ' MG sz', size(u, 1)-Ngh*2
  end subroutine mg_od_relax

  !-------------------------------------------------------------------------
  ! tau correction for nonlinear multigrid
  ! tau = L R Uf - R L Uf
  !-------------------------------------------------------------------------
  subroutine mg_od_tau(mglev, jtau, ju)
    use mg_data
    integer,intent(IN) :: mglev, jtau, ju
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau, u
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: gid, i, j, k, is, ie, js, je, ks, ke, mglevf, amrlev, m
    real(kind=DBL_KIND) :: hi

    mglevf = mglev-1 ! fine level

    if ( mglev > MG_LevelMax ) print *, '*** error: invarid mglev in mg_od_tau'
    if ( mglevf < MG_LevelMin ) print *, '*** error: invarid mglevf in mg_od_tau'

    ! tau = - R L uf
    call mg_get_gridsize(mglevf, is,js,ks,ie,je,ke)
    call mg_od_flux(mglevf, ju) ! in fine grids.
    tau => mg_get_arrp(mglevf, jtau)
    u   => mg_get_arrp(mglevf, ju)
    f   => mg_get_fp  (mglevf)
    hi = 1.d0/mg_get_h(mglevf)
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
    call mg_rstrct(mglev, jtau, jtau)

    ! tau = tau + L R uf = -R L uf + L R uf
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    call mg_od_flux(mglev, ju) ! uf was already restricted (R uf)
    tau => mg_get_arrp(mglev, jtau)
    u   => mg_get_arrp(mglev, ju)
    f   => mg_get_fp  (mglev)
    hi = 1.d0/mg_get_h(mglev)
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
       end do
    end do

  end subroutine mg_od_tau
