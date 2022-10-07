  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine mg_diff_resid(mglev, jres, ju, jrhs)
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: res, u, rhs, eta
    integer :: i, j, k, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: hi
    real(kind=DBL_KIND) :: fxl, fxr, fyl, fyr, fzl, fzr
    real(kind=DBL_KIND) :: etaxl, etaxr, etayl, etayr, etazl, etazr
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    call mg_arrp(mglev,jres, res)
    call mg_arrp(mglev,ju, u)
    call mg_arrp(mglev,jrhs, rhs)
    call mg_arrp(mglev,IETA, eta)
    hi = 1.d0/mg_get_h( mglev )
    call mg_boundary_u(mglev,ju)
    do k = ks, ke
       do j = js, je
          do i = is, ie
#ifdef ETA_LOCAL_CONSTANT
             etaxl = eta(i,j,k)
             etaxr = eta(i,j,k)
             etayl = eta(i,j,k)
             etayr = eta(i,j,k)
             etazl = eta(i,j,k)
             etazr = eta(i,j,k)
#else
             etaxl = (eta(i-1,j,k)+eta(i,j,k))*0.5d0
             etaxr = (eta(i+1,j,k)+eta(i,j,k))*0.5d0
             etayl = (eta(i,j-1,k)+eta(i,j,k))*0.5d0
             etayr = (eta(i,j+1,k)+eta(i,j,k))*0.5d0
             etazl = (eta(i,j,k-1)+eta(i,j,k))*0.5d0
             etazr = (eta(i,j,k+1)+eta(i,j,k))*0.5d0
#endif
             fxl = (u(i,j,k) - u(i-1,j,k))*hi*etaxl
             fyl = (u(i,j,k) - u(i,j-1,k))*hi*etayl
             fzl = (u(i,j,k) - u(i,j,k-1))*hi*etazl
             fxr = (u(i+1,j,k) - u(i,j,k))*hi*etaxr
             fyr = (u(i,j+1,k) - u(i,j,k))*hi*etayr
             fzr = (u(i,j,k+1) - u(i,j,k))*hi*etazr
             res(i,j,k)= rhs(i,j,k) + u(i,j,k) - &
                  (fxr-fxl+fyr-fyl+fzr-fzl)*hi
          enddo
       enddo
    enddo
  end subroutine mg_diff_resid
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine mg_diff_relax(mglev, ju, jrhs)
    integer,intent(IN) :: mglev, ju, jrhs
    real(kind=DBL_KIND) :: res, resmax
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, rhs, eta
    integer :: i, j, k, ipass, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: hi, hi2, alpha6
    real(kind=DBL_KIND) :: fxl, fxr, fyl, fyr, fzl, fzr
    real(kind=DBL_KIND) :: etaxl, etaxr, etayl, etayr, etazl, etazr
    hi = 1.d0/mg_get_h( mglev )
    hi2 = hi**2
    call mg_arrp(mglev, ju, u)
    call mg_arrp(mglev, jrhs, rhs)
    call mg_arrp(mglev,IETA, eta)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    ! Red-black Gauss-Seidel iteration
    resmax = 0.d0
    do ipass=1,2              ! 1-red, 2-black
       call mg_boundary_u(mglev,ju)
       do k = ks, ke
          do j = js, je
             do i = is + mod(j+k+ipass,2), ie, 2
#ifdef ETA_LOCAL_CONSTANT
                etaxl = eta(i,j,k)
                etaxr = eta(i,j,k)
                etayl = eta(i,j,k)
                etayr = eta(i,j,k)
                etazl = eta(i,j,k)
                etazr = eta(i,j,k)
#else
                etaxl = (eta(i-1,j,k)+eta(i,j,k))*0.5d0
                etaxr = (eta(i+1,j,k)+eta(i,j,k))*0.5d0
                etayl = (eta(i,j-1,k)+eta(i,j,k))*0.5d0
                etayr = (eta(i,j+1,k)+eta(i,j,k))*0.5d0
                etazl = (eta(i,j,k-1)+eta(i,j,k))*0.5d0
                etazr = (eta(i,j,k+1)+eta(i,j,k))*0.5d0
#endif
                fxl = (u(i,j,k) - u(i-1,j,k))*hi*etaxl
                fyl = (u(i,j,k) - u(i,j-1,k))*hi*etayl
                fzl = (u(i,j,k) - u(i,j,k-1))*hi*etazl
                fxr = (u(i+1,j,k) - u(i,j,k))*hi*etaxr
                fyr = (u(i,j+1,k) - u(i,j,k))*hi*etayr
                fzr = (u(i,j,k+1) - u(i,j,k))*hi*etazr
                alpha6 = (etaxl+etaxr+etayl+etayr+etazl+etazr)*hi2
                res= (rhs(i,j,k) + u(i,j,k) - &
                     (fxr-fxl+fyr-fyl+fzr-fzl)*hi)
                u(i,j,k)=u(i,j,k) - res/(1.d0 + alpha6)
                resmax = max(resmax, abs(res))
             enddo
          enddo
       enddo
    enddo
!!$    print *, mglev, 'resmax', resmax/h2  ! /maxval(rhs)
    Resmaxg = resmax    !! /maxval(rhs(SZ))
!!$    print *, 'Resmaxg', Resmaxg
  end subroutine mg_diff_relax
