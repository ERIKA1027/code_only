#include "debug_fmg.h"
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine mg_poisson_resid(mglev, jres, ju, jrhs)
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: res, u, rhs
    integer :: i, j, k, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: h2i

    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    call mg_arrp(mglev,jres, res)
    call mg_arrp(mglev,ju, u)
    call mg_arrp(mglev,jrhs, rhs)
    h2i = 1.d0/mg_get_h( mglev ) ** 2
    call mg_boundary_u(mglev,ju)
    do k = ks, ke
       do j = js, je
          do i = is, ie
             res(i,j,k)=-h2i*( &
                   u(i+1,j,k)+u(i-1,j,k) &
                  +u(i,j+1,k)+u(i,j-1,k) &
                  +u(i,j,k+1)+u(i,j,k-1) &
                  -6.d0*u(i,j,k))+rhs(i,j,k)
          enddo
       enddo
    enddo
  end subroutine mg_poisson_resid
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ----------------------------------------------------------------
  subroutine mg_poisson_relax(mglev, ju, jrhs)
    integer,intent(IN) :: mglev, ju, jrhs
    real(kind=DBL_KIND),parameter :: sixth = 1.d0/6.d0
!!$    real(kind=DBL_KIND),parameter :: sixth = 1.d0/4.0d0
    real(kind=DBL_KIND) :: h2, resh2, resh2max
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, rhs
    integer :: i, j, k, ipass, is, ie, js, je, ks, ke
#ifdef DEBUG_VMG_OUTPUT_RES
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dbg !debug 最細グリッドの残差を出力する
#endif !DEBUG_VMG_OUTPUT_RES

    h2 = mg_get_h( mglev ) ** 2
    call mg_arrp(mglev, ju, u)
    call mg_arrp(mglev, jrhs, rhs)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
#ifdef DEBUG_VMG_OUTPUT_RES
    if (mglev == MG_LevelMin) then !debug 最細グリッドの残差を出力する
       call mg_alloc_arr(mglev, IDBG)
       call mg_arrp(mglev, IDBG, dbg)
    end if
#endif !DEBUG_VMG_OUTPUT_RES
    ! Red-black Gauss-Seidel iteration
    resh2max = 0.d0
    do ipass=1,2              ! 1-red, 2-black
       call mg_boundary_u(mglev,ju)
       do k = ks, ke
          do j = js, je
             do i = is + mod(j+k+ipass,2), ie, 2
                resh2 = &
                      u(i+1,j,k) + u(i-1,j,k) &
                     +u(i,j+1,k) + u(i,j-1,k) &
                     +u(i,j,k+1) + u(i,j,k-1) &
                     - 6.d0*u(i,j,k) &
                     - h2*rhs(i,j,k)
                u(i,j,k)=u(i,j,k) + sixth * resh2
                resh2max = max(resh2max, abs(resh2))
#ifdef DEBUG_VMG_OUTPUT_RES
                if (mglev == MG_LevelMin) dbg(i,j,k) = resh2/h2 !debug 最細グリッドの残差を出力する
#endif !DEBUG_VMG_OUTPUT_RES
             enddo
          enddo
       enddo
    enddo
!!$    print *, mglev, 'resmax', resh2max/h2  ! /maxval(rhs)
#define SZ is:ie,js:je,ks:ke
    Resmaxg = resh2max/h2    !! /maxval(rhs(SZ))
#undef SZ
  end subroutine mg_poisson_relax
