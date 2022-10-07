#include "config.h"
#include "debug_fmg.h"
  !-------------------------------------------------------------------------
  subroutine fmg_psi2g(ju)
    use grid, only : Gidmin, GidListMax, GidList, Imin, Imax, Jmin, Jmax, Kmin, Kmax, get_Ucomp
    use fmg_data, only : FMG_LevelMin, AMR_LevelMin, AMR_LevelMax, fmg_get_h, fmg_arrp, fmg_fp, fmg_skip_grid, Mmin
    integer,intent(IN) :: ju
    integer :: amrlev, fmglev, i, j, k, n, gid
    real(kind=DBL_KIND) :: hi
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: gx, gy, gz, u

    fmglev = FMG_LevelMin
    call fmg_poisson_flux(fmglev, ju)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) ! grid ID in AMR data
          call fmg_arrp(amrlev, fmglev, n, ju, u)
          call fmg_fp(amrlev, fmglev, n, f)
          gx => get_Ucomp(MGX, gid)
          gy => get_Ucomp(MGY, gid)
          gz => get_Ucomp(MGZ, gid)
          if ( fmg_skip_grid(n, amrlev, fmglev) ) then
             ! overwrapped region
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i = Imin, Imax
                      gx(i,j,k) = -(u(i+1,j,k) - u(i-1,j,k))*hi*0.5d0
                      gy(i,j,k) = -(u(i,j+1,k) - u(i,j-1,k))*hi*0.5d0
                      gz(i,j,k) = -(u(i,j,k+1) - u(i,j,k-1))*hi*0.5d0
                   enddo
                enddo
             enddo
          else
             ! non overwrapped region
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i = Imin, Imax
                      gx(i,j,k) = -(f(i,j,k,MX)+f(i-1,j,k,MX))*0.5d0
                      gy(i,j,k) = -(f(i,j,k,MY)+f(i,j-1,k,MY))*0.5d0
                      gz(i,j,k) = -(f(i,j,k,MZ)+f(i,j,k-1,MZ))*0.5d0
                   enddo
                enddo
             enddo
          endif
       enddo
    enddo
  end subroutine fmg_psi2g
  ! ----------------------------------------------------------------
  ! solve flux for one composit grid (for one FMG level)
  ! ----------------------------------------------------------------
  subroutine fmg_poisson_flux(fmglev, ju, boundary_fill0)
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
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u
    integer :: amrlev, gid, ndir, lr, i, j, k, imin,jmin,kmin,imax,jmax,kmax
    !KS DEBUG
    ! if(get_myrank() == 231) &
    !      print '(/,A,2I4)', "(Begining of fmg_poisson_flux) myrank, fmglev = ", get_myrank(), fmglev
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
          do k = kmin-1, kmax
             do j = jmin-1, jmax
!VECT
                do i = imin-1, imax
                   f(i,j,k,MX)= (u(i+1,j,k)-u(i,j,k))*hi
                   f(i,j,k,MY)= (u(i,j+1,k)-u(i,j,k))*hi
                   f(i,j,k,MZ)= (u(i,j,k+1)-u(i,j,k))*hi
                      !--------------------------------KS DEBUG ---------------------------------!
!                       if (gid == 0 .and. amrlev==13 .and. i==1 .and. j==1 .and. k==1 .and. myrank==231) then
!                          print '(A,7I5, 1P8E15.7)',"(FMG_POISSON_FLUX, KS DEBUG) ",myrank, fmglev, amrlev, gid, i,j,k,&
! f(i,j,k,MX), f(i,j,k,MY), f(i,j,k,MZ), u(i+1,j,k),u(i,j+1,k),u(i,j,k+1), u(i,j,k), hi
!                       end if
                      !--------------------------------KS DEBUG ---------------------------------!                                    
                end do
             end do
          end do
       end do
    end do
    ! -------------------------------------
    ! フラックス保存 (子のフラックスを信用)
    ! -------------------------------------
    call fmg_fluxcorrection(fmglev)
  end subroutine fmg_poisson_flux
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ju = a code of unknown variable
  ! jrhs = a code of right hand side
  ! ----------------------------------------------------------------
  subroutine fmg_poisson_relax(fmglev, ju, jrhs, resh2maxg, boundary_fill0)
    use mpilib
    use fmg_data
    integer,intent(IN) :: fmglev, ju, jrhs
    real(kind=DBL_KIND),intent(OUT) :: resh2maxg
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=DBL_KIND),parameter :: sixth = 1.d0/6.0d0
    real(kind=DBL_KIND) :: h, dv, h2idv, resh2, resh2max, resh2max_g, ressum
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: u, rhs
    integer :: amrlev, gid, ndir, i, j, k, ipass, is,js,ks,ie,je,ke
#ifdef DEBUG_FMG_OUTPUT_RES
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dbg ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
    !KS DEBUG
    ! if(get_myrank() == 231) &
    !      print '(/,A,2I4)', "(Begining of fmg_poisson_relax) myrank, fmglev = ", get_myrank(), fmglev

    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif

    myrank = get_myrank()
    resh2max = 0.d0
    ressum = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    ! Red-black Gauss-Seidel iteration
    do ipass=1,2              ! 1-red, 2-black
       !KS DEBUG
       ! if(get_myrank() == 231) &
       !      print '(/,A,3I4)', "(in fmg_poisson_relax) myrank, fmglev, ipass = ", get_myrank(), fmglev,  ipass

       call fmg_poisson_flux(fmglev, ju, bool_boundary_fill0)
!!$       call fmg_boundary_minmax(fmglev, ju) ! debug
       do amrlev = AMR_LevelMin, AMR_LevelMax
          h = fmg_get_h( amrlev, fmglev )
          dv = fmg_get_dv( amrlev, fmglev )
          h2idv = dv/h**2
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             ! 最細レベルか子供がいなければ,以下に進む
             call fmg_fp(amrlev, fmglev, gid, f)
             call fmg_arrp(amrlev, fmglev, gid, ju, u)
             call fmg_arrp(amrlev, fmglev, gid, jrhs, rhs)
#ifdef DEBUG_FMG_OUTPUT_RES
             if (fmglev == FMG_LevelMin) call fmg_arrp(amrlev, fmglev, gid, IDBG, dbg) ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
             do k = ks, ke
                do j = js, je
                   do i= is + mod(j+k+ipass,2), ie, 2
                      resh2 = &
                           (f(i,j,k,MX)-f(i-1,j,k,MX) &
                           +f(i,j,k,MY)-f(i,j-1,k,MY) &
                           +f(i,j,k,MZ)-f(i,j,k-1,MZ) &
                           -rhs(i,j,k) * h &
                           )* h
                      u(i,j,k) = u(i,j,k) + sixth * resh2
                      resh2max = max(resh2max,abs(resh2))
                      !--------------------------------KS DEBUG ---------------------------------!
                      ! if (gid == 0 .and. amrlev==13 .and. i==1 .and. j==1 .and. k==1 .and. myrank==231) then
                      ! ! if (gid == 0 .and. i==1 .and. j==1 .and. k==1 .and. (myrank==231 .or. myrank==232)) then
                      ! ! if (resh2 /= resh2) then
                      ! ! if (abs(resh2) > 1d0 .or. resh2 /= resh2) then
                      !    ! if ((i==1 .and. j==1 .and. k==1 .and. myrank==231) .or. resh2 /= resh2) then

                      !    print '(A,7I5, 1P9E15.7)',"(FMG_POISSON_RELAX, KS DEBUG) ",myrank, fmglev, amrlev, gid, i,j,k,&
                      !         resh2, f(i,j,k,MX),f(i-1,j,k,MX), f(i,j,k,MY),f(i,j-1,k,MY),f(i,j,k,MZ),f(i,j,k-1,MZ),rhs(i,j,k),h
                      ! end if
                      !--------------------------------KS DEBUG ---------------------------------!                      
#ifdef DEBUG_FMG_OUTPUT_RES
                      if (fmglev == FMG_LevelMin) dbg(i,j,k) = resh2/h**2 ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
!!$                      ressum = ressum + abs(resh2) * h2idv
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    call mpi_allreduce(resh2max, Resh2maxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  end subroutine fmg_poisson_relax
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine fmg_poisson_resid(fmglev, jres, ju, jrhs, boundary_fill0)
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
    call fmg_poisson_flux(fmglev, ju, bool_boundary_fill0)
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
!VECT
                do i = is, ie
                   res(i,j,k) = &
                        rhs(i,j,k) &
                       -(f(i,j,k,MX)-f(i-1,j,k,MX) &
                        +f(i,j,k,MY)-f(i,j-1,k,MY) &
                        +f(i,j,k,MZ)-f(i,j,k-1,MZ))*hi
                   !--------------------------------KS DEBUG ---------------------------------!
                   ! if (gid == 0 .and. amrlev==13 .and. i==1 .and. j==1 .and. k==1 .and. myrank==231) then
                   !    print '(A,7I5, 1P9E15.7)',"(FMG_POISSON_RESID, KS DEBUG) ",myrank, fmglev, amrlev, gid, i,j,k,&
                   !         res(i,j,k), rhs(i,j,k), &
                   !         f(i,j,k,MX), f(i-1,j,k,MX), f(i,j,k,MY), f(i,j-1,k,MY), f(i,j,k,MZ), f(i,j,k-1,MZ), hi
                   !    end if
                      !--------------------------------KS DEBUG ---------------------------------!                                    
                   
                end do
             enddo
          enddo
       enddo
    enddo

  end subroutine fmg_poisson_resid

