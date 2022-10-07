#include "config.h"
#include "barotropic.h"
#include "debug_fmg.h"
  !-------------------------------------------------------------------------
  ! ambipolar diffusion
  !-------------------------------------------------------------------------
  subroutine fmg_ad_prepare_diffcoeff(jdod, jdhe, jdad)
    use unit
    use fmg_data
#ifndef MP
    use eos, only : Rhocr, Gamma, Cs, Kappa
#endif !MP
    integer,intent(IN) :: jdod, jdhe, jdad  ! IN rho, OUT diffusion coefficient
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: dod, dhe, dad
    real(kind=DBL_KIND) :: xe, rho, nh2, temp, flag, dt
    integer :: fmglev, amrlev, gid
    fmglev = FMG_LevelMin
    dt = fmg_get_dtime()
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          dod => fmg_get_arrp(amrlev,fmglev,gid,jdod)
          dhe => fmg_get_arrp(amrlev,fmglev,gid,jdhe)
          dad => fmg_get_arrp(amrlev,fmglev,gid,jdad)
          dod = dod * dt
          dhe = dhe * dt
          dad = dad * dt
       enddo
    end do
  end subroutine fmg_ad_prepare_diffcoeff
  !-------------------------------------------------------------------------
  ! prepare diffusion coefficient eta
  !-------------------------------------------------------------------------
!   subroutine fmg_ad_prepare_eta(jeta)
!     use unit
!     use fmg_data
!     use eos, only : Rhocr, Gamma, Cs, Kappa
!     integer,intent(IN) :: jeta  ! IN rho, OUT eta
!     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: eta
!     real(kind=DBL_KIND) :: xe, rho, nh2, temp, flag, dt
!     real(kind=DBL_KIND) :: temp_init = ModelParam_temp ! initial temperature in Kelvin
!     integer :: fmglev, amrlev, gid, i,j,k,is,js,ks,ie,je,ke
! !!$    real(kind=DBL_KIND) :: rhomax, etamax

!     fmglev = FMG_LevelMin
!     dt = fmg_get_dtime()
! !!$    rhomax = 0.d0
! !!$    etamax = 0.d0
!     call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
!     do amrlev = AMR_LevelMin, AMR_LevelMax
!        do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
!           eta => fmg_get_arrp(amrlev,fmglev,gid,jeta)
! !!$          rhomax = max(rhomax, maxval(eta(is:ie,js:je,ks:ke,Mmin)))
!           do k = ks, ke
!              do j = js, je
!                 do i= is, ie
!                    nh2 = eta(i,j,k,Mmin) * Unit_n      ! number dinsity in cgs
!                    xe = 5.7d-4 / nh2              ! ionozation degree
!                    GETTEMP( temp, eta(i,j,k,Mmin) )
!                    eta(i,j,k,Mmin) = 740.d0/xe*sqrt(temp/temp_init)*(1.d0-tanh(nh2/1.D15)) ! in cgs [cm^2 s^{-1}]
!                    eta(i,j,k,Mmin) = eta(i,j,k,Mmin) * Unit_t / Unit_l**2 ! in non-dimensional
! !!$                   eta(i,j,k,Mmin) = 1.d0
!                 end do
!              end do
!           end do
! !!$          etamax = max(etamax, maxval(eta(is:ie,js:je,ks:ke,Mmin)))
!           eta = eta * dt
!        enddo
!     end do
! !!$    print "('rhomax etamax =', 1P2E12.5)", rhomax, etamax
! !!$    print "('rhomax etamax =', 1P2E12.5)", rhomax*Unit_n, etamax*Unit_l**2/Unit_t
!   end subroutine fmg_ad_prepare_eta
  !-------------------------------------------------------------------------
  ! prepare source term
  !-------------------------------------------------------------------------
  subroutine fmg_ad_prepare_source(jsrc)
    use fmg_data
    integer,intent(IN) :: jsrc  ! IN B, OUT SRC
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: src
    real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: f
    real(kind=DBL_KIND) :: hi
    integer:: i, j, k, m, fmglev, amrlev, gid, is,js,ks,ie,je,ke
    fmglev = FMG_LevelMin
    call fmg_alloc_f(fmglev)
    call fmg_ad_flux(fmglev, jsrc)
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          src => fmg_get_arrp(amrlev,fmglev,gid,jsrc)
          f => fmg_get_fp(amrlev, fmglev, gid)
          do m = Mmin, Mmax
             do k = ks, ke
                do j = js, je
!VECT
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
  end subroutine fmg_ad_prepare_source
  ! ----------------------------------------------------------------
  ! solve flux for one composit grid (for one FMG level)
  ! ----------------------------------------------------------------
  subroutine fmg_ad_flux(fmglev, ju)
    use fmg_data
    use fmg_boundary_phys, only : fmg_boundary_u
    use fmg_converge
    use fmg_reflux
    use fmg_ghostcell
    use fmg_boundary
    integer,intent(IN) :: fmglev, ju
    real(kind=DBL_KIND) :: hi, h2i
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dod, dhe, dad
    integer :: n, nx, ny, nz
    integer,parameter :: DIMOFFSET = 1-MX
    integer :: amrlev, gid
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: bh, rotBh
    real(kind=DBL_KIND),dimension(:,:,:),allocatable :: dodh, dheh, dadh
    ! ------------
    ! 袖の値の設定
    ! ------------
    call fmg_converge_c2p(fmglev,ju)
    call fmg_boundary_u(fmglev, ju)
    call fmg_ghostcell_fix(fmglev,ju, cubic=.TRUE.)
    ! ----------
    ! フラックス
    ! ----------
    myrank = get_myrank()
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       h2i = hi * 0.5d0
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          ! 最細レベルか子供がいなければ,以下に進む
          call fmg_fp  (amrlev,fmglev,gid, f)
          call fmg_arrp(amrlev,fmglev,gid,ju,   u)
          call fmg_arrp(amrlev,fmglev,gid,IDOD, dod)
          call fmg_arrp(amrlev,fmglev,gid,IDHE, dhe)
          call fmg_arrp(amrlev,fmglev,gid,IDAD, dad)
          allocate( bh(ARRAYSIZE4(u)), rotBh(ARRAYSIZE4(u)), dodh(ARRAYSIZE3(u)), dheh(ARRAYSIZE3(u)), dadh(ARRAYSIZE3(u)))

#define SHFTR( A, NDIM ) cshift((A),  1, (NDIM) + DIMOFFSET)
#define SHFTL( A, NDIM ) cshift((A), -1, (NDIM) + DIMOFFSET)
#define bc u
          do n = MX, MZ
             call cyclecomp( n, nx, ny, nz )
             bh   = (SHFTR(bc, nx) + bc )*0.5d0 !b_(i+1/2,j,k)
             ! rotaion B at (i+1/2, j, k)
             rotBh(:,:,:,nx) &                       ! dBz/dy - dBy/dz
                  = (SHFTR(bh(:,:,:,nz),ny) - SHFTL(bh(:,:,:,nz),ny))*h2i &
                  - (SHFTR(bh(:,:,:,ny),nz) - SHFTL(bh(:,:,:,ny),nz))*h2i
             rotBh(:,:,:,ny) &                       ! dBx/dz - dBz/dx
                  = (SHFTR(bh(:,:,:,nx),nz) - SHFTL(bh(:,:,:,nx),nz))*h2i &
                  - (SHFTR(bc(:,:,:,nz),nx) - bc(:,:,:,nz))*hi
             rotBh(:,:,:,nz) &                       ! dBy/dx - dBx/dy
                  = (SHFTR(bc(:,:,:,ny),nx) - bc(:,:,:,ny))*hi &
                  - (SHFTR(bh(:,:,:,nx),ny) - SHFTL(bh(:,:,:,nx),ny))*h2i
#define Bx bh(:,:,:,nx)
#define By bh(:,:,:,ny)
#define Bz bh(:,:,:,nz)
#define rotBx rotBh(:,:,:,nx)
#define rotBy rotBh(:,:,:,ny)
#define rotBz rotBh(:,:,:,nz)
             ! all flux: OD + HE + AD
             ! Ohmic dissipation
             dodh = (SHFTR(dod, nx) + dod)*0.5d0 !dod_{i+1/2}
             f(:,:,:,nx,nx) = 0.d0
             f(:,:,:,nx,ny) = -rotBz*dodh
             f(:,:,:,nx,nz) =  rotBy*dodh
             ! Hall effect
             dheh = (SHFTR(dhe, nx) + dhe)*0.5d0 !dhe_{i+1/2}
             f(:,:,:,nx,nx) = f(:,:,:,nx,nx) + 0.d0
             f(:,:,:,nx,ny) = f(:,:,:,nx,ny) + (Bx*rotBy - By*rotBx)*dheh
             f(:,:,:,nx,nz) = f(:,:,:,nx,nz) + (Bx*rotBz - Bz*rotBx)*dheh
             ! Ambipolar diffusion
             dadh = (SHFTR(dad, nx) + dad)*0.5d0 !dad_{i+1/2}
             f(:,:,:,nx,nx) = f(:,:,:,nx,nx) + 0.d0
             f(:,:,:,nx,ny) = f(:,:,:,nx,ny) + (Bx*(Bz*rotBx - Bx*rotBz) - By*(By*rotBz - Bz*rotBy))*dadh
             f(:,:,:,nx,nz) = f(:,:,:,nx,nz) + (Bx*(Bx*rotBy - By*rotBx) - Bz*(By*rotBz - Bz*rotBy))*dadh
          end do
          deallocate( bh, rotBh, dodh, dheh, dadh )
       end do
    end do
    ! -------------------------------------
    ! フラックス保存 (子のフラックスを信用)
    ! -------------------------------------
    call fmg_fluxcorrection(fmglev)
#undef SHFTL
#undef SHFTR
#undef bc
#undef Bx
#undef By
#undef Bz
#undef rotBx
#undef rotBy
#undef rotBz
  end subroutine fmg_ad_flux
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  subroutine cyclecomp(ncrd, nx, ny, nz, invert)
    integer,intent(IN) :: ncrd
    integer,intent(OUT) :: nx, ny, nz
    integer,intent(IN),optional :: invert
    integer,dimension(MX:MZ) :: mcycle
    integer,dimension(MX:MZ),parameter :: mcycleV = (/ MX, MY, MZ /)
    integer :: n
    do n = MX, MZ
       mcycle(n) = n
    enddo
    if ( present( invert ) ) then
       mcycle(MX:MX+size(mcycleV)-1) = cshift( mcycleV, -ncrd)
    else
       mcycle(MX:MX+size(mcycleV)-1) = cshift( mcycleV, ncrd )
    endif
    nx = mcycle(MX)
    ny = mcycle(MY)
    nz = mcycle(MZ)
  end subroutine cyclecomp
  ! ----------------------------------------------------------------
  ! smoothing operator
  ! ju = a code of unknown variable
  ! jrhs = a code of right hand side
  ! ----------------------------------------------------------------
  subroutine fmg_ad_relax(fmglev, ju, jrhs, resh2maxg)
    use mpilib
    use fmg_data
    integer,intent(IN) :: fmglev, ju, jrhs
    real(kind=DBL_KIND),intent(OUT) :: resh2maxg
    real(kind=DBL_KIND) :: hi, resh2, resh2max
    real(kind=DBL_KIND) :: ma, mb, mc, mdet, cod, che, cad
    real(kind=DBL_KIND),dimension(MX:MZ) :: db, res
    real(kind=DBL_KIND),dimension(MX:MZ,MX:MZ) :: mi
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, rhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dod, dhe, dad
    integer :: amrlev, gid, i, j, k, m, ipass, is,js,ks,ie,je,ke
#ifdef DEBUG_FMG_OUTPUT_RES
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: dbg ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
    myrank = get_myrank()
    resh2max = 0.d0
    Resmaxl = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    ! Red-black Gauss-Seidel iteration
    do ipass=1,2              ! 1-red, 2-black
       call fmg_ad_flux(fmglev, ju)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          hi = 1.d0/fmg_get_h( amrlev, fmglev )
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             ! 最細レベルか子供がいなければ,以下に進む
             call fmg_fp(amrlev, fmglev, gid, f)
             call fmg_arrp(amrlev, fmglev, gid, ju,   u)
             call fmg_arrp(amrlev, fmglev, gid, jrhs, rhs)
             call fmg_arrp(amrlev, fmglev, gid, IDOD, dod)
             call fmg_arrp(amrlev, fmglev, gid, IDHE, dhe)
             call fmg_arrp(amrlev, fmglev, gid, IDAD, dad)
#ifdef DEBUG_FMG_OUTPUT_RES
             if (fmglev == FMG_LevelMin) call fmg_arrp(amrlev, fmglev, gid, IDBG, dbg) ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
             do k = ks, ke
                do j = js, je
                   do i= is + mod(j+k+ipass,2), ie, 2
                      do m = Mmin, Mmax
                         res(m) = rhs(i,j,k,m) - u(i,j,k,m) &
                              -( f(i,j,k,MX,m)-f(i-1,j,k,MX,m) &
                              +  f(i,j,k,MY,m)-f(i,j-1,k,MY,m) &
                              +  f(i,j,k,MZ,m)-f(i,j,k-1,MZ,m) ) * hi * LAMBDA
                      end do
                      cod = dod(i,j,k) * 4*LAMBDA*hi**2
                      che = dhe(i,j,k)
                      cad = dad(i,j,k)
                      ma = 2*u(i,j,k,MX)*LAMBDA*hi**2
                      mb = 2*u(i,j,k,MY)*LAMBDA*hi**2
                      mc = 2*u(i,j,k,MZ)*LAMBDA*hi**2
                      mi(MX,MX) = cod**2+((3*mc**2+3*mb**2+2*ma**2)*cad+2)*cod+ma**2*che**2+(2*mc**4+(4*mb**2+3*ma**2)*mc**2+2*mb**4+3*ma**2*mb**2+ma**4)*cad**2+(3*mc**2+3*mb**2+2*ma**2)*cad+1
                      mi(MX,MY) = (-mc*che-ma*mb*cad)*cod+ma*mb*che**2+(((-2*mb**2-2*ma**2)*mc-2*mc**3)*cad-mc)*che+(-ma*mb*mc**2-ma*mb**3-ma**3*mb)*cad**2-ma*mb*cad
                      mi(MX,MZ) = (mb*che-ma*mc*cad)*cod+ma*mc*che**2+((2*mb*mc**2+2*mb**3+2*ma**2*mb)*cad+mb)*che+((-ma*mb**2-ma**3)*mc-ma*mc**3)*cad**2-ma*mc*cad
                      mi(MY,MX) = (mc*che-ma*mb*cad)*cod+ma*mb*che**2+((2*mc**3+(2*mb**2+2*ma**2)*mc)*cad+mc)*che+(-ma*mb*mc**2-ma*mb**3-ma**3*mb)*cad**2-ma*mb*cad
                      mi(MY,MY) = cod**2+((3*mc**2+2*mb**2+3*ma**2)*cad+2)*cod+mb**2*che**2+(2*mc**4+(3*mb**2+4*ma**2)*mc**2+mb**4+3*ma**2*mb**2+2*ma**4)*cad**2+(3*mc**2+2*mb**2+3*ma**2)*cad+1
                      mi(MY,MZ) = (-ma*che-mb*mc*cad)*cod+mb*mc*che**2+((-2*ma*mc**2-2*ma*mb**2-2*ma**3)*cad-ma)*che+((-mb**3-ma**2*mb)*mc-mb*mc**3)*cad**2-mb*mc*cad
                      mi(MZ,MX) = (-mb*che-ma*mc*cad)*cod+ma*mc*che**2+((-2*mb*mc**2-2*mb**3-2*ma**2*mb)*cad-mb)*che+((-ma*mb**2-ma**3)*mc-ma*mc**3)*cad**2-ma*mc*cad
                      mi(MZ,MY) = (ma*che-mb*mc*cad)*cod+mb*mc*che**2+((2*ma*mc**2+2*ma*mb**2+2*ma**3)*cad+ma)*che+((-mb**3-ma**2*mb)*mc-mb*mc**3)*cad**2-mb*mc*cad
                      mi(MZ,MZ) = cod**2+((2*mc**2+3*mb**2+3*ma**2)*cad+2)*cod+mc**2*che**2+(mc**4+(3*mb**2+3*ma**2)*mc**2+2*mb**4+4*ma**2*mb**2+2*ma**4)*cad**2+(2*mc**2+3*mb**2+3*ma**2)*cad+1
                      mdet = (cod**3+((4*mc**2+4*mb**2+4*ma**2)*cad+3)*cod**2+((mc**2+mb**2+ma**2)*che**2+(5*mc**4+(10*mb**2+10*ma**2)*mc**2+5*mb**4+10*ma**2*mb**2+5*ma**4)*cad**2+(8*mc**2+8*mb**2+8*ma**2)*cad+3)*cod+((2*mc**4+(4*mb**2+4*ma**2)*mc**2+2*mb**4+4*ma**2*mb**2+2*ma**4)*cad+mc**2+mb**2+ma**2)*che**2+(2*mc**6+(6*mb**2+6*ma**2)*mc**4+(6*mb**4+12*ma**2*mb**2+6*ma**4)*mc**2+2*mb**6+6*ma**2*mb**4+6*ma**4*mb**2+2*ma**6)*cad**3+(5*mc**4+(10*mb**2+10*ma**2)*mc**2+5*mb**4+10*ma**2*mb**2+5*ma**4)*cad**2+(4*mc**2+4*mb**2+4*ma**2)*cad+1)
                      mi = mi / mdet
                      db = matmul(mi, res)
                      u(i,j,k,:) = u(i,j,k,:) + db
                      Resmaxl(amrlev) = max(Resmaxl(amrlev), maxval(abs(res)))
#ifdef DEBUG_FMG_OUTPUT_RES
                      if (fmglev == FMG_LevelMin) dbg(i,j,k,:) = res ! debug 残差を出力
#endif !DEBUG_FMG_OUTPUT_RES
!!$                      u(i,j,k,:) = db
!!$                      resh2max = max(resh2max,maxval(abs(res)))
!!$                      if (amrlev == AMR_LevelMax) print *, maxval(abs(db)), maxval(abs(mi)), maxval(abs(res))
                   enddo
                end do
             enddo
!!$             if (amrlev == AMR_LevelMax) print *, maxval(u(is:ie,js:je,ks:ke,0))
          enddo
       enddo
    enddo
    call mpi_allreduce(Resmaxl, Resmaxg, size(Resmaxg), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
!!$    call mpi_allreduce(resh2max, resh2maxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
!!$    resh2maxg = 1.d5
!!$    if (myrank == PRIMARY_RANK) print *, 'Resmaxg', Resmaxg
!!$    resh2maxg = 1.d-5
  end subroutine fmg_ad_relax
  !-------------------------------------------------------------------------
  ! residual
  !-------------------------------------------------------------------------
  subroutine fmg_ad_resid(fmglev, jres, ju, jrhs)
    use fmg_data
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: res, u, rhs
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: amrlev, gid, i, j, k, m, is, ie, js, je, ks, ke
    real(kind=DBL_KIND) :: hi
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    call fmg_ad_flux(fmglev, ju)
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
  end subroutine fmg_ad_resid
  !-------------------------------------------------------------------------
  ! tau correction for nonlinear multigrid
  !-------------------------------------------------------------------------
  subroutine fmg_ad_tau(fmglev, jtau, ju)
    integer,intent(IN) :: fmglev, jtau, ju
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau, u
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
    integer :: gid, i, j, k, is, ie, js, je, ks, ke, fmglevf, amrlev, m
    real(kind=DBL_KIND) :: hi
#ifdef DEBUG_FMG_OUTPUT_TAU
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: dbg
#endif !DEBUG_FMG_OUTPUT_TAU

    fmglevf = fmglev-1 ! fine level

    if ( fmglev > FMG_LevelMax ) print *, '*** error: invarid fmglev in fmg_tau'
    if ( fmglevf < FMG_LevelMin ) print *, '*** error: invarid fmglevf in fmg_tau'
    myrank = get_myrank()

    ! tau = - R L uf
    call fmg_get_gridsize(fmglevf, is,js,ks,ie,je,ke)
    call fmg_ad_flux(fmglevf, ju) ! in fine grids.
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
    call fmg_ad_flux(fmglev, ju) ! uf was already restricted (R uf)
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

#ifdef DEBUG_FMG_OUTPUT_TAU
    if (fmglev == FMG_LevelMin+1) then
       call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             tau => fmg_get_arrp(amrlev, fmglev, gid, jtau)
             dbg => fmg_get_arrp(amrlev, fmglev, gid, IDBG)
             dbg = tau
          end do
       end do
!!$       print *, '*** tau is stored in IDBG'
    end if
#endif !DEBUG_FMG_OUTPUT_TAU
  end subroutine fmg_ad_tau

