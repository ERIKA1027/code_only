#include "debug_fmg.h"
!-------------------------------------------------------------------------
! tau correction
! tau = L R uf - R L uf
!-------------------------------------------------------------------------
subroutine vmg_ad_tau(amrlev, jtau, ju)
  use fmg_data
  integer,intent(IN) :: amrlev, jtau, ju
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: tau, u
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
  integer :: gid, i, j, k, amrlevf, m
  real(kind=DBL_KIND) :: hi
#ifdef DEBUG_VMG_OUTPUT_TAU
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: dbg
#endif !DEBUG_VMG_OUTPUT_TAU

  amrlevf = amrlev+1 ! fine level
  if ( amrlevf > AMR_LevelMax ) return
  myrank = get_myrank()

  ! tau = - R L uf
  call vmg_ad_flux(amrlevf, ju) ! in child (fine) grids.
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
  call vmg_ad_flux(amrlev, ju) ! uf was already restricted (R uf)
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

#ifdef DEBUG_VMG_OUTPUT_TAU
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
     if (.not. fmg_have_child(gid, amrlev)) cycle
     tau => fmg_get_arrp(amrlev, FMG_Level, gid, jtau)
     dbg => fmg_get_arrp(amrlev, FMG_Level, gid, IDBG)
     dbg = tau
  end do
#endif !DEBUG_VMG_OUTPUT_TAU
!!$  print *, 'tau2', maxval(abs(tau(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:))), amrlev

end subroutine vmg_ad_tau
! ----------------------------------------------------------------
! solve flux for one composit grid (for one VMG level)
! ----------------------------------------------------------------
subroutine vmg_ad_flux(amrlev, ju)
  use fmg_boundary_phys
  integer,intent(IN) :: amrlev, ju
  real(kind=DBL_KIND) :: hi, h2i
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u
  real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dod, dhe, dad
  integer :: n, nx, ny, nz
  integer,parameter :: DIMOFFSET = 1-MX
  integer :: fmglev, gid
  real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: bh, rotBh
  real(kind=DBL_KIND),dimension(:,:,:),allocatable :: dodh, dheh, dadh
  ! --------
  ! 境界条件
  ! --------
  call vmg_ghostcell(amrlev, ju)
  call vmg_boundary_u(amrlev, FMG_Level, ju)
  ! ----------
  ! フラックス
  ! ----------
  fmglev = FMG_Level
  hi = 1.d0/fmg_get_h( amrlev, fmglev )
  h2i = hi * 0.5d0
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
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
#undef SHFTL
#undef SHFTR
#undef bc
#undef Bx
#undef By
#undef Bz
#undef rotBx
#undef rotBy
#undef rotBz
end subroutine vmg_ad_flux
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
  ! ----------------------------------------------------------------
subroutine vmg_ad_relax(amrlev, ju, jrhs)
  use mpilib
  integer,intent(IN) :: amrlev, ju, jrhs
  real(kind=DBL_KIND) :: hi
  real(kind=DBL_KIND) :: ma, mb, mc, mdet, cod, che, cad
  real(kind=DBL_KIND),dimension(MX:MZ) :: db, res
  real(kind=DBL_KIND),dimension(MX:MZ,MX:MZ) :: mi
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: f
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, rhs
  real(kind=DBL_KIND),pointer,dimension(:,:,:) :: dod, dhe, dad
  integer :: gid, i, j, k, m, ipass, is,js,ks,ie,je,ke, fmglev
  real(kind=DBL_KIND) :: resmax, resmax_g
#ifdef DEBUG_VMG_OUTPUT_RES
  real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: dbg ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
  myrank = get_myrank()
  fmglev = FMG_Level
  resmax = 0.d0
  hi = 1.d0/fmg_get_h( amrlev, fmglev )
  call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
  ! Red-black Gauss-Seidel iteration
  do ipass=1,2              ! 1-red, 2-black
     call vmg_ad_flux(amrlev, ju)
     do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
        call fmg_fp(amrlev, fmglev, gid, f)
        call fmg_arrp(amrlev, fmglev, gid, ju,   u)
        call fmg_arrp(amrlev, fmglev, gid, jrhs, rhs)
        call fmg_arrp(amrlev, fmglev, gid, IDOD, dod)
        call fmg_arrp(amrlev, fmglev, gid, IDHE, dhe)
        call fmg_arrp(amrlev, fmglev, gid, IDAD, dad)
#ifdef DEBUG_VMG_OUTPUT_RES
        call fmg_arrp(amrlev, fmglev, gid, IDBG, dbg) ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
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
                 resmax = max(resmax, maxval(abs(res)))
#ifdef DEBUG_VMG_OUTPUT_RES
                 dbg(i,j,k,:) = res ! debug 残差を出力
#endif !DEBUG_VMG_OUTPUT_RES
              enddo
           end do
        enddo
     enddo
  enddo
  call mpi_allreduce(resmax, resmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  Resmaxg(amrlev) = resmax_g
end subroutine vmg_ad_relax
