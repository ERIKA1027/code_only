
!
! This module contains subroutines related to utilities of a uniform grid.
!
#include "config.h"
#define ALLOCATESAFE(A) \
  allocate(A, stat=ierr); \
  if (ierr /= 0) then; \
     print *, '** Allocation fail',  #A, " errorcode =", ierr; \
     stop; \
  endif

module uniformgrid
  use grid, only : Undefi, Imin, Jmin, Kmin, Imax, Jmax, Kmax, Mmin, Mmax, get_Up, get_Xp, get_Yp, get_Zp
  use overBlockCoordinates
  use mpilib
  use string, only : CHARLEN
  implicit none
  private
  character(len=2),parameter :: PrefixDefault = 'ug' ! Prefix of filename for uniform grid
  character(len=CHARLEN) :: ThisPrefix
  integer,save :: LevelFinest = Undefi        ! Finest level of uniform grid
  integer,save :: LevelResolution = Undefi    ! level specifies resolution
  public :: uniformgrid_write
contains
  !---------------------------------------------------------------------
  ! ug_gather
  ! This routine gathers data to the primary rank and packs it to a
  ! uniform grid.
  ! INPUT:
  !  (xmin, ymin, zmin, xmax, ymax, zmax)
  !        ......... vertex of a uniform grid in physical coordinates.
  !  level ........ level specifing resolution of the uniform grid
  !  interpolate .. specify whether interpolate or not (optional). Default is interplate=.false.
  !  prefix ....... prefix of file name (optional). The defalult value is given by PrefixDefalut.
  !---------------------------------------------------------------------
  subroutine uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate, prefix)
    use grid, only : Step, Time, LevelMax, Lmin
    real(KIND=DBL_KIND),intent(IN) :: xmin, ymin, zmin, xmax, ymax, zmax
    integer,intent(IN) :: level
    logical,intent(IN),optional :: interpolate
    character(len=*),intent(IN),optional :: prefix
    type(t_obRectPhys) :: rectPhys, rectPhysIn, rectCompBox, rectPhysAnd
    type(t_obRect) :: rect
    type(t_obPoint) :: pL, pR, pfL, pfR
    type(t_obPointPhys) :: ppL
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: uc, uf, u
    real(kind=DBL_KIND),dimension(:),allocatable :: x, y, z, xc, yc, zc
    logical,dimension(:,:,:),allocatable :: boolmap
    integer :: i, j, k, m, lev, ie, je, ke

#ifdef FORBIT_WRITEDATA
    return
#endif

    if ( present(prefix) ) then
       ThisPrefix = prefix
    else
       ThisPrefix = PrefixDefault
    endif


    if ( level > LevelMax ) then
       print *, '*** uniform_write: level is larger than LevelMax: level =',level, ' LevelMax =', LevelMax
       print *, 'LevelFinest is trancated to ',LevelMax
    endif
    LevelFinest = min(level, LevelMax)
    LevelFinest = max(LevelFinest, Lmin)
    LevelResolution = level

    call ob_assignCoordPhysToRectPhys( (/ xmin, ymin, zmin, xmax, ymax, zmax /), rectPhysIn )
    call ob_computationBoxOfRectPhys( rectCompBox )
    call ob_rectPhysAnd( rectPhysIn, rectCompBox, rectPhysAnd )
    call ob_rectPhysTrim( rectPhysAnd, CellWidth(:,LevelFinest)/2.d0, rectPhys )

    ! -------------------
    ! physical variables
    ! -------------------
    myrank = get_myrank()
    call ob_RectPhys2RectOb(rectPhys, Lmin, rect )
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    if ( myrank == PRIMARY_RANK ) &
         ALLOCATESAFE( u( 0:pR%p(MX)-pL%p(MX), 0:pR%p(MY)-pL%p(MY), 0:pR%p(MZ)-pL%p(MZ), Mmin:Mmax) )
    do lev = Lmin, LevelFinest
       if ( myrank == PRIMARY_RANK ) then
          ALLOCATESAFE( uc(ARRAYSIZE4(u)) )
          ALLOCATESAFE( boolmap(ARRAYSIZE3(u)) )
       endif
       call ug_gather( uc, boolmap, rect )
       if ( myrank == PRIMARY_RANK ) then
          if (lev == Lmin) boolmap(:,:,:) = .true.
          do m = Mmin, Mmax
             where (boolmap) u(:,:,:,m) = uc(:,:,:,m)
          enddo
          deallocate( uc, boolmap )
       endif
       if ( lev == LevelFinest ) exit
       if ( myrank == PRIMARY_RANK ) then
          ALLOCATESAFE( uf( 0:(pR%p(MX)-pL%p(MX)+1)*2-1, 0:(pR%p(MY)-pL%p(MY)+1)*2-1, 0:(pR%p(MZ)-pL%p(MZ)+1)*2-1, Mmin:Mmax) )
          call ug_interp(u, uf, interpolate)
          deallocate( u )

          ! for next level
          call ob_RectPhys2RectOb(rectPhys, lev+1, rect)
          call ob_extractPointFromRect(pfL, rect, 'L')
          call ob_extractPointFromRect(pfR, rect, 'R')
          ALLOCATESAFE( u( 0:pfR%p(MX)-pfL%p(MX), 0:pfR%p(MY)-pfL%p(MY), 0:pfR%p(MZ)-pfL%p(MZ), Mmin:Mmax) )
       ! corresponding fine index = pL%p(MX)*2 : pR%p(MX)*2+1
          u(:,:,:,:) = uf( &
               pfL%p(MX)-pL%p(MX)*2 : pfL%p(MX)-pL%p(MX)*2+pfR%p(MX)-pfL%p(MX), &
               pfL%p(MY)-pL%p(MY)*2 : pfL%p(MY)-pL%p(MY)*2+pfR%p(MY)-pfL%p(MY), &
               pfL%p(MZ)-pL%p(MZ)*2 : pfL%p(MZ)-pL%p(MZ)*2+pfR%p(MZ)-pfL%p(MZ), &
               Mmin:Mmax)
          pL = pfL
          pR = pfR
          deallocate( uf )
       endif
    enddo
    if ( myrank /= PRIMARY_RANK ) return
    ! -----------
    ! coordinates
    ! -----------
    call ob_RectPhys2RectOb(rectPhys, LevelFinest, rect)
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    ALLOCATESAFE( x(0:pR%p(MX)-pL%p(MX)) )
    ALLOCATESAFE( y(0:pR%p(MY)-pL%p(MY)) )
    ALLOCATESAFE( z(0:pR%p(MZ)-pL%p(MZ)) )
    call ob_RectOb2RectPhys( rect, rectPhys )
    call ob_extractPointPhysFromRectPhys(ppL, rectPhys, 'L')
    do i = lbound(x,1), ubound(x,1)
       x(i) = i * CellWidth(MX, LevelFinest) + ppL%p(MX)
    enddo
    do j = lbound(y,1), ubound(y,1)
       y(j) = j * CellWidth(MY, LevelFinest) + ppL%p(MY)
    enddo
    do k = lbound(z,1), ubound(z,1)
       z(k) = k * CellWidth(MZ, LevelFinest) + ppL%p(MZ)
    enddo
    ! --------------------
    ! downsize by sampling
    ! --------------------
    if ( level < Lmin ) then
       call sampling( Lmin - level, u, x, y, z, ie, je, ke )
       ! swap u and trancate u
       ALLOCATESAFE( uc(0:ie,0:je,0:ke,Mmin:Mmax) )
       ALLOCATESAFE( xc(0:ie) )
       ALLOCATESAFE( yc(0:je) )
       ALLOCATESAFE( zc(0:ke) )
       uc = u(0:ie,0:je,0:ke,Mmin:Mmax)
       xc = x(0:ie)
       yc = y(0:je)
       zc = z(0:ke)
       deallocate( u, x, y, z )
       ALLOCATESAFE( u(0:ie,0:je,0:ke,Mmin:Mmax) )
       ALLOCATESAFE( x(0:ie) )
       ALLOCATESAFE( y(0:je) )
       ALLOCATESAFE( z(0:ke) )
       u = uc
       x = xc
       y = yc
       z = zc
       deallocate( uc, xc, yc, zc )
    endif
    ! ----------
    ! write data
    ! ----------
    call ug_write(u, x, y, z, Time(LevelFinest), Step(LevelFinest) )
    call ug_fld(u, x, y, z, Time(LevelFinest), Step(LevelFinest) )
    deallocate( u, x, y, z )
  end subroutine uniformgrid_write
  !---------------------------------------------------------------------
  ! Sampling of u, x, y, z
  !---------------------------------------------------------------------
  subroutine sampling( dlevel, u, x, y, z, ie, je, ke )
    use eos, only : w2u, u2w
    integer,intent(IN) :: dlevel
    real(kind=DBL_KIND),intent(INOUT),dimension(0:,0:,0:,Mmin:) ::  u
    real(kind=DBL_KIND),intent(INOUT),dimension(0:) :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable ::  uc, wc, w
    real(kind=DBL_KIND),dimension(:),allocatable :: xc, yc, zc
    integer,intent(OUT) :: ie, je, ke
    real(kind=DBL_KIND),parameter :: dvh = 1, dvc = 8
    integer :: dlev, is, js, ks, if, jf, kf, m, ic, jc, kc, szi, szj, szk

    is = lbound(u,1);   js = lbound(u,2);   ks = lbound(u,3)
    ie = ubound(u,1);   je = ubound(u,2);   ke = ubound(u,3)
    szi = size(u,1);    szj = size(u,2);    szk = size(u,3)
    do dlev = 1, dlevel
       ALLOCATESAFE( w(is:ie,js:je,ks:ke,Mmin:Mmax) )
       ie = is + szi/2 -1
       je = js + szj/2 -1
       ke = ks + szk/2 -1
       ALLOCATESAFE( wc(is:ie,js:je,ks:ke,Mmin:Mmax) )
       ALLOCATESAFE( uc(is:ie,js:je,ks:ke,Mmin:Mmax) )
       ALLOCATESAFE( xc(is:ie) )
       ALLOCATESAFE( yc(js:je) )
       ALLOCATESAFE( zc(ks:ke) )
       call u2w(u, w, dvh)
       do m=Mmin,Mmax
          do kc=ks,ke
             do jc=js,je
                do ic=is,ie
                   if = IJKF(ic,is)
                   jf = IJKF(jc,js)
                   kf = IJKF(kc,ks)
                   wc(ic,jc,kc,m) = &
                         w(if,  jf,  kf,  m)+w(if+1,jf,  kf,  m) &
                        +w(if,  jf+1,kf,  m)+w(if,  jf,  kf+1,m) &
                        +w(if+1,jf+1,kf,  m)+w(if+1,jf,  kf+1,m) &
                        +w(if,  jf+1,kf+1,m)+w(if+1,jf+1,kf+1,m)
                enddo
             enddo
          enddo
       enddo
       call w2u( wc, uc, dvc)
       u(is:ie,js:je,ks:ke,Mmin:Mmax) = uc(is:ie,js:je,ks:ke,Mmin:Mmax)
       ! x, y, z
       do ic = is, ie
          if = IJKF(ic, is)
          xc(ic) = (x(if) + x(if+1))/2.d0
       enddo
       x(is:ie) = xc(is:ie)
       do jc = js, je
          jf = IJKF(jc, js)
          yc(jc) = (y(jf) + y(jf+1))/2.d0
       enddo
       y(js:je) = yc(js:je)
       do kc = ks, ke
          kf = IJKF(kc, ks)
          zc(kc) = (z(kf) + z(kf+1))/2.d0
       enddo
       z(ks:ke) = zc(ks:ke)
       deallocate(w, wc, uc, xc, yc, zc)
    enddo
  end subroutine sampling
  !---------------------------------------------------------------------
  ! gather data from every node by mpi_allreduce
  ! INPUT:
  !   rect = Rectangle region in over-block coordinates
  !          Indexes specifining rectablel region of interest.
  !          This region is gathered into uGlobal
  !
  ! OUTPUT:
  !   uGlobal = gathered data
  !   boolmap = .ture. where cell has data.
  !---------------------------------------------------------------------
  subroutine ug_gather(uGlobal, boolmap, rect)
    use mpilib
    type(t_obRect),intent(IN) :: rect
    real(kind=DBL_KIND),dimension(0:,0:,0:,Mmin:),intent(OUT) :: uGlobal
    logical,dimension(0:,0:,0:),intent(OUT) :: boolmap
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: uLocal
    logical,dimension(:,:,:),allocatable :: bmap
    type(t_obPoint) :: point, pL, pR
    integer :: gid, rank, i, j, k, is, js, ks, ie, je, ke, ius, jus, kus, iue, jue, kue
    integer :: level
    integer(kind=LLONG_KIND) :: iob, job, kob
    integer :: isz = Imax-Imin+1, jsz = Jmax-Jmin+1, ksz = Kmax-Kmin+1
    logical :: ilast = .false. , jlast = .false. , klast = .false.
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: buf(6)

    if ( get_myrank() == PRIMARY_RANK ) &
         buf = (/lbound(uGlobal,1),lbound(uGlobal,2),lbound(uGlobal,3), ubound(uGlobal,1),ubound(uGlobal,2),ubound(uGlobal,3) /)
    call mpi_bcast(buf, size(buf), MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr )
    is = buf(1) ; js = buf(2) ; ks = buf(3)
    ie = buf(4) ; je = buf(5) ; ke = buf(6)

    ALLOCATESAFE( uLocal(is:ie,js:je,ks:ke,Mmin:Mmax) )
    ALLOCATESAFE( bmap(is:ie,js:je,ks:ke) )

    ! initialize of returned arguments
    if ( get_myrank() == PRIMARY_RANK ) then
       uGlobal(:,:,:,:) = 0.d0
       boolmap(:,:,:) = .false.
    endif
    uLocal(:,:,:,:) = 0.d0
    bmap(:,:,:) = .false.

    myrank = get_myrank()
    call ob_extractLevelFromRect(level, rect)           ! level
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    do kob = pL%p(MZ), pR%p(MZ)+ksz, ksz
       klast = ( kob/ksz == pR%p(MZ)/ksz )
       if ( kob/ksz > pR%p(MZ)/ksz ) exit
       do job = pL%p(MY), pR%p(MY)+jsz, jsz
          jlast =  ( job/jsz == pR%p(MY)/jsz )
          if ( job/jsz > pR%p(MY)/jsz ) exit
          do iob = pL%p(MX), pR%p(MX)+isz, isz
             ilast = ( iob/isz == pR%p(MX)/isz )
             if ( iob/isz > pR%p(MX)/isz ) exit
             call ob_assignCoordToPoint(point, (/iob, job, kob/), level)
             call ob_getBlockFromPoint(point, gid, rank, i,j,k)
             if ( gid == Undefi ) cycle
             if ( rank /= myrank ) cycle
             if ( rank == MPI_PROC_NULL ) cycle
             ! region of block
             is = Imin; js = Jmin; ks = Kmin
             ie = Imax; je = Jmax; ke = Kmax
             if ( iob == pL%p(MX) ) is = i
             if ( job == pL%p(MY) ) js = j
             if ( kob == pL%p(MZ) ) ks = k
             if ( ilast ) ie = mod(pR%p(MX), isz)
             if ( jlast ) je = mod(pR%p(MY), jsz)
             if ( klast ) ke = mod(pR%p(MZ), ksz)
             ! region of ugLocal
             ius = -mod(pL%p(MX), isz) + iob - pL%p(MX) + is
             iue = ius + ie - is
             jus = -mod(pL%p(MY), jsz) + job - pL%p(MY) + js
             jue = jus + je - js
             kus = -mod(pL%p(MZ), ksz) + kob - pL%p(MZ) + ks
             kue = kus + ke - ks
             ! now assign
             u => get_Up(gid)
             uLocal(ius:iue,jus:jue,kus:kue,:) = u(is:ie,js:je,ks:ke,:)
             bmap(ius:iue,jus:jue,kus:kue) = .true.
             if ( ilast ) exit
          enddo
          if ( jlast ) exit
       enddo
       if ( klast ) exit
    enddo
    call mpi_reduce( uLocal, uGlobal, size(uLocal), MPI_DOUBLE_PRECISION, MPI_SUM, PRIMARY_RANK, MPI_COMM_WORLD, ierr )
    call mpi_reduce( bmap, boolmap, size(bmap), MPI_LOGICAL, MPI_LOR, PRIMARY_RANK, MPI_COMM_WORLD, ierr )
    deallocate( uLocal, bmap )
  end subroutine ug_gather
  !---------------------------------------------------------------------
  ! interpolation of u, x, y, z
  !---------------------------------------------------------------------
  subroutine ug_interp(uc, uf, interpolate)
    use eos, only : w2u, u2w
    real(kind=DBL_KIND),dimension(0:,0:,0:,Mmin:),intent(IN) :: uc
    real(kind=DBL_KIND),dimension(0:,0:,0:,Mmin:),intent(OUT) :: uf
    logical,intent(IN),optional :: interpolate
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: grad
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: wc
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: wf
    real(kind=DBL_KIND) :: di, dj, dk, hf, hc, hfi, hci, dvc, dvf
    integer :: i, j, k, m, if, jf, kf, ic, jc, kc
    integer :: ics, ice, jcs, jce, kcs, kce, ifs, ife, jfs, jfe, kfs, kfe
    ! minmod limiter
#define MINMOD(a,b) sign(1.d0,(a))*max(0.d0,min(abs(a),sign(1.d0,(a))*(b)))

    ALLOCATESAFE(grad(ARRAYSIZE3(uc),MX:MZ))
    ALLOCATESAFE(wc(ARRAYSIZE4(uc)))
    ALLOCATESAFE(wf(ARRAYSIZE4(uf)))
    hf = 1.d0                   ! fine cell width
    hc = 2* hf                  ! coarse cell width
    dvc = 1.d0                  ! dummy volume
    dvf = 1.d0                  ! dummy volume
    hfi = 1.d0/hf
    hci = 1.d0/hc
    ! shape of indexs
    ics = lbound(uc, 1)
    ice = ubound(uc, 1)
    jcs = lbound(uc, 2)
    jce = ubound(uc, 2)
    kcs = lbound(uc, 3)
    kce = ubound(uc, 3)
    ifs = lbound(uf, 1)
    ife = ubound(uf, 1)
    jfs = lbound(uf, 2)
    jfe = ubound(uf, 2)
    kfs = lbound(uf, 3)
    kfe = ubound(uf, 3)
    call u2w(uc, wc, dvc)
    do m = Mmin, Mmax
       grad(:,:,:,:) = 0.d0
       if (present( interpolate )) then
          do ic = ics+1, ice-1
             grad(ic,:,:,MX) = ( MINMOD(wc(ic+1,:,:,m)-wc(ic,:,:,m), wc(ic,:,:,m)-wc(ic-1,:,:,m)) )*hci
          enddo
          do jc = jcs+1, jce-1
             grad(:,jc,:,MY) = ( MINMOD(wc(:,jc+1,:,m)-wc(:,jc,:,m), wc(:,jc,:,m)-wc(:,jc-1,:,m)) )*hci
          enddo
          do kc = kcs+1, kce-1
             grad(:,:,kc,MZ) = ( MINMOD(wc(:,:,kc+1,m)-wc(:,:,kc,m), wc(:,:,kc,m)-wc(:,:,kc-1,m)) )*hci
          enddo
          if ( ics /= ice ) then
             grad(ics,:,:,MX) = ( wc(ics+1,:,:,m)-wc(ics,:,:,m) )*hci
             grad(ice,:,:,MX) = ( wc(ice,:,:,m)-wc(ice-1,:,:,m) )*hci
          endif
          if ( jcs /= jce ) then
             grad(:,jcs,:,MY) = ( wc(:,jcs+1,:,m)-wc(:,jcs,:,m) )*hci
             grad(:,jce,:,MY) = ( wc(:,jce,:,m)-wc(:,jce-1,:,m) )*hci
          endif
          if ( kcs /= kce ) then
             grad(:,:,kcs,MZ) = ( wc(:,:,kcs+1,m)-wc(:,:,kcs,m) )*hci
             grad(:,:,kce,MZ) = ( wc(:,:,kce,m)-wc(:,:,kce-1,m) )*hci
          endif
       endif
       do kf = kfs, kfe
          do jf = jfs, jfe
             do if = ifs, ife
                ic = IJKC(if,0)
                jc = IJKC(jf,0)
                kc = IJKC(kf,0)
                di = ( modulo(if,2) - 0.5d0 )*hf
                dj = ( modulo(jf,2) - 0.5d0 )*hf
                dk = ( modulo(kf,2) - 0.5d0 )*hf
                wf(if,jf,kf,m) = wc(ic,jc,kc,m)+grad(ic,jc,kc,MX)*di+grad(ic,jc,kc,MY)*dj+grad(ic,jc,kc,MZ)*dk
             enddo
          enddo
       enddo
    enddo
    call w2u(wf, uf, dvf)
#undef MINMOD
    deallocate( grad, wc, wf)
  end subroutine ug_interp
  !---------------------------------------------------------------------
  ! Write data in uniform grid
  !---------------------------------------------------------------------
#define FH 11
  subroutine ug_write(u,x,y,z,time,step)
    use mpilib
    use string
    use io_util, only : wchar
    real(kind=DBL_KIND),dimension(:,:,:,:),intent(IN) :: u
    real(kind=DBL_KIND),dimension(:),intent(IN) :: x, y, z
    real(kind=DBL_KIND),intent(IN) :: time
    integer(kind=LLONG_KIND),intent(IN) :: step
    character(len=CHARLEN) :: fn, dir, suffix
    real,dimension(:,:,:,:),allocatable :: u4
    real,dimension(:),allocatable :: x4, y4, z4

    if ( get_myrank() /= PRIMARY_RANK ) return
    ALLOCATESAFE( u4(ARRAYSIZE4(u)) )
    ALLOCATESAFE( x4(size(x)) )
    ALLOCATESAFE( y4(size(y)) )
    ALLOCATESAFE( z4(size(z)) )
    u4 = u
    x4 = x
    y4 = y
    z4 = z
    call read_env('DIR', dir)
    call read_env('SUFFIX', suffix)
    fn = get_filename(dir, ThisPrefix, suffix)
    call wchar(6,'uniform grid data = '//fn)
    open(FH, file=fn, form='unformatted')
    write(FH) shape(u), step, time  ! size of u
    write(FH) x4
    write(FH) y4
    write(FH) z4
    write(FH) u4
    call flush(FH)
    close(FH)
    deallocate( u4, x4, y4, z4 )
  end subroutine ug_write
  !---------------------------------------------------------------------
  ! Write fld file
  !---------------------------------------------------------------------
  subroutine ug_fld(u,x,y,z,time,step)
    use mpilib
    use io_util, only : wchar
    use string
    real(kind=DBL_KIND),dimension(:,:,:,:),intent(IN) :: u
    real(kind=DBL_KIND),dimension(:),intent(IN) :: x, y, z
    real(kind=DBL_KIND),intent(IN) :: time
    integer(kind=LLONG_KIND),intent(IN) :: step
    character(len=CHARLEN) :: fn, dir, suffix
    character(len=3),parameter :: Suffixfld = 'fld'
    integer :: ni, nj, nk, nm, ndim, nspace, idepthv, idepthc, hbyte, fbyte, ivol, iskip, m
    if ( get_myrank() /= PRIMARY_RANK ) return
    ni = size(u,1)
    nj = size(u,2)
    nk = size(u,3)
    nm = size(u,4)
    ndim = 3
    nspace = 3
    call read_env('DIR', dir)
    ! filename
    fn = get_filename(dir, ThisPrefix, Suffixfld)
    call wchar(6,'uniform fld file = '//fn)
    open(FH, file=fn)
    write(FH, '(A)') "# AVS"
    call wchar(FH,concat('# Ni =',strcmprs(num2char(ni))))
    call wchar(FH,concat('# Nj =',strcmprs(num2char(nj))))
    call wchar(FH,concat('# Nk =',strcmprs(num2char(nk))))
    call wchar(FH,concat('# Nm =',strcmprs(num2char(nm))))
    call wchar(FH,concat('# lstep=',strcmprs(num2char(step))))
    call wchar(FH,concat('# time =',strcmprs(num2char(time))))

    call wchar(FH,concat('ndim=',strcmprs(num2char(ndim))))
    call wchar(FH,concat('nspace=',strcmprs(num2char(nspace))))
    call wchar(FH,concat('dim1=',strcmprs(num2char(ni))))
    call wchar(FH,concat('dim2=',strcmprs(num2char(nj))))
    call wchar(FH,concat('dim3=',strcmprs(num2char(nk))))
    call wchar(FH,concat('veclen=',strcmprs(num2char(nm))))

    write(FH, '(A)') 'data=float'
    write(FH, '(A)') 'field=rectilinear'
    write(FH, '(A)') '# coordinate'

    !

    idepthv = 4               !byte per point for variables
    idepthc = 4               !byte per point for coordinates
    hbyte = 4                 !header for fortran unformatted format
    fbyte = 4                 !footer for fortran unformatted format
    ivol=ni*nj*nk*idepthv
    iskip = 0
    iskip = iskip + hbyte  !skip header
    iskip = iskip + 4*4    !skip for Ni,Nj,Nk,Nm (int)
    iskip = iskip + 8      !skip for lstep(l) (8byte int)
    iskip = iskip + 8      !skip for time(l) (dbl)
    iskip = iskip + fbyte  !skip footer

    call read_env('SUFFIX', suffix)
    fn = get_filename('', ThisPrefix, suffix)
    ! x
    iskip = iskip + hbyte
    call wchar(FH,concat( &
         concat('coord 1 file=',fn), &
         concat(' filetype=binary skip=',strcmprs(num2char(iskip)))))
    iskip = iskip + ni*idepthc+fbyte
    ! y
    iskip = iskip + hbyte
    call wchar(FH,concat( &
         concat('coord 2 file=',fn), &
         concat(' filetype=binary skip=',strcmprs(num2char(iskip)))))
    iskip = iskip + nj*idepthc+fbyte
    ! z
    iskip = iskip + hbyte
    call wchar(FH,concat( &
         concat('coord 3 file=',fn), &
         concat(' filetype=binary skip=',strcmprs(num2char(iskip)))))
    iskip = iskip + nk*idepthc+fbyte

    iskip = iskip + hbyte
    write(FH,'(A)') '# variable'
    do m=1, Nm
       call wchar(FH,concat( &
            concat('variable '//strcmprs(num2char(m)), &
            concat(' file=',fn)), &
            concat(' filetype=binary skip=', &
            strcmprs(num2char(iskip)))))
       iskip = iskip + ivol
    enddo
    iskip =  iskip + fbyte
    ! check sum
    call wchar(FH,concat( &
         concat('# size of file, ',fn), &
         concat(' = ',strcmprs(num2char(iskip)))))
    call flush(FH)
    close(FH)
  end subroutine ug_fld
  !-------------------------------------------------------------------------
  ! make file name for IO of uniform grid
  ! file name = dir + prefix + stage + . + level + . + suffix
  !-------------------------------------------------------------------------
  function get_filename(dir, prefix, suffix) result(fn)
    use string, only : concat, CHARLEN, num2char
    use grid, only : Step, LevelMax
    character(len=*),intent(IN) :: dir, prefix, suffix
    character(len=CHARLEN) :: fn
    character(len=1) :: dot = '.'
    fn = dir                               ! /dir/
    fn = concat(fn,prefix)                 ! /dir/st
    fn = concat(fn,num2char(Step(LevelMax))) ! /dir/st12000
    fn = concat(fn,dot)                    ! /dir/st12000.
    fn = concat(fn,num2char(LevelResolution))    ! /dir/st12000.0
    fn = concat(fn,dot)                    ! /dir/st12000.0.
    fn = concat(fn,suffix)                 ! /dir/st12000.0.d
  end function get_filename
end module uniformgrid
