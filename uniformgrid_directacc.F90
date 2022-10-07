
!
! This module contains subroutines related to utilities of a uniform grid.
!
! When compiler supports binary IO, use uniformgrid.F90
!
#include "config.h"
#include "packarr.h"
#include "overBlockCoordinates.h"
#include "recl.h"
module uniformgrid
  use grid, only : Undefi, Imin, Jmin, Kmin, Imax, Jmax, Kmax, Mmin, Mmax, get_Up, get_Xp, get_Yp, get_Zp
  use overBlockCoordinates
  use mpilib
  use string, only : CHARLEN
  implicit none
  private
  type(t_obRectPhys),save :: rectCompBoxPhys
  type(t_obRect),save :: rectCompBox
  character(len=2),parameter :: PrefixDefault = 'ug' ! Prefix of filename for uniform grid
  character(len=CHARLEN) :: ThisPrefix
  integer,save :: LevelFinest = Undefi        ! Finest level of uniform grid
  integer,save :: LevelResolution = Undefi    ! level specifies resolution
  integer,save :: LevelFinalResolution = Undefi    ! level of final resolution
  integer,parameter :: U2W = 1
  integer,parameter :: W2U = 2
  character(len=CHARLEN),parameter :: UNIX_CMD_CAT = 'cat'
  type t_ugArray
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: arr => null()
     integer(kind=LLONG_KIND),dimension(MX:MZ) :: offsets
     integer :: nrank = MPI_PROC_NULL
     logical :: allocated = .false.
  end type t_ugArray
  public :: uniformgrid_write
contains
  !---------------------------------------------------------------------
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
    use grid, only : Step, Time
    real(KIND=DBL_KIND),intent(IN) :: xmin, ymin, zmin, xmax, ymax, zmax
    integer,intent(IN) :: level
    logical,intent(IN),optional :: interpolate
    character(len=*),intent(IN),optional :: prefix
    real(kind=DBL_KIND),dimension(OB_COORDS_SZ) ::  coordphys
    type(t_obRectPhys) :: rectPhysIn
    type(t_obRect),dimension(0:NPE-1) :: rect
    type(t_obRect) :: rectg
    type(t_ugArray),pointer :: u
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
#ifdef FORBIT_WRITEDATA
    return
#endif
    if ( present(prefix) ) then
       ThisPrefix = prefix
    else
       ThisPrefix = PrefixDefault
    endif
    if ( level > LevelMax ) then
       if (get_myrank() == PRIMARY_RANK) then
          print *, '*** uniform_write: level is larger than LevelMax: level =',level, ' LevelMax =', LevelMax
          print *, 'LevelFinest is truncated to ',LevelMax
          call flush(6)
       end if
    endif
    LevelFinest = max(min(level, LevelMax), Lmin)
    LevelResolution = level
    LevelFinalResolution = LevelFinest - max(Lmin-LevelResolution, 0)
    coordphys = (/ xmin, ymin, zmin, xmax, ymax, zmax /)
    call ob_assignCoordPhysToRectPhys( coordphys, rectPhysIn )
    call ug_make(rectPhysIn, u, rect, interpolate)
    call refile(u%arr)  ! write components of u into files separately
    call deallocateU(u)
    ! ----------
    ! downsizing
    ! ----------
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call rectGather(rect, rectg)
    call ug_DownSize(Lmin-LevelResolution, rectg)
    call makeCoordinates(x, y, z, rectg)
    ! ----------
    ! write data
    ! ----------
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call ug_write(x, y, z, Time(LevelFinest), Step(LevelFinest) )
    call ug_fld(x, y, z, Time(LevelFinest), Step(LevelFinest) )
    call clearTempFiles
    if ( get_myrank() == PRIMARY_RANK ) deallocate(x, y, z)
  end subroutine uniformgrid_write
  !---------------------------------------------------------------------
  ! make uniform grid
  !---------------------------------------------------------------------
  subroutine ug_make(rectPhysIn, u, rectReturn, interpolate)
    use grid, only : Lmin
    type(t_obRectPhys),intent(IN) :: rectPhysIn
    type(t_ugArray),pointer :: u
    type(t_obRect),dimension(0:NPE-1),intent(OUT) :: rectReturn
    logical,intent(IN),optional :: interpolate
    type(t_obRectPhys) :: rectPhys, rectCompBox, rectPhysAnd
    type(t_obRect),dimension(0:NPE-1) :: rect, rectExt, rectf, rectfExt
    type(t_ugArray),pointer :: uc, uf
    logical,dimension(:,:,:),pointer :: boolmap
    integer :: lev, m
    real(kind=DBL_KIND),dimension(MX:MZ) ::  h

    call ob_computationBoxOfRectPhys( rectCompBox )
    call ob_rectPhysAnd( rectPhysIn, rectCompBox, rectPhysAnd )
    h = CellWidth(:,LevelFinest)/2.d0
    call ob_rectPhysTrim( rectPhysAnd, h, rectPhys )
    ! -------------------
    ! physical variables
    ! -------------------
    myrank = get_myrank()
    call defineRects( rectPhys, Lmin, rect, rectExt )
    call allocateUfromRect( u,  rectExt )
    do lev = Lmin, LevelFinest
       if (myrank == PRIMARY_RANK) print *, 'level/LevelFinest', lev, LevelFinest
       call flush(6)
       call defineRects( rectPhys, lev, rect, rectExt)
       call allocateUfromRect( uc, rectExt )
       if ( allocatedU( u ) ) &
            allocate( boolmap(ARRAYSIZE3(u%arr)) )
       call ug_gather( uc, boolmap, rect )
       if ( allocatedU( uc ) ) then
          if (lev == Lmin) boolmap(:,:,:) = .true.
          do m = Mmin, Mmax
             where (boolmap) u%arr(:,:,:,m) = uc%arr(:,:,:,m)
          enddo
          deallocate( boolmap )
       endif
       call deallocateU( uc )

       if ( lev == LevelFinest ) exit

       call defineRects( rectPhys, lev+1, rectf, rectfExt )
       call allocateUfromRect( uf, rectfExt )
       call ug_interp( uf, u, rectf, rect, interpolate )
       call deallocateU( u )
       u => uf
       nullify( uf )
    enddo
    ! triming for output
    call defineRects( rectPhys, levelFinest, rectf, rectfExt )
    call allocateUfromRect( uf, rectf )
    call transferU( u, uf, rect, rectf )
    call deallocateU( u )
    u => uf
    nullify( uf )
    rectReturn = rectf
  end subroutine ug_make
  !---------------------------------------------------------------------
  ! define rects, (rect, rectg, rectExt), given by rectPhys and level
  !   rect(0:NPE-1) .... distributed rectangles without ghostcell
  !   rectExt .......... distributed rectangles with ghostcell
  !---------------------------------------------------------------------
  subroutine defineRects( rectPhys, level, rect, rectExt )
    type(t_obRectPhys),intent(IN) :: rectPhys
    integer,intent(IN) :: level
    type(t_obRect),dimension(0:NPE-1),intent(OUT) :: rect, rectExt
    type(t_obRect) :: rectE, rectg, rectCompBox
    integer :: rank, nrank, h(MX:MZ)
    nrank = 0
    h = (/1,1,1/)
    call ob_computationBoxOfRect( rectCompBox, level ) ! computation box
    call ob_RectPhys2RectOb(rectPhys, level, rectg ) ! rectg: node global rect
    call rectLocal(rectg, rect)                      ! rect:  node local rect
    do rank = 0, NPE-1          ! rectExt: extended rect for ghost cell
       if ( ob_definedRect(rect(rank)) ) then
          call ob_rectExtend(rect(rank), h, rectE)
          call ob_rectAnd(rectE, rectCompBox, rectExt(rank))
          nrank = nrank + 1
       else
          call ob_undefRect(rectExt(rank))
       endif
    enddo

  end subroutine defineRects
  !-------------------------------------------------------------------------
  ! rect を分割する
  !-------------------------------------------------------------------------
  subroutine rectLocal( rectg, rect )
    type(t_obRect),intent(IN) :: rectg
    type(t_obRect),dimension(0:NPE-1),intent(OUT) :: rect
    type(t_obPoint) :: pLg, pRg, pL, pR
    integer,dimension(0:NPE-1) :: ncells
    integer :: resid, rank
    integer(kind=LLONG_KIND) :: k, sz
    call ob_extractPointFromRect(pLg, rectg, 'L')
    call ob_extractPointFromRect(pRg, rectg, 'R')
    sz = pRg%p(MZ) - pLg%p(MZ) + 1
    ncells(:) = (sz)/(NPE)
    resid = mod( sz , (NPE) )
    ncells(0:resid-1) =  ncells(0:resid-1) + 1
    k = 0
    do rank = 0, (NPE)-1
       pL = pLg
       pR = pRg
       pL%p(MZ) = pLg%p(MZ) + k
       pR%p(MZ) = pLg%p(MZ) + k + ncells(rank) - 1
       if ( k <= sz -1 ) then
          call ob_assignPointToRect( pL, rect(rank), 'L' )
          call ob_assignPointToRect( pR, rect(rank), 'R' )
       else
          call ob_undefRect( rect(rank) )
       end if
       k = k + ncells(rank)
    enddo

  end subroutine rectLocal
  !---------------------------------------------------------------------
  ! gather data from every node
  ! INPUT:
  !   rect = Rectangle region in over-block coordinates
  !          Indexes specifining rectablel region of interest.
  !          This region is gathered into uGlobal
  !
  ! OUTPUT:
  !   uc      = gathered data
  !   boolmap = .ture. where cell has data.
  !---------------------------------------------------------------------
  subroutine ug_gather(uc, boolmap, rect)
    use mpilib
    use packarr
    type(t_obRect),dimension(0:NPE-1),intent(IN) :: rect
    type(t_ugArray),pointer :: uc
    logical,dimension(:,:,:),pointer :: boolmap
    type(t_obPoint) :: point, pL, pR
    integer :: gid, ranks, rankd, i, j, k, is, js, ks, ie, je, ke, ius, jus, kus, iue, jue, kue
    integer :: level
    integer(kind=LLONG_KIND) :: iob, job, kob
    integer :: isz = Imax-Imin+1, jsz = Jmax-Jmin+1, ksz = Kmax-Kmin+1
    logical :: ilast = .false. , jlast = .false. , klast = .false.
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: up, u
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: buf
    integer(kind=LLONG_KIND),dimension(MX:MZ) :: ijkob

    ! initialize of returned arguments

    myrank = get_myrank()
    if ( ob_definedRect( rect(myrank) ) ) then
       if ( .not. allocatedU(uc) ) print *,'*** error uc'
       u => uc%arr
       u(:,:,:,:) = 0.d0
       boolmap(:,:,:) = .false.
    endif

    call pkar_reset
    do rankd = 0, NPE-1
       if ( .not. ob_definedRect( rect(rankd) ) ) cycle
       call ob_extractLevelFromRect(level, rect(rankd) )         ! level
       call ob_extractPointFromRect(pL, rect(rankd), 'L')
       call ob_extractPointFromRect(pR, rect(rankd), 'R')
       do kob = pL%p(MZ), pR%p(MZ)+ksz, ksz
          klast = ( kob/ksz == pR%p(MZ)/ksz )
          if ( kob/ksz > pR%p(MZ)/ksz ) exit
          do job = pL%p(MY), pR%p(MY)+jsz, jsz
             jlast =  ( job/jsz == pR%p(MY)/jsz )
             if ( job/jsz > pR%p(MY)/jsz ) exit
             do iob = pL%p(MX), pR%p(MX)+isz, isz
                ilast = ( iob/isz == pR%p(MX)/isz )
                if ( iob/isz > pR%p(MX)/isz ) exit
                ijkob = (/iob, job, kob/)
                call ob_assignCoordToPoint(point, ijkob, level)
                call ob_getBlockFromPoint(point, gid, ranks, i,j,k)
                if ( gid == Undefi ) cycle
                ! region of block: is:ie,js:je,ks:ke
                is = Imin; js = Jmin; ks = Kmin
                ie = Imax; je = Jmax; ke = Kmax
                if ( iob == pL%p(MX) ) is = i
                if ( job == pL%p(MY) ) js = j
                if ( kob == pL%p(MZ) ) ks = k
                if ( ilast ) ie = mod(pR%p(MX), isz)
                if ( jlast ) je = mod(pR%p(MY), jsz)
                if ( klast ) ke = mod(pR%p(MZ), ksz)
                ! now transfer
                if ( ranks == myrank ) then
                   up => get_Up(gid)
                   allocate(buf(is:ie,js:je,ks:ke,Mmin:Mmax))
                   buf = up(is:ie,js:je,ks:ke,Mmin:Mmax)
                   call pkar_push(buf(PTF4(buf)), size(buf), kind(buf), rankd)
                   deallocate(buf)
                endif
                if ( rankd == myrank ) then
                   call pkar_recvlen((ie-is+1)*(je-js+1)*(ke-ks+1)*(Mmax-Mmin+1), kind(buf), ranks)
                endif
                if ( ilast ) exit
             enddo
             if ( jlast ) exit
          enddo
          if ( klast ) exit
       enddo
    enddo
    call pkar_sendrecv()
    do rankd = 0, NPE-1
       if ( rankd /= myrank ) cycle
       if ( .not. ob_definedRect( rect(rankd) ) ) cycle
       call ob_extractLevelFromRect(level, rect(rankd) )         ! level
       call ob_extractPointFromRect(pL, rect(rankd), 'L')
       call ob_extractPointFromRect(pR, rect(rankd), 'R')
       do kob = pL%p(MZ), pR%p(MZ)+ksz, ksz
          klast = ( kob/ksz == pR%p(MZ)/ksz )
          if ( kob/ksz > pR%p(MZ)/ksz ) exit
          do job = pL%p(MY), pR%p(MY)+jsz, jsz
             jlast =  ( job/jsz == pR%p(MY)/jsz )
             if ( job/jsz > pR%p(MY)/jsz ) exit
             do iob = pL%p(MX), pR%p(MX)+isz, isz
                ilast = ( iob/isz == pR%p(MX)/isz )
                if ( iob/isz > pR%p(MX)/isz ) exit
                ijkob = (/iob, job, kob/)
                call ob_assignCoordToPoint(point, ijkob, level)
                call ob_getBlockFromPoint(point, gid, ranks, i,j,k)
                if ( gid == Undefi ) cycle
                ! region of block: is:ie,js:je,ks:ke
                is = Imin; js = Jmin; ks = Kmin
                ie = Imax; je = Jmax; ke = Kmax
                if ( iob == pL%p(MX) ) is = i
                if ( job == pL%p(MY) ) js = j
                if ( kob == pL%p(MZ) ) ks = k
                if ( ilast ) ie = mod(pR%p(MX), isz)
                if ( jlast ) je = mod(pR%p(MY), jsz)
                if ( klast ) ke = mod(pR%p(MZ), ksz)
                ! region of u: ius:iue,jus:jue,kus:kue
                ius = -mod(pL%p(MX), isz) + iob - uc%offsets(MX) + is
                iue = ius + ie - is
                jus = -mod(pL%p(MY), jsz) + job - uc%offsets(MY) + js
                jue = jus + je - js
                kus = -mod(pL%p(MZ), ksz) + kob - uc%offsets(MZ) + ks
                kue = kus + ke - ks
                ! now assign
                allocate( buf(is:ie,js:je,ks:ke,Mmin:Mmax) )
                UNPACK_RECV4( buf, myrank, ranks, rankd )
                u(ius:iue,jus:jue,kus:kue,:) = buf
                deallocate(buf)
                boolmap(ius:iue,jus:jue,kus:kue) = .true.
                if ( ilast ) exit
             enddo
             if ( jlast ) exit
          enddo
          if ( klast ) exit
       enddo
    enddo
  end subroutine ug_gather
  !---------------------------------------------------------------------
  ! allocate array from local rect
  !---------------------------------------------------------------------
  subroutine allocateUfromRect( u, rect )
    type(t_ugArray),pointer :: u
    type(t_obRect),intent(IN),dimension(0:NPE-1) :: rect
    type(t_obPoint) :: pL, pR
    integer ::  nrank, rank, imax, jmax, kmax
    ! check
    if ( associated(u) ) then
       print *, '*** error: pointer u is already associated.', u%allocated, associated(u%arr)
       HALT
    endif

    !  define nrank
    nrank = 0
    do rank = 0, NPE-1
       if ( ob_definedRect(rect(rank)) ) nrank = nrank + 1
    enddo

    allocate( u )
    u%nrank = nrank
    myrank = get_myrank()
    if ( ob_definedRect(rect(myrank)) ) then
       call ob_extractPointFromRect(pL, rect(myrank), 'L')
       call ob_extractPointFromRect(pR, rect(myrank), 'R')
       imax = pR%p(MX)-pL%p(MX)
       jmax = pR%p(MY)-pL%p(MY)
       kmax = pR%p(MZ)-pL%p(MZ)
       allocate( u%arr( 0:imax, 0:jmax, 0:kmax, Mmin:Mmax) )
       u%arr = 1                ! initialize by safety value
       u%offsets(:) = pL%p(:)
       u%allocated = .true.
    else
       u%allocated = .false.
       nullify( u%arr )
    endif
  end subroutine allocateUfromRect
  !---------------------------------------------------------------------
  ! free arrary
  !---------------------------------------------------------------------
  subroutine deallocateU( u )
    type(t_ugArray),pointer :: u
    if (associated(u%arr)) deallocate( u%arr )
    nullify( u%arr )
    if (associated(u)) deallocate( u )
    nullify( u )
  end subroutine deallocateU
  !---------------------------------------------------------------------
  ! check allocated u
  !---------------------------------------------------------------------
  function allocatedU( u ) result(bool)
    type(t_ugArray),pointer :: u
    logical :: bool
    bool = u%allocated
  end function allocatedU
  !---------------------------------------------------------------------
  ! Write data in uniform grid (for PRIMARY_RANK)
  !---------------------------------------------------------------------
  subroutine ug_write(x,y,z,time,step)
    use systemcall
    use mpilib
    use string
    use io_util, only : wchar, read_env
    real(kind=DBL_KIND),dimension(:),intent(IN) :: x, y, z
    real(kind=DBL_KIND),intent(IN) :: time
    integer(kind=LLONG_KIND),intent(IN) :: step
    integer,parameter :: FH = 11, FHT = 12
    character(len=CHARLEN) :: fn, dir, suffix
    real,dimension(:,:),allocatable :: u4
    real(kind=DBL_KIND),dimension(:,:),allocatable :: u
    real,dimension(:),allocatable :: x4, y4, z4
    integer :: header, lenint, m, k, sx, sy, sz, recr, recw
    if( get_myrank() /= PRIMARY_RANK ) return
    lenint = kind(lenint)
    allocate(u(size(x),size(y)), u4(size(x),size(y)), x4(size(x)), y4(size(y)), z4(size(z)))
    x4 = x
    y4 = y
    z4 = z
    call read_env('DIR', dir)
    call read_env('SUFFIX', suffix)
    fn = get_filename(dir, ThisPrefix, suffix)
    call wchar(6,'uniform grid data = '//fn)

    ! grid size
    sx = size(x4); sy = size(y4); sz = size(z4)
    open(FH, file=trim(fn)//'a', form='unformatted')
    write(FH) sx, sy, sz, Mmax-Mmin+1, step, time  ! size of u
    write(FH) x4
    write(FH) y4
    write(FH) z4
    call flush(FH)
    close(FH)

    ! header
    open(FH, file=trim(fn)//'h', form='unformatted', access='direct', recl=kind(header)/RECL_UNIT)
    header = size(x4)*size(y4)*size(z4)*(Mmax-Mmin+1)*kind(u4)
    write(FH, rec=1) header
    call flush(FH)
    close(FH)
    ! write u
    open(FH, file=trim(fn)//'b', form='unformatted', access='direct', recl=size(u4)*kind(u4)/RECL_UNIT)
    recw = 1
    do m = Mmin, Mmax
       open(unit=FHT, file=get_tempfilename(m), form='unformatted', access='direct', recl=size(u)*kind(u)/RECL_UNIT)
       recr = 1
       do k = 1, size(z4)
          read(FHT, rec=recr) u
          recr = recr + 1
          u4 = u
          write(FH, rec=recw) u4
          recw = recw + 1
       end do
       call flush(FHT)
       close(FHT)
    end do
    call flush(FH)
    close(FH)

    ! header
    open(FH, file=trim(fn)//'f', form='unformatted', access='direct', recl=kind(header)/RECL_UNIT)
    write(FH, rec=1) header
    call flush(FH)
    close(FH)

    ! merge file
    call systemcall_command(trim(UNIX_CMD_CAT)//' '//trim(fn)//'a '//trim(fn)//'h '//trim(fn)//'b '//trim(fn)//'f >'//trim(fn))
    call systemcall_unlink(trim(fn)//'a')
    call systemcall_unlink(trim(fn)//'h')
    call systemcall_unlink(trim(fn)//'b')
    call systemcall_unlink(trim(fn)//'f')

    deallocate( u, u4, x4, y4, z4 )
  end subroutine ug_write
  !---------------------------------------------------------------------
  ! Write fld file (for PRIMARY_RANK)
  !---------------------------------------------------------------------
  subroutine ug_fld(x,y,z,time,step)
    use mpilib
    use io_util, only : wchar, read_env
    use string
    real(kind=DBL_KIND),dimension(:),intent(IN) :: x, y, z
    real(kind=DBL_KIND),intent(IN) :: time
    integer(kind=LLONG_KIND),intent(IN) :: step
    character(len=CHARLEN) :: fn, dir, suffix
    character(len=3),parameter :: Suffixfld = 'fld'
    integer,parameter :: FH = 11
    integer :: ni, nj, nk, nm, ndim, nspace, idepthv, idepthc, hbyte, fbyte, ivol, iskip, m
    if( get_myrank() /= PRIMARY_RANK ) return
    ni = size(x,1)
    nj = size(y,1)
    nk = size(z,1)
    nm = Mmax-Mmin+1
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
  ! make tmp file
  !-------------------------------------------------------------------------
  function get_tempfilename(rank) result(fn)
    use io_util, only : read_env
    use string, only : CHARLEN, concat, num2char
    integer,intent(IN) :: rank
    character(len=CHARLEN) :: fn, dir, suffix
    character(len=1),parameter :: dot = '.'
    call read_env('DIR', dir)
    call read_env('SUFFIX', suffix)
    fn = concat(dir, 'tempUG'//num2char(rank))
    fn = concat(fn, dot)
    fn = concat(fn, suffix)
  end function get_tempfilename
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  function get_tempfilename2(rank) result(fn)
    integer,intent(IN) :: rank
    character(len=CHARLEN) :: fn
    fn = get_tempfilename(rank)
    fn = trim(fn)//trim('d')
  end function get_tempfilename2
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
    fn = concat(fn,num2char(LevelFinalResolution))    ! /dir/st12000.0
    fn = concat(fn,dot)                    ! /dir/st12000.0.
    fn = concat(fn,suffix)                 ! /dir/st12000.0.d
  end function get_filename
  !---------------------------------------------------------------------
  ! interpolation of u
  !---------------------------------------------------------------------
  subroutine ug_interp(uf, uc, rectf, rect, interpolate)
    type(t_ugArray),pointer :: uf, uc
    type(t_obRect),dimension(0:NPE-1),intent(IN) :: rectf, rect
    logical,intent(IN),optional :: interpolate
    type(t_ugArray),pointer :: ut
    type(t_obRect),dimension(0:NPE-1) :: rectt
    integer :: level, rank

    ! define rectt (rectt has no ghostcell)
    do rank = 0, NPE-1
       if ( ob_definedRect( rect(rank) ) ) then
          call ob_extractLevelFromRect(level, rect(rank))
          call ob_rect2RectByLevel(rect(rank), rectt(rank), level+1)
       else
          call ob_undefRect(rectt(rank))
       endif
    enddo
    ! interp uc -> ut
    call overlapFixU( uc )
    call allocateUfromRect( ut,  rectt )
    call interpLocal( ut, uc, interpolate )
    call transferU( ut, uf, rectt, rectf )
    call deallocateU( ut )
  end subroutine ug_interp
  !---------------------------------------------------------------------
  ! interpolation of u
  !---------------------------------------------------------------------
  subroutine interpLocal(ufa, uca, interpolate)
    type(t_ugArray),pointer :: ufa, uca
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: uf, uc
    logical,intent(IN),optional :: interpolate
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: grad
    real(kind=DBL_KIND) :: di, dj, dk, hf, hc, hfi, hci, dvc, dvf
    integer :: m, if, jf, kf, ic, jc, kc
    integer :: ics, ice, jcs, jce, kcs, kce, ifs, ife, jfs, jfe, kfs, kfe
    logical :: bool_interpolate
    ! minmod limiter
#define MINMOD(a,b) sign(1.d0,(a))*max(0.d0,min(abs(a),sign(1.d0,(a))*(b)))

    if (.not. allocatedU( uca ) ) return

    bool_interpolate = .false.
    if ( present( interpolate) ) then
       bool_interpolate = interpolate
    endif

    uf => ufa%arr
    uc => uca%arr

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

    call conv_u2w_smallMem(uc, dvc, U2W)
    allocate( grad(ARRAYSIZE3(uc),MX:MZ) )
    do m = Mmin, Mmax
       grad(:,:,:,:) = 0.d0
       if ( bool_interpolate ) then
          do ic = ics+1, ice-1
             grad(ic,:,:,MX) = ( MINMOD(uc(ic+1,:,:,m)-uc(ic,:,:,m), uc(ic,:,:,m)-uc(ic-1,:,:,m)) )*hci
          enddo
          do jc = jcs+1, jce-1
             grad(:,jc,:,MY) = ( MINMOD(uc(:,jc+1,:,m)-uc(:,jc,:,m), uc(:,jc,:,m)-uc(:,jc-1,:,m)) )*hci
          enddo
          do kc = kcs+1, kce-1
             grad(:,:,kc,MZ) = ( MINMOD(uc(:,:,kc+1,m)-uc(:,:,kc,m), uc(:,:,kc,m)-uc(:,:,kc-1,m)) )*hci
          enddo
          if ( ics /= ice ) then
             grad(ics,:,:,MX) = ( uc(ics+1,:,:,m)-uc(ics,:,:,m) )*hci
             grad(ice,:,:,MX) = ( uc(ice,:,:,m)-uc(ice-1,:,:,m) )*hci
          endif
          if ( jcs /= jce ) then
             grad(:,jcs,:,MY) = ( uc(:,jcs+1,:,m)-uc(:,jcs,:,m) )*hci
             grad(:,jce,:,MY) = ( uc(:,jce,:,m)-uc(:,jce-1,:,m) )*hci
          endif
          if ( kcs /= kce ) then
             grad(:,:,kcs,MZ) = ( uc(:,:,kcs+1,m)-uc(:,:,kcs,m) )*hci
             grad(:,:,kce,MZ) = ( uc(:,:,kce,m)-uc(:,:,kce-1,m) )*hci
          endif
       endif
       do kf = kfs, kfe
          do jf = jfs, jfe
             do if = ifs, ife
                ic = IJKC(if+ufa%offsets(MX),0)-uca%offsets(MX)
                jc = IJKC(jf+ufa%offsets(MY),0)-uca%offsets(MY)
                kc = IJKC(kf+ufa%offsets(MZ),0)-uca%offsets(MZ)
                di = ( modulo(if+ufa%offsets(MX),int(2,LLONG_KIND)) - 0.5d0 )*hf
                dj = ( modulo(jf+ufa%offsets(MY),int(2,LLONG_KIND)) - 0.5d0 )*hf
                dk = ( modulo(kf+ufa%offsets(MZ),int(2,LLONG_KIND)) - 0.5d0 )*hf
                uf(if,jf,kf,m) = uc(ic,jc,kc,m) + grad(ic,jc,kc,MX)*di+grad(ic,jc,kc,MY)*dj+grad(ic,jc,kc,MZ)*dk
             enddo
          enddo
       enddo
    enddo
    deallocate( grad )
    call conv_u2w_smallMem(uf, dvf, W2U)
#undef MINMOD
  end subroutine interpLocal
  !---------------------------------------------------------------------
  ! transfer u between ranks
  !---------------------------------------------------------------------
  subroutine transferU( us, ud, rects, rectd )
    type(t_ugArray),pointer :: us, ud
    type(t_obRect),dimension(0:NPE-1),intent(IN) :: rects, rectd
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: buf
    integer :: ranks, rankd, test
    type(t_obRect) :: rectand, recttest
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ) :: coord
    integer :: iss, ies, jss, jes, kss, kes ! index of send array
    integer :: isd, ied, jsd, jed, ksd, ked ! index of recv array
#define SZS iss:ies,jss:jes,kss:kes,Mmin:Mmax
#define SZD isd:ied,jsd:jed,ksd:ked,Mmin:Mmax
    myrank = get_myrank()
    ! check consistency of u and rect
    if (ob_definedRect( rects(myrank) ) ) then
       call u2rect( us, rects(myrank)%level, recttest )
       test = ob_rectsComp( rects(myrank), recttest )
       if ( test /= OB_RECT_EQUIV .and. test /= OB_RECT_INNER ) then
          print *, '*** error in transferU: us and rects are not consistent'
       end if
    endif
    if (ob_definedRect( rectd(myrank) ) ) then
       call u2rect( ud, rectd(myrank)%level, recttest )
       test = ob_rectsComp( rectd(myrank), recttest )
       if ( test /= OB_RECT_EQUIV .and. test /= OB_RECT_INNER ) then
          print *, '*** error in transferU: ud and rectd are not consistent'
       end if
    end if

    do ranks = 0, NPE-1
       if ( .not. ob_definedRect( rects(ranks) ) ) cycle
       do rankd = 0, NPE-1
          if ( .not. ob_definedRect( rectd(rankd) ) ) cycle
          call ob_rectAnd( rects(ranks), rectd(rankd), rectand )
          if ( ob_testRect(rectand) /= OB_RECT_VALID ) cycle
          call ob_extractCoordFromRect(coord, rectand)

          iss = coord(0)-us%offsets(MX); jss = coord(1)-us%offsets(MY); kss = coord(2)-us%offsets(MZ);
          ies = coord(3)-us%offsets(MX); jes = coord(4)-us%offsets(MY); kes = coord(5)-us%offsets(MZ);
          isd = coord(0)-ud%offsets(MX); jsd = coord(1)-ud%offsets(MY); ksd = coord(2)-ud%offsets(MZ);
          ied = coord(3)-ud%offsets(MX); jed = coord(4)-ud%offsets(MY); ked = coord(5)-ud%offsets(MZ);

          if ( myrank == ranks .and. myrank == rankd ) then
             allocate( buf(SZS) )
             buf = us%arr(SZS)
             ud%arr(SZD) = buf
             deallocate( buf )
!!$             ud%arr(SZD) = us%arr(SZS)
          else if ( myrank == ranks ) then
             allocate( buf(SZS) )
             buf = us%arr(SZS)
             call mpi_send(buf, size(buf), MPI_DOUBLE_PRECISION, rankd, 1, MPI_COMM_WORLD, ierr )
             deallocate( buf )
          else if ( myrank == rankd ) then
             allocate( buf(SZD) )
             call mpi_recv(buf, size(buf), MPI_DOUBLE_PRECISION, ranks, 1, MPI_COMM_WORLD, status, ierr )
             ud%arr(SZD) = buf
             deallocate( buf )
          end if
       enddo
    enddo
#undef SZS
#undef SZD
  end subroutine transferU
  !---------------------------------------------------------------------
  ! make rect form u
  !---------------------------------------------------------------------
  subroutine u2rect( u, level, rect )
    type(t_ugArray),intent(IN) :: u
    integer,intent(IN) :: level
    type(t_obRect),intent(OUT) :: rect
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ) :: coords

    coords(MX:MZ) = u%offsets(MX:MZ)
    coords(MZ+1:OB_COORDS_MAX) = u%offsets(MX:MZ) + (/ubound(u%arr,1),ubound(u%arr,2),ubound(u%arr,3)/)

    call ob_assignCoordToRect( coords, rect, level )

  end subroutine u2rect
  !---------------------------------------------------------------------
  ! overlap fix of u
  !---------------------------------------------------------------------
  subroutine overlapFixU( ua )
    type(t_ugArray),pointer :: ua
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: bufs, bufd
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: rankR, rankL, nrank
    integer,parameter :: TAG = 1

    myrank = get_myrank()
    if ( .not. allocatedU(ua) ) return
    u => ua%arr
    nrank = ua%nrank

    allocate(bufs(lbound(u,1):ubound(u,1), lbound(u,2):ubound(u,2), 1, Mmin:Mmax))
    allocate(bufd(lbound(u,1):ubound(u,1), lbound(u,2):ubound(u,2), 1, Mmin:Mmax))
    rankR = myrank+1
    rankL = myrank-1
    if ( myrank == nrank - 1 )    rankR = MPI_PROC_NULL
    if ( myrank == PRIMARY_RANK ) rankL = MPI_PROC_NULL
    ! overlap fix of u (left ghost cells)
    bufs = u(:,:,ubound(u,3)-1:ubound(u,3)-1,:)
    call mpi_sendrecv( &
         bufs, size(bufs), MPI_DOUBLE_PRECISION, rankR, TAG, &
         bufd, size(bufd), MPI_DOUBLE_PRECISION, rankL, TAG, MPI_COMM_WORLD, status, ierr )
    if (myrank /= PRIMARY_RANK) u(:,:,lbound(u,3):lbound(u,3),:) = bufd
    ! overlap fix of u (right ghost cells)
    bufs = u(:,:,lbound(u,3)+1:lbound(u,3)+1,:)
    call mpi_sendrecv( &
         bufs, size(bufs), MPI_DOUBLE_PRECISION, rankL, TAG, &
         bufd, size(bufd), MPI_DOUBLE_PRECISION, rankR, TAG, MPI_COMM_WORLD, status, ierr )
    if (myrank /= nrank-1) u(:,:,ubound(u,3):ubound(u,3),:) = bufd
    deallocate( bufs, bufd )
  end subroutine overlapFixU
  !---------------------------------------------------------------------
  ! convert between u to w by small memory
  !---------------------------------------------------------------------
  subroutine conv_u2w_smallMem(u, dv, code)
    use eos, only: conv_u2w, conv_w2u
    real(kind=DBL_KIND),dimension(:,:,:,:),intent(INOUT) :: u
    real(kind=DBL_KIND),intent(IN) :: dv
    integer,intent(IN) :: code
    integer,parameter :: szi = 8, szj = 8, szk = 8 ! size of buffer
    real(kind=DBL_KIND),dimension(0:szi-1,0:szj-1,0:szk-1, Mmin:Mmax) :: buffer
    integer :: is, js, ks, ie, je, ke, ii, jj, kk, i, j, k
    is = lbound(u,1) ; js = lbound(u,2) ; ks = lbound(u,3)
    ie = ubound(u,1) ; je = ubound(u,2) ; ke = ubound(u,3)
    do k = ks, ke, szk
       do j = js, je, szj
          do i = is, ie, szi
             buffer = 1.d0      ! initialize by safety value
             ii = min(i+szi-1, ie)
             jj = min(j+szj-1, je)
             kk = min(k+szk-1, ke)
             buffer(0:ii-i, 0:jj-j, 0:kk-k, :) = u(i:ii, j:jj, k:kk, :)
             if ( code == U2W ) then
                call conv_u2w(buffer, dv)
             elseif ( code == W2U ) then
                call conv_w2u(buffer, dv)
             endif
             u(i:ii, j:jj, k:kk, :) = buffer(0:ii-i, 0:jj-j, 0:kk-k, :)
          enddo
       enddo
    enddo
  end subroutine conv_u2w_smallMem
  !---------------------------------------------------------------------
  ! Down-size data of disk
  ! dlev : level of downsize
  ! for PRIMARY_RANK
  !---------------------------------------------------------------------
  subroutine ug_DownSize(dlev, rectg)
    use string, only : CHARLEN, concat
    integer,intent(IN) :: dlev
    type(t_obRect),intent(INOUT) :: rectg
    type(t_obRect) :: rectc
    type(t_obPoint) :: pL, pR
    character(len=CHARLEN) :: fn, fn2
    integer(kind=LLONG_KIND) :: ics, ice, jcs, jce, kcs, kce, iob, job, kob
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ) :: coords
    integer :: szi, szj, szk, level, if, jf, kf, ic, jc, kc, m, lev, rec, rec2
    integer,parameter :: FH=11, FH2 = 12
    real(kind=DBL_KIND),dimension(:,:,:),allocatable :: uf
    real(kind=DBL_KIND),dimension(:,:),allocatable :: uc, ufs
    if (dlev <= 0) return
    myrank = get_myrank()
    if ( myrank /= PRIMARY_RANK ) return
    call conv_u2w_tempFiles(U2W, rectg)
    ! ------------------------------------------
    ! downsize U in i, j, k-directions
    ! -----------------------------------------
    call clearTempFiles2
    do level = 1, dlev
       call ob_extractPointFromRect(pL, rectg, 'L')
       call ob_extractPointFromRect(pR, rectg, 'R')
       call ob_extractLevelFromRect(lev, rectg)
       szi = pR%p(MX) - pL%p(MX) + 1
       szj = pR%p(MY) - pL%p(MY) + 1
       szk = pR%p(MZ) - pL%p(MZ) + 1
       ics = pL%p(MX)/2 + mod(pL%p(MX),2)
       jcs = pL%p(MY)/2 + mod(pL%p(MY),2)
       kcs = pL%p(MZ)/2 + mod(pL%p(MZ),2)
       ice = pR%p(MX)/2 + mod(pR%p(MX),2)-1
       jce = pR%p(MY)/2 + mod(pR%p(MY),2)-1
       kce = pR%p(MZ)/2 + mod(pR%p(MZ),2)-1
       coords = (/ ics, jcs, kcs, ice, jce, kce /)
       call ob_assignCoordToRect( coords, rectc, lev-1 )
       allocate(uf(0:szi-1, 0:szj-1, 0:1))
       allocate(ufs(0:szi-1, 0:szj-1))
       allocate(uc(0:ice-ics, 0:jce-jcs))
       do m = Mmin, Mmax
          fn =  get_tempfilename(m)
          fn2 = get_tempfilename2(m)
          open(unit=FH,  file=fn,  form='unformatted', access='direct', recl=size(ufs)*kind(ufs)/RECL_UNIT)
          open(unit=FH2, file=fn2, form='unformatted', access='direct', recl=size(uc)*kind(uc)/RECL_UNIT)
          rec = 1
          rec2 = 1
          if ( mod( pL%p(MZ), 2 ) /= 0 ) then
             rec = rec+1
          end if
          do kob = kcs, kce
             read(FH, rec=rec) ufs  ; rec=rec+1
             uf(:,:,0) = ufs
             read(FH, rec=rec) ufs  ; rec=rec+1
             uf(:,:,1) = ufs
             do job = jcs, jce
                do iob = ics, ice
                   ! (iob, job, kob) = overBlockCoordinates
                   ic = iob - ics
                   jc = job - jcs
                   kc = kob - kcs
                   if = iob * 2 - pL%p(MX)
                   jf = job * 2 - pL%p(MY)
                   kf = kob * 2 - pL%p(MZ)
                   uc(ic,jc) = &
                        (uf(if, jf,0) + uf(if+1,jf+1,0) + &
                        uf(if+1,jf,0) + uf(if,  jf+1,0) + &
                        uf(if,  jf,1) + uf(if+1,jf+1,1) + &
                        uf(if+1,jf,1) + uf(if,  jf+1,1) )/8.d0
                end do
             end do
             write(FH2, rec=rec2) uc ; rec2=rec2+1
          end do
          call flush(FH)
          call flush(FH2)
          close(unit=FH)
          close(unit=FH2)
       end do
       deallocate(uf, uc, ufs)
       rectg = rectc
       call swapTempFiles
    end do
    call conv_u2w_tempFiles(W2U, rectg)
  end subroutine ug_DownSize
  !---------------------------------------------------------------------
  ! refile u. write data in tempfile gathering in k-direction.
  !---------------------------------------------------------------------
  subroutine refile(u)
    use string, only : CHARLEN
    real(KIND=DBL_KIND),dimension(:,:,:,:),pointer :: u  ! IN
    integer :: m, rank, rec(Mmin:Mmax), k, sz
    integer,parameter :: FH = 11
    myrank = get_myrank()
    call clearTempFiles
    rec(:) = 1
    sz = size(u,1)*size(u,2)
    do rank = 0, NPE-1          ! each component is written in each file separately.
       if (rank == myrank .and. associated(u)) then
          do m = Mmin, Mmax
             open(unit=FH, file=get_tempfilename(m), form='unformatted', access='direct', recl=sz*kind(u)/RECL_UNIT)
             do k=lbound(u,3),ubound(u,3)
                write(FH, rec=rec(m)) u(:,:,k,m)
                rec(m) = rec(m) + 1
             enddo
             call flush(FH)
             close(FH)
          end do
       end if
       call mpi_bcast(rec, size(rec), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr) ! update rec
    end do
  end subroutine refile
  !---------------------------------------------------------------------
  ! clear temp files
  !---------------------------------------------------------------------
  subroutine clearTempFiles
    use systemcall, only : systemcall_unlink
    integer :: m
    character(len=CHARLEN) :: fn
    logical :: exist
    myrank = get_myrank()
    if ( myrank /= PRIMARY_RANK ) return
    do m = Mmin, Mmax
       fn = get_tempfilename(m)
       inquire(file=fn, exist=exist)
       if (exist) call systemcall_unlink(fn)
    end do
  end subroutine clearTempFiles
  !---------------------------------------------------------------------
  ! clear temp files
  !---------------------------------------------------------------------
  subroutine clearTempFiles2
    use systemcall, only : systemcall_unlink
    integer :: m
    character(len=CHARLEN) :: fn
    logical :: exist
    myrank = get_myrank()
    if ( myrank /= PRIMARY_RANK ) return
    do m = Mmin, Mmax
       fn = get_tempfilename2(m)
       inquire(file=fn, exist=exist)
       if (exist) call systemcall_unlink(fn)
    end do
  end subroutine clearTempFiles2
  !---------------------------------------------------------------------
  ! Make coordinates in primary rank
  !---------------------------------------------------------------------
  subroutine makeCoordinates(x, y, z, rectg)
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z ! OUT
    type(t_obRect),intent(IN) :: rectg
    type(t_obPoint) :: pL, pR
    type(t_obRectPhys) :: rectCompBox
    type(t_obPointPhys) ::  pcL
    integer :: imax, jmax, kmax, level, i, j, k
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    if( get_myrank() /= PRIMARY_RANK ) return
    call ob_extractPointFromRect(pL, rectg, 'L')
    call ob_extractPointFromRect(pR, rectg, 'R')
    ! define coordinates, x, y, z
    imax = pR%p(MX)-pL%p(MX)
    jmax = pR%p(MY)-pL%p(MY)
    kmax = pR%p(MZ)-pL%p(MZ)
    allocate( x(0:imax), y(0:jmax), z(0:kmax) )
    call ob_extractLevelFromRect(level, rectg)
    call ob_computationBoxOfRectPhys( rectCompBox )
    call ob_extractPointPhysFromRectPhys(pcL, rectCompBox, 'L')
    if (level >= Lmin) then
       h(:) = CellWidth(:,level)
    else
       h(:) = CellWidth(:,Lmin) * 2.d0 ** (Lmin-level)
    end if
    do i = 0, imax
       x(i) = (i + pL%p(MX) + 0.5d0) * h(MX) + pcL%p(MX)
    end do
    do j = 0, jmax
       y(j) = (j + pL%p(MY) + 0.5d0) * h(MY) + pcL%p(MY)
    end do
    do k = 0, kmax
       z(k) = (k + pL%p(MZ) + 0.5d0) * h(MZ) + pcL%p(MZ)
    end do

  end subroutine makeCoordinates
  !---------------------------------------------------------------------
  ! make rect in global from rects distributed in nodes
  !---------------------------------------------------------------------
  subroutine rectGather(rect, rectg)
    type(t_obRect),dimension(0:NPE-1),intent(IN) :: rect
    type(t_obRect),intent(OUT) :: rectg
    type(t_obRect) :: rectswap
    integer :: rank
    rectg = rect(PRIMARY_RANK)
    do rank = 1, (NPE) - 1
       if (.not. ob_definedRect(rect(rank))) cycle
       call ob_rectOr(rect(rank), rectg, rectswap)
       rectg = rectswap
    end do
  end subroutine rectGather
  !---------------------------------------------------------------------
  ! convert u and w in temp files
  !---------------------------------------------------------------------
  subroutine conv_u2w_tempFiles(code, rect)
    integer,intent(IN) :: code
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: pL, pR
    real(kind=DBL_KIND),dimension(:,:),allocatable :: buf
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: u
    real(kind=DBL_KIND),parameter :: dv = 1.d0
    integer :: m
    integer(kind=LLONG_KIND) :: kob
    integer,parameter :: FH=11, FH2 = FH + Mmax-Mmin + 1
    integer :: recr(Mmin:Mmax), recw(Mmin:Mmax)
    call clearTempFiles2
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    allocate( buf(0:pR%p(MX)-pL%p(MX), 0:pR%p(MY)-pL%p(MY)) )
    allocate(   u(0:pR%p(MX)-pL%p(MX), 0:pR%p(MY)-pL%p(MY), 1, Mmin:Mmax) )
    do m = Mmin, Mmax
       open(unit=FH+m, file=get_tempfilename(m), form='unformatted', access='direct', recl=size(buf)*kind(buf)/RECL_UNIT)
       open(unit=FH2+m, file=get_tempfilename2(m), form='unformatted', access='direct', recl=size(buf)*kind(buf)/RECL_UNIT)
    end do
    recr(:) = 1
    recw(:) = 1
    do kob = pL%p(MZ), pR%p(MZ)
       do m = Mmin, Mmax
          read(unit=FH+m, rec=recr(m)) buf
          recr(m) = recr(m) + 1
          u(:,:,1,m) = buf
       end do
       call conv_u2w_smallMem(u, dv, code)
       do m = Mmin, Mmax
          buf = u(:,:,1,m)
          write(unit=FH2+m, rec=recw(m)) buf
          recw(m) = recw(m) + 1
       end do
    end do
    do m = Mmin, Mmax
       call flush(FH+m)
       call flush(FH2+m)
       close(FH+m)
       close(FH2+m)
    end do
    deallocate( buf, u )
    call swapTempFiles
  end subroutine conv_u2w_tempFiles
  !---------------------------------------------------------------------
  ! swap temp files
  !---------------------------------------------------------------------
  subroutine swapTempFiles
    use systemcall, only : systemcall_unlink, systemcall_rename
    use string, only : CHARLEN
    integer :: m
    character(len=CHARLEN) :: fn, fn2
    logical :: exist, exist2
    do m = Mmin, Mmax
       fn = get_tempfilename(m)
       fn2 = get_tempfilename2(m)
       inquire(file=fn, exist=exist)
       inquire(file=fn2, exist=exist2)
       if ( exist .and. exist2 ) then
          call systemcall_unlink(fn)
       end if
       if ( exist2 ) then
          call systemcall_rename(fn2, fn)
       end if
    end do
  end subroutine swapTempFiles
  !---------------------------------------------------------------------
  subroutine printrect(rect)
    type(t_obRect),dimension(0:NPE-1),intent(IN) :: rect
    integer :: rank
    integer(kind=LLONG_KIND) :: is, ie, js, je, ks, ke
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ) :: coords
    is = HUGE(is)
    js = HUGE(js)
    ks = HUGE(ks)
    ie = -HUGE(ie)
    je = -HUGE(je)
    ke = -HUGE(ke)
    do rank = 0, NPE-1
       if ( .not. ob_definedRect( rect(rank) )) cycle
       call ob_extractCoordFromRect(coords, rect(rank))
       is = min(is,coords(0))
       js = min(js,coords(1))
       ks = min(ks,coords(2))
       ie = max(ie,coords(3))
       je = max(je,coords(4))
       ke = max(ke,coords(5))
    enddo
   if (get_myrank() == PRIMARY_RANK) print *, 'rect', is,ie, js,je, ks,ke
  end subroutine printrect
  !---------------------------------------------------------------------
  subroutine printUwidth( u )
    type(t_ugArray),pointer :: u
    integer(kind=LLONG_KIND) :: is, ie, js, je, ks, ke
    integer :: buf(MX:MZ)
    integer :: nrank, szk, szkr
    myrank = get_myrank()
    if (myrank == PRIMARY_RANK) then
       nrank = u%nrank
    end if
    call mpi_bcast(nrank, 1, MPI_INTEGER, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    if ( myrank == nrank-1) then
       buf(MX) = u%offsets(MX) + ubound(u%arr,1)
       buf(MY) = u%offsets(MY) + ubound(u%arr,2)
       buf(MZ) = u%offsets(MZ) + ubound(u%arr,3)
       call mpi_send(buf(MX), size(buf), MPI_INTEGER, PRIMARY_RANK, 1, MPI_COMM_WORLD, ierr )
    endif
    if (myrank == PRIMARY_RANK) then
       call mpi_recv(buf(MX), size(buf), MPI_INTEGER, nrank-1, 1, MPI_COMM_WORLD, status, ierr )
       ie = buf(MX)
       je = buf(MY)
       ke = buf(MZ)
       is = u%offsets(MX)
       js = u%offsets(MY)
       ks = u%offsets(MZ)
       print *, 'u size ', is,ie,js,je,ks,ke
    endif

    if ( allocatedU(u) ) then
       szk = size(u%arr,3)
    else
       szk = 0
    endif
    call mpi_allreduce(szk, szkr, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
    if (myrank == PRIMARY_RANK) &
         print *, 'u size*', size(u%arr,1),size(u%arr,2), szkr-2*nrank+2
  end subroutine printUwidth
  !---------------------------------------------------------------------
  subroutine checkbound(u, is, ie, js,je,ks,ke )
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer,intent(IN) :: is, ie, js,je,ks,ke

    if (&
         is < lbound(u,1) .or. js < lbound(u,2) .or. ks < lbound(u,3) .or. &
         ie > ubound(u,1) .or. je > ubound(u,2) .or. ke > ubound(u,3)  ) then
       print *, '**** error'
    endif

  end subroutine checkbound
end module uniformgrid
!---------------------------------------------------------------------

