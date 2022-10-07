module uniformgrid
  use grid, only : Undefi, Imin, Jmin, Kmin, Imax, Jmax, Kmax, Mmin, Mmax, get_Up, get_Xp, get_Yp, get_Zp
  use overBlockCoordinates
  use mpilib
  use string, only : CHARLEN
  implicit none
  private
  type(t_obRectPhys),save :: rectCompBoxPhys
  type(t_obRect),save :: rectCompBox
  character(len=2),parameter :: PrefixDefault = 'ug' 
  character(len=CHARLEN) :: ThisPrefix
  integer,save :: LevelFinest = Undefi 
  integer,save :: LevelResolution = Undefi 
  integer,parameter :: U2W = 1
  integer,parameter :: W2U = 2
  type t_ugArray
     real(kind=8),dimension(:,:,:,:),pointer :: arr => null()
     integer(kind=8),dimension(0:2) :: offsets
     integer :: nrank = MPI_PROC_NULL
     logical :: allocated = .false.
  end type t_ugArray
  public :: uniformgrid_write
contains
  subroutine uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate, prefix)
    use grid, only : Step, Time
    use string, only : num2char
    use io_util, only : print_msg
    real(KIND=8),intent(IN) :: xmin, ymin, zmin, xmax, ymax, zmax
    integer,intent(IN) :: level
    logical,intent(IN),optional :: interpolate
    character(len=*),intent(IN),optional :: prefix
    real(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX) :: coordphys
    type(t_obRectPhys) :: rectPhysIn
    type(t_obRect),dimension(0:400 -1) :: rect
    type(t_obRect) :: rectg
    type(t_ugArray),save,pointer :: u => null()
    real(kind=8),dimension(:),pointer :: x, y, z
    if ( present(prefix) ) then
       ThisPrefix = prefix
    else
       ThisPrefix = PrefixDefault
    endif
    LevelFinest = max(min(level, LevelMax), Lmin)
    LevelResolution = min(level, LevelMax)
    if ( level > LevelMax ) then
       call print_msg('*** uniform_write: level is larger than LevelMax: level = '// trim(num2char(level))//' LevelMax = ' // trim(&
&num2char(LevelMax)))
       call print_msg('*** LevelFinest is truncated to '//trim(num2char(LevelFinest)))
       call print_msg('*** LevelResolution is truncated to '//trim(num2char(LevelResolution)))
    endif
    coordphys = (/ xmin, ymin, zmin, xmax, ymax, zmax /)
    call ob_assignCoordPhysToRectPhys( coordphys, rectPhysIn )
    call ug_make(rectPhysIn, u, rect, interpolate)
    call refile(u%arr) 
    call deallocateU(u)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call rectGather(rect, rectg)
    call ug_DownSize(Lmin-LevelResolution, rectg)
    call makeCoordinates(x, y, z, rectg)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call ug_write(x, y, z, Time(LevelFinest), Step(LevelFinest) )
    if ( get_myrank() == 0 ) &
         print '(3(A,I4))', 'Ni x Nj x Nk =', size(x), ' x ', size(y), ' x ', size(z)
    call clearTempFiles
    if ( get_myrank() == 0 ) deallocate(x, y, z)
  end subroutine uniformgrid_write
  subroutine ug_make(rectPhysIn, u, rectReturn, interpolate)
    use grid, only : Lmin
    type(t_obRectPhys),intent(IN) :: rectPhysIn
    type(t_ugArray),pointer :: u
    type(t_obRect),dimension(0:400 -1),intent(OUT) :: rectReturn
    logical,intent(IN),optional :: interpolate
    type(t_obRectPhys) :: rectPhys, rectCompBox, rectPhysAnd
    type(t_obRect),dimension(0:400 -1) :: rect, rectExt, rectf, rectfExt
    type(t_ugArray),save,pointer :: uc => null(), uf => null()
    logical,dimension(:,:,:),pointer :: boolmap
    integer :: lev, m
    real(kind=8),dimension(0:2) :: h
    call ob_computationBoxOfRectPhys( rectCompBox )
    call ob_rectPhysAnd( rectPhysIn, rectCompBox, rectPhysAnd )
    h = CellWidth(:,LevelFinest)/2.d0
    call ob_rectPhysTrim( rectPhysAnd, h, rectPhys )
    myrank = get_myrank()
    call defineRects( rectPhys, Lmin, rect, rectExt )
    call allocateUfromRect( u, rectExt )
    do lev = Lmin, LevelFinest
       if (myrank == 0) print *, 'level/LevelFinest', lev, LevelFinest
       call flush(6)
       call defineRects( rectPhys, lev, rect, rectExt)
       call allocateUfromRect( uc, rectExt )
       if ( allocatedU( u ) ) &
            allocate( boolmap(lbound(u%arr,1):ubound(u%arr,1),lbound(u%arr,2):ubound(u%arr,2),lbound(u%arr,3):ubound(u%arr,3)) )
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
    call defineRects( rectPhys, levelFinest, rectf, rectfExt )
    call allocateUfromRect( uf, rectf )
    call transferU( u, uf, rect, rectf )
    call deallocateU( u )
    u => uf
    nullify( uf )
    rectReturn = rectf
  end subroutine ug_make
  subroutine defineRects( rectPhys, level, rect, rectExt )
    type(t_obRectPhys),intent(IN) :: rectPhys
    integer,intent(IN) :: level
    type(t_obRect),dimension(0:400 -1),intent(OUT) :: rect, rectExt
    type(t_obRect) :: rectE, rectg, rectCompBox
    integer :: rank, nrank, h(0:2)
    nrank = 0
    h = (/1,1,1/)
    call ob_computationBoxOfRect( rectCompBox, level ) 
    call ob_RectPhys2RectOb(rectPhys, level, rectg ) 
    call rectLocal(rectg, rect) 
    do rank = 0, 400 -1 
       if ( ob_definedRect(rect(rank)) ) then
          call ob_rectExtend(rect(rank), h, rectE)
          call ob_rectAnd(rectE, rectCompBox, rectExt(rank))
          nrank = nrank + 1
       else
          call ob_undefRect(rectExt(rank))
       endif
    enddo
  end subroutine defineRects
  subroutine rectLocal( rectg, rect )
    type(t_obRect),intent(IN) :: rectg
    type(t_obRect),dimension(0:400 -1),intent(OUT) :: rect
    type(t_obPoint) :: pLg, pRg, pL, pR
    integer,dimension(0:400 -1) :: ncells
    integer :: resid, rank
    integer(kind=8) :: k, sz
    call ob_extractPointFromRect(pLg, rectg, 'L')
    call ob_extractPointFromRect(pRg, rectg, 'R')
    sz = pRg%p(2) - pLg%p(2) + 1
    ncells(:) = (sz)/(400)
    resid = mod( sz , (400) )
    ncells(0:resid-1) = ncells(0:resid-1) + 1
    k = 0
    do rank = 0, (400)-1
       pL = pLg
       pR = pRg
       pL%p(2) = pLg%p(2) + k
       pR%p(2) = pLg%p(2) + k + ncells(rank) - 1
       if ( k <= sz -1 ) then
          call ob_assignPointToRect( pL, rect(rank), 'L' )
          call ob_assignPointToRect( pR, rect(rank), 'R' )
       else
          call ob_undefRect( rect(rank) )
       end if
       k = k + ncells(rank)
    enddo
  end subroutine rectLocal
  subroutine ug_gather(uc, boolmap, rect)
    use mpilib
    use packarr
    type(t_obRect),dimension(0:400 -1),intent(IN) :: rect
    type(t_ugArray),pointer :: uc
    logical,dimension(:,:,:),pointer :: boolmap
    type(t_obPoint) :: point, pL, pR
    integer :: gid, ranks, rankd, i, j, k, is, js, ks, ie, je, ke, ius, jus, kus, iue, jue, kue
    integer :: level
    integer(kind=8) :: iob, job, kob
    integer :: isz = Imax-Imin+1, jsz = Jmax-Jmin+1, ksz = Kmax-Kmin+1
    logical :: ilast = .false. , jlast = .false. , klast = .false.
    real(kind=8),dimension(:,:,:,:),pointer :: up, u
    real(kind=8),dimension(:,:,:,:),allocatable :: buf
    integer(kind=8),dimension(0:2) :: ijkob
    myrank = get_myrank()
    if ( ob_definedRect( rect(myrank) ) ) then
       if ( .not. allocatedU(uc) ) print *,'*** error uc'
       u => uc%arr
       u(:,:,:,:) = 0.d0
       boolmap(:,:,:) = .false.
    endif
    call pkar_reset
    do rankd = 0, 400 -1
       if ( .not. ob_definedRect( rect(rankd) ) ) cycle
       call ob_extractLevelFromRect(level, rect(rankd) ) 
       call ob_extractPointFromRect(pL, rect(rankd), 'L')
       call ob_extractPointFromRect(pR, rect(rankd), 'R')
       do kob = pL%p(2), pR%p(2)+ksz, ksz
          klast = ( kob/ksz == pR%p(2)/ksz )
          if ( kob/ksz > pR%p(2)/ksz ) exit
          do job = pL%p(1), pR%p(1)+jsz, jsz
             jlast = ( job/jsz == pR%p(1)/jsz )
             if ( job/jsz > pR%p(1)/jsz ) exit
             do iob = pL%p(0), pR%p(0)+isz, isz
                ilast = ( iob/isz == pR%p(0)/isz )
                if ( iob/isz > pR%p(0)/isz ) exit
                ijkob = (/iob, job, kob/)
                call ob_assignCoordToPoint(point, ijkob, level)
                call ob_getBlockFromPoint(point, gid, ranks, i,j,k)
                if ( gid == Undefi ) cycle
                is = Imin
 js = Jmin
 ks = Kmin
                ie = Imax
 je = Jmax
 ke = Kmax
                if ( iob == pL%p(0) ) is = i
                if ( job == pL%p(1) ) js = j
                if ( kob == pL%p(2) ) ks = k
                if ( ilast ) ie = mod(pR%p(0), isz)
                if ( jlast ) je = mod(pR%p(1), jsz)
                if ( klast ) ke = mod(pR%p(2), ksz)
                if ( ranks == myrank ) then
                   up => get_Up(gid)
                   allocate(buf(is:ie,js:je,ks:ke,Mmin:Mmax))
                   buf = up(is:ie,js:je,ks:ke,Mmin:Mmax)
                   call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), kind(buf), rankd)
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
    do rankd = 0, 400 -1
       if ( rankd /= myrank ) cycle
       if ( .not. ob_definedRect( rect(rankd) ) ) cycle
       call ob_extractLevelFromRect(level, rect(rankd) ) 
       call ob_extractPointFromRect(pL, rect(rankd), 'L')
       call ob_extractPointFromRect(pR, rect(rankd), 'R')
       do kob = pL%p(2), pR%p(2)+ksz, ksz
          klast = ( kob/ksz == pR%p(2)/ksz )
          if ( kob/ksz > pR%p(2)/ksz ) exit
          do job = pL%p(1), pR%p(1)+jsz, jsz
             jlast = ( job/jsz == pR%p(1)/jsz )
             if ( job/jsz > pR%p(1)/jsz ) exit
             do iob = pL%p(0), pR%p(0)+isz, isz
                ilast = ( iob/isz == pR%p(0)/isz )
                if ( iob/isz > pR%p(0)/isz ) exit
                ijkob = (/iob, job, kob/)
                call ob_assignCoordToPoint(point, ijkob, level)
                call ob_getBlockFromPoint(point, gid, ranks, i,j,k)
                if ( gid == Undefi ) cycle
                is = Imin
 js = Jmin
 ks = Kmin
                ie = Imax
 je = Jmax
 ke = Kmax
                if ( iob == pL%p(0) ) is = i
                if ( job == pL%p(1) ) js = j
                if ( kob == pL%p(2) ) ks = k
                if ( ilast ) ie = mod(pR%p(0), isz)
                if ( jlast ) je = mod(pR%p(1), jsz)
                if ( klast ) ke = mod(pR%p(2), ksz)
                ius = -mod(pL%p(0), isz) + iob - uc%offsets(0) + is
                iue = ius + ie - is
                jus = -mod(pL%p(1), jsz) + job - uc%offsets(1) + js
                jue = jus + je - js
                kus = -mod(pL%p(2), ksz) + kob - uc%offsets(2) + ks
                kue = kus + ke - ks
                allocate( buf(is:ie,js:je,ks:ke,Mmin:Mmax) )
                if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), kin&
&d(buf), ranks)
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
  subroutine allocateUfromRect( u, rect )
    type(t_ugArray),pointer :: u
    type(t_obRect),intent(IN),dimension(0:400 -1) :: rect
    type(t_obPoint) :: pL, pR
    integer :: nrank, rank, imax, jmax, kmax
    if ( associated(u) ) then
       print *, '*** error: pointer u is already associated.', u%allocated, associated(u%arr)
       write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

    endif
    nrank = 0
    do rank = 0, 400 -1
       if ( ob_definedRect(rect(rank)) ) nrank = nrank + 1
    enddo
    allocate( u )
    u%nrank = nrank
    myrank = get_myrank()
    if ( ob_definedRect(rect(myrank)) ) then
       call ob_extractPointFromRect(pL, rect(myrank), 'L')
       call ob_extractPointFromRect(pR, rect(myrank), 'R')
       imax = pR%p(0)-pL%p(0)
       jmax = pR%p(1)-pL%p(1)
       kmax = pR%p(2)-pL%p(2)
       allocate( u%arr( 0:imax, 0:jmax, 0:kmax, Mmin:Mmax) )
       u%arr = 1 
       u%offsets(:) = pL%p(:)
       u%allocated = .true.
    else
       u%allocated = .false.
       nullify( u%arr )
    endif
  end subroutine allocateUfromRect
  subroutine deallocateU( u )
    type(t_ugArray),pointer :: u
    if (associated(u%arr)) deallocate( u%arr )
    nullify( u%arr )
    if (associated(u)) deallocate( u )
    nullify( u )
  end subroutine deallocateU
  function allocatedU( u ) result(bool)
    type(t_ugArray),pointer :: u
    logical :: bool
    bool = u%allocated
  end function allocatedU
  subroutine ug_write(x,y,z,time,step)
    use mpilib
    use string
    use io_util, only : wchar, read_env
    real(kind=8),dimension(:),intent(IN) :: x, y, z
    real(kind=8),intent(IN) :: time
    integer(kind=8),intent(IN) :: step
    integer,parameter :: FH = 11, FHT = 12
    character(len=CHARLEN) :: fn, dir, suffix
    real,dimension(:,:),allocatable :: u4
    real(kind=8),dimension(:,:),allocatable :: u
    real,dimension(:),allocatable :: x4, y4, z4
    integer :: header, lenint, m, k
    if( get_myrank() /= 0 ) return
    lenint = kind(lenint)
    allocate(u(size(x),size(y)), u4(size(x),size(y)), x4(size(x)), y4(size(y)), z4(size(z)))
    x4 = x
    y4 = y
    z4 = z
    call read_env('DIR', dir)
    call read_env('SUFFIX', suffix)
    fn = get_filename(dir, ThisPrefix, suffix)
    call wchar(6,'uniform grid data = '//fn)
    open(FH, file=fn, form='unformatted', access='stream') 
    header = lenint*4 + kind(step) + kind(time)
    write(FH) header
    write(FH) size(x4), size(y4), size(z4), Mmax-Mmin+1, step, time 
    write(FH) header
    header = size(x4)*kind(x4)
    write(FH) header
    write(FH) x4
    write(FH) header
    header = size(y4)*kind(y4)
    write(FH) header
    write(FH) y4
    write(FH) header
    header = size(z4)*kind(z4)
    write(FH) header
    write(FH) z4
    write(FH) header
    header = size(x4)*size(y4)*size(z4)*(Mmax-Mmin+1)*kind(u4)
    write(FH) header
    do m = Mmin, Mmax
       open(unit=FHT, file=get_tempfilename(m), form='unformatted', access='stream')
       do k = 1, size(z4)
          read(FHT) u
          u4 = u
          write(FH) u4
       end do
       call flush(FHT)
       close(FHT)
    end do
    write(FH) header
    call flush(FH)
    close(FH)
    deallocate( u, u4, x4, y4, z4 )
  end subroutine ug_write
  subroutine ug_fld(x,y,z,time,step)
    use mpilib
    use io_util, only : wchar, read_env
    use string
    real(kind=8),dimension(:),intent(IN) :: x, y, z
    real(kind=8),intent(IN) :: time
    integer(kind=8),intent(IN) :: step
    character(len=CHARLEN) :: fn, dir, suffix
    character(len=3),parameter :: Suffixfld = 'fld'
    integer,parameter :: FH = 11
    integer :: ni, nj, nk, nm, ndim, nspace, idepthv, idepthc, hbyte, fbyte, m
    integer(kind=8) :: ivol, iskip
    if( get_myrank() /= 0 ) return
    ni = size(x,1)
    nj = size(y,1)
    nk = size(z,1)
    nm = Mmax-Mmin+1
    ndim = 3
    nspace = 3
    call read_env('DIR', dir)
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
    idepthv = 4 
    idepthc = 4 
    hbyte = 4 
    fbyte = 4 
    ivol=int(ni, 8)*int(nj, 8)*int(nk, 8)*int(idepthv, 8)
    iskip = 0
    iskip = iskip + hbyte 
    iskip = iskip + 4*4 
    iskip = iskip + 8 
    iskip = iskip + 8 
    iskip = iskip + fbyte 
    call read_env('SUFFIX', suffix)
    fn = get_filename('', ThisPrefix, suffix)
    iskip = iskip + hbyte
    call wchar(FH,concat( &
         concat('coord 1 file=',fn), &
         concat(' filetype=binary skip=',strcmprs(num2char(iskip)))))
    iskip = iskip + ni*idepthc+fbyte
    iskip = iskip + hbyte
    call wchar(FH,concat( &
         concat('coord 2 file=',fn), &
         concat(' filetype=binary skip=',strcmprs(num2char(iskip)))))
    iskip = iskip + nj*idepthc+fbyte
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
    iskip = iskip + fbyte
    call wchar(FH,concat( &
         concat('# size of file, ',fn), &
         concat(' = ',strcmprs(num2char(iskip)))))
    call flush(FH)
    close(FH)
  end subroutine ug_fld
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
  function get_tempfilename2(rank) result(fn)
    integer,intent(IN) :: rank
    character(len=CHARLEN) :: fn
    fn = get_tempfilename(rank)
    fn = trim(fn)//trim('d')
  end function get_tempfilename2
  function get_filename(dir, prefix, suffix) result(fn)
    use string, only : concat, CHARLEN, num2char
    use grid, only : Step, LevelMax
    character(len=*),intent(IN) :: dir, prefix, suffix
    character(len=CHARLEN) :: fn
    character(len=1) :: dot = '.'
    fn = dir 
    fn = concat(fn,prefix) 
    fn = concat(fn,num2char(Step(LevelMax))) 
    fn = concat(fn,dot) 
    fn = concat(fn,num2char(LevelResolution)) 
    fn = concat(fn,dot) 
    fn = concat(fn,suffix) 
  end function get_filename
  subroutine ug_interp(uf, uc, rectf, rect, interpolate)
    type(t_ugArray),pointer :: uf, uc
    type(t_obRect),dimension(0:400 -1),intent(IN) :: rectf, rect
    logical,intent(IN),optional :: interpolate
    type(t_ugArray),save,pointer :: ut => null()
    type(t_obRect),dimension(0:400 -1) :: rectt
    integer :: level, rank
    do rank = 0, 400 -1
       if ( ob_definedRect( rect(rank) ) ) then
          call ob_extractLevelFromRect(level, rect(rank))
          call ob_rect2RectByLevel(rect(rank), rectt(rank), level+1)
       else
          call ob_undefRect(rectt(rank))
       endif
    enddo
    call overlapFixU( uc )
    call allocateUfromRect( ut, rectt )
    call interpLocal( ut, uc, interpolate )
    call transferU( ut, uf, rectt, rectf ) 
    call deallocateU( ut )
  end subroutine ug_interp
  subroutine interpLocal(ufa, uca, interpolate)
    type(t_ugArray),pointer :: ufa, uca
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc
    logical,intent(IN),optional :: interpolate
    real(kind=8),dimension(:,:,:,:),allocatable :: grad
    real(kind=8) :: di, dj, dk, hf, hc, hfi, hci, dvc, dvf
    integer :: m, if, jf, kf, ic, jc, kc
    integer :: ics, ice, jcs, jce, kcs, kce, ifs, ife, jfs, jfe, kfs, kfe
    logical :: bool_interpolate
    if (.not. allocatedU( uca ) ) return
    bool_interpolate = .false.
    if ( present( interpolate) ) then
       bool_interpolate = interpolate
    endif
    uf => ufa%arr
    uc => uca%arr
    hf = 1.d0 
    hc = 2* hf 
    dvc = 1.d0 
    dvf = 1.d0 
    hfi = 1.d0/hf
    hci = 1.d0/hc
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
    allocate( grad(lbound(uc,1):ubound(uc,1),lbound(uc,2):ubound(uc,2),lbound(uc,3):ubound(uc,3),0:2) )
    do m = Mmin, Mmax
       grad(:,:,:,:) = 0.d0
       if ( bool_interpolate ) then
          do ic = ics+1, ice-1
             grad(ic,:,:,0) = ( sign(1.d0,(uc(ic+1,:,:,m)-uc(ic,:,:,m)))*max(0.d0,min(abs(uc(ic+1,:,:,m)-uc(ic,:,:,m)),sign(1.d0,(u&
&c(ic+1,:,:,m)-uc(ic,:,:,m)))*(uc(ic,:,:,m)-uc(ic-1,:,:,m)))) )*hci
          enddo
          do jc = jcs+1, jce-1
             grad(:,jc,:,1) = ( sign(1.d0,(uc(:,jc+1,:,m)-uc(:,jc,:,m)))*max(0.d0,min(abs(uc(:,jc+1,:,m)-uc(:,jc,:,m)),sign(1.d0,(u&
&c(:,jc+1,:,m)-uc(:,jc,:,m)))*(uc(:,jc,:,m)-uc(:,jc-1,:,m)))) )*hci
          enddo
          do kc = kcs+1, kce-1
             grad(:,:,kc,2) = ( sign(1.d0,(uc(:,:,kc+1,m)-uc(:,:,kc,m)))*max(0.d0,min(abs(uc(:,:,kc+1,m)-uc(:,:,kc,m)),sign(1.d0,(u&
&c(:,:,kc+1,m)-uc(:,:,kc,m)))*(uc(:,:,kc,m)-uc(:,:,kc-1,m)))) )*hci
          enddo
          if ( ics /= ice ) then
             grad(ics,:,:,0) = ( uc(ics+1,:,:,m)-uc(ics,:,:,m) )*hci
             grad(ice,:,:,0) = ( uc(ice,:,:,m)-uc(ice-1,:,:,m) )*hci
          endif
          if ( jcs /= jce ) then
             grad(:,jcs,:,1) = ( uc(:,jcs+1,:,m)-uc(:,jcs,:,m) )*hci
             grad(:,jce,:,1) = ( uc(:,jce,:,m)-uc(:,jce-1,:,m) )*hci
          endif
          if ( kcs /= kce ) then
             grad(:,:,kcs,2) = ( uc(:,:,kcs+1,m)-uc(:,:,kcs,m) )*hci
             grad(:,:,kce,2) = ( uc(:,:,kce,m)-uc(:,:,kce-1,m) )*hci
          endif
       endif
       do kf = kfs, kfe
          do jf = jfs, jfe
             do if = ifs, ife
                ic = ((if+ufa%offsets(0))-(0))/2 + mod(min((if+ufa%offsets(0))-(0),0),2) + (0)-uca%offsets(0)
                jc = ((jf+ufa%offsets(1))-(0))/2 + mod(min((jf+ufa%offsets(1))-(0),0),2) + (0)-uca%offsets(1)
                kc = ((kf+ufa%offsets(2))-(0))/2 + mod(min((kf+ufa%offsets(2))-(0),0),2) + (0)-uca%offsets(2)
                di = ( modulo(if+ufa%offsets(0),int(2,8)) - 0.5d0 )*hf
                dj = ( modulo(jf+ufa%offsets(1),int(2,8)) - 0.5d0 )*hf
                dk = ( modulo(kf+ufa%offsets(2),int(2,8)) - 0.5d0 )*hf
                uf(if,jf,kf,m) = uc(ic,jc,kc,m) + grad(ic,jc,kc,0)*di+grad(ic,jc,kc,1)*dj+grad(ic,jc,kc,2)*dk
             enddo
          enddo
       enddo
    enddo
    deallocate( grad )
    call conv_u2w_smallMem(uf, dvf, W2U)
  end subroutine interpLocal
  subroutine transferU( us, ud, rects, rectd )
    type(t_ugArray),pointer :: us, ud
    type(t_obRect),dimension(0:400 -1),intent(IN) :: rects, rectd
    real(kind=8),dimension(:,:,:,:),allocatable :: buf
    integer :: ranks, rankd, test
    type(t_obRect) :: rectand, recttest
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX) :: coord
    integer :: iss, ies, jss, jes, kss, kes 
    integer :: isd, ied, jsd, jed, ksd, ked 
    myrank = get_myrank()
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
    do ranks = 0, 400 -1
       if ( .not. ob_definedRect( rects(ranks) ) ) cycle
       do rankd = 0, 400 -1
          if ( .not. ob_definedRect( rectd(rankd) ) ) cycle
          call ob_rectAnd( rects(ranks), rectd(rankd), rectand )
          if ( ob_testRect(rectand) /= OB_RECT_VALID ) cycle
          call ob_extractCoordFromRect(coord, rectand)
          iss = coord(0)-us%offsets(0)
 jss = coord(1)-us%offsets(1)
 kss = coord(2)-us%offsets(2)

          ies = coord(3)-us%offsets(0)
 jes = coord(4)-us%offsets(1)
 kes = coord(5)-us%offsets(2)

          isd = coord(0)-ud%offsets(0)
 jsd = coord(1)-ud%offsets(1)
 ksd = coord(2)-ud%offsets(2)

          ied = coord(3)-ud%offsets(0)
 jed = coord(4)-ud%offsets(1)
 ked = coord(5)-ud%offsets(2)

          if ( myrank == ranks .and. myrank == rankd ) then
             ud%arr(isd:ied,jsd:jed,ksd:ked,Mmin:Mmax) = us%arr(iss:ies,jss:jes,kss:kes,Mmin:Mmax)
          else if ( myrank == ranks ) then
             allocate( buf(iss:ies,jss:jes,kss:kes,Mmin:Mmax) )
             buf = us%arr(iss:ies,jss:jes,kss:kes,Mmin:Mmax)
             call mpi_send(buf, size(buf), MPI_DOUBLE_PRECISION, rankd, 1, MPI_COMM_WORLD, ierr )
             deallocate( buf )
          else if ( myrank == rankd ) then
             allocate( buf(isd:ied,jsd:jed,ksd:ked,Mmin:Mmax) )
             call mpi_recv(buf, size(buf), MPI_DOUBLE_PRECISION, ranks, 1, MPI_COMM_WORLD, status, ierr )
             ud%arr(isd:ied,jsd:jed,ksd:ked,Mmin:Mmax) = buf
             deallocate( buf )
          end if
       enddo
    enddo
  end subroutine transferU
  subroutine u2rect( u, level, rect )
    type(t_ugArray),intent(IN) :: u
    integer,intent(IN) :: level
    type(t_obRect),intent(OUT) :: rect
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX) :: coords
    coords(0:2) = u%offsets(0:2)
    coords(2 +1:OB_COORDS_MAX) = u%offsets(0:2) + (/ubound(u%arr,1),ubound(u%arr,2),ubound(u%arr,3)/)
    call ob_assignCoordToRect( coords, rect, level )
  end subroutine u2rect
  subroutine overlapFixU( ua )
    type(t_ugArray),pointer :: ua
    real(kind=8),dimension(:,:,:,:),allocatable :: bufs, bufd
    real(kind=8),dimension(:,:,:,:),pointer :: u
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
    if ( myrank == nrank - 1 ) rankR = MPI_PROC_NULL
    if ( myrank == 0 ) rankL = MPI_PROC_NULL
    bufs = u(:,:,ubound(u,3)-1:ubound(u,3)-1,:)
    call mpi_sendrecv( &
         bufs, size(bufs), MPI_DOUBLE_PRECISION, rankR, TAG, &
         bufd, size(bufd), MPI_DOUBLE_PRECISION, rankL, TAG, MPI_COMM_WORLD, status, ierr )
    if (myrank /= 0) u(:,:,lbound(u,3):lbound(u,3),:) = bufd
    bufs = u(:,:,lbound(u,3)+1:lbound(u,3)+1,:)
    call mpi_sendrecv( &
         bufs, size(bufs), MPI_DOUBLE_PRECISION, rankL, TAG, &
         bufd, size(bufd), MPI_DOUBLE_PRECISION, rankR, TAG, MPI_COMM_WORLD, status, ierr )
    if (myrank /= nrank-1) u(:,:,ubound(u,3):ubound(u,3),:) = bufd
    deallocate( bufs, bufd )
  end subroutine overlapFixU
  subroutine conv_u2w_smallMem(u, dv, code)
    use eos, only: conv_u2w, conv_w2u
    real(kind=8),dimension(:,:,:,:),intent(INOUT) :: u
    real(kind=8),intent(IN) :: dv
    integer,intent(IN) :: code
    integer,parameter :: szi = 8, szj = 8, szk = 8 
    real(kind=8),dimension(0:szi-1,0:szj-1,0:szk-1, Mmin:Mmax) :: buffer
    integer :: is, js, ks, ie, je, ke, ii, jj, kk, i, j, k
    is = lbound(u,1) 
 js = lbound(u,2) 
 ks = lbound(u,3)
    ie = ubound(u,1) 
 je = ubound(u,2) 
 ke = ubound(u,3)
    do k = ks, ke, szk
       do j = js, je, szj
          do i = is, ie, szi
             buffer = 1.d0 
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
  subroutine ug_DownSize(dlev, rectg)
    use string, only : CHARLEN, concat
    integer,intent(IN) :: dlev
    type(t_obRect),intent(INOUT) :: rectg
    type(t_obRect) :: rectc
    type(t_obPoint) :: pL, pR
    character(len=CHARLEN) :: fn, fn2
    integer(kind=8) :: ics, ice, jcs, jce, kcs, kce, iob, job, kob
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX) :: coords
    integer :: szi, szj, szk, level, if, jf, kf, ic, jc, kc, m, lev
    integer,parameter :: FH=11, FH2 = 12
    real(kind=8),dimension(:,:,:),allocatable :: uf
    real(kind=8),dimension(:,:),allocatable :: uc, dummy
    logical :: exist
    if (dlev <= 0) return
    myrank = get_myrank()
    if ( myrank /= 0 ) return
    call conv_u2w_tempFiles(U2W, rectg)
    call clearTempFiles2
    do level = 1, dlev
       call ob_extractPointFromRect(pL, rectg, 'L')
       call ob_extractPointFromRect(pR, rectg, 'R')
       call ob_extractLevelFromRect(lev, rectg)
       szi = pR%p(0) - pL%p(0) + 1
       szj = pR%p(1) - pL%p(1) + 1
       szk = pR%p(2) - pL%p(2) + 1
       ics = pL%p(0)/2 + mod(pL%p(0),2)
       jcs = pL%p(1)/2 + mod(pL%p(1),2)
       kcs = pL%p(2)/2 + mod(pL%p(2),2)
       ice = pR%p(0)/2 + mod(pR%p(0),2)-1
       jce = pR%p(1)/2 + mod(pR%p(1),2)-1
       kce = pR%p(2)/2 + mod(pR%p(2),2)-1
       coords = (/ ics, jcs, kcs, ice, jce, kce /)
       call ob_assignCoordToRect( coords, rectc, lev-1 )
       allocate(uf(0:szi-1, 0:szj-1, 0:1))
       allocate(uc(0:ice-ics, 0:jce-jcs))
       do m = Mmin, Mmax
          fn = get_tempfilename(m)
          fn2 = get_tempfilename2(m)
          inquire(file=fn, exist=exist)
          if (.not. exist) then
             print *, '** Error in ug_DownSize. fn not exist', m, level
             call flush(6)
             stop
          endif
          open(unit=FH, file=fn, form='unformatted', access='stream')
          open(unit=FH2, file=fn2, form='unformatted', access='stream')
          if ( mod( pL%p(2), 2 ) /= 0 ) then
             allocate(dummy(lbound(uf,1):ubound(uf,1),lbound(uf,2):ubound(uf,2)))
             read(unit=FH) dummy
             deallocate(dummy)
          end if
          do kob = kcs, kce
             read(FH) uf
             do job = jcs, jce
                do iob = ics, ice
                   ic = iob - ics
                   jc = job - jcs
                   kc = kob - kcs
                   if = iob * 2 - pL%p(0)
                   jf = job * 2 - pL%p(1)
                   kf = kob * 2 - pL%p(2)
                   uc(ic,jc) = &
                        (uf(if, jf,0) + uf(if+1,jf+1,0) + &
                        uf(if+1,jf,0) + uf(if, jf+1,0) + &
                        uf(if, jf,1) + uf(if+1,jf+1,1) + &
                        uf(if+1,jf,1) + uf(if, jf+1,1) )/8.d0
                end do
             end do
             write(FH2) uc
          end do
          call flush(FH)
          call flush(FH2)
          close(unit=FH)
          close(unit=FH2)
       end do
       deallocate(uf, uc)
       rectg = rectc
       call swapTempFiles
    end do
    call conv_u2w_tempFiles(W2U, rectg)
  end subroutine ug_DownSize
  subroutine refile(u)
    use string, only : CHARLEN
    real(KIND=8),dimension(:,:,:,:),pointer :: u 
    integer :: m, rank
    integer,parameter :: FH = 11
    myrank = get_myrank()
    call clearTempFiles
    do rank = 0, 400 -1 
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       if (rank /= myrank) cycle
       if (.not. associated(u)) cycle
       do m = Mmin, Mmax
          open(unit=FH, file=get_tempfilename(m), form='unformatted', access='stream', position='append')
          write(FH) u(:,:,:,m)
          call flush(FH)
          close(FH)
       end do
    end do
  end subroutine refile
  subroutine clearTempFiles
    use systemcall
    integer :: m
    character(len=CHARLEN) :: fn
    logical :: exist
    myrank = get_myrank()
    if ( myrank /= 0 ) return
    do m = Mmin, Mmax
       fn = get_tempfilename(m)
       inquire(file=fn, exist=exist)
       if (exist) call systemcall_unlink(fn)
    end do
  end subroutine clearTempFiles
  subroutine clearTempFiles2
    use systemcall
    integer :: m
    character(len=CHARLEN) :: fn
    logical :: exist
    myrank = get_myrank()
    if ( myrank /= 0 ) return
    do m = Mmin, Mmax
       fn = get_tempfilename2(m)
       inquire(file=fn, exist=exist)
       if (exist) call systemcall_unlink(fn)
    end do
  end subroutine clearTempFiles2
  subroutine makeCoordinates(x, y, z, rectg)
    real(kind=8),dimension(:),pointer :: x, y, z 
    type(t_obRect),intent(IN) :: rectg
    type(t_obPoint) :: pL, pR
    type(t_obRectPhys) :: rectCompBox
    type(t_obPointPhys) :: pcL
    integer :: imax, jmax, kmax, level, i, j, k
    real(kind=8),dimension(0:2) :: h
    if( get_myrank() /= 0 ) return
    call ob_extractPointFromRect(pL, rectg, 'L')
    call ob_extractPointFromRect(pR, rectg, 'R')
    imax = pR%p(0)-pL%p(0)
    jmax = pR%p(1)-pL%p(1)
    kmax = pR%p(2)-pL%p(2)
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
       x(i) = (i + pL%p(0) + 0.5d0) * h(0) + pcL%p(0)
    end do
    do j = 0, jmax
       y(j) = (j + pL%p(1) + 0.5d0) * h(1) + pcL%p(1)
    end do
    do k = 0, kmax
       z(k) = (k + pL%p(2) + 0.5d0) * h(2) + pcL%p(2)
    end do
  end subroutine makeCoordinates
  subroutine rectGather(rect, rectg)
    type(t_obRect),dimension(0:400 -1),intent(IN) :: rect
    type(t_obRect),intent(OUT) :: rectg
    type(t_obRect) :: rectswap
    integer :: rank
    rectg = rect(0)
    do rank = 1, (400) - 1
       if (.not. ob_definedRect(rect(rank))) cycle
       call ob_rectOr(rect(rank), rectg, rectswap)
       rectg = rectswap
    end do
  end subroutine rectGather
  subroutine conv_u2w_tempFiles(code, rect)
    integer,intent(IN) :: code
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: pL, pR
    real(kind=8),dimension(:,:),allocatable :: buf
    real(kind=8),dimension(:,:,:,:),allocatable :: u
    real(kind=8),parameter :: dv = 1.d0
    integer :: m
    integer(kind=8) :: kob
    integer,parameter :: FH=11, FH2 = FH + Mmax-Mmin + 1
    logical :: exist
    call clearTempFiles2
    do m = Mmin, Mmax
       inquire(file=get_tempfilename(m), exist=exist)
       if (.not. exist) then
          print *, 'conv_u2w_tempFiles: FN not exist', m
          call flush(6)
       endif
       open(unit=FH+m, file=get_tempfilename(m), form='unformatted', access='stream')
       open(unit=FH2+m, file=get_tempfilename2(m), form='unformatted', access='stream')
    end do
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    allocate( buf(0:pR%p(0)-pL%p(0), 0:pR%p(1)-pL%p(1)) )
    allocate( u(0:pR%p(0)-pL%p(0), 0:pR%p(1)-pL%p(1), 1, Mmin:Mmax) )
    do kob = pL%p(2), pR%p(2)
       do m = Mmin, Mmax
          read(unit=FH+m) buf
          u(:,:,1,m) = buf
       end do
       call conv_u2w_smallMem(u, dv, code)
       do m = Mmin, Mmax
          buf = u(:,:,1,m)
          write(unit=FH2+m) buf
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
  subroutine swapTempFiles
    use systemcall
    use string, only : CHARLEN
    integer :: m
    character(len=CHARLEN) :: fn, fn2
    logical :: exist, exist2
    do m = Mmin, Mmax
       fn = get_tempfilename(m)
       fn2 = get_tempfilename2(m)
       inquire(file=fn, exist=exist)
       inquire(file=fn2, exist=exist2)
       if (.not. exist) print *, 'swapTempFiles: FN not exist', m
       if (.not. exist2) print *, 'swapTempFiles: FN2 not exist', m
       if ( exist .and. exist2 ) then
          call systemcall_unlink(fn)
       end if
       if ( exist2 ) then
          call systemcall_rename(fn2, fn)
       end if
       inquire(file=fn, exist=exist)
       if (.not. exist) then
          print *, '** Error in end of swapTempFiles: FN not exist', m
       endif
    end do
    call flush(6)
  end subroutine swapTempFiles
  subroutine printrect(rect)
    type(t_obRect),dimension(0:400 -1),intent(IN) :: rect
    integer :: rank
    integer(kind=8) :: is, ie, js, je, ks, ke
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX) :: coords
    is = HUGE(is)
    js = HUGE(js)
    ks = HUGE(ks)
    ie = -HUGE(ie)
    je = -HUGE(je)
    ke = -HUGE(ke)
    do rank = 0, 400 -1
       if ( .not. ob_definedRect( rect(rank) )) cycle
       call ob_extractCoordFromRect(coords, rect(rank))
       is = min(is,coords(0))
       js = min(js,coords(1))
       ks = min(ks,coords(2))
       ie = max(ie,coords(3))
       je = max(je,coords(4))
       ke = max(ke,coords(5))
    enddo
   if (get_myrank() == 0) print *, 'rect', is,ie, js,je, ks,ke
  end subroutine printrect
  subroutine printUwidth( u )
    type(t_ugArray),pointer :: u
    integer(kind=8) :: is, ie, js, je, ks, ke
    integer :: buf(0:2)
    integer :: nrank, szk, szkr
    myrank = get_myrank()
    if (myrank == 0) then
       nrank = u%nrank
    end if
    call mpi_bcast(nrank, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if ( myrank == nrank-1) then
       buf(0) = u%offsets(0) + ubound(u%arr,1)
       buf(1) = u%offsets(1) + ubound(u%arr,2)
       buf(2) = u%offsets(2) + ubound(u%arr,3)
       call mpi_send(buf(0), size(buf), MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr )
    endif
    if (myrank == 0) then
       call mpi_recv(buf(0), size(buf), MPI_INTEGER, nrank-1, 1, MPI_COMM_WORLD, status, ierr )
       ie = buf(0)
       je = buf(1)
       ke = buf(2)
       is = u%offsets(0)
       js = u%offsets(1)
       ks = u%offsets(2)
       print *, 'u size ', is,ie,js,je,ks,ke
    endif
    if ( allocatedU(u) ) then
       szk = size(u%arr,3)
    else
       szk = 0
    endif
    call mpi_allreduce(szk, szkr, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
    if (myrank == 0) &
         print *, 'u size*', size(u%arr,1),size(u%arr,2), szkr-2*nrank+2
  end subroutine printUwidth
  subroutine checkbound(u, is, ie, js,je,ks,ke )
    real(kind=8),dimension(:,:,:,:),pointer :: u
    integer,intent(IN) :: is, ie, js,je,ks,ke
    if (&
         is < lbound(u,1) .or. js < lbound(u,2) .or. ks < lbound(u,3) .or. &
         ie > ubound(u,1) .or. je > ubound(u,2) .or. ke > ubound(u,3) ) then
       print *, '**** error'
    endif
  end subroutine checkbound
end module uniformgrid
