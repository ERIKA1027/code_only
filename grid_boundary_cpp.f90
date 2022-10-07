module grid_boundary
  use grid , only : Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh, Mmin, Mmax, Lmin, Lmax, Ngh, Gidmin, Gidmax, Left, Right, Imin,&
& Imax, Jmin, Jmax, Kmin, Kmax, COMPLETE, PREDICTOR, CORRECTOR, &
       U_StepNumber, U_StepNumberGhostCell, U1_StepNumber, U1_StepNumberGhostCell, U2_StepNumber, U2_StepNumberGhostCell
  implicit none
  private
  integer,parameter :: Ev = 0, Od = 1 
  integer,parameter :: Send = 0, Recv = 1 
  integer,save :: STEP_MODE 
  integer,save :: CurrentLevel 
  integer,save :: LevelUpto 
  logical,save,dimension(Lmin:Lmax) :: Bool_fixed
  logical,save :: Bool_fix_current, Bool_fix_1order, Bool_fix_2order
  integer,save,dimension(Left:Right,0:2,Send:Recv) :: Ins, Ine, Jns, Jne, Kns, Kne
  integer,save,dimension(Left:Right,0:2,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(:),allocatable :: FLS
  type(t_buf),save,dimension(:),allocatable :: FLP
  public :: boundary_grid
contains
  subroutine boundary_grid( level, mode )
    integer,intent(IN) :: level, mode 
    integer :: ndir
    if ( level < Lmin ) return
    STEP_MODE = mode
    LevelUpto = level
    call grid_boundary_init
    call gridsize_init_parent
    call gridsize_init_samelev
    Bool_fix_current = .true.
    Bool_fix_1order = .true.
    Bool_fix_2order = .true.
    CurrentLevel = LevelUpto-1
    do ndir = 0, 2
       call fix_from_samelev(ndir)
    enddo
    Bool_fix_current = .true.
    Bool_fix_1order = .false.
    Bool_fix_2order = .false.
    CurrentLevel = LevelUpto
    do ndir = 0, 2
       call fix_from_parent(ndir)
       call fix_from_samelev(ndir)
    enddo
    if (Bool_fix_current) U_StepNumberGhostCell(level) = U_StepNumber(level)
    if (Bool_fix_1order) U1_StepNumberGhostCell(level) = U1_StepNumber(level)
    if (Bool_fix_2order) U2_StepNumberGhostCell(level) = U2_StepNumber(level)
  end subroutine boundary_grid
  subroutine grid_boundary_init
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary' )
    Bool_fixed(:) = .false.
    bool_initialized = .true.
  end subroutine grid_boundary_init
  subroutine gridsize_init_samelev
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    integer :: lr, ndir
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary: samelev' )
    Ins(:,:,:) = Imingh
    Ine(:,:,:) = Imaxgh
    Jns(:,:,:) = Jmingh
    Jne(:,:,:) = Jmaxgh
    Kns(:,:,:) = Kmingh
    Kne(:,:,:) = Kmaxgh
    Ins(Left ,0,Send) = Imax-Ngh+1 
    Ine(Left ,0,Send) = Imax
    Ins(Right,0,Send) = Imin
    Ine(Right,0,Send) = Imin+Ngh-1
    Jns(Left ,1,Send) = Jmax-Ngh+1
    Jne(Left ,1,Send) = Jmax
    Jns(Right,1,Send) = Jmin
    Jne(Right,1,Send) = Jmin+Ngh-1
    Kns(Left ,2,Send) = Kmax-Ngh+1
    Kne(Left ,2,Send) = Kmax
    Kns(Right,2,Send) = Kmin
    Kne(Right,2,Send) = Kmin+Ngh-1
    Ins(Left ,0,Recv) = Imin-Ngh 
    Ine(Left ,0,Recv) = Imin-1
    Ins(Right,0,Recv) = Imax+1
    Ine(Right,0,Recv) = Imax+Ngh
    Jns(Left ,1,Recv) = Jmin-Ngh
    Jne(Left ,1,Recv) = Jmin-1
    Jns(Right,1,Recv) = Jmax+1
    Jne(Right,1,Recv) = Jmax+Ngh
    Kns(Left ,2,Recv) = Kmin-Ngh
    Kne(Left ,2,Recv) = Kmin-1
    Kns(Right,2,Recv) = Kmax+1
    Kne(Right,2,Recv) = Kmax+Ngh
    if ( .not. allocated(FLS) ) allocate( FLS(0:2) )
    if ( .not. associated(FLS(0)%Bufc) ) then
       lr = Left
       do ndir = 0, 2
          allocate( FLS(ndir)%Bufc( Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(l&
&r,ndir,Send),Mmin:Mmax ) )
       end do
    endif
    bool_initialized = .true.
  end subroutine gridsize_init_samelev
  subroutine gridsize_init_parent
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    integer :: lr, ndir, ieo, jeo, keo
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary: parent' )
    Ics(:,:,Ev) = Imingh
    Ice(:,:,Ev) = (Imax+1)/2-1+Ngh
    Ics(:,:,Od) = (Imax+1)/2-Ngh
    Ice(:,:,Od) = Imaxgh
    Jcs(:,:,Ev) = Jmingh
    Jce(:,:,Ev) = (Jmax+1)/2-1+Ngh
    Jcs(:,:,Od) = (Jmax+1)/2-Ngh
    Jce(:,:,Od) = Jmaxgh
    Kcs(:,:,Ev) = Kmingh
    Kce(:,:,Ev) = (Kmax+1)/2-1+Ngh
    Kcs(:,:,Od) = (Kmax+1)/2-Ngh
    Kce(:,:,Od) = Kmaxgh
    Ics(Left ,0,:) = Imax-Ngh/2
    Ice(Left ,0,:) = Imax+Ngh/2
    Ics(Right,0,:) = Imin-Ngh/2
    Ice(Right,0,:) = Imin+Ngh/2
    Jcs(Left ,1,:) = Jmax-Ngh/2
    Jce(Left ,1,:) = Jmax+Ngh/2
    Jcs(Right,1,:) = Jmin-Ngh/2
    Jce(Right,1,:) = Jmin+Ngh/2
    Kcs(Left ,2,:) = Kmax-Ngh/2
    Kce(Left ,2,:) = Kmax+Ngh/2
    Kcs(Right,2,:) = Kmin-Ngh/2
    Kce(Right,2,:) = Kmin+Ngh/2
    if ( .not. allocated(FLP) ) allocate( FLP(0:2) )
    if ( .not. associated(FLP(0)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = 0, 2
          allocate( FLP(ndir)%Bufc( Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndi&
&r,keo),Mmin:Mmax ) )
       end do
    endif
    bool_initialized = .true.
  end subroutine gridsize_init_parent
  subroutine fix_from_samelev(ndir)
    use mpilib
    use packarr
    use grid
    integer,intent(IN) :: ndir
    integer :: n, lr, rank, rankd, ranks, gid, gidd, gids
    real(kind=8),pointer,dimension(:,:,:,:) :: us, ud, buf
    if ( CurrentLevel < Lmin ) return
    myrank = get_myrank()
    buf => FLS(ndir)%Bufc
    call pkar_reset
    do rank = 0, 400 -1 
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             gids = NeighborGid(lr, ndir, gidd, rankd)
             ranks = NeighborRank(lr, ndir, gidd, rankd)
             if ( gids == Undefi ) cycle 
             if ( myrank == rankd .and. myrank == ranks ) then 
                if (Bool_fix_current ) then
                   us => get_Up(gids)
                   ud => get_Up(gidd)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndir,Send),Mmin:Mm&
&ax)
                endif
                if (Bool_fix_1order) then
                   us => get_U1orderp(gids)
                   ud => get_U1orderp(gidd)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndir,Send),Mmin:Mm&
&ax)
                endif
                if (Bool_fix_2order) then
                   us => get_U2orderp(gids)
                   ud => get_U2orderp(gidd)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndir,Send),Mmin:Mm&
&ax)
                endif
             else if ( myrank == ranks .or. myrank == rankd ) then 
                if (Bool_fix_current) then
                   if (myrank == ranks) then
                      us => get_Up(gids)
                      buf = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndi&
&r,Send),Mmin:Mmax)
                   endif
                   if ((myrank) == (ranks)) call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf),&
& kind(buf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(buf), kind(buf), ranks)
                endif
                if (Bool_fix_1order) then
                   if (myrank == ranks) then
                      us => get_U1orderp(gids)
                      buf = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndi&
&r,Send),Mmin:Mmax)
                   endif
                   if ((myrank) == (ranks)) call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf),&
& kind(buf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(buf), kind(buf), ranks)
                endif
                if (Bool_fix_2order) then
                   if (myrank == ranks) then
                      us => get_U2orderp(gids)
                      buf = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndi&
&r,Send),Mmin:Mmax)
                   endif
                   if ((myrank) == (ranks)) call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf),&
& kind(buf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(buf), kind(buf), ranks)
                endif
             endif
          end do
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1 
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             gids = NeighborGid(lr, ndir, gidd, rankd)
             ranks = NeighborRank(lr, ndir, gidd, rankd)
             if ( gids == Undefi ) cycle 
             if ( myrank == ranks ) cycle
             if ( myrank == rankd ) then 
                if (Bool_fix_current) then
                   ud => get_Up(gidd)
                   if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), &
&kind(buf), ranks)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = buf
                endif
                if (Bool_fix_1order) then
                   ud => get_U1orderp(gidd)
                   if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), &
&kind(buf), ranks)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = buf
                endif
                if (Bool_fix_2order) then
                   ud => get_U2orderp(gidd)
                   if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), &
&kind(buf), ranks)
                   ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),M&
&min:Mmax) = buf
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine fix_from_samelev
  subroutine fix_from_parent(ndir)
    use mpilib
    use packarr
    use grid
    integer,intent(IN) :: ndir
    integer :: levelc, n, pgid, prank
    integer :: gid, gids, gidd, lr, rank, rankd, ranks
    real(kind=8),pointer,dimension(:,:,:,:) :: ud, udbuf
    if ( CurrentLevel <= Lmin ) return
    call pkar_reset
    udbuf => FLP(ndir)%Bufc
    levelc = CurrentLevel -1
    myrank = get_myrank()
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             if ( NeighborGid(lr, ndir, gidd, rankd) /= Undefi ) cycle 
             pgid = ParentGid(gidd,rankd)
             prank = ParentRank(gidd,rankd)
             gids = NeighborGid(lr, ndir, pgid, prank)
             ranks = NeighborRank(lr, ndir, pgid, prank)
             if ( gids == Undefi ) cycle 
             if (myrank == ranks) call interp_time(udbuf, lr, ndir, gids, gidd, rankd)
             if ((myrank) == (ranks)) call pkar_push(udbuf(lbound(udbuf,1),lbound(udbuf,2),lbound(udbuf,3),lbound(udbuf,4)), size(u&
&dbuf), kind(udbuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(udbuf), kind(udbuf), ranks)
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             if ( NeighborGid(lr, ndir, gidd, rankd) /= Undefi ) cycle 
             pgid = ParentGid(gidd,rankd)
             prank = ParentRank(gidd,rankd)
             gids = NeighborGid(lr, ndir, pgid, prank)
             ranks = NeighborRank(lr, ndir, pgid, prank)
             if ( gids == Undefi ) cycle 
             if (myrank == rankd) then
                if ((myrank) == (rankd)) call pkar_pop(udbuf(lbound(udbuf,1),lbound(udbuf,2),lbound(udbuf,3),lbound(udbuf,4)), size&
&(udbuf), kind(udbuf), ranks)
                ud => get_Up(gidd)
                call interp_space(ud, udbuf, ndir, lr)
             end if
          enddo
       enddo
    enddo
  end subroutine fix_from_parent
  subroutine interp_space( ud, buf, ndir, lr)
    use eos
    real(kind=8),pointer,dimension(:,:,:,:) :: ud, buf
    integer,intent(IN) :: ndir, lr
    real(kind=8),dimension(lbound(buf,1):ubound(buf,1),lbound(buf,2):ubound(buf,2),lbound(buf,3):ubound(buf,3),0:2) :: grad
    real(kind=8) :: hf, hc, di, dj, dk, dv, dxl, dxr, dyl, dyr, dzl, dzr
    integer :: i,j,k,m, ic,jc,kc,is,ie,js,je,ks,ke,ms,me,ibs,ibe,jbs,jbe,kbs,kbe,mbs,mbe
    integer :: ig0, jg0, kg0 
    real(kind=8) :: a, b, minmod
    minmod(a,b) = sign(1.d0,a)*max(0.d0,min(abs(a),sign(1.d0,a)*b))
    hf = 1.d0 
    hc = 2* hf 
    dv = 1.d0 
    is = Imingh
    ie = Imaxgh
    js = Jmingh
    je = Jmaxgh
    ks = Kmingh
    ke = Kmaxgh
    if (ndir == 0 .and. lr == Left ) ie = Imin-1
    if (ndir == 0 .and. lr == Right) is = Imax+1
    if (ndir == 1 .and. lr == Left ) je = Jmin-1
    if (ndir == 1 .and. lr == Right) js = Jmax+1
    if (ndir == 2 .and. lr == Left ) ke = Kmin-1
    if (ndir == 2 .and. lr == Right) ks = Kmax+1
    ibs = lbound(buf,1) 
 ibe = ubound(buf,1)
    jbs = lbound(buf,2) 
 jbe = ubound(buf,2)
    kbs = lbound(buf,3) 
 kbe = ubound(buf,3)
    mbs = lbound(buf,4) 
 mbe = ubound(buf,4)
    do m = mbs, mbe
       do k = kbs+1, kbe-1
          do j = jbs+1, jbe-1
             do i = ibs+1, ibe-1
                dxr = (buf(i+1,j,k,m)-buf(i,j,k,m))/hc
                dxl = (buf(i,j,k,m)-buf(i-1,j,k,m))/hc
                dyr = (buf(i,j+1,k,m)-buf(i,j,k,m))/hc
                dyl = (buf(i,j,k,m)-buf(i,j-1,k,m))/hc
                dzr = (buf(i,j,k+1,m)-buf(i,j,k,m))/hc
                dzl = (buf(i,j,k,m)-buf(i,j,k-1,m))/hc
                grad(i,j,k,0) = minmod(dxr, dxl)
                grad(i,j,k,1) = minmod(dyr, dyl)
                grad(i,j,k,2) = minmod(dzr, dzl)
             enddo
          enddo
       enddo
       do k = ks, ke
          do j = js, je
             do i = is, ie
                ic = ((i)-(is))/2 + mod(min((i)-(is),0),2) + (is) - is + ibs+1
                jc = ((j)-(js))/2 + mod(min((j)-(js),0),2) + (js) - js + jbs+1
                kc = ((k)-(ks))/2 + mod(min((k)-(ks),0),2) + (ks) - ks + kbs+1
                di = ( modulo(i,2) - 0.5d0 )*hf
                dj = ( modulo(j,2) - 0.5d0 )*hf
                dk = ( modulo(k,2) - 0.5d0 )*hf
                ud(i,j,k,m) = buf(ic,jc,kc,m)+grad(ic,jc,kc,0)*di+grad(ic,jc,kc,1)*dj+grad(ic,jc,kc,2)*dk
             enddo
          end do
       end do
    end do
  end subroutine interp_space
  subroutine interp_time(udbuf, lr, ndir, gids, gidd, rankd)
    use mpilib
    use grid
    use eos
    integer,intent(IN) :: lr, ndir, gids, gidd, rankd
    real(kind=8),pointer,dimension(:,:,:,:) :: u, udbuf
    integer :: ieo, jeo, keo
    call left_or_right(gidd, rankd, ieo, jeo, keo)
    select case ( STEP_MODE )
    case ( PREDICTOR )
       u => get_U2orderp(gids)
    case ( CORRECTOR )
       u => get_U1orderp(gids)
    case default
       u => get_Up(gids)
    end select
    udbuf = u(Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmax)
  end subroutine interp_time
end module grid_boundary
