module fmg_ghostcell
  use fmg_data
  use mpilib
  implicit none
  private
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(:,:),allocatable :: FLSV 
  type(t_buf),save,dimension(:,:),allocatable :: FLPV 
  type(t_buf),save,dimension(:,:),allocatable :: FLSS 
  type(t_buf),save,dimension(:,:),allocatable :: FLPS 
  type(t_buf),save,dimension(:,:,:),allocatable :: BLV, BLS 
  type(t_buf),save,dimension(:,:,:),allocatable :: CIXV, CIXYV, CIXYZV 
  type(t_buf),save,dimension(:,:,:),allocatable :: CIXS, CIXYS, CIXYZS
  integer,parameter :: Send = 0, Recv = 1 
  integer,parameter :: Ev = 0, Od = 1 
  integer,save,dimension(Left:Right,0:2,Send:Recv) :: Ins, Ine, Jns, Jne, Kns, Kne
  integer,save,dimension(Left:Right,0:2,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  integer,save :: NghCg 
  integer,parameter :: NghFg = 1 
  public :: fmg_ghostcell_fix, fmg_ghfix_samelev, fmg_ghfix_parentlev, fmg_ghfix_samelev_init, fmg_ghfix_parentlev_init, fmg_ghostc&
&ell_finalize
contains
  subroutine fmg_ghostcell_fix(fmglev,icode, cubic, tricubic)
    integer,intent(IN) :: fmglev, icode
    logical,intent(IN),optional :: cubic, tricubic
    integer :: ndir, amrlev
    call fmg_ghfix_samelev_init(fmglev)
    call fmg_ghfix_parentlev_init(fmglev)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do ndir = 0, 2
          call fmg_ghfix_parentlev(amrlev, fmglev, ndir, icode, cubic, tricubic)
          call fmg_ghfix_samelev(amrlev, fmglev,ndir,icode)
       end do
    enddo
  end subroutine fmg_ghostcell_fix
  subroutine fmg_ghfix_samelev_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir, lr
    Ins(:,:,:) = fmg_get_imingh(fmglev)
    Ine(:,:,:) = fmg_get_imaxgh(fmglev)
    Jns(:,:,:) = fmg_get_jmingh(fmglev)
    Jne(:,:,:) = fmg_get_jmaxgh(fmglev)
    Kns(:,:,:) = fmg_get_kmingh(fmglev)
    Kne(:,:,:) = fmg_get_kmaxgh(fmglev)
    Ins(Left ,0,Send) = GridSize(fmglev)%Imax-Ngh+1 
    Ine(Left ,0,Send) = GridSize(fmglev)%Imax
    Ins(Right,0,Send) = GridSize(fmglev)%Imin
    Ine(Right,0,Send) = GridSize(fmglev)%Imin+Ngh-1
    Jns(Left ,1,Send) = GridSize(fmglev)%Jmax-Ngh+1
    Jne(Left ,1,Send) = GridSize(fmglev)%Jmax
    Jns(Right,1,Send) = GridSize(fmglev)%Jmin
    Jne(Right,1,Send) = GridSize(fmglev)%Jmin+Ngh-1
    Kns(Left ,2,Send) = GridSize(fmglev)%Kmax-Ngh+1
    Kne(Left ,2,Send) = GridSize(fmglev)%Kmax
    Kns(Right,2,Send) = GridSize(fmglev)%Kmin
    Kne(Right,2,Send) = GridSize(fmglev)%Kmin+Ngh-1
    Ins(Left ,0,Recv) = GridSize(fmglev)%Imin-Ngh 
    Ine(Left ,0,Recv) = GridSize(fmglev)%Imin-1
    Ins(Right,0,Recv) = GridSize(fmglev)%Imax+1
    Ine(Right,0,Recv) = GridSize(fmglev)%Imax+Ngh
    Jns(Left ,1,Recv) = GridSize(fmglev)%Jmin-Ngh
    Jne(Left ,1,Recv) = GridSize(fmglev)%Jmin-1
    Jns(Right,1,Recv) = GridSize(fmglev)%Jmax+1
    Jne(Right,1,Recv) = GridSize(fmglev)%Jmax+Ngh
    Kns(Left ,2,Recv) = GridSize(fmglev)%Kmin-Ngh
    Kne(Left ,2,Recv) = GridSize(fmglev)%Kmin-1
    Kns(Right,2,Recv) = GridSize(fmglev)%Kmax+1
    Kne(Right,2,Recv) = GridSize(fmglev)%Kmax+Ngh
    if ( .not. allocated(FLSV) ) allocate( FLSV(0:2,FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLSV(0,fmglev)%Bufc) ) then
       lr = Left
       do ndir = 0, 2
          allocate( FLSV(ndir,fmglev)%Bufc( Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Sen&
&d):Kne(lr,ndir,Send),Mmin:Mmax ) )
       end do
    endif
    if ( .not. allocated(FLSS) ) allocate( FLSS(0:2,FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLSS(0,fmglev)%Bufc) ) then
       lr = Left
       do ndir = 0, 2
          allocate( FLSS(ndir,fmglev)%Bufc( Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Sen&
&d):Kne(lr,ndir,Send),Mmin:Mmin ) )
       end do
    endif
  end subroutine fmg_ghfix_samelev_init
  subroutine fmg_ghfix_samelev(amrlev, fmglev,ndir,icode)
    use packarr
    integer,intent(IN) :: amrlev, fmglev, ndir, icode
    integer :: lr, rank, rankd, ranks, gid, gidd, gids
    real(kind=8),pointer,dimension(:,:,:,:) :: us, ud, buf
    myrank = get_myrank()
    call pkar_reset
    if ( fmg_isVector(icode) ) then
       buf => FLSV(ndir,fmglev)%Bufc
    else
       buf => FLSS(ndir,fmglev)%Bufc
    endif
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborSameLevel(lr,ndir) ) cycle 
             gidd = gid
             rankd = rank
             gids = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             if ( myrank == rankd .and. myrank == ranks ) then 
                us => fmg_get_arrp(amrlev, fmglev, gids, icode)
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),:) =&
& us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndir,Send),:)
             else if ( myrank == ranks ) then 
                us => fmg_get_arrp(amrlev, fmglev, gids, icode)
                buf = us(Ins(lr,ndir,Send):Ine(lr,ndir,Send),Jns(lr,ndir,Send):Jne(lr,ndir,Send),Kns(lr,ndir,Send):Kne(lr,ndir,Send&
&),:)
             endif
             if ( (myrank == ranks .or. myrank == rankd) .and. rankd /= ranks) then
                if ((myrank) == (ranks)) call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), ki&
&nd(buf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(buf), kind(buf), ranks)
             endif
          end do
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1 
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborSameLevel(lr,ndir) ) cycle 
             gidd = gid
             rankd = rank
             gids = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             if ( myrank == ranks ) cycle
             if ( myrank == rankd ) then 
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), kin&
&d(buf), ranks)
                ud(Ins(lr,ndir,Recv):Ine(lr,ndir,Recv),Jns(lr,ndir,Recv):Jne(lr,ndir,Recv),Kns(lr,ndir,Recv):Kne(lr,ndir,Recv),:) =&
& buf
             endif
          enddo
       enddo
    enddo
  end subroutine fmg_ghfix_samelev
  subroutine fmg_ghfix_parentlev_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir, lr, ieo, jeo, keo
    integer :: ifs, ife, jfs, jfe, kfs, kfe, ibs, ibe, jbs, jbe, kbs, kbe
    real(kind=8),dimension(:,:,:,:),pointer :: buf
    NghCg = Ngh
    if (Ngh > 2) then
       print *, '*** fmg_ghfix_parentlev_init: error large Ngh is not support. Ngh =', Ngh
       stop
    end if
    Ics(:,:,Ev) = fmg_get_imin(fmglev)-Ngh
    Ice(:,:,Ev) = (GridSize(fmglev)%Imax+1)/2-1+Ngh
    Ics(:,:,Od) = (GridSize(fmglev)%Imax+1)/2-Ngh
    Ice(:,:,Od) = fmg_get_imax(fmglev)+Ngh
    Jcs(:,:,Ev) = fmg_get_jmin(fmglev)-Ngh
    Jce(:,:,Ev) = (GridSize(fmglev)%Jmax+1)/2-1+Ngh
    Jcs(:,:,Od) = (GridSize(fmglev)%Jmax+1)/2-Ngh
    Jce(:,:,Od) = fmg_get_jmax(fmglev)+Ngh
    Kcs(:,:,Ev) = fmg_get_kmin(fmglev)-Ngh
    Kce(:,:,Ev) = (GridSize(fmglev)%Kmax+1)/2-1+Ngh
    Kcs(:,:,Od) = (GridSize(fmglev)%Kmax+1)/2-Ngh
    Kce(:,:,Od) = fmg_get_kmax(fmglev)+Ngh
    Ics(Left, 0,:) = GridSize(fmglev)%Imax-NghCg+1
    Ice(Left, 0,:) = GridSize(fmglev)%Imax
    Ics(Right,0,:) = GridSize(fmglev)%Imin
    Ice(Right,0,:) = GridSize(fmglev)%Imin+NghCg-1
    Jcs(Left, 1,:) = GridSize(fmglev)%Jmax-NghCg+1
    Jce(Left, 1,:) = GridSize(fmglev)%Jmax
    Jcs(Right,1,:) = GridSize(fmglev)%Jmin
    Jce(Right,1,:) = GridSize(fmglev)%Jmin+NghCg-1
    Kcs(Left, 2,:) = GridSize(fmglev)%Kmax-NghCg+1
    Kce(Left, 2,:) = GridSize(fmglev)%Kmax
    Kcs(Right,2,:) = GridSize(fmglev)%Kmin
    Kce(Right,2,:) = GridSize(fmglev)%Kmin+NghCg-1
    if ( .not. allocated(FLPV) ) allocate( FLPV(0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLPV(0, fmglev)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = 0, 2
          allocate( FLPV(ndir, fmglev)%Bufc( Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):K&
&ce(lr,ndir,keo),Mmin:Mmax ) )
       end do
    endif
    if ( .not. allocated(FLPS) ) allocate( FLPS(0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLPS(0, fmglev)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = 0, 2
          allocate( FLPS(ndir, fmglev)%Bufc( Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):K&
&ce(lr,ndir,keo),Mmin:Mmin ) )
       end do
    endif
    if ( .not. allocated(BLS) ) allocate( BLS (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(BLV) ) allocate( BLV (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXV) ) allocate( CIXV (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYV) ) allocate( CIXYV (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYZV) ) allocate( CIXYZV(Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXS) ) allocate( CIXS (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYS) ) allocate( CIXYS (Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYZS) ) allocate( CIXYZS(Left:Right, 0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(CIXV(Left, 0, fmglev)%Bufc )) then
       do ndir = 0, 2
          do lr = Left, Right
             buf => FLPV(ndir, fmglev)%Bufc
             ibs = lbound(buf,1)
             ibe = ubound(buf,1)
             jbs = lbound(buf,2)
             jbe = ubound(buf,2)
             kbs = lbound(buf,3)
             kbe = ubound(buf,3)
             call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe)
             allocate( BLS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmin))
             allocate( BLV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmax))
             allocate( CIXV(lr, ndir, fmglev)%Bufc(ifs:ife,jbs:jbe,kbs:kbe,Mmin:Mmax))
             allocate( CIXYV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kbs:kbe,Mmin:Mmax))
             allocate(CIXYZV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmax))
             allocate( CIXS(lr, ndir, fmglev)%Bufc(ifs:ife,jbs:jbe,kbs:kbe,Mmin:Mmin))
             allocate( CIXYS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kbs:kbe,Mmin:Mmin))
             allocate(CIXYZS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmin))
          end do
       end do
    end if
  end subroutine fmg_ghfix_parentlev_init
  subroutine fmg_ghfix_parentlev(amrlev, fmglev, ndir, icode, cubic, tricubic)
    use packarr
    integer,intent(IN) :: amrlev, fmglev, ndir, icode
    logical,intent(IN),optional :: cubic, tricubic
    integer :: gid, gids, gidd, lr, rank, rankd, ranks, ieo, jeo, keo
    real(kind=8),pointer,dimension(:,:,:,:) :: us, ud, udbuf
    logical :: bool_cubic
    if ( amrlev <= AMR_LevelMin ) return
    bool_cubic = .FALSE. 
 if (present(cubic)) bool_cubic = cubic
    if (present(tricubic)) then
       if (tricubic) bool_cubic = .TRUE.
    end if
    call pkar_reset
    if ( fmg_isVector(icode) ) then
       udbuf => FLPV(ndir, fmglev)%Bufc
    else
       udbuf => FLPS(ndir, fmglev)%Bufc
    endif
    myrank = get_myrank()
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
             gids = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             gidd = gid
             rankd = rank
             if (myrank == ranks) then
                call fmg_left_or_right(gidd, rankd, amrlev, ieo, jeo, keo)
                us => fmg_get_arrp(amrlev-1, fmglev, gids, icode)
                udbuf = us(Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),:)
             end if
             if ((myrank) == (ranks)) call pkar_push(udbuf(lbound(udbuf,1),lbound(udbuf,2),lbound(udbuf,3),lbound(udbuf,4)), size(u&
&dbuf), kind(udbuf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(udbuf), kind(udbuf), ranks)
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
             gids = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             gidd = gid
             rankd = rank
             if (myrank == rankd) then
                call fmg_left_or_right(gidd, rankd, amrlev, ieo, jeo, keo)
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                if ((myrank) == (rankd)) call pkar_pop(udbuf(lbound(udbuf,1),lbound(udbuf,2),lbound(udbuf,3),lbound(udbuf,4)), size&
&(udbuf), kind(udbuf), ranks)
                if (bool_cubic .and. interpCubic_checkSize(udbuf, ndir)) then
                   call interpCubic(ud, udbuf, ndir, lr, fmglev, tricubic)
                else
                   call interp(ud, udbuf, ndir, lr, fmglev)
                end if
             end if
          enddo
       enddo
    enddo
  contains
    subroutine interp(ud, buf, ndir, lr, fmglev)
      real(kind=8),dimension(:,:,:,:),pointer :: buf
      real(kind=8),dimension(:,:,:,:),pointer :: ud
      integer,intent(IN) :: ndir, lr, fmglev
      integer :: ifs, ife, jfs, jfe, kfs, kfe, i, j, k, ic, jc, kc, ic1, jc1, kc1, m
      integer :: ibool, jbool, kbool, ixbool, jxbool, kxbool, inner
      integer :: iba, jba, kba
      real(kind=8) :: uH, u0, uN, uHH
      real(kind=8),dimension(:,:,:,:),pointer :: uintp
      if ( fmg_isVector(icode) ) then
         uintp => BLV(lr, ndir, fmglev)%Bufc
      else
         uintp => BLS(lr, ndir, fmglev)%Bufc
      endif
      ibool = 0
      jbool = 0
      kbool = 0
      if (ndir == 0) ibool = 1
      if (ndir == 1) jbool = 1
      if (ndir == 2) kbool = 1
      ixbool = 1-ibool
      jxbool = 1-jbool
      kxbool = 1-kbool
      if (lr == Left) inner = 1
      if (lr == Right) inner = -1
      iba = lbound(buf,1)
      jba = lbound(buf,2)
      kba = lbound(buf,3)
      call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe) 
      do m = lbound(buf, 4), ubound(buf, 4)
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  ic = (((i)-(GridSize(fmglev)%Imin))/2 + mod(min((i)-(GridSize(fmglev)%Imin),0),2) + (GridSize(fmglev)%Imin))*ixbo&
&ol + (i-ifs+iba)*ibool
                  jc = (((j)-(GridSize(fmglev)%Jmin))/2 + mod(min((j)-(GridSize(fmglev)%Jmin),0),2) + (GridSize(fmglev)%Jmin))*jxbo&
&ol + (j-jfs+jba)*jbool
                  kc = (((k)-(GridSize(fmglev)%Kmin))/2 + mod(min((k)-(GridSize(fmglev)%Kmin),0),2) + (GridSize(fmglev)%Kmin))*kxbo&
&ol + (k-kfs+kba)*kbool
                  ic1 = ( (ic)+2*modulo((i),2)-1 ) * ixbool + ic * ibool
                  jc1 = ( (jc)+2*modulo((j),2)-1 ) * jxbool + jc * jbool
                  kc1 = ( (kc)+2*modulo((k),2)-1 ) * kxbool + kc * kbool
                  uintp(i,j,k,m) &
                       = ibool * (9*buf(ic,jc,kc,m) + buf(ic,jc1,kc1,m) + 3*buf(ic,jc,kc1,m) + 3*buf(ic,jc1,kc,m) )/16 &
                       + jbool * (9*buf(ic,jc,kc,m) + buf(ic1,jc,kc1,m) + 3*buf(ic,jc,kc1,m) + 3*buf(ic1,jc,kc,m) )/16 &
                       + kbool * (9*buf(ic,jc,kc,m) + buf(ic1,jc1,kc,m) + 3*buf(ic,jc1,kc,m) + 3*buf(ic1,jc,kc,m) )/16
               end do
            end do
         end do
      end do
      if (ndir == 0 .and. lr == Left) ifs = ife
      if (ndir == 0 .and. lr == Right) ife = ifs
      if (ndir == 1 .and. lr == Left) jfs = jfe
      if (ndir == 1 .and. lr == Right) jfe = jfs
      if (ndir == 2 .and. lr == Left) kfs = kfe
      if (ndir == 2 .and. lr == Right) kfe = kfs
      do m = lbound(buf,4), ubound(buf,4)
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  uH = uintp(i,j,k,m)
                  u0 = ud(i+inner*ibool, j+inner*jbool, k+inner*kbool, m) 
                  uN = ud(i+inner*ibool*2, j+inner*jbool*2, k+inner*kbool*2, m) 
                  ud(i,j,k, m) = ( 10 * u0 + 8 * uH - 3 * uN ) / 15 
                  if (Ngh == 2) then
                     uHH = uintp(i-inner*ibool, j-inner*jbool, k-inner*kbool, m)
                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) &
                          = 2.d0/15.0*uN - 3.d0/7.d0*u0 + 6.d0/5.d0*uH + 2.d0/21.d0*uHH
                  end if
               enddo
            end do
         end do
      end do
    end subroutine interp
    subroutine interpCubic(ud, buf, ndir, lr, fmglev, tricubic)
      use interpolation
      real(kind=8),dimension(:,:,:,:),pointer :: buf
      real(kind=8),dimension(:,:,:,:),pointer :: ud
      integer,intent(IN) :: ndir, lr, fmglev
      logical,intent(IN),optional :: tricubic
      integer :: ifs, ife, jfs, jfe, kfs, kfe, ibs, ibe, jbs, jbe, kbs, kbe, mbs, mbe,&
           i, j, k, m, ic, jc, kc, ic0, jc0, kc0, ic1, jc1, kc1, ic2, jc2, kc2, ic3, jc3, kc3, &
           icP, jcP, kcP, icN, jcN, kcN, nshift
      integer :: ibool, jbool, kbool, inner
      real(kind=8) :: uH, u0, uN, uHH, xc, yc, zc, a0, a1, a2, a3
      real(kind=8),dimension(:,:,:,:),pointer :: bufx, bufxy, bufxyz
      logical :: bool_tricubic
      bool_tricubic = .FALSE. 
 if (present(tricubic)) bool_tricubic = tricubic
      if ( fmg_isVector(icode) ) then
         bufx => CIXV (lr, ndir, fmglev)%Bufc
         bufxy => CIXYV (lr, ndir, fmglev)%Bufc
         bufxyz => CIXYZV(lr, ndir, fmglev)%Bufc
      else
         bufx => CIXS (lr, ndir, fmglev)%Bufc
         bufxy => CIXYS (lr, ndir, fmglev)%Bufc
         bufxyz => CIXYZS(lr, ndir, fmglev)%Bufc
      endif
      ibs = lbound(buf,1)
      ibe = ubound(buf,1)
      jbs = lbound(buf,2)
      jbe = ubound(buf,2)
      kbs = lbound(buf,3)
      kbe = ubound(buf,3)
      mbs = lbound(buf,4)
      mbe = ubound(buf,4)
      ibool = 0
      jbool = 0
      kbool = 0
      if (ndir == 0) ibool = 1
      if (ndir == 1) jbool = 1
      if (ndir == 2) kbool = 1
      if (lr == Left) inner = 1
      if (lr == Right) inner = -1
      call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe) 
      if (ndir == 0) then
         bufx = buf
      else
         do m = mbs, mbe
            do kc = kbs, kbe
               do jc = jbs, jbe
                  do i = ifs, ife
                     icP = ((i)-(GridSize(fmglev)%Imin))/2 + mod(min((i)-(GridSize(fmglev)%Imin),0),2) + (GridSize(fmglev)%Imin)
                     icN = ( (icP)+2*modulo((i),2)-1 )
                     ic1 = min(icP, icN)
                     ic2 = max(icP, icN)
                     ic0 = 2*ic1 - ic2
                     ic3 = 2*ic2 - ic1
                     nshift = 0
                     if (ic3 > ubound(buf,1)) nshift = ic3 - ubound(buf,1)
                     if (ic0 < lbound(buf,1)) nshift = ic0 - lbound(buf,1)
                     ic0 = ic0 - nshift
                     ic1 = ic1 - nshift
                     ic2 = ic2 - nshift
                     ic3 = ic3 - nshift
                     xc = icP +0.25d0*(2*modulo(i,2)-1)
                     a0 = (xc - ic1)*(xc - ic2)*(xc - ic3)/((ic0 - ic1)*(ic0 - ic2)*(ic0 - ic3))
                     a1 = (xc - ic0)*(xc - ic2)*(xc - ic3)/((ic1 - ic0)*(ic1 - ic2)*(ic1 - ic3))
                     a2 = (xc - ic0)*(xc - ic1)*(xc - ic3)/((ic2 - ic0)*(ic2 - ic1)*(ic2 - ic3))
                     a3 = (xc - ic0)*(xc - ic1)*(xc - ic2)/((ic3 - ic0)*(ic3 - ic1)*(ic3 - ic2))
                     bufx(i, jc, kc, m) = &
                          a0*buf(ic0, jc, kc, m) + &
                          a1*buf(ic1, jc, kc, m) + &
                          a2*buf(ic2, jc, kc, m) + &
                          a3*buf(ic3, jc, kc, m)
                  end do
               end do
            end do
         end do
      end if
      if (ndir == 1) then
         bufxy = bufx
      else
         do m = mbs, mbe
            do kc = kbs, kbe
               do j = jfs, jfe
                  jcP = ((j)-(GridSize(fmglev)%Jmin))/2 + mod(min((j)-(GridSize(fmglev)%Jmin),0),2) + (GridSize(fmglev)%Jmin)
                  jcN = ( (jcP)+2*modulo((j),2)-1 )
                  jc1 = min(jcP, jcN)
                  jc2 = max(jcP, jcN)
                  jc0 = 2*jc1 - jc2
                  jc3 = 2*jc2 - jc1
                  nshift = 0
                  if (jc3 > ubound(buf,2)) nshift = jc3 - ubound(buf,2)
                  if (jc0 < lbound(buf,2)) nshift = jc0 - lbound(buf,2)
                  jc0 = jc0 - nshift
                  jc1 = jc1 - nshift
                  jc2 = jc2 - nshift
                  jc3 = jc3 - nshift
                  yc = jcP +0.25d0*(2*modulo(j,2)-1)
                  a0 = (yc - jc1)*(yc - jc2)*(yc - jc3)/((jc0 - jc1)*(jc0 - jc2)*(jc0 - jc3))
                  a1 = (yc - jc0)*(yc - jc2)*(yc - jc3)/((jc1 - jc0)*(jc1 - jc2)*(jc1 - jc3))
                  a2 = (yc - jc0)*(yc - jc1)*(yc - jc3)/((jc2 - jc0)*(jc2 - jc1)*(jc2 - jc3))
                  a3 = (yc - jc0)*(yc - jc1)*(yc - jc2)/((jc3 - jc0)*(jc3 - jc1)*(jc3 - jc2))
                  do i = ifs, ife
                     bufxy(i, j, kc, m) = &
                          a0*bufx(i, jc0, kc, m) + &
                          a1*bufx(i, jc1, kc, m) + &
                          a2*bufx(i, jc2, kc, m) + &
                          a3*bufx(i, jc3, kc, m)
                  end do
               end do
            end do
         end do
      end if
      if (ndir == 2) then
         bufxyz = bufxy
      else
         do m = mbs, mbe
            do k = kfs, kfe
               kcP = ((k)-(GridSize(fmglev)%Kmin))/2 + mod(min((k)-(GridSize(fmglev)%Kmin),0),2) + (GridSize(fmglev)%Kmin)
               kcN = ( (kcP)+2*modulo((k),2)-1 )
               kc1 = min(kcP, kcN)
               kc2 = max(kcP, kcN)
               kc0 = 2*kc1 - kc2
               kc3 = 2*kc2 - kc1
               nshift = 0
               if (kc3 > ubound(buf,3)) nshift = kc3 - ubound(buf,3)
               if (kc0 < lbound(buf,3)) nshift = kc0 - lbound(buf,3)
               kc0 = kc0 - nshift
               kc1 = kc1 - nshift
               kc2 = kc2 - nshift
               kc3 = kc3 - nshift
               zc = kcP +0.25d0*(2*modulo(k,2)-1)
               a0 = (zc - kc1)*(zc - kc2)*(zc - kc3)/((kc0 - kc1)*(kc0 - kc2)*(kc0 - kc3))
               a1 = (zc - kc0)*(zc - kc2)*(zc - kc3)/((kc1 - kc0)*(kc1 - kc2)*(kc1 - kc3))
               a2 = (zc - kc0)*(zc - kc1)*(zc - kc3)/((kc2 - kc0)*(kc2 - kc1)*(kc2 - kc3))
               a3 = (zc - kc0)*(zc - kc1)*(zc - kc2)/((kc3 - kc0)*(kc3 - kc1)*(kc3 - kc2))
               do j = jfs, jfe
                  do i = ifs, ife
                     bufxyz(i, j, k, m) = &
                          a0*bufxy(i, j, kc0, m) + &
                          a1*bufxy(i, j, kc1, m) + &
                          a2*bufxy(i, j, kc2, m) + &
                          a3*bufxy(i, j, kc3, m)
                  end do
               end do
            end do
         end do
      end if
      if (ndir == 0 .and. lr == Left) ifs = ife
      if (ndir == 0 .and. lr == Right) ife = ifs
      if (ndir == 1 .and. lr == Left) jfs = jfe
      if (ndir == 1 .and. lr == Right) jfe = jfs
      if (ndir == 2 .and. lr == Left) kfs = kfe
      if (ndir == 2 .and. lr == Right) kfe = kfs
      do m = mbs, mbe
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  uH = bufxyz(i,j,k,m)
                  u0 = ud(i+inner*ibool, j+inner*jbool, k+inner*kbool, m) 
                  uN = ud(i+inner*ibool*2, j+inner*jbool*2, k+inner*kbool*2, m) 
                  ud(i,j,k, m) = ( 10 * u0 + 8 * uH - 3 * uN ) / 15 
                  if (Ngh == 2) then
                     uHH = bufxyz(i-inner*ibool, j-inner*jbool, k-inner*kbool, m)
                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) &
                          = 2.d0/15.0*uN - 3.d0/7.d0*u0 + 6.d0/5.d0*uH + 2.d0/21.d0*uHH
                     if (bool_tricubic) &
                          ud(i,j,k,m) = -1.d0/9.d0*uN + 10.d0/21.d0*u0 + 2.d0/3.d0*uH -2.d0/63.d0*uHH
                  end if
               enddo
            end do
         end do
      end do
    end subroutine interpCubic
    function interpCubic_checkSize(buf, ndir) result(validSize)
      real(kind=8),dimension(:,:,:,:),pointer :: buf
      integer,intent(IN) :: ndir
      logical :: validSize
      integer,parameter :: MINBLOCKSIZE = 4
      if (ndir == 0) then
         validSize = size(buf,2) >= MINBLOCKSIZE .and. size(buf,3) >= MINBLOCKSIZE
      elseif (ndir == 1) then
         validSize = size(buf,3) >= MINBLOCKSIZE .and. size(buf,1) >= MINBLOCKSIZE
      elseif (ndir == 2) then
         validSize = size(buf,1) >= MINBLOCKSIZE .and. size(buf,2) >= MINBLOCKSIZE
      else
         print *, '**** error in fmg_ghostcell::interpCubic_checkSize', ndir
         stop
      end if
    end function interpCubic_checkSize
  end subroutine fmg_ghfix_parentlev
  subroutine fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe)
    integer,intent(IN) :: ndir, lr, fmglev
    integer,intent(OUT) :: ifs, jfs, kfs, ife, jfe, kfe
    call fmg_get_gridsizeGh(fmglev, ifs, jfs, kfs, ife, jfe, kfe)
    if (ndir == 0 .and. lr == Left) ife = GridSize(fmglev)%Imin-1
    if (ndir == 0 .and. lr == Right) ifs = GridSize(fmglev)%Imax+1
    if (ndir == 1 .and. lr == Left) jfe = GridSize(fmglev)%Jmin-1
    if (ndir == 1 .and. lr == Right) jfs = GridSize(fmglev)%Jmax+1
    if (ndir == 2 .and. lr == Left) kfe = GridSize(fmglev)%Kmin-1
    if (ndir == 2 .and. lr == Right) kfs = GridSize(fmglev)%Kmax+1
  end subroutine fmg_ghfix_parentlev_get_fineindex
  subroutine fmg_ghostcell_finalize
    integer :: lr, ndir, fmglev
    if (allocated(FLSV)) then 
 do fmglev = lbound(FLSV, 2), ubound(FLSV, 2) 
 do ndir = lbound(FLSV, 1), ubound(FLSV, 1) 
 if (associated(FLSV(ndir,fmglev)%Bufc)) deallocate(FLSV(ndir,fmglev)%Bufc) 
 nullify(FLSV(ndir,fmglev)%Bufc) 
 end do 
 end do 
 deallocate(FLSV) 
 endif
    if (allocated(FLPV)) then 
 do fmglev = lbound(FLPV, 2), ubound(FLPV, 2) 
 do ndir = lbound(FLPV, 1), ubound(FLPV, 1) 
 if (associated(FLPV(ndir,fmglev)%Bufc)) deallocate(FLPV(ndir,fmglev)%Bufc) 
 nullify(FLPV(ndir,fmglev)%Bufc) 
 end do 
 end do 
 deallocate(FLPV) 
 endif
    if (allocated(FLSS)) then 
 do fmglev = lbound(FLSS, 2), ubound(FLSS, 2) 
 do ndir = lbound(FLSS, 1), ubound(FLSS, 1) 
 if (associated(FLSS(ndir,fmglev)%Bufc)) deallocate(FLSS(ndir,fmglev)%Bufc) 
 nullify(FLSS(ndir,fmglev)%Bufc) 
 end do 
 end do 
 deallocate(FLSS) 
 endif
    if (allocated(FLPS)) then 
 do fmglev = lbound(FLPS, 2), ubound(FLPS, 2) 
 do ndir = lbound(FLPS, 1), ubound(FLPS, 1) 
 if (associated(FLPS(ndir,fmglev)%Bufc)) deallocate(FLPS(ndir,fmglev)%Bufc) 
 nullify(FLPS(ndir,fmglev)%Bufc) 
 end do 
 end do 
 deallocate(FLPS) 
 endif
    if (allocated(BLV)) then 
 do fmglev = lbound(BLV, 3), ubound(BLV, 3) 
 do ndir = lbound(BLV, 2), ubound(BLV, 2) 
 do lr = lbound(BLV, 1), ubound(BLV, 1) 
 if (associated(BLV(lr, ndir,fmglev)%Bufc)) deallocate(BLV(lr, ndir,fmglev)%Bufc) 
 nullify(BLV(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(BLV) 
 endif
    if (allocated(BLS)) then 
 do fmglev = lbound(BLS, 3), ubound(BLS, 3) 
 do ndir = lbound(BLS, 2), ubound(BLS, 2) 
 do lr = lbound(BLS, 1), ubound(BLS, 1) 
 if (associated(BLS(lr, ndir,fmglev)%Bufc)) deallocate(BLS(lr, ndir,fmglev)%Bufc) 
 nullify(BLS(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(BLS) 
 endif
    if (allocated(CIXV)) then 
 do fmglev = lbound(CIXV, 3), ubound(CIXV, 3) 
 do ndir = lbound(CIXV, 2), ubound(CIXV, 2) 
 do lr = lbound(CIXV, 1), ubound(CIXV, 1) 
 if (associated(CIXV(lr, ndir,fmglev)%Bufc)) deallocate(CIXV(lr, ndir,fmglev)%Bufc) 
 nullify(CIXV(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXV) 
 endif
    if (allocated(CIXYV)) then 
 do fmglev = lbound(CIXYV, 3), ubound(CIXYV, 3) 
 do ndir = lbound(CIXYV, 2), ubound(CIXYV, 2) 
 do lr = lbound(CIXYV, 1), ubound(CIXYV, 1) 
 if (associated(CIXYV(lr, ndir,fmglev)%Bufc)) deallocate(CIXYV(lr, ndir,fmglev)%Bufc) 
 nullify(CIXYV(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXYV) 
 endif
    if (allocated(CIXYZV)) then 
 do fmglev = lbound(CIXYZV, 3), ubound(CIXYZV, 3) 
 do ndir = lbound(CIXYZV, 2), ubound(CIXYZV, 2) 
 do lr = lbound(CIXYZV, 1), ubound(CIXYZV, 1) 
 if (associated(CIXYZV(lr, ndir,fmglev)%Bufc)) deallocate(CIXYZV(lr, ndir,fmglev)%Bufc) 
 nullify(CIXYZV(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXYZV) 
 endif
    if (allocated(CIXS)) then 
 do fmglev = lbound(CIXS, 3), ubound(CIXS, 3) 
 do ndir = lbound(CIXS, 2), ubound(CIXS, 2) 
 do lr = lbound(CIXS, 1), ubound(CIXS, 1) 
 if (associated(CIXS(lr, ndir,fmglev)%Bufc)) deallocate(CIXS(lr, ndir,fmglev)%Bufc) 
 nullify(CIXS(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXS) 
 endif
    if (allocated(CIXYS)) then 
 do fmglev = lbound(CIXYS, 3), ubound(CIXYS, 3) 
 do ndir = lbound(CIXYS, 2), ubound(CIXYS, 2) 
 do lr = lbound(CIXYS, 1), ubound(CIXYS, 1) 
 if (associated(CIXYS(lr, ndir,fmglev)%Bufc)) deallocate(CIXYS(lr, ndir,fmglev)%Bufc) 
 nullify(CIXYS(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXYS) 
 endif
    if (allocated(CIXYZS)) then 
 do fmglev = lbound(CIXYZS, 3), ubound(CIXYZS, 3) 
 do ndir = lbound(CIXYZS, 2), ubound(CIXYZS, 2) 
 do lr = lbound(CIXYZS, 1), ubound(CIXYZS, 1) 
 if (associated(CIXYZS(lr, ndir,fmglev)%Bufc)) deallocate(CIXYZS(lr, ndir,fmglev)%Bufc) 
 nullify(CIXYZS(lr, ndir,fmglev)%Bufc) 
 end do
 end do
 end do 
 deallocate(CIXYZS) 
 endif
  end subroutine fmg_ghostcell_finalize
end module fmg_ghostcell
