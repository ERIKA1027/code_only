module fmg_reflux
  use fmg_data
  use mpilib
  implicit none
  private
  integer,parameter :: Send = 0, Recv = 1 
  integer,parameter :: Ev = 0, Od = 1 
  integer,save,dimension(Left:Right,0:2,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  integer,save,dimension(Left:Right,0:2) :: Ifs, Ife, Jfs, Jfe, Kfs, Kfe
  type type_surf
     real(kind=8),dimension(:,:,:,:),pointer :: ff => null(), fc => null()
  end type type_surf
  type(type_surf),save,dimension(:,:),allocatable :: Buf
  public :: fmg_fluxcorrection, fmg_reflux_finalize
contains
  subroutine fmg_fluxcorrection(fmglev)
    integer,intent(IN) :: fmglev
    call fmg_reflux_init(fmglev)
    call fmg_reflux_do(fmglev)
  end subroutine fmg_fluxcorrection
  subroutine fmg_reflux_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir
    Ics(:,:,Ev) = GridSize(fmglev)%Imin
    Ice(:,:,Ev) = (GridSize(fmglev)%Imax+1)/2-1
    Ics(:,:,Od) = (GridSize(fmglev)%Imax+1)/2
    Ice(:,:,Od) = GridSize(fmglev)%Imax
    Jcs(:,:,Ev) = GridSize(fmglev)%Jmin
    Jce(:,:,Ev) = (GridSize(fmglev)%Jmax+1)/2-1
    Jcs(:,:,Od) = (GridSize(fmglev)%Jmax+1)/2
    Jce(:,:,Od) = GridSize(fmglev)%Jmax
    Kcs(:,:,Ev) = GridSize(fmglev)%Kmin
    Kce(:,:,Ev) = (GridSize(fmglev)%Kmax+1)/2-1
    Kcs(:,:,Od) = (GridSize(fmglev)%Kmax+1)/2
    Kce(:,:,Od) = GridSize(fmglev)%Kmax
    Ics(Left, 0,:) = GridSize(fmglev)%Imax
    Ice(Left, 0,:) = GridSize(fmglev)%Imax
    Ics(Right,0,:) = GridSize(fmglev)%Imin-1
    Ice(Right,0,:) = GridSize(fmglev)%Imin-1
    Jcs(Left, 1,:) = GridSize(fmglev)%Jmax
    Jce(Left, 1,:) = GridSize(fmglev)%Jmax
    Jcs(Right,1,:) = GridSize(fmglev)%Jmin-1
    Jce(Right,1,:) = GridSize(fmglev)%Jmin-1
    Kcs(Left, 2,:) = GridSize(fmglev)%Kmax
    Kce(Left, 2,:) = GridSize(fmglev)%Kmax
    Kcs(Right,2,:) = GridSize(fmglev)%Kmin-1
    Kce(Right,2,:) = GridSize(fmglev)%Kmin-1
    Ifs(:,:) = GridSize(fmglev)%Imin
    Ife(:,:) = GridSize(fmglev)%Imax
    Jfs(:,:) = GridSize(fmglev)%Jmin
    Jfe(:,:) = GridSize(fmglev)%Jmax
    Kfs(:,:) = GridSize(fmglev)%Kmin
    Kfe(:,:) = GridSize(fmglev)%Kmax
    Ifs(Left,0) = GridSize(fmglev)%Imin-1
    Ife(Left,0) = GridSize(fmglev)%Imin-1
    Ifs(Right,0) = GridSize(fmglev)%Imax
    Ife(Right,0) = GridSize(fmglev)%Imax
    Jfs(Left,1) = GridSize(fmglev)%Jmin-1
    Jfe(Left,1) = GridSize(fmglev)%Jmin-1
    Jfs(Right,1) = GridSize(fmglev)%Jmax
    Jfe(Right,1) = GridSize(fmglev)%Jmax
    Kfs(Left,2) = GridSize(fmglev)%Kmin-1
    Kfe(Left,2) = GridSize(fmglev)%Kmin-1
    Kfs(Right,2) = GridSize(fmglev)%Kmax
    Kfe(Right,2) = GridSize(fmglev)%Kmax
    if ( .not. allocated(Buf) ) allocate( Buf(0:2, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(Buf(0, fmglev)%ff) ) then
       do ndir = 0, 2
          allocate( Buf(ndir, fmglev)%ff(Ifs(Left,ndir):Ife(Left,ndir),Jfs(Left,ndir):Jfe(Left,ndir),Kfs(Left,ndir):Kfe(Left,ndir),&
&Mmin:Mmax) )
       enddo
    end if
    if ( .not. associated(Buf(0, fmglev)%fc) ) then
       do ndir = 0, 2
          allocate( Buf(ndir, fmglev)%fc(Ics(Left,ndir,Ev):Ice(Left,ndir,Ev),Jcs(Left,ndir,Ev):Jce(Left,ndir,Ev),Kcs(Left,ndir,Ev):&
&Kce(Left,ndir,Ev),Mmin:Mmax) )
       enddo
    end if
  end subroutine fmg_reflux_init
  subroutine fmg_reflux_finalize
    integer :: ndir, fmglev
    if (allocated(Buf)) then
       do fmglev = lbound(Buf,2), ubound(Buf,2)
          do ndir = lbound(Buf,1), ubound(Buf,1)
             if (associated(Buf(ndir, fmglev)%ff)) deallocate(Buf(ndir, fmglev)%ff)
             if (associated(Buf(ndir, fmglev)%fc)) deallocate(Buf(ndir, fmglev)%fc)
             nullify(Buf(ndir, fmglev)%ff, Buf(ndir, fmglev)%fc)
          end do
       enddo
       deallocate( Buf )
    endif
  end subroutine fmg_reflux_finalize
  subroutine fmg_reflux_do(fmglev)
    use packarr
    integer,intent(IN) :: fmglev
    integer :: la, gid, gids, gidd, lr, ndir, rank, rankd, ranks, ieo, jeo, keo
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fs, fd
    call pkar_reset
    myrank = get_myrank()
    do rank = 0, 400 -1 
       do la = AMR_LevelMin+1, AMR_LevelMax
          do gid = fmg_get_gidmin_rank(la, rank), fmg_get_gidmax_rank(la, rank)
             do ndir = 0, 2
                do lr = Left, Right
                   if ( .not. Ancestry(la,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
                   gids = gid
                   ranks = rank
                   gidd = Ancestry(la,rank)%Block(gid)%NeighborGid(lr,ndir)
                   rankd = Ancestry(la,rank)%Block(gid)%NeighborRank(lr,ndir)
                   if (myrank == ranks) then
                      fs => fmg_get_fp(la, fmglev, gids)
                      Buf(ndir,fmglev)%ff = fs(Ifs(lr,ndir):Ife(lr,ndir),Jfs(lr,ndir):Jfe(lr,ndir),Kfs(lr,ndir):Kfe(lr,ndir),ndir,:&
&)
                      call restrict(Buf(ndir,fmglev)%ff, Buf(ndir,fmglev)%fc, ndir)
                   end if
                   if ((myrank) == (ranks)) call pkar_push(Buf(ndir,fmglev)%fc(lbound(Buf(ndir,fmglev)%fc,1),lbound(Buf(ndir,fmglev&
&)%fc,2),lbound(Buf(ndir,fmglev)%fc,3),lbound(Buf(ndir,fmglev)%fc,4)), size(Buf(ndir,fmglev)%fc), kind(Buf(ndir,fmglev)%fc), rankd)&
& 
 if ((myrank) == (rankd)) call pkar_recvlen(size(Buf(ndir,fmglev)%fc), kind(Buf(ndir,fmglev)%fc), ranks)
                enddo
             enddo
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do la = AMR_LevelMin+1, AMR_LevelMax
          do gid = fmg_get_gidmin_rank(la, rank), fmg_get_gidmax_rank(la, rank)
             do ndir = 0, 2
                do lr = Left, Right
                   if ( .not. Ancestry(la,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
                   gids = gid
                   ranks = rank
                   gidd = Ancestry(la,rank)%Block(gid)%NeighborGid(lr,ndir)
                   rankd = Ancestry(la,rank)%Block(gid)%NeighborRank(lr,ndir)
                   if (myrank == rankd) then
                      call fmg_left_or_right(gids, ranks, la, ieo, jeo, keo)
                      fd => fmg_get_fp(la-1, fmglev, gidd)
                      if ((myrank) == (rankd)) call pkar_pop(Buf(ndir,fmglev)%fc(lbound(Buf(ndir,fmglev)%fc,1),lbound(Buf(ndir,fmgl&
&ev)%fc,2),lbound(Buf(ndir,fmglev)%fc,3),lbound(Buf(ndir,fmglev)%fc,4)), size(Buf(ndir,fmglev)%fc), kind(Buf(ndir,fmglev)%fc), rank&
&s)
                      fd(Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),ndir&
&,:) = Buf(ndir,fmglev)%fc
                  end if
                end do
             enddo
          enddo
       enddo
    enddo
  end subroutine fmg_reflux_do
  subroutine restrict(ff,fc,ndir)
    real(kind=8),dimension(:,:,:,:),pointer :: ff, fc
    integer,intent(IN) :: ndir
    integer :: ic, jc, kc, if, jf, kf, is, js, ks, ie, je, ke, if0, jf0, kf0, m
    is = lbound(fc, 1)
    js = lbound(fc, 2)
    ks = lbound(fc, 3)
    ie = ubound(fc, 1)
    je = ubound(fc, 2)
    ke = ubound(fc, 3)
    if0 = lbound(ff, 1)
    jf0 = lbound(ff, 2)
    kf0 = lbound(ff, 3)
    if (ndir == 0) then
       if = if0
       ic = is
       do m = lbound(fc,4), ubound(fc,4)
          do kc=ks,ke
             do jc=js,je
                kf = ((kc) - (ks)) * 2 + (ks) - ks + kf0
                jf = ((jc) - (js)) * 2 + (js) - js + jf0
                fc(ic,jc,kc,m) = &
                     (ff(if, jf, kf,m) + ff(if, jf, kf+1,m) &
                     +ff(if, jf+1,kf,m) + ff(if, jf+1,kf+1,m) )*0.25d0
             end do
          end do
       end do
    elseif (ndir == 1) then
       jf = jf0
       jc = js
       do m = lbound(fc,4), ubound(fc,4)
          do kc=ks,ke
             do ic=is,ie
                kf = ((kc) - (ks)) * 2 + (ks) - ks + kf0
                if = ((ic) - (is)) * 2 + (is) - is + if0
                fc(ic,jc,kc,m) = &
                     (ff(if, jf, kf,m) + ff(if, jf, kf+1,m) &
                     +ff(if+1,jf, kf,m) + ff(if+1,jf, kf+1,m) )*0.25d0
             end do
          end do
       end do
    elseif (ndir == 2) then
       kf = kf0
       kc = ks
       do m = lbound(fc,4), ubound(fc,4)
          do jc=js,je
             do ic=is,ie
                jf = ((jc) - (js)) * 2 + (js) - js + jf0
                if = ((ic) - (is)) * 2 + (is) - is + if0
                fc(ic,jc,kc,m) = &
                     (ff(if, jf, kf,m) + ff(if, jf+1, kf,m) &
                     +ff(if+1,jf, kf,m) + ff(if+1,jf+1, kf,m) )*0.25d0
             end do
          end do
       end do
    endif
  end subroutine restrict
end module fmg_reflux
