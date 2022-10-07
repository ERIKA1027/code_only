
! module for flux conservation in FMG cycle
#include "config.h"
#include "packarr.h"
!-------------------------------------------------------------------------
! Module for refluxing of FMG
!-------------------------------------------------------------------------
module fmg_reflux
  use fmg_data
  use mpilib
  implicit none
  private
  integer,parameter :: Send = 0, Recv = 1 ! send reciev code
  integer,parameter :: Ev = 0, Od = 1     ! parity of cell index
  integer,save,dimension(Left:Right,MX:MZ,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  integer,save,dimension(Left:Right,MX:MZ) :: Ifs, Ife, Jfs, Jfe, Kfs, Kfe

  type type_surf
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ff => null(), fc => null()
  end type type_surf
  type(type_surf),save,dimension(:,:),allocatable :: Buf

  public :: fmg_fluxcorrection, fmg_reflux_finalize
contains
  !-------------------------------------------------------------------------
  ! main routine
  !-------------------------------------------------------------------------
  subroutine fmg_fluxcorrection(fmglev)
    integer,intent(IN) :: fmglev
    call fmg_reflux_init(fmglev)
    call fmg_reflux_do(fmglev)
  end subroutine fmg_fluxcorrection
  !-------------------------------------------------------------------------
  ! initialize module
  !-------------------------------------------------------------------------
  subroutine fmg_reflux_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir

    ! ics(LR surf, direction, parity of this grid)
    !     LR,      MX,MY,MZ,  Ev,Od
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

    Ics(Left, MX,:) = GridSize(fmglev)%Imax
    Ice(Left, MX,:) = GridSize(fmglev)%Imax
    Ics(Right,MX,:) = GridSize(fmglev)%Imin-1
    Ice(Right,MX,:) = GridSize(fmglev)%Imin-1

    Jcs(Left, MY,:) = GridSize(fmglev)%Jmax
    Jce(Left, MY,:) = GridSize(fmglev)%Jmax
    Jcs(Right,MY,:) = GridSize(fmglev)%Jmin-1
    Jce(Right,MY,:) = GridSize(fmglev)%Jmin-1

    Kcs(Left, MZ,:) = GridSize(fmglev)%Kmax
    Kce(Left, MZ,:) = GridSize(fmglev)%Kmax
    Kcs(Right,MZ,:) = GridSize(fmglev)%Kmin-1
    Kce(Right,MZ,:) = GridSize(fmglev)%Kmin-1

    ! ifs(LR,ndir), ife
    Ifs(:,:) = GridSize(fmglev)%Imin
    Ife(:,:) = GridSize(fmglev)%Imax
    Jfs(:,:) = GridSize(fmglev)%Jmin
    Jfe(:,:) = GridSize(fmglev)%Jmax
    Kfs(:,:) = GridSize(fmglev)%Kmin
    Kfe(:,:) = GridSize(fmglev)%Kmax

    Ifs(Left,MX)  = GridSize(fmglev)%Imin-1
    Ife(Left,MX)  = GridSize(fmglev)%Imin-1
    Ifs(Right,MX) = GridSize(fmglev)%Imax
    Ife(Right,MX) = GridSize(fmglev)%Imax

    Jfs(Left,MY)  = GridSize(fmglev)%Jmin-1
    Jfe(Left,MY)  = GridSize(fmglev)%Jmin-1
    Jfs(Right,MY) = GridSize(fmglev)%Jmax
    Jfe(Right,MY) = GridSize(fmglev)%Jmax

    Kfs(Left,MZ)  = GridSize(fmglev)%Kmin-1
    Kfe(Left,MZ)  = GridSize(fmglev)%Kmin-1
    Kfs(Right,MZ) = GridSize(fmglev)%Kmax
    Kfe(Right,MZ) = GridSize(fmglev)%Kmax

    ! buffer
    ! とりあえず、大きさだけ正確。添字のオフセットは不定。偶奇は正しい。
#define SZ Ifs(Left,ndir):Ife(Left,ndir),Jfs(Left,ndir):Jfe(Left,ndir),Kfs(Left,ndir):Kfe(Left,ndir),Mmin:Mmax
    if ( .not. allocated(Buf) ) allocate( Buf(MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(Buf(MX, fmglev)%ff) ) then
       do ndir = MX, MZ
          allocate( Buf(ndir, fmglev)%ff(SZ) )
       enddo
    end if
#undef SZ
#define SZ Ics(Left,ndir,Ev):Ice(Left,ndir,Ev),Jcs(Left,ndir,Ev):Jce(Left,ndir,Ev),Kcs(Left,ndir,Ev):Kce(Left,ndir,Ev),Mmin:Mmax
    if ( .not. associated(Buf(MX, fmglev)%fc) ) then
       do ndir = MX, MZ
          allocate( Buf(ndir, fmglev)%fc(SZ) )
       enddo
    end if
#undef SZ

  end subroutine fmg_reflux_init
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! fix ghost cell by neighbor grid at the parent level
  ! 親の袖にアクセスする可能性あり
  !-------------------------------------------------------------------------
  subroutine fmg_reflux_do(fmglev)
    use packarr
    integer,intent(IN) :: fmglev
    integer :: la, gid, gids, gidd, lr, ndir, rank, rankd, ranks, ieo, jeo, keo
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fs, fd

#define SZ Ifs(lr,ndir):Ife(lr,ndir),Jfs(lr,ndir):Jfe(lr,ndir),Kfs(lr,ndir):Kfe(lr,ndir)

    ! 初期化
    call pkar_reset
    myrank = get_myrank()
    do rank = 0, NPE -1         !子セル source
       do la = AMR_LevelMin+1, AMR_LevelMax
          do gid = fmg_get_gidmin_rank(la, rank), fmg_get_gidmax_rank(la, rank)
             do ndir = MX, MZ
                do lr = Left, Right
                   if ( .not. Ancestry(la,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
                   gids  = gid
                   ranks = rank
                   gidd  = Ancestry(la,rank)%Block(gid)%NeighborGid(lr,ndir)
                   rankd = Ancestry(la,rank)%Block(gid)%NeighborRank(lr,ndir)
                   if (myrank == ranks) then
                      fs => fmg_get_fp(la, fmglev, gids)
                      Buf(ndir,fmglev)%ff = fs(SZ,ndir,:)
                      call restrict(Buf(ndir,fmglev)%ff, Buf(ndir,fmglev)%fc, ndir)
                   end if
                   PACK_SEND4( Buf(ndir,fmglev)%fc, myrank, ranks, rankd )
                enddo
             enddo
          enddo
       enddo
    enddo
#undef SZ

    call pkar_sendrecv()

#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo)
    ! 受信バッファから荷ほどき
    do rank = 0, NPE -1
       do la = AMR_LevelMin+1, AMR_LevelMax
          do gid = fmg_get_gidmin_rank(la, rank), fmg_get_gidmax_rank(la, rank)
             do ndir = MX, MZ
                do lr = Left, Right
                   if ( .not. Ancestry(la,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
                   gids  = gid
                   ranks = rank
                   gidd  = Ancestry(la,rank)%Block(gid)%NeighborGid(lr,ndir)
                   rankd = Ancestry(la,rank)%Block(gid)%NeighborRank(lr,ndir)
                   if (myrank == rankd) then
                      call fmg_left_or_right(gids, ranks, la, ieo, jeo, keo)
                      fd => fmg_get_fp(la-1, fmglev, gidd)
                      UNPACK_RECV4(  Buf(ndir,fmglev)%fc, myrank, ranks, rankd )
                      fd(SZ,ndir,:) = Buf(ndir,fmglev)%fc
                  end if
                end do
             enddo
          enddo
       enddo
    enddo
#undef SZ
!!$    ! deallocate arrays
!!$    do rank = 0, NPE -1
!!$       if ( associated( bufs(rank)%pkg ) ) deallocate( bufs(rank)%pkg )
!!$       if ( associated( bufd(rank)%pkg ) ) deallocate( bufd(rank)%pkg )
!!$    enddo
  end subroutine fmg_reflux_do
  !-------------------------------------------------------------------------
  ! restrict flux between grid boundary
  !-------------------------------------------------------------------------
  subroutine restrict(ff,fc,ndir)
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ff, fc
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

    if (ndir == MX) then
       if = if0
       ic = is
       do m = lbound(fc,4), ubound(fc,4)
          do kc=ks,ke
             do jc=js,je
                kf = IJKF(kc, ks) - ks + kf0
                jf = IJKF(jc, js) - js + jf0
                fc(ic,jc,kc,m) = &
                     (ff(if,  jf,  kf,m) + ff(if,  jf,  kf+1,m) &
                     +ff(if,  jf+1,kf,m) + ff(if,  jf+1,kf+1,m) )*0.25d0
             end do
          end do
       end do
    elseif (ndir == MY) then
       jf = jf0
       jc = js
       do m = lbound(fc,4), ubound(fc,4)
          do kc=ks,ke
             do ic=is,ie
                kf = IJKF(kc, ks) - ks + kf0
                if = IJKF(ic, is) - is + if0
                fc(ic,jc,kc,m) = &
                     (ff(if,  jf,  kf,m) + ff(if,  jf,  kf+1,m) &
                     +ff(if+1,jf,  kf,m) + ff(if+1,jf,  kf+1,m) )*0.25d0
             end do
          end do
       end do
    elseif (ndir == MZ) then
       kf = kf0
       kc = ks
       do m = lbound(fc,4), ubound(fc,4)
          do jc=js,je
             do ic=is,ie
                jf = IJKF(jc, js) - js + jf0
                if = IJKF(ic, is) - is + if0
                fc(ic,jc,kc,m) = &
                     (ff(if,  jf,  kf,m) + ff(if,  jf+1,  kf,m) &
                     +ff(if+1,jf,  kf,m) + ff(if+1,jf+1,  kf,m) )*0.25d0
             end do
          end do
       end do
    endif
  end subroutine restrict

end module fmg_reflux
