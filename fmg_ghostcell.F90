!
! FMG にて ghostcell(袖)に値を代入する。
! 親グリッドから子グリッドの袖へ線形補間(bilinear interpolation)を行う。
!
! fix level 0
! fix level 0-1x, 1x, 0-1y, ly, 0-1z, 1z
! fix level 1-2x, 2x, 1-2y, 2y, 1-2z, 2z
! fix level 2-3x, 3x, 2-3y, 3y, 2-3z, 3z
#include "config.h"
#include "packarr.h"
! #define VERBOSE
!-------------------------------------------------------------------------
! Module for ghost cell including inter-node-communications
!-------------------------------------------------------------------------
module fmg_ghostcell
  use fmg_data
  use mpilib
  implicit none
  private
  ! FLSV(fmglev)%Bufc(i,j,k), FLPV(fmglev)%Bufc(i,j,k), FLSS(fmglev)%Bufc(i,j,k), FLPS(fmglev)%Bufc(i,j,k)
  type t_buf
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(:,:),allocatable :: FLSV ! for vector
  type(t_buf),save,dimension(:,:),allocatable :: FLPV ! for vector
  type(t_buf),save,dimension(:,:),allocatable :: FLSS ! for scalar
  type(t_buf),save,dimension(:,:),allocatable :: FLPS ! for scalar
  type(t_buf),save,dimension(:,:,:),allocatable :: BLV, BLS ! buffer for bilinear interpolation
  type(t_buf),save,dimension(:,:,:),allocatable :: CIXV, CIXYV, CIXYZV ! for bicubic interpolation
  type(t_buf),save,dimension(:,:,:),allocatable :: CIXS, CIXYS, CIXYZS

  integer,parameter :: Send = 0, Recv = 1 ! send reciev code
  integer,parameter :: Ev = 0, Od = 1     ! parity of cell index
  integer,save,dimension(Left:Right,MX:MZ,Send:Recv) :: Ins, Ine, Jns, Jne, Kns, Kne
  integer,save,dimension(Left:Right,MX:MZ,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  integer,save :: NghCg   ! buffer size in coarse grid in normal direction
  integer,parameter :: NghFg = 1 ! buffer size in fine grid
  public :: fmg_ghostcell_fix, fmg_ghfix_samelev, fmg_ghfix_parentlev, fmg_ghfix_samelev_init, fmg_ghfix_parentlev_init, fmg_ghostcell_finalize
contains
  !-------------------------------------------------------------------------
  ! fix the boundary value of grids given by fmglev
  !-------------------------------------------------------------------------
  subroutine fmg_ghostcell_fix(fmglev,icode, cubic, tricubic)
    integer,intent(IN) :: fmglev, icode
    logical,intent(IN),optional :: cubic, tricubic
    integer :: ndir, amrlev
    call fmg_ghfix_samelev_init(fmglev)
    call fmg_ghfix_parentlev_init(fmglev)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do ndir = MX, MZ
          call fmg_ghfix_parentlev(amrlev, fmglev, ndir, icode, cubic, tricubic)
          call fmg_ghfix_samelev(amrlev, fmglev,ndir,icode)
       end do
    enddo
  end subroutine fmg_ghostcell_fix
  !-------------------------------------------------------------------------
  ! Initialize fmg_ghostcell
  !-------------------------------------------------------------------------
  subroutine fmg_ghfix_samelev_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir, lr

    Ins(:,:,:) = fmg_get_imingh(fmglev)
    Ine(:,:,:) = fmg_get_imaxgh(fmglev)

    Jns(:,:,:) = fmg_get_jmingh(fmglev)
    Jne(:,:,:) = fmg_get_jmaxgh(fmglev)

    Kns(:,:,:) = fmg_get_kmingh(fmglev)
    Kne(:,:,:) = fmg_get_kmaxgh(fmglev)

    Ins(Left ,MX,Send) = GridSize(fmglev)%Imax-Ngh+1 !x方向左境界を送るセル
    Ine(Left ,MX,Send) = GridSize(fmglev)%Imax
    Ins(Right,MX,Send) = GridSize(fmglev)%Imin
    Ine(Right,MX,Send) = GridSize(fmglev)%Imin+Ngh-1

    Jns(Left ,MY,Send) = GridSize(fmglev)%Jmax-Ngh+1
    Jne(Left ,MY,Send) = GridSize(fmglev)%Jmax
    Jns(Right,MY,Send) = GridSize(fmglev)%Jmin
    Jne(Right,MY,Send) = GridSize(fmglev)%Jmin+Ngh-1

    Kns(Left ,MZ,Send) = GridSize(fmglev)%Kmax-Ngh+1
    Kne(Left ,MZ,Send) = GridSize(fmglev)%Kmax
    Kns(Right,MZ,Send) = GridSize(fmglev)%Kmin
    Kne(Right,MZ,Send) = GridSize(fmglev)%Kmin+Ngh-1

    Ins(Left ,MX,Recv) = GridSize(fmglev)%Imin-Ngh !x方向左境界を受けるセル
    Ine(Left ,MX,Recv) = GridSize(fmglev)%Imin-1
    Ins(Right,MX,Recv) = GridSize(fmglev)%Imax+1
    Ine(Right,MX,Recv) = GridSize(fmglev)%Imax+Ngh

    Jns(Left ,MY,Recv) = GridSize(fmglev)%Jmin-Ngh
    Jne(Left ,MY,Recv) = GridSize(fmglev)%Jmin-1
    Jns(Right,MY,Recv) = GridSize(fmglev)%Jmax+1
    Jne(Right,MY,Recv) = GridSize(fmglev)%Jmax+Ngh

    Kns(Left ,MZ,Recv) = GridSize(fmglev)%Kmin-Ngh
    Kne(Left ,MZ,Recv) = GridSize(fmglev)%Kmin-1
    Kns(Right,MZ,Recv) = GridSize(fmglev)%Kmax+1
    Kne(Right,MZ,Recv) = GridSize(fmglev)%Kmax+Ngh

#define SZ(SENDRECV) Ins(lr,ndir,SENDRECV):Ine(lr,ndir,SENDRECV),Jns(lr,ndir,SENDRECV):Jne(lr,ndir,SENDRECV),Kns(lr,ndir,SENDRECV):Kne(lr,ndir,SENDRECV),Mmin:Mmax
    ! prepare FLSV(fmglev)%Bufc
    if ( .not. allocated(FLSV) ) allocate( FLSV(MX:MZ,FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLSV(MX,fmglev)%Bufc) ) then
       lr = Left
       do ndir = MX, MZ
          allocate( FLSV(ndir,fmglev)%Bufc( SZ(Send) ) )
       end do
    endif
#undef SZ
#define SZ(SENDRECV) Ins(lr,ndir,SENDRECV):Ine(lr,ndir,SENDRECV),Jns(lr,ndir,SENDRECV):Jne(lr,ndir,SENDRECV),Kns(lr,ndir,SENDRECV):Kne(lr,ndir,SENDRECV),Mmin:Mmin
    ! prepare FLSS(fmglev)%Bufc
    if ( .not. allocated(FLSS) ) allocate( FLSS(MX:MZ,FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLSS(MX,fmglev)%Bufc) ) then
       lr = Left
       do ndir = MX, MZ
          allocate( FLSS(ndir,fmglev)%Bufc( SZ(Send) ) )
       end do
    endif
#undef SZ

  end subroutine fmg_ghfix_samelev_init
  !-------------------------------------------------------------------------
  ! fix ghost cell by neighbor grid at the same level
  !-------------------------------------------------------------------------
  subroutine fmg_ghfix_samelev(amrlev, fmglev,ndir,icode)
    use packarr
    integer,intent(IN) :: amrlev, fmglev, ndir, icode
    integer :: lr, rank, rankd, ranks, gid, gidd, gids
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: us, ud, buf
#define SZ(SENDRECV) Ins(lr,ndir,SENDRECV):Ine(lr,ndir,SENDRECV),Jns(lr,ndir,SENDRECV):Jne(lr,ndir,SENDRECV),Kns(lr,ndir,SENDRECV):Kne(lr,ndir,SENDRECV),:
    myrank = get_myrank()

!!$    call mpi_barrier(MPI_COMM_WORLD, ierr) ! 必要。方向別に転送するので

    ! 初期化
    call pkar_reset
    if ( fmg_isVector(icode) ) then
       buf => FLSV(ndir,fmglev)%Bufc
    else
       buf => FLSS(ndir,fmglev)%Bufc
    endif
    ! 送信準備と同じノードは代入していまう。
    do rank = 0, NPE-1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborSameLevel(lr,ndir) ) cycle !隣がいないのは除外
             gidd  = gid
             rankd = rank
             gids  = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             if ( myrank == rankd .and. myrank == ranks ) then !受信と送信が同じnode
                us => fmg_get_arrp(amrlev, fmglev, gids, icode)
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                ud(SZ(Recv)) = us(SZ(Send))
             else if ( myrank == ranks ) then    !送信バッファに詰め込む
                us => fmg_get_arrp(amrlev, fmglev, gids, icode)
                buf = us(SZ(Send))
             endif
             if ( (myrank == ranks .or. myrank == rankd) .and. rankd /= ranks) then
                PACK_SEND4( buf, myrank, ranks, rankd )
             endif
          end do
       enddo
    enddo

    call pkar_sendrecv()

    ! pop
    do rank = 0, NPE -1         ! 送信 rank
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborSameLevel(lr,ndir) ) cycle !隣がいない
             gidd  = gid
             rankd = rank
             gids  = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             if ( myrank == ranks ) cycle
             if ( myrank == rankd ) then !受信バッファから荷解き
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                UNPACK_RECV4( buf,  myrank, ranks, rankd )
                ud(SZ(Recv)) = buf
             endif
          enddo
       enddo
    enddo
    ! check
!!$    do rank = 0, NPE-1
!!$       if (.not. associated(bufd(rank)%pkg) ) cycle
!!$       if (ubound(bufd(rank)%pkg,1) /= bufd(rank)%position-1 ) &
!!$            print *,'*** Error. A buffer has unsed region.'
!!$    end do
!!$    ! deallocate arrays
!!$    do rank = 0, NPE -1
!!$       if ( associated( bufs(rank)%pkg ) ) deallocate( bufs(rank)%pkg )
!!$       if ( associated( bufd(rank)%pkg ) ) deallocate( bufd(rank)%pkg )
!!$    enddo
  end subroutine fmg_ghfix_samelev
  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  subroutine fmg_ghfix_parentlev_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: ndir, lr, ieo, jeo, keo
    integer :: ifs, ife, jfs, jfe, kfs, kfe, ibs, ibe, jbs, jbe, kbs, kbe
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: buf

    NghCg = Ngh
    if (Ngh > 2) then
       print *, '*** fmg_ghfix_parentlev_init: error large Ngh is not support. Ngh =', Ngh
       stop
    end if

    ! ics(LR surf, direction, parity of this grid)
    !     LR,      MX,MY,MZ,  Ev,Od
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

    Ics(Left, MX,:) = GridSize(fmglev)%Imax-NghCg+1
    Ice(Left, MX,:) = GridSize(fmglev)%Imax
    Ics(Right,MX,:) = GridSize(fmglev)%Imin
    Ice(Right,MX,:) = GridSize(fmglev)%Imin+NghCg-1

    Jcs(Left, MY,:) = GridSize(fmglev)%Jmax-NghCg+1
    Jce(Left, MY,:) = GridSize(fmglev)%Jmax
    Jcs(Right,MY,:) = GridSize(fmglev)%Jmin
    Jce(Right,MY,:) = GridSize(fmglev)%Jmin+NghCg-1

    Kcs(Left, MZ,:) = GridSize(fmglev)%Kmax-NghCg+1
    Kce(Left, MZ,:) = GridSize(fmglev)%Kmax
    Kcs(Right,MZ,:) = GridSize(fmglev)%Kmin
    Kce(Right,MZ,:) = GridSize(fmglev)%Kmin+NghCg-1

#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmax
    ! prepare FLPV(fmglev)%Bufc
    if ( .not. allocated(FLPV) ) allocate( FLPV(MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLPV(MX, fmglev)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = MX, MZ
          allocate( FLPV(ndir, fmglev)%Bufc( SZ ) )
       end do
    endif
#undef SZ
#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmin
    ! prepare FLPS(fmglev)%Bufc
    if ( .not. allocated(FLPS) ) allocate( FLPS(MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLPS(MX, fmglev)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = MX, MZ
          allocate( FLPS(ndir, fmglev)%Bufc( SZ ) )
       end do
    endif
#undef SZ

    ! buffer for cubic interpolation
    if ( .not. allocated(BLS)    ) allocate( BLS   (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(BLV)    ) allocate( BLV   (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXV)   ) allocate( CIXV  (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYV)  ) allocate( CIXYV (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYZV) ) allocate( CIXYZV(Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXS)   ) allocate( CIXS  (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYS)  ) allocate( CIXYS (Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. allocated(CIXYZS) ) allocate( CIXYZS(Left:Right, MX:MZ, FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(CIXV(Left, MX, fmglev)%Bufc )) then
       do ndir = MX, MZ
          do lr = Left, Right
             buf => FLPV(ndir, fmglev)%Bufc
             ibs = lbound(buf,1)
             ibe = ubound(buf,1)
             jbs = lbound(buf,2)
             jbe = ubound(buf,2)
             kbs = lbound(buf,3)
             kbe = ubound(buf,3)
             call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe)
             allocate(   BLS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmin))
             allocate(   BLV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmax))
             allocate(  CIXV(lr, ndir, fmglev)%Bufc(ifs:ife,jbs:jbe,kbs:kbe,Mmin:Mmax))
             allocate( CIXYV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kbs:kbe,Mmin:Mmax))
             allocate(CIXYZV(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmax))
             allocate(  CIXS(lr, ndir, fmglev)%Bufc(ifs:ife,jbs:jbe,kbs:kbe,Mmin:Mmin))
             allocate( CIXYS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kbs:kbe,Mmin:Mmin))
             allocate(CIXYZS(lr, ndir, fmglev)%Bufc(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmin))
          end do
       end do
    end if
  end subroutine fmg_ghfix_parentlev_init
  !-------------------------------------------------------------------------
  ! fix ghost cell by neighbor grid at the parent level
  ! 親の袖にアクセスする可能性あり
  ! cubic == .TRUE. のとき、袖に平行な面で３次補間, 垂直な方向に2次補間を行う。
  ! tricubic == .TRUE. のとき、袖を3次元3次補間を行う。Ngh = 2 のときのみ。
  !-------------------------------------------------------------------------
  subroutine fmg_ghfix_parentlev(amrlev, fmglev, ndir, icode, cubic, tricubic)
    use packarr
    integer,intent(IN) :: amrlev, fmglev, ndir, icode
    logical,intent(IN),optional :: cubic, tricubic
    integer :: gid, gids, gidd, lr, rank, rankd, ranks, ieo, jeo, keo
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: us, ud, udbuf
    logical :: bool_cubic
    ! labeled by (lr, ndir, ieo, jeo, keo)
    ! size is labeled by (fmglev)
#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),:

    if ( amrlev <= AMR_LevelMin ) return

!!$    call mpi_barrier(MPI_COMM_WORLD, ierr) ! 必要。方向別に転送するので

    bool_cubic = .FALSE. ; if (present(cubic)) bool_cubic = cubic
    if (present(tricubic)) then
       if (tricubic) bool_cubic = .TRUE.
    end if

    ! 初期化
    call pkar_reset
    if ( fmg_isVector(icode) ) then
       udbuf => FLPV(ndir, fmglev)%Bufc
    else
       udbuf => FLPS(ndir, fmglev)%Bufc
    endif
    myrank = get_myrank()
    do rank = 0, NPE -1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
             gids  = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             gidd  = gid
             rankd = rank
             if (myrank == ranks) then
                call fmg_left_or_right(gidd, rankd, amrlev, ieo, jeo, keo)
                us => fmg_get_arrp(amrlev-1, fmglev, gids, icode)
                udbuf = us(SZ)
             end if
             PACK_SEND4( udbuf, myrank, ranks, rankd )
          enddo
       enddo
    enddo

    call pkar_sendrecv()

    ! 受信バッファから荷ほどき
    do rank = 0, NPE -1
       do gid = fmg_get_gidmin_rank(amrlev, rank), fmg_get_gidmax_rank(amrlev, rank)
          do lr = Left, Right
             if ( .not. Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) cycle
             gids  = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
             ranks = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
             gidd  = gid
             rankd = rank
             if (myrank == rankd) then
                call fmg_left_or_right(gidd, rankd, amrlev, ieo, jeo, keo)
                ud => fmg_get_arrp(amrlev, fmglev, gidd, icode)
                UNPACK_RECV4( udbuf, myrank, ranks, rankd )
                if (bool_cubic .and. interpCubic_checkSize(udbuf, ndir)) then
                   call interpCubic(ud, udbuf, ndir, lr, fmglev, tricubic)
                else
                   call interp(ud, udbuf, ndir, lr, fmglev)
                end if
             end if
          enddo
       enddo
    enddo
!!$    ! deallocate arrays
!!$    do rank = 0, NPE -1
!!$       if ( associated( bufs(rank)%pkg ) ) deallocate( bufs(rank)%pkg )
!!$       if ( associated( bufd(rank)%pkg ) ) deallocate( bufd(rank)%pkg )
!!$    enddo

  contains
    !-------------------------------------------------------------------------
    ! interpolation of buf to ud
    ! ud ....... (out) 全グリッド
    ! buf ...... (in)  袖つき親のグリッド一部（オーバーラップ部）
    ! ndir ..... (in)  x, y, z 方向インデックス
    ! lr ....... (in)  Left, Right インデックス
    ! fmglev ... (in)  fmg level
    !-------------------------------------------------------------------------
    subroutine interp(ud, buf, ndir, lr, fmglev)
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: buf
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ud
      integer,intent(IN) :: ndir, lr, fmglev
      !
      integer :: ifs, ife, jfs, jfe, kfs, kfe, i, j, k, ic, jc, kc, ic1, jc1, kc1, m
      integer :: ibool, jbool, kbool, ixbool, jxbool, kxbool, inner
      integer :: iba, jba, kba
      real(kind=DBL_KIND) :: uH, u0, uN, uHH
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: uintp

      if ( fmg_isVector(icode) ) then
         uintp => BLV(lr, ndir, fmglev)%Bufc
      else
         uintp => BLS(lr, ndir, fmglev)%Bufc
      endif

      ! 方向別 boolen, 0 or 1
      ibool = 0
      jbool = 0
      kbool = 0
      if (ndir == MX) ibool = 1
      if (ndir == MY) jbool = 1
      if (ndir == MZ) kbool = 1
      ! 排他的方向別 boolen, 0 or 1
      ixbool = 1-ibool
      jxbool = 1-jbool
      kxbool = 1-kbool

      ! offset of inner cell index
      if (lr == Left) inner = 1
      if (lr == Right) inner = -1

      ! offset of buffer
      iba = lbound(buf,1)
      jba = lbound(buf,2)
      kba = lbound(buf,3)

      ! ----------------------
      ! 空間的に補間 (2次補間)
      ! ----------------------
      ! 2次補間すると、数値流速が保存する flux を求めるための袖の値が求まる。
      ! see MULTIGRID U.Trottenberg et al. (2001) p. 371
      ! 斜め方向の袖は、3方向の補間を順に繰り返すと、適正な値が入るハズ。
      call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe) ! 変更する ud の範囲 (袖)
#define NEAR(n, nc)  ( (nc)+2*modulo((n),2)-1 )
      do m = lbound(buf, 4), ubound(buf, 4)
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  ! indexes of coarse cell
                  ! nearest cell
                  ic = (IJKC( i, GridSize(fmglev)%Imin ))*ixbool + (i-ifs+iba)*ibool
                  jc = (IJKC( j, GridSize(fmglev)%Jmin ))*jxbool + (j-jfs+jba)*jbool
                  kc = (IJKC( k, GridSize(fmglev)%Kmin ))*kxbool + (k-kfs+kba)*kbool
                  ! ic == ic1 if ndir == MX
                  ic1 = NEAR( i, ic ) * ixbool + ic * ibool
                  jc1 = NEAR( j, jc ) * jxbool + jc * jbool
                  kc1 = NEAR( k, kc ) * kxbool + kc * kbool

                  ! ghost point (south cell), bi-linear interpolation
                  uintp(i,j,k,m) &
                       = ibool * (9*buf(ic,jc,kc,m) + buf(ic,jc1,kc1,m) + 3*buf(ic,jc,kc1,m) + 3*buf(ic,jc1,kc,m) )/16 &
                       + jbool * (9*buf(ic,jc,kc,m) + buf(ic1,jc,kc1,m) + 3*buf(ic,jc,kc1,m) + 3*buf(ic1,jc,kc,m) )/16 &
                       + kbool * (9*buf(ic,jc,kc,m) + buf(ic1,jc1,kc,m) + 3*buf(ic,jc1,kc,m) + 3*buf(ic1,jc,kc,m) )/16
               end do
            end do
         end do
      end do
      ! 最内側の袖をなめるループ
      if (ndir == MX .and. lr == Left)  ifs = ife
      if (ndir == MX .and. lr == Right) ife = ifs
      if (ndir == MY .and. lr == Left)  jfs = jfe
      if (ndir == MY .and. lr == Right) jfe = jfs
      if (ndir == MZ .and. lr == Left)  kfs = kfe
      if (ndir == MZ .and. lr == Right) kfe = kfs
      do m = lbound(buf,4), ubound(buf,4)
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  ! ghost point (south cell), quadratic interpolation
                  uH = uintp(i,j,k,m)
                  u0 = ud(i+inner*ibool,   j+inner*jbool,   k+inner*kbool,   m) ! central cell
                  uN = ud(i+inner*ibool*2, j+inner*jbool*2, k+inner*kbool*2, m) ! north cell
                  ud(i,j,k, m) = ( 10 * u0 + 8 * uH - 3 * uN ) / 15 ! quadratic interpolation
                  ! uN,U0,  uH,  uHH
                  ! ud(i) ..... quadratic interpolatin (conservative interpolation)
                  ! ud(i-1) ... cubic interpolation for fmg_interp
                  ! weight: 2/15. -3/7, 6/5, 2/21
                  if (Ngh == 2) then
                     uHH = uintp(i-inner*ibool, j-inner*jbool, k-inner*kbool, m)
                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) &
                          = 2.d0/15.0*uN - 3.d0/7.d0*u0 + 6.d0/5.d0*uH + 2.d0/21.d0*uHH
!!$                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) = 2.d0*ud(i,j,k,m) - uH
                  end if
               enddo
            end do
         end do
      end do
    end subroutine interp
#undef NEAR
    !-------------------------------------------------------------------------
    ! interpolation of buf to ud
    ! ud ....... (out) 全グリッド
    ! buf ...... (in)  袖つき親のグリッド一部（オーバーラップ部）
    ! ndir ..... (in)  x, y, z 方向インデックス
    ! lr ....... (in)  Left, Right インデックス
    ! fmglev ... (in)  fmg level
    !-------------------------------------------------------------------------
    subroutine interpCubic(ud, buf, ndir, lr, fmglev, tricubic)
      use interpolation
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: buf
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ud
      integer,intent(IN) :: ndir, lr, fmglev
      logical,intent(IN),optional :: tricubic
      !
      integer :: ifs, ife, jfs, jfe, kfs, kfe, ibs, ibe, jbs, jbe, kbs, kbe, mbs, mbe,&
           i, j, k, m, ic, jc, kc, ic0, jc0, kc0, ic1, jc1, kc1, ic2, jc2, kc2, ic3, jc3, kc3, &
           icP, jcP, kcP, icN, jcN, kcN, nshift
      integer :: ibool, jbool, kbool, inner
      real(kind=DBL_KIND) :: uH, u0, uN, uHH, xc, yc, zc, a0, a1, a2, a3
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: bufx, bufxy, bufxyz
      logical :: bool_tricubic
      bool_tricubic = .FALSE. ; if (present(tricubic)) bool_tricubic = tricubic

      if ( fmg_isVector(icode) ) then
         bufx   => CIXV  (lr, ndir, fmglev)%Bufc
         bufxy  => CIXYV (lr, ndir, fmglev)%Bufc
         bufxyz => CIXYZV(lr, ndir, fmglev)%Bufc
      else
         bufx   => CIXS  (lr, ndir, fmglev)%Bufc
         bufxy  => CIXYS (lr, ndir, fmglev)%Bufc
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

!!$      allocate(  bufx(ifs:ife,jbs:jbe,kbs:kbe,Mmin:Mmax))
!!$      allocate( bufxy(ifs:ife,jfs:jfe,kbs:kbe,Mmin:Mmax))
!!$      allocate(bufxyz(ifs:ife,jfs:jfe,kfs:kfe,Mmin:Mmax))

      ! 方向別 boolen, 0 or 1
      ibool = 0
      jbool = 0
      kbool = 0
      if (ndir == MX) ibool = 1
      if (ndir == MY) jbool = 1
      if (ndir == MZ) kbool = 1

      ! offset of inner cell index
      if (lr == Left)  inner = 1
      if (lr == Right) inner = -1

      ! ----------------------
      ! 空間的に補間 (2次補間)
      ! ----------------------
      ! 2次補間すると、数値流速が保存する flux を求めるための袖の値が求まる。
      ! see MULTIGRID U.Trottenberg et al. (2001) p. 371
      ! 斜め方向の袖は、3方向の補間を順に繰り返すと、適正な値が入るハズ。

      call fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe) ! 変更する ud の範囲 (袖)

#define NEAR(n, nc)  ( (nc)+2*modulo((n),2)-1 )
      ! x方向に補間
      if (ndir == MX) then
         ! coorespoinding parent cell
         bufx = buf
      else
         do m = mbs, mbe
            do kc = kbs, kbe
               do jc = jbs, jbe
                  do i = ifs, ife
                     ! coorespoinding parent cell
                     icP = IJKC( i, GridSize(fmglev)%Imin )
                     ! neighbor
                     icN = NEAR( i, icP )
                     ! ic0 < ic1 < ic2
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

      ! y方向に補間
      if (ndir == MY) then
         bufxy = bufx
      else
         do m = mbs, mbe
            do kc = kbs, kbe
               do j = jfs, jfe
                  ! coorespoinding parent cell
                  jcP = IJKC( j, GridSize(fmglev)%Jmin )
                  ! neighbor
                  jcN = NEAR( j, jcP )
                  ! jc0 < jc1 < jc2 < jc3
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

      ! z方向に補間
      if (ndir == MZ) then
         bufxyz = bufxy
      else
         do m = mbs, mbe
            do k = kfs, kfe
               ! coorespoinding parent cell
               kcP = IJKC( k, GridSize(fmglev)%Kmin )
               ! neighbor
               kcN = NEAR( k, kcP )
               ! kc0 < kc1 < kc2 < kc3
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

      ! 最内側の袖をなめるループ
      if (ndir == MX .and. lr == Left)  ifs = ife
      if (ndir == MX .and. lr == Right) ife = ifs
      if (ndir == MY .and. lr == Left)  jfs = jfe
      if (ndir == MY .and. lr == Right) jfe = jfs
      if (ndir == MZ .and. lr == Left)  kfs = kfe
      if (ndir == MZ .and. lr == Right) kfe = kfs
      do m = mbs, mbe
         do k = kfs, kfe
            do j = jfs, jfe
               do i = ifs, ife
                  ! ghost point (south cell), quadratic interpolation
                  uH = bufxyz(i,j,k,m)
                  u0 = ud(i+inner*ibool,   j+inner*jbool,   k+inner*kbool,   m) ! central cell
                  uN = ud(i+inner*ibool*2, j+inner*jbool*2, k+inner*kbool*2, m) ! north cell
                  ud(i,j,k, m) = ( 10 * u0 + 8 * uH - 3 * uN ) / 15 ! quadratic interpolation

                  ! uN,U0,  uH,  uHH
                  ! ud(i) ..... quadratic interpolatin (conservative interpolation)
                  ! ud(i-1) ... cubic interpolation for fmg_interp
                  ! weight: 2/15. -3/7, 6/5, 2/21
                  if (Ngh == 2) then
                     uHH = bufxyz(i-inner*ibool, j-inner*jbool, k-inner*kbool, m)
                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) &
                          = 2.d0/15.0*uN - 3.d0/7.d0*u0 + 6.d0/5.d0*uH + 2.d0/21.d0*uHH
                     if (bool_tricubic) &
                          ud(i,j,k,m) = -1.d0/9.d0*uN + 10.d0/21.d0*u0 + 2.d0/3.d0*uH -2.d0/63.d0*uHH

!!$                     ud(i,j,k,m) = -1.d0/9.d0*uN + 10.d0/21.d0*u0 + 2.d0/3.d0*uH -2.d0/63.d0*uHH
!!$                     ud(i-inner*ibool, j-inner*jbool, k-inner*kbool, m) = 2.d0*ud(i,j,k,m) - uH
                  end if
               enddo
            end do
         end do
      end do

    end subroutine interpCubic
#undef NEAR
    !-------------------------------------------------------------------------
    ! check size of buffer for cubic interpolation
    !-------------------------------------------------------------------------
    function interpCubic_checkSize(buf, ndir) result(validSize)
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: buf
      integer,intent(IN) :: ndir
      logical :: validSize
      integer,parameter :: MINBLOCKSIZE = 4
      if (ndir == MX) then
         validSize = size(buf,2) >= MINBLOCKSIZE .and. size(buf,3) >= MINBLOCKSIZE
      elseif (ndir == MY) then
         validSize = size(buf,3) >= MINBLOCKSIZE .and. size(buf,1) >= MINBLOCKSIZE
      elseif (ndir == MZ) then
         validSize = size(buf,1) >= MINBLOCKSIZE .and. size(buf,2) >= MINBLOCKSIZE
      else
         print *, '**** error in fmg_ghostcell::interpCubic_checkSize', ndir
         stop
      end if
#ifdef VERBOSE
      if (.not. validSize) print *, '**** fmg_ghostcell: too small buffer size'
#endif !VERBOSE
    end function interpCubic_checkSize
  end subroutine fmg_ghfix_parentlev
  !-------------------------------------------------------------------------
  ! 変更される袖の範囲を得る
  !-------------------------------------------------------------------------
  subroutine fmg_ghfix_parentlev_get_fineindex(ndir, lr, fmglev, ifs, jfs, kfs, ife, jfe, kfe)
    integer,intent(IN) :: ndir, lr, fmglev
    integer,intent(OUT) :: ifs, jfs, kfs, ife, jfe, kfe
    call fmg_get_gridsizeGh(fmglev, ifs, jfs, kfs, ife, jfe, kfe)
    if (ndir == MX .and. lr == Left)  ife = GridSize(fmglev)%Imin-1
    if (ndir == MX .and. lr == Right) ifs = GridSize(fmglev)%Imax+1
    if (ndir == MY .and. lr == Left)  jfe = GridSize(fmglev)%Jmin-1
    if (ndir == MY .and. lr == Right) jfs = GridSize(fmglev)%Jmax+1
    if (ndir == MZ .and. lr == Left)  kfe = GridSize(fmglev)%Kmin-1
    if (ndir == MZ .and. lr == Right) kfs = GridSize(fmglev)%Kmax+1
  end subroutine fmg_ghfix_parentlev_get_fineindex
  !-------------------------------------------------------------------------
  ! finalize status of ghost cells
  !-------------------------------------------------------------------------
  subroutine fmg_ghostcell_finalize
    integer :: lr, ndir, fmglev
#define FREE_BUF(BLIST) \
    if (allocated(BLIST)) then ;\
       do fmglev = lbound(BLIST, 2), ubound(BLIST, 2) ;\
          do ndir = lbound(BLIST, 1), ubound(BLIST, 1) ;\
             if (associated(BLIST(ndir,fmglev)%Bufc)) deallocate(BLIST(ndir,fmglev)%Bufc) ;\
             nullify(BLIST(ndir,fmglev)%Bufc) ;\
          end do ;\
       end do ;\
       deallocate(BLIST) ;\
    endif
    FREE_BUF(FLSV)
    FREE_BUF(FLPV)
    FREE_BUF(FLSS)
    FREE_BUF(FLPS)
#undef FREE_BUF

#define FREE_BUF(BLIST) \
    if (allocated(BLIST)) then ;\
       do fmglev = lbound(BLIST, 3), ubound(BLIST, 3) ;\
          do ndir = lbound(BLIST, 2), ubound(BLIST, 2) ;\
             do lr = lbound(BLIST, 1), ubound(BLIST, 1) ;\
                if (associated(BLIST(lr, ndir,fmglev)%Bufc)) deallocate(BLIST(lr, ndir,fmglev)%Bufc) ;\
                nullify(BLIST(lr, ndir,fmglev)%Bufc) ;\
             end do;\
          end do;\
       end do ;\
       deallocate(BLIST) ;\
    endif
    FREE_BUF(BLV)
    FREE_BUF(BLS)
    FREE_BUF(CIXV)
    FREE_BUF(CIXYV)
    FREE_BUF(CIXYZV)
    FREE_BUF(CIXS)
    FREE_BUF(CIXYS)
    FREE_BUF(CIXYZS)

  end subroutine fmg_ghostcell_finalize
end module fmg_ghostcell
