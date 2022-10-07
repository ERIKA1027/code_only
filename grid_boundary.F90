#include "config.h"
#include "packarr.h"

!activating debug output
!#define KS_DEBUG

!-------------------------------------------------------------------------
! Module for grid boundary condition
!
! 粗いグリッドから順に境界条件をfixしてゆく。
! 既に正しい境界条件があるレベルは、パスする。
! 2008/02/29 計算内部と袖のステップ数(U*_StepNumber*)が異なるときだけ、袖の値を更新する。
!-------------------------------------------------------------------------
module grid_boundary
#ifndef KS_DEBUG
  use grid , only : Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh, Mmin, Mmax, Lmin, Lmax, Ngh, Gidmin, Gidmax, Left, Right, Imin, Imax, Jmin, Jmax, Kmin, Kmax, COMPLETE, PREDICTOR, CORRECTOR, &
       U_StepNumber, U_StepNumberGhostCell, U1_StepNumber, U1_StepNumberGhostCell, U2_StepNumber, U2_StepNumberGhostCell
#else !KS DEBUG
  use grid
#endif

  implicit none
  private
  integer,parameter :: Ev = 0, Od = 1     ! even and odd
  integer,parameter :: Send = 0, Recv = 1 ! send reciev code
  integer,save :: STEP_MODE          ! predictor or corrector, complete?
  integer,save :: CurrentLevel       ! 現在のレベル
  integer,save :: LevelUpto          ! このレベルまで
  ! 境界が更新されているかどうかのフラッグ
  logical,save,dimension(Lmin:Lmax) :: Bool_fixed
  ! 境界を更新するかどうか
  logical,save :: Bool_fix_current, Bool_fix_1order, Bool_fix_2order
  integer,save,dimension(Left:Right,MX:MZ,Send:Recv) :: Ins, Ine, Jns, Jne, Kns, Kne
  integer,save,dimension(Left:Right,MX:MZ,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce
  ! 袖のバッファ FLS(ndir)%Bufc(i,j,k), FLP(ndir)%Bufc(i,j,k)

#ifdef KS_DEBUG
  integer,save :: mygid
#endif !KS DEBUG
  
  type t_buf
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(:),allocatable :: FLS
  type(t_buf),save,dimension(:),allocatable :: FLP
  public :: boundary_grid
contains
  !-----------------------------------------------------------------------
  ! set grid boundary condition of grid (level, id)
  !
  ! 境界条件は以下の３つの操作で必要である。
  !
  ! (1) 格子のリファインメント
  ! (2) 流体時間積分 predictor-step
  ! (3) 流体時間積分 corrector-step
  !
  ! 制限事項：
  !   親グリッドに内包されていること。斜めに祖グリッドと隣接する場合は不可。
  !
  ! 注意：
  !   グリッドを細分化するとき、マージン確保のために粗い格子が生成されることがある。
  !   この場合, 生成された粗い格子袖の値は未確定のままであるが、直後に粗い格子を
  !   時間ステップを推進するので、安全に正しい値に更新される。
  !-----------------------------------------------------------------------
  subroutine boundary_grid( level, mode )
    integer,intent(IN) :: level, mode ! 境界条件のレベル, モード
    integer :: ndir
    if ( level < Lmin ) return
    STEP_MODE = mode
    LevelUpto = level
    call grid_boundary_init
    call gridsize_init_parent
    call gridsize_init_samelev

    ! ------------------------------
    ! fix ghostcell on parent level
    ! ------------------------------
    Bool_fix_current = .true.
    Bool_fix_1order = .true.
    Bool_fix_2order = .true.
    CurrentLevel = LevelUpto-1
    do ndir = MX, MZ
       call fix_from_samelev(ndir)
    enddo

    ! ------------------------------
    ! fix ghostcell on this level
    ! ------------------------------
    Bool_fix_current = .true.
    Bool_fix_1order = .false.
    Bool_fix_2order = .false.
    CurrentLevel = LevelUpto
    do ndir = MX, MZ
       call fix_from_parent(ndir)
       call fix_from_samelev(ndir)
    enddo

    ! --------------------------------
    !  update StepNumber in ghost cell
    ! --------------------------------
    if (Bool_fix_current) U_StepNumberGhostCell(level) = U_StepNumber(level)
    if (Bool_fix_1order)  U1_StepNumberGhostCell(level) = U1_StepNumber(level)
    if (Bool_fix_2order)  U2_StepNumberGhostCell(level) = U2_StepNumber(level)

  end subroutine boundary_grid
  !-----------------------------------------------------------------------
  ! grid_boundary を初期化
  !-----------------------------------------------------------------------
  subroutine grid_boundary_init
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary' )
    Bool_fixed(:) = .false.
    bool_initialized = .true.
  end subroutine grid_boundary_init
  !-----------------------------------------------------------------------
  ! Initialize gridsize for transfer (same level)
  !-----------------------------------------------------------------------
  subroutine gridsize_init_samelev
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    integer :: lr, ndir
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary: samelev' )
    ! Ins(LR surf, direction)
    !     LR,      MX,MY,MZ
    Ins(:,:,:) = Imingh
    Ine(:,:,:) = Imaxgh

    Jns(:,:,:) = Jmingh
    Jne(:,:,:) = Jmaxgh

    Kns(:,:,:) = Kmingh
    Kne(:,:,:) = Kmaxgh

    Ins(Left ,MX,Send) = Imax-Ngh+1  !x方向左境界を送るセル
    Ine(Left ,MX,Send) = Imax
    Ins(Right,MX,Send) = Imin
    Ine(Right,MX,Send) = Imin+Ngh-1

    Jns(Left ,MY,Send) = Jmax-Ngh+1
    Jne(Left ,MY,Send) = Jmax
    Jns(Right,MY,Send) = Jmin
    Jne(Right,MY,Send) = Jmin+Ngh-1

    Kns(Left ,MZ,Send) = Kmax-Ngh+1
    Kne(Left ,MZ,Send) = Kmax
    Kns(Right,MZ,Send) = Kmin
    Kne(Right,MZ,Send) = Kmin+Ngh-1

    Ins(Left ,MX,Recv) = Imin-Ngh !x方向左境界を受けるセル
    Ine(Left ,MX,Recv) = Imin-1
    Ins(Right,MX,Recv) = Imax+1
    Ine(Right,MX,Recv) = Imax+Ngh

    Jns(Left ,MY,Recv) = Jmin-Ngh
    Jne(Left ,MY,Recv) = Jmin-1
    Jns(Right,MY,Recv) = Jmax+1
    Jne(Right,MY,Recv) = Jmax+Ngh

    Kns(Left ,MZ,Recv) = Kmin-Ngh
    Kne(Left ,MZ,Recv) = Kmin-1
    Kns(Right,MZ,Recv) = Kmax+1
    Kne(Right,MZ,Recv) = Kmax+Ngh

#define SZ(SENDRECV) Ins(lr,ndir,SENDRECV):Ine(lr,ndir,SENDRECV),Jns(lr,ndir,SENDRECV):Jne(lr,ndir,SENDRECV),Kns(lr,ndir,SENDRECV):Kne(lr,ndir,SENDRECV),Mmin:Mmax
    ! prepare FLS(ndir)%Bufc
    if ( .not. allocated(FLS) ) allocate( FLS(MX:MZ) )
    if ( .not. associated(FLS(MX)%Bufc) ) then
       lr = Left
       do ndir = MX, MZ
          allocate( FLS(ndir)%Bufc( SZ(Send) ) )
       end do
    endif
#undef SZ
    bool_initialized = .true.
  end subroutine gridsize_init_samelev
  !-----------------------------------------------------------------------
  ! Initialize gridsize for transfer (parent level)
  !-----------------------------------------------------------------------
  subroutine gridsize_init_parent
    use io_util, only : print_msg
    logical,save :: bool_initialized = .false.
    integer :: lr, ndir, ieo, jeo, keo
    if ( bool_initialized ) return
    call print_msg( 'initialize grid-boundary: parent' )
    ! Ics(LR surf, direction, parity of this grid)
    !     LR,      MX,MY,MZ,  Ev,Od
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

    Ics(Left ,MX,:) = Imax-Ngh/2
    Ice(Left ,MX,:) = Imax+Ngh/2
    Ics(Right,MX,:) = Imin-Ngh/2
    Ice(Right,MX,:) = Imin+Ngh/2

    Jcs(Left ,MY,:) = Jmax-Ngh/2
    Jce(Left ,MY,:) = Jmax+Ngh/2
    Jcs(Right,MY,:) = Jmin-Ngh/2
    Jce(Right,MY,:) = Jmin+Ngh/2

    Kcs(Left ,MZ,:) = Kmax-Ngh/2
    Kce(Left ,MZ,:) = Kmax+Ngh/2
    Kcs(Right,MZ,:) = Kmin-Ngh/2
    Kce(Right,MZ,:) = Kmin+Ngh/2

#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmax
    ! prepare FLP(ndir)%Bufc
    if ( .not. allocated(FLP) ) allocate( FLP(MX:MZ) )
    if ( .not. associated(FLP(MX)%Bufc) ) then
       lr = Left
       ieo = Ev
       jeo = Ev
       keo = Ev
       do ndir = MX, MZ
          allocate( FLP(ndir)%Bufc( SZ ) )
       end do
    endif
#undef SZ
    bool_initialized = .true.
  end subroutine gridsize_init_parent
  !-----------------------------------------------------------------------
  ! 同じレベルから転送
  !-----------------------------------------------------------------------
  subroutine fix_from_samelev(ndir)
    use mpilib
    use packarr
    use grid
    integer,intent(IN) :: ndir
    integer :: n, lr, rank, rankd, ranks, gid, gidd, gids
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: us, ud, buf
#define SZ(SENDRECV) Ins(lr,ndir,SENDRECV):Ine(lr,ndir,SENDRECV),Jns(lr,ndir,SENDRECV):Jne(lr,ndir,SENDRECV),Kns(lr,ndir,SENDRECV):Kne(lr,ndir,SENDRECV),Mmin:Mmax

    if ( CurrentLevel < Lmin ) return

    myrank = get_myrank()
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr) ! 必要。方向別に転送するので

    ! 初期化
    buf => FLS(ndir)%Bufc
    call pkar_reset

    ! 送信準備と同じノードは代入していまう。
    do rank = 0, NPE-1          ! 受信RANK
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd  = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             gids  = NeighborGid(lr, ndir, gidd, rankd)
             ranks = NeighborRank(lr, ndir, gidd, rankd)
             if ( gids == Undefi ) cycle ! 隣兄弟がいないと無視

             if ( myrank == rankd .and. myrank == ranks ) then !受信と送信が同じnode
                if (Bool_fix_current ) then
                   us => get_Up(gids)
                   ud => get_Up(gidd)
                   ud(SZ(Recv)) = us(SZ(Send))
                endif
                if (Bool_fix_1order) then
                   us => get_U1orderp(gids)
                   ud => get_U1orderp(gidd)
                   ud(SZ(Recv)) = us(SZ(Send))
                endif
                if (Bool_fix_2order) then
                   us => get_U2orderp(gids)
                   ud => get_U2orderp(gidd)
                   ud(SZ(Recv)) = us(SZ(Send))
                endif
             else if ( myrank == ranks .or. myrank == rankd ) then    !送信バッファに詰め込む
                if (Bool_fix_current) then
                   if (myrank == ranks) then
                      us => get_Up(gids)
                      buf = us(SZ(Send))
                   endif
                   PACK_SEND4(buf, myrank, ranks, rankd)
                endif
                if (Bool_fix_1order) then
                   if (myrank == ranks) then
                      us => get_U1orderp(gids)
                      buf = us(SZ(Send))
                   endif
                   PACK_SEND4(buf, myrank, ranks, rankd)
                endif
                if (Bool_fix_2order) then
                   if (myrank == ranks) then
                      us => get_U2orderp(gids)
                      buf = us(SZ(Send))
                   endif
                   PACK_SEND4(buf, myrank, ranks, rankd)
                endif
             endif
          end do
       enddo
    enddo

    call pkar_sendrecv()

    ! pop
    do rank = 0, NPE -1         ! 受信RANK
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd  = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             gids  = NeighborGid(lr, ndir, gidd, rankd)
             ranks = NeighborRank(lr, ndir, gidd, rankd)
             if ( gids == Undefi ) cycle ! 隣兄弟がいないと無視
             if ( myrank == ranks ) cycle
             if ( myrank == rankd ) then !受信バッファから荷解き
                if (Bool_fix_current) then
                   ud => get_Up(gidd)
                   UNPACK_RECV4(buf, myrank, ranks, rankd)
                   ud(SZ(Recv)) = buf
                endif
                if (Bool_fix_1order) then
                   ud => get_U1orderp(gidd)
                   UNPACK_RECV4(buf, myrank, ranks, rankd)
                   ud(SZ(Recv)) = buf
                endif
                if (Bool_fix_2order) then
                   ud => get_U2orderp(gidd)
                   UNPACK_RECV4(buf, myrank, ranks, rankd)
                   ud(SZ(Recv)) = buf
                endif
             endif
          enddo
       enddo
    enddo
#undef SZ
  end subroutine fix_from_samelev
  !-----------------------------------------------------------------------
  ! 親から転送
  ! 時間的に補間 → 転送 → 空間補間
  !-----------------------------------------------------------------------
  subroutine fix_from_parent(ndir)
    use mpilib
    use packarr
    use grid
    integer,intent(IN) :: ndir
    integer :: levelc, n, pgid, prank
    integer :: gid, gids, gidd, lr, rank, rankd, ranks
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ud, udbuf
    if ( CurrentLevel <= Lmin ) return
!!$    call mpi_barrier(MPI_COMM_WORLD, ierr) ! 必要。方向別に転送するので

    ! 初期化
    call pkar_reset
    udbuf => FLP(ndir)%Bufc
    levelc = CurrentLevel -1
    myrank = get_myrank()
    do rank = 0, NPE -1
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd  = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             if ( NeighborGid(lr, ndir, gidd, rankd) /= Undefi ) cycle ! 隣がいるので、除外
             pgid = ParentGid(gidd,rankd)
             prank = ParentRank(gidd,rankd)
             gids = NeighborGid(lr, ndir, pgid, prank)
             ranks = NeighborRank(lr, ndir, pgid, prank)
             if ( gids == Undefi ) cycle ! 親に隣がいないこともある（計算領域の境界など）
             if (myrank == ranks) call interp_time(udbuf, lr, ndir, gids, gidd, rankd)
             PACK_SEND4(udbuf, myrank, ranks, rankd)
          enddo
       enddo
    enddo

    call pkar_sendrecv()

    ! 受信バッファから荷ほどき
    do rank = 0, NPE -1
       do n = Gidmin, GidListNodeMax(CurrentLevel, rank)
          do lr = Left, Right
             gidd = GidListNode(n, CurrentLevel, rank)
             rankd = rank
             if ( NeighborGid(lr, ndir, gidd, rankd) /= Undefi ) cycle ! 隣がいるので、除外
             pgid = ParentGid(gidd,rankd)
             prank = ParentRank(gidd,rankd)
             gids = NeighborGid(lr, ndir, pgid, prank)
             ranks = NeighborRank(lr, ndir, pgid, prank)
             if ( gids == Undefi ) cycle ! 親に隣がいないこともある（計算領域の境界など）
             if (myrank == rankd) then
                UNPACK_RECV4( udbuf, myrank, ranks, rankd )
                ud => get_Up(gidd)
#ifdef KS_DEBUG
                mygid = gidd
#endif !KS DEBUG
                call interp_space(ud, udbuf, ndir, lr)
             end if
          enddo
       enddo
    enddo
  end subroutine fix_from_parent
  !-----------------------------------------------------------------------
  ! 空間方向に補間し、袖に値を入れる
  !-----------------------------------------------------------------------
  subroutine interp_space( ud, buf, ndir, lr)
    use eos
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: ud, buf
    integer,intent(IN) :: ndir, lr
    real(kind=DBL_KIND),dimension(ARRAYSIZE3(buf),MX:MZ) :: grad
    real(kind=DBL_KIND) :: hf, hc, di, dj, dk, dv, dxl, dxr, dyl, dyr, dzl, dzr
    integer :: i,j,k,m, ic,jc,kc,is,ie,js,je,ks,ke,ms,me,ibs,ibe,jbs,jbe,kbs,kbe,mbs,mbe
    integer :: ig0, jg0, kg0                ! origin of grid coordinates
    real(kind=DBL_KIND) :: a, b, minmod
    minmod(a,b) = sign(1.d0,a)*max(0.d0,min(abs(a),sign(1.d0,a)*b))

    hf = 1.d0                   ! fine cell width
    hc = 2* hf                  ! coarse cell width
    dv = 1.d0                   ! volume for dummy

    is = Imingh
    ie = Imaxgh
    js = Jmingh
    je = Jmaxgh
    ks = Kmingh
    ke = Kmaxgh
    if (ndir == MX .and. lr == Left ) ie = Imin-1
    if (ndir == MX .and. lr == Right) is = Imax+1
    if (ndir == MY .and. lr == Left ) je = Jmin-1
    if (ndir == MY .and. lr == Right) js = Jmax+1
    if (ndir == MZ .and. lr == Left ) ke = Kmin-1
    if (ndir == MZ .and. lr == Right) ks = Kmax+1
    ibs = lbound(buf,1) ; ibe = ubound(buf,1)
    jbs = lbound(buf,2) ; jbe = ubound(buf,2)
    kbs = lbound(buf,3) ; kbe = ubound(buf,3)
    mbs = lbound(buf,4) ; mbe = ubound(buf,4)

!!$    call conv_u2w( buf, dv )

    ! -----------
    ! 空間的に補間
    ! -----------
    do m = mbs, mbe
       ! ------------
       ! 勾配を求める
       ! ------------
       do k = kbs+1, kbe-1
          do j = jbs+1, jbe-1
             do i = ibs+1, ibe-1
                dxr = (buf(i+1,j,k,m)-buf(i,j,k,m))/hc
                dxl = (buf(i,j,k,m)-buf(i-1,j,k,m))/hc
                dyr = (buf(i,j+1,k,m)-buf(i,j,k,m))/hc
                dyl = (buf(i,j,k,m)-buf(i,j-1,k,m))/hc
                dzr = (buf(i,j,k+1,m)-buf(i,j,k,m))/hc
                dzl = (buf(i,j,k,m)-buf(i,j,k-1,m))/hc

                !------ WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING -------!
                !                                                                                  !
                !   commented out by KS to avoid a lot of call for rescue around boundaries        !
                !   between coarse and fine grids                                                  !
                !                                                                                  !
                !    -> always use minmod                                                          !
                !                                                                                  !
                !----------------------------------------------------------------------------------!
                ! if (ndir == MX .and. lr == Left) then
                !    dxr = dxl
                ! elseif (ndir == MX .and. lr == Right) then
                !    dxl = dxr
                ! elseif (ndir == MY .and. lr == Left) then
                !    dyr = dyl
                ! elseif (ndir == MY .and. lr == Right) then
                !    dyl = dyr
                ! elseif (ndir == MZ .and. lr == Left) then
                !    dzr = dzl
                ! elseif (ndir == MZ .and. lr == Right) then
                !    dzl = dzr
                ! end if
                !------ WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING -------!                

                grad(i,j,k,MX) = minmod(dxr, dxl)
                grad(i,j,k,MY) = minmod(dyr, dyl)
                grad(i,j,k,MZ) = minmod(dzr, dzl)

!!$                grad(i,j,k,MX) = minmod(buf(i+1,j,k,m)-buf(i,j,k,m), buf(i,j,k,m)-buf(i-1,j,k,m))/hc
!!$                grad(i,j,k,MY) = minmod(buf(i,j+1,k,m)-buf(i,j,k,m), buf(i,j,k,m)-buf(i,j-1,k,m))/hc
!!$                grad(i,j,k,MZ) = minmod(buf(i,j,k+1,m)-buf(i,j,k,m), buf(i,j,k,m)-buf(i,j,k-1,m))/hc
!!$                ! for debug
!!$                grad(i,j,k,MX) = (buf(i+1,j,k,m,0)-buf(i-1,j,k,m,0))/hc/2
!!$                grad(i,j,k,MY) = (buf(i,j+1,k,m,0)-buf(i,j-1,k,m,0))/hc/2
!!$                grad(i,j,k,MZ) = (buf(i,j,k+1,m,0)-buf(i,j,k-1,m,0))/hc/2
             enddo
          enddo
       enddo
       ! This code assumes that odd and even indexes are left and right cells
!OCL NOVREC
       do k = ks, ke
          do j = js, je
             do i = is, ie
                ! indexes of coarse cell
                ic = IJKC( i, is ) - is + ibs+1
                jc = IJKC( j, js ) - js + jbs+1
                kc = IJKC( k, ks ) - ks + kbs+1
                ! deviation of fine grid from the coase cell
                di = ( modulo(i,2) - 0.5d0 )*hf
                dj = ( modulo(j,2) - 0.5d0 )*hf
                dk = ( modulo(k,2) - 0.5d0 )*hf
                ud(i,j,k,m) = buf(ic,jc,kc,m)+grad(ic,jc,kc,MX)*di+grad(ic,jc,kc,MY)*dj+grad(ic,jc,kc,MZ)*dk

#ifdef KS_DEBUG
             if (globdbg_myrank==globdbg_rank .and. mygid==globdbg_gid .and. &
                  i==globdbg_i+lbound(ud,1) .and. j==globdbg_j+lbound(ud,2) .and.k==globdbg_k+lbound(ud,3)) then
                if (m == MRHO .or. m == MP) then
                   print '(A,8I8,/,1P5E15.7,/,1P7E15.7,/,1P7E15.7,/,1P7E15.7)', "(grd_bnd, KS DEBUG)",&
                        globdbg_myrank, mygid, i-lbound(ud,1), j-lbound(ud,2), k-lbound(ud,3),ndir,lr,m,&
                        ud(i,j,k,m), buf(ic,jc,kc,m), grad(ic,jc,kc,MX)*di, grad(ic,jc,kc,MY)*dj, grad(ic,jc,kc,MZ)*dk,&
                        grad(ic,jc,kc,MX)*di, (buf(ic,jc,kc,m)-buf(ic-1,jc,kc,m))/hc, (buf(ic+1,jc,kc,m)-buf(ic,jc,kc,m))/hc, di, buf(ic-1:ic+1,jc,kc,m),&
                        grad(ic,jc,kc,MY)*dj, (buf(ic,jc,kc,m)-buf(ic,jc-1,kc,m))/hc, (buf(ic,jc+1,kc,m)-buf(ic,jc,kc,m))/hc, dj, buf(ic,jc-1:jc+1,kc,m),&
                        grad(ic,jc,kc,MZ)*dk, (buf(ic,jc,kc,m)-buf(ic,jc,kc-1,m))/hc, (buf(ic,jc,kc+1,m)-buf(ic,jc,kc,m))/hc, dk, buf(ic,jc,kc-1:kc+1,m)
                end if
             end if
#endif !KS DEBUG                
             enddo
          end do
       end do
    end do
!!$    call conv_w2u( ud(is:ie,js:je,ks:ke,:), dv )
  end subroutine interp_space
  !------------------------------------------------------------------
  ! 袖を時間方向に補間
  !------------------------------------------------------------------
#ifndef SINGLE_STEP
#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmax
  subroutine interp_time(udbuf, lr, ndir, gids, gidd, rankd)
    use mpilib
    use grid
    use eos
    integer,intent(IN) :: lr, ndir, gids, gidd, rankd
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, u1, u2, udbuf, b0, b1, b2
    real(kind=DBL_KIND) :: rdstp, stp, dstpl, dstpr, dv
    integer :: lstep, lstepc, ldstep, ldstepc, ieo, jeo, keo
    type(t_buf),dimension(:),allocatable,save :: ub0, ub1, ub2

    call interp_time_init
    call left_or_right(gidd, rankd, ieo, jeo, keo)
    ! -----------------------------
    ! 時間的に補間
    ! u(:,:,:,:) .. 未来2次精度(親の現在)
    ! u1(:,:,:,:) .. 未来1次精度
    ! u2(:,:,:,:) .. 過去2次精度
    ! ----------------------------

    lstep = Step( CurrentLevel )
    lstepc = Step( CurrentLevel - 1 )
    ldstep = Dstep( CurrentLevel )
    ldstepc = Dstep( CurrentLevel - 1 )

    if (CurrentLevel == LevelUpto) then
       select case ( STEP_MODE )   !1stepのうち、どのくらい進んでいるか？
       case ( PREDICTOR )
          rdstp = 0.d0
       case ( CORRECTOR )
          rdstp = 0.5d0
       case default
          rdstp = 1.d0
       end select
    else
       rdstp = 1.d0
    endif

    stp = lstep - ( 1.d0 - rdstp ) * ldstep ! 注目するステップ数
    dstpl = stp - ( lstepc - ldstepc )      ! ステップ幅（左）
    dstpr = lstepc - stp                    ! ステップ幅（右）
    ! 2次精度補間
    u  => get_Up(gids)
    u1 => get_U1orderp(gids)
    u2 => get_U2orderp(gids)
    b0 => ub0(ndir)%Bufc
    b1 => ub1(ndir)%Bufc
    b2 => ub2(ndir)%Bufc
    b0 = u(SZ)
    b1 = u1(SZ)
    b2 = u2(SZ)
    dv = 1.d0                   ! volume for dummy
    call conv_u2w( b0, dv )
    call conv_u2w( b1, dv )
    call conv_u2w( b2, dv )
    b1 = 2*b1 - b2              ! 1次精度＠半ステップを全ステップまで進める。
    udbuf = (b0 - b1)*(dstpl/ldstepc)**2 + (b2*dstpr + b1*dstpl)/ldstepc
    call conv_w2u(udbuf, dv)
  contains
  !------------------------------------------------------------------
#undef SZ
#define SZ Ics(l,n,i):Ice(l,n,i),Jcs(l,n,j):Jce(l,n,j),Kcs(l,n,k):Kce(l,n,k),Mmin:Mmax
    subroutine interp_time_init
      logical,save :: bool_init = .false.
      integer :: l, i, j, k, n
      if ( bool_init ) return
      bool_init = .true.
      l = Left
      i = Ev
      j = Ev
      k = Ev
      allocate( ub0(MX:MZ) )
      allocate( ub1(MX:MZ) )
      allocate( ub2(MX:MZ) )
      do n = MX, MZ
         allocate( ub0(n)%Bufc( SZ ) )
         allocate( ub1(n)%Bufc( SZ ) )
         allocate( ub2(n)%Bufc( SZ ) )
      end do
    end subroutine interp_time_init
#undef SZ
  end subroutine interp_time
#else     ! SINGLE_STEP
  !------------------------------------------------------------------
  ! interp in time for SINGLE_STEP mode
  !------------------------------------------------------------------
#define SZ Ics(lr,ndir,ieo):Ice(lr,ndir,ieo),Jcs(lr,ndir,jeo):Jce(lr,ndir,jeo),Kcs(lr,ndir,keo):Kce(lr,ndir,keo),Mmin:Mmax
  subroutine interp_time(udbuf, lr, ndir, gids, gidd, rankd)
    use mpilib
    use grid
    use eos
    integer,intent(IN) :: lr, ndir, gids, gidd, rankd
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, udbuf
    integer :: ieo, jeo, keo
    call left_or_right(gidd, rankd, ieo, jeo, keo)
    select case ( STEP_MODE )
    case ( PREDICTOR )
       u  => get_U2orderp(gids)
    case ( CORRECTOR )
       u  => get_U1orderp(gids)
    case default
       u  => get_Up(gids)
    end select
    udbuf = u(SZ)
  end subroutine interp_time
#undef SZ
#endif    ! SINGLE_STEP
end module grid_boundary
