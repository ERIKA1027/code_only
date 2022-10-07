#include "config.h"
#include "packarr.h"
!-----------------------------------------------------------------------
! refluxing in grid-boundary
!-----------------------------------------------------------------------
module reflux
  use grid
  implicit none
  private
  integer,parameter :: NBOUNDARY =  2*(MZ-MX+1) ! 境界面の個数（３次元では６面）
  integer,parameter :: NIL = 0
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  integer,parameter :: L = Left, R = Right, Ev = Left, Od = Right
#define SZ_FXL Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax
#define SZ_FXR Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax
#define SZ_FYL Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax
#define SZ_FYR Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax
#define SZ_FZL Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax
#define SZ_FZR Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax
  ! buffer for flux in the grid interfaces.
  type t_flux_boundary
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: xl_bnd => null() ! flux at the interfaces
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: xr_bnd => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: yl_bnd => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: yr_bnd => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: zl_bnd => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: zr_bnd => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: xl_sum => null() ! cummulative flux
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: xr_sum => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: yl_sum => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: yr_sum => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: zl_sum => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: zr_sum => null()
  end type t_flux_boundary
  type(t_flux_boundary),save,dimension(Gidmin:Gidmax) :: Fb

  ! array boundaries
  integer,save,dimension(L:R,MX:MZ,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce, Ias, Iae, Jas, Jae, Kas, Kae

  integer,save :: Levelf             ! fine grid level
  integer,save :: Levelc             ! coarse grid level
  ! 子グリッドを間引くためのバッファ FLP(ndir)%Buff(i,j,k,m)
  type t_buf
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: Buff => null()
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: Bufw => null()
  end type t_buf
  type(t_buf),save,dimension(MX:MZ) :: FLP

  public :: fluxcorrection, save_flux, reflux_write, reflux_read
contains

  !-------------------------------------------------------------------------
  ! gid ,surfcode が与えられたとき、flux buffer の pointer を返す。ex) ptr => get_fb_bnd(gid, NIL)
  !-------------------------------------------------------------------------
#define ALLOC_SAFE_GETP(A, SZ) \
       if (.not. associated(A)) then ;\
          allocate(A(SZ)) ;\
       endif ;\
       ptr => A

  function get_fb_bnd(gid,surfcode) result(ptr)
    integer,intent(IN) :: gid, surfcode
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptr
    select case( surfcode )
    case (NIL)
       ALLOC_SAFE_GETP(Fb(gid)%xl_bnd, SZ_FXL)
    case (NIR)
       ALLOC_SAFE_GETP(Fb(gid)%xr_bnd, SZ_FXR)
    case (NJL)
       ALLOC_SAFE_GETP(Fb(gid)%yl_bnd, SZ_FYL)
    case (NJR)
       ALLOC_SAFE_GETP(Fb(gid)%yr_bnd, SZ_FYR)
    case (NKL)
       ALLOC_SAFE_GETP(Fb(gid)%zl_bnd, SZ_FZL)
    case (NKR)
       ALLOC_SAFE_GETP(Fb(gid)%zr_bnd, SZ_FZR)
    end select
  end function get_fb_bnd
  !-------------------------------------------------------------------------
  ! gid ,surfcode が与えられたとき、flux buffer の pointer を返す。ex) ptr => get_fb_sum(gid, NIL)
  !-------------------------------------------------------------------------
#define ALLOC_SAFE_FILL0_GETP(A, SZ) \
       if (.not. associated(A)) then ;\
          allocate(A(SZ)) ;\
          A = 0.d0 ;\
       endif ;\
       ptr => A
  
  function get_fb_sum(gid,surfcode) result(ptr)
    integer,intent(IN) :: gid, surfcode
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ptr
    select case( surfcode )
    case (NIL)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%xl_sum, SZ_FXL)
    case (NIR)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%xr_sum, SZ_FXR)
    case (NJL)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%yl_sum, SZ_FYL)
    case (NJR)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%yr_sum, SZ_FYR)
    case (NKL)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%zl_sum, SZ_FZL)
    case (NKR)
       ALLOC_SAFE_FILL0_GETP(Fb(gid)%zr_sum, SZ_FZR)
    end select
  end function get_fb_sum
  ! -----------------------------------------------------------------
  ! 数値流速補正(子供の数値流速を信用する)
  ! -----------------------------------------------------------------
  subroutine fluxcorrection( level )
    integer,intent(IN) :: level
    integer :: n
    if ( level >= LevelMax )  return
    if ( level < 0 ) return
    Levelf = level + 1
    Levelc = level
    call reflux_init
    call com_restrict
    call dealloc_flux_buffer      ! 不要なバッファをクリアする
  end subroutine fluxcorrection
  !-----------------------------------------------------------------------
  ! Deallocate flux buffer
  !-----------------------------------------------------------------------
#define DEALLOC_SAFE(A) \
  if (associated(A)) then ;\
     deallocate(A) ;\
     nullify(A) ;\
  end if

  subroutine dealloc_flux_buffer
    integer :: gid, n
    ! free fb_sum
    do n = Gidmin, GidListMax( Levelf )
       gid = GidList( n , Levelf )
       DEALLOC_SAFE( Fb(gid)%xl_sum )
       DEALLOC_SAFE( Fb(gid)%xr_sum )
       DEALLOC_SAFE( Fb(gid)%yl_sum )
       DEALLOC_SAFE( Fb(gid)%yr_sum )
       DEALLOC_SAFE( Fb(gid)%zl_sum )
       DEALLOC_SAFE( Fb(gid)%zr_sum )
    enddo
    ! free fb_bnd
    do n = Gidmin, GidListMax( Levelc )
       gid = GidList( n , Levelc )
       DEALLOC_SAFE( Fb(gid)%xl_bnd )
       DEALLOC_SAFE( Fb(gid)%xr_bnd )
       DEALLOC_SAFE( Fb(gid)%yl_bnd )
       DEALLOC_SAFE( Fb(gid)%yr_bnd )
       DEALLOC_SAFE( Fb(gid)%zl_bnd )
       DEALLOC_SAFE( Fb(gid)%zr_bnd )
    end do
  end subroutine dealloc_flux_buffer
  !-----------------------------------------------------------------------
  ! 外側の共有領域(袖)の数値流束をバッファに保存する。保存される数値流速は累積である。
  ! ----------------------------------------------------------------------
  subroutine save_flux( id, f )
    integer,intent(IN) :: id
    real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: f
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: fbuf
    real(kind=DBL_KIND) :: rstp
    integer :: n
    integer(kind=LLONG_KIND) :: dstepc, dstepf
    Levelf = levels(id)
    Levelc = Levelf - 1

    if ( Levelf > LevelMax ) return
    ! save flux directly
    do n = 0, NBOUNDARY - 1
       fbuf => get_fb_bnd( id, n )
       fbuf(:,:,:,:) = f(ARRAYSIZE4(fbuf), decode_direction( n ))
    enddo

    if ( Levelc < 0 ) return
    ! save cummulative flux
    dstepf = Dstep( Levelf )
    dstepc = Dstep( Levelc )
    rstp = dble(dstepf)/dble(dstepc)
    do n = 0, NBOUNDARY - 1
       fbuf => get_fb_sum( id, n )
       if ( Step(Levelc)-Dstep(Levelc) == Step(Levelf)-Dstep(Levelf) ) &
            fbuf(:,:,:,:) = 0 !前のステップで同期してたら、クリアする
       fbuf(:,:,:,:) = fbuf(:,:,:,:) + f(ARRAYSIZE4(fbuf), decode_direction( n ) ) * rstp
    enddo

  end subroutine save_flux
  !-----------------------------------------------------------------------
  ! 子を restriction し、転送
  !-----------------------------------------------------------------------
  subroutine com_restrict
    use mpilib
    use packarr
    use grid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: fp
    integer :: n, ndir, rank, gids, gidd, ranks, rankd, pgid, prank, m, lr
    ! 初期化
    myrank = get_myrank()
    call pkar_reset

    do rank = 0, NPE -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          do m = 0, NBOUNDARY -1 !子からみた面
             ! 送受信 gid, rank
             gids  = GidListNode(n, Levelf, rank)
             ranks = rank
             pgid  = ParentGid(gids,ranks)
             prank = ParentRank(gids,ranks)
             ndir = decode_direction( m )
             lr = decode_LR( m )
             gidd  = NeighborGid(lr, ndir, pgid, prank)
             rankd = NeighborRank(lr, ndir, pgid, prank)
             ! reflux が必要なグリッド (1) 自分隣がいない (2) 親に隣(送信先)がいる
             if (.not. ( &
                  NeighborGid(lr, ndir, gids, ranks) == Undefi .and. &
                  NeighborGid(lr, ndir, pgid, prank) /= Undefi ) ) cycle
             if ( myrank == ranks ) then ! 送信 rank
                fp => get_fb_sum( gids, m )
                call restrict_f( fp, FLP(ndir)%Buff, ndir ) ! 送信前に restriction
             endif
             PACK_SEND4( FLP(ndir)%Buff, myrank, ranks, rankd )
          enddo
       enddo
    enddo

    call pkar_sendrecv()

    ! 受信バッファから荷ほどき
    do rank = 0, NPE -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          do m = 0, NBOUNDARY -1 !子からみた面
             ! 送受信 gid, rank
             gids  = GidListNode(n, Levelf, rank)
             ranks = rank
             pgid  = ParentGid(gids,ranks)
             prank = ParentRank(gids,ranks)
             ndir = decode_direction( m )
             lr = decode_LR( m )
             if (.not. ( &
                  NeighborGid(lr, ndir, gids, ranks) == Undefi .and. &
                  NeighborGid(lr, ndir, pgid, prank) /= Undefi ) ) cycle
             gidd  = NeighborGid(lr, ndir, pgid, prank)
             rankd = NeighborRank(lr, ndir, pgid, prank)
             if (myrank /= rankd) cycle
             UNPACK_RECV4( FLP(ndir)%Buff, myrank, ranks, rankd )
             call modify_u_by_reflux(m)
          end do
       enddo
    enddo

  contains
    !-----------------------------------------------------------------------
    ! modify primitive variables u by re-fluxing
    !-----------------------------------------------------------------------
    subroutine modify_u_by_reflux(m)
      use eos
      integer,intent(IN) :: m
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: fc
      real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
      real(kind=DBL_KIND),dimension(MX:MZ) :: ds
      real(kind=DBL_KIND) :: dt, dv
      integer :: lr, lrc, ndir, lri, lrj, lrk

      ds = get_ds( Levelc )
      dv = get_dv( Levelc )
      lr = decode_LR(m)
      lrc = swap_LR( lr )
      ndir = decode_direction(m)
      call left_or_right(gids, ranks, lri, lrj, lrk)
      ! 左右を入れ替える
      fc => get_fb_bnd( gidd, encode_surf( lrc, ndir ) ) !親の数値流速
      u => get_Up( gidd )           !親境界のセル
      dt = Dtime( Levelc )

#define SZ  Ics(lrc,ndir,lri):Ice(lrc,ndir,lri),Jcs(lrc,ndir,lrj):Jce(lrc,ndir,lrj),Kcs(lrc,ndir,lrk):Kce(lrc,ndir,lrk),Mmin:Mmax
#define SZA Ias(lrc,ndir,lri):Iae(lrc,ndir,lri),Jas(lrc,ndir,lrj):Jae(lrc,ndir,lrj),Kas(lrc,ndir,lrk):Kae(lrc,ndir,lrk),Mmin:Mmax

      FLP(ndir)%Bufw = u(SZ)

      call conv_u2w(FLP(ndir)%Bufw, dv)
      if ( lr == Left ) then
         FLP(ndir)%Bufw  = FLP(ndir)%Bufw - (FLP(ndir)%Buff - fc(SZA) )*dt*ds(ndir)
      else
         FLP(ndir)%Bufw  = FLP(ndir)%Bufw + (FLP(ndir)%Buff - fc(SZA) )*dt*ds(ndir)
      endif
      call conv_w2u(FLP(ndir)%Bufw, dv)
      u(SZ) = FLP(ndir)%Bufw

#undef SZ
#undef SZA
    end subroutine modify_u_by_reflux

  end subroutine com_restrict

  !-----------------------------------------------------------------------
  ! restrict numerical flux (f)  from fine to coarse grid
  !-----------------------------------------------------------------------
  subroutine restrict_f( ff, fc, ncrd )
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ff !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: fc !(OUT)
    integer,intent(IN) :: ncrd                           !(IN)
    integer :: if,jf,kf,ic,jc,kc,ic0,jc0,kc0,if0,jf0,kf0,m
    real(kind=DBL_KIND),dimension(MX:MZ) :: dsf, dsc
    ! check size of array
    if ( ncrd == MX .and. &
         ( &
         size(ff,2) /= size(fc,2)*2 .or. &
         size(ff,3) /= size(fc,3)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NX', size(ff), size(fc)
       stop
    endif
    if ( ncrd == MY .and. &
         ( &
         size(ff,1) /= size(fc,1)*2 .or. &
         size(ff,3) /= size(fc,3)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NY', size(ff), size(fc)
       stop
    endif
    if ( ncrd == MZ .and. &
         ( &
         size(ff,1) /= size(fc,1)*2 .or. &
         size(ff,2) /= size(fc,2)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NZ', size(ff), size(fc)
       stop
    endif
!!$    dsf = 1.d0                  ! volume for fine grid
!!$    dsc = dsf*4                 ! volume for coarse grid
    dsf = get_ds( Levelf )
    dsc = get_ds( Levelc )

    ic0 = lbound(fc,1)
    jc0 = lbound(fc,2)
    kc0 = lbound(fc,3)
    if0 = lbound(ff,1)
    jf0 = lbound(ff,2)
    kf0 = lbound(ff,3)

    if ( ncrd == MX ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if,  jf,  kf,  m) + ff(if,  jf,  kf+1,m) &
                        +ff(if,  jf+1,kf,  m) + ff(if,  jf+1,kf+1,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
    if ( ncrd == MY ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if,  jf,  kf,  m) + ff(if,  jf,  kf+1,m) &
                        +ff(if+1,jf,  kf,  m) + ff(if+1,jf,  kf+1,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
    if ( ncrd == MZ ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if,  jf,  kf,  m) + ff(if,  jf+1,  kf,m) &
                        +ff(if+1,jf,  kf,  m) + ff(if+1,jf+1,  kf,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
  end subroutine restrict_f
  !-----------------------------------------------------------------------
  ! write flux buffer
  !-----------------------------------------------------------------------
#define FB_WRITE(A, SZ) \
  if (.not. associated(A)) then ;\
     allocate(A(SZ)) ;\
     A = 0.d0 ;\
  end if;\
  write(unit) A

  subroutine reflux_write(unit, gid)
    integer,intent(IN) :: unit, gid
    FB_WRITE(Fb(gid)%xl_bnd, SZ_FXL)
    FB_WRITE(Fb(gid)%xr_bnd, SZ_FXR)
    FB_WRITE(Fb(gid)%yl_bnd, SZ_FYL)
    FB_WRITE(Fb(gid)%yr_bnd, SZ_FYR)
    FB_WRITE(Fb(gid)%zl_bnd, SZ_FZL)
    FB_WRITE(Fb(gid)%zr_bnd, SZ_FZR)
    FB_WRITE(Fb(gid)%xl_sum, SZ_FXL)
    FB_WRITE(Fb(gid)%xr_sum, SZ_FXR)
    FB_WRITE(Fb(gid)%yl_sum, SZ_FYL)
    FB_WRITE(Fb(gid)%yr_sum, SZ_FYR)
    FB_WRITE(Fb(gid)%zl_sum, SZ_FZL)
    FB_WRITE(Fb(gid)%zr_sum, SZ_FZR)
  end subroutine reflux_write
  !-----------------------------------------------------------------------
  ! read flux buffer
  !-----------------------------------------------------------------------
#define FB_READ(A, SZ) \
  if (.not. associated(A)) then ;\
     allocate(A(SZ)) ;\
  end if ;\
  read(unit, iostat=eof) A

  subroutine reflux_read(unit, gid, eof )
    integer,intent(IN) :: unit, gid
    integer,intent(OUT) :: eof
    FB_READ(Fb(gid)%xl_bnd, SZ_FXL)
    FB_READ(Fb(gid)%xr_bnd, SZ_FXR)
    FB_READ(Fb(gid)%yl_bnd, SZ_FYL)
    FB_READ(Fb(gid)%yr_bnd, SZ_FYR)
    FB_READ(Fb(gid)%zl_bnd, SZ_FZL)
    FB_READ(Fb(gid)%zr_bnd, SZ_FZR)
    FB_READ(Fb(gid)%xl_sum, SZ_FXL)
    FB_READ(Fb(gid)%xr_sum, SZ_FXR)
    FB_READ(Fb(gid)%yl_sum, SZ_FYL)
    FB_READ(Fb(gid)%yr_sum, SZ_FYR)
    FB_READ(Fb(gid)%zl_sum, SZ_FZL)
    FB_READ(Fb(gid)%zr_sum, SZ_FZR)
  end subroutine reflux_read
  !-----------------------------------------------------------------------
  ! initialize
  !-----------------------------------------------------------------------
  subroutine reflux_init
    use io_util, only : print_msg
    integer :: lr, ndir, lri, lrj, lrk
    logical,save :: bool_initialized = .false.
    if ( bool_initialized ) return
    bool_initialized = .true.
    call print_msg( 'initialize reflux' )

    ! make array for a cell center
    ! ics(LR surf, direction, parity of this grid)
    !     LR,      MX,MY,MZ,  Ev,Od
    Ics(:,:,Ev) = Imin
    Ice(:,:,Ev) = (Imax-Imin+1)/2+Imin-1
    Ics(:,:,Od) = (Imax-Imin+1)/2+Imin
    Ice(:,:,Od) = Imax

    Jcs(:,:,Ev) = Jmin
    Jce(:,:,Ev) = (Jmax-Jmin+1)/2+Jmin-1
    Jcs(:,:,Od) = (Jmax-Jmin+1)/2+Jmin
    Jce(:,:,Od) = Jmax

    Kcs(:,:,Ev) = Kmin
    Kce(:,:,Ev) = (Kmax-Kmin+1)/2+Kmin-1
    Kcs(:,:,Od) = (Kmax-Kmin+1)/2+Kmin
    Kce(:,:,Od) = Kmax

    Ics(L,MX,:) = Imin
    Ice(L,MX,:) = Imin
    Ics(R,MX,:) = Imax
    Ice(R,MX,:) = Imax

    Jcs(L,MY,:) = Jmin
    Jce(L,MY,:) = Jmin
    Jcs(R,MY,:) = Jmax
    Jce(R,MY,:) = Jmax

    Kcs(L,MZ,:) = Kmin
    Kce(L,MZ,:) = Kmin
    Kcs(R,MZ,:) = Kmax
    Kce(R,MZ,:) = Kmax

    ! make array for a cell boundary
    Ias(:,:,Ev) = Imin
    Iae(:,:,Ev) = (Imax-Imin+1)/2+Imin-1
    Ias(:,:,Od) = (Imax-Imin+1)/2+Imin
    Iae(:,:,Od) = Imax

    Jas(:,:,Ev) = Jmin
    Jae(:,:,Ev) = (Jmax-Jmin+1)/2+Jmin-1
    Jas(:,:,Od) = (Jmax-Jmin+1)/2+Jmin
    Jae(:,:,Od) = Jmax

    Kas(:,:,Ev) = Kmin
    Kae(:,:,Ev) = (Kmax-Kmin+1)/2+Kmin-1
    Kas(:,:,Od) = (Kmax-Kmin+1)/2+Kmin
    Kae(:,:,Od) = Kmax

    Ias(L,MX,:) = Imin-1
    Iae(L,MX,:) = Imin-1
    Ias(R,MX,:) = Imax
    Iae(R,MX,:) = Imax

    Jas(L,MY,:) = Jmin-1
    Jae(L,MY,:) = Jmin-1
    Jas(R,MY,:) = Jmax
    Jae(R,MY,:) = Jmax

    Kas(L,MZ,:) = Kmin-1
    Kae(L,MZ,:) = Kmin-1
    Kas(R,MZ,:) = Kmax
    Kae(R,MZ,:) = Kmax

#define SZA Ias(lr,ndir,lri):Iae(lr,ndir,lri),Jas(lr,ndir,lrj):Jae(lr,ndir,lrj),Kas(lr,ndir,lrk):Kae(lr,ndir,lrk),Mmin:Mmax
#define SZ  Ics(lr,ndir,lri):Ice(lr,ndir,lri),Jcs(lr,ndir,lrj):Jce(lr,ndir,lrj),Kcs(lr,ndir,lrk):Kce(lr,ndir,lrk),Mmin:Mmax
    lr = Left
    lri = Ev
    lrj = Ev
    lrk = Ev
    if ( .not. associated(FLP(MX)%Buff) ) then
       do ndir = MX, MZ
          allocate( FLP(ndir)%Buff( SZA ) )
          allocate( FLP(ndir)%Bufw( SZ ) )
       enddo
    endif
#undef SZ
  end subroutine reflux_init
end module reflux
