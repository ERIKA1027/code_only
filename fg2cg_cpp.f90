module fg2cg
  use grid, only : Gidmin, Left, Right, Imin, Imax, Jmin, Jmax, Kmin, Kmax, Mmin, Mmax
  implicit none
  private
  integer,save :: Levelf
  integer,save :: Levelc
  integer,save,dimension(Left:Right) :: Ifs, Ife, Jfs, Jfe, Kfs, Kfe
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(Left:Right,Left:Right,Left:Right) :: FLP
  public :: fg2cg_u
contains
  subroutine fg2cg_u( level )
    integer,intent(IN) :: level 
    if (level < 1) return
    call fg2cg_init( level )
    call com_restrict
  end subroutine fg2cg_u
  subroutine fg2cg_init( level )
    use io_util, only : print_msg
    integer,intent(IN) :: level
    logical,save :: bool_initialized = .false.
    integer :: lri, lrj, lrk
    Levelf = level
    Levelc = level-1
    if ( bool_initialized ) return
    bool_initialized = .true.
    call print_msg( 'initialize fg2cg' )
    Ifs(Left) = Imin
    Ife(Left) = (Imax-Imin+1)/2+Imin-1
    Ifs(Right) = (Imax-Imin+1)/2+Imin
    Ife(Right) = Imax
    Jfs(Left) = Jmin
    Jfe(Left) = (Jmax-Jmin+1)/2+Jmin-1
    Jfs(Right) = (Jmax-Jmin+1)/2+Jmin
    Jfe(Right) = Jmax
    Kfs(Left) = Kmin
    Kfe(Left) = (Kmax-Kmin+1)/2+Kmin-1
    Kfs(Right) = (Kmax-Kmin+1)/2+Kmin
    Kfe(Right) = Kmax
    if ( .not. associated(FLP(Left,Left,Left)%Bufc) ) then
       do lrk = Left, Right
          do lrj = Left, Right
             do lri = Left, Right
                allocate( FLP(lri,lrj,lrk)%Bufc( ifs(lri):ife(lri),jfs(lrj):jfe(lrj),kfs(lrk):kfe(lrk),Mmin:Mmax ) )
             end do
          enddo
       enddo
    endif
  end subroutine fg2cg_init
  subroutine com_restrict
    use mpilib
    use packarr
    use grid
    real(kind=8),dimension(:,:,:,:),pointer :: up
    real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax),target :: ubuf
    integer :: n, rank, lri, lrj, lrk, i, j, k, ranks, rankd, gids, gidd
    call pkar_reset
    myrank = get_myrank()
    lri = -1
    lrj = -1
    lrk = -1
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          gids = GidListNode(n, Levelf, rank)
          ranks = rank
          gidd = ParentGid(gids,ranks)
          rankd = ParentRank(gids,ranks)
          if (myrank == ranks .or. myrank == rankd) then
             call left_or_right(gids, ranks, lri, lrj, lrk)
          endif
          if ( myrank == ranks ) then 
             up => get_Up(gids)
             ubuf = up(Imin:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)
             up => ubuf 
             call restrict_u( up, FLP(lri,lrj,lrk)%Bufc ) 
          endif
          if ((myrank) == (ranks)) call pkar_push(FLP(lri,lrj,lrk)%Bufc(lbound(FLP(lri,lrj,lrk)%Bufc,1),lbound(FLP(lri,lrj,lrk)%Buf&
&c,2),lbound(FLP(lri,lrj,lrk)%Bufc,3),lbound(FLP(lri,lrj,lrk)%Bufc,4)), size(FLP(lri,lrj,lrk)%Bufc), kind(FLP(lri,lrj,lrk)%Bufc), r&
&ankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(FLP(lri,lrj,lrk)%Bufc), kind(FLP(lri,lrj,lrk)%Bufc), ranks)
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          gids = GidListNode(n, Levelf, rank)
          ranks = rank
          gidd = ParentGid(gids,ranks)
          rankd = ParentRank(gids,ranks)
          if (myrank == rankd) then
             call left_or_right(gids, ranks, lri, lrj, lrk)
             if ((myrank) == (rankd)) call pkar_pop(FLP(lri,lrj,lrk)%Bufc(lbound(FLP(lri,lrj,lrk)%Bufc,1),lbound(FLP(lri,lrj,lrk)%B&
&ufc,2),lbound(FLP(lri,lrj,lrk)%Bufc,3),lbound(FLP(lri,lrj,lrk)%Bufc,4)), size(FLP(lri,lrj,lrk)%Bufc), kind(FLP(lri,lrj,lrk)%Bufc),&
& ranks)
             up => get_Up(gidd)
             up(ifs(lri):ife(lri),jfs(lrj):jfe(lrj),kfs(lrk):kfe(lrk),Mmin:Mmax) = FLP(lri,lrj,lrk)%Bufc
          end if
       enddo
    enddo
  end subroutine com_restrict
  subroutine restrict_u( uf, uc )
    use eos
    real(kind=8),dimension(:,:,:,:),pointer :: uf 
    real(kind=8),dimension(:,:,:,:),pointer :: uc 
    real(kind=8),dimension(:,:,:,:),pointer,save :: wf
    real(kind=8) :: dvf, dvc
    integer :: if,jf,kf,ic,jc,kc,ic0,jc0,kc0,if0,jf0,kf0,m
    logical,save :: bool_initialized = .false.
    if ( .not. bool_initialized ) then
       bool_initialized = .true.
       allocate( wf(lbound(uf,1):ubound(uf,1),lbound(uf,2):ubound(uf,2),lbound(uf,3):ubound(uf,3),lbound(uf,4):ubound(uf,4)) )
    endif
    if ( &
         size(uf,1) /= size(uc,1)*2 .or. &
         size(uf,2) /= size(uc,2)*2 .or. &
         size(uf,3) /= size(uc,3)*2 .or. &
         size(uf,4) /= size(uc,4) ) then
       write(*,*) 'restrict_u: uf and uc are incompatible size ', size(uf), size(uc)
       stop
    endif
    dvf = 1.d0 
    dvc = dvf*8 
    call u2w( uf, wf, dvf )
    ic0 = lbound(uc,1)
    jc0 = lbound(uc,2)
    kc0 = lbound(uc,3)
    if0 = lbound(uf,1)
    jf0 = lbound(uf,2)
    kf0 = lbound(uf,3)
    do m= lbound(uc,4), ubound(uc,4)
       do kc = lbound(uc,3), ubound(uc,3)
          do jc = lbound(uc,2), ubound(uc,2)
             do ic = lbound(uc,1), ubound(uc,1)
                if = (ic - ic0) * 2 + if0
                jf = (jc - jc0) * 2 + jf0
                kf = (kc - kc0) * 2 + kf0
                uc(ic,jc,kc,m) = &
                      wf(if, jf, kf, m)+wf(if+1,jf, kf, m) &
                     +wf(if, jf+1,kf, m)+wf(if, jf, kf+1,m) &
                     +wf(if+1,jf+1,kf, m)+wf(if+1,jf, kf+1,m) &
                     +wf(if, jf+1,kf+1,m)+wf(if+1,jf+1,kf+1,m)
             enddo
          enddo
       enddo
    enddo
    call conv_w2u( uc, dvc ) 
  end subroutine restrict_u
end module fg2cg
