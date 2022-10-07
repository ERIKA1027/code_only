module fmg_converge
  use mpilib
  use fmg_data
  implicit none
  private
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: Bufc => null()
  end type t_buf
  type(t_buf),save,dimension(:),allocatable :: FLV 
  type(t_buf),save,dimension(:),allocatable :: FLS 
  integer,parameter :: STENCIL_MIN = 0
  integer,parameter :: STENCIL_MAX = 3
  real(kind=8),save,dimension(STENCIL_MIN:STENCIL_MAX,STENCIL_MIN:STENCIL_MAX,STENCIL_MIN:STENCIL_MAX) :: St 
  integer,save :: FMG_Level = Undefi
  public :: fmg_converge_c2p, fmg_converge_c2p_lev, fmg_converge_init, fmg_converge_finalize
contains
  subroutine fmg_converge_init(fmglev)
    integer,intent(IN) :: fmglev
    if ( .not. allocated(FLV) ) allocate( FLV(FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLV(fmglev)%Bufc) ) &
         allocate( FLV(fmglev)%Bufc( &
         GridSize(fmglev)%Imin: (GridSize(fmglev)%Imax-GridSize(fmglev)%Imin+1)/2+GridSize(fmglev)%Imin-1, &
         GridSize(fmglev)%Jmin: (GridSize(fmglev)%Jmax-GridSize(fmglev)%Jmin+1)/2+GridSize(fmglev)%Jmin-1, &
         GridSize(fmglev)%Kmin: (GridSize(fmglev)%Kmax-GridSize(fmglev)%Kmin+1)/2+GridSize(fmglev)%Kmin-1, &
         Mmin:Mmax) )
    if ( .not. allocated(FLS) ) allocate( FLS(FMG_LevelMin:FMG_LevelMax) )
    if ( .not. associated(FLS(fmglev)%Bufc) ) &
         allocate( FLS(fmglev)%Bufc( &
         GridSize(fmglev)%Imin: (GridSize(fmglev)%Imax-GridSize(fmglev)%Imin+1)/2+GridSize(fmglev)%Imin-1, &
         GridSize(fmglev)%Jmin: (GridSize(fmglev)%Jmax-GridSize(fmglev)%Jmin+1)/2+GridSize(fmglev)%Jmin-1, &
         GridSize(fmglev)%Kmin: (GridSize(fmglev)%Kmax-GridSize(fmglev)%Kmin+1)/2+GridSize(fmglev)%Kmin-1, &
         Mmin:Mmin) )
    FMG_Level = fmglev
  end subroutine fmg_converge_init
  subroutine fmg_converge_finalize
    integer :: fmglev
    if (allocated(FLV)) then
       do fmglev = lbound(FLV,1), ubound(FLV,1)
          if (associated(FLV(fmglev)%Bufc)) deallocate(FLV(fmglev)%Bufc)
          nullify(FLV(fmglev)%Bufc)
       end do
       deallocate(FLV)
    endif
    if (allocated(FLS)) then
       do fmglev = lbound(FLS,1), ubound(FLS,1)
          if (associated(FLS(fmglev)%Bufc)) deallocate(FLS(fmglev)%Bufc)
          nullify(FLS(fmglev)%Bufc)
       end do
       deallocate(FLS)
    endif
  end subroutine fmg_converge_finalize
  subroutine fmg_converge_c2p(fmglev,icode,cubic)
    integer,intent(IN) :: fmglev, icode
    logical,intent(IN),optional :: cubic
    integer :: amrlevc
    call fmg_converge_init(fmglev)
    do amrlevc = AMR_LevelMax-1, AMR_LevelMin, -1
       call fmg_converge_c2p_lev(amrlevc,fmglev,icode,cubic)
    end do
  end subroutine fmg_converge_c2p
  subroutine fmg_converge_c2p_lev(amrlevc,fmglev,icode, cubic)
    use restriction
    use mpilib
    use packarr
    integer,intent(IN) :: amrlevc,fmglev, icode
    logical,intent(IN),optional :: cubic
    integer :: gid, rank, gidd, gids, ranks, rankd, amrlevf
    real(kind=8),dimension(:,:,:,:),pointer :: bufc, uf
    logical :: bool_cubic
    bool_cubic = .FALSE. 
    if (present(cubic)) bool_cubic = cubic
    amrlevf = amrlevc + 1 
    myrank = get_myrank()
    if ( fmg_isVector(icode) ) then
       bufc => FLV(fmglev)%Bufc
    else
       bufc => FLS(fmglev)%Bufc
    endif
    call pkar_reset
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlevf,rank), fmg_get_gidmax_rank(amrlevf,rank)
          gidd = Ancestry(amrlevf,rank)%Block(gid)%ParentGid
          rankd = Ancestry(amrlevf,rank)%Block(gid)%ParentRank
          gids = gid
          ranks = rank
          if (ranks == myrank) then
             uf => fmg_get_arrp(amrlevf, fmglev, gids, icode)
             call rstrct(bufc, uf, &
                  int(lbound(bufc,1)), int(lbound(bufc,2)), int(lbound(bufc,3)), &
                  int(ubound(bufc,1)), int(ubound(bufc,2)), int(ubound(bufc,3)), &
                  GridSize(fmglev)%Imin, GridSize(fmglev)%Jmin, GridSize(fmglev)%Kmin, cubic)
             if (rankd == myrank) call put_buf2uc 
          endif
          if ( (myrank == ranks .or. myrank == rankd) .and. rankd /= ranks) then
             if ((myrank) == (ranks)) call pkar_push(bufc(lbound(bufc,1),lbound(bufc,2),lbound(bufc,3),lbound(bufc,4)), size(bufc),&
& kind(bufc), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(bufc), kind(bufc), ranks)
          end if
       end do
    end do
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlevf,rank), fmg_get_gidmax_rank(amrlevf,rank)
          gidd = Ancestry(amrlevf,rank)%Block(gid)%ParentGid
          rankd = Ancestry(amrlevf,rank)%Block(gid)%ParentRank
          gids = gid
          ranks = rank
          if (rankd /= myrank) cycle
          if (ranks == myrank) cycle 
          if ((myrank) == (rankd)) call pkar_pop(bufc(lbound(bufc,1),lbound(bufc,2),lbound(bufc,3),lbound(bufc,4)), size(bufc), kin&
&d(bufc), ranks)
          call put_buf2uc
       end do
    end do
  contains
    subroutine put_buf2uc
      integer :: is, ie, js, je, ks, ke, lri, lrj, lrk
      real(kind=8),dimension(:,:,:,:),pointer :: uc
      call fmg_left_or_right(gids, ranks, amrlevf, lri, lrj, lrk)
      is = GridSize(fmglev)%Imin
      ie = GridSize(fmglev)%Imax
      js = GridSize(fmglev)%Jmin
      je = GridSize(fmglev)%Jmax
      ks = GridSize(fmglev)%Kmin
      ke = GridSize(fmglev)%Kmax
      if (lri == Left) ie = (ie-is+1)/2-1 + is
      if (lrj == Left) je = (je-js+1)/2-1 + js
      if (lrk == Left) ke = (ke-ks+1)/2-1 + ks
      if (lri == Right) is = (ie-is+1)/2 + is
      if (lrj == Right) js = (je-js+1)/2 + js
      if (lrk == Right) ks = (ke-ks+1)/2 + ks
      uc => fmg_get_arrp(amrlevc, fmglev, gidd, icode)
      uc(is:ie,js:je,ks:ke,:) = bufc
    end subroutine put_buf2uc
  end subroutine fmg_converge_c2p_lev
end module fmg_converge
