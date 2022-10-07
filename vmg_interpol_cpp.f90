module vmg_interpol
  use fmg_data
  implicit none
  private
  real(kind=8),save,dimension(:,:,:,:),allocatable,target :: BufcS, BufcV
  real(kind=8),save,dimension(:,:,:,:),pointer :: bufx => null(), bufxy => null()
  integer,save :: FMG_Level = Undefi
  logical,save :: Initialized = .FALSE.
  integer,save,dimension(Left:Right) :: Is, Ie, Js, Je, Ks, Ke 
  integer,save :: Imin,Jmin,Kmin,Imax,Jmax,Kmax
  public :: vmg_interp_p2c_lev, vmg_interp_init, vmg_interp_finalize
contains
  subroutine vmg_interp_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: szi, szj, szk
    integer :: jcmin, jcmax, kcmin, kcmax
    if ( Initialized ) return
    Initialized = .TRUE.
    call fmg_get_gridsize(fmglev, Imin,Jmin,Kmin,Imax,Jmax,Kmax)
    szi = (Imax - Imin + 1)/2
    szj = (Jmax - Jmin + 1)/2
    szk = (Kmax - Kmin + 1)/2
    Is = (/ Imin-Ngh, Imin+szi-Ngh /)
    Ie = (/ Imin+szi-1+Ngh, Imax+Ngh /)
    Js = (/ Jmin-Ngh, Jmin+szj-Ngh /)
    Je = (/ Jmin+szj-1+Ngh, Jmax+Ngh /)
    Ks = (/ Kmin-Ngh, Kmin+szk-Ngh /)
    Ke = (/ Kmin+szk-1+Ngh, Kmax+Ngh /)
    if ( .not. allocated(BufcS) ) &
         allocate( BufcS( Is(Left):Ie(Left), Js(Left):Je(Left), Ks(Left):Ke(Left), Mmin:Mmin ) )
    if ( .not. allocated(BufcV) ) &
         allocate( BufcV( Is(Left):Ie(Left), Js(Left):Je(Left), Ks(Left):Ke(Left), Mmin:Mmax ) )
    jcmin = ((Jmin)-(Jmin))/2 + mod(min((Jmin)-(Jmin),0),2) + (Jmin)-Ngh
    jcmax = ((Jmax)-(Jmin))/2 + mod(min((Jmax)-(Jmin),0),2) + (Jmin)+Ngh
    kcmin = ((Kmin)-(Kmin))/2 + mod(min((Kmin)-(Kmin),0),2) + (Kmin)-Ngh
    kcmax = ((Kmax)-(Kmin))/2 + mod(min((Kmax)-(Kmin),0),2) + (Kmin)+Ngh
    allocate( bufx (Imin:Imax, jcmin:jcmax, kcmin:kcmax, Mmin:Mmax) )
    allocate( bufxy(Imin:Imax, Jmin:Jmax, kcmin:kcmax, Mmin:Mmax) )
    FMG_Level = fmglev
  end subroutine vmg_interp_init
  subroutine vmg_interp_finalize
    Initialized = .FALSE.
    if (allocated(BufcS)) deallocate( BufcS )
    if (allocated(BufcV)) deallocate( BufcV )
    if (associated(bufx)) then
       deallocate(bufx)
       nullify(bufx)
    end if
   if (associated(bufxy)) then
      deallocate(bufxy)
      nullify(bufxy)
   end if
  end subroutine vmg_interp_finalize
  subroutine vmg_interp_p2c_lev(amrlevf,juf,juc, cubic)
    use mpilib
    use packarr
    use interpolation
    integer,intent(IN) :: amrlevf, juf, juc
    logical,intent(IN),optional :: cubic
    integer :: fmglev, gid, rank, lri, lrj, lrk, &
         gidd, gids, ranks, rankd, amrlevc
    real(kind=8),dimension(:,:,:,:),pointer :: uc, uf, buf
    logical :: bool_cubic
    bool_cubic = .false. 
    if (present(cubic)) bool_cubic = cubic
    amrlevc = amrlevf - 1 
    fmglev = FMG_Level
    myrank = get_myrank()
    if ( fmg_isVector(juf) ) then 
       buf => BufcV
    else
       buf => BufcS
    end if
    call pkar_reset
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlevc,rank), fmg_get_gidmax_rank(amrlevc,rank)
          do lrk = Left, Right
             do lrj = Left, Right
                do lri = Left, Right
                   gidd = Ancestry(amrlevc,rank)%Block(gid)%ChildGid(lri,lrj,lrk)
                   rankd = Ancestry(amrlevc,rank)%Block(gid)%ChildRank(lri,lrj,lrk)
                   gids = gid
                   ranks = rank
                   if (gidd == Undefi) cycle
                   if (ranks == myrank) then
                      uc => fmg_get_arrp(amrlevc, fmglev, gids, juc)
                      buf = uc( Is(lri):Ie(lri), Js(lrj):Je(lrj), Ks(lrk):Ke(lrk),: )
                   endif
                   if ((myrank) == (ranks)) call pkar_push(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf),&
& kind(buf), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(buf), kind(buf), ranks)
                enddo
             enddo
          enddo
       end do
    end do
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do gid = fmg_get_gidmin_rank(amrlevc,rank), fmg_get_gidmax_rank(amrlevc,rank)
          do lrk = Left, Right
             do lrj = Left, Right
                do lri = Left, Right
                   gidd = Ancestry(amrlevc,rank)%Block(gid)%ChildGid(lri,lrj,lrk)
                   rankd = Ancestry(amrlevc,rank)%Block(gid)%ChildRank(lri,lrj,lrk)
                   gids = gid
                   ranks = rank
                   if (gidd == Undefi) cycle
                   if (rankd /= myrank) cycle
                   if ((myrank) == (rankd)) call pkar_pop(buf(lbound(buf,1),lbound(buf,2),lbound(buf,3),lbound(buf,4)), size(buf), &
&kind(buf), ranks)
                   uf => fmg_get_arrp(amrlevf, fmglev, gidd, juf)
                   if (bool_cubic) then
                      call interp_tricubic(uf, buf, bufx, bufxy, &
                           Imin, Jmin, Kmin, Imax, Jmax, Kmax, &
                           Imin, Jmin, Kmin)
                   else
                      call interp_trilinear(uf, buf, &
                           Imin, Jmin, Kmin, Imax, Jmax, Kmax, &
                           Imin, Jmin, Kmin)
                   endif
                enddo
             enddo
          enddo
       end do
    end do
  end subroutine vmg_interp_p2c_lev
end module vmg_interpol
