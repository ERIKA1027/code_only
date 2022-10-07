#include "config.h"
! #define VERBOSE
!-------------------------------------------------------------------------
! Module for convergence between parent and child grids.
!-------------------------------------------------------------------------
module vmg_interpol
  use fmg_data
  implicit none
  private
  real(kind=DBL_KIND),save,dimension(:,:,:,:),allocatable,target :: BufcS, BufcV
  real(kind=DBL_KIND),save,dimension(:,:,:,:),pointer :: bufx => null(), bufxy => null()
  integer,save :: FMG_Level = Undefi
  logical,save :: Initialized = .FALSE.
  integer,save,dimension(Left:Right) :: Is, Ie, Js, Je, Ks, Ke ! lower and higher bound of sub-block
  integer,save :: Imin,Jmin,Kmin,Imax,Jmax,Kmax

  public :: vmg_interp_p2c_lev, vmg_interp_init, vmg_interp_finalize
contains
  !-------------------------------------------------------------------------
  ! initialize
  !-------------------------------------------------------------------------
  subroutine vmg_interp_init(fmglev)
    integer,intent(IN) :: fmglev
    integer :: szi, szj, szk
    integer :: jcmin, jcmax, kcmin, kcmax
    if ( Initialized ) return
    Initialized = .TRUE.

    ! index of boundary of parent block
    call fmg_get_gridsize(fmglev, Imin,Jmin,Kmin,Imax,Jmax,Kmax)
    szi = (Imax - Imin + 1)/2
    szj = (Jmax - Jmin + 1)/2
    szk = (Kmax - Kmin + 1)/2
    Is = (/ Imin-Ngh,     Imin+szi-Ngh /)
    Ie = (/ Imin+szi-1+Ngh, Imax+Ngh /)
    Js = (/ Jmin-Ngh,     Jmin+szj-Ngh /)
    Je = (/ Jmin+szj-1+Ngh, Jmax+Ngh /)
    Ks = (/ Kmin-Ngh,     Kmin+szk-Ngh /)
    Ke = (/ Kmin+szk-1+Ngh, Kmax+Ngh /)

    ! buffer for trilinear interpolation
    if ( .not. allocated(BufcS) ) &
         allocate( BufcS( Is(Left):Ie(Left), Js(Left):Je(Left), Ks(Left):Ke(Left), Mmin:Mmin ) )
    if ( .not. allocated(BufcV) ) &
         allocate( BufcV( Is(Left):Ie(Left), Js(Left):Je(Left), Ks(Left):Ke(Left), Mmin:Mmax ) )

    ! buffer for tricubic interpolation
    jcmin = IJKC(Jmin, Jmin)-Ngh
    jcmax = IJKC(Jmax, Jmin)+Ngh
    kcmin = IJKC(Kmin, Kmin)-Ngh
    kcmax = IJKC(Kmax, Kmin)+Ngh
    allocate( bufx (Imin:Imax, jcmin:jcmax, kcmin:kcmax, Mmin:Mmax) )
    allocate( bufxy(Imin:Imax,  Jmin:Jmax,  kcmin:kcmax, Mmin:Mmax) )

    FMG_Level = fmglev
#ifdef VERBOSE
    print *, '*** vmg_interpol: initialized'
#endif !VERBOSE
  end subroutine vmg_interp_init
  !-------------------------------------------------------------------------
  ! finalize
  !-------------------------------------------------------------------------
  subroutine vmg_interp_finalize
    Initialized = .FALSE.
    if (allocated(BufcS)) deallocate( BufcS )
    if (allocated(BufcV)) deallocate( BufcV )
    if (associated(bufx))  then
       deallocate(bufx)
       nullify(bufx)
    end if
   if (associated(bufxy)) then
      deallocate(bufxy)
      nullify(bufxy)
   end if
#ifdef VERBOSE
    print *, '*** vmg_interpol: finalize'
#endif !VERBOSE
  end subroutine vmg_interp_finalize
  !-------------------------------------------------------------------------
  ! interp child grid to parentgrid (fg2cg)
  !-------------------------------------------------------------------------
  subroutine vmg_interp_p2c_lev(amrlevf,juf,juc, cubic)
    use mpilib
    use packarr
    use interpolation
#include "packarr.h"
    integer,intent(IN) :: amrlevf, juf, juc
    logical,intent(IN),optional :: cubic
    integer :: fmglev, gid, rank, lri, lrj, lrk, &
         gidd, gids, ranks, rankd, amrlevc
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: uc, uf, buf
    logical :: bool_cubic

    bool_cubic = .false.        ! cubic interpolation?
    if (present(cubic)) bool_cubic = cubic

#ifdef VERBOSE
    if (bool_cubic .and. get_myrank() == PRIMARY_RANK) &
         print *, '*** vmg_interpol: cubic interpolation'
#endif !VERBOSE

    amrlevc = amrlevf - 1       ! 子グリッドレベル
    fmglev = FMG_Level
    myrank = get_myrank()

    if ( fmg_isVector(juf) ) then ! select buffer
       buf => BufcV
    else
       buf => BufcS
    end if
    call pkar_reset
    ! 親グリッドを荷造りする。
    do rank = 0, NPE-1
       do gid = fmg_get_gidmin_rank(amrlevc,rank), fmg_get_gidmax_rank(amrlevc,rank)
          do lrk = Left, Right
             do lrj = Left, Right
                do lri = Left, Right
                   gidd  = Ancestry(amrlevc,rank)%Block(gid)%ChildGid(lri,lrj,lrk)
                   rankd = Ancestry(amrlevc,rank)%Block(gid)%ChildRank(lri,lrj,lrk)
                   gids = gid
                   ranks = rank
                   if (gidd == Undefi) cycle
                   if (ranks == myrank) then
                      uc => fmg_get_arrp(amrlevc, fmglev, gids, juc)
                      buf = uc( Is(lri):Ie(lri), Js(lrj):Je(lrj), Ks(lrk):Ke(lrk),: )
                   endif
                   PACK_SEND4( buf, myrank, ranks, rankd )
                enddo
             enddo
          enddo
       end do
    end do
    call pkar_sendrecv()
    ! pop
    do rank = 0, NPE-1
       do gid = fmg_get_gidmin_rank(amrlevc,rank), fmg_get_gidmax_rank(amrlevc,rank)
          do lrk = Left, Right
             do lrj = Left, Right
                do lri = Left, Right
                   gidd  = Ancestry(amrlevc,rank)%Block(gid)%ChildGid(lri,lrj,lrk)
                   rankd = Ancestry(amrlevc,rank)%Block(gid)%ChildRank(lri,lrj,lrk)
                   gids = gid
                   ranks = rank
                   if (gidd == Undefi) cycle
                   if (rankd /= myrank) cycle
                   UNPACK_RECV4( buf, myrank, ranks, rankd )
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
  !-------------------------------------------------------------------------

end module vmg_interpol
