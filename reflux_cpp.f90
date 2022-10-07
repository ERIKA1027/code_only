module reflux
  use grid
  implicit none
  private
  integer,parameter :: NBOUNDARY = 2*(2 -0 +1) 
  integer,parameter :: NIL = 0
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  integer,parameter :: L = Left, R = Right, Ev = Left, Od = Right
  type t_flux_boundary
     real(kind=8),dimension(:,:,:,:),pointer :: xl_bnd => null() 
     real(kind=8),dimension(:,:,:,:),pointer :: xr_bnd => null()
     real(kind=8),dimension(:,:,:,:),pointer :: yl_bnd => null()
     real(kind=8),dimension(:,:,:,:),pointer :: yr_bnd => null()
     real(kind=8),dimension(:,:,:,:),pointer :: zl_bnd => null()
     real(kind=8),dimension(:,:,:,:),pointer :: zr_bnd => null()
     real(kind=8),dimension(:,:,:,:),pointer :: xl_sum => null() 
     real(kind=8),dimension(:,:,:,:),pointer :: xr_sum => null()
     real(kind=8),dimension(:,:,:,:),pointer :: yl_sum => null()
     real(kind=8),dimension(:,:,:,:),pointer :: yr_sum => null()
     real(kind=8),dimension(:,:,:,:),pointer :: zl_sum => null()
     real(kind=8),dimension(:,:,:,:),pointer :: zr_sum => null()
  end type t_flux_boundary
  type(t_flux_boundary),save,dimension(Gidmin:Gidmax) :: Fb
  integer,save,dimension(L:R,0:2,Ev:Od) :: Ics, Ice, Jcs, Jce, Kcs, Kce, Ias, Iae, Jas, Jae, Kas, Kae
  integer,save :: Levelf 
  integer,save :: Levelc 
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: Buff => null()
     real(kind=8),dimension(:,:,:,:),pointer :: Bufw => null()
  end type t_buf
  type(t_buf),save,dimension(0:2) :: FLP
  public :: fluxcorrection, save_flux, reflux_write, reflux_read
contains
  function get_fb_bnd(gid,surfcode) result(ptr)
    integer,intent(IN) :: gid, surfcode
    real(kind=8),dimension(:,:,:,:),pointer :: ptr
    select case( surfcode )
    case (NIL)
       if (.not. associated(Fb(gid)%xl_bnd)) then 
 allocate(Fb(gid)%xl_bnd(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%xl_bnd
    case (NIR)
       if (.not. associated(Fb(gid)%xr_bnd)) then 
 allocate(Fb(gid)%xr_bnd(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%xr_bnd
    case (NJL)
       if (.not. associated(Fb(gid)%yl_bnd)) then 
 allocate(Fb(gid)%yl_bnd(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%yl_bnd
    case (NJR)
       if (.not. associated(Fb(gid)%yr_bnd)) then 
 allocate(Fb(gid)%yr_bnd(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%yr_bnd
    case (NKL)
       if (.not. associated(Fb(gid)%zl_bnd)) then 
 allocate(Fb(gid)%zl_bnd(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%zl_bnd
    case (NKR)
       if (.not. associated(Fb(gid)%zr_bnd)) then 
 allocate(Fb(gid)%zr_bnd(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 endif 
 ptr => Fb(gid)%zr_bnd
    end select
  end function get_fb_bnd
  function get_fb_sum(gid,surfcode) result(ptr)
    integer,intent(IN) :: gid, surfcode
    real(kind=8),dimension(:,:,:,:),pointer :: ptr
    select case( surfcode )
    case (NIL)
       if (.not. associated(Fb(gid)%xl_sum)) then 
 allocate(Fb(gid)%xl_sum(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xl_sum = 0.d0 
 endif 
 ptr => Fb(gid)%xl_sum
    case (NIR)
       if (.not. associated(Fb(gid)%xr_sum)) then 
 allocate(Fb(gid)%xr_sum(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xr_sum = 0.d0 
 endif 
 ptr => Fb(gid)%xr_sum
    case (NJL)
       if (.not. associated(Fb(gid)%yl_sum)) then 
 allocate(Fb(gid)%yl_sum(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yl_sum = 0.d0 
 endif 
 ptr => Fb(gid)%yl_sum
    case (NJR)
       if (.not. associated(Fb(gid)%yr_sum)) then 
 allocate(Fb(gid)%yr_sum(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yr_sum = 0.d0 
 endif 
 ptr => Fb(gid)%yr_sum
    case (NKL)
       if (.not. associated(Fb(gid)%zl_sum)) then 
 allocate(Fb(gid)%zl_sum(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 Fb(gid)%zl_sum = 0.d0 
 endif 
 ptr => Fb(gid)%zl_sum
    case (NKR)
       if (.not. associated(Fb(gid)%zr_sum)) then 
 allocate(Fb(gid)%zr_sum(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 Fb(gid)%zr_sum = 0.d0 
 endif 
 ptr => Fb(gid)%zr_sum
    end select
  end function get_fb_sum
  subroutine fluxcorrection( level )
    integer,intent(IN) :: level
    integer :: n
    if ( level >= LevelMax ) return
    if ( level < 0 ) return
    Levelf = level + 1
    Levelc = level
    call reflux_init
    call com_restrict
    call dealloc_flux_buffer 
  end subroutine fluxcorrection
  subroutine dealloc_flux_buffer
    integer :: gid, n
    do n = Gidmin, GidListMax( Levelf )
       gid = GidList( n , Levelf )
       if (associated(Fb(gid)%xl_sum)) then 
 deallocate(Fb(gid)%xl_sum) 
 nullify(Fb(gid)%xl_sum) 
 end if
       if (associated(Fb(gid)%xr_sum)) then 
 deallocate(Fb(gid)%xr_sum) 
 nullify(Fb(gid)%xr_sum) 
 end if
       if (associated(Fb(gid)%yl_sum)) then 
 deallocate(Fb(gid)%yl_sum) 
 nullify(Fb(gid)%yl_sum) 
 end if
       if (associated(Fb(gid)%yr_sum)) then 
 deallocate(Fb(gid)%yr_sum) 
 nullify(Fb(gid)%yr_sum) 
 end if
       if (associated(Fb(gid)%zl_sum)) then 
 deallocate(Fb(gid)%zl_sum) 
 nullify(Fb(gid)%zl_sum) 
 end if
       if (associated(Fb(gid)%zr_sum)) then 
 deallocate(Fb(gid)%zr_sum) 
 nullify(Fb(gid)%zr_sum) 
 end if
    enddo
    do n = Gidmin, GidListMax( Levelc )
       gid = GidList( n , Levelc )
       if (associated(Fb(gid)%xl_bnd)) then 
 deallocate(Fb(gid)%xl_bnd) 
 nullify(Fb(gid)%xl_bnd) 
 end if
       if (associated(Fb(gid)%xr_bnd)) then 
 deallocate(Fb(gid)%xr_bnd) 
 nullify(Fb(gid)%xr_bnd) 
 end if
       if (associated(Fb(gid)%yl_bnd)) then 
 deallocate(Fb(gid)%yl_bnd) 
 nullify(Fb(gid)%yl_bnd) 
 end if
       if (associated(Fb(gid)%yr_bnd)) then 
 deallocate(Fb(gid)%yr_bnd) 
 nullify(Fb(gid)%yr_bnd) 
 end if
       if (associated(Fb(gid)%zl_bnd)) then 
 deallocate(Fb(gid)%zl_bnd) 
 nullify(Fb(gid)%zl_bnd) 
 end if
       if (associated(Fb(gid)%zr_bnd)) then 
 deallocate(Fb(gid)%zr_bnd) 
 nullify(Fb(gid)%zr_bnd) 
 end if
    end do
  end subroutine dealloc_flux_buffer
  subroutine save_flux( id, f )
    integer,intent(IN) :: id
    real(kind=8),dimension(:,:,:,:,:),pointer :: f
    real(kind=8),dimension(:,:,:,:),pointer :: fbuf
    real(kind=8) :: rstp
    integer :: n
    integer(kind=8) :: dstepc, dstepf
    Levelf = levels(id)
    Levelc = Levelf - 1
    if ( Levelf > LevelMax ) return
    do n = 0, NBOUNDARY - 1
       fbuf => get_fb_bnd( id, n )
       fbuf(:,:,:,:) = f(lbound(fbuf,1):ubound(fbuf,1),lbound(fbuf,2):ubound(fbuf,2),lbound(fbuf,3):ubound(fbuf,3),lbound(fbuf,4):u&
&bound(fbuf,4), decode_direction( n ))
    enddo
    if ( Levelc < 0 ) return
    dstepf = Dstep( Levelf )
    dstepc = Dstep( Levelc )
    rstp = dble(dstepf)/dble(dstepc)
    do n = 0, NBOUNDARY - 1
       fbuf => get_fb_sum( id, n )
       if ( Step(Levelc)-Dstep(Levelc) == Step(Levelf)-Dstep(Levelf) ) &
            fbuf(:,:,:,:) = 0 
       fbuf(:,:,:,:) = fbuf(:,:,:,:) + f(lbound(fbuf,1):ubound(fbuf,1),lbound(fbuf,2):ubound(fbuf,2),lbound(fbuf,3):ubound(fbuf,3),&
&lbound(fbuf,4):ubound(fbuf,4), decode_direction( n ) ) * rstp
    enddo
  end subroutine save_flux
  subroutine com_restrict
    use mpilib
    use packarr
    use grid
    real(kind=8),dimension(:,:,:,:),pointer :: fp
    integer :: n, ndir, rank, gids, gidd, ranks, rankd, pgid, prank, m, lr
    myrank = get_myrank()
    call pkar_reset
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          do m = 0, NBOUNDARY -1 
             gids = GidListNode(n, Levelf, rank)
             ranks = rank
             pgid = ParentGid(gids,ranks)
             prank = ParentRank(gids,ranks)
             ndir = decode_direction( m )
             lr = decode_LR( m )
             gidd = NeighborGid(lr, ndir, pgid, prank)
             rankd = NeighborRank(lr, ndir, pgid, prank)
             if (.not. ( &
                  NeighborGid(lr, ndir, gids, ranks) == Undefi .and. &
                  NeighborGid(lr, ndir, pgid, prank) /= Undefi ) ) cycle
             if ( myrank == ranks ) then 
                fp => get_fb_sum( gids, m )
                call restrict_f( fp, FLP(ndir)%Buff, ndir ) 
             endif
             if ((myrank) == (ranks)) call pkar_push(FLP(ndir)%Buff(lbound(FLP(ndir)%Buff,1),lbound(FLP(ndir)%Buff,2),lbound(FLP(nd&
&ir)%Buff,3),lbound(FLP(ndir)%Buff,4)), size(FLP(ndir)%Buff), kind(FLP(ndir)%Buff), rankd) 
 if ((myrank) == (rankd)) call pkar_recvlen(size(FLP(ndir)%Buff), kind(FLP(ndir)%Buff), ranks)
          enddo
       enddo
    enddo
    call pkar_sendrecv()
    do rank = 0, 400 -1
       do n = Gidmin, GidListNodeMax(Levelf, rank)
          do m = 0, NBOUNDARY -1 
             gids = GidListNode(n, Levelf, rank)
             ranks = rank
             pgid = ParentGid(gids,ranks)
             prank = ParentRank(gids,ranks)
             ndir = decode_direction( m )
             lr = decode_LR( m )
             if (.not. ( &
                  NeighborGid(lr, ndir, gids, ranks) == Undefi .and. &
                  NeighborGid(lr, ndir, pgid, prank) /= Undefi ) ) cycle
             gidd = NeighborGid(lr, ndir, pgid, prank)
             rankd = NeighborRank(lr, ndir, pgid, prank)
             if (myrank /= rankd) cycle
             if ((myrank) == (rankd)) call pkar_pop(FLP(ndir)%Buff(lbound(FLP(ndir)%Buff,1),lbound(FLP(ndir)%Buff,2),lbound(FLP(ndi&
&r)%Buff,3),lbound(FLP(ndir)%Buff,4)), size(FLP(ndir)%Buff), kind(FLP(ndir)%Buff), ranks)
             call modify_u_by_reflux(m)
          end do
       enddo
    enddo
  contains
    subroutine modify_u_by_reflux(m)
      use eos
      integer,intent(IN) :: m
      real(kind=8),dimension(:,:,:,:),pointer :: fc
      real(kind=8),dimension(:,:,:,:),pointer :: u
      real(kind=8),dimension(0:2) :: ds
      real(kind=8) :: dt, dv
      integer :: lr, lrc, ndir, lri, lrj, lrk
      ds = get_ds( Levelc )
      dv = get_dv( Levelc )
      lr = decode_LR(m)
      lrc = swap_LR( lr )
      ndir = decode_direction(m)
      call left_or_right(gids, ranks, lri, lrj, lrk)
      fc => get_fb_bnd( gidd, encode_surf( lrc, ndir ) ) 
      u => get_Up( gidd ) 
      dt = Dtime( Levelc )
      FLP(ndir)%Bufw = u(Ics(lrc,ndir,lri):Ice(lrc,ndir,lri),Jcs(lrc,ndir,lrj):Jce(lrc,ndir,lrj),Kcs(lrc,ndir,lrk):Kce(lrc,ndir,lrk&
&),Mmin:Mmax)
      call conv_u2w(FLP(ndir)%Bufw, dv)
      if ( lr == Left ) then
         FLP(ndir)%Bufw = FLP(ndir)%Bufw - (FLP(ndir)%Buff - fc(Ias(lrc,ndir,lri):Iae(lrc,ndir,lri),Jas(lrc,ndir,lrj):Jae(lrc,ndir,&
&lrj),Kas(lrc,ndir,lrk):Kae(lrc,ndir,lrk),Mmin:Mmax) )*dt*ds(ndir)
      else
         FLP(ndir)%Bufw = FLP(ndir)%Bufw + (FLP(ndir)%Buff - fc(Ias(lrc,ndir,lri):Iae(lrc,ndir,lri),Jas(lrc,ndir,lrj):Jae(lrc,ndir,&
&lrj),Kas(lrc,ndir,lrk):Kae(lrc,ndir,lrk),Mmin:Mmax) )*dt*ds(ndir)
      endif
      call conv_w2u(FLP(ndir)%Bufw, dv)
      u(Ics(lrc,ndir,lri):Ice(lrc,ndir,lri),Jcs(lrc,ndir,lrj):Jce(lrc,ndir,lrj),Kcs(lrc,ndir,lrk):Kce(lrc,ndir,lrk),Mmin:Mmax) = FL&
&P(ndir)%Bufw
    end subroutine modify_u_by_reflux
  end subroutine com_restrict
  subroutine restrict_f( ff, fc, ncrd )
    real(kind=8),dimension(:,:,:,:),pointer :: ff 
    real(kind=8),dimension(:,:,:,:),pointer :: fc 
    integer,intent(IN) :: ncrd 
    integer :: if,jf,kf,ic,jc,kc,ic0,jc0,kc0,if0,jf0,kf0,m
    real(kind=8),dimension(0:2) :: dsf, dsc
    if ( ncrd == 0 .and. &
         ( &
         size(ff,2) /= size(fc,2)*2 .or. &
         size(ff,3) /= size(fc,3)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NX', size(ff), size(fc)
       stop
    endif
    if ( ncrd == 1 .and. &
         ( &
         size(ff,1) /= size(fc,1)*2 .or. &
         size(ff,3) /= size(fc,3)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NY', size(ff), size(fc)
       stop
    endif
    if ( ncrd == 2 .and. &
         ( &
         size(ff,1) /= size(fc,1)*2 .or. &
         size(ff,2) /= size(fc,2)*2 ) ) then
       write(*,*) 'restrict_f: ff and fc are incompatible size NZ', size(ff), size(fc)
       stop
    endif
    dsf = get_ds( Levelf )
    dsc = get_ds( Levelc )
    ic0 = lbound(fc,1)
    jc0 = lbound(fc,2)
    kc0 = lbound(fc,3)
    if0 = lbound(ff,1)
    jf0 = lbound(ff,2)
    kf0 = lbound(ff,3)
    if ( ncrd == 0 ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if, jf, kf, m) + ff(if, jf, kf+1,m) &
                        +ff(if, jf+1,kf, m) + ff(if, jf+1,kf+1,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
    if ( ncrd == 1 ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if, jf, kf, m) + ff(if, jf, kf+1,m) &
                        +ff(if+1,jf, kf, m) + ff(if+1,jf, kf+1,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
    if ( ncrd == 2 ) then
       do m=MMIN,MMAX
          do kc = lbound(fc,3), ubound(fc,3)
             do jc = lbound(fc,2), ubound(fc,2)
                do ic = lbound(fc,1), ubound(fc,1)
                   if = (ic - ic0) * 2 + if0
                   jf = (jc - jc0) * 2 + jf0
                   kf = (kc - kc0) * 2 + kf0
                   fc(ic,jc,kc,m) = &
                        (ff(if, jf, kf, m) + ff(if, jf+1, kf,m) &
                        +ff(if+1,jf, kf, m) + ff(if+1,jf+1, kf,m) )*dsf(ncrd)/dsc(ncrd)
                enddo
             enddo
          enddo
       enddo
    endif
  end subroutine restrict_f
  subroutine reflux_write(unit, gid)
    integer,intent(IN) :: unit, gid
    if (.not. associated(Fb(gid)%xl_bnd)) then 
 allocate(Fb(gid)%xl_bnd(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xl_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%xl_bnd
    if (.not. associated(Fb(gid)%xr_bnd)) then 
 allocate(Fb(gid)%xr_bnd(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xr_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%xr_bnd
    if (.not. associated(Fb(gid)%yl_bnd)) then 
 allocate(Fb(gid)%yl_bnd(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yl_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%yl_bnd
    if (.not. associated(Fb(gid)%yr_bnd)) then 
 allocate(Fb(gid)%yr_bnd(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yr_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%yr_bnd
    if (.not. associated(Fb(gid)%zl_bnd)) then 
 allocate(Fb(gid)%zl_bnd(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 Fb(gid)%zl_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%zl_bnd
    if (.not. associated(Fb(gid)%zr_bnd)) then 
 allocate(Fb(gid)%zr_bnd(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 Fb(gid)%zr_bnd = 0.d0 
 end if
 write(unit) Fb(gid)%zr_bnd
    if (.not. associated(Fb(gid)%xl_sum)) then 
 allocate(Fb(gid)%xl_sum(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xl_sum = 0.d0 
 end if
 write(unit) Fb(gid)%xl_sum
    if (.not. associated(Fb(gid)%xr_sum)) then 
 allocate(Fb(gid)%xr_sum(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%xr_sum = 0.d0 
 end if
 write(unit) Fb(gid)%xr_sum
    if (.not. associated(Fb(gid)%yl_sum)) then 
 allocate(Fb(gid)%yl_sum(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yl_sum = 0.d0 
 end if
 write(unit) Fb(gid)%yl_sum
    if (.not. associated(Fb(gid)%yr_sum)) then 
 allocate(Fb(gid)%yr_sum(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 Fb(gid)%yr_sum = 0.d0 
 end if
 write(unit) Fb(gid)%yr_sum
    if (.not. associated(Fb(gid)%zl_sum)) then 
 allocate(Fb(gid)%zl_sum(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 Fb(gid)%zl_sum = 0.d0 
 end if
 write(unit) Fb(gid)%zl_sum
    if (.not. associated(Fb(gid)%zr_sum)) then 
 allocate(Fb(gid)%zr_sum(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 Fb(gid)%zr_sum = 0.d0 
 end if
 write(unit) Fb(gid)%zr_sum
  end subroutine reflux_write
  subroutine reflux_read(unit, gid, eof )
    integer,intent(IN) :: unit, gid
    integer,intent(OUT) :: eof
    if (.not. associated(Fb(gid)%xl_bnd)) then 
 allocate(Fb(gid)%xl_bnd(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%xl_bnd
    if (.not. associated(Fb(gid)%xr_bnd)) then 
 allocate(Fb(gid)%xr_bnd(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%xr_bnd
    if (.not. associated(Fb(gid)%yl_bnd)) then 
 allocate(Fb(gid)%yl_bnd(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%yl_bnd
    if (.not. associated(Fb(gid)%yr_bnd)) then 
 allocate(Fb(gid)%yr_bnd(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%yr_bnd
    if (.not. associated(Fb(gid)%zl_bnd)) then 
 allocate(Fb(gid)%zl_bnd(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%zl_bnd
    if (.not. associated(Fb(gid)%zr_bnd)) then 
 allocate(Fb(gid)%zr_bnd(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%zr_bnd
    if (.not. associated(Fb(gid)%xl_sum)) then 
 allocate(Fb(gid)%xl_sum(Imin-1:Imin-1,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%xl_sum
    if (.not. associated(Fb(gid)%xr_sum)) then 
 allocate(Fb(gid)%xr_sum(Imax:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%xr_sum
    if (.not. associated(Fb(gid)%yl_sum)) then 
 allocate(Fb(gid)%yl_sum(Imin:Imax,Jmin-1:Jmin-1,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%yl_sum
    if (.not. associated(Fb(gid)%yr_sum)) then 
 allocate(Fb(gid)%yr_sum(Imin:Imax,Jmax:Jmax,Kmin:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%yr_sum
    if (.not. associated(Fb(gid)%zl_sum)) then 
 allocate(Fb(gid)%zl_sum(Imin:Imax,Jmin:Jmax,Kmin-1:Kmin-1,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%zl_sum
    if (.not. associated(Fb(gid)%zr_sum)) then 
 allocate(Fb(gid)%zr_sum(Imin:Imax,Jmin:Jmax,Kmax:Kmax,Mmin:Mmax)) 
 end if 
 read(unit, iostat=eof) Fb(gid)%zr_sum
  end subroutine reflux_read
  subroutine reflux_init
    use io_util, only : print_msg
    integer :: lr, ndir, lri, lrj, lrk
    logical,save :: bool_initialized = .false.
    if ( bool_initialized ) return
    bool_initialized = .true.
    call print_msg( 'initialize reflux' )
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
    Ics(L,0,:) = Imin
    Ice(L,0,:) = Imin
    Ics(R,0,:) = Imax
    Ice(R,0,:) = Imax
    Jcs(L,1,:) = Jmin
    Jce(L,1,:) = Jmin
    Jcs(R,1,:) = Jmax
    Jce(R,1,:) = Jmax
    Kcs(L,2,:) = Kmin
    Kce(L,2,:) = Kmin
    Kcs(R,2,:) = Kmax
    Kce(R,2,:) = Kmax
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
    Ias(L,0,:) = Imin-1
    Iae(L,0,:) = Imin-1
    Ias(R,0,:) = Imax
    Iae(R,0,:) = Imax
    Jas(L,1,:) = Jmin-1
    Jae(L,1,:) = Jmin-1
    Jas(R,1,:) = Jmax
    Jae(R,1,:) = Jmax
    Kas(L,2,:) = Kmin-1
    Kae(L,2,:) = Kmin-1
    Kas(R,2,:) = Kmax
    Kae(R,2,:) = Kmax
    lr = Left
    lri = Ev
    lrj = Ev
    lrk = Ev
    if ( .not. associated(FLP(0)%Buff) ) then
       do ndir = 0, 2
          allocate( FLP(ndir)%Buff( Ias(lr,ndir,lri):Iae(lr,ndir,lri),Jas(lr,ndir,lrj):Jae(lr,ndir,lrj),Kas(lr,ndir,lrk):Kae(lr,ndi&
&r,lrk),Mmin:Mmax ) )
          allocate( FLP(ndir)%Bufw( Ics(lr,ndir,lri):Ice(lr,ndir,lri),Jcs(lr,ndir,lrj):Jce(lr,ndir,lrj),Kcs(lr,ndir,lrk):Kce(lr,ndi&
&r,lrk),Mmin:Mmax ) )
       enddo
    endif
  end subroutine reflux_init
end module reflux
