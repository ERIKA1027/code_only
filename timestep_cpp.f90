module timestep
  use grid
  implicit none
  private
  integer,save,private :: STEP_MODE 
  integer,save,private :: CurrentLevel 
  integer,save,private :: CurrentIndex 
  real(kind=8),save :: Ctime 
  real(kind=8),save :: Cdtime 
  integer(kind=8),save :: Cstep 
  integer(kind=8),save :: Cdstep 
  real(kind=8),save :: Ptime 
  real(kind=8),save :: Pdtime 
  integer(kind=8),save :: Pstep 
  integer(kind=8),save :: Pdstep 
  real(kind=8),save,dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax,0:2),target :: F
  real(kind=8),save,dimension(:,:,:,:,:),pointer :: Fp
  integer,save,dimension(Gidmin:Gidmax) :: Wlist, Wnlist, Ulist
  integer,save :: ListMax
  real(kind=8),save,dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: gam
  public :: step_all_grid
contains
  subroutine timestep_init
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
       Wlist(n) = alloc_block()
       Wnlist(n) = alloc_block()
       Ulist(n) = GidList(n, CurrentLevel )
    enddo
    ListMax = GidListMax( CurrentLevel )
  end subroutine timestep_init
  subroutine timestep_finalize
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
       call dealloc_block(Wlist(n))
       call dealloc_block(Wnlist(n))
    enddo
  end subroutine timestep_finalize
  function get_wp( n ) result( w )
    integer,intent(IN) :: n
    real(kind=8),dimension(:,:,:,:),pointer :: w
    w => get_Up(Wlist(n))
  end function get_wp
  function get_wnp( n ) result( wn )
    integer,intent(IN) :: n
    real(kind=8),dimension(:,:,:,:),pointer :: wn
    wn => get_Up(Wnlist(n))
  end function get_wnp
  subroutine step_all_grid( level )
    use reflux
    use rescue
    use eos 
    integer,intent(IN) :: level
    integer :: id
    CurrentLevel = level
    call timestep_init
    call update_timestep
    STEP_MODE = PREDICTOR
    call boundary_cond
    call rescueLev( CurrentLevel ) 
    do CurrentIndex = Gidmin, ListMax 
       globdbg_mygid = Ulist(CurrentIndex)
       call get_gamma(Ulist(CurrentIndex),gam) 
       call backup_u_2order
       call convert_u2w 
    enddo
    do CurrentIndex = Gidmin, ListMax 
       globdbg_mygid = Ulist(CurrentIndex)
       call get_gamma(Ulist(CurrentIndex),gam) 
       call get_flux 
       call w_update(Cdtime/2) 
       call convert_w2u 
    enddo
    call rescueLev( CurrentLevel )
    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1
    STEP_MODE = CORRECTOR
    call boundary_cond
    call rescueLev( CurrentLevel ) 
    do CurrentIndex = Gidmin, ListMax 
       call backup_u_1order
    enddo
    do CurrentIndex = Gidmin, ListMax 
       globdbg_mygid = Ulist(CurrentIndex)
       call get_gamma(Ulist(CurrentIndex),gam) 
       call get_flux 
       call w_update(Cdtime) 
       call convert_w2u 
       Fp => F
       call save_flux( Ulist(CurrentIndex), Fp )
    enddo
    U_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel) + 1
    STEP_MODE = COMPLETE
    call boundary_cond 
    call rescueLev( CurrentLevel )
    call timestep_finalize
  end subroutine step_all_grid
  subroutine update_timestep
    Ctime = Time( CurrentLevel )
    Cdtime = Dtime( CurrentLevel )
    Cstep = Step( CurrentLevel )
    Cdstep = Dstep( CurrentLevel )
    if ( CurrentLevel > 0 ) then
       Ptime = Time( CurrentLevel-1 )
       Pdtime = Dtime( CurrentLevel-1 )
       Pstep = Step( CurrentLevel-1 )
       Pdstep = Dstep( CurrentLevel-1 )
    endif
    Time( CurrentLevel ) = Ctime + Cdtime
    Step( CurrentLevel ) = Cstep + Cdstep
  end subroutine update_timestep
  subroutine backup_u_2order
    real(kind=8),dimension(:,:,:,:),pointer :: u, u2
    u => get_up( Ulist(CurrentIndex ) )
    u2 => get_u2orderp( Ulist(CurrentIndex) )
    u2 = u
    U2_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel)
    U2_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_2order
  subroutine backup_u_1order
    real(kind=8),dimension(:,:,:,:),pointer :: u, u1
    u => get_up( Ulist(CurrentIndex) )
    u1 => get_u1orderp( Ulist(CurrentIndex) )
    u1 = u
    U1_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel)
    U1_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_1order
  subroutine convert_u2w
    use eos
    use grid
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(:,:,:,:),pointer :: w
    u => get_up( Ulist(CurrentIndex) )
    w => get_wp( CurrentIndex )
    call u2w_withgam( u, w, get_dv(CurrentLevel), gam ) 
  end subroutine convert_u2w
  subroutine convert_w2u
    use eos
    use grid
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(:,:,:,:),pointer :: wn
    u => get_up( Ulist(CurrentIndex) )
    wn => get_wnp( CurrentIndex )
    call w2u_withgam( wn, u, get_dv(CurrentLevel), gam ) 
  end subroutine convert_w2u
  subroutine boundary_cond
    use grid_boundary
    use boundary
    integer :: n
    call boundary_grid( CurrentLevel, STEP_MODE )
    do n = Gidmin, ListMax 
       call boundary_u( Ulist(n), STEP_MODE)
    enddo
  end subroutine boundary_cond
  subroutine get_flux
    use eos
    use util, only : util_arroffset
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax,0:2) :: ul, ur
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,0:2) :: pratio
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax) :: f1d, ul1d, ur1d
    real(kind=8),dimension(:,:,:,:),pointer :: u
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: n, io,jo,ko,i2,j2,k2,i,j,k,m
    u => get_up( Ulist(CurrentIndex) )
    ul = Undefd
    ur = Undefd
    do n = 0, 2
       call util_arroffset(n,io,jo,ko)
       i2 = io*2
       j2 = jo*2
       k2 = ko*2
       do m = Mmin, Mmax
          do k = Kmin-ko, Kmax
             do j = Jmin-jo, Jmax
                do i = Imin-io, Imax
                   ul(i,j,k,m,n) = u(i,j,k,m) &
                        + ((max(0.d0,min((u(i,j,k,m)-u(i-io,j-jo,k-ko,m))*sign(1.d0,(u(i+io,j+jo,k+ko,m)-u(i,j,k,m))),abs(u(i+io,j+&
&jo,k+ko,m)-u(i,j,k,m))))*sign(1.d0,(u(i+io,j+jo,k+ko,m)-u(i,j,k,m)))))*0.5d0
                   ur(i,j,k,m,n) = u(i+io,j+jo,k+ko,m) &
                        - ((max(0.d0,min((u(i+i2,j+j2,k+k2,m)-u(i+io,j+jo,k+ko,m))*sign(1.d0,(u(i+io,j+jo,k+ko,m)-u(i,j,k,m))),abs(&
&u(i+io,j+jo,k+ko,m)-u(i,j,k,m))))*sign(1.d0,(u(i+io,j+jo,k+ko,m)-u(i,j,k,m)))))*0.5d0
                enddo
             enddo
          enddo
       enddo
    enddo
    pratio = min(ul(:,:,:,4,:)/ur(:,:,:,4,:),ur(:,:,:,4,:)/ul(:,:,:,4,:))
    do n = 0, 2
       mcycle = cyclecomp( n )
       ul1d = ul(:,:,:,mcycle,n)
       ur1d = ur(:,:,:,mcycle,n)
       call flux(ul1d, ur1d, pratio, f1d, n, gam) 
       F(:,:,:,mcycle,n) = f1d(:,:,:,:)
    enddo
  end subroutine get_flux
  subroutine w_update( dt )
    use grid
    use util, only : util_arroffset
    use externalForce
    use mpilib 
    real(kind=8),intent(IN) :: dt
    real(kind=8),dimension(0:2) :: ds
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(:,:,:,:),pointer :: w
    real(kind=8),dimension(:,:,:,:),pointer :: wn
    integer :: is,ie,js,je,ks,ke, i,j,k,m,n, io,jo,ko
    w => get_wp( CurrentIndex )
    wn => get_wnp( CurrentIndex )
    ds = get_ds( CurrentLevel )
    wn = w
    do n=0,2
       call util_arroffset(n,io,jo,ko)
       do m = Mmin, Mmax
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i = Imin, Imax
                   wn(i,j,k,m) = wn(i,j,k,m) - dt*ds(n)*(F(i,j,k,m,n)-F(i-io,j-jo,k-ko,m,n))
                enddo
             enddo
          enddo
       enddo
    enddo
    Fp => F
    call source_externalForce(wn, dt, Ulist(CurrentIndex), Fp)
  end subroutine w_update
end module timestep
