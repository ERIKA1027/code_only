module refineCond
  implicit none
  private
  public :: refineCond_eval
contains
  function refineCond_eval( gid ) result( bool )
    integer,intent(IN) :: gid
    logical :: bool
    bool = refineCond_nest( gid ) .or. refineCond_jeans( gid )
  end function refineCond_eval
  function refineCond_nest( gid ) result( bool )
    use grid
    use parameter
    use modelParameter, only : MP_Lmax0, MP_Boxsize
    integer,intent(IN) :: gid
    logical :: bool
    real(kind=8),dimension(0:2) :: h
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: xedgeo, yedgeo, zedgeo
    real(kind=8) :: rmax
    integer :: i,j,k
    bool = .FALSE. 
    if ( gid == Undefi ) return
    if (get_level(gid) >= MP_Lmax0) return 
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    h = CellWidth(:,get_level(gid))
    xedgeo = max(abs(x(Imin)), abs(x(Imax)))
    yedgeo = max(abs(y(Jmin)), abs(y(Jmax)))
    zedgeo = max(abs(z(Kmin)), abs(z(Kmax)))
    rmax = MP_Boxsize / 2**(get_level(gid)+1)
    if ( xedgeo < rmax-0.5*h(0)*(8) .and. & 
         yedgeo < rmax-0.5*h(1)*(8) .and. &
         zedgeo < rmax-0.5*h(2)*(8) ) then
       bool = .TRUE.
    else
       bool = .FALSE.
    end if
  end function refineCond_nest
  function refineCond_jeans( gid ) result( bool )
    use grid
    use parameter, only : Pi
    use modelParameter, only : MP_JeansConst, MP_Gconst
    use unit 
    integer,intent(IN) :: gid
    logical :: bool
    real(kind=8) :: jlength, hmax
    real(kind=8),dimension(0:2) :: h
    real(kind=8),dimension(:,:,:),pointer :: rho
    real(kind=8),dimension(:,:,:),pointer :: p
    real(kind=8) :: flag, csp
    integer :: i,j,k
    bool = .FALSE. 
    if ( gid == Undefi ) return
    rho => get_Ucomp(0,gid)
    p => get_Ucomp(4,gid)
    h = CellWidth(:,get_level(gid))
    hmax = maxval(h)
    jlength = HUGE(jlength)
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             csp = sqrt(p(i,j,k)/rho(i,j,k))
             jlength = min(jlength, csp * sqrt(Pi / (MP_Gconst * rho(i,j,k))))
          enddo
       enddo
    enddo
    if ( jlength < MP_JeansConst * hmax ) bool = .true.
  end function refineCond_jeans
end module refineCond
