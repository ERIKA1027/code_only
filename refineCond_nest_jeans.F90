#include "config.h"
!-------------------------------------------------------------------------
! This module provides refinement condition due to Jeans condition
!-------------------------------------------------------------------------
module refineCond
  implicit none
  private
  public :: refineCond_eval
contains

  !-------------------------------------------------------------------------
  ! evaluate for refinement condition by given by grid ID
  !-------------------------------------------------------------------------
  function refineCond_eval( gid ) result( bool )
    integer,intent(IN) :: gid
    logical :: bool
    bool = refineCond_nest( gid ) .or. refineCond_jeans( gid )
  end function refineCond_eval
  
  !-------------------------------------------------------------------------
  ! evaluate for refinement condition by given by grid ID
  !-------------------------------------------------------------------------
  function refineCond_nest( gid ) result( bool )
    use grid
    use parameter
    use modelParameter, only : MP_Lmax0, MP_Boxsize
    integer,intent(IN) :: gid
    logical :: bool 
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: xedgeo, yedgeo, zedgeo
!!$    real(kind=DBL_KIND),parameter :: rmax = 5.d0
    real(kind=DBL_KIND) :: rmax
    integer :: i,j,k
    bool = .FALSE.              ! defalut value
    if ( gid == Undefi ) return
    if (get_level(gid) >= MP_Lmax0) return ! MP_Lmax0までnested gridを作成
    
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    h = CellWidth(:,get_level(gid))

    ! outer edges of block (defined with cell center)
    xedgeo = max(abs(x(Imin)), abs(x(Imax)))
    yedgeo = max(abs(y(Jmin)), abs(y(Jmax)))
    zedgeo = max(abs(z(Kmin)), abs(z(Kmax)))

    rmax = MP_Boxsize / 2**(get_level(gid)+1)

    ! if ( xedgeo < rmax-h(MX)*(NI) .and. & ! h*(N) is margin
    !      yedgeo < rmax-h(MY)*(NJ) .and. &
    !      zedgeo < rmax-h(MZ)*(NK) ) then
    if ( xedgeo < rmax-0.5*h(MX)*(NI) .and. & ! h*(N) is margin (等式をTrueにしたいのでちょっと減らす)
         yedgeo < rmax-0.5*h(MY)*(NJ) .and. &
         zedgeo < rmax-0.5*h(MZ)*(NK) ) then        
       bool = .TRUE.
    else
       bool = .FALSE.
    end if
  end function refineCond_nest

  !-------------------------------------------------------------------------
  ! evaluate for refinement condition by given by grid ID
  !-------------------------------------------------------------------------
  function refineCond_jeans( gid ) result( bool )
    use grid
    use parameter, only : Pi
    use modelParameter, only : MP_JeansConst, MP_Gconst
    use unit ! KS DEBUG
    integer,intent(IN) :: gid
    logical :: bool
    real(kind=DBL_KIND) :: jlength, hmax
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: p
    real(kind=DBL_KIND) :: flag, csp
    integer :: i,j,k
    bool = .FALSE.              ! defalut value
    if ( gid == Undefi ) return

    ! Jeans condition
    rho => get_Ucomp(MRHO,gid)
    p => get_Ucomp(MP,gid)

    h = CellWidth(:,get_level(gid))

    hmax = maxval(h)
    jlength = HUGE(jlength)
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             ! isothermal sound speed
             csp = sqrt(p(i,j,k)/rho(i,j,k))
             jlength = min(jlength, csp * sqrt(Pi / (MP_Gconst * rho(i,j,k))))
             ! if ( csp * sqrt(Pi / (MP_Gconst * rho(i,j,k))) <  JeansConst * hmax ) then
             !    print '(/,A,I0,10(1P1E9.2))', "KS DEBUG, j_cond: gid, rho, p, cs, lam_j, hmax: ", &
             !         gid, rho(i,j,k)*Unit_rho,p(i,j,k)*Unit_e, csp*Unit_v, jlength*Unit_au, hmax*Unit_au
             ! end if
          enddo
       enddo
    enddo
    if ( jlength < MP_JeansConst * hmax ) bool = .true.
    
  end function refineCond_jeans
end module refineCond
