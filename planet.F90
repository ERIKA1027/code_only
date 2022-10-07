! last updated : 2017/06/06 22:04:44
!
! Get positions of planets.
! 
!
! HAE_JD2000 -> HCI
! References: https://ssd.jpl.nasa.gov/?planet_pos
! 
! 
#include "config.h"
#include "planet.h"

! #define USE_MAIN_PROGRAM


module planet
  implicit none
  private
  logical,save :: Initialized = .FALSE.

  integer,parameter :: LEN_PLANETNAME = 10 ! character lenght of a planet name

  real(kind=DBL_KIND),save :: Time ! in number of centuries past J2000.0

  type t_keplerianElements      ! all the components should be in au and radian
     real(kind=DBL_KIND) :: a   ! semi-major axis in au (a)
     real(kind=DBL_KIND) :: ec  ! eccentircity (e)
     real(kind=DBL_KIND) :: in  ! inclination (i)
     real(kind=DBL_KIND) :: lm  ! mean longitude (lambda or L)
     real(kind=DBL_KIND) :: pi  ! logitude of perihelion (pi)
     real(kind=DBL_KIND) :: om  ! longitude of the ascending node (Omega)
  end type t_keplerianElements
  type t_planet
     ! ke0 .... Keplerian elements at JD2000.0
     ! ke1 .... change rate per century
     character(len=LEN_PLANETNAME) :: name
     type(t_keplerianElements) :: ke0
     type(t_keplerianElements) :: ke1
     real(kind=DBL_KIND) :: x, y, z     ! position
  end type t_planet
  type(t_planet) :: LIST_OF_PLANETS

  real(kind=DBL_KIND) :: DAYS_IN_CENTURY = 36525.d0

  public :: LIST_OF_PLANETS, t_planet
  public :: &
       planet_set_pos, &
       planet_get_pos, &
       planet_set_time  

contains
  !-------------------------------------------------------------------------
  ! set planet
  !-------------------------------------------------------------------------
  subroutine planet_set_pos(planet)
    use parameter, only : Pi, Pi2
    type(t_planet),intent(INOUT) :: planet
    type(t_keplerianElements) :: ke
    real(kind=DBL_KIND) :: w, ma, ea
    real(kind=DBL_KIND) :: xhc, yhc, zhc
    real(kind=DBL_KIND) :: xec, yec, zec
    real(kind=DBL_KIND) :: xhg, yhg, zhg
    real(kind=DBL_KIND) :: sinw, cosw
    real(kind=DBL_KIND) :: sino, coso
    real(kind=DBL_KIND) :: sini, cosi
    real(kind=DBL_KIND) :: inc0, om0
    call planet_init
!!$    print *, planet%name

    ke = get_keplerian_elements(planet % ke0, planet % ke1)
    w = ke%pi - ke%om               ! argument of perihelion (w)
    ma = ke%lm - ke%pi              ! mean anomaly
    ma = modulo(ma + Pi, Pi2) - Pi  ! ma in [-pi, pi]
    ea = solve_kepler_eq(ma, ke%ec) ! eccentric anomaly
    ! heliocentric coordinates in the orbintal plane
    xhc = ke%a * (cos(ea) - ke%ec)
    yhc = ke%a * sqrt(1.d0 - ke%ec**2) * sin(ea)
    zhc = 0.d0
    ! ecliptic plane of J2000
    cosw = cos(w);     sinw = sin(w)
    coso = cos(ke%om); sino = sin(ke%om)
    cosi = cos(ke%in); sini = sin(ke%in)
    xec = (cosw*coso - sinw*sino*cosi)*xhc + (-sinw*coso - cosw*sino*cosi)*yhc
    yec = (cosw*sino + sinw*coso*cosi)*xhc + (-sinw*sino + cosw*coso*cosi)*yhc
    zec = (sinw*sini)*xhc + (cosw*sini)*yhc
!!$    planet%x=xec
!!$    planet%y=yec
!!$    planet%z=zec
!!$    return

    ! transform to HCI coordinates (Franz & Harper 2002)
    inc0 = 7.25d0*Pi/180.d0     ! eq (14)
    om0 = 75.76d0*Pi/180.d0     ! eq (14)
    cosi = cos(inc0); sini = sin(inc0)
    coso = cos(om0);  sino = sin(om0)
    ! Rot(om0, inc0, 0) (Franz & Harper 2002, p221)
    xhg =  coso * xec + sino * yec
    yhg = -sino * cosi * xec + coso * cosi * yec + sini * zec
    zhg =  sino * sini * xec - coso * sini * yec + cosi * zec

    planet%x=xhg
    planet%y=yhg
    planet%z=zhg
    
  end subroutine planet_set_pos
  !-------------------------------------------------------------------------
  ! solve the Keplerian equation,  ma = ea - ec sin(ea), by Newton method.
  !-------------------------------------------------------------------------
  function solve_kepler_eq( ma, ec ) result(ea)
    real(kind=DBL_KIND),intent(IN) :: ma, ec
    real(kind=DBL_KIND) :: ea   ! eccentric anomaly
    real(kind=DBL_KIND) :: dm, de
    real(kind=DBL_KIND),parameter :: DEMIN = 1.e-7
    integer :: i
    integer,parameter :: ITRMAX = 10
    ea = ma + ec * sin(ma)
    do i = 1, ITRMAX
       dm = ma - ( ea - ec * sin(ea))
       de = dm/(1.d0 - ec * cos(ea))
       ea = ea + de
       if (abs(de) < DEMIN) exit
    enddo
    if (abs(de) >= DEMIN) print *, '*** ERROR in solve kepelr eq.'
  end function solve_kepler_eq
  !-------------------------------------------------------------------------
  ! get keplarian elements due to the time
  !-------------------------------------------------------------------------
  function get_keplerian_elements(ke0, ke1) result(ke)
    type(t_keplerianElements),intent(IN) :: ke0, ke1
    type(t_keplerianElements) :: ke ! returned value
#define GET_KE(ELEM) ke% ##ELEM## = get_keplerian_element( ke0% ##ELEM## , ke1% ##ELEM## )
    GET_KE(a)
    GET_KE(ec)
    GET_KE(in)
    GET_KE(lm)
    GET_KE(pi)
    GET_KE(om)
#undef GET_KE
  end function get_keplerian_elements
  !-------------------------------------------------------------------------
  function get_keplerian_element(ke0, ke1) result(ke)
    real(kind=DBL_KIND),intent(IN) :: ke0, ke1
    real(kind=DBL_KIND) :: ke   ! returned value
    ke = ke0 + ke1 * Time 
  end function get_keplerian_element
  !-------------------------------------------------------------------------
  ! set time (a module variable)
  ! time inputed should be in JD2000
  !-------------------------------------------------------------------------
  subroutine planet_set_time(time_jd2000)
    real(kind=DBL_KIND),intent(IN) :: time_jd2000
    Time = time_jd2000 / DAYS_IN_CENTURY
  end subroutine planet_set_time
  !-------------------------------------------------------------------------
  ! Read keplerian elements
  !-------------------------------------------------------------------------
  subroutine read_kepler_elements
    use parameter, only : Pi
    use string, only : CHARLEN
    use io_util, only : readenv
    character(len=*),parameter :: FILENAME_DEFAULT = 'p_elem_t1.txt'
    character(len=CHARLEN) :: filename
    integer,parameter :: LUN = 11
    character(len=81) :: linebuf
    character(len=LEN_PLANETNAME) :: planetName
    real(kind=DBL_KIND),dimension(6) :: keList

    if (.not. readenv('PLANET_TABLE', filename)) filename = FILENAME_DEFAULT
    print *, 'reading ' // trim(filename)
    open(LUN, file=filename)
    do
       read(LUN, "(A)", end=999) linebuf
       if (linebuf(1:3) == '---') exit
    end do
    do
#define DEFKE(PLANET, KE) \
      PLANET%KE% a   = keList(1) ;\
      PLANET%KE% ec  = keList(2) ;\
      PLANET%KE% in  = keList(3) * Pi / 180.d0 ;\
      PLANET%KE% lm  = keList(4) * Pi / 180.d0 ;\
      PLANET%KE% pi  = keList(5) * Pi / 180.d0 ;\
      PLANET%KE% om  = keList(6) * Pi / 180.d0

#define READKE(PLANET) \
      read(LUN, *, end=999) planetName, keList ;\
      PLANET%name = planetName ;\
      DEFKE(PLANET, ke0) ;\
      read(LUN, *, end=999) keList ;\
      DEFKE(PLANET, ke1) ;\
      if ( #PLANET /= trim(planetName) ) print *, '*** ERROR'

      READKE(Mercury)
      READKE(Venus)
      READKE(Earth)
      READKE(Mars)
      READKE(Jupiter)
      READKE(Saturn)
      READKE(Uranus)
      READKE(Neptune)
      READKE(Pluto)
   enddo
999 continue
    close(LUN)

  end subroutine read_kepler_elements
  !-------------------------------------------------------------------------
  ! 
  !-------------------------------------------------------------------------
  subroutine planet_get_pos(planet, x, y, z)
    type(t_planet),intent(IN) :: planet
    real(kind=DBL_KIND),intent(OUT) :: x, y, z
    x = planet%x
    y = planet%y
    z = planet%z
  end subroutine planet_get_pos
  !-------------------------------------------------------------------------
  ! Initialize this module
  !-------------------------------------------------------------------------
  subroutine planet_init
    use parameter, only : parameter_init
    implicit none
    if (Initialized) return
    Initialized = .TRUE.

    call parameter_init
    call read_kepler_elements

  end subroutine planet_init
end module planet

#ifdef USE_MAIN_PROGRAM
program main
  use dates, only : dates_datetoJd2000
  use planet, only : planet_set_time, planet_set_pos, planet_get_pos
  use planet, only : LIST_OF_PLANETS
  integer :: year, month, day, hour, mnt, sec
  real(kind=DBL_KIND) :: jd2000, x, y, z
  year = 2008
  month = 10
  day = 29
  hour = 0
  mnt = 0
  sec = 0
!!$  year = 2000
!!$  month = 3
!!$  day = 20
!!$  hour = 12
!!$  mnt = 0
!!$  sec = 0
  call dates_datetoJd2000(year, month, day, hour, mnt, sec, jd2000)
  call planet_set_time(jd2000)
  call planet_set_pos(Earth)
  call planet_get_pos(Earth, x, y, z)
  print *, Earth%name
  print *, x, y, z

  call planet_set_pos(Mars)
  call planet_get_pos(Mars, x, y, z)
  print *, Mars%name
  print *, x, y, z

  call planet_set_pos(Neptune)
  call planet_get_pos(Neptune, x, y, z)
  print *, Neptune%name
  print *, x, y, z


end program main
#endif !USE_MAIN_PROGRAM
