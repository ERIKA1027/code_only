#include "config.h"

program main
  use mpilib
  use primordial,only : ProstFit2, rad_others  

  implicit none
  real(kind=DBL_KIND):: Mass        ! mass in [M_sun]
  real(kind=DBL_KIND):: Mdot        ! mdot in [M_sun/yr]
  real(kind=DBL_KIND):: tage        ! age of star in [yr] 
  real(kind=DBL_KIND):: Radi        ! radius in [R_sun]
  real(kind=DBL_KIND):: Lum        ! luminosity in [L_sun]
  real(kind=DBL_KIND):: Trad       ! Effective temeprature [K]

  real(kind=DBL_KIND) :: tage_pre, dmass
  type(rad_others) :: rad

  call mpilib_init

  Mass = 1.d-2
  tage = 0.d0
  Mdot = 2.d-6
  dmass= 1.d-2

  do while (tage < 1.d7) 
  
    call ProstFit2(Mass, Mdot, tage, Radi, Lum, Trad, rad)

    write(11, *) tage, Mass, Radi, Lum, Trad, rad%xeuv, rad%xfuv, rad%alpha_euv, rad%heat_euv, rad%hhm &
      , rad%lumeuv, rad%lumfuv, rad%sig_euv, rad%sig_fuv, rad%rOII

    if (Mass < 1.d0) then
      tage  = tage + (Mass*dmass)/Mdot
      Mass  = Mass + Mass*dmass
    else
      Mdot = 0.d0
      tage = tage*1.05d0
    endif

  enddo

  !Mass = 100.e0
  !tage = 1.d2
  !Mdot = 0.d0
  !dmass= 1.d-2

  !do while (tage < 1.d7) 
  !
  !  call ProstFit2(Mass, Mdot, tage, Radi, Lum, Trad, rad)

  !  write(11, *) tage, Mass, Radi, Lum, Trad, rad%xeuv, rad%xfuv, rad%alpha_euv, rad%heat_euv, rad%hhm &
  !    , rad%lumeuv, rad%lumfuv, rad%sig_euv, rad%sig_fuv, rad%rOII

  !  tage  = tage*1.05d0
  !  !Mass  = Mass + Mass*dmass


  !enddo


end program main
