#include "config.h"
!-----------------------------------------------------------------------
! subroutine for external forces and source terms
!-----------------------------------------------------------------------
module externalForce
  implicit none
  private
  real(kind=DBL_KIND),save :: Gm_star
  real(kind=DBL_KIND),save :: Rstar
  public :: source_externalForce
contains
  subroutine source_externalForce(w, dt, gid, f)
    use parameter
    use unit
    use modelParameter
    use grid
!    use sinkParticle, only: sp_SinkRadius, sp_SofteningRadius
#ifdef DM_NFW_PROFILE
    use kinzoku, only: gravDMh
#endif
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: w  ! (OUT)
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: f
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: r, dv, egas, egas_cool, ekin, dtdvdrho

    !real(kind=DBL_KIND),intent(IN)  :: x, y, z
    real(kind=DBL_KIND) :: ghx, ghy, ghz
    real(kind=DBL_KIND) :: g_DMh
    real(kind=DBL_KIND) :: radius, radius2, radius3, Gm_star
    real(kind=DBL_KIND) :: hmax, sp_SinkRadius, sp_SofteningRadius, sp_SofteningRadius2 
   ! integer :: level


#ifdef DM_NFW_PROFILE
    real(kind=DBL_KIND) :: ghx, ghy, ghz
#endif
    integer :: level, i, j, k

!!$    return

    level = get_level(gid)
    dv = get_dv(level)
    x => get_xp(gid)
    y => get_yp(gid)
    z => get_zp(gid)
    rho => get_Ucomp(MRHO, gid)
    vx =>  get_Ucomp(MVX,  gid)
    vy =>  get_Ucomp(MVY,  gid)
    vz =>  get_Ucomp(MVZ,  gid)

    !fx_hpi => get_Ucomp(MXPI,gid)
    !fy_hpi => get_Ucomp(MYPI,gid)
    !fz_hpi => get_Ucomp(MZPI,gid)

    ! Cooling time.
    ! You should specify this value in non-dimension.
    ! A dimensional value can be converted into non-dimensional value by using Unit_* 

    ! w-> MV cm s^-1 g                = g cm s^-1
    ! w-> MP cm^3 g cm^-3 cm^2 s^{-2} = g cm^2  s^-2
    ! w-> MXP cm s^-2 g               = g cm s^-2 
    
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
      
!              w(i,j,k,MVX) = w(i,j,k,MVX) + w(i,j,k,MXPI)*dt ![g cm s^-1]
!              w(i,j,k,MVY) = w(i,j,k,MVY) + w(i,j,k,MYPI)*dt
!              w(i,j,k,MVZ) = w(i,j,k,MVZ) + w(i,j,k,MZPI)*dt  

              dtdvdrho = dt*dv*rho(i,j,k)
!#ifdef MP
!              w(i,j,k,MP)  = w(i,j,k,MP) &
!               +dtdvdrho*(vx(i,j,k)*fx_hpi(i,j,k)+vy(i,j,k)*fy_hpi(i,j,k)+vz(i,j,k)*fz_hpi(i,j,k))
              ! g cm^2 s^-2  ! s g cm s^-1 cm s^-3 = g cm^2 s^-2
!#endif

              radius2= x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
              radius = dsqrt(radius2)

!              do level = Lmin, Lmax
!                 hmax = maxval(CellWidth(:,Lmin))/2.d0**(level-Lmin)
!                 if(1.d-1*MP_Bondi_radius >= MP_spRadius_cell*hmax) then
!                   sp_SinkRadius = MP_spRadius_cell*hmax
!                 end if
!              end do

!              sp_SofteningRadius = sp_SinkRadius

!              if(radius <= sp_SinkRadius) then

!              sp_SofteningRadius2 = sp_SofteningRadius*sp_SofteningRadius
!              radius = dsqrt(radius2+sp_SofteningRadius2)              

!              end if

             
              radius3=radius*radius*radius
              Gm_star = MP_Gconst * MP_Mstar
             ! g_DMh = - Gm_star /radius2


             ! ghx = g_DMh*x(i)/radius
             ! ghy = g_DMh*y(j)/radius
             ! ghz = g_DMh*z(k)/radius
              ghx = -Gm_star*x(i)/radius3
              ghy = -Gm_star*y(j)/radius3
              ghz = -Gm_star*z(k)/radius3




!#ifdef DM_NFW_PROFILE
!              call gravDMh(x(i),y(j),z(k),ghx,ghy,ghz)
              w(i,j,k,MVX) = w(i,j,k,MVX) + ghx*dtdvdrho ![cm s^-2* s*g] => ![g cm s^-1] => [noD] 
              w(i,j,k,MVY) = w(i,j,k,MVY) + ghy*dtdvdrho ![cm s^-2* s*g] => ![g cm s^-1] => [noD] 
              w(i,j,k,MVZ) = w(i,j,k,MVZ) + ghz*dtdvdrho ![cm s^-2* s*g] => ![g cm s^-1] => [noD] 

  #ifdef MP
              w(i,j,k,MP)=w(i,j,k,MP)+dtdvdrho*(vx(i,j,k)*ghx+vy(i,j,k)*ghy+vz(i,j,k)*ghz) ! [ g cm^2 s^-2]
  #endif
!#endif

          end do
       end do
    end do
  end subroutine source_externalForce


!#ifdef DM_NFW_PROFILE
!  subroutine gravDMh(x, y, z, gx, gy, gz)
!    real(kind=DBL_KIND),intent(IN)  :: x, y, z
!    real(kind=DBL_KIND),intent(OUT) :: gx, gy, gz
!    real(kind=DBL_KIND) :: g_DMh
!    real(kind=DBL_KIND) :: Mhalo, radius, radius2, x_r, tap, radi_xy
!
!    radius2= x*x+y*y+z*z
!    radius = dsqrt(radius2)
!    cx_r   = MP_Chalo*radius/MP_Rh
!    tap    = dlog(1.d0+cx_r)-cx_r/(1.d0+cx_r)
!
!    Mhalo = 4.d0*Pi*MP_rhobar*MP_delch*MP_Rs**3.d0*tap !included mass [noD]
!
!    g_DMh = -MP_Gconst*Mhalo/radius2 ! gravitational accerelation [cm s^-2] => [noD]
!
!    ! xyz components ----
!    gx = g_DMh*x/radius
!    gy = g_DMh*y/radius
!    gz = g_DMh*z/radius
!    ! -------------------
!
!  end subroutine gravDMh
!#endif

end module externalForce


















