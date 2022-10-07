#include "config.h"
#define BINARY_GRAVITY
!-----------------------------------------------------------------------
! external force for the rotating frame.
!-----------------------------------------------------------------------
module externalForce
  implicit none
  private
  public :: source_externalForce
contains
  subroutine source_externalForce(w, dt, gid, f)
    use modelParameter
    use grid
#ifdef BINARY_GRAVITY
    use sinkParticle
#endif !BINARY_GRAVITY
    real(kind=DBL_KIND),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: w  ! (OUT)
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: f
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND) :: r, dv, dtdvrho, gx, gy, gz, omega, vx, vy, vz
#ifdef BINARY_GRAVITY
    real(kind=DBL_KIND) :: gxp, gyp, gzp
    real(kind=DBL_KIND),dimension(MX:MZ) :: dr
#endif !BINARY_GRAVITY
    integer :: level, i, j, k
    level = get_level(gid)
    dv = get_dv(level)
    x => get_xp(gid)
    y => get_yp(gid)
    z => get_zp(gid)
    rho => get_Ucomp(MRHO, gid)
!!$    vx =>  get_Ucomp(MVX,  gid)
!!$    vy =>  get_Ucomp(MVY,  gid)
!!$    vz =>  get_Ucomp(MVZ,  gid)
    omega = MP_Omega
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
             ! centrigugal force and corioli force
             vx = (f(i,j,k,MRHO,MX)+f(i-1,j,k,MRHO,MX))/rho(i,j,k)*0.5d0
             vy = (f(i,j,k,MRHO,MY)+f(i,j-1,k,MRHO,MY))/rho(i,j,k)*0.5d0
             vz = (f(i,j,k,MRHO,MZ)+f(i,j,k-1,MRHO,MZ))/rho(i,j,k)*0.5d0
             gx = omega**2 * x(i) + 2.d0 * omega * vy
             gy = omega**2 * y(j) - 2.d0 * omega * vx
             gz = 0.d0
#ifdef BINARY_GRAVITY
             dr = (/x(i)+MP_a1, y(j), z(k)/)
             call sp_gravityOfParticle(dr, MP_ms1, gxp, gyp, gzp)
             gx = gx + gxp
             gy = gy + gyp
             gz = gz + gzp
             dr = (/x(i)-MP_a2, y(j), z(k)/)
             call sp_gravityOfParticle(dr, MP_ms2, gxp, gyp, gzp)
             gx = gx + gxp
             gy = gy + gyp
             gz = gz + gzp
#endif !BINARY_GRAVITY
             dtdvrho = dt*dv*rho(i,j,k)
             w(i,j,k,MVX) = w(i,j,k,MVX) + gx*dtdvrho
             w(i,j,k,MVY) = w(i,j,k,MVY) + gy*dtdvrho
             w(i,j,k,MVZ) = w(i,j,k,MVZ) + gz*dtdvrho
#ifdef MP
             ! Note that a work due to corioli force is canceled out.
             w(i,j,k,MP)  = w(i,j,k,MP) + dtdvrho*( vx*gx+vy*gy+vz*gz )
#endif !MP
          end do
       end do
    end do
  end subroutine source_externalForce
end module externalForce
