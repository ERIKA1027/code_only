module externalForce
  implicit none
  private
  real(kind=8),save :: Gm_star
  real(kind=8),save :: Rstar
  public :: source_externalForce
contains
  subroutine source_externalForce(w, dt, gid, f)
    use parameter
    use unit
    use modelParameter
    use grid
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: w 
    real(kind=8),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=8),dimension(:,:,:,:,:),pointer :: f
    real(kind=8),dimension(:,:,:),pointer :: rho, vx, vy, vz
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: r, dv, egas, egas_cool, ekin, dtdvdrho
    real(kind=8) :: ghx, ghy, ghz
    real(kind=8) :: g_DMh
    real(kind=8) :: radius, radius2, radius3, Gm_star
    real(kind=8) :: hmax, sp_SinkRadius, sp_SofteningRadius, sp_SofteningRadius2
    integer :: level, i, j, k
    level = get_level(gid)
    dv = get_dv(level)
    x => get_xp(gid)
    y => get_yp(gid)
    z => get_zp(gid)
    rho => get_Ucomp(0, gid)
    vx => get_Ucomp(1, gid)
    vy => get_Ucomp(2, gid)
    vz => get_Ucomp(3, gid)
    do k = Kmin, Kmax
       do j = Jmin, Jmax
          do i = Imin, Imax
              dtdvdrho = dt*dv*rho(i,j,k)
              radius2= x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
              radius = dsqrt(radius2)
              radius3=radius*radius*radius
              Gm_star = MP_Gconst * MP_Mstar
              ghx = -Gm_star*x(i)/radius3
              ghy = -Gm_star*y(j)/radius3
              ghz = -Gm_star*z(k)/radius3
              w(i,j,k,1) = w(i,j,k,1) + ghx*dtdvdrho 
              w(i,j,k,2) = w(i,j,k,2) + ghy*dtdvdrho 
              w(i,j,k,3) = w(i,j,k,3) + ghz*dtdvdrho 
              w(i,j,k,4)=w(i,j,k,4)+dtdvdrho*(vx(i,j,k)*ghx+vy(i,j,k)*ghy+vz(i,j,k)*ghz) 
          end do
       end do
    end do
  end subroutine source_externalForce
end module externalForce
