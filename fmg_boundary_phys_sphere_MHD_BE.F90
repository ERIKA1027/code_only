#include "config.h"
! #define DEBUG
!-------------------------------------------------------------------------
! Module for physical boundary condition of multigrid iteration.
! Spherical boundary condition.
!-------------------------------------------------------------------------
module fmg_boundary_phys
  implicit none
  private
  public :: fmg_boundary_u, fmg_boundary_physical_grid, vmg_boundary_u, mg_boundary_u
contains
  ! ----------------------------------------------------------------
  ! ¶­³¦¾ò·ï (AMR FMG cycle)
  ! ----------------------------------------------------------------
  subroutine fmg_boundary_u(fmglev,ju)
    use fmg_data
    integer,intent(IN) :: fmglev, ju
    integer gid, amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call vmg_boundary_u(amrlev, fmglev,ju)
    end do
  end subroutine fmg_boundary_u
  ! ----------------------------------------------------------------
  ! ¶­³¦¾ò·ï (AMR FAS V cycle)
  ! ----------------------------------------------------------------
  subroutine vmg_boundary_u(amrlev,fmglev,ju)
    use fmg_data
    integer,intent(IN) :: amrlev, fmglev, ju
  end subroutine vmg_boundary_u
  ! ----------------------------------------------------------------
  ! ¥Ý¥Æ¥ó¥·¥ã¥ë»Äº¹¤Î¶­³¦¾ò·ï (Base Grid FMG cycle)
  ! ----------------------------------------------------------------
  subroutine mg_boundary_u(mglev, ju)
    use mg_data
    integer,intent(IN) :: mglev, ju
  end subroutine mg_boundary_u
  ! ----------------------------------------------------------------
  ! ÊªÍý¶­³¦¤ËÀÜ¤¹¤ë¥°¥ê¥Ã¥É¤Ë¶­³¦¾ò·ï¤òÍ¿¤¨¤ë
  ! ¶­³¦¾ò·ï¤ÏºÇÁÆ¥°¥ê¥Ã¥É¤Î Â¿½Å¶ËÅ¸³«¤Ë¤è¤ê¡¢µá¤á¤ë¡£
  ! ----------------------------------------------------------------
  subroutine fmg_boundary_physical_grid(ju, jrho)
    use fmg_data
    use grid, only : get_Xp, get_Yp, get_Zp, GidList, GidListMax, get_dv, &
         Imin, Jmin, Kmin, Imax, Jmax, Kmax
    use boundary, only : Rcloud
    integer,intent(IN) :: ju, jrho
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, u
    real(kind=DBL_KIND),dimension(16) :: abuf, abufd
    integer :: i,j,k, amrgid, amrlev, fmglev, gid, la
    real(kind=DBL_KIND) :: pi,ub,rad2, &
         a00,a10,a11r,a11i,a20,a21r,a21i,a22r,a22i,a30,a31r,a31i, &
         a32r,a32i,a33r,a33i, &
         r2,r3,r,rc,cost,sint,cosp,sinp,ri,r2i,r3i,r4i,rci,dv,rhodv
    if (FMG_PDE_TYPE /= FMG_PDE_TYPE_POISSON_EQUATION) then
       print *, '**** error in fmg_boundary_phys_sphere_MHD_BE. invalid FMG_PDE_TYPE', FMG_PDE_TYPE
    endif
    ! --------------
    ! define levels
    ! --------------
    fmglev = FMG_LevelMin       ! real grid only
!!$    print *, 'Rcloud', Rcloud
    ! --------
    ! µå¥Þ¥¹¥¯
    ! --------
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = lbound(GidList,1), GidListMax( amrlev ) ! FMG gid
          amrgid = GidList(gid, amrlev)                 ! AMR gid
          call fmg_arrp(amrlev, fmglev, gid, jrho, rho)
          x => get_Xp(amrgid)
          y => get_Yp(amrgid)
          z => get_Zp(amrgid)
          do k=lbound(rho,3),ubound(rho,3)
             do j=lbound(rho,2),ubound(rho,2)
                do i=lbound(rho,1),ubound(rho,1)
                   if (x(i)**2+y(j)**2+z(k)**2 .ge. Rcloud**2) then
                      rho(i,j,k) = 0.d0
                   endif
                end do
             enddo
          enddo
       enddo
    enddo
!!$#endif
    ! --------
    ! Å¸³«·¸¿ô
    ! --------
    amrlev = AMR_LevelMin       ! coarsest grid only
    pi = 4.d0*atan(1.d0)
    dv = get_dv(amrlev)
    a00 = 0
    a10 = 0
    a11r= 0
    a11i= 0
    a20 = 0
    a21r= 0
    a21i= 0
    a22r= 0
    a22i= 0
    a30 = 0
    a31r= 0
    a31i= 0
    a32r= 0
    a32i= 0
    a33r= 0
    a33i= 0
    do gid = lbound(GidList,1), GidListMax( amrlev ) ! FMG gid
       amrgid = GidList(gid, amrlev)                 ! AMR gid
       call fmg_arrp(amrlev, fmglev, gid, jrho, rho)
       x => get_Xp(amrgid)
       y => get_Yp(amrgid)
       z => get_Zp(amrgid)
       do k=Kmin,Kmax
          do j=Jmin,Jmax
             do i=Imin,Imax
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                r = sqrt(r2)
                r3 = r2*r
                rc = sqrt(x(i)**2 + y(j)**2)
                ri = 1/r
                rci = 1/rc
                cost = z(k)*ri
                sint = rc*ri
                cosp = x(i)*rci
                sinp = y(j)*rci
                rhodv = rho(i,j,k)*dv
                a00 = a00 +  rhodv
                a10 = a10 +  rhodv * r*cost
                a11r= a11r + rhodv * r*sint*cosp
                a11i= a11i + rhodv * r*sint*sinp
                a20 = a20  + rhodv * r2*(3*cost**2-1)
                a21r= a21r + rhodv * r2*sint*cost*cosp
                a21i= a21i + rhodv * r2*sint*cost*sinp
                a22r= a22r + rhodv * r2*sint**2*(cosp**2-sinp**2)
                a22i= a22i + rhodv * r2*sint**2*2*sinp*cosp
                a30 = a30  + rhodv * r3*cost*(5*cost**2-3)
                a31r= a31r + rhodv * r3*sint*(5*cost**2-1)*cosp
                a31i= a31i + rhodv * r3*sint*(5*cost**2-1)*sinp
                a32r= a32r + rhodv * r3*cost*sint**2*(cosp**2-sinp**2)
                a32i= a32i + rhodv * r3*cost*sint**2*2*sinp*cosp
                a33r= a33r + rhodv * r3*sint**3*cosp*(cosp**2-3*sinp**2)
                a33i= a33i + rhodv * r3*sint**3*sinp*(3*cosp**2-sinp**2)
             enddo
          enddo
       enddo
    enddo
    abuf =(/ &
         a00 , &
         a10 , &
         a11r, &
         a11i, &
         a20 , &
         a21r, &
         a21i, &
         a22r, &
         a22i, &
         a30 , &
         a31r, &
         a31i, &
         a32r, &
         a32i, &
         a33r, &
         a33i  &
         /)
    call mpi_allreduce(abuf, abufd, size(abuf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    a00  = abufd(1)
    a10  = abufd(2)
    a11r = abufd(3)
    a11i = abufd(4)
    a20  = abufd(5)
    a21r = abufd(6)
    a21i = abufd(7)
    a22r = abufd(8)
    a22i = abufd(9)
    a30  = abufd(10)
    a31r = abufd(11)
    a31i = abufd(12)
    a32r = abufd(13)
    a32i = abufd(14)
    a33r = abufd(15)
    a33i = abufd(16)

    a00 = -a00/(4.d0 * pi)
    a10 = -a10/(4.d0 * pi)
    a11r= -a11r/(8.d0 * pi)
    a11i= -a11i/(8.d0 * pi)
    a20 = -a20/(16.d0 * pi)
    a21r= -a21r/(8.d0 * pi)*3.d0
    a21i= -a21i/(8.d0 * pi)*3.d0
    a22r= -a22r/(32.d0 * pi)*3.d0
    a22i= -a22i/(32.d0 * pi)*3.d0
    a30 = -a30/(16.d0 * pi)
    a31r= -a31r/(64.d0 * pi)*3.d0
    a31i= -a31i/(64.d0 * pi)*3.d0
    a32r= -a32r/(32.d0 * pi)*15.d0
    a32i= -a32i/(32.d0 * pi)*15.d0
    a33r= -a33r/(64.d0 * pi)*5.d0
    a33i= -a33i/(64.d0 * pi)*5.d0

#ifdef DEBUG
    if ( myrank == 0 ) then
       write(*,*) 'a00 =', a00
       write(*,*) 'a10 =', a10
       write(*,*) 'a11r=', a11r
       write(*,*) 'a11i=', a11i
       write(*,*) 'a20 =', a20
       write(*,*) 'a21r=', a21r
       write(*,*) 'a21i=', a21i
       write(*,*) 'a22r=', a22r
       write(*,*) 'a22i=', a22i
       write(*,*) 'a30 =', a30
       write(*,*) 'a31r=', a31r
       write(*,*) 'a31i=', a31i
       write(*,*) 'a32r=', a32r
       write(*,*) 'a32i=', a32i
       write(*,*) 'a33r=', a33r
       write(*,*) 'a33i=', a33i
    endif
#endif
    ! ---------------------
    ! m=¡ÞL ¤Î¾ì¹ç¤Ï2ÇÜ¤¹¤ë
    ! ---------------------
#define U00 ri*a00
#define U10 r2i*cost*a10
#define U11 r2i*sint*(cosp*a11r + sinp*a11i) *2
#define U20 r3i*(3*cost**2-1)*a20
#define U21 r3i*sint*cost*(cosp*a21r + sinp*a21i) *2
#define U22 r3i*sint**2*((cosp**2-sinp**2)*a22r+2*sinp*cosp*a22i) *2
#define U30 r4i*cost*(5*cost**2-3)*a30
#define U31 r4i*sint*(5*cost**2-1)*(cosp*a31r+sinp*a31i) *2
#define U32 r4i*cost*sint**2*((cosp**2-sinp**2)*a32r+2*sinp*cosp*a32i) *2
#define U33 r4i*sint**3*(cosp*(cosp**2-3*sinp**2)*a33r+sinp*(3*cosp**2-sinp**2)*a33i) *2

    ! --------------------------
    ! ¥°¥ê¥Ã¥É¤Ë¶­³¦¾ò·ï¤òÍ¿¤¨¤ë
    ! --------------------------
#define UB_(I_,IVAL_,UBV_) \
    I_=IVAL_         ;\
    r2i = 1/(x(i)**2 + y(j)**2 + z(k)**2) ;\
    ri = sqrt(r2i)  ;\
    r3i = r2i*ri    ;\
    r4i = r2i*r2i   ;\
    rc = sqrt(x(i)**2 + y(j)**2) ;\
    rci = 1/rc      ;\
    cost = z(k)*ri  ;\
    sint = rc*ri    ;\
    cosp = x(i)*rci ;\
    sinp = y(j)*rci ;\
    UBV_ = U00 + U10 + U11 + U20 + U21 + U22 + U30 + U31 + U32 + U33
    do la = AMR_LevelMin, AMR_LevelMax  ! for all AMR level
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)  ! FMG gid
          amrgid = GidList(gid, la)                     ! AMR gid
          call fmg_arrp(la, fmglev, gid, ju, u)
          x => get_Xp(amrgid)
          y => get_Yp(amrgid)
          z => get_Zp(amrgid)
          if (Geom(la)%Block(gid)%TouchBoundary(Left,MX)) then
             do k=lbound(u,3),ubound(u,3)
                do j=lbound(u,2),ubound(u,2)
                   UB_(i,Imin-1,ub)
                   u(Imin-1,j,k) = ub
                enddo
             enddo
          endif
          if (Geom(la)%Block(gid)%TouchBoundary(Right,MX)) then
             do k=lbound(u,3),ubound(u,3)
                do j=lbound(u,2),ubound(u,2)
                   UB_(i,Imax+1,ub)
                   u(Imax+1,j,k) = ub
                enddo
             enddo
          endif

          if (Geom(la)%Block(gid)%TouchBoundary(Left,MY)) then
             do k=lbound(u,3),ubound(u,3)
                do i=lbound(u,1),ubound(u,1)
                   UB_(j,Jmin-1,ub)
                   u(i,Jmin-1,k) = ub
                enddo
             enddo
          endif
          if (Geom(la)%Block(gid)%TouchBoundary(Right,MY)) then
             do k=lbound(u,3),ubound(u,3)
                do i=lbound(u,1),ubound(u,1)
                   UB_(j,Jmax+1,ub)
                   u(i,Jmax+1,k) = ub
                enddo
             enddo
          endif

          if (Geom(la)%Block(gid)%TouchBoundary(Left,MZ)) then
             do j=lbound(u,2),ubound(u,2)
                do i=lbound(u,1),ubound(u,1)
                   UB_(k,Kmin-1,ub)
                   u(i,j,Kmin-1) = ub
                enddo
             enddo
          endif
          if (Geom(la)%Block(gid)%TouchBoundary(Right,MZ)) then
             do j=lbound(u,2),ubound(u,2)
                do i=lbound(u,1),ubound(u,1)
                   UB_(k,Kmax+1,ub)
                   u(i,j,Kmax+1) = ub
                enddo
             enddo
          end if
       end do
    enddo

  end subroutine fmg_boundary_physical_grid
end module fmg_boundary_phys
