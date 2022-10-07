
#include "config.h"
!-------------------------------------------------------------------------
! module depending on equation of state. (Barotropic equation of state)
!-------------------------------------------------------------------------
module diffuse
  implicit none
  private
  public :: diffuseComp
contains
  subroutine diffuseComp(ncomp, gid)
    use grid
    use mpilib
    integer,intent(IN) :: gid, ncomp
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: i, j, k, m, ii, jj, kk, ncell
    real(kind=DBL_KIND) :: rhomin, rhmin, zero, rhoave
    logical :: bool_msg, anyIsNotFinite, isNotFinite, isFinite
    rho => get_Ucomp(ncomp,gid)
    myrank = get_myrank()
    if (ChildGid(Left,Left,Left,gid, myrank) == Undefi) then
       bool_msg = .TRUE.
    else
       bool_msg = .FALSE.
    endif
    if ( anyIsNotFinite( rho, size(rho) ) ) then
       if (bool_msg) print *, '*** diffusion process sets-in due to NaN or Inf is detected. (m, gid) =', ncomp, gid
       ! NaNを救済する
       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                if ( isNotFinite(rho(i,j,k)) ) then
                   print *, 'rho is ', rho(i,j,k)
                   ! NaNではない値で埋める
                   rhoave = 0.d0
                   ncell = 0
                   do kk = k-1,k+1
                      do jj = j-1,j+1
                         do ii = i-1,i+1
                            if ( isFinite(rho(ii,jj,kk)) ) then
                               rhoave = rhoave + rho(ii,jj,kk)
                               ncell = ncell + 1
                            endif
                         enddo
                      enddo
                   enddo
                   rho(i,j,k) = rhoave/ncell
                endif
             enddo
          enddo
       enddo
    endif

    ! 負の値を救済する
    zero = abs(TINY(zero))
    do
       rhomin = MINVAL(rho(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
       if (rhomin > zero) return
       if (bool_msg) print *, '*** diffusion process for (m, gid,rhomin) =',ncomp, gid, rhomin

       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                if (rho(i,j,k) <= zero) then
                   ! 正で最小値を埋める
                   rhmin = HUGE(rhmin)
                   do kk = k-1,k+1
                      do jj = j-1,j+1
                         do ii = i-1,i+1
                            if ( rho(ii,jj,kk) > zero ) then
                               rhmin = min(rhmin, rho(ii,jj,kk))
                            endif
                         enddo
                      enddo
                   enddo
                   rho(i,j,k) = rhmin
                   if (bool_msg) print *, 'rhmin',rhmin
                endif
             enddo
          enddo
       enddo
    enddo

  end subroutine diffuseComp
end module diffuse
