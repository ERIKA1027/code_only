#include "config.h"
#include "barotropic.h"
!-------------------------------------------------------------------------
! Module for analysis for using sink particles.
! This routine should be called in at the end of subroutine
! step_all_level of multi_timestep.F90
! -------------------------------------------------------------------------
module analysis
  use mpilib
  use grid
  implicit none
  private
  real(kind=DBL_KIND),save :: RhoMax, RhoMin
  real(kind=DBL_KIND),save,dimension(Lmin:Lmax) :: RhoMaxLev, RhoMinLev
  !
  public :: RhoMax, RhoMin
  public :: analysis_keyparam
contains
  !-----------------------------------------------------------------------
  ! Front end of analysis
  !-----------------------------------------------------------------------
  subroutine analysis_keyparam
    use sinkParticle
    use unit
    use modelParameter, only: MP_mu
#if MODEL_ART > 0
    use primordial
#endif !MODEL_ART

    integer :: nparticle, dummy, n
    integer,dimension(:),allocatable :: pid
    real(kind=DBL_KIND),dimension(:),allocatable :: pmass,pmdot,pdm_disk
    real(kind=DBL_KIND),dimension(:,:),allocatable :: pr 
    

    !get RhoMax and RhoMin
    call an_rhomax              !bool_outputの判定に必要

    if (level_sync() /= Lmin)  return ! 全レベルの更新が終わったとき以外は表示をskip


if (get_myrank() == PRIMARY_RANK) print '(A,2((1P1E15.7)))', 'nH_max , nH_min = ', &
#if MODEL_ART > 0         
         RhoMax*Unit_rho/(cgs_amu*MP_mu), RhoMin*Unit_rho/(cgs_amu*MP_mu)
#else !MODEL_ART
    RhoMax*Unit_rho/cgs_mh, RhoMin*Unit_rho/cgs_mh
#endif !MODEL_ART
    
    !sinkParticle
    nparticle = sp_getNparticle()

    if (get_myrank() == PRIMARY_RANK) &
         print *, "nparticle =",nparticle

    if (nparticle > 0) then     ! if sink particles exist
       ! allocate(pid(nparticle), pmass(nparticle), pmdot(nparticle), pr(MX:MZ, nparticle))
       ! call sp_sinkdata2array(dummy, pmass, pmdot=pmdot, pr=pr, pid=pid)
       !---------------- KS DEBUG---------------!
       allocate(pid(nparticle), pmass(nparticle), pmdot(nparticle), pr(MX:MZ, nparticle),pdm_disk(nparticle))
       call sp_sinkdata2array(dummy, pmass, pmdot=pmdot, pr=pr, pid=pid, pdm_disk=pdm_disk)
       !---------------- KS DEBUG---------------!       
       if (get_myrank() == PRIMARY_RANK) then
          do n = 1, nparticle
             print '(A,I0,1P5E12.4)', 'pid, pmass, pmdot, pr(:) = ',&
                  pid(n), pmass(n)*Unit_msun, pmdot(n)*Unit_msun/Unit_yr, pr(:,n)*Unit_au
          end do
       end if
    end if
          
    

  end subroutine analysis_keyparam
  !-----------------------------------------------------------------------
  ! find rhomax
  !-----------------------------------------------------------------------
  subroutine an_rhomax
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    real(kind=DBL_KIND) :: rhmx, rhmn
    integer :: n, level, gid
    myrank = get_myrank()
    ! find rhomax
    do level = Lmin, LevelMax
       rhmx = TINY(rhmx)
       rhmn = HUGE(rhmn)
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) ! gid for U
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( MRHO, gid )
          rhmx = max(rhmx, MAXVAL(rho, mask=GridMask))
          rhmn = min(rhmn, MINVAL(rho, mask=GridMask))
       enddo
       call mpi_allreduce(rhmx, RhoMaxLev(level), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(rhmn, RhoMinLev(level), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    enddo
    RhoMax = MAXVAL(RhoMaxLev(Lmin:LevelMax))
    RhoMin = MINVAL(RhoMinLev(Lmin:LevelMax))

  end subroutine an_rhomax

  
end module analysis
