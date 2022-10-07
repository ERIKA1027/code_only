module analysis
  use mpilib
  use grid
  implicit none
  private
  real(kind=8),save :: RhoMax, RhoMin
  real(kind=8),save,dimension(Lmin:Lmax) :: RhoMaxLev, RhoMinLev
  public :: RhoMax, RhoMin
  public :: analysis_keyparam
contains
  subroutine analysis_keyparam
    use sinkParticle
    use unit
    use modelParameter, only: MP_mu
    use primordial
    integer :: nparticle, dummy, n
    integer,dimension(:),allocatable :: pid
    real(kind=8),dimension(:),allocatable :: pmass,pmdot,pdm_disk
    real(kind=8),dimension(:,:),allocatable :: pr
    call an_rhomax 
    if (level_sync() /= Lmin) return 
if (get_myrank() == 0) print '(A,2((1P1E15.7)))', 'nH_max , nH_min = ', &
         RhoMax*Unit_rho/(cgs_amu*MP_mu), RhoMin*Unit_rho/(cgs_amu*MP_mu)
    nparticle = sp_getNparticle()
    if (get_myrank() == 0) &
         print *, "nparticle =",nparticle
    if (nparticle > 0) then 
       allocate(pid(nparticle), pmass(nparticle), pmdot(nparticle), pr(0:2, nparticle),pdm_disk(nparticle))
       call sp_sinkdata2array(dummy, pmass, pmdot=pmdot, pr=pr, pid=pid, pdm_disk=pdm_disk)
       if (get_myrank() == 0) then
          do n = 1, nparticle
             print '(A,I0,1P5E12.4)', 'pid, pmass, pmdot, pr(:) = ',&
                  pid(n), pmass(n)*Unit_msun, pmdot(n)*Unit_msun/Unit_yr, pr(:,n)*Unit_au
          end do
       end if
    end if
  end subroutine analysis_keyparam
  subroutine an_rhomax
    real(kind=8),dimension(:,:,:),pointer :: rho
    real(kind=8) :: rhmx, rhmn
    integer :: n, level, gid
    myrank = get_myrank()
    do level = Lmin, LevelMax
       rhmx = TINY(rhmx)
       rhmn = HUGE(rhmn)
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level) 
          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
          rho => get_Ucomp( 0, gid )
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
