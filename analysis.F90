#include "config.h"
#define LOG "log.txt"
!-------------------------------------------------------------------------
! Module for analysis.
! This routine should be called in at the end of subroutine
! step_all_level of multi_timestep.F90
! -------------------------------------------------------------------------
module analysis
  implicit none
  integer,parameter :: IntervalStep = 10
  integer,parameter,private :: CHARLEN = 1000 ! length of character
  private
  public :: analysis_keyparam
contains
  !-----------------------------------------------------------------------
  ! get key parameters
  !-----------------------------------------------------------------------
#ifdef EMULATE_2DIM
#define SZ Imin:Imax,Jmin:Jmax,Kmin:Kmin
#else  !EMULATE_2DIM
#define SZ ARRAYSIZE_IJK
#endif !EMULATE_2DIM
  subroutine analysis_keyparam
    use string
    use memstat
    use grid
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, bx, by, bz, p, db, dbp
    real(kind=DBL_KIND) :: dv, ek, et, em, L1, expt, s_ek, s_et, s_em, s_L1
    real(kind=DBL_KIND),dimension(SZ) :: sw
    integer :: n, level, gid
    real(kind=DBL_KIND) :: ds(MX:MZ)
    character(len=CHARLEN) :: fn

    call memstat_print

    if (.not. bool_analysis()) return ! skip

!!$    ek = 0.d0
!!$    et = 0.d0
!!$    em = 0.d0
!!$    L1 = 0.d0
!!$    myrank = get_myrank()
!!$    do level = Lmin, LevelMax
!!$       dv = get_dv( level )
!!$#ifdef EMULATE_2DIM
!!$       ds = get_ds( level )
!!$       dv = ds(MZ)
!!$#endif !EMULATE_2DIM
!!$       do n = Gidmin, GidListMax( level )
!!$          gid = GidList(n, level) ! gid for U
!!$          if ( ChildGid(Left,Left,Left,gid, myrank) /= Undefi ) cycle
!!$          rho => get_Ucomp( MRHO, gid )
!!$          vx => get_Ucomp( MVX, gid )
!!$          vy => get_Ucomp( MVY, gid )
!!$          vz => get_Ucomp( MVZ, gid )
!!$          bx => get_Ucomp( MBX, gid )
!!$          by => get_Ucomp( MBY, gid )
!!$          bz => get_Ucomp( MBZ, gid )
!!$          p => get_Ucomp( MP, gid )
!!$          db => get_Ucomp( MDB, gid )
!!$          dbp => get_U2comp( MDB, gid )
!!$
!!$          sw = rho(SZ)*(vx(SZ)**2+vy(SZ)**2+vz(SZ)**2)
!!$          ek = ek + SUM(sw) * dv
!!$          sw = p(SZ)/(Gamma-1)
!!$          et = et + SUM(sw) * dv
!!$          sw = (bx(SZ)**2+by(SZ)**2+bz(SZ)**2)
!!$          em = em + SUM(sw) * pi8i * dv
!!$          if ( Step(Lmin) == 0 ) then
!!$             L1 = 0.d0
!!$          else
!!$             expt = exp(Dtime(level)*Ch/Cr)
!!$             sw = abs(dbp(SZ) - db(SZ) * expt)
!!$             L1 = L1 + SUM(sw) /(Ch**2*Dtime(level)) * dv
!!$          endif
!!$       enddo
!!$    enddo
!!$    call mpi_allreduce(ek, s_ek, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$    call mpi_allreduce(et, s_et, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$    call mpi_allreduce(em, s_em, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$    call mpi_allreduce(L1, s_L1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$    if ( get_myrank() == PRIMARY_RANK ) then
!!$#define FH 11
!!$#define LOGFILE
!!$       fn = get_logfile()
!!$       open(FH, file=fn, position='append')
!!$       write(FH, '(I14, 5(1PE17.8))') Step(Lmin), Time(Lmin), s_ek, s_et, s_em, s_L1
!!$       close(FH)
!!$    endif
  end subroutine analysis_keyparam
  !-----------------------------------------------------------------------
  ! return ture when analysis shoudl be executed.
  !-----------------------------------------------------------------------
  function bool_analysis() result(bool)
    use grid, only : level_sync, Step, Lmin
    logical :: bool
    bool = .FALSE.
    if ( level_sync() == Lmin .and. mod(Step(Lmin), IntervalStep) == 0 ) bool = .TRUE.
  end function bool_analysis
  !-----------------------------------------------------------------------
  ! make file name for dumpfile
  !-------------------------------------------------------------------------
  function get_logfile() result(fn)
    use io_util, only : read_env
    use string
    integer :: istat
    character(len=CHARLEN) :: fn, dir
    call read_env('DIR', dir)
    fn = concat(dir,LOG)
  end function get_logfile
end module analysis
