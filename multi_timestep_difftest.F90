
#include "config.h"
!-------------------------------------------------------------------------
! Module for multi-timestep
!
!-------------------------------------------------------------------------
module multi_timestep
  implicit none
  private
  integer,save :: CurrentLevel       ! ���ߤΥ�٥�
  public :: multi_timestep_init, step_all_level
contains
  !-----------------------------------------------------------------------
  ! initialize module timestep
  !-----------------------------------------------------------------------
  subroutine multi_timestep_init
  end subroutine multi_timestep_init
  !-----------------------------------------------------------------------
  ! muti timestep �ǻ��ֿ��
  !-----------------------------------------------------------------------
  subroutine step_all_level
    use grid, only : Lmin, Step, level_sync
    use timestep
    use reflux
    use analysis
    use rescue
    use eos, only : Rhocr
    use unit
    use outputdata
#ifdef WITH_SELFGRAVITY
    use fmg
#endif !WITH_SELFGRAVITY
#if defined(SINKPARTICLE) && defined(SINGLE_STEP)
    use sinkParticle
#endif
    do
       call findCurrentLevel
!!$       call grid_refinement

!!$       ! ----------------
!!$       ! CFL condition
!!$       ! ----------------
!!$#ifdef SINGLE_STEP
       call cfl_condition_singlestep
!!$#else
!!$       call cfl_condition
!!$#endif ! SINGLE_STEP
!!$
!!$       ! ------------------------------------------
!!$       ! Multi grid for gravity for multi-timestep
!!$       ! ------------------------------------------
!!$#if defined(WITH_SELFGRAVITY) && !defined(SINGLE_STEP)
!!$       if ( bool_fmg() ) call fmg_multigrid
!!$#endif
!!$
!!$       ! ----------------------
!!$       ! Solve Hydrodynamics
!!$       ! ----------------------
       call step_all_grid( CurrentLevel )
!!$       ! ------------------------------------------
!!$       ! Consistency between fine and coarse grids
!!$       ! ------------------------------------------
!!$       if ( bool_sync() ) then
!!$          call fluxcorrection( CurrentLevel - 1 ) ! for parent grid
!!$          call converge
!!$          call rescueLev( CurrentLevel - 1 )
!!$       endif
!!$       ! ------------------------------------------
!!$       ! Ohmic dissipation
!!$       ! ------------------------------------------
!!$       if ( level_sync() == Lmin ) call fmg_ambipolar_diffusion
       if ( level_sync() == Lmin ) call fmg_ohmic_dissipation
!!$       if ( level_sync() == Lmin ) call fmg_poisson
!!$       ! ------------------------------------------
!!$       ! Multi grid for gravity for single-timestep
!!$       ! ------------------------------------------
!!$#if defined(WITH_SELFGRAVITY) && defined(SINGLE_STEP)
!!$       if ( bool_fmg() ) then
!!$          if (Step(Lmin) == 1) call fmg_multigrid
!!$          call fmg_multigrid
!!$          call source_g_all_level
!!$          call rescueAllLev
!!$       endif
!!$#endif
!!$       ! -------------
!!$       ! sink cell
!!$       ! -------------
!!$#if defined(SINKPARTICLE) && defined(SINGLE_STEP)
!!$       if ( level_sync() == Lmin ) then
!!$          call sp_update
!!$          call converge
!!$       endif
!!$#endif ! SINKPARTICLE
!!$       if ( level_sync() == Lmin ) call analysis_keyparam
       ! -------------------
       ! Output data to file
       ! -------------------
!!$       call output_data

       if ( bool_halt() ) then
!!$          call writeSnap_denseRegion(100./Unit_au)
!!$          call writeSnap_clusters(100./Unit_au)
          return
       endif
    end do
  end subroutine step_all_level
  !-----------------------------------------------------------------------
  ! �����ʤ��٤���٥�򸫤Ĥ��롣CurrentLevel �򥻥å�
  !-----------------------------------------------------------------------
  subroutine findCurrentLevel
    use grid
    integer:: level
    integer(kind=LLONG_KIND) :: nstepmin

    ! find minmum time step
    nstepmin = huge( nstepmin ) ! start value
    do level = Lmin, LevelMax
       nstepmin = min( nstepmin, Step(level))
    enddo
    ! search level which step = nstepmin
    do level = Lmin, LevelMax
       if ( Step(level) == nstepmin ) then
          CurrentLevel = level
          return
       endif
    enddo
  end subroutine findCurrentLevel
  !-----------------------------------------------------------------------
  ! grid refinement using all synchronzied grids
  !-----------------------------------------------------------------------
  subroutine grid_refinement
    use grid
    use io_util
    use refine
    use string
    integer :: lev, slevel, levmax
    logical :: bool

#ifdef NOT_REFINEMENT
    return
#endif !NOT_REFINEMENT

    if ( .not. bool_refinement() ) return
    slevel = level_sync()
    levmax = LevelMax

    do lev = levmax+1, slevel+1, -1 ! �٤����Ƥ�
       call refineLevel(lev, bool)
    enddo

!!$    do lev = slevel+1, Lmax     ! �Ƥ���٤�
!!$       call refineLevel(lev, bool)
!!$       if ( .not. bool ) exit
!!$    enddo


!!$    print *, 'sync lev, this lev', slevel, CurrentLevel
!!$    if ( LevelMax > levmax ) &
!!$         call print_msg( 'new level = '// num2char(LevelMax))
  end subroutine grid_refinement
  !-----------------------------------------------------------------------
  ! �ҥ���åɤ� refinement ���뤿��ξ��(�襹�ƥå�)
  !-----------------------------------------------------------------------
  function bool_refinement_BAK2() result(bool)
    use grid
    logical :: bool
    bool = .true.
  end function bool_refinement_BAK2
#ifdef SINGLE_STEP
  !-----------------------------------------------------------------------
  ! �ҥ���åɤ� refinement ���뤿��ξ��(���ƥ���åɤ������� refinement)
  !-----------------------------------------------------------------------
  function bool_refinement() result(bool)
    use grid
    logical :: bool
    bool = .false.
    if ( level_sync() == Lmin ) bool = .true.
  end function bool_refinement
#else !SINGLE_STEP
  !-----------------------------------------------------------------------
  ! �ҥ���åɤ� refinement ���뤿��ξ��(���Ϲ�����Ʊ�������ߥ�)
  !-----------------------------------------------------------------------
  function bool_refinement() result(bool)
    use grid
    logical :: bool
    integer :: stp
    integer,save :: laststep=-1 ! �Ƕ�refiment�������ƥå�
    bool = .false.
    if ( Step(CurrentLevel) /= laststep ) then
       bool = .TRUE.
       laststep = Step(CurrentLevel)
    endif
  end function bool_refinement
#endif !SINGLE_STEP
  !-----------------------------------------------------------------------
  ! find dtime and dstep with multi-timestep
  !-----------------------------------------------------------------------
  subroutine cfl_condition
    use eos
    use grid
    use io_util
    use string
    use mpilib
    use modelParameter, only : MP_T_Last
#ifdef SINKPARTICLE
    use sinkParticle
#endif ! SINKPARTICLE
    real(kind=DBL_KIND),dimension(Lmin:Lmax) :: t, dt ! �������
    integer(kind=LLONG_KIND),dimension(Lmin:Lmax) :: st, dst ! �������
    integer :: n, gid, l
    real(kind=DBL_KIND) :: dtlocal, thisdt_try, thisdt, dtbuf
    integer(kind=LLONG_KIND) :: thisds_try, syncstep, ds, thisds
    ! ------------------------------
    ! ��ʬ����˻�����(dtlocal)�����
    ! -----------------------------
    dtlocal = huge( dtlocal )
    do n = Gidmin, GidListMax( CurrentLevel )
       gid = GidList( n, CurrentLevel )
       dtlocal = min( dtlocal , get_dt_by_cflcond( gid ) )
    enddo
    call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    dtlocal = dtbuf

#ifdef SINKPARTICLE
    call sp_restrictCFL(dtlocal, CurrentLevel)
#endif ! SINKPARTICLE

    if (dtlocal < 0) then
       write(*,*) 'negative dtlocal', dtlocal
       stop
    endif

    ! ----------------
    ! ����������������
    ! ----------------
    t(Lmin:LevelMax) = Time(Lmin:LevelMax)
    dt(Lmin:LevelMax) = Dtime(Lmin:LevelMax)
    st(Lmin:LevelMax) = Step(Lmin:LevelMax)
    dst(Lmin:LevelMax) = Dstep(Lmin:LevelMax)
    ! --------------------------------------------------
    ! �ƥ�٥��Ĵ�����롣
    ! �ƥ�٥�� step, dstep ���ѹ�����롣
    ! ����٥�� thisdt, thisds ����롣
    ! --------------------------------------------------
    if ( CurrentLevel == Lmin ) then
       thisdt = dtlocal
       thisds = 1
    else
       n=0
       do
          ! ��ԥ��ƥå�
          thisds_try = dst( CurrentLevel - 1 ) / 2**n
          thisdt_try = dt( CurrentLevel - 1 ) / 2**n
          ! ����ʬ��ǽ������
          if (thisds_try == 0) then
             syncstep = st(Lmin)-dst(Lmin) ! �Ǹ�˥������Х�Ʊ���������ƥå�
             do l=Lmin,LevelMax
                ds = st(l)-syncstep  ! �Ǹ�Υ������Х�Ʊ���������ʬ
                st(l) = syncstep + ds*2
                dst(l) = dst(l)*2
             enddo
             cycle
          endif
          ! �ƥ���åɤȤδط���CFL�����������褦��
          ds = st(CurrentLevel)-st(CurrentLevel-1)+dst(CurrentLevel-1) ! �ƤȤ�Ʊ���������ʬ
          if ((mod(ds,thisds_try) == 0) .and. (thisdt_try <= dtlocal)) exit
          n=n+1
       enddo
       thisdt = thisdt_try
       thisds = thisds_try
    endif
    dt( CurrentLevel ) = thisdt
    dst( CurrentLevel ) = thisds

    ! -----------------------
    ! Ʊ����٥�ǻ���򤽤����롣
    ! ����δݤ����������
    ! -----------------------
    t( level_sync() : LevelMax-1 ) = t( LevelMax )

    ! -------------------------
    ! ���礦�ɤλ���˽�λ����
    ! -------------------------
    if (CurrentLevel == Lmin) &
         dt(Lmin) = min(dt(Lmin), MP_T_Last - t(Lmin))
    ! ----------------------------------------------------------
    ! global variable �򥢥åץǡ��Ȥ��� (�ޤ����ֿ�ʤϤ��ʤ�)
    ! ----------------------------------------------------------
    Time(Lmin:LevelMax) = t(Lmin:LevelMax)
    Dtime(Lmin:LevelMax) = dt(Lmin:LevelMax)
    Step(Lmin:LevelMax) = st(Lmin:LevelMax)
    Dstep(Lmin:LevelMax) = dst(Lmin:LevelMax)
    ! -----------------------------------------------
    ! ����ȥ��ƥåפ�ɽ��
    ! ɽ����������ǡ�global variable �Ϲ�������ʤ�
    ! -----------------------------------------------
    if ( get_myrank() == PRIMARY_RANK ) then
       t(CurrentLevel)=t(CurrentLevel)+dt(CurrentLevel)
       st(CurrentLevel)=st(CurrentLevel)+dst(CurrentLevel)
       write(*,*) (st(l),l=Lmin,LevelMax)
       write(*,*) (dst(l),l=Lmin,LevelMax)
       write(*,*) (t(l),l=Lmin,LevelMax)
       write(*,*) (dt(l),l=Lmin,LevelMax)
    endif

!!$    call print_msg( &
!!$         'level = ' // trim( num2char(CurrentLevel) )     // &
!!$         ' t = '    // trim( num2char(t(CurrentLevel)) )  // &
!!$         ' st = '   // trim( num2char(int(st(CurrentLevel))) )  &
!!$         )
  end subroutine cfl_condition
  !-----------------------------------------------------------------------
  ! find dtime and dstep without multi-timestep
  !-----------------------------------------------------------------------
  subroutine cfl_condition_singlestep
    use eos
    use grid
    use io_util
    use string
    use mpilib
    use modelParameter, only : MP_T_Last
#ifdef SINKPARTICLE
    use sinkParticle
#endif ! SINKPARTICLE
    real(kind=DBL_KIND),dimension(Lmin:Lmax) :: t, dt ! �������
    integer(kind=LLONG_KIND),dimension(Lmin:Lmax) :: st, dst ! �������
    integer :: n, gid, l, lev
    real(kind=DBL_KIND) :: dtlocal, dtbuf
    ! ----------------------------------
    ! ����٥�ǻ�����(dtlocal)�����
    ! ----------------------------------
!!$    if ( CurrentLevel == Lmin ) then
!!$       dtlocal = huge( dtlocal )
!!$       do lev = Lmin, LevelMax
!!$          do n = Gidmin, GidListMax( lev )
!!$             gid = GidList( n, lev )
!!$             dtlocal = min( dtlocal , get_dt_by_cflcond( gid ) )
!!$             dtlocal = 1.d0
!!$          enddo
!!$#ifdef SINKPARTICLE
!!$          call sp_restrictCFL(dtlocal, lev)
!!$#endif ! SINKPARTICLE
!!$       enddo
!!$       call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
!!$       dtlocal = dtbuf
!!$       if (dtlocal < 0) then
!!$          write(*,*) 'negative dtlocal', dtlocal
!!$          stop
!!$       endif
!!$    endif
    call read_env('DTIME', dtlocal)

    ! ----------------------
    ! ���ƥå׿�������κ�ʬ
    ! ----------------------
    if ( CurrentLevel == Lmin ) then
       Dstep(CurrentLevel) = 1
       Dtime(CurrentLevel) = dtlocal
    else
       Dstep(CurrentLevel) = Dstep(CurrentLevel-1)
       Dtime(CurrentLevel) = Dtime(CurrentLevel-1)
    endif
    ! -------------------------
    ! ���礦�ɤλ���˽�λ����
    ! -------------------------
    if (CurrentLevel == Lmin) &
         Dtime(Lmin) = min(Dtime(Lmin), MP_T_Last - Time(Lmin))
!!$    ! -----------------------------------------------------
!!$    ! ���������ʥ����ॹ�ƥåפ���(Predictor Step�Τ���)
!!$    ! -----------------------------------------------------
!!$    if (Step(Lmin) == 0 .and. CurrentLevel == Lmin) &
!!$         Dtime(Lmin) = Dtime(Lmin) * 1.D-6
    ! ------------------
    ! ����ȥ��ƥåפ�ɽ��
    ! ------------------
    if ( get_myrank() == PRIMARY_RANK ) then
       t(Lmin:LevelMax) = Time(Lmin:LevelMax)
       dt(Lmin:LevelMax) = Dtime(Lmin:LevelMax)
       st(Lmin:LevelMax) = Step(Lmin:LevelMax)
       dst(Lmin:LevelMax) = Dstep(Lmin:LevelMax)
       t(CurrentLevel)=t(CurrentLevel)+dt(CurrentLevel)
       st(CurrentLevel)=st(CurrentLevel)+dst(CurrentLevel)
       write(*,*) (st(l),l=Lmin,LevelMax)
       write(*,*) (dst(l),l=Lmin,LevelMax)
       write(*,*) (t(l),l=Lmin,LevelMax)
       write(*,*) (dt(l),l=Lmin,LevelMax)
    endif
  end subroutine cfl_condition_singlestep
  !-----------------------------------------------------------------------
  ! �ƥ���åɤ�Ʊ�����Ƥ��뤫�� level = 0 �ΤȤ��� FALSE ���֤�
  !-----------------------------------------------------------------------
  function bool_sync() result(bool)
    use grid
    use io_util
    logical :: bool
    bool = .FALSE.
    if ( CurrentLevel == Lmin ) return
    ! �ƥ���å�Ʊ�Τ�Ʊ�����Ƥ���Ȳ���
    if ( Step(CurrentLevel) /= Step(CurrentLevel-1) ) return
    bool = .TRUE.
    call print_msg('sync to parent level')
  end function bool_sync
  ! -----------------------------------------------------------------
  ! bool halt
  ! -----------------------------------------------------------------
  function bool_halt() result(bool)
    use mpilib
    use grid
    use io_util, only : readenv, print_msg
    use modelParameter, only : MP_T_Last
    logical :: bool
    integer(kind=LLONG_KIND),save :: Hz, Clock0, Clock1
    real(kind=8),save :: ElapseLimit = -1.d0 ! Limit of elapsetime in hour

!!$    bool = .false.
!!$    if ( level_sync() == Lmin ) then
!!$       bool = .true.
!!$       return
!!$    end if

    ! -----------
    ! �����
    ! -----------
    if ( ElapseLimit < 0.e0 ) then
       call system_clock(count=Clock0, count_rate=Hz)
       if ( .not. readenv('ELAPSELIMIT', ElapseLimit) ) then
          call print_msg( '*** error: environment variable ELAPSELIMIT is required.')
          stop
       endif
       if ( Hz == 0 ) then
          call print_msg('*** system_clock returns 0 Hz via count_rate. 1000 Hz is adopted here.')
          Hz = 1000
       endif
       call print_msg('initialize bool_halt')
    endif

    bool = .false.

#ifdef SINGLE_STEP
    if ( level_sync() /= Lmin ) return
#endif

!!$    ! for debug
!!$    if (Step(Lmin) >= 2444) then
!!$       bool = .true.
!!$       return
!!$    endif

    ! for profiling
!!$    if ( Time(Lmin) >= 0.2d0 ) then
!!$    if ( Step(Lmin) == 500 ) then
!!$       bool = .true.
!!$       return
!!$    end if


    ! --------------------
    ! MP_T_LAST �ǽ�λ����
    ! --------------------
    if ( level_sync() == Lmin .and. &
         MP_T_Last > Time(Lmin) - Dtime(Lmin)*0.001 .and. &
         MP_T_Last < Time(Lmin) + Dtime(Lmin)*0.001 ) then
       call print_msg('time reaches T_LAST')
       bool =  .true.
       return
    endif

    ! ---------------------
    ! �в���֤ˤ�꽪λ���롣
    ! ---------------------
    if ( get_myrank() == PRIMARY_RANK ) then
       call system_clock(count=Clock1, count_rate=Hz)
       if ( ( Clock1 - Clock0 )/dble(Hz) >= ElapseLimit*3600.d0 ) bool = .true.
       print *, 'clock in sec:', ( Clock1 - Clock0 )/dble(Hz),  ElapseLimit*3600.d0
    endif
    call mpi_bcast( bool, 1, MPI_LOGICAL, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    if ( bool ) then
       call print_msg('ElapseTime exceeds ElapseLimit')
    endif
  end function bool_halt
  ! -----------------------------------------------------------------
  ! ���Ϥ򹹿����뤫?
  ! -----------------------------------------------------------------
  function bool_fmg() result(bool)
    use grid
    logical :: bool
    integer :: stp
    integer,save :: fmg_laststep=-1 ! �Ƕ���Ϥ�׻��������ƥå�
    integer,save :: fmg_levmax=-1   ! ���ΤȤ��� LevelMax

    bool = .false.
#ifdef SINGLE_STEP
    if ( level_sync() == Lmin ) bool = .TRUE.
#else
    if ( Step(CurrentLevel) /= fmg_laststep .or. LevelMax /= fmg_levmax ) then
       bool = .TRUE.
       fmg_laststep = Step(CurrentLevel)
       fmg_levmax = LevelMax
    endif
#endif
  end function bool_fmg
  ! -----------------------------------------------------------------
  ! ��«���ƥå׻Ҥ���Ƥ�
  ! -----------------------------------------------------------------
  subroutine converge
    use grid
    use fg2cg
    integer:: n
    if ( CurrentLevel == Lmin ) return
    call fg2cg_u( CurrentLevel )
  end subroutine converge
end module multi_timestep