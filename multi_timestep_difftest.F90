
#include "config.h"
!-------------------------------------------------------------------------
! Module for multi-timestep
!
!-------------------------------------------------------------------------
module multi_timestep
  implicit none
  private
  integer,save :: CurrentLevel       ! 現在のレベル
  public :: multi_timestep_init, step_all_level
contains
  !-----------------------------------------------------------------------
  ! initialize module timestep
  !-----------------------------------------------------------------------
  subroutine multi_timestep_init
  end subroutine multi_timestep_init
  !-----------------------------------------------------------------------
  ! muti timestep で時間推進
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
  ! 時刻を進めるべきレベルを見つける。CurrentLevel をセット
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

    do lev = levmax+1, slevel+1, -1 ! 細から粗へ
       call refineLevel(lev, bool)
    enddo

!!$    do lev = slevel+1, Lmax     ! 粗から細へ
!!$       call refineLevel(lev, bool)
!!$       if ( .not. bool ) exit
!!$    enddo


!!$    print *, 'sync lev, this lev', slevel, CurrentLevel
!!$    if ( LevelMax > levmax ) &
!!$         call print_msg( 'new level = '// num2char(LevelMax))
  end subroutine grid_refinement
  !-----------------------------------------------------------------------
  ! 子グリッドを refinement するための条件(毎ステップ)
  !-----------------------------------------------------------------------
  function bool_refinement_BAK2() result(bool)
    use grid
    logical :: bool
    bool = .true.
  end function bool_refinement_BAK2
#ifdef SINGLE_STEP
  !-----------------------------------------------------------------------
  ! 子グリッドを refinement するための条件(最粗グリッドから全て refinement)
  !-----------------------------------------------------------------------
  function bool_refinement() result(bool)
    use grid
    logical :: bool
    bool = .false.
    if ( level_sync() == Lmin ) bool = .true.
  end function bool_refinement
#else !SINGLE_STEP
  !-----------------------------------------------------------------------
  ! 子グリッドを refinement するための条件(重力更新と同じタイミング)
  !-----------------------------------------------------------------------
  function bool_refinement() result(bool)
    use grid
    logical :: bool
    integer :: stp
    integer,save :: laststep=-1 ! 最近refimentしたステップ
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
    real(kind=DBL_KIND),dimension(Lmin:Lmax) :: t, dt ! 作業配列
    integer(kind=LLONG_KIND),dimension(Lmin:Lmax) :: st, dst ! 作業配列
    integer :: n, gid, l
    real(kind=DBL_KIND) :: dtlocal, thisdt_try, thisdt, dtbuf
    integer(kind=LLONG_KIND) :: thisds_try, syncstep, ds, thisds
    ! ------------------------------
    ! 自分勝手に時間幅(dtlocal)を求める
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
    ! 時刻情報を代入する
    ! ----------------
    t(Lmin:LevelMax) = Time(Lmin:LevelMax)
    dt(Lmin:LevelMax) = Dtime(Lmin:LevelMax)
    st(Lmin:LevelMax) = Step(Lmin:LevelMax)
    dst(Lmin:LevelMax) = Dstep(Lmin:LevelMax)
    ! --------------------------------------------------
    ! 各レベルで調整する。
    ! 各レベルの step, dstep が変更される。
    ! 現レベルの thisdt, thisds を求める。
    ! --------------------------------------------------
    if ( CurrentLevel == Lmin ) then
       thisdt = dtlocal
       thisds = 1
    else
       n=0
       do
          ! 試行ステップ
          thisds_try = dst( CurrentLevel - 1 ) / 2**n
          thisdt_try = dt( CurrentLevel - 1 ) / 2**n
          ! 時間分解能を増加
          if (thisds_try == 0) then
             syncstep = st(Lmin)-dst(Lmin) ! 最後にグローバル同期したステップ
             do l=Lmin,LevelMax
                ds = st(l)-syncstep  ! 最後のグローバル同期からの増分
                st(l) = syncstep + ds*2
                dst(l) = dst(l)*2
             enddo
             cycle
          endif
          ! 親グリッドとの関係とCFL条件を満たすように
          ds = st(CurrentLevel)-st(CurrentLevel-1)+dst(CurrentLevel-1) ! 親との同期からの増分
          if ((mod(ds,thisds_try) == 0) .and. (thisdt_try <= dtlocal)) exit
          n=n+1
       enddo
       thisdt = thisdt_try
       thisds = thisds_try
    endif
    dt( CurrentLevel ) = thisdt
    dst( CurrentLevel ) = thisds

    ! -----------------------
    ! 同期レベルで時刻をそろえる。
    ! 時刻の丸め誤差補正。
    ! -----------------------
    t( level_sync() : LevelMax-1 ) = t( LevelMax )

    ! -------------------------
    ! ちょうどの時刻に終了する
    ! -------------------------
    if (CurrentLevel == Lmin) &
         dt(Lmin) = min(dt(Lmin), MP_T_Last - t(Lmin))
    ! ----------------------------------------------------------
    ! global variable をアップデートする (まだ時間推進はしない)
    ! ----------------------------------------------------------
    Time(Lmin:LevelMax) = t(Lmin:LevelMax)
    Dtime(Lmin:LevelMax) = dt(Lmin:LevelMax)
    Step(Lmin:LevelMax) = st(Lmin:LevelMax)
    Dstep(Lmin:LevelMax) = dst(Lmin:LevelMax)
    ! -----------------------------------------------
    ! 時刻とステップを表示
    ! 表示するだけで、global variable は更新されない
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
    real(kind=DBL_KIND),dimension(Lmin:Lmax) :: t, dt ! 作業配列
    integer(kind=LLONG_KIND),dimension(Lmin:Lmax) :: st, dst ! 作業配列
    integer :: n, gid, l, lev
    real(kind=DBL_KIND) :: dtlocal, dtbuf
    ! ----------------------------------
    ! 全レベルで時間幅(dtlocal)を求める
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
    ! ステップ数・時刻の差分
    ! ----------------------
    if ( CurrentLevel == Lmin ) then
       Dstep(CurrentLevel) = 1
       Dtime(CurrentLevel) = dtlocal
    else
       Dstep(CurrentLevel) = Dstep(CurrentLevel-1)
       Dtime(CurrentLevel) = Dtime(CurrentLevel-1)
    endif
    ! -------------------------
    ! ちょうどの時刻に終了する
    ! -------------------------
    if (CurrentLevel == Lmin) &
         Dtime(Lmin) = min(Dtime(Lmin), MP_T_Last - Time(Lmin))
!!$    ! -----------------------------------------------------
!!$    ! 初回は微少なタイムステップを刻む(Predictor Stepのため)
!!$    ! -----------------------------------------------------
!!$    if (Step(Lmin) == 0 .and. CurrentLevel == Lmin) &
!!$         Dtime(Lmin) = Dtime(Lmin) * 1.D-6
    ! ------------------
    ! 時刻とステップを表示
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
  ! 親グリッドと同期しているか？ level = 0 のときは FALSE を返す
  !-----------------------------------------------------------------------
  function bool_sync() result(bool)
    use grid
    use io_util
    logical :: bool
    bool = .FALSE.
    if ( CurrentLevel == Lmin ) return
    ! 親グリッド同士は同期していると仮定
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
    ! 初期化
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
    ! MP_T_LAST で終了する
    ! --------------------
    if ( level_sync() == Lmin .and. &
         MP_T_Last > Time(Lmin) - Dtime(Lmin)*0.001 .and. &
         MP_T_Last < Time(Lmin) + Dtime(Lmin)*0.001 ) then
       call print_msg('time reaches T_LAST')
       bool =  .true.
       return
    endif

    ! ---------------------
    ! 経過時間により終了する。
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
  ! 重力を更新するか?
  ! -----------------------------------------------------------------
  function bool_fmg() result(bool)
    use grid
    logical :: bool
    integer :: stp
    integer,save :: fmg_laststep=-1 ! 最近重力を計算したステップ
    integer,save :: fmg_levmax=-1   ! そのときの LevelMax

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
  ! 収束ステップ子から親へ
  ! -----------------------------------------------------------------
  subroutine converge
    use grid
    use fg2cg
    integer:: n
    if ( CurrentLevel == Lmin ) return
    call fg2cg_u( CurrentLevel )
  end subroutine converge
end module multi_timestep
