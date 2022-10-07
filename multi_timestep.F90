#include "config.h"
! #define USE_RESCUE_RESTRICTCFL
!-------------------------------------------------------------------------
! Module for multi-timestep
!
!-------------------------------------------------------------------------
module multi_timestep
  implicit none
  private
  integer,save :: CurrentLevel       ! 現在のレベル

  !---------------- KS DEBUG (BEGIN)---------------!
  integer :: time_ini, time_prev, time_cur,time_rat
  integer :: time_hydro, time_grav, time_sp, time_ref, time_chem, time_other, time_tot
  !----------------- KS DEBUG (END) ----------------!

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
    ! use grid, only : Lmin, Step, LevelMax, level_sync
    use grid, only : Lmin, Step, LevelMax, level_sync, globdbg_myrank
    use timestep
    use reflux
    use analysis
    use rescue
    use unit
    use outputdata
#if defined(WITH_SELFGRAVITY) || ( defined(FMG_OHMIC_DISSIPATION) && defined(SINGLE_STEP) )
    use fmg
#endif
#if defined(SINKPARTICLE)
    use sinkParticle
#endif
    use chemistry       ! KS ADDED
#if MODEL_ART > 0
    use radiationSource ! KS ADDED
    use grid, only : Time,Dtime,Step,Dstep  ! KS ADDED
#endif
    use mpilib  !KS DEBUG
#ifdef RADTR_M1closer 
    use radtr
#endif

    integer:: n
    real(kind=DBL_KIND),dimension(:,:),pointer :: pos
    real(kind=DBL_KIND),dimension(:),pointer :: lum    

    do
       !------- KS DEBUG -------!
       if (level_sync() == Lmin) then
          call system_clock(time_ini) ! 時間計測開始
          time_hydro=0
          time_grav=0
          time_sp=0
          time_ref=0
          time_chem=0
       endif
       !------- KS DEBUG -------!


       call findCurrentLevel

       call system_clock(time_prev) ! 時間計測 (KS DEBUG)

       call grid_refinement

       !------- KS DEBUG -------!
       call system_clock(time_cur) ! 時間計測
       time_ref = time_ref+ (time_cur - time_prev)
       !------- KS DEBUG -------!

       globdbg_myrank = get_myrank() 
       ! if( .not. bool_refinement()  .and. get_myrank() == PRIMARY_RANK) &
       !      print '(A)', "(MT, KS DEBUG) *** WARNING *** call grid_refinement 2 times"       
       ! call grid_refinement
       !------- KS DEBUG -------!       

       ! ----------------
       ! CFL condition
       ! ----------------
#ifdef SINGLE_STEP
       call cfl_condition_singlestep
#else
       call cfl_condition
#endif ! SINGLE_STEP

       ! ------------------------------------------
       ! Multi grid for gravity for multi-timestep
       ! ------------------------------------------
#if defined(WITH_SELFGRAVITY) && !defined(SINGLE_STEP)
       if ( bool_fmg() ) call fmg_poisson
#endif

       ! ----------------------
       ! Solve Hydrodynamics
       ! ----------------------
#ifndef SKIP_HD_UPDATE
       call system_clock(time_prev) ! 時間計測 (KS DEBUG)

       call step_all_grid( CurrentLevel )

       ! ------------------------------------------
       ! Consistency between fine and coarse grids
       ! ------------------------------------------
       if ( bool_sync() ) then
          call fluxcorrection( CurrentLevel - 1 ) ! for parent grid
          call converge
          call rescueLev( CurrentLevel - 1 )
       endif

       !------- KS DEBUG -------!
       call system_clock(time_cur) ! 時間計測
       time_hydro = time_hydro+ (time_cur - time_prev)
       !------- KS DEBUG -------!
#else
       !hydroを解く代わりにtimeとstepだけ進める
       Time( CurrentLevel ) = Time( CurrentLevel ) + Dtime( CurrentLevel )
       Step( CurrentLevel ) = Step( CurrentLevel ) + Dstep( CurrentLevel )
#endif

       ! ------------------------------------------
       ! Ohmic dissipation
       ! ------------------------------------------
#if defined(FMG_OHMIC_DISSIPATION) && defined(SINGLE_STEP)
       if ( level_sync() == Lmin ) call fmg_ohmic_dissipation
#endif
       ! ------------------------------------------
       ! Multi grid for gravity for single-timestep
       ! ------------------------------------------
#if defined(WITH_SELFGRAVITY) && defined(SINGLE_STEP)
       call system_clock(time_prev) ! 時間計測 (KS DEBUG)

       if ( bool_fmg() ) then
          call rescue_rhopsi_for_poisson ! KS ADDED
          if (Step(Lmin) == 1) call fmg_poisson
          call fmg_poisson
          call source_g_all_level
          call rescueAllLev
       endif

       !------- KS DEBUG -------!
       call system_clock(time_cur) ! 時間計測
       time_grav = time_grav + (time_cur - time_prev)
       !------- KS DEBUG -------!
#endif

       ! -------------
       ! sink cell
       ! -------------
#if defined(SINKPARTICLE) && defined(SINGLE_STEP)
       if ( level_sync() == Lmin ) then
#ifndef CREATE_RADSOURCE_BY_HAND
          call system_clock(time_prev) ! 時間計測 (KS DEBUG)

          call sp_update
          call rescueAllLev ! KS ADDED

          !------- KS DEBUG -------!
          call system_clock(time_cur) ! 時間計測
          time_sp = time_sp + (time_cur - time_prev)
          !------- KS DEBUG -------!
#endif ! CREATE_RADSOURCE_BY_HAND          

#if MODEL_ART > 0
          call system_clock(time_prev) ! 時間計測 (KS DEBUG)

          ! -------------
          ! radiation source
          ! -------------
#ifdef  STOCHASTIC_STELLAR_MODEL
          call stochastic_radiation_source 
#else
          call radiation_source  !KS ADDED
#endif


          ! -------------
          ! ART + Chemistry
          ! -------------
#ifdef RADTR_M1_ONLY
          call radtr_moment(Dtime(Lmin))
#else
          call ch_artchem(Dtime(Lmin))  !KS ADDED
#endif
          ! ------------------    ! HF ADDED
          !  radiation force
          ! ------------------
#ifdef EXTERNALFORCE
          call include_radforce(Dtime(Lmin)) 
#endif

          !------- KS DEBUG -------!
          call system_clock(time_cur) ! 時間計測
          time_chem = time_chem + (time_cur - time_prev)
          !------- KS DEBUG -------!       
#endif ! MODEL_ART


#ifndef RADTR_M1_ONLY
          call converge_alllevel !KS MODIFIED
#endif

       endif
#endif ! SINKPARTICLE

#ifdef ISOTHERMAL
       call ch_reset_temperature()
#endif !ISOTHERMAL


#if defined(SINKPARTICLE) && !defined(SINGLE_STEP)
       if ( CurrentLevel == LevelMax ) then
          call sp_update
       endif
#endif ! SINKPARTICLE


       call analysis_keyparam
       call output_data

       !------- KS DEBUG -------!
       if (level_sync() == Lmin) then
          call system_clock(time_cur,time_rat) ! 時間計測
          time_tot = time_cur - time_ini
          time_other = time_tot - (time_hydro+time_grav+time_sp+time_ref+time_chem)
          if(get_myrank() == PRIMARY_RANK) &
               print '(/,A,/,7(1P1E12.4),/)', &
               "*** time_hydro,   time_grav,   time_sp,   time_ref,   time_chem,   time_other,   time_tot [s] ***", &
               time_hydro/dble(time_rat), time_grav/dble(time_rat), time_sp/dble(time_rat), time_ref/dble(time_rat), &
               time_chem/dble(time_rat), time_other/dble(time_rat), time_tot/dble(time_rat)
       end if
       !------- KS DEBUG -------!

       if ( bool_halt() )  return
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

    !KS ADDED
    integer :: itr_refine, itr_max=10
    logical :: bool_update

#ifdef NOT_REFINEMENT
    return
#endif !NOT_REFINEMENT

    if ( .not. bool_refinement() ) return
    slevel = level_sync()
    levmax = LevelMax

    bool_update = .False. !KS DEBUG

    do lev = levmax+1, slevel+1, -1 ! 細から粗へ
       call refineLevel(lev, bool)
       if (globdbg_gridupdate) bool_update = .True. !gridの張り替えがあったらTrue (KS DEBUG)
    enddo

    !------- KS DEBUG (multiple refinement) -------!
    ! do itr_refine = 0, itr_max-1
    !    if (.not. bool_update) exit !gridの張り替えが無くなったら終了
    !    bool_update = .False.
    !    if (get_myrank() == PRIMARY_RANK) &
    !         print '(/,A,/)', "(grid_refinement, KS DEBUG) call refineLevel again"

    !    do lev = levmax+1, slevel+1, -1 ! 細から粗へ
    !       call refineLevel(lev, bool) 
    !       if (globdbg_gridupdate) then    ! gridの張り替えがあったか？
    !          bool_update = .True.         ! あればTrue
    !          if (get_myrank() == PRIMARY_RANK) &
    !               print '(A,I0,A,I4)', "*** WARNING *** grid refinement DONE in ",itr_refine+2,"th trial ", lev
    !       end if
    !    enddo

    !    if (itr_refine==itr_max-1) then  !収束しなかったら止める
    !       if (get_myrank() == PRIMARY_RANK) &
    !            print '(A,/,A)', "(grid_refinement, KS DEBUG) refinement not converged","stopping..."
    !       stop
    !    end if
    ! end do
    !------- KS DEBUG -------!       


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
    use mpilib ! KS DEBUG
    logical :: bool
    bool = .false.
    if ( level_sync() == Lmin ) bool = .true.
    !---- WARNING ---- WARNING ---- (KS DEBUG) ---- WARNING ---- WARNING ----!
    ! if (level_sync() == Lmin .and. Step(Lmin) == 11000-1) then
    !    bool = .false.
    !    if(get_myrank() == PRIMARY_RANK) &
    !         print '(A)', "(MT, KS DEBUG) *** WARNING *** step = 11000 => skip refinement"
    ! end if
    !---- WARNING ---- WARNING ---- (KS DEBUG) ---- WARNING ---- WARNING ----!    

  end function bool_refinement
#else !SINGLE_STEP
  !-----------------------------------------------------------------------
  ! 子グリッドを refinement するための条件(重力更新と同じタイミング)
  !-----------------------------------------------------------------------
  function bool_refinement() result(bool)
    use mpilib
    use grid
    logical :: bool
    ! Interval of step for trying refinement. It should be 2**N.
    ! STEP_INTERVAL = 2 for the most frequent refinement.
#ifdef WITH_SELFGRAVITY
    integer,parameter :: STEP_INTERVAL = 2
#else !WITH_SELFGRAVITY
    integer,parameter :: STEP_INTERVAL = NI
#endif !WITH_SELFGRAVITY
    integer,save :: laststep=-1 ! 最近テストしたステップ
    integer,save :: count = 0
    integer :: countmax         ! maximum count
    integer :: levsync_req ! このレベルかそれ以下のレベルまで同期すると細分化を検討する。
    integer :: levsync     ! 実際に同期したレベル
    integer :: l
    bool = .FALSE.

    if ( Step(CurrentLevel) == laststep ) return
    laststep = Step(CurrentLevel)

    levsync = level_sync()

    ! grid level have to be synchronized with Lmin and takes count
    countmax = STEP_INTERVAL / max(Dstep(Lmin),1)
    if ( countmax > 0 ) then
       if (levsync == Lmin) count = count + 1
       if (count >= countmax) then
          bool = .TRUE.
          count = 0
       end if
       return
    end if

    ! look for levsync_req
    levsync_req = Lmin - 1      ! initial value
    do l = Lmin, LevelMax
       if ( Dstep(l)/max(Dstep(LevelMax),1) <= STEP_INTERVAL ) then
          levsync_req = l
          exit
       end if
    end do
    if ( levsync <= levsync_req ) bool = .TRUE. ! 同期レベルが十分深い

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
#ifdef USE_RESCUE_RESTRICTCFL
    use rescue, only: rescue_restrictCFL
#endif ! USE_RESCUE_RESTRICTCFL
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
    call sp_restrictCFL(dtlocal)
#endif ! SINKPARTICLE
#ifdef USE_RESCUE_RESTRICTCFL
    call rescue_restrictCFL(dtlocal, CurrentLevel)
#endif !USE_RESCUE_RESTRICTCFL
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
    if ( CurrentLevel == Lmin ) then
       dtlocal = huge( dtlocal )
       do lev = Lmin, LevelMax
          do n = Gidmin, GidListMax( lev )
             gid = GidList( n, lev )
             dtlocal = min( dtlocal , get_dt_by_cflcond( gid ) )
          enddo
       enddo
       call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
       dtlocal = dtbuf
#ifdef SINKPARTICLE
#ifndef CREATE_RADSOURCE_BY_HAND
       call sp_restrictCFL(dtlocal)
#endif ! CREATE_RADSOURCE_BY_HAND
#endif ! SINKPARTICLE
       if (dtlocal < 0) then
          write(*,*) 'negative dtlocal', dtlocal
          stop
       endif
    endif
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
    ! -----------------------------------------------------
    ! 初回は微少なタイムステップを刻む(Predictor Stepのため)
    ! -----------------------------------------------------
    if (Step(Lmin) == 0 .and. CurrentLevel == Lmin) &
         Dtime(Lmin) = Dtime(Lmin) * 1.D-6
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
       ! write(*,*) (st(l),l=Lmin,LevelMax)
       ! write(*,*) (dst(l),l=Lmin,LevelMax)
       ! write(*,*) (t(l),l=Lmin,LevelMax)
       ! write(*,*) (dt(l),l=Lmin,LevelMax)
       write(*,*) &
         'level = ' // trim( num2char(CurrentLevel) )     // &
         ',    step = '   // trim( num2char(int(st(CurrentLevel))) ) // &
         ',    dt = '    // trim( num2char(dt(CurrentLevel)) )  // &
         ',    t = '    // trim( num2char(t(CurrentLevel)) ) 
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
    use string
    !use modelParameter, only : MP_T_Last
    use modelParameter, only : MP_T_Last,  MP_ElapseLimit, MP_CONNECTION_RUN ! KS MODIFIED
    use unit ! KS ADDED    
    logical :: bool
    integer,save :: Hz, Clock0, Clock1, ClockMax
    integer,save :: ClockPrev = 0
    integer(kind=LLONG_KIND),save :: ClockOffset = 0
    integer(kind=LLONG_KIND) :: clockNow
    real(kind=8),save :: ElapseLimit = -1.d0 ! Limit of elapsetime in hour

!!$    bool = .false.
!!$    if ( level_sync() == Lmin ) then
!!$       bool = .true.
!!$       return
!!$    end if

    !----------------- stop because of connection run mode ----------------!
    if (level_sync() == Lmin .and. MP_CONNECTION_RUN > 1) then
       call print_msg('*** THIS IS CONNECTION RUN ***')       
       if(get_myrank() == PRIMARY_RANK) &
            print '(A,I4,A,/)', "transitional peirod continues for another ", MP_CONNECTION_RUN-1, " steps"
       MP_CONNECTION_RUN = MP_CONNECTION_RUN -1       
    else if (level_sync() == Lmin .and. MP_CONNECTION_RUN == 1) then
       call print_msg('*** THIS IS CONNECTION RUN ***')
       call print_msg('stopping as transitional peirod ends')
       bool =  .true.
       return    
    end if
    !----------------------------------------------------------------------!

    ! -----------
    ! 初期化
    ! -----------
    if ( ElapseLimit < 0.e0 ) then
       call system_clock(count=Clock0, count_rate=Hz, count_max=ClockMax)
       !---------- KS MODIFIED -----------!
       ! if ( .not. readenv('ELAPSELIMIT', ElapseLimit) ) then
       !    call print_msg( '*** error: environment variable ELAPSELIMIT is required.')
       !    stop
       ! endif
       ElapseLimit = MP_ElapseLimit
       !---------- KS MODIFIED -----------!
       if ( Hz == 0 ) then
          call print_msg('*** system_clock returns 0 Hz via count_rate. 1000 Hz is adopted here.')
          Hz = 1000
       endif
       call print_msg('initialize bool_halt')
    endif

    bool = .false.

#ifndef HALT_WHENEVER_WO_SYNC
    if ( level_sync() /= Lmin ) return
#endif !HALT_WHENEVER_WO_SYNC

    ! for debug
!!$    if (Step(Lmin) >= 2999) then
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
    if (level_sync() == Lmin .and.  get_myrank() == PRIMARY_RANK ) then !KS ADDED
       print '(A,4((1P1E15.7),A))', 'simulation time (yr): ', Time(Lmin)*Unit_yr, &
            ' / ',MP_T_Last*Unit_yr
    end if
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
       if ( Clock1 < ClockPrev ) ClockOffset = ClockOffset + ClockMax
       clockNow = Clock1 + ClockOffset
       if ( ( clockNow - Clock0 )/dble(Hz) >= ElapseLimit*3600.d0 ) bool = .true.
       !KS MODIFIED
       call print_msg( 'real time (sec): ' // &
            trim(num2char(( clockNow - Clock0 )/dble(Hz))) // ' / ' // &
            trim(num2char(ElapseLimit*3600.d0)) )
       ClockPrev = Clock1
    endif
    call mpi_bcast( bool, 1, MPI_LOGICAL, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    if ( bool ) then
       call print_msg('ElapseTime exceeds ElapseLimit')
       return
    endif

    ! --------------------------
    ! 終了ファイルがあると終了する。
    ! --------------------------
    bool = bool_halt_byfile()
    if ( bool ) then
       call print_msg('Halt file exists')
       return
    endif

  end function bool_halt
  ! -----------------------------------------------------------------
  ! 終了ファイルがあると終了する。
  ! -----------------------------------------------------------------
  function bool_halt_byfile() result(bool)
    use mpilib
    use io_util, only : read_env
    use string, only : CHARLEN, concat
    character(len=CHARLEN),parameter :: FILENAME='HALT'
    character(len=CHARLEN) :: file, dir
    logical :: bool
    if ( get_myrank() == PRIMARY_RANK ) then
       call read_env('DIR', dir)
       file = concat(dir,FILENAME)
       inquire(file=file, exist=bool)
    endif
    call mpi_bcast( bool, 1, MPI_LOGICAL, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
  end function bool_halt_byfile
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

  ! -----------------------------------------------------------------
  ! convergeを全レベルに対して実行 (化学更新後に使うことを想定, KS ADDED)
  ! -----------------------------------------------------------------
  subroutine converge_alllevel
    use grid
    use fg2cg
    integer :: level
    do level = LevelMax, Lmin+1, -1
       call fg2cg_u( level )
    end do
  end subroutine converge_alllevel


  ! ------------------------------------------------------------------
  !   include radiation force
  ! ------------------------------------------------------------------
#ifdef EXTERNALFORCE
 
#define M_MRHO  0
#define M_MVX   1
#define M_MVY   2
#define M_MVZ   3
#define M_MP    4

  subroutine include_radforce(dt)

    use overBlockCoordinates
    use eos, only : w2u_4, u2w_4

    integer :: level, n, gid

    !info for radiation source
    integer :: i, j, k
    real(kind=DBL_KIND), intent(IN) :: dt
    real(kind=DBL_KIND) :: dv, dtdvrho 
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, w


    ! include radforce -------------------------------
    do level = Lmin, Lmax
      dv = get_dv(level)
      do n = Gidmin, GidListMax( level )  
        gid = GidList( n, level )

        u => get_Up(gid)
        allocate( w(ARRAYSIZE4_4(u)) )
        call u2w_4(u, w, dv)

        do k = Kmin, Kmax
          do j = Jmin, Jmax
            do i = Imin, Imax    

              dtdvrho = dt*w(i,j,k,M_MRHO) ! dt*dv*rho      

              w(i,j,k,M_MVX) = w(i,j,k,M_MVX) + u(i,j,k,MXPI)*dtdvrho ![g cm s^-1]
              w(i,j,k,M_MVY) = w(i,j,k,M_MVY) + u(i,j,k,MYPI)*dtdvrho
              w(i,j,k,M_MVZ) = w(i,j,k,M_MVZ) + u(i,j,k,MZPI)*dtdvrho
#ifdef MP
              w(i,j,k,M_MP)  = w(i,j,k,M_MP) &
                +dtdvrho*(u(i,j,k,MVX)*u(i,j,k,MXPI)+u(i,j,k,MVY)*u(i,j,k,MYPI)+u(i,j,k,MVZ)*u(i,j,k,MZPI))
#endif
            enddo
          enddo
        enddo
    
        call w2u_4(w, u, dv)
        deallocate(w)

      enddo
    enddo
  end subroutine include_radforce
#endif


end module multi_timestep
