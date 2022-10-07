#include "config.h"
! #define USE_RESCUE_RESTRICTCFL
!-------------------------------------------------------------------------
! Module for multi-timestep
!
!-------------------------------------------------------------------------
module multi_timestep
  implicit none
  private
  integer,save :: CurrentLevel       ! ���ߤΥ�٥�

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
  ! muti timestep �ǻ��ֿ��
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
          call system_clock(time_ini) ! ���ַ�¬����
          time_hydro=0
          time_grav=0
          time_sp=0
          time_ref=0
          time_chem=0
       endif
       !------- KS DEBUG -------!


       call findCurrentLevel

       call system_clock(time_prev) ! ���ַ�¬ (KS DEBUG)

       call grid_refinement

       !------- KS DEBUG -------!
       call system_clock(time_cur) ! ���ַ�¬
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
       call system_clock(time_prev) ! ���ַ�¬ (KS DEBUG)

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
       call system_clock(time_cur) ! ���ַ�¬
       time_hydro = time_hydro+ (time_cur - time_prev)
       !------- KS DEBUG -------!
#else
       !hydro��������time��step�����ʤ��
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
       call system_clock(time_prev) ! ���ַ�¬ (KS DEBUG)

       if ( bool_fmg() ) then
          call rescue_rhopsi_for_poisson ! KS ADDED
          if (Step(Lmin) == 1) call fmg_poisson
          call fmg_poisson
          call source_g_all_level
          call rescueAllLev
       endif

       !------- KS DEBUG -------!
       call system_clock(time_cur) ! ���ַ�¬
       time_grav = time_grav + (time_cur - time_prev)
       !------- KS DEBUG -------!
#endif

       ! -------------
       ! sink cell
       ! -------------
#if defined(SINKPARTICLE) && defined(SINGLE_STEP)
       if ( level_sync() == Lmin ) then
#ifndef CREATE_RADSOURCE_BY_HAND
          call system_clock(time_prev) ! ���ַ�¬ (KS DEBUG)

          call sp_update
          call rescueAllLev ! KS ADDED

          !------- KS DEBUG -------!
          call system_clock(time_cur) ! ���ַ�¬
          time_sp = time_sp + (time_cur - time_prev)
          !------- KS DEBUG -------!
#endif ! CREATE_RADSOURCE_BY_HAND          

#if MODEL_ART > 0
          call system_clock(time_prev) ! ���ַ�¬ (KS DEBUG)

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
          call system_clock(time_cur) ! ���ַ�¬
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
          call system_clock(time_cur,time_rat) ! ���ַ�¬
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

    do lev = levmax+1, slevel+1, -1 ! �٤����Ƥ�
       call refineLevel(lev, bool)
       if (globdbg_gridupdate) bool_update = .True. !grid��ĥ���ؤ������ä���True (KS DEBUG)
    enddo

    !------- KS DEBUG (multiple refinement) -------!
    ! do itr_refine = 0, itr_max-1
    !    if (.not. bool_update) exit !grid��ĥ���ؤ���̵���ʤä��齪λ
    !    bool_update = .False.
    !    if (get_myrank() == PRIMARY_RANK) &
    !         print '(/,A,/)', "(grid_refinement, KS DEBUG) call refineLevel again"

    !    do lev = levmax+1, slevel+1, -1 ! �٤����Ƥ�
    !       call refineLevel(lev, bool) 
    !       if (globdbg_gridupdate) then    ! grid��ĥ���ؤ������ä�����
    !          bool_update = .True.         ! �����True
    !          if (get_myrank() == PRIMARY_RANK) &
    !               print '(A,I0,A,I4)', "*** WARNING *** grid refinement DONE in ",itr_refine+2,"th trial ", lev
    !       end if
    !    enddo

    !    if (itr_refine==itr_max-1) then  !��«���ʤ��ä���ߤ��
    !       if (get_myrank() == PRIMARY_RANK) &
    !            print '(A,/,A)', "(grid_refinement, KS DEBUG) refinement not converged","stopping..."
    !       stop
    !    end if
    ! end do
    !------- KS DEBUG -------!       


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
  ! �ҥ���åɤ� refinement ���뤿��ξ��(���Ϲ�����Ʊ�������ߥ�)
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
    integer,save :: laststep=-1 ! �Ƕ�ƥ��Ȥ������ƥå�
    integer,save :: count = 0
    integer :: countmax         ! maximum count
    integer :: levsync_req ! ���Υ�٥뤫����ʲ��Υ�٥�ޤ�Ʊ������Ⱥ�ʬ����Ƥ���롣
    integer :: levsync     ! �ºݤ�Ʊ��������٥�
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
    if ( levsync <= levsync_req ) bool = .TRUE. ! Ʊ����٥뤬��ʬ����

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
    ! -----------------------------------------------------
    ! ���������ʥ����ॹ�ƥåפ���(Predictor Step�Τ���)
    ! -----------------------------------------------------
    if (Step(Lmin) == 0 .and. CurrentLevel == Lmin) &
         Dtime(Lmin) = Dtime(Lmin) * 1.D-6
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
    ! �����
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
    ! MP_T_LAST �ǽ�λ����
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
    ! �в���֤ˤ�꽪λ���롣
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
    ! ��λ�ե����뤬����Ƚ�λ���롣
    ! --------------------------
    bool = bool_halt_byfile()
    if ( bool ) then
       call print_msg('Halt file exists')
       return
    endif

  end function bool_halt
  ! -----------------------------------------------------------------
  ! ��λ�ե����뤬����Ƚ�λ���롣
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

  ! -----------------------------------------------------------------
  ! converge������٥���Ф��Ƽ¹� (���ع�����˻Ȥ����Ȥ�����, KS ADDED)
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