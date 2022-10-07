module multi_timestep
  implicit none
  private
  integer,save :: CurrentLevel 
  integer :: time_ini, time_prev, time_cur,time_rat
  integer :: time_hydro, time_grav, time_sp, time_ref, time_chem, time_other, time_tot
  public :: multi_timestep_init, step_all_level
contains
  subroutine multi_timestep_init
  end subroutine multi_timestep_init
  subroutine step_all_level
    use grid, only : Lmin, Step, LevelMax, level_sync, globdbg_myrank
    use timestep
    use reflux
    use analysis
    use rescue
    use unit
    use outputdata
    use sinkParticle
    use chemistry 
    use radiationSource 
    use grid, only : Time,Dtime,Step,Dstep 
    use mpilib 
    use radtr
    integer:: n
    real(kind=8),dimension(:,:),pointer :: pos
    real(kind=8),dimension(:),pointer :: lum
    do
       if (level_sync() == Lmin) then
          call system_clock(time_ini) 
          time_hydro=0
          time_grav=0
          time_sp=0
          time_ref=0
          time_chem=0
       endif
       call findCurrentLevel
       call system_clock(time_prev) 
       call grid_refinement
       call system_clock(time_cur) 
       time_ref = time_ref+ (time_cur - time_prev)
       globdbg_myrank = get_myrank()
       call cfl_condition_singlestep
       call system_clock(time_prev) 
       call step_all_grid( CurrentLevel )
       if ( bool_sync() ) then
          call fluxcorrection( CurrentLevel - 1 ) 
          call converge
          call rescueLev( CurrentLevel - 1 )
       endif
       call system_clock(time_cur) 
       time_hydro = time_hydro+ (time_cur - time_prev)
       if ( level_sync() == Lmin ) then
          call system_clock(time_prev) 
          call sp_update
          call rescueAllLev 
          call system_clock(time_cur) 
          time_sp = time_sp + (time_cur - time_prev)
          call system_clock(time_prev) 
          call radiation_source 
          call radtr_moment(Dtime(Lmin))
          call include_radforce(Dtime(Lmin))
          call system_clock(time_cur) 
          time_chem = time_chem + (time_cur - time_prev)
       endif
       call analysis_keyparam
       call output_data
       if (level_sync() == Lmin) then
          call system_clock(time_cur,time_rat) 
          time_tot = time_cur - time_ini
          time_other = time_tot - (time_hydro+time_grav+time_sp+time_ref+time_chem)
          if(get_myrank() == 0) &
               print '(/,A,/,7(1P1E12.4),/)', &
               "*** time_hydro,   time_grav,   time_sp,   time_ref,   time_chem,   time_other,   time_tot [s] ***", &
               time_hydro/dble(time_rat), time_grav/dble(time_rat), time_sp/dble(time_rat), time_ref/dble(time_rat), &
               time_chem/dble(time_rat), time_other/dble(time_rat), time_tot/dble(time_rat)
       end if
       if ( bool_halt() ) return
    end do
  end subroutine step_all_level
  subroutine findCurrentLevel
    use grid
    integer:: level
    integer(kind=8) :: nstepmin
    nstepmin = huge( nstepmin ) 
    do level = Lmin, LevelMax
       nstepmin = min( nstepmin, Step(level))
    enddo
    do level = Lmin, LevelMax
       if ( Step(level) == nstepmin ) then
          CurrentLevel = level
          return
       endif
    enddo
  end subroutine findCurrentLevel
  subroutine grid_refinement
    use grid
    use io_util
    use refine
    use string
    integer :: lev, slevel, levmax
    logical :: bool
    integer :: itr_refine, itr_max=10
    logical :: bool_update
    if ( .not. bool_refinement() ) return
    slevel = level_sync()
    levmax = LevelMax
    bool_update = .False. 
    do lev = levmax+1, slevel+1, -1 
       call refineLevel(lev, bool)
       if (globdbg_gridupdate) bool_update = .True. 
    enddo
  end subroutine grid_refinement
  function bool_refinement_BAK2() result(bool)
    use grid
    logical :: bool
    bool = .true.
  end function bool_refinement_BAK2
  function bool_refinement() result(bool)
    use grid
    use mpilib 
    logical :: bool
    bool = .false.
    if ( level_sync() == Lmin ) bool = .true.
  end function bool_refinement
  subroutine cfl_condition
    use eos
    use grid
    use io_util
    use string
    use mpilib
    use modelParameter, only : MP_T_Last
    use sinkParticle
    real(kind=8),dimension(Lmin:Lmax) :: t, dt 
    integer(kind=8),dimension(Lmin:Lmax) :: st, dst 
    integer :: n, gid, l
    real(kind=8) :: dtlocal, thisdt_try, thisdt, dtbuf
    integer(kind=8) :: thisds_try, syncstep, ds, thisds
    dtlocal = huge( dtlocal )
    do n = Gidmin, GidListMax( CurrentLevel )
       gid = GidList( n, CurrentLevel )
       dtlocal = min( dtlocal , get_dt_by_cflcond( gid ) )
    enddo
    call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    dtlocal = dtbuf
    call sp_restrictCFL(dtlocal)
    if (dtlocal < 0) then
       write(*,*) 'negative dtlocal', dtlocal
       stop
    endif
    t(Lmin:LevelMax) = Time(Lmin:LevelMax)
    dt(Lmin:LevelMax) = Dtime(Lmin:LevelMax)
    st(Lmin:LevelMax) = Step(Lmin:LevelMax)
    dst(Lmin:LevelMax) = Dstep(Lmin:LevelMax)
    if ( CurrentLevel == Lmin ) then
       thisdt = dtlocal
       thisds = 1
    else
       n=0
       do
          thisds_try = dst( CurrentLevel - 1 ) / 2**n
          thisdt_try = dt( CurrentLevel - 1 ) / 2**n
          if (thisds_try == 0) then
             syncstep = st(Lmin)-dst(Lmin) 
             do l=Lmin,LevelMax
                ds = st(l)-syncstep 
                st(l) = syncstep + ds*2
                dst(l) = dst(l)*2
             enddo
             cycle
          endif
          ds = st(CurrentLevel)-st(CurrentLevel-1)+dst(CurrentLevel-1) 
          if ((mod(ds,thisds_try) == 0) .and. (thisdt_try <= dtlocal)) exit
          n=n+1
       enddo
       thisdt = thisdt_try
       thisds = thisds_try
    endif
    dt( CurrentLevel ) = thisdt
    dst( CurrentLevel ) = thisds
    t( level_sync() : LevelMax-1 ) = t( LevelMax )
    if (CurrentLevel == Lmin) &
         dt(Lmin) = min(dt(Lmin), MP_T_Last - t(Lmin))
    Time(Lmin:LevelMax) = t(Lmin:LevelMax)
    Dtime(Lmin:LevelMax) = dt(Lmin:LevelMax)
    Step(Lmin:LevelMax) = st(Lmin:LevelMax)
    Dstep(Lmin:LevelMax) = dst(Lmin:LevelMax)
    if ( get_myrank() == 0 ) then
       t(CurrentLevel)=t(CurrentLevel)+dt(CurrentLevel)
       st(CurrentLevel)=st(CurrentLevel)+dst(CurrentLevel)
       write(*,*) (st(l),l=Lmin,LevelMax)
       write(*,*) (dst(l),l=Lmin,LevelMax)
       write(*,*) (t(l),l=Lmin,LevelMax)
       write(*,*) (dt(l),l=Lmin,LevelMax)
    endif
  end subroutine cfl_condition
  subroutine cfl_condition_singlestep
    use eos
    use grid
    use io_util
    use string
    use mpilib
    use modelParameter, only : MP_T_Last
    use sinkParticle
    real(kind=8),dimension(Lmin:Lmax) :: t, dt 
    integer(kind=8),dimension(Lmin:Lmax) :: st, dst 
    integer :: n, gid, l, lev
    real(kind=8) :: dtlocal, dtbuf
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
       call sp_restrictCFL(dtlocal)
       if (dtlocal < 0) then
          write(*,*) 'negative dtlocal', dtlocal
          stop
       endif
    endif
    if ( CurrentLevel == Lmin ) then
       Dstep(CurrentLevel) = 1
       Dtime(CurrentLevel) = dtlocal
    else
       Dstep(CurrentLevel) = Dstep(CurrentLevel-1)
       Dtime(CurrentLevel) = Dtime(CurrentLevel-1)
    endif
    if (CurrentLevel == Lmin) &
         Dtime(Lmin) = min(Dtime(Lmin), MP_T_Last - Time(Lmin))
    if (Step(Lmin) == 0 .and. CurrentLevel == Lmin) &
         Dtime(Lmin) = Dtime(Lmin) * 1.D-6
    if ( get_myrank() == 0 ) then
       t(Lmin:LevelMax) = Time(Lmin:LevelMax)
       dt(Lmin:LevelMax) = Dtime(Lmin:LevelMax)
       st(Lmin:LevelMax) = Step(Lmin:LevelMax)
       dst(Lmin:LevelMax) = Dstep(Lmin:LevelMax)
       t(CurrentLevel)=t(CurrentLevel)+dt(CurrentLevel)
       st(CurrentLevel)=st(CurrentLevel)+dst(CurrentLevel)
       write(*,*) &
         'level = ' // trim( num2char(CurrentLevel) ) // &
         ',    step = ' // trim( num2char(int(st(CurrentLevel))) ) // &
         ',    dt = ' // trim( num2char(dt(CurrentLevel)) ) // &
         ',    t = ' // trim( num2char(t(CurrentLevel)) )
    endif
  end subroutine cfl_condition_singlestep
  function bool_sync() result(bool)
    use grid
    use io_util
    logical :: bool
    bool = .FALSE.
    if ( CurrentLevel == Lmin ) return
    if ( Step(CurrentLevel) /= Step(CurrentLevel-1) ) return
    bool = .TRUE.
    call print_msg('sync to parent level')
  end function bool_sync
  function bool_halt() result(bool)
    use mpilib
    use grid
    use io_util, only : readenv, print_msg
    use string
    use modelParameter, only : MP_T_Last, MP_ElapseLimit, MP_CONNECTION_RUN 
    use unit 
    logical :: bool
    integer,save :: Hz, Clock0, Clock1, ClockMax
    integer,save :: ClockPrev = 0
    integer(kind=8),save :: ClockOffset = 0
    integer(kind=8) :: clockNow
    real(kind=8),save :: ElapseLimit = -1.d0 
    if (level_sync() == Lmin .and. MP_CONNECTION_RUN > 1) then
       call print_msg('*** THIS IS CONNECTION RUN ***')
       if(get_myrank() == 0) &
            print '(A,I4,A,/)', "transitional peirod continues for another ", MP_CONNECTION_RUN-1, " steps"
       MP_CONNECTION_RUN = MP_CONNECTION_RUN -1
    else if (level_sync() == Lmin .and. MP_CONNECTION_RUN == 1) then
       call print_msg('*** THIS IS CONNECTION RUN ***')
       call print_msg('stopping as transitional peirod ends')
       bool = .true.
       return
    end if
    if ( ElapseLimit < 0.e0 ) then
       call system_clock(count=Clock0, count_rate=Hz, count_max=ClockMax)
       ElapseLimit = MP_ElapseLimit
       if ( Hz == 0 ) then
          call print_msg('*** system_clock returns 0 Hz via count_rate. 1000 Hz is adopted here.')
          Hz = 1000
       endif
       call print_msg('initialize bool_halt')
    endif
    bool = .false.
    if ( level_sync() /= Lmin ) return
    if (level_sync() == Lmin .and. get_myrank() == 0 ) then 
       print '(A,4((1P1E15.7),A))', 'simulation time (yr): ', Time(Lmin)*Unit_yr, &
            ' / ',MP_T_Last*Unit_yr
    end if
    if ( level_sync() == Lmin .and. &
         MP_T_Last > Time(Lmin) - Dtime(Lmin)*0.001 .and. &
         MP_T_Last < Time(Lmin) + Dtime(Lmin)*0.001 ) then
       call print_msg('time reaches T_LAST')
       bool = .true.
       return
    endif
    if ( get_myrank() == 0 ) then
       call system_clock(count=Clock1, count_rate=Hz)
       if ( Clock1 < ClockPrev ) ClockOffset = ClockOffset + ClockMax
       clockNow = Clock1 + ClockOffset
       if ( ( clockNow - Clock0 )/dble(Hz) >= ElapseLimit*3600.d0 ) bool = .true.
       call print_msg( 'real time (sec): ' // &
            trim(num2char(( clockNow - Clock0 )/dble(Hz))) // ' / ' // &
            trim(num2char(ElapseLimit*3600.d0)) )
       ClockPrev = Clock1
    endif
    call mpi_bcast( bool, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if ( bool ) then
       call print_msg('ElapseTime exceeds ElapseLimit')
       return
    endif
    bool = bool_halt_byfile()
    if ( bool ) then
       call print_msg('Halt file exists')
       return
    endif
  end function bool_halt
  function bool_halt_byfile() result(bool)
    use mpilib
    use io_util, only : read_env
    use string, only : CHARLEN, concat
    character(len=CHARLEN),parameter :: FILENAME='HALT'
    character(len=CHARLEN) :: file, dir
    logical :: bool
    if ( get_myrank() == 0 ) then
       call read_env('DIR', dir)
       file = concat(dir,FILENAME)
       inquire(file=file, exist=bool)
    endif
    call mpi_bcast( bool, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  end function bool_halt_byfile
  function bool_fmg() result(bool)
    use grid
    logical :: bool
    integer :: stp
    integer,save :: fmg_laststep=-1 
    integer,save :: fmg_levmax=-1 
    bool = .false.
    if ( level_sync() == Lmin ) bool = .TRUE.
  end function bool_fmg
  subroutine converge
    use grid
    use fg2cg
    integer:: n
    if ( CurrentLevel == Lmin ) return
    call fg2cg_u( CurrentLevel )
  end subroutine converge
  subroutine converge_alllevel
    use grid
    use fg2cg
    integer :: level
    do level = LevelMax, Lmin+1, -1
       call fg2cg_u( level )
    end do
  end subroutine converge_alllevel
  subroutine include_radforce(dt)
    use overBlockCoordinates
    use eos, only : w2u_4, u2w_4
    integer :: level, n, gid
    integer :: i, j, k
    real(kind=8), intent(IN) :: dt
    real(kind=8) :: dv, dtdvrho
    real(kind=8),dimension(:,:,:,:),pointer :: u, w
    do level = Lmin, Lmax
      dv = get_dv(level)
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        u => get_Up(gid)
        allocate( w(lbound(u,1):ubound(u,1),lbound(u,2):ubound(u,2),lbound(u,3):ubound(u,3),0:4) )
        call u2w_4(u, w, dv)
        do k = Kmin, Kmax
          do j = Jmin, Jmax
            do i = Imin, Imax
              dtdvrho = dt*w(i,j,k,0) 
              w(i,j,k,1) = w(i,j,k,1) + u(i,j,k,14)*dtdvrho 
              w(i,j,k,2) = w(i,j,k,2) + u(i,j,k,15)*dtdvrho
              w(i,j,k,3) = w(i,j,k,3) + u(i,j,k,16)*dtdvrho
              w(i,j,k,4) = w(i,j,k,4) &
                +dtdvrho*(u(i,j,k,1)*u(i,j,k,14)+u(i,j,k,2)*u(i,j,k,15)+u(i,j,k,3)*u(i,j,k,16))
            enddo
          enddo
        enddo
        call w2u_4(w, u, dv)
        deallocate(w)
      enddo
    enddo
  end subroutine include_radforce
end module multi_timestep
