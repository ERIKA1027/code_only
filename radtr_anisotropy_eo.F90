#include "config.h"

#define DEBUG_RADTR   NO
#define NO_CHEMISTRY  NO
#define SET_RADSOURCE NO
#define CHECK_NAN     NO
! parameters for calculation

!----------------------------------------------------------------------------------------
! Module for radiation transfer using M1 closer 
!
!----------------------------------------------------------------------------------------
module radtr
  use unit
  use grid
  use modelParameter
  use radtr_chem
  use parameter, only :  Pi, Pi4, Pi4i
  use radiationSource
  use mpilib
  implicit none
  private

  ! label 
  integer, parameter :: NER = 0, NFX = 1, NFY = 2, NFZ = 3

  real(kind=DBL_KIND), save :: tsum, tend, dt
  integer, save :: CurrentLevel
  integer, save :: CurrentIndex       ! CurrentGridId = Wlist(CurrentIndex)
  logical, save :: loop_end

  integer,save,private :: STEP_MODE
  ! 作業変数リスト for timestep
  integer,save,dimension(Gidmin:Gidmax) :: Ulist
  integer,save :: ListMax, num_step
  ! 流速
  real(kind=DBL_KIND), save,dimension(MX:MZ,NER:NFZ,ARRAYSIZE_IJKGH) :: F ! flux
  real(kind=DBL_KIND), save :: cl_til, cl_til2, cl_th

  real(kind=DBL_KIND), save, dimension(:,:), allocatable :: lambda1, lambda4
  type(t_rs_info),pointer, save :: rs_info

#ifdef EXTERNALFORCE
  real(kind=DBL_KIND), save,dimension(MX:MZ,ARRAYSIZE_IJKGH) :: Frad ! 輻射圧
#endif
  ! number of components
  integer, parameter :: ncom = NUM_COMP_TRANSFER
  integer, dimension(NER:NFZ,ncom), save :: nrads


    !----------EO_added----------!
    integer, save :: loop_countmax

    real(kind=DBL_KIND), parameter :: pi_circ  = 3.14159265
    real(kind=DBL_KIND), save :: del_V
    real(kind=DBL_KIND), save :: norm_fun
    real(kind=DBL_KIND), save, dimension(:), allocatable :: log_cosdisk, log_posx, log_posy, log_posz
    !----------EO_added----------!



  ! -----------------------------------
  ! index of each radiation components
  ! -----------------------------------
#ifdef M1CLOSER_EUV_TRANSFER 
    integer, save :: index_neuv
#endif
#ifdef M1CLOSER_FUV_TRANSFER
    integer, save :: index_nfuv
  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    integer, save :: index_nfuvco, index_nfuvd
  #endif
#endif
#ifdef M1CLOSER_IR_TRANSFER
    integer, save :: index_nir
#endif

#ifdef M1CLOSER_EUV_TRANSFER
  type(rad_mean), save :: rsm
#endif

  ! --------------------- for boundary ---------------------------!
  integer,parameter :: NBOUNDARY =  2*(MZ-MX+1) ! 境界面の個数（３次元では６面）
  integer,parameter :: NIL = 0                  ! 方向コード
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  logical,save :: BoolInitialized = .FALSE.     ! このモジュールが初期化されているか
  ! グリッドの配置
  integer,save :: Igmin(Lmin:Lmax), Jgmin(Lmin:Lmax), Kgmin(Lmin:Lmax)
  integer,save :: Igmax(Lmin:Lmax), Jgmax(Lmin:Lmax), Kgmax(Lmin:Lmax)
  !---------------------------------------------------------------!

  public :: radtr_moment 

contains

  !---------------------------------------------------------------
  ! M1 closerに関連する物理量等を初期化
  !---------------------------------------------------------------

  subroutine radtr_init(dt_hyd)
   
    logical, save :: bool_radinit = .false.
    real(kind=DBL_KIND),intent(IN) :: dt_hyd
    integer :: n, ierr

#ifdef SOLVER_RADTR_HLL
    integer,parameter :: FH = 23
    integer :: i, j, ii, kk, err
    character(100) :: path2hll = "./rt/hll_evals.list"
    real(kind=DBL_KIND) :: dummy
#endif
    integer :: nsource_glob, isrc
    real(kind=DBL_KIND) :: mskr, mskr2
#ifdef M1CLOSER_EUV_TRANSFER
    real(kind=DBL_KIND) :: alpha_EUV, heat_EUV, sig_EUV, sig_FUV, tSion, tlum_euv &
      , tlum_fuv, Sion, Sfuv, lumeuv, lumfuv, alpha_OII, tSfuv
    real(kind=DBL_KIND),parameter :: chi_H = 13.6d0!Hの電離エネルギー (in eV)
#endif

    if(.not. bool_radinit) then
      cl_til = MP_Ctil / Unit_v ! reduced speed of light
      cl_til2= cl_til*cl_til
      cl_th  = cl_til*0.5d0
      bool_radinit = .true.



      !---------------------
      ! read table for hll
      !---------------------
#ifdef SOLVER_RADTR_HLL
      allocate(lambda1(0:100,0:100)) ; allocate(lambda4(0:100,0:100))

      if(get_myrank() == PRIMARY_RANK) then
        open(FH, file=path2hll,status='old',iostat=err)
        read(FH, *) 
        do i = 0, 100
          do j = 0, 100
            read(FH,*,iostat=err)ii,kk,lambda1(ii,kk),dummy,dummy,lambda4(ii,kk)
          enddo
        enddo
      endif

      call mpi_bcast(lambda1, size(lambda1), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
      call mpi_bcast(lambda4, size(lambda4), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
#endif

      !---------------------
      ! set components num
      !---------------------

      n = 0
#ifdef M1CLOSER_EUV_TRANSFER 
      n = n+1
      index_neuv   = n
      nrads(NER,n) = MER 
      nrads(NFX,n) = MFRX
      nrads(NFY,n) = MFRY
      nrads(NFZ,n) = MFRZ
#endif
#ifdef M1CLOSER_FUV_TRANSFER
      n = n+1
      index_nfuv   = n
      nrads(NER,n) = MEF
      nrads(NFX,n) = MFRFX
      nrads(NFY,n) = MFRFY
      nrads(NFZ,n) = MFRFZ

  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
      n = n+1
      index_nfuvco = n
      nrads(NER,n) = MECOF
      nrads(NFX,n) = MFCORFX
      nrads(NFY,n) = MFCORFY
      nrads(NFZ,n) = MFCORFZ

      n = n+1
      index_nfuvd  = n
      nrads(NER,n) = MEDUSTF
      nrads(NFX,n) = MFDUSTRFX
      nrads(NFY,n) = MFDUSTRFY
      nrads(NFZ,n) = MFDUSTRFZ
  #endif

#endif
#ifdef M1CLOSER_IR_TRANSFER
      n = n+1
      index_nir    = n
      nrads(NER,n) = MEIR 
      nrads(NFX,n) = MFIRX
      nrads(NFY,n) = MFIRY
      nrads(NFZ,n) = MFIRZ
#endif

#if DEBUG_RADTR == YES
      ! -----------------------------------------------
      if (.not. (n == ncom)) then
        if(get_myrank() == PRIMARY_RANK) then
          print *, "ncom not equal total components"
        endif
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        stop
      endif
      ! -----------------------------------------------
#endif          


    endif

    tsum     = 0.d0 ! 輻射輸送における時間推進の合計値
    num_step = 0
    tend     = dt_hyd
    loop_end = .false.

    !----------------------------
    ! get information of RS
    !----------------------------
    call rs_GetSourceInfo(rs_info)

    !--- get mskr --------------
#ifdef SUB_GRID_MODEL_DIRECTLIGHT
    nsource_glob = rs_info%nsource
    do isrc = 0, nsource_glob -1
      mskr = rs_info%mskr(isrc)
      call mpi_allreduce( mskr, mskr2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      rs_info%mskr = mskr2
    enddo
    !if(get_myrank() == PRIMARY_RANK) write(*,'(a, 1p1e13.5)'), "mskr", rs_info%mskr
#endif
    
    !------------------------------
    ! get mean value of rad source 
    !------------------------------
#ifdef M1CLOSER_EUV_TRANSFER
    nsource_glob = rs_info%nsource
    rsm%num_rad = nsource_glob

    ! 初期化 ---------
    tSion     = 0.d0
    tSfuv     = 0.d0
    alpha_EUV = 0.d0
    alpha_OII = 0.d0
    heat_EUV  = 0.d0
    tlum_euv  = 0.d0
    tlum_fuv  = 0.d0
    sig_EUV   = 0.d0
    sig_FUV   = 0.d0
    ! ----------------
    
    do isrc = 0, nsource_glob -1
      Sion      = rs_info%x_euv(isrc) *rs_info%lum(isrc)
      Sfuv      = rs_info%x_fuv(isrc) *rs_info%lum(isrc)
      

      lumeuv    = rs_info%lumeuv(isrc)*rs_info%lum(isrc) 
      lumfuv    = rs_info%lumfuv(isrc)*rs_info%lum(isrc) 

      tSion     = tSion + Sion
      tSfuv     = tSfuv + Sfuv

      alpha_EUV = alpha_EUV + rs_info%alpha_euv(isrc)*Sion
      heat_EUV  = heat_EUV  + rs_info%heat_euv(isrc) *Sion 

      alpha_OII = alpha_OII + rs_info%alpha_euv(isrc)*Sion*rs_info%rOII(isrc)

      tlum_euv  = tlum_euv + lumeuv 
      sig_EUV   = sig_EUV  + rs_info%sig_euv(isrc)*lumeuv  

      tlum_fuv  = tlum_fuv + lumfuv
      sig_FUV   = sig_FUV  + rs_info%sig_fuv(isrc)*lumfuv  
    enddo

    if(tSion > 0.d0) then
      rsm%rOII      = alpha_OII / alpha_EUV

      rsm%alpha_EUV = alpha_EUV / tSion
      rsm%heat_EUV  = heat_EUV  / tSion

      rsm%sig_EUV   = sig_EUV   / tlum_euv
      rsm%sig_FUV   = sig_FUV   / tlum_fuv   
    else
      rsm%rOII      = 0.d0
      rsm%alpha_EUV = 0.d0
      rsm%heat_EUV  = 0.d0
      rsm%sig_EUV   = 0.d0
      rsm%sig_FUV   = 0.d0  
    endif

    !rsm%erg_EUV = rsm%heat_EUV+chi_H*cgs_eV ![erg]

    ! energy of one EUV photon
    rsm%erg_EUV = max(tlum_euv/max(tSion, 1.d-30), 13.6d0*cgs_ev)
    rsm%erg_FUV = max(tlum_fuv/max(tSfuv, 1.d-30), 11.174d0*cgs_ev)   ! erg s^-1 / s^-1 = erg   energy per one photon



      !if(get_myrank() == PRIMARY_RANK) write(*,'(a, 1p5e13.5)')  &
      !  , "alpha, heat, sig_EUV, sig_FUV, rOII", rsm%alpha_EUV, rsm%heat_EUV, rsm%sig_EUV &
      !  , rsm%sig_FUV, rsm%rOII

#endif




  end subroutine radtr_init




  !---------------------------------------------------------------
  ! 非等方輻射場の各セルの角度 (EO_added)
  !---------------------------------------------------------------
  subroutine anisotropy_init_new

      use overBlockCoordinates
      use modelParameter
      logical, save :: bool_aniinit = .false.

      integer :: level
      integer :: pos_loopmax
      integer :: loop_count
      integer :: pos_iloop, pos_jloop, pos_kloop

      real(kind=DBL_KIND) :: norm_r2cos
      real(kind=DBL_KIND) :: pos_R2, pos_R
      real(kind=DBL_KIND) :: disk_inc, cos_disk_cell
      real(kind=DBL_KIND), dimension(MX:MZ) :: pos_inj, hmax, disk_axis, pos_vec





      if(.not. bool_aniinit) then
        bool_aniinit = .true.

        level   = MP_Lmax0
        hmax    = CellWidth(:,level)
 
        disk_inc      = pi_circ / 2.d0 !edge-on
      !  disk_inc      = 0.0             !pole-on
        disk_axis(MX) = cos(disk_inc)
        disk_axis(MY) = sin(disk_inc)
        disk_axis(MZ) = 0.d0

      !  h(MX) = 2.d0*MP_Boxsize/(64.0* 2.d0**level) 
      !  h(MY) = 2.d0*MP_Boxsize/(64.0* 2.d0**level) 
      !  h(MZ) = 2.d0*MP_Boxsize/(64.0* 2.d0**level)

        pos_loopmax = 14


        !----------for allocate-----------!
        loop_count = 0
        do pos_iloop=0, pos_loopmax-1
        do pos_jloop=0, pos_loopmax-1
        do pos_kloop=0, pos_loopmax-1

           if (pos_iloop<(pos_loopmax/2)) then
              pos_inj(MX) = 0.d0 + ((pos_iloop+1.d0)*hmax(MX)-hmax(MX)*0.5)
           else
              pos_inj(MX) = 0.d0 - ((pos_iloop+1.d0-(pos_loopmax/2))*hmax(MX)-hmax(MX)*0.5)
           end if

           if (pos_jloop<(pos_loopmax/2)) then
              pos_inj(MY) = 0.d0 + ((pos_jloop+1.d0)*hmax(MY)-hmax(MY)*0.5)
           else
              pos_inj(MY) = 0.d0 - ((pos_jloop+1.d0-(pos_loopmax/2))*hmax(MY)-hmax(MY)*0.5)
           end if

           if (pos_kloop<(pos_loopmax/2)) then
              pos_inj(MZ) = 0.d0 + ((pos_kloop+1.d0)*hmax(MZ)-hmax(MZ)*0.5)
           else
              pos_inj(MZ) = 0.d0 - ((pos_kloop+1.d0-(pos_loopmax/2))*hmax(MZ)-hmax(MZ)*0.5)
           end if

           pos_R2 = pos_inj(MX)**2.d0 + pos_inj(MY)**2.d0 + pos_inj(MZ)**2.d0
           pos_R  = sqrt(pos_R2)


           if(pos_R<(hmax(MX)*7.0-hmax(MX)*0.5)) then
              loop_count = loop_count+1
           end if

        enddo
        enddo
        enddo

        loop_countmax = loop_count
 

        allocate(log_posx(0:loop_countmax))
        allocate(log_posy(0:loop_countmax))
        allocate(log_posz(0:loop_countmax))
        allocate(log_cosdisk(0:loop_countmax))

        log_posx(:)    =0.d0
        log_posy(:)    =0.d0
        log_posz(:)    =0.d0
        log_cosdisk(:) =0.d0


      !-----------for aloocate-------------!


        norm_r2cos = 0.d0
        loop_count = 0

        do pos_iloop=0, pos_loopmax-1
        do pos_jloop=0, pos_loopmax-1
        do pos_kloop=0, pos_loopmax-1

           if (pos_iloop<(pos_loopmax/2)) then
              pos_inj(MX) = 0.d0 + ((pos_iloop+1.d0)*hmax(MX)-hmax(MX)*0.5)
           else
              pos_inj(MX) = 0.d0 - ((pos_iloop+1.d0-(pos_loopmax/2))*hmax(MX)-hmax(MX)*0.5)
           end if

           if (pos_jloop<(pos_loopmax/2)) then
              pos_inj(MY) = 0.d0 + ((pos_jloop+1.d0)*hmax(MY)-hmax(MY)*0.5)
           else
              pos_inj(MY) = 0.d0 - ((pos_jloop+1.d0-(pos_loopmax/2))*hmax(MY)-hmax(MY)*0.5)
           end if

           if (pos_kloop<(pos_loopmax/2)) then
              pos_inj(MZ) = 0.d0 + ((pos_kloop+1.d0)*hmax(MZ)-hmax(MZ)*0.5)
           else
              pos_inj(MZ) = 0.d0 - ((pos_kloop+1.d0-(pos_loopmax/2))*hmax(MZ)-hmax(MZ)*0.5)
           end if

           pos_R2 = pos_inj(MX)**2.d0 + pos_inj(MY)**2.d0 + pos_inj(MZ)**2.d0
           pos_R  = sqrt(pos_R2)



           if(pos_R<(hmax(MX)*7.0-hmax(MX)*0.5)) then
              loop_count = loop_count+1

              pos_vec(MX)    = pos_inj(MX)/pos_R
              pos_vec(MY)    = pos_inj(MY)/pos_R
              pos_vec(MZ)    = pos_inj(MZ)/pos_R
              cos_disk_cell  = dot_product(pos_vec,disk_axis)

              log_posx(loop_count)    = pos_inj(MX)
              log_posy(loop_count)    = pos_inj(MY)
              log_posz(loop_count)    = pos_inj(MZ)

              !----------disk equatorial plane----------!
              if ((pos_inj(MY) < 0.d0 .and. pos_jloop==7 ) .or. (pos_inj(MY) > 0.d0 .and. pos_jloop==0 )) then !edge-on
           !   if ((pos_inj(MX) < 0.d0 .and. pos_iloop==7 ) .or. (pos_inj(MX) > 0.d0 .and. pos_iloop==0 )) then !pole-on

                 cos_disk_cell = 0.d0
              end if
              !---------disk equatorial plane-----------!

              log_cosdisk(loop_count) = abs(cos_disk_cell)
              norm_r2cos = norm_r2cos + (2.0*abs(cos_disk_cell) / pos_R2)


           end if

        enddo
        enddo
        enddo

        norm_fun = 1.d0 / norm_r2cos

        end if

  end subroutine anisotropy_init_new
 

  !---------------------------------------------------------------
  ! 流体の時間ステップまで輻射輸送を時間推進
  !---------------------------------------------------------------
  subroutine radtr_moment(dt_hyd)
    use mpilib
    real(kind=DBL_KIND),intent(IN) :: dt_hyd ! time step lengh of hydro
    real(kind=DBL_KIND) ::  time_radtr
    integer :: time_prev, time_cur, time_rat, time1, time2, time_chem, time_mpi

#ifndef RADTR_DIRECT
       call system_clock(time_prev)
       if(get_myrank() == PRIMARY_RANK) print '(/, A, /, A)', "start RADTR" &
         , "(RADTR ) ---------------------------------------------------------------------------"
#endif
    
    ! ---------------------
    !    initialize time
    ! ---------------------
    time_chem = 0
    time_mpi  = 0

    ! ----------------------
    ! 輻射輸送について初期化
    ! ----------------------
    call radtr_init(dt_hyd)    

    call anisotropy_init_new

    ! ---------------------
    ! 輻射圧について初期化
    ! ---------------------
#ifdef EXTERNALFORCE
  #ifndef RADTR_DIRECT
    call radforce_init
  #endif
#endif

    ! ---------------------
    ! 輻射について時間推進
    ! ---------------------
    do   

      ! cfl condition for photons
      call cfl_condition_radrt
!#if DEBUG_RADTR == YES
      if(get_myrank() == PRIMARY_RANK) write(*,'(a,I5,1p4e13.5)') "step, tsum, dt, dt_hyd", num_step, dt_hyd/dt, tsum, dt, dt_hyd
!#endif

#ifdef SKIP_RADTR_BEFORE_STARFORM  
      if (rs_info%nsource > 0 ) then
#endif
        !-------------------
        ! injection step
        !-------------------
        call injection_step      
        !call injection_step_test


        !------------------
        ! transport step
        !------------------
        call timestep_radtr

#ifdef SKIP_RADTR_BEFORE_STARFORM  
      else
        call set_Erad_nRS  ! set EIR = CMB energy density
      endif
#endif



      !------------------
      ! chemistry step
      !------------------
#ifdef M1CLOSER_EUV_TRANSFER
  #if NO_CHEMISTRY == NO
      call system_clock(time1)
      call ch_radtrchem(dt, dt_hyd, rsm)
      call system_clock(time2)
      time_chem = time_chem + time2-time1
      

  #endif
#else
  #ifdef M1CLOSER_IR_TRANSFER
      call radtr_IRdust(dt, dt_hyd)
  #endif
#endif


#ifdef SKIP_RADTR_BEFORE_STARFORM  
      if (rs_info%nsource > 0 ) then
        call set_Erad_nRS  ! set EIR = CMB energy density
      endif
#endif

      !-------------------
      ! converge all level
      !-------------------
      call system_clock(time1)
      call converge_alllevel
      call system_clock(time2)
      time_mpi = time_mpi + time2-time1



      ! --------------
      ! loop end
      ! --------------
      if (loop_end) exit


    enddo

    !---------------------
    ! no rad force @ sink 
    !---------------------
!#ifdef EXTERNALFORCE
!    call no_radforce_inside_sink
!#endif

#ifndef RADTR_DIRECT
       call system_clock(time_cur, time_rat)
       time_radtr = (time_cur - time_prev)/dble(time_rat)
       if(get_myrank() == PRIMARY_RANK) print '(A,/, A, (1PE12.4), A, A, (1PE12.4), A, A, (1PE12.4), A/)' &
         , "----------------------------------------------------------------------------------", &
         "TIME for RADTR:", time_radtr, "[s]", "CHEM: ", time_chem/dble(time_rat), "[s]", "MPI:", time_mpi/dble(time_rat), "[s]"
#endif


  end subroutine radtr_moment

  
  ! ----------------------------------------------------------------------------------------
  ! find dtime for radation transfer
  ! ----------------------------------------------------------------------------------------

  subroutine cfl_condition_radrt 
    use grid
    use mpilib

    real(kind=DBL_KIND) :: dtlocal, dtcfl, dtbuf
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    integer :: lev, n, gid

#ifdef SKIP_RADTR_BEFORE_STARFORM  
    if (rs_info%nsource > 0 ) then
#endif
      !-------------------------------------
      ! 全レベルで時間幅(dtlocal)を求める
      !-------------------------------------
      dtlocal = huge(dtlocal)

      lev = LevelMax ! 最高levelでのcfl条件
      h   = CellWidth( :, lev )

      do n=MX,MZ
        dtcfl = CFLfac_radtr * h(n) / (3.d0*cl_til)
        dtlocal = min(dtlocal, dtcfl)
      enddo
      
      call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
      dtlocal = dtbuf

      !-------------------------------
      ! step sum & final step
      !-------------------------------
      num_step = num_step + 1
      if(int(MAX_STEP_RADTR_TO_HYDRO)==num_step) then
        dtlocal = tend - tsum
        loop_end = .true.
      else
        if(tsum + dtlocal > tend) then
          dtlocal  = tend - tsum
          loop_end = .true.
        endif
      endif

#ifdef SKIP_RADTR_BEFORE_STARFORM  
    else
      !-------------------------------
      ! step sum & final step
      !-------------------------------
      num_step = num_step + 1
      dtlocal = tend - tsum
      loop_end = .true.
    endif
#endif

    !print *, dtlocal*Unit_yr, tend*Unit_yr

    !-----------------
    ! sum of timestep
    !-----------------
    dt   = dtlocal
    tsum = tsum + dtlocal

    if (dtlocal < 0.d0) then
      write(*,*) 'negative dtlocal'
      stop
    endif

  
  end subroutine cfl_condition_radrt


  !------------------------------------------------------------------------------------------
  ! injection step 
  !------------------------------------------------------------------------------------------

  subroutine injection_step

    use radiationSource
    use overBlockCoordinates
    use modelParameter, only : MP_spRadius_cell

    integer :: level, n, gid

    !info for radiation source
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(MX:MZ) :: ijkg
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, hmax
    real(kind=DBL_KIND) :: Ndot_ion, IR_rate, absrate_Dir

    !----------EO_added----------!
    integer :: pos_loop
    real(kind=DBL_KIND):: R_polor
    real(kind=DBL_KIND):: cos_disk
    !----------EO_added----------!


#ifdef M1CLOSER_EUV_TRANSFER 
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad_euv
    !----------EO_added----------!
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_x_euv
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_y_euv
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_z_euv
    !----------EO_added----------!
#endif
#ifdef M1CLOSER_IR_TRANSFER 
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad_ir
    !----------EO_added----------!
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_x_ir
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_y_ir
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_z_ir
    !----------EO_added----------!
#endif
#ifdef M1CLOSER_FUV_TRANSFER
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad_fuv
    !----------EO_added----------!
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_x_fuv
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_y_fuv
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_z_fuv
    !----------EO_added----------! 
  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad_fuv_CO
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad_fuv_DUST
    !----------EO_added----------!
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_x_fuv_CO
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_y_fuv_CO
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nflux_z_fuv_CO

    real(kind=DBL_KINd), dimension(:,:,:), pointer :: Nflux_x_fuv_DUST
    real(kind=DBL_KINd), dimension(:,:,:), pointer :: Nflux_y_fuv_DUST
    real(kind=DBL_KINd), dimension(:,:,:), pointer :: Nflux_z_fuv_DUST
    !----------EO_added----------! 
  #endif
#endif
    logical :: isNotFinite

   
       level = MP_Lmax0
       isrc  = 0
       hmax    = CellWidth(:,level)
       del_V   = hmax(MX)*hmax(MY)*hmax(MZ)



       do pos_loop =1, loop_countmax

         pos(MX)  = log_posx(pos_loop)
         pos(MY)  = log_posy(pos_loop)
         pos(MZ)  = log_posz(pos_loop)
         cos_disk = log_cosdisk(pos_loop)


         call ob_getIjkgridFromCoordPhys(ijkg, level, pos)
         call get_gid_from_ijkgrid(ijkg(MX),ijkg(MY),ijkg(MZ),level,gid,rank)

         if (gid == Undefi) cycle
         if (rank /= get_myrank() ) cycle

         x => get_Xp(gid)
         y => get_Yp(gid)
         z => get_Zp(gid)

        ! point of cell
        k = int((pos(MZ)-z(Kmingh))/hmax(MZ) + 0.5d0)+Kmingh
        j = int((pos(MY)-y(Jmingh))/hmax(MY) + 0.5d0)+Jmingh
        i = int((pos(MX)-x(Imingh))/hmax(MX) + 0.5d0)+Imingh

        !----------EO_added----------!
        k = min(max(k,Kmin),Kmax)
        j = min(max(j,Jmin),Jmax)
        i = min(max(i,Imin),Imax)
        !----------EO_added----------!

        R_polor = sqrt(x(i)**2.d0 + y(j)**2.d0 + z(k)**2.d0)



        ! -----------------------------------------
        !       injection of EUV photons 
        ! -----------------------------------------
#ifdef M1CLOSER_EUV_TRANSFER 
        Nrad_euv => get_Ucomp(MER, gid)
        Ndot_ion = rs_info%x_euv(isrc)*rs_info%lum(isrc)*Unit_t/MP_PHON ! s^{-1} => noD
        Nrad_euv(i,j,k) = Nrad_euv(i,j,k) + (Ndot_ion * dt * norm_fun * 2.d0*cos_disk / (del_V * R_polor**2.d0))

        Nflux_x_euv => get_Ucomp(MFRX,gid)
        Nflux_y_euv => get_Ucomp(MFRY,gid)
        Nflux_z_euv => get_Ucomp(MFRZ,gid)
        Nflux_x_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * x(i)/R_polor
        Nflux_y_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * y(j)/R_polor
        Nflux_z_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * z(k)/R_polor


  #if CHECK_NAN == YES
        if(isNotFinite(Nrad_euv(i,j,k))) then
          print *, "EUV becomes nan at injection step", isrc, Ndot_ion, dt, del_V
          print *, rs_info%x_euv(isrc), rs_info%lum(isrc)
          stop
        endif
  #endif
#endif




        ! -----------------------------------------
        !       injection of FUV photons 
        ! -----------------------------------------
#ifdef M1CLOSER_FUV_TRANSFER
        Nrad_fuv => get_Ucomp(MEF, gid)
        Ndot_ion = rs_info%x_fuv(isrc)*rs_info%lum(isrc)*Unit_t/MP_PHON ! s^{-1} => noD
        Nrad_fuv(i,j,k) = Nrad_fuv(i,j,k) + (Ndot_ion * dt * norm_fun * 2.d0*cos_disk / (del_V * R_polor**2.d0))


        Nflux_x_fuv => get_Ucomp(MFRFX,gid)
        Nflux_y_fuv => get_Ucomp(MFRFY,gid)
        Nflux_z_fuv => get_Ucomp(MFRFZ,gid)
        Nflux_x_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *x(i)/R_polor
        Nflux_y_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *y(j)/R_polor
        Nflux_z_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *z(k)/R_polor

  #if CHECK_NAN == YES
        if(isNotFinite(Nrad_fuv(i,j,k))) then
          print *, "FUV becomes nan at injection step", isrc, Ndot_ion, dt, del_V
          print *, rs_info%x_fuv(isrc), rs_info%lum(isrc)
          stop
        endif
  #endif


  #ifdef M1CLOSER_SEPARATE_FUV_TRANS
        Nrad_fuv_CO => get_Ucomp(MECOF, gid)
        Nrad_fuv_DUST => get_Ucomp(MEDUSTF, gid)
        Nrad_fuv_CO(i,j,k) = Nrad_fuv_CO(i,j,k) + (Ndot_ion*dt*norm_fun*2.d0*cos_disk/(del_V *R_polor**2.d0))
        Nrad_fuv_DUST(i,j,k) = Nrad_fuv_DUST(i,j,k) + (Ndot_ion*dt*norm_fun*2.d0*cos_disk/(del_V *R_polor**2.d0))

        Nflux_x_fuv_CO => get_Ucomp(MFCORFX,gid)
        Nflux_y_fuv_CO => get_Ucomp(MFCORFY,gid)
        Nflux_z_fuv_CO => get_Ucomp(MFCORFZ,gid) 
        Nflux_x_fuv_DUST => get_Ucomp(MFDUSTRFX,gid)
        Nflux_y_fuv_DUST => get_Ucomp(MFDUSTRFY,gid)
        Nflux_z_fuv_DUST => get_Ucomp(MFDUSTRFZ,gid)

        Nflux_x_fuv_CO(i,j,k) = Nrad_fuv_CO(i,j,k) * cl_til *x(i)/R_polor
        Nflux_y_fuv_CO(i,j,k) = Nrad_fuv_CO(i,j,k) * cl_til *y(j)/R_polor
        Nflux_z_fuv_CO(i,j,k) = Nrad_fuv_CO(i,j,k) * cl_til *z(k)/R_polor
        Nflux_x_fuv_DUST(i,j,k) = Nrad_fuv_DUST(i,j,k) * cl_til *x(i)/R_polor
        Nflux_y_fuv_DUST(i,j,k) = Nrad_fuv_DUST(i,j,k) * cl_til *y(j)/R_polor
        Nflux_z_fuv_DUST(i,j,k) = Nrad_fuv_DUST(i,j,k) * cl_til *z(k)/R_polor
    #if CHECK_NAN == YES
          if(isNotFinite(Nrad_fuv_CO(i,j,k))) then
            print *, "FUV_CO becomes nan at injection step", isrc, Ndot_ion, dt, del_V
            stop
          endif
    #endif
  #endif


#endif



        ! -----------------------------------------
        !       injection of IR photons 
        ! -----------------------------------------
#ifdef M1CLOSER_IR_TRANSFER 

        Nrad_ir => get_Ucomp(MEIR, gid)
  #ifdef RADTR_DIRECT
    #ifdef SUB_GRID_MODEL_DIRECTLIGHT
        absrate_Dir = max(min(1.d0, 1.d0-rs_info%mskr(isrc)),0.d0)
        IR_rate  =  1.d0-rs_info%lumeuv(isrc)*absrate_Dir-rs_info%lumfuv(isrc)


    #else
        IR_rate  =  1.d0-rs_info%lumeuv(isrc)-rs_info%lumfuv(isrc)
    #endif

  #else

    #if defined(M1CLOSER_FUV_TRANSFER) && defined(M1CLOSER_SEPARATE_FUV_TRANS)
        IR_rate  =  1.d0-rs_info%lumeuv(isrc)  ! 直接光が含まれていない場合はFUV もIRとして計算
    #elif defined(M1CLOSER_FUV_TRANSFER)
        IR_rate  =  1.d0-rs_info%lumeuv(isrc)-rs_info%lumfuv(isrc)
    #else
        IR_rate  =  1.d0-rs_info%lumeuv(isrc)  ! 直接光が含まれていない場合はFUV もIRとして計算
    #endif

  #endif

        IR_rate  =  min(max(IR_rate, 0.d0), 1.d0)
        Ndot_ion =  rs_info%lum(isrc)*IR_rate/MP_hnu_IR*Unit_t/MP_PHON ! s^{-1} => noD
        Nrad_ir(i,j,k) = Nrad_ir(i,j,k) + (Ndot_ion*dt*norm_fun*2.d0*cos_disk/(del_V *R_polor**2.d0))


         Nflux_x_ir => get_Ucomp(MFIRX,gid)
         Nflux_y_ir => get_Ucomp(MFIRY,gid)
         Nflux_z_ir => get_Ucomp(MFIRZ,gid)

         Nflux_x_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *x(i)/R_polor
         Nflux_y_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *y(j)/R_polor
         Nflux_z_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *z(k)/R_polor

  #if CHECK_NAN == YES
        if(isNotFinite(Nrad_ir(i,j,k))) then
          print *, "IR becomes nan at injection step", isrc, Ndot_ion, dt, del_V
          print *, rs_info%lum(isrc), IR_rate
          stop
        endif
  #endif
#endif


       enddo




  end subroutine injection_step

  !------------------------------------------------------------------------------------------
  ! injection step IR 
  !------------------------------------------------------------------------------------------

!!#ifdef M1CLOSER_IR_TRANSFER 
!!  subroutine injection_step_IR
!!
!!    use radiationSource
!!    use overBlockCoordinates
!!
!!    integer :: level, n, gid
!!
!!    !info for radiation source
!!    integer :: i, j, k, isrc, rank
!!    integer :: nsource_glob
!!    integer,dimension(MX:MZ) :: ijkg
!!    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
!!    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, h
!!    real(kind=DBL_KIND) :: Ndot_ion, del_V, IR_rate, absrate_Dir
!!    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad
!!
!!
!!    do level = Lmin, Lmax
!!
!!      ! number of source
!!      nsource_glob = rs_info%nsource
!!
!!      do isrc = 0, nsource_glob -1
!!
!!        ! position of radiation source
!!        pos = rs_info%spos(:,isrc) 
!!
!!        call ob_getIjkgridFromCoordPhys(ijkg, level, pos)
!!        call get_gid_from_ijkgrid(ijkg(MX),ijkg(MY),ijkg(MZ),level,gid,rank)
!!
!!        if (gid == Undefi) cycle
!!        if (rank /= get_myrank() ) cycle
!!
!!        !cordinates
!!        x => get_Xp(gid)
!!        y => get_Yp(gid)
!!        z => get_Zp(gid)
!!        h = CellWidth(:,level)
!!
!!        Nrad => get_Ucomp(MEIR, gid)
!!
!!        k = int((pos(MZ)-z(Kmin))/h(MZ) + 0.5d0)+Kmin
!!        j = int((pos(MY)-y(Jmin))/h(MY) + 0.5d0)+Jmin
!!        i = int((pos(MX)-x(Imin))/h(MX) + 0.5d0)+Imin
!!
!!        ! volume of cell
!!        del_V = h(MX)*h(MY)*h(MZ)
!!        !Ndot_ion = (1.d0-rs_info%lumeuv(isrc))*rs_info%lum(isrc)/hnu_IR*Unit_t/MP_PHON ! s^{-1} => noD
!!
!!
!!#ifdef RADTR_DIRECT
!!  #ifdef SUB_GRID_MODEL_DIRECTLIGHT
!!        absrate_Dir = max(min(1.d0, 1.d0-rs_info%mskr(isrc)),0.d0)
!!        IR_rate  =  1.d0-rs_info%lumeuv(isrc)*absrate_Dir-rs_info%lumfuv(isrc)
!!        !print *, "IR_rate", IR_rate, rs_info%lumeuv(isrc)*absrate_Dir, absrate_Dir
!!  #else
!!        IR_rate  =  1.d0-rs_info%lumeuv(isrc)-rs_info%lumfuv(isrc)
!!  #endif
!!
!!#else
!!        IR_rate  =  1.d0-rs_info%lumeuv(isrc)  ! 直接光が含まれていない場合はFUV もIRとして計算
!!#endif
!!
!!        IR_rate  =  min(max(IR_rate, 0.d0), 1.d0)
!!        Ndot_ion =  rs_info%lum(isrc)*IR_rate/MP_hnu_IR*Unit_t/MP_PHON ! s^{-1} => noD
!!
!!        Nrad(i,j,k) = Nrad(i,j,k) + Ndot_ion/del_V*dt
!!
!!      enddo
!!
!!    enddo
!!
!!  end subroutine injection_step_IR
!!#endif

  !------------------------------------------------------------------------------------------
  ! injection step test 
  !------------------------------------------------------------------------------------------

#ifdef M1CLOSER_EUV_TRANSFER 
  subroutine injection_step_test

    use radiationSource
    use overBlockCoordinates

    integer :: level, n, gid

    !info for radiation source
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(MX:MZ) :: ijkg
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, h
    real(kind=DBL_KIND) :: Ndot_ion, del_V
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad, Frx, Fry, Frz
    real(kind=DBL_KIND) :: x_start, N_const, radius_c, radius
    
    x_start = 0.03*cgs_pc/Unit_l
    !x_start = 0.d0*cgs_au/Unit_l
    N_const = 1.d49*Unit_t/MP_PHON
    radius_c= 0.02*cgs_pc/Unit_l

    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList(n, level) ! gid for U


        !cordinates
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)

        Nrad => get_Ucomp(MER, gid)
        Frx  => get_Ucomp(MFRX,gid)
        Fry  => get_Ucomp(MFRY,gid)
        Frz  => get_Ucomp(MFRZ,gid)

        h   = CellWidth( :, level)

        del_V = h(MX)*h(MY)*h(MZ)

        do k = Kmin, Kmax
          if(z(k)-h(MZ) .le. 0.d0 .and. 0.d0 < z(k)+h(MZ)) then
            do j = Jmin, Jmax

                do i = Imin, Imax

                  if(y(j)-h(MY)*0.5d0 .le. -x_start .and. -x_start < y(j)+h(MY)*0.5d0 ) then
                    if( x(i)-h(MX)*0.5d0 .le. 0.d0 .and. 0.d0 < x(i)+h(MX)*0.5d0) then
                      Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
                      Frx(i,j,k)  = 0.d0 !cl_til * Nrad(i,j,k)/dsqrt(2.d0) !(N_const/del_V*dt)
                      Fry(i,j,k)  = cl_til * Nrad(i,j,k)
                      Frz(i,j,k)  = 0.d0
                    endif
                  endif

                  if(y(j)-h(MY)*0.5d0 .le. 0.d0 .and. 0.d0 < y(j)+h(MY)*0.5d0 ) then
                    if( x(i)-h(MX)*0.5d0 .le. -x_start .and. -x_start < x(i)+h(MX)*0.5d0) then
                      Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
                      Frx(i,j,k)  = cl_til * Nrad(i,j,k) !(N_const/del_V*dt)
                      Fry(i,j,k)  = 0.d0 !cl_til * Nrad(i,j,k)/dsqrt(2.d0) 
                      Frz(i,j,k)  = 0.d0
                    endif
                  endif


                enddo
            enddo
          endif
        enddo
      
        !do k = Kmin, Kmax
        !  if(z(k)-h(MZ) .le. 0.d0 .and. 0.d0 < z(k)+h(MZ)) then
        !    do j = Jmin, Jmax

        !      if(y(j)-h(MY)*0.5d0 .le. -x_start .and. -x_start < y(j)+h(MY)*0.5d0 ) then
        !        do i = Imin, Imax

        !          if( x(i)-h(MX)*0.5d0 .le. x_start .and. x_start < x(i)+h(MX)*0.5d0) then
        !            Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
        !            Frx(i,j,k)  = -cl_til * Nrad(i,j,k)/dsqrt(2.d0) !(N_const/del_V*dt)
        !            Fry(i,j,k)  = cl_til * Nrad(i,j,k)/dsqrt(2.d0) 
        !            Frz(i,j,k)  = 0.d0
        !          endif
        !          if( x(i)-h(MX)*0.5d0 .le. -x_start .and. -x_start < x(i)+h(MX)*0.5d0) then
        !            Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
        !            Frx(i,j,k)  = cl_til * Nrad(i,j,k)/dsqrt(2.d0) !(N_const/del_V*dt)
        !            Fry(i,j,k)  = cl_til * Nrad(i,j,k)/dsqrt(2.d0) 
        !            Frz(i,j,k)  = 0.d0
        !          endif


        !        enddo
        !      endif
        !    enddo
        !  endif
        !enddo

      enddo
    enddo

  end subroutine injection_step_test
#endif


  !------------------------------------------------------------------------------------------
  ! transport step 
  !------------------------------------------------------------------------------------------
  subroutine timestep_radtr 

    integer :: level 

    do level = Lmin, Lmax
      call step_all_grid_radtr(level)
    enddo

  end subroutine timestep_radtr 

  !--------------------------------
  ! initialize timestep subroutine
  !--------------------------------
  subroutine timestep_radtr_init
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
      Ulist(n) = GidList(n,  CurrentLevel )
    enddo
    ListMax = GidListMax( CurrentLevel )  
  end subroutine timestep_radtr_init

  !-------------------------
  ! times step at level 
  !-------------------------
  subroutine step_all_grid_radtr(level)
    integer, intent(IN) :: level
    integer :: nr, ier, ifx, ify, ifz

    CurrentLevel = level
    
    call timestep_radtr_init


    STEP_MODE = PREDICTOR
    call boundary_cond
    !call rescueLev( CurrentLevel ) !とりあえずなしで
    call rescueLev_radtr

    do CurrentIndex = Gidmin, ListMax !各グリッド
      call backup_u_2order
    enddo
    
    ! update
    do CurrentIndex = Gidmin, ListMax !各グリッド
      do nr = 1, ncom
        ier = nrads(NER,nr)
        ifx = nrads(NFX,nr)
        ify = nrads(NFY,nr)
        ifz = nrads(NFZ,nr)
        call get_flux_radtr(ier, ifx, ify, ifz)
      enddo
    enddo

    call rescueLev_radtr

  end subroutine step_all_grid_radtr

  !------------------------
  ! boundary condition
  !-----------------------
  subroutine boundary_cond
    use grid_boundary
    use boundary
    integer :: n
    call boundary_grid( CurrentLevel, STEP_MODE )
    do n = Gidmin, ListMax !各グリッド
      !call boundary_u( Ulist(n), STEP_MODE)
      call boundary_u_radtr(Ulist(n))
    enddo
  end subroutine boundary_cond


  !--------------------------------------------------
  !         boundary condition for radtr
  !               輻射場だけ設定
  !--------------------------------------------------

  subroutine boundary_u_radtr(id)
    use grid
    integer, intent(IN) :: id
    integer :: i,j,k,m
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
#ifdef EXTERNAL_FUV_RAD
    real(kind=DBL_KIND) :: urad_fuv_bg, nrad_fuv_bg
#endif

    call boundary_radtr_init
    call touch_boundary_radtr( id, bool_touch )

#ifdef EXTERNAL_FUV_RAD
    urad_fuv_bg = EXTERNAL_FUV_G0*1.6d-3/MP_Ctil  ! [ erg cm^-3 ]
    nrad_fuv_bg = urad_fuv_bg / rsm%erg_FUV ! number density [cm^-3]
    nrad_fuv_bg = nrad_fuv_bg / MP_PHON * Unit_l3 ! [noD]
#endif

    if (.not. any(bool_touch) ) return
    u => get_Up(id)

    ! i=imin, 自由境界
    if ( bool_touch(NIL) ) then
       do i = Imingh,Imin-1

#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER

          !----------EO_added----------!
         ! if(u(i,-1,-1,MER)/=0.d0) then
         !   print*, "u(MER) is not zero"
         ! end if
          !----------EO_added----------!

          u(i,:,:,MER)  = 0.d0
          u(i,:,:,MFRX) = u(Imin,:,:,MFRX)  ! zero gradient at the boundries
          u(i,:,:,MFRY) = u(Imin,:,:,MFRY) 
          u(i,:,:,MFRZ) = u(Imin,:,:,MFRZ) 
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER

    #ifdef EXTERNAL_FUV_RAD
          u(i,:,:,MEF) = nrad_fuv_bg
    #else
          u(i,:,:,MEF) = 0.d0
    #endif
          u(i,:,:,MFRFX) = u(Imin,:,:,MFRFX) 
          u(i,:,:,MFRFY) = u(Imin,:,:,MFRFY) 
          u(i,:,:,MFRFZ) = u(Imin,:,:,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(i,:,:,MECOF)   = 0.d0
          u(i,:,:,MFCORFX) = u(Imin,:,:,MFCORFX) 
          u(i,:,:,MFCORFY) = u(Imin,:,:,MFCORFY) 
          u(i,:,:,MFCORFZ) = u(Imin,:,:,MFCORFZ) 

          u(i,:,:,MEDUSTF)   = 0.d0
          u(i,:,:,MFDUSTRFX) = u(Imin,:,:,MFDUSTRFX) 
          u(i,:,:,MFDUSTRFY) = u(Imin,:,:,MFDUSTRFY) 
          u(i,:,:,MFDUSTRFZ) = u(Imin,:,:,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(i,:,:,MEIR)   = 0.d0
          u(i,:,:,MFIRX)  = u(Imin,:,:,MFIRX) 
          u(i,:,:,MFIRY)  = u(Imin,:,:,MFIRY) 
          u(i,:,:,MFIRZ)  = u(Imin,:,:,MFIRZ) 
  #endif
#endif
       enddo
    endif
    ! i=imax, 自由境界
    if ( bool_touch(NIR) ) then
       do i = Imax+1,Imaxgh

#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER
          u(i,:,:,MER)  = 0.d0
          u(i,:,:,MFRX) = u(Imax,:,:,MFRX)  ! zero gradient at the boundries
          u(i,:,:,MFRY) = u(Imax,:,:,MFRY)  ! zero gradient at the boundries
          u(i,:,:,MFRZ) = u(Imax,:,:,MFRZ)  ! zero gradient at the boundries
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER

    #ifdef EXTERNAL_FUV_RAD
          u(i,:,:,MEF) = nrad_fuv_bg 
    #else
          u(i,:,:,MEF) = 0.d0
    #endif
          u(i,:,:,MFRFX) = u(Imax,:,:,MFRFX) 
          u(i,:,:,MFRFY) = u(Imax,:,:,MFRFY) 
          u(i,:,:,MFRFZ) = u(Imax,:,:,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(i,:,:,MECOF)   = 0.d0
          u(i,:,:,MFCORFX) = u(Imax,:,:,MFCORFX) 
          u(i,:,:,MFCORFY) = u(Imax,:,:,MFCORFY) 
          u(i,:,:,MFCORFZ) = u(Imax,:,:,MFCORFZ) 

          u(i,:,:,MEDUSTF)   = 0.d0
          u(i,:,:,MFDUSTRFX) = u(Imax,:,:,MFDUSTRFX) 
          u(i,:,:,MFDUSTRFY) = u(Imax,:,:,MFDUSTRFY) 
          u(i,:,:,MFDUSTRFZ) = u(Imax,:,:,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(i,:,:,MEIR)   = 0.d0
          u(i,:,:,MFIRX)  = u(Imax,:,:,MFIRX) 
          u(i,:,:,MFIRY)  = u(Imax,:,:,MFIRY) 
          u(i,:,:,MFIRZ)  = u(Imax,:,:,MFIRZ) 
  #endif
#endif

       enddo
    endif

    ! j=jmin, 自由境界
    if ( bool_touch(NJL) ) then
       do j = Jmingh,Jmin-1


#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER
          u(:,j,:,MER) = 0.d0
          u(:,j,:,MFRX) = u(:,Jmin,:,MFRX)  ! zero gradient at the boundries
          u(:,j,:,MFRY) = u(:,Jmin,:,MFRY)  ! zero gradient at the boundries
          u(:,j,:,MFRZ) = u(:,Jmin,:,MFRZ)  ! zero gradient at the boundries
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER

    #ifdef EXTERNAL_FUV_RAD
          u(:,j,:,MEF) = nrad_fuv_bg  
    #else
          u(:,j,:,MEF) = 0.d0
    #endif
          u(:,j,:,MFRFX) = u(:,Jmin,:,MFRFX) 
          u(:,j,:,MFRFY) = u(:,Jmin,:,MFRFY) 
          u(:,j,:,MFRFZ) = u(:,Jmin,:,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(:,j,:,MECOF)   = 0.d0
          u(:,j,:,MFCORFX) = u(:,Jmin,:,MFCORFX) 
          u(:,j,:,MFCORFY) = u(:,Jmin,:,MFCORFY) 
          u(:,j,:,MFCORFZ) = u(:,Jmin,:,MFCORFZ) 

          u(:,j,:,MEDUSTF) = 0.d0
          u(:,j,:,MFDUSTRFX) = u(:,Jmin,:,MFDUSTRFX) 
          u(:,j,:,MFDUSTRFY) = u(:,Jmin,:,MFDUSTRFY) 
          u(:,j,:,MFDUSTRFZ) = u(:,Jmin,:,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(:,j,:,MEIR) = 0.d0
          u(:,j,:,MFIRX)  = u(:,Jmin,:,MFIRX) 
          u(:,j,:,MFIRY)  = u(:,Jmin,:,MFIRY) 
          u(:,j,:,MFIRZ)  = u(:,Jmin,:,MFIRZ) 
  #endif
#endif


       enddo
    endif
    ! j=jmax, 自由境界
    if ( bool_touch(NJR) ) then
       do j = Jmax+1,Jmaxgh

          !u(:,j,:,:) = u(:,Jmax,:,:)

#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER
          u(:,j,:,MER) = 0.d0
          u(:,j,:,MFRX) = u(:,Jmax,:,MFRX)  
          u(:,j,:,MFRY) = u(:,Jmax,:,MFRY)  
          u(:,j,:,MFRZ) = u(:,Jmax,:,MFRZ)  
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER

    #ifdef EXTERNAL_FUV_RAD
          u(:,j,:,MEF) = nrad_fuv_bg
    #else
          u(:,j,:,MEF) = 0.d0
    #endif
          u(:,j,:,MFRFX) = u(:,Jmax,:,MFRFX) 
          u(:,j,:,MFRFY) = u(:,Jmax,:,MFRFY) 
          u(:,j,:,MFRFZ) = u(:,Jmax,:,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(:,j,:,MECOF)   = 0.d0
          u(:,j,:,MFCORFX) = u(:,Jmax,:,MFCORFX) 
          u(:,j,:,MFCORFY) = u(:,Jmax,:,MFCORFY) 
          u(:,j,:,MFCORFZ) = u(:,Jmax,:,MFCORFZ) 

          u(:,j,:,MEDUSTF) = 0.d0
          u(:,j,:,MFDUSTRFX) = u(:,Jmax,:,MFDUSTRFX) 
          u(:,j,:,MFDUSTRFY) = u(:,Jmax,:,MFDUSTRFY) 
          u(:,j,:,MFDUSTRFZ) = u(:,Jmax,:,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(:,j,:,MEIR) = 0.d0
          u(:,j,:,MFIRX)  = u(:,Jmax,:,MFIRX) 
          u(:,j,:,MFIRY)  = u(:,Jmax,:,MFIRY) 
          u(:,j,:,MFIRZ)  = u(:,Jmax,:,MFIRZ) 
  #endif
#endif
       enddo
    endif

    ! k=kmin, 自由境界
    if ( bool_touch(NKL) ) then
       do k = Kmingh,Kmin-1

          !u(:,:,k,:) = u(:,:,Kmin,:)

#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER
          u(:,:,k,MER) = 0.d0
          u(:,:,k,MFRX) = u(:,:,Kmin,MFRX)  
          u(:,:,k,MFRY) = u(:,:,Kmin,MFRY)  
          u(:,:,k,MFRZ) = u(:,:,Kmin,MFRZ)  
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER

    #ifdef EXTERNAL_FUV_RAD
          u(:,:,k,MEF) = nrad_fuv_bg 
    #else
          u(:,:,k,MEF) = 0.d0
    #endif
          u(:,:,k,MFRFX) = u(:,:,Kmin,MFRFX) 
          u(:,:,k,MFRFY) = u(:,:,Kmin,MFRFY) 
          u(:,:,k,MFRFZ) = u(:,:,Kmin,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(:,:,k,MECOF)   = 0.d0
          u(:,:,k,MFCORFX) = u(:,:,Kmin,MFCORFX) 
          u(:,:,k,MFCORFY) = u(:,:,Kmin,MFCORFY) 
          u(:,:,k,MFCORFZ) = u(:,:,Kmin,MFCORFZ) 

          u(:,:,k,MEDUSTF) = 0.d0
          u(:,:,k,MFDUSTRFX) = u(:,:,Kmin,MFDUSTRFX) 
          u(:,:,k,MFDUSTRFY) = u(:,:,Kmin,MFDUSTRFY) 
          u(:,:,k,MFDUSTRFZ) = u(:,:,Kmin,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(:,:,k,MEIR) = 0.d0
          u(:,:,k,MFIRX)  = u(:,:,Kmin,MFIRX) 
          u(:,:,k,MFIRY)  = u(:,:,Kmin,MFIRY) 
          u(:,:,k,MFIRZ)  = u(:,:,Kmin,MFIRZ) 
  #endif
#endif

       enddo
    endif
    ! k=kmax, 自由境界
    if ( bool_touch(NKR) ) then
       do k = Kmax+1,Kmaxgh

          !u(:,:,k,:) = u(:,:,Kmax,:)

#ifdef RADTR_M1closer  
  #ifdef M1CLOSER_EUV_TRANSFER
          u(:,:,k,MER) = 0.d0
          u(:,:,k,MFRX) = u(:,:,Kmax,MFRX)  
          u(:,:,k,MFRY) = u(:,:,Kmax,MFRY)  
          u(:,:,k,MFRZ) = u(:,:,Kmax,MFRZ)  
  #endif
  #ifdef M1CLOSER_FUV_TRANSFER
  
    #ifdef EXTERNAL_FUV_RAD
          u(:,:,k,MEF) = nrad_fuv_bg  
    #else
          u(:,:,k,MEF) = 0.d0
    #endif
          u(:,:,k,MFRFX) = u(:,:,Kmax,MFRFX) 
          u(:,:,k,MFRFY) = u(:,:,Kmax,MFRFY) 
          u(:,:,k,MFRFZ) = u(:,:,Kmax,MFRFZ) 

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
          u(:,:,k,MECOF)   = 0.d0
          u(:,:,k,MFCORFX) = u(:,:,Kmax,MFCORFX) 
          u(:,:,k,MFCORFY) = u(:,:,Kmax,MFCORFY) 
          u(:,:,k,MFCORFZ) = u(:,:,Kmax,MFCORFZ) 

          u(:,:,k,MEDUSTF)   = 0.d0
          u(:,:,k,MFDUSTRFX) = u(:,:,Kmax,MFDUSTRFX) 
          u(:,:,k,MFDUSTRFY) = u(:,:,Kmax,MFDUSTRFY) 
          u(:,:,k,MFDUSTRFZ) = u(:,:,Kmax,MFDUSTRFZ) 
    #endif
  #endif
  #ifdef M1CLOSER_IR_TRANSFER
          u(:,:,k,MEIR)   = 0.d0
          u(:,:,k,MFIRX)  = u(:,:,Kmax,MFIRX) 
          u(:,:,k,MFIRY)  = u(:,:,Kmax,MFIRY) 
          u(:,:,k,MFIRZ)  = u(:,:,Kmax,MFIRZ) 
  #endif
#endif
       enddo
    endif


  end subroutine boundary_u_radtr

  
  !--------------------------------------------------------------------
  ! グリッドが物理境界に接していれば真。そでなければ偽
  !--------------------------------------------------------------------
  subroutine touch_boundary_radtr( id, bool_touch )
    use grid
    integer,intent(IN) :: id
    logical,dimension(0:NBOUNDARY-1),intent(OUT) :: bool_touch
    integer :: level, ig, jg, kg

    !call boundary_radtr_init
    bool_touch(:) = .FALSE.
    level = get_level(id)
    if ( Igrid(id) == Igmin( level ) ) bool_touch(NIL) = .TRUE.
    if ( Igrid(id) == Igmax( level ) ) bool_touch(NIR) = .TRUE.
    if ( Jgrid(id) == Jgmin( level ) ) bool_touch(NJL) = .TRUE.
    if ( Jgrid(id) == Jgmax( level ) ) bool_touch(NJR) = .TRUE.
    if ( Kgrid(id) == Kgmin( level ) ) bool_touch(NKL) = .TRUE.
    if ( Kgrid(id) == Kgmax( level ) ) bool_touch(NKR) = .TRUE.
  end subroutine touch_boundary_radtr

  !--------------------------------
  ! initialized boundary subroutine
  !--------------------------------
  subroutine boundary_radtr_init
    use grid
    use io_util
    integer,parameter :: ni0 = NGI_BASE, nj0 = NGJ_BASE, nk0 = NGK_BASE
    integer :: level
    logical, save :: BoolInitialized = .False.

    if(BoolInitialized) return
    call print_msg( 'initialize boundary_radtr' )
    BoolInitialized = .TRUE.
    Igmin(:) =  0
    Jgmin(:) =  0
    Kgmin(:) =  0
    do level = Lmin, Lmax
       Igmax(level) = ni0*2**(level-Lmin) -1
#ifdef Emulate_1Dim
       Jgmax(level) = Jgmin(level)
#else !Emulate_1Dim
       Jgmax(level) = nj0*2**(level-Lmin) -1
#endif !Emulate_1Dim
#ifdef EMULATE_2DIM
       Kgmax(level) = Kgmin(level)
#else !EMULATE_2DIM
       Kgmax(level) = nk0*2**(level-Lmin) -1
#endif !EMULATE_2DIM
    enddo
  end subroutine boundary_radtr_init 

  !-------------------------
  ! back up second order
  !------------------------
  subroutine backup_u_2order  
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u, u2
    u => get_up( Ulist(CurrentIndex ) )
    u2 => get_u2orderp( Ulist(CurrentIndex) )
    u2 = u
    U2_StepNumber(CurrentLevel)          = U_StepNumber(CurrentLevel)
    U2_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_2order

  !------------------------
  ! rescue N & Frad
  !------------------------
  subroutine rescueLev_radtr
#ifdef METAL
    use kinzoku, only : fdust_solar
#endif
    integer :: n, nr, gid, i, j, k
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
#ifdef DUST_NOTCONSTANT
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
#endif
    real(kind=DBL_KIND) :: f2, reduc_f2, f2_sq
    real(kind=DBL_KIND),parameter :: Nrad_floor = 1.d-20 ! 100pc cubic に一つ=1.e-8
    logical :: isNotFinite
    

    do n = Gidmin, GidListMax( CurrentLevel )

      gid = GidList(n, CurrentLevel) ! gid for U

    
      !no-radiation components ------------------------------
!#ifdef DUST_NOTCONSTANT
!      rhod => get_Ucomp(MDRHO,gid)
!#endif
!      do k = Kmingh, Kmaxgh
!        do j = Jmingh, Jmaxgh
!          do i = Imingh, Imaxgh
!              if(isNotFinite(rhod(i,j,k))) then
!                rhod(i,j,k) = fdust_solar*MP_Metallicity
!              endif
!          enddo
!        enddo
!      enddo
      ! ------------------------------------------------------

#ifdef DUST_NOTCONSTANT
      rhod => get_Ucomp(MDRHO,gid)

      do k = Kmingh, Kmaxgh; do j = Jmingh, Jmaxgh; do i = Imingh, Imaxgh
        rhod(i,j,k) = max(rhod(i,j,k), 0.d0)
      enddo; enddo; enddo

#endif

      !radiation components ----------------------------------
      do nr = 1, ncom

        Nrad => get_Ucomp(nrads(NER,nr),gid)
        Frx  => get_Ucomp(nrads(NFX,nr),gid)
        Fry  => get_Ucomp(nrads(NFY,nr),gid)
        Frz  => get_Ucomp(nrads(NFZ,nr),gid)

        
        do k = Kmingh, Kmaxgh
          do j = Jmingh, Jmaxgh
            do i = Imingh, Imaxgh


#if CHECK_NAN == YES
              if(isNotFinite(Nrad(i,j,k))) then
                print *, "nrad becomes nan:", nr
                stop
              endif

              if(isNotFinite(Frx(i,j,k))) then
                print *, "frx becomes nan:", nr
                stop
              endif
              if(isNotFinite(Fry(i,j,k))) then
                print *, "fry becomes nan:", nr
                stop
              endif
              if(isNotFinite(Frz(i,j,k))) then
                print *, "frz becomes nan:", nr
                stop
              endif
#endif
              
              ! flux おかしくなるのに対処 (いる?)
              f2    = Frx(i,j,k)**2 + Fry(i,j,k)**2 + Frz(i,j,k)**2
              f2_sq = dsqrt(f2)
              if(f2_sq > (cl_til*Nrad(i,j,k)) .and. f2_sq > 0.d0) then
                 Frx(i,j,k) =  Frx(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
                 Fry(i,j,k) =  Fry(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
                 Frz(i,j,k) =  Frz(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
              endif

              if(Nrad(i,j,k) < 0.0) then ! ゴミ排除
                  !print '(A, 1P4E15.7)', "photon density becomes negative", Nrad(i,j,k),Frx(i,j,k),Fry(i,j,k),Frz(i,j,k)
                  !stop
                  Nrad(i,j,k) = 0.d0
                  Frx(i,j,k)  = 0.d0
                  Fry(i,j,k)  = 0.d0
                  Frz(i,j,k)  = 0.d0
              endif



            enddo
          enddo
        enddo
      enddo
      ! ------------------------------------------------------

    enddo
  end subroutine rescueLev_radtr

  ! ------------------------------------------------------------
  !   set_Erad_nRS: set E=0,F=0 when there is no raiation source 
  !                 in the simulation box
  ! ------------------------------------------------------------
  subroutine set_Erad_nRS  
    
    implicit none
    integer :: level, n, gid, nr
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
    real(kind=DBL_KIND) :: nrad_cmb

    !----------------
    nrad_cmb = cgs_asb*DEF_TCMB**4.d0/MP_hnu_IR/MP_PHON*Unit_l3/MP_Crd ![noD]

    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList(n, level) ! gid for U

        !radiation components ----------------------------------
        do nr = 1, ncom

          Nrad => get_Ucomp(nrads(NER,nr),gid)
          Frx  => get_Ucomp(nrads(NFX,nr),gid)
          Fry  => get_Ucomp(nrads(NFY,nr),gid)
          Frz  => get_Ucomp(nrads(NFZ,nr),gid)

#ifdef M1CLOSER_IR_TRANSFER
          if (nr == index_nir) then
            Nrad(:,:,:) = nrad_cmb 
          else
            Nrad(:,:,:) = 0.d0
          endif
#else
          Nrad(:,:,:) = 0.d0
#endif
          Frx(:,:,:) = 0.d0  
          Fry(:,:,:) = 0.d0  
          Frz(:,:,:) = 0.d0  
        enddo
      enddo
    enddo


  end subroutine set_Erad_nRS  



  !-------------------------
  ! get flux 
  !-------------------------
  subroutine get_flux_radtr(ier, ifx, ify, ifz)

    integer, intent(IN) :: ier, ifx, ify, ifz
    integer :: i,j,k,n
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND), dimension(MX:MZ, MX:MZ, ARRAYSIZE_IJKGH) :: PP ! P_ij
    real(kind=DBL_KIND), dimension(MX:MZ) :: fi
    real(kind=DBL_KIND) :: inv_cn, f2, f2_sq, ff, ff2, chi, fchi1, fchi2, c2Nr
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
#ifdef SOLVER_RADTR_HLL
    real(kind=DBL_KIND), dimension(MX:MZ, ARRAYSIZE_IJKGH) :: lam_max, lam_min
    real(kind=DBL_KIND), dimension(MX:MZ) :: costheta, lam_p, lam_m, lcr_mp, lsa_mp
#endif
!#if DEBUG_RADTR == YES
!    real(kind=DBL_KIND) :: u_old
!#endif

    real(kind=DBL_KIND) :: tap
    logical :: isNotFinite
    


    u => get_up( Ulist(CurrentIndex) )

    h   = CellWidth( :, CurrentLevel)

    !------------------------------
    ! set eddington tensor
    !------------------------------
    do k = Kmingh, Kmaxgh
      do j = Jmingh, Jmaxgh
        do i = Imingh, Imaxgh

        !----------EO_added----------!
       ! if(u(i,j,k,ier)/=0.d0 .and. ier==MER) then
       !   print*, "******u_ier is not zero*******", u(i,j,k,ier)
       ! end if
        !----------EO_added----------!

          if(u(i,j,k,ier) .le. 0.d0) then
              PP(:, :, i, j, k) = 0.d0
#ifdef SOLVER_RADTR_HLL
              ff = 0.d0

              fi(MX) = u(i,j,k,ifx)
              fi(MY) = u(i,j,k,ify)
              fi(MZ) = u(i,j,k,ifz)
              f2 = fi(MX)**2.d0+fi(MY)**2.d0+fi(MZ)**2.d0 
              
              if(f2 > 0.d0) then
                f2_sq = dsqrt(f2)

                do n = MX, MZ
                  costheta(n) = fi(n) / f2_sq
                enddo

              else
                costheta = 0.d0
              endif
#endif
          else
          
            ! fi
            fi(MX) = u(i,j,k,ifx)
            fi(MY) = u(i,j,k,ify)
            fi(MZ) = u(i,j,k,ifz)

            f2 = fi(MX)**2.d0+fi(MY)**2.d0+fi(MZ)**2.d0 
            
            f2_sq = dsqrt(f2)
            ff    = f2_sq/(cl_til*u(i,j,k,ier))
            ff2   = ff**2.d0

            !chi
            chi = dmax1(4.d0-3.d0*ff2, 0.d0)
            chi = (3.d0+4.d0*ff2)/(5.d0+2.d0*dsqrt(chi))
            chi = dmax1(dmin1(1.d0, chi), 1.d0/3.d0)
            
            fchi1 = (3.d0*chi-1.d0)/2d0
            fchi2 = (1.d0-chi)/2d0

            c2Nr  = cl_til2*u(i,j,k,ier)

            !print *,CurrentIndex, i,j,k, u(i,j,k,ier), fi(MX)**2.d0/f2, fi(MY)**2.d0/f2, fi(MZ)**2.d0/f2

            if(f2 > 0.d0) then
              ! eddington tensor
              PP(MX, MX, i, j, k) = (fchi1*fi(MX)**2.d0 /f2 + fchi2)*c2Nr 
              PP(MY, MX, i, j, k) = (fchi1*fi(MY)*fi(MX)/f2        )*c2Nr
              PP(MZ, MX, i, j, k) = (fchi1*fi(MZ)*fi(MX)/f2        )*c2Nr

              PP(MX, MY, i, j, k) = PP(MY, MX, i, j, k)
              PP(MY, MY, i, j, k) = (fchi1*fi(MY)**2.d0 /f2 + fchi2)*c2Nr
              PP(MZ, MY, i, j, k) = (fchi1*fi(MZ)*fi(MY)/f2        )*c2Nr

              PP(MX, MZ, i, j, k) = PP(MZ, MX, i, j, k)
              PP(MY, MZ, i, j, k) = PP(MZ, MY, i, j, k)
              PP(MZ, MZ, i, j, k) = (fchi1*fi(MZ)**2.d0 /f2 + fchi2)*c2Nr
              
              !if(u(i,j,k,ier) > 0.d0 )then
              !  print '(1P10E15.7)', u(i,j,k,ier), PP(MX, MY, i, j, k), PP(MY, MY, i, j, k), PP(MZ, MY, i, j, k) &
              !    , fi(MX), fi(MY), fi(MZ), 1.d0-chi, chi , fchi2
              !endif

#ifdef SOLVER_RADTR_HLL
              do n = MX, MZ
                costheta(n) = fi(n) / f2_sq
              enddo
#endif
            else
              PP(MX, MX, i, j, k) = fchi2*c2Nr 
              PP(MY, MX, i, j, k) = 0.d0
              PP(MZ, MX, i, j, k) = 0.d0

              PP(MX, MY, i, j, k) = PP(MY, MX, i, j, k)
              PP(MY, MY, i, j, k) = fchi2*c2Nr
              PP(MZ, MY, i, j, k) = 0.d0 

              PP(MX, MZ, i, j, k) = PP(MZ, MX, i, j, k)
              PP(MY, MZ, i, j, k) = PP(MZ, MY, i, j, k)
              PP(MZ, MZ, i, j, k) = fchi2*c2Nr

#ifdef SOLVER_RADTR_HLL
              costheta = 0.d0
              ff = 0.d0
#endif
            endif

          endif


          !-------------------
          ! set eigenvals
          !-------------------
#ifdef SOLVER_RADTR_HLL
          ff = dmax1(dmin1(ff,1d0),0d0)
          do n = MX, MZ
            costheta(n)=dmax1(dmin1(costheta(n),1.d0),-1.d0)
          enddo
          do n = MX, MZ
            call get_eigenvals(ff, costheta(n), lam_min(n, i, j, k), lam_max(n, i, j, k)) 
          enddo

          !if(u(i,j,k,ier) > 1.d9*Unit_l**3.d0/1.d49) then
          !  print *, ff, costheta, lam_max(:, i, j, k), lam_min(:, i, j, k)
          !endif
#endif


        enddo
      enddo
    enddo
    !------------------------------------
    ! calculate flux (GLF function)
    !------------------------------------
#ifdef SOLVER_RADTR_HLL
    do k=Kmingh, Kmaxgh-1 
      do j = Jmingh, Jmaxgh-1
        do i = Imingh, Imaxgh-1

          lam_p(MX) = dmax1(0.d0, lam_max(MX, i, j, k), lam_max(MX, i+1, j, k))
          lam_p(MY) = dmax1(0.d0, lam_max(MY, i, j, k), lam_max(MY, i, j+1, k))
          lam_p(MZ) = dmax1(0.d0, lam_max(MZ, i, j, k), lam_max(MZ, i, j, k+1))

          lam_m(MX) = dmin1(0.d0, lam_min(MX, i, j, k), lam_min(MX, i+1, j, k))
          lam_m(MY) = dmin1(0.d0, lam_min(MY, i, j, k), lam_min(MY, i, j+1, k))
          lam_m(MZ) = dmin1(0.d0, lam_min(MZ, i, j, k), lam_min(MZ, i, j, k+1))

          !FLG test
          !lam_p = 1.d0
          !lam_m = -1.d0

          !if(u(i,j,k,ier) > 1.d9*Unit_l**3.d0/1.d49) then
          !  print *, lam_p, lam_m
          !endif

          do n = MX, MZ
            lcr_mp(n) = lam_p(n)*lam_m(n)*cl_til
            lsa_mp(n) = lam_p(n)-lam_m(n)
          enddo        

#define GETFLUX(Fh, Fl, Flp, Nl, Nlp, M) \
  Fh=(lam_p(M)*Fl-lam_m(M)*Flp+lcr_mp(M)*(Nlp-Nl))/lsa_mp(M)

          ! Flux ---------------------------------------------------------------------
          GETFLUX(F(MX,NER,i,j,k), u(i,j,k,ifx),u(i+1,j,k,ifx),u(i,j,k,ier),u(i+1,j,k,ier),MX)
          GETFLUX(F(MY,NER,i,j,k), u(i,j,k,ify),u(i,j+1,k,ify),u(i,j,k,ier),u(i,j+1,k,ier),MY)
          GETFLUX(F(MZ,NER,i,j,k), u(i,j,k,ifz),u(i,j,k+1,ifz),u(i,j,k,ier),u(i,j,k+1,ier),MZ)
          
          GETFLUX(F(MX,NFX,i,j,k),PP(MX,MX,i,j,k),PP(MX,MX,i+1,j,k),u(i,j,k,ifx),u(i+1,j,k,ifx),MX)
          GETFLUX(F(MY,NFX,i,j,k),PP(MY,MX,i,j,k),PP(MY,MX,i,j+1,k),u(i,j,k,ifx),u(i,j+1,k,ifx),MY)
          GETFLUX(F(MZ,NFX,i,j,k),PP(MZ,MX,i,j,k),PP(MZ,MX,i,j,k+1),u(i,j,k,ifx),u(i,j,k+1,ifx),MZ)

          GETFLUX(F(MX,NFY,i,j,k),PP(MX,MY,i,j,k),PP(MX,MY,i+1,j,k),u(i,j,k,ify),u(i+1,j,k,ify),MX)
          GETFLUX(F(MY,NFY,i,j,k),PP(MY,MY,i,j,k),PP(MY,MY,i,j+1,k),u(i,j,k,ify),u(i,j+1,k,ify),MY)
          GETFLUX(F(MZ,NFY,i,j,k),PP(MZ,MY,i,j,k),PP(MZ,MY,i,j,k+1),u(i,j,k,ify),u(i,j,k+1,ify),MZ)

          GETFLUX(F(MX,NFZ,i,j,k),PP(MX,MZ,i,j,k),PP(MX,MZ,i+1,j,k),u(i,j,k,ifz),u(i+1,j,k,ifz),MX)
          GETFLUX(F(MY,NFZ,i,j,k),PP(MY,MZ,i,j,k),PP(MY,MZ,i,j+1,k),u(i,j,k,ifz),u(i,j+1,k,ifz),MY)
          GETFLUX(F(MZ,NFZ,i,j,k),PP(MZ,MZ,i,j,k),PP(MZ,MZ,i,j,k+1),u(i,j,k,ifz),u(i,j,k+1,ifz),MZ)
          ! --------------------------------------------------------------------------
#undef GETFLUX

        enddo
      enddo
    enddo

#else
    do k=Kmingh+1, Kmaxgh-1
      do j = Jmingh+1, Jmaxgh-1
        do i = Imingh+1, Imaxgh-1

          F(MX,NER,i,j,k) = ((u(i-1,j,k,ifx)-u(i+1,j,k,ifx))*0.5d0 &
            +cl_th*( u(i+1,j,k,ier)-2.d0*u(i,j,k,ier)+u(i-1,j,k,ier)))/h(MX)
          F(MY,NER,i,j,k) = ((u(i,j-1,k,ify)-u(i,j+1,k,ify))*0.5d0 &
            +cl_th*( u(i,j+1,k,ier)-2.d0*u(i,j,k,ier)+u(i,j-1,k,ier)))/h(MY)
          F(MZ,NER,i,j,k) = ((u(i,j,k-1,ifz)-u(i,j,k+1,ifz))*0.5d0 &
            +cl_th*( u(i,j,k+1,ier)-2.d0*u(i,j,k,ier)+u(i,j,k-1,ier) ))/h(MZ)

!#if DEBUG_RADTR == YES
!          if(isNotFinite(F(MX,NER,i,j,k)) .or. isNotFinite(F(MY,NER,i,j,k)) .or. isNotFinite(F(MZ,NER,i,j,k))) then
!            print *, "nan appear at flux"
!            print *, F(MX,NER,i,j,k), F(MY,NER,i,j,k), F(MZ,NER,i,j,k)
!            print *, u(i-1,j,k,ifx), u(i+1,j,k,ifx), u(i+1,j,k,ier), u(i,j,k,ier), u(i-1,j,k,ier), h(MX)
!            print *, u(i,j-1,k,ify), u(i,j+1,k,ify), u(i,j+1,k,ier), u(i,j,k,ier), u(i,j-1,k,ier), h(MY)
!            print *, u(i,j,k-1,ifz), u(i,j,k+1,ifz), u(i,j,k+1,ier), u(i,j,k,ier), u(i,j,k-1,ier), h(MZ)
!            print *, cl_th
!            stop
!          endif
!#endif

          F(MX,NFX,i,j,k) = ((PP(MX,MX,i-1,j,k)-PP(MX,MX,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ifx)-2.d0*u(i,j,k,ifx)+u(i-1,j,k,ifx)))/h(MX)
          F(MY,NFX,i,j,k) = ((PP(MY,MX,i,j-1,k)-PP(MY,MX,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ifx)-2.d0*u(i,j,k,ifx)+u(i,j-1,k,ifx)))/h(MY)
          F(MZ,NFX,i,j,k) = ((PP(MZ,MX,i,j,k-1)-PP(MZ,MX,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ifx)-2.d0*u(i,j,k,ifx)+u(i,j,k-1,ifx)))/h(MZ)

          F(MX,NFY,i,j,k) = ((PP(MX,MY,i-1,j,k)-PP(MX,MY,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ify)-2.d0*u(i,j,k,ify)+u(i-1,j,k,ify)))/h(MX)
          F(MY,NFY,i,j,k) = ((PP(MY,MY,i,j-1,k)-PP(MY,MY,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ify)-2.d0*u(i,j,k,ify)+u(i,j-1,k,ify)))/h(MY)
          F(MZ,NFY,i,j,k) = ((PP(MZ,MY,i,j,k-1)-PP(MZ,MY,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ify)-2.d0*u(i,j,k,ify)+u(i,j,k-1,ify)))/h(MZ)

          F(MX,NFZ,i,j,k) = ((PP(MX,MZ,i-1,j,k)-PP(MX,MZ,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ifz)-2.d0*u(i,j,k,ifz)+u(i-1,j,k,ifz)))/h(MX)
          F(MY,NFZ,i,j,k) = ((PP(MY,MZ,i,j-1,k)-PP(MY,MZ,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ifz)-2.d0*u(i,j,k,ifz)+u(i,j-1,k,ifz)))/h(MY)
          F(MZ,NFZ,i,j,k) = ((PP(MZ,MZ,i,j,k-1)-PP(MZ,MZ,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ifz)-2.d0*u(i,j,k,ifz)+u(i,j,k-1,ifz)))/h(MZ)

        enddo
      enddo
    enddo
#endif

    !------------
    ! update u 
    !------------

    do k=Kmin, Kmax 
      do j = Jmin, Jmax
        do i = Imin, Imax

!#if DEBUG_RADTR == YES
!          u_old = u(i,j,k,ier)
!          if(isNotFinite(u(i,j,k,ier))) then
!            print *, "Nan find from previous step", CurrentIndex, ier
!          endif
!#endif


          
#ifdef SOLVER_RADTR_HLL
          
#define CAL_FLUX(M, NM) \
          u(i,j,k,M) = u(i,j,k,M) \
           + ((F(MX,NM,i-1,j,k)-F(MX,NM,i,j,k))/h(MX) \
           +  (F(MY,NM,i,j-1,k)-F(MY,NM,i,j,k))/h(MY) \
           +  (F(MZ,NM,i,j,k-1)-F(MZ,NM,i,j,k))/h(MZ))*dt

          CAL_FLUX(ier,  NER)
          CAL_FLUX(ifx, NFX)
          CAL_FLUX(ify, NFY)
          CAL_FLUX(ifz, NFZ)
          
#else
          !MER
          u(i,j,k,ier) = u(i,j,k,ier) + (F(MX,NER,i,j,k)+F(MY,NER,i,j,k)+F(MZ,NER,i,j,k))*dt


          !MFR
          u(i,j,k,ifx)= u(i,j,k,ifx)+(F(MX,NFX,i,j,k)+F(MY,NFX,i,j,k)+F(MZ,NFX,i,j,k))*dt  
          u(i,j,k,ify)= u(i,j,k,ify)+(F(MX,NFY,i,j,k)+F(MY,NFY,i,j,k)+F(MZ,NFY,i,j,k))*dt
          u(i,j,k,ifz)= u(i,j,k,ifz)+(F(MX,NFZ,i,j,k)+F(MY,NFZ,i,j,k)+F(MZ,NFZ,i,j,k))*dt
#endif

          !if(dabs(u(i,j,k,ier)) > 0.d0) then
          !  print *,  CurrentLevel, CurrentIndex, i,j,k, u(i,j,k,ier), u(i,j,k,ifx), u(i,j,k,ify), u(i,j,k,ifz)
          !endif
!#if DEBUG_RADTR == YES
!          if(isNotFinite(u(i,j,k,ier)) .or. isNotFinite(u(i,j,k,ifx))  .or. isNotFinite(u(i,j,k,ify)) .or. &
!          isNotFinite(u(i,j,k,ifz)))then
!              print *, "Nan is found at  get_flux_radtr"
!              print *, "u_old:", u_old
!              print *, CurrentIndex, ier, u(i,j,k,ier), u(i,j,k,ifx), u(i,j,k,ify), u(i,j,k,ifz), dt
!              print *, F(MX,NER,i,j,k), F(MY,NER,i,j,k), F(MZ,NER,i,j,k)
!              print *, F(MX,NFX,i,j,k), F(MY,NFX,i,j,k), F(MZ,NFX,i,j,k)
!              print *, F(MX,NFY,i,j,k), F(MY,NFY,i,j,k), F(MZ,NFY,i,j,k)
!              print *, F(MX,NFZ,i,j,k), F(MY,NFZ,i,j,k), F(MZ,NFZ,i,j,k)
!          endif
!#endif

        !----------EO_added----------!
       ! if(u(i,j,k,ier)/=0.d0) then
       !   print*, "*****u(i,j,k,ier) is zero after get_flux******"
       ! end if
        !----------EO_added----------!

        enddo
      enddo
    enddo



  end subroutine get_flux_radtr

  !--------------------------------------------
  ! get eigenvals of HLL
  !--------------------------------------------
  subroutine get_eigenvals(f2, omega, lam_min, lam_max)

    real(kind=DBL_KIND), intent(in) :: f2, omega
    real(kind=DBL_KIND) :: lam_min, lam_max
    real(kind=DBL_KIND) :: theta,dd1,dd2,de1,de2,lff,ltt
    integer::ii,jj

    theta=ACOS(omega)
    lff  = f2*1.d2
    ltt  = theta/Pi*1.d2


    ii = min(int(lff), 99)
    jj = min(int(ltt), 99)
    dd1 = lff-dble(ii)
    dd2 = ltt-dble(jj)
    de1 = 1.d0-dd1
    de2 = 1.d0-dd2
    
    lam_min = 0.d0
    lam_min = lam_min + de1*de2*lambda1(ii,jj)
    lam_min = lam_min + dd1*de2*lambda1(ii+1,jj)
    lam_min = lam_min + de1*dd2*lambda1(ii,jj+1)
    lam_min = lam_min + dd1*dd2*lambda1(ii+1,jj+1)

    lam_max = 0.d0
    lam_max = lam_max + de1*de2*lambda4(ii,jj)
    lam_max = lam_max + dd1*de2*lambda4(ii+1,jj)
    lam_max = lam_max + de1*dd2*lambda4(ii,jj+1)
    lam_max = lam_max + dd1*dd2*lambda4(ii+1,jj+1)

  end subroutine get_eigenvals

#ifdef EXTERNALFORCE
  subroutine radforce_init
    
    integer :: level, n, gid
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi

    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )

        fx_hpi => get_Ucomp(MXPI,gid)
        fy_hpi => get_Ucomp(MYPI,gid)
        fz_hpi => get_Ucomp(MZPI,gid)

        fx_hpi(:,:,:) = 0.d0
        fy_hpi(:,:,:) = 0.d0
        fz_hpi(:,:,:) = 0.d0

      enddo
    enddo

  end subroutine radforce_init

  subroutine no_radforce_inside_sink

    use overBlockCoordinates
    use sinkParticle, only : sp_getSinkRadius

    integer :: level, n, gid

    !info for radiation source
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(MX:MZ) :: ijkg
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    real(kind=DBL_KIND) :: SinkRadius, r2 

    SinkRadius = sp_getSinkRadius()


    do level = Lmin, Lmax

      ! number of source
      nsource_glob = rs_info%nsource

      do isrc = 0, nsource_glob -1

        ! position of radiation source
        pos = rs_info%spos(:,isrc) 

        call ob_getIjkgridFromCoordPhys(ijkg, level, pos)
        call get_gid_from_ijkgrid(ijkg(MX),ijkg(MY),ijkg(MZ),level,gid,rank)

        if (gid == Undefi) cycle
        if (rank /= get_myrank() ) cycle

        !cordinates
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)

        fx_hpi => get_Ucomp(MXPI,gid)
        fy_hpi => get_Ucomp(MYPI,gid)
        fz_hpi => get_Ucomp(MZPI,gid)

        do k=Kmingh, Kmaxgh 
          do j = Jmingh, Jmaxgh
            do i = Imingh, Imaxgh
                r2 = (x(i)-pos(MX))**2 + (y(j)-pos(MY))**2 + (z(k)-pos(MZ))**2
                if (r2 < SinkRadius**2) then
                    fx_hpi(i,j,k) = 0.d0; fy_hpi(i,j,k) = 0.d0; fz_hpi(i,j,k) = 0.d0
                endif
            enddo
          enddo
        enddo
      enddo

    enddo

  end subroutine no_radforce_inside_sink

#endif

!  !---------------------------------------------
!  ! update N & F
!  !---------------------------------------------
!
!  subroutine u_update_radtr
!
!    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
!    integer :: i,j,k
!    logical :: isNotFinite
!    real(kind=DBL_KIND),dimension(MX:MZ) :: h
!   
!    
!    u => get_up( Ulist(CurrentIndex) )
!    h   = CellWidth( :, CurrentLevel)
!
!    do k=Kmin, Kmax 
!      do j = Jmin, Jmax
!        do i = Imin, Imax
!          
!#ifdef SOLVER_RADTR_HLL
!          
!#define CAL_FLUX(M, NM) \
!          u(i,j,k,M) = u(i,j,k,M) \
!           + ((F(MX,NM,i-1,j,k)-F(MX,NM,i,j,k))/h(MX) \
!           +  (F(MY,NM,i,j-1,k)-F(MY,NM,i,j,k))/h(MY) \
!           +  (F(MZ,NM,i,j,k-1)-F(MZ,NM,i,j,k))/h(MZ))*dt
!
!          CAL_FLUX(MER,  NER)
!          CAL_FLUX(MFRX, NFX)
!          CAL_FLUX(MFRY, NFY)
!          CAL_FLUX(MFRZ, NFZ)
!          
!#else
!          !MER
!          u(i,j,k,MER) = u(i,j,k,MER) + (F(MX,NER,i,j,k)+F(MY,NER,i,j,k)+F(MZ,NER,i,j,k))*dt
!
!          !MFR
!          u(i,j,k,MFRX)= u(i,j,k,MFRX)+(F(MX,NFX,i,j,k)+F(MY,NFX,i,j,k)+F(MZ,NFX,i,j,k))*dt  
!          u(i,j,k,MFRY)= u(i,j,k,MFRY)+(F(MX,NFY,i,j,k)+F(MY,NFY,i,j,k)+F(MZ,NFY,i,j,k))*dt
!          u(i,j,k,MFRZ)= u(i,j,k,MFRZ)+(F(MX,NFZ,i,j,k)+F(MY,NFZ,i,j,k)+F(MZ,NFZ,i,j,k))*dt
!#endif
!
!          !if(dabs(u(i,j,k,MER)) > 0.d0) then
!          !  print *,  CurrentLevel, CurrentIndex, i,j,k, u(i,j,k,MER), u(i,j,k,MFRX), u(i,j,k,MFRY), u(i,j,k,MFRZ)
!          !endif
!#if DEBUG_RADTR == YES
!          if(isNotFinite(u(i,j,k,MER)) .or. isNotFinite(u(i,j,k,MFRX))  .or. isNotFinite(u(i,j,k,MFRY)) .or. &
!          isNotFinite(u(i,j,k,MFRZ)))then
!              print *, CurrentIndex, u(i,j,k,MER), u(i,j,k,MFRX), u(i,j,k,MFRY), u(i,j,k,MFRZ)
!          endif
!#endif
!
!        enddo
!      enddo
!    enddo
!
!
!  end subroutine u_update_radtr


  ! -----------------------------------------------------------------
  ! convergeを全レベルに対して実行 (化学更新後に使うことを想定, KS ADDED)
  ! -----------------------------------------------------------------
  subroutine converge_alllevel
    use fg2cg
    integer :: level
    do level = LevelMax, Lmin+1, -1
       call fg2cg_u( level )
    end do
  end subroutine converge_alllevel


end module radtr
