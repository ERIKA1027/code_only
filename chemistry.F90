#include "config.h"
#include "chemistry_label.h"

!control debug messages
!#define CH_DEBUG

!-----------------------------------------------------------------------
! subroutine for chemistry
!-----------------------------------------------------------------------
module chemistry
  use grid
  use parameter,only : Pi 
  use unit
  use mpilib
#ifdef METAL_TRANSFER
  use kinzoku, only: fdust_solar
#else
  use kinzoku, only: fdust_solar
#endif

#if MODEL_ART > 0
  use radiationSource
  use primordial,only: CoolSolverExplicit, CoolSolverImplicit, adjust_abundance &
    , yHe, c_H2, dbg_flg_prm, chem_vals, get_xmu
#endif !MODEL_ART    


#if defined(M1CLOSER_IR_TRANSFER ) && defined(METAL_TRANSFER)
  use modelParameter, only: MP_PHON, MP_mu, MP_frac_COsum
#elif defined(M1CLOSER_IR_TRANSFER) 
  use modelParameter, only: MP_PHON, MP_Metallicity, MP_mu, MP_frac_COsum
#elif defined(METAL_TRANSFER)

#else
  use modelParameter, only: MP_Metallicity, MP_mu, MP_frac_COsum
#endif

#if defined(M1CLOSER_IR_TRANSFER) && defined(RADTR_DIRECT)
  use radtr, only: radtr_moment
#endif


  implicit none
  private

#if MODEL_ART > 0
  !parameters
  real(kind=DBL_KIND),parameter :: fdt = 1d0 !dt_sub < fdt * t_chem
  integer,parameter :: n_substep_max = 3    !maximum number of sub-steps

  !variables
  real(kind=DBL_KIND),save :: t_chemcool_local,t_chemcool_glob !min(t_chem, t_cool)
  type(t_rs_info),save,pointer :: rs_info    !info for radiation source

  integer,save :: ifirst = 0

  !---------------- KS DEBUG (BEGIN)---------------!
  ! time_ini: initial time, time_prev: previous time,, time_cur: current time, time_rat: time rate 
  integer :: time_ini,time_prev, time_cur,time_rat 
  integer :: num_cell ! num_cell: number of cells for chem. update
  integer :: num_implicit ! number of cells with implicit chem. update
  !----------------- KS DEBUG (END) ----------------!

  public :: ch_artchem
#endif !MODEL_ART
#ifdef ISOTHERMAL
  public :: ch_reset_temperature
#endif !ISOTHERMAL
#ifdef RADTR_M1closer
  public :: ch_CellChemCool
#endif
  

contains
#if MODEL_ART > 0
  !proceed art and chemistry step
  subroutine ch_artchem(dt_code)
    use art
    real(kind=DBL_KIND),intent(IN) :: dt_code !dt of hydro step in code unit

    integer :: n_substep !流体のstep当たりのsubstep数
    integer :: sstep
    real(kind=DBL_KIND):: dt, dt_sub, t_sub

    !for debug
    real(kind=DBL_KIND):: time_chem, time_chem_min,time_chem_mean,time_chem_max,&
         time_cell,time_cell_min,time_cell_mean,time_cell_max,&
         frac_imp,frac_imp_min,frac_imp_mean,frac_imp_max, dt_code_sub, time_radtr
    

    call system_clock(time_ini) ! 時間計測開始 (KS DEBUG)

#ifdef NO_RADIATION
    n_substep = 1       !輻射無しの時はsubstep切らない
#else !NO_RADIATION
    n_substep = 1       !輻射ありの時はsubstep2回切る
#endif !NO_RADIATION


    !initial time step
    dt = dt_code*Unit_t    
    t_sub=0d0

    dt_sub = dt/n_substep    

    dt_code_sub = dt_code/n_substep

    ! if (ifirst == 0) then
    !    dt_sub = dt/n_substep_max !t_chemcool分からないので小さくとっておく
    !    t_chemcool_local = HUGE(1d0)
    !    ifirst = 1
    ! else
    !    !       dt_sub = MIN(dt,MAX(fdt * t_chemcool_glob,dt/n_substep_max))
    !    dt_sub = dt/n_substep_max ! HD updateで化学タイムスケールが変わることもを考慮して念のため (KS TODO)
    !    t_chemcool_local = HUGE(1d0)
    ! end if


    !光源の情報を取得
    call rs_GetSourceInfo(rs_info)

    ! n_substep = 1                    !ひとまずsubstep切らない (for debug)
    ! n_substep = 10                    !substep10回切ってみる


    !begin loop for ART/Chemistry
    call ch_check_all_cells() !chemical abundanceに変な値が入らないように調整

    ! do sstep = 0, n_substep_max - 1
    do sstep = 0, n_substep - 1
       if(get_myrank() == PRIMARY_RANK) &
            print '(/,A,4((1P1E15.7)))', "(ART_CHEM) dt, dt_sub, t_sub (sec) = ", dt, dt_sub, t_sub

#ifdef RADTR_DIRECT

       !--------------------------- BEGIN ART  ------------------------------!       
       if(get_myrank() == PRIMARY_RANK) &
            print *, 'Begin ART'       

       call system_clock(time_prev) ! 時間計測開始 (KS DEBUG)

       call mpi_barrier(MPI_COMM_WORLD, ierr) !RayTracingする前に、全部のセルの更新が終わっている必要あり

       call art_RayTracing

       !KS DEBUG
       call system_clock(time_cur,time_rat) ! 時間計測
       if(get_myrank() == PRIMARY_RANK) then
            if (rs_info%nsource > 0) then
               print '(A,I4,A,2(1PE12.4,A))', "TIME for ART   =>    myrank: ", myrank, ", time: ",(time_cur-time_prev)/dble(time_rat), &
                      "[s], time/source: ",(time_cur-time_prev)/dble(time_rat)/dble(rs_info%nsource), "[s]"
              end if
         end if
       !------------------------------ END ART  ------------------------------!       

       !-------------------------  BEGIN M1 CLOSER ---------------------------!
  #ifdef M1CLOSER_IR_TRANSFER
       call system_clock(time_prev) ! 時間計測開始 (KS DEBUG)

       if(get_myrank() == PRIMARY_RANK) print '(/, A, /, A)', "start RADTR" &
         , "(RADTR IR)------------------------------------------------------------------------"
       call radtr_moment(dt_code_sub)   

       call system_clock(time_cur, time_rat)
       time_radtr = (time_cur - time_prev)/dble(time_rat)
       if(get_myrank() == PRIMARY_RANK) print '(A,/, A, (1PE12.4), A, /)' &
         , "----------------------------------------------------------------------------------", &
         "TIME for RADTR:", time_radtr, "[s]"
  #endif
#endif
       !----------------------------------------------------------------------!
       
       !--------------------------- BEGIN CHEM  ------------------------------!       
       if(get_myrank() == PRIMARY_RANK) &
            print *, 'Begin chem update'

       ! ---------- initialize (KS DEBUG) ----------!
       call system_clock(time_prev) ! 時間計測開始 (KS DEBUG)
       num_cell = 0
       num_implicit = 0 
       ! --------------- (KS DEBUG) ----------------!       

       call ch_PrimChemistry(dt_sub)

       ! ------------------------------ (KS DEBUG, BEG) ------------------------------- !       
       call system_clock(time_cur,time_rat) ! 時間計測

       time_chem = (time_cur-time_prev)/dble(time_rat) !chemistryにかかった時間
       time_cell = time_chem/dble(num_cell)            !1セル当たりchemistryにかかった時間
       frac_imp = dble(num_implicit)/dble(num_cell)    !implicit coolingが呼ばれた割合
       
       call mpi_allreduce(time_chem, time_chem_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)       

       time_chem_mean = time_chem_mean/NPE
       time_cell_mean = time_cell_mean/NPE
       frac_imp_mean = frac_imp_mean/NPE       
       
       if(get_myrank() == PRIMARY_RANK) then
          print '(A,3(A,1P3E10.2))', "TIME for CHEM (min,mean,max)  => ", &
               "tot time: ",time_chem_min,time_chem_mean,time_chem_max, &
               " [s], time/cell: ",time_cell_min,time_cell_mean,time_cell_max, &
               " [s], frac_imp: ", frac_imp_min,frac_imp_mean,frac_imp_max
       end if
       ! ------------------------------ (KS DEBUG, END) ------------------------------- !       

       !--------------------------- END CHEM  ------------------------------!              

       !update sub-timestep
       t_sub = t_sub + dt_sub
       if(dt-t_sub < TINY(1d0)) exit


       !------------------- KS TODO (substep is currently detemined by n_substep) ---------------!
       ! !new sub-timestep
       ! call ch_GetDtSub(dt,t_sub,dt_sub)
       !-----------------------------------------------------------------------------------------!
    end do

    !KS DEBUG
    call system_clock(time_cur,time_rat) ! 時間計測
    if(get_myrank() == PRIMARY_RANK) then
       print '(A,I4,A,2(1PE12.4,A))',"TIME for ARTCHEM   =>    myrank: ", myrank, ", tot time: ",(time_cur-time_ini)/dble(time_rat),"[s]"
    end if

    
  end subroutine ch_artchem


  !-------------------------------------------------------------------------
  ! Get n_substep from dt and t_chemcool
  !-------------------------------------------------------------------------


  subroutine ch_GetDtSub(dt,t_sub,dt_sub)
    real(kind=DBL_KIND),intent(IN) :: dt,t_sub
    real(kind=DBL_KIND),intent(OUT) :: dt_sub


    call mpi_allreduce( t_chemcool_local, t_chemcool_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)


    if (dt-t_sub < fdt * t_chemcool_glob) then
       dt_sub = dt-t_sub
    else
       dt_sub = MAX(fdt * t_chemcool_glob,dt/n_substep_max)
    end if

    ! !KS DEBUG
    ! print '(A, (1P10E15.7))', &
    !      "dt, t_sub/dt, dt_sub, t_chemcool_glob,t_chemcool_local:",dt, t_sub/dt, dt_sub,t_chemcool_glob,t_chemcool_local

    t_chemcool_local = HUGE(1d0)
  end subroutine ch_GetDtSub




  !-------------------------------------------------------------------------
  ! Primordial-gas chemistry based on Hosokawa-san's code (cf. Hosokawa+16)
  !
  ! SPECIES                                                         
  !    0 : H      1 : H2     2 : e     3 : H+     4 : H-    5 : H2+ 
  !-------------------------------------------------------------------------
  subroutine ch_PrimChemistry(dt)
    use modelParameter,only : MP_Tmin,MP_Tmax

    real(kind=DBL_KIND),intent(IN) :: dt

    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
#ifdef METAL
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Tdust
#endif
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: mu, khmpd
#if MODEL_ART == 1
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: kh2pd
#elif MODEL_ART == 2
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kh2pd
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: kh2pd_dbg
#endif
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem0
    real(kind=DBL_KIND) :: xmu, xNc_H2,  cs, T_K0
    integer :: level, n, gid, i, j, k, ic

    logical :: isNotFinite ! defined in naninf.F90
    logical :: naninf_flg
    integer :: count_output=0, max_count=100 !KS DEBUG
    real(kind=DBL_KIND) :: probe_radius, yco0!KS DEBUG

#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: yco
#endif

#ifdef METAL
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kgfuv, kdph, kdco, krOII
  #ifdef DUST_NOTCONSTANT
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
  #endif
#endif
#ifdef EXTERNALFORCE
    real(kind=DBL_KIND) :: dt_code
#endif
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: vx,vy,vz
#ifdef M1CLOSER_IR_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad_IR
#endif
#ifdef METAL_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
#endif

    type(chem_vals) :: xch 

    ! ----------------------------------------------------------
    type chemsp 
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(NCEHM_MIN:NCEHM_MAX)
    integer :: ichem
    ! ----------------------------------------------------------


    ! insert timestep
    xch%dt = dt


    naninf_flg = .False.

    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )

          ! 子グリッドを持つ場合はchemistryは解かなくてよい (あとでconvergeを呼ぶ)
          if (has_child_grid(gid))   cycle

          !座標
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)


          !変数の取り出し
          rho => get_Ucomp(MRHO,gid)
          p => get_Ucomp(MP,gid)

          ! get chemistry pointer
          do ichem = NCEHM_MIN, NCEHM_MAX
            chemary3(ichem)%y => get_Ucomp(ichem,gid)
          enddo

#ifdef CHEM_MODEL_HF2020
          yco  => get_Ucomp(MCO, gid)
#endif

          khpi => get_Ucomp(MKPI,gid)
          heat_hpi => get_Ucomp(MHPI,gid)

#ifdef METAL
          Tdust => get_Ucomp(MTD,gid)
#endif


#ifdef EXTERNALFORCE
          fx_hpi => get_Ucomp(MXPI,gid)
          fy_hpi => get_Ucomp(MYPI,gid)
          fz_hpi => get_Ucomp(MZPI,gid)
#endif
          vx     => get_Ucomp(MVX,gid)
          vy     => get_Ucomp(MVY,gid)
          vz     => get_Ucomp(MVZ,gid)

          !get kh2pd
#if MODEL_ART == 1
          call GetKH2pdThin(x,y,z,kh2pd) !w/o self-shielding
#elif MODEL_ART == 2
  #ifdef RADTR_DIRECT
          kh2pd => get_Ucomp(MKPD,gid) !w/ self-shielding
    #ifdef METAL
          kgfuv => get_Ucomp(MGFUV ,gid)
          kdph  => get_Ucomp(MDPH  ,gid)  
          kdco  => get_Ucomp(MDPCO ,gid)
          krOII => get_Ucomp(MKPOII,gid)
    #endif
  #endif
#endif

#ifdef M1CLOSER_IR_TRANSFER
          Nrad_IR => get_Ucomp(MEIR, gid)
#endif

#ifdef METAL
  #ifdef DUST_NOTCONSTANT
         rhod => get_Ucomp(MDRHO,gid)
  #endif
#endif

#ifdef METAL_TRANSFER
         mmetal => get_Ucomp(MMET, gid)
#endif

         call GetKHmpdThin(x,y,z,khmpd) !get khmpd

          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax

                   ! ------------------------
                   !      metallicity 
                   ! ------------------------
#ifdef METAL_TRANSFER
                   xch%metal = mmetal(i,j,k)
#else
                   xch%metal = MP_Metallicity
#endif

#ifdef METAL
                   xch%Td = Tdust(i,j,k) 
#endif

                   ! --------------------------
                   !   化学組成変数の詰め替え
                   ! --------------------------
                   do ichem = 0, NCHEM-1
                      xch%ychem(ichem) = chemary3(ichem+NCEHM_MIN)%y(i,j,k)
                   enddo 
#ifdef CHEM_MODEL_HF2020
                   xch%yco = yco(i,j,k)
#endif

                   !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
                   call adjust_abundance(xch%ychem &
#ifdef CHEM_MODEL_HF2020
                      , xch%yco &
#endif
                      , xch%metal)

                   !derived variables
                   xch%nH = rho(i,j,k)*Unit_rho/(MP_mu*cgs_amu)                !水素原子核の数密度
                   xmu    = get_xmu(xch%ychem) ! 平均分子量
                   xch%Tg = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb      !温度 [K]

                   !impose T floor (before chem_update)
                   if (xch%Tg < MP_Tmin*0.999d0 .or. xch%Tg > MP_Tmax*1.001d0) then !表示を抑制するため係数をかける
                      !output warning
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** before chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         !suppress warning
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)


#ifdef NO_IONIZATION
                   xch%heat = 0.d0
                   xch%rHpi = 0.d0
#else
                   if (khpi(i,j,k) > 0.d0) then
                      xch%heat = heat_hpi(i,j,k)/khpi(i,j,k) !heat_hpi/khpi = energy deposit per ionization (KS NOTE)
                      xch%rHpi = khpi(i,j,k)
                   else
                      xch%heat = 0.d0
                      xch%rHpi = 0.d0
                   end if
#endif

                   !---------------------- H2 LINE TRAPPING -----------------------!
                   cs  = sqrt(cgs_kb*xch%Tg/(xmu*cgs_mh))          !等温音速 [cm/s]
                   xch%xlmbdj = cs*sqrt(Pi/(cgs_gc*rho(i,j,k)*Unit_rho))   !ジーンズ長さ [cm]
                   xch%xNcH   = xch%nH*0.5d0*xch%xlmbdj   ! column density of hydrogen 
                   xNc_H2 = xch%ychem(1)*xch%nH*0.5*xch%xlmbdj             !水素分子の柱密度の概算 [cm^-2]
                   !---------------------------------------------------------------!

                   !---------------------- velocity gradient ----------------------!
                   xch%dvdr(MX) = dabs(vx(i+1,j,k)-vx(i-1,j,k))/(2.d0*CellWidth(MX,level))/Unit_t ! [s^{-1}]
                   xch%dvdr(MY) = dabs(vy(i,j+1,k)-vy(i,j-1,k))/(2.d0*CellWidth(MY,level))/Unit_t ! [s^{-1}]
                   xch%dvdr(MZ) = dabs(vz(i,j,k+1)-vz(i,j,k-1))/(2.d0*CellWidth(MZ,level))/Unit_t ! [s^{-1}]
                   !---------------------------------------------------------------!

#ifdef METAL                      

  #ifdef RADTR_DIRECT
                   xch%rgfuv = kgfuv(i,j,k)
                   xch%rdph  = kdph(i,j,k) / (rho(i,j,k)*Unit_rho) ! [erg g^-1 s^-1]
                   xch%rcopd = kdco(i,j,k) 
                   xch%rOII  = krOII(i,j,k)
  #else
                   xch%rgfuv = 0.d0
                   xch%rdph  = 0.d0
                   xch%rcopd = 0.d0
                   xch%rOII  = 0.d0
  #endif


  #ifdef DUST_NOTCONSTANT
                   xch%fd    = rhod(i,j,k) / fdust_solar ! dust abundance including sublimation & metallicity    
  #else
    #ifdef METAL_TRANSFER
                   xch%fd    = mmetal(i,j,k)
    #else
                   xch%fd    = MP_Metallicity
    #endif
  #endif

#else
                   xch%rgfuv = 0.d0
                   xch%rdph  = 0.d0
                   xch%rcopd = 0.d0
                   xch%rOII  = 0.d0
                   xch%fd    = 0.d0
#endif


                   ! set --------------------
#ifdef RADTR_DIRECT
                   xch%rH2pd = kh2pd(i,j,k) ! H2 photodissocation rate
                   xch%rHmpd = khmpd(i,j,k) ! H- photo-dissociation rate
#else
                   xch%rH2pd = 0.d0 ! H2 photodissocation rate
                   xch%rHmpd = 0.d0 ! H- photo-dissociation rate
#endif
                   ! ------------------------

                   ! IR photon density ------
#ifdef M1CLOSER_IR_TRANSFER
                   xch%EradIR = Nrad_IR(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR ! [erg cm^-3]
                   xch%EradIR = max(xch%EradIR, 0.d0)
#else
                   xch%EradIR = 0.d0
#endif
                   ! ------------------------


                   
                   !----------------------------------------- KS DEBUG -----------------------------------------!



                   !KS DEBUG
                   ! probe_radius = -8e1/Unit_au !KS DEBUG

                   ! if (abs(x(i)-probe_radius) < abs(x(5)-x(4))*0.5 .and. x(i) > 0. .and. &
                   !      abs(y(j)) < abs(y(5)-y(4)) .and. y(j) > 0. .and. &
                   !      abs(z(k)) < abs(z(5)-z(4)) .and. z(k) > 0. ) then ! x ~ probe_rad au の値を表示

                   ! if (abs(xch%nH- 3.3438341E+07)/xch%nH < 1d-4) then

                   ! if (abs(xch%Tg-1.7395534E+02)/1.7395534E+02 < 1d-7) then
                   ! if (abs(x(i)+4.1424206d2/Unit_au) < abs(x(5)-x(4))*0.5 .and. &
                   !      abs(y(j)-2.0120329d2/Unit_au) < abs(y(5)-y(4))*0.5  .and. &
                   !      abs(z(k)+1.3019036d2/Unit_au) < abs(z(5)-z(4))*0.5 ) then 


                   !    print '(/,A, (1P10E15.7),/)', "*** KS DEBUG ***  x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au                      
                   !    print '(A, 6I6)', "rank,gid,i,j,k,level = ", get_myrank(), gid, i-lbound(rho,1), j-lbound(rho,2), k-lbound(rho,3), level
                   !    print '(A, (1P10E15.7))', "xch%nH, xch%Tg, rHpi, rH2pd, rHmpd, dt:",&
                   !         xnH,xch%Tg,khpi(i,j,k),kh2pd(i,j,k),khmpd(i,j,k), dt
                   !    print '(A, (1P10E15.7))', "ychem: ", ychem(:)
                   !    print '(A,(1P2E15.7))',    "heat, xNc_H2: ", heat, xNc_H2
                   !    ! print '(A,(1P10E15.7))', 'KS DEBUG 2', p(i,j,k), yhn(i,j,k),yh2(i,j,k),yel(i,j,k),yhp(i,j,k)
                   !    ! print '(A,(1P10E15.7))', 'KS DEBUG 2A', &
                   !    !      p(i,j,k)*Unit_e*cgs_mh*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb, p(i,j,k),Unit_e,cgs_mh,xmu,rho(i,j,k),Unit_rho,cgs_kb
                   !    dbg_flg_prm = 1
                   ! end if
                   !----------------------------------------- KS DEBUG -----------------------------------------!

                   !-------------        check Nan/Inf     -------------!
                   ychem0(:) = xch%ychem(:)
                   yco0 = xch%yco
                   T_K0 = xch%Tg
                   if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                      naninf_flg = .True.
                   end if
                   do ic=0,NCHEM-1
                      if (isNotFinite(xch%ychem(ic))) then
                         naninf_flg = .True.
                      end if
                   end do
                   if (isNotFinite(xch%yco)) then  
                      naninf_flg = .True.     
                   endif
                   if (naninf_flg) then
                      print '(A, (1P8E15.7))', "(chemistry) NaN/Inf found before chem update: ", xch%nH,xch%Tg,xch%ychem(:)
                      print *, "stopping..."
                      stop
                   end if

                   if(isNotFinite(xch%rdph)) then ! 安全装置
                      print *, "dph is Nan and dph = 0.d0 for safety at ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au
                      xch%rdph = 0.d0
                   endif



                   !print '(A, (1P3E15.7))', "nH, Tg, Td", xch%nH, xch%Tg, xch%Td
                   !print '(A, (1P6E15.7))', "ychem", xch%ychem
                   !print '(A, (1P6E15.7))', "yco rHpi rH2pd rHmpd heat rgfuv", xch%yco, xch%rHpi, xch%rH2pd, xch%rHmpd, xch%heat,xch%rgfuv
                   !print '(A, (1P7E15.7))', "rdph, rcopd, rOII, xlmbdj, dt, fd, eradir", xch%rdph, xch%rcopd &
                   !  , xch%rOII, xch%xlmbdj, xch%dt, xch%fd, xch%EradIR


                   !-----------    セルの時間発展 (KS NOTE: heat_hpi/khpi = energy deposit per ionization) -------!
                   call ch_CellChemCool(xch)
                   num_cell = num_cell + 1 ! KS DEBUG
                   !--------------------------------------------------------------------------------------------!

                   !-------------        check Nan/Inf     -------------!
                   if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                      naninf_flg = .True.
                   end if
                   if (isNotFinite(xch%yco)) then  
                      naninf_flg = .True.     
                   endif
                   if(isNotFinite(xch%Td))then
                      naninf_flg = .True.     
                   endif

                   do ic=0,NCHEM-1
                      if (isNotFinite(xch%ychem(ic))) then
                         naninf_flg = .True.
                      end if
                   end do
                   if (naninf_flg) then
                      print '(A, (1P6E15.7))', "(chemistry) NaN/Inf found after chem update: ", xch%ychem(:)
                      print '((1P10E15.7))',  xch%nH,T_K0,khpi(i,j,k),xch%heat,x(i),y(j),z(k), xch%dt
                      print '((1P10E15.7))',  xch%ychem(:)
                      print *, xch%rdph, xch%rgfuv
                      print *, "Td", xch%Td
                      print *, "stopping..."
                      stop
                   end if

                   if(dbg_flg_prm == 1) then
                      print '(A, (1P10E15.7))', "(KS DEBUG) x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au
                      print '(A, (1P10E15.7))', "xnH, T_K (after):", xch%nH,xch%Tg
                      print '(A, (1P10E15.7))', "ychem (after): ", xch%ychem(:)
                      dbg_flg_prm = 0
                      !stop
                   end if

                   !impose T floor (after chem_update)
                   if (xch%Tg < MP_Tmin .or. xch%Tg > MP_Tmax) then
                      !output warning
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** after chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         !suppress warning
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)


                   ! ----------------------------------
                   !     packing updated chemistry 
                   ! ----------------------------------
                   do ichem = 0, NCHEM-1
                      chemary3(ichem+NCEHM_MIN)%y(i,j,k) = xch%ychem(ichem)
                   enddo 
#ifdef CHEM_MODEL_HF2020
                   yco(i,j,k) = xch%yco
#endif

                   xmu = get_xmu(xch%ychem) ! 平均分子量
                   p(i,j,k) = (cgs_kb*xch%Tg)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e      !温度 [K]
#ifdef METAL
                   Tdust(i,j,k) = xch%Td
#endif

                end do
             end do
          end do
       end do
    end do

  contains
    !-------------------------------------------------------------------------
    ! Get H2 photo-dissociation rate [s^-1] assuming optically thin LW flux
    !-------------------------------------------------------------------------
    subroutine GetKH2pdThin(x,y,z,kh2pd)
      real(kind=DBL_KIND),dimension(:),pointer,intent(IN) :: x, y, z
      real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK),intent(OUT) :: kh2pd

      real(kind=DBL_KIND),parameter :: pdis = 0.15

      integer :: isrc,i,j,k
      real(kind=DBL_KIND) :: r2,lum_fuv,fuvflux

      kh2pd(:,:,:) = 0.d0
      do isrc=0, rs_info%nsource-1
         lum_fuv = rs_info%x_fuv(isrc) * rs_info%lum(isrc)
         do i = Imin, Imax
            do j = Jmin, Jmax
               do k = Kmin, Kmax
                  r2 = (x(i) - rs_info%spos(MX,isrc))**2 + &
                       (y(j) - rs_info%spos(MY,isrc))**2 + &
                       (z(k) - rs_info%spos(MZ,isrc))**2
                  fuvflux = lum_fuv / (4.*Pi*r2*Unit_l**2)
                  kh2pd(i,j,k)=kh2pd(i,j,k)+pdis*3.4d-10*fuvflux/1.21d+7 !KS+2014の式(13)と対応（誤差10パーセント以下）
                  ! if (i==Imax .and. j==Jmax .and. k==Kmax) then
                  !    print '(A, (1P10E15.7))', "KS DEBUG G: ", sqrt(r2)*Unit_l/cgs_au, fuvflux, kh2pd(i,j,k)
                  ! endif
               end do
            end do
         end do
      end do
    end subroutine GetKH2pdThin



    !-------------------------------------------------------------------------
    ! Get Hm photo-detachemet rate [s^-1] assuming optically thin NIR flux
    !-------------------------------------------------------------------------
    subroutine GetKHmpdThin(x,y,z,khmpd)
      real(kind=DBL_KIND),dimension(:),pointer,intent(IN) :: x, y, z
      real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK),intent(OUT) :: khmpd

      integer :: isrc,i,j,k
      real(kind=DBL_KIND) :: r2,r_star

      khmpd(:,:,:)=0.d0
      do isrc=0, rs_info%nsource-1
         r_star = sqrt(rs_info%lum(isrc)/(4.*Pi * cgs_sigma*rs_info%Trad(isrc)**4)) !radius
         do i = Imin, Imax
            do j = Jmin, Jmax
               do k = Kmin, Kmax
                  r2 = (x(i) - rs_info%spos(MX,isrc))**2 + &
                       (y(j) - rs_info%spos(MY,isrc))**2 + &
                       (z(k) - rs_info%spos(MZ,isrc))**2
                  khmpd(i,j,k)=khmpd(i,j,k)+rs_info%hhm(isrc)*r_star**2/(r2*Unit_l**2) 
               end do
            end do
         end do
      end do

    end subroutine GetKHmpdThin


  end subroutine ch_PrimChemistry


  !-------------------------------------------------------------------------
  ! update chemistry and temperature in a cell
  !
  ! INPUT
  !  - heat:  energy deposited at ionization [erg]
  !
  !-------------------------------------------------------------------------
  subroutine ch_CellChemCool(xch)
    type(chem_vals) :: xch 
    real(kind=DBL_KIND),parameter :: eps_imp = 1d-1
    real(kind=DBL_KIND) :: T_o, y_o(0:NCHEM-1), t_chemcool, t_cool, yco_o
    integer :: myrank ! KS DEBUG

    !初期値を保存しておく
    T_o = xch%Tg
    y_o(:) = xch%ychem(:)
    yco_o  = xch%yco


    !ひとまずexplicitに計算してみる
    call CoolSolverExplicit(xch,t_chemcool,t_cool)

    !dt << t_coolだった場合にはimplicit法で計算し直す
    if (xch%dt > eps_imp * t_cool) then
       xch%Tg       = T_o
       xch%ychem(:) = y_o
       xch%yco      = yco_o
       call CoolSolverImplicit(xch) !HFADDED

       num_implicit = num_implicit + 1 ! KS DEBUG

       ! myrank = get_myrank() 
       ! print '(A,I4, 8(1P1E15.7))', "implicit_cooling called: myrank,xnH,T_K,y=", myrank, xnH,T_K,y(:)
       
    end if

    !chemistry/cooling timescale (node local)
    t_chemcool_local = MIN(t_chemcool_local,t_chemcool)

  end subroutine ch_CellChemCool

  !-------------------------------------------------------------------------
  ! Evaluate radiation force by EUV for given fractional abundance
  !
  ! fx[yz]_hpi obtained in art.F90 is force per Hn particle
  !    ->   convert it to force per volume [g s^-2 cm^-2]
  !-------------------------------------------------------------------------
  subroutine ch_EvalRadForce(xnH,yHn,fx_hpi,fy_hpi,fz_hpi)
    real(kind=DBL_KIND),intent(IN) :: xnH,yHn
    real(kind=DBL_KIND),intent(INOUT) :: fx_hpi,fy_hpi,fz_hpi

    fx_hpi = fx_hpi * xnH * yHn
    fy_hpi = fy_hpi * xnH * yHn
    fz_hpi = fz_hpi * xnH * yHn

  end subroutine ch_EvalRadForce
#endif !MODEL_ART 
  !-------------------------------------------------------------------------
  ! reset temperature to MP_T0 for isotermal calculation
  !-------------------------------------------------------------------------
#ifdef ISOTHERMAL
  subroutine ch_reset_temperature
    use modelParameter,only : MP_T0

    integer :: level, n, gid, i, j, k
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp
    real(kind=DBL_KIND) :: xmu
    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )

          ! 子グリッドを持つ場合はchemistryは解かなくてよい (あとでconvergeを呼ぶ)
          !if (has_child_grid(gid))   cycle


          !変数の取り出し
          rho => get_Ucomp(MRHO,gid)
          p => get_Ucomp(MP,gid)
#if MODEL_ART > 0
          yhn => get_Ucomp(MHN,gid)
          yh2 => get_Ucomp(MH2,gid)
          yel => get_Ucomp(MEL,gid)
          yhp => get_Ucomp(MHP,gid)
#endif !MODEL_ART

          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax                   
#if MODEL_ART > 0
                   xmu = (1.d0 + 4.d0*yHe) /(yhn(i,j,k)+yh2(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe) !平均分子量 (Hm, H2pは無視)
#else !MODEL_ART
                   xmu = 1.d0 !MODEL_ART == 0 -> 化学考慮せず
#endif !MODEL_ART
                   !set pressure to the value corresponding to MP_T0
                   p(i,j,k) = (cgs_kb*MP_T0)*(rho(i,j,k)*Unit_rho)/(cgs_mh*xmu) / Unit_e
                end do
             end do
          end do
       end do
    end do

  end subroutine ch_reset_temperature
#endif !ISOTHERMAL


  !-------------------------------------------------------------------------
  ! check all cells have proper values of abundance
  ! temperature floor is currently not imposed here (KS TODO)
  !-------------------------------------------------------------------------
  subroutine ch_check_all_cells
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
#ifdef CHEM_MODEL_HF2020
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: yco
#endif
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
    real(kind=DBL_KIND) :: yco_l
#ifdef METAL_TRANSFER
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
#endif
    integer :: level, n, gid, i, j, k, ic

    ! ----------------------------------------------------------
    type chemsp 
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(NCEHM_MIN:NCEHM_MAX)
    integer :: ichem
    ! ----------------------------------------------------------

    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )

          ! get chemistry pointer
          do ichem = NCEHM_MIN, NCEHM_MAX
            chemary3(ichem)%y => get_Ucomp(ichem,gid)
          enddo
#ifdef CHEM_MODEL_HF2020
          yco  => get_Ucomp(MCO, gid)
#endif


#ifdef METAL_TRANSFER
          mmetal => get_Ucomp(MMET,gid)
#endif

          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax

                   ! --------------------------
                   !   化学組成変数の詰め替え
                   ! --------------------------
                   do ichem = 0, NCHEM-1
                      ychem(ichem) = chemary3(ichem+NCEHM_MIN)%y(i,j,k)
                   enddo 
#ifdef CHEM_MODEL_HF2020
                   yco_l = yco(i,j,k)
#endif


                   !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
                   call adjust_abundance(ychem &
#ifdef CHEM_MODEL_HF2020
                      , yco_l &
#endif
#ifdef METAL_TRANSFER
                      , mmetal(i,j,k))
#else
                      , MP_Metallicity)
#endif

                   ! --------------------------
                   !   化学組成変数の詰め直し
                   ! --------------------------
                   do ichem = 0, NCHEM-1
                      chemary3(ichem+NCEHM_MIN)%y(i,j,k) = ychem(ichem)
                   enddo 
#ifdef CHEM_MODEL_HF2020
                   yco(i,j,k) = yco_l
#endif

                end do
             end do
          end do
       end do
    end do
  end subroutine ch_check_all_cells


  ! ----------------------------------------------------------------------------
  !     update velocity field
  !     all valuables are noD
  ! ---------------------------------------------------------------------------

  subroutine ch_EvalVelocity(vx, vy, vz,fx_hpi,fy_hpi,fz_hpi, dt)
        
        implicit none

        real(kind=DBL_KIND),intent(IN) :: dt
        real(kind=DBL_KIND),intent(IN) :: fx_hpi,fy_hpi,fz_hpi
        real(kind=DBL_KIND),intent(INOUT) :: vx, vy, vz
        real(kind=DBL_KIND),parameter :: V_UPPER_LIMIT = 1.d2 ! for NL=5
        real(kind=DBL_KIND) :: v2

        ! update velocity field

        vx = vx + fx_hpi * dt / Unit_v  ! [noD]
        vy = vy + fy_hpi * dt / Unit_v
        vz = vz + fz_hpi * dt / Unit_v

        v2 = (vx*vx + vy*vy + vz*vz)*Unit_kms*Unit_kms

        if(.not. v2 <= V_UPPER_LIMIT*V_UPPER_LIMIT) then
            vx = vx * V_UPPER_LIMIT / sqrt(v2)
            vy = vy * V_UPPER_LIMIT / sqrt(v2)
            vz = vz * V_UPPER_LIMIT / sqrt(v2)
        endif


  end subroutine ch_EvalVelocity



end module chemistry
