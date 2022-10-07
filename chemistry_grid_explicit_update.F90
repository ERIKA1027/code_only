#include "config.h"

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
#if MODEL_ART > 0
  use radiationSource
  use primordial,only: CoolSolverExplicit, CoolSolverImplicit, adjust_abundance, yhe, c_H2, dbg_flg_prm, prim_GridChemCoolExplicit
#endif !MODEL_ART    

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
  integer :: num_cell ! num_cell: number of cells for chem update
  integer :: num_implicit ! number of cells with implicit temp update
  integer :: num_expl2 ! number of cells with explicit chem/temp update  
  !----------------- KS DEBUG (END) ----------------!

  public :: ch_artchem
#endif !MODEL_ART
#ifdef ISOTHERMAL
  public :: ch_reset_temperature
#endif !ISOTHERMAL
  

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
    real(kind=DBL_KIND):: time_chem, time_chem_min,time_chem_mean,time_chem_max, &
         time_cell,time_cell_min,time_cell_mean,time_cell_max, &
         frac_imp,frac_imp_min,frac_imp_mean,frac_imp_max, &
         frac_expl2,frac_expl2_min,frac_expl2_mean,frac_expl2_max    
    

    call system_clock(time_ini) ! 時間計測開始 (KS DEBUG)

#ifdef NO_RADIATION
    n_substep = 1       !輻射無しの時はsubstep切らない
#else !NO_RADIATION
    n_substep = 2       !輻射ありの時はsubstep2回切る
#endif !NO_RADIATION


    !initial time step
    dt = dt_code*Unit_t    
    t_sub=0d0

    dt_sub = dt/n_substep    

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
       
       !--------------------------- BEGIN CHEM  ------------------------------!       
       if(get_myrank() == PRIMARY_RANK) &
            print *, 'Begin chem update'

       ! ---------- initialize (KS DEBUG) ----------!
       call system_clock(time_prev) ! 時間計測開始 (KS DEBUG)
       num_cell = 0
       num_implicit = 0
       num_expl2 = 0        
       ! --------------- (KS DEBUG) ----------------!       

       call ch_PrimChemistry(dt_sub)

       ! ------------------------------ (KS DEBUG, BEG) ------------------------------- !       
       call system_clock(time_cur,time_rat) ! 時間計測

       time_chem = (time_cur-time_prev)/dble(time_rat) !chemistryにかかった時間
       time_cell = time_chem/dble(num_cell)            !1セル当たりchemistryにかかった時間
       frac_imp = dble(num_implicit)/dble(num_cell)    !implicit coolingが呼ばれた割合
       frac_expl2 = dble(num_expl2)/dble(num_cell)     !explicit cooling/chemが呼ばれた割合       
       
       call mpi_allreduce(time_chem, time_chem_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)       
       call mpi_allreduce(frac_expl2, frac_expl2_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_expl2, frac_expl2_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_expl2, frac_expl2_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)       

       time_chem_mean = time_chem_mean/NPE
       time_cell_mean = time_cell_mean/NPE
       frac_imp_mean = frac_imp_mean/NPE
       frac_expl2_mean = frac_expl2_mean/NPE              
       
       if(get_myrank() == PRIMARY_RANK) then
          print '(A,2(A,1P3E10.2,A),/,2(A,1P3E10.2))', "TIME for CHEM (min,mean,max)  => ", &
               "tot time: ",time_chem_min,time_chem_mean,time_chem_max, " [s]",&
               " time/cell: ",time_cell_min,time_cell_mean,time_cell_max, " [s]",&
               " frac_imp: ", frac_imp_min,frac_imp_mean,frac_imp_max, &
               " frac_expl2: ", frac_expl2_min,frac_expl2_mean,frac_expl2_max
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
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp, yhm, yh2p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: khmpd, kh2pd_nogh, khpi_nogh !ghost cell含まない変数
#if MODEL_ART == 1
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: kh2pd
#elif MODEL_ART == 2
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kh2pd
#endif
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK,0:NCHEM-1) :: ychem, ychem_o
    real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK)  :: T_K, xNc_H2, xnH, heat, rHpi,rH2pd,rHmpd, T_K_o
    real(kind=DBL_KIND) :: xmu, cs, xlmbdj
    logical,dimension(ARRAYSIZE_IJK) :: flg_recalc

    
    integer :: level, n, gid, i, j, k, ic

    logical :: isNotFinite ! defined in naninf.F90
    logical :: naninf_flg
    integer :: count_output=0, max_count=100 !KS DEBUG
    real(kind=DBL_KIND) :: probe_radius !KS DEBUG

    naninf_flg = .False.

    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )

          globdbg_mygid = gid ! for debug

          ! 子グリッドを持つ場合はchemistryは解かなくてよい (あとでconvergeを呼ぶ)
          if (has_child_grid(gid))   cycle

          !座標
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)


          !変数の取り出し
          rho => get_Ucomp(MRHO,gid)
          p => get_Ucomp(MP,gid)
          yhn => get_Ucomp(MHN,gid)
          yh2 => get_Ucomp(MH2,gid)
          yel => get_Ucomp(MEL,gid)
          yhp => get_Ucomp(MHP,gid)
          yhm => get_Ucomp(MHM,gid)
          yh2p => get_Ucomp(MH2P,gid)
          khpi => get_Ucomp(MKPI,gid)
          heat_hpi => get_Ucomp(MHPI,gid)

#ifndef NO_RADFORCE
          fx_hpi => get_Ucomp(MXPI,gid)
          fy_hpi => get_Ucomp(MYPI,gid)
          fz_hpi => get_Ucomp(MZPI,gid)
#endif !NO_RADFORCE

          !get kh2pd
#if MODEL_ART == 1
          call GetKH2pdThin(x,y,z,kh2pd) !w/o self-shielding
#elif MODEL_ART == 2
          kh2pd => get_Ucomp(MKPD,gid) !w/ self-shielding
#endif
          !get khmpd
          call GetKHmpdThin(x,y,z,khmpd)


          !--------------------- pre-calculation before update ------------------!
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax
                   !化学組成変数の詰め替え
                   ychem(i,j,k,0) = yhn(i,j,k)
                   ychem(i,j,k,1) = yh2(i,j,k)
                   ychem(i,j,k,2) = yel(i,j,k)
                   ychem(i,j,k,3) = yhp(i,j,k)
                   ychem(i,j,k,4) = yhm(i,j,k)
                   ychem(i,j,k,5) = yh2p(i,j,k)


                   !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
                   call adjust_abundance(ychem(i,j,k,:))

                   !derived variables
                   xnH(i,j,k) = rho(i,j,k)*Unit_rho/((1.d0 + 4.d0*yHe)*cgs_mh)                !水素原子核の数密度
                   xmu = (1.d0 + 4.d0*yHe) /(ychem(i,j,k,0)+ychem(i,j,k,1)+ychem(i,j,k,2)+ychem(i,j,k,3)+yHe)  !平均分子量 (Hm, H2pは無視)
                   T_K(i,j,k) = p(i,j,k)*Unit_e*cgs_mh*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb      !温度 [K]

                   !impose T floor (before chem_update)
                   if (T_K(i,j,k) < MP_Tmin*0.999d0 .or. T_K(i,j,k) > MP_Tmax*1.001d0) then !表示を抑制するため係数をかける
                      !output warning
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** before chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              T_K(i,j,k), MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         !suppress warning
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   T_K(i,j,k) = min(max(T_K(i,j,k),MP_Tmin),MP_Tmax)

                   !輻射場も変数を詰め替え
                   if (khpi(i,j,k) > 0.d0) then
                      heat(i,j,k) = heat_hpi(i,j,k)/khpi(i,j,k) !heat_hpi/khpi = energy deposit per ionization (KS NOTE)
                   else
                      heat(i,j,k) = 0.d0
                   end if
                   !ゴーストセル含まない配列を定義 (不要? KS TODO)
                   khpi_nogh(i,j,k) = khpi(i,j,k)
                   kh2pd_nogh(i,j,k) = kh2pd(i,j,k)                   

                   !---------------------- H2 LINE TRAPPING -----------------------!
                   cs  = sqrt(cgs_kb*T_K(i,j,k)/(xmu*cgs_mh))          !等温音速 [cm/s]
                   xlmbdj = cs*sqrt(Pi/(cgs_gc*rho(i,j,k)*Unit_rho))   !ジーンズ長さ [cm]
                   xNc_H2 = ychem(i,j,k,1)*xnH(i,j,k)*0.5*xlmbdj             !水素分子の柱密度の概算 [cm^-2]
                   !---------------------------------------------------------------!

                   !-------------        check Nan/Inf     -------------!
                   if (isNotFinite(xnH(i,j,k)) .or. isNotFinite(T_K(i,j,k))) then
                      naninf_flg = .True.
                   end if
                   do ic=0,NCHEM-1
                      if (isNotFinite(ychem(i,j,k,ic))) then
                         naninf_flg = .True.
                      end if
                   end do
                   if (naninf_flg) then
                      print '(A, (1P8E15.7))', "(chemistry) NaN/Inf found before chem update: ", xnH(i,j,k),T_K(i,j,k),ychem(i,j,k,:)
                      print *, "stopping..."
                      stop
                   end if

                   !save original values
                   ychem_o(i,j,k,:) = ychem(i,j,k,:)
                   T_K_o(i,j,k) = T_K(i,j,k)                   

                   !----------------------------------------- KS DEBUG -----------------------------------------!                   
                   ! probe_radius = 1.5e4/Unit_au !KS DEBUG
                   probe_radius = 3e4/Unit_au !KS DEBUG                   
                   if (abs(x(i)-probe_radius) < abs(x(5)-x(4))*0.5 .and. x(i) > 0. .and. &
                       abs(y(j)) < abs(y(5)-y(4)) .and. y(j) > 0. .and. &
                       abs(z(k)) < abs(z(5)-z(4)) .and. z(k) > 0. ) then ! x ~ probe_rad au の値を表示
                      print '(/,A, 1P3E15.7)', "*** KS DEBUG (INI) ***  x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au                      
                      print '(A, 6I6)', "rank,gid,i,j,k,level = ", get_myrank(), gid, i-lbound(rho,1), j-lbound(rho,2), k-lbound(rho,3), level
                      print '(A, (1P10E15.7))', "xnH, T_K, kHpi, kH2pd, kHmpd, dt:",&
                           xnH(i,j,k),T_K(i,j,k),khpi_nogh(i,j,k),kh2pd_nogh(i,j,k),khmpd(i,j,k), dt
                      print '(A, (1P10E15.7))', "ychem: ", ychem(i,j,k,:)
                      print '(A,(1P2E15.7))',    "heat, xNc_H2: ", heat(i,j,k), xNc_H2(i,j,k)
                   end if
                   !----------------------------------------- KS DEBUG -----------------------------------------!                   
                end do
             end do
          end do
          !--------------------- pre-calculation before update (END) ------------------!          

          !----------------------------------------------------------------------------------!
          !                       try explicit chem/temp update                              !
          !----------------------------------------------------------------------------------!

          call prim_GridChemCoolExplicit(xnH,T_K,ychem,kHpi_nogh,kH2pd_nogh,kHmpd,heat,xNc_H2,dt,flg_recalc)

          !----------------------------------------- KS DEBUG -----------------------------------------!                   
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax
                   if (abs(x(i)-probe_radius) < abs(x(5)-x(4))*0.5 .and. x(i) > 0. .and. &
                       abs(y(j)) < abs(y(5)-y(4)) .and. y(j) > 0. .and. &
                       abs(z(k)) < abs(z(5)-z(4)) .and. z(k) > 0. ) then ! x ~ probe_rad au の値を表示          
                      print '(/,A, (1P10E15.7),/)', "*** KS DEBUG (AFTER EXPL UPDATE) ***  x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au                      
                      print '(A, 6I6)', "rank,gid,i,j,k,level = ", get_myrank(), gid, i-lbound(rho,1), j-lbound(rho,2), k-lbound(rho,3), level
                      print '(A, (1P10E15.7))', "xnH, T_K, rHpi, rH2pd, rHmpd, dt:",&
                           xnH(i,j,k),T_K(i,j,k),khpi(i,j,k),kh2pd(i,j,k),khmpd(i,j,k), dt
                      print '(A, (1P10E15.7))', "ychem: ", ychem(i,j,k,:)
                      print '(A,(1P2E15.7))',    "heat, xNc_H2: ", heat(i,j,k), xNc_H2(i,j,k)
                      print '(A,L)', "flg_recalc: " , flg_recalc(i,j,k)
                   end if
                end do
             end do
          end do
          !----------------------------------------- KS DEBUG -----------------------------------------!                             

          !----------------------------------------------------------------------------------!
          !                implicit chem/temp update                                         !
          !                        for cells w/ flag_recalc = .TRUE.                         !
          !----------------------------------------------------------------------------------!

          ! !KS DEBUG
          ! flg_recalc(:,:,:) = .True. !強制的に再計算

          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax

                   !再計算が必要か
                   if (.not. flg_recalc(i,j,k)) then
                      num_expl2 = num_expl2 + 1 ! KS DEBUG
                      cycle ! 再計算が不要なセルはスキップ
                   end if

                   !再計算の前に値をリセット
                   ychem(i,j,k,:) = ychem_o(i,j,k,:)
                   T_K(i,j,k) = T_K_o(i,j,k)

                   !----------------------------------------- KS DEBUG -----------------------------------------!
                   if (abs(x(i)-probe_radius) < abs(x(5)-x(4))*0.5 .and. x(i) > 0. .and. &
                       abs(y(j)) < abs(y(5)-y(4)) .and. y(j) > 0. .and. &
                       abs(z(k)) < abs(z(5)-z(4)) .and. z(k) > 0. ) then ! x ~ probe_rad au の値を表示

                   ! if (khpi(i,j,k) >1d-10) then
                   
                   
                   ! if (abs(xnH- 3.3438341E+07)/xnH < 1d-4) then

                   ! if (abs(T_K-1.7395534E+02)/1.7395534E+02 < 1d-7) then
                   ! if (abs(x(i)+4.1424206d2/Unit_au) < abs(x(5)-x(4))*0.5 .and. &
                   !      abs(y(j)-2.0120329d2/Unit_au) < abs(y(5)-y(4))*0.5  .and. &
                   !      abs(z(k)+1.3019036d2/Unit_au) < abs(z(5)-z(4))*0.5 ) then 


                      print '(/,A, (1P10E15.7),/)', "*** KS DEBUG (BEFORE IMPL UPDATE)***  x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au                      
                      print '(A, 6I6)', "rank,gid,i,j,k,level = ", get_myrank(), gid, i-lbound(rho,1), j-lbound(rho,2), k-lbound(rho,3), level
                      print '(A, (1P10E15.7))', "xnH, T_K, rHpi, rH2pd, rHmpd, dt:",&
                           xnH(i,j,k),T_K(i,j,k),khpi(i,j,k),kh2pd(i,j,k),khmpd(i,j,k), dt
                      print '(A, (1P10E15.7))', "ychem: ", ychem(i,j,k,:)
                      print '(A,(1P2E15.7))',    "heat, xNc_H2: ", heat(i,j,k), xNc_H2(i,j,k)
                      ! print '(A,(1P10E15.7))', 'KS DEBUG 2', p(i,j,k), yhn(i,j,k),yh2(i,j,k),yel(i,j,k),yhp(i,j,k)
                      ! print '(A,(1P10E15.7))', 'KS DEBUG 2A', &
                      !      p(i,j,k)*Unit_e*cgs_mh*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb, p(i,j,k),Unit_e,cgs_mh,xmu,rho(i,j,k),Unit_rho,cgs_kb
                      dbg_flg_prm = 1
                   end if
                   !----------------------------------------- KS DEBUG -----------------------------------------!

                   !-----------    セルの時間発展 (KS NOTE: heat_hpi/khpi = energy deposit per ionization) -------!
                   call ch_CellChemCool(xnH(i,j,k),T_K(i,j,k),ychem(i,j,k,:),&
                        khpi(i,j,k),kh2pd(i,j,k),khmpd(i,j,k),heat(i,j,k),xNc_H2(i,j,k),dt)
                   num_cell = num_cell + 1 ! KS DEBUG
                   !--------------------------------------------------------------------------------------------!

                   !----------------------------------------- KS DEBUG -----------------------------------------!
                   if(dbg_flg_prm == 1) then
                      print '(A, (1P10E15.7))', "(KS DEBUG) x, y, z [au]: ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au
                      print '(A, (1P10E15.7))', "xnH, T_K (after):", xnH(i,j,k),T_K(i,j,k)
                      print '(A, (1P6E15.7),/)', "ychem (after): ", ychem(i,j,k,:)
                      dbg_flg_prm = 0
                      ! stop
                   end if
                   !----------------------------------------- KS DEBUG -----------------------------------------!
                end do
             end do
          end do

          !----------------------- post-calculation after update -----------------------!          
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax
                   !-------------      check Nan/Inf     -------------!
                   if (isNotFinite(xnH(i,j,k)) .or. isNotFinite(T_K(i,j,k))) then
                      naninf_flg = .True.
                   end if
                   do ic=0,NCHEM-1
                      if (isNotFinite(ychem(i,j,k,ic))) then
                         naninf_flg = .True.
                      end if
                   end do
                   if (naninf_flg) then
                      print '(A, (1P6E15.7))', "(chemistry) NaN/Inf found after chem update: ", ychem(i,j,k,:)
                      print '((1P10E15.7))',  xnH(i,j,k),T_K(i,j,k),khpi(i,j,k),heat(i,j,k),x(i),y(j),z(k), dt
                      print *, "stopping..."
                      stop
                   end if

                   !impose T floor (after chem_update)
                   if (T_K(i,j,k) < MP_Tmin .or. T_K(i,j,k) > MP_Tmax) then
                      !output warning
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** after chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              T_K(i,j,k), MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         !suppress warning
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   T_K(i,j,k) = min(max(T_K(i,j,k),MP_Tmin),MP_Tmax)

                   !更新結果をセル変数に詰め直し
                   yhn(i,j,k) = ychem(i,j,k,0)
                   yh2(i,j,k) =  ychem(i,j,k,1)
                   yel(i,j,k) = ychem(i,j,k,2)
                   yhp(i,j,k) = ychem(i,j,k,3)
                   yhm(i,j,k) = ychem(i,j,k,4)
                   yh2p(i,j,k) = ychem(i,j,k,5)
                   xmu = (1.d0 + 4.d0*yHe) /(ychem(i,j,k,0)+ychem(i,j,k,1)+ychem(i,j,k,2)+ychem(i,j,k,3)+yHe)  !平均分子量 (Hm, H2pは無視)
                   p(i,j,k) = (cgs_kb*T_K(i,j,k))*(rho(i,j,k)*Unit_rho)/(cgs_mh*xmu) / Unit_e      !温度 [K]
                end do
             end do
          end do

          !--------------------- post-calculation before update (END) ------------------!

#ifndef NO_RADFORCE
          !------------------------------- evaluate radiation force ----------------------------------!
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax                   
                   call ch_EvalRadForce(xnH,yhn(i,j,k),fx_hpi(i,j,k),fy_hpi(i,j,k),fz_hpi(i,j,k))
                end do
             end do
          end do
#endif !NO_RADFORCE
                   
          
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
                  khmpd(i,j,k)=khmpd(i,j,k)+2.5d-11*rs_info%hhm(isrc)*r_star**2/(r2*Unit_l**2) !式の意味は細川さんに聞く (KS TODO)
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
  subroutine ch_CellChemCool(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt)
    real(kind=DBL_KIND),intent(IN) :: xnH,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt
    real(kind=DBL_KIND),intent(INOUT) :: T_K,y(0:NCHEM-1)

    real(kind=DBL_KIND),parameter :: eps_imp = 1d-1

    real(kind=DBL_KIND) :: T_o, y_o(0:NCHEM-1), t_chemcool, t_cool

    integer :: myrank ! KS DEBUG

    !初期値を保存しておく
    T_o = T_K
    y_o(:) = y(:)


    !ひとまずexplicitに計算してみる
    call CoolSolverExplicit(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt,t_chemcool,t_cool)

    !dt << t_coolだった場合にはimplicit法で計算し直す
    if (dt > eps_imp * t_cool) then
       T_K = T_o
       y(:) = y_o
       call CoolSolverImplicit(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt)

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
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp, yhm, yh2p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem
    integer :: level, n, gid, i, j, k, ic

    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )

          !変数の取り出し
          yhn => get_Ucomp(MHN,gid)
          yh2 => get_Ucomp(MH2,gid)
          yel => get_Ucomp(MEL,gid)
          yhp => get_Ucomp(MHP,gid)
          yhm => get_Ucomp(MHM,gid)
          yh2p => get_Ucomp(MH2P,gid)

          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax

                   !化学組成変数の詰め替え
                   ychem(0) = yhn(i,j,k)
                   ychem(1) = yh2(i,j,k)
                   ychem(2) = yel(i,j,k)
                   ychem(3) = yhp(i,j,k)
                   ychem(4) = yhm(i,j,k)
                   ychem(5) = yh2p(i,j,k)

                   !流体アップデート時に保存則が破れる可能性が（内挿の際に）あるため、abundanceをadjustしておく
                   call adjust_abundance(ychem)

                   !更新結果をセル変数に詰め直し
                   yhn(i,j,k) = ychem(0)
                   yh2(i,j,k) =  ychem(1)
                   yel(i,j,k) = ychem(2)
                   yhp(i,j,k) = ychem(3)
                   yhm(i,j,k) = ychem(4)
                   yh2p(i,j,k) = ychem(5)
                end do
             end do
          end do
       end do
    end do
  end subroutine ch_check_all_cells

end module chemistry
