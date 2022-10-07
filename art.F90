#include "config.h"

!HFADDED -------------------
!number of frequency bins

#if MODEL_ART == 1
  !EUV RT + thin FUV
#ifdef METAL 
#define NUM_FREQ     2
#else
#define NUM_FREQ 1
#endif

#elif MODEL_ART == 2
 !EUV RT + FUV RT
#ifdef METAL 
#define NUM_FREQ 5
#else
#define NUM_FREQ 2
#endif
#endif

! number of dust bin
#define BIN_EUV           0 
#define BIN_H2           1
#ifdef METAL   
#define BIN_DEUV         2
#define BIN_DFUV         3
#define BIN_CO           4
#endif

!HFADDED -------------------


!control debug messages
!#define ART_DEBUG


!-----------------------------------------------------------------------
! subroutine for Adaptive Ray Tracing (ART)
!-----------------------------------------------------------------------
module art
  use overBlockCoordinates
  use grid
  use unit
  use primordial
  use pix_tools
  use mpilib
  use radiationSource
  use sinkParticle, only : sp_getSinkRadius
#ifdef METAL ! HFADDED
  use kinzoku
#endif

  implicit none
  private

  real(kind=DBL_KIND),parameter :: PI=atan(1d0)*4d0

  integer,parameter :: MTH = 0, MPH = 1 ! direction (theta, phi)

  ! 面の方向コード
  integer,parameter :: NIL = 0                  
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5

  ! direction specified with HEALPix 
  type t_dirHpix
     integer :: res      ! HEALPix level (N_side = 2**res)
     integer(kind=LLONG_KIND) :: ipix         ! ipix = (0, N_pix-1); there are N_pix = 12*N_side**2 patches in total 
  end type t_dirHpix

  ! information for each ray
  type t_ray
     integer :: isrc                               ! index of radiaiton source in rs_info
     type(t_dirHpix) :: Hpix                       ! ray direction specified with HEALPix <-- (level, id)
     real(kind=DBL_KIND) :: dist                   ! distance from the origin     
     real(kind=DBL_KIND),dimension(0:NUM_FREQ-1) :: tau    ! optical depth from the origin (0 -> EUV, 1 -> FUV (NcH2), 2 -> DUST)
     integer :: gid                                ! gid for the grid the ray is entering
     type(t_ray),pointer :: next => null()         ! used to make linklist
 end type t_ray

  ! pointer for t_ray (needed to create array for pointers)
  type t_pray
     type(t_ray),pointer :: ray_list => null()
  end type t_pray

  ! link list for rays
  type(t_ray),save,pointer :: my_ray_list => null()         ! list of rays belong to the current rank (初期化もしておく、番兵ノード)
  type(t_pray),save,dimension(0:NPE-1) :: send_ray_list     ! list of rays to be sent to other ranks
  integer,save,dimension(0:NPE-1) :: N_send                 ! 各rankに送信するrayの本数

  !info for radiation source
  type(t_rs_info),save,pointer :: rs_info

  !data size/type for MPI communication
  integer,parameter :: RAY_BYTE = (4*3+8*(2+NUM_FREQ)) !isrc,res,ipix,dist,tau(0:NUM_FREQ-1),gid
  integer,parameter :: NSEND_BYTE = 4 ! integer  N_send(rank)
  integer,parameter :: NCOMP_BYTE = 8 ! long integer  N_comp_local
  integer,parameter :: N_SEND_MAX = 1024 ! maximum number of rays for each message
  integer,parameter :: BUFSIZE_RAYS = (N_SEND_MAX*RAY_BYTE+NSEND_BYTE+3)/4 !rayの送受信に必要なバッファサイズ (SendRaysのinteger換算のサイズに対応、切り上げ注意)
  integer,parameter :: BUFSIZE_NCOMP = (NCOMP_BYTE+3)/4   !rayの送受信に必要なバッファサイズ (SendRaysのinteger換算のサイズに対応、切り上げ注意)  
  integer,save,dimension(0:BUFSIZE_RAYS-1) :: bufr_rays   !rayの受信用buffer
  integer(kind=LLONG_KIND),save :: bufr_Ncomp             !Ncompの受信用buffer
  integer,parameter :: MPI_TAG_MAX_MACHINE = 2097151  !マシン毎に決まっているMPIで使えるtagの最大値 (cray-mpich/7.4.4 だと最大が2097151らしい)

  ! information for sent message of rays
  type t_msg_rays
     integer :: req                                     ! request number of MPI message
     integer :: tag                                     ! tag of MPI message
     integer,pointer,dimension(:) :: bufs               ! send buffer (message毎に作成、動的に確保)
     type(t_msg_rays),pointer :: next => null()         ! used to make linklist
  end type t_msg_rays

  ! information for sent message of Ncomp
  type t_msg_Ncomp
     integer :: req                                     ! request number of MPI message
     integer :: tag                                     ! tag of MPI message
     integer(kind=LLONG_KIND) :: bufs                   ! send buffer (message毎に作成、サイズはlong integerの8バイト)
     type(t_msg_Ncomp),pointer :: next => null()        ! used to make linklist
  end type t_msg_Ncomp

  ! link list for mpi message info
  type(t_msg_rays),save,pointer :: msg_rays_list => null()         ! list of sent messages for rays (初期化もしておく、番兵ノード)
  type(t_msg_Ncomp),save,pointer :: msg_Ncomp_list => null()         ! list of sent messages for Ncomp (初期化もしておく、番兵ノード)

  ! work completion index
  !  integer,parameter :: HEALPIX_MIN_LEVEL = 2
  integer,parameter :: HEALPIX_MIN_LEVEL = 3      !subgrid disk方向の解像度が低いとartificialな吸収が生じやすくなりそうなので、ひとまず少し解像度を上げてみる (KS TODO)
  integer,parameter :: HEALPIX_MAX_LEVEL = 12     !res = HEALPIX_MAX_LEVELまで許す？  
  integer,parameter :: MIN_RAY_PER_CELL = 3
  !  integer,parameter :: MIN_RAY_PER_CELL = 10   !KS DEBUG
  
  !  integer,parameter :: N_source = 1 ! temporal (KS DEBUG)
  !  integer(kind=LLONG_KIND),parameter :: N_comp_total = (2_8**(HEALPIX_MAX_LEVEL-HEALPIX_MIN_LEVEL))**2 ! temporal (KS DEBUG) 
  !  integer(kind=LLONG_KIND),parameter :: N_comp_total = N_source * 12 * (2_8**(HEALPIX_MAX_LEVEL))**2
  integer(kind=LLONG_KIND) :: N_comp_total ! depends on number of sources
  integer(kind=LLONG_KIND),save,dimension(0:NPE-1) :: N_comp_local ! N_comp_global = \sum_{rank} N_comp_local(rank)
  integer(kind=LLONG_KIND) ::  N_comp_global ! work completed when N_comp_global == N_source * 12 * 2**HEALPIX_MAX_LEVEL
  logical :: N_comp_local_update !N_comp_localにupdateがあったかどうか
  real(kind=DBL_KIND),parameter :: TAU_ENDRAY = 100.d0 !tau がこれ以上大きくなったらrayを終了
  !  real(kind=DBL_KIND),parameter :: TAU_ENDRAY = 1.d4 !KS DEBUG

  ! matrix of random rotation
  real(kind=DBL_KIND),dimension(MX:MZ,MX:MZ) :: m_rot

  ! raduis of sinkparticles
  real(kind=DBL_KIND) :: SinkRadius

  ! parameters controlling the timing of MPI communications
  integer,parameter :: N_SEND_GO = 200   ! MPI_issend is called if maxval(N_send) >= N_SEND_GO
  integer,parameter :: N_GRIDCOUNT_GO = 2000   ! MPI_issend is called if grid_count >= N_GRIDCOUNT_GO 

  ! 今までに使用したtagの最大値 (tagの値は再利用しないことにする)
  integer,save :: tag_max_ray, tag_max_Ncomp

  ! variables for debug
  integer(kind=LLONG_KIND) :: loop_count, grid_count, grid_count_tot !outer loop/cross grid (total)
  integer::t_ini,t_cur,t_rat ! t_ini: initial time, t_cur: current time, t_rat: time rate
  logical :: glob_debug_flag_art = .FALSE. ! global variables for debug message
  integer :: N_myrays ! number of my_rays

  !integer,parameter :: Left = 0, Right = 1
#ifdef RADTR_DIRECT
  public :: art_RayTracing
contains
  !-------------------------------------------------------------------------
  ! art_RayTracing
  ! Tracing rays until all of them are vanished
  !-------------------------------------------------------------------------
  subroutine art_RayTracing
    logical :: flg_alldone 
    type(t_ray),pointer :: tmp_ray
    character(100) :: s

    integer :: level, n, gid,rank,i,j,k
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    integer(kind=LLONG_KIND) :: nsource_glob
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi !EUV RT
    
#if MODEL_ART == 2 
  real(kind=DBL_KIND),dimension(:,:,:),pointer ::kh2pd !FUV RT
#endif

#ifdef METAL
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kgfuv, kdph, kdrco, krOII ! FUV lumi and dust heating rate
#endif

    !---- KS DEBUG --  KS DEBUG --  KS DEBUG --  KS DEBUG ---- !
    type(t_msg_Ncomp),pointer :: msg_Ncomp,tmp_msg_Ncomp
    logical :: flg
    integer :: req
    integer,parameter :: maximum_time_debug = 600  ! time when stop run and output debug message in sec
    !---- KS DEBUG --  KS DEBUG --  KS DEBUG --  KS DEBUG ---- !

    !-----------------------  前回のゴミが残っていないかをチェック --------------------!
    call art_TestRecvd() ! 前回のmsgリストがちゃんと受け取り済みであることを確認
    if (associated(my_ray_list) .or. associated(msg_rays_list) .or. associated(msg_Ncomp_list)) then

       !--------- KS DEBUG ------------#
       do i = 0, 100 !100回testしてみて、それでもまだmsg_Ncomp_listが空にならなければエラーを吐いて止める
          if (.not. associated(msg_Ncomp_list)) exit
          msg_Ncomp => msg_Ncomp_list            ! プリント用のポインタ
          do while (associated(msg_Ncomp)) 
             !print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)',&
             !print *,&
             !     "(KS DEBUG) msg_Ncomp w/ myrank=",myrank,", (req,tag,bufs) = (",&
             !     msg_Ncomp%req,",",msg_Ncomp%tag,",",msg_Ncomp%bufs,&
             !     ") not received at the beginning of art_RayTracing"
             msg_Ncomp => msg_Ncomp%next         ! 次のポインタ
          end do
          call art_TestRecvd() ! もう一度TestRecvdしてみる
          !print '(/,A,A,/)',"** WARNING ** call art_TestRecvd() again", &
          !     " because msg_Ncomp_list remains at the beging of art_RayTracing"
       end do
       !--------- KS DEBUG ------------#       

       if (associated(my_ray_list) .or. associated(msg_rays_list) .or. associated(msg_Ncomp_list)) then
          !print '(A,/,A,/,3L3,/,A)', "*** ERROR *** the last call of ART did not end properly.", &
          !     "Either my_ray_list/msg_rays_list/msg_Ncomp_list exists at the begining of ART...", &
          !     associated(my_ray_list),  associated(msg_rays_list), associated(msg_Ncomp_list), &
          !     "stopping..."
          print *, "*** ERROR *** the last call of ART did not end properly.", &
               "Either my_ray_list/msg_rays_list/msg_Ncomp_list exists at the begining of ART...", &
               associated(my_ray_list),  associated(msg_rays_list), associated(msg_Ncomp_list), &
               "stopping..."
          stop
       end if
    end if
    !------------------------------------------------------------------------------!


    !初期化（Fortranだと宣言時の代入は一度しか効かないので注意）
    N_comp_local(:)=0; N_comp_global=0; loop_count=0; grid_count=0; grid_count_tot=0
    N_send=0; tag_max_ray=-1; tag_max_Ncomp=-2
    N_comp_local_update=.FALSE.; flg_alldone=.FALSE.
    N_myrays = 0

    SinkRadius = sp_getSinkRadius() !sink半径 = sinkの臨界密度に対応するジーンズ長の半分

    call mpi_barrier(MPI_COMM_WORLD, ierr) !RayTracingする前に、全部のセルの更新が終わっていることを要請

#ifdef ART_DEBUG
    call system_clock(t_ini) ! 時間計測開始
#endif
    !KS DEBUG
    call system_clock(t_ini) ! 時間計測開始

#ifndef ART_NO_RANDOM_ROTAION
    call art_random_rot_matrix()
#endif ! ART_NO_RANDOM_ROTAION        

    !輻射の影響を保持するセル変数の初期化
    ! do level = Lmin, Lmax 
    do level = Lmin, LevelMax !LevelMax現存するグリッドの最高レベル
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )
          khpi => get_Ucomp(MKPI,gid); khpi(:,:,:) = 0.d0
          heat_hpi => get_Ucomp(MHPI,gid); heat_hpi(:,:,:) = 0.d0
#ifdef EXTERNALFORCE
          fx_hpi => get_Ucomp(MXPI,gid); fx_hpi(:,:,:) = 0.d0
          fy_hpi => get_Ucomp(MYPI,gid); fy_hpi(:,:,:) = 0.d0
          fz_hpi => get_Ucomp(MZPI,gid); fz_hpi(:,:,:) = 0.d0
#endif

#if MODEL_ART == 2

          kh2pd => get_Ucomp(MKPD,gid); kh2pd(:,:,:) = 0.d0 
#endif !MODEL_ART
#ifdef METAL
          kgfuv => get_Ucomp(MGFUV,gid); kgfuv(:,:,:) = 0.d0
          kdph  => get_Ucomp(MDPH,gid);  kdph(:,:,:)  = 0.d0
          kdrco => get_Ucomp(MDPCO,gid); kdrco(:,:,:) = 0.d0
          krOII => get_Ucomp(MKPOII,gid);krOII(:,:,:) = 0.d0
#endif
       end do
    end do

    !overBlockCoordinates関数の初期化
    call ob_init()                         !全rankから呼んで初期化する必要あり
    call mpi_barrier(MPI_COMM_WORLD, ierr) !全員呼ぶまで待たないと、次のob_***の呼び出しでこける    

    !光源の情報を取得
    call rs_GetSourceInfo(rs_info)

    !全仕事量の算出
    nsource_glob = rs_info%nsource
    !N_comp_total = nsource_glob  * 12 * (2_8**(HEALPIX_MAX_LEVEL))**2
    N_comp_total = nsource_glob  * 12_8 * 4_8**HEALPIX_MAX_LEVEL


    !光源リストに基づくrayの割り当て
    call art_RadSources()

    do while (.not. flg_alldone)
       !自分のrankが担当の仕事をとりあえず終わらせる
       do while (associated(my_ray_list))
          !--------------------- KS DEBUG ---------------------!
          !if (mod(grid_count_tot,100)==0) then  ! 呼び出しすぎて時間計測に時間がかかるのを防止 (KS TODO)
          !   call system_clock(t_cur,t_rat) ! 時間計測
          !   if ((t_cur-t_ini)/dble(t_rat) > maximum_time_debug) then
          !      !print '(I4,A,I4,2I18,3I6)', " sec passed, glob_debug_flag_art is set to .TRUE.", maximum_time_debug, &
          !      !     maximum_time_debug, myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
          !      print *, " sec passed, glob_debug_flag_art is set to .TRUE.", maximum_time_debug, &
          !           maximum_time_debug, myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
          !      glob_debug_flag_art = .TRUE.             
          !   end if
          !end if
          !--------------------- KS DEBUG ---------------------!

          tmp_ray => my_ray_list          !rayの取り出し
          my_ray_list => my_ray_list%next !ポインタを進める
          
          call art_CrossGrid(tmp_ray) !グリッド内を前進させる

          grid_count = grid_count + 1                ! count crossgrid
          grid_count_tot = grid_count_tot + 1        ! count crossgrid

          if (maxval(N_send) >=  N_SEND_GO) exit  !ある程度送るrayが溜まったら、いったんループを抜けて送信
          if (grid_count >=  N_GRIDCOUNT_GO) exit  !ある程度 CrossGridを実行したら、いったんループを抜けて送信          

          !--------------------- KS DEBUG ---------------------!
          if (mod(grid_count_tot,100)==0) then ! 呼び出しすぎて時間計測に時間がかかるのを防止 (KS TODO)
             ! 時間計測
             call system_clock(t_cur,t_rat)
             if ((t_cur-t_ini)/dble(t_rat) > maximum_time_debug+5) then
                !print '(I4,A,I4,5I12)', " sec passed, exiting CrossGrid loop", maximum_time_debug+5, &
                !     myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
                !print *, " sec passed, exiting CrossGrid loop", maximum_time_debug+5, &
                !     myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
                exit
             end if

             ! tagの値も大きくなりすぎないようにチェック
             if (tag_max_ray > MPI_TAG_MAX_MACHINE*0.5 .or. tag_max_Ncomp > MPI_TAG_MAX_MACHINE*0.5) then ! 0.5をかけているのは安全係数
                !print '(A,I4,5I12)', "tag value is unexpectedly large, exiting CrossGrid loop", myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
                !print *, "tag value is unexpectedly large, exiting CrossGrid loop", myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
                !write (s,'(A,I0,A,I0,A)') "(A,I4,A,",NPE,"(I12),A,I12,A,I12)"
                !print *, "(tag value check) myrank=", myrank,&
                !     ", N_comp_local=",N_comp_local,", N_comp_global=",N_comp_global,", N_comp_total=",N_comp_total
                exit
             end if
          end if
          !--------------------- KS DEBUG ---------------------!          
       end do

       !reset grid_count
       grid_count = 0

       !仕事の受け渡しと現況報告
       call art_RecvAll()   ! <- RecvAllを先に呼んでおくことで、rayの受信があった際にNcompの送信を思いとどまらせることが可能？
       call art_SendAll()
       call art_CheckAllDone(flg_alldone) !仕事が全部終わったか？

       loop_count = loop_count + 1 ! count loop

       !--------------------- KS DEBUG ---------------------!

       if (mod(loop_count,100)==0) then  ! 呼び出しすぎて時間計測に時間がかかるのを防止 (KS TODO)
          !時間がかかりすぎていたらループを抜けてDebugメッセージを出力する
          call system_clock(t_cur,t_rat) ! 時間計測
          if ((t_cur-t_ini)/dble(t_rat) > maximum_time_debug+10) then
             !print '(I4,A,I4,5I12)', " sec passed, exiting whole loop", maximum_time_debug+10, &
             !myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
             !write (s,'(A,I0,A,I0,A)') "(A,I4,A,",NPE,"(I12),A,I12,A,I12)"
             !print s, "(timeout check) myrank=", myrank,&
             !     ", N_comp_local=",N_comp_local,", N_comp_global=",N_comp_global,", N_comp_total=",N_comp_total
             exit
          end if

          ! tagの値が大きくなったときもループを抜けてDebugメッセージを出力する
          if (tag_max_ray > MPI_TAG_MAX_MACHINE*0.5 .or. tag_max_Ncomp > MPI_TAG_MAX_MACHINE*0.5) then ! 0.5をかけているのは安全係数
             !print '(A,I4,5I12)', "tag value is unexpectedly large, exiting whole loop", myrank, grid_count_tot, loop_count, N_myrays, tag_max_ray, tag_max_Ncomp
             !write (s,'(A,I0,A,I0,A)') "(A,I4,A,",NPE,"(I12),A,I12,A,I12)"
             !print s, "(tag value check) myrank=", myrank,&
             !     ", N_comp_local=",N_comp_local,", N_comp_global=",N_comp_global,", N_comp_total=",N_comp_total
             exit
          end if
       end if
       !--------------------- KS DEBUG ---------------------!
    end do

#ifdef ART_DEBUG
    call system_clock(t_cur,t_rat) ! 時間計測完了
#endif    

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    ! !msg_Ncompはループを出た後に残っていても間違いでは無いが、次のARTに向けてdeallocateしておく
    ! do while (associated(msg_Ncomp_list))       ! msg_Ncomp_listが残っていたら中身を一つずつdeallocate
    !    msg_Ncomp => msg_Ncomp_list              ! deallocate用のポインタ
    !    msg_Ncomp_list => msg_Ncomp_list%next    ! msg_Ncomp_listの指す先は次のmsgに移動
    !    deallocate(msg_Ncomp) 
    ! end do
       
    !------------------------- CHECK PRPOER ENDING -----------------------!
    ! if (associated(my_ray_list) .or. associated(msg_rays_list) .or. associated(msg_Ncomp_list)) then ! rayが残っていないか、メッセージが残っていないかチェック
    !    print '(A,/,A,/,3L3,/,A)', "*** ERROR *** ART does not end properly.", &
    !         "Either my_ray_list or msg_rays_list or msg_Ncomp_list remains...", &
    !         associated(my_ray_list),  associated(msg_rays_list), associated(msg_Ncomp_list), &
    !         "stopping..."
    !    stop
    ! end if
    if (associated(my_ray_list) ) then ! rayが残っていないかチェック
       print '(A,/,A)', "*** ERROR ***  my_ray_list remains after ART", &
            "stopping..."
       stop
    end if    
    !---------------------------------------------------------------------!

#ifdef ART_DEBUG
    print '(/,A,I4,2(A,I0),A,(1P1E10.2),A,/)', "END OF KS DEBUG: myrank=", myrank,&
         ", loop_count=",loop_count, ", grid_count_tot=",grid_count_tot, &
         ", elapsed_time=",  (t_cur-t_ini)/dble(t_rat),"[s]"
    call flush(6) ! 標準出力 (unit=6) についてflush
#endif        

  end subroutine art_RayTracing

  !-------------------------------------------------------------------------
  ! art_RadSources
  ! injecting rays according to the information of radiation sources
  !-------------------------------------------------------------------------
  subroutine art_RadSources()
    !本当は光源リストについてループを回す
    real(kind=DBL_KIND),dimension(MX:MZ) :: r0 
    type(t_dirHpix) :: Hpix         
    real(kind=DBL_KIND) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1) :: tau
    integer(kind=LLONG_KIND) :: ipix, N_pix, idx
    integer :: i,j,k,rank,gid,level
    integer :: nsource_glob, isrc
    real(kind=DBL_KIND),dimension(:,:),pointer :: spos
    real(kind=DBL_KIND),dimension(:),pointer:: slum

    integer :: l, N_base, N_ref, offset_min, offset_max, offset_skip, i_base, i_offset 

!    call rs_GetSourceData(nsource, spos, slum)

    ! !ソース情報を適当に仮定
    ! ! r0(MX) = -20.d0
    ! ! r0(MY) = -15.d0
    ! ! r0(MZ) = -10.d0
    ! r0(MX) = 0.d0
    ! r0(MY) = 2.d0
    ! r0(MZ) = 3.d0    
    ! tau = 0.
    ! d = 0.
    ! l0 = 1.
    ! nsource = 1
    ! if (myrank /= 1) return
    ! Hpix%res = 3
    ! Hpix%ipix = 451
    ! call art_StartRay(r0, Hpix, d, tau, l0)
    ! return

    !source直近でのray本数 (Hpix res = HEALPIX_MIN_LEVELに分割してスタート)
    ! N_pix = 12 * (2_8**HEALPIX_MIN_LEVEL)**2
    N_pix = 12 * 4_8**HEALPIX_MIN_LEVEL

    !計算領域全体の光源の数
    nsource_glob = rs_info%nsource

    !rayを追加
    do isrc = 0, nsource_glob -1

       r0 = rs_info%spos(:,isrc) !光源の位置

#ifdef SUB_GRID_MODEL_DIRECTLIGHT
       rs_info%mskr(isrc) = 0.d0 ! 初期化
#endif

       ! 光源の所属する最高levelのグリッドを求める
       do level = LevelMax, Lmin, -1 ! 一番上のlevelから探査
          call ob_getBlockFromCoordPhys(r0, level, gid, rank, i,j,k)    
          if (gid /= Undefi) exit    ! gidが見つかったら終わり
          if (gid == Undefi .and. level == Lmin) then    ! gidがLminでも見つからなかったらエラーを吐いて止める
             print *, "(RadSource) gid not found, stopping..."
             stop
          end if
       end do
!       call ob_getBlockFromCoordPhys(r0, 0, gid, rank, i,j,k) ! level=0で光源が所属するrankを取得
       if (rank /= myrank) cycle !自分のrankでなければスキップ

       
       !いろいろな方向に均等にrayが積まれるように工夫
       do idx = 0, N_pix-1
          call ipix_wide_dist(idx, HEALPIX_MIN_LEVEL, ipix)
          Hpix%res = HEALPIX_MIN_LEVEL
          Hpix%ipix = ipix
          dist = 0.d0               ! sink内の通過はsinkセルをoptically thinとすることでひとまず対応 (KS TODO)
          tau(:) = 0.d0
          call art_StartRay(isrc, Hpix, dist, tau)
       end do
          
       ! do ipix = 0, N_pix-1
       ! ! do ipix = (N_pix/12)*0, (N_pix/12)*0 ! KS DEBUG
       !    Hpix%res = HEALPIX_MIN_LEVEL
       !    Hpix%ipix = ipix
       !    dist = 0.d0               ! sink内の通過はsinkセルをoptically thinとすることでひとまず対応 (KS TODO)
       !    tau(:) = 0.d0
       !    call art_StartRay(isrc, Hpix, dist, tau)
       ! end do
    end do

  contains    

    ! Healpixでのpixel idは (12,4,4,...,4,4)進数での値に対応
    !
    ! 4*4*,,,,*4*12 のi番目
    !    -> (4,4,...,4,4,12)進数での (a,b,...,c,d,x)
    !    -> (12,4,4,...,4,4)進数での (x,d,c,...,b,a) の値
    ! という関数を定義
    !
    subroutine ipix_wide_dist(idx, res, ipix)
      integer(kind=LLONG_KIND), intent(IN)  :: idx
      integer, intent(IN) :: res
      integer(kind=LLONG_KIND), intent(OUT)  :: ipix

      integer :: n_dig,n_dig_inv
      integer(kind=LLONG_KIND) :: idx_rest      
      integer, dimension(0:res) :: i_each_digit    !(4,4,...,4,4,12)進数での各桁の値
      integer, dimension(0:res) :: unit_dig !(4,4,...,4,4,12)進数での各桁の単位  (10進数の3桁目なら100)
      integer, dimension(0:res) :: unit_dig_inv !(12,4,4,...,4,4)進数での各桁の単位

   
      !(4,4,...,4,4,12)進数、(12,4,4,...,4,4)進数の定義
      unit_dig(0) = 1
      unit_dig_inv(0) = 1      
      do n_dig = 1,res                           !桁についてdoループ
         unit_dig(n_dig) = 12_8 * 4_8**(n_dig-1)           !(4,4,...,4,4,12)進数でのn_dig桁の単位
         unit_dig_inv(n_dig) = 4_8**n_dig           !(4,4,...,4,4,12)進数でのn_dig桁の単位
      end do

      !idx -> (4,4,...,4,4,12)進数での (a,b,...,c,d,x)
      idx_rest = idx
      do n_dig = res,0,-1                           !桁についてdoループ (一番上の桁から)
         i_each_digit(n_dig) = int(idx_rest/unit_dig(n_dig))   !現在の桁の値
         idx_rest = idx_rest &
              - i_each_digit(n_dig) * unit_dig(n_dig)     !残り
      end do

      !(12,4,4,...,4,4)進数での (x,d,c,...,a,b) の値 = ipix
      ipix = 0
      do n_dig = 0,res
         n_dig_inv = res-n_dig
         ipix = ipix &
              + i_each_digit(n_dig_inv) * unit_dig_inv(n_dig)
      end do

      ! print *, "unit_dig", unit_dig(:)
      ! print *, "unit_dig_inv", unit_dig_inv(:)
      ! print *, "i_each_digit", i_each_digit(:)
      ! print *, "idx, res, ipix", idx, res, ipix
      
    end subroutine ipix_wide_dist

  end subroutine art_RadSources

  !-------------------------------------------------------------------------
  ! art_StartRay
  ! add a ray to ray_list (gid unkonw)
  ! the ray_list can be either my_ray_list or send_ray_list 
  ! depending on the rank containing the ray
  ! INPUT: isrc, Hpix, d, tau(:)
  !-------------------------------------------------------------------------
  subroutine art_StartRay(isrc, Hpix, dist, tau)
    integer,intent(IN) :: isrc                                     ! ray origin
    type(t_dirHpix),intent(IN) :: Hpix                             ! ray direction
    real(kind=DBL_KIND),intent(IN) :: dist                         ! distance from the origin
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1),intent(IN) :: tau  ! distance from the origin


    real(kind=DBL_KIND) :: theta, phi
    real(kind=DBL_KIND),dimension(MX:MZ) :: ndir !direction vector
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    integer :: gid, level, rank,i,j,k

    logical :: flg

    type(t_ray),pointer :: ray

    !HEALPix <-> 角度 の変換
    call pix2ang_nest(2**(Hpix%res), Hpix%ipix, theta, phi)
    call art_ang2vec(theta,phi,ndir)    
    pos(:) = rs_info%spos(:,isrc) + dist * ndir

    ! rayの所属する最高levelのグリッドを求める
    do level = LevelMax, Lmin, -1 ! 一番上のlevelから探査
       call ob_getBlockFromCoordPhys(pos, level, gid, rank, i,j,k)    
       if (gid /= Undefi) exit    ! gidが見つかったら終わり
       if (gid == Undefi .and. level == Lmin) then    ! gidがLminでも見つからなかったらエラーを吐いて止める
          print *, "(StartRay) gid not found", pos
          print *, "stopping..."
          stop
       end if
    end do

    !rayの終了条件の判定
    call art_EndCond(isrc,Hpix,dist,tau,flg)
    if (flg) then                            !終了条件を満たす場合は
#ifdef ART_DEBUG
       print '(/,A,I4,A,/)',"***    ray-ending condition is satisfied.  (myrank=",myrank," in art_StartRay)  ***"
#endif
       call art_RecordComp(Hpix%res)         !完了したことを記録して
       return                                !そのまま終了
    endif

    !rayの作成
    call art_CreateRay(isrc, Hpix, dist, tau, gid, ray)


    if (rank == myrank) then
       ! call art_AddMyRay(isrc, Hpix, dist, tau, gid)
       call art_AddMyRay(ray)       
    else 
#ifdef ART_DEBUG
       print '(/,4(A,I4),/)', "(StartRay) myrank=", myrank,", rank=", rank, ", gid=", gid, ", level=", level
#endif
       !check rank/gid
       if (rank < 0 .or. rank > NPE-1 .or. gid == Undefi) then
          print *, "(StartRay) rank/gid is starnge before callling art_AddSendRay: ", rank, gid
          print *, "myrank, pos, rs_info%spos(:,isrc), dist", myrank, pos, rs_info%spos(:,isrc), dist
          print *, "stopping..."
          stop
       end if
       ! call art_AddSendRay(isrc, Hpix, dist, tau, gid, rank)
       call art_AddSendRay(ray, rank)       
    end if

  end subroutine art_StartRay

  !-------------------------------------------------------------------------
  ! art_CreateRay
  ! create ray from ray info
  ! INPUT: isrc, Hpix, dist, tau, gid
  !-------------------------------------------------------------------------
  subroutine art_CreateRay(isrc, Hpix, dist, tau, gid, ray)
    integer,intent(IN) :: isrc                                     ! ray origin
    type(t_dirHpix),intent(IN) :: Hpix                             ! ray direction
    real(kind=DBL_KIND),intent(IN) :: dist                         ! distance from the origin
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1),intent(IN) :: tau  ! distance from the origin
    integer,intent(IN) :: gid                                      ! gid for the grid the ray is entering
    type(t_ray),pointer,intent(OUT) :: ray

    !rayの作成
    allocate(ray)
    ray%isrc = isrc
    ray%Hpix = Hpix
    ray%dist = dist
    ray%tau = tau
    ray%gid = gid
  end subroutine art_CreateRay

  !-------------------------------------------------------------------------
  ! art_AddMyRay
  ! add ray to my_ray_list
  ! INPUT/OUTPUT: ray
  !-------------------------------------------------------------------------
  subroutine art_AddMyRay(ray)
    type(t_ray),pointer,intent(INOUT) :: ray
    logical :: flg
    character(100) :: s


#ifdef ART_DEBUG
    write (s,'(A,I0,A)') '(A,/,A,I4,I8,I8,I18,(1P',NUM_FREQ+1,'E15.5),I8)'
    print s, "(AddMyRay)  myrank   isrc     res            ipix        dist             tau(:)                    gid",&
    "            ", myrank, ray%isrc, ray%Hpix%res, ray%Hpix%ipix,  ray%dist, ray%tau(:), ray%gid
#endif

    !my_ray_listの一番上にrayを格納
    ray%next => my_ray_list
    my_ray_list => ray

    N_myrays = N_myrays + 1 ! for debug
  end subroutine art_AddMyRay

  !-------------------------------------------------------------------------
  ! art_AddSendRay
  ! add rays to send_ray_list
  ! INPUT: isrc, Hpix, dist, tau, gid, rank
  !-------------------------------------------------------------------------
  subroutine   art_AddSendRay(ray, rank)       
    type(t_ray),pointer,intent(INOUT) :: ray
    integer,intent(IN) :: rank                                     ! rank for the grid the ray is entering

    type(t_ray),pointer :: new_ray
    logical :: flg
    character(100) :: s

    !check rank/gid
    if (rank == myrank .or. rank < 0 .or. rank > NPE-1 .or. ray%gid == Undefi) then
       print *, "(AddSendRay) input rank/gid is strange", rank, ray%gid, myrank
       print *, "stopping..."
       stop
    end if

#ifdef ART_DEBUG
    ! 作ったrayの情報をprint
    write (s,'(A,I0,A)') '(A,/,A,I4,I8,I8,I18,(1P',NUM_FREQ+1,'E15.5),I8)'
    print s, "(AddSendRay) myrank  isrc     res            ipix        dist             tau(:)                    gid",&
    "            ", myrank, ray%isrc, ray%Hpix%res, ray%Hpix%ipix,  ray%dist, ray%tau(:), ray%gid
#endif

    !send_ray_list(rank)の一番上にrayを格納
    ray%next => send_ray_list(rank)%ray_list
    send_ray_list(rank)%ray_list => ray

    !送信待ちのrayが増えたのでカウント
    N_send(rank) = N_send(rank) + 1    
  end subroutine art_AddSendRay

  !-------------------------------------------------------------------------
  ! art_CrossGrid
  ! advancing a ray across a grid
  !
  ! ** relation of indices for cells and their boundaries **
  !           x(Imin)          x(Imin+1)            x(Imax)
  ! xbnd(Imin) <--> xbnd(Imin+1) <-->        ...      <--> xbnd(Imax+1)
  !-------------------------------------------------------------------------
  subroutine art_CrossGrid(ray)
!    use modelParameter, only : MP_SubDiskAbs
#ifdef METAL_TRANSFER
    use modelParameter, only : MP_SubDiskAbs, MP_Boxsize, MP_SubStrmAbs, MP_Metallicity ! for debug
#else
    use modelParameter, only : MP_SubDiskAbs, MP_Boxsize, MP_SubStrmAbs ! for debug
#endif

    type(t_ray),pointer,intent(INOUT) :: ray
    integer :: level,i,j,k,gid, surf
    logical :: flg_no_ijk
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    real(kind=DBL_KIND),dimension(MX:MZ) :: ndir !direction vector
    integer :: isrc
    type(t_dirHpix) :: Hpix         
    real(kind=DBL_KIND) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1) :: tau
    real(kind=DBL_KIND) :: theta, phi
    integer :: cell_count ! 通過したセル数
    integer :: lr, mi
    integer :: next_rank, next_gid ! 次のグリッドのrank、gid
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(Imin:Imax+1) :: xbnd ! cell boundaries in x-direction
    real(kind=DBL_KIND),dimension(Jmin:Jmax+1) :: ybnd ! cell boundaries in y-direction
    real(kind=DBL_KIND),dimension(Kmin:Kmax+1) :: zbnd ! cell boundaries in z-direction
    type(t_obRectPhys) :: cell_bnd !cell boundaryの構造体
    logical :: flg_end, flg_split             !終了・分割条件のflag
    real(kind=DBL_KIND) :: tempK, nH, dl, dv
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yhp, khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi, yco !EUV RT    
    real(kind=DBL_KIND) :: rho_dust, dtau_euv, dtau_fuv
#if MODEL_ART == 2
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kh2pd !FUV RT
    real(kind=DBL_KIND) :: dNcH2 ! column density of H2
#ifdef METAL
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: kgfuv, kdph, Td, kdrco, krOII
    real(kind=DBL_KIND) :: dNcCO, dtau_de, dtau_df
#ifdef DUST_NOTCONSTANT
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
#endif
#endif
#endif
#ifdef SUB_GRID_MODEL_DIRECTLIGHT
       integer :: submask_euv 
#endif
!HFADDED

    integer :: sourceCell_flag ! 光源シンクセルか否か
    integer :: sourceSkin_flag ! 光源セルの皮の部分か否か
    integer :: subDisk_flag    ! サブグリッド円盤吸収か否か
    integer :: subStrm_flag    ! サブグリッド電離光子吸収か否か    

    character(100) :: s
#ifdef METAL_TRANSFER
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: mmetal
#endif


    !ray情報の取り出し
    isrc = ray%isrc
    Hpix = ray%Hpix
    dist = ray%dist
    tau = ray%tau
    gid = ray%gid

    !HEALPix <-> 角度 の変換
    call pix2ang_nest(2**(Hpix%res), Hpix%ipix, theta, phi)
    call art_ang2vec(theta,phi,ndir)    

    !current position
    pos(:) = rs_info%spos(:,isrc) + dist * ndir(:)

    !grid/cell情報
    ! if(gid==-1) then
    !    print *, "(art_CrossGrid) gid = -1 (grid not exists). stopping..."
    !    stop
    ! end if
    level = get_level(gid)
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)    
    h(:) = CellWidth(:,level)
    xbnd(:) = x(Imin:Imax+1) - 5d-1 * h(MX) !ゴーストセルがあるのでこれでOK
    ybnd(:) = y(Jmin:Jmax+1) - 5d-1 * h(MY)
    zbnd(:) = z(Kmin:Kmax+1) - 5d-1 * h(MZ)

    !変数の取り出し
    rho => get_Ucomp(MRHO,gid)
    p => get_Ucomp(MP,gid)
    yhn => get_Ucomp(MHN,gid)
    yh2 => get_Ucomp(MH2,gid)    
    yhp => get_Ucomp(MHP,gid)
    yco => get_Ucomp(MCO,gid)
    khpi => get_Ucomp(MKPI,gid)
    heat_hpi => get_Ucomp(MHPI,gid)
#ifdef EXTERNALFORCE
    fx_hpi => get_Ucomp(MXPI,gid)
    fy_hpi => get_Ucomp(MYPI,gid)
    fz_hpi => get_Ucomp(MZPI,gid)    
#endif

#if MODEL_ART == 2 
  kh2pd => get_Ucomp(MKPD,gid)
  #ifdef METAL
    kgfuv => get_Ucomp(MGFUV,gid)
    kdph  => get_Ucomp(MDPH,gid)
    Td    => get_Ucomp(MTD,gid)
    kdrco => get_Ucomp(MDPCO, gid)
    krOII => get_Ucomp(MKPOII, gid)
    #ifdef DUST_NOTCONSTANT
      rhod => get_Ucomp(MDRHO,gid)
    #endif
  #endif
#endif

#ifdef METAL_TRANSFER
    mmetal => get_Ucomp(MMET, gid)
#endif



#ifdef ART_DEBUG
    write (s,'(A,I0,A)') '(A,/,A,I4,A,I4,A,I18,/,A,(1P',NUM_FREQ+1,'E10.2),A,I4,I18,/,A,(1P5E10.2),/,A,(1P6E10.2),/,A)'
    print s,&
         "************************************************************************", &
         "myrank=", myrank,", gid=", gid, ", grid_count_tot=", grid_count_tot,  &
         "dist, tau(:):", dist, tau(:), ",   Hpix:", Hpix%res, Hpix%ipix, &
         "theta, phi, dxdt, dydt, dzdt:", theta, phi, sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta), &
         "grid_geom:", xbnd(Imin), xbnd(Imax+1), ybnd(Jmin), ybnd(Jmax+1), zbnd(Kmin), zbnd(Kmax+1), &
         "************************************************************************"
#endif

    !--------------------- KS DEBUG ----------------------!
    !if (glob_debug_flag_art) then
    !   !print '(A,I4,A,I4,A,I12,A,1P3E10.2,A,1P6E10.2)', &
    !   print *, &
    !     "myrank=", myrank,", gid=", gid, ", grid_count_tot=", grid_count_tot,  &
    !     ", pos(:):", pos(:), ", grid_geom:", xbnd(Imin), xbnd(Imax+1), ybnd(Jmin), ybnd(Jmax+1), zbnd(Kmin), zbnd(Kmax+1)
    !end if
    !--------------------- KS DEBUG ----------------------!    

    !rayがいるセルの初期位置 (gridは既知)
    call FindInitialCell(pos, xbnd,ybnd,zbnd, i,j,k)

    !-------------------------------------------------------------------------
    !セル通過を繰り返すことでグリッドを通過
    !-------------------------------------------------------------------------
    cell_count = 0
    do

#ifdef ART_DEBUG       
       call art_WriteTrajectory(ray)
#endif

       ! print *, ""
       ! print *, "N_cell=",cell_count
       cell_bnd%left(MX) = xbnd(i)
       cell_bnd%right(MX) = xbnd(i+1)
       cell_bnd%left(MY) = ybnd(j)
       cell_bnd%right(MY) = ybnd(j+1)
       cell_bnd%left(MZ) = zbnd(k)
       cell_bnd%right(MZ) = zbnd(k+1)

       call art_FindIntersection(pos,ndir,cell_bnd,surf)
       
       !update dist
       dist = sqrt(dot_product(pos-rs_info%spos(:,isrc),pos-rs_info%spos(:,isrc)))
       dl =  dist - ray%dist !進んだ距離
       ray%dist = dist       !new dist
       

       ! check wheter the current cell is radiation source cell
       call check_is_sourceCell(x(i),y(j),z(k),rs_info%spos(:,isrc),sourceCell_flag)

       ! check wheter the current cell is skin of radiation source cell
       call check_is_sourceSkin(x(i),y(j),z(k),rs_info%spos(:,isrc),sourceSkin_flag)

       !subgrid diskによる吸収を考慮するか否か
       if (MP_SubDiskAbs >= 1) then
          ! check the effect of absorption by subgrid disk
          call check_abs_subDisk(x(i),y(j),z(k),rs_info%spos(:,isrc),rho(i,j,k),yhn(i,j,k),yh2(i,j,k),subDisk_flag)
          ! subgrid diskによりrayは完全に吸収されると仮定
          if (subDisk_flag == 1) then
             call art_EndRay(ray)
             exit
          end if
       end if

       !------------------------------------- KS DEVELOPING -----------------------------------!
       !strmgren半径を用いた判定によるサブグリッド吸収を考慮するか否か
       !if (MP_SubStrmAbs == 1) then
       !   !check absoprtion of ionizaing photons within sink region (Stromgren check)
       !   call check_strm_in_sink(x(i),y(j),z(k),&
       !        rs_info%spos(:,isrc),rs_info%lum(isrc),rs_info%x_euv(isrc),rho(i,j,k),subStrm_flag)
       !   if (subStrm_flag == 1) then
       !      tau(BIN_EUV) = tau(BIN_EUV) + 1000.  !電離光子のtauを十分大きくして、シンク内で全部吸収されていることを実現
       !   end if
       !end if

       !------------------------------------- KS DEVELOPING -----------------------------------!       

       !------------------------------------- HF DEVELOPING -----------------------------------!
#ifdef SUB_GRID_MODEL_DIRECTLIGHT
        if(sourceSkin_flag == 1 .and. tau(BIN_EUV) < 1000.) then ! 光源セルの皮でsubgridモデルのmaskingを行うか判定
          call check_masksub_euv(Hpix, rho(i,j,k), rs_info%lum(isrc), rs_info%Trad(isrc),rs_info%x_euv(isrc), rs_info%mskr(isrc), submask_euv)
          if (submask_euv == 1) then
            tau(BIN_EUV) = tau(BIN_EUV) + 1000.  !電離光子のtauを十分大きくして、シンク内で全部吸収されていることを実現
          endif
        endif
#endif
       !---------------------------------------------------------------------------------------!




       !get dtau_euv/dNcH2
       !if (sourceCell_flag == 0) then !光源セルでなければ通常のray tracing
       if (sourceCell_flag == 0 .and. sourceSkin_flag == 0) then !光源セル（皮の部分含む）でなければ通常のray tracing          

          call art_GetDtau(rho(i,j,k),yhn(i,j,k),dl,isrc,dtau_euv)   !EUV RT

          ! Dust density
#ifdef DUST_NOTCONSTANT  
          rho_dust = rhod(i,j,k)
#else
  #ifdef METAL_TRANSFER
          rho_dust = fdust_solar*mmetal(i,j,k)
  #else
          rho_dust = fdust_solar*MP_Metallicity
  #endif
#endif


#if MODEL_ART == 2
        call art_GetDNcH2(rho(i,j,k),yh2(i,j,k),dl,dNcH2) !FUV RT
  #ifdef METAL  
        call art_GetDNcnH(rho(i,j,k), rho_dust , dl, isrc, dtau_de, dtau_df)  ! dust opacity for deuv dfuv
        call art_GettauCO(rho(i,j,k), yco(i,j,k),dl,dNcCO)              ! !!cuation!! dNcnHはダストの融解率と柱密度の積
  #endif
#endif ! MODEL_ART
       else               !光源セル（皮の部分含む）ではoptically thinを仮定
          dtau_euv=0d0
#if MODEL_ART == 2
        dNcH2=0d0
  #ifdef METAL
        dtau_de = 0d0
        dtau_df = 0d0
        dNcCO = 0d0
  #endif
#endif ! MODEL_ART
       end if
       
       !update rad_field
       dv = h(MX)*h(MY)*h(MZ) !cell volume

       !store radiation effect
       if (sourceCell_flag == 0) then !光源セルでなければ光反応率を格納（光源セルでの光反応は無視, 皮の部分では光反応を考慮, KS TODO）

          call art_StoreRadEffect(isrc,Hpix,tau,ndir,dl,dv,khpi(i,j,k),heat_hpi(i,j,k) &
            , rho(i,j,k), yhn(i,j,k), dtau_euv &
#ifdef EXTERNALFORCE
               ,fx_hpi(i,j,k),fy_hpi(i,j,k),fz_hpi(i,j,k) &
#endif !EXTERNALFORCE
#if MODEL_ART == 2  
             ,kh2pd(i,j,k) &
  #ifdef METAL
             , kgfuv(i,j,k), kdph(i,j,k), kdrco(i,j,k), krOII(i,j,k), dtau_de, dtau_df &
  #endif
#endif !MODEL_ART
               )

       end if

       !update tau(0) <- tau of EUV
       tau(BIN_EUV)     = tau(BIN_EUV) + dtau_euv !new tau
       ray%tau(BIN_EUV) = tau(BIN_EUV) !ここでrayのtauも更新しておかないと、art_SplitRayの際にうまく情報が伝わらないので注意

       ! if (x(i) > 0 .and. &
       !      abs(y(j)) < abs(y(5)-y(4)) .and. y(j) > 0. .and. &
       !      abs(z(k)) < abs(z(5)-z(4)) .and. z(k) > 0. ) then
       !    print '(A, (1P10E15.7))', "KS DEBUG B", x(i),y(j),z(k),dist, dtau_euv, tau(0)
       ! endif


#if MODEL_ART == 2
      !update tau(1) <- H2 column density for FUV RT
      tau(BIN_H2)     = tau(BIN_H2) + dNcH2 !new column density
      ray%tau(BIN_H2) = tau(BIN_H2)

  #ifdef METAL
      tau(BIN_DEUV)     = tau(BIN_DEUV) + dtau_de 
      ray%tau(BIN_DEUV) = tau(BIN_DEUV)

      tau(BIN_DFUV)     = tau(BIN_DFUV) + dtau_df 
      ray%tau(BIN_DFUV) = tau(BIN_DFUV)

      tau(BIN_CO) = tau(BIN_CO) + dNcCO
      ray%tau(BIN_CO) = tau(BIN_CO)
  #endif
#endif


       ! !------------------- KS DEBUG ---------------!
       ! if (isrc == 1 .and. dist > SinkRadius-maxval(CellWidth(:,LevelMax)) .and. dist < SinkRadius+8.*maxval(CellWidth(:,LevelMax))) then
       !    if (ndir(0) < 0 .and. abs(ndir(1)) < 0.1 .and. abs(ndir(2)) < 0.1) then
       !    write (s,'(A,I0,A)') '(A,/,A,I4,A,I4,A,I18,/,A,(1P',NUM_FREQ+1,'E10.2),A,I4,I18,A,1P3E10.2,/,4(A,1PE10.2),2(A,L3),/,3(A,1PE10.2),/,A)'
       !    print s,&
       !   "************************************************************************", &
       !   "myrank=", myrank,", gid=", gid, ", grid_count_tot=", grid_count_tot,  &
       !   "dist, tau(:):", dist, tau(:), ",   Hpix:", Hpix%res, Hpix%ipix, ", ndir:", ndir(:),&
       !   "nH=", rho(i,j,k)*Unit_rho/(cgs_mh*(1.+4.*yHe)), ", yhn=", yhn(i,j,k),", dl=", dl*Unit_au, ", dist=", dist*Unit_au,", srcCell=",sourceCell_flag,", srcSkin=",sourceSkin_flag,&
       !   "dtau_euv=", dtau_euv,", tau_euv_before=", tau(0)-dtau_euv,", khpi=",khpi(i,j,k),&
       !   "************************************************************************"
       ! end if
       ! end if
       ! !------------------- KS DEBUG ---------------!       


       if (surf == NIR) then
          i = i+1
       else if (surf == NIL) then
          i = i-1
       else if (surf == NJR) then
          j = j+1
       else if (surf == NJL) then
          j = j-1
       else if (surf == NKR) then
          k = k+1
       else if (surf == NKL) then
          k = k-1
       end if

       ! print *, "i,j,k: ", i,j,k
       cell_count = cell_count + 1

       !rayの終了条件をチェック
       call art_EndCond(isrc,Hpix,dist,tau,flg_end,ndir_in=ndir)    
       if (flg_end) then
#ifdef ART_DEBUG
          print '(/,A,I4,A,/)',"***    ray-ending condition is satisfied.  (myrank=",myrank,", CrossGrid)  ***"
#endif
          call art_EndRay(ray)
          ! deallocate(ray)
          exit
       end if

       !rayの分割条件をチェック
       call art_SplitCond(ray, minval(h), flg_split)
       if (flg_split) then
#ifdef ART_DEBUG
          print '(/,A,I4,A,/)',"***    ray-splitting condition is satisfied.  (myrank=",myrank,", CrossGrid)  ***"
#endif
          call art_SplitRay(ray)          
          exit
       end if


       !グリッドの端まで到達したかをチェック、端まで来てたら次のグリッドを取得
       if (i < Imin .or. i > Imax .or. j < Jmin .or. j > Jmax .or. k < Kmin .or. k > kmax) then
          lr = s2lr(surf)
          mi = s2dir(surf)
          call GetNextGrid(lr,mi,gid,myrank,pos,next_gid,next_rank)          

          ray%gid = next_gid !次のグリッドのgid

#ifdef ART_DEBUG
          print '(/,4(A,I4),/)', "(CrossGrid) myrank=", myrank,", gid=", gid,", next_rank=", next_rank,", next_gid=", next_gid
#endif
          N_myrays = N_myrays - 1 ! for debug
          if (next_rank == myrank) then
             call art_AddMyRay(ray) !rayは再利用
          else
             !check rank/gid
             if (next_rank < 0 .or. next_rank > NPE-1 .or. next_gid == Undefi) then
                print *, "(art_CrossGrid) rank/gid is starnge before callling art_AddSendRay: ", next_rank, next_gid
                print *, "myrank, pos, rs_info%spos(:,isrc), dist", myrank, pos, rs_info%spos(:,isrc), dist
                print *, "stopping..."
                stop
             end if
             call art_AddSendRay(ray, next_rank)
          end if
          exit
       end if
    end do

  contains
    
    !find initial cell
    subroutine FindInitialCell(pos, xbnd,ybnd,zbnd, i,j,k)
      real(kind=DBL_KIND),dimension(MX:MZ),intent(INOUT) :: pos     ! start point
      real(kind=DBL_KIND),dimension(Imin:Imax+1),intent(IN) :: xbnd ! cell boundaries in x-direction
      real(kind=DBL_KIND),dimension(Jmin:Jmax+1),intent(IN) :: ybnd ! cell boundaries in y-direction
      real(kind=DBL_KIND),dimension(Kmin:Kmax+1),intent(IN) :: zbnd ! cell boundaries in z-direction
      integer,intent(OUT) :: i,j,k                                  ! indices for the initial cell
      integer :: ii,jj,kk
      real(kind=DBL_KIND),parameter :: eps = 1d-8                  !small number (margin for inside/outside check)

      !check wheter the pos is included in the assumed grid
      if ( pos(MX) < xbnd(Imin) - eps .or. pos(MX) > xbnd(Imax+1) + eps .or. &
           pos(MY) < ybnd(Jmin) - eps .or. pos(MY) > ybnd(Jmax+1) + eps .or. &
           pos(MZ) < zbnd(Kmin) - eps .or. pos(MZ) > zbnd(Kmax+1) + eps ) then
       !print '(/,A,I4,A,/,A,(1P3E20.12),/,A,(1P3E20.12),/,A,(1P3E20.12),/,A)', &
       print *, &
            "(FindInitialCell) pos is not included in the grid (myrank=", myrank,"). something is wrong...", &
            "pos:",pos, "grid_bnd (left):",xbnd(Imin),ybnd(Jmin),zbnd(Kmin), &
            "grid_bnd (right):",xbnd(Imax+1),ybnd(Jmax+1),zbnd(Kmax+1), "stopping..."
         stop
      end if
      !x-direction
      i = Imax
      do ii = Imin, Imax-1
         if (pos(MX) < xbnd(ii+1)) then
            i = ii
            exit
         end if
      end do
      !y-direction
      j = Jmax
      do jj = Jmin, Jmax-1
         if (pos(MY) < ybnd(jj+1)) then
            j = jj
            exit
         end if
      end do
      !z-direction
      k = Kmax
      do kk = Kmin, Kmax-1
         if (pos(MZ) < zbnd(kk+1)) then
            k = kk
            exit
         end if
      end do

      ! print *, xbnd(i), pos(MX), xbnd(i+1)
      ! print *, ybnd(j), pos(MY), ybnd(j+1)
      ! print *, zbnd(k), pos(MZ), zbnd(k+1)
      
      ! print *, i,j,k

    end subroutine FindInitialCell

    !from surface code to direction
    integer function s2dir(surf)
      integer, intent(IN) :: surf
      s2dir = surf/2 !MX/MY/MZ = 0/1/2
    end function s2dir
    !from surface code to left/right
    integer function s2lr(surf)
      integer, intent(IN) :: surf
      s2lr = mod(surf,2) !LEFT/RIGHT = 0/1  
    end function s2lr

    !get next grid (おじさんグリッドしか存在しない場合、おいグリッドが存在する場合に対応)
    subroutine GetNextGrid(lr,mi,gid,rank,pos,next_gid,next_rank)
      integer, intent(IN) :: lr !Left or Right
      integer, intent(IN) :: mi !MX, MY, or MZ
      integer, intent(IN) :: gid, rank !gid and rank for current grid
      real(kind=DBL_KIND),dimension(MX:MZ),intent(INOUT) :: pos ! current position (used to know the position of child grid in the parent grid)
      integer, intent(OUT) :: next_gid, next_rank !gid and rank for next grid
      integer :: tmp_gid, tmp_rank !gid and rank of tmp grid
      

      integer ::  lri, lrj, lrk, sgn
      real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
      real(kind=DBL_KIND),dimension(MX:MZ) :: rcen,rwidth
      real(kind=DBL_KIND) :: xcen, ycen, zcen
      
      !隣に同世代のグリッドが無い場合
      if (NeighborGid(lr,mi,gid,myrank) == Undefi) then
         tmp_gid = ParentGid(gid,myrank)                !親の世代で隣のグリッドを探査 (隣とのレベル差は1以下なので存在は保証される)
         tmp_rank = ParentRank(gid,myrank)
         !rayが計算領域を超えようとしたときにエラーを吐いて止める
         if (tmp_gid == Undefi .or. tmp_rank == Undefi) then
            !print '(A,5I6,/,A,/,A)', "(GetNextGrid: error) gid, myrank, get_level(gid), tmp_gid, tmp_rank: ",&
            print *, "(GetNextGrid: error) gid, myrank, get_level(gid), tmp_gid, tmp_rank: ",&
                 gid, myrank, get_level(gid), tmp_gid, tmp_rank,&
                 "This error may occur when a ray try to go beyond the region where grids exist.",&
                 "stopping..."
            stop
         end if
         next_rank = NeighborRank(lr,mi,tmp_gid,tmp_rank)
         next_gid = NeighborGid(lr,mi,tmp_gid,tmp_rank)
         return
      end if

      !隣に同世代のグリッドがある場合
      next_rank = NeighborRank(lr,mi,gid,myrank)
      next_gid = NeighborGid(lr,mi,gid,myrank)

      !隣のグリッドに子グリッド存在する場合
      if (ChildGid(Left, Left, Left, next_gid, next_rank) /= Undefi) then

         tmp_gid = next_gid     !親グリッドのgid
         tmp_rank = next_rank   !親グリッドのrank
         
         !tmp_rankが別rankの可能性もあるので、tmp_gidのセル情報は使わない
         x => get_Xp(gid)
         y => get_Yp(gid)
         z => get_Zp(gid)

         rcen(MX) = (x(Imin)+x(Imax))/2d0
         rcen(MY) = (y(Jmin)+y(Jmax))/2d0
         rcen(MZ) = (z(Kmin)+z(Kmax))/2d0
         rwidth(MX) = x(Imax+1)-x(Imin) !グリッド幅は両端のセル中心の距離より大きいことに注意
         rwidth(MY) = y(Jmax+1)-y(Jmin)
         rwidth(MZ) = z(Kmax+1)-z(Kmin)

         if (lr==Left) then
            sgn = -1
         else
            sgn = 1
         end if
         rcen(mi) = rcen(mi) + sgn*rwidth(mi) !隣のグリッドの中心

         !posがどのchild_gridに含まれるか
         if(pos(MX)<rcen(MX)) then
            lri = Left
         else
            lri = Right
         end if
         if(pos(MY)<rcen(MY)) then
            lrj = Left
         else
            lrj = Right
         end if
         if(pos(MZ)<rcen(MZ)) then
            lrk = Left
         else
            lrk = Right
         end if
         next_rank = ChildRank(lri, lrj, lrk, tmp_gid, tmp_rank)
         next_gid = ChildGid(lri, lrj, lrk, tmp_gid, tmp_rank)
      end if
   
    end subroutine GetNextGrid

    !
    !現在通過中のセルが光源シンクセルかを判定（他のシンクセルは普通のセルと同様に扱う）
    !sinkセルと判定されたらsourceCell_flag = 1, それ以外は = 0 を返す
    !
    subroutine check_is_sourceCell(xi,yj,zk,spos,sourceCell_flag)
      real(kind=DBL_KIND),intent(IN) :: xi,yj,zk ! coordinates for cell center
      real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: spos ! source position
      integer, intent(OUT) :: sourceCell_flag
      real(kind=DBL_KIND) :: r2
      r2 = (xi-spos(MX))**2 + (yj-spos(MY))**2 + (zk-spos(MZ))**2
      if (r2 < SinkRadius**2) then
         sourceCell_flag = 1
      else
         sourceCell_flag = 0
      end if
      
    end subroutine check_is_sourceCell

    !
    !現在通過中のセルが光源シンクセルの皮の部分かを判定
    !sinkセルの皮の部分と判定されたらsinkCell_flag = 1, それ以外は = 0 を返す
    !皮の部分ではoptically thinを仮定
    ! <- シンクセル内は光が届かないと仮定しているので、シンク内の中性ガスが漏れ出すと電離光子が抜け出せなくなる
    !
    subroutine check_is_sourceSkin(xi,yj,zk,spos,sourceSkin_flag)
      use modelParameter, only : MP_SubAbsCell

      real(kind=DBL_KIND),intent(IN) :: xi,yj,zk ! coordinates for cell center
      real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: spos ! source position
      integer, intent(OUT) :: sourceSkin_flag
      real(kind=DBL_KIND) :: r2
      real(kind=DBL_KIND), save :: r_in, r_out ! バッファ領域の内縁と外縁
      logical,save :: first_call = .True.

      ! 初期化
      if (first_call) then
         r_in = SinkRadius
         r_out = SinkRadius + MP_SubAbsCell*maxval(CellWidth(:,LevelMax))
         first_call = .False.
      end if

      !光源からの距離
      r2 = (xi-spos(MX))**2 + (yj-spos(MY))**2 + (zk-spos(MZ))**2

      !光源セルの皮の部分だったらsourceSkin_flag=1、それ以外はsourceSkin_flag=0
      if (r_in**2 < r2 .and. r2 < r_out**2 ) then
         sourceSkin_flag = 1
      else
         sourceSkin_flag = 0
      end if
    end subroutine check_is_sourceSkin
    

    ! -------------------------------------------------------------------!
    !      subgrid model for EUV photons developed by HF
    ! -------------------------------------------------------------------!

#ifdef SUB_GRID_MODEL_DIRECTLIGHT
    subroutine check_masksub_euv(Hpix, rho, lum, Trad, x_euv, mskr, submask_euv)

      real(kind=DBL_KIND), intent(IN) :: rho, lum, Trad, x_euv
      real(kind=DBL_KIND), intent(INOUT) :: mskr
      type(t_dirHpix),intent(IN) :: Hpix
      integer :: res
      integer(kind=LLONG_KIND) :: N_pix
      integer, intent(INOUT) :: submask_euv

      real(kind=DBL_KIND), save :: r_in3, r_in
      real(kind=DBL_KIND) :: Scr, Sion, nH, r_star
      real(kind=DBL_KIND), parameter :: alpha_B = 2.6d-13 ! cm^3 s^-1
      logical,save :: first_call = .True.

      ! 初期化
      if (first_call) then
         r_in  = SinkRadius*Unit_l
         r_in3 = SinkRadius*SinkRadius*SinkRadius*Unit_l3 ! cm^3
         first_call = .False.
      end if

      ! number ionizing photons
      Sion = lum*x_euv

      ! crtical number density
      nH  = rho*Unit_rhomu 
      r_star = sqrt(lum/(4.*Pi * cgs_sigma*Trad**4))
      !Scr = 4.d0*Pi*alpha_B*nH*nH*r_in3*dlog(r_in/r_star)
      Scr = 4.0*Pi/3.0*alpha_B*nH*nH*r_in3


      if(Sion < Scr) then
          submask_euv = 1
          res = Hpix%res
          N_pix = 12_8 * 4_8**res
          mskr  = mskr + 1.d0/N_pix
      else
          submask_euv = 0
      endif
      !print *,"subgrid", submask_euv, Scr, Sion, mskr
    end subroutine check_masksub_euv
#endif


    !
    !現在通過中のセルの情報からrayがサブグリッド円盤によって吸収済みかを判定
    !サブグリッド円盤によって吸収済みと判定されたらsourceCell_flag = 1, それ以外は = 0 を返す
    !
    !シンク領域に隣接したバッファ領域にある程度密度の高い電離していないガスが存在するとき、
    !その内側に中性水素 or 水素分子円盤があってrayをその内側で吸収済みと仮定
    !具体的な条件としては、nH > MP_SubDiskNcr cm^-3、かつ、y(H)+2y(H2) > 5e-1 をひとまず採用 (KS TODO)
    !シンク半径からMP_SubAbsCell (defualt 2) セルの厚みに入るセルで、サブグリッド円盤による吸収の効果を考慮
    !
    subroutine check_abs_subDisk(xi,yj,zk,spos,rho,yhn,yh2,subDisk_flag)
      use modelParameter, only : MP_SubAbsCell, MP_SubDiskNcr

      real(kind=DBL_KIND),intent(IN) :: xi,yj,zk ! coordinates for cell center
      real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: spos ! source position
      real(kind=DBL_KIND),intent(IN) :: rho, yhn, yh2
      integer, intent(OUT) :: subDisk_flag
      real(kind=DBL_KIND) :: r2
      real(kind=DBL_KIND), save :: r_in, r_out ! バッファ領域の内縁と外縁
      real(kind=DBL_KIND),save :: n_cr, y_cr !criticalな密度と中性+分子割合
      logical,save :: first_call = .True.

      ! 初期化
      if (first_call) then
         n_cr = MP_SubDiskNcr ! criticalな密度はシンク密度の1/10 (MP_SubDiskAbs=1)、1/100 (MP_SubDiskAbs>1)
         y_cr= 5d-1       ! criticalな中性+分子割合
         r_in = SinkRadius
         r_out = SinkRadius + MP_SubAbsCell*maxval(CellWidth(:,LevelMax))
         first_call = .False.
      end if
      
      !何も無ければ、subDisk_flag=0
      subDisk_flag = 0

      !光源からの距離
      r2 = (xi-spos(MX))**2 + (yj-spos(MY))**2 + (zk-spos(MZ))**2

      !bufer領域のセルにおいて、
      if (r_in**2 < r2 .and. r2 < r_out**2 ) then
         ! 密度が高く、主成分が中性水素または水素分子のとき
         if ((rho*Unit_rhomu) > n_cr .and. yhn+2.d0*yh2 > y_cr) then
            subDisk_flag = 1 ! サブグリッド円盤に吸収済みのrayと判定
#ifdef ART_DEBUG
       print '(A,I4,5(A,1P1E10.2))', "(ART subDisk) ray end: myrank=", &
            myrank,", r=", sqrt(r2),", SinkRadius=", SinkRadius,&
            ", nH=", rho*Unit_rhomu,", n_cr=", n_cr,", yhn=", yhn,", yh2=", yh2
#endif
         end if
      end if
    end subroutine check_abs_subDisk

    !
    ! Stromgren半径がシンク半径より小さいかを判定
    ! 小さい場合は、シンク内で電離光子は完全に吸収されて出てこれないと仮定
    ! シンク粒子内での密度は、シンクに隣接するセルの密度と等しいと仮定して見積もり
    !
    subroutine check_strm_in_sink(xi,yj,zk,spos,lum,x_euv,rho,subStrm_flag)
      use modelParameter, only : MP_SubAbsCell

      real(kind=DBL_KIND),intent(IN) :: xi,yj,zk ! coordinates for cell center
      real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: spos ! source position
      real(kind=DBL_KIND),intent(IN) :: lum,x_euv
      real(kind=DBL_KIND),intent(IN) :: rho
      integer, intent(OUT) :: subStrm_flag

      real(kind=DBL_KIND) :: Ndot_ion,nH,r_Strm,r2
      real(kind=DBL_KIND),parameter :: alphaB_HII = 8.13d-14 !alpha_B [cm^3 s^-1] for T_HII=4e4K gas 
      real(kind=DBL_KIND), save :: r_in, r_out ! バッファ領域の内縁と外縁

      logical,save :: first_call = .True.

      ! 初期化
      if (first_call) then
         r_in = SinkRadius
         r_out = SinkRadius + MP_SubAbsCell*maxval(CellWidth(:,LevelMax))
         first_call = .False.
      end if
      
      !何も無ければ、subDisk_flag=0
      subStrm_flag = 0

      !光源からの距離
      r2 = (xi-spos(MX))**2 + (yj-spos(MY))**2 + (zk-spos(MZ))**2

      !bufer領域のセルにおいて、
      if (r_in**2 < r2 .and. r2 < r_out**2 ) then
         Ndot_ion = lum * x_euv
         nH = rho*Unit_rhomu
         r_Strm = (Ndot_ion / ((4.*PI/3.) * nH**2 * alphaB_HII))**(1./3) !nHとNdot_ionに対応するStromgren半径[cm]

         ! Stromgren半径がSink半径より小さければ
         if (r_Strm < SinkRadius * Unit_l) then
            subStrm_flag = 1 ! サブグリッド円盤に吸収済みのrayと判定
#ifdef ART_DEBUG
       print '(A,I4,4(A,1P1E10.2))', "(ART subStrm) EUV absorbed: myrank=", &
            myrank,", r=", sqrt(r2)*Unit_au,", SinkRadius=", SinkRadius*Unit_au,&
            ", nH=", rho*Unit_rhomu,", r_Strm=", r_Strm/cgs_au
#endif
         end if
      end if

    end subroutine check_strm_in_sink
      
  end subroutine art_CrossGrid

  !-------------------------------------------------------------------------
  ! find intersection of ray and cell surface
  ! INPUT:
  !     pos = current position
  !     drdt = direction vector
  !     gid = grid id
  !     i,j,k = cell id within the grid
  ! OUTPUT:
  !     pos = position of the intersection
  !     surf = which surface
  !-------------------------------------------------------------------------
  subroutine art_FindIntersection(pos,drdt,cell_bnd,surf)
    use grid
    use unit
    ! ray trajectory を pos(t) = pos(0) + (sin_th*cos_ph, sin_th*sin_ph, cos_th) * t と記述
    real(kind=DBL_KIND),dimension(MX:MZ),intent(INOUT) :: pos ! start/intersection point
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN)    :: drdt ! normalized direction vector
    type(t_obRectPhys),intent(IN) :: cell_bnd !left & right surfs in MX/MY/MZ directions
    integer,intent(OUT) :: surf ! direction of intersection surface (NIL/NIR/NJL/NJR/NKL/NKR)
    integer :: m
    real(kind=DBL_KIND) :: theta, phi
    real(kind=DBL_KIND),dimension(MX:MZ) :: t_cand ! candidate for t@intersection
    real(kind=DBL_KIND) :: t_int ! t@intersection
!    real(kind=DBL_KIND),parameter :: t_LARGE = Huge(1.d0)  ! Large number (used in the smallest t search)
    real(kind=DBL_KIND),parameter :: eps = 1d-8     ! small number (margin for inside/outside check, normalized by cell size)
    real(kind=DBL_KIND),parameter :: eps_parallel = 1d-10 ! small number (used for parallel check)
    real(kind=DBL_KIND) :: h  !cell size

    h = minval(cell_bnd%right(:) - cell_bnd%left(:))

    !check wheter the pos is included in the assumed cell
    if ( pos(MX) < cell_bnd%left(MX)- eps*h .or. pos(MX) > cell_bnd%right(MX) + eps*h .or. &
         pos(MY) < cell_bnd%left(MY)- eps*h .or. pos(MY) > cell_bnd%right(MY) + eps*h .or. &
         pos(MZ) < cell_bnd%left(MZ)- eps*h .or. pos(MZ) > cell_bnd%right(MZ) + eps*h ) then
       !print '(/,A,I4,A,/,A,(1P3E20.12),/,A,(1P3E20.12),/,A,(1P3E20.12),/,A,(1P3E20.12),/,A)', &
       print *, &
            "(FindIntersection) pos is not included in the cell (myrank=", myrank,"). something is wrong...", &
            "pos:",pos, "cell_bnd%left:",cell_bnd%left,"cell_bnd%right:",cell_bnd%right, &
            "eps, h, eps*h", eps, h, eps*h,  "stopping..."
       stop
    end if



    ! search intersection（posに浮動小数点誤差が入る可能性があることに注意）
    do m = MX, MZ
       if (abs(drdt(m)) < eps_parallel) then !面と平行に進む場合は交わらないとする（0割り防止）
          t_cand(m) = Huge(1.d0)
       else if (drdt(m) > 0) then !正の方向に進む場合
          t_cand(m) = (cell_bnd%right(m) - pos(m)) / drdt(m)
       else if (drdt(m) < 0) then !負の方向に進む場合
          t_cand(m) = (cell_bnd%left(m) - pos(m)) / drdt(m)
       end if
    end do

     t_int = minval(t_cand)

    if (t_int < 0) then
       if (abs(t_int) > eps*h) then  !!注意
          !print '(A,1P1E20.12,A,/,6(A,1P3E20.12,/),A)',&
          print *,&
               "(FindIntersection) t_int =",t_int," is negative.    something is wrong.", &
               "eps, h, eps*h: ", eps, h, eps*h, "t_cand: ", t_cand, &
               "drdt: ",drdt,"pos: ",pos,"cell_bnd%left: ",cell_bnd%left,"cell_bnd%right: ",cell_bnd%right, "stopping..."
          stop
       else
!#ifdef ART_DEBUG
          !print '(A,1P1E20.12,A,/,6(A,1P3E20.12,/))', "(FindIntersection) t_int =",t_int, &
          !print * , "(FindIntersection) t_int =",t_int, &
          !     " is negative but small. Somenthing strange might happen but keep going...", &
          !     "eps, h, eps*h: ", eps, h, eps*h, "t_cand: ", t_cand, &
          !     "drdt: ",drdt,"pos: ",pos,"cell_bnd%left: ",cell_bnd%left,"cell_bnd%right: ",cell_bnd%right
               !"modifying t_int to 0 and going ahead..."
!#endif
          !t_int = 0
       end if
    end if


    pos(:) =  pos(:) + drdt(:)*t_int ! position of the intersection

    ! direction of intersection surface (warning: minloc return not index but position that starts from 1)
    if (minloc(t_cand,1) == 1 .and. drdt(MX) > 0) then
       surf = NIR
    else if (minloc(t_cand,1) == 1 .and. drdt(MX) < 0) then
       surf = NIL
    else if (minloc(t_cand,1) == 2 .and. drdt(MY) > 0) then
       surf = NJR
    else if (minloc(t_cand,1) == 2 .and. drdt(MY) < 0) then
       surf = NJL
    else if (minloc(t_cand,1) == 3 .and. drdt(MZ) > 0) then
       surf = NKR
    else if (minloc(t_cand,1) == 3 .and. drdt(MZ) < 0) then
       surf = NKL
    else
       print *, "(FindIntersection)  surf not able to be determined", minloc(t_cand,1), drdt(:)
       print *, "stopping..."
       stop
    end if

    ! print *, "t_int", t_int
    ! print *, "cell_bnd%left/right:", cell_bnd%left, cell_bnd%right
    ! print *, "newpos:", pos
    ! print *, "surf:", surf
  end subroutine art_FindIntersection

  !-------------------------------------------------------------------------
  ! get dtau_euv from cell data and dl
  !-------------------------------------------------------------------------
  subroutine art_GetDtau(rho,yhn,dl,isrc,dtau)
    real(kind=DBL_KIND),intent(IN) :: rho,yhn,dl
    integer ,intent(IN) :: isrc
    real(kind=DBL_KIND),intent(OUT) :: dtau
    real(kind=DBL_KIND) :: nH,dl_cgs,alpha_euv

    !nH = rho*Unit_rho/(cgs_mh*(1.+4.*yHe)) !in cm^-3
    alpha_euv = rs_info%alpha_euv(isrc)    !in cm^2
    !dl_cgs=dl*Unit_l
    !dtau = nH*yhn*dl_cgs*alpha_euv !注目しているray segmentに対応するtau (ひとまず水素分子は電離しないと仮定; KS TODO)
    
    dtau = rho*yhn*dl*alpha_euv*Unit_lrhomu !注目しているray segmentに対応するtau (ひとまず水素分子は電離しないと仮定; KS TODO)

  end subroutine art_GetDtau

  !-------------------------------------------------------------------------
  ! get dNcH2 from cell data and dl
  !-------------------------------------------------------------------------
  subroutine art_GetDNcH2(rho,yh2,dl,dNcH2)
    real(kind=DBL_KIND),intent(IN) :: rho,yh2,dl
    real(kind=DBL_KIND),intent(OUT) :: dNcH2
    real(kind=DBL_KIND) :: nH,dl_cgs

    !nH = rho*Unit_rho/(cgs_mh*(1.+4.*yHe)) !in cm^-3
    !dl_cgs=dl*Unit_l
    !dNcH2 = nH*yh2*dl_cgs              !ray segmentのH2柱密度

    dNcH2 = rho*yh2*dl*Unit_lrhomu              !ray segmentのH2柱密度
  end subroutine art_GetDNcH2


#ifdef METAL ! HFADDED

  !-------------------------------------------------------------------------
  ! get column local optical depth of dust grain
  !-------------------------------------------------------------------------
  subroutine art_GetDNcnH(rho, rho_dust, dl, isrc, dtau_de, dtau_df)
    real(kind=DBL_KIND),intent(IN) :: rho, dl, rho_dust
    integer ,intent(IN) :: isrc
    real(kind=DBL_KIND),intent(OUT) :: dtau_de, dtau_df
    real(kind=DBL_KIND) :: nH,dl_cgs, eps_d, sigd_euv, sigd_fuv, dNcnH

    ! dust sublimation rate x metallicity          
    eps_d = rho_dust / fdust_solar

    ! dust grain cross section
    sigd_euv = rs_info%sig_euv(isrc) ! [cm^2]
    sigd_fuv = rs_info%sig_fuv(isrc) ! [cm^2]

    dNcnH = rho*dl*Unit_lrhomu       ! [cm^-2] 柱密度

    dtau_de = dNcnH * sigd_euv * eps_d       ! [noD]
    dtau_df = dNcnH * sigd_fuv * eps_d       ! [noD]

  end subroutine art_GetDNcnH

  subroutine art_GettauCO(rho, yco,dl,dNcCO) ! !!cuation!! dNcnHはダストの融解率と柱密度の積
    real(kind=DBL_KIND),intent(IN) :: rho, yco,dl
    real(kind=DBL_KIND),intent(OUT):: dNcCO

    dNcCO = rho*yco*dl*Unit_lrhomu  ! COの柱密度
    
  end subroutine art_GettauCO 

#endif




  !-------------------------------------------------------------------------
  ! store ray contribution to a cell
  !
  ! get optically thin rate here and adjust it in chemistry module to
  ! take into account actual chemical abundance (Tielens&Hollenbach?)
  !-------------------------------------------------------------------------
  subroutine art_StoreRadEffect(isrc,Hpix,tau,ndir,dl,dv,khpi,heat_hpi &
      , rho, yhn, dtau_euv &
#ifdef EXTERNALFORCE
       ,fx_hpi,fy_hpi,fz_hpi&
#endif !EXTERNALFORCE
#if MODEL_ART == 2  
     ,kh2pd &
  #ifdef METAL
     ,kgfuv, kdph, krdco, rOII, dtau_de, dtau_df &
  #endif
#endif ! MODEL_ART   
    )

    ! type(t_ray),pointer,intent(IN) :: ray
    integer,intent(IN) :: isrc
    type(t_dirHpix),intent(IN) :: Hpix
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1),intent(IN) :: tau
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: ndir !direction vector
    real(kind=DBL_KIND),intent(IN) :: dl,dv, rho, yhn, dtau_euv 
    real(kind=DBL_KIND),intent(INOUT) :: khpi,heat_hpi
    real(kind=DBL_KIND) :: lum,dF,dF_EUV,df_hpi,dl_cgs,phi,theta,alpha_euv,heat_euv, xnH, gfuv2, NcCO
    real(kind=DBL_KIND) :: Av, rcod, heat_dust 
    real(kind=DBL_KIND),dimension(MX:MZ) :: drdt ! direction of ray
    ! integer :: isrc,res
    integer :: res    
    integer(kind=LLONG_KIND) :: N_pix
    real(kind=DBL_KIND),parameter :: chi_H=13.6d0!Hの電離エネルギー (in eV)
    real(kind=DBL_KIND),parameter :: tau_max_euv=1d2 !tau(0) > tau_max_euv のときは電離光子の寄与はゼロとする

#ifdef EXTERNALFORCE
    real(kind=DBL_KIND),intent(INOUT) :: fx_hpi,fy_hpi,fz_hpi
#endif !EXTERNALFORCE   

#if MODEL_ART == 2  
    real(kind=DBL_KIND),intent(INOUT) :: kh2pd
    real(kind=DBL_KIND) :: fss,NcH2,dF_FUV, exp_uv
    real(kind=DBL_KIND),parameter :: pdis = 0.15
    real(kind=DBL_KIND),parameter :: alpha_fuv = pdis*3.4d-10/1.21d+7 !conversion factor from F_FUV -> kh2pd
  #ifdef METAL
    real(kind=DBL_KIND),intent(INOUT) :: kgfuv, kdph, krdco, rOII
    real(kind=DBL_KIND),intent(IN)    :: dtau_de, dtau_df
    real(kind=DBL_KIND) :: tau_deuv, tau_dfuv, ddph,test, test2, tap
    real(kind=DBL_KIND) :: exp_dfuv, exp_fuv
  #endif
#endif
    real(kind=DBL_KIND) :: tau_uv, dtau_uv, dtau_fuv, dkhpi, lum_euv, lum_deuv, lum_fuv, mass_cl

    xnH = rho*Unit_rhomu !in cm^-3 

    !ray info
    res = Hpix%res    
    N_pix = 12_8 * 4_8**res ! _8 means long integer in Fortran           
    lum = rs_info%lum(isrc)/N_pix             !luminosity of the ray
    dF = lum * (dl/dv) / Unit_l**2            !sum of dF is F averaged in a cell (attenuation within the cell is neglected)


    ! optical depth of dust grain ------------------------------------------------------
#ifdef METAL
        ! opcaity of dust grain
        tau_deuv = tau(BIN_DEUV) ! optical depth of dust grain for EUV
        tau_dfuv = tau(BIN_DFUV) ! optical depth of dust grain for FUV

  #ifndef NO_IONIZATION
    #ifndef SET_NODUST_ATTENUATION
        tau_uv  = tau(BIN_EUV) + tau_deuv            !total optical depth for euv
    #else  ! ifndef SET_NODUST_ATTENUATION
        tau_uv  = tau(BIN_EUV)                       !total optical depth for euv
    #endif
  #else  ! ifndef NO_IONIZATION
        tau_uv  = tau_deuv
  #endif ! ifndef NO_IONIZATION

#else !METAL
        tau_uv   = tau(BIN_EUV)                   !total optical depth for uv
        tau_dfuv = 0.d0
#endif!METAL

    if (tau_uv < 1.d2) then !tau > tau_max_euvのときは、すでに十分弱くなっているのでEUVの寄与はゼロとする
        exp_uv = dexp(-tau_uv)
    else
        exp_uv = 0.d0
    endif

    if (tau_dfuv < 1.d2) then !tau > tau_max_euvのときは、すでに十分弱くなっているのでEUVの寄与はゼロとする
        exp_dfuv = dexp(-tau_dfuv)
    else
        exp_dfuv = 0.d0
    endif

    !local optical depth ----------------------------------------
#ifdef METAL
  #ifndef NO_IONIZATION
    #ifndef SET_NODUST_ATTENUATION
      dtau_uv = dtau_euv + dtau_de
      dtau_fuv= dtau_df
    #else  ! ifndef SET_NODUST_ATTENUATION
      dtau_uv = dtau_euv
      dtau_fuv= 0.d0
    #endif
  #else  ! ifndef NO_IONIZATION
      dtau_uv = dtau_de
      dtau_fuv= dtau_df
  #endif ! ifndef NO_IONIZATION
#else !METAL
      dtau_uv = dtau_euv
      dtau_fuv= 0.d0
#endif !METAL
    ! -----------------------------------------------------------


    ! EUV photon number density flux
    dF_EUV = dF*rs_info%x_euv(isrc)*exp_uv   
    
    !radiation source info
    alpha_euv = rs_info%alpha_euv(isrc) ! cross section of hydrogen [cm^2]
    heat_euv  = rs_info%heat_euv(isrc)  ! heating rate per ionization event [erg]

    !store radiation effects on a cell
    dkhpi    = dF_EUV*alpha_euv  
    khpi     = khpi     + dkhpi                    !ionization rate per H particle [s^-1]
    heat_hpi = heat_hpi + dkhpi*heat_euv           !heating rate per H particle [erg s^-1]

    !OII ionization rate
#ifdef METAL
    rOII = rOII + dkhpi*rs_info%rOII(isrc)          !OII ionization rate
#endif

    ! radiation pressure ------------------------------------------------
    if(dtau_uv > 1.d-5) then
      lum_euv  = lum*rs_info%lumeuv(isrc)*exp_uv*(1.d0-dexp(-dtau_uv))  ! cellで吸収されるenergy [erg s^-1]
    else
      lum_euv  = lum*rs_info%lumeuv(isrc)*exp_uv*dtau_uv ! cellで吸収されるenergy [erg s^-1]
    endif

    if(dtau_fuv > 1.d-5 )then
      lum_fuv  = lum*rs_info%lumfuv(isrc)*exp_dfuv*(1.d0-dexp(-dtau_fuv))
    else
      lum_fuv  = lum*rs_info%lumfuv(isrc)*exp_dfuv*dtau_fuv
    endif
    ! -------------------------------------------------------------------


#ifdef EXTERNALFORCE
    ! luminosity 
    mass_cl  = dv*rho*Unit_m*cgs_c ! [g cm s^-1]
    df_hpi   = (lum_euv+lum_fuv)/mass_cl/Unit_acc ![cm s^-2] => [noD]
   !df_hpi_i = (dF_EUV*(alpha_euv*(heat_euv+chi_H*cgs_ev)) /cgs_c)*yhn/((1.d0 + 4.d0*yHe)*cgs_mh)  ! [ cm s^-2  ]

    fx_hpi = fx_hpi + df_hpi*ndir(MX)    ! [cm s^{-2}] => [noD]              
    fy_hpi = fy_hpi + df_hpi*ndir(MY)    ! [cm s^{-2}] => [noD]                
    fz_hpi = fz_hpi + df_hpi*ndir(MZ)    ! [cm s^{-2}] => [noD]                
#else
    df_hpi = 0.d0
#endif !EXTERNALFORCE
    ! --------------------------------------------------------------------

    ! heating rate of dust grain ---------------------------------------
#ifdef METAL
    if(dtau_de > 1.d-5) then
      lum_deuv  = lum*rs_info%lumeuv(isrc)*exp_uv*(1.d0-dexp(-dtau_de))  ! cellで吸収されるenergy [erg s^-1]
    else
      lum_deuv  = lum*rs_info%lumeuv(isrc)*exp_uv*dtau_de ! cellで吸収されるenergy [erg s^-1]
    endif

    heat_dust = (lum_deuv+lum_fuv) / (dv*Unit_l3) ! [erg s^-1 cm^-3]
    kdph = kdph + max(heat_dust, 0.d0) ! [erg s^-1 cm^-3] 
#endif
    ! -------------------------------------------------------------------
    

    !effect of FUV -------------------------------------------------------
#if MODEL_ART == 2

#ifdef METAL
    Av   = 5.34d-22/rs_info%sig_fuv(isrc)*tau(BIN_DFUV) 
    gfuv2   = dF*rs_info%lumfuv(isrc)/MP_FISRF  ! [noD] 

    if(Av > 40.d0) then
      exp_fuv = 0.d0
    else
      
    kgfuv   = kgfuv + gfuv2*dexp(-1.8 * Av) ! from Nakatani 2018

  #ifndef SET_NODUST_ATTENUATION
      exp_fuv = dexp(-2.5 * Av)
  #else
      exp_fuv = 1.d0
  #endif
    endif

    NcH2 = tau(BIN_H2)
    fss =  fss_DB98v2(NcH2)
    dF_FUV = dF * rs_info%x_fuv(isrc) * fss * exp_fuv     !sum of dF_FUV is F_FUV averaged in a cell (self-shielding within cell is neglected)

    kh2pd = kh2pd + dF_FUV*alpha_fuv                            ! photo-dissociation rate per H2 particle [s^-1] from Hosokawa-san's code

    NcCO = tau(BIN_CO)
    call COdissociation_rate(NcCO, NcH2, Av, gfuv2, rcod)
    krdco = krdco + rcod

#else ! METAL
    NcH2 = tau(BIN_H2)
    fss =  fss_DB98v2(NcH2)
    dF_FUV = dF * rs_info%x_fuv(isrc) * fss     !sum of dF_FUV is F_FUV averaged in a cell (self-shielding within cell is neglected)
    kh2pd = kh2pd + dF_FUV*alpha_fuv                            ! photo-dissociation rate per H2 particle [s^-1] from Hosokawa-san's code
#endif !METAL
                                                               ! (KS+2014の式(13)と対応、誤差10パーセント以下)
#endif ! MODEL_ART


  contains
    ! !self-shielding factorの表式 (Draine&Bertoldi 1998のsimple版)
    ! function fss_DB98v1(NcH2)
    !   real(kind=DBL_KIND) :: fss_DB98v1
    !   real(kind=DBL_KIND),intent(IN) :: NcH2
    !   real(kind=DBL_KIND),parameter :: NcH2_cr = 1d14
    !   fss_DB98v1 = min(1., (NcH2/NcH2_cr)**(-0.75))
    ! end function fss_DB98v1

    !self-shielding factorの表式 (Draine&Bertoldi 1998のsimple版 + Hosokawa+12のcap)
    function fss_DB98v2(NcH2)
      real(kind=DBL_KIND) :: fss_DB98v2
      real(kind=DBL_KIND),intent(IN) :: NcH2
      real(kind=DBL_KIND),parameter :: NcH2_cr = 1d14, NcH2_max = 1d22

      real(kind=DBL_KIND) :: x,y,z ! temp vars

      if (NcH2 > NcH2_max) then
         fss_DB98v2 = 0d0 ! NcH2 > 1d22 cm^-2のときに強制的に光解離率をゼロにする (Hosokawa+12)
      else if (NcH2 < NcH2_cr) then
         fss_DB98v2 = 1d0 ! NcH2 < 1d14 cm^-2のときはself-shielding効かない
      else
         !-0.75乗をちょっと速くなるように書き換え（2倍くらい変わるらしい）
         x = sqrt(NcH2_cr/NcH2)
         y = sqrt(x)
         fss_DB98v2 = x*y
         ! fss_DB98v2 = (NcH2/NcH2_cr)**(-0.75)         
      end if
    end function fss_DB98v2


  end subroutine art_StoreRadEffect

  !-------------------------------------------------------------------------
  ! split a ray into four sub-rays
  ! INPUT:
  !   ray
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_SplitRay(ray)
    type(t_ray),pointer,intent(INOUT) :: ray

    integer :: isrc
    type(t_dirHpix) :: Hpix, Hpix_new         
    real(kind=DBL_KIND) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1) :: tau 
    integer(kind=LLONG_KIND) :: ipix
  
    !ray情報を取り出して破棄
    isrc = ray%isrc
    Hpix = ray%Hpix
    dist = ray%dist
    tau = ray%tau
    deallocate(ray)

    ! new sub-rays
    !
    ! solid-angle refinement in the nested regime of HEALPix
    !  ipix @res -> 4*ipix, 4*ipix+1, 4*ipix+2, 4*ipix+3 @res+1
    !
    Hpix_new%res = Hpix%res+1

    !rayを4分割
    do ipix = 4_8*Hpix%ipix, 4_8*Hpix%ipix + 3_8 ! _8 means long integer in Fortran
       Hpix_new%ipix = ipix
       call art_StartRay(isrc, Hpix_new, dist, tau) 
    end do
  end subroutine art_SplitRay


  !-------------------------------------------------------------------------
  ! condition for ray splitting
  ! INPUT:
  !   ray
  !   h = cell width
  ! OUTPUT:
  !   bool = .TRUE. if condition is satisfied
  !-------------------------------------------------------------------------
  subroutine art_SplitCond(ray,h,bool)
    type(t_ray),pointer,intent(IN) :: ray
    real(kind=DBL_KIND), intent(IN) :: h
    logical,intent(INOUT) :: bool

    integer(kind=LLONG_KIND) :: N_pix

    !ray variables
    integer :: res
    real(kind=DBL_KIND) :: dist, ds

    integer :: count_output=0 ! warning出力の抑制用 (Fortranではデフォルトでsave属性?)    

    !ray情報の取り出し
    dist = ray%dist
    res = ray%Hpix%res
    ! N_pix = 12_8 * (2_8**res)**2 ! _8 means long integer in Fortran
    N_pix = 12_8 * 4_8**res ! _8 means long integer in Fortran    

    !1本のrayが担当する面積
    ds = dist**2 * 4*PI / N_pix
    !    print *, "h, dist, h**2/ds", h, dist, h**2/ds

    ! 1つのセル面あたりMIN_RAY_PER_CELL本以上のrayが割り当てられているか？
    if (h**2/ds < MIN_RAY_PER_CELL) then 
       bool = .TRUE. ! 足りなければsplit
    else
       bool = .FALSE. ! 足りてたらそのまま
    end if

    !resがHEALPIX_MAX_LEVELを超えないように要請
    if (bool .and. res+1 > HEALPIX_MAX_LEVEL) then
       bool = .FALSE.

       !Warningの出力
       if (count_output < 100) then
          !print '(A,I4,I4,I4,/,A)', " **WARNING** (SplitCond) Hpix res exceeding HEALPIX_MAX_LEVEL", &
          print *, " **WARNING** (SplitCond) Hpix res exceeding HEALPIX_MAX_LEVEL", &
               myrank, res+1, HEALPIX_MAX_LEVEL,"   -> NOT split ray"
       else if (count_output == 100) then
          !print '(A,I4)', "(SplitCond) warning is outputted 100 times. Suppress further waring...", myrank
          print *, "(SplitCond) warning is outputted 100 times. Suppress further waring...", myrank
       end if
       count_output = count_output+1              

       ! stop
    end if
    

  end subroutine art_SplitCond



  !-------------------------------------------------------------------------
  ! execute ray ending and record number of completed rays
  ! INPUT:
  !   ray
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_EndRay(ray)
    type(t_ray),pointer,intent(INOUT) :: ray

    !完了した仕事の記録
    call art_RecordComp(ray%Hpix%res)

    !rayをdeallocate
    deallocate(ray)
    
  end subroutine art_EndRay

  !-------------------------------------------------------------------------
  ! record number of completed rays
  ! INPUT:
  !   res
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_RecordComp(res)
    integer,intent(IN) :: res

    !resの値をチェック
    if (res > HEALPIX_MAX_LEVEL) then
       !print '(A,I4,I4,/,A)', "(RecordComp) res > HEALPIX_MAX_LEVEL", res, HEALPIX_MAX_LEVEL,&
       print *, "(RecordComp) res > HEALPIX_MAX_LEVEL", res, HEALPIX_MAX_LEVEL,&
            "something is strange. stopping..."
       stop
    end if

    !完了した仕事の記録
    N_comp_local(myrank) = N_comp_local(myrank) + (2_8**(HEALPIX_MAX_LEVEL-res))**2 ! 太いrayは貢献大
    N_comp_local_update = .TRUE.
  end subroutine art_RecordComp


  !-------------------------------------------------------------------------
  ! condition for the end of a ray
  ! INPUT:
  !   ray
  ! OUTPUT:
  !   bool = .TRUE. if condition is satisfied
  !-------------------------------------------------------------------------
  subroutine art_EndCond(isrc,Hpix,dist,tau,bool,ndir_in)    
    use modelParameter, only : MP_Boxsize 
    integer,intent(IN) :: isrc
    type(t_dirHpix),intent(IN) :: Hpix         
    real(kind=DBL_KIND),intent(IN) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1),intent(IN) :: tau
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN),optional :: ndir_in
    logical,intent(INOUT) :: bool
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    real(kind=DBL_KIND),dimension(MX:MZ) :: r0 
    real(kind=DBL_KIND) :: theta, phi
    real(kind=DBL_KIND),dimension(MX:MZ) :: ndir !direction vector

    real(kind=DBL_KIND) :: c_eV,hnu_ion,dNiondt
    
    real(kind=DBL_KIND),parameter :: NcH2_max = 1d22 ! maximum H2 column density for FUV
    
    !初期値はFalse
    bool = .FALSE.

    !HEALPix <-> 角度 の変換
    if (.not. present(ndir_in)) then
       call pix2ang_nest(2**(Hpix%res), Hpix%ipix, theta, phi)
       call art_ang2vec(theta,phi,ndir)
    else
       ndir = ndir_in
   end if
    
    !rayの先端
    pos(:) = rs_info%spos(:,isrc) + dist * ndir(:)

    ! FUV はthinとする場合は、tau_euvがthresholdを超えたら終わり
    ! FUV RTを考慮する場合は、基本的には全計算領域にrayを飛ばす必要あり

#if MODEL_ART == 1
    if(tau(BIN_EUV) > TAU_ENDRAY) then
        #ifdef METAL ! HFADDED
            if(tau(BIN_DEUV) > NcnH_max) bool = .TRUE. 
        #else
            bool = .TRUE.
        #endif
#ifdef ART_DEBUG
       print '(A,I4,3(A,1P1E10.2))', "(EndCond) ray end due to tau > TAU_ENDRAY: myrank=", &
            myrank,", dist=", dist,", tau=", tau(BIN_EUV),", TAU_ENDRAY=", TAU_ENDRAY
#endif
    end if
#elif MODEL_ART == 2

  #ifdef METAL
    #ifndef NO_IONIZATION
      #ifndef SET_NODUST_ATTENUATION
        if( (tau(BIN_EUV)+tau(BIN_DEUV)) > TAU_ENDRAY) then
            if(tau(BIN_DFUV) > TAU_ENDRAY .or. tau(BIN_H2) > NcH2_max) bool = .TRUE.
      #else
        if(tau(BIN_EUV) > TAU_ENDRAY .and. tau(BIN_H2) > NcH2_max) then
            bool = .TRUE.
      #endif
    #else

        if( (tau(BIN_DEUV)) > TAU_ENDRAY) then
            if(tau(BIN_DFUV) > TAU_ENDRAY .or. tau(BIN_H2) > NcH2_max) bool = .TRUE.
    #endif
  #else
        if(tau(BIN_EUV) > TAU_ENDRAY .and. tau(BIN_H2) > NcH2_max) then
            bool = .TRUE.
  #endif

  #ifdef ART_DEBUG
          print '(A,I4,5(A,1P1E10.2))', "(EndCond) ray end due to tau > TAU_ENDRAY: myrank=", &
            myrank,", dist=", dist,", tau_EUV=", tau(BIN_EUV),", TAU_ENDRAY=", TAU_ENDRAY, &
            ", NcH2=", tau(1),", NcH2_max=", NcH2_max
  #endif
        end if
#endif ! MODEL_ART

    !壁@x,y,z=\pm MP_Boxsize (in code unit) の付近にきたらrayを終わらせる
    if (maxval(abs(pos(:))) > MP_Boxsize*0.99 ) then
       bool = .TRUE.
#ifdef ART_DEBUG
       print '(A,I4,2(A,1P1E10.2))', "(EndCond) ray end due to box boundary: myrank=", &
            myrank,", dist=", dist,", tau=", tau(BIN_EUV)
#endif
    end if

  end subroutine art_EndCond

  !-------------------------------------------------------------------------
  ! issend MPI messages (rays and N_comp)
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_SendAll()

    call art_TestRecvd() ! mpi_testを実行し、msgのリストから受け取り済みのものを取り除く

    if (maxval(N_send) > 0) then             !送るrayがあれば
       call art_SendRays()                   !rayを送信
    else if (N_comp_local_update .and. (.not. associated(my_ray_list))) then !未完了、および、送るrayが無くて、かつ、N_comp_localに更新あれば
       call art_SendNcomp()                  !N_compを送信
    end if
  end subroutine art_SendAll


  !-------------------------------------------------------------------------
  ! test sent messhage has received
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_TestRecvd()
    integer :: req, msg
    logical :: flg

    type(t_msg_rays),pointer :: msg_rays     !今調べているmsg_rays
    type(t_msg_rays),pointer :: msg_rays0    !最後に調べた未着msg_rays
    type(t_msg_rays),pointer :: tmp_msg_rays !msg_raysの削除の際に一時的に利用
    type(t_msg_Ncomp),pointer :: msg_Ncomp  !今調べているmsg_Ncomp
    type(t_msg_Ncomp),pointer :: msg_Ncomp0 !最後に調べた未着msg_Ncomp
    type(t_msg_Ncomp),pointer :: tmp_msg_Ncomp !msg_Ncompの削除の際に一時的に利用

    !---------------------------- test msg_rays --------------------------------!

    !msg_raysについて、requestが届いているかを一つずつチェック
    msg_rays => msg_rays_list
    nullify(msg_rays0)             ! 初期化
    
    do while (associated(msg_rays))
       req = msg_rays%req
       call mpi_test(req, flg, MPI_STATUS_IGNORE, ierr)  ! メッセージが届いているかを確認

       if (flg) then !メッセージ受け取り済みなら       
#ifdef ART_DEBUG
          print '(A,I0,A,I0,A,I0,A,I0,A)',&
               "(TestRecved, rays) myrank=",myrank,", loop_count=",loop_count,", (req,tag) = (",msg_rays%req,",",msg_rays%tag,") has been received."
#endif

          tmp_msg_rays => msg_rays          !今回調べたmsg (deallocate用)
          msg_rays => msg_rays%next         !次に調べるmsg

          ! 今回調べたmsgを削除
          deallocate(tmp_msg_rays%bufs)          ! 付随の動的bufsを削除
          deallocate(tmp_msg_rays)                ! msg本体を削除

          ! リンクリストのつなぎ替え or 先頭の移動
          if (associated(msg_rays0)) then 
             msg_rays0%next => msg_rays  !つなぎ替え
          else
             msg_rays_list => msg_rays    !先頭の移動
          end if

       else  !メッセージが未着なら
          msg_rays0 => msg_rays           ! 最後に調べた未着msg
          msg_rays => msg_rays%next       ! 次に調べるmsg 
       end if
    end do

    !---------------------------- test msg_Ncomp --------------------------------!

    !msg_Ncompについて、requestが届いているかを一つずつチェック
    msg_Ncomp => msg_Ncomp_list
    nullify(msg_Ncomp0)             ! 初期化
    
    do while (associated(msg_Ncomp))
       req = msg_Ncomp%req
       call mpi_test(req, flg, MPI_STATUS_IGNORE, ierr)  ! メッセージが届いているかを確認

       if (flg) then !メッセージ受け取り済みなら       
#ifdef ART_DEBUG
          print '(A,I0,A,I0,A,I0,A,I0,A)',&
               "(TestRecved, Ncomp) myrank=",myrank,", loop_count=",loop_count,", (req,tag) = (",msg_Ncomp%req,",",msg_Ncomp%tag,") has been received."
#endif

          tmp_msg_Ncomp => msg_Ncomp          !今回調べたmsg (deallocate用)
          msg_Ncomp => msg_Ncomp%next         !次に調べるmsg

          ! 今回調べたmsgを削除
          deallocate(tmp_msg_Ncomp)                ! msg本体を削除 (Ncompの場合はraysの場合と異なり、bufsは静的に確保)

          ! リンクリストのつなぎ替え or 先頭の移動
          if (associated(msg_Ncomp0)) then 
             msg_Ncomp0%next => msg_Ncomp  !つなぎ替え
          else
             msg_Ncomp_list => msg_Ncomp    !先頭の移動
          end if

       else  !メッセージが未着なら
          msg_Ncomp0 => msg_Ncomp           ! 最後に調べた未着msg
          msg_Ncomp => msg_Ncomp%next       ! 次に調べるmsg 
       end if
    end do


  end subroutine art_TestRecvd

  !-------------------------------------------------------------------------
  ! pack and send rays in send_ray_list (each rank will send N_send rays) 
  !-------------------------------------------------------------------------
  subroutine art_SendRays()
    integer :: rank, position
    integer :: n
    character(100) :: s

    !ray info
    type(t_ray),pointer :: ray
    integer,dimension(0:N_SEND_MAX-1) :: isrc
    integer,dimension(0:N_SEND_MAX-1) :: res
    integer(kind=LLONG_KIND),dimension(0:N_SEND_MAX-1) :: ipix
    real(kind=DBL_KIND),dimension(0:N_SEND_MAX-1) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1, 0:N_SEND_MAX-1) :: tau
    integer,dimension(0:N_SEND_MAX-1) :: gid
    integer :: msg, tag, req
    logical :: tag_not_found
    integer :: n_packs, i_pack, N_send_each
    type(t_msg_rays),pointer :: msg_rays
    integer :: bufs_size_byte, bufs_size_int

#ifdef ART_DEBUG
    write (s,'(A,I0,A)') "(A,I4,A,",NPE,"(I4))"
    print s, "(SendRays: N_send) myrank=",myrank, ", N_send=",N_send
#endif

    do rank = 0, NPE-1
       if (rank == myrank) cycle
       if (N_send(rank) <= 0) cycle

       n_packs = N_send(rank)/N_SEND_MAX + 1 ! number of sent packages

#ifdef ART_DEBUG
       if (N_send(rank) > N_SEND_MAX) then       
          print '(A,I8,A,/,4(A,I8))', "(SendRays) N_send(rank) > N_SEND_MAX.            send ", &
               N_send(rank)/N_SEND_MAX, " packages.",&
               "myrank=",myrank, ", rank=",rank, ", N_send(rank)=",N_send(rank), ", N_SEND_MAX=",N_SEND_MAX 
       end if
#endif

       do i_pack = 0, n_packs-1
          if (i_pack == n_packs-1) then 
             N_send_each = N_send(rank) - (n_packs-1)*N_SEND_MAX !最後は残りのray全部
          else
             N_send_each = N_SEND_MAX
          end if

          ! get new tag
          !--------- KS MODIFIED ----------!
          ! tag_max_ray = tag_max_ray + 2
          ! tag = tag_max_ray
          call get_new_tag_rays(tag)
          !--------- KS MODIFIED ----------!

          !ray情報の取り出し
          do n = 0, N_send_each-1
             ray => send_ray_list(rank)%ray_list
             send_ray_list(rank)%ray_list => send_ray_list(rank)%ray_list%next
             isrc(n) = ray%isrc
             res(n) = ray%Hpix%res
             ipix(n) = ray%Hpix%ipix
             dist(n) = ray%dist
             tau(:,n) = ray%tau(:)
             gid(n) = ray%gid
             deallocate(ray)
#ifdef ART_DEBUG
             write (s,'(A,I0,A)') '(A,I4,A,I4,A,I8,A,I8,I8,I18,(1P',NUM_FREQ+1,'E15.5),I8)'
             print s, "(SendRays: ray) myrank=", myrank,&
                  ", rank=",rank,", tag=",tag, ",    ray_info = ", isrc(n),res(n),ipix(n),dist(n), tau(:,n),  gid(n)             
#endif
          end do


          ! allocate msg_rays
          allocate(msg_rays)
          bufs_size_int = (N_send_each*RAY_BYTE+NSEND_BYTE+3)/4 ! size in integer
          allocate(msg_rays%bufs(0:bufs_size_int-1))
          bufs_size_byte = sizeof(msg_rays%bufs)           ! size in byte

          !start packing data
          position = 0
          !pack N_send
          call mpi_pack(N_send_each, 1, MPI_INTEGER, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)
          !pack isrc
          call mpi_pack(isrc(0), N_send_each, MPI_INTEGER, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)
          !pack res
          call mpi_pack(res(0), N_send_each, MPI_INTEGER, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)
          !pack ipix
          call mpi_pack(ipix(0), N_send_each, MPI_INTEGER8, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)          
          !pack d
          call mpi_pack(dist(0), N_send_each, MPI_DOUBLE, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)
          !pack tau
          call mpi_pack(tau(0,0), NUM_FREQ*N_send_each, MPI_DOUBLE, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)
          !pack gid
          call mpi_pack(gid(0), N_send_each, MPI_INTEGER, msg_rays%bufs(0), bufs_size_byte, position, MPI_COMM_WORLD, ierr)                    

          !send data
          call mpi_issend(msg_rays%bufs(0), position, MPI_PACKED, rank, tag, MPI_COMM_WORLD, req, ierr)  ! MPI_PACKED -> 2番目引数はByte数


          ! tagとreqの情報をmsg_raysに格納
          msg_rays%tag = tag
          msg_rays%req = req

          !msg_rays_listの一番上にmsg_raysを格納
          msg_rays%next => msg_rays_list
          msg_rays_list => msg_rays
          
#ifdef ART_DEBUG
          print '(A,I0,A,I0,A,I0,A,I0,A)', "(SendRays: msg) myrank=", myrank,&
               ", rank=",rank,",        message with (req,tag) = (",req,",",tag, ")   has been sent"          
#endif
       end do !i_pack
    end do !rank

    ! reset N_send
    N_send(:) = 0

  contains

    !--------------------------------------------------------
    ! get_new_tag_rays
    !--------------------------------------------------------
    !
    ! 既存のtagと重複しない一番小さい奇数の値を取得
    !    
    !--------------------------------------------------------
    subroutine get_new_tag_rays(tag)
      integer, intent(OUT) :: tag
      integer :: i, tag_cand
      type(t_msg_rays),pointer :: msg_rays
      logical :: flg

      ! tag = -1
      ! do i = 0, MPI_TAG_MAX_MACHINE/2 -1
      !    tag_cand = 2*i+1 !tagの値の候補
      !    flg = .True.     !既存のtagと重複がないときtrue

      !    !msg_rays_listについて既存tagと重複しないかをチェック
      !    msg_rays => msg_rays_list    
      !    do while (associated(msg_rays))     !全てのmsg_raysについてのループ
      !       if (tag_cand == msg_rays%tag) then
      !          flg = .False. !重複あり
      !          exit
      !       end if
      !       msg_rays => msg_rays%next         ! 次へ
      !    end do
         
      !    !重複無かった場合は
      !    if (flg) then
      !       tag = tag_cand !候補の値を代入
      !       exit           !tagの値を見つけたのでループから抜ける
      !    end if
      ! end do

      ! if (tag < 0) then
      !    print '(A,/,A)',  "(get_new_tag_rays) tag cannot be found", "stopping..."
      !    stop
      ! end if
      
      !---------------------------------------------------!
      !   tag の値をいちいち振る必要がない気もする (KS TODO)   !
      !---------------------------------------------------!
      tag = 1 !奇数であればよい？
      !---------------------------------------------------!
      !---------------------------------------------------!
      
    end subroutine get_new_tag_rays
            
  end subroutine art_SendRays

  !-------------------------------------------------------------------------
  ! send N_comp_local
  !-------------------------------------------------------------------------
  subroutine art_SendNcomp()
    integer :: rank, tag, n, msg, req
    logical :: tag_not_found
    character(100) :: s
    integer :: position
    type(t_msg_Ncomp),pointer :: msg_Ncomp    

    do rank = 0, NPE-1
       if (rank == myrank) cycle

       ! get new tag
       !--------- KS MODIFIED ----------!
       ! tag_max_Ncomp = tag_max_Ncomp + 2
       ! tag = tag_max_Ncomp
       call get_new_tag_Ncomp(tag)
       !--------- KS MODIFIED ----------!       

       ! allocate msg_Ncomp
       allocate(msg_Ncomp)

       ! store data
       msg_Ncomp%bufs = N_comp_local(myrank)

       !send data (普通にMPI_INTEGER8で送る)
       call mpi_issend(msg_Ncomp%bufs, 1, MPI_INTEGER8, rank, tag, MPI_COMM_WORLD, req, ierr)

       ! store tag & req
       msg_Ncomp%tag = tag
       msg_Ncomp%req = req

       !msg_Ncomp_listの一番上にmsg_Ncompを格納
       msg_Ncomp%next => msg_Ncomp_list
       msg_Ncomp_list => msg_Ncomp
       
       
#ifdef ART_DEBUG
       print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', "(SendNcomp: msg) myrank=", myrank,&
            ", rank=",rank,", N_comp_local(myrank)=",N_comp_local(myrank),",    message with (req,tag) = (",req,",",tag, ")   has been sent"       
#endif

    end do

    N_comp_local_update = .FALSE.       ! 送信できたらflagをリセット

  contains

    !--------------------------------------------------------
    ! get_new_tag_Ncomp
    !--------------------------------------------------------
    !
    ! 既存のtagと重複しない一番小さい偶数の値を取得
    !    
    !--------------------------------------------------------
    subroutine get_new_tag_Ncomp(tag)
      integer, intent(OUT) :: tag
      integer :: i, tag_cand
      type(t_msg_Ncomp),pointer :: msg_Ncomp
      logical :: flg

      ! tag = -1
      ! do i = 0, MPI_TAG_MAX_MACHINE/2
      !    tag_cand = 2*i !tagの値の候補
      !    flg = .True.     !既存のtagと重複がないときtrue

      !    !msg_Ncomp_listについて既存tagと重複しないかをチェック
      !    msg_Ncomp => msg_Ncomp_list    
      !    do while (associated(msg_Ncomp))     !全てのmsg_Ncompについてのループ
      !       if (tag_cand == msg_Ncomp%tag) then
      !          flg = .False. !重複あり
      !          exit
      !       end if
      !       msg_Ncomp => msg_Ncomp%next         ! 次へ
      !    end do
         
      !    !重複無かった場合は
      !    if (flg) then
      !       tag = tag_cand !候補の値を代入
      !       exit           !tagの値を見つけたのでループから抜ける
      !    end if
      ! end do

      ! if (tag < 0) then
      !    print '(A,/,A)',  "(get_new_tag_Ncomp) tag cannot be found", "stopping..."
      !    stop
      ! end if
      
      !---------------------------------------------------!
      !   tag の値をいちいち振る必要がない気もする (KS TODO)   !
      !---------------------------------------------------!
      tag = 0 !偶数であればよい？
      !---------------------------------------------------!
      !---------------------------------------------------!
      
    end subroutine get_new_tag_Ncomp
            
  end subroutine art_SendNcomp


  !-------------------------------------------------------------------------
  ! probe and recv MPI messages (rays and N_comp)
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_RecvAll()
    logical :: flg
    integer, dimension(0:NPE-1) :: N_onetime !一度に受け取るメッセージの数 (rank毎)
    integer, parameter :: N_MAX_ONETIME = 64 !一度に受け取れるメッセージの数の最大値 (マシン依存あるかも、超えるとエラー吐かずにデータ壊れるので注意)
    !integer :: status(MPI_STATUS_SIZE)  ! 変数statusはmpilib.F90で定義済み

    N_onetime = 0
    !送信中データを全部受信
    do
       call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flg, status, ierr);   !送信中データのチェック
       if (.not. flg) exit ! データ無ければループを抜ける

#ifdef ART_DEBUG
       print '(A,I0,A,I0,A,I0,A,I0,A)',&
            "(RecvAll: iprobe) myrank=",myrank,", rank=",status(MPI_SOURCE),", loop_count=",loop_count,&
            "            message with tag=", status(MPI_TAG),"  is going to be received."
#endif
       if (mod(status(MPI_TAG),2) == 1) then           !奇数のtagはray用
          call art_RecvRays(status(MPI_SOURCE), status(MPI_TAG))     !rayを受け取る
       else                                            !偶数のtagはN_comp用
          call art_RecvNcomp(status(MPI_SOURCE), status(MPI_TAG))    !N_compを受け取る
       end if
    end do
  end subroutine art_RecvAll

  !-------------------------------------------------------------------------
  ! recv, unpack and store rays
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_RecvRays(rank, tag)
    integer, intent(IN) :: rank, tag
    integer :: n, position !, rank

    !ray info
    type(t_ray),pointer :: ray
    integer,dimension(0:N_SEND_MAX-1) :: isrc
    integer,dimension(0:N_SEND_MAX-1) :: res
    integer(kind=LLONG_KIND),dimension(0:N_SEND_MAX-1) :: ipix
    real(kind=DBL_KIND),dimension(0:N_SEND_MAX-1) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1, 0:N_SEND_MAX-1) :: tau
    integer,dimension(0:N_SEND_MAX-1) :: gid
    type(t_dirHpix) :: Hpix
    integer :: N_recv
    character(100) :: s
    integer :: bufr_size_byte
    integer,save,dimension(0:BUFSIZE_RAYS-1) :: bufr !受信用buffer

    !size of bufr
    bufr_size_byte = sizeof(bufr)

    !データが送られてきていれば受信 
    call mpi_recv(bufr(0), bufr_size_byte, MPI_PACKED, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! MPI_PACKED -> 2番目引数はByte数

    !start unpacking data
    position = 0
    !unpack N_recv
    call mpi_unpack(bufr(0), bufr_size_byte, position, N_recv, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !check the size of bufr
    if ((N_recv*RAY_BYTE+NSEND_BYTE+3)/4 > BUFSIZE_RAYS) then !切り上げ注意
       !print '(A,/,4(A,I8),/,A)', "(RecvRays: error) size of bufr is not enough. larger BUFSIZE_RAYS may be needed.",&
       print *, "(RecvRays: error) size of bufr is not enough. larger BUFSIZE_RAYS may be needed.",&
            "myrank= ", myrank, ", N_recv= ", N_recv, ", send_size=",(N_recv*RAY_BYTE+NSEND_BYTE+3)/4, &
            ", BUFSIZE_RAYS=",BUFSIZE_RAYS,  "stopping..."
       stop
    end if

    !unpack isrc
    call mpi_unpack(bufr(0), bufr_size_byte, position, isrc(0), N_recv, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    !unpack res
    call mpi_unpack(bufr(0), bufr_size_byte, position, res(0),N_recv, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    !pack ipix
    call mpi_unpack(bufr(0), bufr_size_byte, position, ipix(0),N_recv, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
    !unpack d
    call mpi_unpack(bufr(0), bufr_size_byte, position, dist(0),N_recv, MPI_DOUBLE, MPI_COMM_WORLD, ierr) 
    !unpack tau
    call mpi_unpack(bufr(0), bufr_size_byte, position, tau(0,0),NUM_FREQ*N_recv, MPI_DOUBLE, MPI_COMM_WORLD, ierr) 
    !unpack gid
    call mpi_unpack(bufr(0), bufr_size_byte, position, gid(0),N_recv, MPI_INTEGER, MPI_COMM_WORLD, ierr)


    if (N_recv > N_SEND_MAX) then
       print *, "N_recv > N_SEND_MAX.    something strange happens.   stopping..."
       stop
    end if

    !受け取ったray情報をmy_ray_listに格納
    do n = 0,N_recv-1
       Hpix%res = res(n)
       Hpix%ipix = ipix(n)

       !rayの作成
       call art_CreateRay(isrc(n), Hpix, dist(n), tau(:,n), gid(n), ray)
       !rayをmy_ray_listに追加
       call art_AddMyRay(ray)

#ifdef ART_DEBUG
       write (s,'(A,I0,A)') '(A,I4,A,I4,A,I8,A,I8,I8,I18,(1P',NUM_FREQ+1,'E15.5),I8)'
       print s, "(RecvRays: ray) myrank=", myrank,&
                  ", rank=",rank,", tag=",tag, ",    ray_info = ", isrc(n),res(n),ipix(n),dist(n), tau(:,n),  gid(n)
#endif       
    end do

  end subroutine art_RecvRays

  !-------------------------------------------------------------------------
  ! recv N_comp_local
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_RecvNcomp(rank, tag)
    integer, intent(IN) :: rank, tag
    integer :: position


    !データが送られてきていれば受信
    call mpi_recv(bufr_Ncomp, 1, MPI_INTEGER8, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! MPI_INTEGER8 -> 2番目引数は個数

    !N_comp_local(rank)に格納
    N_comp_local(rank)  = bufr_Ncomp


#ifdef ART_DEBUG
       print '(A,I4,A,I4,A,I4,A,I18)', "(RecvNcomp: Ncomp) myrank=", myrank,&
            ", rank=",rank,", tag=",tag, ",    N_comp_local(rank) = ", N_comp_local(rank)
#endif       

  end subroutine art_RecvNcomp


  !-------------------------------------------------------------------------
  ! check whether all ray tracing has been finished
  ! INPUT:
  ! OUTPUT:
  !-------------------------------------------------------------------------
  subroutine art_CheckAllDone(bool)
    logical,intent(OUT) :: bool
    character(100) :: s

    N_comp_global = SUM(N_comp_local)

    ! 仕事が完了したかをチェック
    if (N_comp_global > N_comp_total) then !完了した仕事 > 全仕事 の場合何かがおかしいのでstop
       !write (s,'(A,I0,A)') "(A,I4,A,",NPE,"(I18),A,I18,A,I18,/,A,/,A)"
       print *, "(CheckAllDone) myrank=", myrank,&
            ", N_comp_local=",N_comp_local,", N_comp_global=",N_comp_global,", N_comp_total=",N_comp_total,&
            "N_comp_global > N_comp_total.        something is wrong.","stopping..."
       stop
    else if (N_comp_global == N_comp_total) then !全仕事を完了したらbool = True
#ifdef ART_DEBUG
       write (s,'(A,I0,A)') "(A,I4,A,",NPE,"(I18),A,I18,A,I18)"
       print s, "(CheckAllDone) myrank=", myrank,&
            ", N_comp_local=",N_comp_local,", N_comp_global=",N_comp_global,", N_comp_total=",N_comp_total
#endif
       bool = .TRUE.
    else
       bool = .FALSE.       
    end if
  end subroutine art_CheckAllDone  
  
  !-------------------------------------------------------------------------
  ! converting ang (theta, phi) to vec (ndir)
  !-------------------------------------------------------------------------
  subroutine art_ang2vec(theta,phi,ndir)
    real(kind=DBL_KIND),intent(IN) :: theta, phi
    real(kind=DBL_KIND),dimension(MX:MZ),intent(OUT) :: ndir !direction vector
    ndir(MX)=sin(theta)*cos(phi)
    ndir(MY)=sin(theta)*sin(phi)
    ndir(MZ)=cos(theta)
#ifndef ART_NO_RANDOM_ROTAION
    ndir(:) = matmul(m_rot(:,:),ndir(:))
#endif ! ART_NO_RANDOM_ROTAION        
  end subroutine art_ang2vec

  !-------------------------------------------------------------------------
  ! getting random rotation matrix
  ! reshapeと(//)を使った値は、A(0,0), A(1,0), ..., A(1,2), A(2,2)の順
  !-------------------------------------------------------------------------
  subroutine art_random_rot_matrix()
    real(kind=DBL_KIND),dimension(MX:MZ,MX:MZ) :: m_rot_phi, m_rot_theta !z軸周りの回転とz軸を傾ける回転
    real(kind=DBL_KIND),dimension(MX:MZ) :: n_newz
    integer,parameter :: n = 5
    real(kind=DBL_KIND) :: rand1,rand2,rand3, phi_z, phi, costh, sinth, cosph, sinph
    integer i, seedsize
    integer,allocatable:: seed(:)
    integer, save :: i_seed=0

    call random_seed(size=seedsize)  ! シードの格納に必要なサイズを取得する
    allocate(seed(seedsize))         ! シード格納領域を確保
    seed(:) = i_seed                 ! シードに値を代入
    call random_seed(put=seed)         ! シードを渡す

  call random_number(rand1)
  call random_number(rand2)
  call random_number(rand3)


  ! matrix generating rotaion around z-axis
  phi_z = 2.d0*Pi*rand1
  m_rot_phi(:,:) = RESHAPE((/&
       cos(phi_z),  sin(phi_z), 0d0, &
       -sin(phi_z), cos(phi_z), 0d0, &
       0d0,       0d0,      1d0  &
       /),(/3,3/))

  ! new z-axis (homogeneous in solid angle 4pi)
  costh = 2.d0*rand3 - 1.d0 
  sinth = sqrt(1d0-costh**2)
  phi = 2.d0*Pi*rand2
  cosph = cos(phi)  
  sinph = sin(phi)

  ! matrix generating the rotaion of z-axis (using Rodrigues' formula w/ n1 = sinph, n2 = -cosph, n3=0)
  m_rot_theta(:,:) = RESHAPE((/ &
       costh+sinph**2*(1d0-costh), -cosph*sinph*(1d0-costh), cosph*sinth, &
       -cosph*sinph*(1d0-costh), costh+cosph**2*(1d0-costh), sinph*sinth, &
       -cosph*sinth,             -sinph*sinth,               costh        &
       /), (/3,3/))

  
  !matrix generating random rotaion
  m_rot = matmul(m_rot_phi(:,:),m_rot_theta(:,:))

  !KS DEBUG
#ifdef ART_DEBUG
  if (get_myrank() == PRIMARY_RANK)  then
     print *, phi_z, phi, costh
     print '(A, (1P10E15.7))', "m_rot_phi", m_rot_phi(:,:)
     print '(A, (1P10E15.7))', "m_rot_theta", m_rot_theta(:,:)
     print '(A, (1P10E15.7))', "m_rot", m_rot(:,:)
  endif
#endif  

  i_seed = i_seed + 1 ! シードを更新

  end subroutine art_random_rot_matrix

  !-------------------------------------------------------------------------
  ! write current position of ray to a file for debug
  !-------------------------------------------------------------------------
  subroutine art_WriteTrajectory(ray)
    use string, only : CHARLEN, concat
    use io_util, only : read_env
    type(t_ray),intent(IN) :: ray
    integer :: isrc
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos
    real(kind=DBL_KIND),dimension(MTH:MPH) :: dir
    real(kind=DBL_KIND),dimension(MX:MZ) :: r0 
    type(t_dirHpix) :: Hpix         
    real(kind=DBL_KIND) :: dist
    real(kind=DBL_KIND),dimension(0:NUM_FREQ-1) :: tau 
    real(kind=DBL_KIND) :: theta, phi
    real(kind=DBL_KIND),dimension(MX:MZ) :: ndir !direction vector
    integer :: gid
    integer,parameter :: FH = 11
    character(len=CHARLEN) :: fn, direc

    integer,parameter :: ray_base = 0, ray_mod = 0
    integer :: res
    integer(kind=LLONG_KIND) :: ipix
    
    call read_env('DIR', direc)
    fn = concat(direc,'ray_traj.txt')

    !ray情報の取り出し
    isrc = ray%isrc
    Hpix = ray%Hpix
    dist = ray%dist
    tau = ray%tau
    gid = ray%gid

    !注目している一連のray以外は書きこまない
    if(get_ipix(Hpix%res,ray_base,ray_mod) /= Hpix%ipix) return

    !HEALPix <-> 角度 の変換
    call pix2ang_nest(2**(Hpix%res), Hpix%ipix, theta, phi)
    call  art_ang2vec(theta,phi,ndir)
    pos(:) = rs_info%spos(:,isrc) + dist * ndir(:)


    open(FH, file=fn,access='append')
    write(FH,*) Hpix%res,Hpix%ipix,dist,pos, theta, phi
    close(FH)
  contains
    !連なるrayのもつipix
    function get_ipix(res,ray_base,ray_mod)
      integer,intent(IN) :: res,ray_base,ray_mod
      integer(kind=LLONG_KIND) :: get_ipix
      integer :: i
      get_ipix = ray_base
      if (res == 0) return
      do i = 1, res
         get_ipix = 4 * get_ipix + ray_mod
      end do
    end function get_ipix
  end subroutine art_WriteTrajectory


#endif !RADTR_DIRECT
end module art
