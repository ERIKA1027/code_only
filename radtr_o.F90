#include "config.h"

#define DEBUG_RADTR   YES
#define NO_CHEMISTRY  NO
#define SET_RADSOURCE NO
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

  public :: radtr_moment 

contains

  !---------------------------------------------------------------
  ! M1 closerに関連する物理量等を初期化
  !---------------------------------------------------------------

  subroutine radtr_init(dt_hyd)
    
    logical, save :: bool_radinit = .false.
    real(kind=DBL_KIND),intent(IN) :: dt_hyd
#ifdef SOLVER_RADTR_HLL
    integer,parameter :: FH = 23
    integer :: i, j, ii, kk, err
    character(100) :: path2hll = "./rt/hll_evals.list"
    real(kind=DBL_KIND) :: dummy
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
        read(FH, *) i
        do i = 0, 100
          do j = 0, 100
            read(FH,*,iostat=err)ii,kk,lambda1(ii,kk),dummy,dummy,lambda4(ii,kk)
          enddo
        enddo
      endif

      call mpi_bcast(lambda1, size(lambda1), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
      call mpi_bcast(lambda4, size(lambda4), MPI_DOUBLE_PRECISION, PRIMARY_RANK, MPI_COMM_WORLD, ierr)
#endif

    endif

    tsum     = 0.d0 ! 輻射輸送における時間推進の合計値
    num_step = 0
    tend     = dt_hyd
    loop_end = .false.


  end subroutine radtr_init


  !---------------------------------------------------------------
  ! 流体の時間ステップまで輻射輸送を時間推進
  !---------------------------------------------------------------

  subroutine radtr_moment(dt_hyd)
    use mpilib
    real(kind=DBL_KIND),intent(IN) :: dt_hyd ! time step lengh of hydro


    ! ----------------------
    ! 輻射輸送について初期化
    ! ----------------------
    call radtr_init(dt_hyd)    


    ! ---------------------
    ! 輻射について時間推進
    ! ---------------------
    do   

      ! cfl condition for photons
      call cfl_condition_radrt
#if DEBUG_RADTR == YES
      if(get_myrank() == PRIMARY_RANK) write(*,'(a,I5,1p4e13.5)') , "step, tsum, dt, dt_hyd", num_step, dt_hyd/dt, tsum, dt, dt_hyd
#endif


      !-------------------
      ! injection step
      !-------------------
#if SET_RADSOURCE == YES
      call injection_step_test
#else
      call injection_step      
#endif


      !------------------
      ! transport step
      !------------------
      call timestep_radtr

      !------------------
      ! chemistry step
      !------------------
#if NO_CHEMISTRY == NO
      call ch_radtrchem(dt)
#endif
      call rescueLev_radtr ! 救済

      !-------------------
      ! converge all level
      !-------------------
      call converge_alllevel



      ! --------------
      ! loop end
      ! --------------
      if (loop_end) exit


    enddo


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

    integer :: level, n, gid

    !info for radiation source
    type(t_rs_info),pointer :: rs_info
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(MX:MZ) :: ijkg
    real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
    real(kind=DBL_KIND),dimension(MX:MZ) :: pos, h
    real(kind=DBL_KIND) :: Ndot_ion, del_V
    real(kind=DBL_KIND), dimension(:,:,:), pointer :: Nrad

    !----------------------------
    ! get information of RS
    !----------------------------
    call rs_GetSourceInfo(rs_info)


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
        h = CellWidth(:,level)

        Nrad => get_Ucomp(MER, gid)

        k = int((pos(MZ)-z(Kmin))/h(MZ) + 0.5d0)+Kmin
        j = int((pos(MY)-y(Jmin))/h(MY) + 0.5d0)+Jmin
        i = int((pos(MX)-x(Imin))/h(MX) + 0.5d0)+Imin

        ! volume of cell
        del_V = h(MX)*h(MY)*h(MZ)
        Ndot_ion = rs_info%x_euv(isrc)*rs_info%lum(isrc)*Unit_t/MP_PHON ! s^{-1} => noD

        Nrad(i,j,k) = Nrad(i,j,k) + Ndot_ion/del_V*dt

      enddo

    enddo

  end subroutine injection_step

  !------------------------------------------------------------------------------------------
  ! injection step test 
  !------------------------------------------------------------------------------------------

  subroutine injection_step_test

    use radiationSource
    use overBlockCoordinates

    integer :: level, n, gid

    !info for radiation source
    type(t_rs_info),pointer :: rs_info
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
      call get_flux_radtr
      call u_update_radtr
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
      call boundary_u( Ulist(n), STEP_MODE)
    enddo
  end subroutine boundary_cond


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
    integer :: n, gid, i, j, k
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
    real(kind=DBL_KIND) :: f2, reduc_f2
    real(kind=DBL_KIND),parameter :: Nrad_floor = 1.d-20 ! 100pc cubic に一つ=1.e-8
    logical :: isNotFinite

    do n = Gidmin, GidListMax( CurrentLevel )

      gid = GidList(n, CurrentLevel) ! gid for U

      Nrad => get_Ucomp(MER,gid)
      Frx  => get_Ucomp(MFRX,gid)
      Fry  => get_Ucomp(MFRY,gid)
      Frz  => get_Ucomp(MFRZ,gid)

      
      do k = Kmingh, Kmaxgh
        do j = Jmingh, Jmaxgh
          do i = Imingh, Imaxgh


            if(isNotFinite(Nrad(i,j,k))) then
              Nrad(i,j,k) = Nrad_floor
            endif

            if(isNotFinite(Frx(i,j,k))) then
              Frx(i,j,k) = 0.d0
            endif
            if(isNotFinite(Fry(i,j,k))) then
              Fry(i,j,k) = 0.d0
            endif
            if(isNotFinite(Frz(i,j,k))) then
              Frz(i,j,k) = 0.d0
            endif

            if(Nrad(i,j,k) < 0.0) then ! ゴミ排除
                !print '(A, 1P4E15.7)', "photon density becomes negative", Nrad(i,j,k),Frx(i,j,k),Fry(i,j,k),Frz(i,j,k)
                !stop
                Nrad(i,j,k) = 0.d0
                Frx(i,j,k)  = 0.d0
                Fry(i,j,k)  = 0.d0
                Frz(i,j,k)  = 0.d0
            endif
#if DEBUG_RADTR == YES
            !if( Nrad(i,j,k) > Nrad_floor*1.d3) then
            !    print *, i, j, k, Nrad(i,j,k), Frx(i,j,k), Fry(i,j,k), Frz(i,j,k)
            !endif
#endif

            !f2 = (Frx(i,j,k)*Frx(i,j,k)+Fry(i,j,k)*Fry(i,j,k)+Frz(i,j,k)*Frz(i,j,k)) &
            !  /(cl_til*Nrad(i,j,k))**2.d0

!            if(f2 > 4.d0/3.d0) then
!                reduc_f2   = dsqrt(4.d0/3.d0/f2)
!                Frx(i,j,k) = Frx(i,j,k)*reduc_f2
!                Fry(i,j,k) = Fry(i,j,k)*reduc_f2
!                Frz(i,j,k) = Frz(i,j,k)*reduc_f2
!
!                !test
!#if DEBUG_RADTR == YES
!                print *, "f2>4.d0/3.d0 at", i, j, k, f2, Nrad(i,j,k), Frx(i,j,k), Fry(i,j,k), Frz(i,j,k)
!#endif
!            endif



          enddo
        enddo
      enddo
    enddo
  end subroutine rescueLev_radtr

  !-------------------------
  ! get flux 
  !-------------------------
  subroutine get_flux_radtr

    integer :: i,j,k,n
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    real(kind=DBL_KIND), dimension(MX:MZ, MX:MZ, ARRAYSIZE_IJKGH) :: PP ! P_ij
    real(kind=DBL_KIND), dimension(MX:MZ) :: fi
    real(kind=DBL_KIND) :: inv_cn, f2, f2_sq, ff, chi, fchi1, fchi2, c2Nr
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
#ifdef SOLVER_RADTR_HLL
    real(kind=DBL_KIND), dimension(MX:MZ, ARRAYSIZE_IJKGH) :: lam_max, lam_min
    real(kind=DBL_KIND), dimension(MX:MZ) :: costheta, lam_p, lam_m, lcr_mp, lsa_mp
#endif
    real(kind=DBL_KIND) :: tap
    


    u => get_up( Ulist(CurrentIndex) )

    h   = CellWidth( :, CurrentLevel)

    !------------------------------
    ! set eddington tensor
    !------------------------------
    do k = Kmin-1, Kmax+1
      do j = Jmin-1, Jmax+1
        do i = Imin-1, Imax+1

          if(u(i,j,k,MER) .le. 0.d0) then
              PP(:, :, i, j, k) = 0.d0
#ifdef SOLVER_RADTR_HLL
              ff = 0.d0

              fi(MX) = u(i,j,k,MFRX)
              fi(MY) = u(i,j,k,MFRY)
              fi(MZ) = u(i,j,k,MFRZ)
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
            fi(MX) = u(i,j,k,MFRX)
            fi(MY) = u(i,j,k,MFRY)
            fi(MZ) = u(i,j,k,MFRZ)

            f2 = fi(MX)**2.d0+fi(MY)**2.d0+fi(MZ)**2.d0 
            
            f2_sq = dsqrt(f2)
            ff    = f2_sq/(cl_til*u(i,j,k,MER))

            !chi
            chi = dmax1(4.d0-3.d0*ff, 0.d0)
            chi = (3.d0+4.d0*ff)/(5.d0+2.d0*dsqrt(chi))
            chi = dmax1(dmin1(1.d0, chi), 1.d0/3.d0)
            
            fchi1 = (3.d0*chi-1.d0)/2d0
            fchi2 = (1.d0-chi)/2d0

            c2Nr  = cl_til2*u(i,j,k,MER)

            !print *,CurrentIndex, i,j,k, u(i,j,k,MER), fi(MX)**2.d0/f2, fi(MY)**2.d0/f2, fi(MZ)**2.d0/f2

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
              
              !if(u(i,j,k,MER) > 0.d0 )then
              !  print '(1P10E15.7)', u(i,j,k,MER), PP(MX, MY, i, j, k), PP(MY, MY, i, j, k), PP(MZ, MY, i, j, k) &
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

          !if(u(i,j,k,MER) > 1.d9*Unit_l**3.d0/1.d49) then
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
    do k=Kmin-1, Kmax 
      do j = Jmin-1, Jmax
        do i = Imin-1, Imax

          lam_p(MX) = dmax1(0.d0, lam_max(MX, i, j, k), lam_max(MX, i+1, j, k))
          lam_p(MY) = dmax1(0.d0, lam_max(MY, i, j, k), lam_max(MY, i, j+1, k))
          lam_p(MZ) = dmax1(0.d0, lam_max(MZ, i, j, k), lam_max(MZ, i, j, k+1))

          lam_m(MX) = dmin1(0.d0, lam_min(MX, i, j, k), lam_min(MX, i+1, j, k))
          lam_m(MY) = dmin1(0.d0, lam_min(MY, i, j, k), lam_min(MY, i, j+1, k))
          lam_m(MZ) = dmin1(0.d0, lam_min(MZ, i, j, k), lam_min(MZ, i, j, k+1))

          !FLG test
          !lam_p = 1.d0
          !lam_m = -1.d0

          !if(u(i,j,k,MER) > 1.d9*Unit_l**3.d0/1.d49) then
          !  print *, lam_p, lam_m
          !endif

          do n = MX, MZ
            lcr_mp(n) = lam_p(n)*lam_m(n)*cl_til
            lsa_mp(n) = lam_p(n)-lam_m(n)
          enddo        

#define GETFLUX(Fh, Fl, Flp, Nl, Nlp, M) \
  Fh=(lam_p(M)*Fl-lam_m(M)*Flp+lcr_mp(M)*(Nlp-Nl))/lsa_mp(M)


          ! Flux ---------------------------------------------------------------------
          GETFLUX(F(MX,NER,i,j,k), u(i,j,k,MFRX),u(i+1,j,k,MFRX),u(i,j,k,MER),u(i+1,j,k,MER),MX)
          GETFLUX(F(MY,NER,i,j,k), u(i,j,k,MFRY),u(i,j+1,k,MFRY),u(i,j,k,MER),u(i,j+1,k,MER),MY)
          GETFLUX(F(MZ,NER,i,j,k), u(i,j,k,MFRZ),u(i,j,k+1,MFRZ),u(i,j,k,MER),u(i,j,k+1,MER),MZ)
          
          GETFLUX(F(MX,NFX,i,j,k),PP(MX,MX,i,j,k),PP(MX,MX,i+1,j,k),u(i,j,k,MFRX),u(i+1,j,k,MFRX),MX)
          GETFLUX(F(MY,NFX,i,j,k),PP(MY,MX,i,j,k),PP(MY,MX,i,j+1,k),u(i,j,k,MFRX),u(i,j+1,k,MFRX),MY)
          GETFLUX(F(MZ,NFX,i,j,k),PP(MZ,MX,i,j,k),PP(MZ,MX,i,j,k+1),u(i,j,k,MFRX),u(i,j,k+1,MFRX),MZ)

          GETFLUX(F(MX,NFY,i,j,k),PP(MX,MY,i,j,k),PP(MX,MY,i+1,j,k),u(i,j,k,MFRY),u(i+1,j,k,MFRY),MX)
          GETFLUX(F(MY,NFY,i,j,k),PP(MY,MY,i,j,k),PP(MY,MY,i,j+1,k),u(i,j,k,MFRY),u(i,j+1,k,MFRY),MY)
          GETFLUX(F(MZ,NFY,i,j,k),PP(MZ,MY,i,j,k),PP(MZ,MY,i,j,k+1),u(i,j,k,MFRY),u(i,j,k+1,MFRY),MZ)

          GETFLUX(F(MX,NFZ,i,j,k),PP(MX,MZ,i,j,k),PP(MX,MZ,i+1,j,k),u(i,j,k,MFRZ),u(i+1,j,k,MFRZ),MX)
          GETFLUX(F(MY,NFZ,i,j,k),PP(MY,MZ,i,j,k),PP(MY,MZ,i,j+1,k),u(i,j,k,MFRZ),u(i,j+1,k,MFRZ),MY)
          GETFLUX(F(MZ,NFZ,i,j,k),PP(MZ,MZ,i,j,k),PP(MZ,MZ,i,j,k+1),u(i,j,k,MFRZ),u(i,j,k+1,MFRZ),MZ)

          !if(dabs(F(MX,NFY,i,j,k)) > 0.d0 )then
          !    print *, F(MX,NFY,i,j,k), PP(MX,MY,i,j,k), PP(MX,MY,i+1,j,k), u(i,j,k,MFRY),u(i,j+1,k,MFRY), lam_p(MY), lam_m(MY)
          !endif

          !------------------------------------------------------------------
          !tap = u(i,j,k,MER)+ ((F(MX,NER,i-1,j,k)-F(MX,NER,i,j,k))/h(MX)+ &
          !    (F(MY,NER,i,j-1,k)-F(MY,NER,i,j,k))/h(MY)+ &
          !  (F(MZ,NER,i,j,k-1)-F(MZ,NER,i,j,k))/h(MZ))*dt
          !if(tap < 0.d0 ) then
          !    print *, u(i,j,k,MER), (F(MX,NER,i-1,j,k)-F(MX,NER,i,j,k))/h(MX)*dt &
          !      , (F(MY,NER,i,j-1,k)-F(MY,NER,i,j,k))/h(MY)*dt &
          !      ,(F(MZ,NER,i,j,k-1)-F(MZ,NER,i,j,k))/h(MZ)*dt, lam_p, lam_m

          !    print *, dmax1(0.d0, lam_max(MX, i-1, j, k), lam_max(MX, i, j, k)), &
          !             !dmax1(0.d0, lam_max(MY, i, j-1, k), lam_max(MY, i, j, k)), &
          !             !dmax1(0.d0, lam_max(MZ, i, j, k-1), lam_max(MZ, i, j, k)), &
          !             dmin1(0.d0, lam_min(MX, i-1, j, k), lam_min(MX, i, j, k))!, &
          !             !dmin1(0.d0, lam_min(MY, i, j-1, k), lam_min(MY, i, j, k)), &
          !             !dmin1(0.d0, lam_min(MZ, i, j, k-1), lam_min(MZ, i, j, k))
          !    print *, F(MX,NER,i-1,j,k), u(i-1,j,k,MFRX), u(i,j,k,MFRX), u(i-1,j,k,MER),u(i,j,k,MER)
          !    print *, F(MX,NER,i,j,k), u(i,j,k,MFRX), u(i+1,j,k,MFRX), u(i,j,k,MER),u(i+1,j,k,MER)
          !stop
          !endif
          !-------------------------------------------------------------------


          ! --------------------------------------------------------------------------
#undef GETFLUX

        enddo
      enddo
    enddo

#else
    do k=Kmin, Kmax 
      do j = Jmin, Jmax
        do i = Imin, Imax

          F(MX,NER,i,j,k) = ((u(i-1,j,k,MFRX)-u(i+1,j,k,MFRX))*0.5d0 &
            +cl_th*( u(i+1,j,k,MER)-2.d0*u(i,j,k,MER)+u(i-1,j,k,MER)))/h(MX)
          F(MY,NER,i,j,k) = ((u(i,j-1,k,MFRY)-u(i,j+1,k,MFRY))*0.5d0 &
            +cl_th*( u(i,j+1,k,MER)-2.d0*u(i,j,k,MER)+u(i,j-1,k,MER)))/h(MY)
          F(MZ,NER,i,j,k) = ((u(i,j,k-1,MFRZ)-u(i,j,k+1,MFRZ))*0.5d0 &
            +cl_th*( u(i,j,k+1,MER)-2.d0*u(i,j,k,MER)+u(i,j,k-1,MER) ))/h(MZ)

          F(MX,NFX,i,j,k) = ((PP(MX,MX,i-1,j,k)-PP(MX,MX,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,MFRX)-2.d0*u(i,j,k,MFRX)+u(i-1,j,k,MFRX)))/h(MX)
          F(MY,NFX,i,j,k) = ((PP(MY,MX,i,j-1,k)-PP(MY,MX,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,MFRX)-2.d0*u(i,j,k,MFRX)+u(i,j-1,k,MFRX)))/h(MY)
          F(MZ,NFX,i,j,k) = ((PP(MZ,MX,i,j,k-1)-PP(MZ,MX,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,MFRX)-2.d0*u(i,j,k,MFRX)+u(i,j,k-1,MFRX)))/h(MZ)

          F(MX,NFY,i,j,k) = ((PP(MX,MY,i-1,j,k)-PP(MX,MY,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,MFRY)-2.d0*u(i,j,k,MFRY)+u(i-1,j,k,MFRY)))/h(MX)
          F(MY,NFY,i,j,k) = ((PP(MY,MY,i,j-1,k)-PP(MY,MY,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,MFRY)-2.d0*u(i,j,k,MFRY)+u(i,j-1,k,MFRY)))/h(MY)
          F(MZ,NFY,i,j,k) = ((PP(MZ,MY,i,j,k-1)-PP(MZ,MY,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,MFRY)-2.d0*u(i,j,k,MFRY)+u(i,j,k-1,MFRY)))/h(MZ)

          F(MX,NFZ,i,j,k) = ((PP(MX,MZ,i-1,j,k)-PP(MX,MZ,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,MFRZ)-2.d0*u(i,j,k,MFRZ)+u(i-1,j,k,MFRZ)))/h(MX)
          F(MY,NFZ,i,j,k) = ((PP(MY,MZ,i,j-1,k)-PP(MY,MZ,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,MFRZ)-2.d0*u(i,j,k,MFRZ)+u(i,j-1,k,MFRZ)))/h(MY)
          F(MZ,NFZ,i,j,k) = ((PP(MZ,MZ,i,j,k-1)-PP(MZ,MZ,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,MFRZ)-2.d0*u(i,j,k,MFRZ)+u(i,j,k-1,MFRZ)))/h(MZ)

        enddo
      enddo
    enddo
#endif



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

  !---------------------------------------------
  ! update N & F
  !---------------------------------------------

  subroutine u_update_radtr

    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u
    integer :: i,j,k
    logical :: isNotFinite
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
   
    
    u => get_up( Ulist(CurrentIndex) )
    h   = CellWidth( :, CurrentLevel)

    do k=Kmin, Kmax 
      do j = Jmin, Jmax
        do i = Imin, Imax
          
#ifdef SOLVER_RADTR_HLL
          
#define CAL_FLUX(M, NM) \
          u(i,j,k,M) = u(i,j,k,M) \
           + ((F(MX,NM,i-1,j,k)-F(MX,NM,i,j,k))/h(MX) \
           +  (F(MY,NM,i,j-1,k)-F(MY,NM,i,j,k))/h(MY) \
           +  (F(MZ,NM,i,j,k-1)-F(MZ,NM,i,j,k))/h(MZ))*dt

          CAL_FLUX(MER,  NER)
          CAL_FLUX(MFRX, NFX)
          CAL_FLUX(MFRY, NFY)
          CAL_FLUX(MFRZ, NFZ)
          
#else
          !MER
          u(i,j,k,MER) = u(i,j,k,MER) + (F(MX,NER,i,j,k)+F(MY,NER,i,j,k)+F(MZ,NER,i,j,k))*dt

          !MFR
          u(i,j,k,MFRX)= u(i,j,k,MFRX)+(F(MX,NFX,i,j,k)+F(MY,NFX,i,j,k)+F(MZ,NFX,i,j,k))*dt  
          u(i,j,k,MFRY)= u(i,j,k,MFRY)+(F(MX,NFY,i,j,k)+F(MY,NFY,i,j,k)+F(MZ,NFY,i,j,k))*dt
          u(i,j,k,MFRZ)= u(i,j,k,MFRZ)+(F(MX,NFZ,i,j,k)+F(MY,NFZ,i,j,k)+F(MZ,NFZ,i,j,k))*dt
#endif

          !if(dabs(u(i,j,k,MER)) > 0.d0) then
          !  print *,  CurrentLevel, CurrentIndex, i,j,k, u(i,j,k,MER), u(i,j,k,MFRX), u(i,j,k,MFRY), u(i,j,k,MFRZ)
          !endif
#if DEBUG_RADTR == YES
          if(isNotFinite(u(i,j,k,MER)) .or. isNotFinite(u(i,j,k,MFRX))  .or. isNotFinite(u(i,j,k,MFRY)) .or. &
          isNotFinite(u(i,j,k,MFRZ)))then
              print *, CurrentIndex, u(i,j,k,MER), u(i,j,k,MFRX), u(i,j,k,MFRY), u(i,j,k,MFRZ)
          endif
#endif

        enddo
      enddo
    enddo


  end subroutine u_update_radtr


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
