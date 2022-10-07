module radtr
  use unit
  use grid
  use modelParameter
  use radtr_chem
  use parameter, only : Pi, Pi4, Pi4i
  use radiationSource
  use mpilib
  implicit none
  private
  integer, parameter :: NER = 0, NFX = 1, NFY = 2, NFZ = 3
  real(kind=8), save :: tsum, tend, dt
  integer, save :: CurrentLevel
  integer, save :: CurrentIndex 
  logical, save :: loop_end
  integer,save,private :: STEP_MODE
  integer,save,dimension(Gidmin:Gidmax) :: Ulist
  integer,save :: ListMax, num_step
  real(kind=8), save,dimension(0:2,NER:NFZ,Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: F 
  real(kind=8), save :: cl_til, cl_til2, cl_th
  real(kind=8), save, dimension(:,:), allocatable :: lambda1, lambda4
  type(t_rs_info),pointer, save :: rs_info
  real(kind=8), save,dimension(0:2,Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: Frad 
  integer, parameter :: ncom = 3
  integer, dimension(NER:NFZ,ncom), save :: nrads
    integer, save :: loop_countmax
    real(kind=8), parameter :: pi_circ = 3.14159265
    real(kind=8), save :: del_V
    real(kind=8), save :: norm_fun
    real(kind=8), save, dimension(:), allocatable :: log_cosdisk, log_posx, log_posy, log_posz
    integer, save :: index_neuv
    integer, save :: index_nfuv
    integer, save :: index_nir
  type(rad_mean), save :: rsm
  integer,parameter :: NBOUNDARY = 2*(2 -0 +1) 
  integer,parameter :: NIL = 0 
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  logical,save :: BoolInitialized = .FALSE. 
  integer,save :: Igmin(Lmin:Lmax), Jgmin(Lmin:Lmax), Kgmin(Lmin:Lmax)
  integer,save :: Igmax(Lmin:Lmax), Jgmax(Lmin:Lmax), Kgmax(Lmin:Lmax)
  public :: radtr_moment
contains
  subroutine radtr_init(dt_hyd)
    logical, save :: bool_radinit = .false.
    real(kind=8),intent(IN) :: dt_hyd
    integer :: n, ierr
    integer :: nsource_glob, isrc
    real(kind=8) :: mskr, mskr2
    real(kind=8) :: alpha_EUV, heat_EUV, sig_EUV, sig_FUV, tSion, tlum_euv &
      , tlum_fuv, Sion, Sfuv, lumeuv, lumfuv, alpha_OII, tSfuv
    real(kind=8),parameter :: chi_H = 13.6d0
    if(.not. bool_radinit) then
      cl_til = MP_Ctil / Unit_v 
      cl_til2= cl_til*cl_til
      cl_th = cl_til*0.5d0
      bool_radinit = .true.
      n = 0
      n = n+1
      index_neuv = n
      nrads(NER,n) = 22
      nrads(NFX,n) = 23
      nrads(NFY,n) = 24
      nrads(NFZ,n) = 25
      n = n+1
      index_nfuv = n
      nrads(NER,n) = 26
      nrads(NFX,n) = 27
      nrads(NFY,n) = 28
      nrads(NFZ,n) = 29
      n = n+1
      index_nir = n
      nrads(NER,n) = 30
      nrads(NFX,n) = 31
      nrads(NFY,n) = 32
      nrads(NFZ,n) = 33
    endif
    tsum = 0.d0 
    num_step = 0
    tend = dt_hyd
    loop_end = .false.
    call rs_GetSourceInfo(rs_info)
    nsource_glob = rs_info%nsource
    rsm%num_rad = nsource_glob
    tSion = 0.d0
    tSfuv = 0.d0
    alpha_EUV = 0.d0
    alpha_OII = 0.d0
    heat_EUV = 0.d0
    tlum_euv = 0.d0
    tlum_fuv = 0.d0
    sig_EUV = 0.d0
    sig_FUV = 0.d0
    do isrc = 0, nsource_glob -1
      Sion = rs_info%x_euv(isrc) *rs_info%lum(isrc)
      Sfuv = rs_info%x_fuv(isrc) *rs_info%lum(isrc)
      lumeuv = rs_info%lumeuv(isrc)*rs_info%lum(isrc)
      lumfuv = rs_info%lumfuv(isrc)*rs_info%lum(isrc)
      tSion = tSion + Sion
      tSfuv = tSfuv + Sfuv
      alpha_EUV = alpha_EUV + rs_info%alpha_euv(isrc)*Sion
      heat_EUV = heat_EUV + rs_info%heat_euv(isrc) *Sion
      alpha_OII = alpha_OII + rs_info%alpha_euv(isrc)*Sion*rs_info%rOII(isrc)
      tlum_euv = tlum_euv + lumeuv
      sig_EUV = sig_EUV + rs_info%sig_euv(isrc)*lumeuv
      tlum_fuv = tlum_fuv + lumfuv
      sig_FUV = sig_FUV + rs_info%sig_fuv(isrc)*lumfuv
    enddo
    if(tSion > 0.d0) then
      rsm%rOII = alpha_OII / alpha_EUV
      rsm%alpha_EUV = alpha_EUV / tSion
      rsm%heat_EUV = heat_EUV / tSion
      rsm%sig_EUV = sig_EUV / tlum_euv
      rsm%sig_FUV = sig_FUV / tlum_fuv
    else
      rsm%rOII = 0.d0
      rsm%alpha_EUV = 0.d0
      rsm%heat_EUV = 0.d0
      rsm%sig_EUV = 0.d0
      rsm%sig_FUV = 0.d0
    endif
    rsm%erg_EUV = max(tlum_euv/max(tSion, 1.d-30), 13.6d0*cgs_ev)
    rsm%erg_FUV = max(tlum_fuv/max(tSfuv, 1.d-30), 11.174d0*cgs_ev) 
  end subroutine radtr_init
  subroutine anisotropy_init_new
      use overBlockCoordinates
      use modelParameter
      logical, save :: bool_aniinit = .false.
      integer :: level
      integer :: pos_loopmax
      integer :: loop_count
      integer :: pos_iloop, pos_jloop, pos_kloop
      real(kind=8) :: norm_r2cos
      real(kind=8) :: pos_R2, pos_R
      real(kind=8) :: disk_inc, cos_disk_cell
      real(kind=8), dimension(0:2) :: pos_inj, hmax, disk_axis, pos_vec
      if(.not. bool_aniinit) then
        bool_aniinit = .true.
        level = MP_Lmax0
        hmax = CellWidth(:,level)
        disk_inc = pi_circ / 2.d0 
        disk_axis(0) = cos(disk_inc)
        disk_axis(1) = sin(disk_inc)
        disk_axis(2) = 0.d0
        pos_loopmax = 14
        loop_count = 0
        do pos_iloop=0, pos_loopmax-1
        do pos_jloop=0, pos_loopmax-1
        do pos_kloop=0, pos_loopmax-1
           if (pos_iloop<(pos_loopmax/2)) then
              pos_inj(0) = 0.d0 + ((pos_iloop+1.d0)*hmax(0)-hmax(0)*0.5)
           else
              pos_inj(0) = 0.d0 - ((pos_iloop+1.d0-(pos_loopmax/2))*hmax(0)-hmax(0)*0.5)
           end if
           if (pos_jloop<(pos_loopmax/2)) then
              pos_inj(1) = 0.d0 + ((pos_jloop+1.d0)*hmax(1)-hmax(1)*0.5)
           else
              pos_inj(1) = 0.d0 - ((pos_jloop+1.d0-(pos_loopmax/2))*hmax(1)-hmax(1)*0.5)
           end if
           if (pos_kloop<(pos_loopmax/2)) then
              pos_inj(2) = 0.d0 + ((pos_kloop+1.d0)*hmax(2)-hmax(2)*0.5)
           else
              pos_inj(2) = 0.d0 - ((pos_kloop+1.d0-(pos_loopmax/2))*hmax(2)-hmax(2)*0.5)
           end if
           pos_R2 = pos_inj(0)**2.d0 + pos_inj(1)**2.d0 + pos_inj(2)**2.d0
           pos_R = sqrt(pos_R2)
           if(pos_R<(hmax(0)*7.0-hmax(0)*0.5)) then
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
        log_posx(:) =0.d0
        log_posy(:) =0.d0
        log_posz(:) =0.d0
        log_cosdisk(:) =0.d0
        norm_r2cos = 0.d0
        loop_count = 0
        do pos_iloop=0, pos_loopmax-1
        do pos_jloop=0, pos_loopmax-1
        do pos_kloop=0, pos_loopmax-1
           if (pos_iloop<(pos_loopmax/2)) then
              pos_inj(0) = 0.d0 + ((pos_iloop+1.d0)*hmax(0)-hmax(0)*0.5)
           else
              pos_inj(0) = 0.d0 - ((pos_iloop+1.d0-(pos_loopmax/2))*hmax(0)-hmax(0)*0.5)
           end if
           if (pos_jloop<(pos_loopmax/2)) then
              pos_inj(1) = 0.d0 + ((pos_jloop+1.d0)*hmax(1)-hmax(1)*0.5)
           else
              pos_inj(1) = 0.d0 - ((pos_jloop+1.d0-(pos_loopmax/2))*hmax(1)-hmax(1)*0.5)
           end if
           if (pos_kloop<(pos_loopmax/2)) then
              pos_inj(2) = 0.d0 + ((pos_kloop+1.d0)*hmax(2)-hmax(2)*0.5)
           else
              pos_inj(2) = 0.d0 - ((pos_kloop+1.d0-(pos_loopmax/2))*hmax(2)-hmax(2)*0.5)
           end if
           pos_R2 = pos_inj(0)**2.d0 + pos_inj(1)**2.d0 + pos_inj(2)**2.d0
           pos_R = sqrt(pos_R2)
           if(pos_R<(hmax(0)*7.0-hmax(0)*0.5)) then
              loop_count = loop_count+1
              pos_vec(0) = pos_inj(0)/pos_R
              pos_vec(1) = pos_inj(1)/pos_R
              pos_vec(2) = pos_inj(2)/pos_R
              cos_disk_cell = dot_product(pos_vec,disk_axis)
              log_posx(loop_count) = pos_inj(0)
              log_posy(loop_count) = pos_inj(1)
              log_posz(loop_count) = pos_inj(2)
              if ((pos_inj(1) < 0.d0 .and. pos_jloop==7 ) .or. (pos_inj(1) > 0.d0 .and. pos_jloop==0 )) then 
                 cos_disk_cell = 0.d0
              end if
              log_cosdisk(loop_count) = abs(cos_disk_cell)
              norm_r2cos = norm_r2cos + (2.0*abs(cos_disk_cell) / pos_R2)
           end if
        enddo
        enddo
        enddo
        norm_fun = 1.d0 / norm_r2cos
        end if
  end subroutine anisotropy_init_new
  subroutine radtr_moment(dt_hyd)
    use mpilib
    real(kind=8),intent(IN) :: dt_hyd 
    real(kind=8) :: time_radtr
    integer :: time_prev, time_cur, time_rat, time1, time2, time_chem, time_mpi
       call system_clock(time_prev)
       if(get_myrank() == 0) print '(/, A, /, A)', "start RADTR" &
         , "(RADTR ) ---------------------------------------------------------------------------"
    time_chem = 0
    time_mpi = 0
    call radtr_init(dt_hyd)
    call anisotropy_init_new
    call radforce_init
    do
      call cfl_condition_radrt
      if(get_myrank() == 0) write(*,'(a,I5,1p4e13.5)') "step, tsum, dt, dt_hyd", num_step, dt_hyd/dt, tsum, dt, dt_hyd
        call injection_step
        call timestep_radtr
      call system_clock(time1)
      call ch_radtrchem(dt, dt_hyd, rsm)
      call system_clock(time2)
      time_chem = time_chem + time2-time1
      call system_clock(time1)
      call converge_alllevel
      call system_clock(time2)
      time_mpi = time_mpi + time2-time1
      if (loop_end) exit
    enddo
       call system_clock(time_cur, time_rat)
       time_radtr = (time_cur - time_prev)/dble(time_rat)
       if(get_myrank() == 0) print '(A,/, A, (1PE12.4), A, A, (1PE12.4), A, A, (1PE12.4), A/)' &
         , "----------------------------------------------------------------------------------", &
         "TIME for RADTR:", time_radtr, "[s]", "CHEM: ", time_chem/dble(time_rat), "[s]", "MPI:", time_mpi/dble(time_rat), "[s]"
  end subroutine radtr_moment
  subroutine cfl_condition_radrt
    use grid
    use mpilib
    real(kind=8) :: dtlocal, dtcfl, dtbuf
    real(kind=8),dimension(0:2) :: h
    integer :: lev, n, gid
      dtlocal = huge(dtlocal)
      lev = LevelMax 
      h = CellWidth( :, lev )
      do n=0,2
        dtcfl = 0.8 * h(n) / (3.d0*cl_til)
        dtlocal = min(dtlocal, dtcfl)
      enddo
      call mpi_allreduce(dtlocal, dtbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
      dtlocal = dtbuf
      num_step = num_step + 1
      if(int(10.0)==num_step) then
        dtlocal = tend - tsum
        loop_end = .true.
      else
        if(tsum + dtlocal > tend) then
          dtlocal = tend - tsum
          loop_end = .true.
        endif
      endif
    dt = dtlocal
    tsum = tsum + dtlocal
    if (dtlocal < 0.d0) then
      write(*,*) 'negative dtlocal'
      stop
    endif
  end subroutine cfl_condition_radrt
  subroutine injection_step
    use radiationSource
    use overBlockCoordinates
    use modelParameter, only : MP_spRadius_cell
    integer :: level, n, gid
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(0:2) :: ijkg
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(0:2) :: pos, hmax
    real(kind=8) :: Ndot_ion, IR_rate, absrate_Dir
    integer :: pos_loop
    real(kind=8):: R_polor
    real(kind=8):: cos_disk
    real(kind=8), dimension(:,:,:), pointer :: Nrad_euv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_x_euv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_y_euv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_z_euv
    real(kind=8), dimension(:,:,:), pointer :: Nrad_ir
    real(kind=8), dimension(:,:,:), pointer :: Nflux_x_ir
    real(kind=8), dimension(:,:,:), pointer :: Nflux_y_ir
    real(kind=8), dimension(:,:,:), pointer :: Nflux_z_ir
    real(kind=8), dimension(:,:,:), pointer :: Nrad_fuv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_x_fuv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_y_fuv
    real(kind=8), dimension(:,:,:), pointer :: Nflux_z_fuv
    logical :: isNotFinite
       level = MP_Lmax0
       isrc = 0
       hmax = CellWidth(:,level)
       del_V = hmax(0)*hmax(1)*hmax(2)
       do pos_loop =1, loop_countmax
         pos(0) = log_posx(pos_loop)
         pos(1) = log_posy(pos_loop)
         pos(2) = log_posz(pos_loop)
         cos_disk = log_cosdisk(pos_loop)
         call ob_getIjkgridFromCoordPhys(ijkg, level, pos)
         call get_gid_from_ijkgrid(ijkg(0),ijkg(1),ijkg(2),level,gid,rank)
         if (gid == Undefi) cycle
         if (rank /= get_myrank() ) cycle
         x => get_Xp(gid)
         y => get_Yp(gid)
         z => get_Zp(gid)
        k = int((pos(2)-z(Kmingh))/hmax(2) + 0.5d0)+Kmingh
        j = int((pos(1)-y(Jmingh))/hmax(1) + 0.5d0)+Jmingh
        i = int((pos(0)-x(Imingh))/hmax(0) + 0.5d0)+Imingh
        k = min(max(k,Kmin),Kmax)
        j = min(max(j,Jmin),Jmax)
        i = min(max(i,Imin),Imax)
        R_polor = sqrt(x(i)**2.d0 + y(j)**2.d0 + z(k)**2.d0)
        Nrad_euv => get_Ucomp(22, gid)
        Ndot_ion = rs_info%x_euv(isrc)*rs_info%lum(isrc)*Unit_t/MP_PHON 
        Nrad_euv(i,j,k) = Nrad_euv(i,j,k) + (Ndot_ion * dt * norm_fun * 2.d0*cos_disk / (del_V * R_polor**2.d0))
        Nflux_x_euv => get_Ucomp(23,gid)
        Nflux_y_euv => get_Ucomp(24,gid)
        Nflux_z_euv => get_Ucomp(25,gid)
        Nflux_x_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * x(i)/R_polor
        Nflux_y_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * y(j)/R_polor
        Nflux_z_euv(i,j,k) = Nrad_euv(i,j,k) * cl_til * z(k)/R_polor
        Nrad_fuv => get_Ucomp(26, gid)
        Ndot_ion = rs_info%x_fuv(isrc)*rs_info%lum(isrc)*Unit_t/MP_PHON 
        Nrad_fuv(i,j,k) = Nrad_fuv(i,j,k) + (Ndot_ion * dt * norm_fun * 2.d0*cos_disk / (del_V * R_polor**2.d0))
        Nflux_x_fuv => get_Ucomp(27,gid)
        Nflux_y_fuv => get_Ucomp(28,gid)
        Nflux_z_fuv => get_Ucomp(29,gid)
        Nflux_x_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *x(i)/R_polor
        Nflux_y_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *y(j)/R_polor
        Nflux_z_fuv(i,j,k) = Nrad_fuv(i,j,k) * cl_til *z(k)/R_polor
        Nrad_ir => get_Ucomp(30, gid)
        IR_rate = 1.d0-rs_info%lumeuv(isrc)-rs_info%lumfuv(isrc)
        IR_rate = min(max(IR_rate, 0.d0), 1.d0)
        Ndot_ion = rs_info%lum(isrc)*IR_rate/MP_hnu_IR*Unit_t/MP_PHON 
        Nrad_ir(i,j,k) = Nrad_ir(i,j,k) + (Ndot_ion*dt*norm_fun*2.d0*cos_disk/(del_V *R_polor**2.d0))
         Nflux_x_ir => get_Ucomp(31,gid)
         Nflux_y_ir => get_Ucomp(32,gid)
         Nflux_z_ir => get_Ucomp(33,gid)
         Nflux_x_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *x(i)/R_polor
         Nflux_y_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *y(j)/R_polor
         Nflux_z_ir(i,j,k) = Nrad_ir(i,j,k) * cl_til *z(k)/R_polor
       enddo
  end subroutine injection_step
  subroutine injection_step_test
    use radiationSource
    use overBlockCoordinates
    integer :: level, n, gid
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(0:2) :: ijkg
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(0:2) :: pos, h
    real(kind=8) :: Ndot_ion, del_V
    real(kind=8), dimension(:,:,:), pointer :: Nrad, Frx, Fry, Frz
    real(kind=8) :: x_start, N_const, radius_c, radius
    x_start = 0.03*cgs_pc/Unit_l
    N_const = 1.d49*Unit_t/MP_PHON
    radius_c= 0.02*cgs_pc/Unit_l
    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList(n, level) 
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)
        Nrad => get_Ucomp(22, gid)
        Frx => get_Ucomp(23,gid)
        Fry => get_Ucomp(24,gid)
        Frz => get_Ucomp(25,gid)
        h = CellWidth( :, level)
        del_V = h(0)*h(1)*h(2)
        do k = Kmin, Kmax
          if(z(k)-h(2) .le. 0.d0 .and. 0.d0 < z(k)+h(2)) then
            do j = Jmin, Jmax
                do i = Imin, Imax
                  if(y(j)-h(1)*0.5d0 .le. -x_start .and. -x_start < y(j)+h(1)*0.5d0 ) then
                    if( x(i)-h(0)*0.5d0 .le. 0.d0 .and. 0.d0 < x(i)+h(0)*0.5d0) then
                      Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
                      Frx(i,j,k) = 0.d0 
                      Fry(i,j,k) = cl_til * Nrad(i,j,k)
                      Frz(i,j,k) = 0.d0
                    endif
                  endif
                  if(y(j)-h(1)*0.5d0 .le. 0.d0 .and. 0.d0 < y(j)+h(1)*0.5d0 ) then
                    if( x(i)-h(0)*0.5d0 .le. -x_start .and. -x_start < x(i)+h(0)*0.5d0) then
                      Nrad(i,j,k) = Nrad(i,j,k) + N_const/del_V*dt
                      Frx(i,j,k) = cl_til * Nrad(i,j,k) 
                      Fry(i,j,k) = 0.d0 
                      Frz(i,j,k) = 0.d0
                    endif
                  endif
                enddo
            enddo
          endif
        enddo
      enddo
    enddo
  end subroutine injection_step_test
  subroutine timestep_radtr
    integer :: level
    do level = Lmin, Lmax
      call step_all_grid_radtr(level)
    enddo
  end subroutine timestep_radtr
  subroutine timestep_radtr_init
    integer :: n
    do n = Gidmin, GidListMax( CurrentLevel )
      Ulist(n) = GidList(n, CurrentLevel )
    enddo
    ListMax = GidListMax( CurrentLevel )
  end subroutine timestep_radtr_init
  subroutine step_all_grid_radtr(level)
    integer, intent(IN) :: level
    integer :: nr, ier, ifx, ify, ifz
    CurrentLevel = level
    call timestep_radtr_init
    STEP_MODE = PREDICTOR
    call boundary_cond
    call rescueLev_radtr
    do CurrentIndex = Gidmin, ListMax 
      call backup_u_2order
    enddo
    do CurrentIndex = Gidmin, ListMax 
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
  subroutine boundary_cond
    use grid_boundary
    use boundary
    integer :: n
    call boundary_grid( CurrentLevel, STEP_MODE )
    do n = Gidmin, ListMax 
      call boundary_u_radtr(Ulist(n))
    enddo
  end subroutine boundary_cond
  subroutine boundary_u_radtr(id)
    use grid
    integer, intent(IN) :: id
    integer :: i,j,k,m
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    real(kind=8),dimension(:,:,:,:),pointer :: u
    call boundary_radtr_init
    call touch_boundary_radtr( id, bool_touch )
    if (.not. any(bool_touch) ) return
    u => get_Up(id)
    if ( bool_touch(NIL) ) then
       do i = Imingh,Imin-1
          u(i,:,:,22) = 0.d0
          u(i,:,:,23) = u(Imin,:,:,23) 
          u(i,:,:,24) = u(Imin,:,:,24)
          u(i,:,:,25) = u(Imin,:,:,25)
          u(i,:,:,26) = 0.d0
          u(i,:,:,27) = u(Imin,:,:,27)
          u(i,:,:,28) = u(Imin,:,:,28)
          u(i,:,:,29) = u(Imin,:,:,29)
          u(i,:,:,30) = 0.d0
          u(i,:,:,31) = u(Imin,:,:,31)
          u(i,:,:,32) = u(Imin,:,:,32)
          u(i,:,:,33) = u(Imin,:,:,33)
       enddo
    endif
    if ( bool_touch(NIR) ) then
       do i = Imax+1,Imaxgh
          u(i,:,:,22) = 0.d0
          u(i,:,:,23) = u(Imax,:,:,23) 
          u(i,:,:,24) = u(Imax,:,:,24) 
          u(i,:,:,25) = u(Imax,:,:,25) 
          u(i,:,:,26) = 0.d0
          u(i,:,:,27) = u(Imax,:,:,27)
          u(i,:,:,28) = u(Imax,:,:,28)
          u(i,:,:,29) = u(Imax,:,:,29)
          u(i,:,:,30) = 0.d0
          u(i,:,:,31) = u(Imax,:,:,31)
          u(i,:,:,32) = u(Imax,:,:,32)
          u(i,:,:,33) = u(Imax,:,:,33)
       enddo
    endif
    if ( bool_touch(NJL) ) then
       do j = Jmingh,Jmin-1
          u(:,j,:,22) = 0.d0
          u(:,j,:,23) = u(:,Jmin,:,23) 
          u(:,j,:,24) = u(:,Jmin,:,24) 
          u(:,j,:,25) = u(:,Jmin,:,25) 
          u(:,j,:,26) = 0.d0
          u(:,j,:,27) = u(:,Jmin,:,27)
          u(:,j,:,28) = u(:,Jmin,:,28)
          u(:,j,:,29) = u(:,Jmin,:,29)
          u(:,j,:,30) = 0.d0
          u(:,j,:,31) = u(:,Jmin,:,31)
          u(:,j,:,32) = u(:,Jmin,:,32)
          u(:,j,:,33) = u(:,Jmin,:,33)
       enddo
    endif
    if ( bool_touch(NJR) ) then
       do j = Jmax+1,Jmaxgh
          u(:,j,:,22) = 0.d0
          u(:,j,:,23) = u(:,Jmax,:,23)
          u(:,j,:,24) = u(:,Jmax,:,24)
          u(:,j,:,25) = u(:,Jmax,:,25)
          u(:,j,:,26) = 0.d0
          u(:,j,:,27) = u(:,Jmax,:,27)
          u(:,j,:,28) = u(:,Jmax,:,28)
          u(:,j,:,29) = u(:,Jmax,:,29)
          u(:,j,:,30) = 0.d0
          u(:,j,:,31) = u(:,Jmax,:,31)
          u(:,j,:,32) = u(:,Jmax,:,32)
          u(:,j,:,33) = u(:,Jmax,:,33)
       enddo
    endif
    if ( bool_touch(NKL) ) then
       do k = Kmingh,Kmin-1
          u(:,:,k,22) = 0.d0
          u(:,:,k,23) = u(:,:,Kmin,23)
          u(:,:,k,24) = u(:,:,Kmin,24)
          u(:,:,k,25) = u(:,:,Kmin,25)
          u(:,:,k,26) = 0.d0
          u(:,:,k,27) = u(:,:,Kmin,27)
          u(:,:,k,28) = u(:,:,Kmin,28)
          u(:,:,k,29) = u(:,:,Kmin,29)
          u(:,:,k,30) = 0.d0
          u(:,:,k,31) = u(:,:,Kmin,31)
          u(:,:,k,32) = u(:,:,Kmin,32)
          u(:,:,k,33) = u(:,:,Kmin,33)
       enddo
    endif
    if ( bool_touch(NKR) ) then
       do k = Kmax+1,Kmaxgh
          u(:,:,k,22) = 0.d0
          u(:,:,k,23) = u(:,:,Kmax,23)
          u(:,:,k,24) = u(:,:,Kmax,24)
          u(:,:,k,25) = u(:,:,Kmax,25)
          u(:,:,k,26) = 0.d0
          u(:,:,k,27) = u(:,:,Kmax,27)
          u(:,:,k,28) = u(:,:,Kmax,28)
          u(:,:,k,29) = u(:,:,Kmax,29)
          u(:,:,k,30) = 0.d0
          u(:,:,k,31) = u(:,:,Kmax,31)
          u(:,:,k,32) = u(:,:,Kmax,32)
          u(:,:,k,33) = u(:,:,Kmax,33)
       enddo
    endif
  end subroutine boundary_u_radtr
  subroutine touch_boundary_radtr( id, bool_touch )
    use grid
    integer,intent(IN) :: id
    logical,dimension(0:NBOUNDARY-1),intent(OUT) :: bool_touch
    integer :: level, ig, jg, kg
    bool_touch(:) = .FALSE.
    level = get_level(id)
    if ( Igrid(id) == Igmin( level ) ) bool_touch(NIL) = .TRUE.
    if ( Igrid(id) == Igmax( level ) ) bool_touch(NIR) = .TRUE.
    if ( Jgrid(id) == Jgmin( level ) ) bool_touch(NJL) = .TRUE.
    if ( Jgrid(id) == Jgmax( level ) ) bool_touch(NJR) = .TRUE.
    if ( Kgrid(id) == Kgmin( level ) ) bool_touch(NKL) = .TRUE.
    if ( Kgrid(id) == Kgmax( level ) ) bool_touch(NKR) = .TRUE.
  end subroutine touch_boundary_radtr
  subroutine boundary_radtr_init
    use grid
    use io_util
    integer,parameter :: ni0 = 8, nj0 = 8, nk0 = 8
    integer :: level
    logical, save :: BoolInitialized = .False.
    if(BoolInitialized) return
    call print_msg( 'initialize boundary_radtr' )
    BoolInitialized = .TRUE.
    Igmin(:) = 0
    Jgmin(:) = 0
    Kgmin(:) = 0
    do level = Lmin, Lmax
       Igmax(level) = ni0*2**(level-Lmin) -1
       Jgmax(level) = nj0*2**(level-Lmin) -1
       Kgmax(level) = nk0*2**(level-Lmin) -1
    enddo
  end subroutine boundary_radtr_init
  subroutine backup_u_2order
    real(kind=8),dimension(:,:,:,:),pointer :: u, u2
    u => get_up( Ulist(CurrentIndex ) )
    u2 => get_u2orderp( Ulist(CurrentIndex) )
    u2 = u
    U2_StepNumber(CurrentLevel) = U_StepNumber(CurrentLevel)
    U2_StepNumberGhostCell(CurrentLevel) = U_StepNumberGhostCell(CurrentLevel)
  end subroutine backup_u_2order
  subroutine rescueLev_radtr
    use kinzoku, only : fdust_solar
    integer :: n, nr, gid, i, j, k
    real(kind=8),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
    real(kind=8) :: f2, reduc_f2, f2_sq
    real(kind=8),parameter :: Nrad_floor = 1.d-20 
    logical :: isNotFinite
    do n = Gidmin, GidListMax( CurrentLevel )
      gid = GidList(n, CurrentLevel) 
      do nr = 1, ncom
        Nrad => get_Ucomp(nrads(NER,nr),gid)
        Frx => get_Ucomp(nrads(NFX,nr),gid)
        Fry => get_Ucomp(nrads(NFY,nr),gid)
        Frz => get_Ucomp(nrads(NFZ,nr),gid)
        do k = Kmingh, Kmaxgh
          do j = Jmingh, Jmaxgh
            do i = Imingh, Imaxgh
              f2 = Frx(i,j,k)**2 + Fry(i,j,k)**2 + Frz(i,j,k)**2
              f2_sq = dsqrt(f2)
              if(f2_sq > (cl_til*Nrad(i,j,k)) .and. f2_sq > 0.d0) then
                 Frx(i,j,k) = Frx(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
                 Fry(i,j,k) = Fry(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
                 Frz(i,j,k) = Frz(i,j,k)*cl_til*Nrad(i,j,k)/f2_sq
              endif
              if(Nrad(i,j,k) < 0.0) then 
                  Nrad(i,j,k) = 0.d0
                  Frx(i,j,k) = 0.d0
                  Fry(i,j,k) = 0.d0
                  Frz(i,j,k) = 0.d0
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine rescueLev_radtr
  subroutine set_Erad_nRS
    implicit none
    integer :: level, n, gid, nr
    real(kind=8),dimension(:,:,:),pointer :: Nrad, Frx, Fry, Frz
    real(kind=8) :: nrad_cmb
    nrad_cmb = cgs_asb*2.73**4.d0/MP_hnu_IR/MP_PHON*Unit_l3/MP_Crd 
    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList(n, level) 
        do nr = 1, ncom
          Nrad => get_Ucomp(nrads(NER,nr),gid)
          Frx => get_Ucomp(nrads(NFX,nr),gid)
          Fry => get_Ucomp(nrads(NFY,nr),gid)
          Frz => get_Ucomp(nrads(NFZ,nr),gid)
          if (nr == index_nir) then
            Nrad(:,:,:) = nrad_cmb
          else
            Nrad(:,:,:) = 0.d0
          endif
          Frx(:,:,:) = 0.d0
          Fry(:,:,:) = 0.d0
          Frz(:,:,:) = 0.d0
        enddo
      enddo
    enddo
  end subroutine set_Erad_nRS
  subroutine get_flux_radtr(ier, ifx, ify, ifz)
    integer, intent(IN) :: ier, ifx, ify, ifz
    integer :: i,j,k,n
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8), dimension(0:2, 0:2, Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: PP 
    real(kind=8), dimension(0:2) :: fi
    real(kind=8) :: inv_cn, f2, f2_sq, ff, ff2, chi, fchi1, fchi2, c2Nr
    real(kind=8),dimension(0:2) :: h
    real(kind=8) :: tap
    logical :: isNotFinite
    u => get_up( Ulist(CurrentIndex) )
    h = CellWidth( :, CurrentLevel)
    do k = Kmingh, Kmaxgh
      do j = Jmingh, Jmaxgh
        do i = Imingh, Imaxgh
          if(u(i,j,k,ier) .le. 0.d0) then
              PP(:, :, i, j, k) = 0.d0
          else
            fi(0) = u(i,j,k,ifx)
            fi(1) = u(i,j,k,ify)
            fi(2) = u(i,j,k,ifz)
            f2 = fi(0)**2.d0+fi(1)**2.d0+fi(2)**2.d0
            f2_sq = dsqrt(f2)
            ff = f2_sq/(cl_til*u(i,j,k,ier))
            ff2 = ff**2.d0
            chi = dmax1(4.d0-3.d0*ff2, 0.d0)
            chi = (3.d0+4.d0*ff2)/(5.d0+2.d0*dsqrt(chi))
            chi = dmax1(dmin1(1.d0, chi), 1.d0/3.d0)
            fchi1 = (3.d0*chi-1.d0)/2d0
            fchi2 = (1.d0-chi)/2d0
            c2Nr = cl_til2*u(i,j,k,ier)
            if(f2 > 0.d0) then
              PP(0, 0, i, j, k) = (fchi1*fi(0)**2.d0 /f2 + fchi2)*c2Nr
              PP(1, 0, i, j, k) = (fchi1*fi(1)*fi(0)/f2 )*c2Nr
              PP(2, 0, i, j, k) = (fchi1*fi(2)*fi(0)/f2 )*c2Nr
              PP(0, 1, i, j, k) = PP(1, 0, i, j, k)
              PP(1, 1, i, j, k) = (fchi1*fi(1)**2.d0 /f2 + fchi2)*c2Nr
              PP(2, 1, i, j, k) = (fchi1*fi(2)*fi(1)/f2 )*c2Nr
              PP(0, 2, i, j, k) = PP(2, 0, i, j, k)
              PP(1, 2, i, j, k) = PP(2, 1, i, j, k)
              PP(2, 2, i, j, k) = (fchi1*fi(2)**2.d0 /f2 + fchi2)*c2Nr
            else
              PP(0, 0, i, j, k) = fchi2*c2Nr
              PP(1, 0, i, j, k) = 0.d0
              PP(2, 0, i, j, k) = 0.d0
              PP(0, 1, i, j, k) = PP(1, 0, i, j, k)
              PP(1, 1, i, j, k) = fchi2*c2Nr
              PP(2, 1, i, j, k) = 0.d0
              PP(0, 2, i, j, k) = PP(2, 0, i, j, k)
              PP(1, 2, i, j, k) = PP(2, 1, i, j, k)
              PP(2, 2, i, j, k) = fchi2*c2Nr
            endif
          endif
        enddo
      enddo
    enddo
    do k=Kmingh+1, Kmaxgh-1
      do j = Jmingh+1, Jmaxgh-1
        do i = Imingh+1, Imaxgh-1
          F(0,NER,i,j,k) = ((u(i-1,j,k,ifx)-u(i+1,j,k,ifx))*0.5d0 &
            +cl_th*( u(i+1,j,k,ier)-2.d0*u(i,j,k,ier)+u(i-1,j,k,ier)))/h(0)
          F(1,NER,i,j,k) = ((u(i,j-1,k,ify)-u(i,j+1,k,ify))*0.5d0 &
            +cl_th*( u(i,j+1,k,ier)-2.d0*u(i,j,k,ier)+u(i,j-1,k,ier)))/h(1)
          F(2,NER,i,j,k) = ((u(i,j,k-1,ifz)-u(i,j,k+1,ifz))*0.5d0 &
            +cl_th*( u(i,j,k+1,ier)-2.d0*u(i,j,k,ier)+u(i,j,k-1,ier) ))/h(2)
          F(0,NFX,i,j,k) = ((PP(0,0,i-1,j,k)-PP(0,0,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ifx)-2.d0*u(i,j,k,ifx)+u(i-1,j,k,ifx)))/h(0)
          F(1,NFX,i,j,k) = ((PP(1,0,i,j-1,k)-PP(1,0,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ifx)-2.d0*u(i,j,k,ifx)+u(i,j-1,k,ifx)))/h(1)
          F(2,NFX,i,j,k) = ((PP(2,0,i,j,k-1)-PP(2,0,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ifx)-2.d0*u(i,j,k,ifx)+u(i,j,k-1,ifx)))/h(2)
          F(0,NFY,i,j,k) = ((PP(0,1,i-1,j,k)-PP(0,1,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ify)-2.d0*u(i,j,k,ify)+u(i-1,j,k,ify)))/h(0)
          F(1,NFY,i,j,k) = ((PP(1,1,i,j-1,k)-PP(1,1,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ify)-2.d0*u(i,j,k,ify)+u(i,j-1,k,ify)))/h(1)
          F(2,NFY,i,j,k) = ((PP(2,1,i,j,k-1)-PP(2,1,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ify)-2.d0*u(i,j,k,ify)+u(i,j,k-1,ify)))/h(2)
          F(0,NFZ,i,j,k) = ((PP(0,2,i-1,j,k)-PP(0,2,i+1,j,k))*0.5d0 &
            +cl_th*(u(i+1,j,k,ifz)-2.d0*u(i,j,k,ifz)+u(i-1,j,k,ifz)))/h(0)
          F(1,NFZ,i,j,k) = ((PP(1,2,i,j-1,k)-PP(1,2,i,j+1,k))*0.5d0 &
            +cl_th*(u(i,j+1,k,ifz)-2.d0*u(i,j,k,ifz)+u(i,j-1,k,ifz)))/h(1)
          F(2,NFZ,i,j,k) = ((PP(2,2,i,j,k-1)-PP(2,2,i,j,k+1))*0.5d0 &
            +cl_th*(u(i,j,k+1,ifz)-2.d0*u(i,j,k,ifz)+u(i,j,k-1,ifz)))/h(2)
        enddo
      enddo
    enddo
    do k=Kmin, Kmax
      do j = Jmin, Jmax
        do i = Imin, Imax
          u(i,j,k,ier) = u(i,j,k,ier) + (F(0,NER,i,j,k)+F(1,NER,i,j,k)+F(2,NER,i,j,k))*dt
          u(i,j,k,ifx)= u(i,j,k,ifx)+(F(0,NFX,i,j,k)+F(1,NFX,i,j,k)+F(2,NFX,i,j,k))*dt
          u(i,j,k,ify)= u(i,j,k,ify)+(F(0,NFY,i,j,k)+F(1,NFY,i,j,k)+F(2,NFY,i,j,k))*dt
          u(i,j,k,ifz)= u(i,j,k,ifz)+(F(0,NFZ,i,j,k)+F(1,NFZ,i,j,k)+F(2,NFZ,i,j,k))*dt
        enddo
      enddo
    enddo
  end subroutine get_flux_radtr
  subroutine get_eigenvals(f2, omega, lam_min, lam_max)
    real(kind=8), intent(in) :: f2, omega
    real(kind=8) :: lam_min, lam_max
    real(kind=8) :: theta,dd1,dd2,de1,de2,lff,ltt
    integer::ii,jj
    theta=ACOS(omega)
    lff = f2*1.d2
    ltt = theta/Pi*1.d2
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
  subroutine radforce_init
    integer :: level, n, gid
    real(kind=8),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    do level = Lmin, Lmax
      do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        fx_hpi => get_Ucomp(14,gid)
        fy_hpi => get_Ucomp(15,gid)
        fz_hpi => get_Ucomp(16,gid)
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
    integer :: i, j, k, isrc, rank
    integer :: nsource_glob
    integer,dimension(0:2) :: ijkg
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(0:2) :: pos
    real(kind=8),dimension(:,:,:),pointer :: fx_hpi, fy_hpi, fz_hpi
    real(kind=8) :: SinkRadius, r2
    SinkRadius = sp_getSinkRadius()
    do level = Lmin, Lmax
      nsource_glob = rs_info%nsource
      do isrc = 0, nsource_glob -1
        pos = rs_info%spos(:,isrc)
        call ob_getIjkgridFromCoordPhys(ijkg, level, pos)
        call get_gid_from_ijkgrid(ijkg(0),ijkg(1),ijkg(2),level,gid,rank)
        if (gid == Undefi) cycle
        if (rank /= get_myrank() ) cycle
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)
        fx_hpi => get_Ucomp(14,gid)
        fy_hpi => get_Ucomp(15,gid)
        fz_hpi => get_Ucomp(16,gid)
        do k=Kmingh, Kmaxgh
          do j = Jmingh, Jmaxgh
            do i = Imingh, Imaxgh
                r2 = (x(i)-pos(0))**2 + (y(j)-pos(1))**2 + (z(k)-pos(2))**2
                if (r2 < SinkRadius**2) then
                    fx_hpi(i,j,k) = 0.d0
 fy_hpi(i,j,k) = 0.d0
 fz_hpi(i,j,k) = 0.d0
                endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine no_radforce_inside_sink
  subroutine converge_alllevel
    use fg2cg
    integer :: level
    do level = LevelMax, Lmin+1, -1
       call fg2cg_u( level )
    end do
  end subroutine converge_alllevel
end module radtr
