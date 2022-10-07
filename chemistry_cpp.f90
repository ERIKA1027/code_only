module chemistry
  use grid
  use parameter,only : Pi
  use unit
  use mpilib
  use kinzoku, only: fdust_solar
  use radiationSource
  use primordial,only: CoolSolverExplicit, CoolSolverImplicit, adjust_abundance &
    , yHe, c_H2, dbg_flg_prm, chem_vals, get_xmu
  use modelParameter, only: MP_PHON, MP_Metallicity, MP_mu, MP_frac_COsum
  implicit none
  private
  real(kind=8),parameter :: fdt = 1d0 
  integer,parameter :: n_substep_max = 3 
  real(kind=8),save :: t_chemcool_local,t_chemcool_glob 
  type(t_rs_info),save,pointer :: rs_info 
  integer,save :: ifirst = 0
  integer :: time_ini,time_prev, time_cur,time_rat
  integer :: num_cell 
  integer :: num_implicit 
  public :: ch_artchem
  public :: ch_CellChemCool
contains
  subroutine ch_artchem(dt_code)
    use art
    real(kind=8),intent(IN) :: dt_code 
    integer :: n_substep 
    integer :: sstep
    real(kind=8):: dt, dt_sub, t_sub
    real(kind=8):: time_chem, time_chem_min,time_chem_mean,time_chem_max,&
         time_cell,time_cell_min,time_cell_mean,time_cell_max,&
         frac_imp,frac_imp_min,frac_imp_mean,frac_imp_max, dt_code_sub, time_radtr
    call system_clock(time_ini) 
    n_substep = 1 
    dt = dt_code*Unit_t
    t_sub=0d0
    dt_sub = dt/n_substep
    dt_code_sub = dt_code/n_substep
    call rs_GetSourceInfo(rs_info)
    call ch_check_all_cells() 
    do sstep = 0, n_substep - 1
       if(get_myrank() == 0) &
            print '(/,A,4((1P1E15.7)))', "(ART_CHEM) dt, dt_sub, t_sub (sec) = ", dt, dt_sub, t_sub
       if(get_myrank() == 0) &
            print *, 'Begin chem update'
       call system_clock(time_prev) 
       num_cell = 0
       num_implicit = 0
       call ch_PrimChemistry(dt_sub)
       call system_clock(time_cur,time_rat) 
       time_chem = (time_cur-time_prev)/dble(time_rat) 
       time_cell = time_chem/dble(num_cell) 
       frac_imp = dble(num_implicit)/dble(num_cell) 
       call mpi_allreduce(time_chem, time_chem_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_chem, time_chem_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(time_cell, time_cell_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(frac_imp, frac_imp_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       time_chem_mean = time_chem_mean/400
       time_cell_mean = time_cell_mean/400
       frac_imp_mean = frac_imp_mean/400
       if(get_myrank() == 0) then
          print '(A,3(A,1P3E10.2))', "TIME for CHEM (min,mean,max)  => ", &
               "tot time: ",time_chem_min,time_chem_mean,time_chem_max, &
               " [s], time/cell: ",time_cell_min,time_cell_mean,time_cell_max, &
               " [s], frac_imp: ", frac_imp_min,frac_imp_mean,frac_imp_max
       end if
       t_sub = t_sub + dt_sub
       if(dt-t_sub < TINY(1d0)) exit
    end do
    call system_clock(time_cur,time_rat) 
    if(get_myrank() == 0) then
       print '(A,I4,A,2(1PE12.4,A))',"TIME for ARTCHEM   =>    myrank: ", myrank, ", tot time: ",(time_cur-time_ini)/dble(time_rat)&
&,"[s]"
    end if
  end subroutine ch_artchem
  subroutine ch_GetDtSub(dt,t_sub,dt_sub)
    real(kind=8),intent(IN) :: dt,t_sub
    real(kind=8),intent(OUT) :: dt_sub
    call mpi_allreduce( t_chemcool_local, t_chemcool_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    if (dt-t_sub < fdt * t_chemcool_glob) then
       dt_sub = dt-t_sub
    else
       dt_sub = MAX(fdt * t_chemcool_glob,dt/n_substep_max)
    end if
    t_chemcool_local = HUGE(1d0)
  end subroutine ch_GetDtSub
  subroutine ch_PrimChemistry(dt)
    use modelParameter,only : MP_Tmin,MP_Tmax
    real(kind=8),intent(IN) :: dt
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:),pointer :: rho, p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
    real(kind=8),dimension(:,:,:),pointer :: Tdust
    real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax) :: mu, khmpd
    real(kind=8),dimension(:,:,:),pointer :: kh2pd
    real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax) :: kh2pd_dbg
    real(kind=8),dimension(0:6 -1) :: ychem0
    real(kind=8) :: xmu, xNc_H2, cs, T_K0
    integer :: level, n, gid, i, j, k, ic
    logical :: isNotFinite 
    logical :: naninf_flg
    integer :: count_output=0, max_count=100 
    real(kind=8) :: probe_radius, yco0
    real(kind=8),dimension(:,:,:),pointer :: yco
    real(kind=8),dimension(:,:,:),pointer :: kgfuv, kdph, kdco, krOII
    real(kind=8) :: dt_code
    real(kind=8),dimension(:,:,:),pointer :: vx,vy,vz
    real(kind=8),dimension(:,:,:),pointer :: Nrad_IR
    type(chem_vals) :: xch
    type chemsp
      real(kind=8),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(5:10)
    integer :: ichem
    xch%dt = dt
    naninf_flg = .False.
    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )
          if (has_child_grid(gid)) cycle
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)
          rho => get_Ucomp(0,gid)
          p => get_Ucomp(4,gid)
          do ichem = 5, 10
            chemary3(ichem)%y => get_Ucomp(ichem,gid)
          enddo
          yco => get_Ucomp(11, gid)
          khpi => get_Ucomp(12,gid)
          heat_hpi => get_Ucomp(13,gid)
          Tdust => get_Ucomp(21,gid)
          fx_hpi => get_Ucomp(14,gid)
          fy_hpi => get_Ucomp(15,gid)
          fz_hpi => get_Ucomp(16,gid)
          vx => get_Ucomp(1,gid)
          vy => get_Ucomp(2,gid)
          vz => get_Ucomp(3,gid)
          Nrad_IR => get_Ucomp(30, gid)
         call GetKHmpdThin(x,y,z,khmpd) 
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax
                   xch%metal = MP_Metallicity
                   xch%Td = Tdust(i,j,k)
                   do ichem = 0, 6 -1
                      xch%ychem(ichem) = chemary3(ichem+5)%y(i,j,k)
                   enddo
                   xch%yco = yco(i,j,k)
                   call adjust_abundance(xch%ychem &
                      , xch%yco &
                      , xch%metal)
                   xch%nH = rho(i,j,k)*Unit_rho/(MP_mu*cgs_amu) 
                   xmu = get_xmu(xch%ychem) 
                   xch%Tg = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb 
                   if (xch%Tg < MP_Tmin*0.999d0 .or. xch%Tg > MP_Tmax*1.001d0) then 
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** before chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)
                   if (khpi(i,j,k) > 0.d0) then
                      xch%heat = heat_hpi(i,j,k)/khpi(i,j,k) 
                      xch%rHpi = khpi(i,j,k)
                   else
                      xch%heat = 0.d0
                      xch%rHpi = 0.d0
                   end if
                   cs = sqrt(cgs_kb*xch%Tg/(xmu*cgs_mh)) 
                   xch%xlmbdj = cs*sqrt(Pi/(cgs_gc*rho(i,j,k)*Unit_rho)) 
                   xch%xNcH = xch%nH*0.5d0*xch%xlmbdj 
                   xNc_H2 = xch%ychem(1)*xch%nH*0.5*xch%xlmbdj 
                   xch%dvdr(0) = dabs(vx(i+1,j,k)-vx(i-1,j,k))/(2.d0*CellWidth(0,level))/Unit_t 
                   xch%dvdr(1) = dabs(vy(i,j+1,k)-vy(i,j-1,k))/(2.d0*CellWidth(1,level))/Unit_t 
                   xch%dvdr(2) = dabs(vz(i,j,k+1)-vz(i,j,k-1))/(2.d0*CellWidth(2,level))/Unit_t 
                   xch%rgfuv = 0.d0
                   xch%rdph = 0.d0
                   xch%rcopd = 0.d0
                   xch%rOII = 0.d0
                   xch%fd = MP_Metallicity
                   xch%rH2pd = 0.d0 
                   xch%rHmpd = 0.d0 
                   xch%EradIR = Nrad_IR(i,j,k)/Unit_l3*MP_PHON*MP_hnu_IR 
                   xch%EradIR = max(xch%EradIR, 0.d0)
                   ychem0(:) = xch%ychem(:)
                   yco0 = xch%yco
                   T_K0 = xch%Tg
                   if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                      naninf_flg = .True.
                   end if
                   do ic=0,6 -1
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
                   if(isNotFinite(xch%rdph)) then 
                      print *, "dph is Nan and dph = 0.d0 for safety at ", x(i)*Unit_au,y(j)*Unit_au,z(k)*Unit_au
                      xch%rdph = 0.d0
                   endif
                   call ch_CellChemCool(xch)
                   num_cell = num_cell + 1 
                   if (isNotFinite(xch%nH) .or. isNotFinite(xch%Tg)) then
                      naninf_flg = .True.
                   end if
                   if (isNotFinite(xch%yco)) then
                      naninf_flg = .True.
                   endif
                   if(isNotFinite(xch%Td))then
                      naninf_flg = .True.
                   endif
                   do ic=0,6 -1
                      if (isNotFinite(xch%ychem(ic))) then
                         naninf_flg = .True.
                      end if
                   end do
                   if (naninf_flg) then
                      print '(A, (1P6E15.7))', "(chemistry) NaN/Inf found after chem update: ", xch%ychem(:)
                      print '((1P10E15.7))', xch%nH,T_K0,khpi(i,j,k),xch%heat,x(i),y(j),z(k), xch%dt
                      print '((1P10E15.7))', xch%ychem(:)
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
                   end if
                   if (xch%Tg < MP_Tmin .or. xch%Tg > MP_Tmax) then
                      if (count_output < max_count) then
                         print '(A, 6(1PE12.4))', &
                              "*** after chem_update: T lower/upper bound imposed *** (T_K, MP_Tmin, MP_Tmax, x, y,z)", &
                              xch%Tg, MP_Tmin, MP_Tmax, x(i),y(j),z(k)
                         count_output = count_output+1
                         if (count_output == max_count) then
                            print *, 'count_output reaches max_count: no more warning will be shown'
                            count_output = count_output+1
                         end if
                      end if
                   end if
                   xch%Tg = min(max(xch%Tg,MP_Tmin),MP_Tmax)
                   do ichem = 0, 6 -1
                      chemary3(ichem+5)%y(i,j,k) = xch%ychem(ichem)
                   enddo
                   yco(i,j,k) = xch%yco
                   xmu = get_xmu(xch%ychem) 
                   p(i,j,k) = (cgs_kb*xch%Tg)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e 
                   Tdust(i,j,k) = xch%Td
                end do
             end do
          end do
       end do
    end do
  contains
    subroutine GetKH2pdThin(x,y,z,kh2pd)
      real(kind=8),dimension(:),pointer,intent(IN) :: x, y, z
      real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax),intent(OUT) :: kh2pd
      real(kind=8),parameter :: pdis = 0.15
      integer :: isrc,i,j,k
      real(kind=8) :: r2,lum_fuv,fuvflux
      kh2pd(:,:,:) = 0.d0
      do isrc=0, rs_info%nsource-1
         lum_fuv = rs_info%x_fuv(isrc) * rs_info%lum(isrc)
         do i = Imin, Imax
            do j = Jmin, Jmax
               do k = Kmin, Kmax
                  r2 = (x(i) - rs_info%spos(0,isrc))**2 + &
                       (y(j) - rs_info%spos(1,isrc))**2 + &
                       (z(k) - rs_info%spos(2,isrc))**2
                  fuvflux = lum_fuv / (4.*Pi*r2*Unit_l**2)
                  kh2pd(i,j,k)=kh2pd(i,j,k)+pdis*3.4d-10*fuvflux/1.21d+7 
               end do
            end do
         end do
      end do
    end subroutine GetKH2pdThin
    subroutine GetKHmpdThin(x,y,z,khmpd)
      real(kind=8),dimension(:),pointer,intent(IN) :: x, y, z
      real(kind=8),dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax),intent(OUT) :: khmpd
      integer :: isrc,i,j,k
      real(kind=8) :: r2,r_star
      khmpd(:,:,:)=0.d0
      do isrc=0, rs_info%nsource-1
         r_star = sqrt(rs_info%lum(isrc)/(4.*Pi * cgs_sigma*rs_info%Trad(isrc)**4)) 
         do i = Imin, Imax
            do j = Jmin, Jmax
               do k = Kmin, Kmax
                  r2 = (x(i) - rs_info%spos(0,isrc))**2 + &
                       (y(j) - rs_info%spos(1,isrc))**2 + &
                       (z(k) - rs_info%spos(2,isrc))**2
                  khmpd(i,j,k)=khmpd(i,j,k)+rs_info%hhm(isrc)*r_star**2/(r2*Unit_l**2)
               end do
            end do
         end do
      end do
    end subroutine GetKHmpdThin
  end subroutine ch_PrimChemistry
  subroutine ch_CellChemCool(xch)
    type(chem_vals) :: xch
    real(kind=8),parameter :: eps_imp = 1d-1
    real(kind=8) :: T_o, y_o(0:6 -1), t_chemcool, t_cool, yco_o
    integer :: myrank 
    T_o = xch%Tg
    y_o(:) = xch%ychem(:)
    yco_o = xch%yco
    call CoolSolverExplicit(xch,t_chemcool,t_cool)
    if (xch%dt > eps_imp * t_cool) then
       xch%Tg = T_o
       xch%ychem(:) = y_o
       xch%yco = yco_o
       call CoolSolverImplicit(xch) 
       num_implicit = num_implicit + 1 
    end if
    t_chemcool_local = MIN(t_chemcool_local,t_chemcool)
  end subroutine ch_CellChemCool
  subroutine ch_EvalRadForce(xnH,yHn,fx_hpi,fy_hpi,fz_hpi)
    real(kind=8),intent(IN) :: xnH,yHn
    real(kind=8),intent(INOUT) :: fx_hpi,fy_hpi,fz_hpi
    fx_hpi = fx_hpi * xnH * yHn
    fy_hpi = fy_hpi * xnH * yHn
    fz_hpi = fz_hpi * xnH * yHn
  end subroutine ch_EvalRadForce
  subroutine ch_check_all_cells
    real(kind=8),dimension(:,:,:),pointer :: rho, p, &
         khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi
    real(kind=8),dimension(:,:,:),pointer :: yco
    real(kind=8),dimension(0:6 -1) :: ychem
    real(kind=8) :: yco_l
    integer :: level, n, gid, i, j, k, ic
    type chemsp
      real(kind=8),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(5:10)
    integer :: ichem
    do level = Lmin, Lmax
       do n = Gidmin, GidListMax( level )
          gid = GidList( n, level )
          do ichem = 5, 10
            chemary3(ichem)%y => get_Ucomp(ichem,gid)
          enddo
          yco => get_Ucomp(11, gid)
          do i = Imin, Imax
             do j = Jmin, Jmax
                do k = Kmin, Kmax
                   do ichem = 0, 6 -1
                      ychem(ichem) = chemary3(ichem+5)%y(i,j,k)
                   enddo
                   yco_l = yco(i,j,k)
                   call adjust_abundance(ychem &
                      , yco_l &
                      , MP_Metallicity)
                   do ichem = 0, 6 -1
                      chemary3(ichem+5)%y(i,j,k) = ychem(ichem)
                   enddo
                   yco(i,j,k) = yco_l
                end do
             end do
          end do
       end do
    end do
  end subroutine ch_check_all_cells
  subroutine ch_EvalVelocity(vx, vy, vz,fx_hpi,fy_hpi,fz_hpi, dt)
        implicit none
        real(kind=8),intent(IN) :: dt
        real(kind=8),intent(IN) :: fx_hpi,fy_hpi,fz_hpi
        real(kind=8),intent(INOUT) :: vx, vy, vz
        real(kind=8),parameter :: V_UPPER_LIMIT = 1.d2 
        real(kind=8) :: v2
        vx = vx + fx_hpi * dt / Unit_v 
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
