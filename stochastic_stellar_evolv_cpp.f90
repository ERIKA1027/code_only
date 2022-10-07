module stochastic_star
  use string, only : concat
  use io_util, only : read_env
  use mt19937
  implicit none
  private
  integer, save :: lnp_tot 
  type t_sp_mass
    integer :: sid 
    integer :: id_ini 
    real(kind=8) :: mass 
    real(kind=8) :: ftime 
    type (t_sp_mass), pointer :: next => null()
  end type t_sp_mass
  type t_radsc
    integer :: np 
    integer :: lnp 
    integer :: nstar 
    real(kind=8) :: mtot 
    type (t_sp_mass), pointer :: sinfo => null()
    type (t_radsc), pointer :: next => null()
  end type t_radsc
  type(t_radsc), save, pointer :: Sparticle => null()
  real(kind=8), save, allocatable, dimension(:) :: dismass, disf
  integer, save :: idis
  type stellar_evol
    integer :: nmodel
    real(kind=8), pointer, dimension(:) :: Age 
    real(kind=8), pointer, dimension(:) :: Rstar 
    real(kind=8), pointer, dimension(:) :: Lumi 
    real(kind=8), pointer, dimension(:) :: Teff 
    real(kind=8), pointer, dimension(:) :: xeuv 
    real(kind=8), pointer, dimension(:) :: xfuv 
    real(kind=8), pointer, dimension(:) :: alpha_euv 
    real(kind=8), pointer, dimension(:) :: heat_euv 
    real(kind=8), pointer, dimension(:) :: lumeuv 
    real(kind=8), pointer, dimension(:) :: lumfuv 
    real(kind=8), pointer, dimension(:) :: sigd_euv 
    real(kind=8), pointer, dimension(:) :: sigd_fuv 
    real(kind=8), pointer, dimension(:) :: rOII 
  end type stellar_evol
  type (stellar_evol), allocatable, dimension(:), save :: StarTable
  integer, save :: nstarmdl
  type st_rs_info
    integer :: nstar
    real(kind=8) :: lum
    real(kind=8) :: Teff
    real(kind=8) :: xeuv
    real(kind=8) :: xfuv
    real(kind=8) :: alpha_euv
    real(kind=8) :: heat_euv
    real(kind=8) :: lumeuv
    real(kind=8) :: lumfuv
    real(kind=8) :: sigd_euv
    real(kind=8) :: sigd_fuv
    real(kind=8) :: rOII
  end type st_rs_info
  real(kind=8), parameter :: cgs_lsun = 3.904d33 
  public :: set_radsc_stchast, st_rs_info
contains
  subroutine set_radsc_stchast(irank, np_sink, lnp, smass, smdot, stime, rs_info)
    implicit none
    integer, intent(IN) :: irank, np_sink, lnp
    real(kind=8), intent(IN) :: smass, smdot, stime
    integer, save :: ifirst = 1
    type(st_rs_info) :: rs_info
    if (ifirst == 1) then
      call radsc_init(irank, stime)
      ifirst = 0
    endif
    call update_sink_mass(irank, np_sink, lnp, smass, smdot, stime)
    call get_rsinfo(irank, np_sink, lnp, smass, smdot, stime, rs_info)
  end subroutine set_radsc_stchast
  subroutine get_rsinfo(irank, np_sink, lnp, smass, smdot, stime, rs_info)
    implicit none
    integer, intent(IN) :: irank, np_sink, lnp
    real(kind=8), intent(IN) :: smass, smdot, stime
    type(st_rs_info) :: rs_info, rsinfo_star
    type(t_radsc), pointer :: sptr
    type (t_sp_mass), pointer :: sinf
    integer :: isp, nstar, n, id_state
    real(kind=8) :: lumi_tot, sion_tot, sfuv_tot, alpSi, Sion, heatSi, lum_e, lum_f &
      ,lsig_e, lsig_f, alpSirOII, inv_lumtot, inv_siontot
    sptr => Sparticle
    isp = 0
    do
      if(.not. associated(sptr)) then
        print *, "error occurs at update_sink_mass"
        stop
      endif
      isp = isp + 1
      if(isp == lnp) then
        nstar = sptr%nstar
        sinf => sptr%sinfo
        lumi_tot = 0.d0
        sion_tot = 0.d0
        sfuv_tot = 0.d0
        alpSi = 0.d0
        heatSi = 0.d0
        lum_e = 0.d0
        lum_f = 0.d0
        lsig_e = 0.d0
        lsig_f = 0.d0
        alpSirOII = 0.d0
        do n=1, nstar-1 
          call get_invidual_star_inf(sinf, rsinfo_star, stime, id_state)
          if(id_state .ne. 1) then 
            lumi_tot = lumi_tot + rsinfo_star%lum
            Sion = rsinfo_star%xeuv*rsinfo_star%lum
            sion_tot = sion_tot + Sion
            sfuv_tot = sfuv_tot + rsinfo_star%xfuv*rsinfo_star%lum
            alpSi = alpSi + Sion*rsinfo_star%alpha_euv
            heatSi = heatSi + Sion*rsinfo_star%heat_euv
            lum_e = lum_e + rsinfo_star%lumeuv*rsinfo_star%lum
            lum_f = lum_f + rsinfo_star%lumfuv*rsinfo_star%lum
            lsig_e = lsig_e + rsinfo_star%sigd_euv*rsinfo_star%lum
            lsig_f = lsig_f + rsinfo_star%sigd_fuv*rsinfo_star%lum
            alpSirOII = alpSirOII+ Sion*rsinfo_star%alpha_euv*rsinfo_star%rOII
          endif
          sinf => sinf%next
        enddo
        inv_lumtot = 1.d0 / max(lumi_tot, 1.d-30)
        inv_siontot= 1.d0 / max(sion_tot, 1.d-30)
        rs_info%nstar = nstar
        rs_info%lum = lumi_tot
        rs_info%Teff = 1.d4 
        rs_info%xeuv = sion_tot * inv_lumtot
        rs_info%xfuv = sfuv_tot * inv_lumtot
        rs_info%alpha_euv = alpSi * inv_siontot
        rs_info%heat_euv = heatSi * inv_siontot
        rs_info%lumeuv = lum_e * inv_lumtot
        rs_info%lumfuv = lum_f * inv_lumtot
        rs_info%sigd_euv = lsig_e * inv_lumtot
        rs_info%sigd_fuv = lsig_f * inv_lumtot
        rs_info%rOII = alpSirOII / max(alpSi , 1.d-30)
        exit
      endif
      sptr => sptr%next
    enddo
  end subroutine get_rsInfo
  subroutine get_invidual_star_inf(sinf, rsinfo_star, stime, id_state)
    implicit none
    real(kind=8), intent(IN) :: stime 
    type (t_sp_mass), pointer :: sinf
    type(st_rs_info) :: rsinfo_star
    integer, intent(OUT) :: id_state 
    integer :: id_ini, i
    real(kind=8) :: lifetime, age_star, age_hi
    id_state = 0
    id_ini = sinf%id_ini
    lifetime = StarTable(id_ini)%Age(StarTable(id_ini)%nmodel)
    age_star = stime-sinf%ftime
    if(age_star > lifetime) then
      id_state = 1 
      rsinfo_star%lum = 0.d0
      rsinfo_star%Teff = 0.d0
      rsinfo_star%xeuv = 0.d0
      rsinfo_star%xfuv = 0.d0
      rsinfo_star%alpha_euv = 0.d0
      rsinfo_star%heat_euv = 0.d0
      rsinfo_star%lumeuv = 0.d0
      rsinfo_star%lumfuv = 0.d0
      rsinfo_star%sigd_euv = 0.d0
      rsinfo_star%sigd_fuv = 0.d0
      rsinfo_star%rOII = 0.d0
      return
    endif
    if(age_star < StarTable(id_ini)%Age(1)) then
      rsinfo_star%lum = StarTable(id_ini)%Lumi(1)
      rsinfo_star%Teff = StarTable(id_ini)%Teff(1)
      rsinfo_star%xeuv = StarTable(id_ini)%xeuv(1)
      rsinfo_star%xfuv = StarTable(id_ini)%xfuv(1)
      rsinfo_star%alpha_euv = StarTable(id_ini)%alpha_euv(1)
      rsinfo_star%heat_euv = StarTable(id_ini)%heat_euv(1)
      rsinfo_star%lumeuv = StarTable(id_ini)%lumeuv(1)
      rsinfo_star%lumfuv = StarTable(id_ini)%lumfuv(1)
      rsinfo_star%sigd_euv = StarTable(id_ini)%sigd_euv(1)
      rsinfo_star%sigd_fuv = StarTable(id_ini)%sigd_fuv(1)
      rsinfo_star%rOII = StarTable(id_ini)%rOII(1)
      return
    endif
    do i=2, StarTable(id_ini)%nmodel
      if(age_star < StarTable(id_ini)%Age(i)) then
        age_hi=(age_star-StarTable(id_ini)%Age(i-1))/(StarTable(id_ini)%Age(i)-StarTable(id_ini)%Age(i-1))
        rsinfo_star%lum =(StarTable(id_ini)%Lumi(i)-StarTable(id_ini)%Lumi(i-1))*age_hi+StarTable(id_ini)%Lumi(i-1)
        rsinfo_star%Teff=(StarTable(id_ini)%Teff(i)-StarTable(id_ini)%Teff(i-1))*age_hi+StarTable(id_ini)%Teff(i-1)
        rsinfo_star%xeuv=(StarTable(id_ini)%xeuv(i)-StarTable(id_ini)%xeuv(i-1))*age_hi+StarTable(id_ini)%xeuv(i-1)
        rsinfo_star%xfuv=(StarTable(id_ini)%xfuv(i)-StarTable(id_ini)%xfuv(i-1))*age_hi+StarTable(id_ini)%xfuv(i-1)
        rsinfo_star%alpha_euv=(StarTable(id_ini)%alpha_euv(i)-StarTable(id_ini)%alpha_euv(i-1)) &
                                *age_hi+StarTable(id_ini)%alpha_euv(i-1)
        rsinfo_star%heat_euv=(StarTable(id_ini)%heat_euv(i)-StarTable(id_ini)%heat_euv(i-1)) &
                                *age_hi+StarTable(id_ini)%heat_euv(i-1)
        rsinfo_star%lumeuv=(StarTable(id_ini)%lumeuv(i)-StarTable(id_ini)%lumeuv(i-1))*age_hi+StarTable(id_ini)%lumeuv(i-1)
        rsinfo_star%lumfuv=(StarTable(id_ini)%lumfuv(i)-StarTable(id_ini)%lumfuv(i-1))*age_hi+StarTable(id_ini)%lumfuv(i-1)
        rsinfo_star%sigd_euv=(StarTable(id_ini)%sigd_euv(i)-StarTable(id_ini)%sigd_euv(i-1)) &
                                *age_hi+StarTable(id_ini)%sigd_euv(i-1)
        rsinfo_star%sigd_fuv=(StarTable(id_ini)%sigd_fuv(i)-StarTable(id_ini)%sigd_fuv(i-1)) &
                                *age_hi+StarTable(id_ini)%sigd_fuv(i-1)
        rsinfo_star%rOII=(StarTable(id_ini)%rOII(i)-StarTable(id_ini)%rOII(i-1))*age_hi+StarTable(id_ini)%rOII(i-1)
        exit
      endif
    enddo
  end subroutine get_invidual_star_inf
  subroutine update_sink_mass(irank, np_sink, lnp, smass, smdot, stime)
    implicit none
    integer, intent(IN) :: irank, np_sink, lnp
    real(kind=8), intent(IN) :: smass, smdot, stime
    integer :: isp
    type(t_radsc), pointer :: sptr
    if(lnp > lnp_tot) then
      lnp_tot = lnp_tot + 1
      call creat_newSparticle(np_sink, lnp)
      if(lnp .ne. lnp_tot) then
        print *, "error occurs at update sink information at update_sink_mass"
        stop
      endif
    endif
    sptr => Sparticle
    isp = 0
    do
      if(.not. associated(sptr)) then
        print *, "error occurs at update_sink_mass"
        stop
      endif
      isp = isp + 1
      if(isp == lnp) then
        call cal_sinkmass_with_dist(irank, np_sink, lnp, smass, smdot, stime, sptr)
        exit
      endif
      sptr => sptr%next
    enddo
  end subroutine update_sink_mass
  subroutine cal_sinkmass_with_dist(irank, np_sink, lnp, smass, smdot, stime, sptr)
    implicit none
    integer, intent(IN) :: irank, np_sink, lnp
    real(kind=8), intent(IN) :: smass, smdot, stime
    type(t_radsc), pointer :: sptr
    real(kind=8) :: newsmass
    integer :: id_ini
    if(sptr%mtot > smass) then 
      return
    endif
    do
      call get_newmass(newsmass, id_ini)
      sptr%mtot = sptr%mtot + newsmass
      sptr%nstar = sptr%nstar + 1
      call add_newstar(sptr, sptr%nstar, newsmass, stime, id_ini)
      call add_log_newstar(irank, np_sink, lnp, sptr%nstar, newsmass, stime, id_ini)
      if(sptr%mtot > smass) then
        exit
      endif
    enddo
  end subroutine cal_sinkmass_with_dist
  subroutine add_log_newstar(irank, np, lnp, sid, smass, stime, id_ini)
    implicit none
    integer, intent(IN) :: irank, np, lnp, sid, id_ini
    real(kind=8), intent(IN) :: smass, stime
    integer :: FH=11
    character(len=128) :: file, filename, dir
    call read_env('DIR', dir)
    write (filename, '("radsc_log_", i0, ".dat")') irank 
    file=concat(dir,filename)
    open(FH, file=file, position='APPEND')
    write(FH, '(I14, I14, I14, 2(1PE17.9), I14)') np, lnp, sid, smass, stime, id_ini
    close(FH)
  end subroutine add_log_newstar
  subroutine get_newmass(newsmass, id_ini)
    implicit none
    real(kind=8) :: r, newsmass
    integer :: i, id_ini
    r = grnd()
    do i=1,idis
      if(r .le. disf(i)) then
        newsmass = dismass(i)
        id_ini = i
        return
      endif
    enddo
    print *, "error occurs at get_newmass"
    stop
  end subroutine get_newmass
  subroutine radsc_init(irank, stime)
    implicit none
    integer, intent(IN) :: irank
    real(kind=8), intent(IN) :: stime
    character(len=128) :: filename, file, file2, dir
    real(kind=8) :: ini_mass, ftime
    integer :: ios, np, lnp, lnp_max, line_num, iline, sid, i, nstar, dnum, id_ini
    integer, allocatable, dimension(:) :: np_line, lnp_line, np_label, sid_line, idini_line
    real(kind=8), allocatable, dimension(:) :: mass_line, time_line
    real(kind=8) :: tot_smass, mass, dfunc, newsmass
    integer :: FH = 11
    logical :: exist
    type(t_radsc), pointer :: sptr
    type (t_sp_mass), pointer :: sinf
    call sgrnd(irank+1)
    call read_stellar_model
    ios = 1
    dnum= 0
    call read_env('DIR', dir)
    write(filename, '("Stellar_Fit/stellar_bin.dat")') 
    file=concat(dir,filename)
    open(FH, file=file)
    do
      read(FH, fmt=*, iostat=ios) mass, dfunc
      if(ios < 0) then
        exit
      endif
      dnum = dnum + 1
    enddo
    close(FH)
    idis = dnum
    allocate(dismass(dnum), disf(dnum))
    ios = 1
    dnum= 0
    open(FH, file=file)
    do
      read(FH, fmt=*, iostat=ios) mass, dfunc
      if(ios < 0) then
        exit
      endif
      dnum = dnum + 1
      dismass(dnum) = mass
      disf(dnum) = dfunc
    enddo
    disf(dnum) = 1.d0
    close(FH)
    write (filename, '("radsc_log_", i0, ".dat")') irank 
    file2=concat(dir,filename)
    inquire(file=file2, exist=exist)
    if(.not. exist) then
      lnp_tot = 0
      return
    endif
    ios = 1
    lnp_max = 0
    line_num = 0
    open(FH, file=file2)
    do while ( .true. )
      read(FH, fmt=*, iostat=ios) np, lnp, sid, ini_mass, ftime, id_ini
      if(ios < 0) then
        exit
      endif
      if (stime >= ftime ) then
        if(lnp_max < lnp) lnp_max = lnp
        line_num = line_num + 1
        call get_newmass(newsmass, id_ini) 
      endif
    enddo
    close(FH)
    lnp_tot = lnp_max
    allocate(np_line(line_num), lnp_line(line_num), mass_line(line_num) &
      , time_line(line_num), sid_line(line_num), idini_line(line_num))
    allocate(np_label(lnp_tot))
    ios = 1
    iline = 0
    open(FH, file=file2)
    do while (.true.)
      read(FH, fmt=*, iostat=ios) np, lnp, sid, ini_mass, ftime, id_ini
      if(ios < 0) then
        exit
      endif
      if (stime >= ftime ) then
        iline = iline + 1
        np_line(iline) = np
        lnp_line(iline) = lnp
        sid_line(iline) = sid
        idini_line(iline) = id_ini
        mass_line(iline) = ini_mass
        time_line(iline) = ftime
        np_label(lnp_line(iline)) = np_line(iline) 
      endif
    enddo
    close(FH)
    do lnp=1, lnp_tot
      call creat_newSparticle(np_label(lnp), lnp)
    enddo
    sptr => Sparticle
    do lnp = 1, lnp_tot
      tot_smass = 0.d0 
      nstar = 0 
      do iline = 1, line_num
        if(lnp_line(iline) == lnp) then
          call add_newstar(sptr, sid_line(iline), mass_line(iline), time_line(iline), idini_line(iline))
          tot_smass = tot_smass + mass_line(iline)
          nstar = nstar + 1
        endif
      enddo
      sptr%nstar = nstar 
      sptr%mtot = tot_smass 
      sptr => sptr%next
    enddo
    deallocate(np_line, lnp_line, mass_line, time_line, sid_line,idini_line)
    deallocate(np_label)
    open(FH, file=file2, status='replace')
    sptr => Sparticle
    do
      if(.not. associated(sptr)) then
        exit
      endif
      sinf => sptr%sinfo
      do
        if(.not. associated(sinf)) then
          exit
        endif
        write(FH, '(I14, I14, I14, 2(1PE17.9), I14)') sptr%np, sptr%lnp, sinf%sid, sinf%mass, sinf%ftime, sinf%id_ini
        sinf => sinf%next
      enddo
      sptr => sptr%next
    enddo
    close(FH)
  end subroutine radsc_init
  subroutine creat_newSparticle(np, lnp)
    implicit none
    integer, intent(IN) :: np, lnp
    type(t_radsc), pointer :: newSparticle
    type(t_radsc), pointer :: sptr
    allocate(newSparticle)
    newSparticle%np = np
    newSparticle%lnp = lnp
    newSparticle%nstar = 0
    newSparticle%mtot = 0.d0
    if(.not. associated(Sparticle)) then
      Sparticle => newSparticle
      return
    endif
    sptr => Sparticle
    do
      if(.not. associated(sptr%next)) then
        sptr%next => newSparticle
        exit
      endif
      sptr => sptr%next
    enddo
  end subroutine creat_newSparticle
  subroutine add_newstar(sptr, sid, smass, ftime, id_ini)
    implicit none
    type (t_radsc), pointer :: sptr
    integer, intent(IN) :: sid, id_ini
    real(kind=8), intent(IN) :: smass, ftime
    type (t_sp_mass), pointer :: newsinf, sinf
    allocate(newsinf)
    newsinf%sid = sid
    newsinf%id_ini= id_ini
    newsinf%mass = smass
    newsinf%ftime = ftime
    if(.not. associated(sptr%sinfo)) then
      sptr%sinfo => newsinf
      return
    endif
    sinf => sptr%sinfo
    do
      if(.not. associated(sinf%next)) then
        sinf%next => newsinf
        exit
      endif
      sinf => sinf%next
    enddo
  end subroutine add_newstar
  subroutine output_data
    implicit none
    type(t_radsc), pointer :: sptr
    type (t_sp_mass), pointer :: sinfo
    sptr => Sparticle
    print *, "---------------- output sp ---------------------"
    print *, "lnp_tot:", lnp_tot
    do
      if(.not. associated(sptr)) then
        exit
      endif
      print *, "np, lnp, nstar, mass:", sptr%np, sptr%lnp, sptr%nstar, sptr%mtot
      sinfo => sptr%sinfo
      do
        if(.not. associated(sinfo)) then
          exit
        endif
        print *, "sid, mass, ftime", sinfo%sid, sinfo%mass, sinfo%ftime
        sinfo => sinfo%next
      enddo
      sptr => sptr%next
    enddo
  end subroutine output_data
  subroutine read_stellar_model
    implicit none
    character(len=128) :: filename, dummy, dir, file
    real(kind=8) :: mass
    real(kind=8) :: Age, Rstar, Lumi, Teff, xeuv, xfuv, alpha_euv, heat_euv, lumeuv &
      , lumfuv, sigd_euv, sigd_fuv, rOII
    integer :: FH = 11
    integer :: i, n, nmodel
    call read_env('DIR', dir)
    filename = trim("Stellar_Fit/tsevolv_0.0.dat")
    file=concat(dir,filename)
    open(FH, file=file)
    read(FH, fmt=*) nstarmdl
    allocate(StarTable(nstarmdl))
    do i=1, nstarmdl
      read(FH, fmt=*) dummy
      read(FH, fmt=*) mass, nmodel
      StarTable(i)%nmodel = nmodel
      allocate(StarTable(i)%Age(nmodel) &
              ,StarTable(i)%Rstar(nmodel) &
              ,StarTable(i)%Lumi(nmodel) &
              ,StarTable(i)%Teff(nmodel) &
              ,StarTable(i)%xeuv(nmodel) &
              ,StarTable(i)%xfuv(nmodel) &
              ,StarTable(i)%alpha_euv(nmodel) &
              ,StarTable(i)%heat_euv(nmodel) &
              ,StarTable(i)%lumeuv(nmodel) &
              ,StarTable(i)%lumfuv(nmodel) &
              ,StarTable(i)%sigd_euv(nmodel) &
              ,StarTable(i)%sigd_fuv(nmodel) &
              ,StarTable(i)%rOII(nmodel))
      do n=1, nmodel
        read(FH, fmt=*) Age, Rstar, Lumi, Teff, xeuv, xfuv, alpha_euv, heat_euv, lumeuv &
          , lumfuv, sigd_euv, sigd_fuv, rOII
        StarTable(i)%Age(n) = Age
        StarTable(i)%Rstar(n) = Rstar
        StarTable(i)%Lumi(n) = Lumi*cgs_lsun 
        StarTable(i)%Teff(n) = Teff
        StarTable(i)%xeuv(n) = xeuv
        StarTable(i)%xfuv(n) = xfuv
        StarTable(i)%alpha_euv(n) = alpha_euv
        StarTable(i)%heat_euv(n) = heat_euv
        StarTable(i)%lumeuv(n) = lumeuv
        StarTable(i)%lumfuv(n) = lumfuv
        StarTable(i)%sigd_euv(n) = sigd_euv
        StarTable(i)%sigd_fuv(n) = sigd_fuv
        StarTable(i)%rOII(n) = rOII
      enddo
    enddo
    close(FH)
  end subroutine read_stellar_model
end module stochastic_star
