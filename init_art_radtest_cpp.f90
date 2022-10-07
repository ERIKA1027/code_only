program main
  use io
  use grid
  use outputdata 
  implicit none
  call init
  call make_basegrid
  call make_physvar
  call dumpdata
  call output_data 
  call finalize
end program main
subroutine init
  use mpilib
  use grid
  use parameter
  use unit
  use kinzoku
  use modelParameter
  implicit none
  call mpilib_init
  call grid_init
  call parameter_init
  call unit_init
  call init_kinzoku
  call modelParameter_init
  call eos_init
end subroutine init
subroutine finalize
  use mpilib
  use grid
  call grid_finalize
  call mpi_finalize(ierr)
end subroutine finalize
subroutine make_basegrid
  use grid
  use mpilib
  use refine, only : refine_get_norder
  implicit none
  integer :: n, i, j, k, in, jn, kn, gcount, gid, rankid, rank, level, gidn, rankn, sz, szb
  integer :: ngrid_total = (8) * (8) * (8)
  integer,dimension(0:400 -1) :: ngrid_node
  integer,dimension(:,:),allocatable :: ig, jg, kg
  integer :: norder
  integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
  integer :: ifc, jfc, kfc
  LevelMax = 0 
  norder = refine_get_norder(0, 0, 0, 8 -1, 8 -1, 8 -1)
  ngrid_node(:) = int(ngrid_total/(400))
  do n = 1, mod( ngrid_total, 400 )
     ngrid_node(400 -n) = ngrid_node(400 -n) + 1
  enddo
  level = 0 
  rankid = 0 
  myrank = get_myrank() 
  gcount = 0 
  gid = -1 
  ifc = 0 
 jfc = 0 
 kfc = 0 
  call assign()
  call fillcells(norder, R,U,L,D,B,F)
  call update_gidlist(level)
  do n = lbound(GidList,1), GidListMax( level ) 
     call alloc_U1order(GidList(n, level))
     call alloc_U2order(GidList(n, level))
  enddo
  call mpi_allreduce( GidListMax( level )+1, sz, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  allocate( &
       ig(Gidmin:sz-1,0:400 -1), &
       jg(Gidmin:sz-1,0:400 -1), &
       kg(Gidmin:sz-1,0:400 -1) )
  ig = -1
  jg = -1
  kg = -1
  call mpi_allgather( Igrid, sz, MPI_INTEGER, ig, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  call mpi_allgather( Jgrid, sz, MPI_INTEGER, jg, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  call mpi_allgather( Kgrid, sz, MPI_INTEGER, kg, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  do k = 0, (8) -1
     do j = 0, (8) -1
        do i = 0, (8) -1
           do kn = 0, (8) -1
              do jn = 0, (8) -1
                 do in = 0, (8) -1
                    gidn = GidBase(in,jn,kn)
                    rankn = RankBase(in,jn,kn)
                    rank = RankBase(i,j,k)
                    gid = GidBase(i,j,k)
                    if ( ig(gid,rank) - 1 == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rank&
&n) ) then
                       NeighborGid(0, 0, gid, rank) = gidn
                       NeighborRank(0, 0, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) + 1 == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rank&
&n) ) then
                       NeighborGid(1, 0, gid, rank) = gidn
                       NeighborRank(1, 0, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) - 1 == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rank&
&n) ) then
                       NeighborGid(0, 1, gid, rank) = gidn
                       NeighborRank(0, 1, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) + 1 == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rank&
&n) ) then
                       NeighborGid(1, 1, gid, rank) = gidn
                       NeighborRank(1, 1, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) - 1 == kg(gidn,rank&
&n) ) then
                       NeighborGid(0, 2, gid, rank) = gidn
                       NeighborRank(0, 2, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) + 1 == kg(gidn,rank&
&n) ) then
                       NeighborGid(1, 2, gid, rank) = gidn
                       NeighborRank(1, 2, gid, rank) = rankn
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate( ig, jg, kg )
contains
  recursive subroutine fillcells(n, right, up, left, down, back, foward)
    integer,intent(IN) :: n, right, up, left, down, back, foward
    if (n == 0) return
    call fillcells(n-1,foward,right,back,left,down,up)
    call connect(up)
    call fillcells(n-1,up,foward,down,back,left,right)
    call connect(right)
    call fillcells(n-1,up,foward,down,back,left,right)
    call connect(down)
    call fillcells(n-1,left,down,right,up,back,foward)
    call connect(foward)
    call fillcells(n-1,left,down,right,up,back,foward)
    call connect(up)
    call fillcells(n-1,up,back,down,foward,right,left)
    call connect(left)
    call fillcells(n-1,up,back,down,foward,right,left)
    call connect(down)
    call fillcells(n-1,back,right,foward,left,up,down)
  end subroutine fillcells
  subroutine connect(dir)
    integer,intent(IN) :: dir
    select case (dir)
    case (R)
       ifc = ifc + 1
    case (L)
       ifc = ifc - 1
    case (F)
       jfc = jfc + 1
    case (B)
       jfc = jfc - 1
    case (U)
       kfc = kfc + 1
    case (D)
       kfc = kfc - 1
    end select
    call assign()
  end subroutine connect
  subroutine assign()
    if ( ifc < 0 .or. ifc > 8 -1 .or. &
         jfc < 0 .or. jfc > 8 -1 .or. &
         kfc < 0 .or. kfc > 8 -1 ) return
    gcount = gcount + 1
    gid = gid + 1
    if ( gcount > ngrid_node(rankid) ) then 
       rankid = rankid + 1
       gid = 0
       gcount = 1
    endif
    GidBase(ifc,jfc,kfc) = gid
    RankBase(ifc,jfc,kfc) = rankid
    if ( rankid == myrank ) then
       Levels(gid) = level
       call alloc_block_by_gid(gid)
       Igrid(gid) = ifc
       Jgrid(gid) = jfc
       Kgrid(gid) = kfc
    endif
  end subroutine assign
end subroutine make_basegrid
subroutine make_physvar
  use mpilib
  use grid
  use refine
  use grid_boundary
  use fg2cg
  use parameter
  use io_util, only : readenv
  use unit
  use modelParameter
  use primordial 
  use kinzoku, only : fdust_solar
  implicit none
  real(kind=8),dimension(:,:,:),pointer :: rho, vx, vy, vz, p, psi, bx, by, bz, db, gx, gy, gz
  real(kind=8),dimension(:,:,:),pointer :: &
       yhn, yh2, yel, yhp, yhm, yh2p,khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi,yco,rdpco, ychem 
  real(kind=8),dimension(:,:,:),pointer :: kh2pd 
  real(kind=8),dimension(:,:,:),pointer :: kgfuv, kdph, krOII
  real(kind=8),dimension(:,:,:),pointer :: Nrad, Frad_x, Frad_y, Frad_z
  real(kind=8),dimension(:),pointer :: x, y, z
  real(kind=8) :: MachNumber
  integer :: level,gid
  integer :: i,j,k, n
  real(kind=8) :: h,ic,jc,kc,ic0,jc0,kc0,ycoi
  real(kind=8) :: r2, x0, y0, z0, a, rho0, p0, vx0, vy0, vz0, psi0, h1d, db0, b0, rboundary, cs0
  logical :: bool
  real(kind=8) :: yhni, yh2i, yeli, yhpi, yhmi, yh2pi
  real(kind=8) :: yhni2, yh2i2, yeli2, yhpi2, yhmi2, yh2pi2, rho02, cs02, radius, hmax
  real(kind=8) :: urad_cmb, nrad_cmb, nrad_fuv_bg, urad_fuv_bg
  integer :: ncell_HI, ncell_height
  integer :: ichem
  if (MP_GasType == 0) then 
     yh2i=1.d-8
     yhpi=1.d-8
     yhmi=1.d-20
     yh2pi=1.d-20
     yhni=1.d0 - (yhpi+yhmi) - 2*(yh2i+yh2pi)
     yeli=yhpi - yhmi + yh2pi
     ycoi = 0.927d-4
  else if (MP_GasType == 1 .or. MP_GasType == 2) then 
     yhni=1.d-8
     yhpi=1.d-8
     yhmi=1.d-20
     yh2pi=1.d-20
     yh2i= (1.d0 - (yhni+yhpi+yhmi) - 2*yh2pi)/2.d0
     yeli=yhpi - yhmi + yh2pi
     ycoi = 0.927d-4
  else
       if (get_myrank() == 0) &
       print '(A,/,A)', "MP_GasType should be either 0, 1 or 2","stopping..."
       stop
  end if
  rho0 = MP_N0 * cgs_amu * MP_mu / Unit_rho 
  cs0 = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
       /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))/ Unit_v 
  rboundary = MP_Boxsize
  call writeparam
  level = 0 
  Time(:) = 0.d0
  Step(:) = 0.d0
  Dtime(:) = 0.d0
  Dstep(:) = 0.d0
  ic0 = ((8) * (8)-1)/2.d0 
  jc0 = ((8) * (8)-1)/2.d0
  kc0 = ((8) * (8)-1)/2.d0
  h = 2*rboundary/min((8) * (8), (8) * (8), (8) * (8) )
  CellWidth(:, level) = (/h, h, h/)
  do n = Gidmin, GidListMax( level )
     gid = GidList( n, level )
     ic = ic0 - Igrid(gid)* 8
     jc = jc0 - Jgrid(gid)* 8
     kc = kc0 - kgrid(gid)* 8
     x => get_Xp(gid)
     y => get_Yp(gid)
     z => get_Zp(gid)
     do i=Imingh, Imaxgh
        x(i) = (i-ic)*h
     enddo
     do j=Jmingh, Jmaxgh
        y(j) = (j-jc)*h
     enddo
     do k=Kmingh, Kmaxgh
        z(k) = (k-kc)*h
     enddo
  enddo
  do level = Lmin, Lmax
     if ( level > Lmin ) then
        call refineLevel(level, bool)
        LevelMax = max(LevelMax, 0)
        if ( .not. bool ) exit
     endif
     do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )
        psi => get_Ucomp(17,gid)
        gx => get_Ucomp(18,gid)
        gy => get_Ucomp(19,gid)
        gz => get_Ucomp(20,gid)
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)
        psi(:,:,:) = 0.d0
        gx(:,:,:) = 0.d0
        gy(:,:,:) = 0.d0
        gz(:,:,:) = 0.d0
        vx => get_Ucomp(1,gid)
        vy => get_Ucomp(2,gid)
        vz => get_Ucomp(3,gid)
        vx(:,:,:) = 2.d6 / Unit_v
        vy(:,:,:) = 0.d0
        vz(:,:,:) = 0.d0
        do ichem = 5, 10
          ychem => get_Ucomp(ichem,gid)
          ychem(:,:,:) = 1.d-15
        enddo
        yhn => get_Ucomp(5,gid)
        yh2 => get_Ucomp(6,gid)
        yel => get_Ucomp(7,gid)
        yhp => get_Ucomp(8,gid)
        yhm => get_Ucomp(9,gid)
        yh2p => get_Ucomp(10,gid)
        yco => get_Ucomp(11,gid)
        yhn(:,:,:) = yhni
        yh2(:,:,:) = yh2i
        yel(:,:,:) = yeli
        yhp(:,:,:) = yhpi
        yhm(:,:,:) = yhmi
        yh2p(:,:,:) = yh2pi
        yco(:,:,:) = ycoi
        Nrad => get_Ucomp(22,gid)
        Frad_x => get_Ucomp(23,gid)
        Frad_y => get_Ucomp(24,gid)
        Frad_z => get_Ucomp(25,gid)
        Nrad(:,:,:) = 0.d0
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
        Nrad => get_Ucomp(26,gid)
        Frad_x => get_Ucomp(27,gid)
        Frad_y => get_Ucomp(28,gid)
        Frad_z => get_Ucomp(29,gid)
        Nrad(:,:,:) = 0.d0
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
        urad_cmb = cgs_asb*2.73**4.d0/MP_Crd 
        nrad_cmb = urad_cmb / MP_hnu_IR 
        nrad_cmb = nrad_cmb / MP_PHON * Unit_l3 
        Nrad => get_Ucomp(30,gid)
        Frad_x => get_Ucomp(31,gid)
        Frad_y => get_Ucomp(32,gid)
        Frad_z => get_Ucomp(33,gid)
        Nrad(:,:,:) = nrad_cmb
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
        rho => get_Ucomp(0,gid)
        p => get_Ucomp(4,gid)
        rho(:,:,:) = rho0
        p(:,:,:) = rho0 * cs0 * cs0
        if (MP_GasType == 1 .or. MP_GasType == 2) then 
           yh2i2=1.d-8
           yhpi2=1.d-8
           yhmi2=1.d-20
           yh2pi2=1.d-20
           yhni2=1.d0 - (yhpi2+yhmi2) - 2*(yh2i2+yh2pi2)
           yeli2=yhpi2 - yhmi2 + yh2pi2
           rho02 = MP_N0 * cgs_amu * MP_mu / Unit_rho 
           cs02 = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
                /(yhni2+yh2i2+yeli2+yhpi2+yhmi2+yh2pi2+yHe+MP_frac_COsum)))/ Unit_v 
           ncell_HI = 8
           hmax = h/2.**MP_Lmax0 
           do k = Kmingh, Kmaxgh
              do j = Jmingh, Jmaxgh
                 do i = Imingh, Imaxgh
                    radius=sqrt(x(i)**2 + y(j)**2 + z(k)**2)
                    if (radius < ncell_HI * hmax) then
                       rho(i,j,k) = rho02
                       p(i,j,k) = rho02 * cs02 * cs02
                       yhn(i,j,k) = yhni2
                       yh2(i,j,k) = yh2i2
                       yel(i,j,k) = yeli2
                       yhp(i,j,k) = yhpi2
                       yhm(i,j,k) = yhmi2
                       yh2p(i,j,k) = yh2pi2
                    end if
                    if (MP_GasType == 2) then
                       ncell_height = 1
                       if (z(k)**2 < (ncell_height * hmax)**2*0.99) then 
                          rho(i,j,k) = 1d6 * rho(i,j,k)
                          p(i,j,k) = 1d6 * p(i,j,k)
                       end if
                    end if
                  end do
              end do
           end do
        end if 
        do k = Kmingh, Kmaxgh
           do j = Jmingh, Jmaxgh
              do i = Imingh, Imaxgh
                yh2(i,j,k) = 1.3d-6
                yhp(i,j,k) = 1.7d-3
                yhn(i,j,k) = 1.d0-yh2(i,j,k)*2.d0-yhp(i,j,k)
                yel(i,j,k) = yhp(i,j,k)
                yhm(i,j,k) = 1.d-5
                yh2p(i,j,k) = 1.d-5
              end do
           end do
        end do
        khpi => get_Ucomp(12,gid)
        heat_hpi => get_Ucomp(13,gid)
        fx_hpi => get_Ucomp(14,gid)
        fy_hpi => get_Ucomp(15,gid)
        fz_hpi => get_Ucomp(16,gid)
        khpi(:,:,:) = 0.d0
        heat_hpi(:,:,:) = 0.d0
        fx_hpi(:,:,:) = 0.d0
        fy_hpi(:,:,:) = 0.d0
        fz_hpi(:,:,:) = 0.d0
     enddo
  end do
contains
  subroutine writeparam
    use string, only : concat, CHARLEN
    use io_util, only : readenv, wchar
    use modelParameter
    use primordial
    character(len=CHARLEN) :: fn, dir
    if (.not. readenv('DIR', dir) ) stop
    if (get_myrank() == 0) then
       fn = concat(dir, "param.txt")
       call wchar(6,'write parameter file = '//fn)
       open(1, file=fn)
       write(1,'(1P3E15.7)') &
            MP_N0,MP_T0,MP_Boxsize
       write(1,'(2I10,1P1E15.7,I10,1P1E15.7)') &
            MP_Lmax0, MP_JeansConst, MP_spNcr, MP_spRadius_cell, MP_spRadius_lamJ
       write(1,'(6I10)') &
            8, 8, 8, 8, 8, 8
       write(1,'(1P7E15.7)') &
            Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, yHe
       close(1)
    end if
  end subroutine writeparam
end subroutine make_physvar
