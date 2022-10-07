! Turbulent velocity with uniform density
#include "config.h"
#define FNP_  "param.txt"
!-------------------------------------------------------------------------
! Make initial condition
!-------------------------------------------------------------------------
program main
  use io
  use grid
  use outputdata ! KS ADDED (c.f. init_dynamo.F90)
  implicit none
  call init
  call make_basegrid
  call make_physvar
  call dumpdata
  call output_data ! KS ADDED
!!$  call restoredata
  call finalize
end program main
!-------------------------------------------------------------------------
! Initialize proguram
!-------------------------------------------------------------------------
subroutine init
  use mpilib
  use grid
  use parameter
  use unit
#ifdef METAL
  use kinzoku
#endif
  use modelParameter
  implicit none
  call mpilib_init
  call grid_init
  call parameter_init
  call unit_init
#ifdef METAL
  call init_kinzoku
#endif
  call modelParameter_init
  call eos_init
end subroutine init
!-------------------------------------------------------------------------
! Initialize proguram
!-------------------------------------------------------------------------
subroutine finalize
  use mpilib
  use grid
  call grid_finalize
  call mpi_finalize(ierr)
end subroutine finalize
!--------------------------------------------------------------------
! partition of grid for parallel computation
!--------------------------------------------------------------------
subroutine make_basegrid
  use grid
  use mpilib
  use refine, only : refine_get_norder
  implicit none
  integer :: n, i, j, k, in, jn, kn, gcount, gid, rankid, rank, level, gidn, rankn, sz, szb
  integer :: ngrid_total = (NGI_BASE) * (NGJ_BASE) * (NGK_BASE)
  integer,dimension(0:NPE-1) :: ngrid_node
  integer,dimension(:,:),allocatable :: ig, jg, kg
  ! Hillbert filling curve
  integer :: norder
  integer,parameter :: R = 0, U = 1, L = 2, D = 3, B = 4, F =5
  integer :: ifc, jfc, kfc
  
  LevelMax = 0                  ! the first level of hilbert filling curve
  ! width of filled space = 2 ** norder
  norder = refine_get_norder(0, 0, 0, NGI_BASE-1, NGJ_BASE-1, NGK_BASE-1)
  ! number of grid (block) per node
  ngrid_node(:) = int(ngrid_total/(NPE))
  do n = 1, mod( ngrid_total, NPE )
     ngrid_node(NPE-n) = ngrid_node(NPE-n) + 1
  enddo
  ! print *, ngrid_node, norder
  ! fill space
  level = 0                     ! grid level
  rankid = 0                    ! rank id
  myrank = get_myrank()         ! myrank
  gcount = 0                    ! counter for grid
  gid = -1                      ! counter of grid id
  ifc = 0 ; jfc = 0 ; kfc = 0   ! point of filling curve
  call assign()
  call fillcells(norder, R,U,L,D,B,F)
  call update_gidlist(level) 

  do n = lbound(GidList,1), GidListMax( level ) ! make associated link list
     call alloc_U1order(GidList(n, level))
     call alloc_U2order(GidList(n, level))
  enddo

#define L 0
#define R 1
  ! --------------------------------------------------------------
  ! 隣を探すために Igrid, Jgrid, Kgrid を global にする (ig,jg,kg)
  ! --------------------------------------------------------------
  call mpi_allreduce( GidListMax( level )+1, sz, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr ) 
  allocate( &
       ig(Gidmin:sz-1,0:NPE-1), &
       jg(Gidmin:sz-1,0:NPE-1), &
       kg(Gidmin:sz-1,0:NPE-1) )
  ig = -1
  jg = -1
  kg = -1
!!$  allocate( &
!!$       ig(Gidmin:GidListMax( level ),0:NPE-1), &
!!$       jg(Gidmin:GidListMax( level ),0:NPE-1), &
!!$       kg(Gidmin:GidListMax( level ),0:NPE-1) )
  call mpi_allgather( Igrid, sz, MPI_INTEGER, ig, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  call mpi_allgather( Jgrid, sz, MPI_INTEGER, jg, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  call mpi_allgather( Kgrid, sz, MPI_INTEGER, kg, sz, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  ! -----------------------------------
  ! 隣を探す
  ! NeighborGid, NeighborRank を決める。
  ! -----------------------------------
  do k = 0, (NGK_BASE) -1
     do j = 0, (NGJ_BASE) -1
        do i = 0, (NGI_BASE) -1
           do kn = 0, (NGK_BASE) -1
              do jn = 0, (NGJ_BASE) -1
                 do in = 0, (NGI_BASE) -1
                    gidn = GidBase(in,jn,kn)
                    rankn = RankBase(in,jn,kn)
                    rank = RankBase(i,j,k)
                    gid = GidBase(i,j,k)
                    if ( ig(gid,rank) - 1 == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rankn) ) then
                       NeighborGid(L, MX, gid, rank) = gidn
                       NeighborRank(L, MX, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) + 1 == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rankn) ) then
                       NeighborGid(R, MX, gid, rank) = gidn
                       NeighborRank(R, MX, gid, rank) = rankn
                    endif

                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) - 1 == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rankn) ) then
                       NeighborGid(L, MY, gid, rank) = gidn
                       NeighborRank(L, MY, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) + 1 == jg(gidn,rankn) .and. kg(gid,rank) == kg(gidn,rankn) ) then
                       NeighborGid(R, MY, gid, rank) = gidn
                       NeighborRank(R, MY, gid, rank) = rankn
                    endif

                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) - 1 == kg(gidn,rankn) ) then
                       NeighborGid(L, MZ, gid, rank) = gidn
                       NeighborRank(L, MZ, gid, rank) = rankn
                    endif
                    if ( ig(gid,rank) == ig(gidn,rankn) .and. jg(gid,rank) == jg(gidn,rankn) .and. kg(gid,rank) + 1 == kg(gidn,rankn) ) then
                       NeighborGid(R, MZ, gid, rank) = gidn
                       NeighborRank(R, MZ, gid, rank) = rankn
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

!----------EO_removed----------!

  ! -------------------------------
  ! for periodic boundary condition
  ! -------------------------------
!  i = 0                         ! left side
!  in = (NGI_BASE) -1            ! right side (neighbor of i)
!  do k = 0, (NGK_BASE) -1
!     do j = 0, (NGJ_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(L, MX, gid, rank) = GidBase(in,j,k)
!        NeighborRank(L, MX, gid, rank) = RankBase(in,j,k)
!     enddo
!  enddo
!  i = (NGI_BASE) -1             ! right side
!  in = 0                        ! left side (neighbor of i)
!  do k = 0, (NGK_BASE) -1
!     do j = 0, (NGJ_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(R, MX, gid, rank) = GidBase(in,j,k)
!        NeighborRank(R, MX, gid, rank) = RankBase(in,j,k)
!     enddo
!  enddo
!  j = 0                         ! left side
!  jn = (NGJ_BASE) -1            ! right side (neighbor of j)
!  do k = 0, (NGK_BASE) -1
!     do i = 0, (NGI_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(L, MY, gid, rank) = GidBase(i,jn,k)
!        NeighborRank(L, MY, gid, rank) = RankBase(i,jn,k)
!     enddo
!  enddo
!  j = (NGJ_BASE) -1             ! right side
!  jn = 0                        ! left side (neighbor of j)
!  do k = 0, (NGK_BASE) -1
!     do i = 0, (NGI_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(R, MY, gid, rank) = GidBase(i,jn,k)
!        NeighborRank(R, MY, gid, rank) = RankBase(i,jn,k)
!     enddo
!  enddo
!  k = 0                         ! left side
!  kn = (NGK_BASE) -1            ! right side (neighbor of k)
!  do j = 0, (NGJ_BASE) -1
!     do i = 0, (NGI_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(L, MZ, gid, rank) = GidBase(i,j,kn)
!        NeighborRank(L, MZ, gid, rank) = RankBase(i,j,kn)
!     enddo
!  enddo
!  k = (NGK_BASE) -1             ! right side
!  kn = 0                        ! left side (neighbor of k)
!  do j = 0, (NGJ_BASE) -1
!     do i = 0, (NGI_BASE) -1
!        rank = RankBase(i,j,k)
!        gid = GidBase(i,j,k)
!        NeighborGid(R, MZ, gid, rank) = GidBase(i,j,kn)
!        NeighborRank(R, MZ, gid, rank) = RankBase(i,j,kn)
!     enddo
!  enddo

!----------EO_removed----------!

  deallocate( ig, jg, kg )
#undef L
#undef R


contains
  !----------------------------------------------------------------------
  ! Peano-Hilbert filling curve
  !----------------------------------------------------------------------
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
    if ( ifc < 0 .or. ifc > NGI_BASE-1 .or. &
         jfc < 0 .or. jfc > NGJ_BASE-1 .or. &
         kfc < 0 .or. kfc > NGK_BASE-1 ) return
    gcount = gcount + 1
    gid = gid + 1
    if ( gcount > ngrid_node(rankid) ) then ! 次のランク
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
!-------------------------------------------------------------------------
! Define coordinates and physical values
!-------------------------------------------------------------------------
subroutine make_physvar
  use mpilib
  use grid
  use refine
  use grid_boundary
  use fg2cg
!!$  use eos, only : Cs, Gamma, Rhocr, Kappa
  use parameter
  use io_util, only : readenv
  use unit
  use modelParameter
  use primordial ! KS ADDED  
#ifdef METAL
  use kinzoku, only : fdust_solar
#endif ! METAL
  implicit none
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho, vx, vy, vz, p, psi, bx, by, bz, db, gx, gy, gz

#ifdef DUST_NOTCONSTANT
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: rhod
#endif !DUST_NOTCONSTANT

  real(kind=DBL_KIND),dimension(:,:,:),pointer :: &
       yhn, yh2, yel, yhp, yhm, yh2p,khpi,heat_hpi,fx_hpi,fy_hpi,fz_hpi,yco,rdpco, ychem !EUV RT + thin FUV
#if MODEL_ART == 2 
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: kh2pd !EUV RT + FUV RT
#endif !MODEL_ART

#ifdef METAL
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: kgfuv, kdph, krOII 
#endif !METAL       

#ifdef RADTR_M1closer
  real(kind=DBL_KIND),dimension(:,:,:),pointer :: Nrad, Frad_x, Frad_y, Frad_z
#endif !RADTR_M1closer

  
  real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
  real(kind=DBL_KIND) :: MachNumber
  integer :: level,gid
  integer :: i,j,k, n
  real(kind=DBL_KIND) :: h,ic,jc,kc,ic0,jc0,kc0,ycoi
  real(kind=DBL_KIND) :: r2, x0, y0, z0, a, rho0, p0, vx0, vy0, vz0, psi0, h1d, db0, b0, rboundary, cs0
  logical :: bool
  real(kind=DBL_KIND) :: yhni, yh2i, yeli, yhpi, yhmi, yh2pi
  real(kind=DBL_KIND) :: yhni2, yh2i2, yeli2, yhpi2, yhmi2, yh2pi2, rho02, cs02, radius, hmax
  real(kind=DBL_KIND) :: urad_cmb, nrad_cmb, nrad_fuv_bg, urad_fuv_bg
  integer :: ncell_HI, ncell_height
  integer :: ichem

  !not using art module
#if MODEL_ART == 0
  if (get_myrank() == PRIMARY_RANK) &
       print '(A,/,A)', "init_art_radtest is only compatible with MODEL_ART>=1","stopping..."
  stop
#endif !MODEL_ART
  
  
  ! gas type毎の初期化学組成
  if (MP_GasType == 0) then ! HI gas の場合
     yh2i=1.d-8
     yhpi=1.d-8
     yhmi=1.d-20
     yh2pi=1.d-20
     yhni=1.d0 - (yhpi+yhmi) - 2*(yh2i+yh2pi)
     yeli=yhpi - yhmi + yh2pi
     ycoi = 0.927d-4
  else if (MP_GasType == 1 .or. MP_GasType == 2) then ! H2 gas の場合
     yhni=1.d-8
     yhpi=1.d-8
     yhmi=1.d-20
     yh2pi=1.d-20
     yh2i= (1.d0 - (yhni+yhpi+yhmi) - 2*yh2pi)/2.d0
     yeli=yhpi - yhmi + yh2pi
     ycoi = 0.927d-4
  else
       if (get_myrank() == PRIMARY_RANK) &
       print '(A,/,A)', "MP_GasType should be either 0, 1 or 2","stopping..."
       stop
  end if

     !----------EO_removed----------!
    ! yhni=1.d-8
    ! yhpi=1.d-8
    ! yhmi=1.d-20
    ! yh2pi=1.d-20
    ! yh2i= (1.d0 - (yhni+yhpi+yhmi) - 2*yh2pi)/2.d0
    ! yeli=yhpi - yhmi + yh2pi
    ! ycoi = 0.927d-4
     !---------EO_removed-----------!

  ! density, sound speed
  rho0 = MP_N0 * cgs_amu * MP_mu / Unit_rho                     ! central density

  !----------EO_added----------!
 ! rho0  = MP_R0 *MP_Tau0 / Unit_rho
 ! MP_N0 = rho0 * Unit_rho / (cgs_amu * MP_mu)
  !----------EO_added----------!


  cs0 = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
       /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))/ Unit_v    ! isothermal sound speed

  ! box size
  rboundary = MP_Boxsize

#ifdef MHD 
  b0  = 1.d0 / Unit_ugauss ! [micro gausss]
  db0 = 0.d0
#endif !MHD


  call writeparam   

  level = 0                             ! base level
  Time(:) = 0.d0
  Step(:) = 0.d0
  Dtime(:) = 0.d0
  Dstep(:) = 0.d0
  !--------------------
  ! define coordintates
  !--------------------
  ic0 = ((NGI_BASE) * (NI)-1)/2.d0  ! origin in cell number
  jc0 = ((NGJ_BASE) * (NJ)-1)/2.d0
  kc0 = ((NGK_BASE) * (NK)-1)/2.d0
  h = 2*rboundary/min((NGI_BASE) * (NI), (NGJ_BASE) * (NJ), (NGK_BASE) * (NK) )
  CellWidth(:, level) = (/h, h, h/)
  !---------------------------------
  ! define coordinates in base grid
  !---------------------------------
  do n = Gidmin, GidListMax( level )
     gid = GidList( n, level )
     ic = ic0 - Igrid(gid)* NI 
     jc = jc0 - Jgrid(gid)* NJ
     kc = kc0 - kgrid(gid)* NK
     ! coordinates
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
  ! ----------------------
  ! define physical values
  ! ----------------------
  do level = Lmin, Lmax
     if ( level > Lmin ) then
        call refineLevel(level, bool)
        LevelMax = max(LevelMax, 0)
        if ( .not. bool ) exit
     endif
     do n = Gidmin, GidListMax( level )
        gid = GidList( n, level )

#ifdef MPSI
        psi => get_Ucomp(MPSI,gid)
        gx  => get_Ucomp(MGX,gid)
        gy  => get_Ucomp(MGY,gid)
        gz  => get_Ucomp(MGZ,gid)
#endif !MPSI
        x => get_Xp(gid)
        y => get_Yp(gid)
        z => get_Zp(gid)
#ifdef MPSI
        psi(:,:,:) = 0.d0
        gx(:,:,:) = 0.d0
        gy(:,:,:) = 0.d0
        gz(:,:,:) = 0.d0
#endif !MPSI

        !set velocity
        vx  => get_Ucomp(MVX,gid)
        vy  => get_Ucomp(MVY,gid)
        vz  => get_Ucomp(MVZ,gid)

        !----------EO_added----------!
        vx(:,:,:) = 0.d0        
       ! vx(Imingh,:,:) = 1.d6 / Unit_v
!        vx(:,:,:) = 2.d6 / Unit_v
        vy(:,:,:) = 0.d0
        vz(:,:,:) = 0.d0
        !----------EO_added----------!

#ifdef MHD
        bx  => get_Ucomp(MBX,gid)
        by  => get_Ucomp(MBY,gid)
        bz  => get_Ucomp(MBZ,gid)
        db => get_Ucomp(MDB,gid)
        bx(:,:,:) = 0.d0
        by(:,:,:) = 0.d0
        bz(:,:,:) = b0
        db(:,:,:) = db0
#endif !MHD


        ! set initial chemical abundance
        do ichem = NCEHM_MIN, NCEHM_MAX
          ychem => get_Ucomp(ichem,gid)
          ychem(:,:,:) = 1.d-15
        enddo

        !set initial abundance
        yhn => get_Ucomp(MHN,gid)
        yh2 => get_Ucomp(MH2,gid)
        yel => get_Ucomp(MEL,gid)
        yhp  => get_Ucomp(MHP,gid)
        yhm => get_Ucomp(MHM,gid)
        yh2p => get_Ucomp(MH2P,gid)
        yco  => get_Ucomp(MCO,gid)
        yhn(:,:,:) = yhni
        yh2(:,:,:) = yh2i
        yel(:,:,:) = yeli
        yhp(:,:,:) = yhpi
        yhm(:,:,:) = yhmi
        yh2p(:,:,:) = yh2pi
        yco(:,:,:) = ycoi
        
#ifdef RADTR_M1closer
  #ifdef M1CLOSER_EUV_TRANSFER
        Nrad   => get_Ucomp(MER,gid)
        Frad_x => get_Ucomp(MFRX,gid) 
        Frad_y => get_Ucomp(MFRY,gid) 
        Frad_z => get_Ucomp(MFRZ,gid) 
        Nrad(:,:,:) = 0.d0
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
  #endif !M1CLOSER_EUV_TRANSFER

  #ifdef M1CLOSER_FUV_TRANSFER
        Nrad   => get_Ucomp(MEF,gid)
        Frad_x => get_Ucomp(MFRFX,gid) 
        Frad_y => get_Ucomp(MFRFY,gid) 
        Frad_z => get_Ucomp(MFRFZ,gid) 
    #ifdef EXTERNAL_FUV_RAD
        urad_fuv_bg = EXTERNAL_FUV_G0*1.6d-3/MP_Ctil   ! [ erg cm^-3 ]
        nrad_fuv_bg = urad_fuv_bg / (11.174d0*cgs_ev) ! number density [cm^-3]
        nrad_fuv_bg = nrad_fuv_bg / MP_PHON * Unit_l3 ! [noD]
        Nrad(:,:,:) = nrad_fuv_bg
    #else
        Nrad(:,:,:) = 0.d0
    #endif !EXTERNAL_FUV_RAD
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0

    #ifdef M1CLOSER_SEPARATE_FUV_TRANS
        Nrad   => get_Ucomp(MECOF,gid)
        Frad_x => get_Ucomp(MFCORFX,gid) 
        Frad_y => get_Ucomp(MFCORFY,gid) 
        Frad_z => get_Ucomp(MFCORFZ,gid) 
        Nrad(:,:,:) = 0.d0
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0

        Nrad   => get_Ucomp(MEDUSTF,gid)
        Frad_x => get_Ucomp(MFDUSTRFX,gid) 
        Frad_y => get_Ucomp(MFDUSTRFY,gid) 
        Frad_z => get_Ucomp(MFDUSTRFZ,gid) 
        Nrad(:,:,:) = 0.d0
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
    #endif !M1CLOSER_SEPARATE_FUV_TRANS

  #endif !M1CLOSER_FUV_TRANSFER

  #ifdef M1CLOSER_IR_TRANSFER
        urad_cmb = cgs_asb*DEF_TCMB**4.d0/MP_Crd ! [ erg cm^-3]
        nrad_cmb = urad_cmb / MP_hnu_IR  ! number density [cm^-3]
        nrad_cmb = nrad_cmb / MP_PHON * Unit_l3 ! [noD]

        Nrad   => get_Ucomp(MEIR,gid)
        Frad_x => get_Ucomp(MFIRX,gid) 
        Frad_y => get_Ucomp(MFIRY,gid) 
        Frad_z => get_Ucomp(MFIRZ,gid) 
        Nrad(:,:,:)   = nrad_cmb
        Frad_x(:,:,:) = 0.d0
        Frad_y(:,:,:) = 0.d0
        Frad_z(:,:,:) = 0.d0
  #endif !M1CLOSER_IR_TRANSFER
#endif !RADTR_M1closer


        ! set density&pressue
        rho => get_Ucomp(MRHO,gid)
        p => get_Ucomp(MP,gid)
        rho(:,:,:) = rho0
        p(:,:,:) =  rho0 * cs0 * cs0

#ifdef DUST_NOTCONSTANT
        ! set dust density
        rhod => get_Ucomp(MDRHO,gid)
        rhod(:,:,:) = fdust_solar*MP_Metallicity
#endif !DUST_NOTCONSTANT

        !H2 gasの問題のとき、中性水素が全く存在しないことによって、
        !計算開始直後に電離光子が全領域に届いてしまうことがあり、
        !その場合にartificialな振る舞いがあるようだったので対応

        if (MP_GasType == 1 .or. MP_GasType == 2) then ! H2 gas の場合
           yh2i2=1.d-8
          ! yhni2=1.d-8
           yhpi2=1.d-8
           yhmi2=1.d-20
           yh2pi2=1.d-20
           yhni2=1.d0 - (yhpi2+yhmi2) - 2*(yh2i2+yh2pi2)
         !  yh2i2=1.d0 - (yhpi2+yhmi2) - 2*(yhni2+yh2pi2)
           yeli2=yhpi2 - yhmi2 + yh2pi2
          
           rho02 = MP_N0 * cgs_amu * MP_mu / Unit_rho                     ! central density

           !----------EO_added----------!
          ! rho02 = MP_R0 * MP_Tau0/ Unit_rho 
           !----------EO_added----------!

           cs02 = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
                /(yhni2+yh2i2+yeli2+yhpi2+yhmi2+yh2pi2+yHe+MP_frac_COsum)))/ Unit_v    ! isothermal sound speed
           ! assume HI gas for inner most 8 cells
           ncell_HI = 8
           hmax = h/2.**MP_Lmax0 !最大レベルでのセルサイズ
           do k = Kmingh, Kmaxgh
              do j = Jmingh, Jmaxgh
                 do i = Imingh, Imaxgh
                    radius=sqrt(x(i)**2 + y(j)**2 + z(k)**2)
                    if (radius < ncell_HI * hmax) then
                       rho(i,j,k) = rho02
                       p(i,j,k) =  rho02 * cs02 * cs02
                       yhn(i,j,k) = yhni2
                       yh2(i,j,k) = yh2i2
                       yel(i,j,k) = yeli2
                       yhp(i,j,k) = yhpi2
                       yhm(i,j,k) = yhmi2
                       yh2p(i,j,k) = yh2pi2
#ifdef DUST_NOTCONSTANT
                       rhod(i,j,k) = fdust_solar*MP_Metallicity
#endif !DUST_NOTCONSTANT
                    end if

                    !MP_GasType == 2のとき、赤道面を挟んで上下2セルの密度を1d6倍下駄を履かせてみる
                    if (MP_GasType == 2) then
                       ncell_height = 1
                       if (z(k)**2 < (ncell_height * hmax)**2*0.99) then !0.99はif文の中身が等号だったときのための対応
                          rho(i,j,k) = 1d6 * rho(i,j,k)
                          p(i,j,k) = 1d6 * p(i,j,k)
                       end if
                    end if


                  end do
              end do
           end do           
        end if !MP_GasType==1 or 2

        do k = Kmingh, Kmaxgh
           do j = Jmingh, Jmaxgh
              do i = Imingh, Imaxgh
                yh2(i,j,k) = 1.3d-6
                yhp(i,j,k) = 1.7d-3
                yhn(i,j,k) = 1.d0-yh2(i,j,k)*2.d0-yhp(i,j,k)
                yel(i,j,k) = yhp(i,j,k)
                yhm(i,j,k) = 1.d-5
                yh2p(i,j,k) = 1.d-5

!                yhn(i,j,k) = yhni
!                yh2(i,j,k) = yh2i
!                yel(i,j,k) = yeli
!                yhp(i,j,k) = yhpi
!                yhm(i,j,k) = yhmi
!                yh2p(i,j,k) = yh2pi
!                yco(i,j,k) = ycoi
 
              end do
           end do
        end do           

        ! set initial radiaiton field
        khpi => get_Ucomp(MKPI,gid)
        heat_hpi => get_Ucomp(MHPI,gid)
        fx_hpi => get_Ucomp(MXPI,gid)
        fy_hpi => get_Ucomp(MYPI,gid)
        fz_hpi => get_Ucomp(MZPI,gid)
        khpi(:,:,:) = 0.d0
        heat_hpi(:,:,:) = 0.d0
        fx_hpi(:,:,:) = 0.d0
        fy_hpi(:,:,:) = 0.d0
        fz_hpi(:,:,:) = 0.d0

#ifdef RADTR_DIRECT
#if MODEL_ART == 2
        kh2pd => get_Ucomp(MKPD,gid)
        kh2pd(:,:,:) = 0.d0
#ifdef METAL
        kgfuv => get_Ucomp(MGFUV,gid)
        kgfuv(:,:,:) = 0.d0
        kdph  => get_Ucomp(MDPH,gid)
        kdph(:,:,:)  = 0.d0
        rdpco => get_Ucomp(MDPCO,gid)
        rdpco(:,:,:)  = 0.d0
        krOII => get_Ucomp(MKPOII, gid)
        krOII(:,:,:)  = 0.d0
#endif !METAL
#endif !MODEL_ART
#endif ! RADTR_DIRECT




        
     enddo
  end do

contains
  ! -----------------------------------------------------------------
  ! write parameter to a file
  ! -----------------------------------------------------------------
  subroutine writeparam
    use string, only : concat, CHARLEN
    use io_util, only : readenv, wchar
    use modelParameter
    use primordial
    
    character(len=CHARLEN) :: fn, dir
    if (.not. readenv('DIR', dir) ) stop
    if (get_myrank() == PRIMARY_RANK) then
       fn = concat(dir, FNP_)
       call wchar(6,'write parameter file = '//fn)
       open(1, file=fn)

       !model parameters
       write(1,'(1P3E15.7)') &
            MP_N0,MP_T0,MP_Boxsize
       write(1,'(2I10,1P1E15.7,I10,1P1E15.7)') &
            MP_Lmax0, MP_JeansConst,  MP_spNcr, MP_spRadius_cell, MP_spRadius_lamJ
       write(1,'(6I10)') &
            NI, NJ, NK, NGI_BASE, NGJ_BASE, NGK_BASE
       write(1,'(1P7E15.7)') &
            Unit_rho, Unit_v, Unit_t, Unit_l, Unit_m, Unit_e, yHe
       
       close(1)

    end if
  end subroutine writeparam


end subroutine make_physvar
