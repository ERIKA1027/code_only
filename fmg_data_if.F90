#include "config.h"
! #define CHECK_LINKLIST
  !-------------------------------------------------------------------------
  ! prepare data  in LoL
  !-------------------------------------------------------------------------
  subroutine fmg_prepare_data
    call fmg_make_linklist
    call fmg_make_physvar
  end subroutine fmg_prepare_data
  !-------------------------------------------------------------------------
  ! make physical variables
  !-------------------------------------------------------------------------
  subroutine fmg_make_physvar
    use grid, only : Gidmin, Gidmax, GidListMax, GidList
    use parameter, only : Pi4
    use modelParameter, only : MP_Gconst
    use fmg_data
    use mpilib
    use fmg_ghostcell
    use fmg_boundary_phys
    real(kind=DBL_KIND) :: mass, massg, vol, volg, dv, rhoave, pi4g
    integer :: fmglev, amrlev, n, gid, is, ie, js, je, ks, ke
    integer :: timesliceSource, timesliceUnknwon, mb
    ! ---------------
    ! define varlue
    ! ---------------
#ifdef SINGLE_STEP
    timesliceSource = 1         ! predictor
#else
    timesliceSource = 3         ! interpolation
#endif
    timesliceUnknwon = 0        ! current
    fmglev = FMG_LevelMin

#ifdef FMG_POISSON
    ! -------------------
    ! for self-gravity
    ! -------------------
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       pi4g = MP_Gconst * Pi4      ! 4 pi G
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_alloc_arr(fmglev, IU)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do n = Gidmin, GidListMax(amrlev)
             gid = GidList(n, amrlev) ! grid ID in AMR data
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%U,   (/MPSI/), gid, amrlev, timesliceUnknwon)

#ifdef DM_POTENTIAL
             call fmg_prepare_interp_dm(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/MRHO/), gid, amrlev, timesliceSource)
#else
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/MRHO/), gid, amrlev, timesliceSource)
#endif
             GridLevel(amrlev,fmglev)%Block(n)%Rho = GridLevel(amrlev,fmglev)%Block(n)%Rho * pi4g ! Rho := 4 Pi G rho
          enddo
       enddo
       ! --------------------------------
       ! for periodic boundary condition
       ! --------------------------------
#ifdef ZERO_AVERAGE_DENSITY
#define SZ is:ie,js:je,ks:ke,Mmin:Mmin
       myrank = get_myrank()
       mass = 0.d0
       vol = 0.d0
       fmglev = FMG_LevelMin

       call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          dv = fmg_get_dv(amrlev,fmglev)
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi) cycle
             mass = mass + SUM(GridLevel(amrlev,fmglev)%Block(gid)%Rho(SZ)) * dv
             vol = vol + (ie-is+1)*(je-js+1)*(ke-ks+1)*dv
          enddo
       enddo
       call mpi_allreduce(mass, massg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       call mpi_allreduce(vol,  volg,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       rhoave = massg / volg
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             GridLevel(amrlev,fmglev)%Block(gid)%Rho = GridLevel(amrlev,fmglev)%Block(gid)%Rho - rhoave
          enddo
       enddo
!!$       print *, 'density average ', rhoave
#endif !ZERO_AVERAGE_DENSITY
    endif
#endif !FMG_POISSON
#ifdef FMG_DIFFUSION
    ! ----------------------
    ! for diffusion equation
    ! ----------------------
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       mb = MBX
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IETA)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do n = Gidmin, GidListMax(amrlev)
             gid = GidList(n, amrlev) ! grid ID in AMR data
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%U,   (/mb/),   gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/mb/),   gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Eta, (/MRHO/), gid, amrlev, timesliceUnknwon)
          enddo
       enddo
       call fmg_diff_prepare_eta(IETA)
       call fmg_ghostcell_fix(fmglev, IETA)
       call fmg_boundary_u(fmglev, IETA)
       call fmg_diff_prepare_source(IRHO, IETA)
    endif
#endif !FMG_DIFFUSION

#ifdef FMG_OHMIC_DISSIPATION
    ! -----------------------
    ! for Ohmic dissipation
    ! -----------------------
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IETA)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do n = Gidmin, GidListMax(amrlev)
             gid = GidList(n, amrlev) ! grid ID in AMR data
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%U, (/MBX,MBY,MBZ/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/MBX,MBY,MBZ/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Eta, (/MRHO/), gid, amrlev, timesliceUnknwon)
          enddo
       enddo
       call fmg_od_prepare_eta(IETA)
       call fmg_ghostcell_fix(fmglev, IETA)
       call fmg_boundary_u(fmglev, IETA)
       call fmg_od_prepare_source(IRHO)
    endif
#endif !FMG_OHMIC_DISSIPATION

#ifdef FMG_AMBIPOLAR_DIFFUSION
    ! -----------------------
    ! for ambipolar diffusion
    ! -----------------------
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IDOD)
       call fmg_alloc_arr(fmglev, IDHE)
       call fmg_alloc_arr(fmglev, IDAD)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do n = Gidmin, GidListMax(amrlev)
             gid = GidList(n, amrlev) ! grid ID in AMR data
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%U, (/MBX,MBY,MBZ/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/MBX,MBY,MBZ/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Dod, (/MRHO/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Dhe, (/MRHO/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Dad, (/MRHO/), gid, amrlev, timesliceUnknwon)
!!$             GridLevel(amrlev,fmglev)%Block(n)%Rho = GridLevel(amrlev,fmglev)%Block(n)%U
!!$             GridLevel(amrlev,fmglev)%Block(n)%Dhe = GridLevel(amrlev,fmglev)%Block(n)%Dod
!!$             GridLevel(amrlev,fmglev)%Block(n)%Dad = GridLevel(amrlev,fmglev)%Block(n)%Dod
          enddo
       enddo
       call fmg_ad_prepare_diffcoeff(IDOD, IDHE, IDAD)
       call fmg_ghostcell_fix(fmglev, IDOD, cubic=.TRUE.)
       call fmg_ghostcell_fix(fmglev, IDHE, cubic=.TRUE.)
       call fmg_ghostcell_fix(fmglev, IDAD, cubic=.TRUE.)
       call fmg_boundary_u(fmglev, IDOD)
       call fmg_boundary_u(fmglev, IDHE)
       call fmg_boundary_u(fmglev, IDAD)
       call fmg_ad_prepare_source(IRHO)
    endif
#endif !FMG_AMBIPOLAR_DIFFUSION

  end subroutine fmg_make_physvar
  !-------------------------------------------------------------------------
  ! interpolate U in time (time slice of U)
  ! timeslice = 0 : current 2nd order
  ! timeslice = 1 : predictor step
  ! timeslice = 2 : previous corrector step
  ! timeslice = 3 : interpolation in time
  !-------------------------------------------------------------------------
  subroutine fmg_prepare_interp(ufmg, icode, gid, amrlev, timeslice)
    use grid, only : get_Ucomp, get_U1comp, get_U2comp, LevelMax, Step, Dstep
    integer,intent(IN),dimension(:) :: icode
    integer,intent(IN) :: gid, amrlev
    integer,intent(IN) :: timeslice
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: u0, u1, u2
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ufmg
    real(kind=DBL_KIND) :: stp, dstp, stp1, stp2
    integer :: m
    if (timeslice == 0) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 1) then
       do m = lbound(icode,1), ubound(icode,1)
          u1 => get_U1comp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u1(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 2) then
       do m = lbound(icode,1), ubound(icode,1)
          u2 => get_U2comp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u2(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 3) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid)  ! 2nd order future (coarser) or currnet (finer)
          u1 => get_U1comp(icode(m), gid) ! 1st order future value (coarser) or undef (finer)
          u2 => get_U2comp(icode(m), gid) ! 2nd order past value (coarser) or undef (finer)
          if (Step(amrlev) == Step(LevelMax) ) then ! for finer grid
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(ARRAYSIZE3(ufmg))
          else                                  ! for coarser grid
             stp = Step(LevelMax)
             dstp = Dstep(amrlev)
             stp1 = stp - ( Step(amrlev) - Dstep(amrlev) )
             stp2 = Step(amrlev) - stp
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = (u0(ARRAYSIZE3(ufmg)) - (2*u1(ARRAYSIZE3(ufmg))-u2(ARRAYSIZE3(ufmg))))*(stp1/dstp)**2 &
                  + (u2(ARRAYSIZE3(ufmg))*stp2 + (2*u1(ARRAYSIZE3(ufmg))-u2(ARRAYSIZE3(ufmg)))*stp1)/dstp
          endif
       end do
    endif
  end subroutine fmg_prepare_interp


#ifdef DM_POTENTIAL
  subroutine fmg_prepare_interp_dm(ufmg, icode, gid, amrlev, timeslice)
    use grid, only : get_Ucomp, get_U1comp, get_U2comp, LevelMax, Step, Dstep
    integer,intent(IN),dimension(:) :: icode
    integer,intent(IN) :: gid, amrlev
    integer,intent(IN) :: timeslice
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: u0, u1, u2, dmdens
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ufmg
    real(kind=DBL_KIND) :: stp, dstp, stp1, stp2
    integer :: m
    if (timeslice == 0) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid)
          dmdens => get_Ucomp(MDMRHO, gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(ARRAYSIZE3(ufmg)) + dmdens(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 1) then
       do m = lbound(icode,1), ubound(icode,1)
          u1 => get_U1comp(icode(m), gid)
          dmdens => get_Ucomp(MDMRHO, gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u1(ARRAYSIZE3(ufmg)) + dmdens(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 2) then
       do m = lbound(icode,1), ubound(icode,1)
          u2 => get_U2comp(icode(m), gid)
          dmdens => get_Ucomp(MDMRHO, gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u2(ARRAYSIZE3(ufmg)) + dmdens(ARRAYSIZE3(ufmg))
       end do
    elseif (timeslice == 3) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid)  ! 2nd order future (coarser) or currnet (finer)
          u1 => get_U1comp(icode(m), gid) ! 1st order future value (coarser) or undef (finer)
          u2 => get_U2comp(icode(m), gid) ! 2nd order past value (coarser) or undef (finer)
          dmdens => get_Ucomp(MDMRHO, gid)
          if (Step(amrlev) == Step(LevelMax) ) then ! for finer grid
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(ARRAYSIZE3(ufmg)) + dmdens(ARRAYSIZE3(ufmg))
          else                                  ! for coarser grid
             stp = Step(LevelMax)
             dstp = Dstep(amrlev)
             stp1 = stp - ( Step(amrlev) - Dstep(amrlev) )
             stp2 = Step(amrlev) - stp
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = (u0(ARRAYSIZE3(ufmg)) - (2*u1(ARRAYSIZE3(ufmg))-u2(ARRAYSIZE3(ufmg))))*(stp1/dstp)**2 &
                  + (u2(ARRAYSIZE3(ufmg))*stp2 + (2*u1(ARRAYSIZE3(ufmg))-u2(ARRAYSIZE3(ufmg)))*stp1)/dstp + dmdens(ARRAYSIZE3(ufmg))
          endif
       end do
    endif
  end subroutine fmg_prepare_interp_dm
#endif

  !-------------------------------------------------------------------------
  ! restore U frm FMG data to AMR data
  ! INPUT:
  !   icode = list of AMR components
  !   ju    = FMG components
  !-------------------------------------------------------------------------
  subroutine fmg_restore_u(icode, ju)
    use fmg_data
    use grid, only : get_Ucomp, GidList, GidListMax, Imin,Imax,Jmin,Jmax,Kmin,Kmax
    integer,intent(IN),dimension(:) :: icode
    integer,intent(IN) :: ju
    integer :: fmglev, amrlev, n, m, gid
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: uamr
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: ufmg
    fmglev = FMG_LevelMin       ! coarest grid
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) ! grid ID in AMR data
          do m = lbound(icode,1), ubound(icode,1)
             uamr => get_Ucomp(icode(m), gid)
             ufmg => fmg_get_arrp( amrlev, fmglev, n, ju)
             uamr(ARRAYSIZE3(ufmg)) = ufmg(:,:,:,Mmin-lbound(icode,1)+m)
          end do
       end do
    end do
  end subroutine fmg_restore_u
  !-------------------------------------------------------------------------
  ! make link list for fmg_data
  !-------------------------------------------------------------------------
  subroutine fmg_make_linklist
    use mpilib
    use fmg_data
    use grid, only : Gidmin, Gidmax, GidListMax, GidList, ParentGid, ParentRank, ChildGid, ChildRank, encode_surf, NeighborGid, NeighborRank
    use boundary, only : touch_boundary, NBOUNDARY
    integer :: amrlev, n, gid, lri, lrj, lrk, lr, ndir, pgid, prank, cgid, crank, ngid, nrank, npgid, nprank, rank, glistmax
    integer,dimension(Gidmin:Gidmax,0:NPE-1) :: revlink
    integer,dimension(Gidmin:Gidmax) :: myrevlink
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    integer,dimension(:),allocatable :: glist

    ! node global な逆リンク n = revlink(gid,rank) を作る
    revlink(:,:) = Undefi
    myrevlink(:) = Undefi
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) ! grid ID in AMR data
          myrevlink(gid) = n    ! reverse link list
       end do
    end do
    call mpi_allgather(myrevlink, size(myrevlink), MPI_INTEGER, revlink, size(myrevlink), MPI_INTEGER, MPI_COMM_WORLD, ierr )

    ! link list を受け継ぐ
    ! Ancestry()%Block%..が初期化されていること (fmg_data.F90の宣言文で初期化済み)
    do rank = 0, NPE-1
       do amrlev = AMR_LevelMin, AMR_LevelMax

          ! GidList を各 node に転送する。結果 = glist
          myrank = get_myrank()
          if (rank == myrank) glistmax = GidListMax( amrlev )
          call mpi_bcast( glistmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
          allocate(glist(Gidmin:glistmax+Gidmin))
          if (rank == myrank) glist = GidList(ARRAYSIZE1(glist), amrlev)
          call mpi_bcast( glist, size(glist), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)

          do n = Gidmin, glistmax !   n = grid ID in fmg data
             gid = glist(n)       ! gid = grid ID in AMR data

             ! 親へのリンク
             if ( amrlev > AMR_LevelMin ) then
                pgid = ParentGid(gid,rank)
                prank = ParentRank(gid,rank)
                Ancestry(amrlev,rank)%Block(n)%ParentGid = revlink(pgid,prank)
                Ancestry(amrlev,rank)%Block(n)%ParentRank = prank
             endif
             ! 子へのリンク
             if ( amrlev < AMR_LevelMax ) then
                if ( ChildGid(Left,Left,Left,gid,rank) /= Undefi ) then
                   do lrk = Left ,Right
                      do lrj = Left ,Right
                         do lri = Left ,Right
                            cgid = ChildGid(lri,lrj,lrk,gid,rank)
                            crank = ChildRank(lri,lrj,lrk,gid,rank)
                            Ancestry(amrlev,rank)%Block(n)%ChildGid(lri,lrj,lrk) = revlink(cgid,crank)
                            Ancestry(amrlev,rank)%Block(n)%ChildRank(lri,lrj,lrk) = crank
                         enddo
                      enddo
                   enddo
                endif
             endif
             ! 隣へのリンク
             do ndir = MX, MZ
                do lr = Left, Right
                   if ( NeighborGid(lr,ndir,gid,rank) /= Undefi ) then
                      ngid = NeighborGid(lr,ndir,gid,rank)
                      nrank= NeighborRank(lr,ndir,gid,rank)
                      Ancestry(amrlev,rank)%Block(n)%NeighborGid(lr,ndir) = revlink(ngid,nrank)
                      Ancestry(amrlev,rank)%Block(n)%NeighborRank(lr,ndir) = nrank
                      Ancestry(amrlev,rank)%Block(n)%NeighborSameLevel(lr,ndir) = .TRUE.
                   endif
                enddo
             enddo

          enddo
          deallocate(glist)
       enddo
    enddo

    ! 隣が親の場合
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1)+1, ubound(Ancestry,1)
          do n = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             pgid = Ancestry(amrlev,rank)%Block(n)%ParentGid
             prank = Ancestry(amrlev,rank)%Block(n)%ParentRank
             do ndir = MX, MZ
                do lr = Left, Right
                   if ( Ancestry(amrlev,rank)%Block(n)%NeighborSameLevel(lr,ndir) ) cycle
                   npgid  = Ancestry(amrlev-1,prank)%Block(pgid)%NeighborGid(lr,ndir)
                   nprank = Ancestry(amrlev-1,prank)%Block(pgid)%NeighborRank(lr,ndir)
                   if ( npgid == Undefi ) cycle
                   Ancestry(amrlev,rank)%Block(n)%NeighborGid(lr,ndir) = npgid
                   Ancestry(amrlev,rank)%Block(n)%NeighborRank(lr,ndir) = nprank
                   Ancestry(amrlev,rank)%Block(n)%NeighborParentLevel(lr,ndir) = .TRUE.
                enddo
             enddo
          enddo
       enddo
    enddo

    ! 境界条件の情報を受け継ぐ
    do amrlev = lbound(Geom,1), ubound(Geom,1)
       do n = lbound(Geom(amrlev)%Block,1), ubound(Geom(amrlev)%Block,1)
          call touch_boundary( GidList(n, amrlev), bool_touch )
          do ndir = MX, MZ
             do lr = Left, Right
                Geom(amrlev)%Block(n)%TouchBoundary(lr,ndir) = bool_touch(encode_surf(lr,ndir))
             enddo
          enddo
       enddo
    enddo
#ifdef CHECK_LINKLIST
    call fmg_check_linklist
#endif ! CHECK_LINKLIST
  end subroutine fmg_make_linklist
  !-------------------------------------------------------------------------
  ! check link list for fmg_data
  !-------------------------------------------------------------------------
  subroutine fmg_check_linklist
    use grid, only : encode_surf, swap_LR, decode_LR
    integer :: rank, amrlev, fmglev, gid, pgid, prank, cgid, crank, ngid, nrank, lri, lrj, lrk, lr, lrs, ndir
    logical :: bool, bl
    myrank = get_myrank()
    ! size check
    do fmglev = FMG_LevelMin, FMG_LevelMax
       do amrlev = AMR_LevelMin, AMR_LevelMax
          if (ubound(Ancestry(amrlev,myrank)%Block,1) /= ubound(GridLevel(amrlev,fmglev)%Block,1)) &
               print *, '**** error'
          if (lbound(Ancestry(amrlev,myrank)%Block,1) /= lbound(GridLevel(amrlev,fmglev)%Block,1)) &
               print *, '**** error'
       end do
    end do
    do rank = 0, NPE -1
       do amrlev = AMR_LevelMin, AMR_LevelMax
          if (ubound(Ancestry(amrlev,myrank)%Block,1) /= ubound(Ancestry(amrlev,rank)%Block,1) ) &
               print *, '**** error'
          if (lbound(Ancestry(amrlev,myrank)%Block,1) /= lbound(Ancestry(amrlev,rank)%Block,1) ) &
               print *, '**** error'
       end do
    end do
    ! parent -> child -> self
    bool = .true.
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1)+1, ubound(Ancestry,1)
          do gid = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             pgid  = Ancestry(amrlev,rank)%Block(gid)%ParentGid
             prank = Ancestry(amrlev,rank)%Block(gid)%ParentRank
             bl = .false.
             do lrk = Left ,Right
                do lrj = Left ,Right
                   do lri = Left ,Right
                      if (&
                           gid == Ancestry(amrlev-1,prank)%Block(pgid)%ChildGid(lri,lrj,lrk) .and. &
                           rank == Ancestry(amrlev-1,prank)%Block(pgid)%ChildRank(lri,lrj,lrk))  &
                           bl = .true.
                   end do
                end do
             end do
             bool = bool .and. bl
          enddo
       enddo
    enddo
    if (bool) then
       print *, myrank, 'OK: parent -> child -> self'
    else
       print *, myrank, '** BAD: parent -> child -> self'
    endif
    ! child -> parent -> self
    bool = .true.
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1), ubound(Ancestry,1)-1
          do gid = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             if ( Ancestry(amrlev,rank)%Block(gid)%ChildGid(Left, Left, Left) == Undefi ) cycle
             do lrk = Left ,Right
                do lrj = Left ,Right
                   do lri = Left ,Right
                      cgid = Ancestry(amrlev,rank)%Block(gid)%ChildGid(lri,lrj,lrk)
                      crank = Ancestry(amrlev,rank)%Block(gid)%ChildRank(lri,lrj,lrk)
                      if (.not. &
                           (gid == Ancestry(amrlev+1,crank)%Block(cgid)%ParentGid .and. &
                           rank == Ancestry(amrlev+1,crank)%Block(cgid)%ParentRank )) &
                           bool = .false.
                   enddo
                end do
             end do
          enddo
       enddo
    enddo
    if (bool) then
       print *, myrank, 'OK: child -> parent -> self'
    else
       print *, myrank, '** BAD: child -> parent -> self'
    end if
    ! neighbor
    bool = .TRUE.
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1), ubound(Ancestry,1)-1
          do gid = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             do ndir = MX, MZ
                do lr = Left, Right
                   ngid = Ancestry(amrlev,rank)%Block(gid)%NeighborGid(lr,ndir)
                   nrank = Ancestry(amrlev,rank)%Block(gid)%NeighborRank(lr,ndir)
                   if ( Ancestry(amrlev,rank)%Block(gid)%NeighborSameLevel(lr,ndir) ) then
                      if (ngid == Undefi) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      endif
                      if (nrank == MPI_PROC_NULL) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      endif
                      lrs = decode_LR(swap_LR(encode_surf(lr, ndir)))
                      if (gid /= Ancestry(amrlev,nrank)%Block(ngid)%NeighborGid(lrs,ndir) ) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      end if
                      if (rank /= Ancestry(amrlev,nrank)%Block(ngid)%NeighborRank(lrs,ndir) ) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      end if
                      if (Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir)) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      end if
                      if (rank == myrank .and. Geom(amrlev)%Block(gid)%TouchBoundary(lr,ndir) ) then
                         print *, '*** error in link list neighbor'
                         bool = .false.
                      end if
                   else         !兄弟なし
                      if ( Ancestry(amrlev,rank)%Block(gid)%NeighborParentLevel(lr,ndir) ) then
                         if (ngid == Undefi) then
                            print *, '*** error in link list neighbor'
                            bool = .false.
                         endif
                         if (nrank == MPI_PROC_NULL) then
                            print *, '*** error in link list neighbor'
                            bool = .false.
                         endif
                      end if
                   endif
                enddo
             enddo
          enddo
       end do
    end do
    if (bool) then
       print *, myrank, 'OK: Neighbor link'
    else
       print *, myrank, '** BAD: Neighbor link'
    end if

    bool = .true.
    do amrlev = lbound(Ancestry,1), ubound(Ancestry,1)
       do gid = lbound(Ancestry(amrlev,myrank)%Block,1), ubound(Ancestry(amrlev,myrank)%Block,1)
          do ndir = MX, MZ
             do lr = Left, Right
                if ( &
                     count( (/ Ancestry(amrlev,myrank)%Block(gid)%NeighborParentLevel(lr,ndir), &
                     Ancestry(amrlev,myrank)%Block(gid)%NeighborSameLevel(lr,ndir), &
                     Geom(amrlev)%Block(gid)%TouchBoundary(lr,ndir) /) ) /= 1 ) then
                   bool = .false.
                end if
             end do
          end do
       end do
    end do
    if (bool) then
       print *, myrank, 'OK: Neighbor linkf'
    else
       print *, myrank, '** BAD: Neighbor linkf'
    end if
  end subroutine fmg_check_linklist

