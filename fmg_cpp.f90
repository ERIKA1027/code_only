module fmg
  use fmg_data
  implicit none
  private
  real(kind=8),save :: Resh2maxg
  real(kind=8),save,allocatable :: Resmaxg(:), Resmaxl(:), Trerr(:) 
  integer,save :: N_fmg_cycle
  public :: fmg_poisson
contains
  subroutine fmg_psi2g(ju)
    use grid, only : Gidmin, GidListMax, GidList, Imin, Imax, Jmin, Jmax, Kmin, Kmax, get_Ucomp
    use fmg_data, only : FMG_LevelMin, AMR_LevelMin, AMR_LevelMax, fmg_get_h, fmg_arrp, fmg_fp, fmg_skip_grid, Mmin
    integer,intent(IN) :: ju
    integer :: amrlev, fmglev, i, j, k, n, gid
    real(kind=8) :: hi
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    real(kind=8),pointer,dimension(:,:,:) :: gx, gy, gz, u
    fmglev = FMG_LevelMin
    call fmg_poisson_flux(fmglev, ju)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) 
          call fmg_arrp(amrlev, fmglev, n, ju, u)
          call fmg_fp(amrlev, fmglev, n, f)
          gx => get_Ucomp(18, gid)
          gy => get_Ucomp(19, gid)
          gz => get_Ucomp(20, gid)
          if ( fmg_skip_grid(n, amrlev, fmglev) ) then
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i = Imin, Imax
                      gx(i,j,k) = -(u(i+1,j,k) - u(i-1,j,k))*hi*0.5d0
                      gy(i,j,k) = -(u(i,j+1,k) - u(i,j-1,k))*hi*0.5d0
                      gz(i,j,k) = -(u(i,j,k+1) - u(i,j,k-1))*hi*0.5d0
                   enddo
                enddo
             enddo
          else
             do k = Kmin, Kmax
                do j = Jmin, Jmax
                   do i = Imin, Imax
                      gx(i,j,k) = -(f(i,j,k,0)+f(i-1,j,k,0))*0.5d0
                      gy(i,j,k) = -(f(i,j,k,1)+f(i,j-1,k,1))*0.5d0
                      gz(i,j,k) = -(f(i,j,k,2)+f(i,j,k-1,2))*0.5d0
                   enddo
                enddo
             enddo
          endif
       enddo
    enddo
  end subroutine fmg_psi2g
  subroutine fmg_poisson_flux(fmglev, ju, boundary_fill0)
    use fmg_data
    use fmg_boundary_phys, only : fmg_boundary_u
    use fmg_converge
    use fmg_reflux
    use fmg_ghostcell
    use fmg_boundary
    integer,intent(IN) :: fmglev, ju
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=8) :: hi
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    real(kind=8),pointer,dimension(:,:,:) :: u
    integer :: amrlev, gid, ndir, lr, i, j, k, imin,jmin,kmin,imax,jmax,kmax
    if (FmgLevel_fill0 /= fmglev) then
       call fmg_converge_c2p(fmglev,ju)
       call fmg_boundary_u(fmglev, ju)
       call fmg_ghostcell_fix(fmglev,ju)
       if (present(boundary_fill0)) then
          bool_boundary_fill0 = boundary_fill0
       else
          bool_boundary_fill0 = .false.
       endif
       if (bool_boundary_fill0) call fmg_boundary_fill0(fmglev, ju)
       call fmg_boundary_u(fmglev, ju)
    endif
    FmgLevel_fill0 = Undefi
    myrank = get_myrank()
    call fmg_get_gridsize(fmglev, imin,jmin,kmin,imax,jmax,kmax)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       hi = 1.d0/fmg_get_h( amrlev, fmglev )
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          call fmg_fp(amrlev, fmglev, gid, f)
          call fmg_arrp(amrlev, fmglev, gid, ju, u)
          do k = kmin-1, kmax
             do j = jmin-1, jmax
                do i = imin-1, imax
                   f(i,j,k,0)= (u(i+1,j,k)-u(i,j,k))*hi
                   f(i,j,k,1)= (u(i,j+1,k)-u(i,j,k))*hi
                   f(i,j,k,2)= (u(i,j,k+1)-u(i,j,k))*hi
                end do
             end do
          end do
       end do
    end do
    call fmg_fluxcorrection(fmglev)
  end subroutine fmg_poisson_flux
  subroutine fmg_poisson_relax(fmglev, ju, jrhs, resh2maxg, boundary_fill0)
    use mpilib
    use fmg_data
    integer,intent(IN) :: fmglev, ju, jrhs
    real(kind=8),intent(OUT) :: resh2maxg
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=8),parameter :: sixth = 1.d0/6.0d0
    real(kind=8) :: h, dv, h2idv, resh2, resh2max, resh2max_g, ressum
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    real(kind=8),pointer,dimension(:,:,:) :: u, rhs
    integer :: amrlev, gid, ndir, i, j, k, ipass, is,js,ks,ie,je,ke
    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif
    myrank = get_myrank()
    resh2max = 0.d0
    ressum = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do ipass=1,2 
       call fmg_poisson_flux(fmglev, ju, bool_boundary_fill0)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          h = fmg_get_h( amrlev, fmglev )
          dv = fmg_get_dv( amrlev, fmglev )
          h2idv = dv/h**2
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
             call fmg_fp(amrlev, fmglev, gid, f)
             call fmg_arrp(amrlev, fmglev, gid, ju, u)
             call fmg_arrp(amrlev, fmglev, gid, jrhs, rhs)
             do k = ks, ke
                do j = js, je
                   do i= is + mod(j+k+ipass,2), ie, 2
                      resh2 = &
                           (f(i,j,k,0)-f(i-1,j,k,0) &
                           +f(i,j,k,1)-f(i,j-1,k,1) &
                           +f(i,j,k,2)-f(i,j,k-1,2) &
                           -rhs(i,j,k) * h &
                           )* h
                      u(i,j,k) = u(i,j,k) + sixth * resh2
                      resh2max = max(resh2max,abs(resh2))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    call mpi_allreduce(resh2max, Resh2maxg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  end subroutine fmg_poisson_relax
  subroutine fmg_poisson_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    use fmg_data
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    real(kind=8),pointer,dimension(:,:,:) :: res, u, rhs
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    integer :: amrlev, gid, i, j, k, is, ie, js, je, ks, ke
    real(kind=8) :: hi
    if (present(boundary_fill0)) then
       bool_boundary_fill0 = boundary_fill0
    else
       bool_boundary_fill0 = .false.
    endif
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    call fmg_poisson_flux(fmglev, ju, bool_boundary_fill0)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          call fmg_arrp(amrlev,fmglev,gid,jres, res)
          call fmg_arrp(amrlev,fmglev,gid,ju, u)
          call fmg_arrp(amrlev,fmglev,gid,jrhs, rhs)
          call fmg_fp(amrlev, fmglev, gid, f)
          hi = 1.d0/fmg_get_h( amrlev, fmglev )
          do k = ks, ke
             do j = js, je
                do i = is, ie
                   res(i,j,k) = &
                        rhs(i,j,k) &
                       -(f(i,j,k,0)-f(i-1,j,k,0) &
                        +f(i,j,k,1)-f(i,j-1,k,1) &
                        +f(i,j,k,2)-f(i,j,k-1,2))*hi
                end do
             enddo
          enddo
       enddo
    enddo
  end subroutine fmg_poisson_resid
  subroutine fmg_prepare_data
    call fmg_make_linklist
    call fmg_make_physvar
  end subroutine fmg_prepare_data
  subroutine fmg_make_physvar
    use grid, only : Gidmin, Gidmax, GidListMax, GidList
    use parameter, only : Pi4
    use modelParameter, only : MP_Gconst
    use fmg_data
    use mpilib
    use fmg_ghostcell
    use fmg_boundary_phys
    real(kind=8) :: mass, massg, vol, volg, dv, rhoave, pi4g
    integer :: fmglev, amrlev, n, gid, is, ie, js, je, ks, ke
    integer :: timesliceSource, timesliceUnknwon, mb
    timesliceSource = 1 
    timesliceUnknwon = 0 
    fmglev = FMG_LevelMin
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       pi4g = MP_Gconst * Pi4 
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_alloc_arr(fmglev, IU)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do n = Gidmin, GidListMax(amrlev)
             gid = GidList(n, amrlev) 
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%U, (/17/), gid, amrlev, timesliceUnknwon)
             call fmg_prepare_interp(GridLevel(amrlev,fmglev)%Block(n)%Rho, (/0/), gid, amrlev, timesliceSource)
             GridLevel(amrlev,fmglev)%Block(n)%Rho = GridLevel(amrlev,fmglev)%Block(n)%Rho * pi4g 
          enddo
       enddo
       myrank = get_myrank()
       mass = 0.d0
       vol = 0.d0
       fmglev = FMG_LevelMin
       call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
       do amrlev = AMR_LevelMin, AMR_LevelMax
          dv = fmg_get_dv(amrlev,fmglev)
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) /= Undefi) cycle
             mass = mass + SUM(GridLevel(amrlev,fmglev)%Block(gid)%Rho(is:ie,js:je,ks:ke,Mmin:Mmin)) * dv
             vol = vol + (ie-is+1)*(je-js+1)*(ke-ks+1)*dv
          enddo
       enddo
       call mpi_allreduce(mass, massg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       call mpi_allreduce(vol, volg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
       rhoave = massg / volg
       do amrlev = AMR_LevelMin, AMR_LevelMax
          do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
             GridLevel(amrlev,fmglev)%Block(gid)%Rho = GridLevel(amrlev,fmglev)%Block(gid)%Rho - rhoave
          enddo
       enddo
    endif
  end subroutine fmg_make_physvar
  subroutine fmg_prepare_interp(ufmg, icode, gid, amrlev, timeslice)
    use grid, only : get_Ucomp, get_U1comp, get_U2comp, LevelMax, Step, Dstep
    integer,intent(IN),dimension(:) :: icode
    integer,intent(IN) :: gid, amrlev
    integer,intent(IN) :: timeslice
    real(kind=8),dimension(:,:,:),pointer :: u0, u1, u2
    real(kind=8),dimension(:,:,:,:),pointer :: ufmg
    real(kind=8) :: stp, dstp, stp1, stp2
    integer :: m
    if (timeslice == 0) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound&
&(ufmg,3))
       end do
    elseif (timeslice == 1) then
       do m = lbound(icode,1), ubound(icode,1)
          u1 => get_U1comp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u1(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound&
&(ufmg,3))
       end do
    elseif (timeslice == 2) then
       do m = lbound(icode,1), ubound(icode,1)
          u2 => get_U2comp(icode(m), gid)
          ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u2(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound&
&(ufmg,3))
       end do
    elseif (timeslice == 3) then
       do m = lbound(icode,1), ubound(icode,1)
          u0 => get_Ucomp(icode(m), gid) 
          u1 => get_U1comp(icode(m), gid) 
          u2 => get_U2comp(icode(m), gid) 
          if (Step(amrlev) == Step(LevelMax) ) then 
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = u0(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubo&
&und(ufmg,3))
          else 
             stp = Step(LevelMax)
             dstp = Dstep(amrlev)
             stp1 = stp - ( Step(amrlev) - Dstep(amrlev) )
             stp2 = Step(amrlev) - stp
             ufmg(:,:,:,Mmin-lbound(icode,1)+m) = (u0(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ub&
&ound(ufmg,3)) - (2*u1(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3))-u2(lbound(ufmg,1)&
&:ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3))))*(stp1/dstp)**2 &
                  + (u2(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3))*stp2 + (2*u1(lbo&
&und(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3))-u2(lbound(ufmg,1):ubound(ufmg,1),lbound(uf&
&mg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3)))*stp1)/dstp
          endif
       end do
    endif
  end subroutine fmg_prepare_interp
  subroutine fmg_restore_u(icode, ju)
    use fmg_data
    use grid, only : get_Ucomp, GidList, GidListMax, Imin,Imax,Jmin,Jmax,Kmin,Kmax
    integer,intent(IN),dimension(:) :: icode
    integer,intent(IN) :: ju
    integer :: fmglev, amrlev, n, m, gid
    real(kind=8),dimension(:,:,:),pointer :: uamr
    real(kind=8),dimension(:,:,:,:),pointer :: ufmg
    fmglev = FMG_LevelMin 
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) 
          do m = lbound(icode,1), ubound(icode,1)
             uamr => get_Ucomp(icode(m), gid)
             ufmg => fmg_get_arrp( amrlev, fmglev, n, ju)
             uamr(lbound(ufmg,1):ubound(ufmg,1),lbound(ufmg,2):ubound(ufmg,2),lbound(ufmg,3):ubound(ufmg,3)) = ufmg(:,:,:,Mmin-lbou&
&nd(icode,1)+m)
          end do
       end do
    end do
  end subroutine fmg_restore_u
  subroutine fmg_make_linklist
    use mpilib
    use fmg_data
    use grid, only : Gidmin, Gidmax, GidListMax, GidList, ParentGid, ParentRank, ChildGid, ChildRank, encode_surf, NeighborGid, Nei&
&ghborRank
    use boundary, only : touch_boundary, NBOUNDARY
    integer :: amrlev, n, gid, lri, lrj, lrk, lr, ndir, pgid, prank, cgid, crank, ngid, nrank, npgid, nprank, rank, glistmax
    integer,dimension(Gidmin:Gidmax,0:400 -1) :: revlink
    integer,dimension(Gidmin:Gidmax) :: myrevlink
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    integer,dimension(:),allocatable :: glist
    revlink(:,:) = Undefi
    myrevlink(:) = Undefi
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do n = Gidmin, GidListMax(amrlev)
          gid = GidList(n, amrlev) 
          myrevlink(gid) = n 
       end do
    end do
    call mpi_allgather(myrevlink, size(myrevlink), MPI_INTEGER, revlink, size(myrevlink), MPI_INTEGER, MPI_COMM_WORLD, ierr )
    do rank = 0, 400 -1
       do amrlev = AMR_LevelMin, AMR_LevelMax
          myrank = get_myrank()
          if (rank == myrank) glistmax = GidListMax( amrlev )
          call mpi_bcast( glistmax, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
          allocate(glist(Gidmin:glistmax+Gidmin))
          if (rank == myrank) glist = GidList(lbound(glist,1):ubound(glist,1), amrlev)
          call mpi_bcast( glist, size(glist), MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)
          do n = Gidmin, glistmax 
             gid = glist(n) 
             if ( amrlev > AMR_LevelMin ) then
                pgid = ParentGid(gid,rank)
                prank = ParentRank(gid,rank)
                Ancestry(amrlev,rank)%Block(n)%ParentGid = revlink(pgid,prank)
                Ancestry(amrlev,rank)%Block(n)%ParentRank = prank
             endif
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
             do ndir = 0, 2
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
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1)+1, ubound(Ancestry,1)
          do n = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             pgid = Ancestry(amrlev,rank)%Block(n)%ParentGid
             prank = Ancestry(amrlev,rank)%Block(n)%ParentRank
             do ndir = 0, 2
                do lr = Left, Right
                   if ( Ancestry(amrlev,rank)%Block(n)%NeighborSameLevel(lr,ndir) ) cycle
                   npgid = Ancestry(amrlev-1,prank)%Block(pgid)%NeighborGid(lr,ndir)
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
    do amrlev = lbound(Geom,1), ubound(Geom,1)
       do n = lbound(Geom(amrlev)%Block,1), ubound(Geom(amrlev)%Block,1)
          call touch_boundary( GidList(n, amrlev), bool_touch )
          do ndir = 0, 2
             do lr = Left, Right
                Geom(amrlev)%Block(n)%TouchBoundary(lr,ndir) = bool_touch(encode_surf(lr,ndir))
             enddo
          enddo
       enddo
    enddo
  end subroutine fmg_make_linklist
  subroutine fmg_check_linklist
    use grid, only : encode_surf, swap_LR, decode_LR
    integer :: rank, amrlev, fmglev, gid, pgid, prank, cgid, crank, ngid, nrank, lri, lrj, lrk, lr, lrs, ndir
    logical :: bool, bl
    myrank = get_myrank()
    do fmglev = FMG_LevelMin, FMG_LevelMax
       do amrlev = AMR_LevelMin, AMR_LevelMax
          if (ubound(Ancestry(amrlev,myrank)%Block,1) /= ubound(GridLevel(amrlev,fmglev)%Block,1)) &
               print *, '**** error'
          if (lbound(Ancestry(amrlev,myrank)%Block,1) /= lbound(GridLevel(amrlev,fmglev)%Block,1)) &
               print *, '**** error'
       end do
    end do
    do rank = 0, 400 -1
       do amrlev = AMR_LevelMin, AMR_LevelMax
          if (ubound(Ancestry(amrlev,myrank)%Block,1) /= ubound(Ancestry(amrlev,rank)%Block,1) ) &
               print *, '**** error'
          if (lbound(Ancestry(amrlev,myrank)%Block,1) /= lbound(Ancestry(amrlev,rank)%Block,1) ) &
               print *, '**** error'
       end do
    end do
    bool = .true.
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1)+1, ubound(Ancestry,1)
          do gid = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             pgid = Ancestry(amrlev,rank)%Block(gid)%ParentGid
             prank = Ancestry(amrlev,rank)%Block(gid)%ParentRank
             bl = .false.
             do lrk = Left ,Right
                do lrj = Left ,Right
                   do lri = Left ,Right
                      if (&
                           gid == Ancestry(amrlev-1,prank)%Block(pgid)%ChildGid(lri,lrj,lrk) .and. &
                           rank == Ancestry(amrlev-1,prank)%Block(pgid)%ChildRank(lri,lrj,lrk)) &
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
    bool = .TRUE.
    do rank = lbound(Ancestry,2), ubound(Ancestry,2)
       do amrlev = lbound(Ancestry,1), ubound(Ancestry,1)-1
          do gid = lbound(Ancestry(amrlev,rank)%Block,1), ubound(Ancestry(amrlev,rank)%Block,1)
             do ndir = 0, 2
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
                   else 
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
          do ndir = 0, 2
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
  subroutine fmg_poisson
    FMG_PDE_TYPE = FMG_PDE_TYPE_POISSON_EQUATION
    call fmg_data_init 
    call fmg_cycle
  end subroutine fmg_poisson
  subroutine fmg_cycle
    use mpilib
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_converge
    use fmg_ghostcell
    use vmg
    use mg
    use io_util
    use io
    use string
    real(kind=8),parameter :: errormax = 1.d-3
    integer :: nfmgmax
    real(kind=8) :: resh2max, ratio_resh2max, ratio_resh2max_prev
    integer :: fmglev, n
    logical :: bool_converge
    logical :: isNotFinite
    call fmg_init
    call fmg_alloc_LoL 
    call fmg_prepare_data 
    call vmg_init(FMG_LevelMax) 
    call mg_init(FMG_LevelMax) 
    fmglev = FMG_LevelMin
    call fmg_alloc_arr(fmglev, IPSI)
    call fmg_alloc_arr(fmglev, ISRC)
    call fmg_alloc_f(fmglev)
    call fmg_copy(fmglev, IPSI, IU) 
    call fmg_copy(fmglev, ISRC, IRHO) 
    call fmg_converge_c2p(fmglev, ISRC)
    call fmg_converge_c2p(fmglev, IPSI)
    call fmg_boundary_physical_grid(IPSI, ISRC)
    call fmg_boundary_u(fmglev,IPSI)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then 
       resh2max = fmg_get_h2absmax(FMG_LevelMin, ISRC)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       resh2max = fmg_get_absmax(FMG_LevelMin, ISRC)
    else
       resh2max = fmg_get_h2absmax(FMG_LevelMin, ISRC)
    endif
    ratio_resh2max_prev = HUGE(ratio_resh2max_prev)
    nfmgmax = 10
    if (abs(resh2max) <= TINY(resh2max)) nfmgmax = 0
    bool_converge = .true.
    do n = 1, nfmgmax
       call fmg_resid(fmglev, IRHO, IPSI, ISRC)
       call fmg_lin 
       call fmg_add(fmglev, IPSI, IU)
       ratio_resh2max = Resh2maxg/resh2max
       call print_msg( 'fmg error = ' // trim(num2char(ratio_resh2max)) // ' (PDE_TYPE = ' // trim(num2char(FMG_PDE_TYPE)) // ')' )
       if ( ratio_resh2max < errormax ) exit
       if (n == nfmgmax) bool_converge = .false.
       if (ratio_resh2max > ratio_resh2max_prev * 0.5d0 ) bool_converge = .false.
       if ( .not. bool_converge ) exit
       ratio_resh2max_prev = ratio_resh2max
    end do
    if ((.not. bool_converge) .and. (FMG_PDE_TYPE /= FMG_PDE_TYPE_POISSON_EQUATION)) then
       Dt_refineRatio = Dt_refineRatio * 2
       N_fmg_cycle = N_fmg_cycle * 2
       call fmg_set_dtime(fmg_get_dtime()/Dt_refineRatio)
       call print_msg( '** fmg can not converge residual.' )
       call print_msg( '** Dtime = ' // trim(num2char(fmg_get_dtime())) // ' (Dt_refiementRatio = ' // trim(num2char(Dt_refineRatio&
&)) // ')' )
       call mg_finalize
       call vmg_finalize
       call fmg_finalize
       return
    endif
    N_fmg_cycle = N_fmg_cycle + 1
    call fmg_copy(fmglev, IU, IPSI, noskip=.true.) 
    call fmg_converge_c2p(fmglev,IU)
    call fmg_boundary_u(fmglev, IU)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call fmg_psi2g(IU)
       call fmg_restore_u((/17/), IU)
    endif
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
  end subroutine fmg_cycle
  subroutine fmg_init
    allocate( Resmaxg(AMR_LevelMin:AMR_LevelMax) )
    allocate( Resmaxl(AMR_LevelMin:AMR_LevelMax) )
    allocate( Trerr(AMR_LevelMin:AMR_LevelMax) )
  end subroutine fmg_init
  subroutine fmg_finalize
    use fmg_interpol_cubic
    use fmg_converge
    use fmg_ghostcell
    use fmg_reflux
    deallocate(Resmaxg)
    deallocate(Resmaxl)
    deallocate(Trerr)
    call fmg_interp_cubic_finalize
    call fmg_converge_finalize
    call fmg_ghostcell_finalize
    call fmg_reflux_finalize
    call fmg_dealloc_LoL
  end subroutine fmg_finalize
  subroutine fmg_lin
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use string
    INTEGER :: Ncycle, Npre, Npost
    integer :: fmglev, jcycle, jpre, jpost, vlevel
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       Ncycle=2
       Npre=2
       Npost=2
    else
       Ncycle=2
       Npre=2
       Npost=2
    endif
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA)
       end if
    enddo
    fmglev = FMG_LevelMax
    call fmg_alloc_arr(fmglev, IU)
    call fmg_alloc_arr(fmglev, IRHS)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_f(fmglev)
    call fmg_fill0(fmglev, IU) 
    call fmg_slvsml(fmglev, IU, IRHO)
    do fmglev = FMG_LevelMax-1, FMG_LevelMin, -1
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_f(fmglev)
       call fmg_interp(fmglev, IU, IU, cubic=.TRUE.) 
       call fmg_copy(fmglev, IRHS, IRHO)
       do jcycle = 1, Ncycle 
          do vlevel = fmglev, FMG_LevelMax-1 
             do jpre=1,Npre
                call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
             enddo
             call fmg_resid(vlevel, IRES, IU, IRHS, boundary_fill0=.true.)
             call fmg_rstrct(vlevel+1, IRHS, IRES)
             call fmg_fill0(vlevel+1, IU)
             FmgLevel_fill0 = vlevel+1
          enddo
          vlevel = FMG_LevelMax
          call fmg_slvsml(vlevel, IU, IRHS)
          do vlevel = FMG_LevelMax-1, fmglev, -1 
             call fmg_addint(vlevel, IU, IU, IRES)
             do jpost = 1, Npost
                call fmg_relax(vlevel, IU, IRHS, boundary_fill0=.true.)
             enddo
          enddo
       enddo
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
  end subroutine fmg_lin
  subroutine fmg_nonlin
    use vmg
    use mg
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use io
    use string
    use mpilib
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    integer :: jcycle, jpre, jpost, vlevel, fmglev
    integer :: icode
    real(kind=8),parameter :: ALPHA = 1.d0/3.d0
    call fmg_init
    call fmg_alloc_LoL 
    call fmg_prepare_data 
    call vmg_init(FMG_LevelMax) 
    call mg_init(FMG_LevelMax) 
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       call fmg_alloc_arr(fmglev, IU) 
       call fmg_rstrct(fmglev, IU, IU)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA, cubic=.TRUE.)
       end if
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
          do icode = IDOD, IDAD
             call fmg_alloc_arr(fmglev, icode)
             call fmg_rstrct(fmglev, icode, icode)
             call fmg_converge_c2p(fmglev, icode)
             call fmg_boundary_u(fmglev, icode)
             call fmg_ghostcell_fix(fmglev, icode, cubic=.TRUE.)
          end do
       end if
    enddo
    fmglev = FMG_LevelMax
    call fmg_alloc_arr(fmglev, IU)
    call fmg_alloc_arr(fmglev, IRHS)
    call fmg_alloc_arr(fmglev, IRES)
    call fmg_alloc_f(fmglev)
    call fmg_slvsml(fmglev, IU, IRHO, fmg=.TRUE.)
    do fmglev = FMG_LevelMax-1, FMG_LevelMin, -1
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_arr(fmglev+1, IRUF)
       call fmg_alloc_f(fmglev)
       call fmg_interp(fmglev, IU, IU, cubic=.TRUE.) 
       call fmg_copy(fmglev, IRHS, IRHO)
       do jcycle = 1, Ncycle 
          do vlevel = fmglev, FMG_LevelMax-1 
             do jpre=1,Npre
                call fmg_relax(vlevel, IU, IRHS)
             enddo
             call fmg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) 
             call fmg_copy(vlevel+1, IRUF, IU) 
             call fmg_rstrct(vlevel+1, IRHS, IRHS) 
             call fmg_tau(vlevel+1, IRES, IU) 
             call fmg_add(vlevel+1, IRHS, IRES) 
             if (vlevel == fmglev) then 
                call fmg_max_forAMRLevel(Trerr, vlevel+1, IRES, absolute=.TRUE., skip=.TRUE.)
                Trerr = ALPHA*Trerr
             end if
          enddo
          vlevel = FMG_LevelMax
          call fmg_slvsml(vlevel, IU, IRHS)
          do vlevel = FMG_LevelMax-1, fmglev, -1 
             call fmg_sub(vlevel+1, IRES, IU, IRUF) 
             call fmg_interp(vlevel, IRES, IRES, cubic=.TRUE.) 
             call fmg_add(vlevel, IU, IRES) 
             do jpost = 1, Npost
                call fmg_relax(vlevel, IU, IRHS)
             enddo
          enddo
          call fmg_print_errlog(jcycle)
          if (all(Resmaxg < Trerr)) exit
       enddo
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

  end subroutine fmg_nonlin
  subroutine fmg_nonlin_vcycle
    use vmg
    use mg
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary_phys
    use io_util
    use io
    use string
    use mpilib
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    integer :: jcycle, jpre, jpost, vlevel, fmglev
    integer :: icode
    real(kind=8),parameter :: ALPHA = 1.d0/3.d0
    call fmg_init
    call fmg_alloc_LoL 
    call fmg_prepare_data 
    call vmg_init(FMG_LevelMax) 
    call mg_init(FMG_LevelMax) 
    do fmglev = FMG_LevelMin+1, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IRHO)
       call fmg_rstrct(fmglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call fmg_alloc_arr(fmglev, IETA)
          call fmg_rstrct(fmglev, IETA, IETA)
          call fmg_converge_c2p(fmglev, IETA)
          call fmg_boundary_u(fmglev, IETA)
          call fmg_ghostcell_fix(fmglev, IETA, tricubic=.TRUE.)
       end if
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
          do icode = IDOD, IDAD
             call fmg_alloc_arr(fmglev, icode)
             call fmg_rstrct(fmglev, icode, icode)
             call fmg_converge_c2p(fmglev, icode)
             call fmg_boundary_u(fmglev, icode)
             call fmg_ghostcell_fix(fmglev, icode, tricubic=.TRUE.)
          end do
       end if
    enddo
    do fmglev = FMG_LevelMin, FMG_LevelMax
       call fmg_alloc_arr(fmglev, IU)
       call fmg_alloc_arr(fmglev, IRHS)
       call fmg_alloc_arr(fmglev, IRUF)
       call fmg_alloc_arr(fmglev, IRES)
       call fmg_alloc_f(fmglev)
    end do
    call fmg_copy(FMG_LevelMin, IRHS, IRHO) 
    do jcycle = 1, Ncycle 
       do vlevel = FMG_LevelMin, FMG_LevelMax-1 
          do jpre=1,Npre
             call fmg_relax(vlevel, IU, IRHS)
          enddo
          call fmg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) 
          call fmg_copy(vlevel+1, IRUF, IU) 
          call fmg_rstrct(vlevel+1, IRHS, IRHS) 
          call fmg_tau(vlevel+1, IRES, IU) 
          call fmg_add(vlevel+1, IRHS, IRES) 
          if (vlevel == FMG_LevelMin) then 
             call fmg_max_forAMRLevel(Trerr, vlevel+1, IRES, absolute=.TRUE., skip=.TRUE.)
             Trerr = ALPHA*Trerr
          end if
       enddo
       vlevel = FMG_LevelMax
       call fmg_slvsml(vlevel, IU, IRHS)
       do vlevel = FMG_LevelMax-1, FMG_LevelMin, -1 
          call fmg_sub(vlevel+1, IRES, IU, IRUF) 
          call fmg_interp(vlevel, IRES, IRES, cubic=.TRUE.) 
          call fmg_add(vlevel, IU, IRES) 
          do jpost = 1, Npost
             call fmg_relax(vlevel, IU, IRHS)
          enddo
       enddo
       if (all(Resmaxg < Trerr)) exit
       if (jcycle == Ncycle) &
            call print_msg('***** fmg_nonlin was not converged')
    enddo
    call fmg_converge_c2p(FMG_LevelMin, IU)
    call mg_finalize
    call vmg_finalize
    call fmg_finalize
    call dumpdata
    write(*,*) '*** OK'
 call mpi_barrier(MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
 stop

  end subroutine fmg_nonlin_vcycle
  subroutine fmg_slvsml(fmglev, ju, jrhs, fmg)
    use vmg
    integer,intent(IN) :: fmglev, ju, jrhs
    logical,optional :: fmg
    logical :: bool_fmg
    integer :: n
    bool_fmg = .FALSE. 
 if (present(fmg)) bool_fmg = fmg
    if (FMG_PDE_LINEAR) then
       call vmg_fas(ju, jrhs)
    else
       if (bool_fmg) then
          call vmg_fas_fmg(ju, jrhs)
       else
          call vmg_fas(ju, jrhs)
       end if
    end if
  end subroutine fmg_slvsml
  subroutine fmg_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    integer,intent(IN) :: fmglev, jres, ju, jrhs
    logical,optional :: boundary_fill0
    logical :: bool_boundary_fill0
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call fmg_poisson_resid(fmglev, jres, ju, jrhs, boundary_fill0)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    else
       print *, '*** fmg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_resid
  subroutine fmg_rstrct(fmglev, iuc, iuf, cubic)
    use fmg_ghostcell
    use restriction
    integer,intent(IN) :: fmglev, iuc, iuf
    logical,intent(IN),optional :: cubic
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, gid, la, if, jf, kf, ic, jc, kc, m, &
         ics, jcs, kcs, ice, jce, kce
    real(kind=8) :: rdv
    logical :: bool_cubic
    bool_cubic = .FALSE.
    if (present(cubic)) bool_cubic = cubic
    lf = fmglev-1
    lc = fmglev
    call fmg_get_gridsize(lc, ics, jcs, kcs, ice, jce, kce)
    if (bool_cubic) call fmg_ghostcell_fix(lf,iuf,tricubic=cubic)
    do la = AMR_LevelMin, AMR_LevelMax
       rdv = fmg_get_dv(la,lf)/fmg_get_dv(la,lc)
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)
          if ( fmg_skip_grid(gid, la, fmglev) ) cycle
          uf => fmg_get_arrp(la,lf,gid,iuf)
          uc => fmg_get_arrp(la,lc,gid,iuc)
          call rstrct( uc, uf, &
               ics, jcs, kcs, ice, jce, kce, &
               ics, jcs, kcs, cubic)
       end do
    end do
  end subroutine fmg_rstrct
  subroutine fmg_interp(fmglev, iuf, iuc, cubic)
    use fmg_converge
    use fmg_ghostcell
    use fmg_boundary
    use fmg_boundary_phys
    use fmg_interpol_cubic
    use io_util
    integer,intent(IN) :: fmglev, iuf, iuc
    logical,intent(IN),optional :: cubic
    logical :: bool_cubic
    integer :: lf, lc
    integer,parameter :: MINBLOCKSIZE = 4
    lf = fmglev
    lc = fmglev+1
    bool_cubic = .FALSE. 
    if (present(cubic)) bool_cubic = cubic
    call fmg_converge_c2p(lc, iuc, cubic) 
    call fmg_boundary_u(lc, iuc)
    call fmg_ghostcell_fix(lc, iuc, tricubic=bool_cubic)
    call fmg_boundary_extrap(lc, iuc)
    call fmg_boundary_u(lc,iuc)
    if (bool_cubic) then
       call fmg_interp_cubic(fmglev, iuf, iuc)
    else
       call fmg_interp_linear(fmglev, iuf, iuc)
    end if
    call fmg_boundary_fill0(lc, iuc)
    call fmg_boundary_fill0(lf, iuf)
    call fmg_boundary_u(lc,iuc)
    call fmg_boundary_u(lf,iuf)
  end subroutine fmg_interp
  subroutine fmg_interp_linear(fmglev, iuf, iuc)
    use interpolation
    integer,intent(IN) :: fmglev, iuf, iuc
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, gid, la, &
         ifs, jfs, kfs, ife, jfe, kfe
    lf = fmglev
    lc = fmglev+1
    call fmg_get_gridsize(lf, ifs, jfs, kfs, ife, jfe, kfe)
    do la = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)
          if ( fmg_skip_grid(gid, la, lf) ) cycle
          uf => fmg_get_arrp(la,lf,gid,iuf)
          uc => fmg_get_arrp(la,lc,gid,iuc)
          call interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, ifs, jfs, kfs)
       end do
    end do
  end subroutine fmg_interp_linear
  subroutine fmg_fill0(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    real(kind=8),parameter :: zero = 0.d0
    integer :: amrlev, gid
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          u => fmg_get_arrp(amrlev, fmglev, gid, icode)
          u(:,:,:,:) = zero
       end do
    end do
  end subroutine fmg_fill0
  subroutine fmg_relax(fmglev, ju, jrhs, boundary_fill0)
    use mpilib
    integer,intent(IN) :: fmglev, ju, jrhs
    logical,optional :: boundary_fill0
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call fmg_poisson_relax(fmglev, ju, jrhs, Resh2maxg, boundary_fill0)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** fmg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_relax
  subroutine fmg_copy(fmglev, iout, iin, noskip)
    integer,intent(IN) :: fmglev, iout, iin
    logical,optional :: noskip
    logical :: bool_skip
    integer :: amrlev, gid
    real(kind=8),pointer,dimension(:,:,:,:) :: ain, aout
    bool_skip = .TRUE.
 if (present(noskip)) bool_skip = .not. noskip
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( bool_skip .and. fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          aout => fmg_get_arrp(amrlev,fmglev,gid,iout)
          ain => fmg_get_arrp(amrlev,fmglev,gid,iin)
          aout = ain
       end do
    end do
  end subroutine fmg_copy
  subroutine fmg_sub(fmglev, iout, ia, ib)
    integer,intent(IN) :: fmglev, iout, ia, ib
    integer :: amrlev, gid
    real(kind=8),pointer,dimension(:,:,:,:) :: out, a, b
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          out => fmg_get_arrp(amrlev,fmglev,gid,iout)
          a => fmg_get_arrp(amrlev,fmglev,gid,ia)
          b => fmg_get_arrp(amrlev,fmglev,gid,ib)
          out = a - b
       end do
    end do
  end subroutine fmg_sub
  subroutine fmg_add(fmglev, ia, ida)
    integer,intent(IN) :: fmglev, ia, ida
    integer :: amrlev, gid
    real(kind=8),pointer,dimension(:,:,:,:) :: a, da
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          a => fmg_get_arrp(amrlev,fmglev,gid,ia)
          da => fmg_get_arrp(amrlev,fmglev,gid,ida)
          a = a + da
       end do
    end do
  end subroutine fmg_add
  subroutine fmg_addint(fmglev, juf, juc, jres)
    integer,intent(IN) :: fmglev, juf, juc, jres
    call fmg_interp(fmglev, jres, juc)
    call fmg_add(fmglev, juf, jres)
  end subroutine fmg_addint
  subroutine fmg_tau(fmglev, jtau, ju)
    integer,intent(IN) :: fmglev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** fmg_tau: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine fmg_tau
  subroutine fmg_shift(fmglev, icodenew, icode)
    integer,intent(IN) :: fmglev, icodenew, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u, unew
    integer :: i,j,k,is,js,ks,ie,je,ke,if,jf,kf,amrlev,gid, m
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          u => fmg_get_arrp(amrlev,fmglev,gid,icode)
          unew => fmg_get_arrp(amrlev,fmglev,gid,icodenew)
          do m = lbound(u,4), ubound(u,4)
             do k = ks, ke
                do j = js, je
                   do i = is, ie
                      if = i-1
                      jf = j-1
                      kf = k-1
                      unew(i,j,k,m) = u(if,jf,kf,m)
                   end do
                enddo
             enddo
          enddo
       end do
    end do
  end subroutine fmg_shift
  subroutine fmg_fill0_ghcell(fmglev, icode)
    integer,intent(IN) :: fmglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u, unew
    integer :: amrlev,gid
    real(kind=8) :: foo = 0.d0
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          u => fmg_get_arrp(amrlev,fmglev,gid,icode)
          u(lbound(u,1),:,:,:) = foo
          u(ubound(u,1),:,:,:) = foo
          u(:,lbound(u,2),:,:) = foo
          u(:,ubound(u,2),:,:) = foo
          u(:,:,lbound(u,3),:) = foo
          u(:,:,ubound(u,3),:) = foo
       end do
    end do
  end subroutine fmg_fill0_ghcell
  function fmg_get_h2absmax(fmglev, icode) result(normmax)
    use mpilib
    integer,intent(IN) :: icode
    real(kind=8) :: normmax
    real(kind=8) :: h2, normmaxl
    real(kind=8),pointer,dimension(:,:,:,:) :: arr
    integer :: amrlev, gid, fmglev
    integer :: is, ie, js, je, ks, ke
    normmaxl = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       h2 = fmg_get_h(amrlev, fmglev) **2
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          arr => fmg_get_arrp(amrlev,fmglev,gid,icode)
          normmaxl = max(normmaxl, maxval(abs(arr(is:ie,js:je,ks:ke,:)))*h2)
       end do
    end do
    call mpi_allreduce(normmaxl, normmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  end function fmg_get_h2absmax
  function fmg_get_absmax(fmglev, icode) result(normmax)
    use mpilib
    integer,intent(IN) :: fmglev, icode
    real(kind=8) :: normmax
    real(kind=8) :: normmaxl
    real(kind=8),pointer,dimension(:,:,:,:) :: arr
    integer :: amrlev, gid
    integer :: is, ie, js, je, ks, ke
    normmaxl = 0.d0
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do amrlev = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          if ( fmg_skip_grid(gid, amrlev, fmglev) ) cycle
          arr => fmg_get_arrp(amrlev,fmglev,gid,icode)
          normmaxl = max(normmaxl, maxval(abs(arr(is:ie,js:je,ks:ke,:))))
       end do
    end do
    call mpi_allreduce(normmaxl, normmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  end function fmg_get_absmax
  subroutine fmg_max_forAMRLevel(arrmax, fmglev, icode, absolute, skip)
    use mpilib
    real(kind=8),dimension(AMR_LevelMin:),intent(OUT) :: arrmax
    integer,intent(IN) :: fmglev, icode
    logical,intent(IN),optional :: absolute, skip
    real(kind=8),dimension(lbound(arrmax,1):ubound(arrmax,1)) :: arrmaxlocal
    integer :: amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call fmg_max(arrmaxlocal(amrlev), amrlev, fmglev, icode, absolute, skip, mpireduce=.FALSE.)
    end do
    call mpi_allreduce(arrmaxlocal, arrmax, size(arrmax), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  end subroutine fmg_max_forAMRLevel
  subroutine fmg_print_errlog(jcycle)
    use io_util
    use string
    integer,intent(IN) :: jcycle
    integer :: amrlev
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call print_msg( 'error in FMG = ' // trim(num2char(Resmaxg(amrlev))) // ' / ' // trim(num2char(Trerr(amrlev))) &
            // ' at amrlevel = '//trim(num2char(amrlev)) &
            // ', jcycle = '//trim(num2char(jcycle)) &
            )
    end do
  end subroutine fmg_print_errlog
  subroutine fmg_print_prepostlog(prepost, vlevel, jcycle)
    use io_util
    use string
    character(len=*),intent(IN) :: prepost
    integer,intent(IN) :: vlevel, jcycle
    integer,dimension(1) :: maxlev
    maxlev = maxloc(Resmaxg)-1+lbound(Resmaxg,1)
    call print_msg( 'FMG error in '//trim(prepost)//'-smooth = '// trim(num2char(maxval(Resmaxg))) &
         // ' max at amrlev = '//trim(num2char(maxlev(1))) &
         // ', vlevel = '//trim(num2char(vlevel)) &
         // ', jcycle = '//trim(num2char(jcycle)) &
         )
  end subroutine fmg_print_prepostlog
end module fmg
