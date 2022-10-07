module vmg
  use mpilib
  use fmg_data
  use vmg_interpol
  use fmg_converge
  implicit none
  private
  integer,save :: FMG_Level
  integer,save :: Imin, Jmin, Kmin, Imax, Jmax, Kmax
  integer,save :: Imingh, Jmingh, Kmingh, Imaxgh, Jmaxgh, Kmaxgh
  integer,save :: Icmin, Jcmin, Kcmin, Icmax, Jcmax, Kcmax
  real(kind=8),save,allocatable :: Resmaxg(:), Trerr(:)
  logical,save :: CubicInterp
  logical,save :: UseSerialMg
  public :: vmg_fas, vmg_fas_fmg, &
       vmg_init, vmg_finalize, vmg_rstrct, vmg_interp, vmg_relax
contains
subroutine vmg_poisson_tau(amrlev, jtau, ju)
  use fmg_data
  integer,intent(IN) :: amrlev, jtau, ju
  real(kind=8),pointer,dimension(:,:,:) :: tau, u
  real(kind=8),pointer,dimension(:,:,:,:) :: f
  integer :: gid, i, j, k, amrlevf, m
  real(kind=8) :: hi
  amrlevf = amrlev+1 
  if ( amrlevf > AMR_LevelMax ) return
  myrank = get_myrank()
  call vmg_poisson_flux(amrlevf, ju) 
  do gid = fmg_get_gidmin(amrlevf), fmg_get_gidmax(amrlevf)
     call fmg_arrp(amrlevf, FMG_Level, gid, jtau, tau)
     call fmg_arrp(amrlevf, FMG_Level, gid, ju, u)
     call fmg_fp (amrlevf, FMG_Level, gid, f)
     hi = 1.d0/fmg_get_h(amrlevf, FMG_Level)
     tau = 0.d0
     do k = Kmin, Kmax
        do j = Jmin, Jmax
           do i = Imin, Imax
              tau(i,j,k) = &
                   -(f(i,j,k,0)-f(i-1,j,k,0) &
                   + f(i,j,k,1)-f(i,j-1,k,1) &
                   + f(i,j,k,2)-f(i,j,k-1,2))*hi
           end do
        enddo
     enddo
  enddo
  call vmg_rstrct(amrlev, jtau)
  call vmg_poisson_flux(amrlev, ju) 
  do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
     if (.not. fmg_have_child(gid, amrlev)) cycle
     call fmg_arrp(amrlev, FMG_Level, gid, jtau, tau)
     call fmg_arrp(amrlev, FMG_Level, gid, ju, u)
     call fmg_fp (amrlev, FMG_Level, gid, f)
     hi = 1.d0/fmg_get_h(amrlev, FMG_Level)
     do k = Kmin, Kmax
        do j = Jmin, Jmax
           do i = Imin, Imax
              tau(i,j,k) = tau(i,j,k) &
                   +(f(i,j,k,0)-f(i-1,j,k,0) &
                   + f(i,j,k,1)-f(i,j-1,k,1) &
                   + f(i,j,k,2)-f(i,j,k-1,2))*hi
           enddo
        enddo
     enddo
  end do
end subroutine vmg_poisson_tau
subroutine vmg_poisson_flux(amrlev, ju)
  use fmg_boundary_phys
    integer,intent(IN) :: amrlev, ju
    real(kind=8) :: hi
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    real(kind=8),pointer,dimension(:,:,:) :: u
    integer :: gid, ndir, lr, i, j, k
    call vmg_ghostcell(amrlev, ju)
    call vmg_boundary_u(amrlev, FMG_Level, ju)
    hi = 1.d0/fmg_get_h( amrlev, FMG_Level )
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       call fmg_fp(amrlev, FMG_Level, gid, f)
       call fmg_arrp(amrlev, FMG_Level, gid, ju, u)
       do k = Kmin-1, Kmax
          do j = Jmin-1, Jmax
             do i = Imin-1, Imax
                f(i,j,k,0)= (u(i+1,j,k)-u(i,j,k))*hi
                f(i,j,k,1)= (u(i,j+1,k)-u(i,j,k))*hi
                f(i,j,k,2)= (u(i,j,k+1)-u(i,j,k))*hi
             end do
          end do
       end do
    end do
end subroutine vmg_poisson_flux
  subroutine vmg_poisson_relax(amrlev, ju, jrhs)
    use mpilib
    integer,intent(IN) :: amrlev, ju, jrhs
    real(kind=8),parameter :: sixth = 1.d0/6.0d0
    real(kind=8) :: h, resh2, resmax, resmax_g
    real(kind=8),pointer,dimension(:,:,:,:) :: f
    real(kind=8),pointer,dimension(:,:,:) :: u, rhs
    integer :: gid, ndir, i, j, k, ipass
    myrank = get_myrank()
    resmax = 0.d0
    h = fmg_get_h( amrlev, FMG_Level )
    do ipass=1,2 
       call vmg_poisson_flux(amrlev, ju)
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          call fmg_fp(amrlev, FMG_Level, gid, f)
          call fmg_arrp(amrlev, FMG_Level, gid, ju, u)
          call fmg_arrp(amrlev, FMG_Level, gid, jrhs, rhs)
          do k = Kmin, Kmax
             do j = Jmin, Jmax
                do i= Imin + mod(j+k+ipass,2), Imax, 2
                   resh2 = &
                        (f(i,j,k,0)-f(i-1,j,k,0) &
                        +f(i,j,k,1)-f(i,j-1,k,1) &
                        +f(i,j,k,2)-f(i,j,k-1,2) &
                        -rhs(i,j,k) * h &
                        )* h
                   u(i,j,k) = u(i,j,k) + sixth * resh2
                   resmax = max(resmax, abs(resh2))
                enddo
             enddo
          enddo
       enddo
    enddo
    resmax = resmax / h**2
    call mpi_allreduce(resmax, resmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    Resmaxg(amrlev) = resmax_g
  end subroutine vmg_poisson_relax
  subroutine vmg_init(fmglev)
    use util
    integer,intent(IN) :: fmglev
    FMG_Level = fmglev
    Imin = GridSize(FMG_Level)%Imin
    Jmin = GridSize(FMG_Level)%Jmin
    kmin = GridSize(FMG_Level)%Kmin
    Imax = GridSize(FMG_Level)%Imax
    Jmax = GridSize(FMG_Level)%Jmax
    Kmax = GridSize(FMG_Level)%Kmax
    Imingh = fmg_get_imingh(FMG_Level)
    Jmingh = fmg_get_jmingh(FMG_Level)
    Kmingh = fmg_get_kmingh(FMG_Level)
    Imaxgh = fmg_get_imaxgh(FMG_Level)
    Jmaxgh = fmg_get_jmaxgh(FMG_Level)
    Kmaxgh = fmg_get_kmaxgh(FMG_Level)
    Icmin = Imin
    Jcmin = Jmin
    Kcmin = Kmin
    Icmax = ((Imax)-(Imin))/2 + mod(min((Imax)-(Imin),0),2) + (Imin)
    Jcmax = ((Jmax)-(Jmin))/2 + mod(min((Jmax)-(Jmin),0),2) + (Jmin)
    Kcmax = ((Kmax)-(Kmin))/2 + mod(min((Kmax)-(Kmin),0),2) + (Kmin)
    if (Ngh >= 2) then
       CubicInterp = .TRUE.
    else
       CubicInterp = .FALSE.
    end if
    allocate( Resmaxg(AMR_LevelMin:AMR_LevelMax) )
    allocate( Trerr(AMR_LevelMin+1:AMR_LevelMax) )
    call vmg_interp_init(FMG_Level)
    call fmg_converge_init(FMG_Level)
    if ( util_isPowerOfTow(8) .and. &
         util_isPowerOfTow(8) .and. &
         util_isPowerOfTow(8) ) then
       UseSerialMg = .TRUE.
    else
       UseSerialMg = .FALSE.
    end if
  end subroutine vmg_init
  subroutine vmg_finalize
    deallocate(Resmaxg)
    deallocate(Trerr)
    call vmg_interp_finalize
    call fmg_converge_finalize
  end subroutine vmg_finalize
  subroutine vmg_fas(ju, jrho)
    use fmg_boundary
    use fmg_boundary_phys
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    integer :: vmglev, jcycle, jpre, jpost, vlevel
    INTEGER :: Ncycle, Npre, Npost
    real(kind=8),parameter :: ALPHA = 1.d0/3.d0
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       Ncycle=1
       Npre=2
       Npost=2
    else
       Ncycle=1
       Npre=2
       Npost=2
    endif
    call fmg_alloc_arr(FMG_Level, IRHO)
    call fmg_alloc_arr(FMG_Level, IRHS)
    call fmg_alloc_arr(FMG_Level, IRES)
    call fmg_alloc_arr(FMG_Level, IU)
    call fmg_alloc_arr(FMG_Level, IPSI)
    call fmg_alloc_f(FMG_Level)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_alloc_arr(FMG_Level, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call fmg_alloc_arr(FMG_Level, IDOD)
       call fmg_alloc_arr(FMG_Level, IDHE)
       call fmg_alloc_arr(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (jrho /= IRHO) then
       do vmglev = AMR_LevelMin, AMR_LevelMax
          call vmg_copy(vmglev, IRHO, jrho)
       enddo
    endif
    if (ju /= IU) then
       print *, '**** error in mg_dataIF'
       call mpi_finalize(ierr)
       stop
    end if
    call fmg_converge_c2p(FMG_Level, IU)
    call fmg_converge_c2p(FMG_Level, IRHO)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_converge_c2p(FMG_Level, IETA)
       call fmg_boundary_u(FMG_Level, IETA)
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_converge_c2p(FMG_Level, IDOD)
       call fmg_boundary_u(FMG_Level, IDOD)
       call fmg_converge_c2p(FMG_Level, IDHE)
       call fmg_boundary_u(FMG_Level, IDHE)
       call fmg_converge_c2p(FMG_Level, IDAD)
       call fmg_boundary_u(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    call fmg_boundary_fill0(FMG_Level, IU)
    call fmg_boundary_u(FMG_Level, IU)
    do vmglev = AMR_LevelMin, AMR_LevelMax 
       call vmg_copy(vmglev, IRHS, IRHO)
    enddo
    do jcycle = 1, Ncycle 
       do vlevel = AMR_LevelMax, AMR_LevelMin+1, -1 
          call vmg_gridboundary(vlevel, IU)
          do jpre = 1, Npre
             call vmg_relax(vlevel, IU, IRHS)
          enddo
          call vmg_rstrct(vlevel-1, IU, cubic=.TRUE.) 
          call vmg_copy(vlevel-1, IPSI, IU) 
          call vmg_rstrct(vlevel-1, IRHS) 
          call vmg_tau(vlevel-1, IRES, IU) 
          call vmg_add_onov(vlevel-1, IRHS, IRES) 
          Trerr(vlevel) = ALPHA* vmg_get_max(vlevel-1, IRES, absolute=.TRUE., xskip=.TRUE.) 
       enddo
       if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) &
            call vmg_zeroAverage(AMR_LevelMin, IRHS)
       call vmg_slvsml(IU, IRHS)
       do vlevel = AMR_LevelMin+1, AMR_LevelMax 
          call vmg_sub(vlevel-1, IRES, IU, IPSI) 
          call vmg_interp(vlevel, IRES, IRES) 
          call vmg_add(vlevel, IU, IRES) 
          call vmg_gridboundary(vlevel, IU)
          do jpost = 1, Npost
             call vmg_relax(vlevel, IU, IRHS)
          enddo
       enddo
       if (all(Resmaxg(AMR_LevelMin+1:AMR_LevelMax) < Trerr(AMR_LevelMin+1:AMR_LevelMax))) exit
    enddo
  end subroutine vmg_fas
  subroutine vmg_fas_fmg(ju, jrho)
    use fmg_boundary
    use fmg_boundary_phys
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    integer :: vmglev, jcycle, jpre, jpost, vlevel
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    real(kind=8),parameter :: ALPHA = 1.d0/3.d0
    call fmg_alloc_arr(FMG_Level, IRHO)
    call fmg_alloc_arr(FMG_Level, IRHS)
    call fmg_alloc_arr(FMG_Level, IRES)
    call fmg_alloc_arr(FMG_Level, IU)
    call fmg_alloc_arr(FMG_Level, IPSI)
    call fmg_alloc_f(FMG_Level)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_alloc_arr(FMG_Level, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call fmg_alloc_arr(FMG_Level, IDOD)
       call fmg_alloc_arr(FMG_Level, IDHE)
       call fmg_alloc_arr(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (jrho /= IRHO) then
       do vmglev = AMR_LevelMin, AMR_LevelMax
          call vmg_copy(vmglev, IRHO, jrho)
       enddo
    endif
    if (ju /= IU) then
       print *, '**** error in mg_dataIF'
       call mpi_finalize(ierr)
       stop
    end if
    call fmg_converge_c2p(FMG_Level, IU)
    call fmg_converge_c2p(FMG_Level, IRHO)
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call fmg_converge_c2p(FMG_Level, IETA)
       call fmg_boundary_u(FMG_Level, IETA)
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call fmg_converge_c2p(FMG_Level, IDOD)
       call fmg_boundary_u(FMG_Level, IDOD)
       call fmg_converge_c2p(FMG_Level, IDHE)
       call fmg_boundary_u(FMG_Level, IDHE)
       call fmg_converge_c2p(FMG_Level, IDAD)
       call fmg_boundary_u(FMG_Level, IDAD)
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    call fmg_boundary_fill0(FMG_Level, IU)
    call fmg_boundary_u(FMG_Level, IU)
    do vmglev = AMR_LevelMin, AMR_LevelMax
       call vmg_copy(vmglev, IRHS, IRHO)
    enddo
    call vmg_slvsml(IU, IRHO,fmg=.TRUE.)
    do vmglev = AMR_LevelMin+1, AMR_LevelMax
       call vmg_interp(vmglev, IU, IU)
       do jcycle = 1, Ncycle 
          do vlevel = vmglev, AMR_LevelMin+1, -1 
             call vmg_gridboundary(vlevel, IU)
             do jpre = 1, Npre
                call vmg_relax(vlevel, IU, IRHS)
             enddo
             call vmg_rstrct(vlevel-1, IU, cubic=.TRUE.) 
             call vmg_copy(vlevel-1, IPSI, IU) 
             call vmg_rstrct(vlevel-1, IRHS) 
             call vmg_tau(vlevel-1, IRES, IU) 
             call vmg_add_onov(vlevel-1, IRHS, IRES) 
             if (vlevel == vmglev) &
                  Trerr(vlevel) = ALPHA* vmg_get_max(vlevel-1, IRES, absolute=.TRUE., xskip=.TRUE.) 
          enddo
          if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) &
               call vmg_zeroAverage(AMR_LevelMin, IRHS)
          call vmg_slvsml(IU, IRHS)
          do vlevel = AMR_LevelMin+1, vmglev 
             call vmg_sub(vlevel-1, IRES, IU, IPSI) 
             call vmg_interp(vlevel, IRES, IRES) 
             call vmg_add(vlevel, IU, IRES) 
             call vmg_gridboundary(vlevel, IU)
             do jpost = 1, Npost
                call vmg_relax(vlevel, IU, IRHS)
             enddo
             call print_msg( 'error in VMG = ' // trim(num2char(Resmaxg(vlevel))) // ' / ' // trim(num2char(Trerr(vlevel))) &
                  // ' at level = '// trim(num2char(vlevel)) // ', jcycle = '//trim(num2char(jcycle)))
          enddo
          if (all(Resmaxg(AMR_LevelMin+1:AMR_LevelMax) < Trerr(AMR_LevelMin+1:AMR_LevelMax))) exit
       enddo
    end do
  end subroutine vmg_fas_fmg
  subroutine vmg_slvsml(ju, jrhs, fmg)
    use mg_data
    use mg
    use fmg_boundary_phys
    integer,intent(IN) :: ju, jrhs
    logical,optional :: fmg
    logical :: bool_fmg
    integer :: n
    integer,parameter :: N_VMGRELAX = 10
    if (.not. UseSerialMg) then
       do n = 1, N_VMGRELAX
          call vmg_relax( AMR_LevelMin, ju, jrhs)
       end do
       return
    endif
    bool_fmg = .FALSE. 
 if (present(fmg)) bool_fmg = fmg
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_dataPrepare((/ISRC, IPSI/), (/jrhs, ju/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
       call mg_dataPrepare((/ISRC, IPSI, IETA/), (/jrhs, ju, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       call mg_dataPrepare((/ISRC, IPSI, IETA/), (/jrhs, ju, IETA/))
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
       call mg_dataPrepare((/ISRC, IPSI, IDOD, IDHE, IDAD/), (/jrhs, ju, IDOD, IDHE, IDAD/))
    else
       print *, '*** vmg_fas: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    if (get_myrank() == 0) then
       if (FMG_PDE_LINEAR) then
          call mg_cycle(IPSI, ISRC)
       else
          if (bool_fmg) then
             call mg_nonlin(IPSI, ISRC)
          else
             call mg_nonlin_vcycle(IPSI, ISRC)
          end if
       endif
    endif
    call mg_dataRestore((/ju/), (/IPSI/))
  end subroutine vmg_slvsml
  subroutine vmg_tau(amrlev, jtau, ju)
    integer,intent(IN) :: amrlev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call vmg_poisson_tau(amrlev, jtau, ju)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** vmg_tau: this type of PDE is not supported', FMG_PDE_TYPE
    end if
    if (.not. CubicInterp) call vmg_tau_crop(amrlev, jtau) 
  end subroutine vmg_tau
  subroutine vmg_tau_crop(amrlev, jtau)
    use mpilib
    integer,intent(IN) :: amrlev, jtau
    real(kind=8),pointer,dimension(:,:,:,:) :: tau
    integer :: ndir, lr, cgid, lc, lp, crank
    real(kind=8),parameter :: zero = 0.D0
    type(t_ablock),pointer :: cblock 
    myrank = get_myrank()
    lp = amrlev 
    lc = lp + 1 
    do crank = 0, (400)-1
       do cgid = lbound(Ancestry(lc,crank)%Block,1), ubound(Ancestry(lc,crank)%Block,1) 
          cblock => Ancestry(lc,crank)%Block(cgid)
          if ( cblock%NeighborParentLevel(Left, 0) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(Imin,:,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Left, 1) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,Jmin,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Left, 2) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,:,Kmin,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Right, 0) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(Imax,:,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Right, 1) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,Jmax,:,:) = zero
          end if
          if ( cblock%NeighborParentLevel(Right, 2) .and. cblock%ParentRank == myrank ) then
             tau => fmg_get_arrp(lp, FMG_Level, cblock%ParentGid, jtau)
             tau(:,:,Kmax,:) = zero
          end if
       end do
    end do
  end subroutine vmg_tau_crop
  subroutine vmg_sub_onov(amrlev, iout, ia, ib)
    integer,intent(IN) :: amrlev, iout, ia, ib
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: out, a, b
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       out => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       a => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       b => fmg_get_arrp(amrlev,FMG_Level,gid,ib)
       out = a - b
    end do
  end subroutine vmg_sub_onov
  subroutine vmg_sub(amrlev, iout, ia, ib)
    integer,intent(IN) :: amrlev, iout, ia, ib
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: out, a, b
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       out => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       a => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       b => fmg_get_arrp(amrlev,FMG_Level,gid,ib)
       out = a - b
    end do
  end subroutine vmg_sub
  subroutine vmg_add_onov(amrlev, ia, ida)
    integer,intent(IN) :: amrlev, ia, ida
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: a, da
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       a => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       da => fmg_get_arrp(amrlev,FMG_Level,gid,ida)
       a = a + da
    end do
  end subroutine vmg_add_onov
  subroutine vmg_add(amrlev, ia, ida)
    integer,intent(IN) :: amrlev, ia, ida
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: a, da
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       a => fmg_get_arrp(amrlev,FMG_Level,gid,ia)
       da => fmg_get_arrp(amrlev,FMG_Level,gid,ida)
       a = a + da
    end do
  end subroutine vmg_add
  subroutine vmg_rstrct(amrlev, icode, cubic) 
    use fmg_ghostcell
    integer,intent(IN) :: amrlev, icode
    logical,intent(IN),optional :: cubic
    integer :: ndir
    if (present(cubic)) then
       if (cubic) then 
          call fmg_ghfix_parentlev_init(FMG_Level)
          call fmg_ghfix_samelev_init(FMG_Level)
          do ndir = 0, 2
             call fmg_ghfix_parentlev(amrlev+1, FMG_Level, ndir, icode, tricubic=cubic)
             call fmg_ghfix_samelev(amrlev+1, FMG_Level, ndir,icode)
          end do
       end if
    end if
    call fmg_converge_c2p_lev(amrlev,FMG_Level,icode, cubic)
  end subroutine vmg_rstrct
  subroutine vmg_interp(amrlev, juf, juc)
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, juf, juc
    call vmg_ghostcell(amrlev-1, juc)
    call vmg_boundary_extrap(amrlev-1, juc) 
    call vmg_boundary_u(amrlev-1, FMG_Level, juc) 
    call vmg_interp_p2c_lev(amrlev, juf, juc, cubic=CubicInterp)
    call vmg_boundary_fill0(amrlev-1, juc) 
    call vmg_boundary_u(amrlev-1, FMG_Level, juc) 
  end subroutine vmg_interp
  subroutine vmg_cgc(amrlev, ju, jruf, jtmp)
    integer,intent(IN) :: amrlev, ju, jruf, jtmp
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: du, uf
    real(kind=8),pointer,dimension(:,:,:,:) :: tmp, ruf, uc
    call vmg_interp(amrlev, jtmp, ju)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       du => fmg_get_arrp(amrlev, FMG_Level, gid, jtmp)
       uf => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
       uf(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) = uf(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) + du(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) 
    enddo
    call vmg_interp(amrlev, jtmp, jruf)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       du => fmg_get_arrp(amrlev, FMG_Level, gid, jtmp)
       uf => fmg_get_arrp(amrlev, FMG_Level, gid, ju)
       uf(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) = uf(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) - du(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) 
    enddo
  end subroutine vmg_cgc
  subroutine vmg_fill0(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       u => fmg_get_arrp(amrlev, FMG_Level, gid, icode)
       u = 0.d0
    enddo
  end subroutine vmg_fill0
  subroutine vmg_relax(amrlev, ju, jrhs)
    use mpilib
    integer,intent(IN) :: amrlev, ju, jrhs
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call vmg_poisson_relax(amrlev, ju, jrhs)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** vmg_relax: this type of PDE is not supported', FMG_PDE_TYPE
    end if
  end subroutine vmg_relax
  subroutine vmg_copy(amrlev, iout, iin)
    integer,intent(IN) :: amrlev, iout, iin
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: ain, aout
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       aout => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       ain => fmg_get_arrp(amrlev,FMG_Level,gid,iin)
       aout = ain
    end do
  end subroutine vmg_copy
  subroutine vmg_copy_onov(amrlev, iout, iin)
    integer,intent(IN) :: amrlev, iout, iin
    integer :: gid
    real(kind=8),pointer,dimension(:,:,:,:) :: ain, aout
    myrank = get_myrank()
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       if (Ancestry(amrlev,myrank)%Block(gid)%ChildGid(Left,Left,Left) == Undefi) cycle
       aout => fmg_get_arrp(amrlev,FMG_Level,gid,iout)
       ain => fmg_get_arrp(amrlev,FMG_Level,gid,iin)
       aout = ain
    end do
  end subroutine vmg_copy_onov
  subroutine vmg_gridboundary(amrlev, icode)
    use fmg_ghostcell
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, icode
    integer :: ndir
    call vmg_boundary_u(amrlev-1, FMG_Level, icode)
    call fmg_ghfix_parentlev_init(FMG_Level)
    do ndir = 0, 2
       call fmg_ghfix_parentlev(amrlev, FMG_Level, ndir, icode, cubic=CubicInterp)
    enddo
    call vmg_boundary_fill0(amrlev, icode)
    call vmg_boundary_u(amrlev, FMG_Level, icode)
  end subroutine vmg_gridboundary
  subroutine vmg_ghostcell(amrlev, icode)
    use fmg_ghostcell
    use fmg_boundary_phys
    integer,intent(IN) :: amrlev, icode
    integer :: ndir
    call vmg_boundary_u(amrlev, FMG_Level, icode)
    call fmg_ghfix_samelev_init(FMG_Level)
    do ndir = 0, 2
       call fmg_ghfix_samelev(amrlev, FMG_Level, ndir,icode)
    enddo
    call vmg_boundary_fill0(amrlev, icode)
    call vmg_boundary_u(amrlev, FMG_Level, icode)
  end subroutine vmg_ghostcell
  subroutine vmg_boundary_fill0(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    real(kind=8),parameter :: zero = 0.d0
    integer :: id, fmglev
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    fmglev = FMG_Level
    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
          up(Imingh:Imin-1,:,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
          up(Imax+1:Imaxgh,:,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
          up(:,Jmingh:Jmin-1,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
          up(:,Jmax+1:Jmaxgh,:,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
          up(:,:,Kmingh:Kmin-1,:) = zero
       endif
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
          up(:,:,Kmax+1:Kmaxgh,:) = zero
       endif
    enddo
  end subroutine vmg_boundary_fill0
  subroutine vmg_boundary_extrap(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: id, fmglev
    integer :: i, j, k
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    fmglev = FMG_Level
    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
          do i = Imin-1, Imingh, -1
             up(i,:,:,:) = -up(i+1,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
          do i = Imax+1, Imaxgh
             up(i,:,:,:) = -up(i-1,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
          do j = Jmin-1, Jmingh, -1
             up(:,j,:,:) = -up(:,j+1,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
          do j = Jmax+1, Jmaxgh
             up(:,j,:,:) = -up(:,j-1,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
          do k = Kmin-1, Kmingh, -1
             up(:,:,k,:) = -up(:,:,k+1,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
          do k = Kmax+1, Kmaxgh
             up(:,:,k,:) = -up(:,:,k-1,:)
          end do
       end if
    enddo
  end subroutine vmg_boundary_extrap
  subroutine vmg_boundary_extrap_BAK(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    integer :: id, fmglev
    integer :: i, j, k
    real(kind=8),pointer,dimension(:,:,:,:) :: up
    fmglev = FMG_Level
    do id = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       up => fmg_get_arrp(amrlev, fmglev, id, icode)
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,0) ) then
          do i = Imin-1, Imingh, -1
             up(i,:,:,:) = 2*up(i+1,:,:,:) - up(i+2,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,0) ) then
          do i = Imax+1, Imaxgh
             up(i,:,:,:) = 2*up(i-1,:,:,:) - up(i-2,:,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,1) ) then
          do j = Jmin-1, Jmingh, -1
             up(:,j,:,:) = 2*up(:,j+1,:,:) - up(:,j+2,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,1) ) then
          do j = Jmax+1, Jmaxgh
             up(:,j,:,:) = 2*up(:,j-1,:,:) - up(:,j-2,:,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Left,2) ) then
          do k = Kmin-1, Kmingh, -1
             up(:,:,k,:) = 2*up(:,:,k+1,:) - up(:,:,k+2,:)
          end do
       end if
       if ( Geom(amrlev)%Block(id)%TouchBoundary(Right,2) ) then
          do k = Kmax+1, Kmaxgh
             up(:,:,k,:) = 2*up(:,:,k-1,:) - up(:,:,k-2,:)
          end do
       end if
    enddo
  end subroutine vmg_boundary_extrap_BAK
  function vmg_get_max(amrlev, icode, absolute, skip, xskip, mpireduce) result(maxarr)
    integer,intent(IN) :: amrlev, icode
    real(kind=8) :: maxarr
    logical,intent(IN),optional :: absolute, skip, xskip, mpireduce
    call fmg_max(maxarr, amrlev, FMG_Level, icode, absolute, skip, xskip, mpireduce)
  end function vmg_get_max
  subroutine vmg_sum(sumarr, amrlev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, icode
    real(kind=8),intent(OUT) :: sumarr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_sum(sumarr, amrlev, FMG_Level, icode, absolute, skip, mpireduce)
  end subroutine vmg_sum
  subroutine vmg_ave(avearr, amrlev, icode, absolute, skip, mpireduce)
    integer,intent(IN) :: amrlev, icode
    real(kind=8),intent(OUT) :: avearr
    logical,intent(IN),optional :: absolute, skip, mpireduce
    call fmg_ave(avearr, amrlev, FMG_Level, icode, absolute, skip, mpireduce)
  end subroutine vmg_ave
  subroutine vmg_zeroAverage(amrlev, icode)
    integer,intent(IN) :: amrlev, icode
    real(kind=8) :: avearr
    real(kind=8),dimension(:,:,:,:),pointer :: arr
    integer :: fmglev, gid, is,js,ks,ie,je,ke
    fmglev = FMG_Level
    call vmg_ave(avearr, amrlev, icode, skip=.FALSE.)
    call fmg_get_gridsize(fmglev, is,js,ks,ie,je,ke)
    do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
       arr => fmg_get_arrp(amrlev, fmglev, gid, icode)
       arr(is:ie,js:je,ks:ke,Mmin:Mmin) = arr(is:ie,js:je,ks:ke,Mmin:Mmin) - avearr
    enddo
  end subroutine vmg_zeroAverage
  function vmg_get_error(jrhs) result(error)
    use mpilib
    integer,intent(IN) :: jrhs
    real(kind=8) :: error
    real(kind=8) :: rhoave, drho, drhog
    real(kind=8),pointer,dimension(:,:,:,:) :: rho
    real(kind=8),parameter :: eps = 1.D-7
    integer :: amrlev, gid
    error = 0.d0
    do amrlev = AMR_LevelMin, AMR_LevelMax
       call vmg_ave(rhoave, amrlev, jrhs) 
       drho = 0.d0
       do gid = fmg_get_gidmin(amrlev), fmg_get_gidmax(amrlev)
          rho => fmg_get_arrp(amrlev,FMG_Level,gid,jrhs)
          drho = max(drho, maxval(ABS(rho(Imin:Imax,Jmin:Jmax,Kmin:Kmax,:) - rhoave)))
       end do
       call mpi_allreduce(drho, drhog, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       error = max(error, Resmaxg(amrlev)/(drhog*(1.d0+eps)))
    end do
  end function vmg_get_error
end module vmg
