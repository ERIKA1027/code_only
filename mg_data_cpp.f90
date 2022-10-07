module mg_data
  use fmg_data, only : Undefi, t_block, t_gridsize, Ngh, IRHO, IRHS, IRES, IFLX, IDBG, IU, ISRC, IPSI, IRUF, IETA, IDOD, IDHE, IDAD&
&, Left, Right, Mmin, Mmax, Nmin, Nmax, fmg_slice3d, fmg_slice3d_f
  implicit none
  integer,save :: FMG_Level = Undefi
  integer,save :: MG_LevelMin = Undefi
  integer,save :: MG_LevelMax = Undefi
  type(t_block),save,dimension(:),pointer :: MGdata => null()
  type(t_gridsize),save,dimension(:),pointer :: GridSize => null()
  logical,save,dimension(Left:Right,0:2) :: TouchBoundary = .FALSE.
  type t_mgbuffer
     real(kind=8),dimension(:,:,:,:,:),pointer :: sbuf => null()
     real(kind=8),dimension(:,:,:,:,:),pointer :: rbuf => null()
     integer,dimension(:),pointer :: rcounts => null()
     integer,dimension(:),pointer :: displs => null()
     integer :: ncomp
     logical :: bool_v
  end type t_mgbuffer
  type(t_mgbuffer),save,private :: BufPrep, BufRest
  integer,save,private :: Bsi, Bsj, Bsk, Is, Js, Ks, Ie, Je, Ke 
  interface mg_arrp
     module procedure mg_arrp_vect, mg_arrp_scalar
  end interface mg_arrp
  interface mg_fp
     module procedure mg_fp_vect, mg_fp_scalar
  end interface mg_fp
contains
  subroutine mg_data_init(fmglev)
    use mpilib
    use fmg_data, only : fmg_get_gridsize, Geom, AMR_levelMin
    integer,intent(IN) :: fmglev
    integer :: boxsize, mglev, id
    logical,dimension(Left:Right, 0:2) :: tBool
    if ( fmglev == FMG_Level ) return
    FMG_Level = fmglev
    MG_LevelMin = 0
    boxsize = min( &
         ((8)*(8))/2**FMG_Level, &
         ((8)*(8))/2**FMG_Level, &
         ((8)*(8))/2**FMG_Level )
    MG_LevelMax = log(dble(boxsize))/log(2.d0) + 0.5 
    allocate( MGdata(MG_LevelMin:MG_LevelMax) )
    allocate( GridSize(MG_LevelMin:MG_LevelMax) )
    do mglev = MG_LevelMin, MG_LevelMax
       GridSize(mglev)%Imin = 0
       GridSize(mglev)%Jmin = 0
       GridSize(mglev)%Kmin = 0
       GridSize(mglev)%Imax = ((8)*(8))/2**FMG_Level /2**mglev -1
       GridSize(mglev)%Jmax = ((8)*(8))/2**FMG_Level /2**mglev -1
       GridSize(mglev)%Kmax = ((8)*(8))/2**FMG_Level /2**mglev -1
    enddo
    tBool(:,:) = .FALSE.
    myrank = get_myrank()
    do id = lbound(Geom(AMR_levelMin)%Block, 1), ubound(Geom(AMR_levelMin)%Block, 1) 
       tBool(:,:) = tBool(:,:) .or. Geom(AMR_levelMin)%Block(id)%TouchBoundary(:,:)
    enddo
    call mpi_allreduce( tBool, TouchBoundary, size(tBool), MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
    call fmg_get_gridsize(FMG_Level, Is,Js,Ks,Ie,Je,Ke)
    Bsi = Ie - Is + 1
    Bsj = Je - Js + 1
    Bsk = Ke - Ks + 1
  end subroutine mg_data_init
  subroutine mg_buffer_init ( ncomp, buf )
    use mpilib
    use grid, only : GidBase, RankBase
    integer,intent(IN) :: ncomp
    type(t_mgbuffer),intent(INOUT) :: buf
    integer :: blocksize, rank, i, j, k
    if ( associated( buf%sbuf )) return
    blocksize = Bsi * Bsj * Bsk * ncomp
    allocate( buf%sbuf(Bsi,Bsj,Bsk, ncomp, count( RankBase == get_myrank() ) ) )
    allocate( buf%rbuf(Bsi,Bsj,Bsk, ncomp, size(GidBase)) ) 
    allocate(buf%rcounts(0:400 -1))
    buf%rcounts(:) = 0
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             buf%rcounts(RankBase(i,j,k)) = buf%rcounts(RankBase(i,j,k)) + blocksize
          enddo
       enddo
    enddo
    allocate(buf%displs(0:400 -1))
    buf%displs(0) = 0
    do rank = 1, 400 -1
       buf%displs(rank) = buf%displs(rank-1) + buf%rcounts(rank-1)
    enddo
    if ( maxval(buf%rcounts) == minval(buf%rcounts) ) then
       buf%bool_v = .false.
    else
       buf%bool_v = .true.
    endif
  end subroutine mg_buffer_init
  subroutine mg_data_finalize
    use mpilib
    integer :: mglev
    FMG_Level = Undefi
    do mglev = MG_LevelMin, MG_levelMax
       if (associated(MGdata(mglev)%Dbg)) then 
 deallocate(MGdata(mglev)%Dbg) 
 nullify(MGdata(mglev)%Dbg) 
 endif
       if (associated(MGdata(mglev)%U)) then 
 deallocate(MGdata(mglev)%U) 
 nullify(MGdata(mglev)%U) 
 endif
       if (associated(MGdata(mglev)%Res)) then 
 deallocate(MGdata(mglev)%Res) 
 nullify(MGdata(mglev)%Res) 
 endif
       if (associated(MGdata(mglev)%Rhs)) then 
 deallocate(MGdata(mglev)%Rhs) 
 nullify(MGdata(mglev)%Rhs) 
 endif
       if (associated(MGdata(mglev)%Rho)) then 
 deallocate(MGdata(mglev)%Rho) 
 nullify(MGdata(mglev)%Rho) 
 endif
       if (associated(MGdata(mglev)%Src)) then 
 deallocate(MGdata(mglev)%Src) 
 nullify(MGdata(mglev)%Src) 
 endif
       if (associated(MGdata(mglev)%Psi)) then 
 deallocate(MGdata(mglev)%Psi) 
 nullify(MGdata(mglev)%Psi) 
 endif
       if (associated(MGdata(mglev)%Ruf)) then 
 deallocate(MGdata(mglev)%Ruf) 
 nullify(MGdata(mglev)%Ruf) 
 endif
       if (associated(MGdata(mglev)%Eta)) then 
 deallocate(MGdata(mglev)%Eta) 
 nullify(MGdata(mglev)%Eta) 
 endif
       if (associated(MGdata(mglev)%Dod)) then 
 deallocate(MGdata(mglev)%Dod) 
 nullify(MGdata(mglev)%Dod) 
 endif
       if (associated(MGdata(mglev)%Dhe)) then 
 deallocate(MGdata(mglev)%Dhe) 
 nullify(MGdata(mglev)%Dhe) 
 endif
       if (associated(MGdata(mglev)%Dad)) then 
 deallocate(MGdata(mglev)%Dad) 
 nullify(MGdata(mglev)%Dad) 
 endif
       if (associated(MGdata(mglev)%F)) then 
 deallocate(MGdata(mglev)%F) 
 nullify(MGdata(mglev)%F) 
 endif
    end do
    if (associated(MGdata)) then 
 deallocate(MGdata) 
 nullify(MGdata) 
 endif
    if (associated(GridSize)) then 
 deallocate(GridSize) 
 nullify(GridSize) 
 endif
    if (associated(BufPrep%rbuf)) then 
 deallocate(BufPrep%rbuf) 
 nullify(BufPrep%rbuf) 
 endif
    if (associated(BufPrep%sbuf)) then 
 deallocate(BufPrep%sbuf) 
 nullify(BufPrep%sbuf) 
 endif
    if (associated(BufPrep%rcounts)) then 
 deallocate(BufPrep%rcounts) 
 nullify(BufPrep%rcounts) 
 endif
    if (associated(BufPrep%displs)) then 
 deallocate(BufPrep%displs) 
 nullify(BufPrep%displs) 
 endif
    if (associated(BufRest%rbuf)) then 
 deallocate(BufRest%rbuf) 
 nullify(BufRest%rbuf) 
 endif
    if (associated(BufRest%sbuf)) then 
 deallocate(BufRest%sbuf) 
 nullify(BufRest%sbuf) 
 endif
    if (associated(BufRest%rcounts)) then 
 deallocate(BufRest%rcounts) 
 nullify(BufRest%rcounts) 
 endif
    if (associated(BufRest%displs)) then 
 deallocate(BufRest%displs) 
 nullify(BufRest%displs) 
 endif
  end subroutine mg_data_finalize
  subroutine mg_dataPrepare(iumg, iufmg)
    use mpilib
    use fmg_data
    use grid, only : GidBase, RankBase
    integer,dimension(:),intent(IN) :: iufmg, iumg
    integer :: i, j, k, rank, n, m, mp, ncomp
    real(kind=8),pointer,dimension(:,:,:,:) :: u, uf
    if (size(iufmg,1) /= size(iumg,1)) then
       print *, '*** mg_dataPrepare: error. iufmg and iumg are not consistent.', size(iufmg,1),size(iumg,1)
       stop
    endif
    ncomp = 0
    do m = lbound(iumg,1), ubound(iumg,1)
       ncomp = ncomp + fmg_getNcomp(iumg(m))
    end do
    call mg_buffer_init( ncomp, BufPrep )
    myrank = get_myrank()
    n = 1
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             if ( RankBase(i,j,k) == myrank ) then
                mp = lbound(BufPrep%sbuf,4) 
                do m = lbound(iufmg,1), ubound(iufmg,1)
                   uf => fmg_get_arrp(AMR_LevelMin, FMG_Level, GidBase(i,j,k), iufmg(m))
                   ncomp = size(uf,4)
                   BufPrep%sbuf(:,:,:,mp:mp+ncomp-1,n) = uf(Is:Ie,Js:Je,Ks:Ke,:)
                   mp = mp + ncomp
                end do
                n = n + 1
             end if
          end do
       end do
    end do
    if ( BufPrep%bool_v ) then
       call mpi_gatherv( &
            BufPrep%sbuf, size(BufPrep%sbuf), MPI_DOUBLE_PRECISION, &
            BufPrep%rbuf, BufPrep%rcounts, BufPrep%displs, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, ierr)
    else
       call mpi_gather(&
            BufPrep%sbuf, size(BufPrep%sbuf), MPI_DOUBLE_PRECISION, &
            BufPrep%rbuf, size(BufPrep%sbuf), MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, ierr)
    endif
    if (myrank == 0) then
       mp = lbound(BufPrep%rbuf,4) 
       do m = lbound(iumg,1), ubound(iumg,1)
          call mg_alloc_arr(MG_LevelMin, iumg(m))
          u => mg_get_arrp(MG_LevelMin, iumg(m))
          ncomp = size(u,4)
          n = 1
          do rank = 0, 400 -1
             do k = lbound(GidBase,3), ubound(GidBase,3)
                do j = lbound(GidBase,2), ubound(GidBase,2)
                   do i = lbound(GidBase,1), ubound(GidBase,1)
                      if ( RankBase(i,j,k) == rank ) then
                         u( Bsi*i:Bsi*(i+1)-1, Bsj*j:Bsj*(j+1)-1, Bsk*k:Bsk*(k+1)-1,: ) = BufPrep%rbuf(:,:,:,mp:mp+ncomp-1,n)
                         n = n + 1
                      endif
                   end do
                end do
             end do
          enddo
          mp = mp + ncomp
       end do
    endif
  end subroutine mg_dataPrepare
  subroutine mg_dataRestore(iufmg, iumg)
    use mpilib
    use fmg_data
    use grid, only : GidBase, RankBase
    integer,dimension(:),intent(IN) :: iufmg, iumg
    integer :: i, j, k, n, rank, mp, ncomp, m
    real(kind=8),pointer,dimension(:,:,:,:) :: u, uf
    if (size(iufmg,1) /= size(iumg,1)) then
       print *, '*** mg_dataRestore: error. iufmg and iumg are not consistent.', size(iufmg,1),size(iumg,1)
       stop
    endif
    ncomp = 0
    do m = lbound(iumg,1), ubound(iumg,1)
       ncomp = ncomp + fmg_getNcomp(iumg(m))
    end do
    call mg_buffer_init( ncomp, BufRest )
    myrank = get_myrank()
    if (myrank == 0) then
       mp = lbound(BufRest%rbuf,4)
       do m = lbound(iumg,1), ubound(iumg,1)
          u => mg_get_arrp(MG_LevelMin, iumg(m))
          ncomp = size(u, 4)
          n = 1
          do rank = 0, 400 -1
             do k = lbound(GidBase,3), ubound(GidBase,3)
                do j = lbound(GidBase,2), ubound(GidBase,2)
                   do i = lbound(GidBase,1), ubound(GidBase,1)
                      if ( RankBase(i,j,k) == rank ) then
                         BufRest%rbuf(:,:,:,mp:mp+ncomp-1,n) = u( Bsi*i:Bsi*(i+1)-1, Bsj*j:Bsj*(j+1)-1, Bsk*k:Bsk*(k+1)-1, : )
                         n = n + 1
                      end if
                   enddo
                end do
             end do
          end do
          mp = mp + ncomp
       enddo
    endif
    if ( BufRest%bool_v ) then
       call mpi_scatterv( &
            BufRest%rbuf, BufRest%rcounts, BufRest%displs, MPI_DOUBLE_PRECISION, &
            BufRest%sbuf, size(BufRest%sbuf), MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, ierr)
    else
       call mpi_scatter( &
            BufRest%rbuf, size(BufRest%sbuf), MPI_DOUBLE_PRECISION, &
            BufRest%sbuf, size(BufRest%sbuf), MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, ierr)
    endif
    n = 1
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             if ( RankBase(i,j,k) == myrank ) then
                mp = lbound(BufRest%sbuf,4) 
                do m = lbound(iufmg,1), ubound(iufmg,1)
                   uf => fmg_get_arrp(AMR_LevelMin, FMG_Level, GidBase(i,j,k), iufmg(m))
                   ncomp = size(uf,4)
                   uf( Is:Ie,Js:Je,Ks:Ke,: ) = BufRest%sbuf(:,:,:,mp:mp+ncomp-1,n)
                   mp = mp + ncomp
                end do
                n = n + 1
             endif
          end do
       end do
    end do
  end subroutine mg_dataRestore
  function mg_get_imingh(mglev) result(imin)
    integer,intent(IN) :: mglev
    integer :: imin
    imin = GridSize(mglev)%Imin - Ngh
  end function mg_get_imingh
  function mg_get_jmingh(mglev) result(jmin)
    integer,intent(IN) :: mglev
    integer :: jmin
    jmin = GridSize(mglev)%Jmin - Ngh
  end function mg_get_jmingh
  function mg_get_kmingh(mglev) result(kmin)
    integer,intent(IN) :: mglev
    integer :: kmin
    kmin = GridSize(mglev)%Kmin - Ngh
  end function mg_get_kmingh
  function mg_get_imaxgh(mglev) result(imax)
    integer,intent(IN) :: mglev
    integer :: imax
    imax = GridSize(mglev)%Imax + Ngh
  end function mg_get_imaxgh
  function mg_get_jmaxgh(mglev) result(jmax)
    integer,intent(IN) :: mglev
    integer :: jmax
    jmax = GridSize(mglev)%Jmax + Ngh
  end function mg_get_jmaxgh
  function mg_get_kmaxgh(mglev) result(kmax)
    integer,intent(IN) :: mglev
    integer :: kmax
    kmax = GridSize(mglev)%Kmax + Ngh
  end function mg_get_kmaxgh
  subroutine mg_get_gridsize(mglevel, imin,jmin,kmin,imax,jmax,kmax)
    integer,intent(IN) :: mglevel
    integer,intent(OUT) :: imin,jmin,kmin,imax,jmax,kmax
    imin = GridSize(mglevel)%Imin
    jmin = GridSize(mglevel)%Jmin
    kmin = GridSize(mglevel)%Kmin
    imax = GridSize(mglevel)%Imax
    jmax = GridSize(mglevel)%Jmax
    kmax = GridSize(mglevel)%Kmax
  end subroutine mg_get_gridsize
  subroutine mg_get_gridsizeGh(mglevel, imin,jmin,kmin,imax,jmax,kmax)
    integer,intent(IN) :: mglevel
    integer,intent(OUT) :: imin,jmin,kmin,imax,jmax,kmax
    imin = mg_get_imingh(mglevel)
    jmin = mg_get_jmingh(mglevel)
    kmin = mg_get_kmingh(mglevel)
    imax = mg_get_imaxgh(mglevel)
    jmax = mg_get_jmaxgh(mglevel)
    kmax = mg_get_kmaxgh(mglevel)
  end subroutine mg_get_gridsizeGh
  subroutine mg_alloc_arr(mglev, icode)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    integer :: imax, jmax, kmax, imin, jmin, kmin
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    arrp => mg_get_arrp(mglev, icode)
    if ( associated( arrp ) ) return
    imin = mg_get_imingh( mglev )
    imax = mg_get_imaxgh( mglev )
    jmin = mg_get_jmingh( mglev )
    jmax = mg_get_jmaxgh( mglev )
    kmin = mg_get_kmingh( mglev )
    kmax = mg_get_kmaxgh( mglev )
    if (fmg_isVector(icode)) then
       allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmax))
    else
       allocate(arrp(imin:imax,jmin:jmax,kmin:kmax,Mmin:Mmin))
    endif
    call mg_set_arrp(arrp, mglev, icode)
  end subroutine mg_alloc_arr
  subroutine mg_alloc_f( mglev )
    integer,intent(IN) :: mglev
    integer :: imax, jmax, kmax, imin, jmin, kmin, n, m, gid
    if (associated ( MGdata(mglev)%F )) return
    imin = mg_get_imingh( mglev )
    imax = mg_get_imaxgh( mglev )
    jmin = mg_get_jmingh( mglev )
    jmax = mg_get_jmaxgh( mglev )
    kmin = mg_get_kmingh( mglev )
    kmax = mg_get_kmaxgh( mglev )
    allocate( MGdata(mglev)%F(imin:imax,jmin:jmax,kmin:kmax,Nmin:Nmax,Mmin:Mmax))
  end subroutine mg_alloc_f
  subroutine mg_set_arrp(arrp, mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       MGdata(mglev)%Dbg => arrp
    case (IU)
       MGdata(mglev)%U => arrp
    case (IRES)
       MGdata(mglev)%Res => arrp
    case (IRHS)
       MGdata(mglev)%Rhs => arrp
    case (IRHO)
       MGdata(mglev)%Rho => arrp
    case (ISRC)
       MGdata(mglev)%Src => arrp
    case (IPSI)
       MGdata(mglev)%Psi => arrp
    case (IRUF)
       MGdata(mglev)%Ruf => arrp
    case (IETA)
       MGdata(mglev)%Eta => arrp
    case (IDOD)
       MGdata(mglev)%Dod => arrp
    case (IDHE)
       MGdata(mglev)%Dhe => arrp
    case (IDAD)
       MGdata(mglev)%Dad => arrp
    case default
       print *, '*** This code is not supported', icode
       stop
    end select
  end subroutine mg_set_arrp
  function mg_get_arrp(mglev, icode) result(arrp)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       arrp => MGdata(mglev)%Dbg
    case (IU)
       arrp => MGdata(mglev)%U
    case (IRES)
       arrp => MGdata(mglev)%Res
    case (IRHS)
       arrp => MGdata(mglev)%Rhs
    case (IRHO)
       arrp => MGdata(mglev)%Rho
    case (ISRC)
       arrp => MGdata(mglev)%Src
    case (IPSI)
       arrp => MGdata(mglev)%Psi
    case (IRUF)
       arrp => MGdata(mglev)%Ruf
    case (IETA)
       arrp => MGdata(mglev)%Eta
    case (IDOD)
       arrp => MGdata(mglev)%Dod
    case (IDHE)
       arrp => MGdata(mglev)%Dhe
    case (IDAD)
       arrp => MGdata(mglev)%Dad
    case default
       print *, '*** This code is not supported', icode
       stop
    end select
  end function mg_get_arrp
  subroutine mg_arrp_vect(mglev, icode, arrp)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: arrp
    arrp => mg_get_arrp(mglev, icode)
    if (.not. fmg_isVector(icode)) then
       print *, '*** error in mg_arrp_vect: icode is scalar', icode
       stop
    endif
  end subroutine mg_arrp_vect
  subroutine mg_arrp_scalar(mglev, icode, arrp)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:) :: arrp
    real(kind=8),pointer,dimension(:,:,:,:) :: arrpv
    arrpv => mg_get_arrp(mglev, icode)
    arrp => fmg_slice3d(arrpv(:,:,:,Mmin))
    if (fmg_isVector(icode)) then
       print *, '*** error in mg_arrp_scalar: icode is vector', icode
       stop
    endif
  end subroutine mg_arrp_scalar
  function mg_get_fp(mglev) result(fp)
    integer,intent(IN) :: mglev
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fp
    fp => MGdata( mglev )%F
  end function mg_get_fp
  subroutine mg_fp_vect(mglev, fp)
    integer,intent(IN) :: mglev
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fp
    fp => mg_get_fp(mglev)
  end subroutine mg_fp_vect
  subroutine mg_fp_scalar(mglev, fp)
    integer,intent(IN) :: mglev
    real(kind=8),pointer,dimension(:,:,:,:) :: fp
    real(kind=8),pointer,dimension(:,:,:,:,:) :: fpv
    fpv => mg_get_fp(mglev)
    fp => fmg_slice3d_f(fpv(:,:,:,:,Mmin))
  end subroutine mg_fp_scalar
  function mg_get_h(mglev) result(h)
    use fmg_data, only : fmg_get_h, AMR_LevelMin
    integer,intent(IN) :: mglev
    real(kind=8) :: h
    h = fmg_get_h(AMR_LevelMin, FMG_Level)*2**mglev
  end function mg_get_h
  function mg_get_ds(mglev) result(ds)
    integer,intent(IN) :: mglev
    real(kind=8) :: ds
    ds = mg_get_h(mglev)**2
  end function mg_get_ds
  function mg_get_dv(mglev) result(dv)
    integer,intent(IN) :: mglev
    real(kind=8) :: dv
    dv = mg_get_h(mglev)**3
  end function mg_get_dv
end module mg_data
