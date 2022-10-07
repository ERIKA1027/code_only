
#include "config.h"
!-------------------------------------------------------------------------
! Module for uniform Multigrid
!-------------------------------------------------------------------------
module mg_data
  use fmg_data, only : Undefi, t_block, t_gridsize, Ngh, IRHO, IRHS, IRES, IFLX, IDBG, IU, ISRC, IPSI, IRUF, IETA, IDOD, IDHE, IDAD, Left, Right, Mmin, Mmax, Nmin, Nmax, fmg_slice3d, fmg_slice3d_f
  implicit none
  integer,save :: FMG_Level = Undefi
  integer,save :: MG_LevelMin = Undefi
  integer,save :: MG_LevelMax = Undefi
  ! MGData(mglev)%U(i,j,k)
  type(t_block),save,dimension(:),pointer :: MGdata => null()
  ! GridSize(mglev)%{Imin,Imax,Jmin,Jmax,Kmin,Kmax}
  type(t_gridsize),save,dimension(:),pointer :: GridSize => null()
  ! TouchBoundary
  logical,save,dimension(Left:Right,MX:MZ) :: TouchBoundary = .FALSE.

  ! for gather (prepare) and scatter (restore)
  ! Buffer (private)
  type t_mgbuffer
     real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: sbuf => null()
     real(kind=DBL_KIND),dimension(:,:,:,:,:),pointer :: rbuf => null()
     integer,dimension(:),pointer :: rcounts => null()
     integer,dimension(:),pointer :: displs => null()
     integer :: ncomp
     logical :: bool_v
  end type t_mgbuffer
  type(t_mgbuffer),save,private :: BufPrep, BufRest
  integer,save,private :: Bsi, Bsj, Bsk, Is, Js, Ks, Ie, Je, Ke ! size of block
#ifdef FMG_LAMBDA
  real(kind=DBL_KIND),parameter :: LAMBDA = FMG_LAMBDA
#endif !FMG_LAMBDA
  ! generic interface
  interface mg_arrp
     module procedure mg_arrp_vect, mg_arrp_scalar
  end interface mg_arrp
  interface mg_fp
     module procedure mg_fp_vect, mg_fp_scalar
  end interface mg_fp
contains
  !-------------------------------------------------------------------------
  ! initialize (call once in a simulation)
  !-------------------------------------------------------------------------
  subroutine mg_data_init(fmglev)
    use mpilib
    use fmg_data, only : fmg_get_gridsize, Geom, AMR_levelMin
    integer,intent(IN) :: fmglev
    integer :: boxsize, mglev, id
    logical,dimension(Left:Right, MX:MZ) :: tBool

    if ( fmglev == FMG_Level ) return

    FMG_Level = fmglev

    ! define min and max of level
    MG_LevelMin = 0
    boxsize = min( &
         ((NGI_BASE)*(NI))/2**FMG_Level, &
         ((NGJ_BASE)*(NJ))/2**FMG_Level, &
         ((NGK_BASE)*(NK))/2**FMG_Level )
    MG_LevelMax = log(dble(boxsize))/log(2.d0) + 0.5 ! boxsize .... 1

    ! データ格納リスト
    allocate( MGdata(MG_LevelMin:MG_LevelMax) )
    ! データサイズ
    allocate( GridSize(MG_LevelMin:MG_LevelMax) )
    do mglev = MG_LevelMin, MG_LevelMax
       GridSize(mglev)%Imin = 0
       GridSize(mglev)%Jmin = 0
       GridSize(mglev)%Kmin = 0
       GridSize(mglev)%Imax = ((NGI_BASE)*(NI))/2**FMG_Level /2**mglev -1
       GridSize(mglev)%Jmax = ((NGJ_BASE)*(NJ))/2**FMG_Level /2**mglev -1
       GridSize(mglev)%Kmax = ((NGK_BASE)*(NK))/2**FMG_Level /2**mglev -1
    enddo

    ! touch physical boundary (FMG_data から継承される)
    tBool(:,:) = .FALSE.
    myrank = get_myrank()
    do id = lbound(Geom(AMR_levelMin)%Block, 1), ubound(Geom(AMR_levelMin)%Block, 1) ! for all base grid
       tBool(:,:) = tBool(:,:) .or. Geom(AMR_levelMin)%Block(id)%TouchBoundary(:,:)
    enddo
    call mpi_allreduce( tBool, TouchBoundary, size(tBool), MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )

    ! blocksize
    call fmg_get_gridsize(FMG_Level, Is,Js,Ks,Ie,Je,Ke)
    Bsi = Ie - Is + 1
    Bsj = Je - Js + 1
    Bsk = Ke - Ks + 1
  end subroutine mg_data_init
  !-------------------------------------------------------------------------
  ! initialize buffer
  ! INPUT : ncomp = number of components to be transfer
  !         buf   = a object of buffer
  !-------------------------------------------------------------------------
  subroutine mg_buffer_init ( ncomp, buf )
    use mpilib
    use grid, only : GidBase, RankBase
    integer,intent(IN) :: ncomp
    type(t_mgbuffer),intent(INOUT) :: buf
    integer :: blocksize, rank, i, j, k
    if ( associated( buf%sbuf )) return
    ! gather, scatter 用の送受信バッファ
    blocksize = Bsi * Bsj * Bsk * ncomp
    allocate( buf%sbuf(Bsi,Bsj,Bsk, ncomp, count( RankBase == get_myrank() ) ) )! 送信バッファ
    allocate( buf%rbuf(Bsi,Bsj,Bsk, ncomp, size(GidBase)) ) ! 受信バッファ
    ! Rcounts
    allocate(buf%rcounts(0:NPE-1))
    buf%rcounts(:) = 0
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             buf%rcounts(RankBase(i,j,k)) = buf%rcounts(RankBase(i,j,k)) + blocksize
          enddo
       enddo
    enddo
    ! Displs
    allocate(buf%displs(0:NPE-1))
    buf%displs(0) = 0
    do rank = 1, NPE -1
       buf%displs(rank) = buf%displs(rank-1) + buf%rcounts(rank-1)
    enddo
    ! MPI_GATHER or MPI_GATHERV
    if ( maxval(buf%rcounts) == minval(buf%rcounts) ) then
       buf%bool_v = .false.
    else
       buf%bool_v = .true.
    endif
  end subroutine mg_buffer_init
  !-------------------------------------------------------------------------
  ! finalize
  ! Base grid が固定の場合には、実行しないほうが良い(オーバーヘッドのため)。
  !-------------------------------------------------------------------------
#define FREE_POINTER_SAFE(ARRAY) \
             if (associated(ARRAY)) then ;\
                deallocate(ARRAY) ;\
                nullify(ARRAY) ;\
             endif

#define FREE_POINTER(ARRAY) \
             deallocate(ARRAY) ;\
             nullify(ARRAY)

  subroutine mg_data_finalize
    use mpilib
    integer :: mglev
    FMG_Level = Undefi
    do mglev = MG_LevelMin, MG_levelMax
       FREE_POINTER_SAFE(MGdata(mglev)%Dbg)
       FREE_POINTER_SAFE(MGdata(mglev)%U)
       FREE_POINTER_SAFE(MGdata(mglev)%Res)
       FREE_POINTER_SAFE(MGdata(mglev)%Rhs)
       FREE_POINTER_SAFE(MGdata(mglev)%Rho)
       FREE_POINTER_SAFE(MGdata(mglev)%Src)
       FREE_POINTER_SAFE(MGdata(mglev)%Psi)
       FREE_POINTER_SAFE(MGdata(mglev)%Ruf)
       FREE_POINTER_SAFE(MGdata(mglev)%Eta)
       FREE_POINTER_SAFE(MGdata(mglev)%Dod)
       FREE_POINTER_SAFE(MGdata(mglev)%Dhe)
       FREE_POINTER_SAFE(MGdata(mglev)%Dad)
       FREE_POINTER_SAFE(MGdata(mglev)%F)
    end do
    FREE_POINTER_SAFE(MGdata)
    FREE_POINTER_SAFE(GridSize)
    FREE_POINTER_SAFE(BufPrep%rbuf)
    FREE_POINTER_SAFE(BufPrep%sbuf)
    FREE_POINTER_SAFE(BufPrep%rcounts)
    FREE_POINTER_SAFE(BufPrep%displs)
    FREE_POINTER_SAFE(BufRest%rbuf)
    FREE_POINTER_SAFE(BufRest%sbuf)
    FREE_POINTER_SAFE(BufRest%rcounts)
    FREE_POINTER_SAFE(BufRest%displs)
  end subroutine mg_data_finalize
  !-------------------------------------------------------------------------
  ! data parepare
  ! iumg (IN) ...... Index of array for MG data (Input data)
  ! iufmg (IN) ..... Index of array for FMG data (Output data)
  !-------------------------------------------------------------------------
  subroutine mg_dataPrepare(iumg, iufmg)
    use mpilib
    use fmg_data
    use grid, only : GidBase, RankBase
    integer,dimension(:),intent(IN) :: iufmg, iumg
    integer :: i, j, k, rank, n, m, mp, ncomp
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, uf

    if (size(iufmg,1) /= size(iumg,1)) then
       print *, '*** mg_dataPrepare: error. iufmg and iumg are not consistent.', size(iufmg,1),size(iumg,1)
       stop
    endif

    ! total number of components of buffer
    ncomp = 0
    do m = lbound(iumg,1), ubound(iumg,1)
       ncomp = ncomp + fmg_getNcomp(iumg(m))
    end do
    call mg_buffer_init( ncomp, BufPrep )

    myrank = get_myrank()
    ! バッファに push する
#define SZ Is:Ie,Js:Je,Ks:Ke,:
    n = 1
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             if ( RankBase(i,j,k) == myrank ) then
                mp = lbound(BufPrep%sbuf,4) ! pointer of m
                do m = lbound(iufmg,1), ubound(iufmg,1)
                   uf => fmg_get_arrp(AMR_LevelMin, FMG_Level, GidBase(i,j,k), iufmg(m))
                   ncomp = size(uf,4)
                   BufPrep%sbuf(:,:,:,mp:mp+ncomp-1,n) = uf(SZ)
                   mp = mp + ncomp
                end do
                n = n + 1
             end if
          end do
       end do
    end do
    if ( BufPrep%bool_v ) then
       call mpi_gatherv( &
            BufPrep%sbuf, size(BufPrep%sbuf),              MPI_DOUBLE_PRECISION, &
            BufPrep%rbuf, BufPrep%rcounts, BufPrep%displs, MPI_DOUBLE_PRECISION, &
            PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    else
       call mpi_gather(&
            BufPrep%sbuf, size(BufPrep%sbuf), MPI_DOUBLE_PRECISION, &
            BufPrep%rbuf, size(BufPrep%sbuf), MPI_DOUBLE_PRECISION, &
            PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    endif
    ! バッファからpopする
    if (myrank == PRIMARY_RANK) then
       mp = lbound(BufPrep%rbuf,4) ! pointer of m
       do m = lbound(iumg,1), ubound(iumg,1)
          call mg_alloc_arr(MG_LevelMin, iumg(m))
          u => mg_get_arrp(MG_LevelMin, iumg(m))
          ncomp = size(u,4)
          n = 1
          do rank = 0, NPE -1
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
  !-------------------------------------------------------------------------
  ! push data back to the FMG data
  !-------------------------------------------------------------------------
  subroutine mg_dataRestore(iufmg, iumg)
    use mpilib
    use fmg_data
    use grid, only : GidBase, RankBase
    integer,dimension(:),intent(IN) :: iufmg, iumg
    integer :: i, j, k, n, rank, mp, ncomp, m
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: u, uf

    if (size(iufmg,1) /= size(iumg,1)) then
       print *, '*** mg_dataRestore: error. iufmg and iumg are not consistent.', size(iufmg,1),size(iumg,1)
       stop
    endif

    ! total number of components of buffer
    ncomp = 0
    do m = lbound(iumg,1), ubound(iumg,1)
       ncomp = ncomp + fmg_getNcomp(iumg(m))
    end do
    call mg_buffer_init( ncomp, BufRest )

    ! データをプライマリランクから各ノードに送る。
    myrank = get_myrank()
    if (myrank == PRIMARY_RANK) then
       mp = lbound(BufRest%rbuf,4)
       do m = lbound(iumg,1), ubound(iumg,1)
          u => mg_get_arrp(MG_LevelMin, iumg(m))
          ncomp = size(u, 4)
          n = 1
          do rank = 0, NPE -1
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
            BufRest%sbuf, size(BufRest%sbuf),      MPI_DOUBLE_PRECISION, &
            PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    else
       call mpi_scatter( &
            BufRest%rbuf, size(BufRest%sbuf), MPI_DOUBLE_PRECISION, &
            BufRest%sbuf, size(BufRest%sbuf), MPI_DOUBLE_PRECISION, &
            PRIMARY_RANK, MPI_COMM_WORLD, ierr)
    endif
    ! ポップする
#define SZ Is:Ie,Js:Je,Ks:Ke,:
    n = 1
    do k = lbound(GidBase,3), ubound(GidBase,3)
       do j = lbound(GidBase,2), ubound(GidBase,2)
          do i = lbound(GidBase,1), ubound(GidBase,1)
             if ( RankBase(i,j,k) == myrank ) then
                mp = lbound(BufRest%sbuf,4) ! pointer of m
                do m = lbound(iufmg,1), ubound(iufmg,1)
                   uf => fmg_get_arrp(AMR_LevelMin, FMG_Level, GidBase(i,j,k), iufmg(m))
                   ncomp = size(uf,4)
                   uf( SZ ) = BufRest%sbuf(:,:,:,mp:mp+ncomp-1,n)
                   mp = mp + ncomp
                end do
                n = n + 1
             endif
          end do
       end do
    end do
  end subroutine mg_dataRestore
  !-------------------------------------------------------------------------
  ! get grid size
  !-------------------------------------------------------------------------
  function mg_get_imingh(mglev) result(imin)
    integer,intent(IN) :: mglev
    integer :: imin
    imin = GridSize(mglev)%Imin - Ngh
  end function mg_get_imingh
  !-------------------------------------------------------------------------
  function mg_get_jmingh(mglev) result(jmin)
    integer,intent(IN) :: mglev
    integer :: jmin
    jmin = GridSize(mglev)%Jmin - Ngh
  end function mg_get_jmingh
  !-------------------------------------------------------------------------
  function mg_get_kmingh(mglev) result(kmin)
    integer,intent(IN) :: mglev
    integer :: kmin
    kmin = GridSize(mglev)%Kmin - Ngh
  end function mg_get_kmingh
  !-------------------------------------------------------------------------
  function mg_get_imaxgh(mglev) result(imax)
    integer,intent(IN) :: mglev
    integer :: imax
    imax = GridSize(mglev)%Imax + Ngh
  end function mg_get_imaxgh
  !-------------------------------------------------------------------------
  function mg_get_jmaxgh(mglev) result(jmax)
    integer,intent(IN) :: mglev
    integer :: jmax
    jmax = GridSize(mglev)%Jmax + Ngh
  end function mg_get_jmaxgh
  !-------------------------------------------------------------------------
  function mg_get_kmaxgh(mglev) result(kmax)
    integer,intent(IN) :: mglev
    integer :: kmax
    kmax = GridSize(mglev)%Kmax + Ngh
  end function mg_get_kmaxgh
  !-------------------------------------------------------------------------
  ! get grid size without ghost cell
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! get grid size with ghost cell
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! allocation
  !-------------------------------------------------------------------------
  subroutine mg_alloc_arr(mglev, icode)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    integer :: imax, jmax, kmax, imin, jmin, kmin
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
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
  !-------------------------------------------------------------------------
  ! allocate F for all block
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! set pointer
  !-------------------------------------------------------------------------
  subroutine mg_set_arrp(arrp, mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    select case(icode)
    case (IDBG)
       MGdata(mglev)%Dbg => arrp
    case (IU)
       MGdata(mglev)%U   => arrp
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
  !-------------------------------------------------------------------------
  ! get pointer
  !-------------------------------------------------------------------------
  function mg_get_arrp(mglev, icode) result(arrp)
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
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
  !-------------------------------------------------------------------------
  subroutine mg_arrp_vect(mglev, icode, arrp)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrp
    arrp => mg_get_arrp(mglev, icode)
    if (.not. fmg_isVector(icode)) then
       print *, '*** error in mg_arrp_vect: icode is scalar', icode
       stop
    endif
  end subroutine mg_arrp_vect
  !-------------------------------------------------------------------------
  subroutine mg_arrp_scalar(mglev, icode, arrp)
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, icode
    real(kind=DBL_KIND),pointer,dimension(:,:,:) :: arrp
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: arrpv
    arrpv => mg_get_arrp(mglev, icode)
    arrp => fmg_slice3d(arrpv(:,:,:,Mmin))
    if (fmg_isVector(icode)) then
       print *, '*** error in mg_arrp_scalar: icode is vector', icode
       stop
    endif
  end subroutine mg_arrp_scalar
  !-------------------------------------------------------------------------
  ! get pointer of flux
  !-------------------------------------------------------------------------
  function mg_get_fp(mglev) result(fp)
    integer,intent(IN) :: mglev
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fp
    fp => MGdata( mglev )%F
  end function mg_get_fp
  !-------------------------------------------------------------------------
  subroutine mg_fp_vect(mglev, fp)
    integer,intent(IN) :: mglev
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fp
    fp => mg_get_fp(mglev)
  end subroutine mg_fp_vect
  !-------------------------------------------------------------------------
  subroutine mg_fp_scalar(mglev, fp)
    integer,intent(IN) :: mglev
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:) :: fp
    real(kind=DBL_KIND),pointer,dimension(:,:,:,:,:) :: fpv
    fpv => mg_get_fp(mglev)
    fp => fmg_slice3d_f(fpv(:,:,:,:,Mmin))
  end subroutine mg_fp_scalar
  !-------------------------------------------------------------------------
  ! get h, ds, dv
  !-------------------------------------------------------------------------
  function mg_get_h(mglev) result(h)
    use fmg_data, only : fmg_get_h, AMR_LevelMin
    integer,intent(IN) ::  mglev
    real(kind=DBL_KIND) :: h
    h = fmg_get_h(AMR_LevelMin, FMG_Level)*2**mglev
  end function mg_get_h
  !-------------------------------------------------------------------------
  function mg_get_ds(mglev) result(ds)
    integer,intent(IN) ::  mglev
    real(kind=DBL_KIND) :: ds
    ds = mg_get_h(mglev)**2
  end function mg_get_ds
  !-------------------------------------------------------------------------
  function mg_get_dv(mglev) result(dv)
    integer,intent(IN) ::  mglev
    real(kind=DBL_KIND) :: dv
    dv = mg_get_h(mglev)**3
  end function mg_get_dv


end module mg_data
