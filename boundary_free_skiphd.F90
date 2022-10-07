
! Free boundary condition
!
#include "config.h"
module boundary
  use grid , only : Lmin, Lmax
  implicit none
  private
  integer,parameter :: NBOUNDARY =  2*(MZ-MX+1) ! 境界面の個数（３次元では６面）
  integer,parameter :: NIL = 0                  ! 方向コード
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  logical,save :: BoolInitialized = .FALSE.     ! このモジュールが初期化されているか
  ! グリッドの配置
  integer,save :: Igmin(Lmin:Lmax), Jgmin(Lmin:Lmax), Kgmin(Lmin:Lmax)
  integer,save :: Igmax(Lmin:Lmax), Jgmax(Lmin:Lmax), Kgmax(Lmin:Lmax)
  public :: boundary_u, boundary_g, touch_boundary, NBOUNDARY
contains
  !--------------------------------------------------------------------
  ! mirror boundary for z=0, and zmax planes
  !--------------------------------------------------------------------
  subroutine boundary_u(id, stepmode)
    use grid
    use modelParameter
    use unit
    integer,intent(IN) :: id, stepmode
    integer :: i,j,k,m
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: u

    real(kind=DBL_KIND) :: yhni, yh2i, yeli, yhpi, yhmi, yh2pi, ycoi, yHe
    real(kind=DBL_KIND) :: rho0_bound, cs0_bound

    call boundary_init
    call touch_boundary( id, bool_touch )



    if (.not. any(bool_touch) ) return
    u => get_Up(id)

    ! i=imin, 自由境界
    if ( bool_touch(NIL) ) then
       do i = Imingh,Imin-1
          u(i,:,:,:) = u(Imin,:,:,:)

       enddo
    endif

    ! i=imax, 自由境界
    if ( bool_touch(NIR) ) then
       do i = Imax+1,Imaxgh
          u(i,:,:,:) = u(Imax,:,:,:)


       enddo
    endif


    ! j=jmin, 自由境界
    if ( bool_touch(NJL) ) then
       do j = Jmingh,Jmin-1
          u(:,j,:,:) = u(:,Jmin,:,:)

       enddo
    endif
    ! j=jmax, 自由境界
    if ( bool_touch(NJR) ) then
       do j = Jmax+1,Jmaxgh
          u(:,j,:,:) = u(:,Jmax,:,:)

       enddo
    endif

    ! k=kmin, 自由境界
    if ( bool_touch(NKL) ) then
       do k = Kmingh,Kmin-1
          u(:,:,k,:) = u(:,:,Kmin,:)

       enddo
    endif
    ! k=kmax, 自由境界
    if ( bool_touch(NKR) ) then
       do k = Kmax+1,Kmaxgh
          u(:,:,k,:) = u(:,:,Kmax,:)

       enddo
    endif
  end subroutine boundary_u
  !--------------------------------------------------------------------
  ! boundary for gravity
  !--------------------------------------------------------------------
  subroutine boundary_g(gid)
    use grid
    integer,intent(IN) :: gid

  end subroutine boundary_g
  !--------------------------------------------------------------------
  ! グリッドが物理境界に接していれば真。そでなければ偽
  !--------------------------------------------------------------------
  subroutine touch_boundary( id, bool_touch )
    use grid
    integer,intent(IN) :: id
    logical,dimension(0:NBOUNDARY-1),intent(OUT) :: bool_touch
    integer :: level, ig, jg, kg
    call boundary_init
    bool_touch(:) = .FALSE.
    level = get_level(id)
    if ( Igrid(id) == Igmin( level ) ) bool_touch(NIL) = .TRUE.
    if ( Igrid(id) == Igmax( level ) ) bool_touch(NIR) = .TRUE.
    if ( Jgrid(id) == Jgmin( level ) ) bool_touch(NJL) = .TRUE.
    if ( Jgrid(id) == Jgmax( level ) ) bool_touch(NJR) = .TRUE.
    if ( Kgrid(id) == Kgmin( level ) ) bool_touch(NKL) = .TRUE.
    if ( Kgrid(id) == Kgmax( level ) ) bool_touch(NKR) = .TRUE.
  end subroutine touch_boundary
  !--------------------------------------------------------------------
  ! initialized boundary module
  !--------------------------------------------------------------------
  subroutine boundary_init
    use grid
    use io_util
    integer,parameter :: ni0 = NGI_BASE, nj0 = NGJ_BASE, nk0 = NGK_BASE
    integer :: level
    if ( BoolInitialized ) return
    call print_msg( 'initialize boundary' )
    BoolInitialized = .TRUE.
    Igmin(:) =  0
    Jgmin(:) =  0
    Kgmin(:) =  0
    do level = Lmin, Lmax
       Igmax(level) = ni0*2**(level-Lmin) -1
#ifdef Emulate_1Dim
       Jgmax(level) = Jgmin(level)
#else !Emulate_1Dim
       Jgmax(level) = nj0*2**(level-Lmin) -1
#endif !Emulate_1Dim
#ifdef EMULATE_2DIM
       Kgmax(level) = Kgmin(level)
#else !EMULATE_2DIM
       Kgmax(level) = nk0*2**(level-Lmin) -1
#endif !EMULATE_2DIM
    enddo
  end subroutine boundary_init

end module boundary
