
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
          u(i,:,:,MVX) = 1.d6 / Unit_v
          u(i,:,:,MVY) = 0.d0
          u(i,:,:,MVZ) =0.d0
          u(i,:,:,MRHO) = MP_N0 * cgs_amu * MP_mu / Unit_rho

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
           end if

          ! yhni=1.d-8
          ! yhpi=1.d-8
          ! yhmi=1.d-20
          ! yh2pi=1.d-20
          ! yh2i= (1.d0 - (yhni+yhpi+yhmi) - 2*yh2pi)/2.d0
          ! yeli=yhpi - yhmi + yh2pi
          ! ycoi = 0.927d-4
           yHe = 9.7222222d-2 

           rho0_bound = MP_N0 * cgs_amu * MP_mu / Unit_rho
           cs0_bound = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
                /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))/ Unit_v
        
           u(i,:,:,MP) = rho0_bound * cs0_bound * cs0_bound

       enddo
    endif

    ! i=imax, 自由境界
    if ( bool_touch(NIR) ) then
       do i = Imax+1,Imaxgh
          u(i,:,:,:) = u(Imax,:,:,:)

          do j = Jmin, Jmax

             do k = Kmin, Kmax

               if ( u(i,j,k,MVX) < 0.d0 ) then

                  u(i,j,k,MVX) = - u(i,j,k,MVX)

               end if

              end do

           end do

       enddo
    endif


    ! j=jmin, 自由境界
    if ( bool_touch(NJL) ) then
       do j = Jmingh,Jmin-1
          u(:,j,:,:) = u(:,Jmin,:,:)

          do i = Imin, Imax

             do k = Kmin, Kmax

               if ( u(i,j,k,MVY) > 0.d0 ) then

                   u(i,j,k,MVY) = -u(i,j,k,MVY)

               end if

              end do

           end do

       enddo
    endif
    ! j=jmax, 自由境界
    if ( bool_touch(NJR) ) then
       do j = Jmax+1,Jmaxgh
          u(:,j,:,:) = u(:,Jmax,:,:)

          do i = Imin, Imax

             do k = Kmin, Kmax

               if ( u(i,j,k,MVY) < 0.d0 ) then

                   u(i,j,k,MVY) = -u(i,j,k,MVY)

               end if

              end do

           end do

       enddo
    endif

    ! k=kmin, 自由境界
    if ( bool_touch(NKL) ) then
       do k = Kmingh,Kmin-1
          u(:,:,k,:) = u(:,:,Kmin,:)

          do i = Imin, Imax

             do j = Jmin, Jmax

               if ( u(i,j,k,MVZ) > 0.d0 ) then

                   u(i,j,k,MVZ) = -u(i,j,k,MVZ)

               end if

              end do

           end do

       enddo
    endif
    ! k=kmax, 自由境界
    if ( bool_touch(NKR) ) then
       do k = Kmax+1,Kmaxgh
          u(:,:,k,:) = u(:,:,Kmax,:)

          do i = Imin, Imax

             do j = Jmin, Jmax

               if ( u(i,j,k,MVZ) < 0.d0 ) then

                   u(i,j,k,MVZ) = -u(i,j,k,MVZ)

               end if

              end do

           end do

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
