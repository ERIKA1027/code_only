module boundary
  use grid , only : Lmin, Lmax
  implicit none
  private
  integer,parameter :: NBOUNDARY = 2*(2 -0 +1) 
  integer,parameter :: NIL = 0 
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  logical,save :: BoolInitialized = .FALSE. 
  integer,save :: Igmin(Lmin:Lmax), Jgmin(Lmin:Lmax), Kgmin(Lmin:Lmax)
  integer,save :: Igmax(Lmin:Lmax), Jgmax(Lmin:Lmax), Kgmax(Lmin:Lmax)
  public :: boundary_u, boundary_g, touch_boundary, NBOUNDARY
contains
  subroutine boundary_u(id, stepmode)
    use grid
    use modelParameter
    use unit
    integer,intent(IN) :: id, stepmode
    integer :: i,j,k,m
    logical,dimension(0:NBOUNDARY-1) :: bool_touch
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8) :: yhni, yh2i, yeli, yhpi, yhmi, yh2pi, ycoi, yHe
    real(kind=8) :: rho0_bound, cs0_bound
    call boundary_init
    call touch_boundary( id, bool_touch )
    if (.not. any(bool_touch) ) return
    u => get_Up(id)
    if ( bool_touch(NIL) ) then
       do i = Imingh,Imin-1
          u(i,:,:,:) = u(Imin,:,:,:)
          u(i,:,:,1) = 2.d6 / Unit_v
          u(i,:,:,2) = 0.d0
          u(i,:,:,3) =0.d0
          u(i,:,:,0) = MP_N0 * cgs_amu * MP_mu / Unit_rho
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
           end if
           yHe = 9.7222222d-2
           rho0_bound = MP_N0 * cgs_amu * MP_mu / Unit_rho
           cs0_bound = sqrt(cgs_kb*MP_T0/(cgs_amu*MP_mu &
                /(yhni+yh2i+yeli+yhpi+yhmi+yh2pi+yHe+MP_frac_COsum)))/ Unit_v
           u(i,:,:,4) = rho0_bound * cs0_bound * cs0_bound
       enddo
    endif
    if ( bool_touch(NIR) ) then
       do i = Imax+1,Imaxgh
          u(i,:,:,:) = u(Imax,:,:,:)
          do j = Jmin, Jmax
             do k = Kmin, Kmax
               if ( u(i,j,k,1) < 0.d0 ) then
                  u(i,j,k,1) = - u(i,j,k,1)
               end if
              end do
           end do
       enddo
    endif
    if ( bool_touch(NJL) ) then
       do j = Jmingh,Jmin-1
          u(:,j,:,:) = u(:,Jmin,:,:)
          do i = Imin, Imax
             do k = Kmin, Kmax
               if ( u(i,j,k,2) > 0.d0 ) then
                   u(i,j,k,2) = -u(i,j,k,2)
               end if
              end do
           end do
       enddo
    endif
    if ( bool_touch(NJR) ) then
       do j = Jmax+1,Jmaxgh
          u(:,j,:,:) = u(:,Jmax,:,:)
          do i = Imin, Imax
             do k = Kmin, Kmax
               if ( u(i,j,k,2) < 0.d0 ) then
                   u(i,j,k,2) = -u(i,j,k,2)
               end if
              end do
           end do
       enddo
    endif
    if ( bool_touch(NKL) ) then
       do k = Kmingh,Kmin-1
          u(:,:,k,:) = u(:,:,Kmin,:)
          do i = Imin, Imax
             do j = Jmin, Jmax
               if ( u(i,j,k,3) > 0.d0 ) then
                   u(i,j,k,3) = -u(i,j,k,3)
               end if
              end do
           end do
       enddo
    endif
    if ( bool_touch(NKR) ) then
       do k = Kmax+1,Kmaxgh
          u(:,:,k,:) = u(:,:,Kmax,:)
          do i = Imin, Imax
             do j = Jmin, Jmax
               if ( u(i,j,k,3) < 0.d0 ) then
                   u(i,j,k,3) = -u(i,j,k,3)
               end if
              end do
           end do
       enddo
    endif
  end subroutine boundary_u
  subroutine boundary_g(gid)
    use grid
    integer,intent(IN) :: gid
  end subroutine boundary_g
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
  subroutine boundary_init
    use grid
    use io_util
    integer,parameter :: ni0 = 8, nj0 = 8, nk0 = 8
    integer :: level
    if ( BoolInitialized ) return
    call print_msg( 'initialize boundary' )
    BoolInitialized = .TRUE.
    Igmin(:) = 0
    Jgmin(:) = 0
    Kgmin(:) = 0
    do level = Lmin, Lmax
       Igmax(level) = ni0*2**(level-Lmin) -1
       Jgmax(level) = nj0*2**(level-Lmin) -1
       Kgmax(level) = nk0*2**(level-Lmin) -1
    enddo
  end subroutine boundary_init
end module boundary
