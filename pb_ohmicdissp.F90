#include "config.h"
! #define VERBOSE
! #define TEST_PROBLEM
!-------------------------------------------------------------------------
! Module for Ohmic dissipation with explicit solver for patchBlock
! The explicit solver utilizes a sub-cycle method.
! -------------------------------------------------------------------------
module pb_ohmicdissp
  use patchBlock
  implicit none
  private
  real(KIND=DBL_KIND),parameter :: CFL_OD = (CFL)/3.d0
  integer,dimension(MX:MZ),parameter :: B_MLIST = (/MBX,MBY,MBZ/)
  ! list of patch blocks
  type t_patchBlockList
     type(t_upatch),pointer :: patch
  end type t_patchBlockList
  type(t_patchBlockList),dimension(:),allocatable,save :: Patch_List
  ! rank_list
  integer,dimension(:),allocatable,save :: Rank_list
  ! eta_list
  type t_eta_sphere
     real(kind=DBL_KIND),dimension(MX:MZ) :: center
     real(kind=DBL_KIND) :: radius
     real(kind=DBL_KIND) :: dtaper
     real(kind=DBL_KIND) :: etamax
  end type t_eta_sphere
  type(t_eta_sphere),dimension(:),allocatable,save :: Eta_list
  real(kind=DBL_KIND),parameter :: Ratio_Taper_Radius = 1.d0/4.d0

  ! Number cells of margin for patch block.
  ! NOTE: If eta(i) > 0 and eta(i+1) = 0, B(i+1) changes, indicating that N_MARGIN = 1.
  ! For moving sink particles, it should be N_MARGIN >= 3.
  integer,parameter :: N_MARGIN = 4
  ! Number of ghost cells for a patch block
  integer,parameter :: N_GHOSTCELL=1
  integer,save :: Is, Ie, Js, Je, Ks, Ke, Ms, Me, Isgh, Iegh, Jsgh, Jegh, Ksgh, Kegh
  ! Macros for array size
#define AIJK Is:Ie, Js:Je, Ks:Ke
#define AIJKM Is:Ie, Js:Je, Ks:Ke, Ms:Me
#define AIJKGH Isgh:Iegh, Jsgh:Jegh, Ksgh:Kegh
#define AIJKMGH Isgh:Iegh, Jsgh:Jegh, Ksgh:Kegh, Ms:Me
#define AIJKMNGH Isgh:Iegh, Jsgh:Jegh, Ksgh:Kegh, Ms:Me, MX:MZ
  public :: pbod_prepare_patchList, pbod_for_patchList, pbod_restore_patchList, pbod_sphere
contains

  !-----------------------------------------------------------------------
  ! Execute Ohmic dissipation for spherical regions. Diffusion coeffs eta are uniform
  !-----------------------------------------------------------------------
  subroutine pbod_sphere(center_list, radius_list, etamax_list, gridlevel)
    use overBlockCoordinates
    real(kind=DBL_KIND),dimension(MX:,:),intent(IN) :: center_list
    real(kind=DBL_KIND),dimension(:),intent(IN) :: radius_list, etamax_list
    integer,intent(IN) :: gridlevel
    integer :: n
    real(kind=DBL_KIND),dimension(MX:MZ) :: posL, posR
    type(t_obPointPhys) :: pL, pR
    type(t_obRect),dimension(:),allocatable :: rect_list
    type(t_obRectPhys) :: rectPhys
    allocate(rect_list(lbound(center_list, 2): ubound(center_list, 2)))
    allocate( Eta_list(lbound(center_list, 2): ubound(center_list, 2)))
    do n = lbound(rect_list, 1), ubound(rect_list, 1)

       ! Taper region [raduis - dtaper, raduis + dtaper]
       Eta_list(n)%center(:) = center_list(:,n)
       Eta_list(n)%dtaper = radius_list(n)*Ratio_Taper_Radius
       Eta_list(n)%radius = radius_list(n)+Eta_list(n)%dtaper
       Eta_list(n)%etamax = etamax_list(n)

       posL(:) = Eta_list(n)%center(:) - Eta_list(n)%radius - Eta_list(n)%dtaper
       posR(:) = Eta_list(n)%center(:) + Eta_list(n)%radius + Eta_list(n)%dtaper
       call ob_assignCoordPhysToPointPhys(posL, pL)
       call ob_assignCoordPhysToPointPhys(posR, pR)
       call ob_assignPointPhysPairToRectPhys(pL, pR, rectPhys)
       call ob_RectPhys2RectOb(rectPhys, gridlevel, rect_list(n))
    end do
    call pbod_prepare_patchList(rect_list)
    call pbod_for_patchList
    call pbod_restore_patchList
    deallocate(rect_list, Eta_list)
  end subroutine pbod_sphere
  !-----------------------------------------------------------------------
  ! Prepare patch list by rank list
  !-----------------------------------------------------------------------
  subroutine pbod_prepare_patchList(rect_list)
    use overBlockCoordinates
    type(t_obRect),dimension(:),intent(IN) :: rect_list
    type(t_obRect) :: rectgh
    integer,dimension(MX:MZ) :: extension = (/1, 1, 1/)*(N_GHOSTCELL + N_MARGIN )
    integer :: n

    allocate(Patch_List(lbound(rect_list, 1):ubound(rect_list, 1)))
    allocate(Rank_list(lbound(rect_list, 1):ubound(rect_list, 1)))
    call pbod_make_ranklist
    do n = lbound(Patch_List, 1), ubound(Patch_List, 1)
       call ob_rectExtend(rect_list(n), extension, rectgh)
       call pb_gather(rectgh, Patch_List(n)%patch, B_MLIST, Rank_list(n))
    end do
  end subroutine pbod_prepare_patchList
  !-----------------------------------------------------------------------
  ! Restore patch list
  !-----------------------------------------------------------------------
  subroutine pbod_restore_patchList
    integer :: n
    do n = lbound(Patch_List, 1), ubound(Patch_List, 1)
       call pb_scatter(Patch_List(n)%patch, Rank_list(n))
    end do
    deallocate(Patch_List)
    deallocate(Rank_list)
  end subroutine pbod_restore_patchList
  !-----------------------------------------------------------------------
  ! execute Ohmic dissipation for all patch blocks
  !-----------------------------------------------------------------------
  subroutine pbod_for_patchList
    use mpilib
    integer :: n
    myrank = get_myrank()
    do n = lbound(Rank_list, 1), ubound(Rank_list, 1)
       if (myrank == Rank_list(n)) then
          call pbod_subcycle_for_patch(Patch_List(n)%patch, Eta_list(n))
       end if
    end do
  end subroutine pbod_for_patchList
  !-----------------------------------------------------------------------
  ! Execute ohmic dissipation (serial process)
  !-----------------------------------------------------------------------
  subroutine pbod_subcycle_for_patch(patch, etasph)
    use grid
    type(t_upatch),pointer :: patch ! (INOUT)
    type(t_eta_sphere),intent(IN) :: etasph ! parameters of spherical eta
    integer :: gridlevel
    real(kind=DBL_KIND) :: hd_dt, dt, localtime, etamax
    logical :: bool_last
    real(kind=DBL_KIND),dimension(:,:,:,:,:),allocatable :: f
    real(kind=DBL_KIND),dimension(:,:,:,:),allocatable :: bn
    real(kind=DBL_KIND),dimension(:,:,:),allocatable :: eta

    Isgh = lbound(patch%u, 1)   ! array size
    Iegh = ubound(patch%u, 1)
    Jsgh = lbound(patch%u, 2)
    Jegh = ubound(patch%u, 2)
    Ksgh = lbound(patch%u, 3)
    Kegh = ubound(patch%u, 3)
    Is = Isgh + N_GHOSTCELL
    Ie = Iegh - N_GHOSTCELL
    Js = Jsgh + N_GHOSTCELL
    Je = Jegh - N_GHOSTCELL
    Ks = Ksgh + N_GHOSTCELL
    Ke = Kegh - N_GHOSTCELL
    Ms = lbound(patch%u, 4)
    Me = ubound(patch%u, 4)
    ! allocate work array
    allocate(f(AIJKMNGH), bn(AIJKMGH), eta(AIJKGH))

    gridlevel = patch%rect%level
    hd_dt = Dtime(gridlevel)

    bool_last = .FALSE.
    localtime = 0.d0

#ifdef TEST_PROBLEM
    call pbod_get_eta(eta, etamax, patch)
#else !TEST_PROBLEM
    call pbod_get_eta(eta, etamax, patch, etasph)
#endif !TEST_PROBLEM
    dt = CFL_OD * minval(patch%cellwidth(:))**2 / etamax
    print *, 'nsubcycle OD =', int(hd_dt/dt)+1
    do
       if (localtime + dt >= hd_dt) then
          bool_last = .TRUE.
          dt = hd_dt - localtime
       endif
       call pbod_t_integ(dt, patch, f, bn, eta)
       localtime = localtime + dt
#ifdef VERBOSE
       print *, localtime
#endif !VERBOSE
       if (bool_last) exit
    end do

    deallocate(f, bn, eta)
  end subroutine pbod_subcycle_for_patch
  !-----------------------------------------------------------------------
  ! add diffusion for given grid (serial process)
  !-----------------------------------------------------------------------
  subroutine pbod_t_integ(dt, patch, f, bn, eta)
    use grid
    real(kind=DBL_KIND),intent(IN) :: dt
    type(t_upatch),pointer :: patch ! (INOUT)
    real(kind=DBL_KIND),dimension(AIJKMNGH),intent(INOUT) :: f
    real(kind=DBL_KIND),dimension(AIJKMGH),intent(INOUT) :: bn
    real(kind=DBL_KIND),dimension(AIJKGH),intent(IN) :: eta
    real(kind=DBL_KIND) :: dxi, dyi, dzi
    integer :: i, j, k, m
    dxi = 1.d0/patch%cellwidth(MX)
    dyi = 1.d0/patch%cellwidth(MY)
    dzi = 1.d0/patch%cellwidth(MZ)
    ! predictor step
    call pbod_get_flux(f, patch%u, eta, patch%cellwidth)
    do m = Ms, Me
       do k = Ks, Ke
          do j = Js, Je
             do i = Is, Ie
                bn(i,j,k,m) = patch%u(i,j,k,m) -  &
                     ( (f(i,j,k,m,MX)-f(i-1,j,k,m,MX))*dxi &
                     + (f(i,j,k,m,MY)-f(i,j-1,k,m,MY))*dyi &
                     + (f(i,j,k,m,MZ)-f(i,j,k-1,m,MZ))*dzi ) *dt * 0.5d0
             end do
          end do
       end do
    end do
    
    ! corrector step
    call pbod_get_flux(f, bn, eta, patch%cellwidth)
    do m = Ms, Me
       do k = Ks, Ke
          do j = Js, Je
             do i = Is, Ie
                patch%u(i,j,k,m) = patch%u(i,j,k,m) -  &
                     ( (f(i,j,k,m,MX)-f(i-1,j,k,m,MX))*dxi &
                     + (f(i,j,k,m,MY)-f(i,j-1,k,m,MY))*dyi &
                     + (f(i,j,k,m,MZ)-f(i,j,k-1,m,MZ))*dzi ) *dt
             end do
          end do
       end do
    end do
#ifdef VERBOSE
    print *, 'bmax', maxval(abs(patch%u)), patch%u((Is+Ie)/2, (Js+Je)/2, (Ks+Ke)/2, MZ)
#endif !VERBOSE
  end subroutine pbod_t_integ
  !-----------------------------------------------------------------------
  ! get flux of Ohmic dissipation (serial process)
  !-----------------------------------------------------------------------
  subroutine pbod_get_flux(f, b, eta, cellwidth)
    real(kind=DBL_KIND),dimension(AIJKMNGH),intent(OUT) :: f
    real(kind=DBL_KIND),dimension(AIJKMGH),intent(IN) :: b
    real(kind=DBL_KIND),dimension(AIJKGH),intent(IN) :: eta
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: cellwidth
    real(kind=DBL_KIND) :: dxi, dyi, dzi, etabar
    integer :: i, j, k, n, nbx, nby, nbz
    dxi = 1.d0/cellwidth(MX)
    dyi = 1.d0/cellwidth(MY)
    dzi = 1.d0/cellwidth(MZ)
    nbx = Ms
    nby = Ms + 1
    nbz = Ms + 2
    ! flux at i+1/2
    n = MX
    do k = ks, ke
       do j = js, je
          do i = is-1, ie
             etabar = (eta(i+1,j,k) + eta(i,j,k))*0.5d0
             f(i,j,k,nbx,n) = 0.d0
             f(i,j,k,nby,n) = ((-b(i+1,j,k,nby)+b(i,j,k,nby))*dxi + (b(i,j+1,k,nbx)-b(i,j-1,k,nbx)+b(i+1,j+1,k,nbx)-b(i+1,j-1,k,nbx))*0.25d0*dyi)*etabar
             f(i,j,k,nbz,n) = ((-b(i+1,j,k,nbz)+b(i,j,k,nbz))*dxi + (b(i,j,k+1,nbx)-b(i,j,k-1,nbx)+b(i+1,j,k+1,nbx)-b(i+1,j,k-1,nbx))*0.25d0*dzi)*etabar
          end do
       end do
    end do
    ! flux at j+1/2
    n = MY
    do k = ks, ke
       do j = js-1, je
          do i = is, ie
             etabar = (eta(i,j+1,k) + eta(i,j,k))*0.5d0
             f(i,j,k,nbx,n) = ((-b(i,j+1,k,nbx)+b(i,j,k,nbx))*dyi + (b(i+1,j,k,nby)-b(i-1,j,k,nby)+b(i+1,j+1,k,nby)-b(i-1,j+1,k,nby))*dxi*0.25d0)*etabar
             f(i,j,k,nby,n) = 0.d0
             f(i,j,k,nbz,n) = ((-b(i,j+1,k,nbz)+b(i,j,k,nbz))*dyi + (b(i,j,k+1,nby)-b(i,j,k-1,nby)+b(i,j+1,k+1,nby)-b(i,j+1,k-1,nby))*dzi*0.25d0)*etabar
          end do
       end do
    end do
    ! flux at k+1/2
    n = MZ
    do k = ks-1, ke
       do j = js, je
          do i = is, ie
             etabar = (eta(i,j,k+1) + eta(i,j,k))*0.5d0
             f(i,j,k,nbx,n) = ((-b(i,j,k+1,nbx)+b(i,j,k,nbx))*dzi + (b(i+1,j,k,nbz)-b(i-1,j,k,nbz)+b(i+1,j,k+1,nbz)-b(i-1,j,k+1,nbz))*dxi*0.25d0)*etabar
             f(i,j,k,nby,n) = ((-b(i,j,k+1,nby)+b(i,j,k,nby))*dzi + (b(i,j+1,k,nbz)-b(i,j-1,k,nbz)+b(i,j+1,k+1,nbz)-b(i,j-1,k+1,nbz))*dyi*0.25d0)*etabar
             f(i,j,k,nbz,n) = 0.d0
          end do
       end do
    end do
  end subroutine pbod_get_flux
  !-----------------------------------------------------------------------
  ! make Rank_list in order to distribute over the ranks
  !-----------------------------------------------------------------------
  subroutine pbod_make_ranklist
    use mpilib
    integer :: nrank, nrect, skip, n, ns, ne
    nrank = NPE
    nrect = size(Rank_list)
    skip = nrank/nrect
    ns = lbound(Rank_list,1)
    ne = ubound(Rank_list,1)
    if (skip > 0) then ! nrank >= nrect
       do n = ns, ne
          Rank_list(n) = (n+1-ns)*skip-1
       end do
    else ! nrank < nrect
       do n = ns, ne
          Rank_list(n) = mod(n-ns, nrank)
       end do
    end if
#ifdef VERBOSE
    if (get_myrank() == PRIMARY_RANK) print *, 'Rank_list OD =', Rank_list
#endif !VERBOSE
  end subroutine pbod_make_ranklist

#ifdef TEST_PROBLEM
  !-----------------------------------------------------------------------
  ! define eta (serial process)
  !-----------------------------------------------------------------------
  subroutine pbod_get_eta(eta, etamax, patch)
    real(kind=DBL_KIND),dimension(AIJKGH),intent(OUT) :: eta
    real(kind=DBL_KIND),intent(OUT) :: etamax
    type(t_upatch),pointer :: patch ! (IN)
    integer :: i, j, k
    real(kind=DBL_KIND) :: r, posz
    etamax = 1.d0
    posz = (patch%z(ubound(patch%z,1)) + patch%z(lbound(patch%z,1)))/2.d0

    do k = Ksgh, Kegh
       do j = Jsgh, Jegh
          do i = Isgh, Iegh
             r = sqrt(patch%x(i)**2 +patch%y(j)**2 + (patch%z(k)-posz)**2)
             eta(i,j,k) = exp(-r**2)*etamax
          end do
       end do
    end do
  end subroutine pbod_get_eta
#else  !TEST_PROBLEM
  !-----------------------------------------------------------------------
  ! define eta for spherical region (serial process)
  !-----------------------------------------------------------------------
  ! Approximate formula for sigmoid (tanh x + 1)/2
  ! SIGMOID(-1) = 0, SIGMOID(0) = 0.5, SIGMOID(1) = 1
#define SIGMOID(X) min(max(0.5d0*((X)+1.d0), 0.d0), 1.d0)
! #define SIGMOID(X) ((tanh(X)+1.d0)*0.5d0)
  subroutine pbod_get_eta(eta, etamax, patch, etasph)
    real(kind=DBL_KIND),dimension(AIJKGH),intent(OUT) :: eta
    real(kind=DBL_KIND),intent(OUT) :: etamax
    type(t_upatch),pointer :: patch ! (IN)
    type(t_eta_sphere),intent(IN) :: etasph ! parameters of spherical eta
    real(kind=DBL_KIND) :: radius, r
    integer :: i, j, k
!!$    ! vk = sqrt(GM/r)
!!$    vk = sqrt(MP_Gconst * MP_ms1 / sp_SinkRadius)
!!$    eta = vk * sp_SinkRadius

    etamax = etasph%etamax
    do k = Ksgh, Kegh
       do j = Jsgh, Jegh
          do i = Isgh, Iegh
             r = sqrt((patch%x(i)-etasph%center(MX))**2 + (patch%y(j)-etasph%center(MY))**2 + (patch%z(k)-etasph%center(MZ))**2)
             eta(i,j,k) = (1.d0 - SIGMOID((r -etasph%radius)/etasph%dtaper))*etasph%etamax
             ! Taper region: [etasph%radius-etasph%dtaper, etasph%radius+etasph%dtaper]
!!$             if (r <= etasph%radius) then
!!$                eta(i,j,k) = 1.d0
!!$             else
!!$                eta(i,j,k) = 0.d0
!!$             endif
!!$             eta(i,j,k) = exp(-r**2)*etamax
          end do
       end do
    end do

  end subroutine pbod_get_eta
#endif !TEST_PROBLEM
end module pb_ohmicdissp
