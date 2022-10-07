#include "config.h"
!-------------------------------------------------------------------------
! Module for FMG cubic interpolation
!-------------------------------------------------------------------------
module fmg_interpol_cubic
  use fmg_data
  implicit none
  private
  ! buffer for tri-cubic interpolation
  type t_buf
     real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: bufx  => null(), bufxy => null()
  end type t_buf
  type(t_buf),save,dimension(:),pointer :: FLV => null() ! for vector
  type(t_buf),save,dimension(:),pointer :: FLS => null() ! for scalar
  logical,save :: Initialized = .FALSE.
  public :: fmg_interp_cubic, fmg_interp_cubic_finalize, fmg_interp_cubic_check
contains
  ! ----------------------------------------------------------------
  ! Initialize
  ! ----------------------------------------------------------------
  subroutine fmg_interp_cubic_init
    integer :: fmglev
    integer :: ifs, jfs, kfs, ife, jfe, kfe
    integer :: ics, jcs, kcs, ice, jce, kce
    if (Initialized) return
    Initialized = .TRUE.
    allocate( FLV(FMG_LevelMin: FMG_LevelMax-1) )
    allocate( FLS(FMG_LevelMin: FMG_LevelMax-1) )
    do fmglev = lbound(FLV, 1), ubound(FLV, 1)
       call fmg_get_gridsize  (fmglev,   ifs, jfs, kfs, ife, jfe, kfe) ! w/o ghost cell
       call fmg_get_gridsizeGh(fmglev+1, ics, jcs, kcs, ice, jce, kce) !with ghost cell
       allocate(FLV(fmglev)%bufx (ifs:ife,jcs:jce,kcs:kce,Mmin:Mmax))
       allocate(FLV(fmglev)%bufxy(ifs:ife,jfs:jfe,kcs:kce,Mmin:Mmax))
       allocate(FLS(fmglev)%bufx (ifs:ife,jcs:jce,kcs:kce,Mmin:Mmin))
       allocate(FLS(fmglev)%bufxy(ifs:ife,jfs:jfe,kcs:kce,Mmin:Mmin))
    end do
  end subroutine fmg_interp_cubic_init
  ! ----------------------------------------------------------------
  ! Finalize
  ! ----------------------------------------------------------------
  subroutine fmg_interp_cubic_finalize
    integer :: fmglev
    if (.not. Initialized) return
    Initialized = .FALSE.
    do fmglev = lbound(FLV,1), ubound(FLV,1)
       deallocate( &
            FLV(fmglev)%bufx, &
            FLV(fmglev)%bufxy, &
            FLS(fmglev)%bufx, &
            FLS(fmglev)%bufxy )
       nullify( &
            FLV(fmglev)%bufx, &
            FLV(fmglev)%bufxy, &
            FLS(fmglev)%bufx, &
            FLS(fmglev)%bufxy )
    end do
    deallocate(FLV, FLS)
    nullify(FLV, FLS)
  end subroutine fmg_interp_cubic_finalize
  ! ----------------------------------------------------------------
  ! check grid size
  ! ----------------------------------------------------------------
  function fmg_interp_cubic_check(fmglev) result(bool_cubic)
    use io_util
    integer,intent(IN) :: fmglev
    logical :: bool_cubic
    integer,parameter :: MINBLOCKSIZE = 4
    integer :: ics,jcs,kcs,ice,jce,kce
    call fmg_get_gridsizeGh(fmglev+1, ics,jcs,kcs,ice,jce,kce)
    if ( &
         ice-ics+1 < MINBLOCKSIZE .or. &
         jce-jcs+1 < MINBLOCKSIZE .or. &
         kce-kcs+1 < MINBLOCKSIZE ) then
       bool_cubic = .FALSE.
       call print_msg('*** fmg_interp_cubic: Block size is too small to perform cubic interpolation.')
       call print_msg('*** fmg_interp_cubic: Instead, linear interpolation is used.')
    else
       bool_cubic = .TRUE.
    end if
  end function fmg_interp_cubic_check
  ! ----------------------------------------------------------------
  ! Tricubic interpolation
  ! ----------------------------------------------------------------
  subroutine fmg_interp_cubic(fmglev, iuf, iuc)
    use interpolation
    integer,intent(IN) :: fmglev, iuf, iuc
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: uf, uc
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: bufx, bufxy
    type(t_buf),dimension(:),pointer :: flp
    integer :: la, gid, lf, lc, ifs, jfs, kfs, ife, jfe, kfe

    call fmg_interp_cubic_init
    ! switch buffer
    if (fmg_isVector(iuc)) then
       flp => FLV
    else
       flp => FLS
    end if
    bufx  => flp(fmglev)%bufx
    bufxy => flp(fmglev)%bufxy

    lf = fmglev
    lc = fmglev+1
    call fmg_get_gridsize(lf, ifs, jfs, kfs, ife, jfe, kfe)
    do la = AMR_LevelMin, AMR_LevelMax
       do gid = fmg_get_gidmin(la), fmg_get_gidmax(la)
          if ( fmg_skip_grid(gid, la, lf) ) cycle
          uf => fmg_get_arrp(la,lf,gid,iuf)
          uc => fmg_get_arrp(la,lc,gid,iuc)
          call interp_tricubic(uf, uc, bufx, bufxy, &
               ifs, jfs, kfs, ife, jfe, kfe, &
               ifs, jfs, kfs)
       end do
    end do
  end subroutine fmg_interp_cubic
end module fmg_interpol_cubic
