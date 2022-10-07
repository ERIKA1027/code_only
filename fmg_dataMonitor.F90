#include "config.h"
#define VERBOSE
!-------------------------------------------------------------------------
! Monitoring data for skipping procedure
!-------------------------------------------------------------------------
module fmg_dataMonitor
  implicit none
  private
  real(kind=DBL_KIND),save,dimension(:,:,:),allocatable :: Sample
  logical,save :: Initialized = .FALSE.
  public : fmg_dataMonitor_init, fmg_dataMonitor_finalize, fmg_dataMonitor_sample
contains
  !-------------------------------------------------------------------------
  ! initialize data monitor
  !-------------------------------------------------------------------------
  subroutine fmg_dataMonitor_init
    if (Initialized) return
    Initialized = .TRUE.
    if (.not. allocated(Sample)) then
       allocate(Sample(AMR_LevelMin:AMR_LevelMax, FMG_LevelMin:FMG_LevelMax, ICODEMIN:ICODEMAX))
       Sample = HUGE(Sample)    ! set some value for initialized
    end if
  end subroutine fmg_dataMonitor_init
  !-------------------------------------------------------------------------
  ! finalize data monitor
  !-------------------------------------------------------------------------
  subroutine fmg_dataMonitor_finalize
    Initialized = .FALSE.
    deallocate(Sample)
  end subroutine fmg_dataMonitor_finalize
  !-------------------------------------------------------------------------
  ! get Gid
  !-------------------------------------------------------------------------
  function getGid(amrlev, fmglev, icode) result(gid)
    integer,intent(IN) :: amrlev, fmglev, icode
    integer :: gid
    gid = lbound(GridLevel(amrlev,fmglev)%Block,1)
  end function getGid
  !-------------------------------------------------------------------------
  ! extract data
  !-------------------------------------------------------------------------
  function extractData(amrlev, fmglev, icode) result(data)
    integer,intent(IN) :: amrlev, fmglev, icode
    real(kind=DBL_KIND),dimension(:,:,:,:),pointer :: arr
    integer :: gid, imin, jmin, kmin
    gid = getGid(amrlev, fmglev, icode)
    arr => fmg_get_arrp(amrlev, fmglev, gid, icode)
    imin = fmg_get_imin(fmglev)
    jmin = fmg_get_jmin(fmglev)
    kmin = fmg_get_kmin(fmglev)
    data = arr(imin,jmin,kmin,lbound(arr,4))
  end function extractData
  !-------------------------------------------------------------------------
  ! sampling data
  !-------------------------------------------------------------------------
  subroutine fmg_dataMonitor_sample(amrlev, fmglev, icode)
    integer,intent(IN) :: amrlev, fmglev, icode
    Sample(amrlev, fmglev, icode) = extractData(amrlev, fmglev, icode)
  end subroutine fmg_dataMonitor_sample
  !-------------------------------------------------------------------------
  ! skip by data unchanged
  !-------------------------------------------------------------------------
  function fmg_dataMonitor_skip(amrlev, fmglev, icode) result(bool_skip)
    integer,intent(IN) :: amrlev, fmglev, icode
    ! data has not been updated
    bool_skip = (extractData(amrlev, fmglev, icode) == Sample(amrlev, fmglev, icode))
#ifdef VERBOSE
    print *, '*** fmg_dataMonitor_skip', amrlev, fmglev, icode
#endif !VERBOSE
  end function fmg_dataMonitor_skip
  !-------------------------------------------------------------------------
  ! skip by comparing with fmg parent data
  !-------------------------------------------------------------------------
  function fmg_dataMonitor_skip_by_fmgParent(amrlev, fmglev, icode) result(bool_skip)
    integer,intent(IN) :: amrlev, fmglev, icode
    integer :: fmgParent
    fmgParent = fmglev-1
    bool_skip = (extractData(amrlev, fmglev, icode) == Sample(amrlev, fmglev, icode))

  end function fmg_dataMonitor_skip


end module fmg_dataMonitor
