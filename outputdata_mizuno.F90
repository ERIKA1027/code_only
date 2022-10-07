#include "config.h"
!-------------------------------------------------------------------------
!
! Module for controlling output data
!
!-------------------------------------------------------------------------
module outputdata
  implicit none
  private
  public :: output_data
contains
  ! -----------------------------------------------------------------
  ! output procedures 
  ! -----------------------------------------------------------------
  subroutine output_data
    use grid, only : Lmin, LevelMax
    use writeSnap
    use mpilib
    integer :: scope_level
    if (.not. bool_output() ) return

    call writeSnap_whole
    do scope_level=Lmin, LevelMax
       call output_centralbox(scope_level) ! output (Lambda_max/2^scope_level)^3 region
    end do
  end subroutine output_data
  ! -----------------------------------------------------------------
  ! Output central box.
  ! The box size is lambda_max / 2^scope_level
  ! File name is cb stage . resolution_level . d
  ! -----------------------------------------------------------------
  subroutine output_centralbox(scope_level)
    use string, only : CHARLEN, num2char
    use grid, only : LevelMax
    use parameter
    use uniformgrid, only : uniformgrid_write
    integer,intent(IN) :: scope_level
    character(len=2),parameter :: prefixDefault = 'cb'
    real(kind=DBL_KIND) :: halfwidth, boxsize, kzmax, lambda_max
    character(len=CHARLEN) :: prefix

    kzmax = 0.285d0                ! wave number of most unstable perturbation
    lambda_max = 2.d0 * PI / kzmax ! wave length of most unstable perturbation

    boxsize = lambda_max/2.d0 ** scope_level
    halfwidth = boxsize * 0.5d0
!!$    prefix = prefixDefault // trim(num2char(scope_level)) // '.'
    prefix = prefixDefault

    if (scope_level > LevelMax) return !skip

    call uniformgrid_write(-halfwidth,-halfwidth,-halfwidth,halfwidth,halfwidth,halfwidth, scope_level, interpolate=.true.,prefix=prefix)
    
  end subroutine output_centralbox
  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming.
  ! -----------------------------------------------------------------
  function bool_output() result(bool)
    use grid
    use eos
    use analysis, only : RhoMax
    logical :: bool
    real(kind=DBL_KIND),parameter :: logrho_skip = 0.2d0 ! interval for output in unit of log(rho)
    real(kind=DBL_KIND),save :: logrhoio
    real(kind=DBL_KIND),save :: rhomax_prev = Huge(rhomax_prev)
    logical,save :: bool_output_initialized = .false.

!!$    ! -------------
!!$    ! for debug
!!$    ! -------------
!!$    bool = .false.
!!$    if (level_sync() == Lmin .and. mod(Step(Lmin),100) == 0 ) then
!!$       bool = .true.
!!$    endif
!!$    return

    ! ------------------------------------------------------
    ! for gravitational collapsing cloud
    ! 中心密度の対数が等間隔になるようにデータを出力する。
    ! 間隔は logrho_skip で設定する。
    ! ------------------------------------------------------
    bool = .false.
    if ( level_sync() /= Lmin ) return
    if ( .not. bool_output_initialized ) then
       logrhoio = (int(log10(RhoMax)/logrho_skip)+1)*logrho_skip
       bool_output_initialized = .true.
    endif

    if ( log10(RhoMax) > logrhoio .and. log10(rhomax_prev) < logrhoio ) then
       logrhoio = logrhoio + logrho_skip
       bool = .TRUE.
       return
    endif
    rhomax_prev = RhoMax

  end function bool_output
end module outputdata
