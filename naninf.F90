!-------------------------------------------------------------------------
! Detection of NaN, Inf.
! This file should be compiled with a option of
!    -Kieee  ..... for PGI fortran
!    -mieee-fp ... for intel fortran
!-------------------------------------------------------------------------
#include "config.h"
!-------------------------------------------------------------------------
! detect NaN.  It returns .true. if a is NaN.
!-------------------------------------------------------------------------
function isnan(a) result(bool)
  implicit none
  logical :: bool
  real(kind=DBL_KIND),intent(IN) :: a
  if (a /= a ) then
     bool = .true.
  else
     bool = .false.
  end if
end function isnan
!-------------------------------------------------------------------------
! detect Inf.  It returns .true. if a is Inf.
!-------------------------------------------------------------------------
function isinf(a) result(bool)
  implicit none
  logical :: bool
  real(kind=DBL_KIND),intent(IN) :: a
  logical :: isnan
  if ( isnan(a) ) return
  if ( (a*0) /= 0 ) then
     bool = .true.
  else
     bool = .false.
  end if
end function isinf
!-------------------------------------------------------------------------
! return true if a is not finite value.
!-------------------------------------------------------------------------
function isNotFinite(a) result(bool)
  implicit none
  logical :: bool
  real(kind=DBL_KIND),intent(IN) :: a
!  if ( (a*0) /= 0 ) then
  if ( (a*0) /= 0 .or. a /= a) then ! KS MODIFIED
     bool = .true.
  else
     bool = .false.
  end if
end function isNotFinite
!-------------------------------------------------------------------------
! return true if a is finite value.
!-------------------------------------------------------------------------
function isFinite(a) result(bool)
  implicit none
  logical :: bool
  real(kind=DBL_KIND),intent(IN) :: a
  !  if ( (a*0) == 0 ) then
  if ( (a*0) == 0 .and. a == a) then ! KS MODIFIED
     bool = .true.
  else
     bool = .false.
  end if
end function isFinite
!-------------------------------------------------------------------------
! return true if any elements in array a are not finite
!-------------------------------------------------------------------------
function anyIsNotFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=DBL_KIND),intent(IN),dimension(sz) :: a
  ! bool = any( (a*0) /= 0)
  bool = any( (a*0) /= 0 .or. a /= a) ! KS MODIFIED
end function anyIsNotFinite
!-------------------------------------------------------------------------
! return true if any elements in array a are finite
!-------------------------------------------------------------------------
function anyIsFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=DBL_KIND),intent(IN),dimension(sz) :: a
  ! bool = any( (a*0) == 0) 
  bool = any( (a*0) == 0 .and. a == a)  ! KS MODIFIED
end function anyIsFinite
!-------------------------------------------------------------------------
! return true if all the elements in array a are not finite
!-------------------------------------------------------------------------
function allIsNotFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=DBL_KIND),intent(IN),dimension(sz) :: a
  ! bool = all((a*0) /= 0)
  bool = all((a*0) /= 0  .or. a /= a) ! KS MODIFIED
end function allIsNotFinite
!-------------------------------------------------------------------------
! return true if all the elements in array a are finite
!-------------------------------------------------------------------------
function allIsFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=DBL_KIND),intent(IN),dimension(sz) :: a
  ! bool = all((a*0) == 0)
    bool = all((a*0) == 0 .and. a == a) ! KS MODIFIED
end function allIsFinite
