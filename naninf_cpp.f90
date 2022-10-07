function isnan(a) result(bool)
  implicit none
  logical :: bool
  real(kind=8),intent(IN) :: a
  if (a /= a ) then
     bool = .true.
  else
     bool = .false.
  end if
end function isnan
function isinf(a) result(bool)
  implicit none
  logical :: bool
  real(kind=8),intent(IN) :: a
  logical :: isnan
  if ( isnan(a) ) return
  if ( (a*0) /= 0 ) then
     bool = .true.
  else
     bool = .false.
  end if
end function isinf
function isNotFinite(a) result(bool)
  implicit none
  logical :: bool
  real(kind=8),intent(IN) :: a
  if ( (a*0) /= 0 .or. a /= a) then 
     bool = .true.
  else
     bool = .false.
  end if
end function isNotFinite
function isFinite(a) result(bool)
  implicit none
  logical :: bool
  real(kind=8),intent(IN) :: a
  if ( (a*0) == 0 .and. a == a) then 
     bool = .true.
  else
     bool = .false.
  end if
end function isFinite
function anyIsNotFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=8),intent(IN),dimension(sz) :: a
  bool = any( (a*0) /= 0 .or. a /= a) 
end function anyIsNotFinite
function anyIsFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=8),intent(IN),dimension(sz) :: a
  bool = any( (a*0) == 0 .and. a == a) 
end function anyIsFinite
function allIsNotFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=8),intent(IN),dimension(sz) :: a
  bool = all((a*0) /= 0 .or. a /= a) 
end function allIsNotFinite
function allIsFinite(a, sz) result(bool)
  implicit none
  logical :: bool
  integer,intent(IN) :: sz
  real(kind=8),intent(IN),dimension(sz) :: a
    bool = all((a*0) == 0 .and. a == a) 
end function allIsFinite
