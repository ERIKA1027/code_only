module overBlockCoordinates
  use grid
  implicit none
  public
  integer,parameter :: OB_NCOORDS = 2*((2)-(0)+1)
  integer,parameter :: OB_COORDS_MIN = 0
  integer,parameter :: OB_COORDS_MAX = OB_COORDS_MIN + OB_NCOORDS - 1
  type t_obPoint
     integer(kind=8),dimension(0:2) :: p
     integer :: level
  end type t_obPoint
  type t_obRect
     integer(kind=8),dimension(0:2) :: left
     integer(kind=8),dimension(0:2) :: right
     integer :: level
  end type t_obRect
  type t_obPointPhys
     real(kind=8),dimension(0:2) :: p
  end type t_obPointPhys
  type t_obRectPhys
     real(kind=8),dimension(0:2) :: left
     real(kind=8),dimension(0:2) :: right
  end type t_obRectPhys
  integer,parameter :: OB_RECT_VALID = 0
  integer,parameter :: OB_RECT_INVALID = -1
  integer,parameter :: OB_RECT_UNDEF = -2
  integer,parameter :: OB_RECT_INNER = 1
  integer,parameter :: OB_RECT_OUTER = 2
  integer,parameter :: OB_RECT_EQUIV = 3
  integer,parameter :: OB_RECT_CROSS = 4
  integer,parameter :: OB_RECT_DETOUCH = 5
  logical,save :: Initialized = .false.
  type(t_obRectPhys),save :: RectComputationBoxPhys
  type(t_obRect),save,dimension(Lmin:Lmax) :: RectComputationBox
  integer(kind=8),parameter :: LLONG = 1
  integer,parameter :: SHORT = 1
  real(kind=8),parameter :: DBL = 1.d0
  integer,save :: MPI_OB_POINT, MPI_OB_POINTPHYS, MPI_OB_RECT, MPI_OB_RECTPHYS
  private :: Initialized, LLONG, SHORT, DBL, RectComputationBoxPhys, RectComputationBox
contains
  subroutine ob_undefPoint(point)
    type(t_obPoint),intent(INOUT) :: point
    point%p(:) = HUGE(LLONG)
    point%level = HUGE(SHORT)
  end subroutine ob_undefPoint
  subroutine ob_undefPointPhys(point)
    type(t_obPointPhys),intent(INOUT) :: point
    point%p(:) = HUGE(DBL)
  end subroutine ob_undefPointPhys
  subroutine ob_undefRect(rect)
    type(t_obRect),intent(INOUT) :: rect
    type(t_obPoint) :: point
    call ob_undefPoint(point)
    call ob_assignPointToRect(point, rect, 'L')
    call ob_assignPointToRect(point, rect, 'R')
  end subroutine ob_undefRect
  subroutine ob_undefRectPhys(rectPhys)
    type(t_obRectPhys),intent(INOUT) :: rectPhys
    type(t_obPointPhys) :: pointPhys
    call ob_undefPointPhys(pointPhys)
    call ob_assignPointPhysToRectPhys(pointPhys, rectPhys, 'L')
    call ob_assignPointPhysToRectPhys(pointPhys, rectPhys, 'R')
  end subroutine ob_undefRectPhys
  function ob_definedPoint(point) result(bool)
    type(t_obPoint),intent(IN) :: point
    logical :: bool
    if ( point%p(0) >= HUGE(LLONG) .and. &
         point%p(1) >= HUGE(LLONG) .and. &
         point%p(2) >= HUGE(LLONG) .and. &
         point%level >= HUGE(SHORT) ) then
       bool = .false.
    else
       bool = .true.
    end if
  end function ob_definedPoint
  function ob_definedPointPhys(pointPhys) result(bool)
    type(t_obPointPhys),intent(IN) :: pointPhys
    logical :: bool
    if ( pointPhys%p(0) >= HUGE(DBL) .and. &
         pointPhys%p(1) >= HUGE(DBL) .and. &
         pointPhys%p(2) >= HUGE(DBL) ) then
       bool = .false.
    else
       bool = .true.
    end if
  end function ob_definedPointPhys
  function ob_definedRect(rect) result(bool)
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: pointL, pointR
    logical :: bool
    call ob_extractPointFromRect(pointL, rect, 'L')
    call ob_extractPointFromRect(pointR, rect, 'R')
    if (ob_definedPoint(pointL) .and. ob_definedPoint(pointR)) then
       bool = .true.
    else
       bool = .false.
    end if
  end function ob_definedRect
  function ob_definedRectPhys(rectPhys) result(bool)
    type(t_obRectPhys),intent(IN) :: rectPhys
    type(t_obPointPhys) :: pointPhysL, pointPhysR
    logical :: bool
    call ob_extractPointPhysFromRectPhys(pointPhysL, rectPhys, 'L')
    call ob_extractPointPhysFromRectPhys(pointPhysR, rectPhys, 'R')
    if (ob_definedPointPhys(pointPhysL) .and. ob_definedPointPhys(pointPhysR)) then
       bool = .true.
    else
       bool = .false.
    end if
  end function ob_definedRectPhys
  subroutine ob_assignCoordToPoint(point, coords, level )
    type(t_obPoint),intent(OUT) :: point
    integer(kind=8),dimension(0:2),intent(IN) :: coords
    integer,intent(IN) :: level
    point%p(:) = coords
    point%level = level
  end subroutine ob_assignCoordToPoint
  subroutine ob_assignPointToRect(point, rect, leftOrRight)
    type(t_obPoint),intent(IN) :: point
    type(t_obRect),intent(OUT) :: rect
    character(len=1),intent(IN) :: leftOrRight
    if ( leftOrRight .eq. 'L' ) then
       rect%left(:) = point%p(:)
    elseif ( leftOrRight .eq. 'R' ) then
       rect%right(:) = point%p(:)
    endif
    rect%level = point%level
  end subroutine ob_assignPointToRect
  subroutine ob_assignPointPairToRect(pointL, pointR, rect)
    type(t_obPoint),intent(IN) :: pointL, pointR
    type(t_obRect),intent(OUT) :: rect
    rect%left(:) = pointL%p(:)
    rect%right(:) = pointR%p(:)
    rect%level = pointL%level
    if (pointL%level /= pointR%level) print *, '*** Error: inconsistent level', pointL%level, pointR%level
  end subroutine ob_assignPointPairToRect
  subroutine ob_assignCoordToRect( coords, rect, level )
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX),intent(IN) :: coords
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: level
    type(t_obPoint) :: point
    integer(kind=8),dimension(0:2) :: cd
    cd = coords(0:2)
    call ob_assignCoordToPoint( point, cd, level)
    call ob_assignPointToRect( point, rect, 'L' )
    cd = coords(2 +1:OB_COORDS_MAX)
    call ob_assignCoordToPoint( point, cd, level )
    call ob_assignPointToRect( point, rect, 'R' )
  end subroutine ob_assignCoordToRect
  subroutine ob_assignPointPhysToRectPhys(point, rect, leftOrRight)
    type(t_obPointPhys),intent(IN) :: point
    type(t_obRectPhys),intent(OUT) :: rect
    character(len=1),intent(IN) :: leftOrRight
    if ( leftOrRight .eq. 'L' ) then
       rect%left(:) = point%p(:)
    elseif ( leftOrRight .eq. 'R' ) then
       rect%right(:) = point%p(:)
    endif
  end subroutine ob_assignPointPhysToRectPhys
  subroutine ob_assignPointPhysPairToRectPhys(pointL, pointR, rect)
    type(t_obPointPhys),intent(IN) :: pointL, pointR
    type(t_obRectPhys),intent(OUT) :: rect
    rect%left(:) = pointL%p(:)
    rect%right(:) = pointR%p(:)
  end subroutine ob_assignPointPhysPairToRectPhys
  subroutine ob_assignCoordPhysToPointPhys(coords, point)
    real(kind=8),dimension(0:2),intent(IN) :: coords
    type(t_obPointPhys),intent(OUT) :: point
    point%p(0:2) = coords(0:2)
  end subroutine ob_assignCoordPhysToPointPhys
  subroutine ob_assignCoordPhysToRectPhys(coords, rect)
    real(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX),intent(IN) :: coords
    type(t_obRectPhys),intent(OUT) :: rect
    type(t_obPointPhys) :: point
    real(kind=8),dimension(0:2) :: cd
    cd = coords(0:2)
    call ob_assignCoordPhysToPointPhys( cd, point )
    call ob_assignPointPhysToRectPhys( point, rect, 'L' )
    cd = coords(2 +1:OB_COORDS_MAX)
    call ob_assignCoordPhysToPointPhys( cd, point )
    call ob_assignPointPhysToRectPhys( point, rect, 'R' )
  end subroutine ob_assignCoordPhysToRectPhys
  subroutine ob_extractPointFromRect(point, rect, leftOrRight)
    type(t_obPoint),intent(OUT) :: point
    type(t_obRect),intent(IN) :: rect
    character(len=1),intent(IN) :: leftOrRight
    if ( leftOrRight .eq. 'L' ) then
       point%p(:) = rect%left(:)
    elseif ( leftOrRight .eq. 'R' ) then
       point%p(:) = rect%right(:)
    endif
    point%level = rect%level
  end subroutine ob_extractPointFromRect
  subroutine ob_extractPointPairFromRect(pointL, pointR, rect)
    type(t_obPoint),intent(OUT) :: pointL, pointR
    type(t_obRect),intent(IN) :: rect
    pointL%p(:) = rect%left(:)
    pointR%p(:) = rect%right(:)
    pointL%level = rect%level
    pointR%level = rect%level
  end subroutine ob_extractPointPairFromRect
  subroutine ob_extractPointPhysFromRectPhys(point, rect, leftOrRight)
    type(t_obPointPhys),intent(OUT) :: point
    type(t_obRectPhys),intent(IN) :: rect
    character(len=1),intent(IN) :: leftOrRight
    if ( leftOrRight .eq. 'L' ) then
       point%p(:) = rect%left(:)
    elseif ( leftOrRight .eq. 'R' ) then
       point%p(:) = rect%right(:)
    endif
  end subroutine ob_extractPointPhysFromRectPhys
  subroutine ob_extractPointPhysPairFromRectPhys(pointL, pointR, rect)
    type(t_obPointPhys),intent(OUT) :: pointL, pointR
    type(t_obRectPhys),intent(IN) :: rect
    pointL%p(:) = rect%left(:)
    pointR%p(:) = rect%right(:)
  end subroutine ob_extractPointPhysPairFromRectPhys
  subroutine ob_extractCoordFromRect(coords, rect)
    integer(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX),intent(OUT) :: coords
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: point
    call ob_extractPointFromRect(point, rect, 'L')
    coords(0:2) = point%p(:)
    call ob_extractPointFromRect(point, rect, 'R')
    coords(2 +1:OB_COORDS_MAX) = point%p(:)
  end subroutine ob_extractCoordFromRect
  subroutine ob_extractLevelFromRect(level, rect)
    integer,intent(OUT) :: level
    type(t_obRect),intent(IN) :: rect
    level = rect%level
  end subroutine ob_extractLevelFromRect
  subroutine ob_extractLevelFromPoint(level, point)
    integer,intent(OUT) :: level
    type(t_obPoint),intent(IN) :: point
    level = point%level
  end subroutine ob_extractLevelFromPoint
  subroutine ob_extractCoordPhysFromRectPhys(coordPhys, rectPhys)
    real(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX),intent(OUT) :: coordPhys
    type(t_obRectPhys),intent(IN) :: rectPhys
    type(t_obPointPhys) :: point
    call ob_extractPointPhysFromRectPhys(point, rectPhys, 'L')
    coordPhys(0:2) = point%p(:)
    call ob_extractPointPhysFromRectPhys(point, rectPhys, 'R')
    coordPhys(2 +1:OB_COORDS_MAX) = point%p(:)
  end subroutine ob_extractCoordPhysFromRectPhys
  subroutine ob_PointPhys2PointOb( pointPhys, level, pointOb )
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN) :: level
    type(t_obPoint),intent(OUT) :: pointOb
    type(t_obPointPhys) :: compbox
    if ( .not. Initialized ) call ob_init
    call ob_extractPointPhysFromRectPhys(compbox, RectComputationBoxPhys, 'L')
    pointOb%p(:) = (pointPhys%p(:)-compbox%p(:))/CellWidth(:,level)
    pointOb%level = level
  end subroutine ob_PointPhys2PointOb
  subroutine ob_RectPhys2RectOb( rectPhys, level, rectOb)
    type(t_obRectPhys),intent(IN) :: rectPhys
    integer,intent(IN) :: level
    type(t_obRect),intent(OUT) :: rectOb
    type(t_obPointPhys) :: pointPhys
    type(t_obPoint) :: pointOb
    call ob_extractPointPhysFromRectPhys( pointPhys, rectPhys, 'L' )
    call ob_PointPhys2PointOb( pointPhys, level, pointOb )
    call ob_assignPointToRect( pointOb, rectOb, 'L' )
    call ob_extractPointPhysFromRectPhys( pointPhys, rectPhys, 'R' )
    call ob_PointPhys2PointOb( pointPhys, level , pointOb )
    call ob_assignPointToRect( pointOb, rectOb, 'R' )
  end subroutine ob_RectPhys2RectOb
  subroutine ob_PointOb2PointPhys( pointOb, pointPhys )
    type(t_obPoint),intent(IN) :: pointOb
    type(t_obPointPhys),intent(OUT) :: pointPhys
    integer :: level
    type(t_obPointPhys) :: compbox
    if ( .not. Initialized ) call ob_init
    call ob_extractPointPhysFromRectPhys(compbox, RectComputationBoxPhys, 'L')
    level = pointOb%level
    pointPhys%p(:) = (pointOb%p(:) + 0.5d0) * CellWidth(:,level) + compbox%p(:)
  end subroutine ob_PointOb2PointPhys
  subroutine ob_RectOb2RectPhys( rectOb, rectPhys )
    type(t_obRect),intent(IN) :: rectOb
    type(t_obRectPhys),intent(OUT) :: rectPhys
    type(t_obPoint) :: point
    type(t_obPointPhys) :: pointPhys
    call ob_extractPointFromRect(point, rectOb, 'L')
    call ob_PointOb2PointPhys( point, pointPhys )
    call ob_assignPointPhysToRectPhys( pointPhys, rectPhys, 'L')
    call ob_extractPointFromRect(point, rectOb, 'R')
    call ob_PointOb2PointPhys( point, pointPhys )
    call ob_assignPointPhysToRectPhys( pointPhys, rectPhys, 'R')
  end subroutine ob_RectOb2RectPhys
  subroutine ob_point2PointByLevel(pointIn, pointOut, levelOut, LR)
    type(t_obPoint),intent(IN) :: pointIn
    type(t_obPoint),intent(OUT) :: pointOut
    integer,intent(IN) :: levelOut
    character(len=1),intent(IN),optional :: LR
    integer :: levelIn, offset
    integer(KIND=8),parameter :: two =2
    call ob_extractLevelFromPoint(levelIn, pointIn)
    if ( levelOut >= levelIn ) then
       if ( present(LR) ) then
          if (LR == 'L') offset = 0
          if (LR == 'R') offset = 1
       else
          offset = 0
       endif
       if ( offset == 0 ) then
          pointOut%p(:) = pointIn%p(:) * two** (levelOut - levelIn)
       else
          pointOut%p(:) = ( pointIn%p(:) + 1 ) * two** (levelOut - levelIn) -1
       endif
    else
       pointOut%p(:) = pointIn%p(:) / two** (levelIn - levelOut)
    end if
    pointOut%level = levelOut
  end subroutine ob_point2PointByLevel
  subroutine ob_rect2RectByLevel(rectIn, rectOut, levelOut)
    type(t_obRect),intent(IN) :: rectIn
    type(t_obRect),intent(OUT) :: rectOut
    integer,intent(IN) :: levelOut
    type(t_obPoint) :: pointIn, pointOut
    call ob_extractPointFromRect(pointIn, rectIn, 'L')
    call ob_point2PointByLevel(pointIn, pointOut, levelOut, 'L')
    call ob_assignPointToRect(pointOut, rectOut, 'L')
    call ob_extractPointFromRect(pointIn, rectIn, 'R')
    call ob_point2PointByLevel(pointIn, pointOut, levelOut, 'R')
    call ob_assignPointToRect(pointOut, rectOut, 'R')
  end subroutine ob_rect2RectByLevel
  function ob_rectsComp( rectA, rectB ) result( comp )
    type(t_obRect),intent(IN) :: rectA, rectB
    integer :: comp
    if ( rectA%level /= rectB%level ) then
       print *, '*** error in ob_rectsComp: levels are not consistent'
       return
    endif
    if ( ob_testRect(rectA) == OB_RECT_INVALID ) then
       comp = OB_RECT_INVALID
       return
    end if
    if ( ob_testRect(rectB) == OB_RECT_INVALID ) then
       comp = OB_RECT_INVALID
       return
    end if
    if (.not. (ob_definedRect(rectA) .and. ob_definedRect(rectB))) then
       comp = OB_RECT_UNDEF
       return
    endif
    if ( rectA%left(0) == rectB%left(0) .and. &
         rectA%left(1) == rectB%left(1) .and. &
         rectA%left(2) == rectB%left(2) .and. &
         rectA%right(0) == rectB%right(0) .and. &
         rectA%right(1) == rectB%right(1) .and. &
         rectA%right(2) == rectB%right(2) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_EQUIV
       return
    end if
    if ( rectA%left(0) <= rectB%left(0) .and. &
         rectA%left(1) <= rectB%left(1) .and. &
         rectA%left(2) <= rectB%left(2) .and. &
         rectA%right(0) >= rectB%right(0) .and. &
         rectA%right(1) >= rectB%right(1) .and. &
         rectA%right(2) >= rectB%right(2) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_OUTER
       return
    end if
    if ( rectA%left(0) >= rectB%left(0) .and. &
         rectA%left(1) >= rectB%left(1) .and. &
         rectA%left(2) >= rectB%left(2) .and. &
         rectA%right(0) <= rectB%right(0) .and. &
         rectA%right(1) <= rectB%right(1) .and. &
         rectA%right(2) <= rectB%right(2) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_INNER
       return
    end if
    if ( max(rectA%left(0),rectB%left(0)) <= min(rectA%right(0),rectB%right(0) ) .and. &
         max(rectA%left(1),rectB%left(1)) <= min(rectA%right(1),rectB%right(1) ) .and. &
         max(rectA%left(2),rectB%left(2)) <= min(rectA%right(2),rectB%right(2) ) ) then
         comp = OB_RECT_CROSS
         return
      end if
      comp = OB_RECT_DETOUCH
  end function ob_rectsComp
  function ob_rectPhysComp( rectA, rectB ) result ( comp )
    type(t_obRectPhys),intent(IN) :: rectA, rectB
    integer :: comp
    if ( ob_testRectPhys(rectA) == OB_RECT_INVALID ) then
       comp = OB_RECT_INVALID
       return
    end if
    if ( ob_testRectPhys(rectB) == OB_RECT_INVALID ) then
       comp = OB_RECT_INVALID
       return
    end if
    if (.not. (ob_definedRectPhys(rectA) .and. ob_definedRectPhys(rectB))) then
       comp = OB_RECT_UNDEF
       return
    endif
    if ( rectA%left(0) == rectB%left(0) .and. &
         rectA%left(1) == rectB%left(1) .and. &
         rectA%left(2) == rectB%left(2) .and. &
         rectA%right(0) == rectB%right(0) .and. &
         rectA%right(1) == rectB%right(1) .and. &
         rectA%right(2) == rectB%right(2) ) then
       comp = OB_RECT_EQUIV
       return
    end if
    if ( rectA%left(0) <= rectB%left(0) .and. &
         rectA%left(1) <= rectB%left(1) .and. &
         rectA%left(2) <= rectB%left(2) .and. &
         rectA%right(0) >= rectB%right(0) .and. &
         rectA%right(1) >= rectB%right(1) .and. &
         rectA%right(2) >= rectB%right(2) ) then
       comp = OB_RECT_OUTER
       return
    end if
    if ( rectA%left(0) >= rectB%left(0) .and. &
         rectA%left(1) >= rectB%left(1) .and. &
         rectA%left(2) >= rectB%left(2) .and. &
         rectA%right(0) <= rectB%right(0) .and. &
         rectA%right(1) <= rectB%right(1) .and. &
         rectA%right(2) <= rectB%right(2) ) then
       comp = OB_RECT_INNER
       return
    end if
    if ( max(rectA%left(0),rectB%left(0)) <= min(rectA%right(0),rectB%right(0) ) .and. &
         max(rectA%left(1),rectB%left(1)) <= min(rectA%right(1),rectB%right(1) ) .and. &
         max(rectA%left(2),rectB%left(2)) <= min(rectA%right(2),rectB%right(2) ) ) then
         comp = OB_RECT_CROSS
         return
      end if
    comp = OB_RECT_DETOUCH
  end function ob_rectPhysComp
  subroutine ob_rectAnd( rectA, rectB, rect )
    type(t_obRect),intent(IN) :: rectA, rectB
    type(t_obRect),intent(OUT) :: rect
    integer :: m
    if ( rectA%level /= rectB%level ) then
       print *, '*** ob_rectAnd : level is inconsistent between rectA and rectB', rectA%level, rectB%level
       stop
    endif
    do m = 0, 2
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
    rect%level = rectA%level
  end subroutine ob_rectAnd
  subroutine ob_rectOr( rectA, rectB, rect )
    type(t_obRect),intent(IN) :: rectA, rectB
    type(t_obRect),intent(OUT) :: rect
    integer :: m
    if ( rectA%level /= rectB%level ) then
       print *, '*** ob_rectAnd : level is inconsistent between rectA and rectB', rectA%level, rectB%level
       stop
    endif
    do m = 0, 2
       rect%left(m) = min(rectA%left(m), rectB%left(m) )
       rect%right(m) = max(rectA%right(m), rectB%right(m) )
    enddo
    rect%level = rectA%level
  end subroutine ob_rectOr
  subroutine ob_rectPhysAnd( rectA, rectB, rect )
    type(t_obRectPhys),intent(IN) :: rectA, rectB
    type(t_obRectPhys),intent(OUT) :: rect
    integer :: m
    do m = 0, 2
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
  end subroutine ob_rectPhysAnd
  subroutine ob_rectPhysOr( rectA, rectB, rect )
    type(t_obRectPhys),intent(IN) :: rectA, rectB
    type(t_obRectPhys),intent(OUT) :: rect
    type(t_obRectPhys) :: rectTmp
    integer :: m
    call ob_rectPhysAnd( rectA, rectB, rectTmp)
    do m = 0, 2
       if ( rectTmp%left(m) > rectTmp%right(m) ) then
          print *, '*** ob_rectPhysOr: rectA and rectB are detouched'
       endif
    enddo
    do m = 0, 2
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
  end subroutine ob_rectPhysOr
  subroutine ob_rectPhysTrim ( rectIn, h, rectOut )
    type(t_obRectPhys),intent(IN) :: rectIn
    type(t_obRectPhys),intent(OUT) :: rectOut
    real(kind=8),dimension(0:2),intent(IN) :: h
    integer :: n
    do n = 0, 2
       rectOut%left(n) = rectIn%left(n) + h(n)
       rectOut%right(n) = rectIn%right(n) - h(n)
    enddo
  end subroutine ob_rectPhysTrim
  subroutine ob_rectTrim ( rectIn, h, rectOut )
    type(t_obRect),intent(IN) :: rectIn
    type(t_obRect),intent(OUT) :: rectOut
    integer,dimension(0:2),intent(IN) :: h
    integer :: n
    do n = 0, 2
       rectOut%left(n) = rectIn%left(n) + h(n)
       rectOut%right(n) = rectIn%right(n) - h(n)
    enddo
    rectOut%level = rectIn%level
  end subroutine ob_rectTrim
  subroutine ob_rectPhysExtend(rectIn, h, rectOut)
    type(t_obRectPhys),intent(IN) :: rectIn
    type(t_obRectPhys),intent(OUT) :: rectOut
    real(kind=8),dimension(0:2),intent(IN) :: h
    real(kind=8),dimension(0:2) :: minush
    minush = -h
    call ob_rectPhysTrim(rectIn, minush, rectOut)
  end subroutine ob_rectPhysExtend
  subroutine ob_rectExtend(rectIn, h, rectOut)
    type(t_obRect),intent(IN) :: rectIn
    type(t_obRect),intent(OUT) :: rectOut
    integer,dimension(0:2),intent(IN) :: h
    integer,dimension(0:2) :: minush
    minush = -h
    call ob_rectTrim(rectIn, minush, rectOut)
  end subroutine ob_rectExtend
  function ob_testRect( rect ) result( code )
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: pointL, pointR
    integer :: code
    call ob_extractPointFromRect(pointL, rect, 'L')
    call ob_extractPointFromRect(pointR, rect, 'R')
    if ( .not. ob_definedRect( rect ) ) then
       code = OB_RECT_UNDEF
    else if ( &
         ( pointL%p(0) <= pointR%p(0) ) .and. &
         ( pointL%p(1) <= pointR%p(1) ) .and. &
         ( pointL%p(2) <= pointR%p(2) ) ) then
       code = OB_RECT_VALID
    else
       code = OB_RECT_INVALID
    endif
  end function ob_testRect
  function ob_testRectPhys( rect ) result( code )
    type(t_obRectPhys),intent(IN) :: rect
    type(t_obPointPhys) :: pointL, pointR
    integer :: code
    call ob_extractPointPhysFromRectPhys(pointL, rect, 'L')
    call ob_extractPointPhysFromRectPhys(pointR, rect, 'R')
    if ( .not. ob_definedRectPhys( rect ) ) then
       code = OB_RECT_UNDEF
    else if ( &
         ( pointL%p(0) <= pointR%p(0) ) .and. &
         ( pointL%p(1) <= pointR%p(1) ) .and. &
         ( pointL%p(2) <= pointR%p(2) ) ) then
       code = OB_RECT_VALID
    else
       code = OB_RECT_INVALID
    endif
  end function ob_testRectPhys
  subroutine ob_computationBoxOfRectPhys( rectPhys )
    type(t_obRectPhys),intent(OUT) :: rectPhys
    if ( .not. Initialized ) call ob_init
    rectPhys = RectComputationBoxPhys
  end subroutine ob_ComputationBoxOfRectPhys
  subroutine ob_computationBoxOfRect( rect, level )
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: level
    if ( .not. Initialized ) call ob_init
    rect = RectComputationBox(level)
  end subroutine ob_computationBoxOfRect
  subroutine ob_computationBoxOfCoordPhys( coordPhys )
    real(kind=8),dimension(OB_COORDS_MIN:OB_COORDS_MAX),intent(OUT) :: coordPhys
    type(t_obRectPhys) :: rectPhys
    call ob_computationBoxOfRectPhys( rectPhys )
    call ob_extractCoordPhysFromRectPhys(coordPhys, rectPhys)
  end subroutine ob_computationBoxOfCoordPhys
  subroutine ob_getIjkgridFromPointPhys(ijkgrid, level, pointPhys)
    integer,dimension(0:2),intent(OUT) :: ijkgrid
    integer,intent(IN) :: level
    type(t_obPointPhys),intent(IN) :: pointPhys
    type(t_obPoint) :: point
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getIjkgridFromPoint(ijkgrid, point)
  end subroutine ob_getIjkgridFromPointPhys
  subroutine ob_getIjkgridFromPoint(ijkgrid, point)
    integer,dimension(0:2),intent(OUT) :: ijkgrid
    type(t_obPoint),intent(IN) :: point
    ijkgrid(:) = (/ &
         point%p(0) / (Imax-Imin+1), &
         point%p(1) / (Jmax-Jmin+1), &
         point%p(2) / (Kmax-Kmin+1) &
         /)
  end subroutine ob_getIjkgridFromPoint
  subroutine ob_getIjkgridFromCoordPhys(ijkgrid, level, coordPhys)
    integer,dimension(0:2),intent(OUT) :: ijkgrid
    integer,intent(IN) :: level
    real(kind=8),dimension(0:2),intent(IN) :: coordPhys
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coordPhys, pointPhys)
    call ob_getIjkgridFromPointPhys(ijkgrid, level, pointPhys)
  end subroutine ob_getIjkgridFromCoordPhys
  subroutine ob_getBlockFromPoint(point, gid, rank, i,j,k)
    type(t_obPoint),intent(IN) :: point
    integer,intent(OUT) :: gid, rank, i,j,k
    integer :: ig,jg,kg, level
    ig = point%p(0) / (Imax-Imin+1)
    jg = point%p(1) / (Jmax-Jmin+1)
    kg = point%p(2) / (Kmax-Kmin+1)
    level = point%level
    i = mod( point%p(0), Imax-Imin+1 )
    j = mod( point%p(1), Jmax-Jmin+1 )
    k = mod( point%p(2), Kmax-Kmin+1 )
    call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
  end subroutine ob_getBlockFromPoint
  subroutine ob_getBlockFromPointPhys(pointPhys, level, gid, rank, i,j,k)
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN) :: level
    integer,intent(OUT) :: gid, rank, i,j,k
    type(t_obPoint) :: point
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getBlockFromPoint(point, gid, rank, i,j,k)
  end subroutine ob_getBlockFromPointPhys
  subroutine ob_getBlockFromCoordPhys(coords, level, gid, rank, i,j,k)
    real(kind=8),dimension(0:2),intent(IN) :: coords
    integer,intent(IN) :: level
    integer,intent(OUT) :: gid, rank, i,j,k
    type(t_obPoint) :: point
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coords, pointPhys)
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getBlockFromPoint(point, gid, rank, i,j,k)
  end subroutine ob_getBlockFromCoordPhys
  subroutine ob_getPointFromIjkgrid(point, ijkgrid, level, i, j, k)
    type(t_obPoint), intent(OUT) :: point
    integer,dimension(0:2),intent(IN) :: ijkgrid
    integer,intent(IN) :: level, i, j, k
    point%p = (/ &
         ijkgrid(0) * (Imax-Imin+1) + i, &
         ijkgrid(1) * (Jmax-Jmin+1) + j, &
         ijkgrid(2) * (Kmax-Kmin+1) + k /)
    point%level = level
  end subroutine ob_getPointFromIjkgrid
  subroutine ob_getPointFromBlockSerial(point, gid, i, j, k, retcode)
    type(t_obPoint), intent(OUT) :: point
    integer,intent(IN) :: gid, i, j, k
    integer,intent(OUT) :: retcode
    integer :: level
    integer(kind=8),dimension(0:2) :: coords
    level = get_level(gid)
    if ( level == Undefi ) then
       retcode = Undefi
       return
    end if
    retcode = 1
    point%p = (/ &
         Igrid(gid) * (Imax-Imin+1) + i, &
         Jgrid(gid) * (Jmax-Jmin+1) + j, &
         Kgrid(gid) * (Kmax-Kmin+1) + k /)
    point%level = level
  end subroutine ob_getPointFromBlockSerial
  subroutine ob_getPointFromBlock(point, gid, rank, i, j, k, retcode)
    use mpilib
    type(t_obPoint), intent(OUT) :: point
    integer,intent(IN) :: gid, rank, i, j, k
    integer,intent(OUT) :: retcode
    if ( rank == get_myrank() ) then
       call ob_getPointFromBlockSerial(point, gid, i, j, k, retcode)
    endif
    call mpi_bcast( retcode, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, ierr )
    if (retcode == Undefi) return
    call mpi_bcast( point, 1, MPI_OB_POINT, rank, MPI_COMM_WORLD, ierr )
  end subroutine ob_getPointFromBlock
  subroutine ob_getRectFromIjkgrid(rect, ijkgrid, level)
    type(t_obRect),intent(OUT) :: rect
    integer,dimension(0:2),intent(IN) :: ijkgrid
    integer,intent(IN) :: level
    type(t_obPoint) :: pointL, pointR
    call ob_getPointFromIjkgrid(pointL, ijkgrid, level, Imin, Jmin, Kmin)
    call ob_getPointFromIjkgrid(pointR, ijkgrid, level, Imax, Jmax, Kmax)
    call ob_assignPointPairToRect(pointL, pointR, rect)
  end subroutine ob_getRectFromIjkgrid
  subroutine ob_getRectFromGid(rect, gid)
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: gid
    integer,dimension(0:2) :: ijkgrid
    ijkgrid = (/Igrid(gid), Jgrid(gid), Kgrid(gid)/)
    call ob_getRectFromIjkgrid(rect, ijkgrid, get_level(gid))
  end subroutine ob_getRectFromGid
  function ob_PointWithinRect(point, rect) result(bool)
    type(t_obPoint),intent(IN) :: point
    type(t_obRect),intent(IN) :: rect
    logical :: bool
    type(t_obPoint) :: pL, pR
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    if ( pL%p(0) <= point%p(0) .and. point%p(0) <= pR%p(0) .and. &
         pL%p(1) <= point%p(1) .and. point%p(1) <= pR%p(1) .and. &
         pL%p(2) <= point%p(2) .and. point%p(2) <= pR%p(2) ) then
       bool = .true.
    else
       bool = .false.
    endif
  end function ob_PointWithinRect
  function ob_PointPhysWithinRectPhys(point, rect) result(bool)
    type(t_obPointPhys),intent(IN) :: point
    type(t_obRectPhys),intent(IN) :: rect
    logical :: bool
    type(t_obPointPhys) :: pL, pR
    call ob_extractPointPhysFromRectPhys(pL, rect, 'L')
    call ob_extractPointPhysFromRectPhys(pR, rect, 'R')
    if ( pL%p(0) <= point%p(0) .and. point%p(0) <= pR%p(0) .and. &
         pL%p(1) <= point%p(1) .and. point%p(1) <= pR%p(1) .and. &
         pL%p(2) <= point%p(2) .and. point%p(2) <= pR%p(2) ) then
       bool = .true.
    else
       bool = .false.
    endif
  end function ob_PointPhysWithinRectPhys
  subroutine ob_init
    use mpilib
    real(kind=8),dimension(:),pointer :: xp, yp, zp
    integer :: gid, rank, level
    real(kind=8) :: xyzmin(0:2), xyzmax(0:2), coords(OB_COORDS_MIN:OB_COORDS_MAX)
    real(kind=8),save :: xmin, ymin, zmin, xmax, ymax, zmax
    type(t_obPoint) :: point
    integer(KIND=8),parameter :: two =2
    integer(KIND=8) :: power
    level = 0
    call get_gid_from_ijkgrid(0,0,0,level,gid,rank)
    myrank = get_myrank()
    if ( myrank == rank ) then
       xp => get_Xp(gid)
       yp => get_Yp(gid)
       zp => get_Zp(gid)
       xmin = xp(Imin) - CellWidth(0,level)/2.d0
       ymin = yp(Jmin) - CellWidth(1,level)/2.d0
       zmin = zp(Kmin) - CellWidth(2,level)/2.d0
       xyzmin(:) = (/xmin, ymin, zmin/)
    endif
    call mpi_bcast(xyzmin, size(xyzmin), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmin = xyzmin(0)
    ymin = xyzmin(1)
    Zmin = xyzmin(2)
    gid = GidBase( ubound(GidBase,1),ubound(GidBase,2),ubound(GidBase,3) )
    rank = RankBase( ubound(GidBase,1),ubound(GidBase,2),ubound(GidBase,3) )
    if ( myrank == rank ) then
       xp => get_Xp(gid)
       yp => get_Yp(gid)
       zp => get_Zp(gid)
       xmax = xp(Imax) + CellWidth(0,level)/2.d0
       ymax = yp(Jmax) + CellWidth(1,level)/2.d0
       zmax = zp(Kmax) + CellWidth(2,level)/2.d0
       xyzmax(:) = (/xmax, ymax, zmax/)
    endif
    call mpi_bcast(xyzmax, size(xyzmax), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    xmax = xyzmax(0)
    ymax = xyzmax(1)
    zmax = xyzmax(2)
    coords = (/xmin, ymin, zmin, xmax, ymax, zmax/)
    call ob_assignCoordPhysToRectPhys(coords, RectComputationBoxPhys)
    do level = Lmin, Lmax
       point%level = level
       point%p(:) = 0
       call ob_assignPointToRect(point, RectComputationBox(level), 'L')
       power = two**(level-Lmin)
       point%p(0) = (8)*(8)*power-1
       point%p(1) = (8)*(8)*power-1
       point%p(2) = (8)*(8)*power-1
       call ob_assignPointToRect(point, RectComputationBox(level), 'R')
    end do
    call ob_Make_MPI_struct_point
    call ob_Make_MPI_struct_rect
    call ob_Make_MPI_struct_pointPhys
    call ob_Make_MPI_struct_rectPhys
    Initialized = .true.
    if ( get_myrank() == 0 ) print *, 'initialize overBlockCoordinates'
  contains
    subroutine ob_Make_MPI_struct_point
      integer,dimension(2) :: blocklen, disp, type
      integer :: sz
      type(t_obPoint) :: point
      sz = size(point%p)
      blocklen = (/sz, 1/)
      disp = (/0, blocklen(1)*(8)/)
      type = (/MPI_LLONG_INTEGER, MPI_INTEGER /)
      sz = size(blocklen)
      call mpi_type_struct(sz, blocklen, disp, type, MPI_OB_POINT, ierr)
      call mpi_type_commit( MPI_OB_POINT, ierr )
    end subroutine ob_Make_MPI_struct_point
    subroutine ob_Make_MPI_struct_rect
      integer,dimension(3) :: blocklen, disp, type
      type(t_obRect) :: rect
      integer :: szl, szr, sz
      szl = size(rect%left)
      szr = size(rect%right)
      blocklen = (/szl, szr, 1/)
      disp = (/0, blocklen(1)*(8), blocklen(1)*(8)+blocklen(2)*(8) /)
      type = (/MPI_LLONG_INTEGER, MPI_LLONG_INTEGER, MPI_INTEGER /)
      sz = size(blocklen)
      call mpi_type_struct(sz, blocklen, disp, type, MPI_OB_RECT, ierr)
      call mpi_type_commit( MPI_OB_RECT, ierr )
    end subroutine ob_Make_MPI_struct_rect
    subroutine ob_Make_MPI_struct_pointPhys
      type(t_obPointPhys) :: point
      integer :: sz
      sz = size(point%p)
      call mpi_type_contiguous(sz, MPI_DOUBLE_PRECISION, MPI_OB_POINTPHYS, ierr)
      call mpi_type_commit( MPI_OB_POINTPHYS, ierr )
    end subroutine ob_Make_MPI_struct_pointPhys
    subroutine ob_Make_MPI_struct_rectPhys
      type(t_obRectPhys) :: rect
      integer :: sz
      sz = size(rect%left)+size(rect%right)
      call mpi_type_contiguous(sz, MPI_DOUBLE_PRECISION, MPI_OB_RECTPHYS, ierr)
      call mpi_type_commit( MPI_OB_RECTPHYS, ierr )
    end subroutine ob_Make_MPI_struct_rectPhys
  end subroutine ob_init
end module overBlockCoordinates
