#include "config.h"
#include "overBlockCoordinates.h"
!-------------------------------------------------------------------------
! Module for overBlockCoordinates
!
! This module defnes two types of coordinates, over-block coordinates
! and physical coordinates.
!
! Definition of Over-Block coordinates:
!    The coordinates is equivalent to cell number of a given grid level.
!    A cell at the lower-left coorner is specified by (0,0,0).
!    Therefore over-block coordinates depend on the grid level.
!
!-------------------------------------------------------------------------
module overBlockCoordinates
  use grid
  implicit none
  public
  ! number of coordinates for rectangle
  integer,parameter :: OB_NCOORDS = 2*((MZ)-(MX)+1)
  integer,parameter :: OB_COORDS_MIN = MX
  integer,parameter :: OB_COORDS_MAX = OB_COORDS_MIN + OB_NCOORDS - 1
  ! a point in the overblock coordinates
  type t_obPoint
     integer(kind=LLONG_KIND),dimension(MX:MZ) :: p
     integer :: level
  end type t_obPoint
  ! a rectangle in the overblock coordinates
  type t_obRect
     integer(kind=LLONG_KIND),dimension(MX:MZ) :: left
     integer(kind=LLONG_KIND),dimension(MX:MZ) :: right
     integer :: level
  end type t_obRect
  ! a point in the physical coordinates
  type t_obPointPhys
     real(kind=DBL_KIND),dimension(MX:MZ) :: p
  end type t_obPointPhys
  ! a rectangle in the physical coordinates
  type t_obRectPhys
     real(kind=DBL_KIND),dimension(MX:MZ) :: left
     real(kind=DBL_KIND),dimension(MX:MZ) :: right
  end type t_obRectPhys
  ! constant parameter for test
  integer,parameter :: OB_RECT_VALID   =  0
  integer,parameter :: OB_RECT_INVALID = -1
  integer,parameter :: OB_RECT_UNDEF   = -2
  integer,parameter :: OB_RECT_INNER   =  1
  integer,parameter :: OB_RECT_OUTER   =  2
  integer,parameter :: OB_RECT_EQUIV   =  3
  integer,parameter :: OB_RECT_CROSS   =  4
  integer,parameter :: OB_RECT_DETOUCH =  5
  ! private
  logical,save :: Initialized = .false.
  type(t_obRectPhys),save :: RectComputationBoxPhys
  type(t_obRect),save,dimension(Lmin:Lmax) :: RectComputationBox
  integer(kind=LLONG_KIND),parameter :: LLONG = 1
  integer,parameter :: SHORT = 1
  real(kind=DBL_KIND),parameter :: DBL = 1.d0
  ! MPI type
  integer,save :: MPI_OB_POINT, MPI_OB_POINTPHYS, MPI_OB_RECT, MPI_OB_RECTPHYS
  private :: Initialized, LLONG, SHORT, DBL, RectComputationBoxPhys, RectComputationBox
contains
  !-------------------------------------------------------------------------
  ! undefine point
  !-------------------------------------------------------------------------
  subroutine ob_undefPoint(point)
    type(t_obPoint),intent(INOUT) :: point
    point%p(:) = HUGE(LLONG)
    point%level = HUGE(SHORT)
  end subroutine ob_undefPoint
  !-------------------------------------------------------------------------
  ! undefine pointphys
  !-------------------------------------------------------------------------
  subroutine ob_undefPointPhys(point)
    type(t_obPointPhys),intent(INOUT) :: point
    point%p(:) = HUGE(DBL)
  end subroutine ob_undefPointPhys
  !-------------------------------------------------------------------------
  ! undefine rect
  !-------------------------------------------------------------------------
  subroutine ob_undefRect(rect)
    type(t_obRect),intent(INOUT) :: rect
    type(t_obPoint) :: point
    call ob_undefPoint(point)
    call ob_assignPointToRect(point, rect, 'L')
    call ob_assignPointToRect(point, rect, 'R')
  end subroutine ob_undefRect
  !-------------------------------------------------------------------------
  ! undefine rectPhys
  !-------------------------------------------------------------------------
  subroutine ob_undefRectPhys(rectPhys)
    type(t_obRectPhys),intent(INOUT) :: rectPhys
    type(t_obPointPhys) :: pointPhys
    call ob_undefPointPhys(pointPhys)
    call ob_assignPointPhysToRectPhys(pointPhys, rectPhys, 'L')
    call ob_assignPointPhysToRectPhys(pointPhys, rectPhys, 'R')
  end subroutine ob_undefRectPhys
  !-------------------------------------------------------------------------
  ! return true if point is defined.
  !-------------------------------------------------------------------------
  function ob_definedPoint(point) result(bool)
    type(t_obPoint),intent(IN) :: point
    logical :: bool
    if ( point%p(MX) >= HUGE(LLONG) .and. &
         point%p(MY) >= HUGE(LLONG) .and. &
         point%p(MZ) >= HUGE(LLONG) .and. &
         point%level >= HUGE(SHORT) ) then
       bool = .false.
    else
       bool = .true.
    end if
  end function ob_definedPoint
  !-------------------------------------------------------------------------
  ! return true if point is defined.
  !-------------------------------------------------------------------------
  function ob_definedPointPhys(pointPhys) result(bool)
    type(t_obPointPhys),intent(IN) :: pointPhys
    logical :: bool
    if ( pointPhys%p(MX) >= HUGE(DBL) .and. &
         pointPhys%p(MY) >= HUGE(DBL) .and. &
         pointPhys%p(MZ) >= HUGE(DBL) ) then
       bool = .false.
    else
       bool = .true.
    end if
  end function ob_definedPointPhys
  !-------------------------------------------------------------------------
  ! return true if rect is defined.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! return true if rectPhys is defined.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! Assign coordinates to a point
  ! INPUT: iob, job, kob, level
  ! OUTPUT: point
  !-------------------------------------------------------------------------
  subroutine ob_assignCoordToPoint(point, coords, level )
    type(t_obPoint),intent(OUT) :: point
    integer(kind=LLONG_KIND),dimension(MX:MZ),intent(IN) :: coords
    integer,intent(IN) :: level
    point%p(:) = coords
    point%level = level
  end subroutine ob_assignCoordToPoint
  !-------------------------------------------------------------------------
  ! Assgin a point to a rectangle.
  ! INPUT:
  !    point = a point in overblock coordinates.
  !    leftOrRight = 'L' or 'R'.
  ! OUTPUT:
  !    rect = rectangle in overblock coordinates.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! Assgin point pair to a rectangle
  !-------------------------------------------------------------------------
  subroutine ob_assignPointPairToRect(pointL, pointR, rect)
    type(t_obPoint),intent(IN) :: pointL, pointR
    type(t_obRect),intent(OUT) :: rect
    rect%left(:)  = pointL%p(:)
    rect%right(:) = pointR%p(:)
    rect%level = pointL%level
    if (pointL%level /= pointR%level) print *, '*** Error: inconsistent level', pointL%level, pointR%level
  end subroutine ob_assignPointPairToRect
  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  subroutine ob_assignCoordToRect( coords, rect, level )
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ),intent(IN) :: coords
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: level
    type(t_obPoint) :: point
    integer(kind=LLONG_KIND),dimension(MX:MZ) :: cd
    cd = coords(MX:MZ)
    call ob_assignCoordToPoint( point, cd, level)
    call ob_assignPointToRect( point, rect, 'L' )
    cd = coords(MZ+1:OB_COORDS_MAX)
    call ob_assignCoordToPoint( point, cd, level )
    call ob_assignPointToRect( point, rect, 'R' )
  end subroutine ob_assignCoordToRect
  !-------------------------------------------------------------------------
  ! Assgin a point to a rectangle in physical coordinates.
  ! INPUT:
  !    point = a point in physical coordinates.
  !    leftOrRight = 'L' or 'R'.
  ! OUTPUT:
  !    rect = rectangle in physical coordinates.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! Assign two pointPhys pair to rectphys
  !-------------------------------------------------------------------------
  subroutine ob_assignPointPhysPairToRectPhys(pointL, pointR, rect)
    type(t_obPointPhys),intent(IN) :: pointL, pointR
    type(t_obRectPhys),intent(OUT) :: rect
    rect%left(:)  = pointL%p(:)
    rect%right(:) = pointR%p(:)
  end subroutine ob_assignPointPhysPairToRectPhys
  !-------------------------------------------------------------------------
  ! Assign physical coordinates to physical point
  ! INPUT:
  !   coords = (/ x, y, z /), a vector of physical coordinates
  ! OUTPUT:
  !   point = a point in the physical coordinates
  !-------------------------------------------------------------------------
  subroutine ob_assignCoordPhysToPointPhys(coords, point)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: coords
    type(t_obPointPhys),intent(OUT) :: point
    point%p(MX:MZ) = coords(MX:MZ)
  end subroutine ob_assignCoordPhysToPointPhys
  !-------------------------------------------------------------------------
  ! Assign physical coordinates to physical rectangle
  ! INPUT:
  !   coords = (/ xmin, ymin, zmin, xmax, ymax, zmax /), a vector of physical coordinates
  ! OUTPUT:
  !   rect = a rectangle in physical coordinates
  !-------------------------------------------------------------------------
  subroutine ob_assignCoordPhysToRectPhys(coords, rect)
    real(kind=DBL_KIND),dimension(OB_COORDS_SZ),intent(IN) :: coords
    type(t_obRectPhys),intent(OUT) :: rect
    type(t_obPointPhys) :: point
    real(kind=DBL_KIND),dimension(MX:MZ) :: cd
    cd = coords(MX:MZ)
    call ob_assignCoordPhysToPointPhys( cd, point )
    call ob_assignPointPhysToRectPhys( point, rect, 'L' )
    cd = coords(MZ+1:OB_COORDS_MAX)
    call ob_assignCoordPhysToPointPhys( cd, point )
    call ob_assignPointPhysToRectPhys( point, rect, 'R' )
  end subroutine ob_assignCoordPhysToRectPhys
  !-------------------------------------------------------------------------
  ! Exgtract a point from a rectangle.
  ! INPUT:
  !    rect = rectangle in overblock coordinates.
  !    leftOrRight = 'L' or 'R'.
  ! OUTPUT:
  !    point = a point in overblock coordinates.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  subroutine ob_extractPointPairFromRect(pointL, pointR, rect)
    type(t_obPoint),intent(OUT) :: pointL, pointR
    type(t_obRect),intent(IN) :: rect
    pointL%p(:) = rect%left(:)
    pointR%p(:) = rect%right(:)
    pointL%level = rect%level
    pointR%level = rect%level
  end subroutine ob_extractPointPairFromRect
  !-------------------------------------------------------------------------
  ! Extract a point from a rectangle in physical coordinates.
  ! INPUT:
  !    rect = rectangle in physical coordinates.
  !    leftOrRight = 'L' or 'R'.
  ! OUTPUT:
  !    point = a point in physical coordinates.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  subroutine ob_extractPointPhysPairFromRectPhys(pointL, pointR, rect)
    type(t_obPointPhys),intent(OUT) :: pointL, pointR
    type(t_obRectPhys),intent(IN) :: rect
    pointL%p(:) = rect%left(:)
    pointR%p(:) = rect%right(:)
  end subroutine ob_extractPointPhysPairFromRectPhys
  !-------------------------------------------------------------------------
  ! Extract a vector, of which elements are over-block coordinates of a given rectangle.
  ! INPUT:
  !    rect = rectangle in over-block coordinates.
  ! OUTPUT:
  !    coords = a vector of (/iobL, iobR, jobL, jobR, kobL, kobR/)
  !-------------------------------------------------------------------------
  subroutine ob_extractCoordFromRect(coords, rect)
    integer(kind=LLONG_KIND),dimension(OB_COORDS_SZ),intent(OUT) :: coords
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: point
    call ob_extractPointFromRect(point, rect, 'L')
    coords(MX:MZ) = point%p(:)
    call ob_extractPointFromRect(point, rect, 'R')
    coords(MZ+1:OB_COORDS_MAX) = point%p(:)
  end subroutine ob_extractCoordFromRect
  !-------------------------------------------------------------------------
  subroutine ob_extractLevelFromRect(level, rect)
    integer,intent(OUT) :: level
    type(t_obRect),intent(IN) :: rect
    level = rect%level
  end subroutine ob_extractLevelFromRect
  !-------------------------------------------------------------------------
  subroutine ob_extractLevelFromPoint(level, point)
    integer,intent(OUT) :: level
    type(t_obPoint),intent(IN) :: point
    level = point%level
  end subroutine ob_extractLevelFromPoint
  !-------------------------------------------------------------------------
  ! Extract a vector from rectangle in physical coordinates.
  ! INPUT:
  !    rectPhys = rectangle in physical coordinates
  ! OUTPUT:
  !    coordPhys = a vector of physical coordinates (/xL, yL, zL, xR, yR, zR/)
  !-------------------------------------------------------------------------
  subroutine ob_extractCoordPhysFromRectPhys(coordPhys, rectPhys)
    real(kind=DBL_KIND),dimension(OB_COORDS_SZ),intent(OUT) :: coordPhys
    type(t_obRectPhys),intent(IN) :: rectPhys
    type(t_obPointPhys) :: point
    call ob_extractPointPhysFromRectPhys(point, rectPhys, 'L')
    coordPhys(MX:MZ) = point%p(:)
    call ob_extractPointPhysFromRectPhys(point, rectPhys, 'R')
    coordPhys(MZ+1:OB_COORDS_MAX) = point%p(:)
  end subroutine ob_extractCoordPhysFromRectPhys
  !-------------------------------------------------------------------------
  ! get level-dependent over-block coordinates point
  ! when given by physical coordinates pointPhys and level.
  ! INPUT:
  !    pointPhys = physical coordinates
  !    level = grid level where over-block coordinates are measured.
  ! OUTPUT:
  !    pointOb  = corresponding index of cells (iob, job, kob).
  !    For example, a point x is included in the cell of iob
  !                 x - dx/2 <= iob <= x + dx/2.
  !-------------------------------------------------------------------------
  subroutine ob_PointPhys2PointOb( pointPhys, level, pointOb )
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN) :: level
    type(t_obPoint),intent(OUT) :: pointOb
    type(t_obPointPhys) :: compbox
    if ( .not. Initialized ) call ob_init
    call ob_extractPointPhysFromRectPhys(compbox, RectComputationBoxPhys, 'L')
    ! over-dependent coordinates assuming the offset is zero.
    pointOb%p(:) = (pointPhys%p(:)-compbox%p(:))/CellWidth(:,level)
    pointOb%level = level
  end subroutine ob_PointPhys2PointOb
  !-------------------------------------------------------------------------
  ! get a rectangle in over-block coordinates from that in physical coordinates
  ! INPUT:
  !   rectPhys = a rectangle in the physical coordinates
  !   level  = grid level in which overblock coordinates are measured
  ! OUTPUT:
  !   rectOb = a rectangle in the over-block coordinates
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! convert pointOb to pointPhys
  ! INPUT:
  !    pointOb = a point in over-block coordinates
  ! OUTPUT:
  !    pointPhys = a point in physical coordinates
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! Convert Rect from over-block coordinates to physical coordinates.
  ! INPUT:
  !   rectOb = rectangle in over-block coordinates.
  ! OUTPUT:
  !   rectPhys = rectangle in physical coordinates.
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! convert point by level
  ! INPUT:
  !    pointIn ....... a point to be converted
  !    levelOut ...... level of pointOut
  !    LR ............  left or right (optional)
  ! OUTPUT:
  !    pointOut ...... a converted point
  !-------------------------------------------------------------------------
  subroutine ob_point2PointByLevel(pointIn, pointOut, levelOut, LR)
    type(t_obPoint),intent(IN) :: pointIn
    type(t_obPoint),intent(OUT) :: pointOut
    integer,intent(IN) :: levelOut
    character(len=1),intent(IN),optional :: LR
    integer :: levelIn, offset
    integer(KIND=LLONG_KIND),parameter :: two =2
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
  !-------------------------------------------------------------------------
  ! convert point by level
  !-------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------
  ! compare rects
  !-------------------------------------------------------------------------
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

    if ( rectA%left(MX) == rectB%left(MX) .and. &
         rectA%left(MY) == rectB%left(MY) .and. &
         rectA%left(MZ) == rectB%left(MZ) .and. &
         rectA%right(MX) == rectB%right(MX) .and. &
         rectA%right(MY) == rectB%right(MY) .and. &
         rectA%right(MZ) == rectB%right(MZ) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_EQUIV
       return
    end if

    if ( rectA%left(MX) <= rectB%left(MX) .and. &
         rectA%left(MY) <= rectB%left(MY) .and. &
         rectA%left(MZ) <= rectB%left(MZ) .and. &
         rectA%right(MX) >= rectB%right(MX) .and. &
         rectA%right(MY) >= rectB%right(MY) .and. &
         rectA%right(MZ) >= rectB%right(MZ) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_OUTER
       return
    end if

    if ( rectA%left(MX) >= rectB%left(MX) .and. &
         rectA%left(MY) >= rectB%left(MY) .and. &
         rectA%left(MZ) >= rectB%left(MZ) .and. &
         rectA%right(MX) <= rectB%right(MX) .and. &
         rectA%right(MY) <= rectB%right(MY) .and. &
         rectA%right(MZ) <= rectB%right(MZ) .and. &
         rectA%level == rectB%level ) then
       comp = OB_RECT_INNER
       return
    end if

    if ( max(rectA%left(MX),rectB%left(MX)) <= min(rectA%right(MX),rectB%right(MX) ) .and. &
         max(rectA%left(MY),rectB%left(MY)) <= min(rectA%right(MY),rectB%right(MY) ) .and. &
         max(rectA%left(MZ),rectB%left(MZ)) <= min(rectA%right(MZ),rectB%right(MZ) ) ) then
         comp = OB_RECT_CROSS
         return
      end if

      comp = OB_RECT_DETOUCH
  end function ob_rectsComp
  !-------------------------------------------------------------------------
  ! compare rectPhys
  !-------------------------------------------------------------------------
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

    if ( rectA%left(MX) == rectB%left(MX) .and. &
         rectA%left(MY) == rectB%left(MY) .and. &
         rectA%left(MZ) == rectB%left(MZ) .and. &
         rectA%right(MX) == rectB%right(MX) .and. &
         rectA%right(MY) == rectB%right(MY) .and. &
         rectA%right(MZ) == rectB%right(MZ) ) then
       comp = OB_RECT_EQUIV
       return
    end if

    if ( rectA%left(MX) <= rectB%left(MX) .and. &
         rectA%left(MY) <= rectB%left(MY) .and. &
         rectA%left(MZ) <= rectB%left(MZ) .and. &
         rectA%right(MX) >= rectB%right(MX) .and. &
         rectA%right(MY) >= rectB%right(MY) .and. &
         rectA%right(MZ) >= rectB%right(MZ) ) then
       comp = OB_RECT_OUTER
       return
    end if

    if ( rectA%left(MX) >= rectB%left(MX) .and. &
         rectA%left(MY) >= rectB%left(MY) .and. &
         rectA%left(MZ) >= rectB%left(MZ) .and. &
         rectA%right(MX) <= rectB%right(MX) .and. &
         rectA%right(MY) <= rectB%right(MY) .and. &
         rectA%right(MZ) <= rectB%right(MZ) ) then
       comp = OB_RECT_INNER
       return
    end if

    if ( max(rectA%left(MX),rectB%left(MX)) <= min(rectA%right(MX),rectB%right(MX) ) .and. &
         max(rectA%left(MY),rectB%left(MY)) <= min(rectA%right(MY),rectB%right(MY) ) .and. &
         max(rectA%left(MZ),rectB%left(MZ)) <= min(rectA%right(MZ),rectB%right(MZ) ) ) then
         comp = OB_RECT_CROSS
         return
      end if

    comp = OB_RECT_DETOUCH
  end function ob_rectPhysComp
  !-------------------------------------------------------------------------
  ! Get logocal AND between rectA and rectB
  !-------------------------------------------------------------------------
  subroutine ob_rectAnd( rectA, rectB, rect )
    type(t_obRect),intent(IN) :: rectA, rectB
    type(t_obRect),intent(OUT) :: rect
    integer :: m
    if ( rectA%level /= rectB%level ) then
       print *, '*** ob_rectAnd : level is inconsistent between rectA and rectB', rectA%level, rectB%level
       stop
    endif
    do m = MX, MZ
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
    rect%level = rectA%level
  end subroutine ob_rectAnd
  !-------------------------------------------------------------------------
  ! Get logocal OR between rectA and rectB
  !-------------------------------------------------------------------------
  subroutine ob_rectOr( rectA, rectB, rect )
    type(t_obRect),intent(IN) :: rectA, rectB
    type(t_obRect),intent(OUT) :: rect
    integer :: m
    if ( rectA%level /= rectB%level ) then
       print *, '*** ob_rectAnd : level is inconsistent between rectA and rectB', rectA%level, rectB%level
       stop
    endif
    do m = MX, MZ
       rect%left(m) = min(rectA%left(m), rectB%left(m) )
       rect%right(m) = max(rectA%right(m), rectB%right(m) )
    enddo
    rect%level = rectA%level
  end subroutine ob_rectOr
  !-------------------------------------------------------------------------
  ! Get logocal AND between rectA and rectB
  !-------------------------------------------------------------------------
  subroutine ob_rectPhysAnd( rectA, rectB, rect )
    type(t_obRectPhys),intent(IN) :: rectA, rectB
    type(t_obRectPhys),intent(OUT) :: rect
    integer :: m
    do m = MX, MZ
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
  end subroutine ob_rectPhysAnd
  !-------------------------------------------------------------------------
  ! Get logocal OR between rectA and rectB
  !-------------------------------------------------------------------------
  subroutine ob_rectPhysOr( rectA, rectB, rect )
    type(t_obRectPhys),intent(IN) :: rectA, rectB
    type(t_obRectPhys),intent(OUT) :: rect
    type(t_obRectPhys) :: rectTmp
    integer :: m
    call ob_rectPhysAnd( rectA, rectB, rectTmp)
    do m = MX, MZ
       if ( rectTmp%left(m) > rectTmp%right(m) ) then
          print *, '*** ob_rectPhysOr: rectA and rectB are detouched'
       endif
    enddo
    do m = MX, MZ
       rect%left(m) = max(rectA%left(m), rectB%left(m) )
       rect%right(m) = min(rectA%right(m), rectB%right(m) )
    enddo
  end subroutine ob_rectPhysOr
  !-------------------------------------------------------------------------
  !  trim rectangle in physical coordinates.
  !-------------------------------------------------------------------------
  subroutine ob_rectPhysTrim ( rectIn, h, rectOut )
    type(t_obRectPhys),intent(IN) :: rectIn
    type(t_obRectPhys),intent(OUT) :: rectOut
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: h
!!$    real(kind=DBL_KIND) :: ave
    integer :: n
    do n = MX, MZ
       rectOut%left(n) = rectIn%left(n) + h(n)
       rectOut%right(n) = rectIn%right(n) - h(n)
    enddo
!!$    ! avoid crossing case
!!$    do n = MX, MZ
!!$       if ( rectOut%left(n) > rectOut%right(n) ) then
!!$          ave = ( rectOut%left(n) + rectOut%right(n) )/2.d0
!!$          rectOut%left(n) = ave
!!$          rectOut%right(n) = ave
!!$       endif
!!$    enddo
  end subroutine ob_rectPhysTrim
  !-------------------------------------------------------------------------
  !  trim rectangle in overbock coordinates
  !-------------------------------------------------------------------------
  subroutine ob_rectTrim ( rectIn, h, rectOut )
    type(t_obRect),intent(IN) :: rectIn
    type(t_obRect),intent(OUT) :: rectOut
    integer,dimension(MX:MZ),intent(IN) :: h
!!$    integer(kind=LLONG_KIND) :: ave
    integer :: n
    do n = MX, MZ
       rectOut%left(n) = rectIn%left(n) + h(n)
       rectOut%right(n) = rectIn%right(n) - h(n)
    enddo
    rectOut%level = rectIn%level
!!$    ! avoid crossing case
!!$    do n = MX, MZ
!!$       if ( rectOut%left(n) > rectOut%right(n) ) then
!!$          ave = ( rectOut%left(n) + rectOut%right(n) )/2.d0
!!$          rectOut%left(n) = ave
!!$          rectOut%right(n) = ave
!!$       endif
!!$    enddo
  end subroutine ob_rectTrim
  !-------------------------------------------------------------------------
  ! extend rectangle in physical coordinates.
  !-------------------------------------------------------------------------
  subroutine ob_rectPhysExtend(rectIn, h, rectOut)
    type(t_obRectPhys),intent(IN) :: rectIn
    type(t_obRectPhys),intent(OUT) :: rectOut
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: h
    real(kind=DBL_KIND),dimension(MX:MZ) :: minush
    minush = -h
    call ob_rectPhysTrim(rectIn, minush, rectOut)
  end subroutine ob_rectPhysExtend
  !-------------------------------------------------------------------------
  ! extend rectangle in overblock coordinates.
  !-------------------------------------------------------------------------
  subroutine ob_rectExtend(rectIn, h, rectOut)
    type(t_obRect),intent(IN) :: rectIn
    type(t_obRect),intent(OUT) :: rectOut
    integer,dimension(MX:MZ),intent(IN) :: h
    integer,dimension(MX:MZ) :: minush
    minush = -h
    call ob_rectTrim(rectIn, minush, rectOut)
  end subroutine ob_rectExtend
  !-------------------------------------------------------------------------
  ! test rect
  !-------------------------------------------------------------------------
  function ob_testRect( rect ) result( code )
    type(t_obRect),intent(IN) :: rect
    type(t_obPoint) :: pointL, pointR
    integer :: code
    call ob_extractPointFromRect(pointL, rect, 'L')
    call ob_extractPointFromRect(pointR, rect, 'R')
    if ( .not. ob_definedRect( rect ) ) then
       code = OB_RECT_UNDEF
    else if ( &
         ( pointL%p(MX) <= pointR%p(MX) ) .and. &
         ( pointL%p(MY) <= pointR%p(MY) ) .and. &
         ( pointL%p(MZ) <= pointR%p(MZ) ) ) then
       code = OB_RECT_VALID
    else
       code = OB_RECT_INVALID
    endif
  end function ob_testRect
  !-------------------------------------------------------------------------
  ! test rectPhys
  !-------------------------------------------------------------------------
  function ob_testRectPhys( rect ) result( code )
    type(t_obRectPhys),intent(IN) :: rect
    type(t_obPointPhys) :: pointL, pointR
    integer :: code
    call ob_extractPointPhysFromRectPhys(pointL, rect, 'L')
    call ob_extractPointPhysFromRectPhys(pointR, rect, 'R')
    if ( .not. ob_definedRectPhys( rect ) ) then
       code = OB_RECT_UNDEF
    else if ( &
         ( pointL%p(MX) <= pointR%p(MX) ) .and. &
         ( pointL%p(MY) <= pointR%p(MY) ) .and. &
         ( pointL%p(MZ) <= pointR%p(MZ) ) ) then
       code = OB_RECT_VALID
    else
       code = OB_RECT_INVALID
    endif
  end function ob_testRectPhys
  !-------------------------------------------------------------------------
  ! get rectangle of computational box in physical coordinates
  !-------------------------------------------------------------------------
  subroutine ob_computationBoxOfRectPhys( rectPhys )
    type(t_obRectPhys),intent(OUT) :: rectPhys
    if ( .not. Initialized ) call ob_init
    rectPhys = RectComputationBoxPhys
  end subroutine ob_ComputationBoxOfRectPhys
  !-------------------------------------------------------------------------
  ! get rectangle of computational box in overblock coordinates
  !-------------------------------------------------------------------------
  subroutine ob_computationBoxOfRect( rect, level )
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: level
    if ( .not. Initialized ) call ob_init
    rect = RectComputationBox(level)
  end subroutine ob_computationBoxOfRect
  !-------------------------------------------------------------------------
  ! get rectangle of computational box in vector of physical coordinates
  !     coordPhys = (/ xmin, ymin, zmin, xmax, ymax, zmax /)
  !-------------------------------------------------------------------------
  subroutine ob_computationBoxOfCoordPhys( coordPhys )
    real(kind=DBL_KIND),dimension(OB_COORDS_SZ),intent(OUT) :: coordPhys
    type(t_obRectPhys) :: rectPhys
    call ob_computationBoxOfRectPhys( rectPhys )
    call ob_extractCoordPhysFromRectPhys(coordPhys, rectPhys)
  end subroutine ob_computationBoxOfCoordPhys
  !-------------------------------------------------------------------------
  ! get (Igrid, Jgrid Kgrid) of a block where a given point of pointPhys exists.
  ! INPUT:
  !   level = grid level on which (Igrid, Jgrid, Kgrid) are measured.
  !   pointPhys = point in a physical coordinates.
  ! OUTPUT:
  !   ijkgrid = (Igrid, Jgrid Kgrid)
  !-------------------------------------------------------------------------
  subroutine ob_getIjkgridFromPointPhys(ijkgrid, level, pointPhys)
    integer,dimension(MX:MZ),intent(OUT) :: ijkgrid
    integer,intent(IN) :: level
    type(t_obPointPhys),intent(IN) :: pointPhys
    type(t_obPoint) :: point
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getIjkgridFromPoint(ijkgrid, point)
  end subroutine ob_getIjkgridFromPointPhys
  !-------------------------------------------------------------------------
  subroutine ob_getIjkgridFromPoint(ijkgrid, point)
    integer,dimension(MX:MZ),intent(OUT) :: ijkgrid
    type(t_obPoint),intent(IN) :: point
    ijkgrid(:) = (/ &
         point%p(MX) / (Imax-Imin+1), &
         point%p(MY) / (Jmax-Jmin+1), &
         point%p(MZ) / (Kmax-Kmin+1) &
         /)
  end subroutine ob_getIjkgridFromPoint
  !-------------------------------------------------------------------------
  ! same as ob_getIjkgridFromPointPhys() but arument is coordPhys instead of pointPhys
  !-------------------------------------------------------------------------
  subroutine ob_getIjkgridFromCoordPhys(ijkgrid, level, coordPhys)
    integer,dimension(MX:MZ),intent(OUT) :: ijkgrid
    integer,intent(IN) :: level
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: coordPhys
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coordPhys, pointPhys)
    call ob_getIjkgridFromPointPhys(ijkgrid, level, pointPhys)
  end subroutine ob_getIjkgridFromCoordPhys
  !-------------------------------------------------------------------------
  ! get a block and location inside the block, given by over-block coordinates.
  ! INPUT:
  !    point = over-block point
  ! OUTPUT:
  !    gid = grid id of block including a cell of 'point'.
  !    rank = rank of the gid
  !    i,j,k = indexes inside the block
  !    if block is undefined, gid is Undefi, rank is MPI_PROC_NULL.
  ! NOTE: this routine calls get_gid_from_ijkgrid (MPI local, but somewaht expensive)
  !-------------------------------------------------------------------------
  subroutine ob_getBlockFromPoint(point, gid, rank, i,j,k)
    type(t_obPoint),intent(IN) :: point
    integer,intent(OUT) :: gid, rank, i,j,k
    integer :: ig,jg,kg, level
    ! indexes of block (ijkgrid)
    ig = point%p(MX) / (Imax-Imin+1)
    jg = point%p(MY) / (Jmax-Jmin+1)
    kg = point%p(MZ) / (Kmax-Kmin+1)
    level = point%level
    ! indexes of cell inside the block
    i =  mod( point%p(MX), Imax-Imin+1 )
    j =  mod( point%p(MY), Jmax-Jmin+1 )
    k =  mod( point%p(MZ), Kmax-Kmin+1 )
    ! get gid and rank
    call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
  end subroutine ob_getBlockFromPoint
  !-------------------------------------------------------------------------
  ! get a block and location inside the block, given by physcical over-block coordinates.
  ! INPUT:
  !    pointPhys = over-block point in physical coordinates
  ! OUTPUT:
  !    gid = grid id of block including a cell of 'point'.
  !    rank = rank of the gid
  !    i,j,k = indexes inside the block
  !    if block is undefined, gid is Undefi, rank is MPI_PROC_NULL.
  ! NOTE: this routine calls get_gid_from_ijkgrid (MPI local, but somewaht expensive)
  !-------------------------------------------------------------------------
  subroutine ob_getBlockFromPointPhys(pointPhys, level, gid, rank, i,j,k)
    type(t_obPointPhys),intent(IN) :: pointPhys
    integer,intent(IN) :: level
    integer,intent(OUT) :: gid, rank, i,j,k
    type(t_obPoint) :: point
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getBlockFromPoint(point, gid, rank, i,j,k)
  end subroutine ob_getBlockFromPointPhys
  !-------------------------------------------------------------------------
  ! get a block and location inside the block, given by vector of physical coordinates
  ! INPUT:
  !   coords = a vector of physical coordinates
  ! OUTPUT:
  !    gid = grid id of block including a cell of 'point'.
  !    rank = rank of the gid
  !    i,j,k = indexes inside the block
  !    if block is undefined, gid is Undefi, rank is MPI_PROC_NULL.
  ! NOTE: this routine calls get_gid_from_ijkgrid (MPI local, but somewaht expensive)
  !-------------------------------------------------------------------------
  subroutine ob_getBlockFromCoordPhys(coords, level, gid, rank, i,j,k)
    real(kind=DBL_KIND),dimension(MX:MZ),intent(IN) :: coords
    integer,intent(IN) :: level
    integer,intent(OUT) :: gid, rank, i,j,k
    type(t_obPoint) :: point
    type(t_obPointPhys) :: pointPhys
    call ob_assignCoordPhysToPointPhys(coords, pointPhys)
    call ob_PointPhys2PointOb(pointPhys, level, point)
    call ob_getBlockFromPoint(point, gid, rank, i,j,k)
  end subroutine ob_getBlockFromCoordPhys
  !-------------------------------------------------------------------------
  ! get Point given by ijkgrid, level, i, j, k
  !-------------------------------------------------------------------------
  subroutine ob_getPointFromIjkgrid(point, ijkgrid, level, i, j, k)
    type(t_obPoint), intent(OUT) :: point
    integer,dimension(MX:MZ),intent(IN) :: ijkgrid
    integer,intent(IN) :: level, i, j, k
    point%p = (/ &
         ijkgrid(MX) * (Imax-Imin+1) + i, &
         ijkgrid(MY) * (Jmax-Jmin+1) + j, &
         ijkgrid(MZ) * (Kmax-Kmin+1) + k /)
    point%level = level
  end subroutine ob_getPointFromIjkgrid
  !-------------------------------------------------------------------------
  ! get Point given by gid, i, j, k
  ! INPUT:
  !    gid, rank, i, j, k
  ! OUTPUT:
  !    point = point included within a block
  !    retcode = Return code of 1 or Undefi.
  !              Return 1 when a block is found.
  !              Return Undefi when a block is not found.
  ! This subroutine is for serial process (MPI local)
  !-------------------------------------------------------------------------
  subroutine ob_getPointFromBlockSerial(point, gid, i, j, k, retcode)
    type(t_obPoint), intent(OUT) :: point
    integer,intent(IN) :: gid, i, j, k
    integer,intent(OUT) :: retcode
    integer :: level
    integer(kind=LLONG_KIND),dimension(MX:MZ) :: coords
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
  !-------------------------------------------------------------------------
  ! Same as ob_getPointFromBlockSerial but use MPI communication (MPI global).
  ! All nodes have the same results.
  !-------------------------------------------------------------------------
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
    call mpi_bcast( point, 1, MPI_OB_POINT,  rank, MPI_COMM_WORLD, ierr )
  end subroutine ob_getPointFromBlock
  !-------------------------------------------------------------------------
  ! get rect of block, given by (ig, jg, kg)
  !-------------------------------------------------------------------------
  subroutine ob_getRectFromIjkgrid(rect, ijkgrid, level)
    type(t_obRect),intent(OUT) :: rect
    integer,dimension(MX:MZ),intent(IN) :: ijkgrid
    integer,intent(IN) :: level
    type(t_obPoint) :: pointL, pointR
    call ob_getPointFromIjkgrid(pointL, ijkgrid, level, Imin, Jmin, Kmin)
    call ob_getPointFromIjkgrid(pointR, ijkgrid, level, Imax, Jmax, Kmax)
    call ob_assignPointPairToRect(pointL, pointR, rect)
  end subroutine ob_getRectFromIjkgrid
  !-------------------------------------------------------------------------
  ! get rect of block, given by gid
  !-------------------------------------------------------------------------
  subroutine ob_getRectFromGid(rect, gid)
    type(t_obRect),intent(OUT) :: rect
    integer,intent(IN) :: gid
    integer,dimension(MX:MZ) :: ijkgrid
    ijkgrid = (/Igrid(gid), Jgrid(gid), Kgrid(gid)/)
    call ob_getRectFromIjkgrid(rect, ijkgrid, get_level(gid))
  end subroutine ob_getRectFromGid
  !-------------------------------------------------------------------------
  ! whether point within rect
  !-------------------------------------------------------------------------
  function ob_PointWithinRect(point, rect) result(bool)
    type(t_obPoint),intent(IN) :: point
    type(t_obRect),intent(IN) :: rect
    logical :: bool
    type(t_obPoint) :: pL, pR
    call ob_extractPointFromRect(pL, rect, 'L')
    call ob_extractPointFromRect(pR, rect, 'R')
    if ( pL%p(MX) <= point%p(MX) .and. point%p(MX) <= pR%p(MX) .and. &
         pL%p(MY) <= point%p(MY) .and. point%p(MY) <= pR%p(MY) .and. &
         pL%p(MZ) <= point%p(MZ) .and. point%p(MZ) <= pR%p(MZ) ) then
       bool = .true.
    else
       bool = .false.
    endif
  end function ob_PointWithinRect
  !-------------------------------------------------------------------------
  ! whether pointPhys within rectPhys
  !-------------------------------------------------------------------------
  function ob_PointPhysWithinRectPhys(point, rect) result(bool)
    type(t_obPointPhys),intent(IN) :: point
    type(t_obRectPhys),intent(IN) :: rect
    logical :: bool
    type(t_obPointPhys) :: pL, pR
    call ob_extractPointPhysFromRectPhys(pL, rect, 'L')
    call ob_extractPointPhysFromRectPhys(pR, rect, 'R')
    if ( pL%p(MX) <= point%p(MX) .and. point%p(MX) <= pR%p(MX) .and. &
         pL%p(MY) <= point%p(MY) .and. point%p(MY) <= pR%p(MY) .and. &
         pL%p(MZ) <= point%p(MZ) .and. point%p(MZ) <= pR%p(MZ) ) then
       bool = .true.
    else
       bool = .false.
    endif
  end function ob_PointPhysWithinRectPhys
  !-------------------------------------------------------------------------
  ! Initialize this module
  !-------------------------------------------------------------------------
  subroutine ob_init
    use mpilib
    real(kind=DBL_KIND),dimension(:),pointer :: xp, yp, zp
    integer :: gid, rank, level
    real(kind=DBL_KIND) :: xyzmin(MX:MZ), xyzmax(MX:MZ), coords(OB_COORDS_SZ)
    real(kind=DBL_KIND),save :: xmin, ymin, zmin, xmax, ymax, zmax
    type(t_obPoint) :: point
    integer(KIND=LLONG_KIND),parameter :: two =2
    integer(KIND=LLONG_KIND) :: power
    level = 0
    ! -----------------------------------
    ! RectComputationBoxPhys
    ! -----------------------------------
    ! xmin, ymin, zmin
    call get_gid_from_ijkgrid(0,0,0,level,gid,rank)
    myrank = get_myrank()
    if ( myrank == rank ) then
       xp => get_Xp(gid)
       yp => get_Yp(gid)
       zp => get_Zp(gid)
       xmin = xp(Imin) - CellWidth(MX,level)/2.d0
       ymin = yp(Jmin) - CellWidth(MY,level)/2.d0
       zmin = zp(Kmin) - CellWidth(MZ,level)/2.d0
       xyzmin(:) = (/xmin, ymin, zmin/)
    endif
    call mpi_bcast(xyzmin, size(xyzmin), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    ! Lower left corner of compunational box
    xmin = xyzmin(MX)
    ymin = xyzmin(MY)
    Zmin = xyzmin(MZ)

    ! xmax, ymax, zmax
    gid = GidBase( ubound(GidBase,1),ubound(GidBase,2),ubound(GidBase,3) )
    rank = RankBase( ubound(GidBase,1),ubound(GidBase,2),ubound(GidBase,3) )
    if ( myrank == rank ) then
       xp => get_Xp(gid)
       yp => get_Yp(gid)
       zp => get_Zp(gid)
       xmax = xp(Imax) + CellWidth(MX,level)/2.d0
       ymax = yp(Jmax) + CellWidth(MY,level)/2.d0
       zmax = zp(Kmax) + CellWidth(MZ,level)/2.d0
       xyzmax(:) = (/xmax, ymax, zmax/)
    endif
    call mpi_bcast(xyzmax, size(xyzmax), MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
    ! Upper right corner of compunational box
    xmax = xyzmax(MX)
    ymax = xyzmax(MY)
    zmax = xyzmax(MZ)

    coords = (/xmin, ymin, zmin, xmax, ymax, zmax/)
    call ob_assignCoordPhysToRectPhys(coords, RectComputationBoxPhys)

    ! -----------------------------------
    ! RectComputationBox
    ! -----------------------------------
    do level = Lmin, Lmax
       point%level = level
       point%p(:) = 0
       call ob_assignPointToRect(point, RectComputationBox(level), 'L')
       power = two**(level-Lmin)
       point%p(MX) = (NI)*(NGI_BASE)*power-1
       point%p(MY) = (NJ)*(NGJ_BASE)*power-1
       point%p(MZ) = (NK)*(NGK_BASE)*power-1
       call ob_assignPointToRect(point, RectComputationBox(level), 'R')
    end do

    ! -----------------------------------------------
    ! define new type of MPI for overBlockCoordinates
    ! -----------------------------------------------
    call ob_Make_MPI_struct_point
    call ob_Make_MPI_struct_rect
    call ob_Make_MPI_struct_pointPhys
    call ob_Make_MPI_struct_rectPhys
    Initialized = .true.

    if ( get_myrank() == PRIMARY_RANK ) print *, 'initialize overBlockCoordinates'

  contains
    subroutine ob_Make_MPI_struct_point
      integer,dimension(2) :: blocklen, disp, type
      integer :: sz
      type(t_obPoint) :: point
      sz = size(point%p)
      blocklen = (/sz, 1/)
      disp = (/0, blocklen(1)*(LLONG_KIND)/)
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
      disp = (/0, blocklen(1)*(LLONG_KIND), blocklen(1)*(LLONG_KIND)+blocklen(2)*(LLONG_KIND) /)
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
