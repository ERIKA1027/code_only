module writeSnap
  use uniformgrid
  implicit none
  private
  public :: writeSnap_whole
contains
  subroutine writeSnap_whole
    use grid, only : Lmin
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    integer,parameter :: SampleRate = 2 
    integer :: baseLevel
    real(KIND=8) :: xmin, ymin, zmin, xmax, ymax, zmax, coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    call ob_computationBoxOfCoordPhys( coordPhys )
    xmin = coordPhys(0)
    ymin = coordPhys(1)
    zmin = coordPhys(2)
    xmax = coordPhys(2 +1+0)
    ymax = coordPhys(2 +1+1)
    zmax = coordPhys(2 +1+2)
    baseLevel = Lmin - int( log10(dble(SampleRate))/log10(2.d0) + 0.5)
    if ( bool_skip( baseLevel ) ) return
    call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, baseLevel,interpolate=.true.)
  end subroutine writeSnap_whole
  function bool_skip( baseLevel ) result( bool )
    use grid, only : LevelMax
    integer,intent(IN) :: baseLevel
    logical :: bool
    integer,parameter :: DlevelMax = 3
    bool = .false. 
    if ( baseLevel - LevelMax > DlevelMax ) bool = .true.
  end function bool_skip
end module writeSnap
