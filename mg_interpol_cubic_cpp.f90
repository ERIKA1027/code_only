module mg_interpol_cubic
  use mg_data
  implicit none
  private
  type t_buf
     real(kind=8),dimension(:,:,:,:),pointer :: bufx => null(), bufxy => null()
  end type t_buf
  type(t_buf),save,dimension(:),pointer :: FLV => null() 
  type(t_buf),save,dimension(:),pointer :: FLS => null() 
  logical,save :: Initialized = .FALSE.
  public :: mg_interp_cubic, mg_interp_cubic_finalize, mg_interp_cubic_check
contains
  subroutine mg_interp_cubic_init
    integer :: mglev
    integer :: ifs, jfs, kfs, ife, jfe, kfe
    integer :: ics, jcs, kcs, ice, jce, kce
    if (Initialized) return
    Initialized = .TRUE.
    allocate( FLV(MG_LevelMin: MG_LevelMax-1) )
    allocate( FLS(MG_LevelMin: MG_LevelMax-1) )
    do mglev = lbound(FLV, 1), ubound(FLV, 1)
       call mg_get_gridsize (mglev, ifs, jfs, kfs, ife, jfe, kfe) 
       call mg_get_gridsizeGh(mglev+1, ics, jcs, kcs, ice, jce, kce) 
       allocate(FLV(mglev)%bufx (ifs:ife,jcs:jce,kcs:kce,Mmin:Mmax))
       allocate(FLV(mglev)%bufxy(ifs:ife,jfs:jfe,kcs:kce,Mmin:Mmax))
       allocate(FLS(mglev)%bufx (ifs:ife,jcs:jce,kcs:kce,Mmin:Mmin))
       allocate(FLS(mglev)%bufxy(ifs:ife,jfs:jfe,kcs:kce,Mmin:Mmin))
    end do
  end subroutine mg_interp_cubic_init
  subroutine mg_interp_cubic_finalize
    integer :: mglev
    if (.not. Initialized) return
    Initialized = .FALSE.
    do mglev = lbound(FLV,1), ubound(FLV,1)
       deallocate( &
            FLV(mglev)%bufx, &
            FLV(mglev)%bufxy, &
            FLS(mglev)%bufx, &
            FLS(mglev)%bufxy )
       nullify( &
            FLV(mglev)%bufx, &
            FLV(mglev)%bufxy, &
            FLS(mglev)%bufx, &
            FLS(mglev)%bufxy )
    end do
    deallocate(FLV, FLS)
    nullify(FLV, FLS)
  end subroutine mg_interp_cubic_finalize
  function mg_interp_cubic_check(mglev) result(bool_cubic)
    use mg_data, only : mg_get_gridsizeGh
    integer,intent(IN) :: mglev
    logical :: bool_cubic
    integer,parameter :: MINBLOCKSIZE = 4
    integer :: ics,jcs,kcs,ice,jce,kce
    call mg_get_gridsizeGh(mglev+1, ics,jcs,kcs,ice,jce,kce)
    if ( &
         ice-ics+1 < MINBLOCKSIZE .or. &
         jce-jcs+1 < MINBLOCKSIZE .or. &
         kce-kcs+1 < MINBLOCKSIZE ) then
       bool_cubic = .FALSE.
    else
       bool_cubic = .TRUE.
    end if
  end function mg_interp_cubic_check
  subroutine mg_interp_cubic(mglev, iuf, iuc)
    use interpolation
    use fmg_data, only : fmg_isVector
    integer,intent(IN) :: mglev, iuf, iuc
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc
    real(kind=8),dimension(:,:,:,:),pointer :: bufx, bufxy
    type(t_buf),dimension(:),pointer :: flp
    integer :: if, jf, kf, ic, jc, kc, m, nshift
    integer :: ic0, ic1, ic2, ic3, icP, icN, ifmin
    integer :: jc0, jc1, jc2, jc3, jcP, jcN, jfmin
    integer :: kc0, kc1, kc2, kc3, kcP, kcN, kfmin
    integer :: lf, lc, ifs, jfs, kfs, ife, jfe, kfe
    integer :: ibs, ibe, jbs, jbe, kbs, kbe
    real(kind=8) :: xc, yc, zc, a0, a1, a2, a3
    call mg_interp_cubic_init
    if (fmg_isVector(iuc)) then
       flp => FLV
    else
       flp => FLS
    end if
    bufx => flp(mglev)%bufx
    bufxy => flp(mglev)%bufxy
    lf = mglev
    lc = mglev+1
    call mg_get_gridsize(lf, ifs, jfs, kfs, ife, jfe, kfe)
    uf => mg_get_arrp(lf,iuf)
    uc => mg_get_arrp(lc,iuc)
    call interp_tricubic(uf, uc, bufx, bufxy, &
         ifs, jfs, kfs, ife, jfe, kfe, &
         ifs, jfs, kfs)
  end subroutine mg_interp_cubic
end module mg_interpol_cubic
