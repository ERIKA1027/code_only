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
    use io
#ifndef MP
    use eos, only : Rhocr
#endif !MP
    use analysis, only : RhoMax
#ifdef SINKPARTICLE
    use sinkParticle
#endif !SINKPARTICLE
    use writeSnap
    use uniformgrid, only : uniformgrid_write
    use grid, only : Lmin, LevelMax
    use modelParameter, only : MP_Boxsize ! KS ADDED
    use mpilib ! KS ADDED

    integer :: level, np
    real(kind=DBL_KIND) :: hw

    if (.not. bool_output() ) return

    ! np ... number of sinkParticle
    np = sp_getNparticle()    

    ! output whole computational domain
    ! if (get_myrank() == PRIMARY_RANK) then
    !    print '(/,A,I0)', "writeSnap_whole: level = ", -1
    ! end if
    ! call writeSnap_whole

    ! output central box
    call output_centralbox

    if (np > 0) then
       ! output regions around sink particles
       call output_sinkParticles ! (comment out for rad_test)
    end if

    if (np == 0) then !sink粒子形成後は密度最大点の周りでプロットする意味はあまり無いのでスキップ
       ! output around density maximum
       call output_denseRegion()
    end if

  end subroutine output_data

  ! -----------------------------------------------------------------
  ! return ture for when a output timing is comming.
  ! -----------------------------------------------------------------
  function bool_output() result(bool)
    use grid
    use eos
    use analysis, only : RhoMax
    use modelParameter, only : MP_Dstep
    use primordial
    logical :: bool
    real(kind=DBL_KIND),parameter :: logrho_skip = 1d0 ! interval for output in unit of log(rho_max)
    real(kind=DBL_KIND),save :: logrhoio
    !real(kind=DBL_KIND),save :: rhomax_prev = Huge(rhomax_prev)
    real(kind=DBL_KIND),save :: rhomax_prev = Huge(rhomax_prev)
    logical,save :: bool_output_initialized = .false.

    
    !initialization
    if (.not. bool_output_initialized) then
       rhomax_prev = RhoMax
       bool_output_initialized = .True.
    end if


    ! -------------
    ! for debug 
    ! -------------
    bool = .false.
    !MP_Dstep毎にデータを出力
    if (level_sync() == Lmin .and. mod(Step(Lmin),MP_Dstep) == 0 ) then !KS MODIFIED
       bool = .true.
    endif

    !密度が(logrho_skip)桁上がる度にデータを出力
    if (level_sync() == Lmin .and. RhoMax > rhomax_prev * 10**logrho_skip) then !KS MODIFIED
       bool = .true.
       rhomax_prev = rhomax_prev * 10**logrho_skip
    endif


    ! !最初の20ファイルはMP_Dstep/5毎にデータを出力 (tentative)
    ! if (Step(Lmin)/MP_Dstep < 20) then !KS MODIFIED
    !    if (level_sync() == Lmin .and. mod(Step(Lmin),MP_Dstep/5) == 0 ) then !KS MODIFIED
    !       bool = .true.
    !    endif
    ! endif   
    
    return
  end function bool_output

  ! -----------------------------------------------------------------
  ! Output central box.
  ! From level = Lmin-1 to LevelMax
  ! File name is cb step . level . d
  ! -----------------------------------------------------------------
  subroutine output_centralbox()
    use grid, only : Lmin, LevelMax
    use parameter
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use mpilib

    
    character(len=2),parameter :: prefix = 'cb'
    real(kind=DBL_KIND) :: halfwidth
    !integer,parameter :: NIug=32 ! data points in each direction (KS MODIFIED, ~ 2 MB/file)
    integer,parameter :: NIug=64 ! KS DEBUG
    integer,parameter :: l_range=20 ! num of levels to be output
    integer :: level
    

    ! do level = Lmin-1,LevelMax
    do level = Lmin-2,LevelMax       !level = -2まで許してみる
       if (level < LevelMax - l_range + 1) cycle
       halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction       
       if (halfwidth > MP_boxsize) cycle
       if (get_myrank() == PRIMARY_RANK) then
          print '(/,A,I0,A,(1P1E9.2),A)', "output_centralbox: level = ", &
               level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
       end if
       call uniformgrid_write(-halfwidth,-halfwidth,-halfwidth,halfwidth,halfwidth,halfwidth, level, interpolate=.false.,prefix=prefix)
    enddo
  end subroutine output_centralbox

  ! -----------------------------------------------------------------
  ! Output regions around sink particles.
  ! From level = LevelMax-5 to LevelMax
  ! File name is sp. pid . step . level . d
  ! -----------------------------------------------------------------
  subroutine output_sinkParticles()
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use unit, only : Unit_msun
    use sinkParticle
    use mpilib
    use string, only : CHARLEN, num2char, concat

    ! integer,parameter :: NIug=32, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 2 MB/file)
    integer,parameter :: NIug=64, NJug=NIug, NKug=NIug ! minimum resolution (KS MODIFIED, ~ 16 MB/file)    
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=DBL_KIND) :: halfwidth
    integer :: level, np, n, dummy
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    character(len=CHARLEN) :: prefix
    integer,dimension(:),allocatable :: pid
    real(kind=DBL_KIND),dimension(:),allocatable :: pmass
    real(kind=DBL_KIND),dimension(:,:),allocatable :: pr

    ! for barycenter output
    real(kind=DBL_KIND) :: sp_mtot, sp_xbc, sp_ybc, sp_zbc

    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )

    ! pr ... location of sinkParticle
    ! np ... number of sinkParticle
    np = sp_getNparticle()

    if (np > 0) then     ! if sink particles exist
       allocate(pid(np), pmass(np),  pr(MX:MZ, np))
       call sp_sinkdata2array(dummy, pmass, pr=pr, pid=pid) ! from sink particle
    else
       return ! sink粒子存在しなければreturn
    end if

    !--------------------- 各シンク粒子周りでの出力 --------------------!
    do n = 1, np
       !massが小さい場合はスキップ (ひとまず1 M_sunを境界にする)
       if (pmass(n)*Unit_msun < 1d0) then
          if (get_myrank() == PRIMARY_RANK) then
             print '(/,A,I0,A,1P1E9.2,A)', "pid = ", pid(n), ", pmass = ",pmass(n)*Unit_msun, &
             ", skip sinkParticle with mass < 1 Msun"
          end if
          cycle
       end if

       do level = LevelMax-5,LevelMax !最大レベルから6レベル分をプロット
          if (level <= -2) cycle 
          halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction in level -1
          if (get_myrank() == PRIMARY_RANK) then
             print '(/,A,I0,A,1P1E9.2,A,I0,A,(1P1E9.2),A)', "output_sinkParticles: pid = ", pid(n), " (pmass = ", &
             pmass(n)*Unit_msun ,"), level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
          end if

          xmin = coordPhys(MX)
          ymin = coordPhys(MY)
          zmin = coordPhys(MZ)
          xmax = coordPhys(MZ+1+MX)
          ymax = coordPhys(MZ+1+MY)
          zmax = coordPhys(MZ+1+MZ)
          ! 領域の AND をとる
          xmin = max(xmin, pr(MX,n)-halfwidth)
          ymin = max(ymin, pr(MY,n)-halfwidth)
          zmin = max(zmin, pr(MZ,n)-halfwidth)
          xmax = min(xmax, pr(MX,n)+halfwidth)
          ymax = min(ymax, pr(MY,n)+halfwidth)
          zmax = min(zmax, pr(MZ,n)+halfwidth)

          ! speicfy prefix
          prefix = 'sp.'
          prefix = concat(concat(prefix, num2char(pid(n))),'.')

          call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)

       end do ! level
    end do ! n

    ! ----------------------- sink粒子の重心周りでの出力 --------------------------- !
    !重心の取得
    sp_mtot = 0d0
    sp_xbc = 0d0
    sp_ybc = 0d0
    sp_zbc = 0d0
    do n = 1, np
       sp_mtot = sp_mtot + pmass(n)
       sp_xbc = sp_xbc + pr(MX,n) * pmass(n)
       sp_ybc = sp_ybc + pr(MY,n) * pmass(n)
       sp_zbc = sp_zbc + pr(MZ,n) * pmass(n)
    end do
    sp_xbc = sp_xbc/sp_mtot
    sp_ybc = sp_ybc/sp_mtot
    sp_zbc = sp_zbc/sp_mtot    

    !データの出力
    do level = Lmin-2,LevelMax !最大レベルから最小レベルまで出力
       if (level <= -2) cycle 
       halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction in level -1
       if (get_myrank() == PRIMARY_RANK) then
          print '(/,A,1P1E9.2,A,I0,A,(1P1E9.2),A)', "output_sinkParticles: barycenter (mtot = ", &
               sp_mtot*Unit_msun ,"), level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
       end if

       xmin = coordPhys(MX)
       ymin = coordPhys(MY)
       zmin = coordPhys(MZ)
       xmax = coordPhys(MZ+1+MX)
       ymax = coordPhys(MZ+1+MY)
       zmax = coordPhys(MZ+1+MZ)
       ! 領域の AND をとる
       xmin = max(xmin, sp_xbc-halfwidth)
       ymin = max(ymin, sp_ybc-halfwidth)
       zmin = max(zmin, sp_zbc-halfwidth)
       xmax = min(xmax, sp_xbc+halfwidth)
       ymax = min(ymax, sp_ybc+halfwidth)
       zmax = min(zmax, sp_zbc+halfwidth)

       ! speicfy prefix
       prefix = 'bc'
       call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)

    end do ! level
       

  end subroutine output_sinkParticles

  
  ! -----------------------------------------------------------------
  ! Output box around density maximum.
  ! From level = LevelMax - 5 to LevelMax
  ! File name is cb step . level . d
  ! -----------------------------------------------------------------
  subroutine output_denseRegion()
    use grid, only : Lmin, LevelMax
    use parameter
    use modelParameter, only : MP_Boxsize
    use mpilib
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write    

    character(len=2),parameter :: prefix = 'dn'
    real(kind=DBL_KIND) :: halfwidth
    !integer,parameter :: NIug=32 ! data points in each direction (KS MODIFIED, ~ 2 MB/file)
    integer,parameter :: NIug=64 ! KS DEBUG
    integer :: level
    real(kind=DBL_KIND) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=DBL_KIND) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)    

    ! 最大密度の位置を取得
    call rhomaxPosition(xp, yp, zp)

    ! *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** !
    !
    ! slightly shifted to avoid occurence of different cell size in x/y/z directions
    !                                       (NEED TO BE CHECKED)
    !
    xp=xp*1.0001d0
    yp=yp*1.0001d0
    zp=zp*1.0001d0
    ! *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** !    

    !    do level = LevelMax-5,LevelMax !最大レベルから6レベル分をプロット
    do level = Lmin-2,LevelMax !最大レベルから最小レベルまでプロット
       if (level <= -2) cycle
       halfwidth = MP_Boxsize / 2.**(level-Lmin) / (NI*NGI_BASE/dble(NIug))  ! output region  with 32 data points/direction in level -1
       if (halfwidth > MP_boxsize) cycle       
       if (get_myrank() == PRIMARY_RANK) then
          print '(/,A,I0,A,(1P1E9.2),A)', "output_rhomaxbox: level = ", level, &
               ", size = ", halfwidth/MP_Boxsize, " x comp. region"
       end if
    ! upper bound of computational box
    call ob_computationBoxOfCoordPhys( coordPhys )
    xmin = coordPhys(MX)
    ymin = coordPhys(MY)
    zmin = coordPhys(MZ)
    xmax = coordPhys(MZ+1+MX)
    ymax = coordPhys(MZ+1+MY)
    zmax = coordPhys(MZ+1+MZ)

    ! 領域の AND をとる
    xmin = max(xmin, xp-halfwidth)
    ymin = max(ymin, yp-halfwidth)
    zmin = max(zmin, zp-halfwidth)
    xmax = min(xmax, xp+halfwidth)
    ymax = min(ymax, yp+halfwidth)
    zmax = min(zmax, zp+halfwidth)

    ! speicfy prefix
    call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level,interpolate=.false., prefix=prefix)
    enddo


  contains
    !---------------------------------------------------------------------
    ! get a rectangle region which a given level covers
    !---------------------------------------------------------------------
    subroutine rhomaxPosition(xp, yp, zp)
      use mpilib
      use grid, only : get_Ucomp, get_Xp, get_Yp, get_Zp, Undefi, Gidmin, GidListMax, GidList, &
           Imin, Jmin, Kmin, Imax, Jmax, Kmax, Time, Step

      use string, only : CHARLEN, num2char, concat          ! KS ADDED
      use unit, only : Unit_yr, Unit_au, Unit_rho, cgs_mh   ! KS ADDED
      use primordial,only: yhe                              ! KS ADDED
      use io_util, only : read_env                          ! KS ADDED
      
      real(kind=DBL_KIND),intent(OUT) :: xp, yp, zp
      real(kind=DBL_KIND),dimension(:),pointer :: x, y, z
      real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
      real(kind=DBL_KIND),dimension(0:NPE-1) :: rhomaxNode, rhomaxNoder
      real(kind=DBL_KIND),dimension(MX:MZ) ::  buf
      real(kind=DBL_KIND) :: rhomax, xmax, ymax, zmax
      integer :: n, gid, gidmax, maxl(MX:MZ), nodemax(1)

      !KS ADDED
      integer,parameter :: LUN = 11
      character(len=CHARLEN) :: file, dir
      character(len=CHARLEN),parameter :: fnlog='dnLog.dat' !KS ADDED
      real(kind=DBL_KIND) :: xnHmax
      
      ! find gid of maximum rho
      myrank = get_myrank()
      gidmax = Undefi
      rhomaxNode(:) = 0.d0
      do n = Gidmin, GidListMax(LevelMax)
         gid = GidList(n, LevelMax)
         rho => get_Ucomp(MRHO, gid)
         rhomax = maxval(rho(ARRAYSIZE_IJK))
         if (rhomax > rhomaxNode(myrank)) then
            rhomaxNode(myrank) = rhomax
            gidmax = gid
         endif
      enddo
      ! coordinates of rhomax for each node
      if ( gidmax /= Undefi ) then
         rho => get_Ucomp(MRHO, gidmax)
         x => get_Xp(gidmax)
         y => get_Yp(gidmax)
         z => get_Zp(gidmax)
         maxl = maxloc(rho(ARRAYSIZE_IJK))
         maxl = maxl - 1 + (/Imin, Jmin, Kmin/)
         xmax = x(maxl(MX))
         ymax = y(maxl(MY))
         zmax = z(maxl(MZ))
      endif
      ! find node of maximum rho
      call mpi_allreduce(rhomaxNode, rhomaxNoder, size(rhomaxNode), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      nodemax = maxloc(rhomaxNoder) - 1
      buf = (/ xmax, ymax, zmax /)
      call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, nodemax(1), MPI_COMM_WORLD, ierr)
      xp = buf(MX) ; yp = buf(MY) ; zp = buf(MZ)

      !-------------- KS ADDED --------------!
      ! write position and rhomax to a file
      !--------------------------------------!      
      if (get_myrank() == PRIMARY_RANK) then
         call read_env('DIR', dir)
         file = concat(dir,fnlog)
         open(LUN, file=file, position='APPEND')
         xnHmax = maxval(rhomaxNoder)*Unit_rho/((1.d0 + 4.d0*yHe)*cgs_mh)      ! 密度最大点での水素原子核の数密度
         write(LUN, '(I10, 5(1PE17.9))') Step(Lmin), Time(Lmin)*Unit_yr,&
              xp*Unit_au, yp*Unit_au, zp*Unit_au, xnHmax
         print '(/,A,/,I10, 5(1PE17.9))', '(rhomaxPosition) Step, Time, pos, nH_max:', Step(Lmin), Time(Lmin)*Unit_yr, &
              xp*Unit_au, yp*Unit_au, zp*Unit_au, xnHmax
         call flush(LUN)
         close(LUN)
      end if
      !--------------------------------------!      
    end subroutine rhomaxPosition

  end subroutine output_denseRegion

end module outputdata
