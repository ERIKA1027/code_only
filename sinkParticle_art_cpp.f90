module sinkParticle
  use overBlockCoordinates
  use grid
  use string, only : CHARLEN
  use unit
  implicit none
  private
  type t_spParticle
     integer :: pid 
     real(kind=8) :: mass 
     real(kind=8) :: dmass 
     real(kind=8) :: r(0:2) 
     real(kind=8) :: dr(0:2) 
     real(kind=8) :: v(0:2) 
     real(kind=8) :: dv(0:2) 
     real(kind=8) :: dp(0:2) 
     real(kind=8) :: ds(0:2) 
     real(kind=8) :: s(0:2) 
     real(kind=8) :: t_prev 
     real(kind=8) :: dm_disk 
     real(kind=8) :: dJ_disk(0:2) 
     real(kind=8) :: mdot_disk 
     real(kind=8) :: J_disk(0:2) 
     real(kind=8) :: t_crt 
     integer :: lev 
     real(kind=8) :: h(0:2) 
     type(t_spParticle),pointer :: next => null()
  end type t_spParticle
  type t_spPos
     real(8),dimension(0:2) :: r
     type(t_spPos),pointer :: next => null()
  end type t_spPos
  real(kind=8),save :: mass_disk=0.0
  real(kind=8),save :: mass_subdisk=0.0
 integer,save :: Nparticle = 0
  integer,save :: PidMax = -1
  integer,save :: sp_Level = Undefi
  logical,save :: Initialized = .false.
  type(t_spParticle),save,pointer :: Particle => null() 
  real(kind=8),save :: sp_Rhocr, sp_RhocrCreate 
  real(kind=8),save :: sp_Cs 
  real(kind=8),save :: sp_SinkRadius 
  real(kind=8),save :: sp_SofteningRadius 
  real(kind=8),save :: sp_MaskRadius 
  real(kind=8),parameter :: sp_f_rAcc = 4.d0
  real(kind=8),save :: sp_DtimeCFL
  real(kind=8),parameter :: sp_CGS = 1.d-1
  real(kind=8),parameter :: sp_CVS = 5.d-1 
  character(len=CHARLEN),parameter :: FILENAME_LOG='logSinkparticle'
  character(len=CHARLEN),parameter :: FILENAME='sinkparticle.d'
  real(kind=8),parameter :: dt_acc = 3d2 
  real(kind=8), dimension(0:2), save :: hh_sp
  public :: sp_update, sp_write, sp_read, sp_writeLog, sp_restrictCFL, sp_refineCond, sp_getLevel, &
       sp_sinkdata2array, sp_getNparticle, sp_getRhocr, sp_newParticle, sp_gravityOfParticle, &
       sp_getSinkRadius, sp_is_inside_sink, sp_refineCond_KS 
contains
  subroutine sp_update
    use parameter
    use modelParameter
    integer, save :: ifirst = 0
    real(kind=8),dimension(0:2) :: pos, v
    real(kind=8) :: a, vr, mass
    call sp_init
    if (Nparticle == 0) return
    call get_plev
    call sp_accretion_init
   call sp_accretion
    call sp_accretion_comm
    call sp_accretionUpdate
    call sp_update_subdisk
    call sp_boundary_eject
    call sp_writeLog
  end subroutine sp_update
  subroutine sp_createGather(list_newPos)
    use mpilib
    type(t_spPos),pointer :: list_newPos, ptra, ptrb, prev
    type t_spw 
       integer :: w 
       type(t_spw),pointer :: next => null()
    end type t_spw
    type(t_spw),pointer :: weight, list_weight, ptrwa, ptrwb, prevw
    integer :: pida, pidb 
    integer,dimension(2) :: pair
    real(kind=8) :: dmin2, dist2
    if (.not. associated(list_newPos)) return
    nullify(list_weight)
    ptra => list_newPos
    do
       if (.not. associated( ptra ) ) exit
       allocate(weight)
       weight%w = 1
       if (associated(list_weight)) then
          ptrwa%next => weight
          ptrwa => ptrwa%next
       else
          list_weight => weight 
          ptrwa => list_weight
       endif
       ptra => ptra%next
    end do
    do
       pair(:) = Undefi
       dmin2 = huge(dmin2)
       ptra => list_newPos
       pida = 0
       do 
          if (.not. associated( ptra ) ) exit
          ptrb => list_newPos
          pidb = 0
          do 
             if (.not. associated( ptrb ) ) exit
             if (pidb >= pida) exit 
             dist2 = sum((ptra%r - ptrb%r)**2)
             if (dist2 < dmin2) then
                dmin2 = dist2
                pair(:) = (/ min(pida, pidb), max(pida, pidb) /)
             endif
             pidb = pidb + 1
             ptrb => ptrb%next
          end do
          pida = pida + 1
          ptra => ptra%next
       end do
       if (dmin2 > sp_MaskRadius**2) then
          call sp_createGather_dealloc_listweight(list_weight)
          return
       end if
       nullify(prev)
       nullify(prevw)
       ptra => list_newPos
       pida = 0
       ptrwa => list_weight
       do 
          if (.not. associated( ptra ) ) exit
          if (pida == pair(1)) then
             ptrb => list_newPos
             pidb = 0
             ptrwb => list_weight
             do 
                if (.not. associated( ptrb ) ) exit
                if (pidb == pair(2)) then
                   ptra%r = (ptra%r * ptrwa%w + ptrb%r * ptrwb%w)/(ptrwa%w + ptrwb%w) 
                   ptrwa%w = ptrwa%w + ptrwb%w
                   if (get_myrank() == 0) print "('merge particles: pos, weight ',3E12.5, I2)", ptra%r, ptrwa%w
                   prev%next => ptrb%next
                   prevw%next => ptrwb%next
                   deallocate( ptrb )
                   deallocate( ptrwb )
                   exit
                end if
                prev => ptrb
                pidb = pidb + 1
                ptrb => ptrb%next
                prevw => ptrwb
                ptrwb => ptrwb%next
             end do
          end if
          pida = pida + 1
          ptra => ptra%next
          ptrwa => ptrwa%next
       end do
    end do
    call sp_createGather_dealloc_listweight(list_weight)
    contains
      subroutine sp_createGather_dealloc_listweight(list_weight)
        type(t_spw),pointer :: list_weight, next, ptr
        ptr => list_weight
        do
           if (.not. associated(ptr)) exit
           next => ptr%next
           deallocate( ptr )
           ptr => next
        end do
      end subroutine sp_createGather_dealloc_listweight
  end subroutine sp_createGather
  subroutine sp_makeParticle(list_newPos)
    type(t_spPos),pointer :: list_newPos, ptr
    ptr => list_newPos
    do
       if (.not. associated( ptr ) ) exit
       call sp_newParticle( 0.d0, ptr%r, (/0.d0, 0.d0, 0.d0/), (/0.d0, 0.d0, 0.d0/), t_prev=Time(sp_Level) &
         , t_crt=Time(sp_Level) ) 
       ptr => ptr%next
    end do
  end subroutine sp_makeParticle
  subroutine sp_gravityP2P(ptri, grav)
    type(t_spParticle),pointer :: ptri
    real(kind=8),dimension(0:2),intent(OUT) :: grav
    type(t_spParticle),pointer :: ptr
    real(kind=8) :: gx, gy, gz
    real(kind=8),dimension(0:2) :: dr
    ptr => Particle
    grav = 0
    do
       if (.not. associated(ptr)) exit
       if (ptr%pid == ptri%pid) then
          ptr => ptr%next
          cycle
       end if
       dr = ptri%r - ptr%r
       call sp_gravityOfParticle( dr, ptr%mass, gx, gy, gz)
       grav(0) = grav(0) + gx
       grav(1) = grav(1) + gy
       grav(2) = grav(2) + gz
       ptr => ptr%next
    end do
  end subroutine sp_gravityP2P
  subroutine sp_gravityP2P_chi(grav)
    type(t_spParticle),pointer :: ptri
    real(kind=8),dimension(0:, :),intent(OUT) :: grav
    real(kind=8),dimension(:,:),allocatable :: grav_r
    type(t_spParticle),pointer :: ptr1, ptr2
    real(kind=8) :: gx, gy, gz
    real(kind=8),dimension(0:2) :: dr
    integer :: np, next_np
    allocate(grav_r(lbound(grav,1):ubound(grav,1), lbound(grav,2):ubound(grav,2)))
    grav = 0.d0
    grav_r = 0.d0
    np = 0
    ptr1 => Particle
    next_np = get_myrank() + 1
    do
      if (.not. associated(ptr1)) exit
      np = np + 1
      if(np == next_np) then
        ptr2 => Particle
        do
           if (.not. associated(ptr2)) exit
           if (ptr2%pid == ptr1%pid) then
              ptr2 => ptr2%next
              cycle
           end if
           dr = ptr1%r - ptr2%r
           call sp_gravityOfParticle( dr, ptr2%mass, gx, gy, gz)
           grav_r(0,np) = grav_r(0,np) + gx
           grav_r(1,np) = grav_r(1,np) + gy
           grav_r(2,np) = grav_r(2,np) + gz
           ptr2 => ptr2%next
        end do
        next_np = next_np + 400
      endif
      ptr1 => ptr1%next
    enddo
    call mpi_allreduce( grav_r, grav, size(grav_r), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    deallocate(grav_r)
  end subroutine sp_gravityP2P_chi
  subroutine sp_accretion_init
    type(t_spParticle),pointer :: ptr
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ptr%dmass = 0.d0
       ptr%dp(:) = 0.d0
       ptr%ds(:) = 0.d0
       ptr => ptr%next
    end do
  end subroutine sp_accretion_init
  subroutine sp_accretion_comm
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(:,:),allocatable :: buf, bufr
    integer :: n
    real(kind=8) :: r_diskout, v_diskr, t_vis, real_t, t_free, v_free
    real(kind=8) :: dmass_sink, dmass_bh, dmass_disk
    logical :: isNotFinite
    allocate(buf(7, Nparticle), bufr(7, Nparticle))
    ptr => Particle
    do n = 1, Nparticle
       if (.not. associated(ptr)) then
          print *, '*** error in sp_accretion: ptr is no associated.'
       end if
       buf(1, n) = ptr%dmass
       buf(2:4, n) = ptr%dp
       buf(5:7, n) = ptr%ds
       ptr => ptr%next
    end do
    call mpi_allreduce( buf, bufr , size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    ptr => Particle
    do n = 1, Nparticle
       if (.not. associated(ptr)) then
          print *, '*** error in sp_accretion: ptr is no associated.'
       end if
       ptr%dmass = buf(1, n)
       ptr%dp = buf(2:4, n)
       ptr%ds = buf(5:7, n)
       ptr => ptr%next
    end do
    deallocate(buf, bufr)
  end subroutine sp_accretion_comm
  subroutine sp_accretion
    use mpilib
    type(t_spParticle),pointer :: ptr
    integer :: level, gid, rank, i,j,k, ig, jg, kg
    integer,dimension(0:2) :: ijkgL, ijkgR
    real(kind=8),dimension(0:2) :: posL, posR
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:),pointer :: rho, vx, vy, vz, pr 
    real(kind=8) :: dv, rhoacc
    real(kind=8),dimension(0:2) :: dr, du
    ptr => Particle
    if (LevelMax < sp_Level) return
    level = sp_Level
    do
       if (.not. associated(ptr)) exit
       if (ptr%lev == sp_Level) then
          posL(:) = ptr%r(:) - sp_SinkRadius
          posR(:) = ptr%r(:) + sp_SinkRadius
          call sp_restrict_within_boundary(posL, posR)
          call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
          call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
          do kg = ijkgL(2), ijkgR(2)
            do jg = ijkgL(1), ijkgR(1)
              do ig = ijkgL(0), ijkgR(0)
                call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
                if (gid == Undefi) cycle
                if (rank /= get_myrank() ) cycle
                x => get_Xp(gid)
                y => get_Yp(gid)
                z => get_Zp(gid)
                rho => get_Ucomp(0, gid)
                pr => get_Ucomp(4, gid) 
                vx => get_Ucomp(1,gid)
                vy => get_Ucomp(2,gid)
                vz => get_Ucomp(3,gid)
                dv = get_dv(level)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         if ( (x(i)-ptr%r(0))**2 + (y(j)-ptr%r(1))**2 + (z(k)-ptr%r(2))**2 < sp_SinkRadius**2 ) then
                            rhoacc = max(rho(i,j,k) - sp_Rhocr, 0.d0) 
                            ptr%dmass = ptr%dmass + rhoacc*dv
                            ptr%dp(0) = ptr%dp(0) + vx(i,j,k)*rhoacc*dv
                            ptr%dp(1) = ptr%dp(1) + vy(i,j,k)*rhoacc*dv
                            ptr%dp(2) = ptr%dp(2) + vz(i,j,k)*rhoacc*dv
                            dr = (/x(i), y(j), z(k)/) - ptr%r
                            du = (/vx(i,j,k), vy(i,j,k), vz(i,j,k)/) - ptr%v
                            ptr%ds(0) = ptr%ds(0) + (dr(1)*du(2)-dr(2)*du(1))*rhoacc*dv
                            ptr%ds(1) = ptr%ds(1) + (dr(2)*du(0)-dr(0)*du(2))*rhoacc*dv
                            ptr%ds(2) = ptr%ds(2) + (dr(0)*du(1)-dr(1)*du(0))*rhoacc*dv
                            pr(i,j,k) = pr(i,j,k) * (rho(i,j,k)-rhoacc)/rho(i,j,k) 
                            rho(i,j,k) = rho(i,j,k) - rhoacc
                         end if
                      end do 
                   end do 
                end do 
              end do 
            end do 
         end do 
       endif 
       ptr => ptr%next
    end do
  end subroutine sp_accretion
  subroutine sp_newParticle(mass, r, v, s, pid, dmass, dr, dv, dp, ds, &
    t_prev, dm_disk, dJ_disk, mdot_disk, J_disk, t_crt)
    real(kind=8),intent(IN) :: mass
    real(kind=8),dimension(0:2),intent(IN) :: r, v, s
    real(kind=8),intent(IN),optional :: dmass
    real(kind=8),dimension(0:2),intent(IN),optional :: dr, dv, dp, ds
    integer,intent(IN),optional :: pid
    real(kind=8),intent(IN),optional :: t_prev, dm_disk, mdot_disk, t_crt
    real(kind=8),dimension(0:2),intent(IN),optional :: dJ_disk, J_disk
    type(t_spParticle),pointer :: newParticle
    allocate(newParticle)
    newParticle%r = r
    newParticle%v = v
    newParticle%s = s
    newParticle%mass = mass
    if (present(dmass)) then
       newParticle%dmass = dmass
    else
       newParticle%dmass = 0.d0
    end if
    if (present(dp)) then
       newParticle%dp = dp
    else
       newParticle%dp = 0.d0
    end if
    if (present(dr)) then
       newParticle%dr = dr
    else
       newParticle%dr = 0.d0
    end if
    if (present(dv)) then
       newParticle%dv = dv
    else
       newParticle%dv = 0.d0
    end if
    if (present(ds)) then
       newParticle%ds = ds
    else
       newParticle%ds = 0.d0
    end if
    if (present(t_prev)) then
       newParticle%t_prev = t_prev
    else
       newParticle%t_prev = 0.d0
    end if
    if (present(dm_disk)) then
       newParticle%dm_disk = dm_disk
    else
       newParticle%dm_disk = 0.d0
    end if
    if (present(dJ_disk)) then
       newParticle%dJ_disk = dJ_disk
    else
       newParticle%dJ_disk = 0.d0
    end if
    if (present(mdot_disk)) then
       newParticle%mdot_disk = mdot_disk
    else
       newParticle%mdot_disk = 0.d0
    end if
    if (present(t_crt)) then
       newParticle%t_crt = t_crt
    else
       newParticle%t_crt = Time(LevelMax)
    end if
    if (present(J_disk)) then
       newParticle%J_disk = J_disk
    else
       newParticle%J_disk = 0.d0
    end if
    if (present(pid)) then
       if (bool_checkPidDuplicated (pid)) &
            print *, '*** error in sp_newParticle: pid is duplicated.', pid
       newParticle%pid = pid
    else
       newParticle%pid = PidMax + 1
    end if
    PidMax = max(PidMax, newParticle%pid)
    Nparticle = Nparticle + 1
    call sp_appendParticleList(Particle, newParticle)
  contains
    function bool_checkPidDuplicated(pid) result(bool)
      integer,intent(IN) :: pid
      logical :: bool
      type(t_spParticle),pointer :: ptr
      bool = .false.
      ptr => Particle
      do
         if (.not. associated(ptr)) exit
         if ( pid == ptr%pid ) then
            bool = .true.
            exit
         end if
         ptr => ptr%next
      end do
    end function bool_checkPidDuplicated
  end subroutine sp_newParticle
  subroutine sp_deleteParticle(pid)
    integer,intent(IN) :: pid
    type(t_spParticle),pointer :: ptr, ptrdel, next
    if (.not. associated(Particle)) then
       print *, '** error in sp_deleteParticle: No particle is assigned', pid
       call flush(6)
       return
    end if
    Nparticle = Nparticle - 1
    ptr => Particle
    next => ptr%next
    if (Particle%pid == pid) then
       deallocate( Particle )
       Particle => next
       return
    end if
    do
       if (.not. associated(ptr%next)) then
          print *, '*** No particle is deleted. Any pid is not match.'
          call flush(6)
          return
       endif
       if (ptr%next%pid == pid) then
          ptrdel => ptr%next
          ptr%next => ptr%next%next
          deallocate( ptrdel )
          return
       end if
       ptr => ptr%next
    end do
  end subroutine sp_deleteParticle
  subroutine sp_deleteAllParticles
    type(t_spParticle),pointer :: ptr, next
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       next => ptr%next
       call sp_deleteParticle(ptr%pid)
       ptr => next
    end do
  end subroutine sp_deleteAllParticles
  subroutine sp_appendParticleList(a, b)
    type(t_spParticle),pointer :: a, b
    type(t_spParticle),pointer :: ptr
    if ((.not. associated(a)) .and. (.not. associated(b))) then
       print *, '*** error in sp_appendParticleList: Both pointers are not associated.'
       return
    end if
    if (.not. associated(a)) then
       a => b
       return
    end if
    if (.not. associated(b)) then
       b => a
       return
    end if
    ptr => a
    do
       if (.not. associated(ptr%next)) then
          ptr%next => b
          exit
       end if
       ptr => ptr%next
    end do
  end subroutine sp_appendParticleList
  subroutine sp_boundary_eject
    type(t_spParticle),pointer :: ptr, next
    type(t_obRectPhys) :: compDomain, particleDomain
    type(t_obPointPhys) :: pL, pR
    real(kind=8),dimension(0:2) :: posL, posR
    call ob_computationBoxOfRectPhys( compDomain )
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       posL(:) = ptr%r(:) - sp_SinkRadius
       posR(:) = ptr%r(:) + sp_SinkRadius
       call ob_assignCoordPhysToPointPhys(posL, pL)
       call ob_assignCoordPhysToPointPhys(posR, pR)
       call ob_assignPointPhysToRectPhys(pL, particleDomain, 'L')
       call ob_assignPointPhysToRectPhys(pR, particleDomain, 'R')
       next => ptr%next
       if ( ob_rectPhysComp( particleDomain, compDomain ) /= OB_RECT_INNER ) then
          call sp_boundary_eject_write(ptr)
          call sp_deleteParticle(ptr%pid)
       end if
       ptr => next
    end do
  contains
    subroutine sp_boundary_eject_write(ptr)
      use mpilib
      use string, only : concat
      use io_util, only : read_env
      type(t_spParticle),pointer :: ptr
      character(len=CHARLEN),parameter :: FILENAME_EJECT='logSinkparticle_ejected'
      character(len=CHARLEN) :: file, dir
      integer,parameter :: LUN = 11
      if (get_myrank() /= 0) return
      if (.not. associated(Particle)) return
      call read_env('DIR', dir)
      file = concat(dir,FILENAME_EJECT)
      open(LUN, file=file, position='APPEND')
      write(LUN, '(I14, E17.9, I7, 10(1PE17.9))') Step(LevelMax), Time(LevelMax), ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%s
      call flush(LUN)
      close(LUN)
    end subroutine sp_boundary_eject_write
  end subroutine sp_boundary_eject
  subroutine sp_merge
    use parameter
    use eos
    use modelParameter, only : MP_CONNECTION_RUN 
    type(t_spParticle),pointer :: ptri, ptrj
    real(kind=8) :: distanceMerge
    integer :: piddel
    real(kind=8) :: rcm(0:2), vcm(0:2) 
    real(kind=8) :: dri(0:2), drj(0:2), dui(0:2), duj(0:2)
    distanceMerge = sp_SinkRadius * 2.d0
    ptri => Particle
    do
       if (.not. associated(ptri)) exit
       ptrj => Particle
       do
          if (.not. associated(ptrj)) exit
          if ( ptri%pid == ptrj%pid ) then
             ptrj => ptrj%next
             cycle
          endif
          if ( sum((ptri%r - ptrj%r)**2) > distanceMerge**2 ) then
             ptrj => ptrj%next
             cycle
          end if
          if (ptri%mass < TINY(ptri%mass) .and. ptrj%mass < TINY(ptrj%mass)) then
             ptrj => ptrj%next
             cycle
          end if
          if (get_myrank() == 0) then
             print *, '*** particles are merged ', ptri%pid, ptrj%pid
             call flush(6)
          end if
          if (MP_CONNECTION_RUN == 0) then
             rcm = (ptri%r*ptri%mass + ptrj%r*ptrj%mass)/(ptri%mass + ptrj%mass)
             vcm = (ptri%v*ptri%mass + ptrj%v*ptrj%mass)/(ptri%mass + ptrj%mass)
             dri(:) = ptri%r(:) - rcm(:)
             drj(:) = ptrj%r(:) - rcm(:)
             dui(:) = ptri%v(:) - vcm(:)
             duj(:) = ptrj%v(:) - vcm(:)
             ptri%dm_disk = ptri%dm_disk + ptrj%mass
             ptri%dJ_disk(0) = ptri%dJ_disk(0) + &
                  (dri(1)*dui(2)-dri(2)*dui(1))*ptri%mass + (drj(1)*duj(2)-drj(2)*duj(1))*ptrj%mass
             ptri%dJ_disk(1) = ptri%dJ_disk(1) + &
                  (dri(2)*dui(0)-dri(0)*dui(2))*ptri%mass + (drj(2)*duj(0)-drj(0)*duj(2))*ptrj%mass
             ptri%dJ_disk(2) = ptri%dJ_disk(2) + &
                  (dri(0)*dui(1)-dri(1)*dui(0))*ptri%mass + (drj(0)*duj(1)-drj(1)*duj(0))*ptrj%mass
          end if
          ptri%mdot_disk = ptri%mdot_disk + ptrj%mdot_disk
          ptri%J_disk = ptri%J_disk + ptrj%J_disk
          ptri%r = (ptri%r*ptri%mass + ptrj%r*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%v = (ptri%v*ptri%mass + ptrj%v*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%s = (ptri%s*ptri%mass + ptrj%s*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dr = (ptri%dr*ptri%mass + ptrj%dr*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dv = (ptri%dv*ptri%mass + ptrj%dv*ptrj%mass)/(ptri%mass + ptrj%mass)
          ptri%dp = ptri%dp + ptrj%dp
          ptri%ds = ptri%ds + ptrj%ds
          ptri%mass = ptri%mass + ptrj%mass
          ptri%dmass = ptri%dmass + ptrj%dmass
          ptri%t_crt = min(ptri%t_crt, ptrj%t_crt)
          piddel = ptrj%pid
          ptrj => ptrj%next
          call sp_deleteParticle(piddel)
       end do
       ptri => ptri%next
    end do
  end subroutine sp_merge
  subroutine sp_write
    use systemcall
    use mpilib
    use string, only : concat, num2char
    use io_util, only : read_env
    integer,parameter :: LUN = 11
    character(len=CHARLEN) :: file, dir, filen
    character(len=2),parameter :: PREFIX = 'sp', SUFFIX = '.d'
    type(t_spParticle),pointer :: ptr
    logical :: exist
    call sp_init
    if (Nparticle == 0) return
    if (get_myrank() /= 0) return
    call read_env('DIR', dir)
    file = concat(dir,FILENAME) 
    filen = trim(dir)//PREFIX//trim(num2char(Step(LevelMax)))//SUFFIX 
    open(LUN, file=filen, form='unformatted')
    write(LUN) Nparticle, PidMax
    write(LUN) sp_DtimeCFL
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       write(LUN) ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds, &
            ptr%t_prev, ptr%dm_disk, ptr%dJ_disk, ptr%mdot_disk, ptr%J_disk, ptr%t_crt
       ptr => ptr%next
    end do
    call flush(LUN)
    close(LUN)
    inquire(file=filen, exist=exist)
    if (.not. exist) print *, '*** ERROR in sp_write'
    inquire(file=file, exist=exist)
    if (exist) call systemcall_unlink(file)
    call systemcall_symlink(filen, file)
  end subroutine sp_write
  subroutine sp_read
    use mpilib
    use string, only : concat
    use io_util, only : read_env
    integer,parameter :: LUN = 11
    character(len=CHARLEN) :: file, dir
    integer :: n, nptcle, pid
    real(kind=8),dimension(:,:),allocatable :: buf
    integer,dimension(:),allocatable :: bufi
    real(kind=8),dimension(0:2) :: r, v, dr, dv, dp, s, ds
    real(kind=8) :: mass, dmass
    logical :: exist
    real(kind=8) :: t_prev, dm_disk, mdot_disk
    real(kind=8),dimension(0:2) :: dJ_disk, J_disk
    real(kind=8) :: t_crt
    call sp_init
    if (get_myrank() == 0 ) then
       call read_env('DIR', dir)
       file = concat(dir,FILENAME)
       inquire(file=file, exist=exist)
    end if
    call mpi_bcast(exist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (.not. exist) return
    if (get_myrank() == 0 ) then
       open(LUN, file=file, form='unformatted')
       read(LUN) nptcle, PidMax
       read(LUN) sp_DtimeCFL
    end if
    call mpi_bcast(nptcle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(PidMax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sp_DtimeCFL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (nptcle == 0) return
    allocate( bufi(nptcle) )
    allocate( buf(1+1+3*7+1*3+3*2+1, nptcle ) ) 
    if (get_myrank() == 0) then
       do n = 1, nptcle
          read(LUN) bufi(n), buf(:, n)
       end do
       close(LUN)
    end if
    call mpi_bcast(bufi, size(bufi), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(buf, size(buf), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    do n = 1, nptcle
       pid = bufi(n)
       mass = buf(1,n)
       r = buf(2:2+2,n)
       v = buf(5:5+2,n)
       dmass = buf(8,n)
       dr = buf(9:9+2,n)
       dv = buf(12:12+2,n)
       dp = buf(15:15+2,n)
       s = buf(18:18+2,n)
       ds = buf(21:21+2,n)
       t_prev = buf(24,n)
       dm_disk = buf(25,n)
       dJ_disk = buf(26:26+2,n)
       mdot_disk = buf(29,n)
       J_disk = buf(30:30+2,n)
       t_crt = buf(33,n)
       call sp_newparticle( mass, r, v, s, pid=pid, dmass=dmass, dr=dr, dv=dv, dp=dp, ds=ds, &
            t_prev=t_prev, dm_disk=dm_disk, dj_disk=dj_disk, mdot_disk=mdot_disk, &
            j_disk=j_disk, t_crt=t_crt)
    end do
    if (nptcle /= Nparticle) print *, '*** error in sp_read: Nparticle is not consistent', Nparticle, nptcle
    deallocate( bufi, buf )
  end subroutine sp_read
  subroutine sp_writeLog
    use mpilib
    use string, only : concat, num2char
    use io_util, only : read_env
    integer,parameter :: LUN = 11
    character(len=CHARLEN) :: file, dir
    character(len=CHARLEN) :: file_each 
    character(len=CHARLEN),parameter :: prefix_log='spLog' 
    character(len=CHARLEN),parameter :: prefix_log2='logsp' 
    integer,parameter :: LUN2 = 12 
    type(t_spParticle),pointer :: ptr
    integer,parameter :: skip = 10
    integer :: slevel
    slevel = LevelMax
    if (mod(Step(slevel), skip) /= 0) return
    call sp_init
    if (Nparticle == 0) return
    if (get_myrank() /= 0) return
    if (.not. associated(Particle)) return
    call read_env('DIR', dir)
    file = concat(dir,FILENAME_LOG)
    open(LUN, file=file, position='APPEND')
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       write(LUN, '(I14, 2(1PE17.9), I7, 24(1PE23.15))') Step(slevel), Time(slevel), Dtime(slevel), &
            ptr%pid, ptr%mass, ptr%r, ptr%v, ptr%dmass, ptr%dr, ptr%dv, ptr%dp, ptr%s, ptr%ds,Time(slevel)-ptr%t_crt
       file_each = concat(dir,prefix_log2) 
       file_each = concat(file_each,num2char(ptr%pid)) 
       file_each = concat(file_each,'.dat') 
       open(LUN2, file=file_each, position='APPEND')
       write(LUN2, '(I10, 18(1PE15.7))') Step(slevel), Time(slevel)*Unit_yr, &
               ptr%mass*Unit_msun, ptr%dmass/Dtime(slevel)*Unit_msun/Unit_yr,&
               ptr%r*Unit_au, ptr%v*Unit_kms, ptr%s, &
               ptr%t_prev*Unit_yr, ptr%mdot_disk*Unit_msun/Unit_yr, &
               ptr%J_disk*Unit_msun*Unit_kms*Unit_au,(Time(slevel)-ptr%t_crt)*Unit_yr 
       call flush(LUN2)
       close(LUN2)
       ptr => ptr%next
    end do
    call flush(LUN)
    close(LUN)
  end subroutine sp_writeLog
  subroutine sp_gravityOfParticle(dr, mass, gx, gy, gz, psi, no_softening)
    use modelParameter, only : MP_Gconst
    real(kind=8),dimension(0:2),intent(IN) :: dr
    real(kind=8),intent(IN) :: mass
    real(kind=8),intent(OUT) :: gx, gy, gz
    real(kind=8),intent(OUT),optional :: psi
    logical,optional :: no_softening
    real(kind=8) :: flagIn, flagOut, gabs, gm, sradii, sradii2, sradii3, drabs
    logical :: no_soft
    if (present(no_softening)) then
       no_soft = no_softening
    else
       no_soft = .false.
    endif
    sradii = sp_SofteningRadius
    sradii2 = sradii**2
    sradii3 = sradii**3
    gm = MP_Gconst * mass
    drabs = sqrt(sum(dr**2))
    if (no_soft) then
       flagOut = 1
    else
       flagOut = min(int( drabs /sradii), 1) 
    endif
    flagIn = 1 - flagOut 
    gabs = - gm * ( flagIn / sradii3 + flagOut / (drabs**3 + flagIn) )
    gx = gabs * dr(0)
    gy = gabs * dr(1)
    gz = gabs * dr(2)
    if (present(psi)) then
       psi = - gm * ( &
            flagIn * (1.5d0 - 0.5d0*drabs**2/sradii2) / sradii + &
            flagOut / (drabs + flagIn) &
         )
    endif
  end subroutine sp_gravityOfParticle
  subroutine sp_gravityOfParticle2fluid(dr, mass, gx, gy, gz, soft_length, psi, no_softening)
    use modelParameter, only : MP_Gconst
    real(kind=8),dimension(0:2),intent(IN) :: dr
    real(kind=8),intent(IN) :: mass, soft_length
    real(kind=8),intent(OUT) :: gx, gy, gz
    real(kind=8),intent(OUT),optional :: psi
    logical,optional :: no_softening
    real(kind=8) :: flagIn, flagOut, gabs, gm, sradii, sradii2, sradii3, drabs
    logical :: no_soft
    if (present(no_softening)) then
       no_soft = no_softening
    else
       no_soft = .false.
    endif
    sradii = soft_length
    sradii2 = sradii**2
    sradii3 = sradii**3
    gm = MP_Gconst * mass
    drabs = sqrt(sum(dr**2))
    if (no_soft) then
       flagOut = 1
    else
       flagOut = min(int( drabs /sradii), 1) 
    endif
    flagIn = 1 - flagOut 
    gabs = - gm * ( flagIn / sradii3 + flagOut / (drabs**3 + flagIn) )
    gx = gabs * dr(0)
    gy = gabs * dr(1)
    gz = gabs * dr(2)
    if (present(psi)) then
       psi = - gm * ( &
            flagIn * (1.5d0 - 0.5d0*drabs**2/sradii2) / sradii + &
            flagOut / (drabs + flagIn) &
         )
    endif
  end subroutine sp_gravityOfParticle2fluid
  subroutine sp_gravityOfParticle_out_chi(dr, mass, gx, gy, gz)
    use modelParameter, only : MP_Gconst
    real(kind=8),dimension(0:2),intent(IN) :: dr
    real(kind=8),intent(IN) :: mass
    real(kind=8),intent(OUT) :: gx, gy, gz
    real(kind=8) :: gabs, gm, drabs2
    gm = MP_Gconst * mass
    drabs2 = sum(dr**2)
    gabs = -gm/(drabs2*sqrt(drabs2))
    gx = gabs * dr(0)
    gy = gabs * dr(1)
    gz = gabs * dr(2)
  end subroutine sp_gravityOfParticle_out_chi
  subroutine sp_gravityFluid2particle(gravF, dt)
    use mpilib
    real(kind=8),dimension(0:,:),intent(OUT) :: gravF
    real(kind=8),intent(IN) :: dt
    real(kind=8),dimension(:,:), allocatable :: gravBlock
    integer :: level, n, gid
    allocate(gravBlock(0:2,Nparticle))
    gravF = 0.d0
    myrank = get_myrank()
    do level = Lmin, LevelMax
       do n = Gidmin, GidListMax( level )
          gid = GidList(n, level)
          call sp_gravityFluid2particleForBlck(gravBlock, dt, gid)
          if (ChildGid(Left, Left, Left, gid, myrank) == Undefi) then
             gravF = gravF + gravBlock
          end if
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, gravF, size(gravF), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    deallocate(gravBlock)
  end subroutine sp_gravityFluid2particle
  subroutine sp_gravityFluid2particleForBlck(grav, dt, gid)
    use eos, only : w2u_4, u2w_4
    use mpilib 
    use modelParameter, only : MP_Gconst
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(0:, :),intent(OUT) :: grav
    real(kind=8),intent(IN) :: dt
    integer,intent(IN) :: gid
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:,:),pointer :: u, w
    integer :: level, i, j, k, ic, jc, kc, ilev, ilgid, rank, levsp
    real(kind=8) :: dv, dvsc, xc, yc, zc, dtdvrho, gxp, gyp, gzp, nsubm1, dtdv, dtdvsc
    real(kind=8),dimension(0:2) :: h, hh, hsc, bl, br, pl, pr, cl, cr, dr
    real(kind=8) :: gm, drabs2, gabs, dryz_2, drz_2, inv_dt
    real(kind=8),dimension(0:2) :: pos
    real(kind=8) :: soft_length
    integer,parameter :: N_SubCell = 8
    integer :: np
    logical :: isNotFinite 
    level = get_level(gid)
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    u => get_Up(gid)
    dv = get_dv(level)
    allocate( w(lbound(u,1):ubound(u,1),lbound(u,2):ubound(u,2),lbound(u,3):ubound(u,3),0:4) )
    call u2w_4(u, w, dv)
    h(:) = CellWidth(:,level)
    hh = 0.5d0*h
    bl = (/x(Imin), y(Jmin), z(Kmin)/) - hh 
    br = (/x(Imax), y(Jmax), z(Kmax)/) + hh 
    hsc(:) = CellWidth(:, level) / N_SubCell 
    dvsc = dv/N_SubCell**3
    grav = 0.d0
    dtdv = dt*dv
    dtdvsc = dt*dvsc
    inv_dt = 1.d0/dt
    nsubm1 = (N_SubCell-1)*0.5d0
    np = 0
    ptr => Particle
    do
      if (.not. associated(ptr)) exit
      np = np + 1
      gm = MP_Gconst* ptr%mass
      soft_length = sp_SofteningRadius*2.d0**(sp_Level-ptr%lev)
      pl = ptr%r - soft_length
      pr = ptr%r + soft_length
      if ( max(pl(0),bl(0)) < min(pr(0),br(0)) .and. &
           max(pl(1),bl(1)) < min(pr(1),br(1)) .and. &
           max(pl(2),bl(2)) < min(pr(2),br(2)) ) then 
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  cl = (/x(i), y(j), z(k)/) - hh
                  cr = (/x(i), y(j), z(k)/) + hh
                  if ( max(pl(0),cl(0)) < min(pr(0),cr(0)) .and. &
                       max(pl(1),cl(1)) < min(pr(1),cr(1)) .and. &
                       max(pl(2),cl(2)) < min(pr(2),cr(2)) ) then 
                     dtdvrho = dtdvsc*u(i,j,k,0)
                     do kc = 0, N_SubCell-1
                        do jc = 0, N_SubCell-1
                           do ic = 0, N_SubCell-1
                              xc = x(i) + ( ic - nsubm1 ) * hsc(0)
                              yc = y(j) + ( jc - nsubm1 ) * hsc(1)
                              zc = z(k) + ( kc - nsubm1 ) * hsc(2)
                              dr = (/ xc-ptr%r(0), yc-ptr%r(1), zc-ptr%r(2) /)
                              call sp_gravityOfParticle2fluid(dr, ptr%mass, gxp, gyp, gzp, soft_length)
                              w(i,j,k,1) = w(i,j,k,1) + gxp*dtdvrho
                              w(i,j,k,2) = w(i,j,k,2) + gyp*dtdvrho
                              w(i,j,k,3) = w(i,j,k,3) + gzp*dtdvrho
                              w(i,j,k,4) = w(i,j,k,4) + dtdvrho*( &
                                    u(i,j,k,1)*gxp &
                                   +u(i,j,k,2)*gyp &
                                   +u(i,j,k,3)*gzp )
                              grav(:,np) = grav(:,np) - (/gxp, gyp, gzp/) * dvsc * u(i,j,k,0)
                           end do
                        end do
                     end do
                  else 
                     dtdvrho = dtdv*u(i,j,k,0)
                     dr = (/ x(i)-ptr%r(0), y(j)-ptr%r(1), z(k)-ptr%r(2) /)
                     call sp_gravityOfParticle_out_chi(dr, ptr%mass, gxp, gyp, gzp)
                     w(i,j,k,1) = w(i,j,k,1) + gxp*dtdvrho
                     w(i,j,k,2) = w(i,j,k,2) + gyp*dtdvrho
                     w(i,j,k,3) = w(i,j,k,3) + gzp*dtdvrho
                     w(i,j,k,4) = w(i,j,k,4) + dtdvrho*( &
                           u(i,j,k,1)*gxp &
                          +u(i,j,k,2)*gyp &
                          +u(i,j,k,3)*gzp )
                     grav(:,np)=grav(:,np)-(/gxp, gyp, gzp/)*dv*u(i,j,k,0)
                  end if
               end do
            end do
         end do
      else 
         do k = Kmin, Kmax
            do j = Jmin, Jmax
               do i = Imin, Imax
                  dtdvrho = dtdv*u(i,j,k,0)
                  dr = (/ x(i)-ptr%r(0), y(j)-ptr%r(1), z(k)-ptr%r(2) /)
                  call sp_gravityOfParticle_out_chi(dr, ptr%mass, gxp, gyp, gzp)
                     w(i,j,k,1) = w(i,j,k,1) + gxp*dtdvrho
                     w(i,j,k,2) = w(i,j,k,2) + gyp*dtdvrho
                     w(i,j,k,3) = w(i,j,k,3) + gzp*dtdvrho
                     w(i,j,k,4) = w(i,j,k,4) + dtdvrho*( &
                           u(i,j,k,1)*gxp &
                          +u(i,j,k,2)*gyp &
                          +u(i,j,k,3)*gzp )
                     grav(:,np) = grav(:,np) - (/gxp, gyp, gzp/) * dv * u(i,j,k,0)
               end do
            end do
         end do
      end if
      if (ptr%mass > TINY(ptr%mass)) then
        grav(:,np) = grav(:,np) / ptr%mass
      else
        grav(:,np) = 0.d0
      end if
      ptr => ptr%next
    enddo 
    call w2u_4(w, u, dv)
    deallocate( w )
  end subroutine sp_gravityFluid2particleForBlck
  subroutine sp_gravityOfFluid(ptr, grav)
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(0:2),intent(OUT) :: grav
    real(kind=8) :: vol
    real(kind=8) :: vol_spsph 
    call sp_sumUinSphere_chi(ptr%r, sp_SofteningRadius, (/18, 19, 20/), grav, vol)
    grav = grav/vol
    vol_spsph = 4.*3.14/3.*sp_SofteningRadius**3 
    if (abs((vol-vol_spsph)/vol_spsph) > 0.1) then
       if (get_myrank() == 0) &
            print '(A,1P10E15.7)', "(sp_gravityOfFluid) ** WARNING ** sp volume is not consistent: ", &
            ptr%r, vol, vol_spsph
    end if
  end subroutine sp_gravityOfFluid
  subroutine sp_gravityOfFluid_chi(grav)
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(0:,:),intent(OUT) :: grav
    real(kind=8),dimension(0:2) :: grav_ptr
    real(kind=8) :: vol
    real(kind=8) :: vol_spsph 
    integer :: np
    np = 0
    ptr => Particle
    do
      if (.not. associated(ptr)) exit
      np = np + 1
      call sp_sumUinSphere_chi(ptr%r, sp_SofteningRadius, (/18, 19, 20/), grav_ptr, vol)
      grav(:,np) = grav_ptr(:)/vol
      vol_spsph = 4.*3.14/3.*sp_SofteningRadius**3 
      if (abs((vol-vol_spsph)/vol_spsph) > 0.1) then
         if (get_myrank() == 0) &
              print '(A,1P10E15.7)', "(sp_gravityOfFluid) ** WARNING ** sp volume is not consistent: ", &
              ptr%r, vol, vol_spsph
      end if
      ptr => ptr%next
    enddo
  end subroutine sp_gravityOfFluid_chi
  subroutine sp_sumUinSphere(point, radius, mlist, usum, vol)
    real(kind=8),dimension(0:2),intent(IN) :: point
    real(kind=8),intent(IN) :: radius
    integer,intent(IN),dimension(:) :: mlist
    real(kind=8),dimension(size(mlist)),intent(OUT) :: usum
    real(kind=8),intent(OUT) :: vol
    integer,parameter :: N_SubCell = 8
    integer :: level, gid, rank, i,j,k, m, ic, jc, kc, ig, jg, kg
    real(kind=8),dimension(0:2) :: posL, posR
    integer,dimension(0:2) :: ijkgL, ijkgR
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:),pointer :: u
    real(kind=8) :: dv, hsc(0:2), xc, yc, zc
    real(kind=8),dimension(size(mlist)+1) :: buf, bufr
    level = sp_Level
    dv = get_dv(level) / N_SubCell**3
    hsc(:) = CellWidth(:, level) / N_SubCell 
    usum(:) = 0.d0
    vol = 0.d0
    posL(:) = point(:) - radius
    posR(:) = point(:) + radius
    call sp_restrict_within_boundary(posL, posR)
    call ob_getIjkgridFromCoordPhys(ijkgL, level, posL)
    call ob_getIjkgridFromCoordPhys(ijkgR, level, posR)
    do kg = ijkgL(2), ijkgR(2)
       do jg = ijkgL(1), ijkgR(1)
          do ig = ijkgL(0), ijkgR(0)
             call get_gid_from_ijkgrid(ig,jg,kg,level,gid,rank)
             if (gid == Undefi) cycle
             if (rank /= get_myrank() ) cycle
             x => get_Xp(gid)
             y => get_Yp(gid)
             z => get_Zp(gid)
             do m = lbound(mlist, 1), ubound(mlist, 1)
                u => get_Ucomp(mlist(m), gid)
                do k = Kmin, Kmax
                   do j = Jmin, Jmax
                      do i = Imin, Imax
                         do kc = 0, N_SubCell-1
                            do jc = 0, N_SubCell-1
                               do ic = 0, N_SubCell-1
                                  xc = x(i) + ( ic - (N_SubCell-1)*0.5d0 ) * hsc(0)
                                  yc = y(j) + ( jc - (N_SubCell-1)*0.5d0 ) * hsc(1)
                                  zc = z(k) + ( kc - (N_SubCell-1)*0.5d0 ) * hsc(2)
                                  if ( (xc-point(0))**2 + (yc-point(1))**2 + (zc-point(2))**2 <= radius**2 ) then
                                     usum(m) = usum(m) + u(i,j,k)*dv
                                     if (m == 1) vol = vol + dv
                                  end if
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    buf(1) = vol
    buf(2:size(mlist)+1) = usum(:)
    call mpi_allreduce( buf, bufr, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    vol = buf(1)
    usum(:) = buf(2:size(mlist)+1)
  end subroutine sp_sumUinSphere
  subroutine sp_sumUinSphere_chi(point, radius, mlist, usum, vol)
    use modelParameter, only : MP_spRadius_cell
    real(kind=8),dimension(0:2),intent(IN) :: point
    real(kind=8),intent(IN) :: radius
    integer,intent(IN),dimension(:) :: mlist
    real(kind=8),dimension(size(mlist)),intent(OUT) :: usum
    real(kind=8),intent(OUT) :: vol
    integer,parameter :: N_SubCell = 8
    integer :: level, gid, rank, i,j,k, m, ic, jc, kc, ig, jg, kg
    integer :: N_spsub
    real(kind=8),dimension(0:2) :: posL, posR, posSP, posbt, bl, br, psb
    integer,dimension(0:2) :: ijkg, ijkgL, ijkgR
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(:,:,:),pointer :: u
    real(kind=8) :: dv, hsc(0:2), xc, yc, zc
    real(kind=8) :: h(0:2), hspace, radius2
    real(kind=8),dimension(size(mlist)+1) :: buf, bufr
    real(kind=8) :: ti, uj, vk, tim, ujm, vkm, ulin
    integer :: srank, slevel, sgid, splevel_max
    radius2 = radius*radius
    hsc(:) = CellWidth(:, Lmin) / (N_SubCell*2**(sp_Level-Lmin)) 
    dv = hsc(0)*hsc(1)*hsc(2)
    usum(:) = 0.d0
    vol = 0.d0
    posL(:) = point(:) - radius
    posR(:) = point(:) + radius
    call sp_restrict_within_boundary(posL, posR)
    N_spsub = 2*MP_spRadius_cell*N_Subcell
    posSP(0)=hsc(0)*dble(nint(point(0)/hsc(0)))
    posSP(1)=hsc(1)*dble(nint(point(1)/hsc(1)))
    posSP(2)=hsc(2)*dble(nint(point(2)/hsc(2)))
    posbt(0)=hsc(0)*(dble(nint(point(0)/hsc(0)))-dble(N_spsub)*0.5d0-0.5d0)
    posbt(1)=hsc(1)*(dble(nint(point(1)/hsc(1)))-dble(N_spsub)*0.5d0-0.5d0)
    posbt(2)=hsc(2)*(dble(nint(point(2)/hsc(2)))-dble(N_spsub)*0.5d0-0.5d0)
    do kc = 1, N_spsub
      zc = posbt(2)+kc*hsc(2)
      do jc = 1, N_spsub
        yc = posbt(1)+jc*hsc(1)
        do ic = 1, N_spsub
          xc = posbt(0)+ic*hsc(0)
          if ( (xc-posSP(0))**2 + (yc-posSP(1))**2 + (zc-posSP(2))**2 <= radius2 ) then
             psb = (/xc, yc, zc/)
             sgid = Undefi
             srank = Undefi
             splevel_max = 0
             do level = LevelMax, Lmin, -1 
                 call ob_getIjkgridFromCoordPhys(ijkg, level, psb)
                 call get_gid_from_ijkgrid(ijkg(0),ijkg(1),ijkg(2),level,gid,rank)
                 if (gid /= Undefi) then
                   sgid = gid
                   slevel = level
                   srank = rank
                   splevel_max = max(splevel_max, slevel)
                   exit 
                 endif
             enddo
             if (srank == get_myrank()) then
                 x => get_Xp(sgid)
                 y => get_Yp(sgid)
                 z => get_Zp(sgid)
                 h = CellWidth(:,slevel)
                 i = int((psb(0)-x(Imingh))/h(0))+Imingh
                 j = int((psb(1)-y(Jmingh))/h(1))+Jmingh
                 k = int((psb(2)-z(Kmingh))/h(2))+Kmingh
                 ti = (psb(0) - x(i))/(x(i+1) - x(i))
                 uj = (psb(1) - y(j))/(y(j+1) - y(j))
                 vk = (psb(2) - z(k))/(z(k+1) - z(k))
                 tim = 1.d0 - ti
                 ujm = 1.d0 - uj
                 vkm = 1.d0 - vk
                 do m = lbound(mlist, 1), ubound(mlist, 1)
                     u => get_Ucomp(mlist(m), sgid)
                     ulin = tim*ujm*vkm*u(i,j,k) + ti *ujm*vkm*u(i+1,j,k) &
                          + tim*uj *vkm*u(i,j+1,k) + ti *uj *vkm*u(i+1,j+1,k) &
                          + tim*ujm*vk *u(i,j,k+1) + ti *ujm*vk *u(i+1,j,k+1) &
                          + tim*uj *vk *u(i,j+1,k+1)+ ti *uj *vk *u(i+1,j+1,k+1)
                     usum(m) = usum(m) + ulin*dv
                     if (m == 1) vol = vol + dv
                 enddo
             endif
          endif
        enddo
      enddo
    enddo
    buf(1) = vol
    buf(2:size(mlist)+1) = usum(:)
    call mpi_allreduce( buf, bufr, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    buf = bufr
    vol = buf(1)
    usum(:) = buf(2:size(mlist)+1)
  end subroutine sp_sumUinSphere_chi
  subroutine sp_accretionUpdate
    use modelParameter, only : MP_CONNECTION_RUN 
    type(t_spParticle),pointer :: ptr
    real(kind=8) :: mass_new
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       mass_new = ptr%mass + ptr%dmass
       ptr%v = (ptr%v * ptr%mass + ptr%dp) / mass_new
       ptr%s = (ptr%s * ptr%mass + ptr%ds) / mass_new
       ptr%mass = mass_new
       if (MP_CONNECTION_RUN == 0) then 
          ptr%dm_disk = ptr%dm_disk + ptr%dmass
          ptr%dJ_disk = ptr%dJ_disk + ptr%ds
       else 
          ptr%dm_disk = ptr%dm_disk + ptr%mdot_disk * Dtime(LevelMax) 
          ptr%dJ_disk = ptr%dJ_disk + ptr%ds 
       end if
       ptr => ptr%next
    end do
  end subroutine sp_accretionUpdate
  subroutine sp_advectionUpdate
    type(t_spParticle),pointer :: ptr
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       ptr%r = ptr%r + ptr%dr
       ptr%v = ptr%v + ptr%dv
       call get_plev_kobetu(ptr)
       ptr => ptr%next
    end do
  end subroutine sp_advectionUpdate
  subroutine sp_restrictCFL(dt)
    real(kind=8),intent(INOUT) :: dt
    if (Nparticle == 0) return
  end subroutine sp_restrictCFL
  subroutine sp_testIO
    use mpilib
    type(t_spParticle),pointer :: ptr
    if (.not. associated(Particle)) return
    call sp_write
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** a', ptr%pid, ptr%mass 
 call flush(6)
       ptr => ptr%next
    end do
    print *, '** a', Nparticle, PidMax 
 call flush(6)
    call sp_deleteAllParticles
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** b', ptr%pid, ptr%mass 
 call flush(6)
       ptr => ptr%next
    end do
    print *, '** b', Nparticle, PidMax 
 call flush(6)
    call sp_read
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       print *, '** c', ptr%pid, ptr%mass 
 call flush(6)
       ptr => ptr%next
    end do
    print *, '** c', Nparticle, PidMax 
 call flush(6)
  end subroutine sp_testIO
  subroutine sp_refineCond( gid, bool )
    integer,intent(IN) :: gid
    logical,intent(INOUT) :: bool
    integer :: thisLevel
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: margin
    real(kind=8),dimension(0:2) :: hh
    call sp_init
    return
    thisLevel = get_level(gid)
    if (thisLevel + 1 > sp_Level) then 
       bool = .false.
       return
    end if
    if ( bool ) return 
    margin = sp_SinkRadius
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    hh = CellWidth(:,get_level(gid)) * 0.5d0 
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       if ( &
            max(ptr%r(0) - margin, x(Imin)-hh(0)) <= min(ptr%r(0) + margin, x(Imax)+hh(0)) .and. &
            max(ptr%r(1) - margin, y(Jmin)-hh(1)) <= min(ptr%r(1) + margin, y(Jmax)+hh(1)) .and. &
            max(ptr%r(2) - margin, z(Kmin)-hh(2)) <= min(ptr%r(2) + margin, z(Kmax)+hh(2))) then
          bool = .true.
          return
       end if
       ptr => ptr%next
    end do
  end subroutine sp_refineCond
  subroutine sp_refineCond_KS( gid, bool )
    integer,intent(IN) :: gid
    logical,intent(INOUT) :: bool
    integer :: thisLevel
    type(t_spParticle),pointer :: ptr
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: margin
    real(kind=8),dimension(0:2) :: hh
    return
    call sp_init
    thisLevel = get_level(gid)
    if (thisLevel + 1 > sp_Level) then 
       bool = .false.
       return
    end if
    if ( bool ) return 
    margin = sp_SinkRadius 
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    hh = CellWidth(:,get_level(gid))*0.5d0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       if ( &
            max(ptr%r(0) - margin, x(Imin)-hh(0)) <= min(ptr%r(0) + margin, x(Imax)+hh(0)) .and. &
            max(ptr%r(1) - margin, y(Jmin)-hh(1)) <= min(ptr%r(1) + margin, y(Jmax)+hh(1)) .and. &
            max(ptr%r(2) - margin, z(Kmin)-hh(2)) <= min(ptr%r(2) + margin, z(Kmax)+hh(2))) then
          bool = .true.
          return
       end if
       ptr => ptr%next
    end do
  end subroutine sp_refineCond_KS
  SUBROUTINE DSYEVC3(A, W)
    real(kind=8),intent(IN) :: A(3,3)
    real(kind=8),intent(OUT) :: W(3)
    real(kind=8),parameter :: SQRT3 = 1.73205080756887729352744634151D0
    DOUBLE PRECISION M, C1, C0
    DOUBLE PRECISION DE, DD, EE, FF
    DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
    DE = A(1,2) * A(2,3)
    DD = A(1,2)**2
    EE = A(2,3)**2
    FF = A(1,3)**2
    M = A(1,1) + A(2,2) + A(3,3)
    C1 = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
    C0 = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0D0 * A(1,3)*DE
    P = M**2 - 3.0D0 * C1
    Q = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
    SQRTP = SQRT(ABS(P))
    PHI = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1) + C0 * (Q + (27.0D0/4.0D0)*C0) )
    PHI = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)
    C = SQRTP * COS(PHI)
    S = (1.0D0/SQRT3) * SQRTP * SIN(PHI)
    W(2) = (1.0D0/3.0D0) * (M - C)
    W(3) = W(2) + S
    W(1) = W(2) + C
    W(2) = W(2) - S
  END SUBROUTINE DSYEVC3
  function sp_getLevel() result(level)
    integer :: level
    call sp_init
    level = sp_Level
  end function sp_getLevel
  function sp_getRhocr() result(rhocr)
    real(kind=8) :: rhocr
    rhocr = sp_Rhocr
  end function sp_getRhocr
  function sp_getNparticle() result(number)
    integer :: number
    number = Nparticle
  end function sp_getNparticle
  function sp_getSinkRadius() result(SinkRadius)
    real(kind=8) :: SinkRadius
    SinkRadius = sp_SinkRadius
  end function sp_getSinkRadius
  subroutine sp_sinkdata2array(np, pmass, pmdot, pr, pv, ps, pid, pmdot_disk, pJ_disk, pmass_prev &
      , pt_prev, pdm_disk, ptcrt)
    integer,intent(OUT) :: np
    real(kind=8),dimension(:),intent(OUT) :: pmass
    real(kind=8),dimension(0:2,size(pmass)),intent(OUT),optional :: pr, pv, ps
    integer,dimension(size(pmass)),intent(OUT),optional :: pid
    real(kind=8),dimension(:),intent(OUT),optional :: pmdot 
    real(kind=8),dimension(size(pmass)),intent(OUT),optional :: pmdot_disk, pmass_prev, pt_prev, pdm_disk
    real(kind=8),dimension(0:2,size(pmass)),intent(OUT),optional :: pJ_disk
    real(kind=8),dimension(size(pmass)),intent(OUT),optional :: ptcrt
    type(t_spParticle),pointer :: ptr
    np = 0
    ptr => Particle
    do
       if (.not. associated(ptr)) exit
       np = np + 1
       pmass(np) = ptr%mass
       if (present(pr)) pr(:, np) = ptr%r(:)
       if (present(pv)) pv(:, np) = ptr%v(:)
       if (present(ps)) ps(:, np) = ptr%s(:)
       if (present(pid)) pid(np) = ptr%pid
       if (present(pmdot)) then
          if (Dtime(LevelMax) == 0.d0) then
             print *, "sp_sinkdata2array: Dtime(LevelMax) is zero, stopping..."
             stop
          end if
          pmdot(np) = ptr%dmass/Dtime(LevelMax)
       end if
       if (present(pmdot_disk)) pmdot_disk(np) = ptr%mdot_disk
       if (present(pJ_disk)) pJ_disk(:, np) = ptr%J_disk(:)
       if (present(pmass_prev)) pmass_prev(np) = ptr%mass - ptr%dm_disk 
       if (present(pt_prev)) pt_prev(np) = ptr%t_prev
       if (present(pdm_disk)) pdm_disk(np) = ptr%dm_disk
       if (present(ptcrt)) ptcrt(np) = ptr%t_crt
       ptr => ptr%next
    end do
  end subroutine sp_sinkdata2array
  function sp_is_inside_sink(pos) result(bool)
    real(kind=8),dimension(0:2),intent(IN) :: pos
    logical :: bool
    integer :: np
    real(kind=8),dimension(0:2) :: pr
    type(t_spParticle),pointer :: ptr
    np = 0
    ptr => Particle
    bool = .False.
    do
       if (.not. associated(ptr)) exit
       if ( (pos(0)-ptr%r(0))**2 + (pos(1)-ptr%r(1))**2 + (pos(2)-ptr%r(2))**2 < sp_SinkRadius**2 ) &
            bool = .True.
       ptr => ptr%next
    end do
  end function sp_is_inside_sink
  subroutine sp_restrict_within_boundary(posL, posR)
    real(kind=8),dimension(0:2),intent(INOUT) :: posL, posR
    type(t_obRectPhys) :: compDomain, particleDomain
    type(t_obPointPhys) :: pL, pR
    integer :: n
    call ob_computationBoxOfRectPhys( compDomain )
    call ob_extractPointPhysFromRectPhys(pL, compDomain, 'L')
    call ob_extractPointPhysFromRectPhys(pR, compDomain, 'R')
    posL = max(posL, pL%p + CellWidth(:, sp_Level)*0.5d0)
    posR = min(posR, pR%p - CellWidth(:, sp_Level)*0.5d0)
  end subroutine sp_restrict_within_boundary
  subroutine sp_init
    use mpilib
    use parameter, only : Pi
    use modelParameter, only : MP_spNcr,MP_spCs,MP_spRadius_cell,MP_spRadius_lamJ &
      , MP_Gconst, MP_mu, MP_Bondi_radius, MP_Lmax0, MP_Boxsize
    use eos
    real(kind=8) :: csp, jlength, hmax, mjeans, h0
    integer :: level
    if (Initialized) return
    if (get_myrank() == 0) then
       print *, 'initialize sink particle'
       call flush(6)
    end if
    sp_Rhocr = MP_spNcr * cgs_amu * MP_mu / Unit_rho 
    sp_RhocrCreate = sp_Rhocr 
      level = MP_Lmax0
      hmax = maxval(CellWidth(:,Lmin))/2.d0**(level-Lmin)
      sp_SinkRadius = MP_spRadius_cell * hmax
      sp_Level = MP_Lmax0
      if (get_myrank() == 0) then
         print '(A,I4,1P4E12.4)', "level, hmax, jlength, jlength/hmax, Mjeans: ", level, hmax, jlength, jlength/hmax, mjeans
      end if
    if (sp_Level == Undefi) then 
       if (get_myrank() == 0) &
            print '(/,A,/,A,/)', "cannot determine sp_Level.", "stopping..."
       stop
    endif
    if (get_myrank() == 0) then
       print *, 'sp_SinkRadius = ', sp_SinkRadius
       print *, 'sp_SinkRadius[pc]= ', sp_SinkRadius*Unit_pc
       print *, 'sp_SinkRadius[au]= ', sp_SinkRadius*Unit_au
       print *, 'sp_Level = ', sp_Level
       print *, 'MP_spRadius_cell = ', MP_spRadius_cell
       print *, 'LevelMax = ', LevelMax
       call flush(6)
    end if
    sp_SofteningRadius = sp_SinkRadius
    if (get_myrank() == 0) then
       print *, 'sp_SofteningRadius = ', sp_SofteningRadius
       call flush(6)
    end if
    sp_MaskRadius = sp_SinkRadius * 2
    if (get_myrank() == 0) then
       print *, 'sp_MaskRadius = ', sp_MaskRadius
       call flush(6)
    end if
    call write_spparam
    hh_sp(:) = CellWidth(:, Lmin) / (2.d0**(sp_Level-Lmin)) 
    Initialized = .true.
contains
  subroutine write_spparam
    use string, only : concat, CHARLEN
    use io_util, only : readenv, wchar
    character(len=CHARLEN) :: fn, dir
    if (.not. readenv('DIR', dir) ) stop
    if (get_myrank() == 0) then
       fn = concat(dir, "sp_param.txt")
       call wchar(6,'write parameter file = '//fn)
       open(1, file=fn)
       write(1,*) sp_SinkRadius
       write(1,*) sp_Level
       write(1,*) sp_SofteningRadius
       write(1,*) sp_MaskRadius
       close(1)
    end if
  end subroutine write_spparam
end subroutine sp_init
subroutine sp_update_subdisk
  type(t_spParticle),pointer :: ptr
  real(kind=8) :: t
  logical :: isNotFinite
  t = Time(LevelMax) 
  ptr => Particle
  do
     if (.not. associated(ptr)) exit
     if (t-ptr%t_prev > dt_acc/Unit_yr) then
        ptr%mdot_disk = ptr%dm_disk/(t-ptr%t_prev) 
        ptr%J_disk = ptr%dJ_disk 
        ptr%t_prev = t 
        ptr%dm_disk = 0d0 
        ptr%dJ_disk = 0d0 
     end if
     ptr => ptr%next
  end do
end subroutine sp_update_subdisk
subroutine get_plev
  real(kind=8),dimension(0:2) :: pos
  type(t_spParticle),pointer :: ptr
  integer :: np, ilev, ilgid, rank, i, j, k
  np = 0
  ptr => Particle
  do
    if (.not. associated(ptr)) exit
    np = np + 1
    pos = ptr%r
    do ilev = LevelMax, Lmin, -1 
      call ob_getBlockFromCoordPhys(pos, ilev, ilgid, rank, i, j, k)
      if (ilgid /= Undefi) then
        ptr%lev = ilev
        exit 
      endif
      if (ilgid == Undefi .and. ilev == Lmin) then 
        ptr%lev = Lmin
        exit
      endif
    enddo
    ptr%h(:) = hh_sp(:)
    ptr => ptr%next
  enddo
end subroutine get_plev
subroutine get_plev_kobetu(ptr)
  real(kind=8),dimension(0:2) :: pos
  type(t_spParticle),pointer :: ptr
  integer :: ilev, ilgid, rank, i, j, k
  pos = ptr%r
  do ilev = LevelMax, Lmin, -1 
    call ob_getBlockFromCoordPhys(pos, ilev, ilgid, rank, i, j, k)
    if (ilgid /= Undefi) then
      ptr%lev = ilev
      exit 
    endif
    if (ilgid == Undefi .and. ilev == Lmin) then 
      ptr%lev = Lmin
      exit
    endif
  enddo
  ptr%h(:) = hh_sp(:)
end subroutine get_plev_kobetu
  subroutine output_spcrt(list_newPos)
    use grid, only : LevelMax, Lmin, CellWidth, Undefi, Step
    use overBlockCoordinates, only : ob_computationBoxOfCoordPhys, OB_COORDS_MIN, OB_COORDS_MAX
    use uniformgrid, only : uniformgrid_write
    use modelParameter, only : MP_Boxsize
    use unit, only : Unit_msun
    use mpilib
    use string, only : CHARLEN, num2char, concat
    type(t_spPos),pointer :: list_newPos 
    type(t_spPos),pointer :: ptr
    integer,parameter :: NIug=64, NJug=NIug, NKug=NIug 
    real(kind=8) :: xmin, ymin, zmin, xmax, ymax, zmax, xp, yp, zp
    real(kind=8) :: halfwidth, output_width
    integer :: level, dummy
    real(kind=8) :: coordPhys(OB_COORDS_MIN:OB_COORDS_MAX)
    character(len=CHARLEN) :: prefix
    integer :: ntot, npid
    output_width = 2.d0*cgs_pc/Unit_l
    call ob_computationBoxOfCoordPhys( coordPhys )
    ptr => list_newPos
    ntot = 0
    do
      if (.not. associated(ptr)) exit
      ntot = ntot + 1
      ptr => ptr%next
    enddo
    if (ntot < 1) return
    npid = Nparticle-ntot
    ptr => list_newPos
    do
      if (.not. associated(ptr)) exit
      do level = Lmin-2,LevelMax 
        halfwidth = MP_Boxsize / 2.**(level-Lmin) / (8*8/dble(NIug))
        if(halfwidth > output_width) cycle
        if (get_myrank() == 0) then
             print '(/,A,I0,A,I0,A,(1P1E9.2),A)', "outputnewparticle: pid = ", npid &
               , " level = ", level, ", size = ", halfwidth/MP_Boxsize, " x comp. region"
        endif
          xmin = coordPhys(0)
          ymin = coordPhys(1)
          zmin = coordPhys(2)
          xmax = coordPhys(2 +1+0)
          ymax = coordPhys(2 +1+1)
          zmax = coordPhys(2 +1+2)
          xmin = max(xmin, ptr%r(0)-halfwidth)
          ymin = max(ymin, ptr%r(1)-halfwidth)
          zmin = max(zmin, ptr%r(2)-halfwidth)
          xmax = min(xmax, ptr%r(0)+halfwidth)
          ymax = min(ymax, ptr%r(1)+halfwidth)
          zmax = min(zmax, ptr%r(2)+halfwidth)
          prefix = 'crsp.'
          prefix = concat(concat(prefix, num2char(npid)),'.')
          call uniformgrid_write(xmin, ymin, zmin, xmax, ymax, zmax, level, interpolate=.false.,prefix=prefix)
      enddo
      npid = npid + 1
      ptr => ptr%next
    enddo
  end subroutine output_spcrt
end module sinkParticle
