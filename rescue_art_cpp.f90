module rescue
  use grid, only: Lmin, Lmax
  implicit none
  private
  real(kind=8),parameter :: V_UPPER_LIMIT = 1.d2 
  public :: rescueLev, rescueAllLev, rescue_rhopsi_for_poisson 
contains
  subroutine rescueLev(level)
    integer :: level
    call rescueNanInf( level )
    call rescueByFloor ( level ) 
  end subroutine rescueLev
  subroutine rescueAllLev
    use grid, only : Lmin, LevelMax
    integer :: level
    do level = Lmin, LevelMax
       call rescueLev( level )
    enddo
  end subroutine rescueAllLev
  subroutine rescueNanInf( level )
    use grid, only : Gidmin, GidListMax, GidList, Mmin, Mmax, get_Ucomp, has_child_grid, get_Xp, get_Yp, get_Zp, get_Up
    use mpilib 
    integer,intent(IN) :: level
    integer :: n, gid, m, ncell, i, j, k, is, js, ks, ie, je, ke, ii, jj, kk, nw
    logical :: allIsFinite, isFinite 
    real(kind=8),dimension(:,:,:),pointer :: a
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8) :: ave
    logical :: bool_realgrid, bool_write
    real(kind=8),dimension(:,:,:,:),pointer :: u 
    logical :: isNotFinite 
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level) 
       bool_realgrid = .not. has_child_grid(gid)
       bool_write = bool_realgrid
       do m = Mmin, Mmax
          a => get_Ucomp(m, gid)
          if ( allIsFinite( a , size(a) ) ) cycle
          ks = lbound(a,3)
          ke = ubound(a,3)
          js = lbound(a,2)
          je = ubound(a,2)
          is = lbound(a,1)
          ie = ubound(a,1)
          x => get_Xp(gid)
          y => get_Yp(gid)
          z => get_Zp(gid)
          do k = ks, ke
             do j = js, je
                do i = is, ie
                   if ( isFinite(a(i,j,k)) ) cycle
                   u => get_Up(gid)
                   print '(A, 5I5, (1P11E15.7))', "(rescueNanInf) KS DEBUG", &
                        get_myrank(),gid,i-lbound(u,1),j-lbound(u,2),k-lbound(u,3), &
                        u(i,j,k,4), u(i,j,k,0), u(i,j,k,1), u(i,j,k,2), u(i,j,k,3),&
                        x(i), y(j), z(k),u(i,j,k,14), u(i,j,k,15),u(i,j,k,14)
                   if (isNotFinite(u(i,j,k,4)) .or. isNotFinite(u(i,j,k,0)) .or. isNotFinite(u(i,j,k,1)) &
                        .or. isNotFinite(u(i,j,k,2)) .or. isNotFinite(u(i,j,k,3))) then
                      print '(A,/,A)', '(rescue KS DEBUG) NaN/Inf found','stopping'
                      stop
                   end if
                   do nw = 1, max(ie-is+1, je-js+1, ke-ks+1) 
                      ave = 0.d0
                      ncell = 0
                      do kk = max(k-nw, ks), min(k+nw, ke)
                         do jj = max(j-nw, js), min(j+nw, je)
                            do ii = max(i-nw, is), min(i+nw, ie)
                               if ( isFinite(a(ii,jj,kk)) ) then
                                  ave = ave + a(ii,jj,kk)
                                  ncell = ncell + 1
                               endif
                            enddo
                         enddo
                      enddo
                      if ( ncell > 0 ) exit
                   enddo
                   if (ncell == 0) then
                      write(*, *) "*** Error in rescue3D, ncell == 0"
                      stop
                   endif
                   a(i,j,k) = ave/ncell
                   if (bool_write) &
                        write(*,'("*** rescueNanInf (realg, level, gid, x, y, z, m, nw) ", L, I3, I4, 3(1PE14.6), I3, I3)') &
                        bool_realgrid, level, gid, x(i), y(j), z(k), m, nw
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine rescueNanInf
  subroutine setUpperLimitVelocity( level )
    use grid, only : Gidmin, GidListMax, GidList, get_Ucomp, has_child_grid
    integer,intent(IN) :: level
    real(kind=8),dimension(:,:,:),pointer :: v
    real(kind=8) :: vmax
    integer :: n, gid, m
    real(kind=8) :: xloc, yloc, zloc
    logical :: bool_realgrid, bool_write
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       bool_realgrid = .not. has_child_grid(gid)
       bool_write = bool_realgrid
       do m = 1, 3
          v => get_Ucomp(m, gid)
          vmax = maxval(abs(v))
          if ( vmax >= V_UPPER_LIMIT ) then
             call gridloc(gid, xloc, yloc, zloc)
             if (bool_write) then
                write(*, '("*** rescue: V limitted (realg, level, gid, x, y, z) ", L, I3, I4, 3(1PE14.6))') &
                     bool_realgrid, level, gid, xloc, yloc, zloc
             endif
             where( v > V_UPPER_LIMIT ) v = V_UPPER_LIMIT
             where( v < -V_UPPER_LIMIT ) v = -V_UPPER_LIMIT
          end if
       end do
    end do
  end subroutine setUpperLimitVelocity
  subroutine setDiffusion( level )
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax, Mmin, Mmax, get_Ucomp, get_Up, Gidmin, GidListMax, GidList, get_Xp, get_Yp&
&, get_Zp, has_child_grid, Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh
    use modelParameter, only : MP_Tmin, MP_Tmax 
    use unit
    use primordial, only : yHe
    integer,intent(IN) :: level
    integer,dimension(3) :: mlistvel = (/1, 2, 3/)
    integer,dimension(1) :: mlistpos_r = (/0/)
    integer,dimension(1) :: mlistpos_p = (/4/)
    integer,parameter :: NPHASE_MIN = 1
    integer,parameter :: NPHASE_MAX = 3
    real(kind=8),dimension(:,:,:,:),pointer :: u
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax) :: w
    real(kind=8),dimension(Mmin:Mmax) :: wave, uave
    integer :: n, gid, i, j, k, ii, jj, kk, is, ie, js, je, ks, ke, nw, m, nphase, np, ncell, ntrouble, ntrouble_r, ntrouble_p, ntr&
&ouble_v
    real(kind=8),dimension(:),pointer :: x, y, z
    logical :: bool_realgrid, bool_write, bool_allok
    logical :: isNotFinite 
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       u => get_Up(gid)
       if (&
            all(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistpos_p(:)) > 0.d0) .and. & 
            all(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistpos_r(:)) > 0.d0) .and. & 
            all(abs(u(Imin:Imax,Jmin:Jmax,Kmin:Kmax,mlistvel(:))) < V_UPPER_LIMIT)) cycle 
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       bool_realgrid = .not. has_child_grid(gid)
       bool_write = bool_realgrid
       w(:,:,:,:) = u(:,:,:,:)
       w(:,:,:,1) = u(:,:,:,1)*u(:,:,:,0)
       w(:,:,:,2) = u(:,:,:,2)*u(:,:,:,0)
       w(:,:,:,3) = u(:,:,:,3)*u(:,:,:,0)
       do k = Kmin, Kmax
          do j = Jmin, Jmax
             do i = Imin, Imax
                ntrouble_p = 0 
                ntrouble_r = 0 
                ntrouble_v = 0 
                if (.not. all(u(i,j,k,mlistpos_p(:)) > 0.d0)) ntrouble_p = 1
                if (.not. all(u(i,j,k,mlistpos_r(:)) > 0.d0)) ntrouble_r = 1
                if (.not. all(abs(u(i,j,k,mlistvel(:))) < V_UPPER_LIMIT)) ntrouble_v = 1
                ntrouble = ntrouble_p + ntrouble_r * 2 + ntrouble_v * 4
                if (ntrouble == 0) cycle
                print '(A, 4I5,(1P8E15.7))', "(rescue) KS DEBUG", &
                     gid, i-lbound(u,1),j-lbound(u,2),k-lbound(u,3),&
                     u(i,j,k,0), u(i,j,k,4), u(i,j,k,1), u(i,j,k,2), u(i,j,k,3),&
                     x(i), y(j), z(k)
                if (isNotFinite(u(i,j,k,4)) .or. isNotFinite(u(i,j,k,0)) .or. isNotFinite(u(i,j,k,1)) &
                     .or. isNotFinite(u(i,j,k,2)) .or. isNotFinite(u(i,j,k,3))) then
                   print '(A,/,A)', '(rescue KS DEBUG) NaN/Inf found','stopping'
                   stop
                end if
                do nw = 1, max(Imax-Imin+1, Jmax-Jmin+1, Kmax-Kmin+1) 
                   do np = NPHASE_MIN, NPHASE_MAX 
                      nphase = np
                      is = max(i-nw, Imin)
                      ie = min(i+nw, Imax)
                      js = max(j-nw, Jmin)
                      je = min(j+nw, Jmax)
                      ks = max(k-nw, Kmin)
                      ke = min(k+nw, Kmax)
                      if (nphase == 1) then 
                         do m = Mmin, Mmax
                            wave(m) = sum( w(is:ie, js:je, ks:ke, m) )
                         end do
                         wave(:) = wave(:) / ((ie-is+1)*(je-js+1)*(ke-ks+1))
                         uave = wave
                         uave(1:3) = wave(1:3)/wave(0)
                      else if (nphase == 2) then 
                         do m = Mmin, Mmax
                            wave(m) = sum( w(is:ie, js:je, ks:ke, m) )
                         end do
                         wave(:) = wave(:) - w(i,j,k,:)
                         wave(:) = wave(:) / ((ie-is+1)*(je-js+1)*(ke-ks+1)-1)
                         uave = wave
                         uave(1:3) = wave(1:3)/wave(0)
                      else if (nphase == 3) then 
                         uave = 0.d0
                         ncell = 0
                         do kk = ks, ke
                            do jj = js, je
                               do ii = is, ie
                                  if (ii == i .and. jj == j .and. kk == k) cycle
                                  if (&
                                       all(u(ii,jj,kk,mlistpos_p(:)) > 0.d0) .and. &
                                       all(u(ii,jj,kk,mlistpos_r(:)) > 0.d0) .and. &
                                       all(abs(u(ii,jj,kk,mlistvel(:))) < V_UPPER_LIMIT) &
                                       ) then
                                     uave(:) = uave(:) + u(ii,jj,kk,:)
                                     ncell = ncell + 1
                                  end if
                               end do
                            end do
                         end do
                         if (ncell > 0) then
                            uave(:) = uave(:)/ncell
                         endif
                      end if
                      bool_allok = ( &
                           all(uave(mlistpos_p(:)) > 0.d0) .and. &
                           all(uave(mlistpos_r(:)) > 0.d0) .and. &
                           all(abs(uave(mlistvel(:))) < V_UPPER_LIMIT))
                      if ( bool_allok ) exit
                   end do 
                   if ( bool_allok ) then
                      if (nphase == 1) then
                         do m = Mmin, Mmax 
                            u(is:ie,js:je,ks:ke,m) = uave(m)
                         end do
                      else 
                         if (ntrouble_p == 1) u(i,j,k,mlistpos_p) = uave(mlistpos_p) 
                         if (ntrouble_r == 1 .or. ntrouble_v == 1) u(i,j,k,:) = uave(:) 
                      endif
                      exit
                   endif
                end do 
                if ( .not. bool_allok ) then 
                   nphase = nphase + 1
                   if (.not. all(uave(mlistpos_p(:)) > 0.d0)) then
                      u(i,j,k,mlistpos_p) = abs(w(i,j,k,mlistpos_p))
                   endif
                   if (.not. all(uave(mlistpos_r(:)) > 0.d0)) then
                      u(i,j,k,mlistpos_r) = abs(w(i,j,k,mlistpos_r))
                      u(i,j,k,1:3) = w(i,j,k,1:3)/u(i,j,k,0)
                   endif
                   if (.not. all(abs(u(i,j,k,mlistvel(:))) < V_UPPER_LIMIT) ) then
                      u(i,j,k,mlistvel) = w(i,j,k,mlistvel)/sqrt(sum(w(i,j,k,mlistvel)**2)) * V_UPPER_LIMIT
                   endif
                   nw = 1
                endif
                if (bool_write) &
                     write(*,'("*** phase",I2," rescue diff (realg, level, gid, x, y, z, nw, ntb) ", L, I3, I4, 3(1PE14.6), I3, I3)&
&') &
                     nphase, bool_realgrid, level, gid, x(i), y(j), z(k), nw, ntrouble
             end do
          end do
       end do
    end do
  end subroutine setDiffusion
  subroutine gridloc(gid, xloc, yloc, zloc)
    use grid
    integer,intent(IN) :: gid
    real(kind=8),intent(OUT) :: xloc, yloc, zloc
    real(kind=8),dimension(:),pointer :: x, y, z
    x => get_Xp(gid)
    y => get_Yp(gid)
    z => get_Zp(gid)
    xloc = (x(Imin)+x(Imax))*0.5d0
    yloc = (y(Jmin)+y(Jmax))*0.5d0
    zloc = (z(Kmin)+z(Kmax))*0.5d0
  end subroutine gridloc
  subroutine rescue_rhopsi_for_poisson
    use grid, only : Lmin, LevelMax
    integer :: level
    do level = Lmin, LevelMax
       call check_rhopsi_for_poisson( level )
    enddo
  end subroutine rescue_rhopsi_for_poisson
  subroutine check_rhopsi_for_poisson( level )
    use grid, only : Imin, Imax, Jmin, Jmax, Kmin, Kmax, get_Ucomp, get_U1comp, get_U2comp, Gidmin, GidListMax, GidList, Imingh, Im&
&axgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh,&
         get_Xp, get_Yp, get_Zp
    use mpilib
    use unit
    integer,intent(IN) :: level
  real(kind=8),dimension(:,:,:),pointer :: rho,rho1,rho2,psi
    real(kind=8),dimension(:),pointer :: x, y, z
    integer :: n, gid, i, j, k
    logical :: isNotFinite 
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level)
       rho1 => get_U1comp(0,gid)
       psi => get_Ucomp(17,gid)
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       do k = Kmingh, Kmaxgh
          do j = Jmingh, Jmaxgh
             do i = Imingh, Imaxgh
                if (.not. (rho1(i,j,k)*Unit_rho/cgs_mh > 1d-2 .and. rho1(i,j,k)*Unit_rho/cgs_mh < 1d16)) then 
                   print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "** WARNING ** (check U1_rho) ", &
                        rho1(i,j,k)*Unit_rho/cgs_mh, i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                end if
                if (.not. (psi(i,j,k) > -1d2 .and. psi(i,j,k) < 1d0)) then 
                   print '(A, (1P1E15.7), 3I4, (1P3E15.7))', "** WARNING ** (check U_psi) ", &
                        psi(i,j,k), i,j,k, x(i)*Unit_au, y(j)*Unit_au, z(k)*Unit_au
                end if
             end do
          end do
       end do
    end do
  end subroutine check_rhopsi_for_poisson
  subroutine rescueByFloor( level )
    use grid, only : Gidmin, GidListMax, GidList, get_Ucomp, has_child_grid, get_Xp, get_Yp, get_Zp, Imin, Imax, Jmin, Jmax, Kmin, &
&Kmax,&
         Imingh, Imaxgh, Jmingh, Jmaxgh, Kmingh, Kmaxgh
    use modelParameter, only : MP_Tmin, MP_Tmax, MP_Nmin, MP_Vmax, MP_metallicity, MP_mu, MP_frac_COsum
    use unit
    use primordial
    use mpilib
    integer,intent(IN) :: level
    integer :: n, gid, i, j, k, is, ie, js, je, ks, ke
    real(kind=8),dimension(:),pointer :: x, y, z
    logical :: bool_realgrid, bool_gh
    real(kind=8),dimension(:,:,:),pointer :: rho, p, vx, vy, vz
    real(kind=8),dimension(:,:,:),pointer :: yco
    real(kind=8) :: T_K, xmu, v2, xnH, yco_l
    real(kind=8),dimension(0:6 -1) :: ychem
    integer :: count_output_NL=0, count_output_TL=0, count_output_TU=0, count_output_VU=0, max_count=100 
    type chemsp
      real(kind=8),dimension(:,:,:),pointer :: y
    end type chemsp
    type(chemsp), dimension(:) :: chemary3(5:10)
    integer :: ichem
    do n = Gidmin, GidListMax( level )
       gid = GidList(n, level) 
       bool_realgrid = .not. has_child_grid(gid)
       x => get_Xp(gid)
       y => get_Yp(gid)
       z => get_Zp(gid)
       rho => get_Ucomp(0,gid)
       p => get_Ucomp(4,gid)
       do ichem = 5, 10
         chemary3(ichem)%y => get_Ucomp(ichem,gid)
       enddo
       yco => get_Ucomp(11, gid)
       vx => get_Ucomp(1,gid)
       vy => get_Ucomp(2,gid)
       vz => get_Ucomp(3,gid)
       do i = Imingh, Imaxgh
          do j = Jmingh, Jmaxgh
             do k = Kmingh, Kmaxgh
                bool_gh = .not. (i<Imin .or. i>Imax .or. j<Jmin .or. j>Jmax .or. k<Kmin .or. k>Kmax)
                do ichem = 0, 6 -1
                   ychem(ichem) = chemary3(ichem+5)%y(i,j,k)
                enddo
                yco_l = yco(i,j,k)
                call adjust_abundance(ychem &
                , yco_l &
                , MP_Metallicity)
                do ichem = 0, 6 -1
                   chemary3(ichem+5)%y(i,j,k) = ychem(ichem)
                enddo
                yco(i,j,k) = yco_l
                xnH = rho(i,j,k)*Unit_rho/(MP_mu*cgs_amu) 
                if (.not. xnH >= MP_Nmin) then
                   if (count_output_NL <= max_count) then
                      print '(A, 5I5, 2L, I3, 5(1PE14.6))', "*** rescueByFloor (N lower floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, xnH, MP_Nmin, x(i),y(j),z(k)
                      if (count_output_NL == max_count) &
                           print *, 'count_output_NL reaches max_count -> no more warning shown'
                      count_output_NL = count_output_NL+1
                   end if
                   xnH = MP_Nmin
                   rho(i,j,k) = xnH * (MP_mu*cgs_amu) / Unit_rho
                   xmu = get_xmu(ychem) 
                   T_K = sqrt(MP_Tmin*MP_Tmax) 
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e
                end if
                xmu = get_xmu(ychem) 
                T_K = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb 
                if (.not. T_K >= MP_Tmin) then
                   if (count_output_TL <= max_count) then
                      print '(A, 5I5, 2L, I3, 8(1PE14.6))', "*** rescueByFloor (T lower floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, p(i,j,k), rho(i,j,k), xmu, T_K, MP_Tmin, x(i),y(j),z(k)
                      if (count_output_TL == max_count) &
                           print *, 'count_output_TL reaches max_count -> no more warning shown'
                      count_output_TL = count_output_TL+1
                   end if
                   T_K = MP_Tmin
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e
                end if
                if (.not. T_K <= MP_Tmax) then
                   if (count_output_TU <= max_count) then
                      print '(A, 5I5, 2L, I3, 8(1PE14.6))', "*** rescueByFloor (T upper floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, p(i,j,k), rho(i,j,k), xmu, T_K, MP_Tmax, x(i),y(j),z(k)
                      if (count_output_TU == max_count) &
                           print *, 'count_output_TU reaches max_count -> no more warning shown'
                      count_output_TU = count_output_TU+1
                   end if
                   T_K = MP_Tmax
                   p(i,j,k) = (cgs_kb*T_K)*(rho(i,j,k)*Unit_rho)/(cgs_amu*xmu) / Unit_e
                end if
                v2 = (vx(i,j,k)*vx(i,j,k) + vy(i,j,k)*vy(i,j,k) + vz(i,j,k)*vz(i,j,k))*Unit_kms*Unit_kms
                if (.not. v2 <= MP_Vmax*MP_Vmax) then
                   if (count_output_VU <= max_count) then
                      print '(A, 5I5, 2L, I3, 7(1PE14.6))', "*** rescueByFloor (V upper floor): ", &
                           get_myrank(),gid,i-lbound(rho,1),j-lbound(rho,2),k-lbound(rho,3), &
                           bool_realgrid, bool_gh, level, vx(i,j,k)*Unit_kms,vy(i,j,k)*Unit_kms,vz(i,j,k)*Unit_kms, MP_Vmax, x(i),y&
&(j),z(k)
                      if (count_output_VU == max_count) &
                           print *, 'count_output_VU reaches max_count -> no more warning shown'
                      count_output_VU = count_output_VU+1
                   end if
                   vx(i,j,k) = MP_Vmax/sqrt(v2) * vx(i,j,k)
                   vy(i,j,k) = MP_Vmax/sqrt(v2) * vy(i,j,k)
                   vz(i,j,k) = MP_Vmax/sqrt(v2) * vz(i,j,k)
                end if
             end do
          end do
       end do
    enddo
  end subroutine rescueByFloor
end module rescue
