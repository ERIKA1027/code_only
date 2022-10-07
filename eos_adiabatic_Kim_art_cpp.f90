module eos
  use modelParameter, only : MP_Tmin, MP_Tmax, MP_Ctil_nD, MP_mu, MP_frac_COsum 
  implicit none
contains
  function cyclecomp(ncrd, invert) result( mcycle )
    use grid , only : Mmin, Mmax
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(Mmin:Mmax) :: mcycle
    integer,dimension(0:2) :: mcycle3 = (/ 1, 2, 3 /)
    integer :: m
    do m = Mmin, Mmax
       mcycle(m) = m
    enddo
    if ( present( invert ) ) then
       mcycle(1:1 +size(mcycle3)-1) = cshift( mcycle3, -ncrd)
    else
       mcycle(1:1 +size(mcycle3)-1) = cshift( mcycle3, ncrd )
    endif
  end function cyclecomp
  subroutine flux( ql, qr, pratio, f, ncrd, gam )
    use grid, only : Imin,Jmin,Kmin,Mmin, Imax,Jmax,Kmax,Mmax, Imingh,Jmingh,Kmingh, Imaxgh,Jmaxgh,Kmaxgh, &
         globdbg_myrank, globdbg_mygid,globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use util, only : util_arroffset
    integer,intent(IN) :: ncrd
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: ql 
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: qr 
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:,0:) :: pratio 
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:,Mmin:) :: f 
    integer :: i,j,k, is,js,ks,ie,je,ke,io,jo,ko
    real(kind=8) :: &
         rhol,ul,vl,wl,pl,hl,el,cl, &
         rhor,ur,vr,wr,pr,hr,er,cr, &
         drho, drhou, drhov, drhow, drhoh, dp, dh, du, dv, dw, &
         sql,sqr,sqa,rhob,ub,vb,wb,hb,qb2,cb2,cb,ub2,vb2,wb2, &
         gm1,mn, b1, b2, b1b2, db12, ff, gg, rbdq, hh, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5
    real(kind=8) :: el4c, a7c, a8c, a9c, w1c
    real(kind=8) :: eps, x1, x2
    integer :: ichem
    real(kind=8),dimension(Imingh:,Jmingh:,Kmingh:),intent(IN) :: gam
    f(:,:,:,:) = 0.d0
    call util_arroffset(ncrd,io,jo,ko)
    ks = Kmin-ko
    ke = Kmax
    js = Jmin-jo
    je = Jmax
    is = Imin-io
    ie = Imax
    do k = ks, ke
       do j = js, je
          do i = is, ie
             gm1 = gam(i,j,k) - 1
             rhol = abs(ql(i,j,k,0)) 
             ul = ql(i,j,k,1)
             vl = ql(i,j,k,2)
             wl = ql(i,j,k,3)
             pl = abs(ql(i,j,k,4))
             el = rhol * (ul**2 + vl**2 + wl**2) / 2 + pl / gm1
             hl = ( el + pl )/rhol
             cl = sqrt( gam(i,j,k) * pl / rhol )
             rhor = abs(qr(i,j,k,0))
             ur = qr(i,j,k,1)
             vr = qr(i,j,k,2)
             wr = qr(i,j,k,3)
             pr = abs(qr(i,j,k,4))
             er = rhor * (ur**2 + vr**2 + wr**2) / 2 + pr / gm1
             hr = ( er + pr )/rhor
             cr = sqrt( gam(i,j,k) * pr / rhor )
             drho = rhor - rhol
             drhou = rhor * ur - rhol * ul
             drhov = rhor * vr - rhol * vl
             drhow = rhor * wr - rhol * wl
             drhoh = rhor * hr - rhol * hl
             du = ur - ul
             dv = vr - vl
             dw = wr - wl
             dp = pr - pl
             dh = hr - hl
             sql = sqrt(rhol)
             sqr = sqrt(rhor)
             sqa = sql + sqr
             rhob = sql * sqr
             ub = (sql * ul + sqr * ur) / sqa
             vb = (sql * vl + sqr * vr) / sqa
             wb = (sql * wl + sqr * wr) / sqa
             hb = (sql * hl + sqr * hr) / sqa
             ub2 = ub**2
             vb2 = vb**2
             wb2 = wb**2
             qb2 = ub2 + vb2 + wb2
             cb2 = gm1*(hb - qb2/2)
             cb = sqrt( cb2 )
             mn = abs(ub/cb) 
             b1 = max(0.d0, ub+cb, ur+cb)
             b2 = min(0.d0, ub-cb, ul-cb)
             b1b2 = b1 * b2
             db12 = b1 - b2
             if (ncrd == 0) then
                hh = 1.d0 - min(pratio(i,j,k,0), &
                     pratio(i ,j-1,k ,1), pratio(i+1,j-1,k ,1), &
                     pratio(i ,j ,k ,1), pratio(i+1,j ,k ,1), &
                     pratio(i ,j ,k-1,2), pratio(i+1,j ,k-1,2), &
                     pratio(i ,j ,k ,2), pratio(i+1,j ,k ,2))
             else if (ncrd == 1) then
                hh = 1.d0 - min(pratio(i,j,k,1), &
                     pratio(i ,j ,k-1,2), pratio(i ,j+1,k-1,2), &
                     pratio(i ,j ,k ,2), pratio(i ,j+1,k ,2), &
                     pratio(i-1,j ,k ,0), pratio(i-1,j+1,k ,0), &
                     pratio(i ,j ,k ,0), pratio(i ,j+1,k ,0))
             else
                hh = 1.d0 - min(pratio(i,j,k,2), &
                     pratio(i-1,j ,k ,0), pratio(i-1,j ,k+1,0), &
                     pratio(i ,j ,k ,0), pratio(i ,j ,k+1,0), &
                     pratio(i ,j-1,k ,1), pratio(i ,j-1,k+1,1), &
                     pratio(i ,j ,k ,1), pratio(i ,j ,k+1,1))
             endif
             if (qb2 == 0.d0 ) then
                ff = 1.d0
             else
                ff = mn**hh
             endif
             if (mn == 0.d0) then
                gg = 1.d0
             else
                gg = mn**(1.d0-pratio(i,j,k,ncrd))
             endif
             rbdq = drho - ff * dp / cb2
             fl1 = rhol * ul
             fl2 = fl1 * ul + pl
             fl3 = fl1 * vl
             fl4 = fl1 * wl
             fl5 = ul * ( el + pl )
             fr1 = rhor * ur
             fr2 = fr1 * ur + pr
             fr3 = fr1 * vr
             fr4 = fr1 * wr
             fr5 = ur * ( er + pr )
             f(i,j,k,0) = (b1*fl1 - b2*fr1)/db12 &
                  + b1b2/db12 * drho &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq
             f(i,j,k,1) = (b1*fl2 - b2*fr2)/db12 &
                  + b1b2/db12 * drhou &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq * ub
             f(i,j,k,2) = (b1*fl3 - b2*fr3)/db12 &
                  + b1b2/db12 * drhov &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * vb + rhob * dv)
             f(i,j,k,3) = (b1*fl4 - b2*fr4)/db12 &
                  + b1b2/db12 * drhow &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * wb + rhob * dw)
             f(i,j,k,4) = (b1*fl5 - b2*fr5)/db12 &
                  + b1b2/db12 * drhoh &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * hb + rhob * dh)
             do ichem = 5, 10
                f(i,j,k,ichem) = 0.5d0 * ub * (ql(i,j,k,ichem)*rhol + qr(i,j,k,ichem)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,ichem)*rhor - ql(i,j,k,ichem)*rhol)
             enddo
             f(i,j,k,11) = 0.5d0 * ub * (ql(i,j,k,11)*rhol + qr(i,j,k,11)*rhor) &
                  - 0.5d0 * abs(ub) * (qr(i,j,k,11)*rhor - ql(i,j,k,11)*rhol)
          enddo
       enddo
    enddo
  end subroutine flux
  function get_dt_by_cflcond( id ) result( dt )
    use grid
    integer,intent(IN) :: id
    real(kind=8) :: dt
    integer :: i,j,k,level
    real(kind=8),dimension(0:2) :: h
    real(kind=8),dimension(:,:,:),pointer :: rho, vx, vy, vz, p, gx, gy, gz
    real(kind=8),dimension(:),pointer :: x, y, z
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: cs
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh) :: gam
    real(kind=8):: dt_radtr
    level = get_level(id)
    h = CellWidth( :, level )
    rho => get_Ucomp(0,id)
    vx => get_Ucomp(1,id)
    vy => get_Ucomp(2,id)
    vz => get_Ucomp(3,id)
    p => get_Ucomp(4,id)
    x => get_Xp(id)
    y => get_Yp(id)
    z => get_Zp(id)
    call get_gamma(id, gam)
    cs = sqrt(abs(gam*p/rho))
    dt = (0.7) / maxval( &
         (abs(vx)+cs)/h(0) + &
         (abs(vy)+cs)/h(1) + &
         (abs(vz)+cs)/h(2), mask=GridMask )
    dt_radtr = minval(h)*0.8/ (3.d0*MP_Ctil_nD)*10.0
    dt = min(dt, dt_radtr)
  end function get_dt_by_cflcond
  subroutine u2w(u,w,dv)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),intent(IN) :: dv
    integer :: ichem
    real(kind=8) :: T_K, xmu, gam
    integer :: i,j,k
    w(:,:,:,0)=u(:,:,:,0)*dv
    w(:,:,:,1)=u(:,:,:,1)*w(:,:,:,0)
    w(:,:,:,2)=u(:,:,:,2)*w(:,:,:,0)
    w(:,:,:,3)=u(:,:,:,3)*w(:,:,:,0)
    do ichem = 5, 10
      w(:,:,:,ichem)=u(:,:,:,ichem)*w(:,:,:,0)
    enddo
    w(:,:,:,11)=u(:,:,:,11)*w(:,:,:,0)
    w(:,:,:,12)=u(:,:,:,12) * dv
    w(:,:,:,13)=u(:,:,:,13) * dv
    w(:,:,:,14)=u(:,:,:,14)*w(:,:,:,0)
    w(:,:,:,15)=u(:,:,:,15)*w(:,:,:,0)
    w(:,:,:,16)=u(:,:,:,16)*w(:,:,:,0)
    w(:,:,:,21 )=u(:,:,:,21 ) * dv
    w(:,:,:,17) = u(:,:,:,17) * dv
    w(:,:,:,18) = u(:,:,:,18) * u(:,:,:,0) * dv
    w(:,:,:,19) = u(:,:,:,19) * u(:,:,:,0) * dv
    w(:,:,:,20) = u(:,:,:,20) * u(:,:,:,0) * dv
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             xmu = MP_mu / &
                  (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
             T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
             T_K = min(max(T_K,MP_Tmin),MP_Tmax) 
             gam = 1.d0+(1.d0+4.d0*yHe) & 
                  /(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
             w(i,j,k,4) =dv*(u(i,j,k,4)/(gam-1.d0) &
                  +u(i,j,k,0)*(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/2.d0)
          end do
       end do
    end do
    w(:,:,:,22) = u(:,:,:,22) *dv
    w(:,:,:,23) = u(:,:,:,23)*dv
    w(:,:,:,24) = u(:,:,:,24)*dv
    w(:,:,:,25) = u(:,:,:,25)*dv
    w(:,:,:,26) = u(:,:,:,26) *dv
    w(:,:,:,27) = u(:,:,:,27)*dv
    w(:,:,:,28) = u(:,:,:,28)*dv
    w(:,:,:,29) = u(:,:,:,29)*dv
    w(:,:,:,30) = u(:,:,:,30) *dv
    w(:,:,:,31) = u(:,:,:,31)*dv
    w(:,:,:,32) = u(:,:,:,32)*dv
    w(:,:,:,33) = u(:,:,:,33)*dv
  end subroutine u2w
  subroutine u2w_4(u,w,dv)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),intent(IN) :: dv
    real(kind=8) :: tap
    real(kind=8) :: T_K, xmu, gam
    integer :: i,j,k
    tap = MP_mu
    w(:,:,:,0)=u(:,:,:,0)*dv
    w(:,:,:,1) =u(:,:,:,1) *w(:,:,:,0)
    w(:,:,:,2) =u(:,:,:,2) *w(:,:,:,0)
    w(:,:,:,3) =u(:,:,:,3) *w(:,:,:,0)
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             xmu = MP_mu / &
                  (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
             T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
             T_K = min(max(T_K,MP_Tmin),MP_Tmax) 
             gam = 1.d0+tap & 
                  /(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
             w(i,j,k,4) =dv*(u(i,j,k,4)/(gam-1.d0) &
                  +u(i,j,k,0)*(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/2.d0)
          end do
       end do
    end do
  end subroutine u2w_4
  subroutine w2u(w,u,dv)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),intent(IN) :: dv
    real(kind=8) :: T_K, xmu, gam, p_gam_1, gam0
    integer :: i,j,k
    integer :: itr, n_itr=30
    real(kind=8) :: err_gam, max_err_gam=1d-10 
    real(kind=8) :: f1, f2, gam1,dif
    integer :: notConverge
    integer :: ichem
    u(:,:,:,0) =w(:,:,:,0)/dv
    u(:,:,:,1) =w(:,:,:,1) /w(:,:,:,0)
    u(:,:,:,2) =w(:,:,:,2) /w(:,:,:,0)
    u(:,:,:,3) =w(:,:,:,3) /w(:,:,:,0)
    do ichem = 5, 10
      u(:,:,:,ichem)=w(:,:,:,ichem)/w(:,:,:,0)
    enddo
    u(:,:,:,11)=w(:,:,:,11)/w(:,:,:,0)
    u(:,:,:,12)=w(:,:,:,12) / dv
    u(:,:,:,13)=w(:,:,:,13) / dv
    u(:,:,:,14)=w(:,:,:,14) / w(:,:,:,0)
    u(:,:,:,15)=w(:,:,:,15) / w(:,:,:,0)
    u(:,:,:,16)=w(:,:,:,16) / w(:,:,:,0)
    u(:,:,:,21 )=w(:,:,:,21 ) / dv
    u(:,:,:,17) = w(:,:,:,17) / dv
    u(:,:,:,18) = w(:,:,:,18) / w(:,:,:,0)
    u(:,:,:,19) = w(:,:,:,19) / w(:,:,:,0)
    u(:,:,:,20) = w(:,:,:,20) / w(:,:,:,0)
    gam = 5d0/3d0 
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             p_gam_1 = w(i,j,k,4)/dv &
                  - 0.5d0*u(i,j,k,0)*(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2) 
             xmu = MP_mu / &
                  (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
             u(i,j,k,4) = p_gam_1*(gam-1) 
             if ( .not. (u(i,j,k,0) > 0d0 .and. u(i,j,k,4) > 0d0 .and. &
                  u(i,j,k,5) > 0d0 .and. u(i,j,k,6) > 0d0 .and. &
                  u(i,j,k,7) > 0d0 .and. u(i,j,k,8) > 0d0 .and. &
                  u(i,j,k,5) < 2d0 .and. u(i,j,k,6) < 1d0 .and. &
                  u(i,j,k,7) < 2d0 .and. u(i,j,k,8) < 2d0 )) cycle
             notConverge = 1
             do itr = 0, n_itr-1
                gam0 = gam
                T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
                gam = 1.d0+(1.d0+4.d0*yHe)/(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
                f1 = gam - gam0
                err_gam = abs(f1)/max(gam,gam0)
                if (err_gam < max_err_gam) then
                  notConverge = 0
                  exit
                endif
                gam = gam0*(1.d0+1d-5)
                gam1= gam
                u(i,j,k,4) = p_gam_1*(gam-1)
                T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb
                gam = 1.d0+(1.d0+4.d0*yHe)/(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
                f2 = gam - gam1
                dif = (f2-f1)/(gam0*1d-5)
                gam = gam0 - f1/dif
                u(i,j,k,4) = p_gam_1*(gam-1)
                if (itr == n_itr-1) then
                   exit
                end if
             enddo
             if(notConverge == 1) then
                gam = 5d0/3d0 
                xmu = MP_mu / &
                    (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
                u(i,j,k,4) = p_gam_1*(gam-1) 
                do itr = 0, n_itr-1 
                   T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
                   gam0 = gam
                   gam = 1.d0+(1.d0+4.d0*yHe) & 
                        /(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2_2(T_K)*u(i,j,k,6)))
                   err_gam = abs(gam-gam0)/max(gam,gam0)
                   u(i,j,k,4) = p_gam_1*(gam-1) 
                   if (err_gam < max_err_gam) exit 
                   if (itr == n_itr-1) then
                      print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam&
&0, err_gam
                      gam = 5d0/3d0 
                      p_gam_1 = (cgs_kb*MP_Tmin)*(u(i,j,k,0)*Unit_rho)/(cgs_amu*xmu) / Unit_e
                      u(i,j,k,4) = p_gam_1*(gam-1)
                   end if
                end do
             endif
          end do
       end do
    end do
    u(:,:,:,22) = w(:,:,:,22) / dv
    u(:,:,:,23) = w(:,:,:,23) / dv
    u(:,:,:,24) = w(:,:,:,24) / dv
    u(:,:,:,25) = w(:,:,:,25) / dv
    u(:,:,:,26) = w(:,:,:,26) / dv
    u(:,:,:,27) = w(:,:,:,27)/ dv
    u(:,:,:,28) = w(:,:,:,28)/ dv
    u(:,:,:,29) = w(:,:,:,29)/ dv
    u(:,:,:,30) = w(:,:,:,30) / dv
    u(:,:,:,31) = w(:,:,:,31) / dv
    u(:,:,:,32) = w(:,:,:,32) / dv
    u(:,:,:,33) = w(:,:,:,33) / dv
  end subroutine w2u
  subroutine w2u_4(w,u,dv)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),intent(IN) :: dv
    real(kind=8) :: T_K, xmu, gam, p_gam_1, gam0
    integer :: i,j,k
    integer :: itr, n_itr=30
    real(kind=8) :: err_gam, max_err_gam=1d-10 
    real(kind=8) :: f1, f2, gam1,dif
    integer :: notConverge
    real(kind=8) :: tap
    tap = MP_mu
    u(:,:,:,1) =w(:,:,:,1)/w(:,:,:,0)
    u(:,:,:,2) =w(:,:,:,2)/w(:,:,:,0)
    u(:,:,:,3) =w(:,:,:,3)/w(:,:,:,0)
    gam = 5d0/3d0 
    do i = lbound(u,1),ubound(u,1)
       do j = lbound(u,2),ubound(u,2)
          do k = lbound(u,3),ubound(u,3)
             p_gam_1 = w(i,j,k,4)/dv &
                  - 0.5d0*u(i,j,k,0)*(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2) 
             xmu = MP_mu / &
                    (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
             u(i,j,k,4) = p_gam_1*(gam-1) 
             if ( .not. (u(i,j,k,0) > 0d0 .and. u(i,j,k,4) > 0d0 .and. &
                  u(i,j,k,5) > 0d0 .and. u(i,j,k,6) > 0d0 .and. &
                  u(i,j,k,7) > 0d0 .and. u(i,j,k,8) > 0d0 .and. &
                  u(i,j,k,5) < 2d0 .and. u(i,j,k,6) < 1d0 .and. &
                  u(i,j,k,7) < 2d0 .and. u(i,j,k,8) < 2d0 )) cycle
             notConverge = 1
             do itr = 0, n_itr-1
                gam0 = gam
                T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
                gam = 1.d0+tap/(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
                f1 = gam - gam0
                err_gam = abs(f1)/max(gam,gam0)
                if (err_gam < max_err_gam) then
                  notConverge = 0
                  exit
                endif
                gam = gam0*(1.d0+1d-5)
                gam1= gam
                u(i,j,k,4) = p_gam_1*(gam-1)
                T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb
                gam = 1.d0+tap/(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2(T_K)*u(i,j,k,6)))
                f2 = gam - gam1
                dif = (f2-f1)/(gam0*1d-5)
                gam = gam0 - f1/dif
                u(i,j,k,4) = p_gam_1*(gam-1)
                if (itr == n_itr-1) then
                   exit
                end if
             enddo
             if(notConverge == 1) then
                gam = 5d0/3d0 
                xmu = MP_mu / &
                    (u(i,j,k,5)+u(i,j,k,6)+u(i,j,k,7)+u(i,j,k,8)+yHe) 
                u(i,j,k,4) = p_gam_1*(gam-1) 
                do itr = 0, n_itr-1 
                   T_K = u(i,j,k,4)*Unit_e*cgs_amu*xmu /(u(i,j,k,0)*Unit_rho)/cgs_kb 
                   gam0 = gam
                   gam = 1.d0+tap & 
                        /(xmu*(1.5d0*(u(i,j,k,5)+u(i,j,k,7)+u(i,j,k,8)+yHe) + c_H2_2(T_K)*u(i,j,k,6)))
                   err_gam = abs(gam-gam0)/max(gam,gam0)
                   u(i,j,k,4) = p_gam_1*(gam-1) 
                   if (err_gam < max_err_gam) exit 
                   if (itr == n_itr-1) then
                      print *, "*** WARNING *** w2u: itr reaches n_itr, iteration to find consistent gamma not converged", gam, gam&
&0, err_gam
                      gam = 5d0/3d0 
                      p_gam_1 = (cgs_kb*MP_Tmin)*(u(i,j,k,0)*Unit_rho)/(cgs_amu*xmu) / Unit_e
                      u(i,j,k,4) = p_gam_1*(gam-1)
                   end if
                end do
             endif
          end do
       end do
    end do
  end subroutine w2u_4
  subroutine u2w_withgam(u,w,dv,gam)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),dimension(:,:,:),intent(IN) :: gam
    real(kind=8),intent(IN) :: dv
    real(kind=8) :: T_K, xmu
    integer :: i,j,k, ichem
    w(:,:,:,0)=u(:,:,:,0)*dv
    w(:,:,:,1)=u(:,:,:,1)*w(:,:,:,0)
    w(:,:,:,2)=u(:,:,:,2)*w(:,:,:,0)
    w(:,:,:,3)=u(:,:,:,3)*w(:,:,:,0)
    do ichem = 5, 10
      w(:,:,:,ichem)=u(:,:,:,ichem)*w(:,:,:,0)
    enddo
    w(:,:,:,11)=u(:,:,:,11)*w(:,:,:,0)
    w(:,:,:,12)=u(:,:,:,12) * dv
    w(:,:,:,13)=u(:,:,:,13) * dv
    w(:,:,:,14)=u(:,:,:,14) * w(:,:,:,0)
    w(:,:,:,15)=u(:,:,:,15) * w(:,:,:,0)
    w(:,:,:,16)=u(:,:,:,16) * w(:,:,:,0)
    w(:,:,:,21 )=u(:,:,:,21 ) * dv
    w(:,:,:,17) = u(:,:,:,17) * dv
    w(:,:,:,18) = u(:,:,:,18) * u(:,:,:,0) * dv
    w(:,:,:,19) = u(:,:,:,19) * u(:,:,:,0) * dv
    w(:,:,:,20) = u(:,:,:,20) * u(:,:,:,0) * dv
    w(:,:,:,4) =dv*(u(:,:,:,4)/(gam(:,:,:)-1.d0) &
                  +u(:,:,:,0)*(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2)/2.d0)
    w(:,:,:,22) = u(:,:,:,22) *dv
    w(:,:,:,23) = u(:,:,:,23)*dv
    w(:,:,:,24) = u(:,:,:,24)*dv
    w(:,:,:,25) = u(:,:,:,25)*dv
    w(:,:,:,26) = u(:,:,:,26) *dv
    w(:,:,:,27) = u(:,:,:,27)*dv
    w(:,:,:,28) = u(:,:,:,28)*dv
    w(:,:,:,29) = u(:,:,:,29)*dv
    w(:,:,:,30) = u(:,:,:,30) *dv
    w(:,:,:,31) = u(:,:,:,31)*dv
    w(:,:,:,32) = u(:,:,:,32)*dv
    w(:,:,:,33) = u(:,:,:,33)*dv
  end subroutine u2w_withgam
  subroutine w2u_withgam(w,u,dv,gam)
    use grid, only : Mmin, globdbg_myrank, globdbg_mygid, globdbg_rank, globdbg_gid, globdbg_i, globdbg_j, globdbg_k 
    use unit
    use primordial
    real(kind=8),dimension(:,:,:,Mmin:) :: w 
    real(kind=8),dimension(:,:,:,Mmin:) :: u 
    real(kind=8),dimension(:,:,:),intent(IN) :: gam
    real(kind=8),intent(IN) :: dv
    real(kind=8) :: T_K, xmu
    integer :: i,j,k,ichem
    u(:,:,:,0) =w(:,:,:,0)/dv
    u(:,:,:,1) =w(:,:,:,1)/w(:,:,:,0)
    u(:,:,:,2) =w(:,:,:,2)/w(:,:,:,0)
    u(:,:,:,3) =w(:,:,:,3)/w(:,:,:,0)
    do ichem = 5, 10
      u(:,:,:,ichem)=w(:,:,:,ichem)/w(:,:,:,0)
    enddo
    u(:,:,:,11)=w(:,:,:,11)/w(:,:,:,0)
    u(:,:,:,12)=w(:,:,:,12) / dv
    u(:,:,:,13)=w(:,:,:,13) / dv
    u(:,:,:,14)=w(:,:,:,14) / w(:,:,:,0)
    u(:,:,:,15)=w(:,:,:,15) / w(:,:,:,0)
    u(:,:,:,16)=w(:,:,:,16) / w(:,:,:,0)
    u(:,:,:,21 )=w(:,:,:,21 ) / dv
    u(:,:,:,17) = w(:,:,:,17) / dv
    u(:,:,:,18) = w(:,:,:,18) / w(:,:,:,0)
    u(:,:,:,19) = w(:,:,:,19) / w(:,:,:,0)
    u(:,:,:,20) = w(:,:,:,20) / w(:,:,:,0)
    u(:,:,:,4) = ( w(:,:,:,4)/dv -0.5d0*u(:,:,:,0) &
         *(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2) )*(gam(:,:,:)-1.0d0)
    u(:,:,:,22) = w(:,:,:,22) / dv
    u(:,:,:,23) = w(:,:,:,23)/ dv
    u(:,:,:,24) = w(:,:,:,24)/ dv
    u(:,:,:,25) = w(:,:,:,25)/ dv
    u(:,:,:,26) = w(:,:,:,26) / dv
    u(:,:,:,27) = w(:,:,:,27)/ dv
    u(:,:,:,28) = w(:,:,:,28)/ dv
    u(:,:,:,29) = w(:,:,:,29)/ dv
    u(:,:,:,30) = w(:,:,:,30) / dv
    u(:,:,:,31) = w(:,:,:,31)/ dv
    u(:,:,:,32) = w(:,:,:,32)/ dv
    u(:,:,:,33) = w(:,:,:,33)/ dv
  end subroutine w2u_withgam
  subroutine conv_u2w(uw, dv)
    real(kind=8),dimension(:,:,:,:) :: uw 
    real(kind=8),intent(IN) :: dv
    real(kind=8),dimension(lbound(uw,1):ubound(uw,1),lbound(uw,2):ubound(uw,2),lbound(uw,3):ubound(uw,3),lbound(uw,4):ubound(uw,4))&
& :: swap
    call u2w(uw,swap,dv)
    uw = swap
  end subroutine conv_u2w
  subroutine conv_w2u(wu, dv)
    real(kind=8),dimension(:,:,:,:) :: wu 
    real(kind=8),intent(IN) :: dv
    real(kind=8),dimension(lbound(wu,1):ubound(wu,1),lbound(wu,2):ubound(wu,2),lbound(wu,3):ubound(wu,3),lbound(wu,4):ubound(wu,4))&
& :: swap
    call w2u(wu,swap,dv)
    wu = swap
  end subroutine conv_w2u
  subroutine get_gamma(gid, gam)
    use grid, only : Imingh,Jmingh,Kmingh,Imaxgh,Jmaxgh,Kmaxgh, get_Ucomp
    use unit
    use primordial
    integer,intent(IN) :: gid
    real(kind=8),dimension(Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh),intent(OUT) :: gam
    real(kind=8),dimension(:,:,:),pointer :: rho, p, yhn, yh2, yel, yhp
    real(kind=8) :: xmu, T_K
    integer :: i,j,k
    rho => get_Ucomp(0,gid)
    p => get_Ucomp(4,gid)
    yhn => get_Ucomp(5,gid)
    yh2 => get_Ucomp(6,gid)
    yel => get_Ucomp(7,gid)
    yhp => get_Ucomp(8,gid)
    do i = Imingh, Imaxgh
       do j = Jmingh, Jmaxgh
          do k = Kmingh, Kmaxgh
             xmu = MP_mu /(yhn(i,j,k)+yh2(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe) 
             T_K = p(i,j,k)*Unit_e*cgs_amu*xmu /(rho(i,j,k)*Unit_rho)/cgs_kb 
             T_K = min(max(T_K,MP_Tmin),MP_Tmax) 
             gam(i,j,k) = 1.d0+(1.d0+4.d0*yHe) & 
                  /(xmu*(1.5d0*(yhn(i,j,k)+yel(i,j,k)+yhp(i,j,k)+yHe) + c_H2(T_K)*yh2(i,j,k)))
          end do
       end do
    end do
  end subroutine get_gamma
end module eos
subroutine eos_init
end subroutine eos_init
