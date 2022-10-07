module interpolation
  implicit none
  private
  public :: interp_tricubic, interp_trilinear, interp_bicubic_xy, interp_bicubic_yz, interp_bicubic_zx
contains
  subroutine interp_bicubic_xy(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0)
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc, bufx, bufxy
    integer,intent(IN) :: ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0
    call interp_tricubic(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0, boolz=.FALSE.)
  end subroutine interp_bicubic_xy
  subroutine interp_bicubic_yz(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0)
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc, bufx, bufxy
    integer,intent(IN) :: ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0
    call interp_tricubic(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0, boolx=.FALSE.)
  end subroutine interp_bicubic_yz
  subroutine interp_bicubic_zx(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0)
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc, bufx, bufxy
    integer,intent(IN) :: ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0
    call interp_tricubic(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0, booly=.FALSE.)
  end subroutine interp_bicubic_zx
  subroutine interp_tricubic(uf, uc, bufx, bufxy, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0, boolx, booly, boolz)
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc, bufx, bufxy
    integer,intent(IN) :: ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0
    logical,intent(IN),optional :: boolx, booly, boolz
    logical :: bool_x, bool_y, bool_z
    integer :: MINGHCELL = 1 
    integer :: MAXGHCELL = 2 
    integer :: if, jf, kf, ic, jc, kc, m, nshift
    integer :: ic0, ic1, ic2, ic3, icP, icN
    integer :: jc0, jc1, jc2, jc3, jcP, jcN
    integer :: kc0, kc1, kc2, kc3, kcP, kcN
    integer :: ics, jcs, kcs, ice, jce, kce
    integer :: ms, me
    integer,parameter :: MINBLOCKSIZE = 4
    real(kind=8) :: xc, yc, zc, a0, a1, a2, a3
    bool_x = .TRUE.
 bool_y = .TRUE.
 bool_z = .TRUE. 
    if (present(boolx)) bool_x = boolx
    if (present(booly)) bool_y = booly
    if (present(boolz)) bool_z = boolz
    ics = ((ifs)-(if0))/2 + mod(min((ifs)-(if0),0),2) + (if0)
    jcs = ((jfs)-(jf0))/2 + mod(min((jfs)-(jf0),0),2) + (jf0)
    kcs = ((kfs)-(kf0))/2 + mod(min((kfs)-(kf0),0),2) + (kf0)
    ice = ((ife)-(if0))/2 + mod(min((ife)-(if0),0),2) + (if0)
    jce = ((jfe)-(jf0))/2 + mod(min((jfe)-(jf0),0),2) + (jf0)
    kce = ((kfe)-(kf0))/2 + mod(min((kfe)-(kf0),0),2) + (kf0)
    if ( &
         ics-MINGHCELL < lbound(uc,1) .or. &
         jcs-MINGHCELL < lbound(uc,2) .or. &
         kcs-MINGHCELL < lbound(uc,3) .or. &
         ice+MINGHCELL > ubound(uc,1) .or. &
         jce+MINGHCELL > ubound(uc,2) .or. &
         kce+MINGHCELL > ubound(uc,3) ) then
       print *, '*** interp_tricubic: too samll array'
       stop
    end if
    ics = max(ics-MAXGHCELL, lbound(uc,1))
    jcs = max(jcs-MAXGHCELL, lbound(uc,2))
    kcs = max(kcs-MAXGHCELL, lbound(uc,3))
    ice = min(ice+MAXGHCELL, ubound(uc,1))
    jce = min(jce+MAXGHCELL, ubound(uc,2))
    kce = min(kce+MAXGHCELL, ubound(uc,3))
    ms = lbound(uc,4)
    me = ubound(uc,4)
    if ( & 
         ice-ics+1 < MINBLOCKSIZE .or. &
         jce-jcs+1 < MINBLOCKSIZE .or. &
         kce-kcs+1 < MINBLOCKSIZE ) then
       call interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0)
       return
    end if
    do m = ms, me
       if (.not. bool_x) exit
       do kc = kcs, kce
          do jc = jcs, jce
             do if = ifs, ife
                icP = ((if)-(if0))/2 + mod(min((if)-(if0),0),2) + (if0)
                icN = ( (icP)+2*modulo((if),2)-1 )
                ic1 = min(icP, icN)
                ic2 = max(icP, icN)
                ic0 = 2*ic1 - ic2
                ic3 = 2*ic2 - ic1
                nshift = 0
                if (ic3 > ubound(uc,1)) nshift = ic3 - ubound(uc,1)
                if (ic0 < lbound(uc,1)) nshift = ic0 - lbound(uc,1)
                ic0 = ic0 - nshift
                ic1 = ic1 - nshift
                ic2 = ic2 - nshift
                ic3 = ic3 - nshift
                xc = icP +0.25d0*(2*modulo(if,2)-1)
                a0 = (xc - ic1)*(xc - ic2)*(xc - ic3)/((ic0 - ic1)*(ic0 - ic2)*(ic0 - ic3))
                a1 = (xc - ic0)*(xc - ic2)*(xc - ic3)/((ic1 - ic0)*(ic1 - ic2)*(ic1 - ic3))
                a2 = (xc - ic0)*(xc - ic1)*(xc - ic3)/((ic2 - ic0)*(ic2 - ic1)*(ic2 - ic3))
                a3 = (xc - ic0)*(xc - ic1)*(xc - ic2)/((ic3 - ic0)*(ic3 - ic1)*(ic3 - ic2))
                bufx(if, jc, kc, m) = &
                     a0*uc(ic0, jc, kc, m) + &
                     a1*uc(ic1, jc, kc, m) + &
                     a2*uc(ic2, jc, kc, m) + &
                     a3*uc(ic3, jc, kc, m)
             end do
          end do
       end do
    end do
    do m = ms, me
       if (.not. bool_y) exit
       do kc = kcs, kce
          do jf = jfs, jfe
             jcP = ((jf)-(jf0))/2 + mod(min((jf)-(jf0),0),2) + (jf0)
             jcN = ( (jcP)+2*modulo((jf),2)-1 )
             jc1 = min(jcP, jcN)
             jc2 = max(jcP, jcN)
             jc0 = 2*jc1 - jc2
             jc3 = 2*jc2 - jc1
             nshift = 0
             if (jc3 > ubound(uc,2)) nshift = jc3 - ubound(uc,2)
             if (jc0 < lbound(uc,2)) nshift = jc0 - lbound(uc,2)
             jc0 = jc0 - nshift
             jc1 = jc1 - nshift
             jc2 = jc2 - nshift
             jc3 = jc3 - nshift
             yc = jcP +0.25d0*(2*modulo(jf,2)-1)
             a0 = (yc - jc1)*(yc - jc2)*(yc - jc3)/((jc0 - jc1)*(jc0 - jc2)*(jc0 - jc3))
             a1 = (yc - jc0)*(yc - jc2)*(yc - jc3)/((jc1 - jc0)*(jc1 - jc2)*(jc1 - jc3))
             a2 = (yc - jc0)*(yc - jc1)*(yc - jc3)/((jc2 - jc0)*(jc2 - jc1)*(jc2 - jc3))
             a3 = (yc - jc0)*(yc - jc1)*(yc - jc2)/((jc3 - jc0)*(jc3 - jc1)*(jc3 - jc2))
             do if = ifs, ife
                bufxy(if, jf, kc, m) = &
                     a0*bufx(if, jc0, kc, m) + &
                     a1*bufx(if, jc1, kc, m) + &
                     a2*bufx(if, jc2, kc, m) + &
                     a3*bufx(if, jc3, kc, m)
             end do
          end do
       end do
    end do
    do m = ms, me
       if (.not. bool_z) exit
       do kf = kfs, kfe
          kcP = ((kf)-(kf0))/2 + mod(min((kf)-(kf0),0),2) + (kf0)
          kcN = ( (kcP)+2*modulo((kf),2)-1 )
          kc1 = min(kcP, kcN)
          kc2 = max(kcP, kcN)
          kc0 = 2*kc1 - kc2
          kc3 = 2*kc2 - kc1
          nshift = 0
          if (kc3 > ubound(uc,3)) nshift = kc3 - ubound(uc,3)
          if (kc0 < lbound(uc,3)) nshift = kc0 - lbound(uc,3)
          kc0 = kc0 - nshift
          kc1 = kc1 - nshift
          kc2 = kc2 - nshift
          kc3 = kc3 - nshift
          zc = kcP +0.25d0*(2*modulo(kf,2)-1)
          a0 = (zc - kc1)*(zc - kc2)*(zc - kc3)/((kc0 - kc1)*(kc0 - kc2)*(kc0 - kc3))
          a1 = (zc - kc0)*(zc - kc2)*(zc - kc3)/((kc1 - kc0)*(kc1 - kc2)*(kc1 - kc3))
          a2 = (zc - kc0)*(zc - kc1)*(zc - kc3)/((kc2 - kc0)*(kc2 - kc1)*(kc2 - kc3))
          a3 = (zc - kc0)*(zc - kc1)*(zc - kc2)/((kc3 - kc0)*(kc3 - kc1)*(kc3 - kc2))
          do jf = jfs, jfe
             do if = ifs, ife
                uf(if, jf, kf, m) = &
                     a0*bufxy(if, jf, kc0, m) + &
                     a1*bufxy(if, jf, kc1, m) + &
                     a2*bufxy(if, jf, kc2, m) + &
                     a3*bufxy(if, jf, kc3, m)
             end do
          end do
       end do
    end do
  end subroutine interp_tricubic
  subroutine interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0)
    real(kind=8),dimension(:,:,:,:),pointer :: uf, uc
    integer,intent(IN) :: ifs, jfs, kfs, ife, jfe, kfe, if0, jf0, kf0
    integer :: if, jf, kf, ic, jc, kc, ic1, jc1, kc1, m
    do m = lbound(uc,4), ubound(uc,4)
       do kf = kfs, kfe
          do jf = jfs, jfe
             do if = ifs, ife
                ic = ((if)-(if0))/2 + mod(min((if)-(if0),0),2) + (if0)
                jc = ((jf)-(jf0))/2 + mod(min((jf)-(jf0),0),2) + (jf0)
                kc = ((kf)-(kf0))/2 + mod(min((kf)-(kf0),0),2) + (kf0)
                ic1 = ( (ic)+2*modulo((if),2)-1 )
                jc1 = ( (jc)+2*modulo((jf),2)-1 )
                kc1 = ( (kc)+2*modulo((kf),2)-1 )
                uf(if,jf,kf,m) = ( &
                     27.d0*uc(ic,jc,kc,m) &
                     +9.0d0*(uc(ic1,jc,kc,m)+uc(ic,jc1,kc,m)+uc(ic,jc,kc1,m)) &
                     +3.0d0*(uc(ic,jc1,kc1,m)+uc(ic1,jc,kc1,m)+uc(ic1,jc1,kc,m)) &
                     +uc(ic1,jc1,kc1,m) &
                     )/64.d0
             end do
          enddo
       enddo
    enddo
  end subroutine interp_trilinear
end module interpolation
