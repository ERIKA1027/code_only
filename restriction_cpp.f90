module restriction
  implicit none
  private
  integer,parameter :: STENCIL_MIN = 0
  integer,parameter :: STENCIL_MAX = 3
  real(kind=8),save,dimension(STENCIL_MIN:STENCIL_MAX,STENCIL_MIN:STENCIL_MAX,STENCIL_MIN:STENCIL_MAX) :: St 
  logical,save :: Initialized = .FALSE.
  public :: rstrct
contains
  subroutine rstrct(uc, uf, ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0, cubic)
    integer,intent(IN) :: ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    logical,intent(IN),optional :: cubic
    logical :: bool_cubic
    bool_cubic = .FALSE. 
    if (present(cubic)) bool_cubic = cubic
    if (bool_cubic) then
       call rstrct_cubic(uc, uf, ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0)
    else
       call rstrct_linear(uc, uf, ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0)
    end if
  end subroutine rstrct
  subroutine rstrct_cubic_init
    integer,dimension(STENCIL_MIN:STENCIL_MAX) :: ii
    real(kind=8),dimension(STENCIL_MIN:STENCIL_MAX,1) :: s1d
    real(kind=8),dimension(STENCIL_MIN:STENCIL_MAX, STENCIL_MIN:STENCIL_MAX) :: s2d
    real(kind=8) :: ic
    integer :: k
    if (Initialized) return
    Initialized = .TRUE.
    ii = (/0, 1, 2, 3/)
    ic = 1.5d0
    s1d(:,1) = (/ & 
         (ic - ii(1))*(ic - ii(2))*(ic - ii(3))/((ii(0) - ii(1))*(ii(0) - ii(2))*(ii(0) - ii(3))), &
         (ic - ii(0))*(ic - ii(2))*(ic - ii(3))/((ii(1) - ii(0))*(ii(1) - ii(2))*(ii(1) - ii(3))), &
         (ic - ii(0))*(ic - ii(1))*(ic - ii(3))/((ii(2) - ii(0))*(ii(2) - ii(1))*(ii(2) - ii(3))), &
         (ic - ii(0))*(ic - ii(1))*(ic - ii(2))/((ii(3) - ii(0))*(ii(3) - ii(1))*(ii(3) - ii(2))) /)
    s2d = matmul(s1d, transpose(s1d))
    do k = STENCIL_MIN, STENCIL_MAX
       St(:,:,k) = s2d * s1d(k,1)
    end do
  end subroutine rstrct_cubic_init
  subroutine rstrct_cubic(uc, uf, ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0)
    integer,intent(IN) :: ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: ic, jc, kc, m, &
         i0, j0, k0, i1, j1, k1, i2, j2, k2, i3, j3, k3
    call rstrct_cubic_init
    do m = lbound(uc,4), ubound(uc,4)
       do kc = kcs, kce
          do jc = jcs, jce
             do ic = ics, ice
                i1 = ((ic) - (ic0)) * 2 + (ic0)
                j1 = ((jc) - (jc0)) * 2 + (jc0)
                k1 = ((kc) - (kc0)) * 2 + (kc0)
                i0 = i1-1
                j0 = j1-1
                k0 = k1-1
                i2 = i1+1
                j2 = j1+1
                k2 = k1+1
                i3 = i1+2
                j3 = j1+2
                k3 = k1+2
                uc(ic,jc,kc,m) = &
                     uf(i0,j0,k0,m)*St(0,0,0)+uf(i1,j0,k0,m)*St(1,0,0)+uf(i2,j0,k0,m)*St(2,0,0)+uf(i3,j0,k0,m)*St(3,0,0)+&
                     uf(i0,j1,k0,m)*St(0,1,0)+uf(i1,j1,k0,m)*St(1,1,0)+uf(i2,j1,k0,m)*St(2,1,0)+uf(i3,j1,k0,m)*St(3,1,0)+&
                     uf(i0,j2,k0,m)*St(0,2,0)+uf(i1,j2,k0,m)*St(1,2,0)+uf(i2,j2,k0,m)*St(2,2,0)+uf(i3,j2,k0,m)*St(3,2,0)+&
                     uf(i0,j3,k0,m)*St(0,3,0)+uf(i1,j3,k0,m)*St(1,3,0)+uf(i2,j3,k0,m)*St(2,3,0)+uf(i3,j3,k0,m)*St(3,3,0)+&
                     uf(i0,j0,k1,m)*St(0,0,1)+uf(i1,j0,k1,m)*St(1,0,1)+uf(i2,j0,k1,m)*St(2,0,1)+uf(i3,j0,k1,m)*St(3,0,1)+&
                     uf(i0,j1,k1,m)*St(0,1,1)+uf(i1,j1,k1,m)*St(1,1,1)+uf(i2,j1,k1,m)*St(2,1,1)+uf(i3,j1,k1,m)*St(3,1,1)+&
                     uf(i0,j2,k1,m)*St(0,2,1)+uf(i1,j2,k1,m)*St(1,2,1)+uf(i2,j2,k1,m)*St(2,2,1)+uf(i3,j2,k1,m)*St(3,2,1)+&
                     uf(i0,j3,k1,m)*St(0,3,1)+uf(i1,j3,k1,m)*St(1,3,1)+uf(i2,j3,k1,m)*St(2,3,1)+uf(i3,j3,k1,m)*St(3,3,1)+&
                     uf(i0,j0,k2,m)*St(0,0,2)+uf(i1,j0,k2,m)*St(1,0,2)+uf(i2,j0,k2,m)*St(2,0,2)+uf(i3,j0,k2,m)*St(3,0,2)+&
                     uf(i0,j1,k2,m)*St(0,1,2)+uf(i1,j1,k2,m)*St(1,1,2)+uf(i2,j1,k2,m)*St(2,1,2)+uf(i3,j1,k2,m)*St(3,1,2)+&
                     uf(i0,j2,k2,m)*St(0,2,2)+uf(i1,j2,k2,m)*St(1,2,2)+uf(i2,j2,k2,m)*St(2,2,2)+uf(i3,j2,k2,m)*St(3,2,2)+&
                     uf(i0,j3,k2,m)*St(0,3,2)+uf(i1,j3,k2,m)*St(1,3,2)+uf(i2,j3,k2,m)*St(2,3,2)+uf(i3,j3,k2,m)*St(3,3,2)+&
                     uf(i0,j0,k3,m)*St(0,0,3)+uf(i1,j0,k3,m)*St(1,0,3)+uf(i2,j0,k3,m)*St(2,0,3)+uf(i3,j0,k3,m)*St(3,0,3)+&
                     uf(i0,j1,k3,m)*St(0,1,3)+uf(i1,j1,k3,m)*St(1,1,3)+uf(i2,j1,k3,m)*St(2,1,3)+uf(i3,j1,k3,m)*St(3,1,3)+&
                     uf(i0,j2,k3,m)*St(0,2,3)+uf(i1,j2,k3,m)*St(1,2,3)+uf(i2,j2,k3,m)*St(2,2,3)+uf(i3,j2,k3,m)*St(3,2,3)+&
                     uf(i0,j3,k3,m)*St(0,3,3)+uf(i1,j3,k3,m)*St(1,3,3)+uf(i2,j3,k3,m)*St(2,3,3)+uf(i3,j3,k3,m)*St(3,3,3)
             end do
          enddo
       end do
    end do
  end subroutine rstrct_cubic
  subroutine rstrct_linear(uc, uf, ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0)
    integer,intent(IN) :: ics, jcs, kcs, ice, jce, kce, ic0, jc0, kc0
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    real(kind=8),parameter :: rdv = 1.d0/8.d0
    integer :: if, jf, kf, ic, jc, kc, m
    do m = lbound(uc,4), ubound(uc,4)
       do kc = kcs, kce
          do jc = jcs, jce
             do ic = ics, ice
                if = ((ic) - (ic0)) * 2 + (ic0)
                jf = ((jc) - (jc0)) * 2 + (jc0)
                kf = ((kc) - (kc0)) * 2 + (kc0)
                uc(ic,jc,kc,m) = &
                     (uf(if, jf, kf,m) +uf(if+1,jf, kf,m) &
                     +uf(if, jf+1,kf,m) +uf(if, jf, kf+1,m) &
                     +uf(if+1,jf+1,kf,m) +uf(if+1,jf, kf+1,m) &
                     +uf(if, jf+1,kf+1,m)+uf(if+1,jf+1,kf+1,m))*rdv
             end do
          enddo
       end do
    end do
  end subroutine rstrct_linear
end module restriction
