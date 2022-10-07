#include "config.h"
#include "chemistry_label.h"
 
 
 
!--------------------------------------------------
!            subroutine for react_rat
!--------------------------------------------------
  subroutine react_rat_highspeed2(xk, xnH, y, r_f_tot)
    implicit none
    real(kind=DBL_KIND),intent(IN) :: xk(0:NREACT-1),xnH,y(0:NCHEM-1)
    real(kind=DBL_KIND),intent(OUT) :: r_f_tot(0:NCHEM-1)
    integer :: isp, ire
    real (kind=DBL_KIND) :: yhmn, rate, rrate
    
    !initialization
    r_f_tot(:) = 0.d0
    
    ! CHEM_H0)   H  +  e   ->   H+  +  2e
    rate = xk(CHEM_H0)*y(X_HI)*y(X_EL)*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + rate


    ! CHEM_H1)   H+  +  e   ->   H
    rate = xk(CHEM_H1)*y(X_HII)*y(X_EL)*xnH
    r_f_tot(X_HII) = r_f_tot(X_HII) - rate
    r_f_tot(X_EL) = r_f_tot(X_EL) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + rate


    ! CHEM_H2)   H-  +  H   ->   H2  +  e
    rate = xk(CHEM_H2)*y(X_Hm)*y(X_HI)*xnH
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate


    ! CHEM_H3)   H2  +  H+   ->   H2+  +  H
    rate = xk(CHEM_H3)*y(X_H2)*y(X_HII)*xnH
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) - rate
    r_f_tot(X_H2p) = r_f_tot(X_H2p) + rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + rate


    ! CHEM_H4)   H2  +  e   ->   2H  +  e
    rate = xk(CHEM_H4)*y(X_H2)*y(X_EL)*xnH
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H5)   H2  +  H   ->   3H
    rate = xk(CHEM_H5)*y(X_H2)*y(X_HI)*xnH
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H6)   3H   ->   H2  +   H
    rate = xk(CHEM_H6)*y(X_HI)*y(X_HI)*y(X_HI)*xnH*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - 2.d0*rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate


    ! CHEM_H7)   2H  +  H2   ->   2H2
    rate = xk(CHEM_H7)*y(X_HI)*y(X_HI)*y(X_H2)*xnH*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - 2.d0*rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate


    ! CHEM_H8)   2H2   ->   2H  +  H2
    rate = xk(CHEM_H8)*y(X_H2)*y(X_H2)*xnH
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H9)   H  +  e   ->   H-
    rate = xk(CHEM_H9)*y(X_HI)*y(X_EL)*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_EL) = r_f_tot(X_EL) - rate
    r_f_tot(X_Hm) = r_f_tot(X_Hm) + rate


    ! CHEM_H10)   H  +  ph.  ->  H+  +  e
    rate = xk(CHEM_H10)*y(X_HI)
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate


    ! CHEM_H2PHI)   H2  +  ph.  ->  2H+  +  2e
    rate = xk(CHEM_H2PHI)*y(X_H2)
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + 2.d0*rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + 2.d0*rate


    ! CHEM_H11)   H2  +  ph.()  ->  2H
    rate = xk(CHEM_H11)*y(X_H2)
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H12)   2H   ->   H+  +  e  +  H
    rate = xk(CHEM_H12)*y(X_HI)*y(X_HI)*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate


    ! CHEM_H13)   H-  +  e   ->   H  +  2e
    rate = xk(CHEM_H13)*y(X_Hm)*y(X_EL)*xnH
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + rate


    ! CHEM_H14)   H-  +  H+   ->   H2+  +  e
    rate = xk(CHEM_H14)*y(X_Hm)*y(X_HII)*xnH
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) - rate
    r_f_tot(X_H2p) = r_f_tot(X_H2p) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate


    ! CHEM_H15)   H-  +  H+   ->   2H
    rate = xk(CHEM_H15)*y(X_Hm)*y(X_HII)*xnH
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H16)   H-  +  ph.  ->  H  +  e
    rate = xk(CHEM_H16)*y(X_Hm)
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate


    ! CHEM_H17)   H  +  H+   ->   H2+
    rate = xk(CHEM_H17)*y(X_HI)*y(X_HII)*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) - rate
    r_f_tot(X_H2p) = r_f_tot(X_H2p) + rate


    ! CHEM_H18)   H2+  +  H   ->   H2  +  H+
    rate = xk(CHEM_H18)*y(X_H2p)*y(X_HI)*xnH
    r_f_tot(X_H2p) = r_f_tot(X_H2p) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + rate


    ! CHEM_H19)   H2+  +  e   ->   2H
    rate = xk(CHEM_H19)*y(X_H2p)*y(X_EL)*xnH
    r_f_tot(X_H2p) = r_f_tot(X_H2p) - rate
    r_f_tot(X_EL) = r_f_tot(X_EL) - rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + 2.d0*rate


    ! CHEM_H20)   H2+  +  H-   ->   H2  +  H
    rate = xk(CHEM_H20)*y(X_H2p)*y(X_Hm)*xnH
    r_f_tot(X_H2p) = r_f_tot(X_H2p) - rate
    r_f_tot(X_Hm) = r_f_tot(X_Hm) - rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate
    r_f_tot(X_HI) = r_f_tot(X_HI) + rate


    ! CHEM_H21)   2H  +  grain  ->  H2
    rate = xk(CHEM_H21)*y(X_HI)*xnH
    r_f_tot(X_HI) = r_f_tot(X_HI) - 2.d0*rate
    r_f_tot(X_H2) = r_f_tot(X_H2) + rate


#ifdef INCLUDE_COSMICRAY
    ! CHEM_H22)   H  +  cr  ->  H+  +  e
    rate = xk(CHEM_H22)*y(X_HI)
    r_f_tot(X_HI) = r_f_tot(X_HI) - rate
    r_f_tot(X_HII) = r_f_tot(X_HII) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate
#endif


#ifdef INCLUDE_COSMICRAY
    ! CHEM_H23)   H2  +  cr  ->  H2+  +  e
    rate = xk(CHEM_H23)*y(X_H2)
    r_f_tot(X_H2) = r_f_tot(X_H2) - rate
    r_f_tot(X_H2p) = r_f_tot(X_H2p) + rate
    r_f_tot(X_EL) = r_f_tot(X_EL) + rate
#endif


 
 
  end subroutine react_rat_highspeed2
 
 
 
 
!--------------------------------------------------
!            subroutine for react_drdt
!--------------------------------------------------
  subroutine react_drdy(xk, xnH, y, dr_fdy)
    real(kind=DBL_KIND),intent(IN) :: xk(0:NREACT-1),xnH,y(0:NCHEM-1)
    real(kind=DBL_KIND),dimension(0:NCHEM-1,0:NCHEM-1),intent(OUT) :: dr_fdy
    real(kind=DBL_KIND) :: xnH_2, ddr_A, ddr_B, ddr_C
 
    dr_fdy(:,:) = 0.d0
    xnH_2 = xnH*xnH
 
 
    ! CHEM_H0)   H  +  e   ->   H+  +  2e
    ddr_A = xk(CHEM_H0)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H0)*y(X_HI)*xnH
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_A
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) + ddr_A
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) + ddr_A
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) - ddr_B
    dr_fdy(X_EL,X_EL) = dr_fdy(X_EL,X_EL) + ddr_B
    dr_fdy(X_HII,X_EL) = dr_fdy(X_HII,X_EL) + ddr_B


    ! CHEM_H1)   H+  +  e   ->   H
    ddr_A = xk(CHEM_H1)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H1)*y(X_HII)*xnH
    dr_fdy(X_HII,X_HII) = dr_fdy(X_HII,X_HII) - ddr_A
    dr_fdy(X_EL,X_HII) = dr_fdy(X_EL,X_HII) - ddr_A
    dr_fdy(X_HI,X_HII) = dr_fdy(X_HI,X_HII) + ddr_A
    dr_fdy(X_HII,X_EL) = dr_fdy(X_HII,X_EL) - ddr_B
    dr_fdy(X_EL,X_EL) = dr_fdy(X_EL,X_EL) - ddr_B
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) + ddr_B


    ! CHEM_H2)   H-  +  H   ->   H2  +  e
    ddr_A = xk(CHEM_H2)*y(X_HI)*xnH
    ddr_B = xk(CHEM_H2)*y(X_Hm)*xnH
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - ddr_A
    dr_fdy(X_HI,X_Hm) = dr_fdy(X_HI,X_Hm) - ddr_A
    dr_fdy(X_H2,X_Hm) = dr_fdy(X_H2,X_Hm) + ddr_A
    dr_fdy(X_EL,X_Hm) = dr_fdy(X_EL,X_Hm) + ddr_A
    dr_fdy(X_Hm,X_HI) = dr_fdy(X_Hm,X_HI) - ddr_B
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_B
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) + ddr_B
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) + ddr_B


    ! CHEM_H3)   H2  +  H+   ->   H2+  +  H
    ddr_A = xk(CHEM_H3)*y(X_HII)*xnH
    ddr_B = xk(CHEM_H3)*y(X_H2)*xnH
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - ddr_A
    dr_fdy(X_HII,X_H2) = dr_fdy(X_HII,X_H2) - ddr_A
    dr_fdy(X_H2p,X_H2) = dr_fdy(X_H2p,X_H2) + ddr_A
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) + ddr_A
    dr_fdy(X_H2,X_HII) = dr_fdy(X_H2,X_HII) - ddr_B
    dr_fdy(X_HII,X_HII) = dr_fdy(X_HII,X_HII) - ddr_B
    dr_fdy(X_H2p,X_HII) = dr_fdy(X_H2p,X_HII) + ddr_B
    dr_fdy(X_HI,X_HII) = dr_fdy(X_HI,X_HII) + ddr_B


    ! CHEM_H4)   H2  +  e   ->   2H  +  e
    ddr_A = xk(CHEM_H4)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H4)*y(X_H2)*xnH
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - ddr_A
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) + 2.d0*ddr_A
    dr_fdy(X_H2,X_EL) = dr_fdy(X_H2,X_EL) - ddr_B
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) + 2.d0*ddr_B


    ! CHEM_H5)   H2  +  H   ->   3H
    ddr_A = xk(CHEM_H5)*y(X_HI)*xnH
    ddr_B = xk(CHEM_H5)*y(X_H2)*xnH
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - ddr_A
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) + 2.d0*ddr_A
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) - ddr_B
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) + 2.d0*ddr_B


    ! CHEM_H6)   3H   ->   H2  +   H
    ddr_A = 3.d0*xk(CHEM_H6)*y(X_HI)*y(X_HI)*xnH_2
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - 2.d0*ddr_A
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) + ddr_A


    ! CHEM_H7)   2H  +  H2   ->   2H2
    ddr_A = 2.d0*xk(CHEM_H7)*y(X_H2)*y(X_HI)*xnH_2
    ddr_B =      xk(CHEM_H7)*y(X_HI)*y(X_HI)*xnH_2
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - 2.d0*ddr_A
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) + ddr_A
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) - 2.d0*ddr_B
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) + ddr_B


    ! CHEM_H8)   2H2   ->   2H  +  H2
    ddr_A = 2.d0*xk(CHEM_H8)*y(X_H2)*xnH
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - ddr_A
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) + 2.d0*ddr_A


    ! CHEM_H9)   H  +  e   ->   H-
    ddr_A = xk(CHEM_H9)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H9)*y(X_HI)*xnH
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_A
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) - ddr_A
    dr_fdy(X_Hm,X_HI) = dr_fdy(X_Hm,X_HI) + ddr_A
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) - ddr_B
    dr_fdy(X_EL,X_EL) = dr_fdy(X_EL,X_EL) - ddr_B
    dr_fdy(X_Hm,X_EL) = dr_fdy(X_Hm,X_EL) + ddr_B


    ! CHEM_H10)   H  +  ph.  ->  H+  +  e
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - xk(CHEM_H10)
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) + xk(CHEM_H10)
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) + xk(CHEM_H10)


    ! CHEM_H2PHI)   H2  +  ph.  ->  2H+  +  2e
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - xk(CHEM_H2PHI)
    dr_fdy(X_HII,X_H2) = dr_fdy(X_HII,X_H2) + 2.d0*xk(CHEM_H2PHI)
    dr_fdy(X_EL,X_H2) = dr_fdy(X_EL,X_H2) + 2.d0*xk(CHEM_H2PHI)


    ! CHEM_H11)   H2  +  ph.  ->  2H
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - xk(CHEM_H11)
    dr_fdy(X_HI,X_H2) = dr_fdy(X_HI,X_H2) + 2.d0*xk(CHEM_H11)


    ! CHEM_H12)   2H   ->   H+  +  e  +  H
    ddr_A = 2.d0*xk(CHEM_H12)*y(X_HI)*xnH
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_A
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) + ddr_A
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) + ddr_A


    ! CHEM_H13)   H-  +  e   ->   H  +  2e
    ddr_A = xk(CHEM_H13)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H13)*y(X_Hm)*xnH
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - ddr_A
    dr_fdy(X_EL,X_Hm) = dr_fdy(X_EL,X_Hm) + ddr_A
    dr_fdy(X_HI,X_Hm) = dr_fdy(X_HI,X_Hm) + ddr_A
    dr_fdy(X_Hm,X_EL) = dr_fdy(X_Hm,X_EL) - ddr_B
    dr_fdy(X_EL,X_EL) = dr_fdy(X_EL,X_EL) + ddr_B
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) + ddr_B


    ! CHEM_H14)   H-  +  H+   ->   H2+  +  e
    ddr_A = xk(CHEM_H14)*y(X_HII)*xnH
    ddr_B = xk(CHEM_H14)*y(X_Hm)*xnH
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - ddr_A
    dr_fdy(X_HII,X_Hm) = dr_fdy(X_HII,X_Hm) - ddr_A
    dr_fdy(X_H2p,X_Hm) = dr_fdy(X_H2p,X_Hm) + ddr_A
    dr_fdy(X_EL,X_Hm) = dr_fdy(X_EL,X_Hm) + ddr_A
    dr_fdy(X_Hm,X_HII) = dr_fdy(X_Hm,X_HII) - ddr_B
    dr_fdy(X_HII,X_HII) = dr_fdy(X_HII,X_HII) - ddr_B
    dr_fdy(X_H2p,X_HII) = dr_fdy(X_H2p,X_HII) + ddr_B
    dr_fdy(X_EL,X_HII) = dr_fdy(X_EL,X_HII) + ddr_B


    ! CHEM_H15)   H-  +  H+   ->   2H
    ddr_A = xk(CHEM_H15)*y(X_HII)*xnH
    ddr_B = xk(CHEM_H15)*y(X_Hm)*xnH
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - ddr_A
    dr_fdy(X_HII,X_Hm) = dr_fdy(X_HII,X_Hm) - ddr_A
    dr_fdy(X_HI,X_Hm) = dr_fdy(X_HI,X_Hm) + 2.d0*ddr_A
    dr_fdy(X_Hm,X_HII) = dr_fdy(X_Hm,X_HII) - ddr_B
    dr_fdy(X_HII,X_HII) = dr_fdy(X_HII,X_HII) - ddr_B
    dr_fdy(X_HI,X_HII) = dr_fdy(X_HI,X_HII) + 2.d0*ddr_B


    ! CHEM_H16)   H-  +  ph.  ->  H  +  e
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - xk(CHEM_H16)
    dr_fdy(X_HI,X_Hm) = dr_fdy(X_HI,X_Hm) + xk(CHEM_H16)
    dr_fdy(X_EL,X_Hm) = dr_fdy(X_EL,X_Hm) + xk(CHEM_H16)


    ! CHEM_H17)   H  +  H+   ->   H2+
    ddr_A = xk(CHEM_H17)*y(X_HII)*xnH
    ddr_B = xk(CHEM_H17)*y(X_HI)*xnH
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_A
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) - ddr_A
    dr_fdy(X_H2p,X_HI) = dr_fdy(X_H2p,X_HI) + ddr_A
    dr_fdy(X_HI,X_HII) = dr_fdy(X_HI,X_HII) - ddr_B
    dr_fdy(X_HII,X_HII) = dr_fdy(X_HII,X_HII) - ddr_B
    dr_fdy(X_H2p,X_HII) = dr_fdy(X_H2p,X_HII) + ddr_B


    ! CHEM_H18)   H2+  +  H   ->   H2  +  H+
    ddr_A = xk(CHEM_H18)*y(X_HI)*xnH
    ddr_B = xk(CHEM_H18)*y(X_H2p)*xnH
    dr_fdy(X_H2p,X_H2p) = dr_fdy(X_H2p,X_H2p) - ddr_A
    dr_fdy(X_HI,X_H2p) = dr_fdy(X_HI,X_H2p) - ddr_A
    dr_fdy(X_H2,X_H2p) = dr_fdy(X_H2,X_H2p) + ddr_A
    dr_fdy(X_HII,X_H2p) = dr_fdy(X_HII,X_H2p) + ddr_A
    dr_fdy(X_H2p,X_HI) = dr_fdy(X_H2p,X_HI) - ddr_B
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - ddr_B
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) + ddr_B
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) + ddr_B


    ! CHEM_H19)   H2+  +  e   ->   2H
    ddr_A = xk(CHEM_H19)*y(X_EL)*xnH
    ddr_B = xk(CHEM_H19)*y(X_H2p)*xnH
    dr_fdy(X_H2p,X_H2p) = dr_fdy(X_H2p,X_H2p) - ddr_A
    dr_fdy(X_EL,X_H2p) = dr_fdy(X_EL,X_H2p) - ddr_A
    dr_fdy(X_HI,X_H2p) = dr_fdy(X_HI,X_H2p) + 2.d0*ddr_A
    dr_fdy(X_H2p,X_EL) = dr_fdy(X_H2p,X_EL) - ddr_B
    dr_fdy(X_EL,X_EL) = dr_fdy(X_EL,X_EL) - ddr_B
    dr_fdy(X_HI,X_EL) = dr_fdy(X_HI,X_EL) + 2.d0*ddr_B


    ! CHEM_H20)   H2+  +  H-   ->   H2  +  H
    ddr_A = xk(CHEM_H20)*y(X_Hm)*xnH
    ddr_B = xk(CHEM_H20)*y(X_H2p)*xnH
    dr_fdy(X_H2p,X_H2p) = dr_fdy(X_H2p,X_H2p) - ddr_A
    dr_fdy(X_Hm,X_H2p) = dr_fdy(X_Hm,X_H2p) - ddr_A
    dr_fdy(X_H2,X_H2p) = dr_fdy(X_H2,X_H2p) + ddr_A
    dr_fdy(X_HI,X_H2p) = dr_fdy(X_HI,X_H2p) + ddr_A
    dr_fdy(X_H2p,X_Hm) = dr_fdy(X_H2p,X_Hm) - ddr_B
    dr_fdy(X_Hm,X_Hm) = dr_fdy(X_Hm,X_Hm) - ddr_B
    dr_fdy(X_H2,X_Hm) = dr_fdy(X_H2,X_Hm) + ddr_B
    dr_fdy(X_HI,X_Hm) = dr_fdy(X_HI,X_Hm) + ddr_B


    ! CHEM_H21)   2H  +  grain  ->  H2
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - 2.d0*xk(CHEM_H21)*xnH
    dr_fdy(X_H2,X_HI) = dr_fdy(X_H2,X_HI) + xk(CHEM_H21)*xnH


#ifdef INCLUDE_COSMICRAY
    ! CHEM_H22)   H  +  cr  ->  H+  +  e
    dr_fdy(X_HI,X_HI) = dr_fdy(X_HI,X_HI) - xk(CHEM_H22)
    dr_fdy(X_HII,X_HI) = dr_fdy(X_HII,X_HI) + xk(CHEM_H22)
    dr_fdy(X_EL,X_HI) = dr_fdy(X_EL,X_HI) + xk(CHEM_H22)
#endif


#ifdef INCLUDE_COSMICRAY
    ! CHEM_H23)   H2  +  cr  ->  H2+  +  e
    dr_fdy(X_H2,X_H2) = dr_fdy(X_H2,X_H2) - xk(CHEM_H23)
    dr_fdy(X_H2p,X_H2) = dr_fdy(X_H2p,X_H2) + xk(CHEM_H23)
    dr_fdy(X_EL,X_H2) = dr_fdy(X_EL,X_H2) + xk(CHEM_H23)
#endif


 
 
  end subroutine react_drdy
 
 
 
 
!--------------------------------------------------
!            subroutine for adujust abundance
!--------------------------------------------------
  subroutine adjust_abundance(y, yco, metal)
    implicit none
    real(kind=DBL_KIND),intent(IN) :: metal
    real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1), yco
    real(kind=DBL_KIND) :: ycII
    real(kind=DBL_KIND) ::  yHtot, yCtot, yOtot, yDtot, chrgtot, nchrg, pchrg
    real(kind=DBL_KIND), parameter :: min_value = 1d-20 !復帰する際の値
    real(kind=DBL_KIND) :: frac_C, frac_O
    
    ! set frac_C, frac_O
    frac_C = metal*MP_frac_C_solar
    frac_O = metal*MP_frac_O_solar
    
    
    y(X_HI) = MAX(MIN(y(X_HI), 1.d0), min_value)
    y(X_H2) = MAX(MIN(y(X_H2), 0.5d0), min_value)
    y(X_HII) = MAX(MIN(y(X_HII), 1.d0), min_value)
    y(X_Hm) = MAX(MIN(y(X_Hm), 1.d0), min_value)
    y(X_H2p) = MAX(MIN(y(X_H2p), 1.d0), min_value)
    y(X_EL) = MAX(min_value, y(X_EL))
    yco = MAX(MIN(yco, frac_C), min_value)
    ycII= MAX(frac_C-yco, min_value)
    
    
    !------------------------- total C-abundance ---------------------!
    !------------------------- total O-abundance ---------------------!
    
    
    !------------------------- total D-abundance ---------------------!
    
    
    !------------------------- total H-abundance ---------------------!
    yHtot = y(X_HI)+2.d0*y(X_H2)+y(X_HII)+y(X_Hm)+2.d0*y(X_H2p)
    if(abs(yHtot-1) > 1.d-10) then
        if(yHtot < 1.d0) then
            y(X_HI) = y(X_HI)+1.d0-yHtot
        else
           y(X_HI)=y(X_HI)/yHtot
           y(X_H2)=y(X_H2)/yHtot
           y(X_HII)=y(X_HII)/yHtot
           y(X_Hm)=y(X_Hm)/yHtot
           y(X_H2p)=y(X_H2p)/yHtot
        endif
    endif
    
    
    !------------------------- total charge ---------------------!
    pchrg = ycII +y(X_HII)+y(X_H2p)
    nchrg = y(X_EL)+y(X_Hm)
    
    
    chrgtot = 2 * (pchrg-nchrg) / (nchrg+pchrg)
    if(abs(chrgtot) > 1.d-10) then
        ! balancing charge with electron
        y(X_EL) = MAX(min_value, y(X_EL)+pchrg-nchrg)
    endif
    
    
  end subroutine adjust_abundance
    
    
!--------------------------------------------------
!            subroutine for setting y_max & y_min
!--------------------------------------------------
  subroutine setting_ymax_ymin(y_max, y_min, metal)
      implicit none
      real(kind=DBL_KIND),dimension(0:NCHEM-1) :: y_max, y_min
      real(kind=DBL_KIND), parameter :: min_value = 1d-20 !復帰する際の値
      real(kind=DBL_KIND) :: frac_C, frac_O, metal
    
      frac_C = metal*MP_frac_C_solar
      frac_O = metal*MP_frac_O_solar
    
    
      y_min(:) = min_value
    
      y_max(X_H2) = 0.5d0
      y_max(X_EL) = 3.d0 ! いい加減な値
      y_max(X_HI) = 1.d0
      y_max(X_HII) = 1.d0
      y_max(X_Hm) = 1.d0
      y_max(X_H2p) = 1.d0
    
    
  end subroutine setting_ymax_ymin
