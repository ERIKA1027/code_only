#include "config.h" 

#define USE_LAPACK_CHEM NO
! --------------------------------------------------------------
!  additional module for line cooling
! --------------------------------------------------------------
module kinzoku2

  implicit none
  private

  ! line3
  type line3
    real(kind=DBL_KIND) :: xnu_10
    real(kind=DBL_KIND) :: xnu_20 
    real(kind=DBL_KIND) :: xnu_21
    real(kind=DBL_KIND) :: Q_10
    real(kind=DBL_KIND) :: Q_20
    real(kind=DBL_KIND) :: Q_21
    real(kind=DBL_KIND) :: A_10 
    real(kind=DBL_KIND) :: A_20
    real(kind=DBL_KIND) :: A_21
    real(kind=DBL_KIND) :: C_10
    real(kind=DBL_KIND) :: C_20
    real(kind=DBL_KIND) :: C_21
    real(kind=DBL_KIND) :: C_01
    real(kind=DBL_KIND) :: C_02
    real(kind=DBL_KIND) :: C_12
    real(kind=DBL_KIND) :: g_0
    real(kind=DBL_KIND) :: g_1
    real(kind=DBL_KIND) :: g_2
    real(kind=DBL_KIND) :: xNclmn
    real(kind=DBL_KIND) :: v_th
    real(kind=DBL_KIND) :: tau_cnt
  end type line3

  ! line2
  type line2
    real(kind=DBL_KIND) :: Q_10
    real(kind=DBL_KIND) :: A_10
    real(kind=DBL_KIND) :: C_10
    real(kind=DBL_KIND) :: C_01
    real(kind=DBL_KIND) :: g_0
    real(kind=DBL_KIND) :: g_1
    real(kind=DBL_KIND) :: tau0
    real(kind=DBL_KIND) :: tau_cnt
  end type line2

  public :: CIcool, CIIcool, HDcool, OHcool, H2Ocool, OIcool

contains

  ! -------------------------------------------------------
  !                   CII cooling
  ! -------------------------------------------------------
  subroutine CIIcool(xnH,T_K,Trad,xNc_CII,y_m,y_a,y_e,esc_10,tau_c,xLd_CII)
    implicit none 
    real(kind=DBL_KIND) :: xnH,T_K,Trad,xNc_CII,y_m,y_a,y_e,esc_10,tau_c,xLd_CII
    real(kind=DBL_KIND), parameter :: xk_B = 1.38d-16
    real(kind=DBL_KIND), parameter :: h_P  = 6.63d-27
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24
    real(kind=DBL_KIND), parameter :: pi   = 3.14159265358979d0
    type(line2) :: l2
    real(kind=DBL_KIND) :: DT_10, DE_10, xnu_10, gamma_e, gamma_H, gamma_H2
    real(kind=DBL_KIND) :: xm_C, v_th, f_0, f_1,S_10

    l2%tau_cnt=tau_c

!   set coefficients
    l2%g_0=2.d0
    l2%g_1=4.d0
!
    DT_10=92.d0
    DE_10=DT_10*xk_B
    xnu_10=DE_10/h_P
!
    l2%Q_10=Q_bg(DT_10, Trad)
!
    l2%A_10=2.4d-6
    gamma_e=2.8d-7/dsqrt(T_K*1.d-2)
    gamma_H=8.0d-10*(T_K*1.d-2)**7.d-2
    gamma_H2=5.d-1*gamma_H
    l2%C_10=xnH*(y_e*gamma_e+y_a*gamma_H+y_m*gamma_H2)
    l2%C_01=(l2%g_1/l2%g_0)*l2%C_10*dexp(-DT_10/T_K)

    xm_C=12.d0*xm_p
    v_th=dsqrt(2.d0*xk_B*T_K/xm_C)
    l2%tau0=(l2%A_10/8.d0/pi)*(3.d10/xnu_10)**3*xNc_CII/v_th

    call CIIpop(f_0,f_1,esc_10, l2)

    f_0 = max(f_0, 1.d-30)
    f_1 = max(f_1, 1.d-30)

    S_10=1.d0/(l2%g_1*f_0/(l2%g_0*f_1)-1.d0)
    xLd_CII=DE_10*l2%A_10*f_1*esc_10*(1.d0-l2%Q_10/S_10) !/xnH

    return
  end subroutine CIIcool

  subroutine CIIpop(f_0,f_1,esc_10, l2)
    implicit none
    type(line2) :: l2
    real(kind=DBL_KIND) :: f_0, f_1, esc_10
    real(kind=DBL_KIND), parameter :: err_eps=1.d-4
    real(kind=DBL_KIND) :: esc(1),error(1),desc(1),esc_f(1),error_f(1),A(1,1)
    real(kind=DBL_KIND) :: err_max, fact
    integer :: itr, i, j, delta_esc,jj
    integer :: ipiv(1:1), info
    integer, parameter :: inc=1, n=1
#if USE_LAPACK_CHEM == NO
    integer, dimension(1) :: indx
    real(kind=DBL_KIND) :: d
#endif


    esc(1)=esc_10

    do itr=1,1000
       call twolevel(f_0,f_1,esc,error,l2)

!   evaluate error
       err_max=0.d0
       do i=1,1
          err_max=max(dabs(error(i)),err_max)
       enddo
       if(err_max.lt.err_eps) exit

       do j=1,1
          if(esc(j).eq.0.d0) then
             delta_esc=1.d-10
          else
             delta_esc=1.d-5*esc(j)
          endif

          do jj=1,1
             if(jj.eq.j) then
                esc_f(jj)=esc(jj)+delta_esc
             else
                esc_f(jj)=esc(jj)
             endif
          enddo

          call twolevel(f_0,f_1,esc_f,error_f,l2)
! ---- set the matrix A
          do i=1,1              
             A(i,j)=(error_f(i)-error(i))/(esc_f(j)-esc(j))
          enddo
       enddo

! ----- set the vector desc
       do i=1,1
          desc(i)=-error(i)
       enddo

! ----- solve linear equations
       !call gaussj(A,1,1,desc,1,1,0)
#if USE_LAPACK_CHEM == YES
       call dgesv(n,inc, A, n, ipiv, desc,n,info)
#else
          call ludcmp(A, indx, d, 1)
          call lubksb(A, indx, desc, 1)
#endif
        

       fact=1.d0
       if(itr.gt.20) then
          if(err_max.gt.1.d0) fact=1.d-2
       endif
       do i=1,1
          if(esc(i)*desc(i).ne.0.d0) then
             fact=min(fact,4.d-1*dabs(esc(i)/desc(i)))
          endif 
       enddo

!     improve esc
       do i=1,1
          esc(i)=esc(i)+fact*desc(i)
       enddo 
      enddo
      
      esc_10=esc(1)
  end subroutine CIIpop

  subroutine twolevel(f_0,f_1,esc,error,l2)
!   solves level population (f_0, f_1) of two-level system
!   for given escape probability (esc)
!   returns the difference of old and new esc as (error). 
    implicit none
    type(line2) :: l2
    real(kind=DBL_KIND) :: f_0, f_1, esc_10, R_01, R_10, tau_10
    real(kind=DBL_KIND) :: esc(1),error(1)
    
    esc_10=esc(1)
!   population of level 0 & 1
    R_01=(l2%g_1/l2%g_0)*l2%A_10*esc_10*l2%Q_10+l2%C_01
    R_10=esc_10*l2%A_10*(1.d0+l2%Q_10)+l2%C_10
    f_1=R_01/(R_10+R_01)
    f_0=R_10/(R_10+R_01)
    tau_10=l2%tau0*(f_0*l2%g_1/l2%g_0-f_1)
    esc(1)=beta_esc(tau_10,l2%tau_cnt)
    error(1)=esc_10-esc(1)
    return
  end subroutine twolevel



  ! -------------------------------------------------------
  !                   OI cooling 
  ! -------------------------------------------------------
  subroutine OIcool(xnH,T_K,Trad,xNc_OI,y_m,y_a,y_e,esc,tau_c,xLd_OI)
    implicit none
    real(kind=DBL_KIND) :: xnH,T_K,Trad,xNc_OI,y_m,y_a,y_e,esc(1:3),tau_c,xLd_OI
    real(kind=DBL_KIND), parameter :: xk_B = 1.38d-16
    real(kind=DBL_KIND), parameter :: h_P  = 6.63d-27
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24
    type(line3) :: l3
    real(kind=DBL_KIND) :: DT_10, DT_20, DT_21, DE_10, DE_20, DE_21 &
      , gamma10_e, gamma20_e, gamma21_e, gamma10_H, gamma20_H, gamma21_H &
      , gamma10_H2, gamma20_H2, gamma21_H2, xm_O, f_0, f_1, f_2, S_10, S_20 &
      , S_21, x_10, x_20, x_21


    l3%tau_cnt=tau_c
    if(T_K < 10.d0) then
       xLd_OI=0.d0
       return
    endif
!   set coefficients
    l3%g_0=5.d0
    l3%g_1=3.d0
    l3%g_2=1.d0
!
    DT_10=2.3d2
    DT_20=3.28d2
    DT_21=9.8d1
!
    DE_10=DT_10*xk_B
    DE_20=DT_20*xk_B
    DE_21=DT_21*xk_B
!   
    l3%xnu_10=DE_10/h_P
    l3%xnu_20=DE_20/h_P
    l3%xnu_21=DE_21/h_P
!
    l3%A_10=9.0d-5
    l3%A_20=1.0d-10
    l3%A_21=1.7d-5
!
    l3%Q_10=Q_bg(DT_10,Trad)
    l3%Q_20=Q_bg(DT_20,Trad)
    l3%Q_21=Q_bg(DT_21,Trad)
!
    gamma10_e=1.4d-8
    gamma20_e=1.4d-8
    gamma21_e=5.0d-9
!
    gamma10_H=9.2d-11*(T_K*1.d-2)**0.67d0
    gamma20_H=4.3d-11*(T_K*1.d-2)**0.80d0
    gamma21_H=1.1d-10*(T_K*1.d-2)**0.44d0
!
    gamma10_H2=5.d-1*gamma10_H
    gamma20_H2=5.d-1*gamma20_H
    gamma21_H2=5.d-2*gamma21_H
!
    l3%C_10=xnH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2)
    l3%C_20=xnH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2)
    l3%C_21=xnH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2)
    l3%C_01=(l3%g_1/l3%g_0)*l3%C_10*dexp(-DT_10/T_K)
    l3%C_02=(l3%g_2/l3%g_0)*l3%C_20*dexp(-DT_20/T_K)
    l3%C_12=(l3%g_2/l3%g_1)*l3%C_21*dexp(-DT_21/T_K)
!

    l3%xNclmn=xNc_OI
    xm_O=16.d0*xm_p
    l3%v_th=dsqrt(2.d0*xk_B*T_K/xm_O)           
    call pop3lev(f_0,f_1,f_2,esc,l3)

    f_0 = max(f_0, 1.d-30)
    f_1 = max(f_1, 1.d-30)
    f_2 = max(f_2, 1.d-30)
!
    S_10=1.d0/(l3%g_1*f_0/(l3%g_0*f_1)-1.d0)
    S_20=1.d0/(l3%g_2*f_0/(l3%g_0*f_2)-1.d0)
    S_21=1.d0/(l3%g_2*f_1/(l3%g_1*f_2)-1.d0)
!
    x_10=esc(1)*(1.d0-l3%Q_10/S_10)
    x_20=esc(2)*(1.d0-l3%Q_20/S_20)
    x_21=esc(3)*(1.d0-l3%Q_21/S_21)
!
    xLd_OI=(DE_10*l3%A_10*f_1*x_10+DE_20*l3%A_20*f_2*x_20+DE_21*l3%A_21*f_2*x_21) ! erg s^-1!/xnH

    !print *, "SS", xLd_OI, S_10, S_20, S_21, f_0, f_1

    end subroutine OIcool






  ! --------------------------------------------------------
  !                   H2O cooling 
  ! --------------------------------------------------------

  subroutine H2Ocool(xnH,T_K,y_H2,xNc_H2O,tau_cnt,xLd_H2O)
    
    !H2O cooling function by Neufeld et al.
    implicit none
    real(kind=DBL_KIND), parameter :: xk_B = 1.38d-16
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24
    real(kind=DBL_KIND) :: xnH,T_K,y_H2,xNc_H2O,tau_cnt,xLd_H2O
    real(kind=DBL_KIND) :: xn_c, xlT, cs, xlN, aL0, aLLTE, xlnh, alpha &
        , xL0inv, xLLTEinv, xn_h, xLinv, xL, xlNo, xLo, xlNp, xLp, aL0p
    ! 100, 200, 400, 1000, 2000, 4000K
    real(kind=DBL_KIND) :: xlTa(1:6)=(/2.000d0, 2.301d0, 2.602d0, 3.000d0, 3.301d0, 3.602d0/)
    real(kind=DBL_KIND) :: xlNa(1:10)=(/10.0d0, 11.0d0, 12.0d0, 13.0d0, 14.0d0, &
          15.0d0, 16.0d0, 17.0d0, 18.0d0, 19.0d0/)  
    real(kind=DBL_KIND) :: aL0a(1:6)=(/24.35d0, 23.87d0, 23.42d0, 22.88d0, 22.50d0, 22.14d0/)
    real(kind=DBL_KIND) :: aLLTEa(1:6,1:10)=reshape((/14.59d0 &
       , 13.85d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
         14.59d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
         14.60d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
         14.68d0, 13.88d0, 13.17d0, 12.32d0, 11.86d0, 11.64d0,&
         14.98d0, 14.05d0, 13.25d0, 12.34d0, 11.87d0, 11.65d0,&
         15.53d0, 14.46d0, 13.53d0, 12.49d0, 11.97d0, 11.72d0,&
         16.22d0, 15.05d0, 14.02d0, 12.87d0, 12.35d0, 12.06d0,&
         17.00d0, 15.74d0, 14.63d0, 13.46d0, 12.97d0, 12.66d0,&
         17.83d0, 16.50d0, 15.32d0, 14.16d0, 13.69d0, 13.36d0,&
         18.70d0, 17.31d0, 16.07d0, 14.94d0, 14.46d0, 14.13d0/),(/6,10/))
     real(kind=DBL_KIND) :: xlnha(1:6,1:10) = reshape((/9.00d0&
        , 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0, &
          8.99d0, 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0, &
          8.96d0, 9.03d0, 9.19d0, 9.50d0, 9.66d0, 9.59d0, &
          8.74d0, 8.89d0, 9.11d0, 9.47d0, 9.65d0, 9.59d0, &
          8.11d0, 8.37d0, 8.73d0, 9.31d0, 9.56d0, 9.53d0, &
          7.20d0, 7.51d0, 7.95d0, 8.74d0, 9.15d0, 9.20d0, &
          6.22d0, 6.53d0, 6.99d0, 7.87d0, 8.38d0, 8.50d0, &
          5.22d0, 5.57d0, 6.03d0, 6.94d0, 7.48d0, 7.64d0, &
          4.24d0, 4.59d0, 5.09d0, 6.02d0, 6.59d0, 6.78d0, &
          3.21d0, 3.58d0, 4.10d0, 5.08d0, 5.69d0, 5.89d0/),(/6,10/))
     real(kind=DBL_KIND) :: alphaa(1:6,1:10) = reshape((/0.43d0 &
        , 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0, &  
          0.43d0, 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0, &  
          0.42d0, 0.41d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0, &  
          0.41d0, 0.39d0, 0.37d0, 0.35d0, 0.33d0, 0.33d0, &  
          0.42d0, 0.38d0, 0.34d0, 0.33d0, 0.32d0, 0.32d0, &  
          0.45d0, 0.38d0, 0.34d0, 0.32d0, 0.30d0, 0.30d0, &  
          0.47d0, 0.40d0, 0.35d0, 0.32d0, 0.29d0, 0.30d0, &  
          0.50d0, 0.42d0, 0.36d0, 0.32d0, 0.28d0, 0.29d0, &  
          0.52d0, 0.44d0, 0.37d0, 0.31d0, 0.27d0, 0.28d0, &  
          0.53d0, 0.45d0, 0.39d0, 0.31d0, 0.27d0, 0.27d0/),(/6,10/))
!     
!     10, 20, 30, 50, 80, 100K
      real(kind=DBL_KIND) :: xlTb(1:6)=(/1.000d0, 1.301d0, 1.477d0, &
          1.699d0, 1.903d0, 2.000d0/)
      real(kind=DBL_KIND) :: aL0bo(1:6)=(/26.81d0, 25.88d0, 25.43d0,&
          24.96d0, 24.58d0, 24.41d0/)
      real(kind=DBL_KIND) :: aLLTEbo(1:6,1:10)=reshape((/17.94d0 &
        , 16.71d0, 16.08d0, 15.41d0, 14.85d0, 14.60d0, & 
          17.96d0, 16.72d0, 16.09d0, 15.42d0, 14.86d0, 14.60d0, &
          18.14d0, 16.86d0, 16.19d0, 15.47d0, 14.88d0, 14.62d0, &
          18.77d0, 17.36d0, 16.58d0, 15.72d0, 15.02d0, 14.73d0, &
          19.70d0, 18.11d0, 17.25d0, 16.27d0, 15.47d0, 15.11d0, &
          20.67d0, 18.93d0, 18.05d0, 17.01d0, 16.12d0, 15.72d0, &
          21.58d0, 19.81d0, 18.91d0, 17.82d0, 16.88d0, 16.45d0, &
          22.53d0, 20.71d0, 19.80d0, 18.69d0, 17.70d0, 17.25d0, &
          23.50d0, 21.64d0, 20.72d0, 19.58d0, 18.56d0, 18.10d0, &
          24.41d0, 22.58d0, 21.65d0, 20.49d0, 19.46d0, 18.99d0/),(/6,10/))
      real(kind=DBL_KIND) :: xlnhbo(1:6,1:10) = reshape((/8.81d0 &
          , 8.90d0, 9.03d0, 9.14d0, 9.19d0, 9.20d0, & 
            8.79d0, 8.88d0, 9.01d0, 9.13d0, 9.20d0, 9.20d0, & 
            8.60d0, 8.73d0, 8.88d0, 9.03d0, 9.13d0, 9.15d0, & 
            7.96d0, 8.14d0, 8.31d0, 8.55d0, 8.78d0, 8.85d0, & 
            7.01d0, 7.21d0, 7.40d0, 7.69d0, 8.00d0, 8.11d0, & 
            6.02d0, 6.20d0, 6.41d0, 6.71d0, 7.04d0, 7.17d0, & 
            5.03d0, 5.21d0, 5.42d0, 5.70d0, 6.06d0, 6.19d0, & 
            4.02d0, 4.21d0, 4.41d0, 4.71d0, 5.05d0, 5.19d0, & 
            3.02d0, 3.22d0, 3.41d0, 3.71d0, 4.06d0, 4.20d0, & 
            2.03d0, 2.21d0, 2.42d0, 2.70d0, 3.06d0, 3.20d0/),(/6,10/))
      real(kind=DBL_KIND) :: alphabo(1:6,1:10) =reshape((/0.71d0 &
         , 0.49d0, 0.48d0, 0.46d0, 0.45d0, 0.45d0, & 
           0.64d0, 0.49d0, 0.48d0, 0.45d0, 0.45d0, 0.45d0, & 
           0.57d0, 0.50d0, 0.47d0, 0.45d0, 0.44d0, 0.44d0, & 
           0.56d0, 0.53d0, 0.48d0, 0.44d0, 0.42d0, 0.41d0, & 
           0.59d0, 0.60d0, 0.53d0, 0.47d0, 0.45d0, 0.43d0, & 
           0.72d0, 0.64d0, 0.58d0, 0.52d0, 0.49d0, 0.47d0, & 
           0.85d0, 0.68d0, 0.61d0, 0.56d0, 0.52d0, 0.50d0, & 
           0.86d0, 0.70d0, 0.63d0, 0.58d0, 0.55d0, 0.53d0, & 
           0.93d0, 0.72d0, 0.66d0, 0.61d0, 0.57d0, 0.55d0, & 
           0.87d0, 0.73d0, 0.67d0, 0.62d0, 0.59d0, 0.56d0/),(/6,10/))
!
      real(kind=DBL_KIND) :: aL0bp(1:6)=(/27.01d0, 25.73d0, 25.24d0, &
           24.75d0, 24.38d0, 24.22d0/)
      real(kind=DBL_KIND) :: aLLTEbp(1:6,1:10)=reshape ((/17.72d0&
         , 16.60d0, 16.12d0, 15.43d0, 14.86d0, 14.60d0, &
           17.76d0, 16.63d0, 16.13d0, 15.43d0, 14.86d0, 14.60d0, & 
           18.07d0, 16.85d0, 16.24d0, 15.48d0, 14.88d0, 14.62d0, & 
           18.83d0, 17.41d0, 16.61d0, 15.72d0, 15.03d0, 14.73d0, &
           19.68d0, 18.11d0, 17.26d0, 16.28d0, 15.47d0, 15.11d0, & 
           20.50d0, 18.94d0, 18.06d0, 17.01d0, 16.12d0, 15.72d0, &
           21.37d0, 19.83d0, 18.93d0, 17.82d0, 16.87d0, 16.45d0, & 
           22.28d0, 20.75d0, 19.81d0, 18.69d0, 17.70d0, 17.25d0, & 
           23.22d0, 21.69d0, 20.73d0, 19.58d0, 18.56d0, 18.10d0, & 
           24.19d0, 22.60d0, 21.65d0, 20.49d0, 19.45d0, 18.99d0/),(/6,10/))
      real(kind=DBL_KIND) :: xlnhbp(1:6,1:10) =reshape((/9.30d0 &
         , 9.11d0, 8.94d0, 8.70d0, 8.55d0, 8.50d0, & 
           9.25d0, 9.06d0, 8.91d0, 8.69d0, 8.54d0, 8.49d0, &
           8.94d0, 8.82d0, 8.71d0, 8.56d0, 8.46d0, 8.43d0, &
           8.17d0, 8.16d0, 8.16d0, 8.12d0, 8.10d0, 8.11d0, &
           7.21d0, 7.24d0, 7.28d0, 7.30d0, 7.32d0, 7.35d0, & 
           6.20d0, 6.24d0, 6.31d0, 6.33d0, 6.35d0, 6.38d0, &
           5.21d0, 5.25d0, 5.30d0, 5.32d0, 5.35d0, 5.39d0, &
           4.22d0, 4.26d0, 4.31d0, 4.33d0, 4.36d0, 4.40d0, &
           3.23d0, 3.24d0, 3.31d0, 3.33d0, 3.37d0, 3.41d0, &
           2.21d0, 2.25d0, 2.30d0, 2.32d0, 2.37d0, 2.41d0/),(/6,10/))
      real(kind=DBL_KIND) :: alphabp(1:6,1:10) =reshape ((/0.49d0 &
         , 0.72d0, 0.69d0, 0.53d0, 0.46d0, 0.44d0, & 
           0.52d0, 0.65d0, 0.66d0, 0.53d0, 0.46d0, 0.44d0, &
           0.49d0, 0.65d0, 0.63d0, 0.51d0, 0.45d0, 0.43d0, &
           0.57d0, 0.73d0, 0.65d0, 0.49d0, 0.43d0, 0.41d0, &
           0.75d0, 0.68d0, 0.64d0, 0.52d0, 0.45d0, 0.42d0, &
           0.76d0, 0.70d0, 0.68d0, 0.55d0, 0.47d0, 0.43d0, &
           0.79d0, 0.73d0, 0.69d0, 0.58d0, 0.48d0, 0.44d0, &
           0.80d0, 0.75d0, 0.71d0, 0.58d0, 0.50d0, 0.46d0, &
           0.82d0, 0.77d0, 0.73d0, 0.60d0, 0.52d0, 0.48d0, &
           0.91d0, 0.80d0, 0.74d0, 0.62d0, 0.53d0, 0.50d0/),(/6,10/))
!
      xn_c=(1.d0-y_H2)*xnH
!

      xlT=dlog10(T_K)
      cs=1.d-5*dsqrt(2.d0*xk_B*T_K/(18.d0*xm_p))
      xNc_H2O=xNc_H2O+1.d-10
      if(xlT > 2.d0) then
         xlN=dlog10(xNc_H2O/cs)
         call linear(xlTa,aL0a,6,xlT,aL0)
         call bilinear(xlTa,xlNa,aLLTEa,6,10,xlT,xlN,aLLTE)
         call bilinear(xlTa,xlNa,xlnha,6,10,xlT,xlN,xlnh)
         call bilinear(xlTa,xlNa,alphaa,6,10,xlT,xlN,alpha)
         xL0inv=10.d0**aL0
         xLLTEinv=10.d0**aLLTE
         xn_h=10.d0**xlnh
         xLinv=xL0inv+xn_c*xLLTEinv+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
         xL=1.d0/xLinv
      else
         xlNo=dlog10(0.75d0*xNc_H2O/cs)
         call linear(xlTb,aL0bo,6,xlT,aL0)
         call bilinear(xlTb,xlNa,aLLTEbo,6,10,xlT,xlNo,aLLTE)
         call bilinear(xlTb,xlNa,xlnhbo,6,10,xlT,xlNo,xlnh)
         call bilinear(xlTb,xlNa,alphabo,6,10,xlT,xlNo,alpha)
         xL0inv=10.d0**aL0
         xLLTEinv=10.d0**aLLTE
         xn_h=10.d0**xlnh
         xLinv=xL0inv+xn_c*xLLTEinv+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
         xLo=1.d0/xLinv
!
         xlNp=dlog10(0.25d0*xNc_H2O/cs)
         call linear(xlTb,aL0bp,6,xlT,aL0p)
         call bilinear(xlTb,xlNa,aLLTEbp,6,10,xlT,xlNp,aLLTE)
         call bilinear(xlTb,xlNa,xlnhbp,6,10,xlT,xlNp,xlnh)
         call bilinear(xlTb,xlNa,alphabp,6,10,xlT,xlNp,alpha)
         xL0inv=10.d0**aL0
         xLLTEinv=10.d0**aLLTE
         xn_h=10.d0**xlnh
         xLinv=xL0inv+xn_c*xLLTEinv+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
         xLp=1.d0/xLinv
!
         xL=0.75d0*xLo+0.25d0*xLp
         xL=1.148d0*xL
      endif
!
      xLd_H2O=(1.d0-y_H2)*xL*dexp(-tau_cnt)
      

   end subroutine H2Ocool


  ! ---------------------------------------------------------
  !                     OH cooling
  ! ---------------------------------------------------------
  subroutine OHcool(xnH,T_K,y_H2,xNc_OH,tau_cnt,xLd_OH)

    !OH cooling function by Leiden data 
    implicit none
    real(kind=DBL_KIND), parameter :: xk_B = 1.38d-16
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24
    real(kind=DBL_KIND) :: xnH, T_K, y_H2, xNc_OH, tau_cnt, xLd_OH
    real(kind=DBL_KIND) :: xn_c, xlT, cd, xlN, aL0, aLLTE, xlnh, alpha
    real(kind=DBL_KIND) :: xL0inv, xLLTEinv, xn_h, xLinv, xL
    real(kind=DBL_KIND) :: xlTa(1:6)=(/1.477d0, 1.699d0, 1.903d0, 2.000d0, 2.477d0, 2.778d0/)
    real(kind=DBL_KIND) :: xlNa(1:9)=(/10.0d0, 11.d0, 12.0d0, 13.d0, 14.d0, 15.d0 &
      , 16.d0, 17.d0, 18.d0/)
    real(kind=DBL_KIND) :: aL0a(1:6)= &
        (/0.2531d2, 0.2453d2, 0.2402d2, 0.2382d2, 0.2316d2, 0.2290d2/)
    real(kind=DBL_KIND) :: aLLTEa(1:6,1:9) = reshape((/ 0.1620d2 &
        , 0.1544d2, 0.1490d2, 0.1467d2, 0.1381d2, 0.1358d2, & 
          0.1620d2, 0.1545d2, 0.1490d2, 0.1467d2, 0.1381d2, 0.1358d2,   &
          0.1622d2, 0.1546d2, 0.1491d2, 0.1467d2, 0.1381d2, 0.1359d2,   &
          0.1636d2, 0.1555d2, 0.1496d2, 0.1471d2, 0.1383d2, 0.1359d2,   &
          0.1702d2, 0.1600d2, 0.1525d2, 0.1495d2, 0.1394d2, 0.1367d2,   &
          0.1774d2, 0.1660d2, 0.1577d2, 0.1545d2, 0.1449d2, 0.1415d2,   &
          0.1858d2, 0.1738d2, 0.1650d2, 0.1616d2, 0.1510d2, 0.1477d2,   &
          0.1942d2, 0.1828d2, 0.1740d2, 0.1705d2, 0.1584d2, 0.1540d2,   &
          0.2038d2, 0.1921d2, 0.1836d2, 0.1802d2, 0.1683d2, 0.1635d2/),(/6,9/))
    real(kind=DBL_KIND) :: xlnha(1:6,1:9) = reshape((/  0.905d1, &
          0.898d1,  0.890d1,  0.888d1,  0.877d1,  0.872d1, &
          0.905d1,  0.898d1,  0.890d1,  0.887d1,  0.877d1,  0.871d1,    &  
          0.903d1,  0.897d1,  0.889d1,  0.887d1,  0.876d1,  0.871d1,    &
          0.888d1,  0.885d1,  0.881d1,  0.880d1,  0.874d1,  0.869d1,    &
          0.820d1,  0.827d1,  0.836d1,  0.840d1,  0.850d1,  0.850d1,    &
          0.725d1,  0.745d1,  0.771d1,  0.779d1,  0.799d1,  0.801d1,    &
          0.626d1,  0.646d1,  0.675d1,  0.685d1,  0.722d1,  0.729d1,    &
          0.526d1,  0.546d1,  0.575d1,  0.585d1,  0.623d1,  0.633d1,    &
          0.426d1,  0.446d1,  0.475d1,  0.485d1,  0.523d1,  0.533d1/),(/6,9/))
    real(kind=DBL_KIND) :: alphaa(1:6,1:9) = reshape((/ 0.353d0      &
        , 0.543d0,  0.536d0,  0.530d0,  0.466d0,  0.434d0,  &
          0.354d0,  0.543d0,  0.536d0,  0.530d0,  0.466d0,  0.434d0,    & 
          0.364d0,  0.544d0,  0.536d0,  0.529d0,  0.466d0,  0.434d0,    &
          0.399d0,  0.556d0,  0.534d0,  0.525d0,  0.463d0,  0.432d0,    &
          0.655d0,  0.586d0,  0.532d0,  0.521d0,  0.447d0,  0.424d0,    &
          0.641d0,  0.538d0,  0.538d0,  0.538d0,  0.448d0,  0.443d0,    &
          0.649d0,  0.572d0,  0.552d0,  0.516d0,  0.381d0,  0.370d0,    &
          0.566d0,  0.618d0,  0.580d0,  0.538d0,  0.381d0,  0.363d0,    &
          0.701d0,  0.645d0,  0.598d0,  0.550d0,  0.383d0,  0.369d0/),(/6,9/))
!
    xn_c=(1.d0-y_H2)*xnH
!
    xlT=dlog10(T_K)
    cd=1.d-5*dsqrt(2.d0*xk_B*T_K/(17.d0*xm_p))
    xNc_OH=xNc_OH+1.d-10
    xlN=dlog10(xNc_OH/cd)
    call linear(xlTa,aL0a,6,xlT,aL0)
    call bilinear(xlTa,xlNa,aLLTEa,6,9,xlT,xlN,aLLTE)
    call bilinear(xlTa,xlNa,xlnha,6,9,xlT,xlN,xlnh)
    call bilinear(xlTa,xlNa,alphaa,6,9,xlT,xlN,alpha)
!
    xL0inv=10.d0**aL0
    xLLTEinv=10.d0**aLLTE
    xn_h=10.d0**xlnh
    xLinv=xL0inv+xn_c*xLLTEinv+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
    xL=1.d0/xLinv
!
    xLd_OH=(1.d0-y_H2)*xL*dexp(-tau_cnt)
  end subroutine OHcool


  subroutine linear(xa,ya,m,x,y)
    implicit none
    integer :: m, ms,i
    real(kind=DBL_KIND), dimension(m) :: xa, ya
    real(kind=DBL_KIND) :: x, y, y1, y2, t
    
    do 11 i=1,m
        if(x-xa(i).le.0.d0) then
            ms=i
            go to 12
        endif
11  continue
    ms=m
12  continue
    if(ms.eq.1) ms=2
    y1=ya(ms-1)
    y2=ya(ms)
    t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
    y=(1.d0-t)*y1+t*y2
    return
  end subroutine linear

  subroutine bilinear(x1a,x2a,ya,m,n,x1,x2,y)
    implicit none
    integer :: m, n, i, ms, ns
    real(kind=DBL_KIND) :: x1, x2, y
    real(kind=DBL_KIND), dimension(m) :: x1a
    real(kind=DBL_KIND), dimension(n) :: x2a
    real(kind=DBL_KIND), dimension(m,n) :: ya
    real(kind=DBL_KIND) :: y1, y2, y3, y4, t, u

    do 11 i=1,m
        if(x1-x1a(i).le.0.d0) then
            ms=i
            go to 12  
        endif
11  continue 
    ms=m
12  continue
    do 13 i=1,n
        if(x2-x2a(i).le.0.d0) then
            ns=i
            go to 14  
        endif
13  continue
    ns=n
14  continue
    if(ms.eq.1) ms=2
    if(ns.eq.1) ns=2
    y1=ya(ms-1,ns-1)
    y2=ya(ms,ns-1)
    y3=ya(ms,ns)
    y4=ya(ms-1,ns)
    t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
    u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
    y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
    return
  end subroutine bilinear
     

  ! ---------------------------------------------------------
  !                     HD cooling
  ! ---------------------------------------------------------
  subroutine HDcool(xnH,T_K, Trad,y_H,y_H2,xNc_HD,tau_cnt,xLd_HD)
    ! HD cooling from Galli and Palla (1998), A&A, 335, 403
    ! cooling is for T < 3000 K, in erg cm^3 s-1
    implicit none
    real(kind=DBL_KIND), parameter :: xk_B = 1.38066d-16
    real(kind=DBL_KIND), parameter :: h_Pl = 6.62618d-27
    real(kind=DBL_KIND), parameter :: pi   = 3.14159265358979d0
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24
    real(kind=DBL_KIND) :: xnH,T_K, Trad,y_H,y_H2,xNc_HD,tau_cnt
    real(kind=DBL_KIND) :: g0, g1, g2, g3
    real(kind=DBL_KIND) :: a10d, a21d, a32d, a20d, a31d
    real(kind=DBL_KIND) :: xnu10d, xnu21d, xnu32d, xnu20d, xnu31d
    real(kind=DBL_KIND) :: DT10d, DT21d, DT32d, DT20d, DT31d
    real(kind=DBL_KIND) :: c10d, c21d, c32d, c20d, c31d, c01d, c12d, c23d, c02d, c13d
    real(kind=DBL_KIND) :: q10d, q21d, q32d, q20d, q31d
    real(kind=DBL_KIND) :: w10d, w21d, w32d, w20d, w31d, w30d, w01d, w12d, w23d, w02d, w13d, w03d
    real(kind=DBL_KIND) :: w0d, w1d, w2d, an0d, an1d, an2d, an3d, antotd
    real(kind=DBL_KIND) :: f0d, f1d, f2d, f3d, xNc0d, xNc1d, xNc2d, xNc3d, v_D, tau10d, tau21d, tau32d, tau20d, tau31d
    real(kind=DBL_KIND) :: esc10d, esc21d, esc32d, esc20d, esc31d, s10d, x10d, s21d, s20d, x21d &
      , x20d, s32d, s31d, x32d, x31d
    real(kind=DBL_KIND) :: xLd10d, xLd21d, xLd32d, xLd20d, xLd31d, xLd_HD

    !     USES Q_bg, beta_esc
    !------------------------------------------------
    !     Various parameters for HD
    !-----------------------------------------------
    g0=1.d0
    g1=3.d0
    g2=5.d0
    g3=7.d0

!     radiative transitions
    a10d = 5.12d-8
    a21d = 4.86d-7
    a32d = 1.72d-6
    a20d = 7.05d-12
    a31d = 1.15d-10

!     frequencies in Hz
    xnu10d = 2.67d12
    xnu21d = 5.32d12
    xnu32d = 7.93d12
    xnu20d = xnu21d+xnu10d
    xnu31d = xnu32d+xnu21d
!     in K
    DT10d = h_Pl*xnu10d/xk_B
    DT21d = h_Pl*xnu21d/xk_B
    DT32d = h_Pl*xnu32d/xk_B
    DT20d = h_Pl*xnu20d/xk_B
    DT31d = h_Pl*xnu31d/xk_B


    !---------------------------------------------------------------
    !Collisional transitions for HD, from Shaefer (1990) [HD,He]
    !with a corrective factor for [HD,H] from Wright & Morton (1979)
    !---------------------------------------------------------------
    c10d = (4.4d-12+3.6d-13*T_K**0.77d0)*y_H*xnH/1.27d0
    c21d = (4.1d-12+2.0d-13*T_K**0.92d0)*y_H*xnH/1.27d0
    c32d = (2.4d-12+8.7d-14*T_K**1.03d0)*y_H*xnH/1.27d0
    c20d = (3.4d-13+1.1d-14*T_K**1.12d0)*y_H*xnH/1.27d0
    c31d = (3.2d-13+1.3d-15*T_K**1.47d0)*y_H*xnH/1.27d0
!
    c01d = (g1/g0)*c10d*dexp(-h_Pl*xnu10d/(xk_B*T_K))
    c12d = (g2/g1)*c21d*dexp(-h_Pl*xnu21d/(xk_B*T_K))
    c23d = (g3/g2)*c32d*dexp(-h_Pl*xnu32d/(xk_B*T_K))
    c02d = (g2/g0)*c20d*dexp(-h_Pl*xnu20d/(xk_B*T_K))
    c13d = (g3/g1)*c31d*dexp(-h_Pl*xnu31d/(xk_B*T_K))
!
!   
    q10d = Q_bg(DT10d, Trad)
    q21d = Q_bg(DT21d, Trad)
    q32d = Q_bg(DT32d, Trad)
    q20d = Q_bg(DT20d, Trad)
    q31d = Q_bg(DT31d, Trad)

    !---------------------------------------------------------------
    ! Level populations for HD
    !---------------------------------------------------------------
      
    w10d=a10d*(1.d0+q10d)+c10d
    w21d=a21d*(1.d0+q21d)+c21d
    w32d=a32d*(1.d0+q32d)+c32d
    w20d=a20d*(1.d0+q20d)+c20d
    w31d=a31d*(1.d0+q31d)+c31d
    w30d=0.d0
!    
    w01d=(g1/g0)*a10d*q10d+c01d
    w12d=(g2/g1)*a21d*q21d+c12d
    w23d=(g3/g2)*a32d*q32d+c23d
    w02d=(g2/g0)*a20d*q20d+c02d
    w13d=(g3/g1)*a31d*q31d+c13d
    w03d=0.d0
!
    w0d = w01d+w02d+w03d
    w1d = w10d+w12d+w13d
    w2d = w20d+w21d+w23d
!
    an0d = -(w10d*w21d*w32d+w30d*w2d*w1d+w20d*w31d*w12d-  &
              w12d*w21d*w30d+w2d*w31d*w10d+w32d*w20d*w1d)
    an1d = -w0d*w21d*w32d+w20d*w31d*w02d-w30d*w2d*w01d-   & 
              w02d*w21d*w30d-w2d*w31d*w0d-w32d*w20d*w01d
    an2d = -(w0d*w1d*w32d+w10d*w31d*w02d+w30d*w12d*w01d+  &
              w02d*w1d*w30d+w12d*w31d*w0d-w32d*w10d*w01d)
    an3d = -w0d*w1d*w2d+w10d*w21d*w02d+w20d*w12d*w01d+    &
              w02d*w1d*w20d+w12d*w21d*w0d+w2d*w10d*w01d
    antotd = an0d+an1d+an2d+an3d

    f0d = an0d/antotd
    f1d = an1d/antotd
    f2d = an2d/antotd
    f3d = an3d/antotd
!
    xNc0d = f0d*xNc_HD
    xNc1d = f1d*xNc_HD
    xNc2d = f2d*xNc_HD
    xNc3d = f3d*xNc_HD      
!
    v_D=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p)) 
!
    tau10d = (a10d/8.d0/pi)*(3.d10/xnu10d)**3*(xNc0d*g1/g0-xNc1d)/v_D
    tau21d = (a21d/8.d0/pi)*(3.d10/xnu21d)**3*(xNc1d*g2/g1-xNc2d)/v_D
    tau32d = (a32d/8.d0/pi)*(3.d10/xnu32d)**3*(xNc2d*g3/g2-xNc3d)/v_D
    tau20d = (a20d/8.d0/pi)*(3.d10/xnu20d)**3*(xNc0d*g2/g0-xNc2d)/v_D
    tau31d = (a31d/8.d0/pi)*(3.d10/xnu31d)**3*(xNc1d*g3/g1-xNc3d)/v_D
!
    esc10d=beta_esc(tau10d,tau_cnt)
    esc21d=beta_esc(tau21d,tau_cnt)
    esc32d=beta_esc(tau32d,tau_cnt)
    esc20d=beta_esc(tau20d,tau_cnt)
    esc31d=beta_esc(tau31d,tau_cnt)

!    
    if(f1d /= 0.d0) then
       s10d=1.d0/(g1*f0d/(g0*f1d)-1.d0)
       x10d=esc10d*(1.d0-q10d/s10d)
    else
       x10d=0.d0
    endif
    if(f2d /= 0.d0) then
       s21d=1.d0/(g2*f1d/(g1*f2d)-1.d0)
       s20d=1.d0/(g2*f0d/(g0*f2d)-1.d0)
       x21d=esc21d*(1.d0-q21d/s21d)
       x20d=esc20d*(1.d0-q20d/s20d)
    else
       x21d=0.d0
       x20d=0.d0
    endif
    if(f3d /= 0.d0) then
       s32d=1.d0/(g3*f2d/(g2*f3d)-1.d0)
       s31d=1.d0/(g3*f1d/(g1*f3d)-1.d0)
       x32d=esc32d*(1.d0-q32d/s32d)
       x31d=esc31d*(1.d0-q31d/s31d)
    else
       x32d=0.d0
       x31d=0.d0
    endif
!     
    !---------------------------------------------------------------
    ! Computes cooling rate for HD
    !---------------------------------------------------------------
    xLd10d = x10d*f1d*a10d*h_Pl*xnu10d  !/xnH 
    xLd21d = x21d*f2d*a21d*h_Pl*xnu21d  !/xnH 
    xLd32d = x32d*f3d*a32d*h_Pl*xnu32d  !/xnH 
    xLd20d = x20d*f2d*a20d*h_Pl*xnu20d  !/xnH 
    xLd31d = x31d*f3d*a31d*h_Pl*xnu31d  !/xnH      
    xLd_HD = xLd10d+xLd21d+xLd32d+xLd20d+xLd31d ! [erg s^-1]
    return
  end subroutine HDcool


  ! ---------------------------------------------------------
  !                     CI cooling
  ! ---------------------------------------------------------
  subroutine CIcool(xnH,T_K, Trad,xNc_CI,y_m,y_a,y_e,esc,tau_c,xLd_CI)
    implicit none
    !implicit real*8(a-h,o-z)
    real(kind=DBL_KIND) :: xnH,T_K, Trad,xNc_CI,y_m,y_a,y_e,tau_c,xLd_CI 
    real(kind=DBL_KIND) :: esc(1:3)

    real(kind=DBL_KIND), parameter :: xk_B = 1.38d-16
    real(kind=DBL_KIND), parameter :: h_P  = 6.63d-27
    real(kind=DBL_KIND), parameter :: xm_p = 1.67d-24

    real(kind=DBL_KIND) :: DT_10, DT_20, DT_21, DE_10, DE_20, DE_21
    real(kind=DBL_KIND) :: gamma10_e, gamma20_e, gamma21_e, gamma10_H, gamma20_H, gamma21_H
    real(kind=DBL_KIND) :: gamma10_H2, gamma20_H2, gamma21_H2
    real(kind=DBL_KIND) :: xm_C, f_0, f_1, f_2, S_10, S_20, S_21, x_10, x_20, x_21

    type(line3) :: l3

    ! ---------------------
    if(T_K < 1.d0) then
       xLd_CI=0.d0
       return
    endif
    l3%tau_cnt=tau_c
    ! ---------------------
!
!   set coefficients
!
    l3%g_0=1.d0
    l3%g_1=3.d0
    l3%g_2=5.d0
!
    DT_10=2.4d1
    DT_20=6.3d1
    DT_21=3.9d1
!   
    DE_10=DT_10*xk_B
    DE_20=DT_20*xk_B
    DE_21=DT_21*xk_B
!   
    l3%xnu_10=DE_10/h_P
    l3%xnu_20=DE_20/h_P
    l3%xnu_21=DE_21/h_P
!
    l3%A_10=7.9d-8
    l3%A_20=2.0d-14
    l3%A_21=2.7d-7
!
    l3%Q_10=Q_bg(DT_10, Trad)
    l3%Q_20=Q_bg(DT_20, Trad)
    l3%Q_21=Q_bg(DT_21, Trad)
!
    gamma10_e=3.0d-9
    gamma20_e=5.0d-9
    gamma21_e=1.5d-8
!
    gamma10_H=1.6d-10*(T_K*1.d-2)**0.14d0
    gamma20_H=9.2d-11*(T_K*1.d-2)**0.26d0
    gamma21_H=2.9d-10*(T_K*1.d-2)**0.26d0
!
    gamma10_H2=5.d-2*gamma10_H
    gamma20_H2=5.d-1*gamma20_H
    gamma21_H2=5.d-1*gamma21_H
!
    l3%C_10=xnH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2)
    l3%C_20=xnH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2)
    l3%C_21=xnH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2)
    l3%C_01=(l3%g_1/l3%g_0)*l3%C_10*dexp(-DT_10/T_K)
    l3%C_02=(l3%g_2/l3%g_0)*l3%C_20*dexp(-DT_20/T_K)
    l3%C_12=(l3%g_2/l3%g_1)*l3%C_21*dexp(-DT_21/T_K)
!
    l3%xNclmn=xNc_CI
    xm_C=12.d0*xm_p
    l3%v_th=dsqrt(2.d0*xk_B*T_K/xm_C)           
    call pop3lev(f_0,f_1,f_2,esc,l3)

    f_0 = max(f_0, 1.d-30)
    f_1 = max(f_1, 1.d-30)
    f_2 = max(f_2, 1.d-30)
!    
    S_10=1.d0/(l3%g_1*f_0/(l3%g_0*f_1)-1.d0)
    S_20=1.d0/(l3%g_2*f_0/(l3%g_0*f_2)-1.d0)
    S_21=1.d0/(l3%g_2*f_1/(l3%g_1*f_2)-1.d0)
!
    x_10=esc(1)*(1.d0-l3%Q_10/S_10)
    x_20=esc(2)*(1.d0-l3%Q_20/S_20)
    x_21=esc(3)*(1.d0-l3%Q_21/S_21)
!
    xLd_CI=(DE_10*l3%A_10*f_1*x_10+DE_20*l3%A_20*f_2*x_20+DE_21*l3%A_21*f_2*x_21) ! erg s^-1    !/xnH
    return
  end subroutine CIcool


  subroutine pop3lev(f_0,f_1,f_2,esc,l3)
    implicit none
    real(kind=DBL_KIND), parameter :: err_eps = 1.d-5
    real(kind=DBL_KIND), dimension(3) :: esc,error,desc,esc_f,error_f
    real(kind=DBL_KIND), dimension(3,3) :: A
    integer :: itr, i, j, jj
    real(kind=DBL_KIND) :: f_0,f_1,f_2
    real(kind=DBL_KIND) :: err_max, delta_esc,fact
    integer, parameter :: inc=1, n=3
    integer :: ipiv(1:3), info
    type(line3) :: l3
#if USE_LAPACK_CHEM == NO
    integer, dimension(3) :: indx
    real(kind=DBL_KIND) :: d
#endif

    do itr=1,100
      
       !determin population under given esc
       call thrlev(f_0,f_1,f_2,esc,error,l3)
       !evaluate error
       err_max=0.d0
       do i=1,3
          err_max=max(dabs(error(i)),err_max)
       enddo
       
       !error small enough ?
       if(err_max.lt.err_eps) exit
       
       !if not, improve guess for esc by Newton-Raphson scheme
       do j=1,3
          if(esc(j).eq.0.d0) then
             delta_esc=1.d-10
          else
             delta_esc=1.d-5*esc(j)
          endif
 
          !do jj=1,3
             !if(jj.eq.j) then
             !   esc_f(jj)=esc(jj)+delta_esc
             !else
             !   esc_f(jj)=esc(jj)
             !endif
          !enddo            
          esc_f(:) = esc(:)
          esc_f(j) = esc_f(j) + delta_esc
          
          call thrlev(f_0,f_1,f_2,esc_f,error_f,l3)
          
          !set the matrix A
          do i=1,3               
             A(i,j)=(error_f(i)-error(i))/(esc_f(j)-esc(j))
          enddo
       enddo
       
       !set the vector desc
       do i=1,3
          desc(i)=-error(i)
       enddo
       
       !solve linear equations
#if USE_LAPACK_CHEM == YES
       call dgesv(n,inc, A, n, ipiv, desc,n,info)
#else
       call ludcmp(A, indx, d, 3)
       call lubksb(A, indx, desc,3)
#endif


       !call gaussj(A,3,3,desc,1,1,0)
!
       fact=1.d0
       if(itr.gt.20) then
          if(err_max.gt.1.d0) fact=1.d-2
       endif
       do i=1,3
          if(esc(i)*desc(i).ne.0.d0) then
             fact=min(fact,4.d-1*dabs(esc(i)/desc(i)))
          endif 
       enddo
       
       !next guess for esc
       do i=1,3
          esc(i)=esc(i)+fact*desc(i)
       enddo 

    enddo

    return
  end subroutine pop3lev


  subroutine thrlev(f_0,f_1,f_2,esc,error,l3)
      
    implicit none
    type(line3) :: l3
    real(kind=DBL_KIND) :: f_0,f_1,f_2
    real(kind=DBL_KIND),dimension(3) :: esc, error
    real(kind=DBL_KIND), parameter :: pi=3.14159265358979d0
    real(kind=DBL_KIND) :: esc_10, esc_20, esc_21
    real(kind=DBL_KIND) :: R_10, R_20, R_21, R_01, R_02, R_12
    real(kind=DBL_KIND) :: xNc_0, xNc_1, xNc_2, tau_10, tau_20, tau_21
    real(kind=DBL_KIND) :: err_10, err_20, err_21

    esc_10=esc(1)
    esc_20=esc(2)
    esc_21=esc(3)

    R_10=esc_10*l3%A_10*(1.d0+l3%Q_10)+l3%C_10
    R_20=esc_20*l3%A_20*(1.d0+l3%Q_20)+l3%C_20
    R_21=esc_21*l3%A_21*(1.d0+l3%Q_21)+l3%C_21
    R_01=(l3%g_1/l3%g_0)*esc_10*l3%A_10*l3%Q_10+l3%C_01
    R_02=(l3%g_2/l3%g_0)*esc_20*l3%A_20*l3%Q_20+l3%C_02
    R_12=(l3%g_2/l3%g_1)*esc_21*l3%A_21*l3%Q_21+l3%C_12
    
    f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) ) &
        /( (R_01+R_02+R_20)*(R_10+R_12+R_21)      &
         -(R_01-R_21)*(R_10-R_20) )
    f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
    f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

    xNc_0=l3%xNclmn*f_0
    xNc_1=l3%xNclmn*f_1
    xNc_2=l3%xNclmn*f_2

    tau_10=(l3%A_10/8.d0/pi)*(3.d10/l3%xnu_10)**3*(xNc_0*l3%g_1/l3%g_0-xNc_1)/l3%v_th 
    tau_20=(l3%A_20/8.d0/pi)*(3.d10/l3%xnu_20)**3*(xNc_0*l3%g_2/l3%g_0-xNc_2)/l3%v_th 
    tau_21=(l3%A_21/8.d0/pi)*(3.d10/l3%xnu_21)**3*(xNc_1*l3%g_2/l3%g_1-xNc_2)/l3%v_th 

    err_10=esc_10-beta_esc(tau_10,l3%tau_cnt)
    err_20=esc_20-beta_esc(tau_20,l3%tau_cnt)
    err_21=esc_21-beta_esc(tau_21,l3%tau_cnt)
    
    error(1)=err_10
    error(2)=err_20
    error(3)=err_21
    
    return
  end subroutine thrlev


  function  Q_bg(T_nu, Trad)
    real(kind=DBL_KIND) :: Q_bg
    real(kind=DBL_KIND),intent(IN) :: T_nu, Trad

    real(kind=DBL_KIND) :: x, Q_bg_CMB 

    x     = T_nu/Trad
    if(x.gt.1.d2) then
       Q_bg_CMB = 0.d0
    else
       Q_bg_CMB = 1.d0/(dexp(x)-1.d0)
    endif
    Q_bg = Q_bg_CMB
    return
  end function Q_bg


  function beta_esc(tau_L,tau_C)
    implicit none
    real(kind=DBL_KIND) :: beta_esc
    real(kind=DBL_KIND) :: tau_L, tau_C

    if(tau_L.lt.0.d0) then
       beta_esc=1.d0
    elseif(tau_L.lt.1.d-5) then 
       beta_esc=dexp(-tau_C)
    else
       beta_esc=dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
    endif
    return
  end function beta_esc

#if USE_LAPACK_CHEM == NO
  ! ----------------------------------------- 
  !             inverse matrix 
  ! -----------------------------------------

  subroutine ludcmp(A, indx, d, ndim)
    
    real(kind=DBL_KIND), dimension(:,:) :: A
    integer, dimension(:) :: indx
    integer, intent(IN) :: ndim
    real(kind=DBL_KIND) :: d
    integer :: i, j, kk, imax_1
    real(kind=DBL_KIND) :: big,dum,summ, tempora
    real(kind=DBL_KIND), dimension(:), allocatable :: vv

    allocate(vv(ndim))

    ! initial -----------
    imax_1 = 1
    tempora= 0.d0
    vv = 0.d0
    ! -------------------

    d = 1.d0

    do i=1, ndim
      
      big = 0.d0

      do j=1, ndim
        tempora = dabs(A(i,j))
        if (tempora > big) big = tempora
      enddo

      if(big == 0.d0) then 
        print *, "Singular matrix in rourine ludcmp"
        print *, "stopping ...."
        stop
      endif

      vv(i) = 1.d0/big
    enddo

    
    do j = 1, ndim
      
      do i=1, j-1
        summ = A(i,j)

        do kk=1, i-1
          summ = summ - A(i,kk)*A(kk,j)
        enddo
        A(i,j) = summ
      enddo

      big = 0.d0

      do i=j, ndim
        
        summ = A(i,j)

        do kk=1, j-1
          summ = summ - A(i, kk)*A(kk,j)
        enddo
        A(i,j) = summ

        dum = vv(i)*dabs(summ)

        if(dum >= big) then
          big    = dum 
          imax_1 = i
        endif
      enddo

      if (j .ne. imax_1) then
        
        do kk = 1, ndim
          
          dum           = A(imax_1, kk)
          A(imax_1, kk) = A(j,kk)
          A(j,kk)       = dum

        enddo
        
        d = -d
        vv(imax_1) = vv(j)

      endif

      indx(j) = imax_1

      if( j .ne. ndim) then
        dum = 1.d0 / (A(j,j))
  
        do i = j+1, ndim
          A(i, j) = A(i, j)*dum
        enddo

      endif

    enddo

  end subroutine ludcmp


  subroutine lubksb(A, indx, b, ndim)

    real(kind=DBL_KIND), dimension(:, :),intent(IN):: A
    integer, dimension(:),intent(IN) :: indx
    real(kind=DBL_KIND), dimension(:), intent(INOUT) :: b
    integer :: i, ii, ip, j, ndim
    real(kind=DBL_KIND) :: summ

    ii = 1

    do i=1, ndim
      ip = indx(i)
      summ = b(ip)
      b(ip)= b(i)

      if (ii .ne. 1) then
        do j=ii-1, i-1
          summ = summ - a(i,j)*b(j)
        enddo
      else if (summ .ne. 0.d0) then
        ii = i+1
      endif

      b(i) = summ
    enddo
    
    do i = ndim, 1, -1
      summ = b(i)

      do j=i+1, ndim
        summ = summ - a(i,j)*b(j)
      enddo
      b(i) = summ/a(i,i)

    enddo

  end subroutine lubksb
#endif

end module kinzoku2
