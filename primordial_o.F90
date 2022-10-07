#include "config.h"

!-- WARNING -- WARNING -- WARNING -- eARNING !
!
! do not use config.h but directly define here
! for debuging purpose
!
!#define NCHEM 6
!#define NREACT 21
!#define DBL_KIND 8
!#define ARRAYSIZE_IJK Imin:Imax,Jmin:Jmax,Kmin:Kmax
!-- WARNING -- WARNING -- WARNING -- WARNING !

!control debug messages
!#define PRIM_DEBUG

!-----------------------------------------------------------------------
! subroutine for primordial-gas proccess (cf. Hosokawa+ 2017)
!-----------------------------------------------------------------------
module primordial
  use kinzoku 
#ifdef DM_NFW_PROFILE
  use modelParameter, only: MP_Tcmb
#endif
  implicit none
  private

  !physical constants
  real(kind=DBL_KIND), parameter :: boltz = 1.380662d-16  ! Boltzman const.  
  real(kind=DBL_KIND), parameter :: evolt = 1.6022E-12      ! eV in erg

  !physical parameters
  real(kind=DBL_KIND), parameter :: yHe      = 9.7222222d-2  ! Helium abundance  
  real(kind=DBL_KIND), parameter :: zred     = 0.d0   ! redshift
#ifndef DM_NFW_PROFILE
  real(kind=DBL_KIND), parameter :: MP_Tcmb     = DEF_TCMB  !2.73d0*(1.d0 + zred) ! CMB temp.
#endif
  real(kind=DBL_KIND), parameter :: CONST_mH = 1.6733e-24 ! hydrogen atom
  real(kind=DBL_KIND), parameter :: CONST_G  = 6.6726e-8
  real(kind=DBL_KIND), parameter :: CONST_amu= 1.66053886e-24
  real(kind=DBL_KIND), parameter :: CONST_kB = 1.3806505e-16
  real(kind=DBL_KIND), parameter :: cgs_sigma = 5.67051D-5

  ! types for chemistry solver
  type chem_vals
    real(kind=DBL_KIND) :: nH                        ! number density
    real(kind=DBL_KIND) :: Tg                         ! temperature
    real(kind=DBL_KIND) :: Td                         ! dust temp
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: ychem ! chemical abundance
    real(kind=DBL_KIND) :: yco                        ! co abundance
    real(kind=DBL_KIND) :: rHpi                       ! H  photo-ionization rate
    real(kind=DBL_KIND) :: rH2pd                      ! H2 photo-dissociation rate
    real(kind=DBL_KIND) :: rHmpd                      ! H- photo-dissociation rate
    real(kind=DBL_KIND) :: heat                       ! heating rate of H ionization
    real(kind=DBL_KIND) :: rgfuv                      ! FUV intensity for photoelectric heating
    real(kind=DBL_KIND) :: rdph                       ! dust heating rate
    real(kind=DBL_KIND) :: rcopd                      ! CO phodissociation rate
    real(kind=DBL_KIND) :: rOII                       ! OII photoionization rate
    real(kind=DBL_KIND) :: xlmbdj                     ! jeans length
    real(kind=DBL_KIND) :: dt                         ! timestep of chemical solver
    real(kind=DBL_KIND) :: fd                         ! dust abundance ratio to solar metallicity 
    real(kind=DBL_KIND) :: EradIR                     ! energy density of IR photons
    real(kind=DBL_KIND) :: chi_d                      ! chi_d for radtr
    real(kind=DBL_KIND) :: metal                      ! metallicity of gas
  end type chem_vals

  type rad_others
    real(kind=DBL_KIND) :: xeuv       ! the emissivity of stellar EUV photons xeuv * lum 
    real(kind=DBL_KIND) :: xfuv       ! the emissivity of stellar FUV photons
    real(kind=DBL_KIND) :: alpha_euv  ! mean cross section for EUV [cm^-2]
    real(kind=DBL_KIND) :: heat_euv   ! mean heating rate per one ionization [erg per one ionization]
    real(kind=DBL_KIND) :: hhm        ! nu averaged H^- cross section at stellar surface [s^-1] 
    real(kind=DBL_KIND) :: lumeuv     ! the raio of EUV luminosity
    real(kind=DBL_KIND) :: lumfuv     ! the raio of EUV luminosity
    real(kind=DBL_KIND) :: sig_euv    ! dust cross section for EUV at Z=Zsun [cm^2]
    real(kind=DBL_KIND) :: sig_fuv    ! dust cross section for FUV at Z=Zsun [cm^2]
    real(kind=DBL_KIND) :: rOII       ! ratio of photodisociation rate of OII and HI
  end type rad_others
  

  integer,save :: dbg_flg_prm = 0 !for debugging

  public :: chemreact, HUVAHHH, yHe, zred, c_H2, c_H2_2, CoolSolverExplicit, CoolSolverImplicit, adjust_abundance, chemreact_adptv
  !public :: prim_GridChemCoolExplicit
  public :: dbg_flg_prm, H2cool, H2cool_HM, H2cool_LTEfit, H2cool_LowD_G15, H2cool_Omukai
  public :: chem_vals
  public :: ProstFit2, rad_others
  external DGESV  
contains

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$             CHEMISTRY          $$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !-----------------------------------------------------------------------
  ! solve the chemical reaction eqs.
  !
  ! tchem is negative if not converging
  !
  ! SPECIES                                                         
  !    0 : H      1 : H2     2 : e     3 : H+     4 : H-    5 : H2+ 
  !-----------------------------------------------------------------------  
  subroutine chemreact(xnH,T_K,Tdust,y,yco,rHpi,rH2pd,rHmpd,fd,metal,dt,tchem,xk_in)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,rHpi,rH2pd,rHmpd,dt,Tdust,fd,metal
    real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1),yco
    real(kind=DBL_KIND),intent(OUT) :: tchem

    real(kind=DBL_KIND),dimension(0:NREACT-1),intent(IN),optional :: xk_in !高速化のため反応係数を外から受け取れるようにする

    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: y_init,y_tmp,dy,ddy
    real(kind=DBL_KIND),dimension(0:NREACT-1) :: xk
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: r_f,r_f_fw,r_f_bw
    real(kind=DBL_KIND),dimension(0:NCHEM-1,0:NCHEM-1) :: dr_fdy,A
    real(kind=DBL_KIND),dimension(0:NCHEM-1,0:NCHEM-1) :: dr_fdy_dbg !for debug
    integer :: isp,jsp,ia,info,itr,iguess
    real(kind=DBL_KIND) :: delta_y,dr_f,err,err_max,r_f_big,tch

    !parameters controlling iteration
    ! real(kind=DBL_KIND),parameter :: eps=1.d-4        !relaitve displacement
    ! real(kind=DBL_KIND),parameter :: eps_y=1.d-10     !absolute displacement
    real(kind=DBL_KIND),parameter :: eps=1.d-5        !relaitve displacement
    real(kind=DBL_KIND),parameter :: eps_y=1.d-10     !absolute displacement    
    real(kind=DBL_KIND),parameter :: eps_conv = 1.d-5 !1.d-5 !1.d-8
    integer,parameter :: itrmax   = 30    !10
    !variables/parameters for dgesv
    integer,parameter :: inc=1, n=NCHEM
    integer :: ipiv(0:NCHEM-1)
    logical :: nan_flg

    !initialize
    y_init(:) = y(:)

    !外からxk_inが与えられていたらそれをxkに代入して終わり
    if (present(xk_in)) then
       xk(:) = xk_in(:)
    else
       !与えられていなければここで計算
       !get xk(:) from xnH, T_K, and radiation
       call react_coef(xnH,T_K,rHpi,rH2pd,rHmpd, Tdust, fd, xk)
    end if
    

    do iguess=0,2
       if (iguess==0) then
          y(:) = y_init(:) ! initial guess of previous abundance
       else if (iguess==1) then
          y(:) = 1d-12 ! initial guess of almost dissociated gas
          y(X_HI) = 1d0
          call adjust_abundance(y,yco,metal)
       else if (iguess==2) then
          y(:) = 1d-12 ! initial guess of almost ionized gas
          y(X_ele) = 1d0
          y(X_HII) = 1d0
          call adjust_abundance(y,yco,metal)          
       end if
       dy(:) = y(:) - y_init(:)             
       nan_flg = .false.

       !begin iteration
       do itr=0,itrmax-1

          !reaction rate@current position (y)
          ! call react_rat(xk,xnH,y,r_f)
          call react_rat_highspeed(xk,xnH,y,r_f)

          !KS DEBUG
          !          if (dbg_flg_prm == 1 .or. itr > itrmax-3) then
          if (dbg_flg_prm == 1) then
             print '(/,A,I0,A,(1P10E15.7))', "itr: ",itr, ", xnH T_K rHpi rH2pd rHmpd dt: ",&
                  xnH,T_K,rHpi,rH2pd,rHmpd,dt
             print '(A,(1P10E15.7))', "y_init: ",y_init(:)
             print '(A,(1P10E15.7))', "y:",y(:)
             print '(A,(1P10E15.7))', "r_f:",r_f(:)
          end if



          !================= KS MODIFIED for SPEED-UP ================!

          !get dr/dy by using analytical expressions
          call react_drdy(xk, xnH, y, dr_fdy)


          ! !get Jacobian dr/dy below
          ! do jsp=0,NCHEM-1

          !    !define delta_y(j)
          !    !             if(dabs(y(jsp)).le.1.d-10) then
          !    if(dabs(y(jsp)) .le. eps_y/eps) then ! KS TODO
          !       delta_y = eps_y
          !    else
          !       delta_y = eps*y(jsp)
          !    endif

          !    !KS MODIFIED
          !    y_tmp(:) = y(:)
          !    ! do isp=0,NCHEM-1
          !    !    if(isp.eq.jsp) then
          !    !       y_tmp(isp) = y(isp) + 5.d-1*delta_y
          !    !    else
          !    !       y_tmp(isp) = y(isp)
          !    !    endif
          !    ! enddo

          !    !reaction rate@forward displacement (y+delta_y)
          !    y_tmp(jsp) = y(jsp) + 5.d-1*delta_y
          !    call react_rat(xk,xnH,y_tmp,r_f_fw)

          !    !reaction rate@backward dispalcement (y-delta_y)
          !    y_tmp(jsp) = y(jsp) - 5.d-1*delta_y
          !    call react_rat(xk,xnH,y_tmp,r_f_bw)

          !    ! !--------------- KS DEBUG ----------------!
          !    ! print *, "jsp", jsp
          !    ! print *, "r_f_fw", r_f_fw
          !    ! print *, "r_f_bw", r_f_bw
          !    ! !--------------- KS DEBUG ----------------!

          !    !derivative of reaction rate
          !    do isp=0,NCHEM-1
          !       dr_fdy(isp,jsp)= (r_f_fw(isp) - r_f_bw(isp))/delta_y

          !       r_f_big = max(dabs(r_f_fw(isp)),dabs(r_f_bw(isp)))
          !       if(r_f_big.ne.0.d0) then
          !          dr_f = dabs(r_f_fw(isp)-r_f_bw(isp))
          !          if(dr_f/r_f_big.lt.1.d-15) then
          !             dr_fdy(isp,jsp) = 0.d0
          !          endif
          !       endif
          !    enddo
          ! enddo

          ! !get dr/dy by using analytical expressions
          ! call react_drdy(xk, xnH, y, dr_fdy_dbg)

 !          print *, "numerical derivative"
 !          do jsp=0,NCHEM-1
 !             print *, dr_fdy(:,jsp)
 !          end do
 !          print *, "analytical derivative"
 !          do jsp=0,NCHEM-1
 !             print *, dr_fdy_dbg(:,jsp)
 !          end do
 !          print *, "(analytical - numerical)/numerical derivative"
 !          do jsp=0,NCHEM-1
 !             print *, (dr_fdy_dbg(:,jsp) - dr_fdy(:,jsp))/abs(dr_fdy(:,jsp))
 !          end do
 !          print *, " numerical and analytical conservation check (Hnuc and charge)"
 !          do jsp=0,NCHEM-1
 !             print *, dr_fdy_dbg(0,jsp)+2*dr_fdy_dbg(1,jsp)+dr_fdy_dbg(3,jsp)+dr_fdy_dbg(4,jsp)+2*dr_fdy_dbg(5,jsp), &
 !                  -dr_fdy_dbg(2,jsp)+dr_fdy_dbg(3,jsp)-dr_fdy_dbg(4,jsp)+dr_fdy_dbg(5,jsp),&
 ! dr_fdy(0,jsp)+2*dr_fdy(1,jsp)+dr_fdy(3,jsp)+dr_fdy(4,jsp)+2*dr_fdy(5,jsp), &
 !                  -dr_fdy(2,jsp)+dr_fdy(3,jsp)-dr_fdy(4,jsp)+dr_fdy(5,jsp)             
 !          end do                    
          ! stop
          !================= KS MODIFIED for SPEED-UP ================!         

          ! -------------- SET matrix A -----------------
          do isp=0,NCHEM-1
             do jsp=0,NCHEM-1
                if(isp.ne.jsp) then
                   A(isp,jsp) = -dt*dr_fdy(isp,jsp)
                else
                   A(isp,jsp) = 1.d0-dt*dr_fdy(isp,jsp)
                endif
             enddo

             !            write(*,'(a,I5,1p4e13.5)') 'Matrix A:', &
             !                         isp, (A(isp,jsp),jsp=1,4) 
          enddo
          ! ----------------------------------------------

          ! -------- SET vector ddy ----------
          do isp=0,NCHEM-1
             ddy(isp) = r_f(isp)*dt-dy(isp)
          enddo
          ! -----------------------------------

          ! !--------------- KS DEBUG ----------------!
          ! do isp=0,NCHEM-1
          !    do jsp=0,NCHEM-1
          !       print *, "A(isp,jsp)=",isp, jsp, A(isp,jsp)
          !    end do
          ! end do
          ! do isp=0,NCHEM-1
          !    print *, "ddy(isp) = ", isp, ddy(isp)
          ! enddo
          ! !--------------- KS DEBUG ----------------!


          ! -------- Solve linear equation ----------
          call dgesv(n,inc,A,n,ipiv,ddy,n,info) !dgesv from LAPCK library
          ! call gaussj(A,N_sp,N_sp,ddy,1,1) !Gauss-Jordan method
          ! ----------------------------------------------

          ! !--------------- KS DEBUG ----------------!
          ! do isp=0,NCHEM-1
          !    print *, "ddy new (isp) = ", isp, ddy(isp)
          ! enddo
          ! !--------------- KS DEBUG ----------------!

          if(itr < 15) then
            do isp=0,NCHEM-1
               dy(isp) = dy(isp) + ddy(isp)
               y(isp)  = y(isp)  + ddy(isp) !!update!!
            enddo
          else
            do isp=0,NCHEM-1
               dy(isp) = dy(isp) + ddy(isp)*0.5d0
               y(isp)  = y(isp)  + ddy(isp)*0.5d0 !!update!!
            enddo
          endif


          !---------- 物理的にありえない値にならないよう調整 ----------
          ! KS TODO: 値の調整をすべきか？BH降着計算では調整してた
          if (dbg_flg_prm == 1) then
             print '(A,(1P10E15.7))', "(Before_adjust) y:",y(:) !KS DEBUG
          end if
          call adjust_abundance2(y,metal)
          dy(:) = y(:) - y_init(:) !yをADJUSTしたらdyにも反映させる
          if (dbg_flg_prm == 1) then
             print '(A,(1P10E15.7))', "(After_adjust) y:",y(:) !KS DEBUG
          end if

          !------------  NaN/Inf/Huge check -> あったら次のinitial guessへ  --------!
          do isp=0,NCHEM-1
             ! if (y(isp)/=y(isp) .or. y(isp)*0/=0) then
             if (y(isp)/=y(isp) .or. y(isp)*0/=0 .or. y(isp)==Huge(1d0)) then !Note: adjustしきれないときHugeが入る
                nan_flg = .true. 
                print *, "(chemreact) NaN/Inf/Huge found, move to next initial guess"
                exit
             end if
          end do
          if (nan_flg) exit ! exit from itr do-loop

          !----------------   iterationが収束したかの判定   -------------------!
          err_max = 0.d0
          do isp=0,NCHEM-1
             if (y(isp) == 0d0) then
                err_max=Huge(1d0) ! zeroがあったら収束させない
                exit
             end if
             !if(y(isp) > 1d-8) then ! 組成比が小さい分子は判定に使わない（まずはこの条件使わないでやってみて様子を見る KS TODO）
             if(y(isp) > 1d-15 .or. y_init(isp) > 1d-15) then ! 組成比が小さい分子は判定に使わない（しばらくこの条件でやって様子を見る、doubleの桁数〜15 KS TODO）
             ! if(y(isp) > 1d-30 .or. y_init(isp) > 1d-30) then ! 組成比が小さい分子は判定に使わない（しばらくこの条件でやって様子を見る KS TODO）
                err = dabs(ddy(isp)/y(isp))
                err_max = max(err,err_max)
             end if   
         enddo
          if (dbg_flg_prm == 1) then
             print '(A,(1P10E15.7))', "(error_y, error_max) y:",ddy(:)/y(:), err_max !KS DEBUG
          end if

          !------------  収束望み薄 -> 次のinitial guessへ  --------!
         if (itr >= 2) then !iteration 3回目になってもerr_maxが1より大きかったら望み薄
            if (err_max > 1d0) then
               if (dbg_flg_prm == 1) then
                  print '(A,1P1E15.7,A,I0,A)', "(chemreact) err_max =", err_max, " (>1) at ",&
                       itr, "th iteration, move to next initial guess"
               end if
               exit !exit from itr do-loop
            end if
         end if
         



          ! ----------- TEST OUTPUT -------------
          !         write(*,'(I5,1p1e13.5)') itr,err_max 
          ! -------------------------------------

          !収束判定とtchemの評価
          if(err_max.lt.eps_conv) then
             tchem = HUGE(1d0)
             do isp=0,NCHEM-1
                if (y(isp) /= y_init(isp)) then 
                   !                if (isp == 4 .or. isp == 5) then !HmとH2pはabundanceが小さいのでt_chemに影響しないようにする (KS MEMO)
                   if ((y(isp)+y_init(isp))/2. < 1d-5) then !少量の分子はt_chemの見積もりに使わない (KS TODO)
                      cycle
                   endif
                   tch = dabs((y(isp) + y_init(isp)) &
                        /(2.d0*(y(isp) - y_init(isp))))*dt  
                   tchem = min(tch,tchem)
                end if
             enddo
             return
          endif
          

          !収束しなかった場合のメッセージ出力
          if(itr.eq.itrmax-1 .and. dbg_flg_prm == 1) then 
             if (iguess==0) then
                print *, "chemreact: iteration not converged (initial guess of previous abundance)"
                print *, "try initial guess of fully dissociated gas"
             else if (iguess==1) then
                print *, "chemreact: iteration not converged (initial guess of fully dissociated gas)"
                print *, "try initial guess of fully ionized gas"
             else if (iguess==2) then
                print *, "chemreact: iteration not converged (initial guess of fully ionized gas)"
                print *, "try another timestep update"
             end if
          end if
       end do!itr
    end do !iguess
    !最後まで収束しなかった場合
    tchem = -1d0 !negative value meaning non-convergence

  end subroutine chemreact

  !-----------------------------------------------------------------------
  ! chemreact wrapper
  !
  ! dt が tchemより大きすぎてimplicit法でも収束しない場合に、
  ! substepを切って一旦収束させてから時間発展を解くことで対応
  !-----------------------------------------------------------------------  
  subroutine chemreact_adptv(xnH,T_K,Td,y,yco,rHpi,rH2pd,rHmpd,fd,metal,dt,tchem,force_substep,xk_in)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,rHpi,rH2pd,rHmpd,dt,Td,fd,metal
    real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1),yco
    real(kind=DBL_KIND),intent(OUT) :: tchem
    logical,intent(IN),optional :: force_substep
    real(kind=DBL_KIND),dimension(0:NREACT-1),intent(IN),optional :: xk_in !高速化のため反応係数を外から受け取れるようにする    

    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: y_tmp
    real(kind=DBL_KIND) :: t_sub, dt_sub
    integer :: i,n
    integer, parameter :: imax = 15 !10 !substepの時間幅を短くする回数のmax
    integer, parameter :: nmax = 100 !20 !substep数のmax
    logical :: bool_force_substep

    !For debug
    real(kind=DBL_KIND) :: y_o(0:NCHEM-1), y_dbg(0:20,0:NCHEM-1), dt_sub_n(0:20)
    integer :: nn, i_itr_n(0:20)

    !KS DEBUG
    y_o(:) = y(:)
    y_dbg(:,:) = -1d10
    dt_sub_n(:) = -1d10
    i_itr_n(:) = -100

    !強制的にsubstepを切るか
    bool_force_substep = .false.
    if ( present( force_substep) ) then
       bool_force_substep = force_substep
    endif


    t_sub=0d0 !substepのt

    do n=0,nmax
       do i=0,imax
          y_tmp(:)=y(:)
          dt_sub = (dt-t_sub)*10.d0**(-i) !substepのdelta_t
          if (n==0 .and. i==0 .and. bool_force_substep) then
             dt_sub = 0.5d0 * dt_sub   !bool_force_substepがtrueのときは、強制的に少なくとも2回以上substepを経るようにする
          end if

          !外からxk_inが与えられていたら、chemreactにも渡す
          if (present(xk_in)) then
             call chemreact(xnH,T_K,Td,y_tmp,yco,rHpi,rH2pd,rHmpd,fd,metal,dt_sub,tchem,xk_in=xk_in) ! Tdは一つ前で計算されているのでよい?HFADDED
          else
             !そうでなければ通常のchemreactを実行
             call chemreact(xnH,T_K,Td,y_tmp,yco,rHpi,rH2pd,rHmpd,fd,metal,dt_sub,tchem)
          end if

          !KS DEBUG
          i_itr_n(n) = i

          if (tchem>0d0) exit !tchem > 0 if chemreact converged

          !収束しなかった場合
          if (i==imax) then
             print *, "(chemreact_adptv), i == imax ", i, imax
             print '(2(A,1P1E15.7))', "(chemreact_adptv), chemreact not converged w/ dt_sub=",dt_sub, ", dt=", dt
             print '(/,A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt_sub: ",&
                  xnH,T_K,rHpi,rH2pd,rHmpd,dt_sub
             print '(A,(1P10E15.7))', "y_o: ",y_o(:)
             print *, "stopping..."
             stop
          end if
       end do !i

       !値を更新
       y(:) = y_tmp(:)
       t_sub = t_sub + dt_sub

       !KS DEBUG
       y_dbg(n,:) = y_tmp(:)
       dt_sub_n(n) = dt_sub

       !dtに到達したら終了
       if ((dt-t_sub) < 1d-3*dt) exit

       !dtまで到達できなかった場合
       if (n==nmax) then
          print *, "(chemreact_adptv), n == nmax ", n, nmax
          print '(3(A,1P1E15.7))', "(chemreact_adptv), max itr. reaches before reaching dt w/ t_sub=",&
               t_sub, ", dt=", dt, ", dt_sub=", dt_sub 
          print '(/,A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt_sub: ",&
               xnH,T_K,rHpi,rH2pd,rHmpd,dt_sub
          print '(A,(1P10E15.7))', "y_o: ",y_o(:)
          do nn=0, nmax
             print '(A,2I4, (1P11E15.7))', "n, i_itr_n, dt_sub_n, y_n: ",nn,i_itr_n(nn),dt_sub_n(nn), y_dbg(nn,:)
          end do
          print *, "stopping..."
          stop
       end if
    end do !n
    
  end subroutine chemreact_adptv

  !-----------------------------------------------------------------------
  ! abundance がおかしな値にならないよう調整する
  !-----------------------------------------------------------------------  
  subroutine adjust_abundance(y, yco, metal)
    real(kind=DBL_KIND),intent(IN) :: metal
    real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1)    
    real(kind=DBL_KIND),intent(INOUT) :: yco
    real(kind=DBL_KIND) ::  yHtot, chrgtot, nchrg, pchrg, pchrg_o, pchrg_c;
    real(kind=DBL_KIND), parameter :: min_value = 1d-50 !復帰する際の値
    real(kind=DBL_KIND) :: frac_C, frac_O

    ! set frac_C, frac_O
    frac_C = metal*frac_C_solar
    frac_O = metal*frac_O_solar


    !物理的にありえない値を取らないように制限
    y(0) = MAX(MIN(1.d0,y(0)),min_value)   !H
    y(1) = MAX(MIN(0.5d0,y(1)),min_value)  !H2
    y(2) = MAX(MIN(1.d0,y(2)),min_value)   !e
    y(3) = MAX(MIN(1.d0,y(3)),min_value)   !H+
    y(4) = MAX(MIN(1.d0,y(4)),min_value)   !H-
    y(5) = MAX(MIN(0.5d0,y(5)),min_value)  !H2+
    yco  = MAX(MIN(frac_C,yco),min_value)  !CO

    ! y(0) = MAX(MIN(1.d0,y(0)),0d0)   !H
    ! y(1) = MAX(MIN(0.5d0,y(1)),0d0)  !H2
    ! y(2) = MAX(MIN(1.d0,y(2)),0d0)   !e
    ! y(3) = MAX(MIN(1.d0,y(3)),0d0)   !H+
    ! y(4) = MAX(MIN(1.d0,y(4)),0d0)   !H-
    ! y(5) = MAX(MIN(0.5d0,y(5)),0d0)  !H2+

    
    !保存則が破れないように調整
    !------------------------- total H-abundance ---------------------!
    yHtot = y(0)+2.0*y(1)+y(3)+y(4)+2.0*y(5);
    if(yHtot == 0.d0) then
       print *, "(adjust_abundance) yHtot is zero -> unable to adjust"
       print *, "*** WARNING ***    Fill Huge(1d0)    *** WARNING ***"
       y(:)=Huge(1d0)
       return
    endif
    if(abs(yHtot-1) > 1.d-10) then
      !水素原子核の合計が保存していなかったら水素原子核を含む粒子の数を均等に増減させる
      y(0)=y(0)/yHtot; y(1)=y(1)/yHtot; y(3)=y(3)/yHtot; y(4)=y(4)/yHtot; y(5)=y(5)/yHtot;
   end if
    !-------------------------  total charge  -----------------------!
    pchrg = y(3) + y(5)
    nchrg = y(2) + y(4)

    ! heavy elements
    pchrg_o = y(3) * ( frac_O - yco ) ! assuming ionization rate is same with HII
    pchrg_c = (frac_C - yco) * 2.d0   ! all atomic carbon becomes CII
    pchrg   = pchrg + pchrg_o + pchrg_c
    

    chrgtot = 2 * (pchrg-nchrg) / (nchrg+pchrg)
    if(abs(chrgtot) > 1.d-10) then
       !電荷の合計が保存していなかったら電子の数を増減させる
       y(2) = MAX(min_value, pchrg - y(4));
       !それでもマイナス電荷が多かったら、HmをHnに置き換える
       y(0) = y(0) + (y(2) + y(4) - pchrg)
       y(4) = pchrg - y(2)
   end if
 end subroutine adjust_abundance

  !-----------------------------------------------------------------------
  ! abundance がおかしな値にならないよう調整する
  !-----------------------------------------------------------------------  
  subroutine adjust_abundance2(y,metal)
    real(kind=DBL_KIND),intent(IN)::metal
    real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1)    
    real(kind=DBL_KIND) ::  yHtot, chrgtot, nchrg, pchrg;
    real(kind=DBL_KIND), parameter :: min_value = 1d-50 !復帰する際の値
    real(kind=DBL_KIND) :: frac_C, frac_O

    ! set frac_C, frac_O
    frac_C = metal*frac_C_solar
    frac_O = metal*frac_O_solar

    !物理的にありえない値を取らないように制限
    y(0) = MAX(MIN(1.d0,y(0)),min_value)   !H
    y(1) = MAX(MIN(0.5d0,y(1)),min_value)  !H2
    y(2) = MAX(MIN(1.d0+frac_C*2.d0+frac_O,y(2)),min_value)   !e
    y(3) = MAX(MIN(1.d0,y(3)),min_value)   !H+
    y(4) = MAX(MIN(1.d0,y(4)),min_value)   !H-
    y(5) = MAX(MIN(0.5d0,y(5)),min_value)  !H2+

    ! y(0) = MAX(MIN(1.d0,y(0)),0d0)   !H
    ! y(1) = MAX(MIN(0.5d0,y(1)),0d0)  !H2
    ! y(2) = MAX(MIN(1.d0,y(2)),0d0)   !e
    ! y(3) = MAX(MIN(1.d0,y(3)),0d0)   !H+
    ! y(4) = MAX(MIN(1.d0,y(4)),0d0)   !H-
    ! y(5) = MAX(MIN(0.5d0,y(5)),0d0)  !H2+

    
    !保存則が破れないように調整
    !------------------------- total H-abundance ---------------------!
    yHtot = y(0)+2.0*y(1)+y(3)+y(4)+2.0*y(5);
    if(yHtot == 0.d0) then
       print *, "(adjust_abundance) yHtot is zero -> unable to adjust"
       print *, "*** WARNING ***    Fill Huge(1d0)    *** WARNING ***"
       y(:)=Huge(1d0)
       return
    endif
    if(abs(yHtot-1) > 1.d-10) then
      !水素原子核の合計が保存していなかったら水素原子核を含む粒子の数を均等に増減させる
      y(0)=y(0)/yHtot; y(1)=y(1)/yHtot; y(3)=y(3)/yHtot; y(4)=y(4)/yHtot; y(5)=y(5)/yHtot;
   end if
    !-------------------------  total charge  -----------------------!
    !pchrg = y(3) + y(5)
    !nchrg = y(2) + y(4)
    !chrgtot = 2 * (pchrg-nchrg) / (nchrg+pchrg)
    !if(abs(chrgtot) > 1.d-10) then
    !   !電荷の合計が保存していなかったら電子の数を増減させる
    !   y(2) = MAX(min_value, pchrg - y(4));
    !   !それでもマイナス電荷が多かったら、HmをHnに置き換える
    !   y(0) = y(0) + (y(2) + y(4) - pchrg)
    !   y(4) = pchrg - y(2)
    !end if
 end subroutine adjust_abundance2
  

  !-----------------------------------------------------------------------
  ! giving chemical reaction coefficients
  !-----------------------------------------------------------------------  
  subroutine react_coef(xnH,T_K,rHpi,rH2pd,rHmpd, Tdust, fd, xk) 
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,rHpi,rH2pd,rHmpd,Tdust, fd
    real(kind=DBL_KIND),intent(OUT) :: xk(0:NREACT-1)

    integer :: i
    real(kind=DBL_KIND) :: T_eV, xlnT_eV, log_T, log_T2, log_T3,&
       xk_L, xk_H, rlgT4, rn_cr, a, xkcid, xkdt;
    real(kind=DBL_KIND) :: f_a, r_a, inv_T, sqrt_T

    !xlnT_eV, log_T
    T_eV     = 8.61735d-5*T_K;
    xlnT_eV = dlog(T_eV);
    log_T  = dlog10(T_K);
    log_T2 = log_T*log_T;
    log_T3 = log_T2*log_T;
    inv_T  = 1.d0/T_K
    sqrt_T = dsqrt(T_K)

    ! 0)   H     +   e     ->   H+    + 2 e   
    !  Abel et al. (1997)
    xk(0) = dexp(-32.71396786 + (13.536556 &
         + (-5.73932875 + (1.56315498 + (-0.2877056 &
         + (3.48255977d-2 + (-2.63197617d-3 &
         + (1.11954395d-4 - 2.03914985d-6*xlnT_eV)*xlnT_eV) &
         *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)

    !----------------- KS DEBUG ---------------------!
    !xk(0)=0d0
    !----------------- KS DEBUG ---------------------!    

    !  1)   H+    +   e     ->   H     +   ph.
    ! Glover and Jappsen (Case A, B; 2007) based on Ferland et al. 1992

#ifdef RECOM_ION_PHOTO
    xk(1) = 1.269d-13 * (315614.d0*inv_T)**1.503 *(1+(604625.d0*inv_T)**0.470)**(-1.923)  !Case A
#else
    xk(1) = 2.753d-14 * (315614.d0*inv_T)**1.5   *(1+(115188.d0*inv_T)**0.407)**(-2.242)  !Case B
#endif


    ! 2)   H-    +   H     ->   H2    +   e  !<8>!            
    ! Kreckel et al. (2010)
    xk(2) = 1.35d-9* (T_K**0.098493 + 0.32852*T_K**0.5561 + 2.771d-7*T_K**2.1826) &
         /(1.0 + 6.191d-3*T_K**1.0461 + 8.9712d-11*T_K**3.0424 + 3.2576d-14*T_K**3.7741)

    !  3)   H2    +   H+    ->   H2+   +   H  !<11>!       
    !     Galli & Palla (1998) (H15)
    if(T_K <= 1.e4) then
       xk(3) = 3.d-10*dexp(-2.1050d4*inv_T);
    else
       xk(3) = 1.5d-10*dexp(-1.4000d4*inv_T);
    end if

    !  4)   H2    +   e     -> 2 H     +   e  !<12>!
    !     Galli & Palla (1998) (H17)
    xk(4) = 4.4d-10*(T_K**3.5d-1)*dexp(-1.02000d5*inv_T)

    !  5)   H2    +   H     -> 3 H            !<13>!
    !     Martin, Schwarz, Mandy (1996)
    xkcid = xkcidM96(xnH, T_K, log_T, log_T2, log_T3)
    xkdt  = xkdtM96 (xnH, T_K, log_T, log_T2, log_T3)
    xk(5) = xkcid + xkdt

    !   6) 3 H               ->   H2    +   H    !<19>!          
    !        Forrey (2013)
    xk(6) = 6.d-32*T_K**(-0.25) + 2.d-31*T_K**(-0.5)

    !   7) 2 H     +   H2    -> 2 H2             !<20>!     
    !     Palla, Salpeter & Stahler (1983)
    xk(7) = xk(6)*0.125d0

    !   8) 2 H2              -> 2 H     +   H2   !<21>!           
    !     Palla, Salpeter & Stahler (1983)
    xk_L = 1.18d-10*dexp(-6.95d+4*inv_T)
    xk_H = 8.125d-8/sqrt_T*dexp(-5.2e+4*inv_T)*(1.0 - dexp(-6.d3*inv_T))
    rlgT4 = log_T - 4.d0  !dlog10(T_K*1.d-4);
    rn_cr = 1.d1**(4.845d0 - 1.3d0*rlgT4 + 1.62d0*rlgT4**2)
    a = 1.d0/(1.d0 + xnH/rn_cr)
    if(a.eq.1.d0) then
       xk(8) = xk_L
    elseif(a.eq.0.d0) then
       xk(8) = xk_H
    else
       xk(8) = xk_H**(1.d0-a)*xk_L**a
    endif

    !  9)   H     +   e     ->   H-    +   ph. !< only reaction rate, 7>!
    !     Galli & Palla (1998) (H3)
    xk(9) = 1.4d-18*(T_K**9.28d-1)*dexp(-T_K/1.62d4)

    ! 10)   H     + ph.(uv) ->  H+   +  e  (photo-ionization)
    xk(10) = rHpi

    ! 11)   H2    + ph.(fuv) ->  2H   (photo-dissociation)
    xk(11) = rH2pd

    ! 12)   2 H               ->   H+    +   e     +   H !<22>! 
    xk(12) = 1.7d-4*xk(0);


    !  13)   H-    +   e     ->   H     + 2 e  !< only reaction rate, 14>! 
    !       Abel et al. (1997)   (14)
    xk(13) = dexp(- 18.01849334d0 + (2.3608522d0 &
         + (-0.28274430d0+(1.62331664d-2+(-3.36501203d-2 &
         + (1.17832978d-2+(-1.65619470d-3+(1.06827520d-4 &
         - 2.63128581d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV) &
         *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)

    ! 14) H-    +   H+    ->   H2+   +   e  !<16>! 
    !      Galli & Palla (1998) (H6)
    if(T_K.le.8.d3) then
       xk(14) = 6.9d-9/(T_K**3.5d-1)
    else
       xk(14) = 9.6d-7/(T_K**9.d-1)
    endif

    !  15)   H-    +   H+    ->  2H  
    !      Galli & Palla (1998) (H7)
    xk(15) = 6.3d-8 + 5.7d-6/sqrt_T - 9.2d-11*sqrt_T  + 4.4d-13*T_K

    !  16)  H- + ph. -> H + e
    xk(16) = rHmpd

    !  17) H     +   H+    ->   H2+   +   ph.    
    !      Galli & Palla (1998) (H8)
    xk(17) = 1.d1**(-19.38d0 - 1.523d0*log_T &
         + 1.118d0*log_T**2 - 0.1269d0*log_T**3)


    !  18)  H2+   +   H     ->   H2    +   H+         
    !       Galli & Palla (1998) (H10)
    xk(18) = 6.4d-10;

    !  19)   H2+   +   e     -> 2 H                
    xk(19) = 2.d-7/sqrt_T;

    !  20)   H2+   +   H-    ->   H2    +   H         
    !        Millar et al. (1991) 
    xk(20) = 39.8371685741d-7/sqrt_T      !sqrt(T_K/300.0);

    !  21)   2 H  +  grain  ->    H2  
#ifdef METAL

#if DUST_GAS_COUPLING == YES
    f_a=1.d0/(1.d0+dexp(7.5d2*(1.d0/75.d0-1.d0/Tdust))) 
    r_a=0.3464101615d-17*sqrt_T*f_a/(1.d0+4.0d-2*dsqrt(T_K+Tdust)+2.0d-3*T_K+ &
      8.0d-6*T_K*T_K)*fd  ! ダスト数密度を考慮する場合はここに加える HFADDED
    xk(21) = r_a
#else
    xk(21) = 0.d0
#endif
#endif

  end subroutine react_coef

  !---------------------------------------------------------------------------------
  ! give chemical reaction rate r_f_tot(i), which is defined as dy(i)/dt = r_f_tot(i)
  !---------------------------------------------------------------------------------
  subroutine react_rat(xk, xnH, y, r_f_tot)
    real(kind=DBL_KIND),intent(IN) :: xk(0:NREACT-1),xnH,y(0:NCHEM-1)
    real(kind=DBL_KIND),intent(OUT) :: r_f_tot(0:NCHEM-1)

    real(kind=DBL_KIND) :: r_f(0:NCHEM-1,0:NREACT-1)
    integer :: isp, ire
    real (kind=DBL_KIND) :: yhmn, rate, rrate

    !initialization
    r_f(:,:) = 0.d0
    r_f_tot(:) = 0.d0    

    !   0)   H     +   e     ->   H+    + 2 e 
    !       (0         2          3       2*2) 
    rate = xk(0)*y(0)*y(2)*xnH
    r_f(0,0) = -rate
    r_f(2,0) =  rate
    r_f(3,0) =  rate

    !   1)   H+    +   e     ->   H     +   ph.
    !       (3         2          0)
    rate = xk(1)*y(3)*y(2)*xnH
    r_f(3,1) = -rate
    r_f(2,1) = -rate
    r_f(0,1) =  rate

    !   2)   H-    +   H     ->   H2    +   e       
    !        (4        0          1         2)
    rate = xk(2)*y(4)*y(0)*xnH
    r_f(4,2) = -rate
    r_f(0,2) = -rate
    r_f(1,2) =  rate
    r_f(2,2) =  rate 

    !  3)   H2    +   H+    ->   H2+   +   H  
    !      (1         3          5         0)
    rate = xk(3)*y(1)*y(3)*xnH
    r_f(1,3) = -rate
    r_f(3,3) = -rate
    r_f(0,3) =  rate
    r_f(5,3) =  rate

    !  4)   H2    +   e     -> 2 H     +   e        
    !      (1         2        2*0         2)
    rate = xk(4)*y(1)*y(2)*xnH
    r_f(1,4) = -rate
    r_f(0,4) = 2.0*rate

    !  5)   H2    +   H     -> 3 H    
    !      (1         0        3*0)
    rate = xk(5)*y(1)*y(0)*xnH
    r_f(1,5) = -rate
    r_f(0,5) =  2.0*rate

    !  6) 3 H           ->   H2    +   H  
    !    (3*0                 1        0)
    rate = xk(6)*y(0)*y(0)*y(0)*xnH*xnH
    r_f(0,6) = -2.0*rate
    r_f(1,6) = rate

    !  7)  2 H     +   H2    -> 2 H2  
    !     (2*0         1        2*1)
    rate = xk(7)*y(0)*y(0)*y(1)*xnH*xnH 
    r_f(0,7) = -2.0*rate
    r_f(1,7) = rate

    !   8) 2 H2              -> 2 H     +   H2        
    !     (2*1                  2*0         1)
    rate = xk(8)*y(1)*y(1)*xnH
    r_f(1,8) = -rate
    r_f(0,8) =  2.0*rate

    !   9)  H     +   e     ->   H-    +   ph. 
    !      (0         2          4)
    rate = xk(9)*y(0)*y(2)*xnH
    r_f(2,9) = -rate
    r_f(0,9) = -rate
    r_f(4,9) =  rate


    !  10)   H     + ph.(uv) ->  H+   +  e  (photo-ionization)   
    !       (0                   3      2)
    !      also  H2    + ph.(uv) -> 2H+   + 2e 
    !        (   1                  2*3       2*2 )
    rate = xk(10)*y(0)
    r_f(0,10) = -rate
    r_f(2,10) =  rate
    r_f(3,10) =  rate

    !tentatively kill H2 photo-ionization (KS TODO)
    ! rate = xk(10)*y(0)
    ! r_f(2,10) = -rate
    ! r_f(2,10) =  2.*rate
    ! r_f(3,10) =  2.*rate
    
    ! rate = xk(10)*y(0)
    ! r_f(0,10) = -xk(10)*y(0)           
    ! r_f(1,10) = -xk(10)*y(1)  
    ! r_f(3,10) = xk(10)*(1.0 - y(3))   
    ! r_f(2,10) = xk(10)*(1.0 - y(2))   

    ! 11)  H2    + ph.(fuv) ->  2 H   (photo-dissociation)   
    !      (1                   2*0  )
    rate = xk(11)*y(1)
    r_f(1,11) = -rate
    r_f(0,11) =  2.0*rate

    ! 12) 2H ->  H+   +   e   +   H   
    !   ( 2*0     3       2       0  )
    rate = xk(12)*y(0)*y(0)*xnH 
    r_f(0,12) = -rate
    r_f(3,12) =  rate
    r_f(2,12) =  rate


    !  13)   H-    +   e     ->   H     + 2 e 
    !      ( 4         2          0       2 )
    rate = xk(13)*y(4)*y(2)*xnH
    r_f(4,13) = -rate
    r_f(0,13) =  rate
    r_f(2,13) =  rate

    !  14)  H-    +   H+    ->   H2+   +   e 
    !     ( 4         3          5         2 )
    rate = xk(14)*y(4)*y(3)*xnH
    r_f(4,14) = -rate
    r_f(3,14) = -rate
    r_f(5,14) =  rate
    r_f(2,14) =  rate

    !  15)   H-    +   H+    ->  2H
    !      ( 4         3          0 )
    rate = xk(15)*y(4)*y(3)*xnH
    r_f(4,15) = -rate
    r_f(3,15) = -rate
    r_f(0,15) = 2.0*rate

    !  16)   H-    +    ph   ->  H   +   e
    !      ( 4                   0       2 )
    rate = xk(16)*y(4)
    r_f(4,16) = -rate
    r_f(0,16) = rate
    r_f(2,16) = rate

    !   17) H  +  H+  ->  H2+  +  ph. 
    !     ( 0     3        5 )
    rate = xk(17)*y(0)*y(3)*xnH
    r_f(0,17) = -rate
    r_f(3,17) = -rate
    r_f(5,17) =  rate

    !  18) H2+  +  H  ->  H2  +  H+   
    !    ( 5       0      1      3 )
    rate = xk(18)*y(5)*y(0)*xnH
    r_f(5,18) = -rate
    r_f(0,18) = -rate
    r_f(1,18) =  rate
    r_f(3,18) =  rate

    !  19)   H2+   +   e     -> 2 H 
    !      ( 5         2          0 )
    rate = xk(19)*y(5)*y(2)*xnH
    r_f(5,19) = -rate
    r_f(2,19) = -rate
    r_f(0,19) = 2.0*rate

    !  20)   H2+   +   H-    ->   H2    +   H    
    !      ( 5         4           1        0 )
    rate = xk(20)*y(5)*y(4)*xnH
    r_f(5,20) = -rate
    r_f(4,20) = -rate
    r_f(1,20) =  rate  
    r_f(0,20) =  rate

    !  21)  2 H    + grain   ->  H2
    !     ( 0                    1  )
#ifdef METAL
    rate = xk(21)*y(0)*xnH
    r_f(0,21) = - 2.d0 * rate
    r_f(1,21) = rate
#endif
    !===================  summing up  =====================
    do isp=0,NCHEM-1
       do ire=0,NREACT-1
          r_f_tot(isp) = r_f_tot(isp) + r_f(isp,ire)
       end do       
    end do
    !=====================================================

    !KS DEBUG
    ! if (dbg_flg_prm == 1) then
    !    isp = 0 !Hn
    !    print '(A,I,/,1P1E15.7,A, (1P21E15.7))', "in react_rat: species = ",isp,&
    !         r_f_tot(isp), " <- ", r_f(isp,:)
    !    isp = 1 !H2
    !    print '(A,I,/,1P1E15.7,A, (1P21E15.7))', "in react_rat: species = ",isp,&
    !         r_f_tot(isp), " <- ", r_f(isp,:)
    !    isp = 4 !Hm
    !    print '(A,I,/,1P1E15.7,A, (1P21E15.7))', "in react_rat: species = ",isp,&
    !         r_f_tot(isp), " <- ", r_f(isp,:)       
    ! end if


  end subroutine react_rat

  !
  ! 高速化されたreact_rat
  ! 結果は同じだが、r_fを経由しない分だけ速い
  ! （r_fはNCHEM*NREACTの配列なので結構重い）
  !
  subroutine react_rat_highspeed(xk, xnH, y, r_f_tot)
    real(kind=DBL_KIND),intent(IN) :: xk(0:NREACT-1),xnH,y(0:NCHEM-1)
    real(kind=DBL_KIND),intent(OUT) :: r_f_tot(0:NCHEM-1)
    integer :: isp, ire
    real (kind=DBL_KIND) :: yhmn, rate, rrate

    !initialization
    r_f_tot(:) = 0.d0    

    !   0)   H     +   e     ->   H+    + 2 e 
    !       (0         2          3       2*2) 
    rate = xk(0)*y(0)*y(2)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(2) = r_f_tot(2) + rate
    r_f_tot(3) = r_f_tot(3) + rate


    !   1)   H+    +   e     ->   H     +   ph.
    !       (3         2          0)
    rate = xk(1)*y(3)*y(2)*xnH
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(0) = r_f_tot(0) + rate

    !   2)   H-    +   H     ->   H2    +   e       
    !        (4        0          1         2)
    rate = xk(2)*y(4)*y(0)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(1) = r_f_tot(1) + rate
    r_f_tot(2) = r_f_tot(2) + rate

    ! print *, "KS DEBUG 2: ", r_f_tot(1)
    

    !  3)   H2    +   H+    ->   H2+   +   H  
    !      (1         3          5         0)
    rate = xk(3)*y(1)*y(3)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(0) = r_f_tot(0) + rate
    r_f_tot(5) = r_f_tot(5) + rate

    ! print *, "KS DEBUG 3: ", r_f_tot(1)    

    !  4)   H2    +   e     -> 2 H     +   e        
    !      (1         2        2*0         2)
    rate = xk(4)*y(1)*y(2)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.0*rate

    ! print *, "KS DEBUG 4: ", r_f_tot(1)

    !  5)   H2    +   H     -> 3 H    
    !      (1         0        3*0)
    rate = xk(5)*y(1)*y(0)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) +  2.0*rate

    ! print *, "KS DEBUG 5: ", r_f_tot(1)

    !  6) 3 H           ->   H2    +   H  
    !    (3*0                 1        0)
    rate = xk(6)*y(0)*y(0)*y(0)*xnH*xnH
    r_f_tot(0) = r_f_tot(0) - 2.0*rate
    r_f_tot(1) = r_f_tot(1) + rate

    ! print *, "KS DEBUG 6: ", r_f_tot(1)

    !  7)  2 H     +   H2    -> 2 H2  
    !     (2*0         1        2*1)
    rate = xk(7)*y(0)*y(0)*y(1)*xnH*xnH 
    r_f_tot(0) = r_f_tot(0) - 2.0*rate
    r_f_tot(1) = r_f_tot(1) + rate

    ! print *, "KS DEBUG 7: ", r_f_tot(1)

    !   8) 2 H2              -> 2 H     +   H2        
    !     (2*1                  2*0         1)
    rate = xk(8)*y(1)*y(1)*xnH
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.0*rate

    ! print *, "KS DEBUG 8: ", r_f_tot(1)

    !   9)  H     +   e     ->   H-    +   ph. 
    !      (0         2          4)
    rate = xk(9)*y(0)*y(2)*xnH
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(4) = r_f_tot(4) + rate


    !  10)   H     + ph.(uv) ->  H+   +  e  (photo-ionization)   
    !       (0                   3      2)
    !      also  H2    + ph.(uv) -> 2H+   + 2e 
    !        (   1                  2*3       2*2 )
    rate = xk(10)*y(0)
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(2) = r_f_tot(2) + rate
    r_f_tot(3) = r_f_tot(3) + rate

    !tentatively kill H2 photo-ionization (KS TODO)
    ! rate = xk(10)*y(0)
    ! r_f_tot(2,10) = -rate
    ! r_f_tot(2,10) =  2.*rate
    ! r_f_tot(3,10) =  2.*rate
    
    ! rate = xk(10)*y(0)
    ! r_f_tot(0,10) = -xk(10)*y(0)           
    ! r_f_tot(1,10) = -xk(10)*y(1)  
    ! r_f_tot(3,10) = xk(10)*(1.0 - y(3))   
    ! r_f_tot(2,10) = xk(10)*(1.0 - y(2))   

    ! 11)  H2    + ph.(fuv) ->  2 H   (photo-dissociation)   
    !      (1                   2*0  )
    rate = xk(11)*y(1)
    r_f_tot(1) = r_f_tot(1) - rate
    r_f_tot(0) = r_f_tot(0) + 2.0*rate

    ! print *, "KS DEBUG 11: ", r_f_tot(1), xk(11), y(1), r_f_tot(1)+xk(11)*y(1) 

    ! 12) 2H ->  H+   +   e   +   H   
    !   ( 2*0     3       2       0  )
    rate = xk(12)*y(0)*y(0)*xnH 
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) + rate
    r_f_tot(2) = r_f_tot(2) + rate


    !  13)   H-    +   e     ->   H     + 2 e 
    !      ( 4         2          0       2 )
    rate = xk(13)*y(4)*y(2)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(0) = r_f_tot(0) + rate
    r_f_tot(2) = r_f_tot(2) + rate

    !  14)  H-    +   H+    ->   H2+   +   e 
    !     ( 4         3          5         2 )
    rate = xk(14)*y(4)*y(3)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(5) = r_f_tot(5) + rate
    r_f_tot(2) = r_f_tot(2) + rate

    !  15)   H-    +   H+    ->  2H
    !      ( 4         3          0 )
    rate = xk(15)*y(4)*y(3)*xnH
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(0) = r_f_tot(0) + 2.0*rate

    !  16)   H-    +    ph   ->  H   +   e
    !      ( 4                   0       2 )
    rate = xk(16)*y(4)
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(0) = r_f_tot(0) + rate
    r_f_tot(2) = r_f_tot(2) + rate

    !   17) H  +  H+  ->  H2+  +  ph. 
    !     ( 0     3        5 )
    rate = xk(17)*y(0)*y(3)*xnH
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(3) = r_f_tot(3) - rate
    r_f_tot(5) = r_f_tot(5) + rate

    !  18) H2+  +  H  ->  H2  +  H+   
    !    ( 5       0      1      3 )
    rate = xk(18)*y(5)*y(0)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(0) = r_f_tot(0) - rate
    r_f_tot(1) = r_f_tot(1) + rate
    r_f_tot(3) = r_f_tot(3) + rate

    ! print *, "KS DEBUG 18: ", r_f_tot(1), xk(18),y(5),y(0),xnH, r_f_tot(1)- xk(18)*y(5)*y(0)*xnH


    !  19)   H2+   +   e     -> 2 H 
    !      ( 5         2          0 )
    rate = xk(19)*y(5)*y(2)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(2) = r_f_tot(2) - rate
    r_f_tot(0) = r_f_tot(0) + 2.0*rate

    !  20)   H2+   +   H-    ->   H2    +   H    
    !      ( 5         4           1        0 )
    rate = xk(20)*y(5)*y(4)*xnH
    r_f_tot(5) = r_f_tot(5) - rate
    r_f_tot(4) = r_f_tot(4) - rate
    r_f_tot(1) = r_f_tot(1) + rate  
    r_f_tot(0) = r_f_tot(0) + rate

    !  21)  2 H    + grain   ->  H2
    !     ( 0                    1  )
#ifdef METAL
    rate = xk(21)*y(0)*xnH
    r_f_tot(0) = r_f_tot(0) - 2.d0*rate
    r_f_tot(1) = r_f_tot(1) + rate
#endif

    ! print *, "KS DEBUG 20: ", r_f_tot(1), xk(20),y(5),y(4),xnH, r_f_tot(1)-xk(20)*y(5)*y(4)*xnH

  end subroutine react_rat_highspeed


  !---------------------------------------------------------------------------------
  ! 計算時間の短縮のため、dr_f/dy を数値的に求めるのではなく直接求めてしまう
  !---------------------------------------------------------------------------------
  subroutine react_drdy(xk, xnH, y, dr_fdy)
    real(kind=DBL_KIND),intent(IN) :: xk(0:NREACT-1),xnH,y(0:NCHEM-1)
    real(kind=DBL_KIND),dimension(0:NCHEM-1,0:NCHEM-1),intent(OUT) :: dr_fdy
    real(kind=DBL_KIND) :: xnH_2, ddr_A, ddr_B, ddr_C

    xnH_2 = xnH*xnH

!A -> B (21)
#define DRDY_1T01_21(i,A,B) \
    dr_fdy(A,A) = dr_fdy(A,A) - 2.d0*xk(i)*xnH ;\
    dr_fdy(B,A) = dr_fdy(B,A) + xk(i)*xnH


!A -> B + C
#define DRDY_1TO2(i,A,B,C) \
    dr_fdy(A,A) = dr_fdy(A,A) - xk(i) ;\
    dr_fdy(B,A) = dr_fdy(B,A) + xk(i) ;\
    dr_fdy(C,A) = dr_fdy(C,A) + xk(i)    



!A+B -> C
#define DRDY_2TO1(i,A,B,C) \
    ddr_A = xk(i)*y(B)*xnH ;\
    ddr_B = xk(i)*y(A)*xnH ;\
    dr_fdy(A,A) = dr_fdy(A,A) - ddr_A ;\
    dr_fdy(B,A) = dr_fdy(B,A) - ddr_A ;\
    dr_fdy(C,A) = dr_fdy(C,A) + ddr_A ;\
    dr_fdy(A,B) = dr_fdy(A,B) - ddr_B ;\
    dr_fdy(B,B) = dr_fdy(B,B) - ddr_B ;\
    dr_fdy(C,B) = dr_fdy(C,B) + ddr_B


!A+B -> C+D
#define DRDY_2TO2(i,A,B,C,D) \
    ddr_A = xk(i)*y(B)*xnH ;\
    ddr_B = xk(i)*y(A)*xnH ;\
    dr_fdy(A,A) = dr_fdy(A,A) - ddr_A ;\
    dr_fdy(B,A) = dr_fdy(B,A) - ddr_A ;\
    dr_fdy(C,A) = dr_fdy(C,A) + ddr_A ;\
    dr_fdy(D,A) = dr_fdy(D,A) + ddr_A ;\
    dr_fdy(A,B) = dr_fdy(A,B) - ddr_B ;\
    dr_fdy(B,B) = dr_fdy(B,B) - ddr_B ;\
    dr_fdy(C,B) = dr_fdy(C,B) + ddr_B ;\
    dr_fdy(D,B) = dr_fdy(D,B) + ddr_B


!A+B -> C+D+E
#define DRDY_2TO3(i,A,B,C,D,E) \
    ddr_A = xk(i)*y(B)*xnH ;\
    ddr_B = xk(i)*y(A)*xnH ;\
    dr_fdy(A,A) = dr_fdy(A,A) - ddr_A ;\
    dr_fdy(B,A) = dr_fdy(B,A) - ddr_A ;\
    dr_fdy(C,A) = dr_fdy(C,A) + ddr_A ;\
    dr_fdy(D,A) = dr_fdy(D,A) + ddr_A ;\
    dr_fdy(E,A) = dr_fdy(E,A) + ddr_A ;\    
    dr_fdy(A,B) = dr_fdy(A,B) - ddr_B ;\
    dr_fdy(B,B) = dr_fdy(B,B) - ddr_B ;\
    dr_fdy(C,B) = dr_fdy(C,B) + ddr_B ;\
    dr_fdy(D,B) = dr_fdy(D,B) + ddr_B ;\    
    dr_fdy(E,B) = dr_fdy(E,B) + ddr_B


!A+B+C -> D+E
#define DRDY_3TO2(i,A,B,C,D,E) \
    ddr_A = xk(i)*y(B)*y(C)*xnH_2 ;\
    ddr_B = xk(i)*y(A)*y(C)*xnH_2 ;\
    ddr_C = xk(i)*y(A)*y(B)*xnH_2 ;\
    dr_fdy(A,A) = dr_fdy(A,A) - ddr_A ;\
    dr_fdy(B,A) = dr_fdy(B,A) - ddr_A ;\
    dr_fdy(C,A) = dr_fdy(C,A) - ddr_A ;\
    dr_fdy(D,A) = dr_fdy(D,A) + ddr_A ;\
    dr_fdy(E,A) = dr_fdy(E,A) + ddr_A ;\    
    dr_fdy(A,B) = dr_fdy(A,B) - ddr_B ;\
    dr_fdy(B,B) = dr_fdy(B,B) - ddr_B ;\
    dr_fdy(C,B) = dr_fdy(C,B) - ddr_B ;\
    dr_fdy(D,B) = dr_fdy(D,B) + ddr_B ;\
    dr_fdy(E,B) = dr_fdy(E,B) + ddr_B ;\    
    dr_fdy(A,C) = dr_fdy(A,C) - ddr_C ;\
    dr_fdy(B,C) = dr_fdy(B,C) - ddr_C ;\
    dr_fdy(C,C) = dr_fdy(C,C) - ddr_C ;\
    dr_fdy(D,C) = dr_fdy(D,C) + ddr_C ;\
    dr_fdy(E,C) = dr_fdy(E,C) + ddr_C
    

    !initialization
    dr_fdy(:,:) = 0.d0

    !   0)   H     +   e     ->   H+    + 2 e 
    !       (0         2          3       2*2) 
    DRDY_2TO3(0,0,2,3,2,2)

    !   1)   H+    +   e     ->   H     +   ph.
    !       (3         2          0)
    DRDY_2TO1(1,3,2,0)

    
    !   2)   H-    +   H     ->   H2    +   e       
    !        (4        0          1         2)
    DRDY_2TO2(2,4,0,1,2)

    !  3)   H2    +   H+    ->   H2+   +   H  
    !      (1         3          5         0)
    DRDY_2TO2(3,1,3,5,0)

    !  4)   H2    +   e     -> 2 H     +   e        
    !      (1         2        2*0         2)
    DRDY_2TO3(4,1,2,0,0,2)

    !  5)   H2    +   H     -> 3 H    
    !      (1         0        3*0)
    DRDY_2TO3(5,1,0,0,0,0)

    !  6) 3 H           ->   H2    +   H  
    !    (3*0                 1        0)
    DRDY_3TO2(6,0,0,0,1,0)

    !  7)  2 H     +   H2    -> 2 H2  
    !     (2*0         1        2*1)
    DRDY_3TO2(7,0,0,1,1,1)

    !   8) 2 H2              -> 2 H     +   H2        
    !     (2*1                  2*0         1)
    DRDY_2TO3(8,1,1,0,0,1)

    !   9)  H     +   e     ->   H-    +   ph. 
    !      (0         2          4)
    DRDY_2TO1(9,0,2,4)


    !  10)   H     + ph.(uv) ->  H+   +  e  (photo-ionization)   
    !       (0                   3      2)
    DRDY_1TO2(10,0,3,2)

    !      also  H2    + ph.(uv) -> 2H+   + 2e  <-  ignore in this code (KS TODO)
    !        (   1                  2*3       2*2 )

    ! 11)  H2    + ph.(fuv) ->  2 H   (photo-dissociation)   
    !      (1                   2*0  )
    DRDY_1TO2(11,1,0,0)

    ! 12) 2H ->  H+   +   e   +   H   
    !   ( 2*0     3       2       0  )
    DRDY_2TO3(12,0,0,3,2,0)

    !  13)   H-    +   e     ->   H     + 2 e 
    !      ( 4         2          0       2*2)
    DRDY_2TO3(13,4,2,0,2,2)

    !  14)  H-    +   H+    ->   H2+   +   e 
    !     ( 4         3          5         2 )
    DRDY_2TO2(14,4,3,5,2)

    !  15)   H-    +   H+    ->  2H
    !      ( 4         3          2*0 )
    DRDY_2TO2(15,4,3,0,0)

    !  16)   H-    +    ph   ->  H   +   e
    !      ( 4                   0       2 )
    DRDY_1TO2(16,4,0,2)

    !   17) H  +  H+  ->  H2+  +  ph. 
    !     ( 0     3        5 )
    DRDY_2TO1(17,0,3,5)

    !  18) H2+  +  H  ->  H2  +  H+   
    !    ( 5       0      1      3 )
    DRDY_2TO2(18,5,0,1,3)

    !  19)   H2+   +   e     -> 2 H 
    !      ( 5         2          2*0 )
    DRDY_2TO2(19,5,2,0,0)

    !  20)   H2+   +   H-    ->   H2    +   H    
    !      ( 5         4           1        0 )
    DRDY_2TO2(20,5,4,1,0)

    !  21)  2 H    + grain   ->  H2
    !     ( 0                    1  )
#ifdef METAL
    DRDY_1T01_21(21,0,1) 
#endif

  end subroutine react_drdy

  !-----------------------------------------------------------------------
  ! Martin+96, collision induced dissociation
  !-----------------------------------------------------------------------  
  function xkcidM96(xnH, T, logT, logT2, logT3)
    real(kind=DBL_KIND) :: xkcidM96
    real(kind=DBL_KIND),intent(IN) :: xnH, T, logT, logT2, logT3
    real(kind=DBL_KIND),dimension(21) :: a
    real(kind=DBL_KIND) :: log_g_h1, log_g_h2, log_g_l1, log_g_l2,&
         log_n_c1, n_c1, log_n_c2, n_c2, p, log_gamma_cd
    data a /-1.784239e2, -6.842243e1, 4.320243e1, -4.633167, 6.970086e1,&
         4.087038e4, -2.370570e4, 1.288953e2, -5.391334e1, 5.315517,&
!         -1.973427e1, 1.678905e4, -2.578611e4, 1.482123e1, -4.890915,&
         -1.973427e1, 1.678095e4, -2.578611e4, 1.482123e1, -4.890915,& !KS MODIFIED
         4.749030e-1, -1.338283e2, -1.164408, 8.227443e-1, 5.864073e-1, -2.056313/

  log_g_h1 = a(1) + a(2)*logT + a(3)*logT2 + a(4)*logT3 + a(5)*log10(1+a(6)/T)
  log_g_h2 = a(7)/T

  log_g_l1 = a(8) + a(9)*logT + a(10)*logT2 + a(11)*log10(1+a(12)/T)
  log_g_l2 = a(13)/T

  log_n_c1 = a(14) + a(15)*logT + a(16)*logT2 + a(17)/T
  n_c1 = 10.0**log_n_c1

  log_n_c2 = a(18) + log_n_c1
  n_c2 = 10.0**log_n_c2

  p = a(19) + a(20)*exp(-T/1850.0) + a(21)*exp(-T/440.0)

  log_gamma_cd =   log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+(xnH/n_c1)**p) &
                 + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+(xnH/n_c2)**p)
  
  xkcidM96 = 10.0**log_gamma_cd

  end function xkcidM96

  !-----------------------------------------------------------------------
  ! Martin+96, dissociative tunneling
  !-----------------------------------------------------------------------  
  function xkdtM96(xnH, T, logT, logT2, logT3)
    real(kind=DBL_KIND) :: xkdtM96
    real(kind=DBL_KIND),intent(IN) :: xnH, T, logT, logT2, logT3
    real(kind=DBL_KIND),dimension(21) :: a
    real(kind=DBL_KIND) :: log_g_h1, log_g_h2, log_g_l1, log_g_l2,&
         log_n_c1, n_c1, log_n_c2, n_c2, p, log_gamma_dt

    data a /-1.427664e2, 4.270741e1, -2.027365, -2.582097e-1,&
         2.136094e1,  2.753531e4, -2.146779e4,  6.034928e1, -2.743096e1,&
         2.676150,   -1.128215e1,  1.425455e4, -2.312520e4,  9.305564,  &
         -2.464009,   1.985955e-1, 7.430600e2, -1.174242,    7.502286e-1,&
         2.358848e-1, 2.937507/

  log_g_h1 = a(1) + a(2)*logT + a(3)*logT2 + a(4)*logT3 + a(5)*log10(1+a(6)/T)
  log_g_h2 = a(7)/T

  log_g_l1 = a(8) + a(9)*logT + a(10)*logT2 + a(11)*log10(1+a(12)/T)
  log_g_l2 = a(13)/T

  log_n_c1 = a(14) + a(15)*logT + a(16)*logT2 + a(17)/T
  n_c1 = 10.0**log_n_c1

  log_n_c2 = a(18) + log_n_c1
  n_c2 = 10.0**log_n_c2

  p = a(19) + a(20)*exp(-T/1850.0) + a(21)*exp(-T/440.0)

  log_gamma_dt =   log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+(xnH/n_c1)**p) &
                 + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+(xnH/n_c2)**p)
  
  xkdtM96 = 10.0**log_gamma_dt

  end function xkdtM96


  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$u
  !$$$$$$$$$$$$$$$$$$$$$$$             COOLING            $$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !-----------------------------------------------------------------------
  ! total cooling rate via all processes
  !-----------------------------------------------------------------------  
  subroutine tot_cool(xch,T_K,xLmbd_tot,xk_in)
    type(chem_vals) :: xch
    real(kind=DBL_KIND),intent(IN) :: T_K
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_tot !cooling rate (erg/cm^3/s)
    real(kind=DBL_KIND),dimension(0:NREACT-1),intent(IN),optional :: xk_in !高速化のため反応係数を外から受け取れるようにする
    real(kind=DBL_KIND) :: xkd, tau_cnt, esc_cont

    real(kind=DBL_KIND) :: xk(0:NREACT-1), chi_d
    real(kind=DBL_KIND) :: xLmbd_chem,xLmbd_line, dust_cool, ph_heat, ne, rOpi
    real(kind=DBL_KIND) :: inside_exp


    !外からxk_inが与えられていたらそれをxkに代入して終わり
    if (present(xk_in)) then
       xk(:) = xk_in(:)
    else
       !与えられていなければここで計算
       !get xk(:) from xnH, T_K, and radiation
       call react_coef(xch%nH,T_K,xch%rHpi,xch%rH2pd,xch%rHmpd, xch%Td, xch%fd, xk)
    end if

#ifdef METAL
    rOpi = xch%rOII ! OII のionization rate [s^{-1}]
#endif

    call chem_cool(xch%nH,T_K,xch%Td,xch%ychem,xk,xch%heat,xLmbd_chem)
    call line_cool(xch%nH,T_K,xch%Td,xch%ychem,xch%yco,xch%xlmbdj,xLmbd_line,xch%fd,xch%metal &
#ifdef METAL
        , rOpi &
#endif
      )

    ! dust gas coupling
#if DUST_GAS_COUPLING == YES
#ifndef METALCOOLOFF
    call dtemp_radtr(xch%nH, T_K, xch%Td, xch%EradIR, chi_d, dph_in=xch%rdph)

    !call plank_dop(Td, xkd)
    !tau_cnt  = xkd * rho * 0.5d0 * xlmbdj
    !esc_cont = min(1.d0, tau_cnt**(-2.d0))
    !dust_cool = 4.d0 * cgs_sigma* xkd * (Td**4.d0 * esc_cont- Tcmb**4.d0)  * rho - rdph ! [erg/cm^3/s]

    dust_cool = 5.83e-8*xch%nH*((1.d0 + 4.d0*yHe)*CONST_mH)*xch%nH*sqrt(T_K*1.d-3) &
        *(1.d0-0.8d0*exp(-75.d0/T_K))*(T_K-xch%Td)*xch%fd  ! [erg/cm^3/s] ダスト数密度の変化考慮すべき
#else
    dust_cool = 0.e0
#endif
#endif
    
#if PHOTOELECTRIC_HEATING == YES
    ne = xch%ychem(2) * xch%nH 
    call PhotoelectricHeating(ne, xch%nH, T_K, xch%rgfuv, ph_heat, xch%fd) 
#else
    ph_heat = 0.e0
#endif
#ifdef NO_RADIATION
    ph_heat = 0.e0
#endif

    xLmbd_tot = xLmbd_chem + xLmbd_line

#if DUST_GAS_COUPLING == YES
    xLmbd_tot = xLmbd_tot  + dust_cool
#endif

#if PHOTOELECTRIC_HEATING == YES
    inside_exp = (xch%ychem(3) - 0.01)*1.d5
    if(inside_exp < 1.d2) then
      ph_heat   = ph_heat * 1.d0 / (1.d0 + dexp(inside_exp))
    else
      ph_heat   = 0.d0
    endif
    xLmbd_tot = xLmbd_tot  - ph_heat  !電離領域内部で温度が急激にあがるのを防ぐ
#endif

    !KS DEBUG
    if (dbg_flg_prm == 1) then
       print '(A,(1P3E15.7))', "xLmbd_chem, xLmbd_line, xLmbd_tot, dust_cool, ph_heat: ",xLmbd_chem&
        , xLmbd_line, xLmbd_tot, dust_cool, ph_heat, xch%rgfuv !KS DEBUG
       print '(A,/)', "#-----------------------------------------------------------#"
    end if    


  end subroutine tot_cool

  !-----------------------------------------------------------------------
  ! total cooling rate via chemical reactions
  !-----------------------------------------------------------------------  
  subroutine chem_cool(xnH,T_K,Td,y,xk,heat,xLmbd_chem)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y(0:NCHEM-1),xk(0:NREACT-1), Td
    real(kind=DBL_KIND),intent(IN) :: heat !mean energy transferred to gas per one photoionizaiton event (erg)
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_chem !cooling rate (erg/cm^3/s)
    real(kind=DBL_KIND) :: xn_cr,crit,rtHm,rt3b,rtdis,rtHci,rtHra,rtHpi,xL_H2,xL_ci,xL_rec,xL_HM,xG_pi
    real(kind=DBL_KIND) :: heat_H2disso, rtH2d
    real(kind=DBL_KIND) :: log_T,log_T2,log_T3,log_T4,log_T5 
    
    real(kind=DBL_KIND),dimension(0:NCHEM-1) :: r_f ! KS DEBUG
#ifdef METAL
    real(kind=DBL_KIND) :: xL_H2d
#endif

    !log10(T_K)
    log_T  = log10(T_K)
    log_T2 = log_T*log_T
    log_T3 = log_T2*log_T
    log_T4 = log_T3*log_T
    log_T5 = log_T4*log_T

    !< H_2 formation/dissociation cooling >!
    !KS NOTE: energy gain at H2 formation from H-: 4.48 - 0.75 = 3.73 eV does not equal to 3.53eV
    xn_cr =  1.d6/dsqrt(T_K) &
         /(1.6d0*y(0)*dexp(-(4.d2/T_K)**2) + 1.4d0*y(1)*dexp(-1.2d4/(T_K+1.2d3)))
    crit = 1.d0/(1.d0 + xn_cr/xnH)
    rtHm = xk(2)*y(4)*y(0)     !< H- + H -> H2 +  e >!
    rt3b = (xk(6)*(y(0)**3) + xk(7)*(y(0)**2)*y(1))*xnH  !< 3 H  -> H2 +  H  >,  < 2H  + H2  -> 2H2 >!
    rtdis = xk(5)*y(1)*y(0) + xk(8)*(y(1)**2)     !< H2 + H  -> 3 H   >,  < 2 H2  -> 2 H + H2 >!
    xL_H2 = - (3.73*rtHm*crit + 4.48*(rt3b*crit-rtdis)) * xnH*xnH*evolt ! <-- 3.53 is corrected to 3.73
#ifdef METAL
#if DUST_GAS_COUPLING == YES
    rtH2d  = xk(21)*y(0)
    xL_H2d = -(0.2d0+4.2d0*crit)*rtH2d*xnH*xnH*evolt
#else
    xL_H2d = 0.d0
#endif
#endif

#ifdef METALCOOLOFF
    xL_H2d = 0.d0
#endif


    !collisional ionization cooling: < H + e -> H+ + 2 e > 
    rtHci = xk(0)*y(0)*y(2) 
    xL_ci = rtHci*13.6*evolt*xnH*xnH


    !radiative recombination cooling: < H+ + e -> H + ph. >
    ! Sugimura+17's fit (< 1% error for 1e2K < T < 1e7K) to Ferland et al. (1992)
#ifdef RECOM_ION_PHOTO
    xL_rec = 1d1 ** (-26.02 + 9.187d-1*log_T - 3.733d-1*log_T2&
         + 1.174d-1*log_T3 - 1.707d-2*log_T4 + 8.119d-4*log_T5)* y(3) * y(2) * xnH * xnH ! Case A
#else
    xL_rec = 1d1 ** (-25.87 + 4.958d-1*log_T - 1.052d-1*log_T2&
         + 4.264d-2*log_T3 - 9.165d-3*log_T4 + 5.491d-4*log_T5) * y(3) * y(2) * xnH * xnH ! Case B
#endif

    !radiative association cooling: < H + e -> H- + ph. > 
    rtHra = xk(9)*y(0)*y(2)
    xL_HM  = boltz*T_K*rtHra*xnH*xnH !<- kb * T_K * xhmrat * xnH の方がよいのでは？（電子の運動エネルギー分の冷却, KS TODO）
    !xL_HM  = rtHra*0.755*evolt*xnH*xnH !<- kb * T_K * xhmrat * xnH の方がよいのでは？（電子の運動エネルギー分の冷却, KS TODO）


    !photo-ionization heating: < H + ph. -> H+ + e >
    rtHpi  = xk(10)*y(0)  
    xG_pi = rtHpi*heat*xnH

    
    ! heating of photodissociation
#ifdef H2PD_HEATING
    rtH2d        = xk(11)*y(1) ! [s^{-1}]
    heat_H2disso = 9.d0*rtH2d*(2.2d0*evolt*crit)*xnH ! [ erg cm^-3 s^-1] (Hollenbach & McKee 1979 ))
#else
    heat_H2disso = 0.d0
#endif



    ! ================ total rate per vol. ==============
    xLmbd_chem =  xL_H2 + xL_ci + xL_rec + xL_HM - xG_pi - heat_H2disso
#ifdef METAL
    xLmbd_chem =  xLmbd_chem + xL_H2d
#endif
    ! ===================================================
    ! print '(A,(1P4E15.7))', "xL_H2: ", rtHm,crit,rt3b,rtdis !KS DEBUG
    
    !KS DEBUG
    if (dbg_flg_prm == 1) then
       print '(/,A)', "#-------------   reaction and cooling rates   --------------#"
       ! call react_rat(xk,xnH,y,r_f)
       call react_rat_highspeed(xk,xnH,y,r_f)       
       !ydot(H): radiative recomb., collisional ion., photo-ion., H2 form/dis
       print '(A,(1P6E16.7E3))', "ydot(H) rr, ci, pi, -2*ydot(H2), sum, full: ",&
            xk(1)*y(3)*y(2)*xnH, xk(0)*y(0)*y(2)*xnH, rtHpi, -2.*r_f(1),  &
            xk(1)*y(3)*y(2)*xnH - xk(0)*y(0)*y(2)*xnH - rtHpi - 2.*r_f(1), &
            r_f(0)
       !ydot(H2): associative detacth..,  H2p-channel, 3-body formatoin, collisional dis., photo-dis.
       print '(A,(1P7E16.7E3))', "ydot(H2) ad, 3b, h2p, cd, pd, sum, full: ",& 
            rtHm*xnH, rt3b*xnH, xk(18)*y(5)*y(0)*xnH, rtdis*xnH, xk(11)*y(1), &
            rtHm*xnH + rt3b*xnH + xk(18)*y(5)*y(0)*xnH - rtdis*xnH - xk(11)*y(1), r_f(1)
       !ydot(Hm): radiative association, associative det., el-collsion det., photo-det.
       print '(A,(1P6E16.7E3))', "ydot(Hm) ra, ad, pd, ecd, sum, full: ",&
            rtHra*xnH, rtHm*xnH, xk(16)*y(4), xk(13)*y(4)*y(2)*xnH, &
            rtHra*xnH - rtHm*xnH - xk(16)*y(4) - xk(13)*y(4)*y(2)*xnH, r_f(4)
       !cooling
       print '(A,(1P5E16.7E3))', "xL_H2,xL_ci,xL_rec,xL_HM,xG_pi: ",&
            xL_H2,xL_ci,xL_rec,xL_HM,xG_pi
    end if


    end subroutine chem_cool

  !-----------------------------------------------------------------------
  ! total cooling rate via line emission (w/o chemical transition)
  ! dimension of xLmdb_line is erg/cm^3/s
  !-----------------------------------------------------------------------  
  subroutine line_cool(xnH,T_K,Td,y,yco,xlmbdj,xLmbd_line,fd,metal &
#ifdef METAL
        , rOpi &
#endif
    )
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y(0:NCHEM-1),Td,yco,fd,metal
#ifdef METAL
    real(kind=DBL_KIND),intent(IN) :: rOpi 
#endif
    real(kind=DBL_KIND),intent(IN) :: xlmbdj ! jeans length
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_line !cooling rate (erg/cm^3/s)

    real(kind=DBL_KIND) :: xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, y_H2, y_H, y_e, y_Hp, T_5, tau_cnt

    real(kind=DBL_KIND) :: gff1, xLmbd_ff, xNc_H2, cooling_CO, heating_CO, xLmbd_CO
#if METAL_LINE == YES
    real(kind=DBL_KIND) :: y_CII, nCII, xNc_CII, ne, nH, cool_rate_CII, xLmbd_CII
    real(kind=DBL_KIND) :: y_OI, nOI, xNc_OI, cool_rate_OI, xLmbd_OI, y_OII, nOII,&
      xNc_OII, cool_rate_OII, xLmbd_OII, y_OIII, nOIII, xNc_OIII, cool_rate_OIII, xLmbd_OIII
    real(kind=DBL_KIND),dimension(0:2) :: dvdr, b_cont
    real(kind=DBL_KIND) :: rho, t_ff, radius, v_bulk, dv_d, kappa
    real(kind=DBL_KIND) :: nCO, xNc_CO
#endif
    integer :: ii
    real(kind=DBL_KIND), parameter :: min_value = 1.d-50
    real(kind=DBL_KIND) :: frac_C, frac_O
    
    real(kind=DBL_KIND) :: pi
    data pi/3.14159265358979d0/


    ! set frac_C, frac_O
    frac_C = metal*frac_C_solar
    frac_O = metal*frac_O_solar

    y_H = y(0)  !same as y_Hn
    y_H2 = y(1) !same as y_H2
    y_e = y(2)  !same as y_el
    y_Hp = y(3) !same as y_Hp

    ! -------- column density each particle ----------
    xNc_H2 = y_H2*xnH*0.5*xlmbdj 

#if METAL_LINE == YES

    ne      = y_e*xnH
    nH      = y_H*xnH

    ! ----------------------------------------
    rho    = ((1.d0 + 4.d0*yHe)*CONST_mH)*xnH      ! [g cm^{-3}]
    t_ff   = sqrt(3.0*pi/(32.0*CONST_G*rho))    ! free fall time
    radius = 0.5d0*xlmbdj 
    v_bulk = radius/3.0/t_ff
    ! ----------------------------------------

    ! ---------------------------
    y_CII  = ( frac_C - yco)
    y_CII  = MAX(MIN(frac_C,y_CII),min_value)
    nCII   = y_CII*xnH
    xNc_CII = y_CII*xnH*radius
    ! ---------------------------

    ! ---------------------------
    nCO    = yco * xnH
    xNc_CO = yco*xnH*radius 
    ! ---------------------------

#if Oxygen_COOL_CAL == YES

    ! --------------------------
    y_OI   = frac_O*(1.d0 - y_Hp - yco)  ! oxgen ionization rate is equal to HI
    y_OI  = MAX(MIN(frac_O,y_OI),min_value)
    nOI    = y_OI*xnH
    xNc_OI = y_OI*xnH*radius
    ! --------------------------

    call OII_OII_ratio(ne, y_Hp, T_K, rOpi, y_OII, y_OIII, metal)

    ! --------------------------
    nOII   = y_OII*xnH
    xNc_OII= y_OII*xnH*radius
    ! --------------------------

    ! --------------------------
    nOIII   = y_OIII*xnH
    xNc_OIII= y_OIII*xnH*radius
    ! --------------------------
#else
    ! --------------------------
    y_OI   = frac_O  ! oxgen ionization rate is equal to HI
    nOI    = y_OI*xnH
    xNc_OI = y_OI*xnH*radius
    ! --------------------------
#endif

    ! --------- optical depth for dust continium --------
    !print *, 'before cool_tot', Td
    !call plank_dop(Td, kappa)
    call find_dop(Td, kappa)  ! opacity at solar metallicity
    tau_cnt = (kappa*fd)*rho*radius                 ! [noD] 

    do ii = 0, 2
        b_cont(ii) = exp(-tau_cnt) 
    enddo
#else
    tau_cnt = 0.d0 !assume opt-thin for continuum
#endif
    ! ------------------ H2 line cooling -------------------------
    call H2cool(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)


    ! ------------------  HI Lya cooling ------------------------
    T_5 = 1.d-5*T_K
    xLmbd_Lya = y_e*y_H*7.50d-19/(1.d0 + dsqrt(T_5)) &
         *dexp(-1.18348d5/T_K) *xnH*xnH   

    
    ! -----------------  Compton cooling --------------------- 
    xLmbd_cpt = 5.65d-36*((1.d0 + zred)**4)*(T_K - MP_Tcmb)*xnH*y_e


    ! -------------------- HeII cooling ------------------------ 
    xLmbd_Hep = 5.54d-17/(1.d0 + sqrt(T_5))/(T_K**0.397) &
         *exp(-4.73638d+5/T_K)*yHe*(y_e*xnH)**2


    ! -------------------- H free-free (GJ07)  ------------------------    (KS TODO)
    if (T_K<3.2d5) then
       gff1 = 0.79464 + 0.1243 * dlog10(T_K)
    else
       gff1 = 2.13164 - 0.1240 * dlog10(T_K)
    endif
    xLmbd_ff = 1.426e-27 * sqrt(T_K) * gff1* y_Hp * y_e * xnH * xnH
    ! -----------------------------------------------------------


#if METAL_LINE == YES
    ! ------------------- Metal line cooling  -----------------------

    ! --------------------
    !   CII line cooling
    ! --------------------
    dv_d     = sqrt(2.0*CONST_kB*T_K/(CONST_AC)/CONST_amu) ! [cm g^{-1}]
    do ii = 0, 2
        dvdr(ii)   = max(dv_d/radius, v_bulk/radius) 
    end do
    call carbonII_cooling(nCII, ne, nH, T_K, MP_Tcmb, dvdr, b_cont, cool_rate_CII) 
    xLmbd_CII = nCII * cool_rate_CII

    ! ---------------------
    !   CO line cooling 
    ! ---------------------
    call COcool(nH,T_K ,y_H2,xNc_CO,tau_cnt,cooling_CO)
    call COcool(nH,MP_Tcmb,y_H2,xNc_CO,tau_cnt,heating_CO)
    xLmbd_CO  = nCO * (cooling_co - heating_co) ! [erg/s/cm^3]

    ! --------------------
    !   OI line cooling
    ! --------------------
    dv_d     = sqrt(2.0*CONST_kB*T_K/(CONST_AO)/CONST_amu) ! [cm g^{-1}]
    do ii = 0, 2
        dvdr(ii)   = max(dv_d/radius, v_bulk/radius) 
    end do

    call oxygenI_cooling(nOI, ne, nH, T_K, MP_Tcmb, dvdr, b_cont, cool_rate_OI)
    xLmbd_OI = nOI * cool_rate_OI                          ! [erg/s/cm^3] 


#if Oxygen_COOL_CAL == YES
    call oxygenII_cooling(nOII, ne,   T_K, MP_Tcmb, dvdr, b_cont, cool_rate_OII)
    xLmbd_OII   = nOII* cool_rate_OII                          ! [erg/s/cm^3] 


    call oxygenIII_cooling(nOIII, ne, nH, T_K, MP_Tcmb, dvdr, b_cont, cool_rate_OIII)
    xLmbd_OIII  = nOIII * cool_rate_OIII

#endif

    !test
    !print *, 'line', xLmbd_H2, xLmbd_Lya, xLmbd_cpt,xLmbd_Hep, xLmbd_ff, xLmbd_CII, xLmbd_OI, xLmbd_OII &
    !     ,xLmbd_OIII 

#endif


    ! ================ total rate per vol. ==============
    ! xLmbd_line =  xLmbd_H2 + xLmbd_Lya + xLmbd_cpt + xLmbd_Hep
    xLmbd_line =  xLmbd_H2 + xLmbd_Lya + xLmbd_cpt + xLmbd_Hep  + xLmbd_ff  

#if METAL_LINE == YES
#ifndef METALCOOLOFF
    xLmbd_line =  xLmbd_line + xLmbd_CII + xLmbd_OI + xLmbd_CO  
#endif
#if Oxygen_COOL_CAL == YES
    xLmbd_line =  xLmbd_line + xLmbd_OII + xLmbd_OIII
#endif
#endif

    ! ===================================================
    ! print '(A,(1P4E15.7))', "xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya: ",&
    !      xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya ! KS DEBUG

    !KS DEBUG
    if (dbg_flg_prm == 1) then
! print '(A,(1P4E15.7))', "xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya: ",&
!      xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya ! KS DEBUG
print '(A,(1P5E15.7))', "xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, xLmbd_ff: ",&
     xLmbd_H2, xLmbd_Hep, xLmbd_cpt, xLmbd_Lya, xLmbd_ff ! KS DEBUG
    end if        

  end subroutine line_cool

  !-----------------------------------------------------------------------
  ! Cooling rate via H2 molecular line emission
  !  using HM formula to interpolate low-dens and LTE(+line esc) formulae
  !-----------------------------------------------------------------------  
  subroutine H2cool(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_H2

    real(kind=DBL_KIND) :: xLdH2_n0,xLdH2_LTE,y_He,T_rad, xLdH2_dbg,xLdH2_dbg2

    y_He = yHe !assume all helium is in the neutral state
    T_rad = 0d0 !assume no absoprition of stellar radiation by H2 lines

    ! call H2cool_LowD_G15(xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He,xLdH2_n0)
    ! call H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau_cnt,T_rad,xLdH2_LTE)

    !KS DEBUG
    ! call H2cool_LTEfit(xnH,T_K,y_H2,0d0,tau_cnt,T_rad,xLdH2_dbg)
    ! call H2cool_HM(xnH,T_K,y_H,y_H2,xLdH2_dbg2)

    !lineがthickになった後に、opt-thin low-density cooling rateが現れないように場合分け
    if(xnH.lt.1.d+8) then
       call H2cool_LowD_G15(xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He,xLdH2_n0)
       call H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau_cnt,T_rad,xLdH2_LTE)
       xLmbd_H2 = xLdH2_LTE/(1.d0+xLdH2_LTE/xLdH2_n0)
    else
       call H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau_cnt,T_rad,xLdH2_LTE)
       xLmbd_H2 = xLdH2_LTE 
    end if

    !KS DEBUG
    ! print *, xLdH2_n0/xLdH2_dbg, xLdH2_LTE/xLdH2_dbg, xLmbd_H2/xLdH2_dbg, xLdH2_dbg2/xLdH2_dbg

  end subroutine H2cool

  !-----------------------------------------------------------------------
  ! H2 cooling rate using the Hollenbach-McKee formula
  !-----------------------------------------------------------------------  
  subroutine H2cool_HM(xnH,T_K,y_H,y_H2,xLmbd_H2)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y_H,y_H2
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_H2

    real(kind=DBL_KIND), parameter :: E_r_20 = 512.0d0*boltz
    real(kind=DBL_KIND), parameter :: E_r_31 = 853.3d0*boltz
    real(kind=DBL_KIND), parameter :: E_v_10 = 5860.0d0*boltz
    real(kind=DBL_KIND), parameter :: E_v_20 = 11720.0d0*boltz

    real(kind=DBL_KIND) :: lgT, sqT, T3,&
         kd_a, kd_m, L_r_H_LTE, L_r_H2_LTE, L_v_H_LTE, L_v_H2_LTE,&
         lg_L_r_H_n0, L_r_H_n0, L_r_H2_n0, L_v_H_n0, L_v_H2_n0,&
         kd_cr_r_H, kd_cr_v_H, kd_cr_r_H2, kd_cr_v_H2, L_vr_H,&
         L_vr_H2, Lda_rv,&
         Tiv, xkTiv, T3iv, kdaiv, kdmiv,&
         L_r_LTE, L_v_LTE, gamma_r_H2_2, gamma_r_H2_3, gamma_v_10_H,&
         gamma_v_20_H, gamma_v_10_H2

    if(y_H2.gt.0.d0 .and. T_K.gt.3.d+1) then

      kd_a = xnH*y_H
      kd_m = xnH*y_H2

      kdaiv = 1.d0/kd_a
      kdmiv = 1.d0/kd_m

      Tiv    = 1.d0/T_K
      xkTiv  = Tiv/boltz

      lgT    = log10(T_K)
      sqT    = sqrt(T_K)
      T3     = 1.d-3*T_K
      T3iv   = 1.d+3*Tiv

      L_r_LTE = ( 9.5d-22*(T3**3.76)/(1.d0 + 0.12d0*(T3**2.1)))*exp(-(0.13*T3iv)**3)&
           + 3.0d-24*exp(-0.51d0*T3iv)       
      L_v_LTE = ( 6.7d-19*exp(-5.86d0*T3iv) + 1.6d-18*exp(-1.17d1*T3iv))

      L_r_H_LTE  = L_r_LTE*kdaiv
      L_r_H2_LTE = L_r_LTE*kdmiv           
      L_v_H_LTE  = L_v_LTE*kdaiv
      L_v_H2_LTE = L_v_LTE*kdmiv

      lg_L_r_H_n0 = - 1.03d2 + 9.759d+1*lgT - 4.805d+1*(lgT*lgT)&  !<= GP98
           + 1.08d+1*(lgT*lgT*lgT) - 0.9032d0*(lgT**4)
      L_r_H_n0  = 10.0**lg_L_r_H_n0

      gamma_r_H2_2 = ( 3.3d-12 + 6.6d-12*T3 )*7.0071412d-1
      gamma_r_H2_3 = ( 3.3d-12 + 6.6d-12*T3 )*1.0041880d+0

      L_r_H2_n0 =   0.25d0*( 5.0d0*gamma_r_H2_2 *exp(-E_r_20*xkTiv)*E_r_20 )&
           + 0.75d0*( 2.333d0*gamma_r_H2_3 *exp(-E_r_31*xkTiv)*E_r_31 )

      gamma_v_10_H  = 1.0d-12*sqT*exp(-1.0d+3*Tiv)
      gamma_v_20_H  = 1.6d-12*sqT*exp(-(4.0d+2*Tiv)**2)
      gamma_v_10_H2 = 1.4d-12*sqT*exp(-1.2d+4/(T_K + 1.2d+3))

      L_v_H_n0  =  gamma_v_10_H*exp(-E_v_10*xkTiv)*E_v_10 + gamma_v_20_H*exp(-E_v_20*xkTiv)*E_v_20         
      L_v_H2_n0 =  gamma_v_10_H2*exp(-E_v_10*xkTiv)*E_v_10
 
      kd_cr_r_H  = kd_a*L_r_H_LTE/L_r_H_n0
      kd_cr_v_H  = kd_a*L_v_H_LTE/L_v_H_n0
      kd_cr_r_H2 = kd_m*L_r_H2_LTE/L_r_H2_n0
      kd_cr_v_H2 = kd_m*L_v_H2_LTE/L_v_H2_n0

      L_vr_H  =   L_r_H_LTE/( 1.d0 + kd_cr_r_H*kdaiv ) + L_v_H_LTE/( 1.d0 + kd_cr_v_H*kdaiv )
      L_vr_H2 =   L_r_H2_LTE/( 1.d0 + kd_cr_r_H2*kdmiv ) + L_v_H2_LTE/( 1.d0 + kd_cr_v_H2*kdmiv )

      Lda_rv = kd_m*( y_H*L_vr_H + y_H2*L_vr_H2 )

      xLmbd_H2 = 0.55d0*Lda_rv*xnH               ! per unit vol.

      else
        xLmbd_H2 = 0.d0
      end if

      return 
  end subroutine H2cool_HM

  !-----------------------------------------------------------------------
  ! background radiation
  !-----------------------------------------------------------------------  
  function  Q_bg(T_nu)
    real(kind=DBL_KIND) :: Q_bg
    real(kind=DBL_KIND),intent(IN) :: T_nu

    real(kind=DBL_KIND) :: T_rad, x, Q_bg_CMB 

    T_rad = MP_Tcmb
    x     = T_nu/T_rad
    
    if(x.gt.1.d2) then
       Q_bg_CMB = 0.d0
    else
       Q_bg_CMB = 1.d0/(dexp(x)-1.d0)
    endif
    
    Q_bg = Q_bg_CMB
    
    return
  end function Q_bg

  !-----------------------------------------------------------------------
  ! escape probability
  !-----------------------------------------------------------------------  
  function  beta_esc(tau_L,tau_C)
    real(kind=DBL_KIND) :: beta_esc
    real(kind=DBL_KIND) :: tau_L, tau_C

    if(tau_L.lt.0.d0) then
       beta_esc=1.d0
    elseif(tau_L.lt.1.d-5) then 
       beta_esc = dexp(-tau_C)
    else
       beta_esc = dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
    endif
    return
  end function beta_esc


  !-----------------------------------------------------------------------
  ! gamma of H2 molecules
  !-----------------------------------------------------------------------  

  function c_H2(T_K)
    real(kind=DBL_KIND),intent(IN) :: T_K
    real(kind=DBL_KIND) :: c_H2, dT
    integer :: ii

    real(kind=DBL_KIND) :: logTk
    real(kind=DBL_KIND),dimension(50) :: ca, xlTa

    data xlTa/ 1.00000d0,1.08163d0,1.16327d0,1.24490d0,1.32653d0,1.40816d0 &
    ,1.48980d0,1.57143d0,1.65306d0,1.73469d0,1.81633d0,1.89796d0,1.97959d0 &
    ,2.06122d0,2.14286d0,2.22449d0,2.30612d0,2.38776d0,2.46939d0,2.55102d0 &
    ,2.63265d0,2.71429d0,2.79592d0,2.87755d0,2.95918d0,3.04082d0,3.12245d0 &
    ,3.20408d0,3.28571d0,3.36735d0,3.44898d0,3.53061d0,3.61224d0,3.69388d0 &
    ,3.77551d0,3.85714d0,3.93878d0,4.02041d0,4.10204d0,4.18367d0,4.26531d0 &
    ,4.34694d0,4.42857d0,4.51020d0,4.59184d0,4.67347d0,4.75510d0,4.83673d0 &
    ,4.91837d0,5.00000d0 /

    data ca/ 1.50000d0,1.50000d0,1.50000d0,1.50000d0,1.50000d0,1.50000d0   &
    ,1.50000d0,1.50000d0,1.50482d0,1.51494d0,1.53443d0,1.59157d0,1.68857d0 &
    ,1.81833d0,1.97800d0,2.14694d0,2.27310d0,2.37931d0,2.44969d0,2.48382d0 &
    ,2.49866d0,2.50000d0,2.50559d0,2.52365d0,2.55735d0,2.62265d0,2.71562d0 &
    ,2.83053d0,2.95857d0,3.07868d0,3.17516d0,3.26061d0,3.32957d0,3.38085d0 &
    ,3.41796d0,3.44314d0,3.46176d0,3.47367d0,3.48120d0,3.48937d0,3.49000d0 &
    ,3.49518d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0,3.50000d0 &
    ,3.50000d0,3.50000d0 /

    logTk = dlog10(T_K)
    ii = int((logTk-1.d0)*49.d0/4.d0) + 1

    if(ii < 1) then
      c_H2 = ca(1)
    else if(ii>49) then
      c_H2 = ca(50)
    else
      c_H2=(ca(ii+1)-ca(ii))/(xlTa(ii+1)-xlTa(ii))*(logTk-xlTa(ii))+ca(ii)
    endif

  end function c_H2

  function c_H2_2(T_K)
    real(kind=DBL_KIND),intent(IN) :: T_K    
    real(kind=DBL_KIND) :: c_H2_2

    real(kind=DBL_KIND),dimension(30) :: xlTa, ca
    real(kind=DBL_KIND) :: xlT
    
      data xlta / 0.150d1, 0.160d1, 0.171d1, 0.181d1, 0.191d1,&
                  0.202d1, 0.212d1, 0.222d1, 0.233d1, 0.243d1,&
                  0.253d1, 0.264d1, 0.274d1, 0.284d1, 0.295d1,&
                  0.305d1, 0.316d1, 0.326d1, 0.336d1, 0.347d1,&
                  0.357d1, 0.367d1, 0.378d1, 0.388d1, 0.398d1,&
                  0.409d1, 0.419d1, 0.429d1, 0.440d1, 0.450d1 /

      data ca / 0.150d1, 0.150d1, 0.151d1, 0.153d1, 0.160d1,&
                0.174d1, 0.193d1, 0.214d1, 0.231d1, 0.243d1,&
                0.248d1, 0.250d1, 0.250d1, 0.251d1, 0.255d1,&
                0.263d1, 0.276d1, 0.292d1, 0.307d1, 0.320d1,&
                0.330d1, 0.337d1, 0.342d1, 0.345d1, 0.347d1,&
                0.348d1, 0.349d1, 0.349d1, 0.350d1, 0.350d1 /

      xlT = dlog10(T_K)
      call linear(xlTa,ca,30,xlT,c_H2_2)
      return
    contains
    !-----------------------------------------------------------------------
    ! liner interpolation
    !-----------------------------------------------------------------------  
    subroutine linear(xa,ya,m,x,y)
      integer,intent(IN) :: m
      real(kind=DBL_KIND),dimension(m),intent(IN) :: xa,ya
      real(kind=DBL_KIND),intent(IN) :: x
      real(kind=DBL_KIND),intent(OUT) :: y

      integer :: ms, i
      real(kind=DBL_KIND) :: t, y1, y2
      
      ms = 0
      do i=1,m
         if(x-xa(i).le.0.d0 .or. i.eq.m) then
            ms = i
            exit
         endif
      end do

      if(ms.eq.0) then
        print *, 'zero ms in linear (primordial)'
        stop
      end if

      if(ms.eq.1) ms=2
      y1 = ya(ms-1)
      y2 = ya(ms)
      t  = (x - xa(ms-1))/(xa(ms) - xa(ms-1))
      y  = (1.d0-t)*y1 + t*y2

      return

    end subroutine linear
  end function c_H2_2

  !-----------------------------------------------------------------------
  ! H2 cooling under the assumption of statistical equilibrium
  !
  !     xLmbd_H2 : H2 cooling rate per unit volume
  !-----------------------------------------------------------------------
  subroutine H2cool_Omukai(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLmbd_H2)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_H2

    real(kind=DBL_KIND) :: A(0:2,0:2,0:22,0:22),ET(0:2,0:22),p(0:22),f(0:2,0:22),f_v(0:2)
    real(kind=DBL_KIND) :: xk_B,h_P,pi,xm_p,c_light
    data xk_B/1.380662d-16/,h_P/6.626176d-27/,pi/3.14159265358979d0/,xm_p/1.67d-24/,c_light/2.99792458d10/

    real(kind=DBL_KIND) :: gamma_h_j2down_omukai(0:18)
    real(kind=DBL_KIND) g_0, g_1, g_2, A_10, A_20, A_21,&
         gamma_10H, gamma_20H, gamma_21H, gamma_10H2, gamma_20H2, gamma_21H2,&
         gamma_10e, gamma_20e, gamma_21e, gamma_10Hp, gamma_20Hp, gamma_21Hp,&
         xn_a, xn_m,&
         C_10, C_20, C_21, DT_10, DT_20, DT_21, C_01, C_02, C_12,&
         Q_10, Q_20, Q_21,&
         R_10, R_20, R_21, R_01, R_02, R_12,&
         f_0, f_1, f_2, f_para, f_ortho, z_para, z_ortho, &
         gamma_H, gamma_H2, gamma_e, gamma_Hp, xj, xj2, xj3,  C_ul, DT, Q_ul, g_u, g_l, r,&
         v_th, DE, xnu_Hz, xnu_Hz3, xNc_u, xNc_l, tau_ul, esc, x, xmeV, xn_e, xn_Hp

    integer :: i,ii,iv,ivf,ivi,j,jf,ji


    data (((A(i,ii,j,j+2),ii=0,i-1),i=1,2),j=0,20)  &
         /8.54d-7 ,3.47d-7 ,1.29d-6   &
         ,4.23d-7 ,1.61d-7 ,6.40d-7   &
         ,2.90d-7 ,1.03d-7 ,4.41d-7   &
         ,2.09d-7 ,6.98d-8 ,3.18d-7   &
         ,1.50d-7 ,4.72d-8 ,2.28d-7   &
         ,1.06d-7 ,3.15d-8 ,1.62d-7   &
         ,7.38d-8 ,2.05d-8 ,1.12d-7   &
         ,4.98d-8 ,1.31d-8 ,7.50d-8   &
         ,3.27d-8 ,8.07d-9 ,4.88d-8   &
         ,2.07d-8 ,4.83d-9 ,3.06d-8   &
         ,1.27d-8 ,2.79d-9 ,1.85d-8   &
         ,7.44d-9 ,1.55d-9 ,1.07d-8   &
         ,4.16d-9 ,8.27d-10,5.84d-9   &
         ,2.19d-9 ,4.21d-10,3.00d-9   &
         ,1.08d-9 ,2.04d-10,1.43d-9   &
         ,4.87d-10,9.45d-11,6.17d-10  &
         ,1.96d-10,4.23d-11,2.33d-10  &
         ,6.71d-11,1.91d-11,7.32d-11  &
         ,1.82d-11,9.63d-12,1.73d-11  &
         ,3.38d-12,6.26d-12,2.46d-12  &
         ,2.96d-13,5.78d-12,1.16d-13/

    data (((A(i,ii,j,j),ii=0,i-1),i=1,2),j=1,20)  &
         /4.29d-7 ,1.94d-7 ,6.37d-7   &
         ,3.03d-7 ,1.38d-7 ,4.50d-7   &
         ,2.78d-7 ,1.29d-7 ,4.12d-7   &
         ,2.65d-7 ,1.25d-7 ,3.91d-7   &
         ,2.55d-7 ,1.23d-7 ,3.74d-7   &
         ,2.45d-7 ,1.21d-7 ,3.58d-7   &
         ,2.34d-7 ,1.20d-7 ,3.40d-7   &
         ,2.23d-7 ,1.18d-7 ,3.22d-7   &
         ,2.12d-7 ,1.17d-7 ,3.03d-7   &
         ,1.99d-7 ,1.15d-7 ,2.84d-7   &
         ,1.87d-7 ,1.13d-7 ,2.63d-7   &
         ,1.74d-7 ,1.11d-7 ,2.42d-7   &
         ,1.61d-7 ,1.08d-7 ,2.21d-7   &
         ,1.47d-7 ,1.05d-7 ,2.01d-7   &
         ,1.34d-7 ,1.02d-7 ,1.80d-7   &
         ,1.21d-7 ,9.88d-8 ,1.60d-7   &
         ,1.08d-7 ,9.50d-8 ,1.41d-7   &
         ,9.61d-8 ,9.07d-8 ,1.22d-7   &
         ,8.43d-8 ,8.62d-8 ,1.05d-7   &
         ,7.32d-8 ,8.14d-8 ,8.85d-8/

    data (((A(i,ii,j,j-2),ii=0,i),i=0,2),j=2,20) &
         /2.94d-11,2.53d-7 ,2.79d-11  &
         ,1.27d-7 ,3.68d-7 ,2.56d-11  &
         ,4.76d-10,3.47d-7,4.50d-10   &
         ,1.90d-7 ,4.98d-7 ,4.12d-10  &
         ,2.76d-9 ,3.98d-7 ,2.59d-9   &
         ,2.38d-7 ,5.60d-7 ,2.37d-9   &
         ,9.84d-9 ,4.21d-7 ,9.21d-9   &
         ,2.77d-7 ,5.77d-7 ,8.37d-9   &
         ,2.64d-8 ,4.19d-7 ,2.46d-8   &
         ,3.07d-7 ,5.57d-7 ,2.22d-8   &
         ,5.88d-8 ,3.96d-7 ,5.44d-8   &
         ,3.28d-7 ,5.05d-7 ,4.88d-8   &
         ,1.14d-7 ,3.54d-7 ,1.05d-7   &
         ,3.39d-7 ,4.30d-7 ,9.33d-8   &
         ,2.00d-7 ,2.98d-7 ,1.82d-7   &
         ,3.40d-7 ,3.38d-7 ,1.61d-7   &
         ,3.24d-7 ,2.34d-7 ,2.92d-7   &
         ,3.30d-7 ,2.41d-7 ,2.55d-7   &
         ,4.90d-7 ,1.68d-7 ,4.38d-7   &
         ,3.12d-7 ,1.49d-7 ,3.79d-7   &
         ,7.03d-7 ,1.05d-7 ,6.21d-7   &
         ,2.85d-7 ,7.20d-8 ,5.32d-7   &
         ,9.64d-7 ,5.30d-8 ,8.42d-7   &
         ,2.53d-7 ,1.96d-8 ,7.13d-7   &
         ,1.27d-6 ,1.65d-8 ,1.10d-6   &
         ,2.16d-7 ,1.49d-11,9.18d-7   &
         ,1.62d-6 ,4.26d-10,1.38d-6   &
         ,1.78d-7 ,1.91d-8 ,1.14d-6   &
         ,2.00d-6 ,8.38d-9 ,1.69d-6   &
         ,1.39d-7 ,8.07d-8 ,1.37d-6   &
         ,2.41d-6 ,4.27d-8 ,2.00d-6   &
         ,1.01d-7 ,1.86d-7 ,1.61d-6   &
         ,2.83d-6 ,1.04d-7 ,2.32d-6   &
         ,6.80d-8 ,3.35d-7 ,1.84d-6   &
         ,3.26d-6 ,1.93d-7 ,2.64d-6   &
         ,3.98d-8 ,5.24d-7 ,2.05d-6   &
         ,3.68d-6 ,3.08d-7 ,2.93d-6   &
         ,1.84d-8 ,7.49d-7 ,2.23d-6/ 

    g_0=1.d0
    g_1=1.d0
    g_2=1.d0

    A_10=8.3d-7
    A_20=4.1d-7
    A_21=1.1d-6

    !     Hollenbach & McKee (1979) corrected with (1989)
    gamma_10H=1.0d-12*dsqrt(T_K)*dexp(-1.d3/T_K)
    gamma_20H=1.6d-12*dsqrt(T_K)*dexp(-(4.d2/T_K)**2)
    gamma_21H=4.5d-12*dsqrt(T_K)*dexp(-(5.d2/T_K)**2)

    gamma_10H2=1.4d-12*dsqrt(T_K)*dexp(-1.81d4/(T_K+1.2d3))
    gamma_20H2=0.d0
    gamma_21H2=gamma_10H2

    do iv=0,2
       do j=0,22
          ET(iv,j)=E_H2_BFM(iv,J)
       end do
    end do

    DT_10=ET(1,0)-ET(0,0)
    DT_20=ET(2,0)-ET(0,0)
    DT_21=ET(2,0)-ET(1,0)

    !     Draine, Reberge, Dalgarno (1983)
    gamma_10e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_10/T_K)
    gamma_20e=2.5d-12*(T_K**0.5d0)/(1.d0+0.5d0*DT_20/T_K)
    gamma_21e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_21/T_K)

    !     Krstic (2002), my fit
    gamma_10Hp=1.4d-4*T_K**(-1.344d0) &
         *(1.d0+4.005d-9*T_K**2.066d0) &
         *dexp((DT_10-9589d0)/T_K)
    gamma_20Hp=4.585d-5*T_K**(-1.291d0) &
         *(1.d0+1.378d-8*T_K**1.903d0) &
         *dexp((DT_20-1.933d4)/T_K)
    gamma_21Hp=5.593d-3*T_K**(-1.770d0) &
         *(1.d0+3.505d-10*T_K**2.406d0) &
         *dexp((DT_21-1.460d4)/T_K)

    xn_a=y_H*xnH
    xn_m=y_H2*xnH
    xn_e=y_e*xnH
    xn_Hp=y_Hp*xnH

    C_10=gamma_10H*xn_a+gamma_10H2*xn_m &
         +gamma_10e*xn_e+gamma_10Hp*xn_Hp
    C_20=gamma_20H*xn_a+gamma_20H2*xn_m &
         +gamma_20e*xn_e+gamma_20Hp*xn_Hp
    C_21=gamma_21H*xn_a+gamma_21H2*xn_m &
         +gamma_21e*xn_e+gamma_21Hp*xn_Hp

    C_01=C_10*dexp(-DT_10/T_K)
    C_02=C_20*dexp(-DT_20/T_K)
    C_12=C_21*dexp(-DT_21/T_K)

    Q_10=Q_bg(DT_10)
    Q_20=Q_bg(DT_20)
    Q_21=Q_bg(DT_21)

    R_10=A_10*(1.d0+Q_10)+C_10
    R_20=A_20*(1.d0+Q_20)+C_20
    R_21=A_21*(1.d0+Q_21)+C_21
    R_01=(g_1/g_0)*A_10*Q_10+C_01
    R_02=(g_2/g_0)*A_20*Q_20+C_02
    R_12=(g_2/g_1)*A_21*Q_21+C_12

    f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) ) &
         /( (R_01+R_02+R_20)*(R_10+R_12+R_21)  &
         -(R_01-R_21)*(R_10-R_20) )
    f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
    f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

    f_v(0)=f_0
    f_v(1)=f_1
    f_v(2)=f_2


    f_para=0.25d0
    f_ortho=0.75d0
    do iv=0,2
       z_para=0.d0
       z_ortho=0.d0
       p(0)=1.d0
       p(1)=1.d0
       do j=0,18
          if(mod(j,2).eq.0) then
             z_para=z_para+p(j)
          else
             z_ortho=z_ortho+p(j)
          endif
          !
          !     collisional exciation/de-excitation
          !     H2-H collision
          !            Galli & Palla (1998)
          if(j == 0) then
             gamma_H=2.93d-14+1.21d-15*T_K &
                  +2.16d-19*T_K**2+1.32d-21*T_K**3
             !               if(iv==0) print *,"gamma_H(0) ", gamma_H*1.d-6
          elseif(j == 1) then
             gamma_H=8.34d-14+5.97d-16*T_K &
                  +7.76d-19*T_K**2+1.72d-21*T_K**3
             !               if(iv==0) print *,"gamma_H(1)", gamma_H
          elseif(j == 2) then
             gamma_H=7.54d-14+2.65d-16*T_K &
                  +2.14d-19*T_K**2+2.65d-21*T_K**3
          elseif(j == 3) then
             gamma_H=2.95d-14+1.21d-16*T_K &
                  +5.53d-19*T_K**2+2.51d-21*T_K**3
          else
             !     HM(1989) (Holenbach and Macee)
             gamma_H=4.6d-12*(2.d0*dble(j+2)-3.d0)*dsqrt(T_K)         &
                  *dsqrt(1.d0+2.d0*85.25d0*(2.d0*dble(j+2)-1.d0)/T_K)  &
                  *dexp(-(1.d1*85.25d0*(2.d0*dble(j+2)-1.d0))          &
                  /(T_K+85.25d0*dble(j+2)*dble(j+3))                   &
                  -0.1187d0*(4.d0*dble(j+2)-2.d0))
          endif
          gamma_H_j2down_omukai(j)=gamma_H

          !     H2-H2 collision (Holenbach and Macee 1979)
          gamma_H2=(3.3d-12+6.6d-15*T_K) &
               *0.276d0*dble((j+2)**2)   &
               *dexp(-(dble(j+2)/3.18d0)**1.7d0)

          !     H2-Hp collision : Gerlich(1990)
          xmeV=11.604d0
          if(j == 0) then 
             !        2-0
             gamma_Hp=2.889d-10*dexp(0.171d0*xmeV/T_K)
          elseif(j == 1) then
             !        3-1
             gamma_Hp=8.202d-10*dexp(0.250d0*xmeV/T_K)
          elseif(j == 2) then
             !        4-2
             gamma_Hp=5.700d-10*dexp(0.159d0*xmeV/T_K)
          elseif(j == 3) then
             !        5-3
             gamma_Hp=8.883d-10*dexp(0.223d0*xmeV/T_K)
          elseif(j == 4) then
             !        6-4
             gamma_Hp=5.209d-10*dexp(0.124d0*xmeV/T_K)
          elseif(j == 5) then
             !        7-5
             gamma_Hp=7.977d-10*dexp(0.177d0*xmeV/T_K)
          elseif(j == 6) then
             !        8-6
             gamma_Hp=4.489d-10*dexp(0.101d0*xmeV/T_K)
          else
             !        9-7
             gamma_Hp=5.915d-10*dexp(0.094d0*xmeV/T_K)
          endif

          !     H2-e collisional rate from Draine(1990)
          xj=dble(j)
          x=(ET(iv,j+2)-ET(iv,j))/T_K
          gamma_e=1.d-10*(xj+2.d0)*(xj+1.d0)/(2.d0*xj+5.d0)*(1.d0+1.5d0/x)

          C_ul=gamma_H*xn_a+gamma_H2*xn_m &
               +gamma_Hp*xn_Hp+gamma_e*xn_e
          DT=ET(iv,j+2)-ET(iv,j)
          Q_ul=Q_bg(DT)
          g_u=2.d0*dble(j)+5.d0
          g_l=2.d0*dble(j)+1.d0
          r=(g_u/g_l)*(A(iv,iv,j+2,j)*Q_ul+C_ul*dexp(-DT/T_K)) &
               /(A(iv,iv,j+2,j)*(1.d0+Q_ul)+C_ul)
          p(j+2)=p(j)*r

       enddo

       do j=0,20
          if(mod(j,2).eq.0) then
             f(iv,j)=(p(j)/z_para)*f_para*f_v(iv)
          else
             f(iv,j)=(p(j)/z_ortho)*f_ortho*f_v(iv)
          endif
       enddo
    enddo

    v_th=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))

    xLmbd_H2=0.d0
    do ji=0,18
       do ivi=1,2
          do ivf=0,ivi-1
             jf=ji+2
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l) &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do

    do ji=1,20
       do ivi=1,2
          do ivf=0,ivi-1
             jf=ji
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l)                     &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do


    do ji=2,20
       do ivi=0,2
          do ivf=0,ivi
             jf=ji-2
             DT=ET(ivi,ji)-ET(ivf,jf)
             DE=DT*xk_B
             xnu_Hz=DE/h_P
             xNc_l=xNc_H2*f(ivf,jf)
             xNc_u=xNc_H2*f(ivi,ji)
             g_l=dble(2*jf+1)
             g_u=dble(2*ji+1)
             Q_ul=Q_bg(DT)
             tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3 &
                  *(xNc_l*g_u/g_l-xNc_u)/v_th
             esc=beta_esc(tau_ul,tau_cnt)
             if(f(ivi,ji).ne.0.d0) then
                xLmbd_H2=xLmbd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc &
                     *(1.d0-Q_ul*((g_u/g_l)    &
                     *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
             endif
          end do
       end do
    end do

    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 ! per unit vol. (KS TODO)

    !!      print *, "T_k", T_k

    !      print *, "f(0,j)", f(0,0:10)
    !      ivi=0;ji=2;ivf=0;jf=0
    !      DT=ET(ivi,ji)-ET(ivf,jf)
    !      DE=DT*xk_B
    !      print *, "Omukai"
    !      print *, "ivi,ji,ivf,jf", ivi,ji,ivf,jf
    !      print *, "f(ivi,ji), A(ivi,ivf,ji,jf),DT", 
    !     &     f(ivi,ji), A(ivi,ivf,ji,jf),DT
    !      print *, "f*A*DE ",f(ivi,ji)*A(ivi,ivf,ji,jf)*DE
    !      print *, "f*A*DE/xnH", xLmbd_H2*xnH
    !      read(*,'()')
    !      print *, "gamma_H_j2down_omukai(0:2)"
    !      print *, gamma_H_j2down_omukai(0:2)*1.d-6
    return 
  end subroutine H2cool_Omukai



  ! Borysow, Frommhold, and Moraldi (1989) ApJ,336,495
  ! Equation (A11) (in K)
  function E_H2_BFM(iv,J)
    integer,intent(IN) :: iv,J
    real(kind=DBL_KIND) :: E_H2_BFM

    real(kind=DBL_KIND) :: vv,rrot,Ev,Ev0,Ev1,Ev2,Ev3,Ev4,Bv,Bv0,Bv1,Bv2,Bv3,Bv4,&
         Dv,Dv0,Dv1,Dv2,Dv3,Dv4,Fv,Fv0,Fv1,Fv2,Fv3,Fv4,Gv,Gv0,Gv1,Gv2,Gv3,Gv4,&
         Hv,Hv0,Hv1,Hv2,Hv3,Hv4,Ov,Ov0,Ov1,Ov2,Ov3,Ov4,Pv,Pv0,Pv1,Pv2,Pv3,Pv4

    Ev0=0.38496d0
    Ev1=-0.04609d0
    Ev2=0.00178d0
    Ev3=-7.d-5
    Ev4=2.9511d-6

    Bv0=54.438d0
    Bv1=4.6063d0
    Bv2=-2.0050d0
    Bv3=0.19260d0
    Bv4=-6.2953d-3

    Dv0=-0.72593d0
    Dv1=5.9990d0
    Dv2=-1.5187d0
    Dv3=0.12721d0
    Dv4=-3.3391d-3

    Fv0=-12.662d0
    Fv1=20.047d0
    Fv2=-4.7873d0
    Fv3=0.36900d0
    Fv4=-9.003d-3

    Gv0=-24.006d0
    Gv1=33.989d0
    Gv2=-7.6841d0
    Gv3=0.52413d0
    Gv4=-1.1297d-2

    Hv0=-22.384d0
    Hv1=31.150d0
    Hv2=-6.6139d0
    Hv3=0.39226d0
    Hv4=-9.496d-3

    Ov0=-10.541d0
    Ov1=14.746d0
    Ov2=-2.9476d0
    Ov3=0.16016d0
    Ov4=-6.5005d-3

    Pv0=-2.0021d0
    Pv1=2.8400d0
    Pv2=-0.54654d0
    Pv3=0.031636d0
    Pv4=-2.4398d-3

    vv=dble(iv)+0.5d0
    Ev=Ev0+Ev1*vv+Ev2*vv**2+Ev3*vv**3+Ev4*vv**4
    Bv=Bv0+Bv1*vv+Bv2*vv**2+Bv3*vv**3+Bv4*vv**4
    Dv=Dv0+Dv1*vv+Dv2*vv**2+Dv3*vv**3+Dv4*vv**4
    Fv=Fv0+Fv1*vv+Fv2*vv**2+Fv3*vv**3+Fv4*vv**4
    Gv=Gv0+Gv1*vv+Gv2*vv**2+Gv3*vv**3+Gv4*vv**4
    Hv=Hv0+Hv1*vv+Hv2*vv**2+Hv3*vv**3+Hv4*vv**4
    Ov=Ov0+Ov1*vv+Ov2*vv**2+Ov3*vv**3+Ov4*vv**4
    Pv=Pv0+Pv1*vv+Pv2*vv**2+Pv3*vv**3+Pv4*vv**4

    rrot=dble(J*(J+1))
    E_H2_BFM=1.43879d0* &
         (-Ev*1.d5+Bv*rrot-Dv*1.d-2*rrot**2+Fv*1.d-5*rrot**3 &
         -Gv*1.d-8*rrot**4+Hv*1.d-11*rrot**5-Ov*1.d-14*rrot**6 &
         +Pv*1.d-17*rrot**7)
    return
  end function E_H2_BFM


  !-----------------------------------------------------------------------
  ! LTE H2 cooling fit by Glover(2015) and fesc fit by Omukai-san
  !
  !     xLmbd_H2: H2 cooling rate per unit volume
  !-----------------------------------------------------------------------
  subroutine H2cool_LTEfit(xnH,T_K,y_H2,xNc_H2,tau,T_rad,xLmbd_H2)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y_H2,xNc_H2,tau,T_rad
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_H2

    real(kind=DBL_KIND) :: fesc,xlgT3,a,xLmbd,ar,xLmbdr

    fesc=0.d0
    if(T_K > 1.d2) then
       xlgT3=dlog10(T_K*1.d-3)
       a=-20.584225d0+5.0194035d0*xlgT3-1.5738805d0*xlgT3**2 &
            -4.7155769d0*xlgT3**3+2.4714161d0*xlgT3**4 &
            +5.4710750d0*xlgT3**5-3.9467356d0*xlgT3**6 &
            -2.2148338d0*xlgT3**7+1.8161874d0*xlgT3**8
       xLmbd=10.d0**a
    else
       xLmbd=0.d0
    endif
    if(T_rad > 1.d2) then 
       xlgT3=dlog10(T_rad*1.d-3)
       ar=-20.584225d0+5.0194035d0*xlgT3-1.5738805d0*xlgT3**2 &
            -4.7155769d0*xlgT3**3+2.4714161d0*xlgT3**4 &
            +5.4710750d0*xlgT3**5-3.9467356d0*xlgT3**6 &
            -2.2148338d0*xlgT3**7+1.8161874d0*xlgT3**8
       xLmbdr=10.d0**ar
    else
       xLmbdr=0.d0
    endif

    xLmbd_H2=(xLmbd-xLmbdr)/xnH

    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 ! per unit vol. (KS TODO)

    call f_fit(xNc_H2,T_K,fesc)
    xLmbd_H2=xLmbd_H2*fesc*dexp(-tau)
  end subroutine H2cool_LTEfit

  !  gives fitting value for fesc
  subroutine f_fit(xNc_H2,T_K,fesc)
    real(kind=DBL_KIND),intent(IN) :: xNc_H2,T_K
    real(kind=DBL_KIND),intent(OUT) :: fesc

    real(kind=DBL_KIND) :: a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5,alpha,beta,xlnf

    if(T_K < 1.d3) then
       a0= 0.978382d0
       a1= 0.000399572d0
       a2= -1.73108d-006
       a3= 1.15363d-009
       a4= 8.24607d-013
       a5= -7.65975d-016
    elseif(T_K <4.d3) then
       a0=0.827472d0         
       a1= 0.000107697d0      
       a2= -8.25123d-008    
       a3= 1.92812d-012     
       a4= 5.7192d-015      
       a5= -7.869d-019      
    else
       a0= 1.11569d0          
       a1= -0.000329302d0     
       a2= 1.01846d-007     
       a3= -1.46666d-011    
       a4= 1.00764d-015     
       a5= -2.68873d-020    
    endif
    alpha=a0+a1*T_K+a2*T_K**2+a3*T_K**3+a4*T_K**4+a5*T_K**5
    b0 = 24.0561d0     
    b1 = 0.00110043d0  
    b2 = -2.87224d-007 
    b3 = 6.11525d-011  
    b4 = -6.55034d-015 
    b5 = 2.54997d-019    
    beta=b0+b1*T_K+b2*T_K**2+b3*T_K**3+b4*T_K**4+b5*T_K**5

    if(beta < 200.d0) then
      xlnf=-alpha*dlog10(1.d0+xNc_H2/10.d0**beta)
      fesc=10.d0**xlnf
    else
      fesc = 1.d0
    endif
  end subroutine f_fit

  !-----------------------------------------------------------------------
  ! low-density H2 cooling fit by Glover(2015) 
  !
  !     xLmbd_H2 : H2 cooling rate per unit volume
  !
  !-----------------------------------------------------------------------
  subroutine H2cool_LowD_G15(xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He,xLmbd_H2)
    real(kind=DBL_KIND),intent(IN) :: xnH,T_K,y_H,y_H2,y_e,y_Hp,y_He
    real(kind=DBL_KIND),intent(OUT) :: xLmbd_H2

    real(kind=DBL_KIND) :: x_o,x_p,a0,a1,a2,a3,a4,a5,a6,a7,a8, &
         T3,xlgT3,xlgT3_2,xlgT3_3,xlgT3_4,xlgT3_5,xlgT3_6,xlgT3_7,xlgT3_8, &
         xLd_pH2H,xlgLd_pH2H,xLd_oH2H,xlgLd_oH2H,xLd_H2H,xlgLd_pH2pH2,xLd_pH2pH2, &
         xlgLd_pH2oH2,xLd_pH2oH2,xlgLd_oH2pH2,xLd_oH2pH2,&
         xlgLd_oH2oH2,xLd_oH2oH2,xLd_H2H2,&
         xlgLd_pH2He,xLd_pH2He,xlgLd_oH2He,xLd_oH2He,xLd_H2He, &
         xlgLd_H2Hp,xLd_H2Hp,xlgLd_H2e,xLd_H2e


    x_o=0.75d0
    x_p=0.25d0

    !     ortho-H2 excited by collisions with H
    if(T_K < 1000.d0) then
       a0=-24.330855d0
       a1=4.4404496d0
       a2=-4.0460989d0
       a3=-1.1390725d0
       a4=9.8094223d0
       a5=8.6273872d0
       a6=0.0d0
       a7=0.0d0
    else
       a0=-24.329086d0
       a1=4.6105087d0
       a2=-3.9505350d0
       a3=12.363818d0
       a4=-32.403165d0
       a5=48.853562d0
       a6=-38.542008d0
       a7=12.066770d0
    endif

    T3=T_K*1.d-3
    xlgT3=dlog10(T3)
    xlgT3_2=xlgT3**2
    xlgT3_3=xlgT3_2*xlgT3
    xlgT3_4=xlgT3_3*xlgT3
    xlgT3_5=xlgT3_4*xlgT3
    xlgT3_6=xlgT3_5*xlgT3
    xlgT3_7=xlgT3_6*xlgT3
    xlgT3_8=xlgT3_7*xlgT3
    !     
    if(T_K < 100.d0) then
       xLd_oH2H=5.09d-27*(T3**0.5d0)*dexp(-852.5d0/T_K)
    else
       xlgLd_oH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7
       xLd_oH2H=10.d0**xlgLd_oH2H
    endif

    !------------------------------------------------------!
    !     para-H2 excited by collisions with H
    if(T_K < 1000.d0) then
       a0=-24.216387d0
       a1=3.3237480d0
       a2=-11.642384d0
       a3=-35.553366d0
       a4=-35.105689d0
       a5=-10.922078d0
       a6=0.0d0
       a7=0.0d0
    else
       a0=-24.216387d0
       a1=4.2046488d0
       a2=-1.3155285d0
       a3=-1.6552763d0
       a4=4.1780102d0
       a5=-0.56949697d0
       a6=-3.3824407d0
       a7=1.0904027d0
    endif

    if(T_K < 100.d0) then
       xLd_pH2H=8.16d-26*(T3**0.5d0)*dexp(-509.85d0/T_K)
    else
       xlgLd_pH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7
       xLd_pH2H=10.d0**xlgLd_pH2H
    endif


    xLd_H2H=(x_o/(x_o+x_p))*xLd_oH2H +(x_p/(x_o+x_p))*xLd_pH2H

    !------------------------------------------------------!
    !     para-H2 with para-H2
    a0=-23.889798d0
    a1=1.8550774d0
    a2=-0.55593388d0
    a3=0.28429361d0
    a4=-0.20581113d0
    a5=0.13112378d0
    xlgLd_pH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2pH2=10.d0**xlgLd_pH2pH2

    !     para-H2 with ortho-H2
    a0=-23.748534d0
    a1=1.76676480d0
    a2=-0.58634325d0
    a3=0.31074159d0
    a4=-0.17455629d0
    a5=0.18530758d0
    xlgLd_pH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2oH2=10.d0**xlgLd_pH2oH2

    !     ortho-H2 with para-H2
    a0=-24.126177d0
    a1=2.3258217d0
    a2=-1.0082491d0
    a3=0.54823768d0
    a4=-0.33679759d0
    a5=0.20771406d0
    xlgLd_oH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2pH2=10.d0**xlgLd_oH2pH2

    !     ortho-H2 with ortho-H2
    a0=-24.020047d0
    a1=2.2687566d0
    a2=-1.0200304d0
    a3=0.83561432d0
    a4=-0.40772247d0
    a5=0.096025713d0
    xlgLd_oH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2oH2=10.d0**xlgLd_oH2oH2


    xLd_H2H2=x_p**2*xLd_pH2pH2+x_p*x_o*xLd_pH2oH2 &
         +x_o*x_p*xLd_oH2pH2+x_o**2*xLd_oH2oH2

    !------------------------------------------------------!
    !     para-H2 with He 
    a0=-23.489029d0
    a1=1.8210825d0
    a2=-0.59110559d0
    a3=0.42280623d0
    a4=-0.30171138d0
    a5=0.12872839d0
    xlgLd_pH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_pH2He=10.d0**xlgLd_pH2He
    !     ortho-H2 with He
    a0=-23.7749d0
    a1=2.40654d0
    a2=-1.23449d0
    a3=0.739874d0
    a4=-0.258940d0
    a5=0.120573d0
    xlgLd_oH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
         +a4*xlgT3_4+a5*xlgT3_5
    xLd_oH2He=10.d0**xlgLd_oH2He


    xLd_H2He=(x_o/(x_o+x_p))*xLd_oH2He &
         +(x_p/(x_o+x_p))*xLd_pH2He


    !------------------------------------------------------!
    !     H2 with proton
    if(10.d0<T_K) then
       a0=-22.089523d0
       a1=1.5714711d0
       a2=0.015391166d0
       a3=-0.23619985d0
       a4=-0.51002221d0
       a5=0.32168730d0
       xlgLd_H2Hp=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5
       xLd_H2Hp=10.d0**xlgLd_H2Hp
    else
       xLd_H2Hp=0.d0
    endif

    !------------------------------------------------------!
    !     H2 with e (KS MODIFIED)
    if(100.d0<T_K .and. T_K<500.d0) then
       a0=-21.928796d0 
       a1=16.815730d0  
       a2=96.743155d0
       a3=343.19180d0
       a4=734.71651d0
       a5=983.67576d0
       a6=801.81247d0
       a7=364.14446d0
       a8=70.609154d0
       xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 &
            +a7*xlgT3_7+a8*xlgT3_8
       xLd_H2e=10.d0**xlgLd_H2e
    else if(500.d0<T_K) then
       a0=-22.921189d0
       a1=1.6802758d0 
       a2=0.93310622d0
       a3=4.0406627d0
       a4=-4.7274036d0
       a5=-8.8077017d0
       a6=8.9167183d0
       a7=6.4380698d0
       a8=-6.3701156d0
       xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 &
            +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 &
            +a7*xlgT3_7+a8*xlgT3_8
       xLd_H2e=10.d0**xlgLd_H2e
    else
       xLd_H2e=0.d0 
    endif


    !------------------------------------------------------!
    ! total
    xLmbd_H2=xLd_H2H*y_H+xLd_H2H2*y_H2+xLd_H2He*y_He &
         +xLd_H2Hp*y_Hp+xLd_H2e*y_e

    xLmbd_H2 = xLmbd_H2 * xnH**2 * y_H2 ! per unit vol. (KS TODO)      
  end subroutine H2cool_LowD_G15


  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$        COOLING SOLVER          $$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !-------------------------------------------------------------------------
  ! explicit update of both chemistry and temperature
  ! グリッド単位で時間発展
  ! 温度は内部エネルギーに着目して時間発展させる
  ! i,j,k配列のサイズはどの変数も同じとする
  !
  !----- WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  ----
  !  解離光子が届いている場合、H2形成と解離が何桁もキャンセルしていることがあり、
  !  解離光子フラックスが10パーセント変わるだけでexplicitな発展が不可能になる
  !  なので以下の関数は作ったモノの使い物にならなかったのでコメントアウト
  !  ただし、explicitな発展では無く、iteration1回のimplicitな発展とかにしたら
  !  うまく回るかもしれないが、そのようなコード書き換えは今後の可能性（KS TODO）
  !----- WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  ----
  !-------------------------------------------------------------------------
  
  ! subroutine prim_GridChemCoolExplicit(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt,flg_recalc)
  !   use grid 
  !   real(kind=DBL_KIND),dimension(Imin:,Jmin:,Kmin:),intent(IN) :: xnH,rHpi,rH2pd,rHmpd,heat,xNc_H2
  !   real(kind=DBL_KIND),intent(IN) :: dt
  !   real(kind=DBL_KIND),dimension(Imin:,Jmin:,Kmin:),intent(INOUT) :: T_K
  !   real(kind=DBL_KIND),dimension(Imin:,Jmin:,Kmin:,0:),intent(INOUT) :: y
  !   logical,dimension(Imin:,Jmin:,Kmin:),intent(OUT) :: flg_recalc

  !   integer :: i,j,k,isp
  !   real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK,0:NCHEM-1) :: r_f          !reaction rate
  !   real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK,0:NREACT-1) :: xk           !rate coefficient
  !   real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: xLmbd_tot
  !   real(kind=DBL_KIND) :: xLmbd_line, xLmbd_chem
  !   real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK,0:NCHEM-1) :: y_new
  !   real(kind=DBL_KIND),dimension(ARRAYSIZE_IJK) :: T_K_new
  !   real(kind=DBL_KIND) :: xmu, en, gamma, pr, delta_y_max, delta_T

  !   !initialization
  !   flg_recalc(:,:,:) = .False.     !再計算が必要かどうかのフラグ

  !   !get ydot and lambda
  !   do i = lbound(xnH,1),ubound(xnH,1)
  !      do j = lbound(xnH,2),ubound(xnH,2)
  !         do k = lbound(xnH,3),ubound(xnH,3)
  !            call react_coef(xnH(i,j,k),T_K(i,j,k),&                                     !get xk
  !                 rHpi(i,j,k),rH2pd(i,j,k),rHmpd(i,j,k),xk(i,j,k,:))  
  !            call react_rat_highspeed(xk(i,j,k,:),xnH(i,j,k),y(i,j,k,:),r_f(i,j,k,:))    !get ydot
  !            call chem_cool(xnH(i,j,k),T_K(i,j,k),y(i,j,k,:),xk(i,j,k,:),heat(i,j,k),xLmbd_chem)
  !            call line_cool(xnH(i,j,k),T_K(i,j,k),y(i,j,k,:),xNc_H2(i,j,k),xLmbd_line)
  !            xLmbd_tot(i,j,k) = xLmbd_line + xLmbd_chem                                  !get lambda
  !         end do
  !      end do
  !   end do

  !   do i = lbound(xnH,1),ubound(xnH,1)
  !      do j = lbound(xnH,2),ubound(xnH,2)
  !         do k = lbound(xnH,3),ubound(xnH,3)
  !            !get en from xnH and T_K
  !            xmu = (1.d0 + 4.d0*yHe) &
  !                 /(y(i,j,k,0)+y(i,j,k,1)+y(i,j,k,2)+y(i,j,k,3)+yHe)  !平均分子量 (Hm, H2pは無視)
  !            pr = xnH(i,j,k)*(1.d0 + 4.d0*yHe)/xmu * boltz * T_K(i,j,k)             !圧力 pr = rho/(mH*xmu) kT
  !            gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
  !                 /(xmu*(1.5d0*(y(i,j,k,0)+y(i,j,k,2)+y(i,j,k,3)+yHe) + c_H2(T_K(i,j,k))*y(i,j,k,1)))
  !            en = pr / (gamma-1.)                                      !内部エネルギー密度 [erg cm^-3]

  !            !update y and en
  !            y_new(i,j,k,:) = y(i,j,k,:) + r_f(i,j,k,:)*dt
  !            en = en + xLmbd_tot(i,j,k)*dt

  !            !get T_K from xnH and en
  !            xmu = (1.d0 + 4.d0*yHe) &
  !                 /(y_new(i,j,k,0)+y_new(i,j,k,1)+y_new(i,j,k,2)+y_new(i,j,k,3)+yHe)   ! 平均分子量 (Hm, H2pは無視)
  !            gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
  !                 /(xmu*(1.5d0*(y_new(i,j,k,0)+y_new(i,j,k,2)+y_new(i,j,k,3)+yHe) + c_H2(T_K(i,j,k))*y_new(i,j,k,1)))    
  !            pr = (gamma-1.) * en
  !            T_K_new(i,j,k) = pr / (xnH(i,j,k)*(1.d0 + 4.d0*yHe)/xmu * boltz)     !new temperature

  !            !get delta_y_max and delta_T_K
  !            delta_y_max = 0d0
  !            do isp = 0, NCHEM-1
  !               if(y(i,j,k,isp) > 1d-15 .or. y_new(i,j,k,isp) > 1d-15) then
  !                  delta_y_max = max(delta_y_max,abs((y_new(i,j,k,isp) - y(i,j,k,isp))/y(i,j,k,isp)))
  !               end if
  !            end do
  !            delta_T = abs((T_K_new(i,j,k) - T_K(i,j,k))/T_K(i,j,k))

  !            !KS DEBUG
  !            if (globdbg_myrank == globdbg_rank .and. globdbg_mygid == globdbg_gid &
  !                 .and. i-Imin+Ngh == globdbg_i .and. j-Jmin+Ngh == globdbg_j &
  !                 .and. k-Kmin+Ngh == globdbg_k) then
  !               print '(/,A)', "*** KS DEBUG (in GridChemCoolExpl) ***"
  !               print '(A, 5I6)', "rank,gid,i,j,k = ", globdbg_myrank,globdbg_mygid, &
  !                    i-Imin+Ngh, j-Jmin+Ngh, k-Kmin+Ngh
  !               print '(A, (1P10E15.7))', "xnH, T_K, rHpi, rH2pd, rHmpd, dt:",&
  !                    xnH(i,j,k),T_K(i,j,k),rHpi(i,j,k),rh2pd(i,j,k),rHmpd(i,j,k), dt
  !               print '(A, (1P10E15.7))', "y: ", y(i,j,k,:)
  !               print '(A,(1P2E15.7))', "heat, xNc_H2: ", heat(i,j,k), xNc_H2(i,j,k)
  !               print '(/,A, (1P10E15.7))', "T_K_new: ", T_K_new(i,j,k)
  !               print '(A, (1P10E15.7))', "r_f: ", r_f(i,j,k,:)
  !               print '(A, (1P10E15.7))', "r_f*dt: ", r_f(i,j,k,:)*dt                
  !               print '(A, (1P10E15.7))', "y_new: ",y_new(i,j,k,:)                
  !            end if


  !            ! update y and T_K
  !            y(i,j,k,:) = y_new(i,j,k,:)
  !            call adjust_abundance(y(i,j,k,:))
  !            T_K(i,j,k) = T_K_new(i,j,k)
             
  !            !再計算が必要かをチェック
  !            if (delta_y_max > 0.1 .or. delta_T > 0.1) then !変化が大きい場合には
  !               flg_recalc(i,j,k) = .True.       !あとでimplicitに再計算するためにフラグを立てておく
  !            end if
  !         end do
  !      end do
  !   end do


  ! end subroutine prim_GridChemCoolExplicit
    

  !-------------------------------------------------------------------------
  ! implicit update of chemistry but explicit update of temperature
  ! 温度は内部エネルギーに着目して時間発展させる
  ! HFADDED: 190424 use xlmbdj
  !-------------------------------------------------------------------------
  subroutine CoolSolverExplicit(xch,t_chemcool,t_cool)
    type(chem_vals) :: xch
    real(kind=DBL_KIND),intent(INOUT) :: t_chemcool, t_cool
    real(kind=DBL_KIND) :: xmu, gamma, c_v,  t_chem, xLmbd_tot, en, pr, cs, T_new, rho, radius, chi_d, Tg
    real(kind=DBL_KIND),dimension(0:NREACT-1) :: xk    

    !get en from xnH and T_K
    xmu = (1.d0 + 4.d0*yHe) /(xch%ychem(0)+xch%ychem(1)+xch%ychem(2)+xch%ychem(3)+yHe)        !平均分子量 (Hm, H2pは無視)
    pr = xch%nH*(1.d0 + 4.d0*yHe)/xmu * boltz * xch%Tg             !圧力 pr = rho/(mH*xmu) kT
    gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))
    en = pr / (gamma-1.)                                      !内部エネルギー密度 [erg cm^-3]

    rho = ((1.d0 + 4.d0*yHe)*CONST_mH)*xch%nH !HFADDED 

!    print '(A,(1P5E15.7))', "(BEFORE) T_K,pr,xmu,gamma,xch%nH: ",T_K,pr,xmu,gamma,xch%nH !KS DEBUG    

    !get dust temp
#ifdef METAL
    radius = 0.5d0*xch%xlmbdj
    !call dust_temp(xch%nH, xch%Tg, Tcmb, xch%Td, radius, rho, xch%rdph)
    call dtemp_radtr(xch%nH, xch%Tg, xch%Td, xch%EradIR, xch%chi_d, dph_in=xch%rdph)
#else
    xch%Td = xch%Tg
#endif

    !pre-calculate xk here to speed-up of the code
    call react_coef(xch%nH,xch%Tg,xch%rHpi,xch%rH2pd,xch%rHmpd, xch%Td, xch%fd, xk)    

    !update y
    !call chemreact(xnH,T_K,y,rHpi,rH2pd,rHmpd,dt,t_chem)
    call chemreact_adptv(xch%nH,xch%Tg,xch%Td,xch%ychem,xch%yco,xch%rHpi,xch%rH2pd,xch%rHmpd &
      ,xch%fd,xch%metal,xch%dt,t_chem,xk_in=xk)

    !update yco
    call update_yco(xch%nH,xch%Tg,xch%ychem,xch%yco,xch%rcopd,xch%metal,xch%dt) 

    !get xLmbd_tot
    Tg = xch%Tg
    call tot_cool(xch,Tg,xLmbd_tot,xk_in=xk)
    


    !update en in response to xLmbd_tot
    en = en - xLmbd_tot * xch%dt


    !get T_K from xnH and en
    xmu = (1.d0 + 4.d0*yHe) /(xch%ychem(0)+xch%ychem(1)+xch%ychem(2)+xch%ychem(3)+yHe)   ! 平均分子量 (Hm, H2pは無視)
    gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))    
    pr = (gamma-1.) * en
    T_new = pr / (xch%nH*(1.d0 + 4.d0*yHe)/xmu * boltz)     !temperature update

    !update T
    xch%Tg = T_new

!    print '(A,(1P5E15.7))', "(AFTER) xch%Tg,pr,xmu,gamma,xnH: ",xch%Tg,pr,xmu,gamma,xnH !KS DEBUG    

    !timescales
    t_cool = abs(en/xLmbd_tot)             !cooling timescale in sec
    t_chemcool = min(t_cool,t_chem)        !cooling and chemical timescale in sec

    ! print '(A,(1P2E15.7))', "t_cool, t_chem: ",t_cool,t_chem !KS DEBUG

    
  end subroutine CoolSolverExplicit

  !-------------------------------------------------------------------------
  ! implicit update of chemistry and temperature
  !
  ! 温度は比熱と加熱冷却率に従って時間発展させる
  !
  !                use IMPLICIT Newton-Rhapson/bisection integrator
  ! HFADDED: 190424 use xlmbdj
  !-------------------------------------------------------------------------
 subroutine CoolSolverImplicit(xch)

    type(chem_vals) :: xch 
    real(kind=DBL_KIND),parameter :: err_eps=1d-2         !iteration ends when err_Tn < err_eps (KS TODO)
    real(kind=DBL_KIND),parameter :: T_min=2.d0, T_max=1d5 !Range of T [K] where T_n is searched (MP_Tmin/MP_Tmaxと揃える？KS TODO)
    integer,parameter ::  k_nr_max = 30                   !maximum iteration for NR method
    integer,parameter ::  k_srch_max = 100               !maximum iteration for lower/upper limit search
    integer,parameter ::  k_bs_max = 100                   !maximum iteration for BS method
    
    integer :: k_nr, k_bs, k_srch    
    real(kind=DBL_KIND) :: G, G0, G1,G_o, err_G, deno,  T_o, T_n, T_n0, T_n1, dT_K, err_Tn, y_o(0:NCHEM-1)
    logical :: force_substep !chemreact_adptvのsubstepを切ったり切らなかったりででcoolingのiterationで収束しないことを防ぐ

    real(kind=DBL_KIND) :: gfuv, dph, yco_o, drcopd
    real(kind=DBL_KIND),dimension(0:NREACT-1) :: xk        

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$            Initialization        $$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    force_substep = .False.

    !save original values
    T_o = xch%Tg
    y_o(:) = xch%ychem(:)
    yco_o  = xch%yco

#ifdef METAL 
    gfuv = xch%rgfuv
    dph  = xch%rdph
    drcopd = xch%rcopd
#else
    gfuv = 0.d0
    dph  = 0.d0
    drcopd = 0.d0
#endif

    !update T_K and y with rates@T_o
    call UpdateTY_with_Tn(xch, T_o)       

    !=============!
    G_o = T_o - xch%Tg
    !=============!

    !KS DEBUG
    if (dbg_flg_prm == 1) then
       print '(A,(1P3E15.7),/)', "(CoolSolverImplicit, beginning) T_o, T_K, G_o",&
            T_o, xch%Tg, G_o
    end if

    !updateしてみて温度変化が大きそうだったら、chemreact_adptvでsubstepを必ず切ることにして計算し直す
    if (abs(G_o/T_o) > 0.5) then
       !KS DEBUG
       if (dbg_flg_prm == 1) print '(A,(1P1E15.7),A)', "(CoolSolverImplicit) abs(G_o/T_o) = ", &
            abs(G_o/T_o), " > 0.5    =>    force_substep"
       force_substep = .True.
       xch%Tg = T_o       
       xch%ychem(:) = y_o(:)
       xch%yco  = yco_o
       call UpdateTY_with_Tn(xch,T_o,force_substep=force_substep)       
       G_o = T_o - xch%Tg
    end if
    
    

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$     Newton-Rapson iteration      $$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    !initialization for NR method
    T_n0 = T_o
    dT_K   = T_o*1.d-3
    T_n  = T_o + dT_K
    G0 = G_o
    G = G_o

    !begin NR iteration
    do k_nr=0, k_nr_max -1

       !upddate T_K and y with rates@T_n
       xch%Tg = T_o       
       xch%ychem(:) = y_o(:)
       xch%yco  = yco_o
       call UpdateTY_with_Tn(xch,T_n,force_substep=force_substep)       

       !=============!
       G = T_n - xch%Tg    !implicit法はこのGをzeroにすることを目指す（T_nから求めたTdotでT_o --> T_nに発展） (KS MEMO)             !=============!
       !=============!

       !数回やっても残差が大きい時はNRじゃなくてbisec使う (KS MODIFIED)
       if (abs(G/T_o) > 1. .and. k_nr >= 2) then
          !KS DEBUG
          if (dbg_flg_prm == 1) &
               print '(A,/,I0,(1P4E15.7),A)',"(CoolSolverImplicit, NR) abs(G/T_o) >1. with k_nr >= 2",&
               k_nr, T_o, T_n, xch%Tg, G, "    =>     move to bisection method"
          exit
       end if

       !NR法での次の位置
       deno = G - G0
       dT_K   = G*(T_n0 - T_n)/deno 

       !if(abs(dT_K) > 0.05*T_n0) dT_K = 0.05*T_n0*dT_K/abs(dT_K) !動きすぎないよう制限 (とりあえずコメントアウト, KS TODO)
       T_n0 = T_n
       G0   = G
       T_n    = T_n + dT_K                  !<--- T_n update
       if(T_n > T_max) T_n = T_max 
       if(T_n < T_min) T_n = T_min 

       !同じ値に張り付いてしまった場合 (e.g., floorに張り付いた場合) もbisec使う (KS MODIFIED)
       if (T_n == T_n0) then
          !KS DEBUG
          if (dbg_flg_prm == 1) &
               print '(A,/,I0,(1P4E15.7),A)',"(CoolSolverImplicit, NR) T_n == T_n0 occurs during NR iteration",&
               k_nr, T_o, T_n, xch%Tg, G, "    =>     move to bisection method"
          exit
       end if


       err_Tn = abs((T_n - T_n0)/T_n)     !<--- convergence of T update
       err_G = abs(G/T_o)                 !<--- residual of T

       !KS DEBUG
       if (dbg_flg_prm == 1) then
          print '(A,I0,(1P4E15.7))', "(CoolSolverImplicit, NR) k_nr, T_o, T_K, T_n0, T_n: ",&
               k_nr, T_o, xch%Tg, T_n0, T_n
          print '(A,(1P3E15.7),/)', "G, err_Tn, err_G: ",G, err_Tn, err_G
       end if

       !------ EXIT -------!
       if (err_Tn < err_eps .and. err_G < err_eps) then
          return
       end if
       !-------------------!
    end do

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$       bisection iteration        $$$$$$$$$$$$$$$$$$$$$$$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


    !--------------------------------------------------------------------------------------
    ! bisection methodで調べるT_nの範囲（下限、上限）を求める ==> その後2分法を用いてT_nを同定 (KS MEMO)
    ! (KS MODIFIED)
    !    ・二分法における上限(T_n0)、下限(T_n1)の決め方を効率化する
    !    ・上限下限が不必要に広がらないようにして、最終的なiterationの数を減らす
    !--------------------------------------------------------------------------------------

    ! print *, "*** Start bisec iteration ***"


    !+--+--+--+--+--+--+--+--+--+ determine range of T  +--+--+--+--+--+--+--+--+--+!
    !----------      冷却が優勢なら下限を探しに行く（上限は途中の値を利用） ----------!
    if (G_o>0.) then 
       ! print *, "*** try to find lower bound ***"
       T_n0 = T_o;
       do k_srch = 0, k_srch_max-1
          T_n1 = T_n0     !上限は途中の値を利用
          T_n0 = T_n0*0.9

          !update T_K and y with rates@T_n0
          xch%ychem(:) = y_o(:)
          xch%yco  = yco_o
          xch%Tg = T_o
          call UpdateTY_with_Tn(xch,T_n0,force_substep=force_substep)       

          !========================!
          G = T_n0 - xch%Tg
          !========================!

       !KS DEBUG
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec_lb) k_srch, T_o, T_n0, T_n1, G: ",k_srch, T_o, T_n0, T_n1, G
       end if


          if (G < 0.) then
             exit           !温度の下限が見つかったら終了
          end if

          if (T_n0 < T_min .or. k_srch == k_srch_max-1) then
             print *, "CLODE_IMP_bisec: lower_bound_not_found", T_o, T_n0, xch%Td, T_min
             print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
             print '(A,(1P10E15.7))', "y_init: ",y_o(:)
             print *, 'gfuv, dph, rcopd, rOII', gfuv, dph, xch%rcopd, xch%rOII
             print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
             print *, 'fd, EradIR, chi_d', xch%fd, xch%EradIR, xch%chi_d
             print *, "stopping..."

             if(k_srch == k_srch_max-1) then ! この場合だけ計算を止める
                stop !stop せずにwarningで進めてもいいかも (KS TODO)
             end if
          end if

       end do
       !--------------------------------------------------------------------------!
    else 
       !--------    加熱が優勢なら上限を探しに行く（下限は途中の値を利用）   -------------!
       ! print *, "*** try to find upper bound ***"
       T_n1 = T_o
       do k_srch = 0, k_srch_max-1
          T_n0 = T_n1     !下限は途中の値を利用
          T_n1 = T_n1*1.1

          !update T_K and y with rates@T_n1
          xch%ychem(:) = y_o(:)
          xch%yco  = yco_o
          xch%Tg = T_o          
          call UpdateTY_with_Tn(xch,T_n1,force_substep=force_substep)       
          !========================!
          G = T_n1 - xch%Tg
          !========================!

       !KS DEBUG
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec_ub) k_srch, T_o, T_n0, T_n1, G: ",k_srch, T_o, T_n0, T_n1, G
       end if

          if (G > 0.) then
             exit           !温度の下限が見つかったら終了
          end if

          if (T_n1 > T_max .or. k_srch == k_srch_max-1) then
             print *, "CLODE_IMP_bisec: upper_bound_not_found", T_o, T_n1, T_max
             print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
             print '(A,(1P10E15.7))', "y_init: ",y_o(:)
             print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
             print *, "stopping..."
             stop !stop せずにwarningで進めてもいいかも (KS TODO)
          end if
       end do
       !--------------------------------------------------------------------------!
    end if

    !+--+--+--+--+--+--+--+--+--+ Bisec iteration +--+--+--+--+--+--+--+--+--+!
    do k_bs = 0, k_bs_max-1
       !中間の値をとってくる
       T_n = 0.5*(T_n1 + T_n0)

       !update T_K and y with rates@T_n
       xch%ychem(:) = y_o(:)
       xch%yco  = yco_o
       xch%Tg = T_o
       call UpdateTY_with_Tn(xch,T_n,force_substep=force_substep)       

       !========================!
       G = T_n - xch%Tg       
       !========================!

       !幅[T_n0,T_n1]を狭めていく
       if (G < 0.) then 
          T_n0 = T_n  !T_nだと温度低すぎるとき、T_n2を新たにT_n0（下限）とする
       else
          T_n1 = T_n  !T_nだと温度高すぎるとき、T_n2を新たにT_n1（上限）とする
       end if

       err_Tn = abs((T_n1-T_n0)/T_n0)     !<--- convegence of T update
       err_G = abs(G/T_o)                 !<--- residual of T    

       !KS DEBUG
       if (dbg_flg_prm == 1) then
          print '(A,I0, (1P4E15.7))', "(CoolSolverImplicit, bisec) k_bs, T_o, T_K, err_Tn, err_G: ", k_bs, T_o, xch%Tg, err_Tn, err_G
          print '(A,(1P5E20.12))', " T_n0, T_n1, T_n, T_K, G: ", T_n0, T_n1, T_n, xch%Tg, G
       end if

       !------ EXIT -------!
       !if (err_Tn < err_eps .and. err_G < err_eps) then
       !err_Gを消すのは諦める、気持ちとしては一旦bisectionで近づけておいて、あとでNRで一致させる (KS TODO)
       !冷却率の温度依存性が非常に鋭く、かつ、化学組成がchemreactのsubstepによって微妙に変わるときなどに収束しなかったので、ひとまず変えてみた
       if (err_Tn < err_eps) then
          xch%Tg = T_n
          return
       end if
       !-------------------!

       !------------------- ITERATION LIMIT -----------------------!
       if (k_bs == k_bs_max-1) then
          print '(A,(1P3E15.7))', "CLODE_IMP_bisec: hitting the iteration limit",T_o, T_n0, T_n1
          print '(A,(1P10E15.7))', "xnH T_K rHpi rH2pd rHmpd dt: ",xch%nH,T_o,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%dt
          print '(A,(1P10E15.7))', "y_init: ",y_o(:)
          print '(A,(1P2E15.7))', "heat, xlmbdj: ", xch%heat, xch%xlmbdj
          print *, "stopping..."
          stop !stop せずにwarningで進めてもいいかも (KS TODO)
       end if
       !----------------------------------------------------------!
    end do
  end subroutine CoolSolverImplicit 

  !-------------------------------------------------------------------------
  ! update T_K and y using the reaction/cooling rates for T_n
  ! HFADDED: use xlmbdj
  !-------------------------------------------------------------------------
subroutine UpdateTY_with_Tn(xch, T_n, force_substep)
    type(chem_vals) :: xch 
    real(kind=DBL_KIND),intent(IN) :: T_n
    logical,intent(IN),optional :: force_substep
    real(kind=DBL_KIND),dimension(0:NREACT-1) :: xk !高速化のため反応係数@T_nを渡せるようにする

    real(kind=DBL_KIND) :: t_chem,xLmbd_tot,xmu,gamma,T_new,pr,en,t_cool, t_chemcool
    real(kind=DBL_KIND) :: rho, radius, Tdust

    logical :: bool_force_substep
    !強制的にsubstepを切るか
    bool_force_substep = .false.
    if ( present( force_substep) ) then
       bool_force_substep = force_substep
    endif

    ! dust temp
#ifdef METAL
    rho    = ((1.d0 + 4.d0*yHe)*CONST_mH)*xch%nH      ! [g cm^{-3}]
    radius = 0.5d0*xch%xlmbdj
    !call dust_temp(xch%nH, xch%Tg, Tcmb, xch%Td, radius, rho, xch%rdph)
    call dtemp_radtr(xch%nH, xch%Tg, xch%Td, xch%EradIR, xch%chi_d, dph_in=xch%rdph)
#else
    xch%Td = xch%Tg
#endif

    !pre-calculate xk@T_n here to speed-up of the code
    call react_coef(xch%nH,T_n,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%Td,xch%fd,xk)    

    !get en from xnH and T_K
    xmu = (1.d0 + 4.d0*yHe) /(xch%ychem(0)+xch%ychem(1)+xch%ychem(2)+xch%ychem(3)+yHe)        !平均分子量 (Hm, H2pは無視)
    pr = xch%nH*(1.d0 + 4.d0*yHe)/xmu * boltz *xch%Tg             !圧力 pr = rho/(mH*xmu) kT
    gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))
    en = pr / (gamma-1.)                                      !内部エネルギー密度 [erg cm^-3]

    !KS DEBUG
    if (dbg_flg_prm == 1) then
       print '(A,(1P7E15.7))', "(BEFORE chemreact) T_K,T_n,en,pr,xmu,gamma,xnH: ",xch%Tg,T_n,en,pr,xmu,gamma,xch%nH !KS DEBUG
    end if

    !update y
    !call chemreact(xnH,T_n,y,rHpi,rH2pd,rHmpd,dt,t_chem)
    if(xch%nH > 1e11 .or. bool_force_substep) then        !密度高いときやCoolSolver側からsubstepを切ることを要請された場合（温度変化が大きい場合に対応）
       !xk@T_nをchemreact_adptvに渡す
       call chemreact_adptv(xch%nH,T_n,xch%Td,xch%ychem,xch%yco,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%fd,xch%metal &
         ,xch%dt,t_chem, force_substep=.True., xk_in=xk)   !強制的にsubstepを切る
    else
       !xk@T_nをchemreact_adptvに渡す
       call chemreact_adptv(xch%nH,T_n,xch%Td,xch%ychem,xch%yco,xch%rHpi,xch%rH2pd,xch%rHmpd,xch%fd,xch%metal &
       ,xch%dt,t_chem,xk_in=xk)
    end if

    ! update co
    call update_yco(xch%nH,T_n,xch%ychem,xch%yco,xch%rcopd,xch%metal,xch%dt)

    !get xLmbd_tot (xk@T_nをtot_coolに渡す)
    call tot_cool(xch,T_n,xLmbd_tot,xk_in=xk)


    !update en in response to xLmbd_tot
    en = en - xLmbd_tot * xch%dt

    !get T_K from xnH and en
    xmu = (1.d0 + 4.d0*yHe) /(xch%ychem(0)+xch%ychem(1)+xch%ychem(2)+xch%ychem(3)+yHe)   ! 平均分子量 (Hm, H2pは無視)
    gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
         /(xmu*(1.5d0*(xch%ychem(0)+xch%ychem(2)+xch%ychem(3)+yHe) + c_H2(xch%Tg)*xch%ychem(1)))    
    pr = (gamma-1.) * en
    T_new = pr / (xch%nH*(1.d0 + 4.d0*yHe)/xmu * boltz)     !temperature update

    !update T
    xch%Tg = T_new

    !KS DEBUG
    if (dbg_flg_prm == 1) then
       print '(A,(1P7E15.7))', "(AFTER chemreact) T_K,T_n,en,pr,xmu,gamma,xnH: ",xch%Tg,T_n,en,pr,xmu,gamma,xch%nH !KS DEBUG
    end if

    !timescales
    t_cool = abs(en/xLmbd_tot)             !cooling timescale in sec
    t_chemcool = min(t_cool,t_chem)        !cooling and chemical timescale in sec

    ! print '(A,(1P2E15.7))', "t_cool, t_chem: ",t_cool,t_chem !KS DEBUG

  end subroutine UpdateTY_With_Tn

  ! !-------------------------------------------------------------------------
  ! ! call chemreact and then get Tdot by calling tot_cool
  ! !-------------------------------------------------------------------------
  ! subroutine chemreact_and_Tdot(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt,Tdot,t_chemcool)
  !   real(kind=DBL_KIND),intent(IN) :: xnH,T_K,rHpi,rH2pd,rHmpd,heat,xNc_H2,dt
  !   real(kind=DBL_KIND),intent(OUT) :: Tdot,t_chemcool
  !   real(kind=DBL_KIND),intent(INOUT) :: y(0:NCHEM-1)    
  !   real(kind=DBL_KIND) :: t_chem,xLmbd_tot,xmu,gamma,T_new,pr,en,t_cool

  !   !get en from xnH and T_K
  !   xmu = (1.d0 + 4.d0*yHe) /(y(0)+y(1)+y(2)+y(3)+yHe)        !平均分子量 (Hm, H2pは無視)
  !   pr = xnH*(1.d0 + 4.d0*yHe)/xmu * boltz * T_K             !圧力 pr = rho/(mH*xmu) kT
  !   gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
  !        /(xmu*(1.5d0*(y(0)+y(2)+y(3)+yHe) + c_H2(T_K)*y(1)))
  !   en = pr / (gamma-1.)                                      !内部エネルギー密度 [erg cm^-3]

  !   ! print '(A,(1P5E15.7))', "(BEFORE) T_K,pr,xmu,gamma,xnH: ",T_K,pr,xmu,gamma,xnH !KS DEBUG    


  !   !update y
  !   call chemreact(xnH,T_K,y,rHpi,rH2pd,rHmpd,dt,t_chem)

  !   !get xLmbd_tot
  !   call tot_cool(xnH,T_K,y,rHpi,rH2pd,rHmpd,heat,xNc_H2,xLmbd_tot)

  !   !update en in response to xLmbd_tot
  !   en = en - xLmbd_tot * dt

  !   print '(A,(1P3E15.7))', 'delta_P, P_o, P_new',xLmbd_tot * dt *(gamma-1.), pr, pr+xLmbd_tot * dt *(gamma-1.)


  !   !get T_K from xnH and en
  !   xmu = (1.d0 + 4.d0*yHe) /(y(0)+y(1)+y(2)+y(3)+yHe)   ! 平均分子量 (Hm, H2pは無視)
  !   gamma = 1.d0+(1.d0+4.d0*yHe) &                            !比熱比 (Hm, H2pは無視)
  !        /(xmu*(1.5d0*(y(0)+y(2)+y(3)+yHe) + c_H2(T_K)*y(1)))    
  !   pr = (gamma-1.) * en
  !   T_new = pr / (xnH*(1.d0 + 4.d0*yHe)/xmu * boltz)     !temperature update

  !   print '(A,(1P4E15.7))', 'T_o, T_new, delta_T, T_new/ T_o',T_K, T_new, T_new-T_K, T_new/ T_K

  !   !get time derivative
  !   Tdot = (T_new-T_K)/dt

  !   ! print '(A,(1P5E15.7))', "(AFTER) T_K,pr,xmu,gamma,xnH: ",T_K,pr,xmu,gamma,xnH !KS DEBUG    

  !   !timescales
  !   t_cool = abs(en/xLmbd_tot)             !cooling timescale in sec
  !   t_chemcool = min(t_cool,t_chem)        !cooling and chemical timescale in sec
  !   print '(A,(1P2E15.7))', "t_cool, t_chem: ",t_cool,t_chem !KS DEBUG

  ! end subroutine chemreact_and_Tdot

  


  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$            RADIATION           $$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  !----------------------------------------------------------------------!
  !                           < HUVAHHH >                                !
  !                                                                      !
  !     input T: stellar effective temperature (10**3 -- 10**5.5K)       !
  !     output                                                           !
  !      - huv: huv*xslum is emissivity of stellar EUV photons           !
  !      - alpha: mean cross section for EUV photons (/cm^2)             !
  !      - heat: mean energy transferred to gas                          !
  !              per one photoionizaiton event (erg)                     !
  !      - hve: hve*xslum is emissivity of stellar FUV photons           !
  !             ( 11.174 -- 13.6 eV )                                    !
  !      - hhm: nu-averaged H- dissociation cross section x pi*Bnu       !
  !                                                                      !
  !    < HUVAH > and < HMDSR > have been merged by KS (2014.7.18)        !
  !----------------------------------------------------------------------! 
  subroutine HUVAHHH(T, huv, alpha, heat, hve, hhm)
    real(kind=DBL_KIND),intent(IN) :: T
    real(kind=DBL_KIND),intent(OUT) :: huv, alpha, heat, hve, hhm

    real(kind=DBL_KIND),dimension(54) :: tbhuv,tbalph,tbheat,tbhved,tbhmds
    real(kind=DBL_KIND) :: tl1,tl2,dtl,tlg
    integer ::nt,it,it1
    
    
    data  tbhuv /-5.21206d+01,-4.52041d+01,-3.90139d+01,-3.34752d+01,&
         -2.85209d+01,-2.40910d+01,-2.01313d+01,-1.65935d+01,-1.34341d+01,&
         -1.06141d+01,-8.09864d+00,-5.85626d+00,-3.85882d+00,-2.08108d+00,&
         -5.00370d-01, 9.03649d-01, 2.14923d+00, 3.25274d+00, 4.22889d+00,&
         5.09087d+00, 5.85053d+00, 6.51851d+00, 7.10437d+00, 7.61669d+00,&
         8.06318d+00, 8.45079d+00, 8.78575d+00, 9.07370d+00, 9.31969d+00,&
         9.52832d+00, 9.70371d+00, 9.84963d+00, 9.96948d+00, 1.00663d+01,&
         1.01430d+01, 1.02020d+01, 1.02455d+01, 1.02757d+01, 1.02943d+01,&
         1.03028d+01, 1.03028d+01, 1.02954d+01, 1.02817d+01, 1.02627d+01,&
         1.02390d+01, 1.02116d+01, 1.01808d+01, 1.01473d+01, 1.01115d+01,&
         1.00737d+01, 1.00344d+01, 9.99361d+00, 9.95174d+00, 9.90893d+00/

    data tbalph /-1.72075d+01,-1.72084d+01,-1.72093d+01,-1.72103d+01,&
         -1.72115d+01,-1.72128d+01,-1.72142d+01,-1.72158d+01,-1.72176d+01,&
         -1.72196d+01,-1.72218d+01,-1.72243d+01,-1.72271d+01,-1.72301d+01,&
         -1.72335d+01,-1.72373d+01,-1.72415d+01,-1.72462d+01,-1.72515d+01,&
         -1.72573d+01,-1.72637d+01,-1.72709d+01,-1.72789d+01,-1.72877d+01,&
         -1.72976d+01,-1.73085d+01,-1.73206d+01,-1.73339d+01,-1.73487d+01,&
         -1.73651d+01,-1.73831d+01,-1.74030d+01,-1.74248d+01,-1.74486d+01,&
         -1.74747d+01,-1.75030d+01,-1.75336d+01,-1.75667d+01,-1.76022d+01,&
         -1.76403d+01,-1.76808d+01,-1.77239d+01,-1.77695d+01,-1.78176d+01,&
         -1.78681d+01,-1.79210d+01,-1.79762d+01,-1.80337d+01,-1.80935d+01,&
         -1.81554d+01,-1.82193d+01,-1.82852d+01,-1.83530d+01,-1.84226d+01/

    data tbheat /-1.28617d+01,-1.28148d+01,-1.27678d+01,-1.27209d+01,&
         -1.26740d+01,-1.26272d+01,-1.25803d+01,-1.25336d+01,-1.24868d+01,&
         -1.24401d+01,-1.23935d+01,-1.23469d+01,-1.23004d+01,-1.22540d+01,&
         -1.22076d+01,-1.21613d+01,-1.21151d+01,-1.20690d+01,-1.20230d+01,&
         -1.19772d+01,-1.19314d+01,-1.18858d+01,-1.18404d+01,-1.17951d+01,&
         -1.17500d+01,-1.17051d+01,-1.16604d+01,-1.16160d+01,-1.15718d+01,&
         -1.15278d+01,-1.14842d+01,-1.14410d+01,-1.13981d+01,-1.13558d+01,&
         -1.13140d+01,-1.12729d+01,-1.12325d+01,-1.11928d+01,-1.11539d+01,&
         -1.11159d+01,-1.10788d+01,-1.10426d+01,-1.10074d+01,-1.09731d+01,&
         -1.09397d+01,-1.09073d+01,-1.08758d+01,-1.08452d+01,-1.08155d+01,&
         -1.07867d+01,-1.07587d+01,-1.07315d+01,-1.07052d+01,-1.06796d+01/

    !< FUV from 11.174eV (1100A) to 13.6eV (912A) >!
    data tbhved /-4.00350E+01,-3.43797E+01,-2.93209E+01,-2.47971E+01,&
         -2.07532E+01,-1.71399E+01,-1.39128E+01,-1.10321E+01,-8.46216E+00,&
         -6.17092E+00,-4.12969E+00,-2.31271E+00,-6.96923E-01, 7.38349E-01,&
         2.01159E+00, 3.13933E+00, 4.13631E+00, 5.01571E+00, 5.78929E+00,&
         6.46755E+00, 7.05991E+00, 7.57479E+00, 8.01978E+00, 8.40170E+00,&
         8.72672E+00, 9.00041E+00, 9.22781E+00, 9.41352E+00, 9.56170E+00,&
         9.67618E+00, 9.76043E+00, 9.81768E+00, 9.85086E+00, 9.86267E+00,&
         9.85558E+00, 9.83184E+00, 9.79345E+00, 9.74225E+00, 9.67985E+00,&
         9.60767E+00, 9.52697E+00, 9.43886E+00, 9.34428E+00, 9.24407E+00,&
         9.13894E+00, 9.02954E+00, 8.91638E+00, 8.79995E+00, 8.68064E+00,&
         8.55880E+00, 8.43475E+00, 8.30875E+00, 8.18102E+00, 8.05178E+00/

    data tbhmds /  1.1613E+01,  1.2113E+01,  1.2571E+01,  1.2992E+01,&
        1.3380E+01,  1.3739E+01,  1.4070E+01,  1.4376E+01,  1.4662E+01,&
        1.4927E+01,  1.5175E+01,  1.5407E+01,  1.5624E+01,  1.5829E+01,&
        1.6023E+01,  1.6205E+01,  1.6379E+01,  1.6544E+01,  1.6701E+01,&
        1.6851E+01,  1.6995E+01,  1.7134E+01,  1.7267E+01,  1.7395E+01,&
        1.7519E+01,  1.7639E+01,  1.7756E+01,  1.7869E+01,  1.7979E+01,&
        1.8087E+01,  1.8191E+01,  1.8293E+01,  1.8391E+01,  1.8487E+01,&
        1.8579E+01,  1.8668E+01,  1.8754E+01,  1.8837E+01,  1.8917E+01,&
        1.8994E+01,  1.9068E+01,  1.9139E+01,  1.9208E+01,  1.9275E+01,&
        1.9340E+01,  1.9404E+01,  1.9465E+01,  1.9525E+01,  1.9584E+01,&
        1.9642E+01,  1.9698E+01,  1.9754E+01,  1.9808E+01,  1.9862E+01/
    
    data tl1,tl2,dtl,nt/ 3.0000d+00, 5.5000d+00, 4.716981132d-02,54/

    tlg = (min(max(tl1,log10(t+1.0e-37)),tl2)-tl1)/dtl
    it  = min(int(tlg)+1,nt)
    it1 = min(it+1,nt)
    tlg = tlg-int(tlg)

    huv = 10**(tbhuv(it)*(1.d0-tlg) + tbhuv(it1)*tlg)
    alpha = 10**(tbalph(it)*(1.d0-tlg) + tbalph(it1)*tlg)
    heat  = 10**(tbheat(it)*(1.d0-tlg) + tbheat(it1)*tlg)
    hve   = 10**(tbhved(it)*(1.d0-tlg)  + tbhved(it1)*tlg)
    hhm = 10**(tbhmds(it)*(1.d0-tlg) + tbhmds(it1)*tlg)
    
  end subroutine HUVAHHH

  !-------------------------------------------------------------------------
  ! get R_* and L_* from M_* and Mdot using a fitting of Prost data
  !-------------------------------------------------------------------------
  subroutine ProstFit2(Mass, Mdot, Rad, Lum, Trad, rda)
    use mpilib

    real(kind=DBL_KIND),intent(IN) :: Mass        ! mass in [M_sun]
    real(kind=DBL_KIND),intent(IN) :: Mdot        ! mdot in [M_sun/yr]
    real(kind=DBL_KIND),intent(Out) :: Rad        ! radius in [R_sun]
    real(kind=DBL_KIND),intent(Out) :: Lum        ! luminosity in [L_sun]
    real(kind=DBL_KIND),intent(Out) :: Trad       ! Effective temeprature [K]

    type(rad_others), intent(Out), optional :: rda

    !data file
    integer,parameter :: nmd=11  !number of mdot in data
    integer,parameter :: npt=250 !number of mass in data
    character(100) :: fname(0:nmd-1) = (/"0e00", "1e-5", "1e-4", "3e-4", "1e-3", "3e-3" &
      , "6e-3", "1e-2", "3e-2", "6e-2", "1e-1"/)
    character(100) :: path2fitdat = "./ProstFit/"//trim(FOL_RADSOURCE)//"/"
    integer,parameter :: FH = 11
    character(len=100) :: ffn 
    integer :: err 

        !data container
    real(kind=DBL_KIND),save :: xm(0:npt-1)                !mass in [M_sun]
    real(kind=DBL_KIND),save :: xr(0:npt-1,0:nmd-1)        !radius in [R_sun]
    real(kind=DBL_KIND),save :: xls(0:npt-1,0:nmd-1)       !luminosity in [L_sun]
    real(kind=DBL_KIND),save :: xtrad(0:npt-1,0:nmd-1)     !Trad [K]
    real(kind=DBL_KIND),save :: xxeuv(0:npt-1,0:nmd-1)     !xeuv
    real(kind=DBL_KIND),save :: xxfuv(0:npt-1,0:nmd-1)     !xfuv
    real(kind=DBL_KIND),save :: xaleuv(0:npt-1,0:nmd-1)    !alpha_euv
    real(kind=DBL_KIND),save :: xheat(0:npt-1,0:nmd-1)     !heat_euv
    real(kind=DBL_KIND),save :: xhhm(0:npt-1,0:nmd-1)      !hhm
    real(kind=DBL_KIND),save :: xlumeuv(0:npt-1,0:nmd-1)   !lumeuv
    real(kind=DBL_KIND),save :: xlumfuv(0:npt-1,0:nmd-1)   !lumfuv
    real(kind=DBL_KIND),save :: xsigeuv(0:npt-1,0:nmd-1)   !sigeuv
    real(kind=DBL_KIND),save :: xsigfuv(0:npt-1,0:nmd-1)   !sigfuv
    real(kind=DBL_KIND),save :: xrOII(0:npt-1,0:nmd-1)     !rOII
    real(kind=DBL_KIND),parameter ::  xmd(0:nmd-1) = (/0.d0, 1.d-5, 1.d-4, 3.d-4, 1.d-3, 3.d-3, 6.d-3,&
         1.d-2, 3.d-2, 6.d-2, 1.d-1/)          !mdot in [M_sun/yr]

    !other variables
    integer,save :: ifirst=0
    integer :: i, j, id, ix, iy
    real(kind=DBL_KIND) :: x_in, y_in, x1, x2, y1, y2, z_y1, z_y2
    real(kind=DBL_KIND) :: tap1, tap2, tap3, tap4

    !initialization (read data)
    if (ifirst == 0) then
       do id = 0, nmd-1
          ffn= trim(path2fitdat)//trim(fname(id))
          if(get_myrank() == PRIMARY_RANK) print *, "PROST_FIT: reading ", ffn
          open(FH, file=ffn,status='old',iostat=err)
          if (err > 0) then
             print '(A,/,A)', "PROST_FIT: file not found","stopping..."
             stop
          end if
          do i = 0, npt-1
             read(FH,*,iostat=err) xm(i), xr(i,id), xls(i,id), xtrad(i,id), xxeuv(i,id), xxfuv(i,id), xaleuv(i,id)&
               , xheat(i,id), xhhm(i,id), xlumeuv(i,id), xlumfuv(i,id), xsigeuv(i,id), xsigfuv(i,id), xrOII(i,id)
             if (err > 0) then
                print *, "PROST_FIT: **WARNING** data size might be inconsistent"
                exit !¿¿¿¿¿¿¿¿¿¿¿¿¿¿
             end if
          end do
          close(FH)
       end do
       ifirst = 1
    end if


    ! linear interpolation
    x_in    = dlog10(Mass)
    x_in    = min(max(x_in, xm(0)),  xm(npt-1))
    y_in    = min(max(Mdot,xmd(0)), xmd(nmd-1))

    ! get index of m & mdot
    ix  = int((x_in+2.0)/0.02d0)
    ix  = max(min(ix, npt-2),0)
    iy = nmd-2
    do j = 0, nmd-2
      if (y_in <= xmd(j+1)) then
        iy = j
        exit
      endif
    enddo


    !interpolation
    x1 = xm(ix)
    x2 = xm(ix+1)
    y1 = xmd(iy)
    y2 = xmd(iy+1)

    tap1 = (x2-x_in)/(x2-x1)
    tap2 = (x_in-x1)/(x2-x1)
    tap3 = (y2-y_in)/(y2-y1)
    tap4 = (y_in-y1)/(y2-y1)

#define INTERPOLATION_LIN(zval, zarr) \
    z_y1 = tap1*zarr(ix,iy)  + tap2*zarr(ix+1,iy)    ;\
    z_y2 = tap1*zarr(ix,iy+1)+ tap2*zarr(ix+1,iy+1)  ;\
    zval = tap3*z_y1 + tap4*z_y2; \
    zval = max(zval, 0.d0)

    INTERPOLATION_LIN( Rad , xr   )
    INTERPOLATION_LIN( Lum , xls  )
    INTERPOLATION_LIN( Trad, xtrad)

    if(present(rda)) then
      INTERPOLATION_LIN( rda%xeuv, xxeuv)
      INTERPOLATION_LIN( rda%xfuv, xxfuv)
      INTERPOLATION_LIN( rda%alpha_euv, xaleuv)
      INTERPOLATION_LIN( rda%heat_euv, xheat)
      INTERPOLATION_LIN( rda%hhm, xhhm)
      INTERPOLATION_LIN( rda%lumeuv, xlumeuv)
      INTERPOLATION_LIN( rda%lumfuv, xlumfuv)
      INTERPOLATION_LIN( rda%sig_euv, xsigeuv)
      INTERPOLATION_LIN( rda%sig_fuv, xsigfuv)
      INTERPOLATION_LIN( rda%rOII, xrOII)
    endif

    !print *, x_in, x1, x2, ix, y_in, iy, xeuv, xxeuv(ix, iy)!, x_in -x1, x2 - x_in
      !(x2-x_in)/(x2-x1)*xxeuv(ix,iy)+(x_in-x1)/(x2-x1)*xxeuv(ix+1,iy), &
      !(x2-x_in)/(x2-x1)*xxeuv(ix,iy+1)+(x_in-x1)/(x2-x1)*xxeuv(ix+1,iy+1)

#undef INTERPOLATION_LIN

  end subroutine ProstFit2

end module primordial
