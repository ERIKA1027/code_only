module mg
  use mg_data
  use fmg_boundary_phys
  implicit none
  real(kind=8),save :: Resmaxg
  private
  public :: mg_lin, mg_cycle, mg_boundary_fill0, mg_fill0, mg_nonlin, mg_nonlin_vcycle, mg_init, mg_finalize
contains
  subroutine mg_poisson_resid(mglev, jres, ju, jrhs)
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=8),pointer,dimension(:,:,:) :: res, u, rhs
    integer :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8) :: h2i
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    call mg_arrp(mglev,jres, res)
    call mg_arrp(mglev,ju, u)
    call mg_arrp(mglev,jrhs, rhs)
    h2i = 1.d0/mg_get_h( mglev ) ** 2
    call mg_boundary_u(mglev,ju)
    do k = ks, ke
       do j = js, je
          do i = is, ie
             res(i,j,k)=-h2i*( &
                   u(i+1,j,k)+u(i-1,j,k) &
                  +u(i,j+1,k)+u(i,j-1,k) &
                  +u(i,j,k+1)+u(i,j,k-1) &
                  -6.d0*u(i,j,k))+rhs(i,j,k)
          enddo
       enddo
    enddo
  end subroutine mg_poisson_resid
  subroutine mg_poisson_relax(mglev, ju, jrhs)
    integer,intent(IN) :: mglev, ju, jrhs
    real(kind=8),parameter :: sixth = 1.d0/6.d0
    real(kind=8) :: h2, resh2, resh2max
    real(kind=8),pointer,dimension(:,:,:) :: u, rhs
    integer :: i, j, k, ipass, is, ie, js, je, ks, ke
    h2 = mg_get_h( mglev ) ** 2
    call mg_arrp(mglev, ju, u)
    call mg_arrp(mglev, jrhs, rhs)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    resh2max = 0.d0
    do ipass=1,2 
       call mg_boundary_u(mglev,ju)
       do k = ks, ke
          do j = js, je
             do i = is + mod(j+k+ipass,2), ie, 2
                resh2 = &
                      u(i+1,j,k) + u(i-1,j,k) &
                     +u(i,j+1,k) + u(i,j-1,k) &
                     +u(i,j,k+1) + u(i,j,k-1) &
                     - 6.d0*u(i,j,k) &
                     - h2*rhs(i,j,k)
                u(i,j,k)=u(i,j,k) + sixth * resh2
                resh2max = max(resh2max, abs(resh2))
             enddo
          enddo
       enddo
    enddo
    Resmaxg = resh2max/h2 
  end subroutine mg_poisson_relax
  subroutine mg_cycle(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: ju, jrho
    real(kind=8),parameter :: errormax = 1.e-2
    integer,parameter :: nmg = 10
    real(kind=8) :: rhomax
    integer :: n, mglev
    mglev = MG_LevelMin
    if (ju /= IPSI) then
       call mg_alloc_arr(mglev, IPSI)
       call mg_copy(mglev, IPSI, ju) 
    end if
    if (jrho /= ISRC) then
       call mg_alloc_arr(mglev, ISRC)
       call mg_copy(mglev, ISRC, jrho) 
    end if
    call mg_boundary_fill0(mglev, IPSI) 
    call mg_boundary_u(mglev, IPSI) 
    call mg_alloc_arr(mglev, IRHO)
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
    elseif ( &
         FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
       call mg_alloc_arr(mglev, IETA)
       call mg_boundary_u(mglev, IETA)
    else
       print *, '*** mg_cycle: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    do n = 1, nmg
       call mg_resid(mglev, IRHO, IPSI, ISRC)
       if (n == 1) rhomax = mg_get_absmax(mglev, IRHO)
       call mg_lin 
       call mg_add(mglev, IPSI, IU)
       if ( Resmaxg / rhomax < errormax ) exit
    end do
    if (ju /= IPSI) call mg_copy(mglev, ju, IPSI) 
  end subroutine mg_cycle
  subroutine mg_lin
    use string
    use io_util
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    INTEGER :: Ncycle, Npre, Npost
    integer :: mglev, jcycle, jpre, jpost, vlevel
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
       Ncycle=2
       Npre=2
       Npost=2
    else
       Ncycle=2
       Npre=2
       Npost=2
    endif
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
    elseif ( &
         FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or.&
         FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    else
       print *, '*** mg_lin: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    do mglev = MG_LevelMin+1, MG_LevelMax
       call mg_alloc_arr(mglev, IRHO)
       call mg_rstrct(mglev, IRHO, IRHO)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       else
          print *, '*** mg_lin: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
    mglev = MG_LevelMax
    call mg_alloc_arr(mglev, IU)
    call mg_alloc_arr(mglev, IRHS)
    call mg_fill0(mglev, IU) 
    call mg_slvsml(mglev, IU, IRHO)
    do mglev = MG_LevelMax-1, MG_LevelMin, -1
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRES)
       call mg_interp(mglev, IU, IU, cubic=.TRUE.)
       call mg_copy(mglev, IRHS, IRHO)
       do jcycle = 1, Ncycle 
          do vlevel = mglev, MG_LevelMax-1 
             do jpre=1,Npre
                call mg_relax(vlevel, IU, IRHS)
             enddo
             call mg_resid(vlevel, IRES, IU, IRHS)
             call mg_rstrct(vlevel+1, IRHS, IRES)
             call mg_fill0(vlevel+1, IU)
          enddo
          vlevel = MG_LevelMax
          call mg_slvsml(vlevel, IU, IRHS)
          do vlevel = MG_LevelMax-1, mglev, -1 
             call mg_addint(vlevel, IU, IU, IRES)
             do jpost = 1, Npost
                call mg_relax(vlevel, IU, IRHS)
             enddo
          enddo
       enddo
    enddo
  end subroutine mg_lin
  subroutine mg_nonlin(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    INTEGER,parameter :: Ncycle=2, Npre=2, Npost=2
    integer :: mglev, jcycle, jpre, jpost, vlevel, icode
    real(kind=8),parameter :: ALPHA = 1.d-20
    real(kind=8) :: trerr
    mglev = MG_LevelMin
    if (ju /= IU) then
       call mg_alloc_arr(mglev, IU)
       call mg_copy(mglev, IU, ju)
    end if
    if (jrho /= IRHO) then
       call mg_alloc_arr(mglev, IRHO)
       call mg_copy(mglev, IRHO, jrho)
    end if
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call mg_boundary_u(MG_LevelMin, IDOD)
       call mg_boundary_u(MG_LevelMin, IDHE)
       call mg_boundary_u(MG_LevelMin, IDAD)
    else
       print *, '*** mg_nonlin: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    do mglev = MG_LevelMin+1, MG_LevelMax
       call mg_alloc_arr(mglev, IRHO)
       call mg_rstrct(mglev, IRHO, IRHO)
       call mg_alloc_arr(mglev, IU) 
       call mg_rstrct(mglev, IU, IU)
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
          do icode = IDOD, IDAD
             call mg_alloc_arr(mglev, icode)
             call mg_rstrct(mglev, icode, icode)
             call mg_boundary_u(mglev, icode)
          end do
       else
          print *, '*** mg_nonlin: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
    mglev = MG_LevelMax
    call mg_slvsml(mglev, IU, IRHO)
    call mg_alloc_arr(mglev, IRHS)
    call mg_alloc_arr(mglev, IRES)
    do mglev = MG_LevelMax-1, MG_LevelMin, -1
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRES)
       call mg_alloc_arr(mglev+1, IRUF)
       call mg_interp(mglev, IU, IU, cubic=.TRUE.)
       call mg_copy(mglev, IRHS, IRHO)
       do jcycle = 1, Ncycle 
          do vlevel = mglev, MG_LevelMax-1 
             do jpre=1,Npre
                call mg_relax(vlevel, IU, IRHS)
             enddo
             call mg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) 
             call mg_copy(vlevel+1, IRUF, IU) 
             call mg_rstrct(vlevel+1, IRHS, IRHS) 
             call mg_tau(vlevel+1, IRES, IU) 
             call mg_add(vlevel+1, IRHS, IRES) 
             if (vlevel == mglev) trerr = ALPHA*mg_get_absmax(vlevel+1, IRES)
          enddo
          vlevel = MG_LevelMax
          call mg_slvsml(vlevel, IU, IRHS)
          do vlevel = MG_LevelMax-1, mglev, -1 
             call mg_sub(vlevel+1, IRES, IU, IRUF) 
             call mg_interp(vlevel, IRES, IRES, cubic=.TRUE.) 
             call mg_add(vlevel, IU, IRES) 
             do jpost = 1, Npost
                call mg_relax(vlevel, IU, IRHS)
             enddo
          enddo
          if (Resmaxg < trerr) exit
          if (trerr <= tiny(trerr)) exit
       enddo
    enddo
    call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
         // ' at level = 0' // ', jcycle = '//trim(num2char(jcycle)) // ' (final)')
    if (ju /= IU) call mg_copy(MG_LevelMin, ju, IU) 
  end subroutine mg_nonlin
  subroutine mg_nonlin_vcycle(ju, jrho)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    use string
    use io_util
    integer,intent(IN) :: ju, jrho
    INTEGER,parameter :: Ncycle=1, Npre=2, Npost=2
    integer :: mglev, jcycle, jpre, jpost, vlevel, icode
    real(kind=8),parameter :: ALPHA = 1.d-20
    real(kind=8) :: trerr
    mglev = MG_LevelMin
    if (ju /= IU) then
       call mg_alloc_arr(mglev, IU)
       call mg_copy(mglev, IU, ju)
    end if
    if (jrho /= IRHO) then
       call mg_alloc_arr(mglev, IRHO)
       call mg_copy(mglev, IRHO, jrho)
    end if
    if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION ) then
    elseif ( FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION )then
       call mg_boundary_u(MG_LevelMin, IETA)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION ) then
       call mg_boundary_u(MG_LevelMin, IDOD)
       call mg_boundary_u(MG_LevelMin, IDHE)
       call mg_boundary_u(MG_LevelMin, IDAD)
    else
       print *, '*** mg_nonlin_vcycle: this type of PDE is not supported', FMG_PDE_TYPE
    endif
    do mglev = MG_LevelMin+1, MG_LevelMax
       if ( FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       elseif ( &
            FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION .or. &
            FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION ) then
          call mg_alloc_arr(mglev, IETA)
          call mg_rstrct(mglev, IETA, IETA)
          call mg_boundary_u(mglev, IETA)
       elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
          do icode = IDOD, IDAD
             call mg_alloc_arr(mglev, icode)
             call mg_rstrct(mglev, icode, icode)
             call mg_boundary_u(mglev, icode)
          end do
       else
          print *, '*** mg_nonlin_vcycle: this type of PDE is not supported', FMG_PDE_TYPE
       end if
    enddo
    do mglev = MG_LevelMin, MG_LevelMax
       call mg_alloc_arr(mglev, IU)
       call mg_alloc_arr(mglev, IRHS)
       call mg_alloc_arr(mglev, IRUF)
       call mg_alloc_arr(mglev, IRES)
    end do
    call mg_copy(MG_LevelMin, IRHS, IRHO) 
    do jcycle = 1, Ncycle 
       do vlevel = MG_LevelMin, MG_LevelMax-1 
          do jpre=1,Npre
             call mg_relax(vlevel, IU, IRHS)
          enddo
          call mg_rstrct(vlevel+1, IU, IU, cubic=.TRUE.) 
          call mg_copy(vlevel+1, IRUF, IU) 
          call mg_rstrct(vlevel+1, IRHS, IRHS) 
          call mg_tau(vlevel+1, IRES, IU) 
          call mg_add(vlevel+1, IRHS, IRES) 
          if (vlevel == MG_LevelMin) trerr = ALPHA*mg_get_absmax(vlevel+1, IRES)
       enddo
       vlevel = MG_LevelMax
       call mg_slvsml(vlevel, IU, IRHS)
       do vlevel = MG_LevelMax-1, MG_LevelMin, -1 
          call mg_sub(vlevel+1, IRES, IU, IRUF) 
          call mg_interp(vlevel, IRES, IRES, cubic=.TRUE.) 
          call mg_add(vlevel, IU, IRES) 
          do jpost = 1, Npost
             call mg_relax(vlevel, IU, IRHS)
          enddo
       enddo
       if (Resmaxg < trerr) exit
       if (trerr <= tiny(trerr)) exit
    enddo
    call print_msg( 'error in MG = ' // trim(num2char(Resmaxg)) // ' / ' // trim(num2char(trerr)) &
         // ' at level = 0' // ', jcycle = '//trim(num2char(jcycle)) // ' (final)')
    if (ju /= IU) call mg_copy(MG_LevelMin, ju, IU) 
  end subroutine mg_nonlin_vcycle
  subroutine mg_init(fmglev)
    integer,intent(IN) :: fmglev
    call mg_data_init(fmglev)
  end subroutine mg_init
  subroutine mg_finalize
    use mg_interpol_cubic
    call mg_data_finalize
    call mg_interp_cubic_finalize
  end subroutine mg_finalize
  subroutine mg_slvsml(mglev, ju, jrhs)
    integer,intent(IN) :: mglev, ju, jrhs
    integer :: n
    do n = 1, 2
       call mg_relax(mglev, ju, jrhs)
    end do
  end subroutine mg_slvsml
  subroutine mg_resid(mglev, jres, ju, jrhs)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, jres, ju, jrhs
    real(kind=8) :: h2i
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_poisson_resid(mglev, jres, ju, jrhs)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    else
       print *, '*** mg_resid: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine mg_resid
  subroutine mg_rstrct(mglev, iuc, iuf, cubic)
    use restriction
    use fmg_data, only : FMG_PDE_LINEAR
    integer,intent(IN) :: mglev, iuc, iuf
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    logical,intent(IN),optional :: cubic
    logical :: bool_cubic
    integer :: ics, jcs, kcs, ice, jce, kce, lf, lc
    bool_cubic = .FALSE. 
    if (present(cubic)) bool_cubic = cubic
    lf = mglev-1
    lc = mglev
    call mg_get_gridsize(lc, ics,jcs,kcs,ice,jce,kce)
    uf => mg_get_arrp(lf,iuf)
    uc => mg_get_arrp(lc,iuc)
    if (bool_cubic) call mg_boundary_u(lf, iuf)
    call rstrct(uc, uf, &
            ics, jcs, kcs, ice, jce, kce, &
            ics, jcs, kcs, cubic)
  end subroutine mg_rstrct
  subroutine mg_interp(mglev, iuf, iuc, cubic)
    use mg_interpol_cubic
    integer,intent(IN) :: mglev, iuf, iuc
    logical,intent(IN),optional :: cubic
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    logical :: bool_cubic
    integer :: lf, lc
    lf = mglev
    lc = mglev+1
    bool_cubic = .FALSE. 
    if (present(cubic)) bool_cubic = cubic
    call mg_boundary_extrap(lc, iuc)
    call mg_boundary_u(lc,iuc)
    if (bool_cubic) then
       call mg_interp_cubic(lf, iuf, iuc)
    else
       call mg_interp_linear(lf, iuf, iuc)
    end if
    call mg_boundary_fill0(lf, iuf)
    call mg_boundary_fill0(lc, iuc)
    call mg_boundary_u(lf,iuf)
    call mg_boundary_u(lc,iuc)
  end subroutine mg_interp
  subroutine mg_interp_linear(mglev, iuf, iuc)
    use interpolation
    integer,intent(IN) :: mglev, iuf, iuc
    real(kind=8),pointer,dimension(:,:,:,:) :: uf, uc
    integer :: lf, lc, &
         if, jf, kf, &
         ic, jc, kc, &
         ic1, jc1, kc1, &
         ifs, jfs, kfs, ife, jfe, kfe, m
    lf = mglev
    lc = mglev+1
    call mg_get_gridsize(lf, ifs,jfs,kfs,ife,jfe,kfe)
    uf => mg_get_arrp(lf,iuf)
    uc => mg_get_arrp(lc,iuc)
    call interp_trilinear(uf, uc, ifs, jfs, kfs, ife, jfe, kfe, ifs, jfs, kfs)
  end subroutine mg_interp_linear
  subroutine mg_fill0(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    real(kind=8),parameter :: zero = 0.d0
    u => mg_get_arrp(mglev, icode)
    u = zero
  end subroutine mg_fill0
  subroutine mg_relax(mglev, ju, jrhs)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, ju, jrhs
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_POISSON_EQUATION) then
       call mg_poisson_relax(mglev, ju, jrhs)
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_DIFFUSION_EQUATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** mg_relax: this type of PDE is not supported', FMG_PDE_TYPE
    endif
  end subroutine mg_relax
  subroutine mg_copy(mglev, iout, iin)
    integer,intent(IN) :: mglev, iout, iin
    real(kind=8),pointer,dimension(:,:,:,:) :: ain, aout
    aout => mg_get_arrp(mglev,iout)
    ain => mg_get_arrp(mglev,iin)
    aout = ain
  end subroutine mg_copy
  subroutine mg_sub(mglev, iout, ia, ib)
    integer,intent(IN) :: mglev, iout, ia, ib
    integer :: amrlev, gid
    real(kind=8),pointer,dimension(:,:,:,:) :: out, a, b
    out => mg_get_arrp(mglev,iout)
    a => mg_get_arrp(mglev,ia)
    b => mg_get_arrp(mglev,ib)
    out = a - b
  end subroutine mg_sub
  subroutine mg_add(mglev, ia, ida)
    integer,intent(IN) :: mglev, ia, ida
    real(kind=8),pointer,dimension(:,:,:,:) :: a, da
    a => mg_get_arrp(mglev,ia)
    da => mg_get_arrp(mglev,ida)
    a = a + da
  end subroutine mg_add
  subroutine mg_addint(mglev, juf, juc, jres)
    integer,intent(IN) :: mglev, juf, juc, jres
    call mg_interp(mglev, jres, juc)
    call mg_add(mglev, juf, jres)
  end subroutine mg_addint
  subroutine mg_boundary_fill0(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    real(kind=8) :: zero = 0.d0
    integer :: is, ie, js, je, ks, ke
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    u => mg_get_arrp(mglev,icode)
    if (TouchBoundary(Left, 0)) u(lbound(u,1):is-1,:,:,:) = zero
    if (TouchBoundary(Right,0)) u(ie+1:ubound(u,1),:,:,:) = zero
    if (TouchBoundary(Left, 1)) u(:,lbound(u,2):js-1,:,:) = zero
    if (TouchBoundary(Right,1)) u(:,je+1:ubound(u,2),:,:) = zero
    if (TouchBoundary(Left, 2)) u(:,:,lbound(u,3):ks-1,:) = zero
    if (TouchBoundary(Right,2)) u(:,:,ke+1:ubound(u,3),:) = zero
  end subroutine mg_boundary_fill0
  subroutine mg_boundary_extrap(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    real(kind=8) :: zero = 0.d0
    integer :: is, ie, js, je, ks, ke
    integer :: isg, ieg, jsg, jeg, ksg, keg
    call mg_get_gridsize (mglev, is,js,ks,ie,je,ke)
    call mg_get_gridsizeGh(mglev, isg,jsg,ksg,ieg,jeg,keg)
    u => mg_get_arrp(mglev,icode)
    if (TouchBoundary(Left, 0)) u(isg:is-1,:,:,:) = -u(is:is-1+Ngh:-1,:,:,:)
    if (TouchBoundary(Right,0)) u(ie+1:ieg,:,:,:) = -u(ie+1-Ngh:ie:-1,:,:,:)
    if (TouchBoundary(Left, 1)) u(:,jsg:js-1,:,:) = -u(:,js:js-1+Ngh:-1,:,:)
    if (TouchBoundary(Right,1)) u(:,je+1:jeg,:,:) = -u(:,je+1-Ngh:je:-1,:,:)
    if (TouchBoundary(Left, 2)) u(:,:,ksg:ks-1,:) = -u(:,:,ks:ks-1+Ngh:-1,:)
    if (TouchBoundary(Right,2)) u(:,:,ke+1:keg,:) = -u(:,:,ke+1-Ngh:ke:-1,:)
  end subroutine mg_boundary_extrap
  function mg_get_absmax(mglev, icode) result(umax)
    integer,intent(IN) :: mglev, icode
    real(kind=8) :: umax
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    umax = maxval(abs(u(is:ie,js:je,ks:ke,:)))
  end function mg_get_absmax
  subroutine mg_minmax(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    print *, 'minmax', minval(u(is:ie,js:je,ks:ke,:)), maxval(u(is:ie,js:je,ks:ke,:))
  end subroutine mg_minmax
  subroutine mg_minmaxall(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    u => mg_get_arrp(mglev,icode)
    print *, 'minmaxall', minval(u), maxval(u)
  end subroutine mg_minmaxall
  subroutine mg_minmaxbnd(mglev, icode)
    integer,intent(IN) :: mglev, icode
    real(kind=8),pointer,dimension(:,:,:,:) :: u
    integer :: is, ie, js, je, ks, ke
    real(kind=8) :: umin, umax
    u => mg_get_arrp(mglev,icode)
    call mg_get_gridsize(mglev, is,js,ks,ie,je,ke)
    umin = minval( u(lbound(u,1):is-1,:,:,:) )
    umin = min(umin, minval(u(ie+1:ubound(u,1),:,:,:)))
    umin = min(umin, minval(u(:,lbound(u,2):js-1,:,:)))
    umin = min(umin, minval(u(:,je+1:ubound(u,2),:,:)))
    umin = min(umin, minval(u(:,:,lbound(u,3):ks-1,:)))
    umin = min(umin, minval(u(:,:,ke+1:ubound(u,3),:)))
    umax = maxval( u(lbound(u,1):is-1,:,:,:) )
    umax = max(umax, maxval(u(ie+1:ubound(u,1),:,:,:)))
    umax = max(umax, maxval(u(:,lbound(u,2):js-1,:,:)))
    umax = max(umax, maxval(u(:,je+1:ubound(u,2),:,:)))
    umax = max(umax, maxval(u(:,:,lbound(u,3):ks-1,:)))
    umax = max(umax, maxval(u(:,:,ke+1:ubound(u,3),:)))
    print *, 'minmaxbnd', umin, umax
  end subroutine mg_minmaxbnd
  subroutine mg_tau(mglev, jtau, ju)
    use fmg_data, only : &
         FMG_PDE_TYPE, &
         FMG_PDE_TYPE_POISSON_EQUATION, &
         FMG_PDE_TYPE_DIFFUSION_EQUATION, &
         FMG_PDE_TYPE_OHMIC_DISSIPATION, &
         FMG_PDE_TYPE_AMBIP_DIFFUSION
    integer,intent(IN) :: mglev, jtau, ju
    if (FMG_PDE_TYPE == FMG_PDE_TYPE_OHMIC_DISSIPATION) then
    elseif (FMG_PDE_TYPE == FMG_PDE_TYPE_AMBIP_DIFFUSION) then
    else
       print *, '*** mg_tau: this type of PDE is not supported'
    endif
  end subroutine mg_tau
end module mg
