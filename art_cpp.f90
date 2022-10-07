module art
  use overBlockCoordinates
  use grid
  use unit
  use primordial
  use pix_tools
  use mpilib
  use radiationSource
  use sinkParticle, only : sp_getSinkRadius
  use kinzoku
  implicit none
  private
  real(kind=8),parameter :: PI=atan(1d0)*4d0
  integer,parameter :: MTH = 0, MPH = 1 
  integer,parameter :: NIL = 0
  integer,parameter :: NIR = 1
  integer,parameter :: NJL = 2
  integer,parameter :: NJR = 3
  integer,parameter :: NKL = 4
  integer,parameter :: NKR = 5
  type t_dirHpix
     integer :: res 
     integer(kind=8) :: ipix 
  end type t_dirHpix
  type t_ray
     integer :: isrc 
     type(t_dirHpix) :: Hpix 
     real(kind=8) :: dist 
     real(kind=8),dimension(0:5 -1) :: tau 
     integer :: gid 
     type(t_ray),pointer :: next => null() 
 end type t_ray
  type t_pray
     type(t_ray),pointer :: ray_list => null()
  end type t_pray
  type(t_ray),save,pointer :: my_ray_list => null() 
  type(t_pray),save,dimension(0:400 -1) :: send_ray_list 
  integer,save,dimension(0:400 -1) :: N_send 
  type(t_rs_info),save,pointer :: rs_info
  integer,parameter :: RAY_BYTE = (4*3+8*(2+5)) 
  integer,parameter :: NSEND_BYTE = 4 
  integer,parameter :: NCOMP_BYTE = 8 
  integer,parameter :: N_SEND_MAX = 1024 
  integer,parameter :: BUFSIZE_RAYS = (N_SEND_MAX*RAY_BYTE+NSEND_BYTE+3)/4 
  integer,parameter :: BUFSIZE_NCOMP = (NCOMP_BYTE+3)/4 
  integer,save,dimension(0:BUFSIZE_RAYS-1) :: bufr_rays 
  integer(kind=8),save :: bufr_Ncomp 
  integer,parameter :: MPI_TAG_MAX_MACHINE = 2097151 
  type t_msg_rays
     integer :: req 
     integer :: tag 
     integer,pointer,dimension(:) :: bufs 
     type(t_msg_rays),pointer :: next => null() 
  end type t_msg_rays
  type t_msg_Ncomp
     integer :: req 
     integer :: tag 
     integer(kind=8) :: bufs 
     type(t_msg_Ncomp),pointer :: next => null() 
  end type t_msg_Ncomp
  type(t_msg_rays),save,pointer :: msg_rays_list => null() 
  type(t_msg_Ncomp),save,pointer :: msg_Ncomp_list => null() 
  integer,parameter :: HEALPIX_MIN_LEVEL = 3 
  integer,parameter :: HEALPIX_MAX_LEVEL = 12 
  integer,parameter :: MIN_RAY_PER_CELL = 3
  integer(kind=8) :: N_comp_total 
  integer(kind=8),save,dimension(0:400 -1) :: N_comp_local 
  integer(kind=8) :: N_comp_global 
  logical :: N_comp_local_update 
  real(kind=8),parameter :: TAU_ENDRAY = 100.d0 
  real(kind=8),dimension(0:2,0:2) :: m_rot
  real(kind=8) :: SinkRadius
  integer,parameter :: N_SEND_GO = 200 
  integer,parameter :: N_GRIDCOUNT_GO = 2000 
  integer,save :: tag_max_ray, tag_max_Ncomp
  integer(kind=8) :: loop_count, grid_count, grid_count_tot 
  integer::t_ini,t_cur,t_rat 
  logical :: glob_debug_flag_art = .FALSE. 
  integer :: N_myrays 
end module art
