! configuration file for AMR

! model of adaptive ray tracing (ART)
! 0 -> no art
! 1 -> EUV + thin FUV
! 2 -> EUV + FUV w/ self-shielding considered
#define MODEL_ART 2

! use numerical flux RoeM2
#define FLUX_ROEM2

#define NPE 400

! gridsize for each grid
#define NI 8
#define NJ 8
#define NK 8

! number of grid in each direction for base grid.
! number of cell is NI*NGI_BASE x NJ*NGJ_BASE x NK*NGK_BASE
#define NGI_BASE 8
#define NGJ_BASE 8
#define NGK_BASE 8

! maximun level
#define NL 16

! maximum grid id
!#define NGID 1500
#define NGID 5000

! skip HD update
!#define SKIP_HD_UPDATE

! isothermal-gas calculation
!#define ISOTHERMAL

! create radiation source by hand (for test)
!#define CREATE_RADSOURCE_BY_HAND

! random rotaion for ART is not performed
!#define ART_NO_RANDOM_ROTAION

! no radiation force (radiation force has been implemented but not been tested yet)
!#define NO_RADFORCE
!#define NO_DUST_RADFORCE
!#define NO_IONIZATION


!YES OR NO
#define YES  0
#define NO   1


! metal
#define METAL
!#define TEST_STELLAR_100


! Label component of hydro
! Self-gravitational adabatic HD (no chemistry)
#if MODEL_ART == 0
#define MRHO 0
#define MVX  1
#define MVY  2
#define MVZ  3
#define MP   4
#define MPSI 5
#define MGX  6
#define MGY  7
#define MGZ  8
#define NM   9

! Self-gravitational adabatic HD (H,H2,El,Hp,Hm,H2p,EUV-RT,FUV-thin)
#elif MODEL_ART == 1
#define MRHO 0
#define MVX  1
#define MVY  2
#define MVZ  3
#define MP   4
#define MHN  5
#define MH2  6
#define MEL  7
#define MHP  8
#define MHM  9
#define MH2P 10
#define MKPI 11
#define MHPI 12
#define MXPI 13
#define MYPI 14
#define MZPI 15
#define MPSI 16
#define MGX  17
#define MGY  18
#define MGZ  19
#ifdef METAL
#define MTD  20
#endif
#ifdef METAL
#define NM   21
#else
#define NM   20
#endif


! Self-gravitational adabatic HD (H,H2,El,Hp,Hm,H2p,EUV-RT,FUV-RT)
#elif MODEL_ART == 2

#include "chemistry_radtr.h"

#endif


! Self-gravitational barotropic MHD
! #define MRHO 0
! #define MVX  1
! #define MVY  2
! #define MVZ  3
! #define MBX  4
! #define MBY  5
! #define MBZ  6
! #define MDB  7
! #define MPSI 8
! #define MGX  9
! #define MGY  10
! #define MGZ  11
! #define NM   12

! label of coordinates
#define MX 0
#define MY 1
#define MZ 2

! number of ghost cell
#define N_GHOST_CELL 2

! CFL number
#define CFL 0.7

! Spatial order of accuracy
! #define RECONSTRUCTION_MUSCL3

! sink cell
#define SINKPARTICLE

! Precision
#define DBL_KIND 8
#define LLONG_KIND 8

! Primary rank
#define PRIMARY_RANK 0

! use self-gravity (multigrid)
#ifdef MPSI
!#define WITH_SELFGRAVITY
#define FMG_POISSON
#endif
#define SINGLE_STEP

! use single timestep rather than multi-timestep of Berger and Colella (1998)
! 最終的に重力を考慮するので、テスト計算でもSINGLE_STEPで計算する
!#ifdef WITH_SELFGRAVITY
!#define SINGLE_STEP
!#endif

! use external force
! #define EXTERNALFORCE

!#define REFINE_NO_MARGIN  

!Turbulence
#define TURBINI
#ifdef TURBINI

! kind of complex
#define CMPLX_KIND 8

! seed of random
#define SEED 3

! filename
#define FOUT 'TurbulentVelocity'

#endif !TURBINI


! MHD
#if defined(MBX) || defined(MBY) || defined(MBZ) || defined(MDB)
#define MHD
! #define MHD_ROE
#endif

! Define this argument when periodic boundary condition
#define ZERO_AVERAGE_DENSITY

! FMG_LAMDA = 1 for complete implicit. FMG_LAMDA = 1/2 for Crank-Nicolson
! #define FMG_OHMIC_DISSIPATION
! #define FMG_LAMBDA 1.d0

! fixed grid
!#define NOT_REFINEMENT

! used in refineCond
! #define BENCH_MARK

! It is usefull for estimate of performance
! #define FORBIT_WRITEDATA

! Emulation of two-dimension
! #define EMULATE_2DIM

!-----------------------------------------------------------------------
! macro for programing
!-----------------------------------------------------------------------
! array size
#define ARRAYSIZE5(A) lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2),lbound(A,3):ubound(A,3),lbound(A,4):ubound(A,4),lbound(A,5):ubound(A,5)
#define ARRAYSIZE4(A) lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2),lbound(A,3):ubound(A,3),lbound(A,4):ubound(A,4)
#define ARRAYSIZE3(A) lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2),lbound(A,3):ubound(A,3)
#define ARRAYSIZE2(A) lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2)
#define ARRAYSIZE1(A) lbound(A,1):ubound(A,1)

#define ARRAYSIZE4_4(A) lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2),lbound(A,3):ubound(A,3),0:4

#define ARRAYSIZE_IJK Imin:Imax,Jmin:Jmax,Kmin:Kmax
#define ARRAYSIZE_IJKM Imin:Imax,Jmin:Jmax,Kmin:Kmax,Mmin:Mmax
#define ARRAYSIZE_IJKGH Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh
#define ARRAYSIZE_IJKMGH Imingh:Imaxgh,Jmingh:Jmaxgh,Kmingh:Kmaxgh,Mmin:Mmax


#define PTF5(A) lbound(A,1),lbound(A,2),lbound(A,3),lbound(A,4),lbound(A,5)
#define PTF4(A) lbound(A,1),lbound(A,2),lbound(A,3),lbound(A,4)
#define PTF3(A) lbound(A,1),lbound(A,2),lbound(A,3)
#define PTF2(A) lbound(A,1),lbound(A,2)
#define PTF1(A) lbound(A,1)

!-----------------------------------------------------------------------
! for debug
!-----------------------------------------------------------------------
#define PRINTV(f) write(*,*) "** " ,  #f , " =" , f
#define PRINTRV(f) write(*,*) "** (", myrank, ") " ,  #f , " =" , f
#define HALT   write(*,*) '*** OK'; call mpi_barrier(MPI_COMM_WORLD, ierr); call mpi_finalize(ierr); stop;
#define POINTOK     print *, '*** OK'; call mpi_barrier(MPI_COMM_WORLD, ierr)


!-----------------------------------------------------------------------
! to access coarser subscripts
!-----------------------------------------------------------------------
! floor(real(i) / 2), max(i/2, abs(i-1)/2)*sign(1,i), i/2 + mod(min(i,0),2)
#define IJKC(ijk,ijk0)  ((ijk)-(ijk0))/2 + mod(min((ijk)-(ijk0),0),2) + (ijk0)
!#define IJKC(ijk,ijk0)  floor(real((ijk)-(ijk0)) / 2)+ (ijk0)
#define IJKF(ijk,ijk0)  ((ijk) - (ijk0)) * 2 + (ijk0)

!-----------------------------------------------------------------------
! number of dimension
!-----------------------------------------------------------------------
#define NDIM(array) ubound(shape(array),1)

!-----------------------------------------------------------------------
! vector loop (intel compiler)
!-----------------------------------------------------------------------
#define VECT  DEC$ VECTOR ALWAYS
