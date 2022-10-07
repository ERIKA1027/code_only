! set chemistry
! variables
 
! radiation source -------------
!#define RADSOURCE_SC
!#define STOCHASTIC_STELLAR_MODEL
#define FOL_RADSOURCE "Z1"
#define MINMASS_RS 0.1d0
! ------------------------------
#define RADTR_M1_ONLY
#define RADTR_M1closer
#define CFLfac_radtr 0.8
#define MAX_STEP_RADTR_TO_HYDRO 10.0
!#define SOLVER_RADTR_HLL
 
#define M1CLOSER_EUV_TRANSFER
#define RECOM_ION_PHOTO
#define M1CLOSER_FUV_TRANSFER
#define FUV_COLUMNDENS_OPTION 2
#define SELF_SHIELDING_H2_FORMULA 2
#define M1CLOSER_IR_TRANSFER
#define NUM_COMP_TRANSFER 3
 
! dust evolution
!#define DUST_NOTCONSTANT
 
! Radiation pressure
#define EXTERNALFORCE
 
! photoelectric heating
#define PHOTOELECTRIC_HEATING YES
 

! cosmic ray
#define INCLUDE_COSMICRAY
#define COSMICRAY_RATE_ION 2.d-16



#define MRHO 0
#define MVX 1
#define MVY 2
#define MVZ 3
#define MP 4
#define MHN 5
#define MH2 6
#define MEL 7
#define MHP 8
#define MHM 9
#define MH2P 10
#define MCO 11
#define MKPI 12
#define MHPI 13
#define MXPI 14
#define MYPI 15
#define MZPI 16
#define MPSI 17
#define MGX 18
#define MGY 19
#define MGZ 20
#define MTD 21
#define MER 22
#define MFRX 23
#define MFRY 24
#define MFRZ 25
#define MEF 26
#define MFRFX 27
#define MFRFY 28
#define MFRFZ 29
#define MEIR 30
#define MFIRX 31
#define MFIRY 32
#define MFIRZ 33
#define NM 34
#define NCEHM_MIN 5
#define NCEHM_MAX 10
 
#define DEF_TCMB 2.73
#define H2PD_HEATING
#define CHEM_MODEL_HF2020
