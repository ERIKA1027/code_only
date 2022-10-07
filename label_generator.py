# -*- coding: utf-8- -*-
# hydro dynamics 用のlabel を生成する

import sys
import os
import shutil

ON = 1
OFF= 0

# ----------------------------------
# 1: DIRECT_LIGHT + IR_TRANSFER
#    DIRECT_LIGHT  = ON
#    M1_CLOSURE    = ON
#    EUV_TRANSFER  = OFF
#    IR_TRANSFER   = ON
# ----------------------------------


# ---------------------
#     setting !!
# --------------------

# ==============================
#             MHD
# ==============================
MHD           = OFF


# ==============================
#             OUTFLOW
# ==============================
OUTFLOW       = OFF
OUTFLOW_FW    = "0.3d0"     # fraction of outflow rate to inflow rate
OUTFLOW_FV    = "0.3d0"     # ratio of outflow velocity to kep velo at stellar surface


# ==============================
#       Radiation sources
# ==============================

STC_STELLEV   = OFF          # stochastic stellar model  
RADSOURCE_SC  = OFF         # radiation source ON: sink particles are treated as mini star cluster 
FILE_RADSOURCE= "Z1"        # folder name for radiation source 
MINMASS_RS    = ON          # minimal mass of radiation source
MINMASS       = "0.1d0"     # [Msun]
USE_TKH_TACC  = OFF         # using TKH & TACC for calculating luminosity
REDUCE_SMASS  = OFF         # reducing stellar mass at radiation_source
REDUCE_SMRATE = "0.3d0"     # reduce rate

# ==============================
#    OPTION FOR SC formation
# ==============================
SKIP_RADTR_BEFORE_SF = OFF  # Skip radtr before creation of radiation sorces



# ==============================
#           RT Solver
# ==============================

DIRECT_LIGHT  = OFF
M1_CLOSURE    = ON


# ==============================
#           M1 closer 
# ==============================

CFL_CONDITION       = 0.8   # CFL condition
MAX_STEP_RAD        = 10.0  # maximum steps in one step of hydrodynamics
SOLVER              = "GLF" # GLF or HLL

EUV_TRANSFER        = ON    # EUV photons
FUV_TRANSFER        = ON    # H2 photodissociation
IR_TRANSFER         = ON    # IR components

RADIATION_PRESSURE = ON     # Radiation pressure caused by EUV, EUV_TRANSFER == ONの場合にしか有効にならない
DUST_ATTENUATION   = ON     # Dust attenuation をON

EXTERNAL_FUV_RAD    = OFF
EXTERNAL_G0         = 1.0   # background UV strength Habing flux G0 = 1.6e-3 erg cm^-2 s^-1

IGNORING_SINK_TRANS_M1 = OFF # ignore radiative processes in sink cell

# ==============================
#    OPTION FOR EUV TRANS
# ==============================
NO_IONIZATION       = OFF   # ignoring photoionization of HI
RECOM_ION_PHOTO     = ON    # Recombination photons 

# ==============================
#    OPTION FOR FUV TRANS
# ==============================
SEPARATE_FUV_TRANS      = OFF       # seperate components of FUV for H2, CO, DUST
PHOTOELECTRIC_HEATING   = ON        # inlcude photoelctric heating
COLUMN_DENSITY_OPTION   = "Sobolev" # How to calculate the column density 
                                    # Jeans:  jeans length
                                    # Sobolev:Sobolev length 
                                    # 今のところ, SEPARATE_FUV_TRANS = OFF の場合しかSobolevに対応していない
                                    # SEPARATE_FUV_TRANS = ON でも使える様にするには, radtr_chem.F90のchemistry計算後の
                                    # self-shieldingを考慮するところを変更しないといけない

#SELF_SHIELDING_H2       = "DB96_SIMPLE" # simple formula of Draine Bertoldi 1996
SELF_SHIELDING_H2       = "WLCOTGREEN19" # modified formula of Wolcott-Green 2019
                                         # # 今のところ, SEPARATE_FUV_TRANS = OFF の場合しか対応していない

# ==============================
#       DIRECT LIGHT 
# ==============================
SUB_GRID_MODEL_DIRECTLIGHT = OFF # direct light のmaskingをするかしないか



# ==============================
#           CHEMISTRY
# ==============================
DEF_TCMB = 2.73         # tempearture of background radiation 

DUST_NOTCONSTANT = OFF          # calculation of dust grain abundance
H2PD_HEATING     = ON           # photodissociation heating
METAL_TRANSFER   = OFF          # metal transport 
COSMIC_RAY       = ON           # include cosmic ray
RATE_COSMICRAY   = "2.d-16"     # ionization of primary cosmic ray [s^-1]

# option for shelf shielding

# =============================
#       sink cells
# ============================
OUTPUTSINK_CRTEPOCH    = OFF      # option for data output at the epoch of sink creation 
OUTPUTSINK_CRT_BOXSIZE = "1.d0"  # boxsize [pc]


# ==============================
#       OTHER OPTIONS
# ==============================
DM_NFW_PROFILE   = OFF
DM_POTENTIAL     = OFF  # include DM matter potential in calculations




# end
# ----------------------------------------------------------------------------------------------------------------------------------

if M1_CLOSURE == OFF:
    EUV_TRANSFER        = OFF    # EUV photons
    FUV_TRANSFER        = OFF    # H2 photodissociation
    IR_TRANSFER         = OFF    # IR components
    RECOM_ION_PHOTO     = OFF    # Recombination photons 

if FUV_TRANSFER == OFF:
    SEPARATE_FUV_TRANS == OFF 
    PHOTOELECTRIC_HEATING == OFF

if NO_IONIZATION == ON:
    RECOM_ION_PHOTO = OFF    

# -----------------------------------
filename = "label_chem2hydro.txt"
with open(filename) as f:
    s= f.read()
    labels_chem = s.split("\n")
labels_chem.remove("")
# ----------------------------------

HF2020_switch = OFF
if any("HF2020" in lab for lab in labels_chem):
    labels_chem.remove("HF2020")
    HF2020_switch = ON

# labels 
labels = ["MRHO"]    
labels += ["MVX"]     
labels += ["MVY"]     
labels += ["MVZ"]     
labels += ["MP"]      
for lab in labels_chem:
    labels += [lab]
if HF2020_switch == ON:
    labels += ["MCO"]
labels += ["MKPI"]    
labels += ["MHPI"]    
labels += ["MXPI"]    
labels += ["MYPI"]    
labels += ["MZPI"]    
labels += ["MPSI"]    
labels += ["MGX"]     
labels += ["MGY"]     
labels += ["MGZ"]     

if MHD == ON:
    labels += ["MBX"]
    labels += ["MBY"]
    labels += ["MBZ"]
    labels += ["MDB"]



if DIRECT_LIGHT == ON:
    labels += ["MKPD"]    
    labels += ["MDPH"]    
    labels += ["MGFUV"]   
    labels += ["MDPCO"]  
    labels += ["MKPOII"]
labels += ["MTD"]     
if DUST_NOTCONSTANT == ON:
    labels += ["MDRHO"]


if EUV_TRANSFER  == ON:
    labels += ["MER"]     
    labels += ["MFRX"]    
    labels += ["MFRY"]    
    labels += ["MFRZ"]    
if FUV_TRANSFER == ON:
    labels += ["MEF"]
    labels += ["MFRFX"]
    labels += ["MFRFY"]
    labels += ["MFRFZ"]

    if SEPARATE_FUV_TRANS == ON:
        labels += ["MECOF"]
        labels += ["MFCORFX"]
        labels += ["MFCORFY"]
        labels += ["MFCORFZ"]
        labels += ["MEDUSTF"]
        labels += ["MFDUSTRFX"]
        labels += ["MFDUSTRFY"]
        labels += ["MFDUSTRFZ"]


if IR_TRANSFER  == ON:
    labels += ["MEIR"]     
    labels += ["MFIRX"]    
    labels += ["MFIRY"]    
    labels += ["MFIRZ"]    

if METAL_TRANSFER == ON:
    labels += ["MMET"]
if DM_POTENTIAL == ON:
    labels += ["MDMRHO"]



filename = "./chemistry_radtr.h"
if(os.path.isfile(filename)):
    os.remove(filename)

lists = ["! set chemistry"]
lists +=["! variables"]
lists += [" "]

lists += ["! radiation source -------------"]
if STC_STELLEV == ON and RADSOURCE_SC == ON:
    print("you cannot define STC_STELLEV == ON and RADSOURCE_SC == ON simultaneously")
    sys.exit()

if RADSOURCE_SC == ON:
    lists += ["#define RADSOURCE_SC"]
else:
    lists += ["!#define RADSOURCE_SC"]
    
    if STC_STELLEV == ON:
        lists += ["#define STOCHASTIC_STELLAR_MODEL"]

        if FILE_RADSOURCE == "Z1":
            lists += ["#define STELLAR_MODEL_FILENAME \"tsevolv_0.0.dat\""]
        elif FILE_RADSOURCE == "Z-1":
            lists += ["#define STELLAR_MODEL_FILENAME \"tsevolv_-1.0.dat\""]
        elif FILE_RADSOURCE == "Z-2":
            lists += ["#define STELLAR_MODEL_FILENAME \"tsevolv_-2.0.dat\""]
        else:
            print ("not defined yet \n")
            sys.exit()

    else:
        lists += ["!#define STOCHASTIC_STELLAR_MODEL"]



if USE_TKH_TACC == ON:
    lists += ["#define USE_TKH_TACC"]

if SKIP_RADTR_BEFORE_SF == ON:
    lists += ["#define SKIP_RADTR_BEFORE_STARFORM"]


if(RADSOURCE_SC == ON):
    if FILE_RADSOURCE == "Z1":
        lists += ["#define FOL_RADSOURCE 0"]
    elif FILE_RADSOURCE == "Z-1":
        lists += ["#define FOL_RADSOURCE 1"]
    elif FILE_RADSOURCE == "Z-2":
        lists += ["#define FOL_RADSOURCE 2"]
    else:
        print ("not defined yet \n")
        sys.exit()
else:
    lists += ["#define FOL_RADSOURCE \"" + FILE_RADSOURCE+ "\""]

if MINMASS_RS == ON:
    lists += ["#define MINMASS_RS " +MINMASS ]


if REDUCE_SMASS == ON:        # reducing stellar mass at radiation_source
    lists += ["#define REDUCE_SMRATE " + REDUCE_SMRATE]


lists += ["! ------------------------------"]

if(os.path.isfile("./Configs/objects.defs")):
    os.remove("./Configs/objects.defs")

if DIRECT_LIGHT == OFF and M1_CLOSURE == ON:
    lists += ["#define RADTR_M1_ONLY"]
    shutil.copy("./Configs/objects_m1only.defs", "./Configs/objects.defs")
else:
    shutil.copy("./Configs/objects_direct.defs", "./Configs/objects.defs")
    


if M1_CLOSURE == ON:
    lists += ["#define RADTR_M1closer"]
if DIRECT_LIGHT == ON:
    lists += ["#define RADTR_DIRECT"]
    if SUB_GRID_MODEL_DIRECTLIGHT == ON:
        lists += ["#define SUB_GRID_MODEL_DIRECTLIGHT"]


lists += ["#define CFLfac_radtr "+ str(CFL_CONDITION)]
lists += ["#define MAX_STEP_RADTR_TO_HYDRO " + str(MAX_STEP_RAD)]

if SOLVER == "HLL":
    lists += ["#define SOLVER_RADTR_HLL"]
elif SOLVER == "GLF":
    lists += ["!#define SOLVER_RADTR_HLL"]
else:
    print("SOLVER not defined\n")
    sys.exit()
lists += [" "]

sum_components = 0
if EUV_TRANSFER  == ON:
    lists += ["#define M1CLOSER_EUV_TRANSFER"]
    sum_components += 1
    
    if RECOM_ION_PHOTO == ON:
        lists += ["#define RECOM_ION_PHOTO"]

if FUV_TRANSFER == ON:
    lists += ["#define M1CLOSER_FUV_TRANSFER"]
    sum_components += 1

    if SEPARATE_FUV_TRANS == ON:
        lists += ["#define M1CLOSER_SEPARATE_FUV_TRANS"]
        sum_components += 2

    if COLUMN_DENSITY_OPTION == "Jeans":
        lists += ["#define FUV_COLUMNDENS_OPTION 1"]
    elif COLUMN_DENSITY_OPTION == "Sobolev":
        lists += ["#define FUV_COLUMNDENS_OPTION 2"]
    else:
        print("not defined yet this option of COLUMN_DENSITY_OPTION")
        sys.exit()
 
    if EXTERNAL_FUV_RAD == ON:
        lists += ["#define EXTERNAL_FUV_RAD"]
        lists += ["#define EXTERNAL_FUV_G0 "+str(EXTERNAL_G0)]



if SELF_SHIELDING_H2 == "DB96_SIMPLE":
    lists += ["#define SELF_SHIELDING_H2_FORMULA 1"]
elif SELF_SHIELDING_H2 == "WLCOTGREEN19":
    lists += ["#define SELF_SHIELDING_H2_FORMULA 2"]
else:
    print("not defined yet this option of SELF_SHIELDING_H2_FORMULA")
    sys.exit()


if IR_TRANSFER  == ON:
    lists += ["#define M1CLOSER_IR_TRANSFER"]
    sum_components += 1

#sum_components = EUV_TRANSFER + IR_TRANSFER + FUV_TRANSFER + CO_FUV_TRANSFER
lists += ["#define NUM_COMP_TRANSFER "+ str(sum_components)]
lists += [" "]

lists += ["! dust evolution"]
if DUST_NOTCONSTANT == ON:
    lists += ["#define DUST_NOTCONSTANT"]
else:
    lists += ["!#define DUST_NOTCONSTANT"]
lists += [" "]
lists += ["! Radiation pressure"]
if RADIATION_PRESSURE == ON:
    lists += ["#define EXTERNALFORCE"]
else:
    lists += ["!#define EXTERNALFORCE"]

if DUST_ATTENUATION == OFF:
    lists += ["#define SET_NODUST_ATTENUATION"]

if NO_IONIZATION == ON:
    lists += ["#define NO_IONIZATION"]

lists += [" "]
lists += ["! photoelectric heating"]
if PHOTOELECTRIC_HEATING == ON:
    lists += ["#define PHOTOELECTRIC_HEATING YES"]
else:
    lists += ["#define PHOTOELECTRIC_HEATING NO"]

lists += [" "]
if IGNORING_SINK_TRANS_M1 == ON:
    lists += ["#define IGNORING_SINK_TRANS_M1"]


lists += [""]

if METAL_TRANSFER == ON:
    lists += ["! metal transfer"]
    lists += ["#define METAL_TRANSFER"]

if DM_NFW_PROFILE == ON:
    lists += ["#define DM_NFW_PROFILE"]

if DM_POTENTIAL == ON:
    lists += ["! DM potential"]
    lists += ["#define DM_POTENTIAL"]

if COSMIC_RAY == ON:
    lists += ["! cosmic ray"]
    lists += ["#define INCLUDE_COSMICRAY"]
    lists += ["#define COSMICRAY_RATE_ION " + RATE_COSMICRAY]



lists += [""]

# OUTFLOW 
if OUTFLOW == ON:
    lists += ["! OUTFLOW PARAMETERS"]
    lists += ["#define OUTFLOW_ON"] 
    lists += ["#define OUTFLOW_FW " + OUTFLOW_FW ]
    lists += ["#define OUTFLOW_FV " + OUTFLOW_FV ]

lists += [""]

# sink output
if OUTPUTSINK_CRTEPOCH == ON:
    lists += ["! OUTPUT cells around sink particles at the SPCERATIONS"]
    lists += ["#define OUTPUT_SP_CRTEPOCH"]
    lists += ["#define OUTPUT_SP_CRTBOXSIZE " + OUTPUTSINK_CRT_BOXSIZE]


lists += [""]




num = 0
chem_num = []
for label in labels:
    write = "#define " + label + " " + str(num)
    lists += [write]
    if label in labels_chem:
        chem_num += [num]
    num   += 1

lists += ["#define NM "+ str(num)]      
lists += ["#define NCEHM_MIN "+str(min(chem_num))]
lists += ["#define NCEHM_MAX "+str(max(chem_num))]


lists += [" "]
lists += ["#define DEF_TCMB "+ str(DEF_TCMB) ]

if H2PD_HEATING == ON:
    lists += ["#define H2PD_HEATING"]

if HF2020_switch == ON: 
    lists += ["#define CHEM_MODEL_HF2020"]


f = open(filename, 'w')
for lis in lists:
    f.write(lis)
    f.write("\n")
f.close()







