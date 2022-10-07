#!/bin/sh
#PBS -N radtest_IR
#PBS -j oe
#PBS -q large-md
#PBS -l nodes=10
#PBS -l walltime=24:00:00

#
# Run shell 
#
set -x
CRAY_OMP_CHECK_AFFINITY=TRUE; export CRAY_OMP_CHECK_AFFINITY
export OMP_NUM_THREADS=1
NPE=400;                            # process number

BINDIR=/work/$USER/sfumato/ART_sublimation/bin;                #executable directory
DIR=/work/$USER/sfumato/ART_sublimation/data/test_eo/;     export DIR #output directory

PREFIX=st;		export PREFIX	   # prefix of filename
SUFFIX=d;		export SUFFIX	   # suffix of filename
DUMPFILE=dump;		export DUMPFILE	   # file name of restart files
ELAPSELIMIT=20.0;        export ELAPSELIMIT # maximum calculation time in hour
T_LAST=3.D8;		export T_LAST	   # maximum time in yr
#T_LAST=4.2020054D2;	export T_LAST	   # maximum time in yr
Dstep=500;              export Dstep       # data output interval
N0=1.d4;                 export N0          # central density of cloud in cm^-3
T0=1.8d2;               export T0          # temperature of cloud in K
Lmax0=10;               export Lmax0        # initial maximum nest level
sp_radius_in_cell=7;    export sp_radius_in_cell # radius of sink particles wrt finest cell size
sp_Ncr=1.d0;            export sp_Ncr
Mstar=1.d4;             export Mstar       # mass of protostar in M_sun
Mdot=1D-3;              export Mdot        # mdot of protostar in M_sun/yr
Boxsize=6.d1;           export Boxsize     # half of side length of computational box in au
GasType=1;              export GasType     # type of ambient gas: 0 -> HI, 1 -> H2, 2 -> H2&disk-like
Metallicity=1D-1;       export Metallicity
Lstar=1.2d0;             export Lstar
Tsubl=1.d3;            export Tsubl       #1200[K] or 1500[K]
#sinksubl=0.5;           export sinksubl    #ratio of sink to sublimation_radi (0~1,0)

#rm -rf $DIR
mkdir -p $DIR
cd $DIR
NUM=`ls -1 stdout*| wc -l`
STDOUT="stdout"$((NUM))

date 2>&1 | tee -a $STDOUT

#create init data
if [ ! -L $DIR/dump.0.d ]; then
ln -sf $BINDIR/Dust_op  $DIR/Dust_op
ln -sf $BINDIR/ProstFit $DIR/ProstFit
ln -sf $BINDIR/static_rho1D $DIR/static_rho1D
ln -sf $BINDIR/set_restartfiles.sh $DIR/set_restartfiles.sh
ln -sf $BINDIR/remove_old_st_files.sh $DIR/remove_old_st_files.sh
ln -sf $BINDIR/rt $DIR/rt
    aprun -n $NPE $BINDIR/init_test_eo 2>&1 | tee -a $STDOUT
./set_restartfiles.sh $DIR $NPE
fi

aprun -n $NPE $BINDIR/main_test_eo 2>&1 | tee -a $STDOUT
 ./set_restartfiles.sh $DIR $NPE
$BINDIR/remove_old_st_files.sh
date 2>&1 | tee -a $STDOUT
