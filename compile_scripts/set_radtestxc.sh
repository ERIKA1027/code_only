#!/bin/sh
#
# usage:
# path_to_this_script/set_radtest_xc.sh [clean]
#
# note: 
# generate code for collapse of rotating BE sphere to accretion test without rad FB

#move to src directory
SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR/..

#tuning parameters
cp config.h.art config.h
sed -i -e "s/^#define MODEL_ART.*/#define MODEL_ART 2/g" config.h
sed -i -e "s/^#define WITH_SELFGRAVITY.*/!#define WITH_SELFGRAVITY/g" config.h
sed -i -e "s/^!#define NOT_REFINEMENT.*/#define NOT_REFINEMENT/g" config.h
sed -i -e "s/^#define SKIP_HD_UPDATE.*/!#define SKIP_HD_UPDATE/g" config.h
sed -i -e "s/^!#define CREATE_RADSOURCE_BY_HAND.*/#define CREATE_RADSOURCE_BY_HAND/g" config.h
sed -i -e "s/^!#define METAL.*/#define METAL/g" config.h

#for xc
sed -i -e "s/^#define NPE .*/#define NPE 400/g" config.h
sed -i -e "s/^#define NI .*/#define NI 8/g" config.h
sed -i -e "s/^#define NJ .*/#define NJ 8/g" config.h
sed -i -e "s/^#define NK .*/#define NK 8/g" config.h
sed -i -e "s/^#define NGI_BASE .*/#define NGI_BASE 8/g" config.h
sed -i -e "s/^#define NGJ_BASE .*/#define NGJ_BASE 8/g" config.h
sed -i -e "s/^#define NGK_BASE .*/#define NGK_BASE 8/g" config.h
sed -i -e "s/^#define NGID .*/#define NGID 5000/g" config.h


cp Makefile.art.xc Makefile
sed -i -e "s/^INITSRC.*/INITSRC = init_art_radtest/g" Makefile
sed -i -e "s/^MODPAR.*/MODPAR = modelParameter_art_radtest/g" Makefile
sed -i -e "s/^REFINECOND.*/REFINECOND = refineCond_nest_jeans/g" Makefile
sed -i -e "s/^OUTPUTD.*/OUTPUTD = outputdata_art_radtest/g" Makefile
sed -i -e "s/^SP.*/SP      = sinkParticle_art/g" Makefile
#sed -i -e "s/^SUFFIX.*/SUFFIX = _radtest/g" Makefile


#make 
if [ "$1" = "clean" ]; then
    make clean
else
    make
    #copy files for static density profile
    cp -a static_rho1D /work/$USER/sfumato/ART_sublimation/bin/
    #copy files for stellar radiation data
    cp -a ProstFit /work/$USER/sfumato/ART_sublimation/bin/
    #copy files for Dust opacity
    cp -a Dust_op /work/$USER/sfumato/ART_sublimation/bin/
fi

