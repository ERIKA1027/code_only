#! /bin/sh
#
# Setup files for restarting the simulation.
# This script should be executed at the last part of a jobscript.
#
# Calling sequence
# ./set_restartfiles.sh $DIR $NPE $STEP
#
#  st100.10.d -> dump.10.d
#  sp100.d -> sinkparticle.d
#
DIR=$1
NPE=$2
STEP=$3

echo $DIR $NPE $STEP

cd $DIR

## dumpfile for restarting
#LASTDATA=`ls st[0-9]*.0.d | sort -n -r -t t -k2,2 | head -1`
#STEP=`echo $LASTDATA | cut -d t  -f2 | cut -d. -f1`
RANKMIN=0
RANKMAX=`expr $NPE - 1`
for RANK in `seq $RANKMIN $RANKMAX`
do
    DATAFILE=st$STEP.$RANK.d
    SYMLFILE=dump.$RANK.d
    #    if [ -f $SYMLFILE ]; then
    /bin/rm -f $SYMLFILE
    #    fi
    if [ -f $DATAFILE ]; then
	ln -s $DATAFILE $SYMLFILE
    else
	echo "$DATAFILE does not exist!" >&2
    fi
done

## sink particle
DATAFILE=sp$STEP.d
SYMLFILE=sinkparticle.d
#if [ -f $SYMLFILE ]; then
/bin/rm -rf $SYMLFILE
#fi
if [ -f $DATAFILE ]; then
    ln -s $DATAFILE $SYMLFILE
fi

