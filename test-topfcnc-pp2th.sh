#!/usr/bin/bash

PYTHON=python3

TAG=test-topfcnc-pp2th
OUTDIR=$TAG
MT=172.5
MYCOEFF=RCtphi


### redirect output to log file
exec >> ${TAG}.log  2>&1

### setup madgraph
VERSION=354
MGLINK=https://github.com/mg5amcnlo/mg5amcnlo/archive/refs/tags/v3.5.4.tar.gz
MGTAR=mg5_v${VERSION}.tar.gz
if [[ ! -f "$MGTAR" ]] ; then

### get MG
date && echo "--- Get madgraph"
wget -O $MGTAR $MGLINK
tar xzf $MGTAR
MGBASE=`tar tzf $MGTAR | grep -m1 'madgraph/$' | sed 's/\/madgraph\///'`
MG="$PYTHON $MGBASE/bin/mg5_aMC"
date && echo "--- Using MG executable: $MG"


### get utilities
date && echo "--- Install utilities"
echo "install ninja
install collier
" > ${OUTDIR}.cmd
time $MG -f ${OUTDIR}.cmd

# end of setup MG
fi



### setup model
MODEL=./TopFCNC
RESTRICTION=$MYCOEFF
if  [[ ! -f TopFCNC.tar.gz ]] ; then

### get model
wget https://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/TopFCNC/TopFCNC.tar.gz
tar xzf TopFCNC.tar.gz

### make restriction
cd $MODEL
cp restrict_onlyu.dat param_card.dat
#python2 ./write_param_card.py
# set all coeff to zero
for COEFF in RCtphi ICtphi RCtG ICtG RCtW ICtW RCtB ICtB RCuphi ICuphi RCuG ICuG RCuW ICuW RCuB ICuB RC1utR IC1utR RC1utL IC1utL RC3utL IC3utL RCtcphi ICtcphi RCtcG ICtcG RCtcW ICtcW RCtcB ICtcB RCctphi ICctphi RCctG ICctG RCctW ICctW RCctB ICctB RC1ctR IC1ctR RC1ctL IC1ctL RC3ctL IC3ctL ; do
sed -i "s/[0-9.e+-]\+ *# *$COEFF/0. # $COEFF/" param_card.dat
done
# set MYCOEFF to one
sed -i "s/[0-9.e+-]\+ *# *$MYCOEFF/1. # $MYCOEFF/" param_card.dat
mv param_card.dat restrict_${RESTRICTION}.dat
cd $OLDPWD

fi


### MG
if [[ "$MG" == "" ]] ; then
MGBASE=`tar tzf $MGTAR | grep -m1 'madgraph/$' | sed 's/\/madgraph\///'`
MG="$PYTHON $MGBASE/bin/mg5_aMC"
date && echo "--- MGBASE: $MGBASE"
fi



### only generate if output directory is absent
if [[ ! -d "${OUTDIR}" ]] ; then

### generate and output
echo "set acknowledged_v3.1_syntax True --global
import model ${MODEL}-${RESTRICTION}
generate p p > t h \$\$ t~ QCD=1 QED=1 NP=2 [QCD]
output $OUTDIR" > ${TAG}.cmd

$MG -f ${TAG}.cmd

fi


### launch
echo "set automatic_html_opening False
launch $OUTDIR
fixed_order=ON
order=NLO
done
set fixed_ren_scale True           # central scales are fixed to MT
set fixed_fac_scale True
set MT              $MT
set muR_ref_fixed   $MT
set muF_ref_fixed   $MT
0
exit
exit" > ${TAG}.cmd

$MG -f ${TAG}.cmd


### results should be similar to
#   --------------------------------------------------------------
#      Final results and run summary:
#      Process p p > t h $$ t~ QCD=1 QED=1 NP=2 [QCD]
#      Run at p-p collider (6500.0 + 6500.0 GeV)
#      Total cross section: 3.190e-01 +- 1.3e-03 pb
#   --------------------------------------------------------------
#      Scale variation (computed from histogram information):
#          Dynamical_scale_choice -2 (envelope of 9 values):
#              3.190e-01 pb  +7.2% -7.2%
#   --------------------------------------------------------------
