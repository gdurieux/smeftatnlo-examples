#!/usr/bin/bash

MYCOEFF=cQq13

### cQq13 results should be similar to:
#order 4                             43.3 pb    -1%  +0.5%       (±0.5%)
#order 20002                        -40.3 pb    +9%   -10%       (±0.4%)
#order 40000                          141 pb   -10%    +8%       (±0.3%)
#order 204                           3.05 pb    -9%   +40%       (±10%)
#order 20202                         4.44 pb  -100%   +90%       (±5%)
#order 40200                        -2.45 pb  +400%  -500%       (±30%)
##
#order 4 + order 20002 + order 40000 + order 204 + order 20202 + order 40200          149 pb    -2%    +1%       (±1%)
#order 40000 + order 40200            139 pb    -3%    +3%       (±0.8%)
#order 20002 + order 20202          -35.9 pb    +3%    -1%       (±1%)
#order 4 + order 204                 46.3 pb    -1%    +3%       (±1%)


### cQq83 results should be similar to:
#order 4                             9.54 pb    -1%  +0.4%       (±0.6%)
#order 20002                            0 pb
#order 40000                          141 pb   -10%    +8%       (±0.3%)
#order 204                         -0.568 pb   +10%   -30%       (±20%)
#order 20202                         2.61 pb   -20%   +20%       (±2%)
#order 40200                        -2.24 pb  +400%  -600%       (±30%)
##
#order 4 + order 20002 + order 40000 + order 204 + order 20202 + order 40200          150 pb    -2%    +3%       (±0.9%)
#order 40000 + order 40200            139 pb    -2%    +3%       (±0.8%)
#order 20002 + order 20202           2.61 pb   -20%   +20%       (±2%)
#order 4 + order 204                 8.97 pb  -0.7%    +1%       (±2%)


#VERSION=310
#MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.1.x/+download/MG5_aMC_v3.1.0.tar.gz

#VERSION=332
#MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.3.x/+download/MG5_aMC_v3.3.2.tar.gz

#VERSION=340
#MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.3.x/+download/MG5_aMC_v3.4.0.RC3.tar.gz

#VERSION=342
#MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.4.x/+download/MG5_aMC_v3.4.2.tar.gz

#VERSION=351
#MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.1.tar.gz

VERSION=353
MGLINK=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.3.tar.gz

OUTDIR=gen_single_top_4f_${MYCOEFF}_v${VERSION}
exec >> ${OUTDIR}.log  2>&1

# need to specify python3.8 for CS8, while CS9 has python3.9 as default already
V=3
PYTHON=python${V}
F2PY=f2py${V}

INIT="#
#set run_mode 2                                  # !!! 1 is for cluster, 2 for multicore  
#set nb_core 10
set f2py_compiler ${F2PY}
set auto_convert_model True                     # convert model to python3 automatically
set acknowledged_v3.1_syntax True --global      # needed for v3.1 and above
save options"                                   # write all that in config file


### only set up madgraph and utilities once, assumes it is all done if the tar.gz file is there
MGTAR=mg5_v${VERSION}.tar.gz
if [[ ! -f "$MGTAR" ]] ; then

### get MG
date
wget -O $MGTAR $MGLINK
tar xzf $MGTAR
MG="$PYTHON `tar tzf $MGTAR | grep mg5_aMC`"



### get utilities
date
echo "$INIT
#install lhapdf6
install ninja
install collier
" > ${OUTDIR}.cmd
time $MG -f ${OUTDIR}.cmd

### end of madgraph setup
fi


### need to point to LHAPDF (required for systematics reweighting, and lhapdf)
MGCONFIG="`tar tzf $MGTAR | grep '/mg5_configuration.txt'`"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`grep -o "[^ ]\+/lhapdf6[^/]\+" ${MGCONFIG}`/lib/
echo "--- Using LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

### MG executable and initialisation
MG="$PYTHON `tar tzf $MGTAR | grep mg5_aMC`"
echo "--- Using MG executable: $MG"


### get model
MODELTARBALL=SMEFTatNLO_v1.0.3.tar.gz
MODEL="./SMEFTatNLO"
if [[ ! -d "$MODEL" ]] ; then
wget https://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/SMEFTatNLO/${MODELTARBALL}
tar xzf ${MODELTARBALL}
fi


### create restriction card for one operator coefficient
RESTRICTION=NLO_$MYCOEFF
COEFFS=($MYCOEFF)
RANDOMS=(0.97410347 0.42856892 0.32560595 0.18687959 0.66700969 0.95601325 0.34684861 0.98139281 0.57624127 0.48050942 0.59998608 0.5559571 0.80537718 0.83778932 0.70840278 0.49491561 0.7044201 0.95969375 0.53159833 0.68032797 0.99473372 0.88755928 0.94969546 0.42232173 0.61368649 0.95255866 0.0356684 0.54192876 0.84393958 0.07092631)

cd ${MODEL}
cp restrict_NLO.dat restrict_${RESTRICTION}.dat

# set all coeffs to zero
for COEFF in c3pl1 c3pl2 c3pl3 cblS3 cdp cG cll1111 cll1122 cll1133 cll1221 cll1331 cll2222 cll2233 cll2332 cll3333 cp cpBB cpd cpDC cpe cpG cpl1 cpl2 cpl3 cpmu cpQ3 cpq3i cpQM cpqMi cpt cpta cpu cpW cpWB cQd1 cQd8 cQe1 cQe2 cQe3 cQl31 cQl32 cQl33 cQlM1 cQlM2 cQlM3 cQQ1 cQq11 cQq13 cQQ8 cQq81 cQq83 cQt1 cQt8 cQu1 cQu8 ctd1 ctd8 cte1 cte2 cte3 ctG ctl1 ctl2 ctl3 ctlS3 ctlT3 ctp ctq1 ctq8 ctt1 ctu1 ctu8 ctW ctZ cWWW ; do
sed -i "s/[0-9.e+-]\+ # ${COEFF} /0. # ${COEFF} /i" restrict_${RESTRICTION}.dat
done

# set COEFFS to random numbers
for ICOEFF in `seq 0 $((${#COEFFS[*]}-1))`; do
sed -i "s/[0-9.e+-]\+ # ${COEFFS[$ICOEFF]} /${RANDOMS[$ICOEFF]} # ${COEFFS[$ICOEFF]} /" restrict_${RESTRICTION}.dat
done

# set scales
SCALE=173.
sed -i "s/[0-9.e+-]\+ # MU_R /${SCALE} # MU_R /" restrict_${RESTRICTION}.dat
sed -i "s/[0-9.e+-]\+ # mueft /${SCALE} # mueft /" restrict_${RESTRICTION}.dat


cd ../



### loop filtering for cQq83, allowing all gluon loops
echo "--- Setting up loop filtering for cQq83, allowing all gluon loops"
LOOPS=`tar tzf $MGTAR | grep madgraph/loop/loop_diagram_generation.py`
if [[ ! -f ${LOOPS}.bkp ]]; then
cp $LOOPS ${LOOPS}.bkp
else
cp ${LOOPS}.bkp $LOOPS
fi
read -d '' PATCH <<EOF
=== modified file '$LOOPS'
--- $LOOPS	2020-03-11 09:28:14 +0000
+++ $LOOPS	2020-04-03 21:08:18 +0000
@@ -384,7 +384,7 @@
         # By default the user filter does nothing if filter is not set, 
         # if you want to turn it on and edit it by hand, then set the 
         # variable edit_filter_manually to True
-        edit_filter_manually = False 
+        edit_filter_manually = True 
         if not edit_filter_manually and filter in [None,'None']:
             return
         if isinstance(filter,str) and  filter.lower() == 'true':
@@ -415,6 +415,10 @@
                     raise InvalidCmd("The user-defined filter '%s' did not"%filter+
                                  " returned the following error:\n       > %s"%str(e))
 
+            # requires a gluon to run in all loops
+            if (21 not in diag.get_loop_lines_pdgs()) and (9000005 not in diag.get_loop_lines_pdgs()):
+                valid_diag = False
+
 #            if any([abs(pdg) not in range(1,7) for pdg in diag.get_loop_lines_pdgs()]):
 #                valid_diag = False
 
@@ -538,7 +542,7 @@
             
             if valid_diag:
                 newloopselection.append(diag)
-        self['loop_diagrams']=newloopselection
+        #self['loop_diagrams']=newloopselection
         # To monitor what are the diagrams filtered, simply comment the line
         # directly above and uncomment the two directly below.
 #        self['loop_diagrams'] = base_objects.DiagramList(
EOF
if [[ "$MYCOEFF" = "cQq83" ]]; then
patch -p0 <<< "$PATCH"
fi


### generate diagrams, and output directory if not already present
echo "import model ${MODEL}-${RESTRICTION}
generate p p > t j \$$ W- NP=2 QCD=0 QED=2 [QCD]
output ${OUTDIR}
y                                 # just in case some installation or overwritting is needed
" > ${OUTDIR}.cmd
if [[ ! -d "${OUTDIR}" ]] ; then
	date
	echo "--- Generate and output"
	time $MG -f ${OUTDIR}.cmd
fi

### revert madgraph patch accepting all gluon loops
if [[ "$MYCOEFF" = "cQq83" ]]; then
patch -R -p0 <<< "$PATCH"
fi


### Lower pole check threshold
#date
#echo "--- Lower IRPoleCheckThreshold"
#sed -i '/#IRPoleCheckThreshold/{N;s/\n.\+/\n1.0d-4/}' ${OUTDIR}/Cards/FKS_params.dat
#echo "--- Grep that"
#grep -A3 "#IRPoleCheckThreshold" ${OUTDIR}/Cards/FKS_params.dat



### patch fixed order analysis to extract all orders
date
echo "--- Patching fixed order analysis"
echo "--- The ordering of couplings is assumed to be NP, QCD, QED so that tag = NP + 100*QCD + 10000* QED"
grep -m1 "the order of the coupling orders is" ${OUTDIR}/SubProcesses/*/orders.inc

ORDERS=(4 20002 40000 204 20202 40200)
BOOKS=""
CASES=""
for ICOEFF in `seq 0 $((${#ORDERS[*]}-1))`; do

if [[ "$BOOKS" != "" ]]; then
BOOKS="$BOOKS
"
fi
BOOKS="$BOOKS+      call HwU_book($((1+$ICOEFF)),'order ${ORDERS[$ICOEFF]}', 1,0.d0,1.d0)"

if [[ "$CASES" != "" ]]; then
CASES="$CASES
"
fi
CASES="$CASES+         case (${ORDERS[$ICOEFF]})
+            call HwU_fill($((1+$ICOEFF)),var,wgts)"

done

echo "--- BOOKS"
echo "$BOOKS"

echo "--- CASES"
echo "$CASES"

echo "--- patching"
cp `tar tzf $MGTAR | grep /FixedOrderAnalysis/analysis_HwU_template.f` ${OUTDIR}/FixedOrderAnalysis/analysis_HwU_template.f
patch -p0 <<< "--- ${OUTDIR}/FixedOrderAnalysis/analysis_HwU_template.f	2023-01-20 16:22:52.000000000 +0100
+++ ${OUTDIR}/FixedOrderAnalysis/analysis_HwU_template.f	2023-11-02 23:43:26.137064453 +0100
@@ -38,8 +38,$(( 8 - 2 + ${#ORDERS[*]} )) @@
 c weights and the information on the weights:
       call HwU_inithist(nwgt,weights_info)
 c declare (i.e. book) the histograms
-      call HwU_book(1,'total rate      ', 5,0.5d0,5.5d0)
-      call HwU_book(2,'total rate Born ', 5,0.5d0,5.5d0)
$BOOKS
       return
       end
 
@@ -104,16 +$(( 104 - 2 + ${#ORDERS[*]} )),$(( 16 + 4 - 5 + 3 + 2 * ${#ORDERS[*]} )) @@
       integer ibody
 c local variable
       double precision var
+c orders
+      integer orders_tag_plot
+      integer iii
+      common /corderstagplot/ orders_tag_plot
 c
 c Fill the histograms here using a call to the HwU_fill()
 c subroutine. The first argument is the histogram label, the second is
 c the numerical value of the variable to plot for the current
 c phase-space point and the final argument is the weight of the current
 c phase-space point.
-      var=1d0
-c always fill the total rate
-      call HwU_fill(1,var,wgts)
-c only fill the total rate for the Born when ibody=3
-      if (ibody.eq.3) call HwU_fill(2,var,wgts)
+      var=.5d0
+      select case (orders_tag_plot)
$CASES
+      end select
       return
       end"


### set COEFFS 1.
SETS=""
for COEFF in ${COEFFS[*]}; do
SETS="$SETS#
set $COEFF 1.  "
done
echo "--- SETS"
echo "$SETS"


### launch
date
echo "launch ${OUTDIR}
fixed_order=ON
order=NLO
done
#
set fixed_ren_scale True
set fixed_fac_scale True
set mur_ref_fixed ${SCALE}
set muf_ref_fixed ${SCALE}
set mueft ${SCALE}
set WW 0.                                # set W width to zero for gauge invariance
$SETS                                    # set coeffs 
0
exit" > ${OUTDIR}.cmd

echo "--- Launch"
time $MG -f ${OUTDIR}.cmd


### get histogram information
echo "--- Get histogram information"
VARIOUS=$(dirname `tar tzf $MGTAR | grep madgraph/various/histograms.py`)
PYCODE="import sys, glob
import numpy as np
here = '$VARIOUS/'
if here not in sys.path:
    sys.path.append(here)
import histograms as hist

def fun(v):
    return(np.round(v,decimals=-int(np.floor(np.log10(np.abs(v))))))

def printout(dat,names):
    datsum = np.sum([dat[name] for name in names],axis=0)
    c = datsum[0]
    e = datsum[1]
    d = np.min(datsum[2:])
    u = np.max(datsum[2:])
    if c==0.:
        print('{:25}\t{:8.3g} pb'.format(
            ' + '.join(names),
            c
            ))
    else:
        print('{:25}\t{:8.3g} pb {:+5.3g}% {:+5.3g}%\t(±{:<.2g}%)'.format(
            ' + '.join(names),
            c,
            fun((d/c-1)*100),
            fun((u/c-1)*100),
            fun(e/np.abs(c)*100)
            ))
for hwu in glob.glob('$OUTDIR/Events/*/MADatNLO.HwU'):
    print(hwu)
    file = hist.HwUList(hwu)
    dat = {}
    dat['entries'] = file.get_wgt_names()
    for name in file.get_hist_names():
        newname = name.replace(', #1','')
        dat[newname] = np.array(
                    [file.get(name).get('central')[0],
                     file.get(name).get('stat_error')[0]]+\
                    [file.get(name).get(entry)[0] for entry in dat['entries'] if 'scale_adv' in entry])
        printout(dat,[newname])
    print('#')
    printout(dat,[name for name in dat if name!='entries'])
    if 'order 40000' in dat and 'order 40200' in dat:
        printout(dat,['order 40000','order 40200'])
    if 'order 20002' in dat and 'order 20202' in dat:
        printout(dat,['order 20002','order 20202'])
    if 'order 4' in dat and 'order 204' in dat:
        printout(dat,['order 4','order 204'])

"
$PYTHON -c "$PYCODE"

