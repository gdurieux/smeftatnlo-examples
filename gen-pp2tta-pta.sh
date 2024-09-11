#!/usr/bin/bash

PYTHON=python3

TAG=gen-pp2tta-pta
OUTDIR=$TAG
MT=172.5
MYCOEFF=ctZ


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



### get model
MODELTARBALL=SMEFTatNLO_v1.0.3.tar.gz
MODEL="./SMEFTatNLO"

if [[ ! -d "$MODEL" ]] ; then
wget https://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/SMEFTatNLO/${MODELTARBALL}
tar xzf ${MODELTARBALL}

### make restriction
RESTRICTION=NLO_${MYCOEFF}
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
generate p p > t t~ a QCD=2 QED=1 NP=2 [QCD]
output $OUTDIR" > ${TAG}.cmd

$MG -f ${TAG}.cmd



### set fixed order analysis
date && echo "--- Patching fixed order analysis"
echo "--- The ordering of couplings is assumed to be NP, QCD, QED so that tag = NP + 100*QCD + 10000* QED"
cat << EOF > ${OUTDIR}/FixedOrderAnalysis/analysis_HwU_template.f
c
c This file contains the default histograms for fixed order runs: it
c only plots the total rate as an example. It can be used as a template
c to make distributions for other observables.
c
c This uses the HwU package and generates histograms in the HwU/GnuPlot
c format. This format is human-readable. After running, the histograms
c can be found in the Events/run_XX/ directory.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called once at the start of each run. Here the
c histograms should be declared. 
c
c Declare the histograms using 'HwU_book'.
c     o) The first argument is an integer that labels the histogram. In
c     the analysis_end and analysis_fill subroutines this label is used
c     to keep track of the histogram. The label should be a number
c     starting at 1 and be increased for each plot.
c     o) The second argument is a string that will apear above the
c     histogram. Do not use brackets "(" or ")" inside this string.
c     o) The third, forth and fifth arguments are the number of bis, the
c     lower edge of the first bin and the upper edge of the last
c     bin.
c     o) When including scale and/or PDF uncertainties, declare a
c     histogram for each weight, and compute the uncertainties from the
c     final set of histograms
c
      implicit none
c When including scale and/or PDF uncertainties the total number of
c weights considered is nwgt
      integer nwgt
c In the weights_info, there is an text string that explains what each
c weight will mean. The size of this array of strings is equal to nwgt.
      character*(*) weights_info(*)
c Initialize the histogramming package (HwU). Pass the number of
c weights and the information on the weights:
      call HwU_inithist(nwgt,weights_info)
c declare (i.e. book) the histograms
      call HwU_book(1,'rate SM NLO',   1,0.d0,1.d0)
      call HwU_book(2,'rate SM LO',    1,0.d0,1.d0)
      call HwU_book(3,'rate lin NLO',  1,0.d0,1.d0)
      call HwU_book(4,'rate lin LO',   1,0.d0,1.d0)
      call HwU_book(5,'rate quad NLO', 1,0.d0,1.d0)
      call HwU_book(6,'rate quad LO',  1,0.d0,1.d0)
      
      call HwU_book( 7,'pta SM NLO',   101, -5.d0, 500.d0)
      call HwU_book( 8,'pta SM LO' ,   101, -5.d0, 500.d0)
      call HwU_book( 9,'pta lin NLO',  101, -5.d0, 500.d0)
      call HwU_book(10,'pta lin LO' ,  101, -5.d0, 500.d0)
      call HwU_book(11,'pta quad NLO', 101, -5.d0, 500.d0)
      call HwU_book(12,'pta quad LO' , 101, -5.d0, 500.d0)
      
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called once at the end of the run. Here the
c histograms are written to disk. Note that this is done for each
c integration channel separately. There is an external script that will
c read the HwU data files in each of the integration channels and
c combines them by summing all the bins in a final single HwU data file
c to be put in the Events/run_XX directory, together with a gnuplot
c file to convert them to a postscript histogram file.
      implicit none
      double precision dummy
      call HwU_write_file
      return                
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called for each n-body and (n+1)-body configuration
c that passes the generation cuts. Here the histrograms are filled.
      implicit none
c This includes the 'nexternal' parameter that labels the number of
c particles in the (n+1)-body process
      include 'nexternal.inc'
c This is an array which is '-1' for initial state and '1' for final
c state particles
      integer istatus(nexternal)
c This is an array with (simplified) PDG codes for the particles. Note
c that channels that are combined (i.e. they have the same matrix
c elements) are given only 1 set of PDG codes. This means, e.g., that
c when using a 5-flavour scheme calculation (massless b quark), no
c b-tagging can be applied.
      integer iPDG(nexternal)
c The array of the momenta and masses of the initial and final state
c particles in the lab frame. The format is "E, px, py, pz, mass", while
c the second dimension loops over the particles in the process. Note
c that these are the (n+1)-body particles; for the n-body there is one
c momenta equal to all zero's (this is not necessarily the last particle
c in the list). If one uses IR-safe obserables only, there should be no
c difficulty in using this.
      double precision p(0:4,nexternal)
c The weight of the current phase-space point is wgts(1). If scale
c and/or PDF uncertainties are included through reweighting, the rest of
c the array contains the list of weights in the same order as described
c by the weigths_info strings in analysis_begin
      double precision wgts(*)
c The ibody variable is:
c     ibody=1 : (n+1)-body contribution
c     ibody=2 : n-body contribution (excluding the Born)
c     ibody=3 : Born contribution
c The histograms need to be filled for all these contribution to get the
c physics NLO results. (Note that the adaptive phase-space integration
c is optimized using the sum of the contributions, therefore plotting
c them separately might lead to larger than expected statistical
c fluctuations).
      integer ibody
c local variable
      double precision var

c orders
      integer orders_tag_plot
      integer iii
      common /corderstagplot/ orders_tag_plot

c photon pt
      integer i
      real*8 pt, pta
      external pt

c find the photons
      pta=-1.d0
      do i=1,nexternal
         if (istatus(i).eq.1 .and. ipdg(i).eq.22 ) then
            pta=pt(p(0,i))
         endif
      enddo


c NP + 100* QCD + 10000* QED
c 2040{0,2,4} is LO
c 2060{0,2,4} is NLO
      var = 0.5d0
      select case (orders_tag_plot)
         case (20400)
            call HwU_fill( 1,var,wgts)
            call HwU_fill( 2,var,wgts)
            call HwU_fill( 7,pta,wgts)
            call HwU_fill( 8,pta,wgts)
         case (20600)
            call HwU_fill( 1,var,wgts)
            call HwU_fill( 7,pta,wgts)
         case (20402)
            call HwU_fill( 3,var,wgts)
            call HwU_fill( 4,var,wgts)
            call HwU_fill( 9,pta,wgts)
            call HwU_fill(10,pta,wgts)
         case (20602)
            call HwU_fill( 3,var,wgts)
            call HwU_fill( 9,pta,wgts)
         case (20404)
            call HwU_fill( 5,var,wgts)
            call HwU_fill( 6,var,wgts)
            call HwU_fill(11,pta,wgts)
            call HwU_fill(12,pta,wgts)
         case (20604)
            call HwU_fill( 5,var,wgts)
            call HwU_fill(11,pta,wgts)
      end select
      

      return
      end


EOF


# end of generate and output
fi




### launch
if [[ ! -d ${OUTDIR}/Events/run_01 ]] ; then
echo "set automatic_html_opening False
launch $OUTDIR
fixed_order=ON
order=NLO
done
set req_acc_fo      0.001          # high precision
set $MYCOEFF        -1.            # !!!! negative value, since interference is negative !!!!
set fixed_ren_scale True           # central scales are fixed to MT
set fixed_fac_scale True
set MT              $MT
set muR_ref_fixed   $MT
set muF_ref_fixed   $MT
set mueft           $MT
#set quarkphreco    True
set ptgmin          20.            # default photon pt cut
set etagamma        2.5            # photon eta cut
0
exit
exit" > ${TAG}.cmd

$MG -f ${TAG}.cmd
fi




### get histogram information
date && echo "--- Get histogram information"
VARIOUS=$MGBASE/madgraph/various
PYCODE="import sys, glob
import numpy as np
import matplotlib.pyplot as plt
here = '$VARIOUS/'
if here not in sys.path:
    sys.path.append(here)
import histograms as hist

def plot(x,y,fmt=None,**args):
    X = x
    Y = y.repeat(2)
    if fmt:
        return plt.plot(X,Y,fmt,**args)
    else:
        return plt.plot(X,Y,**args)

def plotfill(x,y,z,**args):
    X = x
    Y = y.repeat(2)
    Z = z.repeat(2)
    return plt.fill_between(X,Y,Z,linewidth=0,**args)

# convert from pb to fb
fac = 1e3

# loop over histogram files
for fname in sorted(glob.glob('${OUTDIR}/Events/run_*/MADatNLO.HwU')):
    print(fname)
    
    ## collect data
    dat = {}
    
    ## rate panel
    plt.subplots(2,1,figsize=(6,4*2))
    plt.subplot(2,1,1)
    
    ## loop over the different histograms for each order
    file = hist.HwUList(fname)
    for name in file.get_hist_names():
    
        # just pick the pT(photon) histograms
        if 'pta' not in name:
            continue
        h = file.get(name)
        
        # group n bins
        h.rebin(5)
        
        # bins, central value, scale up/down uncertainties, and statistical error
        x = np.array([bb for b in h.bins for bb in b.boundaries])
        y = np.array( h.get('central') )*fac
        scales  = np.array([h.get(entry) for entry in file.get_wgt_names()
            if entry[0] == 'scale_adv' and entry[2] in [0.5,1,2] and entry[3] in [0.5,1,2] ])
        u = np.max( scales, axis=0 )*fac
        d = np.min( scales, axis=0 )*fac
        e = np.array(h.get('stat_error'))*fac
        
        dat[name] = [x,y,u,d,e]
        
            
        # the interference is negative, so take an absolute value for this log-scale plot
        def fun(z):
            return np.abs(z)
        
        # central value with nice colour and line style
        args = {}
        if ' LO,' in name:
            args['color'] = p[0].get_color()
            args['fmt'] = '--'
        p = plot(x,fun(y),
            label=name.replace('pta ','').replace(', #1','').replace('lin','-lin'),
            **args)
        
        # scale and statistical bands
        if ' LO' not in name:
            plotfill(x,fun(u),fun(d), color=p[0].get_color(), alpha=.3)
            plotfill(x,fun(y+e),fun(y-e), color='black', alpha=.4)
        
    
    plt.yscale('log')
    plt.xlim(20,370)
    plt.ylim(.1,2e3)
    plt.legend(loc='upper right')
    plt.xlabel('\$p_T(\\gamma)\$ [GeV]')
    plt.title('$MYCOEFF dependence of \$p p \\\to t \\\bar{t} \\\gamma\$ [fb] at \$\\sqrt{s}=13\$ TeV')
    
    print(repr(dat))
    
    ## k-factor panel
    plt.subplot(2,1,2)
    for tag in ['SM', 'lin', 'quad']:
        nlo = dat['pta {:} NLO, #1'.format(tag)]
        lo  = dat['pta {:} LO, #1'.format(tag)]
        p = plot(nlo[0], nlo[1]/lo[1], label=tag)
        plotfill(nlo[0],nlo[2]/lo[1],nlo[3]/lo[1], color=p[0].get_color(), alpha=.3)
        plotfill(nlo[0],(nlo[1]+nlo[4])/lo[1],(nlo[1]-nlo[4])/lo[1], color='black', alpha=.4)
    plt.ylabel('K factor')
    plt.xlim(20,370)
    plt.ylim(0,2)
    
    plt.savefig('${OUTDIR}_{:}.pdf'.format(fname.split('/')[-2]))
    #plt.show()
    
    
    
    ##### specific coefficient values
    
    nlo0 = dat['pta SM NLO, #1']
    lo0  = dat['pta SM LO, #1']
    nlo1 = dat['pta lin NLO, #1']
    lo1  = dat['pta lin LO, #1']
    nlo2 = dat['pta quad NLO, #1']
    lo2  = dat['pta quad LO, #1']
    
    plt.subplots(2,1,figsize=(6,4*2))
    plt.subplot(2,1,1)
    c = 0
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/nlo0[1], 'k-', label='SM')
    
    c = 0.5
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/nlo0[1], label='ctZ/Λ²={:}/TeV²'.format(c))
    plot(lo0[0], (lo0[1]+c*lo1[1]+c*c*lo2[1])/lo0[1], '--', color=p[0].get_color())
    
    c = -0.5
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/nlo0[1], label='ctZ/Λ²={:}/TeV²'.format(c))
    plot(lo0[0], (lo0[1]+c*lo1[1]+c*c*lo2[1])/lo0[1], '--', color=p[0].get_color())
    
    plt.xlim(20,370)
    plt.ylim(.95,1.5)
    plt.ylabel('EFT / SM')
    plt.legend(loc='upper left')
    plt.xlabel('\$p_T(\\gamma)\$ [GeV]')
    
    
    plt.subplot(2,1,2)
    c = 0
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/(lo0[1]+c*lo1[1]+c*c*lo2[1]), 'k-', label='SM')
    c = 0.5
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/(lo0[1]+c*lo1[1]+c*c*lo2[1]), label='ctZ/Λ²={:}/TeV²'.format(c))
    c = -0.5
    p = plot(nlo0[0], (nlo0[1]+c*nlo1[1]+c*c*nlo2[1])/(lo0[1]+c*lo1[1]+c*c*lo2[1]), label='ctZ/Λ²={:}/TeV²'.format(c))
    
    
#    plt.legend(loc='upper right')
    plt.ylabel('K factor')
    plt.xlim(20,370)
    plt.ylim(1.25,1.55)
    
    plt.savefig('${OUTDIR}_{:}_ctZ05.pdf'.format(fname.split('/')[-2]))
    
"
$PYTHON -c "$PYCODE"



### resulting histogram file should look like the following:
cat << EOF > MADatNLO_forcomparison.HwU
##& xmin & xmax & central value & dy & delta_mu_cen -2 @aux & delta_mu_min -2 @aux & delta_mu_max -2 @aux & dyn=-2 muR= 1.000 muF= 1.000 & dyn=-2 muR= 2.000 muF= 1.000 & dyn=-2 muR= 0.500 muF= 1.000 & dyn=-2 muR= 1.000 muF= 2.000 & dyn=-2 muR= 2.000 muF= 2.000 & dyn=-2 muR= 0.500 muF= 2.000 & dyn=-2 muR= 1.000 muF= 0.500 & dyn=-2 muR= 2.000 muF= 0.500 & dyn=-2 muR= 0.500 muF= 0.500

<histogram> 1 "rate SM NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   +2.4252143e+00   +2.3390015e-03  +2.4252143e+00   +2.0974929e+00   +2.7223627e+00   +2.4252143e+00   +2.1918154e+00   +2.6629431e+00   +2.3376837e+00   +2.0974929e+00   +2.5939510e+00   +2.5060998e+00   +2.2814089e+00   +2.7223627e+00
<\histogram>


<histogram> 1 "rate SM LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   +1.6376173e+00   +4.9546769e-04  +1.6376173e+00   +1.2553741e+00   +2.1656522e+00   +1.6376173e+00   +1.3643410e+00   +2.0037476e+00   +1.5068244e+00   +1.2553741e+00   +1.8437127e+00   +1.7699383e+00   +1.4745811e+00   +2.1656522e+00
<\histogram>


<histogram> 1 "rate lin NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   -1.0889546e-01   +2.3323399e-04  -1.0889546e-01   -1.1254058e-01   -9.9342947e-02   -1.0889546e-01   -1.0256158e-01   -1.1218035e-01   -1.0691070e-01   -9.9342947e-02   -1.1254058e-01   -1.1014683e-01   -1.0516334e-01   -1.1093294e-01
<\histogram>


<histogram> 1 "rate lin LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   -9.3065188e-02   +3.4494256e-05  -9.3065188e-02   -1.2338522e-01   -7.0825467e-02   -9.3065188e-02   -7.7535001e-02   -1.1387224e-01   -8.5011745e-02   -7.0825467e-02   -1.0401825e-01   -1.0083993e-01   -8.4012337e-02   -1.2338522e-01
<\histogram>


<histogram> 1 "rate quad NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   +3.5483928e-01   +2.1115552e-04  +3.5483928e-01   +3.1004640e-01   +3.7966619e-01   +3.5483928e-01   +3.2692455e-01   +3.7851042e-01   +3.4151902e-01   +3.1004640e-01   +3.7251036e-01   +3.6532270e-01   +3.4220913e-01   +3.7966619e-01
<\histogram>


<histogram> 1 "rate quad LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  +0.0000000e+00   +1.0000000e+00   +2.6897715e-01   +7.5605440e-05  +2.6897715e-01   +1.9760094e-01   +3.7126957e-01   +2.6897715e-01   +2.2409180e-01   +3.2911378e-01   +2.3718023e-01   +1.9760094e-01   +2.9020784e-01   +3.0343016e-01   +2.5279546e-01   +3.7126957e-01
<\histogram>


<histogram> 101 "pta SM NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   +3.4983048e-01   +1.4823981e-03  +3.4983048e-01   +3.0393101e-01   +3.9251364e-01   +3.4983048e-01   +3.1604120e-01   +3.8433983e-01   +3.3841305e-01   +3.0393101e-01   +3.7499634e-01   +3.6025093e-01   +3.2729246e-01   +3.9251364e-01
  +2.5000000e+01   +3.0000000e+01   +2.6957949e-01   +1.4183699e-03  +2.6957949e-01   +2.3395910e-01   +3.0266130e-01   +2.6957949e-01   +2.4346494e-01   +2.9630879e-01   +2.6062718e-01   +2.3395910e-01   +2.8900108e-01   +2.7773981e-01   +2.5230302e-01   +3.0266130e-01
  +3.0000000e+01   +3.5000000e+01   +2.1616475e-01   +1.1162894e-03  +2.1616475e-01   +1.8728388e-01   +2.4314015e-01   +2.1616475e-01   +1.9505126e-01   +2.3790670e-01   +2.0884002e-01   +1.8728388e-01   +2.3190937e-01   +2.2286830e-01   +2.0230322e-01   +2.4314015e-01
  +3.5000000e+01   +4.0000000e+01   +1.7805809e-01   +7.5449503e-04  +1.7805809e-01   +1.5400229e-01   +2.0068074e-01   +1.7805809e-01   +1.6062384e-01   +1.9604352e-01   +1.7177131e-01   +1.5400229e-01   +1.9081549e-01   +1.8387094e-01   +1.6685660e-01   +2.0068074e-01
  +4.0000000e+01   +4.5000000e+01   +1.4802237e-01   +8.4395484e-04  +1.4802237e-01   +1.2814033e-01   +1.6649058e-01   +1.4802237e-01   +1.3361386e-01   +1.6282271e-01   +1.4288429e-01   +1.2814033e-01   +1.5866014e-01   +1.5277187e-01   +1.3877437e-01   +1.6649058e-01
  +4.5000000e+01   +5.0000000e+01   +1.2558755e-01   +6.5481171e-04  +1.2558755e-01   +1.0858180e-01   +1.4137184e-01   +1.2558755e-01   +1.1333715e-01   +1.3819050e-01   +1.2111361e-01   +1.0858180e-01   +1.3454696e-01   +1.2970701e-01   +1.1781309e-01   +1.4137184e-01
  +5.0000000e+01   +5.5000000e+01   +1.0840208e-01   +7.0604112e-04  +1.0840208e-01   +9.3711521e-02   +1.2199669e-01   +1.0840208e-01   +9.7849521e-02   +1.1924213e-01   +1.0452049e-01   +9.3711521e-02   +1.1610284e-01   +1.1199350e-01   +1.0176238e-01   +1.2199669e-01
  +5.5000000e+01   +6.0000000e+01   +9.3898199e-02   +7.5327456e-04  +9.3898199e-02   +8.1224512e-02   +1.0550969e-01   +9.3898199e-02   +8.4838417e-02   +1.0314374e-01   +9.0525629e-02   +8.1224512e-02   +1.0044919e-01   +9.7027647e-02   +8.8267193e-02   +1.0550969e-01
  +6.0000000e+01   +6.5000000e+01   +8.3401859e-02   +6.1759318e-04  +8.3401859e-02   +7.1896650e-02   +9.4129176e-02   +8.3401859e-02   +7.5171449e-02   +9.1940710e-02   +8.0324257e-02   +7.1896650e-02   +8.9440785e-02   +8.6253655e-02   +7.8277815e-02   +9.4129176e-02
  +6.5000000e+01   +7.0000000e+01   +7.2743083e-02   +6.1238245e-04  +7.2743083e-02   +6.2901003e-02   +8.1650239e-02   +7.2743083e-02   +6.5718820e-02   +7.9915667e-02   +7.0140813e-02   +6.2901003e-02   +7.7888820e-02   +7.5137282e-02   +6.8384385e-02   +8.1650239e-02
  +7.0000000e+01   +7.5000000e+01   +6.4939065e-02   +5.1022961e-04  +6.4939065e-02   +5.6053944e-02   +7.3035789e-02   +6.4939065e-02   +5.8647938e-02   +7.1378584e-02   +6.2523133e-02   +5.6053944e-02   +6.9457591e-02   +6.7170132e-02   +6.1108981e-02   +7.3035789e-02
  +7.5000000e+01   +8.0000000e+01   +5.7233052e-02   +4.5075661e-04  +5.7233052e-02   +4.9567817e-02   +6.4004313e-02   +5.7233052e-02   +5.1809592e-02   +6.2692534e-02   +5.5187311e-02   +4.9567817e-02   +6.1146546e-02   +5.9117862e-02   +5.3938187e-02   +6.4004313e-02
  +8.0000000e+01   +8.5000000e+01   +5.2537353e-02   +4.2650888e-04  +5.2537353e-02   +4.5244562e-02   +5.9243891e-02   +5.2537353e-02   +4.7373712e-02   +5.7878911e-02   +5.0556694e-02   +4.5244562e-02   +5.6308480e-02   +5.4379922e-02   +4.9408259e-02   +5.9243891e-02
  +8.5000000e+01   +9.0000000e+01   +4.7011569e-02   +4.2535632e-04  +4.7011569e-02   +4.0564972e-02   +5.2816652e-02   +4.7011569e-02   +4.2480451e-02   +5.1631940e-02   +4.5244385e-02   +4.0564972e-02   +5.0258940e-02   +4.8654642e-02   +4.4312938e-02   +5.2816652e-02
  +9.0000000e+01   +9.5000000e+01   +4.2366657e-02   +3.7653164e-04  +4.2366657e-02   +3.6573903e-02   +4.7528406e-02   +4.2366657e-02   +3.8316773e-02   +4.6470727e-02   +4.0765126e-02   +3.6573903e-02   +4.5238829e-02   +4.3855453e-02   +3.9986061e-02   +4.7528406e-02
  +9.5000000e+01   +1.0000000e+02   +3.8080418e-02   +3.7222932e-04  +3.8080418e-02   +3.2901730e-02   +4.2627785e-02   +3.8080418e-02   +3.4482340e-02   +4.1694289e-02   +3.6634838e-02   +3.2901730e-02   +4.0595617e-02   +3.9420634e-02   +3.5995484e-02   +4.2627785e-02
  +1.0000000e+02   +1.0500000e+02   +3.5488017e-02   +5.4230891e-04  +3.5488017e-02   +3.0573739e-02   +3.9881705e-02   +3.5488017e-02   +3.2022676e-02   +3.9055860e-02   +3.4166150e-02   +3.0573739e-02   +3.8057621e-02   +3.6699628e-02   +3.3400758e-02   +3.9881705e-02
  +1.0500000e+02   +1.1000000e+02   +3.2161871e-02   +5.3114699e-04  +3.2161871e-02   +2.7650413e-02   +3.6264581e-02   +3.2161871e-02   +2.9069863e-02   +3.5308797e-02   +3.0835883e-02   +2.7650413e-02   +3.4246847e-02   +3.3426954e-02   +3.0456311e-02   +3.6264581e-02
  +1.1000000e+02   +1.1500000e+02   +2.9180458e-02   +3.5907425e-04  +2.9180458e-02   +2.5175580e-02   +3.2682190e-02   +2.9180458e-02   +2.6413035e-02   +3.1968010e-02   +2.8054008e-02   +2.5175580e-02   +3.1122202e-02   +3.0227368e-02   +2.7603462e-02   +3.2682190e-02
  +1.1500000e+02   +1.2000000e+02   +2.8141185e-02   +8.6958985e-04  +2.8141185e-02   +2.4010492e-02   +3.2083030e-02   +2.8141185e-02   +2.5281201e-02   +3.1170095e-02   +2.6933963e-02   +2.4010492e-02   +3.0164754e-02   +2.9297868e-02   +2.6526948e-02   +3.2083030e-02
  +1.2000000e+02   +1.2500000e+02   +2.4485867e-02   +4.0439904e-04  +2.4485867e-02   +2.1139552e-02   +2.7355935e-02   +2.4485867e-02   +2.2178522e-02   +2.6798492e-02   +2.3549991e-02   +2.1139552e-02   +2.6115149e-02   +2.5346493e-02   +2.3173761e-02   +2.7355935e-02
  +1.2500000e+02   +1.3000000e+02   +2.2463957e-02   +2.7231974e-04  +2.2463957e-02   +1.9353292e-02   +2.5174382e-02   +2.2463957e-02   +2.0339689e-02   +2.4598895e-02   +2.1564733e-02   +1.9353292e-02   +2.3921132e-02   +2.3305192e-02   +2.1295334e-02   +2.5174382e-02
  +1.3000000e+02   +1.3500000e+02   +2.1297425e-02   +3.7929568e-04  +2.1297425e-02   +1.8239722e-02   +2.4090056e-02   +2.1297425e-02   +1.9204050e-02   +2.3463045e-02   +2.0401851e-02   +1.8239722e-02   +2.2755835e-02   +2.2151391e-02   +2.0150134e-02   +2.4090056e-02
  +1.3500000e+02   +1.4000000e+02   +1.8945659e-02   +2.7822487e-04  +1.8945659e-02   +1.6357342e-02   +2.1132681e-02   +1.8945659e-02   +1.7195223e-02   +2.0672918e-02   +1.8190846e-02   +1.6357342e-02   +2.0121636e-02   +1.9650117e-02   +1.8007921e-02   +2.1132681e-02
  +1.4000000e+02   +1.4500000e+02   +1.7455566e-02   +3.7312107e-04  +1.7455566e-02   +1.5098230e-02   +1.9393428e-02   +1.7455566e-02   +1.5859619e-02   +1.9017004e-02   +1.6779694e-02   +1.5098230e-02   +1.8543225e-02   +1.8076883e-02   +1.6592718e-02   +1.9393428e-02
  +1.4500000e+02   +1.5000000e+02   +1.6591402e-02   +3.6906014e-04  +1.6591402e-02   +1.4262782e-02   +1.8620372e-02   +1.6591402e-02   +1.5004244e-02   +1.8200694e-02   +1.5920212e-02   +1.4262782e-02   +1.7704080e-02   +1.7222760e-02   +1.5728333e-02   +1.8620372e-02
  +1.5000000e+02   +1.5500000e+02   +1.5608381e-02   +2.3943446e-04  +1.5608381e-02   +1.3349907e-02   +1.7652514e-02   +1.5608381e-02   +1.4081063e-02   +1.7183272e-02   +1.4930337e-02   +1.3349907e-02   +1.6649726e-02   +1.6255679e-02   +1.4801569e-02   +1.7652514e-02
  +1.5500000e+02   +1.6000000e+02   +1.4088011e-02   +2.5057555e-04  +1.4088011e-02   +1.2133910e-02   +1.5740937e-02   +1.4088011e-02   +1.2774003e-02   +1.5394471e-02   +1.3512759e-02   +1.2133910e-02   +1.4977066e-02   +1.4627295e-02   +1.3399233e-02   +1.5740937e-02
  +1.6000000e+02   +1.6500000e+02   +1.3108379e-02   +2.3757009e-04  +1.3108379e-02   +1.1288072e-02   +1.4650417e-02   +1.3108379e-02   +1.1883649e-02   +1.4327724e-02   +1.2574528e-02   +1.1288072e-02   +1.3943125e-02   +1.3612530e-02   +1.2468821e-02   +1.4650417e-02
  +1.6500000e+02   +1.7000000e+02   +1.1823386e-02   +2.1796357e-04  +1.1823386e-02   +1.0207554e-02   +1.3139366e-02   +1.1823386e-02   +1.0764569e-02   +1.2841474e-02   +1.1323789e-02   +1.0207554e-02   +1.2480936e-02   +1.2292270e-02   +1.1310056e-02   +1.3139366e-02
  +1.7000000e+02   +1.7500000e+02   +1.1333507e-02   +2.4977663e-04  +1.1333507e-02   +9.7476328e-03   +1.2672709e-02   +1.1333507e-02   +1.0265405e-02   +1.2404151e-02   +1.0871412e-02   +9.7476328e-03   +1.2075251e-02   +1.1766400e-02   +1.0772649e-02   +1.2672709e-02
  +1.7500000e+02   +1.8000000e+02   +1.0375195e-02   +2.5805413e-04  +1.0375195e-02   +8.9461362e-03   +1.1542482e-02   +1.0375195e-02   +9.4289318e-03   +1.1299123e-02   +9.9462273e-03   +8.9461362e-03   +1.0997623e-02   +1.0776238e-02   +9.9019087e-03   +1.1542482e-02
  +1.8000000e+02   +1.8500000e+02   +9.6485002e-03   +1.9966547e-04  +9.6485002e-03   +8.2920641e-03   +1.0793289e-02   +9.6485002e-03   +8.7660799e-03   +1.0512055e-02   +9.2177063e-03   +8.2920641e-03   +1.0189958e-02   +1.0063562e-02   +9.2391111e-03   +1.0793289e-02
  +1.8500000e+02   +1.9000000e+02   +9.0380750e-03   +2.0553048e-04  +9.0380750e-03   +7.7826900e-03   +1.0068849e-02   +9.0380750e-03   +8.2063041e-03   +9.8562239e-03   +8.6619014e-03   +7.7826900e-03   +9.5922510e-03   +9.3904680e-03   +8.6225850e-03   +1.0068849e-02
  +1.9000000e+02   +1.9500000e+02   +8.3401407e-03   +1.9288446e-04  +8.3401407e-03   +7.1952027e-03   +9.2509912e-03   +8.3401407e-03   +7.5881358e-03   +9.0674224e-03   +7.9942179e-03   +7.1952027e-03   +8.8307112e-03   +8.6606444e-03   +7.9722529e-03   +9.2509912e-03
  +1.9500000e+02   +2.0000000e+02   +7.7369216e-03   +2.1224759e-04  +7.7369216e-03   +6.6792808e-03   +8.5718828e-03   +7.7369216e-03   +7.0501105e-03   +8.3923463e-03   +7.4100938e-03   +6.6792808e-03   +8.1679650e-03   +8.0440153e-03   +7.4161084e-03   +8.5718828e-03
  +2.0000000e+02   +2.0500000e+02   +7.3342825e-03   +2.1341816e-04  +7.3342825e-03   +6.3108220e-03   +8.1678693e-03   +7.3342825e-03   +6.6649232e-03   +7.9881979e-03   +7.0196491e-03   +6.3108220e-03   +7.7670335e-03   +7.6312824e-03   +7.0155050e-03   +8.1678693e-03
  +2.0500000e+02   +2.1000000e+02   +6.5309586e-03   +1.8226893e-04  +6.5309586e-03   +5.6802544e-03   +7.1317379e-03   +6.5309586e-03   +5.9816459e-03   +7.0299615e-03   +6.2769236e-03   +5.6802544e-03   +6.8790008e-03   +6.7613549e-03   +6.2747062e-03   +7.1317379e-03
  +2.1000000e+02   +2.1500000e+02   +6.6544422e-03   +1.9475372e-04  +6.6544422e-03   +5.6827670e-03   +7.5000462e-03   +6.6544422e-03   +6.0036294e-03   +7.3252711e-03   +6.3658260e-03   +5.6827670e-03   +7.1153409e-03   +6.9292921e-03   +6.3232935e-03   +7.5000462e-03
  +2.1500000e+02   +2.2000000e+02   +5.9293556e-03   +1.7238816e-04  +5.9293556e-03   +5.1003853e-03   +6.5936591e-03   +5.9293556e-03   +5.3968436e-03   +6.4426322e-03   +5.6655843e-03   +5.1003853e-03   +6.2565065e-03   +6.1791029e-03   +5.6916815e-03   +6.5936591e-03
  +2.2000000e+02   +2.2500000e+02   +5.4298410e-03   +1.3459154e-04  +5.4298410e-03   +4.6981428e-03   +5.9798237e-03   +5.4298410e-03   +4.9607555e-03   +5.8667882e-03   +5.2032376e-03   +4.6981428e-03   +5.7210194e-03   +5.6437188e-03   +5.2223937e-03   +5.9798237e-03
  +2.2500000e+02   +2.3000000e+02   +5.1526915e-03   +1.4346413e-04  +5.1526915e-03   +4.4321247e-03   +5.7308274e-03   +5.1526915e-03   +4.6956732e-03   +5.5885027e-03   +4.9182271e-03   +4.4321247e-03   +5.4231073e-03   +5.3812740e-03   +4.9632262e-03   +5.7308274e-03
  +2.3000000e+02   +2.3500000e+02   +5.0873324e-03   +1.5021935e-04  +5.0873324e-03   +4.3414080e-03   +5.7230113e-03   +5.0873324e-03   +4.5932509e-03   +5.5940050e-03   +4.8610283e-03   +4.3414080e-03   +5.4298512e-03   +5.2982758e-03   +4.8414678e-03   +5.7230113e-03
  +2.3500000e+02   +2.4000000e+02   +4.5599708e-03   +1.4474169e-04  +4.5599708e-03   +3.9229330e-03   +5.0640593e-03   +4.5599708e-03   +4.1519837e-03   +4.9519596e-03   +4.3582854e-03   +3.9229330e-03   +4.8138709e-03   +4.7527925e-03   +4.3821443e-03   +5.0640593e-03
  +2.4000000e+02   +2.4500000e+02   +4.0726488e-03   +1.6847052e-04  +4.0726488e-03   +3.5379580e-03   +4.4397602e-03   +4.0726488e-03   +3.7351448e-03   +4.3748346e-03   +3.9071522e-03   +3.5379580e-03   +4.2779834e-03   +4.2246123e-03   +3.9296766e-03   +4.4397602e-03
  +2.4500000e+02   +2.5000000e+02   +4.0989120e-03   +1.8072890e-04  +4.0989120e-03   +3.5088387e-03   +4.5838978e-03   +4.0989120e-03   +3.7200526e-03   +4.4728747e-03   +3.9100572e-03   +3.5088387e-03   +4.3377509e-03   +4.2796892e-03   +3.9325117e-03   +4.5838978e-03
  +2.5000000e+02   +2.5500000e+02   +3.9102659e-03   +1.6483237e-04  +3.9102659e-03   +3.3339358e-03   +4.3999321e-03   +3.9102659e-03   +3.5386052e-03   +4.2852640e-03   +3.7251779e-03   +3.3339358e-03   +4.1486793e-03   +4.0885563e-03   +3.7452401e-03   +4.3999321e-03
  +2.5500000e+02   +2.6000000e+02   +3.4107605e-03   +1.1041811e-04  +3.4107605e-03   +2.9529977e-03   +3.7375564e-03   +3.4107605e-03   +3.1230557e-03   +3.6728401e-03   +3.2662130e-03   +2.9529977e-03   +3.5843830e-03   +3.5467535e-03   +3.2934313e-03   +3.7375564e-03
  +2.6000000e+02   +2.6500000e+02   +3.2095748e-03   +1.0794036e-04  +3.2095748e-03   +2.7853261e-03   +3.4982683e-03   +3.2095748e-03   +2.9385874e-03   +3.4566452e-03   +3.0829768e-03   +2.7853261e-03   +3.3868747e-03   +3.3240945e-03   +3.0892853e-03   +3.4982683e-03
  +2.6500000e+02   +2.7000000e+02   +3.0797597e-03   +1.4499401e-04  +3.0797597e-03   +2.6496056e-03   +3.4081080e-03   +3.0797597e-03   +2.8091698e-03   +3.3356634e-03   +2.9403734e-03   +2.6496056e-03   +3.2424854e-03   +3.2121808e-03   +2.9697817e-03   +3.4081080e-03
  +2.7000000e+02   +2.7500000e+02   +3.0152958e-03   +1.3400856e-04  +3.0152958e-03   +2.5748613e-03   +3.3803656e-03   +3.0152958e-03   +2.7311661e-03   +3.3000696e-03   +2.8762403e-03   +2.5748613e-03   +3.2019785e-03   +3.1492414e-03   +2.8896898e-03   +3.3803656e-03
  +2.7500000e+02   +2.8000000e+02   +2.5319738e-03   +9.9641157e-05  +2.5319738e-03   +2.2068638e-03   +2.7359688e-03   +2.5319738e-03   +2.3386195e-03   +2.6904844e-03   +2.4218318e-03   +2.2068638e-03   +2.6269592e-03   +2.6365561e-03   +2.4720474e-03   +2.7359688e-03
  +2.8000000e+02   +2.8500000e+02   +2.6529784e-03   +1.0623288e-04  +2.6529784e-03   +2.2615586e-03   +2.9795840e-03   +2.6529784e-03   +2.4034160e-03   +2.9027750e-03   +2.5256552e-03   +2.2615586e-03   +2.8107136e-03   +2.7764333e-03   +2.5479507e-03   +2.9795840e-03
  +2.8500000e+02   +2.9000000e+02   +2.1489778e-03   +1.6811391e-04  +2.1489778e-03   +1.8935207e-03   +2.2679597e-03   +2.1489778e-03   +1.9987609e-03   +2.2587561e-03   +2.0668469e-03   +1.8935207e-03   +2.2238493e-03   +2.2192902e-03   +2.1004601e-03   +2.2679597e-03
  +2.9000000e+02   +2.9500000e+02   +2.3295454e-03   +1.7982884e-04  +2.3295454e-03   +2.0004749e-03   +2.5827644e-03   +2.3295454e-03   +2.1150310e-03   +2.5406479e-03   +2.2318750e-03   +2.0004749e-03   +2.4802523e-03   +2.4205662e-03   +2.2297376e-03   +2.5827644e-03
  +2.9500000e+02   +3.0000000e+02   +2.2179982e-03   +1.3334978e-04  +2.2179982e-03   +1.8871198e-03   +2.4953630e-03   +2.2179982e-03   +2.0087853e-03   +2.4278519e-03   +2.1082311e-03   +1.8871198e-03   +2.3473586e-03   +2.3247905e-03   +2.1332127e-03   +2.4953630e-03
  +3.0000000e+02   +3.0500000e+02   +1.8621398e-03   +8.2520291e-05  +1.8621398e-03   +1.6336313e-03   +1.9815582e-03   +1.8621398e-03   +1.7214235e-03   +1.9760677e-03   +1.7945793e-03   +1.6336313e-03   +1.9495322e-03   +1.9205806e-03   +1.8071634e-03   +1.9815582e-03
  +3.0500000e+02   +3.1000000e+02   +1.9859365e-03   +1.5875672e-04  +1.9859365e-03   +1.6839632e-03   +2.2443958e-03   +1.9859365e-03   +1.7966241e-03   +2.1773775e-03   +1.8831909e-03   +1.6839632e-03   +2.0998559e-03   +2.0871302e-03   +1.9128198e-03   +2.2443958e-03
  +3.1000000e+02   +3.1500000e+02   +1.7889944e-03   +1.0636413e-04  +1.7889944e-03   +1.5350489e-03   +1.9831561e-03   +1.7889944e-03   +1.6295827e-03   +1.9416192e-03   +1.7072248e-03   +1.5350489e-03   +1.8886092e-03   +1.8681597e-03   +1.7265957e-03   +1.9831561e-03
  +3.1500000e+02   +3.2000000e+02   +1.5253831e-03   +8.9595194e-05  +1.5253831e-03   +1.3376279e-03   +1.6203112e-03   +1.5253831e-03   +1.4168624e-03   +1.6066815e-03   +1.4620750e-03   +1.3376279e-03   +1.5764119e-03   +1.5821404e-03   +1.4954769e-03   +1.6203112e-03
  +3.2000000e+02   +3.2500000e+02   +1.6442884e-03   +7.6498207e-05  +1.6442884e-03   +1.3990049e-03   +1.8452468e-03   +1.6442884e-03   +1.4889740e-03   +1.8002470e-03   +1.5641318e-03   +1.3990049e-03   +1.7434706e-03   +1.7211954e-03   +1.5806116e-03   +1.8452468e-03
  +3.2500000e+02   +3.3000000e+02   +1.3814445e-03   +7.1881994e-05  +1.3814445e-03   +1.2130043e-03   +1.4632726e-03   +1.3814445e-03   +1.2809031e-03   +1.4591009e-03   +1.3290141e-03   +1.2130043e-03   +1.4380934e-03   +1.4267812e-03   +1.3474667e-03   +1.4632726e-03
  +3.3000000e+02   +3.3500000e+02   +1.3199858e-03   +7.6699463e-05  +1.3199858e-03   +1.1511441e-03   +1.4140424e-03   +1.3199858e-03   +1.2180561e-03   +1.4046338e-03   +1.2665609e-03   +1.1511441e-03   +1.3791734e-03   +1.3666262e-03   +1.2836629e-03   +1.4140424e-03
  +3.3500000e+02   +3.4000000e+02   +1.3761107e-03   +8.2568238e-05  +1.3761107e-03   +1.1708486e-03   +1.5439491e-03   +1.3761107e-03   +1.2482902e-03   +1.5027781e-03   +1.3070358e-03   +1.1708486e-03   +1.4536924e-03   +1.4442301e-03   +1.3287206e-03   +1.5439491e-03
  +3.4000000e+02   +3.4500000e+02   +1.2150757e-03   +1.1581445e-04  +1.2150757e-03   +1.0467050e-03   +1.3296661e-03   +1.2150757e-03   +1.1151647e-03   +1.3038368e-03   +1.1561193e-03   +1.0467050e-03   +1.2661514e-03   +1.2699484e-03   +1.1840699e-03   +1.3296661e-03
  +3.4500000e+02   +3.5000000e+02   +1.0284697e-03   +1.2609003e-04  +1.0284697e-03   +9.1866003e-04   +1.0725509e-03   +1.0284697e-03   +9.6235917e-04   +1.0707045e-03   +1.0003989e-03   +9.1866003e-04   +1.0725509e-03   +1.0496690e-03   +1.0042315e-03   +1.0535002e-03
  +3.5000000e+02   +3.5500000e+02   +1.2231007e-03   +7.9303007e-05  +1.2231007e-03   +1.0267957e-03   +1.4010370e-03   +1.2231007e-03   +1.1015502e-03   +1.3498398e-03   +1.1529753e-03   +1.0267957e-03   +1.2931191e-03   +1.2935211e-03   +1.1798522e-03   +1.4010370e-03
  +3.5500000e+02   +3.6000000e+02   +1.0203761e-03   +6.0561882e-05  +1.0203761e-03   +8.7837169e-04   +1.1205260e-03   +1.0203761e-03   +9.3625392e-04   +1.0953077e-03   +9.7054584e-04   +8.7837169e-04   +1.0634908e-03   +1.0694352e-03   +9.9666625e-04   +1.1205260e-03
  +3.6000000e+02   +3.6500000e+02   +9.7374089e-04   +7.3357124e-05  +9.7374089e-04   +8.4320681e-04   +1.0549712e-03   +9.7374089e-04   +8.9274481e-04   +1.0465285e-03   +9.3384213e-04   +8.4320681e-04   +1.0267421e-03   +1.0092726e-03   +9.4201017e-04   +1.0549712e-03
  +3.6500000e+02   +3.7000000e+02   +9.1700201e-04   +6.9271380e-05  +9.1700201e-04   +7.9073578e-04   +1.0008181e-03   +9.1700201e-04   +8.4342045e-04   +9.8074516e-04   +8.7205752e-04   +7.9073578e-04   +9.5289929e-04   +9.6014228e-04   +8.9772413e-04   +1.0008181e-03
  +3.7000000e+02   +3.7500000e+02   +9.7644695e-04   +6.8996934e-05  +9.7644695e-04   +8.2596661e-04   +1.1029561e-03   +9.7644695e-04   +8.7851753e-04   +1.0792158e-03   +9.3008546e-04   +8.2596661e-04   +1.0472932e-03   +1.0204546e-03   +9.3208181e-04   +1.1029561e-03
  +3.7500000e+02   +3.8000000e+02   +7.8954179e-04   +7.1105735e-05  +7.8954179e-04   +6.8881306e-04   +8.4229947e-04   +7.8954179e-04   +7.2895513e-04   +8.3949252e-04   +7.5826392e-04   +6.8881306e-04   +8.2631229e-04   +8.1669527e-04   +7.6864912e-04   +8.4229947e-04
  +3.8000000e+02   +3.8500000e+02   +8.0844249e-04   +7.8355385e-05  +8.0844249e-04   +6.9578431e-04   +8.8368880e-04   +8.0844249e-04   +7.3662044e-04   +8.7702890e-04   +7.7530331e-04   +6.9578431e-04   +8.6004350e-04   +8.3736092e-04   +7.7683498e-04   +8.8368880e-04
  +3.8500000e+02   +3.9000000e+02   +6.9008523e-04   +7.4529624e-05  +6.9008523e-04   +6.0056333e-04   +7.3885559e-04   +6.9008523e-04   +6.4082778e-04   +7.2715432e-04   +6.5667967e-04   +6.0056333e-04   +7.0842867e-04   +7.2158784e-04   +6.8214131e-04   +7.3885559e-04
  +3.9000000e+02   +3.9500000e+02   +6.8132119e-04   +5.7790173e-05  +6.8132119e-04   +5.9101914e-04   +7.3423746e-04   +6.8132119e-04   +6.2632478e-04   +7.2926293e-04   +6.5326073e-04   +5.9101914e-04   +7.1617791e-04   +7.0652369e-04   +6.6183730e-04   +7.3423746e-04
  +3.9500000e+02   +4.0000000e+02   +6.7718969e-04   +5.1581229e-05  +6.7718969e-04   +5.8465973e-04   +7.3717945e-04   +6.7718969e-04   +6.1987890e-04   +7.2955999e-04   +6.4859415e-04   +5.8465973e-04   +7.1487149e-04   +7.0358899e-04   +6.5572638e-04   +7.3717945e-04
  +4.0000000e+02   +4.0500000e+02   +6.0587218e-04   +4.6338582e-05  +6.0587218e-04   +5.2728847e-04   +6.4860548e-04   +6.0587218e-04   +5.5934880e-04   +6.4425810e-04   +5.8045213e-04   +5.2728847e-04   +6.3254140e-04   +6.2872869e-04   +5.9164775e-04   +6.4860548e-04
  +4.0500000e+02   +4.1000000e+02   +5.1402244e-04   +5.5810038e-05  +5.1402244e-04   +4.5726243e-04   +5.3087917e-04   +5.1402244e-04   +4.8336638e-04   +5.3087917e-04   +4.9554347e-04   +4.5726243e-04   +5.2734850e-04   +5.2861075e-04   +5.0861236e-04   +5.2540030e-04
  +4.1000000e+02   +4.1500000e+02   +5.4839408e-04   +5.5153711e-05  +5.4839408e-04   +4.7733540e-04   +5.8637324e-04   +5.4839408e-04   +5.0582204e-04   +5.8396228e-04   +5.2611542e-04   +4.7733540e-04   +5.7438553e-04   +5.6804161e-04   +5.3433035e-04   +5.8637324e-04
  +4.1500000e+02   +4.2000000e+02   +5.4327915e-04   +4.6077695e-05  +5.4327915e-04   +4.6495771e-04   +5.9931053e-04   +5.4327915e-04   +4.9662369e-04   +5.8650113e-04   +5.1602654e-04   +4.6495771e-04   +5.6911748e-04   +5.6972308e-04   +5.2962709e-04   +5.9931053e-04
  +4.2000000e+02   +4.2500000e+02   +5.1942465e-04   +5.3414690e-05  +5.1942465e-04   +4.4320906e-04   +5.7458926e-04   +5.1942465e-04   +4.7467014e-04   +5.6101191e-04   +4.9186212e-04   +4.4320906e-04   +5.4242364e-04   +5.4577253e-04   +5.0709670e-04   +5.7458926e-04
  +4.2500000e+02   +4.3000000e+02   +5.5762653e-04   +1.6508570e-04  +5.5762653e-04   +4.7016920e-04   +6.3293230e-04   +5.5762653e-04   +5.0086239e-04   +6.1781070e-04   +5.3045745e-04   +4.7016920e-04   +5.9891926e-04   +5.8439450e-04   +5.3306087e-04   +6.3293230e-04
  +4.3000000e+02   +4.3500000e+02   +3.4976728e-04   +1.6210370e-04  +3.4976728e-04   +3.1724439e-04   +3.6291928e-04   +3.4976728e-04   +3.3831623e-04   +3.4446804e-04   +3.3449259e-04   +3.1724439e-04   +3.4064823e-04   +3.6291928e-04   +3.5949136e-04   +3.4235390e-04
  +4.3500000e+02   +4.4000000e+02   +4.0391790e-04   +4.4527931e-05  +4.0391790e-04   +3.5443452e-04   +4.2649698e-04   +4.0391790e-04   +3.7459187e-04   +4.2649698e-04   +3.8921870e-04   +3.5443452e-04   +4.2260620e-04   +4.1620183e-04   +3.9457370e-04   +4.2415821e-04
  +4.4000000e+02   +4.4500000e+02   +3.6682104e-04   +4.0954228e-05  +3.6682104e-04   +3.2814281e-04   +3.8333902e-04   +3.6682104e-04   +3.4387938e-04   +3.8074796e-04   +3.5742347e-04   +3.2814281e-04   +3.8333902e-04   +3.7195843e-04   +3.5797122e-04   +3.6954834e-04
  +4.4500000e+02   +4.5000000e+02   +4.2136070e-04   +6.0602588e-05  +4.2136070e-04   +3.5641162e-04   +4.7348635e-04   +4.2136070e-04   +3.8113510e-04   +4.6208401e-04   +3.9959948e-04   +3.5641162e-04   +4.4720175e-04   +4.4257647e-04   +4.0698196e-04   +4.7348635e-04
  +4.5000000e+02   +4.5500000e+02   +2.6775085e-04   +5.7646490e-05  +2.6775085e-04   +2.4455615e-04   +2.7087086e-04   +2.6775085e-04   +2.6108195e-04   +2.5995674e-04   +2.6284761e-04   +2.4960655e-04   +2.6712744e-04   +2.6850616e-04   +2.7087086e-04   +2.4455615e-04
  +4.5500000e+02   +4.6000000e+02   +3.7537472e-04   +3.7388388e-05  +3.7537472e-04   +3.1645089e-04   +4.2400628e-04   +3.7537472e-04   +3.3989304e-04   +4.1102296e-04   +3.5418634e-04   +3.1645089e-04   +3.9540877e-04   +3.9646568e-04   +3.6466320e-04   +4.2400628e-04
  +4.6000000e+02   +4.6500000e+02   +3.8308007e-04   +4.3880096e-05  +3.8308007e-04   +3.1987697e-04   +4.4022915e-04   +3.8308007e-04   +3.4435093e-04   +4.2394986e-04   +3.6026525e-04   +3.1987697e-04   +4.0576897e-04   +4.0629670e-04   +3.7050246e-04   +4.4022915e-04
  +4.6500000e+02   +4.7000000e+02   +3.1606994e-04   +3.9414274e-05  +3.1606994e-04   +2.7415073e-04   +3.3823247e-04   +3.1606994e-04   +2.9060045e-04   +3.3823247e-04   +3.0333128e-04   +2.7415073e-04   +3.3304429e-04   +3.2654729e-04   +3.0663272e-04   +3.3803929e-04
  +4.7000000e+02   +4.7500000e+02   +3.2839435e-04   +3.6246901e-05  +3.2839435e-04   +2.7897832e-04   +3.6594617e-04   +3.2839435e-04   +2.9853135e-04   +3.5748160e-04   +3.1123663e-04   +2.7897832e-04   +3.4585408e-04   +3.4488093e-04   +3.1883897e-04   +3.6594617e-04
  +4.7500000e+02   +4.8000000e+02   +2.4899578e-04   +3.9727990e-05  +2.4899578e-04   +2.2492866e-04   +2.5947516e-04   +2.4899578e-04   +2.3507801e-04   +2.5549980e-04   +2.4378544e-04   +2.2492866e-04   +2.5947516e-04   +2.5088220e-04   +2.4391323e-04   +2.4486288e-04
  +4.8000000e+02   +4.8500000e+02   +3.0147582e-04   +3.6009373e-05  +3.0147582e-04   +2.5310485e-04   +3.4287935e-04   +3.0147582e-04   +2.7422969e-04   +3.2787756e-04   +2.8168472e-04   +2.5310485e-04   +3.1191822e-04   +3.2215326e-04   +2.9723888e-04   +3.4287935e-04
  +4.8500000e+02   +4.9000000e+02   +2.2047997e-04   +3.4245687e-05  +2.2047997e-04   +1.9783262e-04   +2.2613294e-04   +2.2047997e-04   +2.0920606e-04   +2.2436788e-04   +2.1294893e-04   +1.9783262e-04   +2.2423864e-04   +2.2613294e-04   +2.2027513e-04   +2.1995228e-04
  +4.9000000e+02   +4.9500000e+02   +2.7864794e-04   +3.5063086e-05  +2.7864794e-04   +2.3316387e-04   +3.1832582e-04   +2.7864794e-04   +2.5093926e-04   +3.0755192e-04   +2.6230051e-04   +2.3316387e-04   +2.9495105e-04   +2.9506073e-04   +2.6983903e-04   +3.1832582e-04
  +4.9500000e+02   +5.0000000e+02   +2.0307631e-04   +3.3183001e-05  +2.0307631e-04   +1.8025894e-04   +2.1063546e-04   +2.0307631e-04   +1.9185033e-04   +2.0815798e-04   +1.9461534e-04   +1.8025894e-04   +2.0589815e-04   +2.1063546e-04   +2.0380510e-04   +2.0732727e-04
<\histogram>


<histogram> 101 "pta SM LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   +2.3564634e-01   +3.0513047e-04  +2.3564634e-01   +1.8286765e-01   +3.0751131e-01   +2.3564634e-01   +1.9632302e-01   +2.8833098e-01   +2.1949587e-01   +1.8286765e-01   +2.6856968e-01   +2.5132199e-01   +2.0938280e-01   +3.0751131e-01
  +2.5000000e+01   +3.0000000e+01   +1.8122869e-01   +2.3170424e-04  +1.8122869e-01   +1.4039539e-01   +2.3692366e-01   +1.8122869e-01   +1.5098627e-01   +2.2174692e-01   +1.6851646e-01   +1.4039539e-01   +2.0619255e-01   +1.9363230e-01   +1.6132004e-01   +2.3692366e-01
  +3.0000000e+01   +3.5000000e+01   +1.4450343e-01   +1.8524956e-04  +1.4450343e-01   +1.1176429e-01   +1.8922989e-01   +1.4450343e-01   +1.2038952e-01   +1.7681080e-01   +1.3415058e-01   +1.1176429e-01   +1.6414330e-01   +1.5465327e-01   +1.2884561e-01   +1.8922989e-01
  +3.5000000e+01   +4.0000000e+01   +1.1882807e-01   +1.5375936e-04  +1.1882807e-01   +9.1773029e-02   +1.5584335e-01   +1.1882807e-01   +9.8998721e-02   +1.4539506e-01   +1.1015508e-01   +9.1773029e-02   +1.3478301e-01   +1.2736721e-01   +1.0611290e-01   +1.5584335e-01
  +4.0000000e+01   +4.5000000e+01   +9.9183305e-02   +1.3076468e-04  +9.9183305e-02   +7.6483985e-02   +1.3028674e-01   +9.9183305e-02   +8.2632158e-02   +1.2135822e-01   +9.1803658e-02   +7.6483985e-02   +1.1232867e-01   +1.0648038e-01   +8.8711538e-02   +1.3028674e-01
  +4.5000000e+01   +5.0000000e+01   +8.4029636e-02   +1.1290907e-04  +8.4029636e-02   +6.4695693e-02   +1.1056062e-01   +8.4029636e-02   +7.0007252e-02   +1.0281657e-01   +7.7654180e-02   +6.4695693e-02   +9.5015718e-02   +9.0358679e-02   +7.5280136e-02   +1.1056062e-01
  +5.0000000e+01   +5.5000000e+01   +7.2632224e-02   +1.0011921e-04  +7.2632224e-02   +5.5855347e-02   +9.5683791e-02   +7.2632224e-02   +6.0511770e-02   +8.8870974e-02   +6.7043123e-02   +5.5855347e-02   +8.2032292e-02   +7.8200183e-02   +6.5150579e-02   +9.5683791e-02
  +5.5000000e+01   +6.0000000e+01   +6.3295273e-02   +8.8231096e-05  +6.3295273e-02   +4.8614317e-02   +8.3493897e-02   +6.3295273e-02   +5.2732921e-02   +7.7446518e-02   +5.8351720e-02   +4.8614317e-02   +7.1397712e-02   +6.8237659e-02   +5.6850548e-02   +8.3493897e-02
  +6.0000000e+01   +6.5000000e+01   +5.5355938e-02   +7.8966712e-05  +5.5355938e-02   +4.2450462e-02   +7.3138061e-02   +5.5355938e-02   +4.6118455e-02   +6.7732145e-02   +5.0953249e-02   +4.2450462e-02   +6.2345122e-02   +5.9774072e-02   +4.9799314e-02   +7.3138061e-02
  +6.5000000e+01   +7.0000000e+01   +4.9008401e-02   +7.2410538e-05  +4.9008401e-02   +3.7537124e-02   +6.4836246e-02   +4.9008401e-02   +4.0830157e-02   +5.9965458e-02   +4.5055775e-02   +3.7537124e-02   +5.5129123e-02   +5.2989186e-02   +4.4146653e-02   +6.4836246e-02
  +7.0000000e+01   +7.5000000e+01   +4.3654397e-02   +6.5597226e-05  +4.3654397e-02   +3.3398891e-02   +5.7821306e-02   +4.3654397e-02   +3.6369597e-02   +5.3414427e-02   +4.0088655e-02   +3.3398891e-02   +4.9051479e-02   +4.7256035e-02   +3.9370215e-02   +5.7821306e-02
  +7.5000000e+01   +8.0000000e+01   +3.9044825e-02   +5.9968169e-05  +3.9044825e-02   +2.9836118e-02   +5.1782597e-02   +3.9044825e-02   +3.2529244e-02   +4.7774274e-02   +3.5812261e-02   +2.9836118e-02   +4.3818988e-02   +4.2320736e-02   +3.5258494e-02   +5.1782597e-02
  +8.0000000e+01   +8.5000000e+01   +3.4969047e-02   +5.4877810e-05  +3.4969047e-02   +2.6688294e-02   +4.6438574e-02   +3.4969047e-02   +2.9133611e-02   +4.2787256e-02   +3.2033933e-02   +2.6688294e-02   +3.9195922e-02   +3.7953190e-02   +3.1619774e-02   +4.6438574e-02
  +8.5000000e+01   +9.0000000e+01   +3.1712362e-02   +5.0827558e-05  +3.1712362e-02   +2.4176579e-02   +4.2163502e-02   +3.1712362e-02   +2.6420379e-02   +3.8802450e-02   +2.9019125e-02   +2.4176579e-02   +3.5507080e-02   +3.4459267e-02   +2.8708905e-02   +4.2163502e-02
  +9.0000000e+01   +9.5000000e+01   +2.8737118e-02   +4.7851520e-05  +2.8737118e-02   +2.1880867e-02   +3.8258568e-02   +2.8737118e-02   +2.3941632e-02   +3.5162023e-02   +2.6263586e-02   +2.1880867e-02   +3.2135468e-02   +3.1267857e-02   +2.6050055e-02   +3.8258568e-02
  +9.5000000e+01   +1.0000000e+02   +2.6028021e-02   +4.4561981e-05  +2.6028021e-02   +1.9795287e-02   +3.4695001e-02   +2.6028021e-02   +2.1684615e-02   +3.1847236e-02   +2.3760266e-02   +1.9795287e-02   +2.9072467e-02   +2.8355437e-02   +2.3623642e-02   +3.4695001e-02
  +1.0000000e+02   +1.0500000e+02   +2.3727441e-02   +4.1255214e-05  +2.3727441e-02   +1.8026180e-02   +3.1665246e-02   +2.3727441e-02   +1.9767939e-02   +2.9032304e-02   +2.1636807e-02   +1.8026180e-02   +2.6474259e-02   +2.5879284e-02   +2.1560695e-02   +3.1665246e-02
  +1.0500000e+02   +1.1000000e+02   +2.1732267e-02   +3.9237428e-05  +2.1732267e-02   +1.6491988e-02   +2.9037814e-02   +2.1732267e-02   +1.8105710e-02   +2.6591061e-02   +1.9795318e-02   +1.6491988e-02   +2.4221058e-02   +2.3731943e-02   +1.9771692e-02   +2.9037814e-02
  +1.1000000e+02   +1.1500000e+02   +1.9896512e-02   +3.6600202e-05  +1.9896512e-02   +1.5081343e-02   +2.6618004e-02   +1.9896512e-02   +1.6576295e-02   +2.4344874e-02   +1.8102122e-02   +1.5081343e-02   +2.2149305e-02   +2.1754288e-02   +1.8124056e-02   +2.6618004e-02
  +1.1500000e+02   +1.2000000e+02   +1.8287420e-02   +3.4822946e-05  +1.8287420e-02   +1.3851343e-02   +2.4486845e-02   +1.8287420e-02   +1.5235719e-02   +2.2376031e-02   +1.6625754e-02   +1.3851343e-02   +2.0342858e-02   +2.0012538e-02   +1.6672960e-02   +2.4486845e-02
  +1.2000000e+02   +1.2500000e+02   +1.6765491e-02   +3.2653781e-05  +1.6765491e-02   +1.2683070e-02   +2.2478409e-02   +1.6765491e-02   +1.3967761e-02   +2.0513837e-02   +1.5223477e-02   +1.2683070e-02   +1.8627066e-02   +1.8371090e-02   +1.5305427e-02   +2.2478409e-02
  +1.2500000e+02   +1.3000000e+02   +1.5345989e-02   +3.0606624e-05  +1.5345989e-02   +1.1597402e-02   +2.0598289e-02   +1.5345989e-02   +1.2785138e-02   +1.8776970e-02   +1.3920350e-02   +1.1597402e-02   +1.7032594e-02   +1.6834511e-02   +1.4025264e-02   +2.0598289e-02
  +1.3000000e+02   +1.3500000e+02   +1.4174908e-02   +2.8650815e-05  +1.4174908e-02   +1.0697480e-02   +1.9054457e-02   +1.4174908e-02   +1.1809480e-02   +1.7344064e-02   +1.2840176e-02   +1.0697480e-02   +1.5710919e-02   +1.5572774e-02   +1.2974078e-02   +1.9054457e-02
  +1.3500000e+02   +1.4000000e+02   +1.3136287e-02   +2.7694035e-05  +1.3136287e-02   +9.9083473e-03   +1.7670086e-02   +1.3136287e-02   +1.0944178e-02   +1.6073233e-02   +1.1892980e-02   +9.9083473e-03   +1.4551954e-02   +1.4441359e-02   +1.2031467e-02   +1.7670086e-02
  +1.4000000e+02   +1.4500000e+02   +1.2182335e-02   +2.6630617e-05  +1.2182335e-02   +9.1782036e-03   +1.6408448e-02   +1.2182335e-02   +1.0149416e-02   +1.4906002e-02   +1.1016590e-02   +9.1782036e-03   +1.3479624e-02   +1.3410251e-02   +1.1172424e-02   +1.6408448e-02
  +1.4500000e+02   +1.5000000e+02   +1.1248386e-02   +2.4635977e-05  +1.1248386e-02   +8.4643515e-03   +1.5169735e-02   +1.1248386e-02   +9.3713196e-03   +1.3763245e-02   +1.0159753e-02   +8.4643515e-03   +1.2431221e-02   +1.2397879e-02   +1.0328991e-02   +1.5169735e-02
  +1.5000000e+02   +1.5500000e+02   +1.0420814e-02   +2.4194445e-05  +1.0420814e-02   +7.8357916e-03   +1.4066609e-02   +1.0420814e-02   +8.6818483e-03   +1.2750649e-02   +9.4052933e-03   +7.8357916e-03   +1.1508082e-02   +1.1496319e-02   +9.5778783e-03   +1.4066609e-02
  +1.5500000e+02   +1.6000000e+02   +9.7098200e-03   +2.3272005e-05  +9.7098200e-03   +7.2940616e-03   +1.3120915e-02   +9.7098200e-03   +8.0895004e-03   +1.1880693e-02   +8.7550553e-03   +7.2940616e-03   +1.0712467e-02   +1.0723425e-02   +8.9339606e-03   +1.3120915e-02
  +1.6000000e+02   +1.6500000e+02   +9.0247661e-03   +2.1904516e-05  +9.0247661e-03   +6.7744830e-03   +1.2205789e-02   +9.0247661e-03   +7.5187643e-03   +1.1042478e-02   +8.1314061e-03   +6.7744830e-03   +9.9493855e-03   +9.9755127e-03   +8.3108560e-03   +1.2205789e-02
  +1.6500000e+02   +1.7000000e+02   +8.3561266e-03   +2.0873886e-05  +8.3561266e-03   +6.2665383e-03   +1.1313486e-02   +8.3561266e-03   +6.9617037e-03   +1.0224348e-02   +7.5217202e-03   +6.2665383e-03   +9.2033890e-03   +9.2462540e-03   +7.7032916e-03   +1.1313486e-02
  +1.7000000e+02   +1.7500000e+02   +7.7594674e-03   +1.9961231e-05  +7.7594674e-03   +5.8115501e-03   +1.0520747e-02   +7.7594674e-03   +6.4646119e-03   +9.4942910e-03   +6.9755980e-03   +5.8115501e-03   +8.5351675e-03   +8.5983669e-03   +7.1635203e-03   +1.0520747e-02
  +1.7500000e+02   +1.8000000e+02   +7.2518896e-03   +1.9035043e-05  +7.2518896e-03   +5.4270780e-03   +9.8417624e-03   +7.2518896e-03   +6.0417359e-03   +8.8732315e-03   +6.5141167e-03   +5.4270780e-03   +7.9705107e-03   +8.0434480e-03   +6.7012030e-03   +9.8417624e-03
  +1.8000000e+02   +1.8500000e+02   +6.7324839e-03   +1.8424869e-05  +6.7324839e-03   +5.0342645e-03   +9.1450916e-03   +6.7324839e-03   +5.6090054e-03   +8.2376995e-03   +6.0426232e-03   +5.0342645e-03   +7.3936032e-03   +7.4740750e-03   +6.2268438e-03   +9.1450916e-03
  +1.8500000e+02   +1.9000000e+02   +6.2821491e-03   +1.7447292e-05  +6.2821491e-03   +4.6938810e-03   +8.5416435e-03   +6.2821491e-03   +5.2338198e-03   +7.6866810e-03   +5.6340610e-03   +4.6938810e-03   +6.8936964e-03   +6.9808904e-03   +5.8159591e-03   +8.5416435e-03
  +1.9000000e+02   +1.9500000e+02   +5.8702266e-03   +1.6726570e-05  +5.8702266e-03   +4.3808306e-03   +7.9920016e-03   +5.8702266e-03   +4.8906366e-03   +7.1826627e-03   +5.2583068e-03   +4.3808306e-03   +6.4339332e-03   +6.5316806e-03   +5.4417108e-03   +7.9920016e-03
  +1.9500000e+02   +2.0000000e+02   +5.4965484e-03   +1.6116464e-05  +5.4965484e-03   +4.0992735e-03   +7.4891143e-03   +5.4965484e-03   +4.5793157e-03   +6.7254393e-03   +4.9203541e-03   +4.0992735e-03   +6.0204226e-03   +6.1206827e-03   +5.0992978e-03   +7.4891143e-03
  +2.0000000e+02   +2.0500000e+02   +5.1243248e-03   +1.5379158e-05  +5.1243248e-03   +3.8184289e-03   +6.9890478e-03   +5.1243248e-03   +4.2692062e-03   +6.2699959e-03   +4.5832565e-03   +3.8184289e-03   +5.6079586e-03   +5.7119891e-03   +4.7588048e-03   +6.9890478e-03
  +2.0500000e+02   +2.1000000e+02   +4.7832317e-03   +1.5043986e-05  +4.7832317e-03   +3.5602729e-03   +6.5320698e-03   +4.7832317e-03   +3.9850331e-03   +5.8526432e-03   +4.2733927e-03   +3.5602729e-03   +5.2288166e-03   +5.3385119e-03   +4.4476513e-03   +6.5320698e-03
  +2.1000000e+02   +2.1500000e+02   +4.4443912e-03   +1.4336069e-05  +4.4443912e-03   +3.3047690e-03   +6.0759648e-03   +4.4443912e-03   +3.7027363e-03   +5.4380461e-03   +3.9667114e-03   +3.3047690e-03   +4.8535691e-03   +4.9657473e-03   +4.1370916e-03   +6.0759648e-03
  +2.1500000e+02   +2.2000000e+02   +4.1833780e-03   +1.4659960e-05  +4.1833780e-03   +3.1089443e-03   +5.7235189e-03   +4.1833780e-03   +3.4852797e-03   +5.1186770e-03   +3.7316631e-03   +3.1089443e-03   +4.5659699e-03   +4.6777015e-03   +3.8971133e-03   +5.7235189e-03
  +2.2000000e+02   +2.2500000e+02   +3.9184213e-03   +1.3204776e-05  +3.9184213e-03   +2.9101123e-03   +5.3651794e-03   +3.9184213e-03   +3.2645375e-03   +4.7944827e-03   +3.4930048e-03   +2.9101123e-03   +4.2739537e-03   +4.3848389e-03   +3.6531217e-03   +5.3651794e-03
  +2.2500000e+02   +2.3000000e+02   +3.6624645e-03   +1.3662542e-05  +3.6624645e-03   +2.7166597e-03   +5.0216341e-03   +3.6624645e-03   +3.0512934e-03   +4.4813002e-03   +3.2608042e-03   +2.7166597e-03   +3.9898386e-03   +4.1040674e-03   +3.4192037e-03   +5.0216341e-03
  +2.3000000e+02   +2.3500000e+02   +3.4140748e-03   +1.2349909e-05  +3.4140748e-03   +2.5313083e-03   +4.6835672e-03   +3.4140748e-03   +2.8443535e-03   +4.1773765e-03   +3.0383272e-03   +2.5313083e-03   +3.7176214e-03   +3.8277726e-03   +3.1890158e-03   +4.6835672e-03
  +2.3500000e+02   +2.4000000e+02   +3.2244876e-03   +1.2262916e-05  +3.2244876e-03   +2.3893390e-03   +4.4269446e-03   +3.2244876e-03   +2.6864032e-03   +3.9454026e-03   +2.8679214e-03   +2.3893390e-03   +3.5091171e-03   +3.6180409e-03   +3.0142829e-03   +4.4269446e-03
  +2.4000000e+02   +2.4500000e+02   +3.0065407e-03   +1.2091507e-05  +3.0065407e-03   +2.2248135e-03   +4.1339886e-03   +3.0065407e-03   +2.5048262e-03   +3.6787282e-03   +2.6704416e-03   +2.2248135e-03   +3.2674859e-03   +3.3786146e-03   +2.8148106e-03   +4.1339886e-03
  +2.4500000e+02   +2.5000000e+02   +2.8413390e-03   +1.2526890e-05  +2.8413390e-03   +2.1018432e-03   +3.9088886e-03   +2.8413390e-03   +2.3671924e-03   +3.4765915e-03   +2.5228406e-03   +2.1018432e-03   +3.0868847e-03   +3.1946454e-03   +2.6615409e-03   +3.9088886e-03
  +2.5000000e+02   +2.5500000e+02   +2.6623355e-03   +1.0873792e-05  +2.6623355e-03   +1.9671571e-03   +3.6671904e-03   +2.6623355e-03   +2.2180603e-03   +3.2575676e-03   +2.3611767e-03   +1.9671571e-03   +2.8890771e-03   +2.9971108e-03   +2.4969700e-03   +3.6671904e-03
  +2.5500000e+02   +2.6000000e+02   +2.4941120e-03   +1.0071714e-05  +2.4941120e-03   +1.8418550e-03   +3.4377540e-03   +2.4941120e-03   +2.0779087e-03   +3.0517332e-03   +2.2107769e-03   +1.8418550e-03   +2.7050516e-03   +2.8095978e-03   +2.3407483e-03   +3.4377540e-03
  +2.6000000e+02   +2.6500000e+02   +2.3458055e-03   +9.7935820e-06  +2.3458055e-03   +1.7306472e-03   +3.2369819e-03   +2.3458055e-03   +1.9543508e-03   +2.8702689e-03   +2.0772942e-03   +1.7306472e-03   +2.5417254e-03   +2.6455118e-03   +2.2040438e-03   +3.2369819e-03
  +2.6500000e+02   +2.7000000e+02   +2.2011593e-03   +9.9826069e-06  +2.2011593e-03   +1.6235708e-03   +3.0386735e-03   +2.2011593e-03   +1.8338425e-03   +2.6932837e-03   +1.9487706e-03   +1.6235708e-03   +2.3844673e-03   +2.4834386e-03   +2.0690161e-03   +3.0386735e-03
  +2.7000000e+02   +2.7500000e+02   +2.0646091e-03   +9.2978740e-06  +2.0646091e-03   +1.5216132e-03   +2.8528715e-03   +2.0646091e-03   +1.7200788e-03   +2.5262039e-03   +1.8263909e-03   +1.5216132e-03   +2.2347262e-03   +2.3315867e-03   +1.9425049e-03   +2.8528715e-03
  +2.7500000e+02   +2.8000000e+02   +1.9467816e-03   +9.6105968e-06  +1.9467816e-03   +1.4335072e-03   +2.6927530e-03   +1.9467816e-03   +1.6219138e-03   +2.3820334e-03   +1.7206373e-03   +1.4335072e-03   +2.1053289e-03   +2.2007255e-03   +1.8334809e-03   +2.6927530e-03
  +2.8000000e+02   +2.8500000e+02   +1.8185331e-03   +8.6122525e-06  +1.8185331e-03   +1.3382923e-03   +2.5171284e-03   +1.8185331e-03   +1.5150666e-03   +2.2251117e-03   +1.6063511e-03   +1.3382923e-03   +1.9654911e-03   +2.0571917e-03   +1.7138993e-03   +2.5171284e-03
  +2.8500000e+02   +2.9000000e+02   +1.7177457e-03   +8.6768277e-06  +1.7177457e-03   +1.2631586e-03   +2.3798322e-03   +1.7177457e-03   +1.4310981e-03   +2.1017909e-03   +1.5161682e-03   +1.2631586e-03   +1.8551456e-03   +1.9449826e-03   +1.6204149e-03   +2.3798322e-03
  +2.9000000e+02   +2.9500000e+02   +1.6186123e-03   +8.8654435e-06  +1.6186123e-03   +1.1903846e-03   +2.2427379e-03   +1.6186123e-03   +1.3485075e-03   +1.9804937e-03   +1.4288176e-03   +1.1903846e-03   +1.7482655e-03   +1.8329386e-03   +1.5270682e-03   +2.2427379e-03
  +2.9500000e+02   +3.0000000e+02   +1.5176887e-03   +7.8043835e-06  +1.5176887e-03   +1.1145056e-03   +2.1061524e-03   +1.5176887e-03   +1.2644254e-03   +1.8570059e-03   +1.3377401e-03   +1.1145056e-03   +1.6368254e-03   +1.7213103e-03   +1.4340678e-03   +2.1061524e-03
  +3.0000000e+02   +3.0500000e+02   +1.4387624e-03   +7.6136610e-06  +1.4387624e-03   +1.0557262e-03   +1.9985154e-03   +1.4387624e-03   +1.1986700e-03   +1.7604338e-03   +1.2671871e-03   +1.0557262e-03   +1.5504985e-03   +1.6333411e-03   +1.3607784e-03   +1.9985154e-03
  +3.0500000e+02   +3.1000000e+02   +1.3495297e-03   +7.1970315e-06  +1.3495297e-03   +9.8879350e-04   +1.8774838e-03   +1.3495297e-03   +1.1243279e-03   +1.6512509e-03   +1.1868479e-03   +9.8879350e-04   +1.4521975e-03   +1.5344247e-03   +1.2783686e-03   +1.8774838e-03
  +3.1000000e+02   +3.1500000e+02   +1.2681227e-03   +7.5129792e-06  +1.2681227e-03   +9.2951314e-04   +1.7638604e-03   +1.2681227e-03   +1.0565056e-03   +1.5516434e-03   +1.1156937e-03   +9.2951314e-04   +1.3651350e-03   +1.4415629e-03   +1.2010030e-03   +1.7638604e-03
  +3.1500000e+02   +3.2000000e+02   +1.2103627e-03   +6.8990006e-06  +1.2103627e-03   +8.8633952e-04   +1.6853553e-03   +1.2103627e-03   +1.0083843e-03   +1.4809696e-03   +1.0638725e-03   +8.8633952e-04   +1.3017279e-03   +1.3774025e-03   +1.1475494e-03   +1.6853553e-03
  +3.2000000e+02   +3.2500000e+02   +1.1241013e-03   +6.8997961e-06  +1.1241013e-03   +8.2262832e-04   +1.5665514e-03   +1.1241013e-03   +9.3651770e-04   +1.3754224e-03   +9.8740003e-04   +8.2262832e-04   +1.2081581e-03   +1.2803068e-03   +1.0666564e-03   +1.5665514e-03
  +3.2500000e+02   +3.3000000e+02   +1.0854982e-03   +7.7589335e-06  +1.0854982e-03   +7.9433768e-04   +1.5131686e-03   +1.0854982e-03   +9.0435641e-04   +1.3281885e-03   +9.5344281e-04   +7.9433768e-04   +1.1666088e-03   +1.2366783e-03   +1.0303084e-03   +1.5131686e-03
  +3.3000000e+02   +3.3500000e+02   +1.0095907e-03   +6.5330780e-06  +1.0095907e-03   +7.3793923e-04   +1.4090590e-03   +1.0095907e-03   +8.4111594e-04   +1.2353100e-03   +8.8574780e-04   +7.3793923e-04   +1.0837789e-03   +1.1515918e-03   +9.5942073e-04   +1.4090590e-03
  +3.3500000e+02   +3.4000000e+02   +9.5095407e-04   +6.2574611e-06  +9.5095407e-04   +6.9446990e-04   +1.3286078e-03   +9.5095407e-04   +7.9226423e-04   +1.1635637e-03   +8.3357157e-04   +6.9446990e-04   +1.0199374e-03   +1.0858409e-03   +9.0464193e-04   +1.3286078e-03
  +3.4000000e+02   +3.4500000e+02   +9.0069291e-04   +6.2267760e-06  +9.0069291e-04   +6.5764938e-04   +1.2587406e-03   +9.0069291e-04   +7.5039037e-04   +1.1020654e-03   +7.8937600e-04   +6.5764938e-04   +9.6586078e-04   +1.0287400e-03   +8.5706972e-04   +1.2587406e-03
  +3.4500000e+02   +3.5000000e+02   +8.4932549e-04   +5.7087983e-06  +8.4932549e-04   +6.1985461e-04   +1.1876665e-03   +8.4932549e-04   +7.0759487e-04   +1.0392135e-03   +7.4401090e-04   +6.1985461e-04   +9.1035321e-04   +9.7065277e-04   +8.0867574e-04   +1.1876665e-03
  +3.5000000e+02   +3.5500000e+02   +8.0779676e-04   +5.9109621e-06  +8.0779676e-04   +5.8888515e-04   +1.1310872e-03   +8.0779676e-04   +6.7299621e-04   +9.8840001e-04   +7.0683830e-04   +5.8888515e-04   +8.6486971e-04   +9.2441178e-04   +7.7015120e-04   +1.1310872e-03
  +3.5500000e+02   +3.6000000e+02   +7.5532941e-04   +6.0568272e-06  +7.5532941e-04   +5.5082201e-04   +1.0574030e-03   +7.5532941e-04   +6.2928431e-04   +9.2420220e-04   +6.6115117e-04   +5.5082201e-04   +8.0896808e-04   +8.6419136e-04   +7.1998000e-04   +1.0574030e-03
  +3.6000000e+02   +3.6500000e+02   +7.1742223e-04   +7.9774383e-06  +7.1742223e-04   +5.2234906e-04   +1.0060620e-03   +7.1742223e-04   +5.9770289e-04   +8.7781998e-04   +6.2697512e-04   +5.2234906e-04   +7.6715111e-04   +8.2223157e-04   +6.8502221e-04   +1.0060620e-03
  +3.6500000e+02   +3.7000000e+02   +6.8831599e-04   +5.3684306e-06  +6.8831599e-04   +5.0080942e-04   +9.6613880e-04   +6.8831599e-04   +5.7345375e-04   +8.4220629e-04   +6.0112110e-04   +5.0080942e-04   +7.3551677e-04   +7.8960320e-04   +6.5783873e-04   +9.6613880e-04
  +3.7000000e+02   +3.7500000e+02   +6.4069567e-04   +4.5888100e-06  +6.4069567e-04   +4.6589041e-04   +8.9980700e-04   +6.4069567e-04   +5.3378005e-04   +7.8393928e-04   +5.5920781e-04   +4.6589041e-04   +6.8423275e-04   +7.3539173e-04   +6.1267373e-04   +8.9980700e-04
  +3.7500000e+02   +3.8000000e+02   +6.0568068e-04   +5.4620454e-06  +6.0568068e-04   +4.4040253e-04   +8.5088876e-04   +6.0568068e-04   +5.0460815e-04   +7.4109580e-04   +5.2861477e-04   +4.4040253e-04   +6.4679984e-04   +6.9541196e-04   +5.7936557e-04   +8.5088876e-04
  +3.8000000e+02   +3.8500000e+02   +5.7407798e-04   +4.6063883e-06  +5.7407798e-04   +4.1690549e-04   +8.0748444e-04   +5.7407798e-04   +4.7827912e-04   +7.0242749e-04   +5.0041127e-04   +4.1690549e-04   +6.1229071e-04   +6.5993859e-04   +5.4981181e-04   +8.0748444e-04
  +3.8500000e+02   +3.9000000e+02   +5.4680402e-04   +4.7711820e-06  +5.4680402e-04   +3.9722157e-04   +7.6911935e-04   +5.4680402e-04   +4.5555648e-04   +6.6905574e-04   +4.7678468e-04   +3.9722157e-04   +5.8338184e-04   +6.2858373e-04   +5.2368924e-04   +7.6911935e-04
  +3.9000000e+02   +3.9500000e+02   +5.0987343e-04   +4.9876716e-06  +5.0987343e-04   +3.6996376e-04   +7.1810909e-04   +5.0987343e-04   +4.2478863e-04   +6.2386841e-04   +4.4406715e-04   +3.6996376e-04   +5.4334950e-04   +5.8689417e-04   +4.8895663e-04   +7.1810909e-04
  +3.9500000e+02   +4.0000000e+02   +4.9430621e-04   +4.9932412e-06  +4.9430621e-04   +3.5893084e-04   +6.9574574e-04   +4.9430621e-04   +4.1181920e-04   +6.0482073e-04   +4.3082434e-04   +3.5893084e-04   +5.2714595e-04   +5.6861713e-04   +4.7372950e-04   +6.9574574e-04
  +4.0000000e+02   +4.0500000e+02   +4.6463676e-04   +4.1892086e-06  +4.6463676e-04   +3.3713296e-04   +6.5451705e-04   +4.6463676e-04   +3.8710082e-04   +5.6851796e-04   +4.0466037e-04   +3.3713296e-04   +4.9513234e-04   +5.3492187e-04   +4.4565713e-04   +6.5451705e-04
  +4.0500000e+02   +4.1000000e+02   +4.3572700e-04   +4.0181135e-06  +4.3572700e-04   +3.1570912e-04   +6.1473069e-04   +4.3572700e-04   +3.6301535e-04   +5.3314471e-04   +3.7894536e-04   +3.1570912e-04   +4.6366811e-04   +5.0240540e-04   +4.1856683e-04   +6.1473069e-04
  +4.1000000e+02   +4.1500000e+02   +4.1837948e-04   +4.1337313e-06  +4.1837948e-04   +3.0324578e-04   +5.9013294e-04   +4.1837948e-04   +3.4856270e-04   +5.1191869e-04   +3.6398561e-04   +3.0324578e-04   +4.4536372e-04   +4.8230221e-04   +4.0181834e-04   +5.9013294e-04
  +4.1500000e+02   +4.2000000e+02   +3.9336751e-04   +3.8428649e-06  +3.9336751e-04   +2.8477490e-04   +5.5564826e-04   +3.9336751e-04   +3.2772459e-04   +4.8131465e-04   +3.4181503e-04   +2.8477490e-04   +4.1823635e-04   +4.5411866e-04   +3.7833793e-04   +5.5564826e-04
  +4.2000000e+02   +4.2500000e+02   +3.7539995e-04   +4.1203487e-06  +3.7539995e-04   +2.7153488e-04   +5.3075919e-04   +3.7539995e-04   +3.1275534e-04   +4.5933004e-04   +3.2592309e-04   +2.7153488e-04   +3.9879136e-04   +4.3377740e-04   +3.6139109e-04   +5.3075919e-04
  +4.2500000e+02   +4.3000000e+02   +3.6193416e-04   +3.9629679e-06  +3.6193416e-04   +2.6215645e-04   +5.1112688e-04   +3.6193416e-04   +3.0153665e-04   +4.4285360e-04   +3.1466613e-04   +2.6215645e-04   +3.8501760e-04   +4.1773235e-04   +3.4802354e-04   +5.1112688e-04
  +4.3000000e+02   +4.3500000e+02   +3.4082085e-04   +4.2920781e-06  +3.4082085e-04   +2.4682820e-04   +4.8143463e-04   +3.4082085e-04   +2.8394662e-04   +4.1701994e-04   +2.9626766e-04   +2.4682820e-04   +3.6250573e-04   +3.9346551e-04   +3.2780624e-04   +4.8143463e-04
  +4.3500000e+02   +4.4000000e+02   +3.1772162e-04   +3.3282908e-06  +3.1772162e-04   +2.2945553e-04   +4.5009748e-04   +3.1772162e-04   +2.6470205e-04   +3.8875628e-04   +2.7541524e-04   +2.2945553e-04   +3.3699124e-04   +3.6785441e-04   +3.0646894e-04   +4.5009748e-04
  +4.4000000e+02   +4.4500000e+02   +3.0593130e-04   +3.8186138e-06  +3.0593130e-04   +2.2116040e-04   +4.3304546e-04   +3.0593130e-04   +2.5487920e-04   +3.7432988e-04   +2.6545865e-04   +2.2116040e-04   +3.2480856e-04   +3.5391819e-04   +2.9485835e-04   +4.3304546e-04
  +4.4500000e+02   +4.5000000e+02   +2.8605580e-04   +3.1213866e-06  +2.8605580e-04   +2.0623039e-04   +4.0600759e-04   +2.8605580e-04   +2.3832043e-04   +3.5001074e-04   +2.4753815e-04   +2.0623039e-04   +3.0288150e-04   +3.3182076e-04   +2.7644840e-04   +4.0600759e-04
  +4.5000000e+02   +4.5500000e+02   +2.7078216e-04   +3.3807251e-06  +2.7078216e-04   +1.9518756e-04   +3.8444236e-04   +2.7078216e-04   +2.2559558e-04   +3.3132233e-04   +2.3428342e-04   +1.9518756e-04   +2.8666339e-04   +3.1419598e-04   +2.6176472e-04   +3.8444236e-04
  +4.5500000e+02   +4.6000000e+02   +2.5650364e-04   +3.4920618e-06  +2.5650364e-04   +1.8492935e-04   +3.6418834e-04   +2.5650364e-04   +2.1369976e-04   +3.1385144e-04   +2.2197052e-04   +1.8492935e-04   +2.7159759e-04   +2.9764283e-04   +2.4797389e-04   +3.6418834e-04
  +4.6000000e+02   +4.6500000e+02   +2.4990018e-04   +3.4547070e-06  +2.4990018e-04   +1.8023252e-04   +3.5475337e-04   +2.4990018e-04   +2.0819826e-04   +3.0577163e-04   +2.1633293e-04   +1.8023252e-04   +2.6469954e-04   +2.8993183e-04   +2.4154965e-04   +3.5475337e-04
  +4.6500000e+02   +4.7000000e+02   +2.3674021e-04   +2.9390320e-06  +2.3674021e-04   +1.7068908e-04   +3.3616178e-04   +2.3674021e-04   +1.9723436e-04   +2.8966943e-04   +2.0487795e-04   +1.7068908e-04   +2.5068355e-04   +2.7473736e-04   +2.2889072e-04   +3.3616178e-04
  +4.7000000e+02   +4.7500000e+02   +2.2995045e-04   +3.7143205e-06  +2.2995045e-04   +1.6604169e-04   +3.2615744e-04   +2.2995045e-04   +1.9157763e-04   +2.8136164e-04   +1.9929968e-04   +1.6604169e-04   +2.4385813e-04   +2.6656103e-04   +2.2207887e-04   +3.2615744e-04
  +4.7500000e+02   +4.8000000e+02   +2.1546109e-04   +2.8128957e-06  +2.1546109e-04   +1.5522006e-04   +3.0629645e-04   +2.1546109e-04   +1.7950617e-04   +2.6363282e-04   +1.8631049e-04   +1.5522006e-04   +2.2796488e-04   +2.5032911e-04   +2.0855560e-04   +3.0629645e-04
  +4.8000000e+02   +4.8500000e+02   +2.1189764e-04   +6.0875692e-06  +2.1189764e-04   +1.5269211e-04   +3.0127052e-04   +2.1189764e-04   +1.7653738e-04   +2.5927269e-04   +1.8327621e-04   +1.5269211e-04   +2.2425219e-04   +2.4622154e-04   +2.0513348e-04   +3.0127052e-04
  +4.8500000e+02   +4.9000000e+02   +1.9573254e-04   +2.9760264e-06  +1.9573254e-04   +1.4090582e-04   +2.7852307e-04   +1.9573254e-04   +1.6306980e-04   +2.3949345e-04   +1.6912913e-04   +1.4090582e-04   +2.0694219e-04   +2.2763055e-04   +1.8964487e-04   +2.7852307e-04
  +4.9000000e+02   +4.9500000e+02   +1.8395329e-04   +2.4373592e-06  +1.8395329e-04   +1.3227928e-04   +2.6208388e-04   +1.8395329e-04   +1.5325622e-04   +2.2508067e-04   +1.5877470e-04   +1.3227928e-04   +1.9427277e-04   +2.1419519e-04   +1.7845151e-04   +2.6208388e-04
  +4.9500000e+02   +5.0000000e+02   +1.7631533e-04   +2.4606509e-06  +1.7631533e-04   +1.2664921e-04   +2.5151427e-04   +1.7631533e-04   +1.4689282e-04   +2.1573503e-04   +1.5201694e-04   +1.2664921e-04   +1.8600413e-04   +2.0555689e-04   +1.7125472e-04   +2.5151427e-04
<\histogram>


<histogram> 101 "pta lin NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   -5.4175987e-03   +6.6737901e-05  -5.4175987e-03   -5.9327464e-03   -4.7630935e-03   -5.4175987e-03   -4.9422373e-03   -5.8666306e-03   -5.2635973e-03   -4.7630935e-03   -5.7687605e-03   -5.5392839e-03   -5.0900837e-03   -5.9327464e-03
  +2.5000000e+01   +3.0000000e+01   -5.6225121e-03   +1.1667155e-04  -5.6225121e-03   -6.1268573e-03   -4.9539831e-03   -5.6225121e-03   -5.1331829e-03   -6.0813778e-03   -5.4730613e-03   -4.9539831e-03   -5.9959327e-03   -5.7336528e-03   -5.2765814e-03   -6.1268573e-03
  +3.0000000e+01   +3.5000000e+01   -5.7206354e-03   +1.1701418e-04  -5.7206354e-03   -6.2502224e-03   -5.0331031e-03   -5.7206354e-03   -5.2257686e-03   -6.1821581e-03   -5.5563252e-03   -5.0331031e-03   -6.0804627e-03   -5.8515473e-03   -5.3865446e-03   -6.2502224e-03
  +3.5000000e+01   +4.0000000e+01   -5.5926452e-03   +6.6261347e-05  -5.5926452e-03   -6.0750886e-03   -4.9360674e-03   -5.5926452e-03   -5.1300179e-03   -6.0061158e-03   -5.4278749e-03   -4.9360674e-03   -5.9054640e-03   -5.7258891e-03   -5.2938279e-03   -6.0750886e-03
  +4.0000000e+01   +4.5000000e+01   -5.4813107e-03   +6.7922388e-05  -5.4813107e-03   -5.9327967e-03   -4.8472201e-03   -5.4813107e-03   -5.0326560e-03   -5.8780611e-03   -5.3274635e-03   -4.8472201e-03   -5.7918224e-03   -5.6043118e-03   -5.1888871e-03   -5.9327967e-03
  +4.5000000e+01   +5.0000000e+01   -5.3662772e-03   +8.6670660e-05  -5.3662772e-03   -5.8235417e-03   -4.7361867e-03   -5.3662772e-03   -4.9245994e-03   -5.7590480e-03   -5.2073999e-03   -4.7361867e-03   -5.6644889e-03   -5.4959535e-03   -5.0855007e-03   -5.8235417e-03
  +5.0000000e+01   +5.5000000e+01   -4.9365812e-03   +8.4896249e-05  -4.9365812e-03   -5.2358227e-03   -4.4147372e-03   -4.9365812e-03   -4.5754824e-03   -5.2173198e-03   -4.8125802e-03   -4.4147372e-03   -5.1679213e-03   -5.0313389e-03   -4.7091172e-03   -5.2358227e-03
  +5.5000000e+01   +6.0000000e+01   -4.8277041e-03   +7.1411068e-05  -4.8277041e-03   -5.1459138e-03   -4.3044256e-03   -4.8277041e-03   -4.4665076e-03   -5.1166192e-03   -4.6999574e-03   -4.3044256e-03   -5.0594579e-03   -4.9289907e-03   -4.6040101e-03   -5.1459138e-03
  +6.0000000e+01   +6.5000000e+01   -4.6328456e-03   +7.4333528e-05  -4.6328456e-03   -4.9309000e-03   -4.1315542e-03   -4.6328456e-03   -4.2857688e-03   -4.9109173e-03   -4.5127799e-03   -4.1315542e-03   -4.8605400e-03   -4.7251548e-03   -4.4148520e-03   -4.9309000e-03
  +6.5000000e+01   +7.0000000e+01   -4.2282909e-03   +7.9385874e-05  -4.2282909e-03   -4.4236337e-03   -3.8124450e-03   -4.2282909e-03   -3.9443158e-03   -4.4236337e-03   -4.1346480e-03   -3.8124450e-03   -4.4050063e-03   -4.2958866e-03   -4.0527813e-03   -4.4134188e-03
  +7.0000000e+01   +7.5000000e+01   -4.2615461e-03   +7.6638164e-05  -4.2615461e-03   -4.5590693e-03   -3.7881784e-03   -4.2615461e-03   -3.9385714e-03   -4.5239523e-03   -4.1408012e-03   -3.7881784e-03   -4.4649237e-03   -4.3603814e-03   -4.0690817e-03   -4.5590693e-03
  +7.5000000e+01   +8.0000000e+01   -3.8030529e-03   +6.9127777e-05  -3.8030529e-03   -3.9606498e-03   -3.4376909e-03   -3.8030529e-03   -3.5577929e-03   -3.9606498e-03   -3.7192071e-03   -3.4376909e-03   -3.9475853e-03   -3.8643408e-03   -3.6578787e-03   -3.9482829e-03
  +8.0000000e+01   +8.5000000e+01   -3.5927956e-03   +6.6820027e-05  -3.5927956e-03   -3.7446835e-03   -3.2438095e-03   -3.5927956e-03   -3.3594094e-03   -3.7446835e-03   -3.5113470e-03   -3.2438095e-03   -3.7300888e-03   -3.6534927e-03   -3.4567961e-03   -3.7355275e-03
  +8.5000000e+01   +9.0000000e+01   -3.4489908e-03   +6.8304375e-05  -3.4489908e-03   -3.6027760e-03   -3.1060940e-03   -3.4489908e-03   -3.2204705e-03   -3.6027760e-03   -3.3661793e-03   -3.1060940e-03   -3.5823090e-03   -3.5129923e-03   -3.3183692e-03   -3.6016585e-03
  +9.0000000e+01   +9.5000000e+01   -3.1820190e-03   +6.7364731e-05  -3.1820190e-03   -3.2975644e-03   -2.8866370e-03   -3.1820190e-03   -2.9859641e-03   -3.2975644e-03   -3.1161113e-03   -2.8866370e-03   -3.2960637e-03   -3.2278801e-03   -3.0682068e-03   -3.2752126e-03
  +9.5000000e+01   +1.0000000e+02   -2.8234028e-03   +7.3697424e-05  -2.8234028e-03   -2.8938170e-03   -2.5967942e-03   -2.8234028e-03   -2.6786685e-03   -2.8738382e-03   -2.7769080e-03   -2.5967942e-03   -2.8938170e-03   -2.8504985e-03   -2.7442295e-03   -2.8303864e-03
  +1.0000000e+02   +1.0500000e+02   -2.8842576e-03   +7.6472109e-05  -2.8842576e-03   -3.0074659e-03   -2.6026965e-03   -2.8842576e-03   -2.6961822e-03   -3.0074659e-03   -2.8196855e-03   -2.6026965e-03   -2.9991731e-03   -2.9326126e-03   -2.7760191e-03   -2.9961575e-03
  +1.0500000e+02   +1.1000000e+02   -2.5159638e-03   +6.4002871e-05  -2.5159638e-03   -2.5705241e-03   -2.3084180e-03   -2.5159638e-03   -2.3870843e-03   -2.5607384e-03   -2.4678162e-03   -2.3084180e-03   -2.5705241e-03   -2.5503431e-03   -2.4541854e-03   -2.5342712e-03
  +1.1000000e+02   +1.1500000e+02   -2.4094730e-03   +6.4554311e-05  -2.4094730e-03   -2.4768076e-03   -2.2051693e-03   -2.4094730e-03   -2.2788242e-03   -2.4652282e-03   -2.3652852e-03   -2.2051693e-03   -2.4768076e-03   -2.4370153e-03   -2.3389443e-03   -2.4326834e-03
  +1.1500000e+02   +1.2000000e+02   -2.2531178e-03   +6.3749109e-05  -2.2531178e-03   -2.3143595e-03   -2.0628847e-03   -2.2531178e-03   -2.1321200e-03   -2.3031646e-03   -2.2116961e-03   -2.0628847e-03   -2.3143595e-03   -2.2810938e-03   -2.1904307e-03   -2.2750189e-03
  +1.2000000e+02   +1.2500000e+02   -2.0095229e-03   +5.8960525e-05  -2.0095229e-03   -2.0482796e-03   -1.8659651e-03   -2.0095229e-03   -1.9211091e-03   -2.0194000e-03   -1.9839014e-03   -1.8659651e-03   -2.0482796e-03   -2.0208340e-03   -1.9648477e-03   -1.9720853e-03
  +1.2500000e+02   +1.3000000e+02   -1.8976853e-03   +6.0063031e-05  -1.8976853e-03   -1.9352785e-03   -1.7528214e-03   -1.8976853e-03   -1.8091017e-03   -1.9160854e-03   -1.8677405e-03   -1.7528214e-03   -1.9352785e-03   -1.9149251e-03   -1.8553542e-03   -1.8803486e-03
  +1.3000000e+02   +1.3500000e+02   -1.7503743e-03   +6.0720266e-05  -1.7503743e-03   -1.7732376e-03   -1.6212796e-03   -1.7503743e-03   -1.6750401e-03   -1.7559873e-03   -1.7213701e-03   -1.6212796e-03   -1.7732376e-03   -1.7688824e-03   -1.7205623e-03   -1.7249886e-03
  +1.3500000e+02   +1.4000000e+02   -1.6196049e-03   +5.5111816e-05  -1.6196049e-03   -1.6406097e-03   -1.5072319e-03   -1.6196049e-03   -1.5539876e-03   -1.6175116e-03   -1.5973685e-03   -1.5072319e-03   -1.6406097e-03   -1.6311452e-03   -1.5924719e-03   -1.5801812e-03
  +1.4000000e+02   +1.4500000e+02   -1.5001832e-03   +5.2344061e-05  -1.5001832e-03   -1.5144950e-03   -1.3995744e-03   -1.5001832e-03   -1.4431374e-03   -1.4915904e-03   -1.4799759e-03   -1.3995744e-03   -1.5144950e-03   -1.5101614e-03   -1.4788941e-03   -1.4548905e-03
  +1.4500000e+02   +1.5000000e+02   -1.3396471e-03   +5.6084663e-05  -1.3396471e-03   -1.3424043e-03   -1.2614402e-03   -1.3396471e-03   -1.3017375e-03   -1.3087476e-03   -1.3295356e-03   -1.2674858e-03   -1.3424043e-03   -1.3397382e-03   -1.3284190e-03   -1.2614402e-03
  +1.5000000e+02   +1.5500000e+02   -1.3677426e-03   +5.4528987e-05  -1.3677426e-03   -1.3913578e-03   -1.2587844e-03   -1.3677426e-03   -1.3032293e-03   -1.3821936e-03   -1.3418826e-03   -1.2587844e-03   -1.3913578e-03   -1.3870267e-03   -1.3428061e-03   -1.3638926e-03
  +1.5500000e+02   +1.6000000e+02   -1.1325995e-03   +5.2331294e-05  -1.1325995e-03   -1.1330944e-03   -1.0477026e-03   -1.1325995e-03   -1.1096121e-03   -1.0903227e-03   -1.1246183e-03   -1.0800326e-03   -1.1214224e-03   -1.1323161e-03   -1.1330944e-03   -1.0477026e-03
  +1.6000000e+02   +1.6500000e+02   -1.0744582e-03   +5.0151160e-05  -1.0744582e-03   -1.0747845e-03   -9.9347485e-04   -1.0744582e-03   -1.0492499e-03   -1.0404130e-03   -1.0707335e-03   -1.0243039e-03   -1.0747845e-03   -1.0692369e-03   -1.0676502e-03   -9.9347485e-04
  +1.6500000e+02   +1.7000000e+02   -9.6646992e-04   +4.3523836e-05  -9.6646992e-04   -9.7285528e-04   -8.8977558e-04   -9.6646992e-04   -9.5032446e-04   -9.2420934e-04   -9.5787245e-04   -9.2314866e-04   -9.4935580e-04   -9.6858256e-04   -9.7285528e-04   -8.8977558e-04
  +1.7000000e+02   +1.7500000e+02   -9.4603129e-04   +4.4866214e-05  -9.4603129e-04   -9.4873409e-04   -8.8806032e-04   -9.4603129e-04   -9.1972061e-04   -9.2339043e-04   -9.3936516e-04   -8.9536941e-04   -9.4873409e-04   -9.4574061e-04   -9.3910164e-04   -8.8806032e-04
  +1.7500000e+02   +1.8000000e+02   -7.8908132e-04   +5.8360931e-05  -7.8908132e-04   -7.9668657e-04   -6.7597823e-04   -7.8908132e-04   -7.8867357e-04   -7.3180978e-04   -7.9668657e-04   -7.7620160e-04   -7.7463958e-04   -7.7235123e-04   -7.9457235e-04   -6.7597823e-04
  +1.8000000e+02   +1.8500000e+02   -7.9075509e-04   +5.6970562e-05  -7.9075509e-04   -7.9853135e-04   -7.3122652e-04   -7.9075509e-04   -7.7746861e-04   -7.5631421e-04   -7.8194825e-04   -7.5360606e-04   -7.7498834e-04   -7.9535206e-04   -7.9853135e-04   -7.3122652e-04
  +1.8500000e+02   +1.9000000e+02   -6.9079118e-04   +4.7214018e-05  -6.9079118e-04   -6.9747608e-04   -6.0048074e-04   -6.9079118e-04   -6.8886009e-04   -6.4345918e-04   -6.9466007e-04   -6.7584481e-04   -6.7713620e-04   -6.8059393e-04   -6.9747608e-04   -6.0048074e-04
  +1.9000000e+02   +1.9500000e+02   -7.3517100e-04   +4.6472541e-05  -7.3517100e-04   -7.4059449e-04   -6.8938151e-04   -7.3517100e-04   -7.1042612e-04   -7.2523778e-04   -7.2699358e-04   -6.8938151e-04   -7.4059449e-04   -7.3825939e-04   -7.2800117e-04   -7.0227969e-04
  +1.9500000e+02   +2.0000000e+02   -5.0247530e-04   +4.8379802e-05  -5.0247530e-04   -5.2532602e-04   -3.6504103e-04   -5.0247530e-04   -5.2532602e-04   -4.2481634e-04   -5.1977182e-04   -5.2435680e-04   -4.7339634e-04   -4.7748448e-04   -5.2088257e-04   -3.6504103e-04
  +2.0000000e+02   +2.0500000e+02   -5.9940997e-04   +5.2466503e-05  -5.9940997e-04   -5.9940997e-04   -5.4427881e-04   -5.9940997e-04   -5.8876337e-04   -5.7432816e-04   -5.9868501e-04   -5.7510003e-04   -5.9671635e-04   -5.9503677e-04   -5.9897676e-04   -5.4427881e-04
  +2.0500000e+02   +2.1000000e+02   -4.5960048e-04   +5.1234893e-05  -4.5960048e-04   -4.7767750e-04   -3.6301320e-04   -4.5960048e-04   -4.7350397e-04   -4.0103937e-04   -4.6615689e-04   -4.6617210e-04   -4.3186682e-04   -4.4830846e-04   -4.7767750e-04   -3.6301320e-04
  +2.1000000e+02   +2.1500000e+02   -5.5988361e-04   +4.7494845e-05  -5.5988361e-04   -5.7045188e-04   -5.2686246e-04   -5.5988361e-04   -5.4030378e-04   -5.5362853e-04   -5.5725035e-04   -5.2686246e-04   -5.7045188e-04   -5.5818239e-04   -5.5088856e-04   -5.3015519e-04
  +2.1500000e+02   +2.2000000e+02   -4.1147499e-04   +3.5091567e-05  -4.1147499e-04   -4.2137818e-04   -3.2446098e-04   -4.1147499e-04   -4.2012893e-04   -3.6580735e-04   -4.2137818e-04   -4.1695873e-04   -3.9828321e-04   -3.9557391e-04   -4.1916034e-04   -3.2446098e-04
  +2.2000000e+02   +2.2500000e+02   -3.5028074e-04   +3.9279560e-05  -3.5028074e-04   -3.6732719e-04   -2.5099139e-04   -3.5028074e-04   -3.6732719e-04   -2.9415316e-04   -3.6254659e-04   -3.6649763e-04   -3.2885718e-04   -3.3235158e-04   -3.6429476e-04   -2.5099139e-04
  +2.2500000e+02   +2.3000000e+02   -2.8656230e-04   +4.8734219e-05  -2.8656230e-04   -3.1098007e-04   -1.8641290e-04   -2.8656230e-04   -3.1098007e-04   -2.2198013e-04   -2.9708841e-04   -3.0961113e-04   -2.5293300e-04   -2.7319420e-04   -3.1061872e-04   -1.8641290e-04
  +2.3000000e+02   +2.3500000e+02   -4.3995169e-04   +4.5931153e-05  -4.3995169e-04   -4.5059075e-04   -4.0909719e-04   -4.3995169e-04   -4.2123456e-04   -4.4097418e-04   -4.3551586e-04   -4.0909719e-04   -4.5059075e-04   -4.4080125e-04   -4.3106223e-04   -4.2575962e-04
  +2.3500000e+02   +2.4000000e+02   -3.0106153e-04   +4.3856170e-05  -3.0106153e-04   -3.1057558e-04   -2.2972252e-04   -3.0106153e-04   -3.1057558e-04   -2.6197638e-04   -3.0941667e-04   -3.0867637e-04   -2.8799367e-04   -2.8863905e-04   -3.0979235e-04   -2.2972252e-04
  +2.4000000e+02   +2.4500000e+02   -2.7833792e-04   +6.6446114e-05  -2.7833792e-04   -2.9191472e-04   -2.2224062e-04   -2.7833792e-04   -2.8705749e-04   -2.4233902e-04   -2.8016395e-04   -2.8077910e-04   -2.5847530e-04   -2.7411411e-04   -2.9191472e-04   -2.2224062e-04
  +2.4500000e+02   +2.5000000e+02   -2.6903656e-04   +4.6947881e-05  -2.6903656e-04   -2.7578007e-04   -2.1144551e-04   -2.6903656e-04   -2.7504192e-04   -2.3855882e-04   -2.7578007e-04   -2.7301593e-04   -2.6043647e-04   -2.5892850e-04   -2.7489196e-04   -2.1144551e-04
  +2.5000000e+02   +2.5500000e+02   -2.3199674e-04   +3.2368010e-05  -2.3199674e-04   -2.4131154e-04   -1.7276924e-04   -2.3199674e-04   -2.4131154e-04   -1.9834310e-04   -2.3821554e-04   -2.3945357e-04   -2.1849997e-04   -2.2224568e-04   -2.4084041e-04   -1.7276924e-04
  +2.5500000e+02   +2.6000000e+02   -2.0443770e-04   +2.9221867e-05  -2.0443770e-04   -2.1682395e-04   -1.4058306e-04   -2.0443770e-04   -2.1682395e-04   -1.6733546e-04   -2.1263588e-04   -2.1679942e-04   -1.8958630e-04   -1.9339451e-04   -2.1505023e-04   -1.4058306e-04
  +2.6000000e+02   +2.6500000e+02   -1.9638263e-04   +2.7166460e-05  -1.9638263e-04   -2.0629446e-04   -1.3854583e-04   -1.9638263e-04   -2.0575772e-04   -1.6523924e-04   -2.0465084e-04   -2.0629446e-04   -1.8667937e-04   -1.8466908e-04   -2.0293206e-04   -1.3854583e-04
  +2.6500000e+02   +2.7000000e+02   -1.8931796e-04   +3.3676433e-05  -1.8931796e-04   -1.9695117e-04   -1.3751004e-04   -1.8931796e-04   -1.9674243e-04   -1.6217036e-04   -1.9666368e-04   -1.9695117e-04   -1.8169639e-04   -1.7860722e-04   -1.9430017e-04   -1.3751004e-04
  +2.7000000e+02   +2.7500000e+02   -1.7060959e-04   +2.7715388e-05  -1.7060959e-04   -1.7868799e-04   -1.2242120e-04   -1.7060959e-04   -1.7822383e-04   -1.4449887e-04   -1.7768949e-04   -1.7868799e-04   -1.6284953e-04   -1.6116613e-04   -1.7625856e-04   -1.2242120e-04
  +2.7500000e+02   +2.8000000e+02   -1.5686250e-04   +2.6847624e-05  -1.5686250e-04   -1.6513627e-04   -1.1067738e-04   -1.5686250e-04   -1.6513627e-04   -1.3058671e-04   -1.6250611e-04   -1.6473344e-04   -1.4659200e-04   -1.4872760e-04   -1.6394322e-04   -1.1067738e-04
  +2.8000000e+02   +2.8500000e+02   -1.0029555e-04   +3.1988876e-05  -1.0029555e-04   -1.2369299e-04   -3.0383066e-05   -1.0029555e-04   -1.1858396e-04   -6.0328492e-05   -1.1274402e-04   -1.2369299e-04   -8.4943010e-05   -8.4292656e-05   -1.1106384e-04   -3.0383066e-05
  +2.8500000e+02   +2.9000000e+02   -1.1134986e-04   +3.5479070e-05  -1.1134986e-04   -1.2654911e-04   -5.6836274e-05   -1.1134986e-04   -1.2490649e-04   -7.9003720e-05   -1.1893719e-04   -1.2654911e-04   -9.6628579e-05   -1.0077084e-04   -1.2126575e-04   -5.6836274e-05
  +2.9000000e+02   +2.9500000e+02   -1.2692581e-04   +3.2687565e-05  -1.2692581e-04   -1.3374691e-04   -8.9903065e-05   -1.2692581e-04   -1.3374691e-04   -1.0543945e-04   -1.3158353e-04   -1.3347126e-04   -1.1854751e-04   -1.2082358e-04   -1.3318963e-04   -8.9903065e-05
  +2.9500000e+02   +3.0000000e+02   -1.1163523e-04   +2.8888037e-05  -1.1163523e-04   -1.2020409e-04   -7.2567317e-05   -1.1163523e-04   -1.1922284e-04   -8.9906627e-05   -1.1758231e-04   -1.2020409e-04   -1.0426714e-04   -1.0381805e-04   -1.1707068e-04   -7.2567317e-05
  +3.0000000e+02   +3.0500000e+02   -1.1656914e-04   +2.3535558e-05  -1.1656914e-04   -1.2149167e-04   -8.8915591e-05   -1.1656914e-04   -1.2083670e-04   -1.0039521e-04   -1.1948342e-04   -1.1973970e-04   -1.1024453e-04   -1.1276220e-04   -1.2149167e-04   -8.8915591e-05
  +3.0500000e+02   +3.1000000e+02   -8.9525086e-05   +2.2821599e-05  -8.9525086e-05   -9.8803782e-05   -5.2819967e-05   -8.9525086e-05   -9.7475023e-05   -6.8775849e-05   -9.5257007e-05   -9.8803782e-05   -8.1934140e-05   -8.2017804e-05   -9.5017635e-05   -5.2819967e-05
  +3.1000000e+02   +3.1500000e+02   -6.6906301e-05   +2.4148046e-05  -6.6906301e-05   -8.3481857e-05   -1.8002080e-05   -6.6906301e-05   -7.8679746e-05   -4.1004916e-05   -7.6958311e-05   -8.3481857e-05   -5.9674974e-05   -5.3926672e-05   -7.1859196e-05   -1.8002080e-05
  +3.1500000e+02   +3.2000000e+02   -1.0100724e-04   +2.4379718e-05  -1.0100724e-04   -1.0516457e-04   -8.3250248e-05   -1.0100724e-04   -1.0310477e-04   -8.9844741e-05   -1.0133391e-04   -1.0067047e-04   -9.5068210e-05   -9.9928050e-05   -1.0516457e-04   -8.3250248e-05
  +3.2000000e+02   +3.2500000e+02   -7.1139631e-05   +2.2554848e-05  -7.1139631e-05   -7.9652001e-05   -3.9992802e-05   -7.1139631e-05   -7.8390436e-05   -5.2987777e-05   -7.6061531e-05   -7.9652001e-05   -6.4071685e-05   -6.5054147e-05   -7.6432661e-05   -3.9992802e-05
  +3.2500000e+02   +3.3000000e+02   -5.0609592e-05   +1.8874666e-05  -5.0609592e-05   -6.4801333e-05   -9.9305797e-06   -5.0609592e-05   -6.1790645e-05   -2.6961904e-05   -5.7688509e-05   -6.4801333e-05   -4.0771242e-05   -4.1394231e-05   -5.7340650e-05   -9.9305797e-06
  +3.3000000e+02   +3.3500000e+02   -5.6946532e-05   +2.0698889e-05  -5.6946532e-05   -6.5734880e-05   -2.7164044e-05   -5.6946532e-05   -6.4256750e-05   -3.9731998e-05   -6.1731056e-05   -6.5734880e-05   -5.0057806e-05   -5.0734632e-05   -6.1867259e-05   -2.7164044e-05
  +3.3500000e+02   +3.4000000e+02   -3.9560825e-05   +1.9897827e-05  -3.9560825e-05   -5.3062382e-05   -3.3482294e-06   -3.9560825e-05   -4.9949117e-05   -1.8138121e-05   -4.6310700e-05   -5.3062382e-05   -3.0873377e-05   -3.1498027e-05   -4.5993258e-05   -3.3482294e-06
  +3.4000000e+02   +3.4500000e+02   -8.3355120e-05   +1.8847300e-05  -8.3355120e-05   -8.4080898e-05   -7.5610202e-05   -8.3355120e-05   -8.2107528e-05   -7.9452022e-05   -8.2828066e-05   -7.9796944e-05   -8.2142459e-05   -8.3235243e-05   -8.4080898e-05   -7.5610202e-05
  +3.4500000e+02   +3.5000000e+02   -1.6472307e-06   +2.5251720e-05  -1.6472307e-06   -2.4875089e-05   +4.9364602e-05   -1.6472307e-06   -1.9018203e-05   +2.9433904e-05   -1.1737436e-05   -2.4875089e-05   +1.2540363e-05   +1.0487918e-05   -1.1757393e-05   +4.9364602e-05
  +3.5000000e+02   +3.5500000e+02   -7.1134645e-05   +2.2848229e-05  -7.1134645e-05   -7.2023590e-05   -6.2374533e-05   -7.1134645e-05   -7.0234749e-05   -6.7510060e-05   -7.2023590e-05   -6.9281890e-05   -7.1616286e-05   -6.9657915e-05   -7.0871784e-05   -6.2374533e-05
  +3.5500000e+02   +3.6000000e+02   +7.1870247e-06   +1.8553295e-05  +7.1870247e-06   -1.5715134e-05   +5.5672985e-05   +7.1870247e-06   -9.8896463e-06   +3.7094320e-05   -2.5240142e-06   -1.5715134e-05   +2.1171892e-05   +1.8682233e-05   -2.8428425e-06   +5.5672985e-05
  +3.6000000e+02   +3.6500000e+02   -5.3251004e-05   +1.8669925e-05  -5.3251004e-05   -5.5043286e-05   -4.3214201e-05   -5.3251004e-05   -5.4364467e-05   -4.7352492e-05   -5.4030810e-05   -5.3567819e-05   -5.0884665e-05   -5.2169907e-05   -5.5043286e-05   -4.3214201e-05
  +3.6500000e+02   +3.7000000e+02   -5.2266549e-05   +1.6275928e-05  -5.2266549e-05   -5.3288163e-05   -4.4538948e-05   -5.2266549e-05   -5.2529222e-05   -4.7956746e-05   -5.2748044e-05   -5.1644826e-05   -5.0837269e-05   -5.1504163e-05   -5.3288163e-05   -4.4538948e-05
  +3.7000000e+02   +3.7500000e+02   -4.8118892e-06   +1.7283145e-05  -4.8118892e-06   -2.1847391e-05   +3.3524744e-05   -4.8118892e-06   -1.6275819e-05   +1.5973823e-05   -1.3819016e-05   -2.1847391e-05   +1.5055826e-06   +6.3257731e-06   -9.1958075e-06   +3.3524744e-05
  +3.7500000e+02   +3.8000000e+02   -4.5481489e-05   +1.7517775e-05  -4.5481489e-05   -4.7220840e-05   -3.9376317e-05   -4.5481489e-05   -4.5898215e-05   -4.1395825e-05   -4.5569048e-05   -4.4831427e-05   -4.3534330e-05   -4.5606193e-05   -4.7220840e-05   -3.9376317e-05
  +3.8000000e+02   +3.8500000e+02   -3.6710564e-05   +1.7452315e-05  -3.6710564e-05   -3.8991525e-05   -2.5302475e-05   -3.6710564e-05   -3.8505781e-05   -3.0812710e-05   -3.8776405e-05   -3.8991525e-05   -3.5542907e-05   -3.4198224e-05   -3.7779088e-05   -2.5302475e-05
  +3.8500000e+02   +3.9000000e+02   +2.1494589e-05   +2.1054022e-05  +2.1494589e-05   -1.6333240e-07   +6.5293902e-05   +2.1494589e-05   +6.2050629e-06   +4.7165067e-05   +1.1812096e-05   -1.6333240e-07   +3.2287572e-05   +3.3469742e-05   +1.4233678e-05   +6.5293902e-05
  +3.9000000e+02   +3.9500000e+02   -3.3188854e-05   +1.7010151e-05  -3.3188854e-05   -3.5315215e-05   -2.2718039e-05   -3.3188854e-05   -3.4527882e-05   -2.8362908e-05   -3.5315215e-05   -3.5213003e-05   -3.2901686e-05   -3.0352776e-05   -3.3384712e-05   -2.2718039e-05
  +3.9500000e+02   +4.0000000e+02   -2.1070554e-05   +1.8107338e-05  -2.1070554e-05   -2.6694651e-05   -5.7665993e-06   -2.1070554e-05   -2.5069996e-05   -1.2393682e-05   -2.4359514e-05   -2.6694651e-05   -1.8407178e-05   -1.7425834e-05   -2.3248922e-05   -5.7665993e-06
  +4.0000000e+02   +4.0500000e+02   -1.2457224e-05   +1.3321732e-05  -1.2457224e-05   -1.9245567e-05   +4.2172831e-06   -1.2457224e-05   -1.8025438e-05   -1.6174400e-06   -1.5186061e-05   -1.9245567e-05   -6.8346619e-06   -9.3764818e-06   -1.6616943e-05   +4.2172831e-06
  +4.0500000e+02   +4.1000000e+02   +1.2761845e-05   +1.6262946e-05  +1.2761845e-05   -2.5516137e-06   +4.3980656e-05   +1.2761845e-05   +1.9344853e-06   +3.1121321e-05   +5.8722355e-06   -2.5516137e-06   +2.0454343e-05   +2.1191451e-05   +7.5310384e-06   +4.3980656e-05
  +4.1000000e+02   +4.1500000e+02   -2.4331083e-05   +1.7674813e-05  -2.4331083e-05   -2.7342997e-05   -1.2922281e-05   -2.4331083e-05   -2.6186979e-05   -1.9234979e-05   -2.6839964e-05   -2.7342997e-05   -2.3970649e-05   -2.0746372e-05   -2.4280704e-05   -1.2922281e-05
  +4.1500000e+02   +4.2000000e+02   -1.3186406e-05   +1.5435383e-05  -1.3186406e-05   -1.7078820e-05   -2.9990712e-06   -1.3186406e-05   -1.6941728e-05   -5.5241331e-06   -1.4255420e-05   -1.7078820e-05   -8.1755266e-06   -1.2288645e-05   -1.6994035e-05   -2.9990712e-06
  +4.2000000e+02   +4.2500000e+02   -4.0449705e-06   +1.7446623e-05  -4.0449705e-06   -1.1671096e-05   +1.3187565e-05   -4.0449705e-06   -9.3539731e-06   +5.7145222e-06   -7.8756915e-06   -1.1671096e-05   -5.3231510e-07   +6.3170760e-07   -6.4391362e-06   +1.3187565e-05
  +4.2500000e+02   +4.3000000e+02   -5.0649150e-05   +5.0370514e-05  -5.0649150e-05   -5.7621360e-05   -4.3345470e-05   -5.0649150e-05   -4.4935115e-05   -5.7110443e-05   -4.9389300e-05   -4.3345470e-05   -5.6531000e-05   -5.1937388e-05   -4.6606374e-05   -5.7621360e-05
  +4.3000000e+02   +4.3500000e+02   -9.6487775e-05   +1.4629831e-04  -9.6487775e-05   -1.2090448e-04   -7.7441759e-05   -9.6487775e-05   -7.9396516e-05   -1.1985750e-04   -9.4676484e-05   -7.7441759e-05   -1.1843497e-04   -9.8088470e-05   -8.1241766e-05   -1.2090448e-04
  +4.3500000e+02   +4.4000000e+02   -2.2595380e-05   +1.6981474e-05  -2.2595380e-05   -2.3296854e-05   -1.9243227e-05   -2.2595380e-05   -2.2414855e-05   -2.1256326e-05   -2.3296854e-05   -2.2445755e-05   -2.3101406e-05   -2.1827293e-05   -2.2376995e-05   -1.9243227e-05
  +4.4000000e+02   +4.4500000e+02   +3.8712017e-06   +1.1579568e-05  +3.8712017e-06   -4.5205802e-06   +2.1208375e-05   +3.8712017e-06   -2.1761516e-06   +1.4364791e-05   +1.4969454e-07   -4.5205802e-06   +8.4624721e-06   +8.2339574e-06   +6.1475543e-07   +2.1208375e-05
  +4.4500000e+02   +4.5000000e+02   -9.1579605e-06   +1.0307812e-05  -9.1579605e-06   -1.2424706e-05   -5.5589032e-07   -9.1579605e-06   -1.1730993e-05   -3.8989557e-06   -1.0687318e-05   -1.2424706e-05   -6.8052608e-06   -7.3659002e-06   -1.0883071e-05   -5.5589032e-07
  +4.5000000e+02   +4.5500000e+02   +8.7792235e-06   +1.0577406e-05  +8.7792235e-06   -4.6457908e-07   +2.7381569e-05   +8.7792235e-06   +2.1037669e-06   +2.0031539e-05   +4.8055943e-06   -4.6457908e-07   +1.3845350e-05   +1.3552510e-05   +5.2344532e-06   +2.7381569e-05
  +4.5500000e+02   +4.6000000e+02   -1.3305630e-05   +1.0805442e-05  -1.3305630e-05   -1.4650806e-05   -1.0105193e-05   -1.3305630e-05   -1.4201997e-05   -1.0730072e-05   -1.3409173e-05   -1.3909767e-05   -1.1531385e-05   -1.3370244e-05   -1.4650806e-05   -1.0105193e-05
  +4.6000000e+02   +4.6500000e+02   -4.4521744e-07   +1.1273575e-05  -4.4521744e-07   -6.3870983e-06   +1.2593582e-05   -4.4521744e-07   -4.1854002e-06   +6.2535962e-06   -3.7955973e-06   -6.3870983e-06   +1.1021829e-06   +3.7532331e-06   -1.3616594e-06   +1.2593582e-05
  +4.6500000e+02   +4.7000000e+02   -2.5487567e-05   +1.0036500e-05  -2.5487567e-05   -3.0927280e-05   -2.1162909e-05   -2.5487567e-05   -2.3215590e-05   -2.7663520e-05   -2.3277936e-05   -2.1162909e-05   -2.5336576e-05   -2.8436182e-05   -2.5865793e-05   -3.0927280e-05
  +4.7000000e+02   +4.7500000e+02   +4.8462357e-06   +9.9832306e-06  +4.8462357e-06   -3.1145781e-06   +2.1525023e-05   +4.8462357e-06   +9.2031592e-08   +1.2963402e-05   +2.1490849e-07   -3.1145781e-06   +6.1332153e-06   +1.0726481e-05   +4.2253193e-06   +2.1525023e-05
  +4.7500000e+02   +4.8000000e+02   -2.5430709e-06   +1.1200857e-05  -2.5430709e-06   -6.1817281e-06   +5.4861005e-06   -2.5430709e-06   -5.2210473e-06   +2.4167643e-06   -4.1982872e-06   -6.1817281e-06   -3.5464703e-07   -6.9804629e-07   -4.1390370e-06   +5.4861005e-06
  +4.8000000e+02   +4.8500000e+02   -1.0672609e-05   +9.9799524e-06  -1.0672609e-05   -1.1591888e-05   -7.1221937e-06   -1.0672609e-05   -1.1389518e-05   -8.6104232e-06   -1.1308376e-05   -1.1591888e-05   -9.9718979e-06   -9.9687881e-06   -1.1154836e-05   -7.1221937e-06
  +4.8500000e+02   +4.9000000e+02   -2.4477094e-06   +1.0665078e-05  -2.4477094e-06   -5.8185297e-06   +4.9675379e-06   -2.4477094e-06   -5.0298286e-06   +2.3342797e-06   -3.8807547e-06   -5.8185297e-06   -1.4183701e-07   -9.3696166e-07   -4.2111948e-06   +4.9675379e-06
  +4.9000000e+02   +4.9500000e+02   -4.3606929e-06   +1.1945315e-05  -4.3606929e-06   -7.1065557e-06   +2.3272789e-06   -4.3606929e-06   -6.0231483e-06   -1.0772010e-06   -6.0919629e-06   -7.1065557e-06   -3.8358947e-06   -2.1720257e-06   -4.6069072e-06   +2.3272789e-06
  +4.9500000e+02   +5.0000000e+02   +1.3424208e-05   +1.3947650e-05  +1.3424208e-05   +5.0741106e-06   +2.9128908e-05   +1.3424208e-05   +7.6253768e-06   +2.2772705e-05   +9.7698599e-06   +5.0741106e-06   +1.7420947e-05   +1.7798375e-05   +1.0707081e-05   +2.9128908e-05
<\histogram>


<histogram> 101 "pta lin LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   -3.8750365e-03   +7.5817423e-06  -3.8750365e-03   -5.0602578e-03   -2.9848902e-03   -3.8750365e-03   -3.2283925e-03   -4.7413979e-03   -3.5827608e-03   -2.9848902e-03   -4.3837766e-03   -4.1356335e-03   -3.4455026e-03   -5.0602578e-03
  +2.5000000e+01   +3.0000000e+01   -4.0405066e-03   +7.5057754e-06  -4.0405066e-03   -5.2833087e-03   -3.1089466e-03   -4.0405066e-03   -3.3662500e-03   -4.9438632e-03   -3.7316656e-03   -3.1089466e-03   -4.5659728e-03   -4.3179279e-03   -3.5973765e-03   -5.2833087e-03
  +3.0000000e+01   +3.5000000e+01   -4.1251657e-03   +7.8136351e-06  -4.1251657e-03   -5.4004000e-03   -3.1709748e-03   -4.1251657e-03   -3.4367816e-03   -5.0474498e-03   -3.8061182e-03   -3.1709748e-03   -4.6570715e-03   -4.4136241e-03   -3.6771036e-03   -5.4004000e-03
  +3.5000000e+01   +4.0000000e+01   -4.1326009e-03   +7.8207584e-06  -4.1326009e-03   -5.4167672e-03   -3.1734995e-03   -4.1326009e-03   -3.4429758e-03   -5.0565473e-03   -3.8091484e-03   -3.1734995e-03   -4.6607790e-03   -4.4270005e-03   -3.6882476e-03   -5.4167672e-03
  +4.0000000e+01   +4.5000000e+01   -4.0727730e-03   +7.8487913e-06  -4.0727730e-03   -5.3447735e-03   -3.1244721e-03   -4.0727730e-03   -3.3931319e-03   -4.9833437e-03   -3.7503010e-03   -3.1244721e-03   -4.5887743e-03   -4.3681615e-03   -3.6392275e-03   -5.3447735e-03
  +4.5000000e+01   +5.0000000e+01   -3.9758089e-03   +7.6750723e-06  -3.9758089e-03   -5.2237908e-03   -3.0470195e-03   -3.9758089e-03   -3.3123483e-03   -4.8647006e-03   -3.6573347e-03   -3.0470195e-03   -4.4750233e-03   -4.2692849e-03   -3.5568510e-03   -5.2237908e-03
  +5.0000000e+01   +5.5000000e+01   -3.8704693e-03   +7.5186955e-06  -3.8704693e-03   -5.0907906e-03   -2.9637571e-03   -3.8704693e-03   -3.2245874e-03   -4.7358095e-03   -3.5573949e-03   -2.9637571e-03   -4.3527396e-03   -4.1605872e-03   -3.4662921e-03   -5.0907906e-03
  +5.5000000e+01   +6.0000000e+01   -3.7471233e-03   +7.4139432e-06  -3.7471233e-03   -4.9335427e-03   -2.8669260e-03   -3.7471233e-03   -3.1218250e-03   -4.5848869e-03   -3.4411688e-03   -2.8669260e-03   -4.2105281e-03   -4.0320719e-03   -3.3592226e-03   -4.9335427e-03
  +6.0000000e+01   +6.5000000e+01   -3.5937168e-03   +7.2083818e-06  -3.5937168e-03   -4.7366577e-03   -2.7470743e-03   -3.5937168e-03   -2.9940178e-03   -4.3971822e-03   -3.2973107e-03   -2.7470743e-03   -4.0345072e-03   -3.8711625e-03   -3.2251649e-03   -4.7366577e-03
  +6.5000000e+01   +7.0000000e+01   -3.4344069e-03   +7.0485324e-06  -3.4344069e-03   -4.5312134e-03   -2.6231852e-03   -3.4344069e-03   -2.8612926e-03   -4.2022545e-03   -3.1486066e-03   -2.6231852e-03   -3.8525566e-03   -3.7032574e-03   -3.0852788e-03   -4.5312134e-03
  +7.0000000e+01   +7.5000000e+01   -3.2881972e-03   +6.8384219e-06  -3.2881972e-03   -4.3424505e-03   -2.5095664e-03   -3.2881972e-03   -2.7394817e-03   -4.0233563e-03   -3.0122303e-03   -2.5095664e-03   -3.6856899e-03   -3.5489858e-03   -2.9567512e-03   -4.3424505e-03
  +7.5000000e+01   +8.0000000e+01   -3.1368575e-03   +6.6248926e-06  -3.1368575e-03   -4.1464855e-03   -2.3922414e-03   -3.1368575e-03   -2.6133967e-03   -3.8381808e-03   -2.8714051e-03   -2.3922414e-03   -3.5133798e-03   -3.3888281e-03   -2.8233196e-03   -4.1464855e-03
  +8.0000000e+01   +8.5000000e+01   -2.9554908e-03   +6.3759054e-06  -2.9554908e-03   -3.9115981e-03   -2.2516552e-03   -2.9554908e-03   -2.4622955e-03   -3.6162648e-03   -2.7026594e-03   -2.2516552e-03   -3.3069069e-03   -3.1968598e-03   -2.6633858e-03   -3.9115981e-03
  +8.5000000e+01   +9.0000000e+01   -2.8161074e-03   +6.1954935e-06  -2.8161074e-03   -3.7295181e-03   -2.1444024e-03   -2.8161074e-03   -2.3461714e-03   -3.4457189e-03   -2.5739241e-03   -2.1444024e-03   -3.1493894e-03   -3.0480499e-03   -2.5394087e-03   -3.7295181e-03
  +9.0000000e+01   +9.5000000e+01   -2.6677433e-03   +5.9956957e-06  -2.6677433e-03   -3.5372524e-03   -2.0294128e-03   -2.6677433e-03   -2.2225655e-03   -3.2641843e-03   -2.4359024e-03   -2.0294128e-03   -2.9805093e-03   -2.8909158e-03   -2.4084962e-03   -3.5372524e-03
  +9.5000000e+01   +1.0000000e+02   -2.5047809e-03   +5.7905988e-06  -2.5047809e-03   -3.3239510e-03   -1.9042044e-03   -2.5047809e-03   -2.0867974e-03   -3.0647876e-03   -2.2856148e-03   -1.9042044e-03   -2.7966214e-03   -2.7165892e-03   -2.2632602e-03   -3.3239510e-03
  +1.0000000e+02   +1.0500000e+02   -2.3692660e-03   +5.5584272e-06  -2.3692660e-03   -3.1472411e-03   -1.7996842e-03   -2.3692660e-03   -1.9738962e-03   -2.8989746e-03   -2.1601597e-03   -1.7996842e-03   -2.6431174e-03   -2.5721685e-03   -2.1429396e-03   -3.1472411e-03
  +1.0500000e+02   +1.1000000e+02   -2.2324826e-03   +5.3548247e-06  -2.2324826e-03   -2.9677224e-03   -1.6948692e-03   -2.2324826e-03   -1.8599385e-03   -2.7316099e-03   -2.0343500e-03   -1.6948692e-03   -2.4891797e-03   -2.4254519e-03   -2.0207061e-03   -2.9677224e-03
  +1.1000000e+02   +1.1500000e+02   -2.1039539e-03   +5.1903805e-06  -2.1039539e-03   -2.8001859e-03   -1.5956409e-03   -2.1039539e-03   -1.7528579e-03   -2.5743452e-03   -1.9152464e-03   -1.5956409e-03   -2.3434476e-03   -2.2885282e-03   -1.9066316e-03   -2.8001859e-03
  +1.1500000e+02   +1.2000000e+02   -1.9729494e-03   +4.9640584e-06  -1.9729494e-03   -2.6275623e-03   -1.4955909e-03   -1.9729494e-03   -1.6437148e-03   -2.4140517e-03   -1.7951565e-03   -1.4955909e-03   -2.1965088e-03   -2.1474468e-03   -1.7890932e-03   -2.6275623e-03
  +1.2000000e+02   +1.2500000e+02   -1.8515299e-03   +4.8191042e-06  -1.8515299e-03   -2.4680435e-03   -1.4025817e-03   -1.8515299e-03   -1.5425571e-03   -2.2654856e-03   -1.6835175e-03   -1.4025817e-03   -2.0599102e-03   -2.0170758e-03   -1.6804777e-03   -2.4680435e-03
  +1.2500000e+02   +1.3000000e+02   -1.7245020e-03   +4.6154115e-06  -1.7245020e-03   -2.3011192e-03   -1.3051968e-03   -1.7245020e-03   -1.4367269e-03   -2.1100578e-03   -1.5666265e-03   -1.3051968e-03   -1.9168852e-03   -1.8806522e-03   -1.5668198e-03   -2.3011192e-03
  +1.3000000e+02   +1.3500000e+02   -1.6206610e-03   +4.4784858e-06  -1.6206610e-03   -2.1642911e-03   -1.2257672e-03   -1.6206610e-03   -1.3502143e-03   -1.9830004e-03   -1.4712872e-03   -1.2257672e-03   -1.8002304e-03   -1.7688255e-03   -1.4736540e-03   -2.1642911e-03
  +1.3500000e+02   +1.4000000e+02   -1.5188456e-03   +4.2979635e-06  -1.5188456e-03   -2.0296874e-03   -1.1482365e-03   -1.5188456e-03   -1.2653893e-03   -1.8584215e-03   -1.3782272e-03   -1.1482365e-03   -1.6863644e-03   -1.6588172e-03   -1.3820032e-03   -2.0296874e-03
  +1.4000000e+02   +1.4500000e+02   -1.4244428e-03   +4.1981497e-06  -1.4244428e-03   -1.9053009e-03   -1.0760616e-03   -1.4244428e-03   -1.1867399e-03   -1.7429126e-03   -1.2915958e-03   -1.0760616e-03   -1.5803643e-03   -1.5571590e-03   -1.2973091e-03   -1.9053009e-03
  +1.4500000e+02   +1.5000000e+02   -1.3334099e-03   +4.0438820e-06  -1.3334099e-03   -1.7849417e-03   -1.0066337e-03   -1.3334099e-03   -1.1108980e-03   -1.6315271e-03   -1.2082614e-03   -1.0066337e-03   -1.4783986e-03   -1.4587922e-03   -1.2153572e-03   -1.7849417e-03
  +1.5000000e+02   +1.5500000e+02   -1.2397781e-03   +3.9127725e-06  -1.2397781e-03   -1.6605772e-03   -9.3562430e-04   -1.2397781e-03   -1.0328909e-03   -1.5169615e-03   -1.1230290e-03   -9.3562430e-04   -1.3741103e-03   -1.3571519e-03   -1.1306781e-03   -1.6605772e-03
  +1.5500000e+02   +1.6000000e+02   -1.1700265e-03   +3.9039789e-06  -1.1700265e-03   -1.5682292e-03   -8.8249612e-04   -1.1700265e-03   -9.7477910e-04   -1.4316152e-03   -1.0592593e-03   -8.8249612e-04   -1.2960832e-03   -1.2816780e-03   -1.0677988e-03   -1.5682292e-03
  +1.6000000e+02   +1.6500000e+02   -1.0939406e-03   +3.6687592e-06  -1.0939406e-03   -1.4674822e-03   -8.2458917e-04   -1.0939406e-03   -9.1139005e-04   -1.3385185e-03   -9.8975364e-04   -8.2458917e-04   -1.2110379e-03   -1.1993398e-03   -9.9920081e-04   -1.4674822e-03
  +1.6500000e+02   +1.7000000e+02   -1.0147568e-03   +3.5338747e-06  -1.0147568e-03   -1.3622559e-03   -7.6441109e-04   -1.0147568e-03   -8.4541994e-04   -1.2416311e-03   -9.1752197e-04   -7.6441109e-04   -1.1226570e-03   -1.1133406e-03   -9.2755271e-04   -1.3622559e-03
  +1.7000000e+02   +1.7500000e+02   -9.4379514e-04   +3.4128226e-06  -9.4379514e-04   -1.2678069e-03   -7.1061368e-04   -9.4379514e-04   -7.8629996e-04   -1.1548042e-03   -8.5294894e-04   -7.1061368e-04   -1.0436471e-03   -1.0361497e-03   -8.6324290e-04   -1.2678069e-03
  +1.7500000e+02   +1.8000000e+02   -8.8869120e-04   +3.3966259e-06  -8.8869120e-04   -1.1947272e-03   -6.6873489e-04   -8.8869120e-04   -7.4039147e-04   -1.0873804e-03   -8.0268188e-04   -6.6873489e-04   -9.8214155e-04   -9.7642335e-04   -8.1348337e-04   -1.1947272e-03
  +1.8000000e+02   +1.8500000e+02   -8.2990268e-04   +3.2051439e-06  -8.2990268e-04   -1.1167189e-03   -6.2403473e-04   -8.2990268e-04   -6.9141324e-04   -1.0154482e-03   -7.4902831e-04   -6.2403473e-04   -9.1649232e-04   -9.1266891e-04   -7.6036791e-04   -1.1167189e-03
  +1.8500000e+02   +1.9000000e+02   -7.7057686e-04   +3.0817380e-06  -7.7057686e-04   -1.0372312e-03   -5.7935101e-04   -7.7057686e-04   -6.4198737e-04   -9.4285860e-04   -6.9539449e-04   -5.7935101e-04   -8.5086733e-04   -8.4770540e-04   -7.0624513e-04   -1.0372312e-03
  +1.9000000e+02   +1.9500000e+02   -7.1318082e-04   +2.9822988e-06  -7.1318082e-04   -9.6040851e-04   -5.3597310e-04   -7.1318082e-04   -5.9416926e-04   -8.7263030e-04   -6.4332801e-04   -5.3597310e-04   -7.8716010e-04   -7.8491998e-04   -6.5393703e-04   -9.6040851e-04
  +1.9500000e+02   +2.0000000e+02   -6.7478886e-04   +2.9067437e-06  -6.7478886e-04   -9.0959413e-04   -5.0675094e-04   -6.7478886e-04   -5.6218396e-04   -8.2565484e-04   -6.0825271e-04   -5.0675094e-04   -7.4424281e-04   -7.4339056e-04   -6.1933777e-04   -9.0959413e-04
  +2.0000000e+02   +2.0500000e+02   -6.2637480e-04   +2.8129837e-06  -6.2637480e-04   -8.4446634e-04   -4.7038111e-04   -6.2637480e-04   -5.2184895e-04   -7.6641660e-04   -5.6459804e-04   -4.7038111e-04   -6.9082806e-04   -6.9016310e-04   -5.7499261e-04   -8.4446634e-04
  +2.0500000e+02   +2.1000000e+02   -5.8424287e-04   +2.7177022e-06  -5.8424287e-04   -7.8842683e-04   -4.3839809e-04   -5.8424287e-04   -4.8674774e-04   -7.1486497e-04   -5.2620887e-04   -4.3839809e-04   -6.4385602e-04   -6.4436328e-04   -5.3683561e-04   -7.8842683e-04
  +2.1000000e+02   +2.1500000e+02   -5.3967423e-04   +2.5997306e-06  -5.3967423e-04   -7.2880902e-04   -4.0471619e-04   -5.3967423e-04   -4.4961645e-04   -6.6033193e-04   -4.8578050e-04   -4.0471619e-04   -5.9438889e-04   -5.9563897e-04   -4.9624212e-04   -7.2880902e-04
  +2.1500000e+02   +2.2000000e+02   -5.0519206e-04   +2.5969833e-06  -5.0519206e-04   -6.8226402e-04   -3.7888461e-04   -5.0519206e-04   -4.2088850e-04   -6.1814041e-04   -4.5477483e-04   -3.7888461e-04   -5.5645114e-04   -5.5759885e-04   -4.6454992e-04   -6.8226402e-04
  +2.2000000e+02   +2.2500000e+02   -4.7566478e-04   +2.4798652e-06  -4.7566478e-04   -6.4312842e-04   -3.5641997e-04   -4.7566478e-04   -3.9628858e-04   -5.8201159e-04   -4.2781053e-04   -3.5641997e-04   -5.2345831e-04   -5.2561423e-04   -4.3790271e-04   -6.4312842e-04
  +2.2500000e+02   +2.3000000e+02   -4.3847757e-04   +2.6147665e-06  -4.3847757e-04   -5.9303011e-04   -3.2851399e-04   -4.3847757e-04   -3.6530692e-04   -5.3651023e-04   -3.9431505e-04   -3.2851399e-04   -4.8247404e-04   -4.8467000e-04   -4.0379101e-04   -5.9303011e-04
  +2.3000000e+02   +2.3500000e+02   -4.0837518e-04   +2.2857406e-06  -4.0837518e-04   -5.5260829e-04   -3.0582643e-04   -4.0837518e-04   -3.4022785e-04   -4.9967766e-04   -3.6708318e-04   -3.0582643e-04   -4.4915382e-04   -4.5163416e-04   -3.7626802e-04   -5.5260829e-04
  +2.3500000e+02   +2.4000000e+02   -3.8462420e-04   +2.2709936e-06  -3.8462420e-04   -5.2055908e-04   -2.8804606e-04   -3.8462420e-04   -3.2044029e-04   -4.7061658e-04   -3.4574141e-04   -2.8804606e-04   -4.2304057e-04   -4.2544109e-04   -3.5444589e-04   -5.2055908e-04
  +2.4000000e+02   +2.4500000e+02   -3.5523349e-04   +2.1824345e-06  -3.5523349e-04   -4.8117233e-04   -2.6585939e-04   -3.5523349e-04   -2.9595415e-04   -4.3465484e-04   -3.1911076e-04   -2.6585939e-04   -3.9045600e-04   -3.9325117e-04   -3.2762764e-04   -4.8117233e-04
  +2.4500000e+02   +2.5000000e+02   -3.3194750e-04   +2.1689920e-06  -3.3194750e-04   -4.4960994e-04   -2.4847264e-04   -3.3194750e-04   -2.7655399e-04   -4.0616270e-04   -2.9824147e-04   -2.4847264e-04   -3.6492085e-04   -3.6745593e-04   -3.0613699e-04   -4.4960994e-04
  +2.5000000e+02   +2.5500000e+02   -3.0573396e-04   +2.0842907e-06  -3.0573396e-04   -4.1411924e-04   -2.2885928e-04   -3.0573396e-04   -2.5471481e-04   -3.7408844e-04   -2.7469962e-04   -2.2885928e-04   -3.3611560e-04   -3.3845024e-04   -2.8197157e-04   -4.1411924e-04
  +2.5500000e+02   +2.6000000e+02   -2.8909992e-04   +2.0665907e-06  -2.8909992e-04   -3.9192155e-04   -2.1628946e-04   -2.8909992e-04   -2.4085656e-04   -3.5373544e-04   -2.5961203e-04   -2.1628946e-04   -3.1765484e-04   -3.2030855e-04   -2.6685728e-04   -3.9192155e-04
  +2.6000000e+02   +2.6500000e+02   -2.6582142e-04   +1.9035693e-06  -2.6582142e-04   -3.6031306e-04   -1.9888972e-04   -2.6582142e-04   -2.2146264e-04   -3.2525245e-04   -2.3872714e-04   -1.9888972e-04   -2.9210058e-04   -2.9447566e-04   -2.4533524e-04   -3.6031306e-04
  +2.6500000e+02   +2.7000000e+02   -2.4865758e-04   +1.8957259e-06  -2.4865758e-04   -3.3712584e-04   -1.8605594e-04   -2.4865758e-04   -2.0716300e-04   -3.0425115e-04   -2.2332277e-04   -1.8605594e-04   -2.7325219e-04   -2.7552526e-04   -2.2954717e-04   -3.3712584e-04
  +2.7000000e+02   +2.7500000e+02   -2.2843604e-04   +1.8358581e-06  -2.2843604e-04   -3.0957758e-04   -1.7100499e-04   -2.2843604e-04   -1.9031592e-04   -2.7950863e-04   -2.0525714e-04   -1.7100499e-04   -2.5114753e-04   -2.5301070e-04   -2.1078971e-04   -3.0957758e-04
  +2.7500000e+02   +2.8000000e+02   -2.1602749e-04   +1.7774600e-06  -2.1602749e-04   -2.9310990e-04   -1.6155118e-04   -2.1602749e-04   -1.7997805e-04   -2.6432583e-04   -1.9390974e-04   -1.6155118e-04   -2.3726313e-04   -2.3955203e-04   -1.9957695e-04   -2.9310990e-04
  +2.8000000e+02   +2.8500000e+02   -1.9936565e-04   +1.6662189e-06  -1.9936565e-04   -2.7074031e-04   -1.4899318e-04   -1.9936565e-04   -1.6609664e-04   -2.4393882e-04   -1.7883638e-04   -1.4899318e-04   -2.1881974e-04   -2.2126989e-04   -1.8434563e-04   -2.7074031e-04
  +2.8500000e+02   +2.9000000e+02   -1.8954868e-04   +1.7471783e-06  -1.8954868e-04   -2.5731518e-04   -1.4171803e-04   -1.8954868e-04   -1.5791788e-04   -2.3192704e-04   -1.7010403e-04   -1.4171803e-04   -2.0813505e-04   -2.1029783e-04   -1.7520454e-04   -2.5731518e-04
  +2.9000000e+02   +2.9500000e+02   -1.7539469e-04   +1.6395082e-06  -1.7539469e-04   -2.3814690e-04   -1.3114096e-04   -1.7539469e-04   -1.4612581e-04   -2.1460857e-04   -1.5740837e-04   -1.3114096e-04   -1.9260096e-04   -1.9463203e-04   -1.6215294e-04   -2.3814690e-04
  +2.9500000e+02   +3.0000000e+02   -1.6174802e-04   +1.5428272e-06  -1.6174802e-04   -2.1977300e-04   -1.2085654e-04   -1.6174802e-04   -1.3475643e-04   -1.9791083e-04   -1.4506399e-04   -1.2085654e-04   -1.7749668e-04   -1.7961547e-04   -1.4964227e-04   -2.1977300e-04
  +3.0000000e+02   +3.0500000e+02   -1.5167445e-04   +1.5160762e-06  -1.5167445e-04   -2.0604989e-04   -1.1335902e-04   -1.5167445e-04   -1.2636388e-04   -1.8558508e-04   -1.3606473e-04   -1.1335902e-04   -1.6648541e-04   -1.6839989e-04   -1.4029827e-04   -2.0604989e-04
  +3.0500000e+02   +3.1000000e+02   -1.3850002e-04   +1.4375046e-06  -1.3850002e-04   -1.8820928e-04   -1.0349436e-04   -1.3850002e-04   -1.1538792e-04   -1.6946517e-04   -1.2422418e-04   -1.0349436e-04   -1.5199761e-04   -1.5381916e-04   -1.2815069e-04   -1.8820928e-04
  +3.1000000e+02   +3.1500000e+02   -1.3098525e-04   +1.4940968e-06  -1.3098525e-04   -1.7785115e-04   -9.7972182e-05   -1.3098525e-04   -1.0912718e-04   -1.6027028e-04   -1.1759592e-04   -9.7972182e-05   -1.4388743e-04   -1.4535369e-04   -1.2109789e-04   -1.7785115e-04
  +3.1500000e+02   +3.2000000e+02   -1.2388620e-04   +1.4283867e-06  -1.2388620e-04   -1.6818947e-04   -9.2682553e-05   -1.2388620e-04   -1.0321277e-04   -1.5158406e-04   -1.1124678e-04   -9.2682553e-05   -1.3611879e-04   -1.3745742e-04   -1.1451930e-04   -1.6818947e-04
  +3.2000000e+02   +3.2500000e+02   -1.1445500e-04   +1.4550183e-06  -1.1445500e-04   -1.5543568e-04   -8.5615910e-05   -1.1445500e-04   -9.5355404e-05   -1.4004429e-04   -1.0276470e-04   -8.5615910e-05   -1.2574033e-04   -1.2703404e-04   -1.0583532e-04   -1.5543568e-04
  +3.2500000e+02   +3.3000000e+02   -1.0980053e-04   +1.4724926e-06  -1.0980053e-04   -1.4909422e-04   -8.2165385e-05   -1.0980053e-04   -9.1477636e-05   -1.3434919e-04   -9.8623035e-05   -8.2165385e-05   -1.2067269e-04   -1.2185131e-04   -1.0151745e-04   -1.4909422e-04
  +3.3000000e+02   +3.3500000e+02   -9.8715586e-05   +1.2729922e-06  -9.8715586e-05   -1.3424273e-04   -7.3763129e-05   -9.8715586e-05   -8.2242494e-05   -1.2078593e-04   -8.8537818e-05   -7.3763129e-05   -1.0833267e-04   -1.0971352e-04   -9.1405151e-05   -1.3424273e-04
  +3.3500000e+02   +3.4000000e+02   -9.3595107e-05   +1.2874936e-06  -9.3595107e-05   -1.2706153e-04   -7.0048996e-05   -9.3595107e-05   -7.7976491e-05   -1.1452064e-04   -8.4079744e-05   -7.0048996e-05   -1.0287788e-04   -1.0384450e-04   -8.6515519e-05   -1.2706153e-04
  +3.4000000e+02   +3.4500000e+02   -8.8202408e-05   +1.2564471e-06  -8.8202408e-05   -1.1982389e-04   -6.5987328e-05   -8.8202408e-05   -7.3483693e-05   -1.0792227e-04   -7.9204535e-05   -6.5987328e-05   -9.6912689e-05   -9.7929331e-05   -8.1587437e-05   -1.1982389e-04
  +3.4500000e+02   +3.5000000e+02   -8.3701735e-05   +1.2120442e-06  -8.3701735e-05   -1.1375627e-04   -6.2605379e-05   -8.3701735e-05   -6.9734063e-05   -1.0241535e-04   -7.5145176e-05   -6.2605379e-05   -9.1945762e-05   -9.2970409e-05   -7.7456039e-05   -1.1375627e-04
  +3.5000000e+02   +3.5500000e+02   -7.6047703e-05   +1.2219614e-06  -7.6047703e-05   -1.0319615e-04   -5.6963677e-05   -7.6047703e-05   -6.3357296e-05   -9.3050075e-05   -6.8373452e-05   -5.6963677e-05   -8.3660049e-05   -8.4339855e-05   -7.0265702e-05   -1.0319615e-04
  +3.5500000e+02   +3.6000000e+02   -7.2344218e-05   +1.3855042e-06  -7.2344218e-05   -9.8287842e-05   -5.4151649e-05   -7.2344218e-05   -6.0271826e-05   -8.8518577e-05   -6.4998175e-05   -5.4151649e-05   -7.9530145e-05   -8.0328415e-05   -6.6923667e-05   -9.8287842e-05
  +3.6000000e+02   +3.6500000e+02   -6.5348833e-05   +1.3936536e-06  -6.5348833e-05   -8.8611767e-05   -4.8989120e-05   -6.5348833e-05   -5.4443793e-05   -7.9959202e-05   -5.8801596e-05   -4.8989120e-05   -7.1948170e-05   -7.2420378e-05   -6.0335277e-05   -8.8611767e-05
  +3.6500000e+02   +3.7000000e+02   -6.0229250e-05   +1.1131922e-06  -6.0229250e-05   -8.1412340e-05   -4.5269844e-05   -6.0229250e-05   -5.0178533e-05   -7.3695007e-05   -5.4337355e-05   -4.5269844e-05   -6.6485832e-05   -6.6536450e-05   -5.5433225e-05   -8.1412340e-05
  +3.7000000e+02   +3.7500000e+02   -5.9442841e-05   +1.0445441e-06  -5.9442841e-05   -8.0743140e-05   -4.4507619e-05   -5.9442841e-05   -4.9523355e-05   -7.2732778e-05   -5.3422455e-05   -4.4507619e-05   -6.5366381e-05   -6.5989529e-05   -5.4977572e-05   -8.0743140e-05
  +3.7500000e+02   +3.8000000e+02   -5.3296983e-05   +1.0299143e-06  -5.3296983e-05   -7.2290874e-05   -3.9954321e-05   -5.3296983e-05   -4.4403087e-05   -6.5212862e-05   -4.7957137e-05   -3.9954321e-05   -5.8679151e-05   -5.9081682e-05   -4.9222470e-05   -7.2290874e-05
  +3.8000000e+02   +3.8500000e+02   -4.9892182e-05   +9.6259486e-07  -4.9892182e-05   -6.7871230e-05   -3.7306767e-05   -4.9892182e-05   -4.1566460e-05   -6.1046831e-05   -4.4779279e-05   -3.7306767e-05   -5.4790809e-05   -5.5469610e-05   -4.6213158e-05   -6.7871230e-05
  +3.8500000e+02   +3.9000000e+02   -4.7775653e-05   +1.0612855e-06  -4.7775653e-05   -6.4670363e-05   -3.5899192e-05   -4.7775653e-05   -3.9803123e-05   -5.8457095e-05   -4.3089765e-05   -3.5899192e-05   -5.2723564e-05   -5.2853616e-05   -4.4033703e-05   -6.4670363e-05
  +3.9000000e+02   +3.9500000e+02   -4.3768030e-05   +9.3331656e-07  -4.3768030e-05   -5.9396115e-05   -3.2806524e-05   -4.3768030e-05   -3.6464270e-05   -5.3553470e-05   -3.9377642e-05   -3.2806524e-05   -4.8181498e-05   -4.8543091e-05   -4.0442494e-05   -5.9396115e-05
  +3.9500000e+02   +4.0000000e+02   -4.2624898e-05   +1.3584708e-06  -4.2624898e-05   -5.7634357e-05   -3.2071789e-05   -4.2624898e-05   -3.5511894e-05   -5.2154760e-05   -3.8495741e-05   -3.2071789e-05   -4.7102426e-05   -4.7103243e-05   -3.9242923e-05   -5.7634357e-05
  +4.0000000e+02   +4.0500000e+02   -4.0294499e-05   +9.1713020e-07  -4.0294499e-05   -5.4688993e-05   -3.0214324e-05   -4.0294499e-05   -3.3570383e-05   -4.9303345e-05   -3.6266228e-05   -3.0214324e-05   -4.4374450e-05   -4.4696069e-05   -3.7237443e-05   -5.4688993e-05
  +4.0500000e+02   +4.1000000e+02   -3.6608733e-05   +8.9089678e-07  -3.6608733e-05   -4.9484180e-05   -2.7543781e-05   -3.6608733e-05   -3.0499676e-05   -4.4793532e-05   -3.3060774e-05   -2.7543781e-05   -4.0452339e-05   -4.0442294e-05   -3.3693512e-05   -4.9484180e-05
  +4.1000000e+02   +4.1500000e+02   -3.6205728e-05   +9.7732100e-07  -3.6205728e-05   -4.9025800e-05   -2.7212787e-05   -3.6205728e-05   -3.0163924e-05   -4.4300427e-05   -3.2663484e-05   -2.7212787e-05   -3.9966223e-05   -4.0067669e-05   -3.3381404e-05   -4.9025800e-05
  +4.1500000e+02   +4.2000000e+02   -3.2576133e-05   +8.4323936e-07  -3.2576133e-05   -4.4096472e-05   -2.4487214e-05   -3.2576133e-05   -2.7140013e-05   -3.9859344e-05   -2.9391979e-05   -2.4487214e-05   -3.5963292e-05   -3.6039039e-05   -3.0025052e-05   -4.4096472e-05
  +4.2000000e+02   +4.2500000e+02   -2.9578590e-05   +8.5459351e-07  -2.9578590e-05   -3.9889764e-05   -2.2303560e-05   -2.9578590e-05   -2.4642684e-05   -3.6191624e-05   -2.6770942e-05   -2.2303560e-05   -3.2756259e-05   -3.2600997e-05   -2.7160728e-05   -3.9889764e-05
  +4.2500000e+02   +4.3000000e+02   -3.0244742e-05   +8.5538549e-07  -3.0244742e-05   -4.0993799e-05   -2.2718451e-05   -3.0244742e-05   -2.5197671e-05   -3.7006712e-05   -2.7268935e-05   -2.2718451e-05   -3.3365590e-05   -3.3503299e-05   -2.7912461e-05   -4.0993799e-05
  +4.3000000e+02   +4.3500000e+02   -2.8377945e-05   +9.7863691e-07  -2.8377945e-05   -3.8343375e-05   -2.1376094e-05   -2.8377945e-05   -2.3642396e-05   -3.4722548e-05   -2.5657705e-05   -2.1376094e-05   -3.1394129e-05   -3.1337163e-05   -2.6107798e-05   -3.8343375e-05
  +4.3500000e+02   +4.4000000e+02   -2.4652213e-05   +8.6905231e-07  -2.4652213e-05   -3.3312724e-05   -1.8565828e-05   -2.4652213e-05   -2.0538392e-05   -3.0163833e-05   -2.2284547e-05   -1.8565828e-05   -2.7266816e-05   -2.7225729e-05   -2.2682454e-05   -3.3312724e-05
  +4.4000000e+02   +4.4500000e+02   -2.4122533e-05   +7.9127804e-07  -2.4122533e-05   -3.2552145e-05   -1.8191173e-05   -2.4122533e-05   -2.0097101e-05   -2.9515729e-05   -2.1834847e-05   -1.8191173e-05   -2.6716575e-05   -2.6604126e-05   -2.2164580e-05   -3.2552145e-05
  +4.4500000e+02   +4.5000000e+02   -2.2459045e-05   +7.0254684e-07  -2.2459045e-05   -3.0448247e-05   -1.6869198e-05   -2.2459045e-05   -1.8711208e-05   -2.7480327e-05   -2.0248084e-05   -1.6869198e-05   -2.4775051e-05   -2.4884658e-05   -2.0732049e-05   -3.0448247e-05
  +4.5000000e+02   +4.5500000e+02   -2.1542289e-05   +7.9405555e-07  -2.1542289e-05   -2.9235938e-05   -1.6167861e-05   -2.1542289e-05   -1.7947434e-05   -2.6358609e-05   -1.9406269e-05   -1.6167861e-05   -2.3745027e-05   -2.3893867e-05   -1.9906594e-05   -2.9235938e-05
  +4.5500000e+02   +4.6000000e+02   -1.9240850e-05   +8.1494673e-07  -1.9240850e-05   -2.5846670e-05   -1.4573917e-05   -1.9240850e-05   -1.6030046e-05   -2.3542624e-05   -1.7493057e-05   -1.4573917e-05   -2.1404069e-05   -2.1123894e-05   -1.7598858e-05   -2.5846670e-05
  +4.6000000e+02   +4.6500000e+02   -1.8124217e-05   +6.9699966e-07  -1.8124217e-05   -2.4303182e-05   -1.3741339e-05   -1.8124217e-05   -1.5099751e-05   -2.2176338e-05   -1.6493717e-05   -1.3741339e-05   -2.0181301e-05   -1.9862437e-05   -1.6547907e-05   -2.4303182e-05
  +4.6500000e+02   +4.7000000e+02   -1.8062731e-05   +6.4058161e-07  -1.8062731e-05   -2.4452867e-05   -1.3586814e-05   -1.8062731e-05   -1.5048524e-05   -2.2101105e-05   -1.6308241e-05   -1.3586814e-05   -1.9954357e-05   -1.9984771e-05   -1.6649827e-05   -2.4452867e-05
  +4.7000000e+02   +4.7500000e+02   -1.6929435e-05   +7.2476442e-07  -1.6929435e-05   -2.2664633e-05   -1.2866934e-05   -1.6929435e-05   -1.4104348e-05   -2.0714434e-05   -1.5444169e-05   -1.2866934e-05   -1.8897099e-05   -1.8523287e-05   -1.5432228e-05   -2.2664633e-05
  +4.7500000e+02   +4.8000000e+02   -1.5487445e-05   +6.2793275e-07  -1.5487445e-05   -2.0800700e-05   -1.1733219e-05   -1.5487445e-05   -1.2902988e-05   -1.8950052e-05   -1.4083372e-05   -1.1733219e-05   -1.7232062e-05   -1.6999938e-05   -1.4163085e-05   -2.0800700e-05
  +4.8000000e+02   +4.8500000e+02   -1.5423548e-05   +1.2685618e-06  -1.5423548e-05   -2.0604370e-05   -1.1746388e-05   -1.5423548e-05   -1.2849754e-05   -1.8871868e-05   -1.4099179e-05   -1.1746388e-05   -1.7251405e-05   -1.6839482e-05   -1.4029405e-05   -2.0604370e-05
  +4.8500000e+02   +4.9000000e+02   -1.4928198e-05   +6.6867202e-07  -1.4928198e-05   -2.0169332e-05   -1.1255389e-05   -1.4928198e-05   -1.2437064e-05   -1.8265770e-05   -1.3509833e-05   -1.1255389e-05   -1.6530295e-05   -1.6483936e-05   -1.3733190e-05   -2.0169332e-05
  +4.9000000e+02   +4.9500000e+02   -1.2754370e-05   +7.2404730e-07  -1.2754370e-05   -1.7036249e-05   -9.7109399e-06   -1.2754370e-05   -1.0625994e-05   -1.5605930e-05   -1.1656031e-05   -9.7109399e-06   -1.4262030e-05   -1.3923338e-05   -1.1599890e-05   -1.7036249e-05
  +4.9500000e+02   +5.0000000e+02   -1.2169446e-05   +6.0291657e-07  -1.2169446e-05   -1.6300272e-05   -9.2450885e-06   -1.2169446e-05   -1.0138677e-05   -1.4890230e-05   -1.1096872e-05   -9.2450885e-06   -1.3577856e-05   -1.3321841e-05   -1.1098769e-05   -1.6300272e-05
<\histogram>


<histogram> 101 "pta quad NLO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   +5.1613695e-03   +3.2963470e-05  +5.1613695e-03   +4.4928587e-03   +5.7178399e-03   +5.1613695e-03   +4.6842093e-03   +5.6324413e-03   +4.9915831e-03   +4.4928587e-03   +5.5135688e-03   +5.3006917e-03   +4.8480412e-03   +5.7178399e-03
  +2.5000000e+01   +3.0000000e+01   +5.6985457e-03   +7.1578096e-05  +5.6985457e-03   +4.9592370e-03   +6.3065946e-03   +5.6985457e-03   +5.1703319e-03   +6.2211262e-03   +5.5127608e-03   +4.9592370e-03   +6.0941078e-03   +5.8478323e-03   +5.3492691e-03   +6.3065946e-03
  +3.0000000e+01   +3.5000000e+01   +6.1112897e-03   +7.4027076e-05  +6.1112897e-03   +5.3139125e-03   +6.7709514e-03   +6.1112897e-03   +5.5475165e-03   +6.6669108e-03   +5.9044135e-03   +5.3139125e-03   +6.5228773e-03   +6.2828846e-03   +5.7499381e-03   +6.7709514e-03
  +3.5000000e+01   +4.0000000e+01   +6.4607463e-03   +4.3313461e-05  +6.4607463e-03   +5.6073396e-03   +7.1731788e-03   +6.4607463e-03   +5.8641102e-03   +7.0492533e-03   +6.2306880e-03   +5.6073396e-03   +6.8837136e-03   +6.6545462e-03   +6.0891229e-03   +7.1731788e-03
  +4.0000000e+01   +4.5000000e+01   +6.6214590e-03   +4.6307863e-05  +6.6214590e-03   +5.7631499e-03   +7.3136782e-03   +6.6214590e-03   +6.0240684e-03   +7.1994983e-03   +6.3918871e-03   +5.7631499e-03   +7.0426689e-03   +6.8143604e-03   +6.2531833e-03   +7.3136782e-03
  +4.5000000e+01   +5.0000000e+01   +6.9139493e-03   +6.3668608e-05  +6.9139493e-03   +5.9960010e-03   +7.6732701e-03   +6.9139493e-03   +6.2750442e-03   +7.5444810e-03   +6.6656196e-03   +5.9960010e-03   +7.3691444e-03   +7.1231450e-03   +6.5207269e-03   +7.6732701e-03
  +5.0000000e+01   +5.5000000e+01   +6.9751672e-03   +6.4695833e-05  +6.9751672e-03   +6.0668150e-03   +7.6991133e-03   +6.9751672e-03   +6.3474615e-03   +7.5812390e-03   +6.7292722e-03   +6.0668150e-03   +7.4153544e-03   +7.1812795e-03   +6.5945819e-03   +7.6991133e-03
  +5.5000000e+01   +6.0000000e+01   +7.0546156e-03   +5.6003235e-05  +7.0546156e-03   +6.1406846e-03   +7.7722905e-03   +7.0546156e-03   +6.4269977e-03   +7.6546911e-03   +6.8051391e-03   +6.1406846e-03   +7.4892049e-03   +7.2648759e-03   +6.6805767e-03   +7.7722905e-03
  +6.0000000e+01   +6.5000000e+01   +7.2358157e-03   +6.5300376e-05  +7.2358157e-03   +6.2784021e-03   +8.0049471e-03   +7.2358157e-03   +6.5759131e-03   +7.8801133e-03   +6.9746563e-03   +6.2784021e-03   +7.7029385e-03   +7.4550001e-03   +6.8389962e-03   +8.0049471e-03
  +6.5000000e+01   +7.0000000e+01   +7.1251441e-03   +7.3317637e-05  +7.1251441e-03   +6.2081071e-03   +7.8250374e-03   +7.1251441e-03   +6.4990182e-03   +7.7173763e-03   +6.8751116e-03   +6.2081071e-03   +7.5585789e-03   +7.3346506e-03   +6.7570280e-03   +7.8250374e-03
  +7.0000000e+01   +7.5000000e+01   +7.4032352e-03   +7.2616932e-05  +7.4032352e-03   +6.3997213e-03   +8.2302135e-03   +7.4032352e-03   +6.7171358e-03   +8.0819173e-03   +7.1207203e-03   +6.3997213e-03   +7.8823690e-03   +7.6456742e-03   +7.0023981e-03   +8.2302135e-03
  +7.5000000e+01   +8.0000000e+01   +7.0774159e-03   +6.8105737e-05  +7.0774159e-03   +6.1794327e-03   +7.7325751e-03   +7.0774159e-03   +6.4712504e-03   +7.6375820e-03   +6.8303811e-03   +6.1794327e-03   +7.4885197e-03   +7.2828646e-03   +6.7301999e-03   +7.7325751e-03
  +8.0000000e+01   +8.5000000e+01   +7.2068028e-03   +7.1501762e-05  +7.2068028e-03   +6.2559149e-03   +7.9463113e-03   +7.2068028e-03   +6.5626973e-03   +7.8250776e-03   +6.9409150e-03   +6.2559149e-03   +7.6516198e-03   +7.4331274e-03   +6.8386527e-03   +7.9463113e-03
  +8.5000000e+01   +9.0000000e+01   +7.1866897e-03   +7.8891870e-05  +7.1866897e-03   +6.2361703e-03   +7.9234139e-03   +7.1866897e-03   +6.5478096e-03   +7.7971306e-03   +6.9158313e-03   +6.2361703e-03   +7.6188641e-03   +7.4186733e-03   +6.8295321e-03   +7.9234139e-03
  +9.0000000e+01   +9.5000000e+01   +7.0713484e-03   +8.0204207e-05  +7.0713484e-03   +6.1493834e-03   +7.7608593e-03   +7.0713484e-03   +6.4548950e-03   +7.6502960e-03   +6.8094808e-03   +6.1493834e-03   +7.4854618e-03   +7.2932519e-03   +6.7301174e-03   +7.7608593e-03
  +9.5000000e+01   +1.0000000e+02   +6.9303478e-03   +8.9699498e-05  +6.9303478e-03   +6.0365206e-03   +7.5773271e-03   +6.9303478e-03   +6.3346760e-03   +7.4826205e-03   +6.6780656e-03   +6.0365206e-03   +7.3306397e-03   +7.1410595e-03   +6.6017853e-03   +7.5773271e-03
  +1.0000000e+02   +1.0500000e+02   +6.9993991e-03   +9.1427096e-05  +6.9993991e-03   +6.0718587e-03   +7.7035793e-03   +6.9993991e-03   +6.3788776e-03   +7.5908846e-03   +6.7357552e-03   +6.0718587e-03   +7.4239258e-03   +7.2242253e-03   +6.6573466e-03   +7.7035793e-03
  +1.0500000e+02   +1.1000000e+02   +6.8513640e-03   +8.3135625e-05  +6.8513640e-03   +5.9466548e-03   +7.5284037e-03   +6.8513640e-03   +6.2552234e-03   +7.4102775e-03   +6.5856207e-03   +5.9466548e-03   +7.2403928e-03   +7.0806845e-03   +6.5374750e-03   +7.5284037e-03
  +1.1000000e+02   +1.1500000e+02   +6.6892417e-03   +9.6370837e-05  +6.6892417e-03   +5.8221436e-03   +7.3077424e-03   +6.6892417e-03   +6.1204245e-03   +7.2113723e-03   +6.4368971e-03   +5.8221436e-03   +7.0594478e-03   +6.9028702e-03   +6.3910346e-03   +7.3077424e-03
  +1.1500000e+02   +1.2000000e+02   +6.7909913e-03   +1.2717692e-04  +6.7909913e-03   +5.8680073e-03   +7.5102433e-03   +6.7909913e-03   +6.1808538e-03   +7.3792882e-03   +6.5185197e-03   +5.8680073e-03   +7.1988033e-03   +7.0297076e-03   +6.4701881e-03   +7.5102433e-03
  +1.2000000e+02   +1.2500000e+02   +6.4415412e-03   +1.2213048e-04  +6.4415412e-03   +5.6146990e-03   +7.0085609e-03   +6.4415412e-03   +5.8996792e-03   +6.9338361e-03   +6.2049262e-03   +5.6146990e-03   +6.8008143e-03   +6.6384525e-03   +6.1570261e-03   +7.0085609e-03
  +1.2500000e+02   +1.3000000e+02   +6.3619518e-03   +9.5607321e-05  +6.3619518e-03   +5.5269691e-03   +6.9579789e-03   +6.3619518e-03   +5.8162979e-03   +6.8668558e-03   +6.1178468e-03   +5.5269691e-03   +6.7212934e-03   +6.5696521e-03   +6.0808411e-03   +6.9579789e-03
  +1.3000000e+02   +1.3500000e+02   +6.2585340e-03   +1.0374642e-04  +6.2585340e-03   +5.4358412e-03   +6.8425820e-03   +6.2585340e-03   +5.7267029e-03   +6.7464028e-03   +6.0122839e-03   +5.4358412e-03   +6.5977627e-03   +6.4703527e-03   +5.9946688e-03   +6.8425820e-03
  +1.3500000e+02   +1.4000000e+02   +6.1498821e-03   +1.0268725e-04  +6.1498821e-03   +5.3368742e-03   +6.7300487e-03   +6.1498821e-03   +5.6253930e-03   +6.6326509e-03   +5.9050485e-03   +5.3368742e-03   +6.4836722e-03   +6.3613158e-03   +5.8920903e-03   +6.7300487e-03
  +1.4000000e+02   +1.4500000e+02   +5.9423593e-03   +9.7108249e-05  +5.9423593e-03   +5.1747854e-03   +6.4551379e-03   +5.9423593e-03   +5.4503031e-03   +6.3825754e-03   +5.7136046e-03   +5.1747854e-03   +6.2539814e-03   +6.1349338e-03   +5.7022727e-03   +6.4551379e-03
  +1.4500000e+02   +1.5000000e+02   +5.8426939e-03   +1.0917540e-04  +5.8426939e-03   +5.0858268e-03   +6.3492297e-03   +5.8426939e-03   +5.3577492e-03   +6.2775617e-03   +5.6174079e-03   +5.0858268e-03   +6.1519539e-03   +6.0337181e-03   +5.6078637e-03   +6.3492297e-03
  +1.5000000e+02   +1.5500000e+02   +5.8814428e-03   +1.1330262e-04  +5.8814428e-03   +5.0718657e-03   +6.4927564e-03   +5.8814428e-03   +5.3637855e-03   +6.3717668e-03   +5.6275684e-03   +5.0718657e-03   +6.2043596e-03   +6.1070153e-03   +5.6387339e-03   +6.4927564e-03
  +1.5500000e+02   +1.6000000e+02   +5.5865692e-03   +1.1160455e-04  +5.5865692e-03   +4.8582370e-03   +6.0727869e-03   +5.5865692e-03   +5.1261399e-03   +5.9965685e-03   +5.3633490e-03   +4.8582370e-03   +5.8693934e-03   +5.7784685e-03   +5.3750313e-03   +6.0727869e-03
  +1.6000000e+02   +1.6500000e+02   +5.4363569e-03   +1.1035244e-04  +5.4363569e-03   +4.7298002e-03   +5.9017492e-03   +5.4363569e-03   +4.9914272e-03   +5.8297721e-03   +5.2194749e-03   +4.7298002e-03   +5.7085805e-03   +5.6235460e-03   +5.2355399e-03   +5.9017492e-03
  +1.6500000e+02   +1.7000000e+02   +5.3164161e-03   +1.0148756e-04  +5.3164161e-03   +4.6185466e-03   +5.7780328e-03   +5.3164161e-03   +4.8780237e-03   +5.7069956e-03   +5.1001065e-03   +4.6185466e-03   +5.5835222e-03   +5.5009114e-03   +5.1185702e-03   +5.7780328e-03
  +1.7000000e+02   +1.7500000e+02   +5.1476204e-03   +9.9842323e-05  +5.1476204e-03   +4.4823197e-03   +5.5697760e-03   +5.1476204e-03   +4.7299926e-03   +5.5135990e-03   +4.9449951e-03   +4.4823197e-03   +5.4061542e-03   +5.3196578e-03   +4.9599419e-03   +5.5697760e-03
  +1.7500000e+02   +1.8000000e+02   +5.0343905e-03   +1.1842851e-04  +5.0343905e-03   +4.3812381e-03   +5.4464048e-03   +5.0343905e-03   +4.6269106e-03   +5.3906054e-03   +4.8327943e-03   +4.3812381e-03   +5.2823835e-03   +5.2048362e-03   +4.8546524e-03   +5.4464048e-03
  +1.8000000e+02   +1.8500000e+02   +5.0763898e-03   +1.2426767e-04  +5.0763898e-03   +4.3691555e-03   +5.6022253e-03   +5.0763898e-03   +4.6368532e-03   +5.4866509e-03   +4.8434049e-03   +4.3691555e-03   +5.3326594e-03   +5.2887500e-03   +4.8947599e-03   +5.6022253e-03
  +1.8500000e+02   +1.9000000e+02   +4.6667237e-03   +1.2323339e-04  +4.6667237e-03   +4.0811093e-03   +4.9962919e-03   +4.6667237e-03   +4.3044582e-03   +4.9693770e-03   +4.4897711e-03   +4.0811093e-03   +4.8881220e-03   +4.8115584e-03   +4.5094822e-03   +4.9962919e-03
  +1.9000000e+02   +1.9500000e+02   +4.6928446e-03   +1.1851325e-04  +4.6928446e-03   +4.0718529e-03   +5.0912766e-03   +4.6928446e-03   +4.3080292e-03   +5.0337673e-03   +4.4974821e-03   +4.0718529e-03   +4.9255065e-03   +4.8589902e-03   +4.5282769e-03   +5.0912766e-03
  +1.9500000e+02   +2.0000000e+02   +4.4165118e-03   +1.1819008e-04  +4.4165118e-03   +3.8671851e-03   +4.7115038e-03   +4.4165118e-03   +4.0806763e-03   +4.6904499e-03   +4.2491086e-03   +3.8671851e-03   +4.6174931e-03   +4.5540276e-03   +4.2778651e-03   +4.7115038e-03
  +2.0000000e+02   +2.0500000e+02   +4.3824078e-03   +1.2727574e-04  +4.3824078e-03   +3.8125736e-03   +4.7247031e-03   +4.3824078e-03   +4.0335711e-03   +4.6820240e-03   +4.2026729e-03   +3.8125736e-03   +4.5890424e-03   +4.5336151e-03   +4.2394354e-03   +4.7247031e-03
  +2.0500000e+02   +2.1000000e+02   +4.0264780e-03   +1.2211983e-04  +4.0264780e-03   +3.5510228e-03   +4.2336862e-03   +4.0264780e-03   +3.7441678e-03   +4.2336862e-03   +3.8812768e-03   +3.5510228e-03   +4.1846069e-03   +4.1410797e-03   +3.9207612e-03   +4.2293774e-03
  +2.1000000e+02   +2.1500000e+02   +4.2319800e-03   +1.3442517e-04  +4.2319800e-03   +3.6679747e-03   +4.5895503e-03   +4.2319800e-03   +3.8789752e-03   +4.5500814e-03   +4.0606039e-03   +3.6679747e-03   +4.4619238e-03   +4.3762664e-03   +4.0761224e-03   +4.5895503e-03
  +2.1500000e+02   +2.2000000e+02   +4.0125577e-03   +1.2826555e-04  +4.0125577e-03   +3.4903337e-03   +4.3181782e-03   +4.0125577e-03   +3.6953082e-03   +4.2830608e-03   +3.8471069e-03   +3.4903337e-03   +4.2002139e-03   +4.1517957e-03   +3.8872238e-03   +4.3181782e-03
  +2.2000000e+02   +2.2500000e+02   +3.8069462e-03   +1.1773337e-04  +3.8069462e-03   +3.3273247e-03   +4.0606670e-03   +3.8069462e-03   +3.5219664e-03   +4.0350490e-03   +3.6531561e-03   +3.3273247e-03   +3.9653707e-03   +3.9373652e-03   +3.7058078e-03   +4.0606670e-03
  +2.2500000e+02   +2.3000000e+02   +3.6726694e-03   +1.2871841e-04  +3.6726694e-03   +3.2038515e-03   +3.9298868e-03   +3.6726694e-03   +3.4009692e-03   +3.8869734e-03   +3.5132809e-03   +3.2038515e-03   +3.8065442e-03   +3.8136129e-03   +3.5911025e-03   +3.9298868e-03
  +2.3000000e+02   +2.3500000e+02   +3.8113122e-03   +1.3253757e-04  +3.8113122e-03   +3.2861210e-03   +4.1558282e-03   +3.8113122e-03   +3.4866463e-03   +4.1098253e-03   +3.6447465e-03   +3.2861210e-03   +4.0160286e-03   +3.9518003e-03   +3.6743547e-03   +4.1558282e-03
  +2.3500000e+02   +2.4000000e+02   +3.4962609e-03   +1.3953602e-04  +3.4962609e-03   +3.0548774e-03   +3.7226881e-03   +3.4962609e-03   +3.2326309e-03   +3.7091472e-03   +3.3578899e-03   +3.0548774e-03   +3.6511379e-03   +3.6111337e-03   +3.3996213e-03   +3.7226881e-03
  +2.4000000e+02   +2.4500000e+02   +3.3018104e-03   +1.5602508e-04  +3.3018104e-03   +2.9037761e-03   +3.4772840e-03   +3.3018104e-03   +3.0671906e-03   +3.4772840e-03   +3.1800682e-03   +2.9037761e-03   +3.4387635e-03   +3.3965286e-03   +3.2172279e-03   +3.4664455e-03
  +2.4500000e+02   +2.5000000e+02   +3.4791537e-03   +1.4882404e-04  +3.4791537e-03   +2.9974512e-03   +3.7948541e-03   +3.4791537e-03   +3.1864601e-03   +3.7450963e-03   +3.3217479e-03   +2.9974512e-03   +3.6555861e-03   +3.6166695e-03   +3.3675442e-03   +3.7948541e-03
  +2.5000000e+02   +2.5500000e+02   +3.2317830e-03   +1.1631739e-04  +3.2317830e-03   +2.8139003e-03   +3.4514224e-03   +3.2317830e-03   +2.9838467e-03   +3.4361388e-03   +3.0975906e-03   +2.8139003e-03   +3.3755327e-03   +3.3420035e-03   +3.1427877e-03   +3.4514224e-03
  +2.5500000e+02   +2.6000000e+02   +3.0743777e-03   +1.0827598e-04  +3.0743777e-03   +2.6861080e-03   +3.2631329e-03   +3.0743777e-03   +2.8480882e-03   +3.2517209e-03   +2.9484529e-03   +2.6861080e-03   +3.1993129e-03   +3.1795109e-03   +3.0014707e-03   +3.2631329e-03
  +2.6000000e+02   +2.6500000e+02   +2.8767701e-03   +1.1818749e-04  +2.8767701e-03   +2.5352477e-03   +3.0114137e-03   +2.8767701e-03   +2.6825881e-03   +3.0114137e-03   +2.7680839e-03   +2.5352477e-03   +2.9796085e-03   +2.9618800e-03   +2.8192200e-03   +2.9984405e-03
  +2.6500000e+02   +2.7000000e+02   +3.0084220e-03   +1.6118885e-04  +3.0084220e-03   +2.5999120e-03   +3.2543486e-03   +3.0084220e-03   +2.7651631e-03   +3.2208580e-03   +2.8735745e-03   +2.5999120e-03   +3.1501095e-03   +3.1251703e-03   +2.9238103e-03   +3.2543486e-03
  +2.7000000e+02   +2.7500000e+02   +2.8382640e-03   +1.5502193e-04  +2.8382640e-03   +2.4627245e-03   +3.0429057e-03   +2.8382640e-03   +2.6202422e-03   +3.0182266e-03   +2.7115883e-03   +2.4627245e-03   +2.9558309e-03   +2.9465114e-03   +2.7709136e-03   +3.0429057e-03
  +2.7500000e+02   +2.8000000e+02   +2.6088321e-03   +1.1814478e-04  +2.6088321e-03   +2.2914317e-03   +2.7317524e-03   +2.6088321e-03   +2.4348490e-03   +2.7271689e-03   +2.4990751e-03   +2.2914317e-03   +2.6854716e-03   +2.7005920e-03   +2.5717368e-03   +2.7317524e-03
  +2.8000000e+02   +2.8500000e+02   +2.6429955e-03   +1.3185568e-04  +2.6429955e-03   +2.2961906e-03   +2.8235991e-03   +2.6429955e-03   +2.4429130e-03   +2.8053373e-03   +2.5262968e-03   +2.2961906e-03   +2.7507273e-03   +2.7421141e-03   +2.5833110e-03   +2.8235991e-03
  +2.8500000e+02   +2.9000000e+02   +2.5006696e-03   +1.3964461e-04  +2.5006696e-03   +2.1915408e-03   +2.6329561e-03   +2.5006696e-03   +2.3233204e-03   +2.6329561e-03   +2.4020513e-03   +2.1915408e-03   +2.6006882e-03   +2.5767558e-03   +2.4450089e-03   +2.6221732e-03
  +2.9000000e+02   +2.9500000e+02   +2.4394060e-03   +1.3438591e-04  +2.4394060e-03   +2.1336642e-03   +2.5795872e-03   +2.4394060e-03   +2.2601535e-03   +2.5795872e-03   +2.3451209e-03   +2.1336642e-03   +2.5496293e-03   +2.5116931e-03   +2.3768403e-03   +2.5674277e-03
  +2.9500000e+02   +3.0000000e+02   +2.4804008e-03   +1.2512390e-04  +2.4804008e-03   +2.1403459e-03   +2.6784989e-03   +2.4804008e-03   +2.2811029e-03   +2.6532956e-03   +2.3660836e-03   +2.1403459e-03   +2.5945047e-03   +2.5790513e-03   +2.4168943e-03   +2.6784989e-03
  +3.0000000e+02   +3.0500000e+02   +2.1001345e-03   +1.0782318e-04  +2.1001345e-03   +1.8762080e-03   +2.1570656e-03   +2.1001345e-03   +1.9816998e-03   +2.1568579e-03   +2.0308013e-03   +1.8762080e-03   +2.1570656e-03   +2.1464041e-03   +2.0766067e-03   +2.1130408e-03
  +3.0500000e+02   +3.1000000e+02   +2.2338615e-03   +1.2491776e-04  +2.2338615e-03   +1.9435156e-03   +2.3682769e-03   +2.2338615e-03   +2.0713008e-03   +2.3594011e-03   +2.1334293e-03   +1.9435156e-03   +2.3151015e-03   +2.3172060e-03   +2.1929905e-03   +2.3682769e-03
  +3.1000000e+02   +3.1500000e+02   +2.0044224e-03   +1.3539852e-04  +2.0044224e-03   +1.7763199e-03   +2.0760222e-03   +2.0044224e-03   +1.8815879e-03   +2.0760222e-03   +1.9310820e-03   +1.7763199e-03   +2.0649766e-03   +2.0608552e-03   +1.9807335e-03   +2.0521810e-03
  +3.1500000e+02   +3.2000000e+02   +1.9969197e-03   +1.2226140e-04  +1.9969197e-03   +1.7584590e-03   +2.0759792e-03   +1.9969197e-03   +1.8702092e-03   +2.0759792e-03   +1.9144458e-03   +1.7584590e-03   +2.0517460e-03   +2.0620997e-03   +1.9756155e-03   +2.0646738e-03
  +3.2000000e+02   +3.2500000e+02   +2.1266043e-03   +1.1515234e-04  +2.1266043e-03   +1.8294267e-03   +2.3004115e-03   +2.1266043e-03   +1.9517589e-03   +2.2819220e-03   +2.0274224e-03   +1.8294267e-03   +2.2312906e-03   +2.2107832e-03   +2.0693228e-03   +2.3004115e-03
  +3.2500000e+02   +3.3000000e+02   +1.6974501e-03   +1.0494598e-04  +1.6974501e-03   +1.5380131e-03   +1.7311354e-03   +1.6974501e-03   +1.6196346e-03   +1.7113755e-03   +1.6510407e-03   +1.5380131e-03   +1.7311354e-03   +1.7202578e-03   +1.6899628e-03   +1.6478118e-03
  +3.3000000e+02   +3.3500000e+02   +1.8503301e-03   +1.2292982e-04  +1.8503301e-03   +1.6233252e-03   +1.9416864e-03   +1.8503301e-03   +1.7227655e-03   +1.9416864e-03   +1.7784134e-03   +1.6233252e-03   +1.9241105e-03   +1.9020323e-03   +1.8133618e-03   +1.9202672e-03
  +3.3500000e+02   +3.4000000e+02   +1.8850337e-03   +1.3957888e-04  +1.8850337e-03   +1.6251293e-03   +2.0285984e-03   +1.8850337e-03   +1.7385760e-03   +2.0075116e-03   +1.7930121e-03   +1.6251293e-03   +1.9604400e-03   +1.9663354e-03   +1.8503101e-03   +2.0285984e-03
  +3.4000000e+02   +3.4500000e+02   +1.8653339e-03   +1.9467151e-04  +1.8653339e-03   +1.6040017e-03   +2.0110631e-03   +1.8653339e-03   +1.7150002e-03   +1.9961681e-03   +1.7751850e-03   +1.6040017e-03   +1.9498057e-03   +1.9406142e-03   +1.8210567e-03   +2.0110631e-03
  +3.4500000e+02   +3.5000000e+02   +1.1693659e-03   +2.1020032e-04  +1.1693659e-03   +9.4201361e-04   +1.1981174e-03   +1.1693659e-03   +1.1759394e-03   +1.0717006e-03   +1.1749127e-03   +1.1389328e-03   +1.1526822e-03   +1.1360820e-03   +1.1981174e-03   +9.4201361e-04
  +3.5000000e+02   +3.5500000e+02   +1.8406519e-03   +1.4345513e-04  +1.8406519e-03   +1.5549594e-03   +2.0508204e-03   +1.8406519e-03   +1.6778116e-03   +1.9955910e-03   +1.7321076e-03   +1.5549594e-03   +1.9205255e-03   +1.9427769e-03   +1.8020370e-03   +2.0508204e-03
  +3.5500000e+02   +3.6000000e+02   +1.5481612e-03   +1.1190704e-04  +1.5481612e-03   +1.3513006e-03   +1.6229878e-03   +1.5481612e-03   +1.4423321e-03   +1.6229878e-03   +1.4791355e-03   +1.3513006e-03   +1.5982563e-03   +1.6058476e-03   +1.5307756e-03   +1.6216164e-03
  +3.6000000e+02   +3.6500000e+02   +1.5274595e-03   +1.5043522e-04  +1.5274595e-03   +1.3385097e-03   +1.6038929e-03   +1.5274595e-03   +1.4215825e-03   +1.6038929e-03   +1.4682163e-03   +1.3385097e-03   +1.5914755e-03   +1.5702889e-03   +1.4982277e-03   +1.5833054e-03
  +3.6500000e+02   +3.7000000e+02   +1.4013081e-03   +1.4595408e-04  +1.4013081e-03   +1.2367433e-03   +1.4481279e-03   +1.4013081e-03   +1.3202090e-03   +1.4428527e-03   +1.3408417e-03   +1.2367433e-03   +1.4278218e-03   +1.4481279e-03   +1.3993848e-03   +1.4285633e-03
  +3.7000000e+02   +3.7500000e+02   +1.4491641e-03   +1.1459202e-04  +1.4491641e-03   +1.2622115e-03   +1.5265071e-03   +1.4491641e-03   +1.3460059e-03   +1.5265071e-03   +1.3863422e-03   +1.2622115e-03   +1.5056815e-03   +1.4978965e-03   +1.4249554e-03   +1.5178016e-03
  +3.7500000e+02   +3.8000000e+02   +1.3063313e-03   +1.1056980e-04  +1.3063313e-03   +1.1594258e-03   +1.3528412e-03   +1.3063313e-03   +1.2276749e-03   +1.3505032e-03   +1.2622888e-03   +1.1594258e-03   +1.3528412e-03   +1.3341316e-03   +1.2893341e-03   +1.3159138e-03
  +3.8000000e+02   +3.8500000e+02   +1.1854966e-03   +9.4536084e-05  +1.1854966e-03   +1.0731300e-03   +1.2155940e-03   +1.1854966e-03   +1.1305316e-03   +1.1963250e-03   +1.1548419e-03   +1.0731300e-03   +1.2155940e-03   +1.1966059e-03   +1.1787525e-03   +1.1404743e-03
  +3.8500000e+02   +3.9000000e+02   +1.2944300e-03   +1.1922088e-04  +1.2944300e-03   +1.1234682e-03   +1.3650756e-03   +1.2944300e-03   +1.2043340e-03   +1.3598647e-03   +1.2314308e-03   +1.1234682e-03   +1.3333401e-03   +1.3490489e-03   +1.2844225e-03   +1.3650756e-03
  +3.9000000e+02   +3.9500000e+02   +1.0522098e-03   +1.3660657e-04  +1.0522098e-03   +9.5835628e-04   +1.0694615e-03   +1.0522098e-03   +1.0104191e-03   +1.0493545e-03   +1.0253772e-03   +9.5835628e-04   +1.0694615e-03   +1.0609896e-03   +1.0541953e-03   +9.9511511e-04
  +3.9500000e+02   +4.0000000e+02   +1.1796409e-03   +1.1732389e-04  +1.1796409e-03   +1.0381134e-03   +1.2282110e-03   +1.1796409e-03   +1.1037412e-03   +1.2282110e-03   +1.1343667e-03   +1.0381134e-03   +1.2225408e-03   +1.2116443e-03   +1.1646786e-03   +1.2062935e-03
  +4.0000000e+02   +4.0500000e+02   +1.1862114e-03   +9.8789654e-05  +1.1862114e-03   +1.0312920e-03   +1.2482877e-03   +1.1862114e-03   +1.1024628e-03   +1.2482877e-03   +1.1325568e-03   +1.0312920e-03   +1.2297960e-03   +1.2285353e-03   +1.1703322e-03   +1.2419718e-03
  +4.0500000e+02   +4.1000000e+02   +1.0393595e-03   +1.1537181e-04  +1.0393595e-03   +9.2786166e-04   +1.0673664e-03   +1.0393595e-03   +9.8379121e-04   +1.0620041e-03   +1.0045393e-03   +9.2786166e-04   +1.0673664e-03   +1.0582615e-03   +1.0328215e-03   +1.0258199e-03
  +4.1000000e+02   +4.1500000e+02   +1.0651901e-03   +1.3037607e-04  +1.0651901e-03   +9.3190316e-04   +1.1079230e-03   +1.0651901e-03   +9.9728562e-04   +1.1079230e-03   +1.0168072e-03   +9.3190316e-04   +1.0933942e-03   +1.1032325e-03   +1.0598929e-03   +1.0993876e-03
  +4.1500000e+02   +4.2000000e+02   +1.0776473e-03   +1.1040293e-04  +1.0776473e-03   +9.3360671e-04   +1.1390625e-03   +1.0776473e-03   +9.9874699e-04   +1.1390625e-03   +1.0282832e-03   +9.3360671e-04   +1.1214442e-03   +1.1171394e-03   +1.0613928e-03   +1.1343852e-03
  +4.2000000e+02   +4.2500000e+02   +1.0072321e-03   +1.2792059e-04  +1.0072321e-03   +8.8318927e-04   +1.0488874e-03   +1.0072321e-03   +9.4232243e-04   +1.0488874e-03   +9.6530656e-04   +8.8318927e-04   +1.0407124e-03   +1.0361216e-03   +9.9653016e-04   +1.0305291e-03
  +4.2500000e+02   +4.3000000e+02   +8.1736860e-04   +1.3008127e-04  +8.1736860e-04   +7.4882787e-04   +8.2768888e-04   +8.1736860e-04   +7.9182550e-04   +8.0281649e-04   +8.0102179e-04   +7.5302566e-04   +8.2768888e-04   +8.1861970e-04   +8.2401870e-04   +7.4882787e-04
  +4.3000000e+02   +4.3500000e+02   +1.0595643e-03   +1.0701035e-04  +1.0595643e-03   +8.9481028e-04   +1.1691289e-03   +1.0595643e-03   +9.6866264e-04   +1.1436981e-03   +9.9573056e-04   +8.9481028e-04   +1.1024127e-03   +1.1178845e-03   +1.0430314e-03   +1.1691289e-03
  +4.3500000e+02   +4.4000000e+02   +8.2527078e-04   +1.1184669e-04  +8.2527078e-04   +7.3813981e-04   +8.4459873e-04   +8.2527078e-04   +7.8585929e-04   +8.3485533e-04   +7.9509763e-04   +7.3813981e-04   +8.3817251e-04   +8.4459873e-04   +8.3001739e-04   +8.0850812e-04
  +4.4000000e+02   +4.4500000e+02   +7.1332687e-04   +1.1309591e-04  +7.1332687e-04   +6.1353589e-04   +7.3530236e-04   +7.1332687e-04   +6.9553343e-04   +6.9260977e-04   +7.1322803e-04   +6.7142938e-04   +7.3530236e-04   +6.9256757e-04   +7.0834911e-04   +6.1353589e-04
  +4.4500000e+02   +4.5000000e+02   +8.3047773e-04   +1.5160126e-04  +8.3047773e-04   +7.2781058e-04   +8.5939033e-04   +8.3047773e-04   +7.8046515e-04   +8.5857368e-04   +7.9216086e-04   +7.2781058e-04   +8.4862704e-04   +8.5939033e-04   +8.3049346e-04   +8.4772708e-04
  +4.5000000e+02   +4.5500000e+02   +5.9268061e-04   +1.4022267e-04  +5.9268061e-04   +4.6915632e-04   +6.0373303e-04   +5.9268061e-04   +5.9366515e-04   +5.4736325e-04   +5.9606081e-04   +5.7362515e-04   +5.9223654e-04   +5.7043335e-04   +6.0373303e-04   +4.6915632e-04
  +4.5500000e+02   +4.6000000e+02   +8.5980999e-04   +9.8751808e-05  +8.5980999e-04   +7.3270091e-04   +9.2958480e-04   +8.5980999e-04   +7.9111105e-04   +9.1905452e-04   +8.1169031e-04   +7.3270091e-04   +8.9281216e-04   +9.0194613e-04   +8.4924283e-04   +9.2958480e-04
  +4.6000000e+02   +4.6500000e+02   +7.5116621e-04   +9.8486604e-05  +7.5116621e-04   +6.6349723e-04   +7.7825673e-04   +7.5116621e-04   +7.0577428e-04   +7.7685636e-04   +7.2386635e-04   +6.6349723e-04   +7.7825673e-04   +7.6697785e-04   +7.4358827e-04   +7.5229130e-04
  +4.6500000e+02   +4.7000000e+02   +6.8274156e-04   +9.5307621e-05  +6.8274156e-04   +6.1428264e-04   +7.0084247e-04   +6.8274156e-04   +6.5034080e-04   +6.9030689e-04   +6.6290594e-04   +6.1428264e-04   +7.0084247e-04   +6.8939758e-04   +6.8054090e-04   +6.5451137e-04
  +4.7000000e+02   +4.7500000e+02   +8.0659087e-04   +1.2867179e-04  +8.0659087e-04   +6.8706096e-04   +8.7145183e-04   +8.0659087e-04   +7.3935988e-04   +8.6713072e-04   +7.6444762e-04   +6.8706096e-04   +8.4618540e-04   +8.4182488e-04   +7.9048535e-04   +8.7145183e-04
  +4.7500000e+02   +4.8000000e+02   +4.6668179e-04   +1.2892325e-04  +4.6668179e-04   +3.2731873e-04   +4.7854799e-04   +4.6668179e-04   +4.7782590e-04   +4.1251818e-04   +4.7854799e-04   +4.6702603e-04   +4.6390976e-04   +4.3512080e-04   +4.7766184e-04   +3.2731873e-04
  +4.8000000e+02   +4.8500000e+02   +7.7641035e-04   +1.0228717e-04  +7.7641035e-04   +6.5358842e-04   +8.5729451e-04   +7.7641035e-04   +7.0907990e-04   +8.3934585e-04   +7.2874624e-04   +6.5358842e-04   +8.0913563e-04   +8.2021609e-04   +7.6558687e-04   +8.5729451e-04
  +4.8500000e+02   +4.9000000e+02   +4.6255342e-04   +9.8229259e-05  +4.6255342e-04   +3.5681938e-04   +4.7600490e-04   +4.6255342e-04   +4.6685718e-04   +4.2088477e-04   +4.6426350e-04   +4.4945182e-04   +4.5653837e-04   +4.4492139e-04   +4.7600490e-04   +3.5681938e-04
  +4.9000000e+02   +4.9500000e+02   +6.2657611e-04   +9.8145621e-05  +6.2657611e-04   +5.4868830e-04   +6.5022364e-04   +6.2657611e-04   +5.8746834e-04   +6.5022364e-04   +5.9938393e-04   +5.4868830e-04   +6.4568286e-04   +6.4515147e-04   +6.2345110e-04   +6.3640852e-04
  +4.9500000e+02   +5.0000000e+02   +4.3439659e-04   +1.3364820e-04  +4.3439659e-04   +3.4473488e-04   +4.5464496e-04   +4.3439659e-04   +4.3911414e-04   +3.9406006e-04   +4.3145956e-04   +4.1902994e-04   +4.2190007e-04   +4.2640604e-04   +4.5464496e-04   +3.4473488e-04
<\histogram>


<histogram> 101 "pta quad LO |X_AXIS@LIN |Y_AXIS@LOG |TYPE@#1"
  -5.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +0.0000000e+00   +5.0000000e+00   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +5.0000000e+00   +1.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.0000000e+01   +1.5000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +1.5000000e+01   +2.0000000e+01   +0.0000000e+00   +0.0000000e+00  +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00   +0.0000000e+00
  +2.0000000e+01   +2.5000000e+01   +3.5773621e-03   +4.9579255e-06  +3.5773621e-03   +2.7360928e-03   +4.7108701e-03   +3.5773621e-03   +2.9803921e-03   +4.3771708e-03   +3.2841295e-03   +2.7360928e-03   +4.0183791e-03   +3.8500867e-03   +3.2076061e-03   +4.7108701e-03
  +2.5000000e+01   +3.0000000e+01   +3.9431202e-03   +5.1260991e-06  +3.9431202e-03   +3.0110670e-03   +5.2017996e-03   +3.9431202e-03   +3.2851149e-03   +4.8247037e-03   +3.6141812e-03   +3.0110670e-03   +4.4222221e-03   +4.2513119e-03   +3.5418772e-03   +5.2017996e-03
  +3.0000000e+01   +3.5000000e+01   +4.2414347e-03   +5.5228776e-06  +4.2414347e-03   +3.2342042e-03   +5.6043986e-03   +4.2414347e-03   +3.5336483e-03   +5.1897135e-03   +3.8820125e-03   +3.2342042e-03   +4.7499341e-03   +4.5803472e-03   +3.8160050e-03   +5.6043986e-03
  +3.5000000e+01   +4.0000000e+01   +4.4810238e-03   +5.8224339e-06  +4.4810238e-03   +3.4120745e-03   +5.9304445e-03   +4.4810238e-03   +3.7332558e-03   +5.4828690e-03   +4.0955099e-03   +3.4120745e-03   +5.0111640e-03   +4.8468172e-03   +4.0380080e-03   +5.9304445e-03
  +4.0000000e+01   +4.5000000e+01   +4.6588602e-03   +6.1400989e-06  +4.6588602e-03   +3.5425015e-03   +6.1756293e-03   +4.6588602e-03   +3.8814159e-03   +5.7004651e-03   +4.2520612e-03   +3.5425015e-03   +5.2027161e-03   +5.0472009e-03   +4.2049528e-03   +6.1756293e-03
  +4.5000000e+01   +5.0000000e+01   +4.7933916e-03   +6.3287190e-06  +4.7933916e-03   +3.6394290e-03   +6.3643180e-03   +4.7933916e-03   +3.9934977e-03   +5.8650748e-03   +4.3684032e-03   +3.6394290e-03   +5.3450690e-03   +5.2014118e-03   +4.3334297e-03   +6.3643180e-03
  +5.0000000e+01   +5.5000000e+01   +4.9152517e-03   +6.5848381e-06  +4.9152517e-03   +3.7273939e-03   +6.5352278e-03   +4.9152517e-03   +4.0950226e-03   +6.0141795e-03   +4.4739872e-03   +3.7273939e-03   +5.4742594e-03   +5.3410927e-03   +4.4498016e-03   +6.5352278e-03
  +5.5000000e+01   +6.0000000e+01   +5.0053363e-03   +6.7814693e-06  +5.0053363e-03   +3.7908936e-03   +6.6645017e-03   +5.0053363e-03   +4.1700743e-03   +6.1244046e-03   +4.5502060e-03   +3.7908936e-03   +5.5675185e-03   +5.4467456e-03   +4.5378235e-03   +6.6645017e-03
  +6.0000000e+01   +6.5000000e+01   +5.0577440e-03   +6.9681436e-06  +5.0577440e-03   +3.8254732e-03   +6.7442950e-03   +5.0577440e-03   +4.2137368e-03   +6.1885296e-03   +4.5917121e-03   +3.8254732e-03   +5.6183043e-03   +5.5119586e-03   +4.5921543e-03   +6.7442950e-03
  +6.5000000e+01   +7.0000000e+01   +5.0919700e-03   +7.1534576e-06  +5.0919700e-03   +3.8466815e-03   +6.7993513e-03   +5.0919700e-03   +4.2422511e-03   +6.2304076e-03   +4.6171681e-03   +3.8466815e-03   +5.6494516e-03   +5.5569550e-03   +4.6296417e-03   +6.7993513e-03
  +7.0000000e+01   +7.5000000e+01   +5.1232806e-03   +7.3417290e-06  +5.1232806e-03   +3.8656914e-03   +6.8503980e-03   +5.1232806e-03   +4.2683368e-03   +6.2687188e-03   +4.6399858e-03   +3.8656914e-03   +5.6773709e-03   +5.5986741e-03   +4.6643993e-03   +6.8503980e-03
  +7.5000000e+01   +8.0000000e+01   +5.1321411e-03   +7.4503271e-06  +5.1321411e-03   +3.8676439e-03   +6.8717122e-03   +5.1321411e-03   +4.2757186e-03   +6.2795598e-03   +4.6423292e-03   +3.8676439e-03   +5.6802383e-03   +5.6160940e-03   +4.6789121e-03   +6.8717122e-03
  +8.0000000e+01   +8.5000000e+01   +5.0994260e-03   +7.5363518e-06  +5.0994260e-03   +3.8379196e-03   +6.8381990e-03   +5.0994260e-03   +4.2484630e-03   +6.2395305e-03   +4.6066513e-03   +3.8379196e-03   +5.6365837e-03   +5.5887041e-03   +4.6560928e-03   +6.8381990e-03
  +8.5000000e+01   +9.0000000e+01   +5.1013424e-03   +7.6708338e-06  +5.1013424e-03   +3.8352894e-03   +6.8489835e-03   +5.1013424e-03   +4.2500595e-03   +6.2418753e-03   +4.6034946e-03   +3.8352894e-03   +5.6327212e-03   +5.5975181e-03   +4.6634361e-03   +6.8489835e-03
  +9.0000000e+01   +9.5000000e+01   +5.0768228e-03   +7.7936194e-06  +5.0768228e-03   +3.8120791e-03   +6.8256999e-03   +5.0768228e-03   +4.2296315e-03   +6.2118737e-03   +4.5756353e-03   +3.8120791e-03   +5.5986331e-03   +5.5784889e-03   +4.6475826e-03   +6.8256999e-03
  +9.5000000e+01   +1.0000000e+02   +5.0155915e-03   +7.8077165e-06  +5.0155915e-03   +3.7613292e-03   +6.7530715e-03   +5.0155915e-03   +4.1786181e-03   +6.1369526e-03   +4.5147197e-03   +3.7613292e-03   +5.5240986e-03   +5.5191312e-03   +4.5981302e-03   +6.7530715e-03
  +1.0000000e+02   +1.0500000e+02   +4.9764507e-03   +7.9065321e-06  +4.9764507e-03   +3.7278407e-03   +6.7087886e-03   +4.9764507e-03   +4.1460091e-03   +6.0890611e-03   +4.4745238e-03   +3.7278407e-03   +5.4749155e-03   +5.4829399e-03   +4.5679781e-03   +6.7087886e-03
  +1.0500000e+02   +1.1000000e+02   +4.9242365e-03   +7.9750595e-06  +4.9242365e-03   +3.6845248e-03   +6.6470058e-03   +4.9242365e-03   +4.1025079e-03   +6.0251726e-03   +4.4225318e-03   +3.6845248e-03   +5.4112994e-03   +5.4324463e-03   +4.5259105e-03   +6.6470058e-03
  +1.1000000e+02   +1.1500000e+02   +4.8699877e-03   +8.0469486e-06  +4.8699877e-03   +3.6397146e-03   +6.5823056e-03   +4.8699877e-03   +4.0573121e-03   +5.9587955e-03   +4.3687462e-03   +3.6397146e-03   +5.3454887e-03   +5.3795681e-03   +4.4818564e-03   +6.5823056e-03
  +1.1500000e+02   +1.2000000e+02   +4.7901506e-03   +8.0852982e-06  +4.7901506e-03   +3.5761051e-03   +6.4825813e-03   +4.7901506e-03   +3.9907975e-03   +5.8611088e-03   +4.2923954e-03   +3.5761051e-03   +5.2520680e-03   +5.2980658e-03   +4.4139549e-03   +6.4825813e-03
  +1.2000000e+02   +1.2500000e+02   +4.7174138e-03   +8.1866015e-06  +4.7174138e-03   +3.5178581e-03   +6.3924465e-03   +4.7174138e-03   +3.9301987e-03   +5.7721099e-03   +4.2224821e-03   +3.5178581e-03   +5.1665238e-03   +5.2244007e-03   +4.3525826e-03   +6.3924465e-03
  +1.2500000e+02   +1.3000000e+02   +4.6097171e-03   +8.1713236e-06  +4.6097171e-03   +3.4334234e-03   +6.2548617e-03   +4.6097171e-03   +3.8404739e-03   +5.6403349e-03   +4.1211349e-03   +3.4334234e-03   +5.0425180e-03   +5.1119559e-03   +4.2589018e-03   +6.2548617e-03
  +1.3000000e+02   +1.3500000e+02   +4.5581181e-03   +8.2922185e-06  +4.5581181e-03   +3.3908200e-03   +6.1934032e-03   +4.5581181e-03   +3.7974852e-03   +5.5771997e-03   +4.0699982e-03   +3.3908200e-03   +4.9799481e-03   +5.0617273e-03   +4.2170552e-03   +6.1934032e-03
  +1.3500000e+02   +1.4000000e+02   +4.4700784e-03   +8.4135953e-06  +4.4700784e-03   +3.3224400e-03   +6.0800529e-03   +4.4700784e-03   +3.7241372e-03   +5.4694764e-03   +3.9879215e-03   +3.3224400e-03   +4.8795213e-03   +4.9690886e-03   +4.1398752e-03   +6.0800529e-03
  +1.4000000e+02   +1.4500000e+02   +4.3886613e-03   +8.4329887e-06  +4.3886613e-03   +3.2576451e-03   +5.9782035e-03   +4.3886613e-03   +3.6563067e-03   +5.3698566e-03   +3.9101484e-03   +3.2576451e-03   +4.7843601e-03   +4.8858494e-03   +4.0705267e-03   +5.9782035e-03
  +1.4500000e+02   +1.5000000e+02   +4.3096762e-03   +8.4593385e-06  +4.3096762e-03   +3.1956016e-03   +5.8776637e-03   +4.3096762e-03   +3.5905019e-03   +5.2732125e-03   +3.8356776e-03   +3.1956016e-03   +4.6932395e-03   +4.8036805e-03   +4.0020697e-03   +5.8776637e-03
  +1.5000000e+02   +1.5500000e+02   +4.1992858e-03   +8.5460887e-06  +4.1992858e-03   +3.1104719e-03   +5.7343042e-03   +4.1992858e-03   +3.4985327e-03   +5.1381410e-03   +3.7334964e-03   +3.1104719e-03   +4.5682131e-03   +4.6865162e-03   +3.9044568e-03   +5.7343042e-03
  +1.5500000e+02   +1.6000000e+02   +4.1360996e-03   +8.5991015e-06  +4.1360996e-03   +3.0606006e-03   +5.6543932e-03   +4.1360996e-03   +3.4458907e-03   +5.0608284e-03   +3.6736362e-03   +3.0606006e-03   +4.4949694e-03   +4.6212068e-03   +3.8500461e-03   +5.6543932e-03
  +1.6000000e+02   +1.6500000e+02   +4.0395852e-03   +8.5787321e-06  +4.0395852e-03   +2.9859072e-03   +5.5293995e-03   +4.0395852e-03   +3.3654825e-03   +4.9427358e-03   +3.5839817e-03   +2.9859072e-03   +4.3852707e-03   +4.5190522e-03   +3.7649382e-03   +5.5293995e-03
  +1.6500000e+02   +1.7000000e+02   +3.9350128e-03   +8.6526833e-06  +3.9350128e-03   +2.9055157e-03   +5.3926947e-03   +3.9350128e-03   +3.2783603e-03   +4.8147834e-03   +3.4874880e-03   +2.9055157e-03   +4.2672032e-03   +4.4073264e-03   +3.6718569e-03   +5.3926947e-03
  +1.7000000e+02   +1.7500000e+02   +3.8423279e-03   +8.6109905e-06  +3.8423279e-03   +2.8337871e-03   +5.2728140e-03   +3.8423279e-03   +3.2011420e-03   +4.7013766e-03   +3.4013919e-03   +2.8337871e-03   +4.1618586e-03   +4.3093507e-03   +3.5902308e-03   +5.2728140e-03
  +1.7500000e+02   +1.8000000e+02   +3.7623395e-03   +8.6334220e-06  +3.7623395e-03   +2.7719263e-03   +5.1691953e-03   +3.7623395e-03   +3.1345017e-03   +4.6035049e-03   +3.3271403e-03   +2.7719263e-03   +4.0710061e-03   +4.2246654e-03   +3.5196771e-03   +5.1691953e-03
  +1.8000000e+02   +1.8500000e+02   +3.6587107e-03   +8.5915579e-06  +3.6587107e-03   +2.6928280e-03   +5.0325472e-03   +3.6587107e-03   +3.0481659e-03   +4.4767070e-03   +3.2321989e-03   +2.6928280e-03   +3.9548380e-03   +4.1129860e-03   +3.4266345e-03   +5.0325472e-03
  +1.8500000e+02   +1.9000000e+02   +3.5603915e-03   +8.5327702e-06  +3.5603915e-03   +2.6177443e-03   +4.9034040e-03   +3.5603915e-03   +2.9662534e-03   +4.3564059e-03   +3.1420763e-03   +2.6177443e-03   +3.8445660e-03   +4.0074407e-03   +3.3387014e-03   +4.9034040e-03
  +1.9000000e+02   +1.9500000e+02   +3.4836323e-03   +8.6104310e-06  +3.4836323e-03   +2.5583954e-03   +4.8038223e-03   +3.4836323e-03   +2.9023037e-03   +4.2624856e-03   +3.0708396e-03   +2.5583954e-03   +3.7574030e-03   +3.9260544e-03   +3.2708967e-03   +4.8038223e-03
  +1.9500000e+02   +2.0000000e+02   +3.4025091e-03   +8.6866168e-06  +3.4025091e-03   +2.4963983e-03   +4.6971398e-03   +3.4025091e-03   +2.8347176e-03   +4.1632251e-03   +2.9964241e-03   +2.4963983e-03   +3.6663497e-03   +3.8388652e-03   +3.1982574e-03   +4.6971398e-03
  +2.0000000e+02   +2.0500000e+02   +3.3027629e-03   +8.5411238e-06  +3.3027629e-03   +2.4206411e-03   +4.5650854e-03   +3.3027629e-03   +2.7516166e-03   +4.0411785e-03   +2.9054933e-03   +2.4206411e-03   +3.5550890e-03   +3.7309402e-03   +3.1083420e-03   +4.5650854e-03
  +2.0500000e+02   +2.1000000e+02   +3.2144729e-03   +8.5754834e-06  +3.2144729e-03   +2.3533323e-03   +4.4487851e-03   +3.2144729e-03   +2.6780599e-03   +3.9331489e-03   +2.8247025e-03   +2.3533323e-03   +3.4562357e-03   +3.6358907e-03   +3.0291539e-03   +4.4487851e-03
  +2.1000000e+02   +2.1500000e+02   +3.1133418e-03   +8.4724596e-06  +3.1133418e-03   +2.2771205e-03   +4.3134000e-03   +3.1133418e-03   +2.5938049e-03   +3.8094072e-03   +2.7332255e-03   +2.2771205e-03   +3.3443066e-03   +3.5252435e-03   +2.9369709e-03   +4.3134000e-03
  +2.1500000e+02   +2.2000000e+02   +3.0341455e-03   +8.5797753e-06  +3.0341455e-03   +2.2171072e-03   +4.2084712e-03   +3.0341455e-03   +2.5278246e-03   +3.7125048e-03   +2.6611915e-03   +2.2171072e-03   +3.2561676e-03   +3.4394875e-03   +2.8655253e-03   +4.2084712e-03
  +2.2000000e+02   +2.2500000e+02   +2.9541141e-03   +8.4227412e-06  +2.9541141e-03   +2.1561855e-03   +4.1026627e-03   +2.9541141e-03   +2.4611481e-03   +3.6145801e-03   +2.5880677e-03   +2.1561855e-03   +3.1666950e-03   +3.3530129e-03   +2.7934810e-03   +4.1026627e-03
  +2.2500000e+02   +2.3000000e+02   +2.8651266e-03   +8.3854722e-06  +2.8651266e-03   +2.0890377e-03   +3.9839130e-03   +2.8651266e-03   +2.3870106e-03   +3.5056973e-03   +2.5074698e-03   +2.0890377e-03   +3.0680775e-03   +3.2559616e-03   +2.7126250e-03   +3.9839130e-03
  +2.3000000e+02   +2.3500000e+02   +2.7720638e-03   +8.3255728e-06  +2.7720638e-03   +2.0195519e-03   +3.8580587e-03   +2.7720638e-03   +2.3094775e-03   +3.3918281e-03   +2.4240662e-03   +2.0195519e-03   +2.9660273e-03   +3.1531033e-03   +2.6269312e-03   +3.8580587e-03
  +2.3500000e+02   +2.4000000e+02   +2.7040425e-03   +8.3702030e-06  +2.7040425e-03   +1.9681078e-03   +3.7676890e-03   +2.7040425e-03   +2.2528074e-03   +3.3085990e-03   +2.3623181e-03   +1.9681078e-03   +2.8904734e-03   +3.0792464e-03   +2.5653992e-03   +3.7676890e-03
  +2.4000000e+02   +2.4500000e+02   +2.6212538e-03   +9.0108116e-06  +2.6212538e-03   +1.9057806e-03   +3.6570121e-03   +2.6212538e-03   +2.1838338e-03   +3.2073007e-03   +2.2875067e-03   +1.9057806e-03   +2.7989361e-03   +2.9887929e-03   +2.4900398e-03   +3.6570121e-03
  +2.4500000e+02   +2.5000000e+02   +2.5478036e-03   +8.4167189e-06  +2.5478036e-03   +1.8505782e-03   +3.5585019e-03   +2.5478036e-03   +2.1226406e-03   +3.1174290e-03   +2.2212472e-03   +1.8505782e-03   +2.7178628e-03   +2.9082822e-03   +2.4229648e-03   +3.5585019e-03
  +2.5000000e+02   +2.5500000e+02   +2.4794728e-03   +8.2861750e-06  +2.4794728e-03   +1.7991867e-03   +3.4669055e-03   +2.4794728e-03   +2.0657128e-03   +3.0338214e-03   +2.1595621e-03   +1.7991867e-03   +2.6423860e-03   +2.8334227e-03   +2.3605974e-03   +3.4669055e-03
  +2.5500000e+02   +2.6000000e+02   +2.4038050e-03   +8.3685086e-06  +2.4038050e-03   +1.7427352e-03   +3.3645732e-03   +2.4038050e-03   +2.0026717e-03   +2.9412358e-03   +2.0918032e-03   +1.7427352e-03   +2.5594784e-03   +2.7497888e-03   +2.2909197e-03   +3.3645732e-03
  +2.6000000e+02   +2.6500000e+02   +2.3320424e-03   +8.1099535e-06  +2.3320424e-03   +1.6889581e-03   +3.2679826e-03   +2.3320424e-03   +1.9428844e-03   +2.8534289e-03   +2.0272547e-03   +1.6889581e-03   +2.4804982e-03   +2.6708477e-03   +2.2251518e-03   +3.2679826e-03
  +2.6500000e+02   +2.7000000e+02   +2.2494109e-03   +8.2796862e-06  +2.2494109e-03   +1.6279012e-03   +3.1550947e-03   +2.2494109e-03   +1.8740420e-03   +2.7523230e-03   +1.9539682e-03   +1.6279012e-03   +2.3908266e-03   +2.5785869e-03   +2.1482872e-03   +3.1550947e-03
  +2.7000000e+02   +2.7500000e+02   +2.1762624e-03   +7.9793659e-06  +2.1762624e-03   +1.5729217e-03   +3.0568790e-03   +2.1762624e-03   +1.8131002e-03   +2.6628203e-03   +1.8879765e-03   +1.5729217e-03   +2.3100811e-03   +2.4983176e-03   +2.0814125e-03   +3.0568790e-03
  +2.7500000e+02   +2.8000000e+02   +2.1247961e-03   +8.3092820e-06  +2.1247961e-03   +1.5348899e-03   +2.9867247e-03   +2.1247961e-03   +1.7702221e-03   +2.5998473e-03   +1.8423268e-03   +1.5348899e-03   +2.2542254e-03   +2.4409820e-03   +2.0336449e-03   +2.9867247e-03
  +2.8000000e+02   +2.8500000e+02   +2.0403896e-03   +7.8579985e-06  +2.0403896e-03   +1.4723157e-03   +2.8714210e-03   +2.0403896e-03   +1.6999011e-03   +2.4965698e-03   +1.7672191e-03   +1.4723157e-03   +2.1623252e-03   +2.3467470e-03   +1.9551352e-03   +2.8714210e-03
  +2.8500000e+02   +2.9000000e+02   +1.9868560e-03   +7.8268596e-06  +1.9868560e-03   +1.4324034e-03   +2.7990329e-03   +1.9868560e-03   +1.6553008e-03   +2.4310675e-03   +1.7193126e-03   +1.4324034e-03   +2.1037078e-03   +2.2875860e-03   +1.9058465e-03   +2.7990329e-03
  +2.9000000e+02   +2.9500000e+02   +1.9087428e-03   +7.8532152e-06  +1.9087428e-03   +1.3751557e-03   +2.6912766e-03   +1.9087428e-03   +1.5902228e-03   +2.3354902e-03   +1.6505982e-03   +1.3751557e-03   +2.0196309e-03   +2.1995192e-03   +1.8324760e-03   +2.6912766e-03
  +2.9500000e+02   +3.0000000e+02   +1.8605686e-03   +7.7129538e-06  +1.8605686e-03   +1.3388104e-03   +2.6268953e-03   +1.8605686e-03   +1.5500874e-03   +2.2765450e-03   +1.6069728e-03   +1.3388104e-03   +1.9662519e-03   +2.1469020e-03   +1.7886388e-03   +2.6268953e-03
  +3.0000000e+02   +3.0500000e+02   +1.8123577e-03   +7.6721118e-06  +1.8123577e-03   +1.3027954e-03   +2.5618408e-03   +1.8123577e-03   +1.5099216e-03   +2.2175553e-03   +1.5637441e-03   +1.3027954e-03   +1.9133584e-03   +2.0937340e-03   +1.7443436e-03   +2.5618408e-03
  +3.0500000e+02   +3.1000000e+02   +1.7553937e-03   +7.7358353e-06  +1.7553937e-03   +1.2606583e-03   +2.4840099e-03   +1.7553937e-03   +1.4624635e-03   +2.1478560e-03   +1.5131669e-03   +1.2606583e-03   +1.8514733e-03   +2.0301248e-03   +1.6913491e-03   +2.4840099e-03
  +3.1000000e+02   +3.1500000e+02   +1.6836013e-03   +7.5019828e-06  +1.6836013e-03   +1.2083641e-03   +2.3841776e-03   +1.6836013e-03   +1.4026516e-03   +2.0600125e-03   +1.4503984e-03   +1.2083641e-03   +1.7746712e-03   +1.9485340e-03   +1.6233737e-03   +2.3841776e-03
  +3.1500000e+02   +3.2000000e+02   +1.6568712e-03   +7.5658482e-06  +1.6568712e-03   +1.1879141e-03   +2.3492177e-03   +1.6568712e-03   +1.3803820e-03   +2.0273063e-03   +1.4258522e-03   +1.1879141e-03   +1.7446372e-03   +1.9199620e-03   +1.5995696e-03   +2.3492177e-03
  +3.2000000e+02   +3.2500000e+02   +1.5764572e-03   +7.5665153e-06  +1.5764572e-03   +1.1292534e-03   +2.2376257e-03   +1.5764572e-03   +1.3133869e-03   +1.9289136e-03   +1.3554418e-03   +1.1292534e-03   +1.6584849e-03   +1.8287603e-03   +1.5235873e-03   +2.2376257e-03
  +3.2500000e+02   +3.3000000e+02   +1.5492364e-03   +7.5192251e-06  +1.5492364e-03   +1.1088532e-03   +2.2010520e-03   +1.5492364e-03   +1.2907087e-03   +1.8956070e-03   +1.3309555e-03   +1.1088532e-03   +1.6285238e-03   +1.7988694e-03   +1.4986844e-03   +2.2010520e-03
  +3.3000000e+02   +3.3500000e+02   +1.4873961e-03   +7.4592054e-06  +1.4873961e-03   +1.0635263e-03   +2.1155608e-03   +1.4873961e-03   +1.2391879e-03   +1.8199405e-03   +1.2765496e-03   +1.0635263e-03   +1.5619542e-03   +1.7289998e-03   +1.4404742e-03   +2.1155608e-03
  +3.3500000e+02   +3.4000000e+02   +1.4375520e-03   +7.3429753e-06  +1.4375520e-03   +1.0270315e-03   +2.0466681e-03   +1.4375520e-03   +1.1976615e-03   +1.7589527e-03   +1.2327450e-03   +1.0270315e-03   +1.5083560e-03   +1.6726953e-03   +1.3935654e-03   +2.0466681e-03
  +3.4000000e+02   +3.4500000e+02   +1.3970554e-03   +7.1815353e-06  +1.3970554e-03   +9.9731331e-04   +1.9907846e-03   +1.3970554e-03   +1.1639227e-03   +1.7094019e-03   +1.1970742e-03   +9.9731331e-04   +1.4647102e-03   +1.6270229e-03   +1.3555146e-03   +1.9907846e-03
  +3.4500000e+02   +3.5000000e+02   +1.3507987e-03   +7.1396889e-06  +1.3507987e-03   +9.6356861e-04   +1.9265667e-03   +1.3507987e-03   +1.1253851e-03   +1.6528036e-03   +1.1565705e-03   +9.6356861e-04   +1.4151508e-03   +1.5745389e-03   +1.3117889e-03   +1.9265667e-03
  +3.5000000e+02   +3.5500000e+02   +1.3102727e-03   +7.1157710e-06  +1.3102727e-03   +9.3339348e-04   +1.8716569e-03   +1.3102727e-03   +1.0916219e-03   +1.6032171e-03   +1.1203514e-03   +9.3339348e-04   +1.3708339e-03   +1.5296626e-03   +1.2744012e-03   +1.8716569e-03
  +3.5500000e+02   +3.6000000e+02   +1.2487536e-03   +6.9471746e-06  +1.2487536e-03   +8.8908421e-04   +1.7849180e-03   +1.2487536e-03   +1.0403687e-03   +1.5279437e-03   +1.0671670e-03   +8.8908421e-04   +1.3057589e-03   +1.4587729e-03   +1.2153412e-03   +1.7849180e-03
  +3.6000000e+02   +3.6500000e+02   +1.2251629e-03   +7.2862465e-06  +1.2251629e-03   +8.7146998e-04   +1.7531673e-03   +1.2251629e-03   +1.0207147e-03   +1.4990788e-03   +1.0460246e-03   +8.7146998e-04   +1.2798896e-03   +1.4328237e-03   +1.1937223e-03   +1.7531673e-03
  +3.6500000e+02   +3.7000000e+02   +1.1995191e-03   +7.1682941e-06  +1.1995191e-03   +8.5221479e-04   +1.7188878e-03   +1.1995191e-03   +9.9935015e-04   +1.4677016e-03   +1.0229126e-03   +8.5221479e-04   +1.2516103e-03   +1.4048078e-03   +1.1703815e-03   +1.7188878e-03
  +3.7000000e+02   +3.7500000e+02   +1.1496020e-03   +6.7611901e-06  +1.1496020e-03   +8.1636823e-04   +1.6481216e-03   +1.1496020e-03   +9.5776295e-04   +1.4066243e-03   +9.7988605e-04   +8.1636823e-04   +1.1989641e-03   +1.3469722e-03   +1.1221971e-03   +1.6481216e-03
  +3.7500000e+02   +3.8000000e+02   +1.1038304e-03   +7.1912606e-06  +1.1038304e-03   +7.8319826e-04   +1.5842025e-03   +1.1038304e-03   +9.1962948e-04   +1.3506194e-03   +9.4007216e-04   +7.8319826e-04   +1.1502488e-03   +1.2947326e-03   +1.0786750e-03   +1.5842025e-03
  +3.8000000e+02   +3.8500000e+02   +1.0790691e-03   +6.7730829e-06  +1.0790691e-03   +7.6518708e-04   +1.5496057e-03   +1.0790691e-03   +8.9900015e-04   +1.3203220e-03   +9.1845335e-04   +7.6518708e-04   +1.1237966e-03   +1.2664573e-03   +1.0551181e-03   +1.5496057e-03
  +3.8500000e+02   +3.9000000e+02   +1.0365015e-03   +6.6659012e-06  +1.0365015e-03   +7.3416493e-04   +1.4904933e-03   +1.0365015e-03   +8.6353605e-04   +1.2682374e-03   +8.8121751e-04   +7.3416493e-04   +1.0782358e-03   +1.2181462e-03   +1.0148689e-03   +1.4904933e-03
  +3.9000000e+02   +3.9500000e+02   +9.9070234e-04   +6.5221905e-06  +9.9070234e-04   +7.0111383e-04   +1.4260718e-03   +9.9070234e-04   +8.2537957e-04   +1.2121987e-03   +8.4154630e-04   +7.0111383e-04   +1.0296950e-03   +1.1654960e-03   +9.7100466e-04   +1.4260718e-03
  +3.9500000e+02   +4.0000000e+02   +9.7382763e-04   +6.5099199e-06  +9.7382763e-04   +6.8885514e-04   +1.4025471e-03   +9.7382763e-04   +8.1132081e-04   +1.1915512e-03   +8.2683221e-04   +6.8885514e-04   +1.0116913e-03   +1.1462698e-03   +9.5498678e-04   +1.4025471e-03
  +4.0000000e+02   +4.0500000e+02   +9.4426341e-04   +6.7093163e-06  +9.4426341e-04   +6.6748188e-04   +1.3610914e-03   +9.4426341e-04   +7.8669011e-04   +1.1553772e-03   +8.0117791e-04   +6.6748188e-04   +9.8030129e-04   +1.1123890e-03   +9.2675984e-04   +1.3610914e-03
  +4.0500000e+02   +4.1000000e+02   +9.1128636e-04   +6.4455321e-06  +9.1128636e-04   +6.4361440e-04   +1.3148145e-03   +9.1128636e-04   +7.5921607e-04   +1.1150273e-03   +7.7252977e-04   +6.4361440e-04   +9.4524817e-04   +1.0745679e-03   +8.9525014e-04   +1.3148145e-03
  +4.1000000e+02   +4.1500000e+02   +8.8231876e-04   +6.3087693e-06  +8.8231876e-04   +6.2285998e-04   +1.2737321e-03   +8.8231876e-04   +7.3508241e-04   +1.0795832e-03   +7.4761827e-04   +6.2285998e-04   +9.1476707e-04   +1.0409922e-03   +8.6727738e-04   +1.2737321e-03
  +4.1500000e+02   +4.2000000e+02   +8.4457237e-04   +6.2096297e-06  +8.4457237e-04   +5.9529157e-04   +1.2214138e-03   +8.4457237e-04   +7.0363491e-04   +1.0333977e-03   +7.1452790e-04   +5.9529157e-04   +8.7427853e-04   +9.9823366e-04   +8.3165411e-04   +1.2214138e-03
  +4.2000000e+02   +4.2500000e+02   +8.3101327e-04   +6.2782624e-06  +8.3101327e-04   +5.8537090e-04   +1.2026591e-03   +8.3101327e-04   +6.9233851e-04   +1.0168071e-03   +7.0262014e-04   +5.8537090e-04   +8.5970850e-04   +9.8290587e-04   +8.1888412e-04   +1.2026591e-03
  +4.2500000e+02   +4.3000000e+02   +8.0219397e-04   +6.2362833e-06  +8.0219397e-04   +5.6482374e-04   +1.1616481e-03   +8.0219397e-04   +6.6832842e-04   +9.8154455e-04   +6.7795741e-04   +5.6482374e-04   +8.2953177e-04   +9.4938850e-04   +7.9095992e-04   +1.1616481e-03
  +4.3000000e+02   +4.3500000e+02   +7.6761926e-04   +6.0964481e-06  +7.6761926e-04   +5.4017310e-04   +1.1122971e-03   +7.6761926e-04   +6.3952337e-04   +9.3923982e-04   +6.4836925e-04   +5.4017310e-04   +7.9332845e-04   +9.0905505e-04   +7.5735711e-04   +1.1122971e-03
  +4.3500000e+02   +4.4000000e+02   +7.4577257e-04   +6.1715920e-06  +7.4577257e-04   +5.2407632e-04   +1.0823423e-03   +7.4577257e-04   +6.2132230e-04   +9.1250873e-04   +6.2904832e-04   +5.2407632e-04   +7.6968781e-04   +8.8457367e-04   +7.3696107e-04   +1.0823423e-03
  +4.4000000e+02   +4.4500000e+02   +7.2127737e-04   +5.8915802e-06  +7.2127737e-04   +5.0659575e-04   +1.0474222e-03   +7.2127737e-04   +6.0091470e-04   +8.8253698e-04   +6.0806640e-04   +5.0659575e-04   +7.4401490e-04   +8.5603426e-04   +7.1318413e-04   +1.0474222e-03
  +4.4500000e+02   +4.5000000e+02   +7.0170293e-04   +5.9360840e-06  +7.0170293e-04   +4.9229766e-04   +1.0202431e-03   +7.0170293e-04   +5.8460672e-04   +8.5858616e-04   +5.9090446e-04   +4.9229766e-04   +7.2301596e-04   +8.3382143e-04   +6.9467807e-04   +1.0202431e-03
  +4.5000000e+02   +4.5500000e+02   +6.7357932e-04   +5.8595541e-06  +6.7357932e-04   +4.7242496e-04   +9.7974214e-04   +6.7357932e-04   +5.6117623e-04   +8.2417488e-04   +5.6705125e-04   +4.7242496e-04   +6.9382976e-04   +8.0072089e-04   +6.6710117e-04   +9.7974214e-04
  +4.5500000e+02   +4.6000000e+02   +6.4676950e-04   +5.7997979e-06  +6.4676950e-04   +4.5319978e-04   +9.4177341e-04   +6.4676950e-04   +5.3884028e-04   +7.9137103e-04   +5.4397530e-04   +4.5319978e-04   +6.6559461e-04   +7.6968994e-04   +6.4124846e-04   +9.4177341e-04
  +4.6000000e+02   +4.6500000e+02   +6.3395643e-04   +5.7743803e-06  +6.3395643e-04   +4.4370575e-04   +9.2436890e-04   +6.3395643e-04   +5.2816539e-04   +7.7569328e-04   +5.3257961e-04   +4.4370575e-04   +6.5165109e-04   +7.5546566e-04   +6.2939785e-04   +9.2436890e-04
  +4.6500000e+02   +4.7000000e+02   +6.1793504e-04   +5.7037929e-06  +6.1793504e-04   +4.3248833e-04   +9.0100934e-04   +6.1793504e-04   +5.1481756e-04   +7.5608991e-04   +5.1911536e-04   +4.3248833e-04   +6.3517659e-04   +7.3637441e-04   +6.1349242e-04   +9.0100934e-04
  +4.7000000e+02   +4.7500000e+02   +5.9361882e-04   +5.7607316e-06  +5.9361882e-04   +4.1506736e-04   +8.6660723e-04   +5.9361882e-04   +4.9455909e-04   +7.2633718e-04   +4.9820498e-04   +4.1506736e-04   +6.0959116e-04   +7.0825833e-04   +5.9006822e-04   +8.6660723e-04
  +4.7500000e+02   +4.8000000e+02   +5.7923466e-04   +5.6646912e-06  +5.7923466e-04   +4.0476363e-04   +8.4614838e-04   +5.7923466e-04   +4.8257527e-04   +7.0873708e-04   +4.8583743e-04   +4.0476363e-04   +5.9445854e-04   +6.9153778e-04   +5.7613788e-04   +8.4614838e-04
  +4.8000000e+02   +4.8500000e+02   +5.5908630e-04   +5.4724774e-06  +5.5908630e-04   +3.9024532e-04   +8.1780068e-04   +5.5908630e-04   +4.6578913e-04   +6.8408407e-04   +4.6841111e-04   +3.9024532e-04   +5.7313614e-04   +6.6836985e-04   +5.5683609e-04   +8.1780068e-04
  +4.8500000e+02   +4.9000000e+02   +5.4234703e-04   +5.5646541e-06  +5.4234703e-04   +3.7841873e-04   +7.9363291e-04   +5.4234703e-04   +4.5184324e-04   +6.6360229e-04   +4.5421566e-04   +3.7841873e-04   +5.5576696e-04   +6.4861809e-04   +5.4038036e-04   +7.9363291e-04
  +4.9000000e+02   +4.9500000e+02   +5.2294282e-04   +5.4102488e-06  +5.2294282e-04   +3.6462273e-04   +7.6586331e-04   +5.2294282e-04   +4.3567712e-04   +6.3985979e-04   +4.3765632e-04   +3.6462273e-04   +5.3550536e-04   +6.2592260e-04   +5.2147222e-04   +7.6586331e-04
  +4.9500000e+02   +5.0000000e+02   +5.1251674e-04   +5.6016573e-06  +5.1251674e-04   +3.5692245e-04   +7.5166427e-04   +5.1251674e-04   +4.2699086e-04   +6.2710269e-04   +4.2841369e-04   +3.5692245e-04   +5.2419632e-04   +6.1431805e-04   +5.1180414e-04   +7.5166427e-04
<\histogram>

EOF
