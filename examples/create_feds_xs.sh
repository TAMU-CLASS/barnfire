#! /usr/bin/env bash
# A. Till Feb 2016

############################
# INITIALIZE VARIABLES
############################

## !!!!!
# You need to set and export the following:
# export SCRATCH_BARN=...
# export ENDF=...
# export NJOY=...
## !!!!!

# This is a resolution parameter used in indicators.py only. I usually just set it at 9.
res=9

# In addition to coarse groups, I use energy penalties.
# See my dissertation, Appendix B, section 1 for more details
penalty=0.96

# Which problem to run, materials-wise.
# This is only used in indicators.py and indicators_clustering.py,
# where there is a mapping from this string to a set of materials defined in
# materials_materials.py.
# This particular m specifies the 2D pincell set of problems (with only U-238 as an indicator)
# (This is a somewhat complicated parameter that I should simplify.)
m=kpin2d

# The specific list of materials to generate NJOY inputs for / read NJOY outputs from
# This is only used in calls to materials.py
mList=(kFUEL kMOD kRFUEL kREFUEL kRMFUEL)

# How many groups to use in the resolved resonance region (RRR).
g=10 #30

# These are the coarse-group boundaries in eV. c also specifies the extent of the RRR
c=(4.21983E+0 9.50002E+0 2.78852E+1 5.17847E+1 1.54176E+2 5.39204E+2 2.08410E+3 9.11881E+3 2.26994E+4)

# The group structure to use outside the RRR in the thermal and fast energy ranges.
# G=shem-361 would point to ../dat/energy_groups/shem-361.txt
G=shem-361

# How many Legendre moments to use for the scattering transfer matrices
L=3

# Which clustering algorithm to use. Recommended: use 'tmg' for MG XS and 'har' for FEDS XS
# Use './indicators_clustering.py -h' for more information.
clusterer=har

# Which apportioning algorithm to use. Recommended: use 'var' for FEDS and 'L1' or 'max' for MG XS
# Use './indicators_clustering.py -h' for more information.
appt=var

# Which reactions to include in the PDT XS file. Recommended: use 'abs' if absorption edits are needed
# Use './materials.py -h' for more information
rxnOpt=abs

############################
# GENERATE CROSS SECTIONS
############################

#NB: In many cases, if a previous step has run successfully, you don't need
# to rerun it to run a later step if you make changes that affect that later
# step only.

# Step 0: Initialize your scratch directory
scriptdir=`pwd`
srcdir=../src
cd $srcdir
srcdir=`pwd`
./Initialize.py $scriptdir $0
# Step 1: Generate the NJOY inputs
./materials.py -m ${mList[*]} -v
# Step 2: Run NJOY to generate the PENDF files
cd $SCRATCH_BARN/xs
./RunPendf.sh
cd $srcdir
# Step 3: Generate the indicators and their energy grid
./indicators.py -R ${c[0]} ${c[*]: -1} -m $m -w flux -r $res -G $G -v
# Step 4: Generate more indicators on the same energy grid but using escape XS
./indicators.py -R ${c[0]} ${c[*]: -1} -m $m -w fluxe -r $res -v
# Step 5: Cluster the indicators to get the 'clust' file, which specifies the energy mesh
./indicators_clustering.py -c ${c[*]} -e $g -m $m -w $clusterer -E $penalty -r $res -a $appt
# Step 6: Generate the weighting spectrum on the subelements
./indicators.py -R ${c[0]} ${c[*]: -1} -m $m -w wgt -r $res -M clust-$g-$res
# Step 7: Regenerate NJOY inputs with a group structure of the subelements
./materials.py -m ${mList[*]} -G clust-$g-$res -L $L -v
# Step 8: Run NJOY to generate GENDF files (on the subelements)
cd $SCRATCH_BARN/xs
./RunGendf.sh
# Step 9: Sum/average over the subelements to get PDT-formatted FEDS XS on the elements
cd $srcdir
./materials.py -b -m ${mList[*]} -p $rxnOpt -M clust-$g-$res -v
cd $scriptdir
# The result should be a list of .data files in $SCRATCH_BARN/xs/pdtxs

