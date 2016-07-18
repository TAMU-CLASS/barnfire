#! /usr/bin/env bash
# A. Till Jul 2016

##################
# SET DIRECTORIES
##################

# You need to set and export the following:
# export SCRATCH_BARN=...
# export ENDF=...
# export NJOY=...

##################################
# SET PROBLEM-DEPENDENT VARIABLES
##################################

# The specific list of materials for which to generate NJOY inputs
# Used in calls to materials.py and when weighting spectra are calculated in the last call to indicators.py
mList=(mtTFUEL_0 mtTFUEL_1 mtTFUEL_2 mtTFUEL_3 mtTFUEL_4 mtTFUEL_5 mtTFUEL_6 mtTFUEL_7 mtH2O_0 mtH2O_1)

# The specific list of materials to use as the indicators when determining the group structure
# indList should be a subset of mList
indList=(mtTFUEL_0 mtTFUEL_7)

# The importances are the relative weightings of the indicators when clustering
# The order corresponds to the materials in indList
impList=(1 2)

# How many groups to use in the resolved resonance region (RRR).
g=30

# Total groups used (sum of fast, RRR, and thermal groups)
#gTotal=65 #35+$g
gTotal=166 #136+$g

# The group structure to use outside the RRR in the thermal and fast energy ranges.
# G=shem-361 would point to ../dat/energy_groups/shem-361.txt
#G=scale-44-tweaked
G=shem-361

#############################
# SET LESS-CHANGED VARIABLES
#############################

# These are the coarse-group boundaries in eV. c also specifies the extent of the RRR
c=(4.21983E+0 9.50002E+0 2.78852E+1 5.17847E+1 1.54176E+2 5.39204E+2 2.08410E+3 9.11881E+3 2.26994E+4)

# How many Legendre moments to use for the scattering transfer matrices
L=1

# Which reactions to include in the PDT XS file. Recommended: use 'abs' if absorption edits are needed
# Use './materials.py -h' for more information
rxnOpt=abs

##############################
# SET RARELY-CHANGED VARIABLES
##############################

# Which clustering algorithm to use. Recommended: use 'tmg' for MG XS and 'har' for FEDS XS
# Use './indicators_clustering.py -h' for more information.
clusterer=har

# Which apportioning algorithm to use. Recommended: use 'var' for FEDS and 'L1' or 'max' for MG XS
# Use './indicators_clustering.py -h' for more information.
appt=var

# This is a resolution parameter used in indicators.py only. I usually just set it at 9.
res=9

# In addition to coarse groups, I use energy penalties.
# See my dissertation, Appendix B, section 1 for more details
penalty=0.96

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

# Step 3: Generate the indicators and their energy grid
# (Add -t option to force temperature dependence)
cd $srcdir
./indicators.py -R ${c[0]} ${c[*]: -1} -m manual -i ${indList[*]} -t -w flux -r $res -G $G -v 
# Step 4: Generate more indicators on the same energy grid but using escape XS
./indicators.py -R ${c[0]} ${c[*]: -1} -m manual -i ${indList[*]} -t -w fluxe -r $res -p -v
# Step 5: Cluster the indicators to get the 'clust' file, which specifies the energy mesh
./indicators_clustering.py -c ${c[*]} -e $g -m manual -i ${indList[*]} -I ${impList[*]} -w $clusterer -E $penalty -r $res -a $appt -s 1
# Step 6: Generate the weighting spectrum on the subelements
./indicators.py -R ${c[0]} ${c[*]: -1} -m manual -i ${mList[*]} -t -w wgt -r $res -M clust-$g-$res
# Step 7: Regenerate NJOY inputs with a group structure of the subelements
./materials.py -m ${mList[*]} -G clust-$g-$res -L $L -v
# Step 8: Run NJOY to generate GENDF files (on the subelements)
cd $SCRATCH_BARN/xs
./RunGendf.sh
# Step 9: Sum/average over the subelements to get PDT-formatted FEDS XS on the elements
cd $srcdir
./materials.py -b -m ${mList[*]} -p $rxnOpt -M clust-$g-$res -v
# Step 10: Create multi-temperature PDT-formatted FEDS XS on the elements
./mergeXS.py -m ${mList[*]} -g $gTotal -v
# The result should be a list of .data files in $SCRATCH_BARN/xs/pdtxs

# Step 11: Run NJOY to generate the ACE files (can be done after step 2)
cd $SCRATCH_BARN/xs
./RunAce.sh
# Step 12: Copy ACE files and create xsdir
cd ace/xdata
./copyAce.sh
cd $scriptdir
# The result should be ACE files and xsdir (tells MCNP where the ACE files are) in $SCRATCH_BARN/xs/ace/xdata
exit

# Addendum:
# If there is an error with step 6 (likely caused by temperature interpolation of PENDF files,
# which is not supported in indicators.py), try the following:
# (No -t because interpolation not currently supported)
#./indicators.py -R ${c[0]} ${c[*]: -1} -m manual -i ${mList[*]} -w wgt -r $res -M clust-$g-$res
# Redo for materials that need correct temperatures
#./indicators.py -R ${c[0]} ${c[*]: -1} -m manual -i ${indList[*]} -t -w wgt -r $res -M clust-$g-$res

