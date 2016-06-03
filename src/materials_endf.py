'''
Andrew Till
Summer 2014

Functions for processing ENDF files

As of 5/31/2016,LANL's T2 nuclear data (http://t2.lanl.gov/nis/data/endf/index.html) had
ENDF/B-VII.1 neutron data but only ENDF/B-VII.0 bound thermal neutron data.
BNL's non-thermal neutron data files can sometimes be bad (will produce garbage when processed with NJOY).
The code uses the most up-to-date data while requiring consistency.
Hence, automatic downloads of the non-thermal neutron data from LANL T2
and automatic downloads of the bound-thermal neutron data from BNL.
A copy of the bound-thermal neutron data from BNL is included locally in case BNL changes its URL's.
'''

#TODO: Implement automatic downloads of the bound-thermal neutron data from BNL

#MINE
from directories import get_common_directories, copy_thermal_endf_xs, make_endf_directory

#STDLIB
import os
import sys
import shutil
from subprocess import check_output
sys.path.insert(1, get_common_directories()['nuclideData'])

#TPL
import nuclide_data as nd

#MINE
from Readgroupr import get_pot_scat_xs
import materials_util as util

###############################################################################
def get_endf_files(ZAList, dirDict, potXSDict, useBXS=False, verbosity=False):
    if verbosity > 1:
        print '------- ENDF Files -------'
    make_endf_directory()
    endfDirr = dirDict['endf']
    for (Z, A) in sorted(ZAList):
        get_endf_file(Z, A, endfDirr, verbosity)
        if useBXS:
            potXSDict[(Z,A)] = get_pot_scat_xs_wrapper(Z, A, endfDirr, verbosity)
        else:
            potXSDict[(Z,A)] = 0.0
    if verbosity > 1:
        print 'potXSDict', potXSDict

def get_endf_file(Z, A, endfDirr, verbosity=False):
    '''Download ENDF file from LANL T2'''
    sym = nd.z2sym[Z]
    outDirr = endfDirr
    metastableStr = ''
    # Metastable isomeric states use the groundstate A + 400
    effA = A % 400
    if A // 400 > 0:
        metastableStr = 'm'
    outname = 'endf_{0:02g}{1:03g}{2}_vii1'.format(Z, effA, metastableStr)
    outPath = os.path.join(outDirr, outname)
    if not(os.path.isfile(outPath)):
        if A == 0:
            effA = 'nat'
        metastableStr = ''
        if A // 400 > 0:
            metastableStr = 'm1'
        webLocation = 'http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/{0}/{1}{2}'.format(sym, effA, metastableStr)
        if verbosity:
            print 'Downloading {0} from {1} to {2}'.format(outname, webLocation, outDirr)
        os.system('wget {0} -O {1}'.format(webLocation, outPath))
    elif verbosity > 2:
        print '{0} already downloaded in {1}'.format(outname, outDirr)

def get_thermal_endf_files(ZASabList, dirDict, verbosity=False):
    if verbosity > 1:
        print '------- Bound Thermal ENDF Files -------'
    # Copy all bound thermal ENDF XS that exist locally (in ../dat) to the ENDF folder
    copy_thermal_endf_xs()
    # Download any missing bound thermal ENDF XS to the ENDF folder
    endfDirr = dirDict['endf']
    ZSabSet = {(Z,Sab) for (Z,A,Sab) in ZASabList}
    for (Z, Sab) in sorted(ZSabSet):
        get_thermal_endf_file(Z, Sab, endfDirr, verbosity)
        
def get_thermal_endf_file(Z, Sab, endfDirr, verbosity=False):
    '''Download bound thermal ENDF file from BNL'''
    if Sab in util.get_non_bound_names():
        return
    Sab2thermalName = util.get_element_thermal_name_to_thermal_name_dict()
    Sab2bnl = util.get_element_thermal_name_to_bnl_id_dict()
    bnlID = Sab2bnl[Sab]
    outDirr = endfDirr
    outname = util.format_thermal_filename()(Z=Z, thermalName=Sab2thermalName[Sab])
    outPath = os.path.join(outDirr, outname)
    tempOutPath = '{0}.tmp'.format(outPath)
    if not(os.path.isfile(outPath)):
        webLocation = 'http://www.nndc.bnl.gov/sigma/getDataset.jsp?evalid={}'.format(bnlID)
        if verbosity:
            print 'Downloading {0} from {1} to {2}'.format(outname, webLocation, outDirr)
        os.system('wget {0} -O {1}'.format(webLocation, tempOutPath))
        # BNL returns an HTML file. Remove first two lines and last line to get the ENDF file
        numLines = int(check_output(["wc", "-l", tempOutPath]).split()[0])
        os.system('head  -n {0} {1} | tail -n {2} > {3}'.format(
            numLines-1, tempOutPath, numLines-3, outPath))
        os.remove('{0}'.format(tempOutPath))
    elif verbosity > 2:
        print '{0} already downloaded in {1}'.format(outname, outDirr)
    
def get_pot_scat_xs_wrapper(Z, A, endfDirr, verbosity=False):
    '''Look in the ENDF file for the potential scattering XS'''
    # Metastable isomeric states use the groundstate A + 400
    effA = A % 400
    metastableStr = ''
    if A // 400 > 0:
        metastableStr = 'm'
    filename = 'endf_{0:02}{1:03}{2}_vii1'.format(Z, effA, metastableStr)
    filePath = os.path.join(endfDirr, filename)
    readgrouprDict = get_pot_scat_xs(filePath)
    potScatXS = readgrouprDict['pot']
    if verbosity > 2:
        print 'Potential cross section for ({0}, {1}) is {2} b'.format(Z, A, potScatXS)
    return potScatXS
