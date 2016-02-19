'''
Andrew Till
Summer 2014

Functions for processing ENDF files
'''

#MINE
from directories import get_common_directories, copy_thermal_endf_xs, make_endf_directory

#STDLIB
import os
import sys
import shutil
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
    for (Z, A) in ZAList:
        get_endf_file(Z, A, endfDirr, verbosity)
        if useBXS:
            potXSDict[(Z,A)] = get_pot_scat_xs_wrapper(Z, A, endfDirr, verbosity)
        else:
            potXSDict[(Z,A)] = 0.0
    copy_thermal_endf_xs()
    if verbosity > 1:
        print 'potXSDict', potXSDict

def get_endf_file(Z, A, endfDirr, verbosity=False):
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
