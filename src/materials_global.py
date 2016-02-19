'''
Andrew Till
Summer 2014

Global parameters for unionizing materials.

'''
#TPL
import numpy as np
#MINE
from directories import get_common_directories
from materials_endf import get_endf_files
import materials_util as util

###############################################################################
def get_union_parameters(materials, globalTDict, globalBXSDict, globalTXSDict, useBXS=False, verbosity=False):
    globalZAList = get_union_nuclide_list(materials)
    globalZList = get_union_element_list(materials)
    #
    check_unique_short_names(materials, verbosity)
    #
    potXSDict = {}
    dirDict = get_common_directories()
    get_endf_files(globalZAList, dirDict, potXSDict, useBXS, verbosity)
    #
    get_global_temperature_dict(globalZAList, materials, globalTDict)
    #
    calc_background_xs(materials, potXSDict, useBXS)
    get_global_background_xs_dict(globalZAList, materials, globalBXSDict)
    #
    get_global_thermal_xs_dict(globalZList, materials, globalTXSDict)
    #
    return globalZAList, globalZList

def check_unique_short_names(materials, verbosity=False):
    shortNameSet = set()
    for material in materials:
        shortName = material.shortName
        if shortName in shortNameSet:
            raise ValueError('All materials should have unique short names')
        else:
            shortNameSet.update([shortName])
    if verbosity > 1:
        print 'Short names contained in problem'
        print shortNameSet

def get_union_nuclide_list(materials):
    globalZAList = set([])
    for material in materials:
        globalZAList |= material.ZAList
    return globalZAList

def get_union_element_list(materials):
    globalZList = set([])
    for material in materials:
        globalZList |= material.ZList
    return globalZList

def get_global_temperature_dict(globalZAList, materials, globalTDict):
    for (Z, A) in globalZAList:
        globalTDict[(Z, A)] = []
    for material in materials:
        for (Z, A) in material.ZAList:
            globalTDict[(Z, A)].append(material.temperature)
    for key in globalTDict:
        globalTDict[key] = set(globalTDict[key])

def get_global_thermal_xs_dict(globalZList, materials, globalTXSDict):
    for Z in globalZList:
        globalTXSDict[Z] = []
    for material in materials:
        for Z in material.ZList:
            globalTXSDict[Z] += material.thermalXSDict[Z]
    for key in globalTXSDict:
        globalTXSDict[key] = set(globalTXSDict[key])

###############################################################################
def get_global_background_xs_dict(globalZAList, materials, globalBXSDict):
    for (Z, A) in globalZAList:
        globalBXSDict[(Z, A)] = [np.inf]
    for material in materials:
        for (Z, A) in material.ZAList:
            globalBXSDict[(Z, A)].append(material.backgroundXSDict[(Z, A)])
    for key in globalBXSDict:
        globalBXSDict[key] = set(globalBXSDict[key])

def calc_background_xs(materials, potXSDict, useBXS=False):
    '''Use the potential scattering XS plus the chord length to determine background XS for each nuclide in each material'''
    for material in materials:
        material.backgroundXSDict = {}
        for (Zto, Ato) in material.ZAList:
            backgroundXS = 0.0
            if useBXS:
                for (Zfrom, Afrom) in material.ZAList:
                    if (Zto, Ato) != (Zfrom, Afrom):
                        backgroundXS += (
                        material.abundanceDict[(Zfrom, Afrom)] *
                        material.elemAtomFracDict[Zfrom] *
                        potXSDict[(Zfrom, Afrom)] )
                backgroundXS += material.chordLength / material.atomDensity
                backgroundXS /= (material.abundanceDict[(Zto, Ato)] * material.elemAtomFracDict[Zto])
                if util.has_bondarenko_iteration(Zto):
                    material.backgroundXSDict[(Zto, Ato)] = backgroundXS
                else:
                    material.backgroundXSDict[(Zto, Ato)] = np.inf
            else:
                material.backgroundXSDict[(Zto, Ato)] = np.inf
        material.check_background_xs_keys_consistency()

###############################################################################
def print_globals(globalZAList, globalZList, globalTDict, globalBXSDict, globalTXSDict, verbosity=False):
    if verbosity:
        print '------- Global Variables -------'
        print 'globalZAList', sorted(globalZAList)
        print 'globalZList', sorted(globalZList)
        print 'globalTDict', globalTDict
        print 'globalBXSDict', globalBXSDict
        print 'globalTXSDict', globalTXSDict



