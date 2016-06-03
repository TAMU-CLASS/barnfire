'''
Andrew Till
Summer 2014

Global parameters for unionizing materials.

A nuclide is uniquely defined by its atomic number (Z), atomic mass (A), and thermal treatment (Sab)

'''
#TPL
import numpy as np
#MINE
from directories import get_common_directories
from materials_endf import get_endf_files, get_thermal_endf_files
import materials_util as util

###############################################################################
def get_union_parameters(materials, globalTDict, globalBXSDict, globalTXSDict, useBXS=False, verbosity=False):
    # Sab is for S(alpha,beta), the free/bound thermal treatment used by each nuclide
    globalZASabList = get_union_sab_nuclide_list(materials)
    globalZAList = get_union_nuclide_list(materials)
    globalZList = get_union_element_list(materials)
    #
    check_unique_short_names(materials, verbosity)
    #
    potXSDict = {}
    dirDict = get_common_directories()
    get_endf_files(globalZAList, dirDict, potXSDict, useBXS, verbosity)
    get_thermal_endf_files(globalZASabList, dirDict, verbosity)
    #
    get_global_temperature_dict(globalZASabList, materials, globalTDict)
    #
    calc_background_xs(materials, potXSDict, useBXS)
    get_global_background_xs_dict(globalZASabList, materials, globalBXSDict)
    #
    get_global_thermal_xs_dict(globalZASabList, materials, globalTXSDict)
    #
    #
    return globalZASabList, globalZAList, globalZList

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

def get_union_sab_nuclide_list(materials):
    globalZASabList = set()
    for material in materials:
        za2sab = material.SabDict
        for (Z, A) in material.ZAList:
            Sab = za2sab[(Z,A)]
            globalZASabList.add((Z,A,Sab))
    return globalZASabList
        
def get_union_nuclide_list(materials):
    globalZAList = set()
    for material in materials:
        globalZAList |= material.ZAList
    return globalZAList

def get_union_element_list(materials):
    globalZList = set()
    for material in materials:
        globalZList |= material.ZList
    return globalZList

def get_global_temperature_dict(globalZASabList, materials, globalTDict):
    for (Z,A,Sab) in globalZASabList:
        globalTDict[(Z,A,Sab)] = set()
    for material in materials:
        za2sab = material.SabDict
        for (Z, A) in material.ZAList:
            Sab = za2sab[(Z,A)]
            globalTDict[(Z,A,Sab)].add(material.temperature)

def get_global_thermal_xs_dict(globalZASabList, materials, globalTXSDict):
    for (Z,A,Sab) in globalZASabList:
        globalTXSDict[(Z,A,Sab)] = set()
    for material in materials:
        za2sab = material.SabDict
        for (Z, A) in material.ZAList:
            Sab = za2sab[(Z,A)]
            globalTXSDict[(Z,A,Sab)].update(material.thermalXSDict[(Z,A)])

###############################################################################
def get_global_background_xs_dict(globalZASabList, materials, globalBXSDict):
    for (Z,A,Sab) in globalZASabList:
        globalBXSDict[(Z,A,Sab)] = set([np.inf])
    for material in materials:
        za2sab = material.SabDict
        for (Z, A) in material.ZAList:
            Sab = za2sab[(Z,A)]
            globalBXSDict[(Z,A,Sab)].add(material.backgroundXSDict[(Z,A)])

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
def print_globals(globalZASabList, globalZAList, globalZList, globalTDict, globalBXSDict, globalTXSDict, verbosity=False):
    if verbosity:
        print '------- Global Variables -------'
        print 'globalZASabList', sorted(globalZASabList)
        print 'globalZAList', sorted(globalZAList)
        print 'globalZList', sorted(globalZList)
        print 'globalTDict', sorted(globalTDict.items())
        print 'globalBXSDict', sorted(globalBXSDict.items())
        print 'globalTXSDict', sorted(globalTXSDict.items())


