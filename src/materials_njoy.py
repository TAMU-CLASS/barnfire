'''
Andrew Till
Summer 2014

NJOY utilities for materials
'''

#MINE
from directories import get_common_directories

#STDLIB
import os
import sys
import shutil
import copy
sys.path.insert(1, get_common_directories()['nuclideData'])

#TPL
import numpy as np
import nuclide_data as nd
#MINE
from Readgroupr import get_mt2short_dict, get_short2mt_dict, get_endf_mt_list
import materials_util as util
import write_njoy as njoy
import makegroups as mg

def create_njoy_decks(inputDict, globalZAList, njoyTDict, njoyBXSDict, globalTXSDict, verbosity=False):
    endfLib = 'vii1'
    legendreOrder = inputDict['legendreorder']
    groupOpt = inputDict['groupopt']
    ignoreTransferMatrices = inputDict['smallscattering']
    #
    groupBdrs = get_group_boundaries(inputDict, verbosity)
    #
    endfName = 'endf/b-{0}'.format(endfLib)
    pendfScriptPaths, gendfScriptPaths, aceScriptPaths = [], [], []
    aceFilesDict = {}
    allowedInelasticThermalMTList = util.get_inelastic_thermal_mt_list()
    short2mtDict = get_short2mt_dict(get_endf_mt_list())
    mt2shortDict = get_mt2short_dict(get_endf_mt_list())
    short2matDict = util.get_thermal_short2mat_dict()
    name2filenameDict = util.get_thermal_name2filename_dict()
    njoyGroupOpt = get_njoy_group_opt(groupOpt)
    tapes = njoy.NJOYTape()
    if verbosity:
        print '------- Creating NJOY decks -------'
    for (Z,A) in sorted(globalZAList):
        # Metastable isomeric states use A of the groundstate plus 400 (like MCNP)
        if A // 400 > 0:
            isExcitedState = True
            metastableStr = 'm'
        else:
            isExcitedState = False
            metastableStr = ''
        Atrue = A % 400
        dat = njoy.NJOYDat(A=Atrue, Z=Z, endfName=endfName, groupBdrs=groupBdrs, groupOpt=njoyGroupOpt,
            legendreOrder=legendreOrder)
        dat.sig0List = sorted(njoyBXSDict[(Z,A)], reverse=True)
        dat.thermList = sorted(njoyTDict[(Z,A)])
        dat.mat = nd.mats[(Z, Atrue, isExcitedState)]
        dat.isFissionable = util.is_fissionable((Z,A))
        dat.includeMF6 = not(ignoreTransferMatrices)
        # Always process free so you can do simple reactions later.
        thermalNameList = sorted(globalTXSDict[Z] | set(['free']) )
        #
        dat.nuclideName = '{0}-{1}{2}'.format(nd.z2sym[Z].lower(), Atrue, metastableStr)
        dat.endfFile = 'endf_{0:02d}{1:03d}{2}_{3}'.format(Z, Atrue, metastableStr, endfLib)
        thermalFilenameList = [name2filenameDict[name] for name in thermalNameList
            if name in name2filenameDict]
        #
        thermalMTList = [short2mtDict[name] for name in thermalNameList]
        inelasticThermalMTList = [mt for mt in thermalMTList if mt in allowedInelasticThermalMTList]
        thermalMATList = [short2matDict[name] for name in thermalNameList if name in short2matDict]
        dat.thermalMTList = thermalMTList
        dat.matThermalList = thermalMATList
        dat.inelasticMTList = inelasticThermalMTList
        dat.endfThermalFileList=thermalFilenameList
        #
        if verbosity:
            print Z, A, dat.sig0List, dat.thermList, inelasticThermalMTList, thermalMATList
        #
        pendfScriptPath, gendfScriptPath, aceScriptPath, aceFiles = njoy.create_njoy_script(dat, tapes)
        pendfScriptPaths.append(pendfScriptPath)
        gendfScriptPaths.append(gendfScriptPath)
        aceScriptPaths.append(aceScriptPath)
        aceFilesDict[dat.nuclideName] = aceFiles
    # TODO Insert loop here over thermal options (that adds to both aceScriptPaths and aceFilesDict)
    njoy.create_njoy_driver(pendfScriptPaths, gendfScriptPaths, aceScriptPaths)
    njoy.create_ace_copier(aceFilesDict)

###############################################################################
def print_njoys(njoyTDict, njoyBXSDict, verbosity=False):
    if verbosity:
        print '------- NJOY -------'
        print 'njoyTDict', njoyTDict
        print 'njoyBXSDict', njoyBXSDict

###############################################################################
def get_group_boundaries(inputDict, verbosity):
    groupOpt = inputDict['groupopt']
    desiredGroups = inputDict['groups']
    energyBounds = inputDict['energybounds']
    resolvedRange = inputDict['resolvedrange']
    if verbosity:
        print '---- Make Groups ----'
    #
    if groupOpt == 'equal':
        logEnergyBounds = np.log10(energyBounds)
        groupBdrs = np.logspace(logEnergyBounds[0], logEnergyBounds[1], desiredGroups + 1)
        return groupBdrs
    elif groupOpt == 'time':
        # This will make 1/v_{g} - 1/v_{g+1} = constant
        groupBdrs = np.zeros(desiredGroups + 1)
        c = (  (1 / (desiredGroups * np.sqrt(energyBounds[0]))) *
            (1 - np.sqrt(energyBounds[0]/energyBounds[1]))  )
        groupBdrs[desiredGroups] = energyBounds[1]
        for g in range(desiredGroups, 0, -1):
            groupBdrs[g-1] = 1 / np.square(c + 1/np.sqrt(groupBdrs[g]))
        return groupBdrs
    elif groupOpt.startswith('njoy'):
        return []
    #
    parserStr = ''
    energyBoundsStr = ' '.join(['{0:.10e}'.format(val) for val in energyBounds])
    resolvedRangeStr = ' '.join(['{0:.10e}'.format(val) for val in resolvedRange])
    parserStr += ' -R {0} -r {1}'.format(energyBoundsStr, resolvedRangeStr)
    #
    if groupOpt == 'default':
        parserStr += ' -e'
    else:
        filename = groupOpt
        if len(groupOpt.rsplit('.', 1)) == 1:
            filename += '.txt'
        if verbosity:
            print "Reading in group structure from '{0}'".format(filename)
        parserStr += ' -I {0}'.format(filename)
        parserStr += ' -e' #TODO: Don't always specify this, probably
    #
    parserStr += ' -g {0}'.format(desiredGroups)
    if verbosity:
        print "Saving group structure to 'energy_out_custom.txt'"
    if verbosity > 1:
        print './makegroups.py', parserStr
    mgInputDat = mg.define_input_parser().parse_args(parserStr.split())
    mg.finish_parsing_inputs(mgInputDat)
    groupBdrs = mg.resolve_group_structure(mgInputDat)
    numGroups = len(groupBdrs) - 1
    if verbosity:
        print numGroups, 'groups will be used'
    return sorted(groupBdrs)

def get_njoy_temperatures(globalZAList, njoyTDict, globalTDict, globalTXSDict):
    '''Create a grid of temperatures per nuclide that will be used as NJOY input. Complicated because the non-free thermal XS are defined at discrete temperatures, and each thermal XS is defined on a different temperature grid.'''
    #To use for NJOY input:
    #for T in njoyTDict[(Z,A)]:
    #    for rxn in globalTXSDict[Z]:
    #        if T in allowedTDict[rxn]:
    #            include rxn
    maxNumTemperatures = 10
    allowedTDict = {}
    get_allowed_njoy_temperatures(allowedTDict)
    for (Z, A) in globalZAList:
        allowedTList = set()
        njoyTDict[(Z, A)] = set()
        for thermalXS in globalTXSDict[Z]:
            if thermalXS != 'free':
                allowedTList |= allowedTDict[thermalXS]
        if allowedTList:
            # Use temperature grid where thermal XS's are defined
            sortedAllowedTList = np.array(sorted(allowedTList))
            for T in globalTDict[(Z, A)]:
                njoyTDict[(Z, A)] |= set(util.get_nearest_points(T, sortedAllowedTList))
        else:
            # No grid specified, free to use actual temperatures
            njoyTDict[(Z,A)] = copy.copy(globalTDict[(Z,A)])
        while len(njoyTDict[(Z,A)]) > maxNumTemperatures:
            # This may unintentially preferentially thin one thermal XS grid over another
            # in the case of multiple non-free thermal XS per nuclide
            itemToRemove = thin_list(njoyTDict[(Z,A)])
            njoyTDict[(Z,A)].remove(itemToRemove)

def get_njoy_background_xs(globalZAList, njoyBXSDict, globalBXSDict, useCommonGrid=False):
    maxNumBXS = 10
    maxBXS = 1.e10
    minBXS = 1.e-1
    for (Z,A) in globalZAList:
        if useCommonGrid:
            njoyBXSDict[(Z,A)] = set(get_common_sigma0_grid())
        else:
            njoyBXSDict[(Z,A)] = copy.copy(globalBXSDict[(Z,A)])
        #
        outOfRangeBXS = [bxs for bxs in njoyBXSDict[(Z,A)] if (bxs > maxBXS) or (bxs < minBXS)]
        for offendingBXS in outOfRangeBXS:
            njoyBXSDict[(Z,A)].remove(offendingBXS)
        #The first background XS must always be inf
        njoyBXSDict[(Z,A)].update([maxBXS])
        #
        while len(njoyBXSDict[(Z,A)]) > maxNumBXS:
            itemToRemove = thin_list(njoyBXSDict[(Z,A)], False)
            njoyBXSDict[(Z,A)].remove(itemToRemove)

def get_common_sigma0_grid():
    return [1.e10, 1.e4, 1.e3, 5.e2, 2.e2, 1.e2, 5.e1, 2.e1, 1.e1, 1.e0]

def get_njoy_group_opt(parserGroupOpt):
    '''Map the groupopt from the parser to groupOpt for NJOY.'''
    if parserGroupOpt.startswith('njoy'):
        return int(parserGroupOpt.split('-')[-1])
    else:
        return 1

def get_allowed_njoy_temperatures(allowedTDict):
    '''S(alpha, beta) tables are only evaluated at certain temperatures and cannot be used at other temperatures. All T in K.'''
    # Tfuel works for both UO2 and ZrH
    Tfuel = set([296., 400., 500., 600., 700., 800., 1000., 1200.])
    Tmod = set([293.6, 350., 400., 450., 500., 550., 600., 650., 800.])
    # Tstructural works for both Fe and Al
    Tstructural = set([20., 80., 293.6, 400., 600., 800.])
    Tgraphite = set([296., 400., 500., 600., 700., 800., 1000., 1200., 1600., 2000.])
    allowedTDict['ouo2inel'] = Tfuel
    allowedTDict['ouo2elas'] = Tfuel
    allowedTDict['uuo2inel'] = Tfuel
    allowedTDict['uuo2elas'] = Tfuel
    allowedTDict['hzrhinel'] = Tfuel
    allowedTDict['hzrhelas'] = Tfuel
    allowedTDict['zrzrhinel'] = Tfuel
    allowedTDict['zrzrhelas'] = Tfuel
    allowedTDict['alinel'] = Tstructural
    allowedTDict['alelas'] = Tstructural
    allowedTDict['feinel'] = Tstructural
    allowedTDict['feelas'] = Tstructural
    allowedTDict['hh2o'] = Tmod
    allowedTDict['graphinel'] = Tgraphite
    allowedTDict['graphelas'] = Tgraphite



