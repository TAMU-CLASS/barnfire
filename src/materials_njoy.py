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
from Readgroupr import get_short2mt_dict, get_endf_mt_list
import materials_util as util
import write_njoy as njoy
import makegroups as mg

def create_njoy_decks(inputDict, globalZASabList, njoyTDict, njoyBXSDict, globalTXSDict, verbosity=False):
    endfLib = 'vii1'
    legendreOrder = inputDict['legendreorder']
    groupOpt = inputDict['groupopt']
    ignoreTransferMatrices = inputDict['smallscattering']
    usePURR = not(inputDict['purroff'])
    #
    groupBdrs = get_group_boundaries(inputDict, verbosity)
    #
    endfName = 'endf/b-{0}'.format(endfLib)
    pendfScriptPaths, gendfScriptPaths, aceScriptPaths = [], [], []
    aceFilesDict = {}
    allowedInelasticThermalMTList = util.get_inelastic_thermal_mt_list()
    # txs means thermal xs name; et means element thermal name
    txs2mtDict = get_short2mt_dict(get_endf_mt_list())
    et2filenameDict = util.get_element_thermal_name_to_endf_filename_dict()
    et2matDict = util.get_element_thermal_name_to_mat_number_dict()
    et2inelmtDict = util.get_element_thermal_name_to_inelastic_mt_number_dict()
    njoyGroupOpt = get_njoy_group_opt(groupOpt)
    tapes = njoy.NJOYTape()
    if verbosity:
        print '------- Creating NJOY decks -------'
        print 'Folder ACE_file temperature(K)'
    #
    # Generate fast (non bound-thermal) ACE deck
    #
    for (Z,A,Sab) in sorted(globalZASabList):
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
        dat.sig0List = sorted(njoyBXSDict[(Z,A,Sab)], reverse=True)
        dat.thermList = sorted(njoyTDict[(Z,A,Sab)])
        dat.mat = nd.mats[(Z, Atrue, isExcitedState)]
        dat.isFissionable = util.is_fissionable((Z,A))
        dat.includeMF6 = not(ignoreTransferMatrices)
        dat.usePURR = usePURR
        dat.aceExt = util.get_ace_extension(Sab)
        #
        # Always process free so you can do simple reactions later. Ignore null thermal treatments.
        etSet = set([Sab, 'free'])
        etSet.discard('none')
        etList = sorted(etSet)
        thermalXSNameList = sorted(globalTXSDict[(Z,A,Sab)] | set(['free']) )
        # The following 3 variables are used by write_njoy in THERMR and must have the same order
        thermalFilenameList = [et2filenameDict[et] for et in etList]
        thermalMATList = [et2matDict[et] for et in etList]
        inelasticThermalMTList = [et2inelmtDict[et] for et in etList]
        dat.endfThermalFileList=thermalFilenameList
        dat.matThermalList = thermalMATList
        dat.inelasticMTList = inelasticThermalMTList
        # write_njoy uses the thermalMTList in GROUPR only
        thermalMTList = sorted([txs2mtDict[txs] for txs in thermalXSNameList])
        dat.thermalMTList = thermalMTList
        #
        dat.nuclideName = util.get_nuclide_dirr(nd.z2sym[Z], Atrue, Sab, metastableStr)
        dat.endfFile = 'endf_{0:02d}{1:03d}{2}_{3}'.format(Z, Atrue, metastableStr, endfLib)
        #
        if verbosity > 1:
            print Z, A, Sab, dat.sig0List, inelasticThermalMTList, thermalMATList
        if verbosity:
            for i,T in enumerate(dat.thermList):
                print '{0} {1:2d}{2:03d}.{3}{4}c {5:.1f}'.format(
                    dat.nuclideName, Z, Atrue, dat.aceExt, i, T)
        #
        pendfScriptPath, gendfScriptPath, aceScriptPath, aceFiles = njoy.create_njoy_script(dat, tapes)
        pendfScriptPaths.append(pendfScriptPath)
        gendfScriptPaths.append(gendfScriptPath)
        aceScriptPaths.append(aceScriptPath)
        aceFilesDict[dat.nuclideName] = aceFiles
    # Generate bound thermal ACE decks
    # (et for element thermal name; txs for thermal xs name)
    et2txsDict = util.get_element_thermal_name_to_thermal_xs_list_dict()
    et2mcnpDict = util.get_element_thermal_name_to_mcnp_thermal_name_dict()
    mcnp2zaDict = util.get_mcnp_thermal_name_to_main_za_dict()
    mcnp2zaidListDict = util.get_mcnp_thermal_name_to_zaid_list_dict()
    nonBoundETs = util.get_non_bound_names()
    # Determine unique bound element thermal names
    uniqueETs = {et for (Z,A,et) in globalZASabList if et not in nonBoundETs}
    #
    # Create an ACE deck for each unique bound thermal name
    #
    for et in sorted(uniqueETs):
        mcnpName = et2mcnpDict[et]
        # Which Z,A to reference for PENDF file
        Zref, Aref = mcnp2zaDict[mcnpName]
        #
        if Aref // 400 > 0:
            isExcitedState = True
            metastableStr = 'm'
        else:
            isExcitedState = False
            metastableStr = ''
        Atrue = Aref % 400
        dat = njoy.NJOYDat(A=Atrue, Z=Zref, endfName=endfName)
        dat.thermList = sorted(njoyTDict[(Zref,Aref,et)])
        dat.mat = nd.mats[(Zref, Atrue, isExcitedState)]
        dat.thermalNuclideName = mcnpName
        dat.thermalFileName = dat.thermalNuclideName.replace('/', '_')
        dat.aceExt = util.get_ace_extension(et)
        #
        dat.nuclideName = util.get_nuclide_dirr(nd.z2sym[Zref], Atrue, et, metastableStr)
        dat.endfFile = 'endf_{0:02d}{1:03d}{2}_{3}'.format(Zref, Atrue, metastableStr, endfLib)
        dat.ZAIDs = mcnp2zaidListDict[mcnpName]
        #
        thermalXSNameList = et2txsDict[et]
        thermalMTList = [txs2mtDict[txs] for txs in thermalXSNameList]
        dat.thermalMTList = thermalMTList
        #
        thermalFilenameList = [et2filenameDict[et]]
        thermalMATList = [et2matDict[et]]
        inelasticThermalMTList = [et2inelmtDict[et]]
        dat.endfThermalFileList=thermalFilenameList
        dat.matThermalList = thermalMATList
        dat.inelasticMTList = inelasticThermalMTList
        #
        if verbosity > 1:
            print et, et2mcnpDict[et], Zref, Aref, inelasticThermalMTList, thermalMATList
        if verbosity:
            for i,T in enumerate(dat.thermList):
                print '{} {}.{}{}t {:.1f}'.format(dat.thermalFileName, mcnpName, dat.aceExt, i, T)
        #
        thermalAceScriptPath, thermalAceFiles = njoy.create_thermal_ace_njoy_script(dat, tapes)
        aceScriptPaths.append(thermalAceScriptPath)
        aceFilesDict[dat.thermalFileName] = thermalAceFiles
    # Create a driver for both fast (non-thermal) and (bound-)thermal NJOY calls
    njoy.create_njoy_driver(pendfScriptPaths, gendfScriptPaths, aceScriptPaths)
    # Create a script to copy ACE files from ace directory to xdata directory
    njoy.create_ace_copier(aceFilesDict)

###############################################################################
def print_njoys(njoyTDict, njoyBXSDict, verbosity=False):
    if verbosity:
        print '------- NJOY -------'
        print 'njoyTDict', sorted(njoyTDict.items())
        print 'njoyBXSDict', sorted(njoyBXSDict.items())

def print_mcnp_material_inputs(materials, njoyTDict, verbosity=False):
    '''Print MCNP-style input for each material'''
    if verbosity:
        Sab2mcnp = util.get_element_thermal_name_to_mcnp_thermal_name_dict()
        for i, material in enumerate(materials):
            print '------- Material {} (ACE) -------'.format(i)
            name = material.shortName
            T = material.temperature
            T_MeV = 8.6173324E-11 * T
            # NB: Some of the abundanceDict's may be unnormalized because they neglect some of the isotopes
            a = material.abundanceDict
            e = material.elemAtomFracDict
            aD = material.atomDensity
            iT = material.temperatureIndex
            #
            thermalName = material.thermalOpt
            thermalStr = ''
            if thermalName not in util.get_non_bound_names():
                thermalStr = '-{}'.format(thermalName)
            #
            mcnpStrings = []
            mcnpThermalStrings = set()
            for (Z,A) in sorted(material.ZAList):
                sym = material.symDict[Z]
                Sab = material.SabDict[(Z,A)]
                #
                Tgrid = np.sort(list(njoyTDict[(Z,A,Sab)]))
                closestTindex = np.argmin(np.abs(T - Tgrid))
                distanceToClosestT = np.abs(T - Tgrid[closestTindex])
                nearestTs = util.get_nearest_points(T, Tgrid)
                #
                aceTindex = closestTindex
                aceSabIndex = util.get_ace_extension(Sab)
                aceExtension = 10*int(aceSabIndex) + int(aceTindex)
                #
                atomFrac = a[(Z,A)] * e[Z]
                #
                # Save MCNP material input for nuclide
                mcnpStrings.append('     {:02}{:03}.{:02}c    {:.8e}'.format(Z, A, aceExtension, atomFrac))
                if Sab in Sab2mcnp:
                    mcnpThermalName = Sab2mcnp[Sab]
                    mcnpThermalStrings.add('     {}.{:02}t'.format(mcnpThermalName, aceExtension))
                #
                # Warn if > 1K difference between desired temperature and NJOY evaluated temperature
                if distanceToClosestT > 1.0:
                    print '    Warning: Desired temperature for material {} is {} K,'.format(name, T)
                    print 'but nearest grid points are {} K and {} K for nuclide {}-{}{}.'.format(nearestTs[0], nearestTs[-1], sym, A, thermalStr)
                    print 'Input deck below uses ACE file corresponding to a temperature of {} K.'.format(Tgrid[closestTindex])
                #
                # Warn if material's temperature index does not match ACE index
                if iT != aceTindex:
                    print '    Warning: Material {} has a temperature index of {},'.format(name, iT)
                    print 'but nearest temperature index in ACE filesfor nuclide {}-{}{} is {}.'.format(
                        sym, A, thermalStr, aceTindex)
            #
            # Print MCNP material input
            print 'MCNP-style material input for {}:'.format(name)
            print '     atom_density (1/b-cm): {:.5f}'.format(aD)
            print '     tmp (MeV): {:.3e}'.format(T_MeV)
            print 'm{}'.format(i)
            for strr in mcnpStrings:
                print strr
            if mcnpThermalStrings:
                print 'mt{}'.format(i)
            for strr in sorted(mcnpThermalStrings):
                print strr
        
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
        #Do not crop group structure's energy bounds inputDict['energybounds']
        #(Do we want this??)
        parserStr += ' -e' 
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

def get_njoy_temperatures(globalZASabList, njoyTDict, globalTDict, globalTXSDict):
    '''Create a grid of temperatures per nuclide that will be used as NJOY input. Complicated because the non-free thermal XS are defined at discrete temperatures, and each thermal XS is defined on a different temperature grid.'''
    #To use for NJOY input:
    #for T in njoyTDict[(Z,A,Sab)]:
    #    for rxn in globalTXSDict[(Z,A,Sab)]:
    #        if T in allowedTDict[rxn]:
    #            include rxn
    maxNumTemperatures = 10
    allowedTDict = get_allowed_njoy_temperatures()
    for (Z,A,Sab) in globalZASabList:
        njoyTDict[(Z,A,Sab)] = set()
        if Sab in allowedTDict:
            # Use temperature grid where thermal XS's are defined
            sortedAllowedTList = np.array(sorted(allowedTDict[Sab]))
            for T in globalTDict[(Z,A,Sab)]:
                njoyTDict[(Z,A,Sab)].update(util.get_nearest_points(T, sortedAllowedTList))
        else:
            # No grid specified, free to use actual temperatures
            njoyTDict[(Z,A,Sab)].update(globalTDict[(Z,A,Sab)])
        #
        while len(njoyTDict[(Z,A,Sab)]) > maxNumTemperatures:
            # This may unintentially preferentially thin one thermal XS grid over another
            # in the case of multiple non-free thermal XS per nuclide
            itemToRemove = util.thin_list(njoyTDict[(Z,A,Sab)])
            njoyTDict[(Z,A,Sab)].remove(itemToRemove)

def get_njoy_background_xs(globalZASabList, njoyBXSDict, globalBXSDict, useCommonGrid=False):
    maxNumBXS = 10
    maxBXS = 1.e10
    minBXS = 1.e-1
    for (Z,A,Sab) in globalZASabList:
        if useCommonGrid:
            njoyBXSDict[(Z,A,Sab)] = set(get_common_sigma0_grid())
        else:
            njoyBXSDict[(Z,A,Sab)] = copy.copy(globalBXSDict[(Z,A,Sab)])
        #
        outOfRangeBXS = [bxs for bxs in njoyBXSDict[(Z,A,Sab)] if (bxs > maxBXS) or (bxs < minBXS)]
        for offendingBXS in outOfRangeBXS:
            njoyBXSDict[(Z,A,Sab)].remove(offendingBXS)
        #The first background XS must always be inf
        njoyBXSDict[(Z,A,Sab)].update([maxBXS])
        #
        while len(njoyBXSDict[(Z,A,Sab)]) > maxNumBXS:
            itemToRemove = util.thin_list(njoyBXSDict[(Z,A,Sab)], False)
            njoyBXSDict[(Z,A,Sab)].remove(itemToRemove)

def get_common_sigma0_grid():
    return [1.e10, 1.e4, 1.e3, 5.e2, 2.e2, 1.e2, 5.e1, 2.e1, 1.e1, 1.e0]

def get_njoy_group_opt(parserGroupOpt):
    '''Map the groupopt from the parser to groupOpt for NJOY.'''
    if parserGroupOpt.startswith('njoy'):
        return int(parserGroupOpt.split('-')[-1])
    else:
        return 1

def get_allowed_njoy_temperatures():
    '''S(alpha, beta) tables are only evaluated at certain temperatures and cannot be used at other temperatures. All T in K.'''
    # Tfuel works for both UO2 and ZrH
    Tfuel = set([296., 400., 500., 600., 700., 800., 1000., 1200.])
    # Tmod works for H2O
    Tmod = set([293.6, 350., 400., 450., 500., 550., 600., 650., 800.])
    # Tstructural works for both Fe and Al
    Tstructural = set([20., 80., 293.6, 400., 600., 800.])
    # Tmod works for graphite
    Tgraphite = set([296., 400., 500., 600., 700., 800., 1000., 1200., 1600., 2000.])
    #
    # Each element thermal name corresponds to a bound thermal ENDF file
    # Use the temperature grid from that file as the allowed NJOY evaluation temperatures for 
    # that element thermal name (aka Sab).
    allowedTDict = {
        'hh2o': Tmod,
        'ouo2': Tfuel,
        'uuo2': Tfuel,
        'hzrh': Tfuel,
        'zrzrh': Tfuel,
        'graph': Tgraphite,
        'al': Tstructural,
        'fe': Tstructural,
        }
    return allowedTDict




