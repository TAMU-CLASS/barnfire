#! /usr/bin/env python

'''
Andrew Till
Fall 2014
PhD Research

Driver for creating the generalized energy mesh using clustering of observations based on similarity indicators. The indicators should be correlated with the flux: particles with energies that have similar indicators should have similar fluxes for all (r,Omega). These fluxes are allowed to depend on (r,Omega).
'''

#MINE
from directories import get_common_directories

#STDLIB
import os
import sys
sys.path.insert(1, get_common_directories()['nuclideData'])

#TPL
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nuclide_data as nd

#MINE
import plotutil as putil
from readxs import read_group_file
from makegroups import write_egrid
import Readgroupr as readgroupr
import materials_materials as mat
from materials_util import calc_chord_length, get_nuclide_dirr
from materials_global import get_union_parameters

mpl.rcParams.update({'font.size': 16, 'lines.linewidth': 2})

def define_defaults():
    '''Specify default parameters'''
    # Main parameters (plotOutput is always by default False)
    verbosity = True
    resolution = 9
    energySpacing = 1.5
    #workOpt = 'sigt'
    #workOpt = 'wgt'
    #workOpt = 'flux'
    workOpt = 'fluxe'

    # Specify range of interest (for output fluxes and plotting)
    #rrrRange = [6E-1, 1.06E3]
    #rrrRange = [3.0, 55.6]
    ##rrrRange = [55.6, 1.06E3]
    ####rrrRange = [3.0E0, 1.06E3]
    #rrrRange = [3.0, 2.5E4]
    rrrRange = [1E-5,2E7]

    # Specify which coarse group structure to use outside the RRR
    groupOpt = 'scale-44'

    # Specify constants for the fast Maxwellian shape
    #WattConstants = [0.88111, 3.4005] #thermal fission from U-238
    #WattConstants = [0.966, 2.842] # thermal fission from Pu-239
    WattConstants = [0.988, 2.2249] #thermal fission from U-235

    # Specify which set of materials to use
    #materialOpt = '3'
    #materialOpt = '4'
    ##materialOpt = '5'
    #materialOpt = 'c5g7'
    #materialOpt = 'graphite'
    #materialOpt = 'iron'
    materialOpt = 'manual'

    # If materialOpt is 'manual', list of materials to use
    materialsList = ['deb']

    # Specify energy mesh for output weights / fluxes
    #meshName = 'clust-{r}'
    meshName = 'res-{r}'

    return {'verbosity': verbosity, 'workopt': workOpt, 'resolution':resolution, 'energyspacing': energySpacing, 'rrr': rrrRange,'groupopt': groupOpt, 'materialopt': materialOpt, 'meshname': meshName, 'fastspectrumparam': WattConstants, 'listmaterials': materialsList}

def do_all(inputDict):
    '''Create Sigma_t for the specified materials on a proper grid for later clustering'''
    # Read in options
    verbosity = inputDict['verbosity']
    plotOutput = inputDict['plot']
    workOpt = inputDict['workopt']
    materialOpt = inputDict['materialopt']
    materialsList = inputDict['listmaterials']
    pwResFactor = inputDict['resolution']
    eSpacing = inputDict['energyspacing']
    rrrRange = inputDict['rrr']
    coarseGroupName = inputDict['groupopt']
    meshName = inputDict['meshname']
    useLowZScat = inputDict['Zlow']
    noScatOpt = inputDict['totalonly']
    temperatureDependence = inputDict['temperaturedependence']
    WattConstants = inputDict['fastspectrumparam']

    # Initialize options
    meshPath = get_mesh_path(meshName, pwResFactor)
    calcFlux = False
    useEscapeXS = False
    if workOpt in ['flux', 'fluxe', 'wgt']:
        calcFlux = True
    if workOpt in ['fluxe', 'wgt']:
        useEscapeXS = True
    normalizeSource = True
    if workOpt == 'wgt':
        normalizeSource = False
    linearTol, maxXSJump, maxFluxJump, maxdEJump = get_tolerances(pwResFactor)

    # Specify materials
    materials = []
    if materialOpt == '3':
        materials.append(mat.get_inner_hot_mox_material())
        materials.append(mat.get_middle_hot_mox_material())
        materials.append(mat.get_outer_hot_mox_material())
        # O in UO2 does not have a temperature at 550 K, so the PENDF lookup fails
        #if workOpt == 'wgt':
        #    materials.append(mat.get_hot_h2o_material())
        temperatureDependence = True
    if materialOpt in ['4','5']:
        materials.append(mat.get_cold_mox_material())
        materials.append(mat.get_cold_uo2_material())
        if workOpt == 'wgt':
            materials.append(mat.get_cold_h2o_material())
        if materialOpt == '5':
            chordLength = calc_chord_length(0.597)
            materials[1].update_chord_length(chordLength)
    if materialOpt == 'graphite':
        materials.append(mat.get_graphite_material())
    if materialOpt == 'iron':
        materials.append(mat.get_iron_material())
        if workOpt == 'wgt':
            materials.append(mat.get_pu_metal_material())
            materials.append(mat.get_thick_iron_material())
    elif materialOpt == 'c5g7':
        materials.append(mat.get_c5g7_uo2_material())
        materials.append(mat.get_c5g7_high_mox_material())
        if workOpt == 'wgt':
            materials.append(mat.get_c5g7_low_mox_material())
            materials.append(mat.get_c5g7_med_mox_material())
            materials.append(mat.get_c5g7_control_rod_material())
            materials.append(mat.get_c5g7_fission_chamber_material())
            materials.append(mat.get_c5g7_guide_tube_material())
            materials.append(mat.get_c5g7_moderator_material())
    elif materialOpt == 'kpin':
        materials.append(mat.get_kord_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_enriched_fuel_material())
            materials.append(mat.get_kord_zirconium_material())
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_rod_fuel_material())
    elif materialOpt == 'kenrichedpin':
        materials.append(mat.get_kord_fuel_material())
        materials.append(mat.get_kord_enriched_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_zirconium_material())
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_rod_fuel_material())
    elif materialOpt == 'kpin2d':
        materials.append(mat.get_kord_rod_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_fuel_material())
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_enriched_rod_fuel_material())
            materials.append(mat.get_kord_mox_rod_fuel_material())
    elif materialOpt == 'kenrichedpin2d':
        materials.append(mat.get_kord_rod_fuel_material())
        materials.append(mat.get_kord_enriched_rod_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_fuel_material())
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_mox_rod_fuel_material())
    elif materialOpt == 'kmoxpin2d':
        materials.append(mat.get_kord_rod_fuel_material())
        materials.append(mat.get_kord_mox_rod_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_fuel_material())
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_enriched_rod_fuel_material())
    elif materialOpt == 'kmoxenrichedpin2d':
        materials.append(mat.get_kord_rod_fuel_material())
        materials.append(mat.get_kord_enriched_rod_fuel_material())
        materials.append(mat.get_kord_mox_rod_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_fuel_material())
            materials.append(mat.get_kord_moderator_material())
    elif materialOpt == 'kcladpin':
        materials.append(mat.get_kord_fuel_material())
        materials.append(mat.get_kord_enriched_fuel_material())
        materials.append(mat.get_kord_zirconium_material())
        if workOpt == 'wgt':
            materials.append(mat.get_kord_moderator_material())
            materials.append(mat.get_kord_rod_fuel_material())
    elif materialOpt == 'trigafuel':
        materials.append(mat.get_triga_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_triga_clad_material())
            materials.append(mat.get_triga_zirconium_material())
            materials.append(mat.get_triga_graphite_material())
            materials.append(mat.get_triga_borated_graphite_material())
            materials.append(mat.get_triga_b4c_material())
            materials.append(mat.get_triga_moderator_material())
            materials.append(mat.get_triga_air_material())
            materials.append(mat.get_triga_grid_plate_material())
            materials.append(mat.get_triga_lead_material())
    elif materialOpt == 'ctrigafuel':
        materials.append(mat.get_depleted_triga_fuel_material())
        if workOpt == 'wgt':
            materials.append(mat.get_triga_clad_material())
            materials.append(mat.get_triga_zirconium_material())
            materials.append(mat.get_triga_graphite_material())
            materials.append(mat.get_triga_borated_graphite_material())
            materials.append(mat.get_triga_b4c_material())
            materials.append(mat.get_triga_moderator_material())
            materials.append(mat.get_triga_air_material())
            materials.append(mat.get_triga_grid_plate_material())
            materials.append(mat.get_triga_lead_material())
    elif materialOpt == 'trigamore':
        materials.append(mat.get_triga_fuel_material())
        materials.append(mat.get_triga_clad_material())
        materials.append(mat.get_triga_zirconium_material())
        if workOpt == 'wgt':
            materials.append(mat.get_triga_graphite_material())
            materials.append(mat.get_triga_borated_graphite_material())
            materials.append(mat.get_triga_b4c_material())
            materials.append(mat.get_triga_moderator_material())
            materials.append(mat.get_triga_air_material())
            materials.append(mat.get_triga_grid_plate_material())
            materials.append(mat.get_triga_lead_material())
    elif materialOpt == 'deb':
        materials.append(mat.get_bruss_enriched_rod_fuel_material())
    else:
        #materialOpt == 'manual'
        materialFunctionDict = mat.get_materials_name2function_dict()
        for materialName in materialsList:
            materials.append(materialFunctionDict[materialName]())

    # For each material, compute total cross sections or infinite-medium fluxes
    if verbosity:
        print 'Using resolution {0}'.format(pwResFactor)
    if calcFlux:
        # Infinite-medium flux
        energyGrid, fluxMat = compute_infinite_medium_flux(materials, linearTol, maxFluxJump, maxdEJump, eSpacing, WattConstants, useEscapeXS, normalizeSource, temperatureDependence, useLowZScat, noScatOpt, plotOutput, verbosity)
        indicatorMat = fluxMat
    else:
        # Total XS
        energyGrid, totalXSMat = compute_total_xs(materials, linearTol, maxXSJump, maxdEJump, temperatureDependence, plotOutput, verbosity)
        indicatorMat = totalXSMat

    # Get desired common group boundaries
    if (not useEscapeXS) and calcFlux:
        coarseGroupBdrs = read_default_group_structure(coarseGroupName)
        groupBdrs = form_reference_group_structure(
                energyGrid, coarseGroupBdrs, rrrRange, meshPath)
        if plotOutput:
            plot_flux(energyGrid, fluxMat, rrrRange)
    else:
        groupBdrs = read_group_file(meshPath)
    if verbosity:
        print 'Using energy mesh {0} with {1} groups'.format(meshPath, len(groupBdrs)-1)

    # Save flux/xs on common group boundaries
    if calcFlux:
        save_material_fluxes(materials, energyGrid, fluxMat, groupBdrs, pwResFactor, workOpt, plotOutput)
    else:
        save_material_xs(materials, energyGrid, totalXSMat, groupBdrs, pwResFactor, plotOutput)

###############################################################################
def compute_infinite_medium_flux(materials, linearTol, maxFluxJump, maxdEJump, eSpacing, WattConstants, useEscapeXS=False, normalizeSource=True, temperatureDependence=False, useLowZScat=False, noScatOpt=False, plotOutput=False, verbosity=False):
    '''Compute Q/Sigma_t for each material on desired energy grid'''
    '''eSpacing is for the final grid. maxdEJump is used before the flux calculation and thinning.'''

    #Parse materials
    numMaterials = len(materials)
    globalTDict, globalBXSDict, globalTXSDict = {}, {}, {}
    globalZASabList, globalZAList, _ = get_union_parameters(
        materials, globalTDict, globalBXSDict, globalTXSDict, False, verbosity)

    nuclideDirDict = {}
    globalZATList = []
    # (Z,A,Sab,T) are needed to get PENDF directory and temperature within the PENDF file.
    globalZASabTList = []
    for (Z,A,Sab) in globalZASabList:
        metastableStr = ''
        if A // 400 > 0:
            metastableStr = 'm'
        Atrue = A % 400
        nuclideDirr = get_nuclide_dirr(nd.z2sym[Z], Atrue, Sab, metastableStr)
        nuclideDirDict[(Z,A,Sab)] = nuclideDirr
        globalZATList.append((Z,A,-1))
        globalZASabTList.append((Z,A,Sab,-1))
    # Parse temperatures, if desired
    if temperatureDependence:
        globalZATList = []
        globalZASabTList = []
        # For now, require that each material temperature be on the pendf temperature grid (no interpolations)
        for material in materials:
            T = material.temperature
            for (Z,A) in material.ZAList:
                Sab = material.SabDict[(Z,A)]
                globalZATList.append((Z,A,T))
                globalZASabTList.append((Z,A,Sab,T))
        globalZATList = sorted(globalZATList)
        globalZASabTList = sorted(globalZASabTList)
        if verbosity:
            print 'globalZATList:', globalZATList
            print 'globalZASabTList:', globalZASabTList
    else:
        for material in materials:
            material.update_temperature(-1)

    # Read in pointwise cross sections from their PENDF files at the desired temperatures
    energyGridDict = {}
    totalXSDict = {}
    scatXSDict = {}
    filename = 'pendf_ascii.txt'
    dirDict = get_common_directories()
    rootDirr = dirDict['pendf']
    mts = [1, 2] # total xs, elastic scattering
    mtsStr = ' '.join([str(mtStr) for mtStr in mts])
    for (Z,A,Sab,T) in globalZASabTList:
        nuclideDirr = nuclideDirDict[(Z,A,Sab)]
        desiredT = T
        inDirr = os.path.join(rootDirr, nuclideDirr)
        #
        parser = readgroupr.define_input_parser()
        parseStr = '-i {i} -I {I} -w pen -T {T} -m {m}'.format(i=inDirr, I=filename, T=desiredT, m=mtsStr)
        if verbosity:
            print 'Looking for ({Z}, {A}, {Sab})... '.format(Z=Z,A=A,Sab=Sab),
            parseStr += ' -v'
        readerDict = vars(parser.parse_args(parseStr.split()))
        readgroupr.finish_parsing(readerDict)
        xsDict = readgroupr.execute_reader(readerDict)
        # No bound thermal XS are taken into account in indicators.py.
        # XS are stored as functions of (Z,A,T) and not Sab.
        energyGridDict[(Z,A,T)] = xsDict['energy']
        totalXSDict[(Z,A,T)] = xsDict[(3,1)]
        scatXSDict[(Z,A,T)] = xsDict[(3,2)]

    # Determine cross sections on union grid
    unionEnergyGrid = np.unique(np.concatenate(energyGridDict.values()))
    if plotOutput > 2:
        plot_dE(unionEnergyGrid, 'de_flux_0_union')

    materialIndexDict = {}
    # We  macro XS together to make it easier to thin / thicken
    unionXSMat = np.zeros((numMaterials, len(unionEnergyGrid)))
    compute_macro_xs(materialIndexDict, unionXSMat, unionEnergyGrid, energyGridDict, totalXSDict, materials, globalZATList)

    alphaDict = {}
    get_scattering_widths(alphaDict, globalZAList)

    # Thin grid to linear tolerance on cross sections
    fullEnergyGrid = unionEnergyGrid
    fullXSMat = unionXSMat
    unionEnergyGrid, unionXSMat = thin_grid(unionEnergyGrid, unionXSMat, linearTol, verbosity)
    if verbosity > 1:
        compute_thinning_error(fullEnergyGrid, fullXSMat, unionEnergyGrid, unionXSMat, True)
    if plotOutput > 2:
        plot_dE(unionEnergyGrid, 'de_flux_1_after_xs_thin')

    # Thicken the grid so the maximum relative jump in energy is bounded (for downscattering)
    unionEnergyGrid, unionXSMat = thicken_grid(unionEnergyGrid, unionXSMat, maxdEJump, 'x', verbosity)
    #plot_dE(unionEnergyGrid, 'de_flux_2_after_dE_thick')

    # Compute escape XS
    escapeXS = np.zeros(numMaterials)
    if useEscapeXS:
        if verbosity:
            print 'Using escape cross sections (1/cm):'
        for material in materials:
            iMat = materialIndexDict[material.shortName]
            escapeXS[iMat] += material.chordLength
            # It may be better to ignore O completely
            #oxygenXS = 3.8883 * material.elemAtomFracDict[8] * material.atomDensity
            #escapeXS[iMat] += oxygenXS
            if verbosity:
                print material.shortName, escapeXS[iMat]

    # Compute Q/Sigma_t
    unionFluxMat = np.zeros((numMaterials, len(unionEnergyGrid)))
    t0 = time.time()
    perform_slowing_down_calc(energyGridDict, scatXSDict, unionFluxMat, materialIndexDict, alphaDict, unionXSMat, unionEnergyGrid, materials, WattConstants, escapeXS, normalizeSource, useLowZScat, noScatOpt, plotOutput)
    if verbosity:
        print 'Performed slowing-down calculation in {0} sec'.format(time.time() - t0)
    if plotOutput > 1:
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionFluxMat, 'flux_2_after_dE_thick', 'Flux(E)/M(E)')

    # Thin grid to linear tolerance on flux
    fullEnergyGrid = unionEnergyGrid
    fullFluxMat = unionFluxMat
    unionEnergyGrid, unionFluxMat = thin_grid(unionEnergyGrid, unionFluxMat, linearTol, verbosity)
    if verbosity > 1:
        compute_thinning_error(fullEnergyGrid, fullFluxMat, unionEnergyGrid, unionFluxMat, True)
    if plotOutput > 2:
        plot_dE(unionEnergyGrid, 'de_flux_3_after_sigma_thin')
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionFluxMat, 'flux_3_after_sigma_thin', 'Flux(E)/M(E)')

    # Thicken the grid so the maximum relative jump in the flux is bounded
    unionEnergyGrid, unionFluxMat = thicken_grid(unionEnergyGrid, unionFluxMat, maxFluxJump, 'y', verbosity)
    if plotOutput > 1:
        plot_dE(unionEnergyGrid, 'de_flux_4_after_sigma_thick')
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionFluxMat, 'flux_4_after_sigma_thick', 'Flux(E)/M(E)')

    # Thicken the grid so the maximum relative jump in energy is bounded (for downscattering)
    unionEnergyGrid, unionFluxMat = thicken_grid(unionEnergyGrid, unionFluxMat, get_min_dE(eSpacing), 'x', verbosity)
    #plot_dE(unionEnergyGrid, 'de_flux_5_after_dE_thick')

    return unionEnergyGrid, unionFluxMat

###############################################################################
def compute_total_xs(materials, linearTol, maxXSJump, maxdEJump, temperatureDependence=False, plotOutput=False, verbosity=False):
    '''Get total cross sections for each material on desired energy grid'''

    #Parse materials
    numMaterials = len(materials)
    globalTDict, globalBXSDict, globalTXSDict = {}, {}, {}
    globalZASabList, globalZAList, _ = get_union_parameters(
        materials, globalTDict, globalBXSDict, globalTXSDict, False, verbosity)

    nuclideDirDict = {}
    globalZATList = []
    # (Z,A,Sab,T) are needed to get PENDF directory and temperature within the PENDF file.
    globalZASabTList = []
    for (Z,A,Sab) in globalZASabList:
        metastableStr = ''
        if A // 400 > 0:
            metastableStr = 'm'
        Atrue = A % 400
        nuclideDirr = get_nuclide_dirr(nd.z2sym[Z], Atrue, Sab, metastableStr)
        nuclideDirDict[(Z,A,Sab)] = nuclideDirr
        globalZATList.append((Z,A,-1))
        globalZASabTList.append((Z,A,Sab,-1))
    # Parse temperatures, if desired
    if temperatureDependence:
        globalZATList = []
        globalZASabTList = []
        # For now, require that each material temperature be on the pendf temperature grid (no interpolations)
        for material in materials:
            T = material.temperature
            for (Z,A) in material.ZAList:
                Sab = material.SabDict[(Z,A)]
                globalZATList.append((Z,A,T))
                globalZASabTList.append((Z,A,Sab,T))
        globalZATList = sorted(globalZATList)
        globalZASabTList = sorted(globalZASabTList)
        if verbosity:
            print 'globalZATList:', globalZATList
            print 'globalZASabTList:', globalZASabTList
    else:
        for material in materials:
            material.update_temperature(-1)

    # Read in pointwise cross sections from their PENDF files at the desired temperatures
    energyGridDict = {}
    totalXSDict = {}
    filename = 'pendf_ascii.txt'
    dirDict = get_common_directories()
    rootDirr = dirDict['pendf']
    mts = [1] # total xs
    mtsStr = ' '.join([str(mtStr) for mtStr in mts])
    for (Z,A,Sab,T) in globalZASabTList:
        nuclideDirr = nuclideDirDict[(Z,A,Sab)]
        desiredT = T
        inDirr = os.path.join(rootDirr, nuclideDirr)
        #
        parser = readgroupr.define_input_parser()
        parseStr = '-i {i} -I {I} -w pen -T {T} -m {m}'.format(i=inDirr, I=filename, T=desiredT, m=mtsStr)
        if verbosity:
            print 'Looking for ({Z}, {A}, {Sab})... '.format(Z=Z,A=A,Sab=Sab),
            parseStr += ' -v'
        readerDict = vars(parser.parse_args(parseStr.split()))
        readgroupr.finish_parsing(readerDict)
        xsDict = readgroupr.execute_reader(readerDict)
        # No bound thermal XS are taken into account in indicators.py.
        # XS are stored as functions of (Z,A,T) and not Sab.
        energyGridDict[(Z,A,T)] = xsDict['energy']
        totalXSDict[(Z,A,T)] = xsDict[(3,1)]

    # Concatenate the individual energy grids into a 1-D, ordered, unionized, numpy array
    unionEnergyGrid = np.unique(np.concatenate(energyGridDict.values()))

    # Compute the macroscopic XS on the union grid by linear interpolation of microscopic XS
    materialIndexDict = {}
    unionMacroTotalXSMat = np.zeros((numMaterials, len(unionEnergyGrid)))
    compute_macro_xs(materialIndexDict, unionMacroTotalXSMat, unionEnergyGrid, energyGridDict, totalXSDict, materials, globalZATList)
    if plotOutput > 1:
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionMacroTotalXSMat, 'xs_0_union')
        plot_dE(unionEnergyGrid, 'de_xs_0_union')

    # Thin the grid using a maximum error from linear approximation
    fullEnergyGrid = unionEnergyGrid
    fullMacroTotalXSMat = unionMacroTotalXSMat
    unionEnergyGrid, unionMacroTotalXSMat = thin_grid(unionEnergyGrid, unionMacroTotalXSMat, linearTol, verbosity)
    if verbosity > 1:
        compute_thinning_error(fullEnergyGrid, fullMacroTotalXSMat, unionEnergyGrid, unionMacroTotalXSMat, True)
    if plotOutput > 1:
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionMacroTotalXSMat, 'xs_1_after_xs_thin')
        plot_dE(unionEnergyGrid, 'de_xs_1_after_xs_thin')

    # Thicken the grid so the maximum relative jump in the total cross section is bounded
    unionEnergyGrid, unionMacroTotalXSMat = thicken_grid(unionEnergyGrid, unionMacroTotalXSMat, maxXSJump, 'y', verbosity)
    if verbosity > 1:
        compute_thinning_error(fullEnergyGrid, fullMacroTotalXSMat, unionEnergyGrid, unionMacroTotalXSMat, True)
    if plotOutput > 1:
        plot_sigma(materials, materialIndexDict, unionEnergyGrid, unionMacroTotalXSMat, 'xs_2_after_xs_thick')
        plot_dE(unionEnergyGrid, 'de_xs_2_after_xs_thick')

    return unionEnergyGrid, unionMacroTotalXSMat

###############################################################################
def perform_slowing_down_calc(energyGridDict, scatXSDict, unionFluxMat, materialIndexDict, alphaDict, unionXSMat, unionEnergyGrid, materials, WattConstants, macroEscapeXS=0., normalizeSource=True, useLowZScat=False, noScatOpt=False, plotOutput=False):
    '''Solve the slowing-down system to get an infinite-medium Q/Sigma_t
    These sources are meant to approximate NJOY's sources. See NJOY 2016 manual, page ~236 ff. (GROUPR)
    '''

    # Set the output plot directory
    figureDirr = get_common_directories()['figures/indicators']

    # Determine source
    #Ethermal = 0.1 #eV (NJOY 4 iwt example)
    #Efast = 820.3E3 #eV (NJOY 4 iwt example)
    #Efast = 1E5 #eV (NJOY 5 iwt approximately)
    Ethermal = 0.1 #eV
    thermalTinK = 300. #K
    Efast = 5E4 #eV
    EveryFast = 10E6 #eV (NJOY 5 iwt)
    #
    iThermal = np.argmin(np.abs(unionEnergyGrid - Ethermal))
    iFast = np.argmin(np.abs(unionEnergyGrid - Efast))
    iVeryFast = np.argmin(np.abs(unionEnergyGrid - EveryFast))
    iAll = len(unionEnergyGrid) - 1
    #
    baseSource = 1. / (unionEnergyGrid + 1E-12)
    #
    #fastSpectrum = get_watt_spectrum(unionEnergyGrid[iFast:], WattConstants)
    #baseSource[iFast:] = fastSpectrum * (baseSource[iFast] / fastSpectrum[0])
    fastSpectrum = get_watt_spectrum(unionEnergyGrid, WattConstants)
    fastSpectrum /= (fastSpectrum[iFast] / baseSource[iFast])
    baseSource += fastSpectrum
    #
    baseSource[iVeryFast:] = baseSource[iVeryFast]
    #
    thermalT = convert_K_to_eV(thermalTinK)
    thermalSpectrum = get_alternate_maxwellian_spectrum(unionEnergyGrid[:iThermal], thermalT)
    baseSource[:iThermal] = thermalSpectrum * (baseSource[iThermal-1] / thermalSpectrum[-1])
    if useLowZScat:
        # If low-Z scattering is included, don't double-count it by adding in 1/E
        saveSource = baseSource.copy()
        baseSource[iThermal:iVeryFast] -= 1. / (unionEnergyGrid[iThermal:iVeryFast] + 1E-12)
        baseSource[iVeryFast:] -= 1. / (unionEnergyGrid[iVeryFast] + 1E-12)
        baseSource[:iThermal] *= baseSource[iThermal] / baseSource[iThermal-1]
    if plotOutput:
        plt.figure(4)
        plt.clf()
        plt.loglog(unionEnergyGrid, baseSource)
        #plt.loglog(unionEnergyGrid, baseSource, 'o')
        if unionEnergyGrid[-1] > 2E7:
            plt.xlim(right=2E7)
        filename = 'p_src.pdf'
        plotPath = os.path.join(figureDirr, filename)
        plt.savefig(plotPath)
        #
    # Compute slowing-dow flux for each material
    unionFluxMat[:, :] = 0.
    for material in materials:
        # Determine inhomogeneous source (which includes incident partial currents)
        T = material.temperature
        source = baseSource.copy()
        iMat = materialIndexDict[material.shortName]
        matlEscapeXS = macroEscapeXS[iMat]
        if matlEscapeXS:
            source *= matlEscapeXS
        unionFluxMat[iMat, iFast:] = source[iFast:] / unionXSMat[iMat, iFast:]
        #
        if noScatOpt:
            # Do not perform down-scattering calculation
            unionFluxMat[iMat, :] = source / (unionXSMat[iMat,:] + matlEscapeXS)
            continue
        # Determine which (Z,A) to include in the downscattering calculation
        if useLowZScat:
            ZAList = list(material.ZAList)
        else:
            # For systems where the important behavior is high-Z, can get away with this
            ZAList = [(Z,A) for (Z,A) in list(material.ZAList) if Z >= 90]
        #ZAList = [(92,235), (92,238), (94,239)]
        # Should probably do this for Zr:
        #ZAList = [(Z,A) for (Z,A) in list(material.ZAList) if Z >= 26]
        #
        # Solve slowing-down equation for a given source on the union energy grid
        for iEnergy in range(iVeryFast-1,-1,-1):
            energy = unionEnergyGrid[iEnergy]
            lhsCoeff = (matlEscapeXS + unionXSMat[iMat, iEnergy])
            scatSrc = 0.
            for (Z,A) in ZAList:
                alpha = alphaDict[(Z,A)]
                #
                highestEnergy = energy / alpha
                iStrt = iEnergy
                iEnd = iStrt
                while iEnd < iAll and unionEnergyGrid[iEnd] < highestEnergy:
                    iEnd += 1
                # iEnd is the highest energy index: ueg[iEnd-1] < highestEnergy <= ueg[iEnd]
                EL = unionEnergyGrid[iStrt:iEnd]
                ER = unionEnergyGrid[iStrt+1:iEnd+1]
                ERinside = ER.copy()
                ERinside[-1] = highestEnergy
                wgtL = EL / (ERinside + EL)
                wgtR = 1. - wgtL
                dE = ERinside - EL
                #
                scatCoeff = np.zeros(iEnd - iStrt + 1)
                scatCoeff[:-1] += wgtL * dE / EL
                scatCoeff[1:]  += wgtR * dE / ER
                #
                atomDensity = material.atomDensity * material.elemAtomFracDict[Z] * material.abundanceDict[(Z,A)]
                scatXS = np.interp(unionEnergyGrid[iStrt:iEnd+1], energyGridDict[(Z,A,T)], scatXSDict[(Z,A,T)])
                scatCoeff *= scatXS * atomDensity / (1 - alpha)
                #
                scatSrc += np.sum(scatCoeff[1:] * unionFluxMat[iMat, iStrt+1:iEnd+1])
                lhsCoeff -= scatCoeff[0]
            unionFluxMat[iMat, iEnergy] = (source[iEnergy] + scatSrc) / lhsCoeff
        unionFluxMat[iMat, iVeryFast:] = unionFluxMat[iMat, iVeryFast-1]
        if useLowZScat and normalizeSource:
            # Total hack to get the long-term energy shape approximately correct
            unionFluxMat[iMat, :] /= saveSource
            # Original
            #unionFluxMat[iMat, iFast:] *= np.sqrt(unionEnergyGrid[iFast:])
            #unionFluxMat[iMat, :iFast] *= np.sqrt(unionEnergyGrid[iFast])
            unionFluxMat[iMat, :] *= np.sqrt(unionEnergyGrid[:])
            unionFluxMat[iMat, iVeryFast:] = unionFluxMat[iMat, iVeryFast-1]
        elif normalizeSource:
            # Not sure this is right (divide by source or baseSource?)
            unionFluxMat[iMat, :] /= source

def get_scattering_widths(alphaDict, globalZAList):
    '''Determine the alpha's for a given nuclide'''
    for (Z,A) in globalZAList:
        if A:
            wgt = nd.weight(Z, A)
        else:
            wgt = nd.weight(Z)
        alpha = ((wgt-1)/(wgt+1))**2
        alphaDict[(Z,A)] = alpha

###############################################################################
def compute_macro_xs(materialIndexDict, unionXSMat, finalEnergyGrid, energyGridDict, xsDict, materials, globalZATList, offset=0):

    # Linearly interpolate the microscopic cross sections onto the final energy grid
    unionXSDict = {}
    for (Z,A,T) in globalZATList:
       unionXSDict[(Z,A,T)] = np.interp(finalEnergyGrid, energyGridDict[(Z,A,T)], xsDict[(Z,A,T)])

    # Form the macroscopic cross sections on the final energy grid
    numMaterials = len(materials)
    unionXSMat[offset:(numMaterials + offset), :] = 0
    for iMat, material in enumerate(materials):
        shortName = material.shortName
        materialIndexDict[shortName] = iMat + offset
        ZAList = list(material.ZAList)
        T = material.temperature
        for (Z,A) in ZAList:
            atomDensity = material.atomDensity * material.elemAtomFracDict[Z] * material.abundanceDict[(Z,A)]
            unionXSMat[iMat+offset, :] += atomDensity * unionXSDict[(Z,A,T)]

###############################################################################
def get_min_dE(resFactor, Amax=250.):
    '''Return the minimum dE/E = 1 - alpha. This is a conservative measure of du. Amax is the heaviest nuclide likely to be encountered. A large resFactor is a low resolution.'''
    alphaMax = ((Amax - 1) / (Amax + 1)) ** 2 # Maximum downscattering
    return resFactor * (1 - alphaMax)

def compute_observations(signal, strt, end):
    '''
    Takes a positive signal that varies over several magnitudes and returns an observation.
    Uses only the part of the signal that is between strt and end.
    Assumes signal is indexed by [dimension, gridPoint].
    Returns an observation indexed by [gridPoint, dimension]
    '''
    return np.log10(np.transpose(signal[:, strt:end]))

def find_range(grid, strtVal, endVal):
    '''Returns strtIndex, endIndex such that
    grid[strtIndex] >= strtVal and grid[endIndex-1] < endVal.
    Assumes grid is sorted ascendingly.'''

    strtIndex = np.argmin(np.abs(grid - strtVal))
    if grid[strtIndex] < strtVal and strtIndex != len(grid)-1:
        strtIndex += 1
    endIndex = np.argmin(np.abs(grid - endVal))
    if grid[endIndex] >= endVal and endIndex != 0:
        endIndex -= 1
    endIndex += 1
    return strtIndex, endIndex

###############################################################################
def get_watt_spectrum(E, constants):
    '''The constants are assumed in MeV. E is assumed in eV.'''
    a = constants[0] * 1E6
    b = constants[1] / 1E6
    c = np.exp(-a * b / 4) / np.sqrt(np.pi * a**3 * b / 4)
    alpha = E/a
    beta = np.sqrt(b*E)
    spectrum = np.zeros(len(E))
    # Split into two pieces to deal with overflows of sinh
    mask = beta < 300
    spectrum[mask] = c * np.exp(-alpha[mask]) * np.sinh(beta[mask])
    mask = np.logical_not(mask)
    spectrum[mask] = (c/2) * (np.exp(beta[mask]-alpha[mask]) - np.exp(-beta[mask]-alpha[mask]))
    return spectrum

def get_maxwellian_spectrum(E, T):
    '''E and T must be in units of eV.
    q(E) ~ E/T * exp(-E/T)
    Use for q(E) if (Sigma_0 + Sigma_t(E)) * f(E) = Sigma_0 * q(E) and Sigma_0 >> Sigma_t(E)
    '''
    c = np.sqrt(8 / (np.pi * T))
    return c * E/T * np.exp(-E/T)

def get_alternate_maxwellian_spectrum(E, T):
    '''E and T must be in units of eV.
    q(E) ~ sqrt(E/T) * exp(-E/T)
    Use for q(E) if Sigma_t(E) * f(E) = q(E) and Sigma_t(E) ~ 1/sqrt(E)
    '''
    c = 2 / (np.sqrt(np.pi) * T)
    return c * np.sqrt(E/T) * np.exp(-E/T)

def convert_K_to_eV(T):
    return 8.6173324E-5 * T

###############################################################################
def get_tolerances(pwResFactor):
    '''The higher pwResFactor, the better the indicators will be, but the larger they will be as well'''
    if pwResFactor == 0:
        linearTol = 3E-1
        maxXSJump = 7E-1
        maxFluxJump = 7E-1
        maxdEJump = get_min_dE(7)
    elif pwResFactor == 1:
        linearTol = 2E-1
        maxXSJump = 7E-1
        maxFluxJump = 7E-1
        maxdEJump = get_min_dE(7)
    elif pwResFactor == 2:
        linearTol = 1E-1
        maxXSJump = 6.05E-1
        maxFluxJump = 6.05E-1
        maxdEJump = get_min_dE(7)
    elif pwResFactor == 3:
        linearTol = 7E-2
        maxXSJump = 6E-1
        maxFluxJump = 6E-1
        maxdEJump = get_min_dE(7)
    elif pwResFactor == 4:
        linearTol = 7E-2
        maxXSJump = 5E-1
        maxFluxJump = 5E-1
        maxdEJump = get_min_dE(5)
    elif pwResFactor == 5:
        linearTol = 5E-2
        maxXSJump = 3E-1
        maxFluxJump = 3E-1
        maxdEJump = get_min_dE(5)
    elif pwResFactor == 6:
        linearTol = 2E-2
        maxXSJump = 2E-1
        maxFluxJump = 2E-1
        maxdEJump = get_min_dE(5)
    elif pwResFactor == 7:
        linearTol = 1E-2
        maxXSJump = 1.5E-1
        maxFluxJump = 1.5E-1
        maxdEJump = get_min_dE(5)
    elif pwResFactor == 8:
        linearTol = 8E-3
        maxXSJump = 1E-1
        maxFluxJump = 1E-1
        maxdEJump = get_min_dE(5)
    elif pwResFactor == 9:
        linearTol = 5E-3
        maxXSJump = 5E-2
        maxFluxJump = 5E-2
        maxdEJump = get_min_dE(0.1)
    elif pwResFactor == 10:
        # Use only for graphite
        linearTol = 1E-3
        maxXSJump = 3E-3
        maxFluxJump = 3E-3
        maxdEJump = get_min_dE(0.1, 12)
    return linearTol, maxXSJump, maxFluxJump, maxdEJump

###############################################################################
def thin_grid(xj, yij, linearTol, verbosity=False):
    '''
    xj is the initial grid. yij are the data points that live on the grid.
    There may be more than one set of these data points (multiple dimensions).
    f_i(x) is linearly interpolated using the xj and yij.
    This function removes x points if the error made by linear interpolation
    on the new grid without those points is less than the input tolerance.
    It does this iteratively, removing only half of the points if several
    sequential points are to be removed.
    '''

    if verbosity > 1:
        print 'Thinning with tolerance {0}'.format(linearTol)
    numDim = len(yij)
    numPoints = len(xj)
    originalNumPoints = numPoints
    pointsRemoved = True
    iteration = 0
    while pointsRemoved:
        iteration += 1
        error = np.zeros(numPoints)
        # Do not remove endpoints
        error[0] = np.inf
        error[-1] = np.inf
        for iDim in range(numDim):
            linearityError = (np.abs(
                yij[iDim, 1:-1] - ( yij[iDim, :-2] +
                ( yij[iDim, 2:]- yij[iDim, :-2]) * (xj[1:-1] - xj[:-2]) / (xj[2:] - xj[:-2]))) /
                yij[iDim, 1:-1])
            error[1:-1] = np.maximum(error[1:-1], linearityError)
        toRemove = error <= linearTol
        # If two sequential points are to be removed, only remove the one with the odd index
        toRemove[:-1] = np.logical_and(toRemove[:-1], np.logical_not(
            np.logical_and(toRemove[1:], np.arange(numPoints-1)%2==1)))
        #toRemove[1:] = np.logical_xor(np.logical_and(np.logical_and(toRemove[1:], toRemove[:-1]), np.arange(numPoints-1)%2==1 ),toRemove[1:])
        numToRemove = np.sum(toRemove)
        if not numToRemove:
            pointsRemoved = False
            break
        numPoints -= numToRemove
        toKeep = np.logical_not(toRemove)
        xj = xj[toKeep]
        yij = yij[:, toKeep]
        if verbosity > 1:
            print 'Removing {0} of {1} points in iteration {2}'.format(
                numToRemove, numPoints+numToRemove, iteration)
    if verbosity:
        print 'Reduced grid from {0} to {1} points in {2} iterations'.format(
            originalNumPoints, numPoints, iteration)
    return xj, yij

def compute_thinning_error(xOriginal, yOriginal, xThinned, yThinned, verbosity=False):
    '''
    xOriginal and yOriginal are inputs to thin_grid
    xThinned and yThinned are outputs from thin_grid
    x is xj; y is yij, where i is the dimension index
    Points are not removed in thin_grid if the error in removing those points
    is larger than the tolerance. The error in removing those points is
    computed with respect to the current grid, which means that the maximum
    error reported here may not be less than the tolerance used in thin_grid,
    especially if multiple iterations of thin_grid are used.
    '''
    numDim, numPoints = yOriginal.shape
    yErr = np.zeros(numPoints)
    for iDim in range(numDim):
        yNewInterp = np.interp(xOriginal, xThinned, yThinned[iDim,:])
        yErr = np.maximum(yErr, np.abs((yOriginal[iDim,:] - yNewInterp) / yOriginal[iDim,:]))
    maxErr = max(yErr)
    if verbosity:
        print 'Maximum fractional error from thinning was {0}'.format(maxErr)
    return maxErr

def thicken_grid(xj, yij, maxRelJump, opt='x', verbosity=False):
    '''
    xj is the initial grid. yij are the data points that live on the grid.
    There may be more than one set of these data points (multiple dimensions).
    f_i(x) is linearly interpolated using the xj and yij.
    This function ensures a minimum grid spacing in x or y, depending on opt:
    if 'opt' == 'x':
        adds a mid-point between x[j] and x[j-1] if abs(x[j] - x[j-1]) / x[j-1/2] > tol
    elif 'opt' == 'y':
        adds a mid-point between x[j] and x[j-1] if abs(y[j] - y[j-1]) / y[j-1/2] > tol
    Assumes x and y are positive and x is sorted.
    Midpoints are defined geometrically.
    It does this iteratively, until there are no points left to add.
    Thickening should be applied before thinning, or else fidelity will be lost.
    Do not take the log of xj or yij before calling this thicken function.
    '''

    if verbosity > 1:
        print 'Thickening with jump tolerance {0} on {1}'.format(maxRelJump, opt)
    numDim = len(yij)
    numPoints = len(xj)
    originalNumPoints = numPoints
    pointsAdded = True
    iteration = 0
    initialX = xj[0]
    if xj[0] == 0:
        xj[0] = xj[1]
    while pointsAdded:
        iteration += 1
        if opt == 'x':
            jump = np.diff(xj) / xj[:-1]
        elif opt == 'y':
            jump = np.max(np.abs(np.diff(yij)) / np.minimum(yij[:, 1:], yij[:, :-1]), axis=0)
        toAdd = np.where(jump > maxRelJump)[0]
        numToAdd = len(toAdd)
        if not numToAdd:
            pointsAdded = False
            break
        addLocations = np.sqrt(xj[toAdd] * xj[toAdd+1])
        # No aliasing is done
        xjOld = xj
        xj = np.union1d(xj, addLocations)
        yijOld = yij
        numPoints = len(xj)
        yij = np.zeros((numDim, numPoints))
        for iDim in range(numDim):
            yij[iDim,:] = np.interp(xj, xjOld, yijOld[iDim, :])
        if verbosity > 1:
            print 'Adding {0} points to existing {1} points in iteration {2}'.format(
                numToAdd, numPoints-numToAdd, iteration)
    if verbosity:
        print 'Enlarged grid from {0} to {1} points in {2} iterations'.format(
            originalNumPoints, numPoints, iteration)
    xj[0] = initialX
    return xj, yij

def compute_x_jump(xj, verbosity=False):
    '''Compute the maximum relative change in the ascending grid'''
    jump = np.max(np.diff(xj) / xj[:-1])
    if verbosity:
        print 'Maximum jump in the grid is {0}'.format(jump)
    return jump

def compute_y_jump(yij, verbosity=False):
    '''Compute the maximum relative change in the data points'''
    jump = np.max(np.abs(np.diff(yij)) / np.minimum(yij[:, 1:], yij[:, :-1]))
    if verbosity:
        print 'Maximum jump in the data points is {0}'.format(jump)
    return jump


###############################################################################
def insert_group_structure(insertedGrid, existingGrid, minSpacing=0.003):
    '''Inserts a grid into an existing grid. Assumes both are ascendingly sorted.
    Remove nearby points in existingGrid if they are too close to the boundaries
    of insertedGrid.'''
    gStart, gEnd = find_range(existingGrid, insertedGrid[0], insertedGrid[-1])
    if gStart > 1 and abs(1 - insertedGrid[0] / existingGrid[gStart - 1]) <= minSpacing:
        gStart -= 1
    if gEnd < len(existingGrid) and abs(1 - insertedGrid[-1] / existingGrid[gEnd]) <= minSpacing:
        gEnd += 1
    unionSize = len(insertedGrid) + (len(existingGrid) - gEnd) + gStart
    unionGrid = np.zeros(unionSize)
    unionGrid[:gStart] = existingGrid[:gStart]
    unionGrid[gStart:gStart+len(insertedGrid)] = insertedGrid
    unionGrid[gStart+len(insertedGrid):] = existingGrid[gEnd:]
    unionGrid = np.union1d(unionGrid, existingGrid)[::-1] # Hack!!!
    return unionGrid

def form_dual(pwGrid, pwData, doIntegral=True):
    '''Form the dual grid and make averages on each dual cell equal to the pw point inside. This differs in behavior from calc_avg_from_pw, which does proper lin-lin trapezoid integration and tends to smear the averages on the dual grid out.'''
    dualGrid = form_dual_of_pw_grid(pwGrid)
    dE = np.diff(dualGrid)
    if doIntegral:
        dualData = pwData * dE[np.newaxis, :]
    else:
        dualData = pwData.copy()
    return dualGrid, dualData

def form_dual_of_pw_grid(pwGrid):
    '''Creates the dual of a pointwise grid. Defines centerpoint between two points as the sqrt of the two points (assumes pointwise grid is logarithmically spaced). Does not extend outer boundaries of pw grid. Returns the boundaries of the dual grid.'''
    dualGrid = np.zeros(len(pwGrid) + 1)
    dualGrid[1:-1] = np.sqrt(pwGrid[1:] * pwGrid[:-1])
    dualGrid[0] = pwGrid[0]
    dualGrid[-1] = pwGrid[-1]
    if pwGrid[0] == 0:
        dualGrid[1] = 0.5 * pwGrid[1]
    return dualGrid

def condense_like_points(dualGrid, labels):
    '''Given the boundaries of a 1-D grid and a labeling for each cell, returns a smaller grid where contiguous elements of the original grid have been combined. Assumes len(labels) + 1 == len(dualGrid).'''
    numOldBoundaries = len(dualGrid)
    keepMask = np.ones(numOldBoundaries, dtype=np.bool)
    # Keep unless the cell before you has the same label. Always keep the first cell. Always keep the outer boundaries.
    keepMask[1:-1] = labels[1:] != labels[:-1]
    labels = labels[keepMask[:-1]]
    dualGrid = dualGrid[keepMask]
    return dualGrid, labels

def calc_avg_from_pw(pwGrid, pwMat, avgGrid, doIntegral=True):
    '''Given data on a pointwise grid, compute the data averaged to an average grid.
    pwMat is indexed by [dimension/point, energy]. pwGrid is indexed by [energy].
    avgGrid should contain the boundaries for the average grid.
    Both grids must be sorted ascendingly.
    Does not extrapolate well, and fails if the entire cell is outside the pw grid.
    If doIntegral is true, the values on the avgGrid will be integrated over each cell.
    If doIntegral is false, the values on the avgGrid will be averaged over each cell.'''

    numDim, numPWEnergies = pwMat.shape
    numAvgGroups = avgGrid.shape[0] - 1
    avgMat = np.zeros((numDim,numAvgGroups))

    dEPW = np.diff(pwGrid)
    dEAvg = np.diff(avgGrid)

    for g in range(numAvgGroups):
        # See find_range for details on how boundaries are handled
        gPWStart, gPWEnd = find_range(pwGrid, avgGrid[g], avgGrid[g+1])
        if avgGrid[g+1] < pwGrid[0] or avgGrid[g] > pwGrid[-1]:
            print 'Error: Do not use for extrapolation.'
            exit(1)
        if gPWStart < gPWEnd:
            # Do a multi-point trapezoid rule to determine the integral inside avg cell g.
            # Start with the contributions from trapezoids internal to the cell g.
            numInternalPoints = gPWEnd - gPWStart
            wgts = np.zeros(numInternalPoints)
            wgts[:-1] += 0.5 * dEPW[gPWStart:gPWEnd-1]
            wgts[1:] += 0.5 * dEPW[gPWStart:gPWEnd-1]
            avgMat[:, g] = np.sum(pwMat[:, gPWStart:gPWEnd] * wgts[np.newaxis, :], axis=1)
            # Handle left boundary
            dEL = (pwGrid[gPWStart] - avgGrid[g])
            if gPWStart > 0:
                wL = 0.5 * dEL / dEPW[gPWStart-1]
                avgMat[:,g] += pwMat[:, gPWStart-1] * wL * dEL
                avgMat[:,g] += pwMat[:, gPWStart] * (1 - wL) * dEL
            else:
                avgMat[:,g] += pwMat[:,0] * dEL
            # Handle right boundary
            dER = (avgGrid[g+1] - pwGrid[gPWEnd-1])
            if gPWEnd < len(pwGrid):
                wR = 0.5 *  dER / dEPW[gPWEnd-1]
                avgMat[:,g] += pwMat[:, gPWEnd] * wR * dER
                avgMat[:,g] += pwMat[:, gPWEnd-1] * (1 - wR) * dER
            else:
                avgMat[:,g] += pwMat[:, -1] * dER
        else:
            #Avg cell g does not contain any pw points. Do one trapezoid for interpolation.
            gL = gPWEnd-1
            gR = gPWEnd
            valL = pwMat[:,gL] + (pwMat[:,gR] - pwMat[:,gL]) * (avgGrid[g] - pwGrid[gL]) / (pwGrid[gR] - pwGrid[gL])
            valR = pwMat[:,gL] + (pwMat[:,gR] - pwMat[:,gL]) * (avgGrid[g+1] - pwGrid[gL]) / (pwGrid[gR] - pwGrid[gL])
            avgMat[:,g] = 0.5 * (valL + valR) * dEAvg[g]
        if not doIntegral:
            avgMat[:,g] /= dEAvg[g]
    return avgMat

###############################################################################
def read_default_group_structure(groupName):
    '''Read in the background MG group structure'''
    groupFilename = '{0}.txt'.format(groupName)
    datDirr = get_common_directories()['dat/energy_groups']
    groupFilePath = os.path.join(datDirr, groupFilename)
    groupBdrs = read_group_file(groupFilePath)
    return groupBdrs

def get_mesh_path(meshName, pwResFactor):
    '''Determine the mesh path from a base filename'''
    meshPath = meshName.format(r=pwResFactor)
    if not meshName.count('.'):
        meshPath += '.txt'
    if not os.path.split(meshPath)[0]:
        datDirr = get_common_directories()['dat/energy_groups']
        meshPath = os.path.join(datDirr, meshPath)
    return meshPath

###############################################################################
def plot_dE(energyGrid, baseName='energy'):
    '''Plot the lethargy width between points in the energy grid and the cumulative points vs energy'''

    # Set the output plot directory
    figureDirr = get_common_directories()['figures/indicators']

    #Plot lethargy spacing (distance between points in the energy grid) as a function of energy
    plt.figure(1)
    plt.clf()
    du = np.log(energyGrid[1:] / energyGrid[:-1])
    du[0] = du[1]
    x, y = putil.get_stairs(energyGrid, du)
    plt.loglog(x, y, rasterized=True)
    #
    # dE/E ~= du for small dE/E
    xLim = [x[0], x[-1]]
    dEScat = get_min_dE(3) * np.array([1,1])
    plt.loglog(xLim, dEScat, '--')
    #
    plt.ylim([1E-7, 1E+1])
    plt.xlabel('Energy (eV)')
    plt.ylabel('Lethargy  width')
    plt.title('{0:,} total energy points'.format(len(energyGrid)))
    figureName = 'p_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=300)
    #
    #Elim = [1.0E0, 1.0E2] #eV
    #Elim = [3.0E0, 3.0E4] #eV
    Elim = [1.0E3, 1.0E6] #eV
    plt.xlim(Elim)
    plt.ylim([3E-6, 3E-1])
    figureName = 'p_zoom_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=300)

    # Plot how many points in energy you need to represent the RRR
    plt.clf()
    Elim = [1.0E-1, 1.0E5] #eV
    indexElim = find_range(energyGrid, Elim[0], Elim[1])
    subsetSize = energyGrid[indexElim[0]:indexElim[1]].size
    indexCDF = np.arange(subsetSize)
    plt.loglog(energyGrid[indexElim[0]:indexElim[1]], indexCDF, rasterized=True)
    #
    # Determine how much of the RRR you can afford to resolve using a given number of energy points
    desiredPointsList = [900, 2000, 4500, 10000]
    for desiredPoints in desiredPointsList:
        if (indexElim[0]+desiredPoints) > subsetSize:
            maxEnergyUsingDesiredPoints = Elim[1]
        else:
            maxEnergyUsingDesiredPoints = energyGrid[indexElim[0]+desiredPoints]
        x = [maxEnergyUsingDesiredPoints, maxEnergyUsingDesiredPoints]
        y = [1, desiredPoints]
        plt.loglog(x, y, '--', color='g')
        x = [Elim[0], maxEnergyUsingDesiredPoints]
        y = [desiredPoints, desiredPoints]
        plt.loglog(x, y, '--', color='g', label='{0:,} points'.format(desiredPoints))
    # Plot Master's points as reference
    refEnergies = [3.0, 55.6, 1000.0] #eV
    for refEnergy in refEnergies:
        x = [refEnergy, refEnergy]
        y = plt.ylim()
        plt.loglog(x, y, ':', color='k', label='{0:.1e} eV'.format(refEnergy))
    #
    # Finish plotting
    plt.xlabel('Energy (eV)')
    plt.ylabel('Cumulative required energy points')
    plt.legend(loc='best')
    figureName = 'p_cdf_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=300)

def plot_sigma(materials, materialIndexDict, energyGrid, sigmaMat, baseName='test', yName='Cross section'):
    # Set the output plot directory
    figureDirr = get_common_directories()['figures/indicators']
    #
    plt.figure(1)
    plt.clf()
    for material in materials:
        shortName = material.shortName
        xsIndex = materialIndexDict[shortName]
        plt.loglog(energyGrid, sigmaMat[xsIndex, :], 'o-', markersize=1, label=shortName, rasterized=True)
    plt.xlabel('Energy (eV)')
    plt.ylabel(yName)
    plt.title('{0:,} total energy points'.format(len(energyGrid)))
    #
    figureName = 'p_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=500)
    #
    leg = plt.legend(loc='best')
    figureName = 'p_legend_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=500)
    #
    #plt.xlim([3.0E0, 3.0E4])
    leg.set_visible(False)
    #plt.xlim([6E-1, 1.E2])
    plt.xlim([1.0E3, 1.0E7])
    figureName = 'p_zoom_{0}.pdf'.format(baseName)
    figurePath = os.path.join(figureDirr, figureName)
    plt.savefig(figurePath, dpi=500)

def plot_flux(energyGrid, fluxMat, rrrRange):
    # Set the output plot directory
    figureDirr = get_common_directories()['figures/indicators']
    #
    numPoints, numGroups = fluxMat.shape
    fluxMat = fluxMat.copy()
    for ip in range(numPoints):
        plt.figure(1)
        plt.clf()
        strt, end = find_range(energyGrid, rrrRange[0], rrrRange[1])
        energy = energyGrid[strt:end]
        flux = fluxMat[ip, strt:end]
        flux /= np.max(flux)
        plt.loglog(energy, flux, rasterized=True)
        plt.xlim(rrrRange[0], rrrRange[1])
        #yMax = np.power(10, np.ceil(np.log10(np.max(flux))+2))
        yMax = 1E2
        plt.ylim(1E-10*yMax, yMax)
        plt.xlabel('Energy (eV)')
        plt.ylabel('E * flux(E)')
        # Add dpi
        plotName = 'p_flux_{0}'.format(ip)
        plotPath = os.path.join(figureDirr, plotName)
        plt.savefig(plotPath, dpi=400)

def plot_output(materials, energyGrid, outputMat, groupBdrs, pwResFactor, outputType):
    '''Plot the output for each material'''
    # Set the output plot directory
    figureDirr = get_common_directories()['figures/indicators']
    #
    doIntegral = False
    avgOutputMat = calc_avg_from_pw(energyGrid, outputMat, groupBdrs, doIntegral)
    plt.figure(5)
    plt.clf()
    for iMat, material in enumerate(materials):
        materialName = material.shortName
        x, y = putil.get_stairs(groupBdrs, avgOutputMat[iMat, :])
        plt.loglog(x, y, marker='', label=materialName, rasterized=True)
    plt.xlabel('Energy (eV)')
    if outputType == 'flux':
        plt.ylabel('Flux (arb./eV)')
    else:
        plt.ylabel('Total cross section (1/cm-eV)')
    plt.xlim(np.min(groupBdrs), np.max(groupBdrs))
    filename = 'p_out.pdf'
    filePath = os.path.join(figureDirr, filename)
    plt.savefig(filePath, dpi=400)
    plt.legend(loc='best')
    filename = 'p_legend_out.pdf'
    filePath = os.path.join(figureDirr, filename)
    plt.savefig(filePath, dpi=300)

###############################################################################
def save_material_fluxes(materials, energyGrid, fluxMat, groupBdrs, pwResFactor, workOpt, plotOutput=False):
    '''Save the MG infinite-medium fluxes for each material. Save for the entire energy range.
    Notice they have already been multiplied by the long-range energy shapes, such as 1/E.'''
    # Set the output dat directory
    datDirr = get_common_directories()['dat/indicators']
    groupBdrs = groupBdrs[::-1]
    doIntegral = True
    mgFluxMat = calc_avg_from_pw(energyGrid, fluxMat, groupBdrs, doIntegral)
    #mgEnergyGrid, mgFluxMat = form_dual(energyGrid, fluxMat)
    # Flip to being ordered descendingly
    mgEnergyGrid = groupBdrs[::-1]
    mgFluxMat = mgFluxMat[:, ::-1]
    for iMat, material in enumerate(materials):
        materialName = material.shortName
        if workOpt == 'fluxe':
            filename = 'inf_flux_{0}_e_{1}.txt'.format(materialName, pwResFactor)
        elif workOpt == 'flux':
            filename = 'inf_flux_{0}_{1}.txt'.format(materialName, pwResFactor)
        else:
            filename = 'wgt_{0}_{1}.txt'.format(materialName, pwResFactor)
        filePath = os.path.join(datDirr, filename)
        np.savetxt(filePath, np.transpose((mgEnergyGrid[:-1], mgFluxMat[iMat, :])), header='upper_energy mg_flux')
    if plotOutput:
        plot_output(materials, energyGrid, fluxMat, groupBdrs, pwResFactor, 'flux')

def save_material_xs(materials, energyGrid, totalXSMat, groupBdrs, pwResFactor, plotOutput=False):
    '''Save the MG XS for each material. Save for the entire energy range.'''
    # Set the output dat directory
    datDirr = get_common_directories()['dat/indicators']
    # Multiply back by E
    groupBdrs = groupBdrs[::-1]
    doIntegral = False
    mgXSMat = calc_avg_from_pw(energyGrid, totalXSMat, groupBdrs, doIntegral)
    # Flip to being ordered descendingly
    mgEnergyGrid = groupBdrs[::-1]
    mgXSMat = mgXSMat[:, ::-1]
    for iMat, material in enumerate(materials):
        materialName = material.shortName
        filename = 'tot_xs_{0}_{1}.txt'.format(materialName, pwResFactor)
        filePath = os.path.join(datDirr, filename)
        np.savetxt(filePath, np.transpose((mgEnergyGrid[:-1], mgXSMat[iMat, :])), header='upper_energy mg_xs')
    if plotOutput:
        plot_output(materials, energyGrid, totalXSMat, groupBdrs, pwResFactor, 'xs')

def form_reference_group_structure(resolvedPWGrid, coarseGroupBdrs, rrrRange, outPath):
    '''Given a set of pointwise energies and boundaries for the RRR, create a MG group structure'''
    resolvedMGBdrs = resolvedPWGrid
    #resolvedMGBdrs = form_dual_of_pw_grid(resolvedPWGrid) # Unnecessary
    strt, end = find_range(resolvedMGBdrs, rrrRange[0], rrrRange[1])
    rrrGrid = np.union1d([rrrRange[0], rrrRange[1]], resolvedMGBdrs[strt:end])
    coarseGroupBdrs = coarseGroupBdrs[::-1]
    #groupBdrs = insert_group_structure(rrrGrid, coarseGroupBdrs)
    groupBdrs = insert_group_structure(rrrGrid, coarseGroupBdrs, 0.0) # Hack!!!
    # TODO: Insert eSpacing here?
    write_egrid(outPath, groupBdrs)
    return groupBdrs

###############################################################################

def define_input_parser():
    import argparse
    #
    parser = argparse.ArgumentParser(description='Calculator of indicators.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    defaults = define_defaults()
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=defaults['verbosity'], choices=[0,1,2,3,4], type=int)
    parser.add_argument('-p', '--plot', help='Make various diagnostic plots depending on desired level of prolificity', nargs='?', const=1, type=int, choices=[0,1,2,3], default=0)
    parser.add_argument('-f', '--fastspectrumparam', help="The 'a' and 'b' parameters to determine the Maxwellian used for the fast spectrum. See MCNP manual for examples.", type=float, nargs=2, default=defaults['fastspectrumparam'])
    parser.add_argument('-T', '--totalonly', help='If specified, do not calculate downscattering.', action='store_true', default=False)
    parser.add_argument('-Z', '--Zlow', help='If specified, add low-Z components to the elastic scattering source and turn off the 1/E fixed source. If not specified, the elastic scattering source comes only from high-Z components and low-Z components are taken into account with a 1/E fixed source.', action='store_true', default=False)
    parser.add_argument('-w', '--workopt', help='What to do. flux means calculate infinite-medium flux. fluxe means calculate infinite-medium flux with escape cross section. sigt means calculate total macroscopic cross section. flux and fluxe store psi*dE/M(E), where M(E)=1/E in the RRR. wgt means do fluxe but store the MG flux (psi*dE).', choices=['flux', 'fluxe', 'sigt', 'wgt'], default=defaults['workopt'])
    parser.add_argument('-R', '--resolvedresrange', dest='rrr', help='Energy range of resolved resonance region.', type=float, nargs=2, default=defaults['rrr'])
    parser.add_argument('-G', '--groupopt', help="Base coarse group structure file head. If the flux work option is used, the output group structure is the base coarse group structure with the RRR overwritten by the hyperfine group structure. Some examples include 'c5g7', 'scale-44', and 'res-N' where N=1,..9. Always looks in default directory of makegroups.", default=defaults['groupopt'])
    parser.add_argument('-r', '--resolution', help='Resolution to use for the pointwise flux calculations', type=int, choices=range(11), default=defaults['resolution'])
    parser.add_argument('-E', '--energyspacing', help='The maximum energy jump in the final grid is equal to this factor multiplied by a jump normalization based on downscattering distance off a heavy nucleus.', type=float, default=defaults['energyspacing'])
    parser.add_argument('-m', '--materialopt', help=" Unless 'manual' is used, specifies a set of materials to use. If 'manual' is used, give a space-separated list of material names in 'listmaterials'.", choices=['3','4','5','c5g7', 'graphite', 'iron', 'kpin', 'kenrichedpin', 'kcladpin', 'kpin2d', 'kenrichedpin2d', 'kmoxpin2d', 'kmoxenrichedpin2d', 'trigafuel', 'ctrigafuel', 'trigamore', 'deb', 'manual'], default=defaults['materialopt'])
    parser.add_argument('-i', '--indicatormaterials', dest='listmaterials', help="When manual 'materialopt' is used, specify the materials to use.", nargs='+', default=defaults['listmaterials'], choices=mat.get_materials_name2function_dict().keys())
    parser.add_argument('-t', '--temperaturedependence', help="If specified, use the temperature in the PENDF file that corresponds to the temperature of the material (as specified in 'materials_materials.py'). If this temperature does not exist, an error is thrown. If not specified, the first temperature in the PENDF file is used. This flag is added for 'materialopt' of '3'.", action='store_true', default=False)
    parser.add_argument('-M', '--meshname', help="The energy mesh on which the output fluxes are constructed. {r} is replaced by the resolution. If the 'flux' workopt is used, the mesh file is written to. Else, is it read from. If no file type is specified, '.txt' is used. If no directory is specified, the 'energyDat' directory from directories is used.", default=defaults['meshname'])
    return parser

###############################################################################
if __name__ == '__main__':
    parser = define_input_parser()
    inputDict = vars(parser.parse_args())
    if inputDict['verbosity'] > 1:
        print 'Summary of inputs:'
        print inputDict
    #test_grid_thinning_and_thickening()
    do_all(inputDict)

