#! /usr/bin/env python

'''
Andrew Till
Fall 2014
PhD Research

This clustering module clusters indicators and uses labels to compute a generalized energy mesh:
    Create observations from subset of indicators, viz. indicators[strt:end]
    Cluster observations and extract labels
    Combine labels from each subset
    Merge labels into existing group structure
    Create generalized energy mesh using the dual of the input grid
'''

#STDLIB
import os
import time

#TPL
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from sklearn import cluster
from sklearn import neighbors
import matplotlib as mpl
#import matplotlib.ticker as mtick
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 18, 'lines.linewidth': 2})

#MINE
import plotutil as putil
from directories import get_common_directories

def define_defaults():
    '''Specify default parameters'''
    # Main parameters
    verbosity = False
    resolution = 9
    #workOpt = 'amg'
    #workOpt = 'tmg'
    #workOpt = 'mg'
    workOpt = 'har'
    numElements = 64

    # Misc parameters
    showNumbers = False
    condenseSubelements = True
    plotOpt = 'half' # none, first, last, firstlast, all, half
    energyPenalty = 3.6E-1
    dpi = 100

    # Specify range of interest (for output fluxes and plotting)
    numCoarseGroups = 0
    #coarseBdrs = [6E-1, 3E0, 54.7, 2E3, 2.5E4]
    coarseBdrs = [3.0, 54.7, 1.06E3]

    # How to assign the number of elements per coarse group
    apportionOpt = 'equal'

    # Specify how many elements to use in each coarse group
    listNumElements = []

    # Specify which set of materials to use
    #materialOpt = '4'
    #materialOpt = 'c5g7'
    materialOpt = 'manual'

    # If materialOpt is 'manual', list of materials to use
    materialsList = ['deb']
    importancesLst = []

    return {'verbosity': verbosity, 'numelements': numElements, 'apportionopt': apportionOpt, 'workopt': workOpt, 'resolution':resolution, 'coarsebdrs': coarseBdrs, 'numcoarsegroups': numCoarseGroups, 'materialopt': materialOpt, 'shownumbers': showNumbers, 'condensesubelements': condenseSubelements, 'plotopt': plotOpt, 'energypenalty': energyPenalty, 'dpi': dpi, 'listnumelements': listNumElements, 'listmaterials': materialsList}

def do_all(inputDict):
    '''Driver function to read indicators, compute observations, perform clustering, write energy mesh'''
    ''' >>>>> NB: See compute_map() for external calls! <<<<<< '''
    # Create timeDict to house timing results
    timeDict = {}
    # Create dataDict to house problem parameters
    # Future functions should be assumed to modify dataDict
    dataDict = {}
    # Initialize dataDict
    copy_inputDict_to_dataDict(inputDict, dataDict)
    populate_directories(dataDict)
    # Determine work options
    parse_work_opt(dataDict)
    # Read in energy mesh and indicators
    indicators = read_energy_and_indicators(dataDict)
    # Define observations and rescale energy; determine va
    observations = compute_observations(dataDict, indicators)
    # Determine coarse group structure and the number of energy elements per coarse group
    globalLabels, thermalOffset = apportion_elements(dataDict, observations)
    # Determine number of neighbors for the clustering
    find_num_neighbors(dataDict, timeDict, observations)
    # Loop over each coarse group
    offset = thermalOffset
    print_timing_header()
    for coarseGroup in range(len(dataDict['coarseBdrs'])-1):
        # Cluster within each coarse group and plot
        offset = cluster_one_coarse_group(
            dataDict, timeDict, observations, globalLabels, offset, coarseGroup)
    # Create one plot over all energy
    plot_summary(dataDict, observations, globalLabels)
    # Write out generalized energy mesh
    write_mesh(dataDict, globalLabels)



###############################################################################
def copy_inputDict_to_dataDict(inputDict, dataDict):
    '''Capitalize entries in dataDict using camelCase'''
    # Inputs / outputs
    dataDict['verbosity'] = inputDict['verbosity']
    dataDict['workOpt'] = inputDict['workopt']
    dataDict['useSigt'] = inputDict['sigt']
    dataDict['resolution'] = inputDict['resolution']
    dataDict['materialOpt'] = inputDict['materialopt']
    dataDict['materialsList'] = inputDict['listmaterials']
    dataDict['importancesList'] = inputDict['listimportances']
    dataDict['numElements'] = inputDict['numelements']
    dataDict['numElementsIsTotal'] = inputDict['numelementsistotal']
    dataDict['numElementsList'] = inputDict['listnumelements']
    dataDict['apportionOpt'] = inputDict['apportionopt']
    dataDict['coarseBdrs'] = inputDict['coarsebdrs']
    dataDict['numCoarseGroups'] = inputDict['numcoarsegroups']
    dataDict['showNumbers'] = inputDict['shownumbers']
    dataDict['condenseSubelements'] = inputDict['condensesubelements']
    dataDict['plotOpt'] = inputDict['plotopt']
    dataDict['energyPenalty'] = inputDict['energypenalty']
    dataDict['dpi'] = inputDict['dpi']

def populate_directories(dataDict):
    '''Add directory entries to dataDict'''
    # Inputs
    dirDict = get_common_directories()
    # Outputs
    dataDict['indicatorsDatDirr'] = dirDict['dat/indicators']
    dataDict['energyDatDirr'] = dirDict['dat/energy_groups']
    dataDict['plotDirr'] = dirDict['figures/clustering']

def parse_work_opt(dataDict):
    '''Populate work option parameters'''
    # Inputs
    workOpt = dataDict['workOpt']

    forceContiguity = False
    useEqualMGSpacing = False
    useEqualIndexSpacing = False
    if workOpt == 'amg':
        forceContiguity = True
    elif workOpt == 'tmg':
        forceContiguity = True
        useEqualIndexSpacing = True
    elif workOpt == 'mg':
        forceContiguity = True
        useEqualMGSpacing = True
    if not forceContiguity:
        clusterName = 'har'
    elif useEqualMGSpacing:
        clusterName = 'mg'
    elif useEqualIndexSpacing:
        clusterName = 'tmg'
    else:
        clusterName = 'amg'

    # Outputs
    dataDict['forceContiguity'] = forceContiguity
    dataDict['useEqualMGSpacing'] = useEqualMGSpacing
    dataDict['useEqualIndexSpacing'] = useEqualIndexSpacing
    dataDict['clusterName'] = clusterName

def read_energy_and_indicators(dataDict):
    '''Read input fine energy mesh and indicators that live on that mesh'''
    # Inputs
    resolution = dataDict['resolution']
    materialOpt = dataDict['materialOpt']
    materialsList = dataDict['materialsList']
    importancesList = dataDict['importancesList']
    useSigt = dataDict['useSigt']
    energyDatDirr = dataDict['energyDatDirr']
    indicatorsDatDirr = dataDict['indicatorsDatDirr']

    # Read in hyperfine MG energy mesh
    energyName = 'res-{0}.txt'.format(resolution)
    energyPath = os.path.join(energyDatDirr, energyName)
    groupBdrs = load_txt(energyPath, skiprows=2, usecols=[1])
    numGroups = len(groupBdrs) - 1
    dE =  - np.diff(groupBdrs)
    groupBdrs[groupBdrs == 0] = 1E-5

    # Read indicator relative importances
    if materialOpt == 'c5g7':
        #importanceDict = {'cUO2': 4, 'clowMOX': 2, 'cmedMOX': 2, 'chighMOX': 2, 'cCR': 1}
        importanceDict = {'cUO2': 1, 'chighMOX': 1}
    elif materialOpt == 'graphite':
        importanceDict = {'graphite': 1}
    elif materialOpt == 'iron':
        importanceDict = {'iron': 1}
    elif materialOpt == 'kpin':
        importanceDict = {'kFUEL': 1}
    elif materialOpt == 'kenrichedpin':
        importanceDict = {'kEFUEL': 4, 'kFUEL': 1}
    elif materialOpt == 'kcladpin':
        importanceDict = {'kEFUEL': 10, 'kFUEL': 2.5, 'kZR': 1}
    elif materialOpt == 'kpin2d':
        importanceDict = {'kRFUEL': 1}
    elif materialOpt == 'kenrichedpin2d':
        importanceDict = {'kREFUEL': 4, 'kRFUEL': 1}
    elif materialOpt == 'kmoxpin2d':
        importanceDict = {'kRMFUEL': 4, 'kRFUEL': 1}
    elif materialOpt == 'kmoxenrichedpin2d':
        importanceDict = {'kRMFUEL': 4, 'kREFUEL': 4, 'kRFUEL': 1}
    elif materialOpt == 'trigafuel':
        importanceDict = {'tFUEL': 1}
    elif materialOpt == 'ctrigafuel':
        importanceDict = {'tcFUEL': 1}
    elif materialOpt == 'ctrigafuel_0':
        importanceDict = {'tdFUEL_0': 1}
    elif materialOpt == 'trigamore':
        importanceDict = {'tFUEL': 10, 'tCLAD': 2, 'tZIRC': 2, 'tIRRADIATIONTUBE': 1}
    elif materialOpt == 'manual':
        if not importancesList:
            importancesList = [1]*len(materialsList)
        importanceDict = {material:importance for material,importance in zip(
            materialsList, importancesList)}

    # Appending to materialNames does not alias onto importanceDict
    materialNames = importanceDict.keys()
    numMaterials = len(materialNames)
    # Use infinite medium flux, infinite medium flux with escape xs, and energy itself
    numIndicators = 2 * numMaterials + 1

    # Read indicators
    indicators = np.zeros((numIndicators, numGroups))
    i = 0
    for materialName in materialNames:
        weight = np.power(10, importanceDict[materialName])
        fluxName = 'inf_flux_{0}_{1}.txt'.format(materialName, resolution)
        fluxPath = os.path.join(indicatorsDatDirr, fluxName)
        indicators[i, :] = weight * load_txt(fluxPath) / dE
        if useSigt:
            fluxName = 'tot_xs_{0}_{1}.txt'.format(materialName, resolution)
            fluxPath = os.path.join(indicatorsDatDirr, fluxName)
            indicators[i+1, :] = weight * load_txt(fluxPath)
        else:
            fluxName = 'inf_flux_{0}_e_{1}.txt'.format(materialName, resolution)
            fluxPath = os.path.join(indicatorsDatDirr, fluxName)
            indicators[i+1, :] = weight * load_txt(fluxPath) / dE
        i += 2
    energyGrid = np.sqrt(groupBdrs[1:] * groupBdrs[:-1])
    indicators[-1, :] = energyGrid
    materialNames += 'E'

    # Outputs
    dataDict['dE'] = dE
    dataDict['groupBdrs'] = groupBdrs
    dataDict['energyGrid'] = energyGrid
    dataDict['materialNames'] = materialNames
    dataDict['numGroups'] = numGroups
    dataDict['numIndicators'] = numIndicators
    return indicators

def compute_observations(dataDict, indicators):
    '''Compute observations from indicators.
    Only uses first and last values of coarseBdrs, which are not changed in apportion_elements()'''
    # Inputs
    energyPenalty = dataDict['energyPenalty']
    groupBdrs = dataDict['groupBdrs']
    coarseBdrs = dataDict['coarseBdrs']

    numIndicators, numGroups = indicators.shape
    observations = np.log10(indicators)
    medians = np.median(observations, axis=1)
    observations -= medians[:, np.newaxis]
    #
    obsRange = np.max(observations[:-1,:]) - np.min(observations[:-1,:])
    strt = np.argmin(np.abs(groupBdrs - coarseBdrs[-1]))
    end = np.argmin(np.abs(groupBdrs - coarseBdrs[0]))
    energyRange = np.max(observations[-1,strt:end]) - np.min(observations[-1,strt:end])
    observations[-1,:] *= np.sqrt(numIndicators) * energyPenalty * (obsRange / energyRange)

    # Transpose to be [numGroups, numPoints]. Copy for later speed
    observations = observations.transpose().copy()

    # Outputs
    return observations

def apportion_elements(dataDict, observations):
    '''Determine the coarse group boundaries and the number of elements in each coarse group.
    The first and last coarse group boundaries determine the extent of the RRR.
    coarseBdrs are replaced with an equal lethargy grid if numCoarseGroups is nonzero.
    numElementsList has precedence over numElements/apportionOpt if both are given.
    numElements refers to the elements in the RRR unless numElementsIsTotal is True.
    '''
    # Inputs
    verbosity = dataDict['verbosity']
    coarseBdrs = dataDict['coarseBdrs']
    numCoarseGroups = dataDict['numCoarseGroups']
    numElementsList = dataDict['numElementsList']
    groupBdrs = dataDict['groupBdrs']
    apportionOpt = dataDict['apportionOpt']
    numElements = dataDict['numElements']
    numElementsIsTotal = dataDict['numElementsIsTotal']
    numGroups = dataDict['numGroups']
    indicatorsDatDirr = dataDict['indicatorsDatDirr']

    # Determine the coarse group boundaries
    if numCoarseGroups != 0:
        # Use an equal lethargy grid with numCoarseGroups groups
        coarseBdrs = np.logspace(np.log10(coarseBdrs[0]), np.log10(coarseBdrs[-1]), numCoarseGroups+1)
    else:
        numCoarseGroups = len(coarseBdrs) - 1

    # Determine the total number of elements
    if numElementsList:
        numElementsIsTotal = False
        numElements = np.sum(numElementsList)
    numElementsIsRRR = not(numElementsIsTotal)

    # Determine the RRR bounding indices on the fine group structure
    strt = np.argmin(np.abs(groupBdrs - coarseBdrs[-1]))
    end = np.argmin(np.abs(groupBdrs - coarseBdrs[0]))

    # Initialize the global labels, which map from subelement index to element index
    globalLabels = np.zeros(numGroups, dtype=np.int)
    numThermal = numGroups - end
    numFast = strt
    numFastAndThermal = numThermal + numFast
    # On the old group structure,
    #   [:strt) are fast, [end:) are thermal, and [strt:end) are resonance
    # On the new element structure,
    #   [:newStrt) are fast, [:newEnd) are thermal, and [newStrt:newEnd) are resonance
    newStrt = strt
    if numElementsIsRRR:
        # numElements is just the number of resonance elements
        newEnd = newStrt + numElements
    else:
        # numElements is the total number of elements:
        newEnd = end + (numElements - numGroups)
    thermalOffset = newEnd
    globalLabels[:strt] = np.arange(numFast)
    globalLabels[end:] = np.arange(numThermal) + thermalOffset

    # Determine the total number of elements
    numElementsTot = numElements
    if numElementsIsRRR:
        numElementsTot += numFastAndThermal
    numElementsRRR = numElementsTot - numFastAndThermal

    # Determine the number of elements per coarse group and store in numClustersList
    if numElementsList:
        numClustersList = numElementsList
        apportionOpt = 'manual'
    elif apportionOpt in ['var', 'max', 'L1', 'birch']:
        numClustersList = auto_apportion(observations, numElementsRRR, groupBdrs, coarseBdrs, apportionOpt, verbosity)
    else: # apportionOpt == 'equal'
        # Apportion as equally as possible. Give high-energy coarse groups the remainder
        numElementsPerCoarseGroup = numElementsRRR // numCoarseGroups
        remainder = numElementsRRR % numCoarseGroups
        numClustersList = numElementsPerCoarseGroup * np.ones(numCoarseGroups, dtype=np.int)
        if remainder:
            numClustersList[-remainder:] += 1

    # Check for validity of the number of clusters list
    if np.any(numClustersList <= 0):
        minNumElements = numFastAndThermal + len(coarseBdrs) - 1
        raise ValueError('{0} clusters specified, but at least {1} should have been used'.format(numElementsTot, minNumElements))

    # Output 0
    # Print the number of elements per coarse group
    print 'final elements per coarse group ({0}):\n'.format(apportionOpt), numClustersList

    # Output 1
    # Save the number of elements per coarse group
    baseName = 'aptn_{0}'.format(apportionOpt)
    #filename = '{0}_{1}_{2}.txt'.format(baseName, numElements, resolution)
    filename = '{0}_e{1}_g{2}.txt'.format(baseName, numElements, numCoarseGroups)
    filePath = os.path.join(indicatorsDatDirr, filename)
    # Output
    with open(filePath, 'w') as fid:
        fid.write('# Energy mesh with {0} coarse groups and {1} resonance elements\n'.format(numCoarseGroups, numElementsRRR))
        fid.write('# coarse_group upper_bound(eV) num_elements\n')
        for coarseGroup in range(numCoarseGroups):
            energy = coarseBdrs[coarseGroup]
            numElem = numClustersList[coarseGroup]
            # Write the same number of digits as NJOY. This prints from coarseGroup, not groupBdrs.
            fid.write('{0:g} {1:.6e} {2}\n'.format(coarseGroup, energy, numElem))
        fid.write('{0:g} {1:.6e} {2}\n'.format(-1, coarseBdrs[-1], 0))

    # Output 2
    dataDict['apportionOpt'] = apportionOpt
    dataDict['numCoarseGroups'] = numCoarseGroups
    dataDict['coarseBdrs'] = coarseBdrs
    dataDict['numElementsIsTotal'] = numElementsIsTotal
    dataDict['numElementsList'] = numElementsList
    dataDict['numElements'] = numElements
    dataDict['numElementsTot'] = numElementsTot
    dataDict['numClustersList'] = numClustersList
    return globalLabels, thermalOffset

def auto_apportion(observations, numElementsRRR, groupBdrs, coarseBdrs, apportionOpt, verbosity):
    '''Assign the number of clusters / elements proportional to the relative variance within a coarse group'''

    numCoarseGroups = len(coarseBdrs) - 1
    metric = np.zeros(numCoarseGroups)
    numFineGroupsPerCoarseGroup = np.zeros(numCoarseGroups, dtype=np.int)
    timeBirch = 0.0
    if numCoarseGroups == 1:
        return np.array([numElementsRRR])
    for coarseGroup in range(len(coarseBdrs)-1):
        # Only look at the fine groups within the current coarse group
        strt = np.argmin(np.abs(groupBdrs - coarseBdrs[coarseGroup+1]))
        end = np.argmin(np.abs(groupBdrs - coarseBdrs[coarseGroup]))
        obs = observations[strt:end, :]
        # Points are groups and dim are materials
        numPoints, numDim = obs.shape
        numFineGroupsPerCoarseGroup[coarseGroup] = numPoints
        if apportionOpt == 'var':
            # For each dimension, compute the average over the points
            means = np.mean(obs, axis=0)
            # The variance metric is the sum of the square errors from the means
            # var is actually a standard deviation
            var = np.sqrt(np.sum(np.square(obs - means[np.newaxis,:])) / (numPoints * numDim))
            metric[coarseGroup] = var
        elif apportionOpt == 'max':
            # The max error metric is the maximum obs range in over all dimensions,
            # where each obs range is the difference between the largest and smallest
            # point values for that dim
            maxErr = np.max(np.max(obs, axis=0) - np.min(obs, axis=0))
            metric[coarseGroup] = maxErr
        elif apportionOpt == 'L1':
            # The L1 error metric is the unnormalized sum of the absolute pointwise differences,
            # where each pointwise difference is a maximum over all dimensions
            L1Err = np.sum(np.max(np.abs(obs[:-1,:] - obs[1:,:]), axis=1))
            metric[coarseGroup] = L1Err
        elif apportionOpt == 'birch':
            # The BIRCH error metric uses the number of clusters required to reduce the variance of an
            # energy element to threshold. A value of 1.0 for threshold would produce clusters
            # that span a factor of 10 in flux / xs space (based on how obs are defined).
            # Do not include energy penalty when looking at clusters.
            threshold = 0.05
            clusterer = cluster.Birch(n_clusters=None, threshold=threshold)
            t0 = time.time()
            numClusters = len(np.unique(clusterer.fit_predict(obs[:,:-1])))
            timeBirch += time.time() - t0
            # Converted to float implicitly
            metric[coarseGroup] = numClusters
    # Normalize the metric
    metric /= np.sum(metric)
    # Each coarse group needs at least one energy element and the number of energy elements
    # per coarse group should be proportional to the normalized metric (relative variance)
    desiredElementsPerCoarseGroup = np.minimum(numFineGroupsPerCoarseGroup,
            np.maximum(1, metric * numElementsRRR))
    sumDesiredElements = np.sum(desiredElementsPerCoarseGroup)
    # Compute a new metric that takes into account the maximum(1,:)
    newMetric = (desiredElementsPerCoarseGroup - 1) / (sumDesiredElements - numCoarseGroups)
    fractions, floored = np.modf(newMetric * (numElementsRRR - numCoarseGroups))
    elementsPerCoarseGroup = np.array(floored, dtype=np.int) + 1
    # Add elements to the largest coarse groups with the largest fractional number of elements
    remainder = numElementsRRR - np.sum(elementsPerCoarseGroup)
    if remainder:
        locLargestFrac = np.argsort(fractions)[-remainder:]
        elementsPerCoarseGroup[locLargestFrac] += 1
    if verbosity and apportionOpt == 'birch':
        print 'Birch clustering took {0} s'.format(timeBirch)
    if verbosity:
        print 'desired elements per coarse group:\n', desiredElementsPerCoarseGroup
        if remainder:
            print 'remainder: {0}; smallest large fractional: {1:.3f}'.format(
                    remainder, fractions[locLargestFrac[0]])

    # Outputs
    return elementsPerCoarseGroup

def find_num_neighbors(dataDict, timeDict, observations):
    '''Use connected components to determine the number of neighbors'''
    # (no inputs from dataDict)

    t0 = time.time()
    numMinNeighbors = find_minimum_num_neighbors(observations[:,:-1])
    numNeighbors = max(15+numMinNeighbors, 200) # Changed from 100
    timeInitialNeighbors = time.time() - t0

    # Outputs
    dataDict['numNeighbors'] = numNeighbors
    timeDict['timeInitialNeighbors'] = timeInitialNeighbors

def cluster_one_coarse_group(dataDict, timeDict, observations, globalLabels, offset, coarseGroup):
    '''Apply clustering to one coarse group and return the labels of the clustering'''
    # Inputs
    groupBdrs = dataDict['groupBdrs']
    coarseBdrs = dataDict['coarseBdrs']
    energyGrid = dataDict['energyGrid']
    numClustersList = dataDict['numClustersList']
    numNeighbors = dataDict['numNeighbors']
    forceContiguity = dataDict['forceContiguity']
    useEqualMGSpacing = dataDict['useEqualMGSpacing']
    useEqualIndexSpacing = dataDict['useEqualIndexSpacing']
    plotOpt = dataDict['plotOpt']

    # Slice arrays for current coarse group
    strt = np.argmin(np.abs(groupBdrs - coarseBdrs[coarseGroup+1]))
    end = np.argmin(np.abs(groupBdrs - coarseBdrs[coarseGroup]))
    obs = observations[strt:end,:]
    eGrid = energyGrid[strt:end]
    gBdr = groupBdrs[strt:end+1]
    numGroups, numPoints = obs.shape
    numClusters = numClustersList[coarseGroup]

    # Create connectivity graph based on number of neighbors
    t0 = time.time()
    useNeighbors = min(numGroups, numNeighbors)
    knnGraph = neighbors.kneighbors_graph(obs[:,:-1], useNeighbors, include_self=True)
    timeNeighbors = time.time() - t0
    connectivityGraph = knnGraph

    # Perform appropriate clustering for the current coarse group
    t0 = time.time()
    if not forceContiguity:
        labels = cluster_using_hierarchical_agglomeration(obs, numClusters, connectivityGraph)
    elif useEqualMGSpacing:
        labels = cluster_using_equal_energy_spacing(gBdr, numClusters)
    elif useEqualIndexSpacing:
        labels = cluster_using_equal_topological_spacing(obs, numClusters)
    else:
        labels = cluster_using_mg_squared_error(obs, numClusters)
    timeCluster = time.time() - t0

    # Reorder labels for maximal downscattering (make labels descendingly sorted)
    temp = np.zeros((numClusters, 1))
    labels, temp = reorder_codebook(labels, temp, eGrid)

    # Output 0: Set global labels based on local labels and global offset (aliased)
    offset -= numClusters
    globalLabels[strt:end] = labels + offset
    uniqueLabels = np.unique(labels)

    # Output 1: Print number of neighbors and required times
    timeDict['timeNeighbors'] = timeNeighbors
    timeDict['timeCluster'] = timeCluster
    print_timing(dataDict, timeDict, numGroups, coarseGroup)

    # Output 2: Plot observations colored by label / cluster
    if plotOpt not in ['none', 'sum']:
        plot_clustering(dataDict, coarseGroup, uniqueLabels, labels, eGrid, obs, numClusters, numPoints, offset)

    # Output3:
    return offset

###############################################################################
def print_timing_header():
    '''Print the header for print_timing()'''
    # Output
    print 'coarseGroup, numNeighbors, fracNeighbors, numGroups, timeInitialNeighbors, timeNeighbors, timeCluster'

def print_timing(dataDict, timeDict, numGroups, coarseGroup):
    '''Print interesting size and time information'''
    # Inputs
    numNeighbors = dataDict['numNeighbors']
    timeInitialNeighbors = timeDict['timeInitialNeighbors']
    timeNeighbors = timeDict['timeNeighbors']
    timeCluster = timeDict['timeCluster']
    # Output
    print coarseGroup, numNeighbors, float(numNeighbors) / numGroups, numGroups, timeInitialNeighbors, timeNeighbors, timeCluster

def plot_summary(dataDict, observations, globalLabels):
    '''Plot the entire energy range range'''
    # Inputs
    groupBdrs = dataDict['groupBdrs']
    coarseBdrs = dataDict['coarseBdrs']
    energyGrid = dataDict['energyGrid']
    numElements = dataDict['numElements']

    # Slice arrays for the RRR
    strt = np.argmin(np.abs(groupBdrs - coarseBdrs[-1]))
    end = np.argmin(np.abs(groupBdrs - coarseBdrs[0]))
    obs = observations[strt:end,:]
    eGrid = energyGrid[strt:end]
    gBdr = groupBdrs[strt:end+1]
    numGroups, numPoints = obs.shape
    labels = globalLabels[strt:end].copy()
    offset = np.min(labels)
    labels -= offset
    uniqueLabels = np.unique(labels)

    numClusters = numElements
    coarseGroup = 'sum'
    # Output
    plot_clustering(dataDict, coarseGroup, uniqueLabels, labels, eGrid, obs, numClusters, numPoints, offset)

def plot_clustering(dataDict, coarseGroup, uniqueLabels, labels, eGrid, obs, numClusters, numPoints, offset):
    '''Plot observations in coarse group'''
    # Inputs
    forceContiguity = dataDict['forceContiguity']
    clusterName = dataDict['clusterName']
    materialNames = dataDict['materialNames']
    numElements = dataDict['numElements']
    numCoarseGroups = dataDict['numCoarseGroups']
    plotOpt = dataDict['plotOpt']
    plotDirr = dataDict['plotDirr']
    showNumbers = dataDict['showNumbers']
    dpi = dataDict['dpi']

    if plotOpt != 'none':
        colors = putil.get_colors(max(uniqueLabels)-min(uniqueLabels)+1)
        #if forceContiguity and not showNumbers:
        if True or numClusters < 1000:
            colors = colors[np.argsort(np.random.random(colors.shape[0]))]
        if plotOpt == 'first':
            pointsToPlot = [0]
        elif plotOpt == 'last':
            pointsToPlot = [numPoints-1]
        elif plotOpt == 'firstlast':
            pointsToPlot = [0, numPoints-1]
        elif plotOpt == 'half':
            pointsToPlot = range(0, numPoints-1, 2)
        else:
            pointsToPlot = range(numPoints)
        avgLabels, avgEGrid, avgObs = average_observations(labels, eGrid, obs)
        for ip in pointsToPlot:
            material = materialNames[ip/2]
            if ip % 2 == 1:
                material += '_e'
            print 'Saving plot for {0}, {1}, {2}'.format(material, numElements, coarseGroup)
            plt.figure(3)
            plt.clf()
            labelColors = [colors[label] for label in labels]
            #plt.semilogx(eGrid, obs[:,ip], '-', color=[0.2, 0.2, 0.2], rasterized=True)
            plt.scatter(eGrid, obs[:,ip], c=labelColors, edgecolor='none', rasterized=True)
            plt.xscale('log')
            if showNumbers:
                for label in uniqueLabels:
                    mask = (avgLabels == label)
                    color = colors[label]
                    marker = r'${0}$'.format(label + offset)
                    sz = 10
                    if len(marker) == 5:
                        sz = 12
                    white = [1, 1, 1, 1.0]
                    plt.semilogx(avgEGrid[mask], avgObs[mask,ip], linestyle='', markersize=sz, color=white, markeredgecolor=white, rasterized=False, marker='o')
                    plt.semilogx(avgEGrid[mask], avgObs[mask,ip], linestyle='', markersize=sz, color=color, markeredgecolor=color, rasterized=False, marker=marker)
            if (eGrid[0] / eGrid[-1]) < 5:
                plt.xscale('linear')
            plt.xlabel('Energy (eV)')
            plt.ylabel('Observation (arb.)')
            plt.title('{0} elements'.format(len(uniqueLabels)))
            plt.xlim(np.min(eGrid), np.max(eGrid))
            # Hack the zoom
            #plt.xlim([1E3,1E7])
            #plt.xlim([7.E6, 1.E7])
            baseName = 'p_obs'
            if forceContiguity:
                baseName += '_{0}'.format(clusterName)
            effCoarseGroups = 'of_{0}'.format(numCoarseGroups - 1)
            if coarseGroup == 'sum':
                effCoarseGroups = '{0}'.format(numCoarseGroups)
            plotName = '{0}_{1}_{2}_{3}_{4}.pdf'.format(
                baseName, numElements, coarseGroup, effCoarseGroups, material)
            plotPath = os.path.join(plotDirr, plotName)
            # Output
            plt.tight_layout()
            plt.savefig(plotPath, dpi=dpi)

def write_mesh(dataDict, globalLabels):
    '''Write the energy mesh by writing subelements'''
    # Inputs
    forceContiguity = dataDict['forceContiguity']
    condenseSubelements = dataDict['condenseSubelements']
    resolution = dataDict['resolution']
    groupBdrs = dataDict['groupBdrs']
    clusterName = dataDict['clusterName']
    numElements = dataDict['numElements']
    numElementsTot = dataDict['numElementsTot']
    energyDatDirr = dataDict['energyDatDirr']

    if condenseSubelements:
        globalLabels, groupBdrs = condense_subelements(globalLabels, groupBdrs)
    numSubelements = len(globalLabels)
    baseName = 'clust'
    if forceContiguity:
        baseName += '-{0}'.format(clusterName)
    filename = '{0}-{1}-{2}.txt'.format(baseName, numElements, resolution)
    filePath = os.path.join(energyDatDirr, filename)
    # Output
    with open(filePath, 'w') as fid:
        fid.write('# Energy mesh with {0} elements and {1} subelements\n'.format(numElementsTot, numSubelements))
        fid.write('# element upper bound region\n')
        for label, energy in zip(globalLabels, groupBdrs[:-1]):
            if energy > 2.5E4:
                energyType = 'fast'
            elif energy <= 3.0:
                energyType = 'thermal'
            else:
                energyType = 'resonance'
            fid.write('{0:g} {1:.8e} {2}\n'.format(label, energy, energyType))
        fid.write('-1 {0:.8e} thermal\n'.format(groupBdrs[-1]))

###############################################################################
def compute_map(inputDict):
    '''Callable interface to use these clustering methods from another script.
    Returns the uncondensed labels'''
    # Create timeDict to house timing results
    timeDict = {}
    # Create dataDict to house problem parameters. Future functions should be assumed to modify dataDict
    dataDict = {}
    # Initialize dataDict, assuming inputDict has the observations
    # indicators sets: dE, groupBdrs, energyGrid (just energyAvg), materialNames, numGroups, numIndicators
    observations = extract_external_inputDict(inputDict, dataDict) #(not complete)
    # Determine work options
    parse_work_opt(dataDict)
    # Determine coarse group structure and the number of energy elements per coarse group
    globalLabels, thermalOffset = apportion_elements(dataDict, observations)
    # Determine number of neighbors for the clustering
    find_num_neighbors(dataDict, timeDict, observations)
    # Loop over each coarse group
    offset = thermalOffset
    print_timing_header()
    for coarseGroup in range(len(dataDict['coarseBdrs'])-1):
        # Cluster within each coarse group and plot
        offset = cluster_one_coarse_group(
            dataDict, timeDict, observations, globalLabels, offset, coarseGroup)
    return globalLabels

###############################################################################
def reorder_codebook(codebook, centroids, energyAvg, positionFunc=np.mean):
    '''Sort so group with highest value of positionFunc is cluster 0'''
    masks = get_masks(codebook)
    numGroups, numClusters = masks.shape
    energyCentroids = np.zeros(numClusters)
    newCodebook = np.zeros(codebook.size, dtype=np.int64)
    newCentroids = np.zeros(centroids.shape)
    for i in range(numClusters):
        energyCentroids[i] = positionFunc(energyAvg[masks[:,i]])
    clusterOrder = np.argsort(-energyCentroids)
    for ic, cluster in enumerate(clusterOrder):
        newCodebook[codebook==cluster] = ic
        newCentroids[ic,:] = centroids[cluster,:]
    return newCodebook, newCentroids

def get_masks(codebook):
    minCode = min(codebook)
    maxCode = max(codebook)
    numCodes = maxCode - minCode + 1
    numGroups = len(codebook)
    masks = np.zeros([numGroups, numCodes], dtype=bool)
    for i in range(minCode, maxCode+1):
        masks[:,i] = codebook==i
    return masks

###############################################################################
def cluster_using_mg_squared_error(observations, numGroups):
    '''This may return fewer groups than numGroups when numGroups is large'''
    '''If this fails to produce the desired number of groups, split up existing groups starting with LOW energies and going to HIGH'''
    # Calculate the max difference in obs over all space/angle points between neighboring energy points
    numPoints, numDim = observations.shape
    obsErr = np.zeros(numPoints)
    obsErr[1:] = np.max(np.abs(observations[:-1,:] - observations[1:,:]) , axis=1)
    # Sum this error
    cumErr = np.cumsum(np.square(obsErr))
    totErr = cumErr[-1]
    # Evenly divide the total error into numGroups groups
    errPerGroup = totErr / numGroups
    labels = np.zeros(numPoints, dtype=np.int)
    groupBdrLower = 0
    for g in range(numGroups):
        desiredErr = (g + 1) * errPerGroup
        groupBdrUpper = np.argmin(np.abs(cumErr - desiredErr)) + 1
        if groupBdrLower >= groupBdrUpper:
            groupBdrUpper = groupBdrLower + 1
        labels[groupBdrLower:groupBdrUpper] = g
        groupBdrLower = groupBdrUpper
    return labels

def cluster_using_equal_topological_spacing(observations, numGroups):
    '''Split observations into numGroups pieces so that an even number of points go into each piece'''
    numPoints, numDim = observations.shape
    pointsPerGroup = numPoints // numGroups
    remainder = numPoints % numGroups
    labels = np.zeros(numPoints, dtype=np.int)
    groupBdrLower = 0
    for g in range(numGroups):
        groupBdrUpper = groupBdrLower + pointsPerGroup
        if g < remainder:
            groupBdrUpper += 1
        labels[groupBdrLower:groupBdrUpper] = g
        groupBdrLower = groupBdrUpper
    return labels

def cluster_using_equal_energy_spacing(fineGroupBdrs, numGroups):
    '''This may return fewer groups than numGroups when fineGroupBdrs has unequal spacing and numGroups is large'''
    '''If this fails to produce the desired number of groups, split up existing groups starting with HIGH energies and going to LOW'''
    # Calculate an energy mesh using equal log spacing with numGroups groups
    numFineGroups = len(fineGroupBdrs) - 1
    desiredGroupBdrs = np.logspace(np.log10(fineGroupBdrs[0]), np.log10(fineGroupBdrs[-1]), numGroups+1)
    # Get as close as possible to this mesh using the fineGroupBdrs
    labels = np.zeros(numFineGroups, dtype=np.int)
    groupBdrLower = 0
    for g in range(numGroups):
        groupBdrUpper = np.argmin(np.abs(fineGroupBdrs - desiredGroupBdrs[g+1]))
        if groupBdrLower >= groupBdrUpper:
            groupBdrUpper = groupBdrLower + 1
        labels[groupBdrLower:groupBdrUpper] = g
        groupBdrLower = groupBdrUpper
    return labels

def cluster_using_hierarchical_agglomeration(obs, numClusters, connectivityGraph):
    '''Use hierarhical agglomerative clustering with a connectivity graph'''
    har = cluster.AgglomerativeClustering(n_clusters=numClusters, connectivity=connectivityGraph)
    har.fit(obs)
    return har.labels_
    # Hack to use Birch:
    #birch = cluster.Birch(n_clusters=numClusters, threshold=0.05)
    #return birch.fit_predict(obs)

###############################################################################
def condense_subelements(labels, energies):
    '''Combine all energy points that have the same labels. If dual, combine all energy groups that have the same labels, and interpret energies as group boundaries. Output a dual mesh'''
    isDual = True
    if len(energies) == len(labels):
        isDual = False
    energies = energies.copy()
    small = 1E-5
    energies[energies == 0] = small
    if isDual:
        numBdrs = len(energies)
        toKeep = np.ones(numBdrs, dtype=bool)
        toKeep[1:-1] = (labels[1:] != labels[:-1])
        labelsOut = labels[toKeep[:-1]]
        energiesOut = energies[toKeep]
    else:
        numEnergies = len(energies)
        toKeep = np.ones(numEnergies, dtype=bool)
        toKeep = (labels[1:] != labels[:-1])
        groupsToKeep = np.sum(toKeep) + 1
        labelsOut = np.zeros(groupsToKeep, dtype=np.int)
        labelsOut[:-1] = labels[toKeep]
        labelsOut[-1] = labels[-1]
        energiesOut = np.zeros(groupsToKeep+1)
        energiesOut[0] = energies[0]
        toKeepIndices = np.where(toKeep)[0]
        energiesOut[1:-1] = np.sqrt(energies[toKeepIndices] * energies[toKeepIndices+1])
        energiesOut[-1] = energies[-1]
    return labelsOut, energiesOut

def average_observations(labels, energies, observations):
    '''Currently only works for primal meshes'''
    energies = energies.copy()
    small = 1E-5
    energies[energies == 0] = small
    #
    numGroups, numPoints = observations.shape
    #
    toKeep = np.ones(numGroups+1, dtype=bool)
    toKeep[1:-1] = (labels[1:] != labels[:-1])
    subelementBdrs = np.where(toKeep)[0]
    numSubelements = len(subelementBdrs) - 1
    observationsOut = np.zeros((numSubelements, numPoints))
    energiesOut = np.zeros(numSubelements)
    labelsOut = np.zeros(numSubelements, dtype=int)
    for i in range(numSubelements):
        strt, end = subelementBdrs[i], subelementBdrs[i+1]
        observationsOut[i, :] = np.mean(observations[strt:end, :], axis=0)
        energiesOut[i] = np.exp(np.mean(np.log(energies[strt:end])))
        labelsOut[i] = labels[strt]
    return labelsOut, energiesOut, observationsOut

###############################################################################
def find_minimum_neighbors_radius(observations):
    '''Find the minimum radius such that the number of connected components is 1.
    First, find a bounding interval. Then use bisection within the interval.'''
    #
    radiusNeighbors = 0.10
    bounds = [0, np.inf]
    foundBounds = [False, False]
    while not np.all(foundBounds):
        numRadComponents, radLabels = sparse.csgraph.connected_components(
            neighbors.radius_neighbors_graph(observations, radiusNeighbors))
        if numRadComponents > 1:
            bounds[0] = radiusNeighbors
            foundBounds[0] = True
            radiusNeighbors *= 2
        else:
            bounds[1] = radiusNeighbors
            foundBounds[1] = True
            radiusNeighbors /= 2
    #
    converged = False
    tol = 1E-2
    its = 0
    maxIts = 10
    while (its < maxIts and not converged):
        its += 1
        radiusNeighbors = 0.5 * (bounds[0] + bounds[1])
        numRadComponents, radLabels = sparse.csgraph.connected_components(
            neighbors.radius_neighbors_graph(observations, radiusNeighbors))
        if numRadComponents > 1:
            bounds[0] = radiusNeighbors
        else:
            bounds[1] = radiusNeighbors
        sz = (bounds[1] - bounds[0]) / radiusNeighbors
        if sz <= tol:
            converged = True
    radiusNeighbors = bounds[1]
    return radiusNeighbors

def get_index_neighbors(length):
    '''Return a neighbors graph that is tridiagonal, viz., the neighbors are determined by index'''
    arr = np.ones(length)
    return sparse.spdiags([arr, arr, arr], [-1, 0, 1], length, length)

def find_minimum_num_neighbors(observations):
    '''Find the minimum number of neighbors such that the number of connected components is 1.
    First, find a bounding interval. Then use bisection within the interval.'''
    #
    numNeighbors = 10
    numPoints = observations.shape[0]
    bounds = [0, numPoints]
    foundBounds = [False, False]
    while not np.all(foundBounds):
        numKnnComponents, knnLabels = sparse.csgraph.connected_components(
            neighbors.kneighbors_graph(observations, numNeighbors, include_self=True))
        if numKnnComponents > 1:
            bounds[0] = numNeighbors
            foundBounds[0] = True
            numNeighbors *= 2
        else:
            bounds[1] = numNeighbors
            foundBounds[1] = True
            numNeighbors /= 2
        if numNeighbors == 0:
            bounds[0] = 0
            foundBounds[0] = True
    #
    converged = False
    minSz = 1
    its = 0
    maxIts = 10
    while (its < maxIts and not converged):
        its += 1
        numNeighbors = int(np.round(0.5 * (bounds[0] + bounds[1])))
        numKnnComponents, knnLabels = sparse.csgraph.connected_components(
            neighbors.kneighbors_graph(observations, numNeighbors, include_self=True))
        if numKnnComponents > 1:
            bounds[0] = numNeighbors
        else:
            bounds[1] = numNeighbors
        sz = bounds[1] - bounds[0]
        if sz <= minSz:
            converged = True
    numNeighbors = bounds[1]
    if not converged:
        print "Not converged!"
    return numNeighbors


###############################################################################
def load_txt(filePath, skiprows=1, usecols=[1]):
    return np.loadtxt(filePath, skiprows=skiprows, usecols=usecols)

###############################################################################
def define_input_parser():
    import argparse
    #
    parser = argparse.ArgumentParser(description='Minimizer of observation errors and creator of generalized energy meshes.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    defaults = define_defaults()
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=defaults['verbosity'], choices=[0,1,2,3,4], type=int)
    parser.add_argument('-e', '--elements', dest='numelements', help="Number of total energy elements to use for the energy mesh. Ignored if the 'numcoarsegroups' option is used.", type=int, default=defaults['numelements'])
    parser.add_argument('-T', '--totalnumelements', dest='numelementsistotal', help="If specified, the number of energy elements is taken to be the total number of elements. If not specified, the number of energy elements is taken to be the number of elements in the resolved resonance region. Ignored if the 'numcoarsegroups' option is used.", action='store_true', default=False)
    parser.add_argument('-a', '--apportionopt', help="Specify how to assign the number of energy elements per coarse group, if not explicitly specified using 'listnumelements'. 'equal' means use an equal number of elements per coarse group. 'var', 'max', 'birch', and 'L1' are four ways that assign elements in proportion to the relative variance within a coarse group. 'L1' is not normalized to the number of fine points per coarse group and is more useful for 'amg' and 'tmg' work options. 'birch' uses the number of Birch clusters, and is approximate.", choices=['equal', 'var', 'max', 'birch', 'L1'], default=defaults['apportionopt'])
    parser.add_argument('-S', '--sigt', help='If specified, use Sigma_t as the other indicator. If not specified, use phi calculated with an escape cross sections as the other indicator. Always use the infinite-medium flux as the first indicator. Fluxes are normalized in shape to be more constant: the Maxwellian-1/E-fission source shape is divided out.', action='store_true', default=False)
    parser.add_argument('-w', '--workopt', help='What to do. har means do hierarchical agglomerative clustering (with restricted connectivity based on a set number of nearest neighbors). mg means do even spacing in lethargy (or as close as possible given the input energy mesh). amg means minimize squared error within each group. tmg means evenly divide in index (topology) space from input grid', choices=['amg','tmg','mg','har'], default=defaults['workopt'])
    parser.add_argument('-c', '--coarsebdrs', help='The resolved resonance range and how it is to be split into coarse groups (one clustering calculation per coarse group).', type=float, nargs='+', default=defaults['coarsebdrs'])
    parser.add_argument('-n', '--numcoarsegroups', help="The number of coarse groups to be used. If nonzero, overwrites the internal members of 'coarsebdrs'", type=int, default=defaults['numcoarsegroups'])
    parser.add_argument('-l', '--listnumelements', help='Number of elements to be put in each coarse boundary. Number of arguments should be one less than the number of coarse boundaries. Takes priority over "elements" if set', type=int, nargs='+', default=defaults['listnumelements'])
    parser.add_argument('-r', '--resolution', help='Resolution to use for the pointwise flux calculations', type=int, choices=range(11), default=defaults['resolution'])
    parser.add_argument('-m', '--materialopt', help="Unless 'manual' is used, specifies a set of materials to use. If 'manual' is used, give a space-separated list of material names in 'listmaterials'.", choices=['4','5','c5g7', 'graphite', 'iron', 'kpin', 'kenrichedpin', 'kcladpin', 'kpin2d', 'kenrichedpin2d', 'kmoxpin2d', 'kmoxenrichedpin2d', 'trigafuel', 'ctrigafuel', 'ctrigafuel_0' 'trigamore', 'manual'], default=defaults['materialopt'])
    parser.add_argument('-i', '--indicatormaterials', dest='listmaterials', help="When manual 'materialopt' is used, specify the materials to use.", nargs='+', default=defaults['listmaterials'])
    parser.add_argument('-I', '--importances', dest='listimportances', help="When manual 'materialopt' is used, specify the weightings (importances) to use when clustering.", nargs='+', type=int, default=[])
    parser.add_argument('-p', '--plotopt', help='Which observations to plot', choices=['none', 'first', 'last', 'firstlast', 'half', 'all', 'sum'], default=defaults['plotopt'])
    parser.add_argument('-s', '--shownumbers', help='If true, show element numbers on the plots', type=int, default=defaults['shownumbers'])
    parser.add_argument('-E', '--energypenalty', help='The energy variable is added to the observations to encourage contiguity for high numbers of elements. A value of 0 will not penalize in energy at all. A very large value will yield equal-lethargy-spaced MG', type=float, default=defaults['energypenalty'])
    parser.add_argument('-C', '--condensesubelements', help='If true, condense contiguous energy ranges before outputting', type=int, default=defaults['condensesubelements'])
    parser.add_argument('-d', '--dpi', help='Resolution to use for output plots', type=int, default=defaults['dpi'])
    return parser

###############################################################################
if __name__ == '__main__':
    parser = define_input_parser()
    inputDict = vars(parser.parse_args())
    if inputDict['verbosity'] > 1:
        print 'Summary of inputs:'
        print inputDict
    do_all(inputDict)

