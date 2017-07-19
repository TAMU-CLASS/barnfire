#! /usr/bin/env python

'''
Andrew Till
Summer 2014

Read in a GROUPR file from NJOY and store in a useful data structure.

See Readgroupr_readme.py for more information.

TODO: Find way to sneak in (3,457) and (3,458) into inputDict
TODO: Add depletion print option, or more generally, allow a list of mt numbers to be printed

TODO: Renormalize MT 2 with the thermal cross section
'''

#STDLIB
import os
import sys
import cPickle as pickle
#import cProfile as profile
from datetime import datetime
#TPL
import numpy as np
from scipy import sparse
#MINE

def execute_reader(inputDict):
    if 'finishedParsing' not in inputDict:
        raise ValueError('Invalid inputDict specified. Be sure to call finish_parsing(inputDict).')
    verbosity = inputDict['verbosity']
    #
    inDirr = inputDict['inputdir']
    gendfFilename = inputDict['inname']
    filePathIn = os.path.join(inDirr, gendfFilename)
    #
    outDirr = inputDict['outputdir']
    pdtXSName = inputDict['outname']
    pdtXSPath = os.path.join(outDirr, pdtXSName)
    pickleName = inputDict['picklename']
    if pickleName.lower() == 'none':
        picklePath = None
    else:
        picklePath = os.path.join(outDirr, pickleName)
    #
    mtsForMF3 = inputDict['mf3list']
    mtsForMF5 = inputDict['mf5list']
    mtsForMF6 = inputDict['mf6list']
    #thermalMTs = [241, 242]
    #thermalMTs = [221]
    #thermalMTs = []
    #thermalMTs = [239, 240]
    #thermalMTs = [222]
    thermalMTs = inputDict['thermallist']
    #
    desiredT = inputDict['temperature']
    desiredSig0 = inputDict['sig0']
    if inputDict['sig0Vec'] is not None:
        desiredSig0 = inputDict['sig0Vec']
    #
    flux = inputDict['flux']
    energyMesh = inputDict['energyMesh']
    doCondensation = False
    # TODO: Generalize in a smart way if flux is not specified
    if flux is not None and energyMesh is not None:
        doCondensation = True
    #
    format = inputDict['format']
    whichXS = inputDict['printopt']
    #
    workOpt = inputDict['workopt']
    outputDict = {}
    if workOpt == 'gendf':
        xsData = read_xs(filePathIn, mtsForMF3, mtsForMF5, mtsForMF6, verbosity)
        if picklePath is not None:
            pickle_xs(xsData, picklePath)
        if doCondensation:
            xsData = condense_xs(xsData, energyMesh, flux, verbosity)
        interpolate_T_sig0_xs(xsData, desiredT, desiredSig0, outputDict, verbosity)
        if whichXS != 'none':
            combine_transfer_matrix(xsData, thermalMTs)
            convert_to_scipy_sparse_scat_matrix(xsData, format)
            write_pdt_xs(pdtXSPath, xsData, desiredT, format, whichXS)
            print_xs_summary(filePathIn, xsData, verbosity)
    elif workOpt == 'pickle':
        xsData = unpickle_xs(picklePath)
        if doCondensation:
            xsData = condense_xs(xsData, energyMesh, flux, verbosity)
        interpolate_T_sig0_xs(xsData, desiredT, desiredSig0, outputDict, verbosity)
        if whichXS != 'none':
            combine_transfer_matrix(xsData, thermalMTs)
            convert_to_scipy_sparse_scat_matrix(xsData, format)
            write_pdt_xs(pdtXSPath, xsData, desiredT, format, whichXS)
            print_xs_summary(filePathIn, xsData, verbosity)
    elif workOpt == 'rxn':
        get_reaction_list(filePathIn, verbosity)
    elif workOpt == 'scat':
        get_scattering_sizes(filePathIn, verbosity)
    elif workOpt == 'pen':
        get_pointwise_xs(filePathIn, mtsForMF3, desiredT, outputDict, verbosity)
    elif workOpt == 'pot':
        get_pot_scat_xs(filePathIn, outputDict, verbosity)
    elif workOpt == 'decay':
        get_decay_parameters(filePathIn, outputDict, verbosity)
    elif workOpt == 'meta':
        get_metastable_branching_ratios(filePathIn, outputDict, verbosity)
    elif workOpt == 'fp':
        get_fission_product_yields(filePathIn, outputDict, verbosity)
    elif workOpt == 'qvalue':
        get_reaction_Q_values(filePathIn, mtsForMF3, outputDict, verbosity)
    return outputDict

####################################################################################
def read_xs(filePath, desiredMTsForMF3, desiredMTsForMF5, desiredMTsForMF6, verbosity):
    '''Read desired xs and transfer matrices from GROUPR file with handle fid'''
    data = {}
    # Loop through entire file once to figure out temperatures and reactions
    mfmtsAvail, numGroups, numLegMoments, thermList, numSig0 = get_reaction_list(filePath, verbosity)
    # Keep intersection of what you desire and what is available
    mfmtsKeep = determine_rxns_to_keep(desiredMTsForMF3, desiredMTsForMF5, desiredMTsForMF6, mfmtsAvail, verbosity)
    # Initialize data structure, which houses everything.
    populate_data_dict(data, mfmtsKeep, numGroups, numLegMoments, thermList, numSig0)
    # Go through file and store reactions in mfmtsKeep
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        line = get_next_line(fid)
        while line != '':
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if not (is_continue(mf, mt, lineNum)):
                # Read header
                xsDict = {}
                read_xs_header(line, xsDict)
            if (mf == 1) and (mt == 451):
                # Read new temperature
                eBdrDict = {}
                read_group_bdr_header(line, eBdrDict)
                read_group_bdr_body(fid, eBdrDict)
                merge_data_ebdrdict(data, eBdrDict)
                line = get_next_line(fid)
                line = get_next_line(fid)
            elif (mf, mt) in mfmtsKeep:
                if mf == 3:
                    # Read new xs
                    read_xs_vector_body(fid, xsDict, verbosity)
                    merge_data_xs_vector(data, xsDict)
                    line = get_next_line(fid)
                elif (mf, mt) == (6, 18):
                    # Read prompt fission transfer matrices, which are special
                    read_prompt_fission_matrix_body(fid, xsDict, verbosity)
                    merge_data_prompt_fission_matrix(data, xsDict)
                    line = get_next_line(fid)
                elif (mf, mt) == (5, 455):
                    # Read delayed fission chi, nu, lambda (decay constant)
                    read_delayed_fission_spectrum(fid, xsDict, verbosity)
                    merge_data_delayed_fission_xs(data, xsDict)
                    line = get_next_line(fid)
                    key = 2055
                elif mf == 6:
                    # Read non-fission matrix
                    read_xs_matrix_body(fid, xsDict, verbosity)
                    #if mt == 2:
                    #    profile.runctx('read_xs_matrix_body(fid, xsDict, verbosity)', globals(), locals())
                    #else:
                    #    read_xs_matrix_body(fid, xsDict, verbosity)
                    merge_data_xs_matrix(data, xsDict)
                    line = get_next_line(fid)
                    key = get_pdt_mt(mf, mt)
                else:
                    line = seek_to_next_useful_entry(fid, mf, mt, verbosity)
            else:
                # Skip reaction
                line = seek_to_next_useful_entry(fid, mf, mt, verbosity)
    unify_scattering_sparsity_patterns(data)
    return data

####################################################################################
def populate_data_dict(data, mfmtsSet, numGroups, numLegMoments, thermList, numSig0):
    '''Populate the data dictionary with sizes and reaction lists. Allocate what you know the size of (xs, flux, csc param; not sparse xfer matrix, fission matrix)'''
    numTherm = len(thermList)
    data['mfmts'] = mfmtsSet
    data['thermList'] = sorted(thermList)
    data['numLegMoments'] = numLegMoments
    data['numTherm'] = numTherm
    data['numSig0'] = numSig0
    data['rxn'] = {}
    data['flux'] = np.zeros((numTherm, numSig0, numGroups))
    for (mf,mt) in mfmtsSet:
        actualNumTherm = lookup_num_therm(numTherm, mf, mt)
        actualNumSig0 = lookup_num_sig0(numSig0, mf, mt)
        data['rxn'][(mf,mt)] = {}
        rxn = data['rxn'][(mf,mt)]
        rxn['numSig0'] = actualNumSig0
        rxn['numTherm'] = actualNumTherm
        # Cannot allocate XS for MF 6 here because we don't know the sparsity patterns yet
        if mf == 3:
            rxn['xs'] = np.zeros((actualNumTherm, actualNumSig0, numGroups))
            rxn['numLegMoments'] = 1
        elif mf == 5 and mt == 455:
            rxn['numLegMoments'] = 6  # Assuming 6 delayed neutron groups
            rxn['numDNGs'] = 6
        elif mf == 6 and mt == 18:
            rxn['numLegMoments'] = 1
        elif mf == 6:
            #'List' is indexed by temperature
            # Cannot set numLegMoments because some transfer matrices have fewer
            rxn['xsList'] = []
            rxn['rowStartMat'] = np.zeros((actualNumTherm, numGroups), dtype=np.int)
            rxn['colSizeMat'] = np.zeros((actualNumTherm, numGroups), dtype=np.int)
            rxn['numLegMoments'] = numLegMoments

def merge_data_ebdrdict(data, eBdrDict):
    '''Merge eBdrDict into data. Overwrite group boundaries, sigma0 values. Make energy run high to low.'''
    data['Z'] = eBdrDict['Z']
    data['A'] = eBdrDict['A']
    data['weight'] = eBdrDict['weight']
    data['numSig0'] = eBdrDict['numSig0']
    data['numGroups'] = eBdrDict['numNeutronGroups']
    #The following data is copied, not aliased
    data['sig0List'] = sorted(eBdrDict['sigma0List'], reverse=True)
    data['groupBdrs'] = np.array(eBdrDict['neutronEnergyBdrs'][::-1])

def merge_data_xs_vector(data, xsDict):
    '''Merge xsDict into data. Overwrite flux, add to xs of (mf,mt). Make energy run high to low.'''
    mf,mt = xsDict['mf'], xsDict['mt']
    rxnDict = data['rxn'][(mf,mt)]
    #
    flux = data['flux']
    xs = rxnDict['xs']
    thermIndex = 0
    if rxnDict['numTherm'] != 1:
        thermIndex = data['thermList'].index(xsDict['temperature'])
    #The data is copied, not aliased
    xs[thermIndex, :, :] = xsDict['xs'][:, ::-1]
    if mt == 1:
        flux[thermIndex, :, :] = xsDict['flux'][:, ::-1]
    #xs and flux are indexed by xs[thermIndex, sig0Index, groupIndex]

def unify_scattering_sparsity_patterns(data):
    '''For each reaction, take a union over temperature of the sparsity patterns, because they may be different. Copy the xs into this global sparse format. This will allow for easy interpolation between temperatures later.'''
    #
    numGroups = data['numGroups']
    identityForMin = numGroups
    mfmtList = [(mf,mt) for (mf,mt) in data['mfmts'] if (mf == 6 and mt != 18)]
    for (mf,mt) in mfmtList:
        rxnDict = data['rxn'][(mf,mt)]
        numTherm = rxnDict['numTherm']
        numSig0 = rxnDict['numSig0']
        numLegMoments = rxnDict['numLegMoments']
        #
        rowStartMat = rxnDict['rowStartMat']
        rowStartMat[rowStartMat == -1] = identityForMin
        minRowStart = np.min(rowStartMat, axis=0)
        rowEndMat = rowStartMat + rxnDict['colSizeMat']
        maxRowEnd = np.max(rowEndMat, axis=0)
        maxColSize = maxRowEnd - minRowStart

        minRowStart, maxColSize = remove_sparse_holes(minRowStart, maxColSize)
        minRowStart[minRowStart == identityForMin] = -1

        unionIndexPtr = np.zeros(numGroups + 1, dtype=np.int)
        unionIndexPtr[1:] = np.cumsum(maxColSize)
        xsSize = unionIndexPtr[-1]
        rxnDict['xs'] = np.zeros((numTherm, numSig0, numLegMoments, xsSize))
        xs = rxnDict['xs']
        #
        for thermIndex in range(numTherm):
            localXS = rxnDict['xsList'][thermIndex]
            localColSize = rxnDict['colSizeMat'][thermIndex, :]
            localIndexPtr = np.zeros(numGroups+1, dtype=np.int)
            localIndexPtr[1:] = np.cumsum(localColSize)
            localRowStart = rxnDict['rowStartMat'][thermIndex, :]
            for fromGroup in range(numGroups):
                localStrt, localEnd = localIndexPtr[fromGroup], localIndexPtr[fromGroup + 1]
                unionStrt = unionIndexPtr[fromGroup] + (localRowStart[fromGroup] - minRowStart[fromGroup])
                unionEnd = unionStrt + localColSize[fromGroup]
                xs[thermIndex, :, :, unionStrt:unionEnd] = localXS[:, :, localStrt:localEnd]
        rxnDict['rowStart'] = minRowStart
        rxnDict['indexPtr'] = unionIndexPtr
        del(rxnDict['xsList'])
        del(rxnDict['rowStartMat'])
        del(rxnDict['colSizeMat'])
    #Scattering from group g may be found via this indexing:
    #xs[thermIndex, sig0Index, legIndex, indexPtr[g]:indexPtr[g+1]]
    #The groups into which g scatters are:
    #rowStart[g]:(rowStart[g]+indexPtr[g+1]-indexPtr[g]-1)

def merge_data_xs_matrix(data, xsDict):
    '''Merge xsDict into data. Append new temperatures to the scattering matrix. Make energy run high to low. Because reactions at different temperatures may have different sparsity patterns (e.g., thermal reactions), we store the scattering kernel as a list of CSC matrices, where the list is indexed by temperature.'''
    #
    mf,mt = xsDict['mf'], xsDict['mt']
    rxnDict = data['rxn'][(mf,mt)]
    #
    numGroups = data['numGroups']
    #
    # The number of Legendre moments may be truncated due to insufficient data in evaluations
    actualNumLegMoments = xsDict['numLegMoments']
    numLegMoments = rxnDict['numLegMoments']
    #
    numSig0 = rxnDict['numSig0']
    numTherm = rxnDict['numTherm']
    #
    high2lowIndexPtr = np.zeros(numGroups+1, dtype=np.int)
    high2lowIndexPtr[1:] = np.cumsum(xsDict['colSize'][::-1])
    #
    # Copy from xs to rxnDict and flip ordering of sparse matrix
    xsSize = high2lowIndexPtr[-1]
    rxnDict['xsList'].append(np.zeros((numSig0, numLegMoments, xsSize)))
    xs = rxnDict['xsList'][-1]
    numRows = xsDict['lastRow'] - xsDict['firstRow'] + 1
    firstRow = xsDict['firstRow']
    for g in range(numRows):
        high2LowGroupFrom = numGroups - (g + firstRow) - 1
        strt = high2lowIndexPtr[high2LowGroupFrom]
        end = high2lowIndexPtr[high2LowGroupFrom + 1]
        for legIndex in range(actualNumLegMoments):
                for sig0Index in range(numSig0):
                    # The data is copied, not aliased
                    # Be careful: the index ordering is different on the left and right
                    xs[sig0Index, legIndex, strt:end] = xsDict['xs'][g][legIndex, sig0Index, :][::-1]
    # Use the full number of Legendre moments for all transfer rxn, even if there is no data.
    xs[:, actualNumLegMoments:, :] = 0.

    # Each temperature in general has a different sparsity pattern, so store CSC indexing for each
    thermIndex = 0
    if numTherm !=1:
        thermIndex = data['thermList'].index(xsDict['temperature'])
    # Copy rowStart to rxnDict and flip ordering
    high2lowRowStart = rxnDict['rowStartMat'][thermIndex, :]
    # The starting group from high to low is the ending group from low to high
    high2lowRowStart[:] = (xsDict['colSize'] + xsDict['rowStart'] - 1)[::-1]
    # Flipping order changes the index of a group and we refer directly to indexes with this array
    high2lowRowStart[:] = numGroups - high2lowRowStart[:] - 1
    unusedRows = (xsDict['rowStart'] == -1)
    high2lowRowStart[unusedRows[::-1]] = -1
    # Copy colSize to rxnDict and flip ordering
    rxnDict['colSizeMat'][thermIndex, :] = xsDict['colSize'][::-1]

def merge_data_prompt_fission_matrix(data, xsDict):
    '''Merge xsDict into data. Overwrite fission spectrum, fission matrix, and group where matrix starts'''
    mf,mt = xsDict['mf'], xsDict['mt']
    rxnDict = data['rxn'][(mf,mt)]
    numGroups = data['numGroups']
    #The data is copied, not aliased
    rxnDict['flux'] = xsDict['flux'][::-1].copy()
    rxnDict['lowEnergySpectrum'] = xsDict['lowEnergySpectrum'][::-1].copy()
    rxnDict['lowEnergyProd'] = xsDict['lowEnergyProd'][::-1].copy()
    rxnDict['highEnergyMatrix'] = xsDict['highEnergyMatrix'][::-1,::-1].copy()
    rxnDict['highestHighEnergyGroup'] = numGroups - xsDict['lowestHighEnergyGroup'] - 1
    rxnDict['promptChi'] = xsDict['promptChi'][::-1].copy()
    rxnDict['promptProd'] = xsDict['promptProd'][::-1].copy()
    rxnDict['FissionMatrix'] = xsDict['FissionMatrix'][::-1, ::-1].copy()
    #lowEnergySpectrum is indexed by (groupTo)
    #highEnergyMatrix is indexed by (groupFrom, groupTo)
    #highEnergyMatrix applies to groups 0 through highestHighEnergyGroup (inclusive)
    #FissionMatrix is indexed by (groupFrom, groupTo)

def merge_data_delayed_fission_xs(data, xsDict):
    '''Merge xsDict into data. Overwrite delayed fission chi, nu, and decay constant'''
    mf,mt = xsDict['mf'], xsDict['mt']
    rxnDict = data['rxn'][(mf,mt)]
    numGroups = data['numGroups']
    numLegMoments = data['numLegMoments']
    #The data is copied, not aliased
    rxnDict['decayConst'] = xsDict['decayConst']
    rxnDict['delayedChi'] = xsDict['delayedChi'][:,::-1].copy()
    rxnDict['numDNGs'] = xsDict['numDNGs']  # For mt 455 numLegMoments is numDelayedNeutronGroups

####################################################################################
def read_group_bdr_header(line, paramDict):
    '''Reads in a CONT record (see NJOY manual, page 2801). Adds Z, A, weight, numSig0, titleLength to paramDict.'''
    ZZAAA = '{0:5g}'.format(int(get_float(line, 1)))
    Z = int(ZZAAA[0:2])
    A = int(ZZAAA[2:5])
    # weight is relative to a neutron
    weight = get_float(line, 2)
    numSig0 = get_int(line, 4)
    titleLength = get_int(line, 6)
    headerDict = {'Z': Z, 'A': A, 'weight': weight, 'numSig0': numSig0, 'titleLength': titleLength}
    paramDict.update(headerDict)

def read_group_bdr_body(fid, paramDict):
    '''Reads in a LIST record. Advances fid. Adds temperature, title, sigma0List, neutronEnergyList, and gammaEnergyList to paramDict'''
    titleLength = paramDict['titleLength']
    numSig0 = paramDict['numSig0']
    line = get_next_line(fid)
    #
    temperature = get_float(line, 1)
    numNeutronGroups = get_int(line, 3)
    numGammaGroups = get_int(line, 4)
    numEntries = get_int(line, 5)
    pos = 6
    titleList, pos, line = get_list(pos, line, fid, titleLength, get_str)
    sigma0List, pos, line = get_list(pos, line, fid, numSig0, get_float)
    neutronEnergyBdrs, pos, line = get_list(pos, line, fid, numNeutronGroups+1, get_float)
    gammaEnergyBdrs, pos, line = get_list(pos, line, fid, numGammaGroups+1, get_float)
    title = ''.join(titleList)
    bodyDict = {'temperature': temperature, 'numNeutronGroups': numNeutronGroups, 'numGammaGroups': numGammaGroups, 'title': title, 'sigma0List': sigma0List, 'neutronEnergyBdrs': neutronEnergyBdrs, 'gammaEnergyBdrs': gammaEnergyBdrs}
    paramDict.update(bodyDict)

####################################################################################
def read_xs_header(line, xsDict):
    '''Reads in a CONT record. Adds Z, A, weight, numLegMoments, numSig0, numGroups, mf, mt, lineNum to xsDict.'''
    ZZAAA = '{0:5g}'.format(int(get_float(line, 1)))
    Z = int(ZZAAA[0:2])
    A = int(ZZAAA[2:5])
    weight = get_float(line, 2)
    numLegMoments = get_int(line, 3)
    numSig0 = get_int(line, 4)
    numGroups = get_int(line, 6)
    mf, mt, lineNum = get_mf_mt_lineNum(line)
    headerDict = {'Z': Z, 'A': A, 'weight': weight, 'numLegMoments': numLegMoments, 'numSig0': numSig0, 'numGroups': numGroups, 'mf': mf, 'mt': mt, 'lineNum': lineNum}
    xsDict.update(headerDict)

def read_xs_vector_body(fid, xsDict, verbosity=0):
    '''Reads in a LIST record. Advances fid. Writes flux and xs to xsDict.'''
    numLegMoments = xsDict['numLegMoments']
    numSig0 = xsDict['numSig0']
    numGroups = xsDict['numGroups']
    mf = xsDict['mf']
    mt = xsDict['mt']
    pdtMT = get_pdt_mt(mf, mt)
    if verbosity > 1:
        print 'Reading  MF {0}, MT {1:3g} and storing in PDT-MT {2:4g}'.format(mf, mt, pdtMT)
    line = get_next_line(fid)
    temperature = get_float(line, 1)
    groupIndex = 0
    flux = np.zeros((numSig0, numGroups))
    xs = np.zeros((numSig0, numGroups))
    # Read cross section from file and store in xsDict
    while groupIndex < numGroups:
        numSecondaryPos = get_int(line, 3)
        indexLowestGroup = get_int(line, 4)
        numEntries = get_int(line, 5)
        groupIndex = get_int(line, 6)
        pos = 0
        #xsForOneGroup, pos, line = get_list(pos, line, fid, numEntries, get_float)
        xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
        strt, end = 0, numSig0 * numLegMoments
        flux[:, groupIndex-1] = xsForOneGroup[strt:end:numLegMoments]
        strt, end = numSig0 * numLegMoments, 2 * numSig0 * numLegMoments
        xs[:, groupIndex-1] = xsForOneGroup[strt:end:numLegMoments]
        line = get_next_line(fid)
    bodyDict = {'flux': flux, 'xs': xs, 'temperature': temperature}
    xsDict.update(bodyDict)

def read_xs_matrix_body(fid, xsDict, verbosity=0):
    '''Reads in a LIST record. Advances fid. Writes xs, colSize, rowStart, firstRow, and LastRow to xsDict.'''
    numLegMoments = xsDict['numLegMoments']
    numSig0 = xsDict['numSig0']
    numGroups = xsDict['numGroups']
    mf = xsDict['mf']
    mt = xsDict['mt']
    pdtMT = get_pdt_mt(mf, mt)
    if verbosity > 1:
        print 'Reading  MF {0}, MT {1:3g} and storing in PDT-MT {2:4g}'.format(mf, mt, pdtMT)
    line = get_next_line(fid)
    temperature = get_float(line, 1)
    groupIndex = 0
    xs = []
    colSize = np.zeros(numGroups, dtype=np.int)
    rowStart = -1 * np.ones(numGroups, dtype=np.int)
    firstGroupFrom = get_int(line, 6) - 1
    lastGroupFrom = firstGroupFrom - 1
    # Read transfer matrix from file and store in CSC-like format in xsDict
    while groupIndex < numGroups:
        numSecondaryPos = get_int(line, 3)
        indexLowestGroup = get_int(line, 4)
        numEntries = get_int(line, 5)
        groupIndex = get_int(line, 6)
        pos = 0
        #xsForOneGroup, pos, line = get_list(pos, line, fid, numEntries, get_float)
        xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
        # Once the thermal data is done, GENDF skips to a null record for the highest group, which we skip
        if mt in range(221, 250): # if mt is thermal process
            if (groupIndex-1) != (lastGroupFrom+1):
                line = get_next_line(fid)
                break
        lastGroupFrom += 1
        colSize[groupIndex-1] = numSecondaryPos - 1
        rowStart[groupIndex-1] = indexLowestGroup - 1
        # When reshaping, the fastest index should go first for Fortran ordering
        xsRearrange = np.reshape(np.asarray(xsForOneGroup), (numLegMoments, numSig0, numSecondaryPos), order='F')
        xs.append(xsRearrange[:,:,1:])
        if np.all(xsRearrange[:,:,1:] == 0.0) and groupIndex == numGroups:
            # If thermal reaction is null, prevent combine_transfer_matrix from overwriting xfer matrix with all zeros
            rowStart[-1] = -1
            colSize[-1] = 0
        line = get_next_line(fid)
    bodyDict = {'xs': xs, 'colSize': colSize, 'rowStart': rowStart, 'firstRow': firstGroupFrom, 'lastRow': lastGroupFrom, 'temperature': temperature}
    xsDict.update(bodyDict)


def read_prompt_fission_matrix_body(fid, xsDict, verbosity=0):
    '''Reads in a LIST record. Advances fid. Writes xs to xsDict. No temperature or sigma0 dependence is expected for the fission spectrum.'''
    numLegMoments = xsDict['numLegMoments']
    numSig0 = xsDict['numSig0']
    numGroups = xsDict['numGroups']
    mf = xsDict['mf']
    mt = xsDict['mt']
    pdtMT = 1018
    if verbosity > 1:
        print 'Reading  MF {0}, MT {1:3g} and storing in PDT-MT {2:4g}'.format(mf, mt, pdtMT)
    flux = np.zeros(numGroups)
    lowEnergySpectrum = np.zeros(numGroups)
    line = get_next_line(fid)
    temperature = get_float(line, 1)
    numSecondaryPos = get_int(line, 3)
    indexLowestGroup = get_int(line, 4)
    numEntries = get_int(line, 5)
    groupIndex = get_int(line, 6)
    if groupIndex == 0:
        # Read piece 1: Low-energy fission spectrum (may not exist)
        pos = 0
        #xsForOneGroup, pos, line = get_list(pos, line, fid, numEntries, get_float)
        xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
        strt, end = (indexLowestGroup - 1), (indexLowestGroup - 1) + numSecondaryPos
        lowEnergySpectrum[strt:end] = xsForOneGroup
        # Read piece 2: Low-energy neutron production cross section
        line = get_next_line(fid)  #read CONT card for first prod section
        indexLowestGroup = get_int(line, 4)
        numEntries = get_int(line, 5)
        indexIncidentGroup = get_int(line, 6)
        lowEnergyProd = np.zeros(0)
        while not(indexLowestGroup) and not(indexIncidentGroup == numGroups):
            pos = 0
            #read data in this prod section
            xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
            flux[indexIncidentGroup - 1] = xsForOneGroup[0]
            lowEnergyProd = np.append(lowEnergyProd, xsForOneGroup[1])  #2nd entry is nuSig_f
            line = get_next_line(fid)  #read CONT card for consecutive prod section
            indexLowestGroup = get_int(line, 4)
            numEntries = get_int(line, 5)
            indexIncidentGroup = get_int(line, 6)
        hasHighEnergyPiece = True
        if not(indexLowestGroup) and (indexIncidentGroup == numGroups):
            hasHighEnergyPiece = False
            pos = 0
            #read data in this prod section
            xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
            lowEnergyProd = np.append(lowEnergyProd, xsForOneGroup[1])  #2nd entry is nuSig_f
            line = get_next_line(fid)  # skip blank line at the end of this section
    else:
        hasHighEnergyPiece = True
    # Read piece 3: High-energy transfer matrix (may not exist)
    if hasHighEnergyPiece:
        groupIndex = get_int(line, 6)
        lowestHighEnergyGroup = groupIndex - 1
        numHighEnergyGroups = numGroups - groupIndex + 1
        highEnergyMatrix = np.zeros((numHighEnergyGroups, numGroups))
        groupIndex = 0
        while groupIndex < numGroups:
            numSecondaryPos = get_int(line, 3)
            indexLowestGroup = get_int(line, 4)
            numEntries = get_int(line, 5)
            groupIndex = get_int(line, 6)
            pos = 0
            #xsForOneGroup, pos, line = get_list(pos, line, fid, numEntries, get_float)
            xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
            strt, end = (indexLowestGroup - 1), (indexLowestGroup - 1) + (numSecondaryPos - 1)
            rowIndex = groupIndex - lowestHighEnergyGroup - 1
            flux[groupIndex - 1] = xsForOneGroup[0]
            highEnergyMatrix[rowIndex, strt:end] = xsForOneGroup[1:]
            line = get_next_line(fid)
    else:
        lowestHighEnergyGroup = numGroups
        highEnergyMatrix = np.zeros((0, numGroups))
    assert len(lowEnergyProd) + numHighEnergyGroups == numGroups
    assert len(lowEnergyProd) == lowestHighEnergyGroup

    # Reconstruct full prompt fission matrix (F)
    FissionMatrix = np.zeros((numGroups, numGroups))
    FissionMatrix[0:lowestHighEnergyGroup, :] = np.outer(lowEnergyProd, lowEnergySpectrum)
    FissionMatrix[lowestHighEnergyGroup:numGroups, :] = highEnergyMatrix
    # Compute prompt fission source = F'*phi
    F_phi = np.dot(flux, FissionMatrix)
    # prompt Chi is normalized promt fission source
    promptChi = F_phi/np.sum(F_phi)
    # prompt nu*sig_f (Prod) is row sum of F
    promptProd = np.sum(FissionMatrix, 1)

    bodyDict = {'flux': flux, 'lowEnergySpectrum': lowEnergySpectrum, 'lowEnergyProd': lowEnergyProd, \
                'lowestHighEnergyGroup': lowestHighEnergyGroup, 'highEnergyMatrix': highEnergyMatrix, \
                'FissionMatrix': FissionMatrix, 'promptChi': promptChi, 'promptProd': promptProd, \
                'temperature': temperature}

    xsDict.update(bodyDict)


def read_delayed_fission_spectrum(fid, xsDict, verbosity=0):
    '''Reads in a LIST record. Advances fid. Writes xs to xsDict. No temperature or sigma0 dependence is expected for the fission spectrum.'''
    numDNGs = xsDict['numLegMoments']  # number of delayed neutron groups
    numSig0 = xsDict['numSig0']
    assert numSig0 == 1               # delayed Chi does not have sigma0 depencence
    numGroups = xsDict['numGroups']
    mf = xsDict['mf']
    mt = xsDict['mt']
    pdtMT = 2055
    if verbosity > 1:
        print 'Reading  MF {0}, MT {1:3g} and storing in PDT-MT {2:4g}'.format(mf, mt, pdtMT)
    nu = np.zeros(numGroups)
    line = get_next_line(fid)
    temperature = get_float(line, 1)
    numSecondaryPos = get_int(line, 3)
    numSecondaryGroups = numSecondaryPos - 1  # First position is for decay constant
    indexLowestGroup = get_int(line, 4) - 1
    numEntries = get_int(line, 5)
    groupIndex = get_int(line, 6)
    decayConst = np.zeros(numDNGs)
    delayedChi = np.zeros((numDNGs, numGroups))
    pos = 0
    xsForOneGroup, pos, line = get_float_list(pos, line, fid, numEntries)
    xsStrt, xsEnd =  0, numDNGs
    decayConst =  xsForOneGroup[xsStrt:xsEnd]
    xsStrt, xsEnd = numDNGs, numDNGs * (numSecondaryGroups + 1)
    chiStrt, chiEnd = indexLowestGroup, indexLowestGroup + numSecondaryGroups
    delayedChi[:, chiStrt:chiEnd] = np.reshape(np.asarray(xsForOneGroup[xsStrt:xsEnd]), (numDNGs, numSecondaryGroups), order='F')
    bodyDict = {'decayConst': decayConst, 'delayedChi': delayedChi, 'numDNGs': numDNGs, 'temperature': temperature}
    line = get_next_line(fid)
    xsDict.update(bodyDict)


####################################################################################
def get_pointwise_xs(filePathIn, mtsForMF3, desiredT, outputDict, verbosity):
    '''Read in a pointwise xs from a pendf file. Reads in the temperature specified by desiredT.
    If this temperature is negative, read in the first temperature.
    Multiple temperatures are represented in a PENDF file by simply repeating the entire PENDF
    structure for each temperature.
    '''
    if 0 in mtsForMF3:
        mtsForMF3 = [1,2]
    mfmtList = [(3,mt) for mt in mtsForMF3]
    foundMFMTDict = {}
    for (mf,mt) in mfmtList:
        foundMFMTDict[(mf,mt)] = False
    foundT = False
    if desiredT < 0:
        foundT = True
    with open(filePathIn, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if (mf, mt) == (1, 451):
                for i in range(3):
                    line = get_next_line(fid)
                temperature = get_float(line, 1)
                if abs(temperature - desiredT) < (1E-8 * desiredT):
                    foundT = True
            elif (foundT and (mf, mt) in mfmtList):
                foundMFMTDict[(mf,mt)] = True
                for i in range(2):
                    line = get_next_line(fid)
                numPoints = get_int(line, 1)
                numEntries = 2 * numPoints
                outputDict['energy'] = np.zeros(numPoints)
                outputDict[(mf,mt)] = np.zeros(numPoints)
                data, pos, line = get_float_list(0, line, fid, numEntries)
                strt, end = 0, numEntries - 1
                outputDict['energy'][:] = data[strt:end:2]
                strt, end = 1, numEntries
                outputDict[(mf,mt)][:] = data[strt:end:2]
            if np.all(foundMFMTDict.values()):
                break
    if not foundT:
        print 'Did not find a temperature at {0} K.'.format(desiredT)
    elif verbosity:
        print 'Found a temperature at {0} K with {1} points.'.format(temperature, numPoints)

####################################################################################
def get_pot_scat_xs(filePath, outputDict, verbosity=False):
    '''Read a potential scattering XS from an ENDF file'''
    '''Resonance extents are given in mf 2, line 3 (first two numbers) for the resolved range.
    The unresolved range is harder to find, but always starts at the end of the resolved range,
    so if you know where the resolved range ends, you can look for the first number in each row
    to match. The second number of the matching row is the end of the unresolved resonance range.
    This location may not occur (e.g., H-1), so make sure you're still in mf 2, mt 151.
    '''
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        # Look through the file for each reaction type
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if (mf, mt) == (2,151):
                for i in range(3):
                    line = get_next_line(fid)
                scatLength = get_float(line, 2)
                potXS = 4 * np.pi * scatLength * scatLength
                outputDict['pot'] = potXS
                if verbosity:
                    print os.path.split(filePath)[-1], potXS
                return


def get_decay_parameters(filePath, outputDict, verbosity=0):
    '''Read the radioactive decay parameters in MF 8, MT 457, skipping the spectrum information)'''
    decayParticleDict = {0: None, 1: 'B-', 2: 'EC', 3: 'IT', 4: 'A', 5: 'N', 6: 'SF', 7: 'P'}
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if (mf,mt) == (8,457):
                line = get_next_line(fid)
                outputDict['halfLife'] = get_float(line, 1)
                #
                numDecayEnergyEntries = get_int(line, 5)
                numDecayLines = numDecayEnergyEntries // 6
                if (numDecayEnergyEntries % 6):
                    numDecayLines += 1
                for i in range(numDecayLines+1):
                    line = get_next_line(fid)
                #
                numDecayModes = get_int(line, 5) // 6
                decayModes = []
                for i in range(numDecayModes):
                    line = get_next_line(fid)
                    reactionType = get_float(line, 1)
                    metastableState = int(get_float(line, 2))
                    branchingRatio = get_float(line, 5)
                    #
                    decayParticles = []
                    while reactionType > 1E-6:
                        # divmod(a,b) returns (a//b, a%b)
                        decayParticle, remainder = divmod(reactionType, 1)
                        decayParticleStr = decayParticleDict[int(decayParticle)]
                        reactionType = 10 * remainder
                        if decayParticle:
                            decayParticles.append(decayParticleStr)
                    if branchingRatio > 1E-4:
                        decayModes.append(DecayMode(decayParticles, metastableState, branchingRatio))
                outputDict['decayModes'] = decayModes
                line, pos = seek_to_next_entry(fid, mf, mt)

class DecayMode():
    def __init__(self, decayParticles, metastableState, branchingRatio):
        self.decayParticles = decayParticles
        self.M = metastableState
        self.branchingRatio = branchingRatio

def get_metastable_branching_ratios(filePath, outputDict, verbosity=0):
    '''Read branching ratios (MF9) and cross sections (MF10) for production of metastable states'''
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if mf in (9,10):
                numMetaStates = get_int(line, 5)
                for i in range(numMetaStates):
                    line = get_next_line(fid)
                    ZZAAA = get_int(line, 3)
                    Z, A = divmod(ZZAAA, 1000)
                    # Metastable state number is specified in file 8, which we don't read.
                    # This assumes groundstate is M == 0 and first metastable state is M == 1.
                    M = get_int(line, 4)
                    if M > 1:
                        M = 1
                    if (mf,mt,M) in outputDict:
                        raise ValueError('({0},{1},{2}) already exists. Too many metastable states?'.format(mf,mt,M))
                    outputDict[(mf,mt,M)] = {}
                    xsData = outputDict[(mf,mt,M)]
                    xsData['Z'] = Z
                    xsData['A'] = A
                    #
                    numEnergyRanges = get_int(line, 5)
                    numEnergyPoints = get_int(line, 6)
                    energyInterpolation, pos, line = get_list(0, line, fid, 2*numEnergyRanges, get_int)
                    strt, end = 0, 2 * numEnergyRanges - 1
                    xsData['numEnergies'] = energyInterpolation[strt:end:2]
                    strt, end = 1, 2 * numEnergyRanges
                    #To interpolate in 'x-y', the interpolation scheme is:
                    #{1: 'histogram', 2: 'lin-lin', 3: 'log-lin', 4: 'lin-log', 5: 'log-log'}
                    xsData['interpolationSchemes'] = energyInterpolation[strt:end:2]
                    #
                    energyBranchingRatios, pos, line = get_float_list(0, line, fid, 2*numEnergyPoints)
                    strt, end = 0, 2 * numEnergyPoints - 1
                    xsData['energies'] = np.array(energyBranchingRatios[strt:end:2])
                    strt, end = 1, 2 * numEnergyPoints
                    if mf == 9:
                        xsData['branchingRatios'] = np.array(energyBranchingRatios[strt:end:2])
                    if mf == 10:
                        xsData['xs'] = np.array(energyBranchingRatios[strt:end:2])
                #line, pos = seek_to_next_entry(fid, mf, mt)

def get_fission_product_yields(filePath, outputDict, verbosity=0):
    '''Read independent (MT454) and cumulative (MT459) fission product yields (MF8) for spontaneous fission or neutron-induced fission'''
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if (mf,mt) in [(8,454), (8,459)]:
                if mt == 454:
                    # MT 454 is independent yield
                    outputDict['indepYield'] = {}
                    yieldData = outputDict['indepYield']
                else:
                    # MT 459 is cumulative yield
                    outputDict['cumulYield'] = {}
                    yieldData = outputDict['cumulYield']
                # Fission product yield data is given on an energy grid
                numEnergies = get_int(line, 3) #??
                yieldData['energies'] = np.zeros(numEnergies)
                energies = yieldData['energies']
                #{1: 'histogram', 2: 'lin-lin', 3: 'log-lin', 4: 'lin-log', 5: 'log-log'}
                yieldData['interpolationSchemes'] = np.zeros(numEnergies-1, dtype=np.int)
                interpSchemes =  yieldData['interpolationSchemes']
                #
                line = get_next_line(fid)
                yieldLists = []
                fpDicts = []
                for energyPoint in range(numEnergies):
                    energies[energyPoint] = get_float(line, 1)
                    numFissionProducts = get_int(line, 6)
                    if energyPoint:
                        interpSchemes[energyPoint-1] = get_int(line, 3)
                    fpYieldData, pos, line = get_float_list(0, line, fid, 4*numFissionProducts)
                    #
                    strt, end = 0, 4 * numFissionProducts - 3
                    ZAList = np.array(fpYieldData[strt:end:4], dtype=np.int)
                    ZList, AList = divmod(ZAList, 1000)
                    #
                    strt, end = 1, 4 * numFissionProducts - 2
                    MList = np.array(fpYieldData[strt:end:4], dtype=np.int)
                    #
                    # Temporarily store the yields as a list of lists
                    strt, end = 2, 4 * numFissionProducts - 1
                    yieldLists.append(fpYieldData[strt:end:4])
                    #
                    # Temporarily store the fission products as a dict of (Z,A,M)
                    # whose value is an index in the yieldMatrix
                    ZAMDict = {}
                    for i in range(numFissionProducts):
                        ZAMDict[(ZList[i], AList[i], MList[i])] = i
                    fpDicts.append(ZAMDict)
                    line = get_next_line(fid)
                # Create yield matrix (energy,nuclide). Different energies may have different fp.
                ZAMSet = set()
                for ZAMDict in fpDicts:
                    ZAMSet.update(ZAMDict.keys())
                numNuclides = len(ZAMSet)
                indexToZAMList = []
                ZAMtoIndexDict = {}
                for i, (Z,A,M) in enumerate(sorted(ZAMSet)):
                    indexToZAMList.append((Z,A,M))
                    ZAMtoIndexDict[(Z,A,M)] = i
                yieldMatrix = np.zeros((numNuclides, numEnergies))
                for energyPoint in range(numEnergies):
                    for (Z,A,M) in fpDicts[energyPoint]:
                        nuclideIndexFrom = fpDicts[energyPoint][(Z,A,M)]
                        nuclideIndexTo = ZAMtoIndexDict[(Z,A,M)]
                        yieldMatrix[nuclideIndexTo, energyPoint] = yieldLists[energyPoint][nuclideIndexFrom]
                yieldData['indexToZAM'] = indexToZAMList
                yieldData['ZAMtoIndex'] = ZAMtoIndexDict
                yieldData['yields'] = yieldMatrix

def get_reaction_Q_values(filePathIn, mtsForMF3, outputDict, verbosity=0):
    '''Read Q values (QI parameter) from ENDF file for each desired MT'''
    if 0 in mtsForMF3:
        mtsForMF3 = [18]
    mfmtList = [(3,mt) for mt in mtsForMF3]
    foundMFMTDict = {}
    for (mf,mt) in mfmtList:
        foundMFMTDict[(mf,mt)] = False
    with open(filePathIn, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if (mf, mt) in mfmtList:
                foundMFMTDict[(mf,mt)] = True
                line = get_next_line(fid)
                # The Q-value, in eV
                outputDict[(mf,mt)] = get_float(line, 2)
            if np.all(foundMFMTDict.values()):
                break

####################################################################################
def get_scattering_sizes(filePath, verbosity=False):
    '''Determine the number of secondary groups (summed over incident group, max'd over temperatures)'''
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        startList = []
        mfList = []
        mtList = []
        # Look through the file for each reaction type
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if not(is_continue(mf, mt, lineNum)) and not(is_section_end(mf, mt, lineNum)):
                if mf == 6:
                    startList.append(pos)
                    mfList.append(mf)
                    mtList.append(mt)
        # Find the number of groups
        fid.seek(0)
        for i in range(3):
            line = get_next_line(fid)
        numGroups = get_int(line, 3)
        # Find the secondary size for each reaction
        secondarySizeDict = {}
        for mf, mt, startPos in zip(mfList, mtList, startList):
            secondarySizeDict[(mf,mt)] = 0
        for mf, mt, startPos in zip(mfList, mtList, startList):
            fid.seek(startPos)
            #Have to look at each group
            line = get_next_line(fid)
            line = get_next_line(fid)
            temperature = get_float(line, 1)
            groupIndex = 0
            secondarySize = 0
            # Read transfer matrix from file and store in CSC-like format in xsDict
            while groupIndex < numGroups:
                numSecondaryPos = get_int(line, 3)
                indexLowestGroup = get_int(line, 4)
                numEntries = get_int(line, 5)
                groupIndex = get_int(line, 6)
                pos = 0
                line = get_next_line(fid)
                #xs, pos, line = get_list(pos, line, fid, numEntries, get_float)
                xs, pos, line = get_float_list(pos, line, fid, numEntries)
                secondarySize += numSecondaryPos - 1
                if verbosity > 2:
                    print '{0} {1:3g} {2:5.1f} {3:3g} {4:3g}'.format(mf, mt, temperature, groupIndex, numSecondaryPos - 1)
            if verbosity > 1:
                print 'total: {0}, {1:3g}, {2:5.1f}, {3:3g}'.format(mf, mt, temperature, secondarySize)
            secondarySizeDict[(mf,mt)] = max(secondarySizeDict[(mf,mt)], secondarySize)
        if verbosity:
            print 'secondarySizeDict (max over temperatures)'
            print secondarySizeDict
        return secondarySizeDict

def get_reaction_list(filePath, verbosity):
    '''Get all the reactions in the file.'''
    print filePath
    with open(filePath, 'r') as fid:
        line = get_next_line(fid)
        print line
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        startList = []
        mfList = []
        mtList = []
        # Look through the file for each reaction type
        while line != '':
            line, pos = seek_to_next_entry(fid, mf, mt)
            if line == '':
                break
            mf, mt, lineNum = get_mf_mt_lineNum(line)
            if not(is_continue(mf, mt, lineNum)) and not(is_section_end(mf, mt, lineNum)):
                startList.append(pos)
                mfList.append(mf)
                mtList.append(mt)
        # Find the number of groups
        fid.seek(0)
        for i in range(3):
            line = get_next_line(fid)
        numGroups = get_int(line, 3)
        # For each reaction, figure out the temperatures, Leg moments, and sig0
        thermList = []
        numLegList = []
        numSig0List = []
        for startPos in startList:
            fid.seek(startPos)
            line = get_next_line(fid)
            numLegList.append(get_int(line, 3))
            numSig0List.append(get_int(line, 4))
            line = get_next_line(fid)
            thermList.append(get_float(line, 1))
        if verbosity:
            print 'File contains the following reactions'
            print '{0:4s} {1:2s} {2:26s} {3} {4} {5}'.format('MF', 'MT', 'Name', 'numLeg', 'numSig0', 'T(K)')
        pdtMTList, combineTransferList, validMTsForMF3, validMTsForMF5, validMTsForMF6 = get_mt_lists()
        mt2long = get_mt2long_dict(pdtMTList)
        for mf, mt, numLeg, numSig0, T in zip(mfList, mtList, numLegList, numSig0List, thermList):
            if mf == 3 or mf == 6:
                mtEff = get_pdt_mt(mf, mt)
                try:
                    name = mt2long[mtEff]
                except KeyError:
                    name = 'unknown'
            elif mf == 1 and mt == 451:
                name = 'energyBounds'
            elif mf == 5 and mt == 455:
                try:
                    name = mt2long[2055]
                except KeyError:
                    name = 'unknown'
            else:
                name = ''
            if verbosity:
                print '{0:2g} {1:4g} {2:30s} {3:2g} {4:2g} {5:8.1f}'.format(mf, mt, name, numLeg, numSig0, T)
        uniqueMFMTs = set([ (mf, mt) for mf,mt in zip(mfList, mtList) if mf in (3,5,6) ])
        uniqueTemperatures = sorted(set(thermList))
        numLegMoments = 1
        numLegMomentsList = [numLeg for mf, numLeg in zip(mfList, numLegList) if mf == 6]
        if numLegMomentsList:
            numLegMoments = max(numLegMomentsList)
        numSig0 = max(numSig0List)
    return uniqueMFMTs, numGroups, numLegMoments, uniqueTemperatures, numSig0

####################################################################################
def get_list(pos, line, fid, listLength, get_method):
    '''Read a list of length listLength from fid. Advances fid.'''
    newList = []
    # Be careful about line wrapping
    numbersPerLine = 6
    for i in range(listLength):
        pos = (pos % numbersPerLine) + 1
        if pos == 1:
            line = get_next_line(fid)
        newList.append(get_method(line, pos))
    return newList, pos, line

def get_float_list(pos, line, fid, listLength):
    '''Optimized version to read a list of floats of length listLength from fid. Advances fid.'''
    # Be careful about line wrapping
    numbersPerLine = 6
    stride = 11
    sz = 0
    numLines = listLength / numbersPerLine
    strrs = []
    # Read full lines
    for i in range(numLines):
        line = fid.readline()[0:80]
        for j in range(numbersPerLine):
            strrs.append(line[stride*j:stride*(j+1)])
    # Read partial lines
    if listLength % numbersPerLine != 0:
        line = fid.readline()[0:80]
        for j in range(listLength % numbersPerLine):
            strrs.append(line[stride*j:stride*(j+1)])
    endPos = ((listLength - 1) % numbersPerLine) + 1
    newList = map(njoy_to_float, strrs)
    return newList, endPos, line

def njoy_to_float(strr):
    # try:
    #     return float(strr[0:9] + 'e' + strr[9:11])
    # except ValueError:
    #     # TODO: FIX!!! (Recover from this error instead of passing it along; interpolate?)
    if strr == '   -inf+***':
        return 0.0
    if strr == '    NaN+***':
        return 0.0
    if not(strr.strip()):
            return 0.0
    numStrr = strr[0]
    for char in strr[1:]:
        if char == '+' or char == '-':
            numStrr += 'e'
        numStrr += char
    num = float(numStrr)
    return num

def get_float(line, loc=1):
    strr = get_str(line, loc)
    #Null case
    if not(strr.strip()):
        return 0.0
    numStrr = ''
    for char in strr:
        if char == '+' or char == '-':
            numStrr += 'e'
        numStrr += char
    num = float(numStrr)
    return num

def get_float_old(line, loc=1):
    strr = get_str(line, loc)
    # Null case
    if not(strr.strip()):
        return 0.0
    # Figure out where the exponent is
    # exponentPos == -1 means not found
    # exponentPos == 0 means sign, not exponent location
    exponentPos = strr.rfind('+')
    if exponentPos > 0:
        numStrr = '{0}e{1}'.format(strr[:exponentPos], strr[exponentPos:])
        return float(numStrr)
    exponentPos = strr.rfind('-')
    if exponentPos > 0:
        numStrr = '{0}e{1}'.format(strr[:exponentPos], strr[exponentPos:])
        return float(numStrr)
    # No exponent and not null, so safe to use value of string
    # (signed zeros still have correct logical comparisons)
    return float(strr)

def get_int(line, loc=1):
    strr = get_str(line, loc)
    # Null case
    if not(strr.strip()):
        return 0
    return int(strr)

def get_str(line, loc=1):
    stride = 11
    strt = stride * (loc - 1)
    end = stride * loc
    return line[strt:end]

####################################################################################
def get_next_line(fid):
    return fid.readline()[0:80]

def get_mf_mt_lineNum(line):
    return int(line[70:72]), int(line[72:75]), int(line[75:80])

def is_continue(mf, mt, lineNum):
    if all((mf == 0, mt == 0, lineNum == 0)):
        return True
    return False

def is_section_end(mf, mt, lineNum):
    if mt == 0 and lineNum == 99999:
        return True
    return False

####################################################################################
def seek_to_next_useful_entry(fid, mf, mt, verbosity=0):
    '''Skip this (MF,MT) pair and advance fid until the next non-trivial one is found'''
    if verbosity > 1:
        print 'Skipping MF {0}, MT {1:3g}'.format(mf, mt)
    line, pos = seek_to_next_entry(fid, mf, mt)
    if line == '':
        return line
    mf, mt, lineNum = get_mf_mt_lineNum(line)
    while is_continue(mf, mt, lineNum) or is_section_end(mf, mt, lineNum):
        line = get_next_line(fid)
        if line == '':
            return line
        mf, mt, lineNum = get_mf_mt_lineNum(line)
    return line

def seek_to_next_entry(fid, currentMf, currentMt):
    '''Skip this (MF,MT) pair and advance fid until the next one is found'''
    finished = False
    while not(finished):
        pos = fid.tell()
        line = get_next_line(fid)
        if line == '':
            finished = True
            break
        mf, mt, lineNum = get_mf_mt_lineNum(line)
        if mf != currentMf or mt != currentMt:
            finished = True
            break
    return line, pos

####################################################################################
def get_endf_mt_list():
    '''Provides a subset of the officially-recognized ENDF MT numbers and their meaning.
    See Tables 4, 21, and 25 in NJOY 2016 manual for MT numbers and names of thermal XS
    References:
    endf6: ENDF-6 Formats Manual, BNL-90365-2009 Rev.2 (Description of the ENDF-6 format)
    njoy2016: The NJOY Nuclear Data Processing System, Version 2016, LA-UR-17-20093 (NJOY 2016 manual)
    '''
    endfMTList = []
    # Neutrons
    endfMTList.append((1, 'ntot', 'nP0Total'))
    endfMTList.append((2, 'nelas', 'nElastic'))
    endfMTList.append((4, 'nineltot', 'nTotalInelastic'))
    endfMTList.append((16, 'n2n', 'n2n'))
    endfMTList.append((17, 'n3n', 'n3n'))
    endfMTList.append((18, 'nf', 'fission'))
    endfMTList.append((22, 'nna', 'nnAlpha'))
    endfMTList.append((28, 'nnp', 'nnProton'))
    endfMTList.append((37, 'n4n', 'n4n'))
    # Inelastic
    for i in range(51, 90+1):
        iLevel = i - 51 + 1
        endfMTList.append((i, 'ninel{0:02g}'.format(iLevel), 'nInelastic{0:02g}'.format(iLevel)))
    endfMTList.append((91, 'ncontinelas', 'nContinuumInelastic'))
    #
    endfMTList.append((102, 'ng', 'nGamma'))
    endfMTList.append((103, 'np', 'nProton'))
    endfMTList.append((104, 'nd', 'nDeuteron'))
    endfMTList.append((105, 'nt', 'nTriton'))
    endfMTList.append((106, 'nhe3', 'nHe3'))
    endfMTList.append((107, 'na', 'nAlpha'))
    # Thermal
    #Benzene and BeO are special cases (see manual)
    endfMTList.append((221, 'free', 'freeThermal'))
    endfMTList.append((222, 'hh2o', 'HinH2OThermal'))
    endfMTList.append((223, 'hch2inel', 'HinCH2ThermalInelastic'))
    endfMTList.append((224, 'hch2elas', 'HinCH2ThermalElastic'))
    endfMTList.append((225, 'hzrhinel', 'HinZrHThermalInelastic'))
    endfMTList.append((226, 'hzrhelas', 'HinZrHThermalElastic'))
    endfMTList.append((227, 'benz', 'benzeneThermal'))
    endfMTList.append((228, 'dd2o', 'DinD2OThermal'))
    endfMTList.append((229, 'graphinel', 'CinGraphiteThermalInelastic'))
    endfMTList.append((230, 'graphelas', 'CinGraphiteThermalElastic'))
    endfMTList.append((231, 'beinel', 'BeThermalInelastic'))
    endfMTList.append((232, 'beelas', 'BeThermalElastic'))
    endfMTList.append((233, 'bebeoinel', 'BeInBeOThermalInelastic'))
    endfMTList.append((234, 'bebeoelas', 'BeInBeOThermalElastic'))
    endfMTList.append((235, 'zrzrhinel', 'ZrinZrHThermalInelastic'))
    endfMTList.append((236, 'zrzrhelas', 'ZrinZrHThermalElastic'))
    endfMTList.append((237, 'obeoinel', 'OInBeOThermalInelastic'))
    endfMTList.append((238, 'obeoelas', 'OInBeOThermalElastic'))
    endfMTList.append((239, 'ouo2inel', 'OinUO2ThermalInelastic'))
    endfMTList.append((240, 'ouo2elas', 'OinUO2ThermalElastic'))
    endfMTList.append((241, 'uuo2inel', 'UinUO2ThermalInelastic'))
    endfMTList.append((242, 'uuo2elas', 'UinUO2ThermalElastic'))
    endfMTList.append((243, 'alinel', 'AlThermalInelastic'))
    endfMTList.append((244, 'alelas', 'AlThermalElastic'))
    endfMTList.append((245, 'feinel', 'FeThermalInelastic'))
    endfMTList.append((246, 'feelas', 'FeThermalElastic'))
    #
    endfMTList.append((259, 'invel', 'inverseVelocity'))
    endfMTList.append((301, 'etot', 'energyReleaseTot'))
    endfMTList.append((318, 'efission', 'energyReleaseFission'))
    endfMTList.append((452, 'nutot', 'totalNu'))
    endfMTList.append((455, 'nudelay', 'delayedNu'))
    endfMTList.append((456, 'nuprompt', 'promptNu'))
    # Gammas
    endfMTList.append((501, 'gtot', 'gammaTotal'))
    return endfMTList

def get_neutron_transfer_list():
    '''
    It makes sense to have the following MT numbers be associated with
    transfer matrices for neutrons to neutrons
    (viz., they must produce neutrons from neutrons
    and they must be a valid MF 6 option in GROUPR).
    '''
    neutronTransferList = [2, 16, 17, 18, 22, 28, 37]
    for i in range(51, 91+1):
        neutronTransferList.append(i)
    for i in range(221, 246+1):
        neutronTransferList.append(i)
    return neutronTransferList

def get_gamma_transfer_list():
    '''Gamma xfer matrices are currently NOT supported because of the inconsistency between PDT and ENDF in the naming schemes'''
    gammaTransferList = []
    return gammaTransferList

def get_mts_for_combining():
    '''
    Use MT 2501 for the sum of all transfer reactions except fission
    (new versions of PDT use MT 2519 instead of MT 2501)
    One (and only one) scattering source should be added to 2501 at the very end.
    Except in the case of Be (see manual), the thermal scattering xfer matrix
       should overwrite the existing 2501 values and the total xs should be updated
       to be consistent with the new scattering xs in the thermal range.
    '''
    combineTransferList = [2, 16, 17, 22, 28, 37]
    for i in range(51, 91+1):
        # Inelasti scatterings
        combineTransferList.append(i)
    return combineTransferList

def get_pdt_mt_list(endfMTList, neutronTransferList, gammaTransferList):
    '''Maps ENDF MT nubmers to PDT-MT's. Determines which MTs to combine and which can be MF 6s.'''
    #There is no concept of MF number in PDT.
    #The MT number for xs (MF 3) stays the same.
    #The MT number of transfer matrices (MF 6) are equal to the MT numbers
    #   of their MF 3 counterparts plus 2500.
    pdtMTList = []
    #Append all xs
    for endfMT in endfMTList:
        pdtMTList.append(endfMT)
    #PN (N>1) Legendre moments of the spectrum will use MT 4000 + N
    #PN (N>1) Legendre moments of the total cross section will use MT 4100 + N
    #This limits the number of Legendre moments to 100.
    #This preserves the notion of a MF 3 (xs) being a vector without dependence on Legendre moment.
    #All transfer matrices may have multiple Legendre moments in one MT number.
    for n in range(1,100+1):
        wgtMT = n + 4000
        pdtMTList.append((wgtMT, 'nwgt{0:02g}'.format(n), 'nP{0:02g}Spectrum'.format(n)))
    for n in range(1,100+1):
        totMT = n + 4100
        pdtMTList.append((totMT, 'ntot{0:02g}'.format(n), 'nP{0:02g}Total'.format(n)))
    #The flux and chi (fission spectrum) are given MT numbers where they had none before.
    pdtMTList.append((1099, 'nwgt', 'nP0Spectrum'))
    pdtMTList.append((1018, 'chi', 'nFissionSpectrum'))
    pdtMTList.append((1452, 'nufission', 'nNuFission'))
    pdtMTList.append((2055, 'chid', 'nDelayedFissionSpectrum'))
    pdtMTList.append((2501, 'xfernofission', 'combinedNonFissionTransfer'))
    #Append valid neutron xfer matrices
    for endfMT in endfMTList:
        mt, shortName, longName = endfMT
        if mt in neutronTransferList:
            mtX = mt + 2500
            shortNameX = shortName + 'xfer'
            longNameX = longName + 'Transfer'
            pdtMTList.append((mtX, shortNameX, longNameX))
    #Append some common gamma xfer matrices.
    #These align with the PDT xs converter nomenclature, not the ENDF
    #(ENDF: coherent xs is 502, incoherent is 504, pair prod is 515,517 or 516)
    pdtMTList.append((3501, 'gcohxfer', 'gammaCoherentTransfer'))
    pdtMTList.append((3502, 'gincohxfer', 'gammaIncoherentTransfer'))
    pdtMTList.append((3503, 'gpairxfer', 'gammaPairProductionTransfer'))
    return pdtMTList

def get_valid_mts_lists(pdtMTList, neutronTransferList, gammaTransferList):
    '''Determine the valid MT numbers for MF's 3, 5, and 6'''
    validMTsForMF3 = []
    validMTsForMF5 = []
    for pdtMT in pdtMTList:
        mt, shortName, longName = pdtMT
        if mt < 1000:
            validMTsForMF3.append(mt)
        if mt == 2055:  #special case for delayedChi
            validMTsForMF5.append(455)
    validMTsForMF6 = np.concatenate((neutronTransferList, gammaTransferList))
    validMTsForMF3 = sorted(validMTsForMF3)
    validMTsForMF5 = sorted(validMTsForMF5)
    validMTsForMF6 = sorted(validMTsForMF6)
    return validMTsForMF3, validMTsForMF5, validMTsForMF6

def get_mt_lists():
    '''Return which MT lists'''
    endfMTList = get_endf_mt_list()
    neutronTransferList = get_neutron_transfer_list()
    gammaTransferList = get_gamma_transfer_list()
    combineTransferList = get_mts_for_combining()
    pdtMTList = get_pdt_mt_list(endfMTList, neutronTransferList, gammaTransferList)
    validMTsForMF3, validMTsForMF5, validMTsForMF6 = get_valid_mts_lists(pdtMTList, neutronTransferList, gammaTransferList)
    return pdtMTList, combineTransferList, validMTsForMF3, validMTsForMF5, validMTsForMF6

####################################################################################
def lookup_num_sig0(maxNumSig0, mf, mt):
    '''Look up the number of background xs the reaction contains (either maxNumSig0 or 1).'''
    if mf == 3:
        # NJOY 99 -- 2016:
        mtsAtMultSig0 = {1, 2, 18, 102, 301, 318}
        # NJOY 2012 -- 2016:
        mtsAtMultSig0 |= {i for i in range(51,91+1)}
        if mt in mtsAtMultSig0:
            return maxNumSig0
        else:
            return 1
    elif mf == 6:
        mtsAtMultSig0 = [2]
        if mt in mtsAtMultSig0:
            return maxNumSig0
        else:
            return 1
    else:
        return 1

def lookup_num_therm(maxNumTherm, mf, mt):
    '''Look up the number of temperatures the reaction contains (either maxNumTherm or 1).'''
    thermalMTs = range(221, 246+1)
    if mf == 3:
        mtsDependOnT = [1, 2, 18, 102]
        if (mt in mtsDependOnT) or (mt in thermalMTs):
            return maxNumTherm
        else:
            return 1
    elif mf == 6:
        mtsDependOnT = [2]
        if (mt in mtsDependOnT) or (mt in thermalMTs):
            return maxNumTherm
        else:
            return 1
    else:
        return 1

def determine_rxns_to_keep(desiredMTsForMF3, desiredMTsForMF5, desiredMTsForMF6, mfmtsAvail, verbosity):
    '''Determine which reactions in the GENDF file to read in and store'''
    #
    mfmtsDesired = []
    for mt in desiredMTsForMF3:
        mfmtsDesired.append((3, mt))
    for mt in desiredMTsForMF5:
        mfmtsDesired.append((5, mt))
    for mt in desiredMTsForMF6:
        mfmtsDesired.append((6, mt))
    mfmtsDesired = set(mfmtsDesired)
    #
    mfmtsValid = []
    pdtMTList, combineTransferList, validMTsForMF3, validMTsForMF5, validMTsForMF6 = get_mt_lists()
    for mt in validMTsForMF3:
        mfmtsValid.append((3, mt))
    for mt in validMTsForMF5:
        mfmtsValid.append((5, mt))
    for mt in validMTsForMF6:
        mfmtsValid.append((6, mt))
    mfmtsValid = set(mfmtsValid)
    #
    availMTsForMF3 = [(mf, mt) for mf,mt in mfmtsAvail if mf == 3]
    availMTsForMF5 = [(mf, mt) for mf,mt in mfmtsAvail if mf == 5]
    availMTsForMF6 = [(mf, mt) for mf,mt in mfmtsAvail if mf == 6]
    #
    mfmtsKeep = list(mfmtsAvail & mfmtsDesired & mfmtsValid)
    if (3, 0) in mfmtsDesired:
        mfmtsKeep += list(set(availMTsForMF3) & mfmtsValid)
    if (5, 0) in mfmtsDesired:
        mfmtsKeep += list(set(availMTsForMF5) & mfmtsValid)
    if (6, 0) in mfmtsDesired:
        mfmtsKeep += list(set(availMTsForMF6) & mfmtsValid)
    mfmtsKeep = set(mfmtsKeep)
    mfmtsKeepSorted = list(sorted(mfmtsKeep))
    #
    if verbosity >=2:
        print_mts(pdtMTList, combineTransferList, validMTsForMF3, validMTsForMF5, validMTsForMF6)
    if verbosity:
        mt2long = get_mt2long_dict(pdtMTList)
        names = [mt2long[get_pdt_mt(mf, mt)] for mf,mt in mfmtsKeepSorted]
        print_kept_mts(mfmtsKeepSorted, names)
    return mfmtsKeep

def get_mt2long_dict(mtList):
    '''Get dictionary for mt number to long name.'''
    mt2long = {}
    for mtDat in mtList:
        mt, shortName, longName = mtDat
        mt2long[mt] = longName
    return mt2long

def get_short2mt_dict(mtList):
    '''Get dictionary for short name to mt number.'''
    short2mt = {}
    for mtDat in mtList:
        mt, shortName, longName = mtDat
        short2mt[shortName] = mt
    return short2mt

def get_mt2short_dict(mtList):
    '''Get dictionary for mt number to short name.'''
    mt2short = {}
    for mtDat in mtList:
        mt, shortName, longName = mtDat
        mt2short[mt] = shortName
    return mt2short

def get_pdt_mt(mf, mt):
    '''PDT does not use MFs. Map ENDF (MF,MT) to PDT-MT.'''
    if mf == 3:
        return mt
    else:
        return 2500 + mt

####################################################################################
def print_mts(pdtMTList, combineTransferList, validMTsForMF3, validMTsForMF5, validMTsForMF6):
    print 'PDT rxns'
    for pdtMT in sorted(pdtMTList):
        print pdtMT
    print 'Combine these rxns into 2501'
    print combineTransferList
    print 'Valid ENDF MT numbers for MF 3'
    print validMTsForMF3
    print 'Valid ENDF MT numbers for MF 5'
    print validMTsForMF5
    print 'Valid ENDF MT numbers for MF 6'
    print validMTsForMF6

def print_kept_mts(mfmts, names):
    for mf2print in [3,6]:
        print 'Keeping these ENDF MT numbers for MF {0}'.format(mf2print)
        for (mf,mt), name in filter(lambda ((mf,mt), name): mf==mf2print, zip(mfmts, names)):
            print '{0:3g} ({1})'.format(mt, name)

def print_xs_summary(filePathIn, data, verbosity):
    filenameOut = '{0:02g}-{1:03g}.data'.format(data['Z'], data['A'])
    filePathOut = os.path.join(os.path.split(filePathIn)[0], filenameOut)
    mfmtsSet = data['mfmts']
    mfmtsSorted = sorted(mfmtsSet)
    if verbosity > 1:
        print 'Summary of XS...'
        print 'MF MT numLeg numSig0 numTherm'
        for (mf, mt) in mfmtsSorted:
            print '{0} {1:3g} {2:2g} {3:2g} {4:2g}'.format(
                mf, mt,
                data['rxn'][(mf,mt)]['numLegMoments'],
                data['rxn'][(mf,mt)]['numSig0'],
                data['rxn'][(mf,mt)]['numTherm'])

####################################################################################
def find_interpolation_indices_fractions(value, sortedArray, monotonicFunction=(lambda x: x)):
    '''Finds the two points in sortedArray closest to value using the distance metric monotonicFunction. Finds a convex list of fractions for interpolation, which is linear wrt monotonicFunction. If value is outside sortedArray, uses closest value (no extrapolation).'''
    fValue = monotonicFunction(value)
    fArray = monotonicFunction(sortedArray)
    nearestIndex = np.argmin(np.abs(fArray - fValue))
    nearestIndices = [nearestIndex]
    lastIndex = len(sortedArray) - 1
    if value < sortedArray[nearestIndex] and nearestIndex != 0:
        nearestIndices.append(nearestIndex-1)
    elif value > sortedArray[nearestIndex] and nearestIndex != lastIndex:
        nearestIndices.append(nearestIndex+1)
    nearestIndices = sorted(nearestIndices)
    if len(nearestIndices) > 1:
        fraction = (fValue - fArray[nearestIndices[0]]) / (fArray[nearestIndices[1]] - fArray[nearestIndices[0]])
        fractions = np.array([1.0-fraction, fraction])
    else:
        fractions = np.array([1.0])
    return nearestIndices, fractions

def find_interpolation_indices_fractions_array(inArray, sortedGrid, monotonicFunction=(lambda x: x), reverse=False):
    '''For each value of inArray, finds the two closest points in sortedGrid using the distance metric monotonicFunction. Finds a convex list of fractions for interpolation, which is linear wrt monotonicFunction. If value is outside sortedArray, uses closest value (no extrapolation). If reverse is True, will flip indices.'''
    inSize = len(inArray)
    gridSize = len(sortedGrid)
    fractionMat = np.zeros((gridSize, inSize))
    fIn = monotonicFunction(inArray)
    fGrid = monotonicFunction(sortedGrid)
    lastIndex = len(sortedGrid) - 1
    for i in range(inSize):
        value = inArray[i]
        fValue = fIn[i]
        nearestIndex = np.argmin(np.abs(fGrid - fValue))
        nearestIndices = [nearestIndex]
        if value < sortedGrid[nearestIndex] and nearestIndex != 0:
            nearestIndices.append(nearestIndex-1)
        elif value > sortedGrid[nearestIndex] and nearestIndex != lastIndex:
            nearestIndices.append(nearestIndex+1)
        nearestIndices = np.array(sorted(nearestIndices))
        if len(nearestIndices) > 1:
            fraction = (fValue - fGrid[nearestIndices[0]]) / (fGrid[nearestIndices[1]] - fGrid[nearestIndices[0]])
            fractions = np.array([1.0-fraction, fraction])
        else:
            fractions = np.array([1.0])
        if reverse:
            nearestIndices = gridSize - nearestIndices - 1
        fractionMat[nearestIndices, i] = fractions
    return fractionMat

####################################################################################
def pickle_xs(data, filePathOut, binaryOpt='b'):
    writeStr = 'w'
    if binaryOpt:
        writeStr = 'wb'
        protocol = pickle.HIGHEST_PROTOCOL
    pickle.dump(data, open(filePathOut, writeStr), protocol)

def unpickle_xs(picklePath, binaryOpt='b'):
    readStr = 'r'
    if binaryOpt:
        readStr = 'rb'
    return pickle.load(open(picklePath, readStr))

def multiline_string(vector, spacing, numberPerLine, decimals):
    outStr = ''
    N = int(np.ceil(len(vector)/float(numberPerLine)))*numberPerLine
    for i in range(numberPerLine, N+1, numberPerLine):
        strs = ['{0:>{1}.{2}g}'.format(vi, spacing, decimals) for vi in vector[i-numberPerLine:i]]
        outStr += ''.join(strs) + '\n'
    return outStr

####################################################################################
def condense_xs(xsDataIn, energyMesh, flux, verbosity):
    '''Condense the cross sections based on the energy mesh, which maps existing groups
    to energy elements, which may be discontiguous. If the flux is specified, use it
    instead of the flux in xsDataIn. Return xsData, which has the xs on the elements.'''
    #
    xsData = {}
    # >>> Determine energy mesh properties <<<
    numGroups = len(energyMesh)
    if numGroups != xsDataIn['numGroups']:
        print 'Error. Input energyMesh should contain the same number of groups as the GENDF file.'
        exit(1)
    energyMesh -= np.min(energyMesh)
    meshElements = np.unique(energyMesh)
    numElements = len(meshElements)
    # np.bincount applies the scatter pattern: out[n] += 1 if in[i] == n
    numSubelementsPerElement = np.bincount(energyMesh)
    largeElements = np.where(numSubelementsPerElement > 1)[0]
    # Nothing needs to be done if the energyMesh is not condensing
    numLargeElements = np.sum(largeElements)
    if not numLargeElements:
        return xsDataIn
    # These are the first/last group indices that are to be condensed
    # np.in1d returns true at pos i if arg1[i] is in arg2
    firstLargeGroupLoc = np.where(np.in1d(energyMesh, largeElements))[0][0]
    lastLargeGroupLoc = np.where(np.in1d(energyMesh, largeElements))[0][-1]
    # The element locations are determined to be consistent with the group locations
    # The element locations correspond to indices in meshElements list
    firstLargeElementLoc = firstLargeGroupLoc
    lastLargeElementLoc = numElements - (numGroups - lastLargeGroupLoc)
    #
    # Fill in element boundaries using the energy mesh.
    # np.bincount applies the scatter pattern: out[n] += weight[i] if in[i] == n
    groupSizes = -np.diff(xsDataIn['groupBdrs'])
    elementSizes = np.bincount(energyMesh, weights=groupSizes)
    elementBdrs = xsDataIn['groupBdrs'][0] * np.ones(numElements + 1)
    elementBdrs[1:] -= np.cumsum(elementSizes)
    # Due to finite precision, use the original boundaries for uncondensed groups, including those
    # that border the condensed region. The "+1"'s convert from group indices to boundary indices.
    elementBdrs[0:firstLargeElementLoc+1] = xsDataIn['groupBdrs'][0:firstLargeGroupLoc+1]
    elementBdrs[lastLargeElementLoc+1:] = xsDataIn['groupBdrs'][lastLargeGroupLoc+1:]
    xsData['groupBdrs'] = elementBdrs
    xsData['numGroups'] = numElements
    #
    # >>> Determine the flux <<<
    # Determine the fine-group sig0- and T-dependent flux
    numThermFlux, numSig0Flux, numGroupFlux = xsDataIn['flux'].shape
    if flux is None:
        flux = xsDataIn['flux']
    if flux.shape[-1] != numGroupFlux:
        print flux.shape[-1]
        print numGroupFlux
        print 'Error. Input flux should contain the same number of groups as the GENDF file.'
        exit(1)
    if len(flux.shape) == 1:
        fluxSmall = flux
        flux = np.zeros((numThermFlux, numSig0Flux, numGroupFlux))
        flux[:,:,:] = fluxSmall[np.newaxis, np.newaxis, :]
    elif flux.shape != xsDataIn['flux'].shape:
        print 'Error. If input flux is an ndarray of dimension 3, it must contain the same number of sigma0 and T values as the GENDF file.'
        exit(1)
    else:
        print 'Error. Input flux must either be an ndarray of dimension 1 or 3.'
        exit(1)
    #
    # Define the elementwise or "coarse-group" sig0- and T-dependent flux
    # NB, the lambda does not copy local variables; energyMesh is referenced at invocation time
    # apply_along_axis applies condenseFunc to each 1D slice of flux. The slice is over groups (index 2)
    condenseFunc = lambda x_flux: np.bincount(energyMesh, weights=x_flux)
    elementFlux = np.apply_along_axis(condenseFunc, 2, flux)
    xsData['flux'] = elementFlux
    #
    # >>> Copy over unchanged parameters <<<
    # Copy over nuclide-wide group-structure-agnostic information; arrays are aliased
    keysToCopy = ['Z', 'A', 'weight', 'mfmts', 'numLegMoments', 'numTherm', 'numSig0', 'sig0List', 'thermList']
    for key in keysToCopy:
        xsData[key] = xsDataIn[key]
    #
    # Copy over group-structure-agnostic information for each reaction
    xsData['rxn'] = {}
    for (mf,mt) in xsData['mfmts']:
        xsData['rxn'][(mf,mt)] = {}
        rxn = xsData['rxn'][(mf,mt)]
        rxnIn = xsDataIn['rxn'][(mf,mt)]
        keysToCopy = ['numTherm', 'numSig0', 'numLegMoments']
        if (mf,mt) == (5,455):  # if delayed fission spectrum, numLegMoments is replaced by numDNGs
            keysToCopy = ['numTherm', 'numSig0', 'numDNGs']
        for key in keysToCopy:
            rxn[key] = rxnIn[key]
    #
    # >>> Condense the cross sections for MF 3 <<<
    # Reactions with numThermRxn == 1 or numSig0Rxn == 1 use the flux at iTherm = 0 or iSig0 = 0.
    # The flux at energies at which these reactions occur is not sensitive its thermal or sig0 index.
    mfmtList = [(mf,mt) for (mf,mt) in xsDataIn['mfmts'] if (mf == 3 and mt != 1018)]
    for (mf,mt) in mfmtList:
        rxn = xsData['rxn'][(mf,mt)]
        numThermRxn = rxn['numTherm']
        numSig0Rxn = rxn['numSig0']
        rxn['xs'] = np.zeros((numThermRxn, numSig0Rxn, numElements))
        xs = rxn['xs']
        xsIn = xsDataIn['rxn'][(mf,mt)]['xs']
        for iTherm in range(numThermRxn):
            for iSig0 in range(numSig0Rxn):
                wgt = flux[iTherm, iSig0, :] * xsIn[iTherm, iSig0, :]
                norm = elementFlux[iTherm, iSig0, :]
                xs[iTherm, iSig0, :] = np.bincount(energyMesh, weights=wgt) / norm
    #
    # >>> Condense the cross sections for MF 5 MT 455 (delayed neutron spectrum) <<<
    # Assuming insensitive to thermal and sig0 index
    mfmtList = [(mf,mt) for (mf,mt) in xsDataIn['mfmts'] if (mf == 5 and mt == 455)]
    for (mf,mt) in mfmtList:
        rxn = xsData['rxn'][(mf,mt)]
        numThermRxn = rxn['numTherm']
        numSig0Rxn = rxn['numSig0']
        numDNGs = rxn['numDNGs']
        delayedChiIn = xsDataIn['rxn'][(mf,mt)]['delayedChi']
        # for iDNG in range(numDNGs):
        #     wgt = xsIn[iDNG, :]
        #     delayedChi[iDNG, :] = np.bincount(energyMesh, weights=wgt)
        rxn['delayedChi'] = np.apply_along_axis(condenseFunc, 1, delayedChiIn)
        rxn['decayConst'] = xsDataIn['rxn'][(mf,mt)]['decayConst']
    #
    # >>> Condense the cross sections for MF 6 (all but fission) <<<
    # The scalar flux is used to weight all Legendre moments.
    mfmtList = [(mf,mt) for (mf,mt) in xsDataIn['mfmts'] if (mf == 6 and mt != 18)]
    for (mf,mt) in mfmtList:
        rxn = xsData['rxn'][(mf,mt)]
        rxnIn = xsDataIn['rxn'][(mf,mt)]
        numThermRxn = rxn['numTherm']
        numSig0Rxn = rxn['numSig0']
        numLegMomentsRxn = rxn['numLegMoments']
        # We have indexing for g <- g' (to group, from group). We determine both
        # intermediate indexing for g <- e' (to group, from element), and
        # final indexing for e <- e' (to element, from element).
        indPtrIn = rxnIn['indexPtr']
        rowStartIn = rxnIn['rowStart']
        rowSizeIn = np.diff(indPtrIn)
        rowEndIn = rowStartIn + rowSizeIn
        # rowStartIntermediate[e'] is the first group g scattered to from element e'
        rowStartIntermediate = np.zeros(numElements, dtype=np.int)
        # rowSizeIntermediate[e'] is the number of groups g scattered to from element e'
        rowSizeIntermediate = np.zeros(numElements, dtype=np.int)
        # rowStart[e'] is the first element e scattered to from element e'
        rowStart = np.zeros(numElements, dtype=np.int)
        # rowSize[e'] is the number of elements e scattered to from element e'
        rowSize = np.zeros(numElements, dtype=np.int)
        # Deal with empty rows
        sentinelVal = np.max(rowStartIn) + 1
        rowStartInSafe = rowStartIn.copy()
        rowStartInSafe[rowStartInSafe == -1] = sentinelVal
        for ie, elemFrom in enumerate(meshElements):
            mask = (energyMesh == elemFrom)
            rowStartIntermediate[ie] = np.min(rowStartInSafe[mask])
            rowEndIntermediate = np.max(rowEndIn[mask])
            if rowStartIntermediate[ie] != sentinelVal:
                rowSizeIntermediate[ie] = rowEndIntermediate - rowStartIntermediate[ie]
                strtG, endG = rowStartIntermediate[ie], rowStartIntermediate[ie] + rowSizeIntermediate[ie]
                groupsTo = np.arange(strtG, endG)
                elementsTo = np.unique(energyMesh[groupsTo])
                strtE, endE = elementsTo[0], elementsTo[-1] + 1
            else:
                rowStartIntermediate[ie] = -1
                rowSizeIntermediate[ie] = 0
                strtE, endE = -1, -1
            rowStart[ie] = strtE
            rowSize[ie] = endE - strtE
        indPtr = np.zeros(numElements + 1, dtype=np.int)
        indPtr[1:] = np.cumsum(rowSize)
        rxn['indexPtr'] = indPtr
        rxn['rowStart'] = rowStart
        # Knowing sizes, we condense the scattering matrix on the elements
        scatMatSize = indPtr[-1]
        rxn['xs'] = np.zeros((numThermRxn, numSig0Rxn, numLegMomentsRxn, scatMatSize))
        xs = rxn['xs']
        xsIn = rxnIn['xs']
        for iTherm in range(numThermRxn):
            for iSig0 in range(numSig0Rxn):
                fluxSlice = flux[iTherm, iSig0, :]
                elementFluxSlice = elementFlux[iTherm, iSig0, :]
                for iLegMoment in range(numLegMomentsRxn):
                    xsSlice = xs[iTherm, iSig0, iLegMoment, :]
                    xsSliceIn = xsIn[iTherm, iSig0, iLegMoment, :]
                    for elemFrom in range(numElements):
                        offsetI = rowStartIntermediate[elemFrom]
                        szI = rowSizeIntermediate[elemFrom]
                        xsCol = np.zeros(szI)
                        groupsFrom = np.where(energyMesh == elemFrom)[0]
                        for groupFrom in groupsFrom:
                            strtG, endG = indPtrIn[groupFrom], indPtrIn[groupFrom + 1]
                            offsetG = rowStartIn[groupFrom] - offsetI
                            szG = rowSizeIn[groupFrom]
                            xsCol[offsetG:offsetG+szG] += xsSliceIn[strtG:endG] * fluxSlice[groupFrom]
                        xsCol /= elementFluxSlice[elemFrom]
                        strtE, endE = indPtr[elemFrom], indPtr[elemFrom + 1]
                        offsetE = rowStart[elemFrom]
                        partialEnergyMesh = energyMesh[offsetI:offsetI+szI] - offsetE
                        xsSlice[strtE:endE] = np.bincount(partialEnergyMesh, weights=xsCol)
    #
    # >>> Condense the cross sections for MF 6, MT 18 (fission) <<<
    mfmtList = []
    if (6,18) in xsDataIn['mfmts']:
        mfmtList = [(6,18)]
    for (mf,mt) in mfmtList:
        rxn = xsData['rxn'][(mf,mt)]
        rxnIn = xsDataIn['rxn'][(mf,mt)]
        numThermRxn = rxn['numTherm']
        numSig0Rxn = rxn['numSig0']
        numLegMomentsRxn = rxn['numLegMoments']
        #
        # First, we condense low-energy spectrum by column to an intermediate form (e <- g').
        lowEnergySpectrumIn = rxnIn['lowEnergySpectrum']
        lowEnergySpectrum = np.apply_along_axis(condenseFunc, 0, lowEnergySpectrumIn)
        norm = np.sum(lowEnergySpectrum)
        if norm:
            lowEnergySpectrum /= norm
        rxn['lowEnergySpectrum'] = lowEnergySpectrum
        #
        # Second, we condense prompt fission spectrum by column to an intermediate form (e <- g').
        promptChiIn = rxnIn['promptChi']
        promptChi = np.apply_along_axis(condenseFunc, 0, promptChiIn)
        norm = np.sum(promptChi)
        if norm:
            promptChi /= norm
        rxn['promptChi'] = promptChi
        #
        # Third, we condense prompt fission nu*sig_f (nu) using flux as weighting
        promptProdIn = rxnIn['promptProd']
        iTherm, iSig0 = 0, 0  # The 0th temperature and sig0 index is used for the weighting flux
        wgt = flux[iTherm, iSig0, :] * promptProdIn
        norm = elementFlux[iTherm, iSig0, :]
        rxn['promptProd'] = np.bincount(energyMesh, weights=wgt) / norm
        #
        # Fourth, we condense the full FissionMatrix
        FissionMatrix = rxnIn['FissionMatrix']
          # condense prompt fission matrix by column to an intermediate form (e <- g').
        FissionMatrixIntermediate = np.apply_along_axis(condenseFunc, 1, FissionMatrix)
        iTherm, iSig0 = 0, 0  # The 0th temperature and sig0 index is used for the weighting flux
        wgt = np.dot(np.diag(flux[iTherm, iSig0, :]), FissionMatrixIntermediate[:,:])
        norm = elementFlux[iTherm, iSig0, :]
          # final condense for e <- e' (to element, from element).
        rxn['FissionMatrix'] = np.dot(np.diag(1./norm), np.apply_along_axis(condenseFunc, 0, wgt))
        #
        # Fifth, we determine which elements are high-energy elements, if any
        lastHighGroup = rxnIn['highestHighEnergyGroup']
        if lastHighGroup == -1:
            # The high-energy portion does not exist.
            rxn['highestHighEnergyGroup'] = rxnIn['highestHighEnergyGroup']
            rxn['highEnergyMatrix'] = rxnIn['highEnergyMatrix']
        else:
            highEnergyGroups = np.arange(0, lastHighGroup + 1)
            highEnergyElements = np.unique(energyMesh[highEnergyGroups])
            # We assume that the elements are intelligently ordered st low elements have higher energies
            lastHighElement = highEnergyElements[-1]
            rxn['highestHighEnergyGroup'] = lastHighElement
            #
            # We condense fission matrix by column to an intermediate form (e <- g').
            highEnergyMatIn = rxnIn['highEnergyMatrix']
            highEnergyMatIntermediate = np.apply_along_axis(condenseFunc, 1, highEnergyMatIn)
            #
            # Finally, we condense the fission matrix by row to the final form (e <- e').
            # The 0th temperature and sig0 index is used for the weighting flux
            rxn['highEnergyMatrix'] = np.zeros((lastHighElement+1, numElements))
            highEnergyMat = rxn['highEnergyMatrix']
            for elem in range(lastHighElement+1):
                mask = (energyMesh == elem)
                groupsInElem = np.where(energyMesh == elem)[0]
                lowEnergyGroupsInElem = groupsInElem[groupsInElem > lastHighGroup]
                highEnergyGroupsInElem = groupsInElem[groupsInElem <= lastHighGroup]
                #
                wgtHigh = flux[0, 0, highEnergyGroupsInElem]
                wgtLow = flux[0, 0, lowEnergyGroupsInElem]
                # The fission matrix in an element with high- and low-energy groups is a
                # flux-weighted average
                xsLow = xsData['rxn'][(3,18)]['xs'][0, 0, elem] * xsData['rxn'][(3,452)]['xs'][0, 0, elem]
                sumLow = 0.
                if xsLow:
                    sumLow = xsLow * np.sum(wgtLow)
                # The intermediate high energy fission matrix is nuSigf_{g'->e}
                highEnergyMat[elem, :] = (
                    np.sum(highEnergyMatIntermediate[highEnergyGroupsInElem, :] *
                        wgtHigh[:, np.newaxis], axis=0) +
                    sumLow * lowEnergySpectrum ) / (np.sum(wgtLow) + np.sum(wgtHigh))

    return xsData

####################################################################################
def remove_sparse_holes(holeyRowStart, holeyColSize):
    '''The following pattern, though contiguous in CSC, is holey in CSR:
    X X X
    X   X
    Given a CSC rowStart / colSize, this function determines the colStart / rowSize for CSR that removes holes.
    Because we are still in CSC, we then convert back to fullRowStart / fullColSize.
    Probably only needs to be done for the combined matrix, because holes do not occur for individual reactions.
    '''
    # Use for individual matrices?
    holeyRowEnd = holeyRowStart + holeyColSize - 1
    numGroups = len(holeyRowEnd)
    identityForMin = numGroups
    colStart = identityForMin * np.ones(numGroups, dtype=np.int)
    colEnd = np.zeros(numGroups, dtype=np.int)
    fullRowStart = identityForMin * np.ones(numGroups, dtype=np.int)
    fullRowEnd = np.zeros(numGroups, dtype=np.int)
    for groupFrom in range(numGroups):
        strt, end = holeyRowStart[groupFrom], holeyRowEnd[groupFrom] + 1
        colStart[strt:end] = np.minimum(groupFrom, colStart[strt:end])
        colEnd[strt:end] = np.maximum(groupFrom, colEnd[strt:end])
    for groupTo in range(numGroups):
        if colStart[groupTo] == numGroups and colEnd[groupTo] == 0:
            colStart[groupTo] = numGroups
            colEnd[groupTo] = numGroups - 1
    for groupTo in range(numGroups):
        strt, end = colStart[groupTo], colEnd[groupTo] + 1
        fullRowStart[strt:end] = np.minimum(groupTo, fullRowStart[strt:end])
        fullRowEnd[strt:end] = np.maximum(groupTo, fullRowEnd[strt:end])
    for groupFrom in range(numGroups):
        if fullRowStart[groupFrom] == numGroups and fullRowEnd[groupFrom] == 0:
            fullRowStart[groupFrom] = numGroups
            fullRowEnd[groupFrom] = numGroups - 1
    fullColSize = fullRowEnd - fullRowStart + 1
    return fullRowStart, fullColSize

def combine_transfer_matrix(data, thermalMTList, thermalMultList=[], verbosity=False):
    '''Create a new transfer matrix that is the sum of all transfer matrices except thermal and fission. Overwrite the thermal portion of this new matrix with the sum of the transfer matrices corresponding to thermal reactions in thermalList. Store in extended union CSC-like structure. The structure is extended so there will be no holes when the conversion to CSR occurs. Call after interpolation over T and sig0.'''
    thermalMFMTList = [(6,mt) for mt in thermalMTList]
    numGroups = data['numGroups']
    numLegMoments = data['numLegMoments']
    identityForMin = numGroups
    maxColSize = np.zeros(numGroups, dtype=np.int)
    minRowStart = identityForMin * np.ones(numGroups, dtype=np.int)
    maxRowEnd = np.zeros(numGroups, dtype=np.int)     #new
    mfmtList = [(mf,mt) for (mf,mt) in sorted(data['mfmts'], reverse=True) if (mf == 6 and mt != 18 and mt < 221)]
    mfmtFullList = mfmtList + thermalMFMTList
    if not mfmtFullList:
        # If there are no transfer matrices to combine, then thermal scattering is not well
        # defined; set xsTherms as the uncorrected total / scat. cross sections and return
        totalMFMT = (3,1)
        scatMFMT = (3,2)
        data['rxn'][totalMFMT]['xsTherm'] = data['rxn'][totalMFMT]['xsOut']
        data['rxn'][scatMFMT]['xsTherm'] = data['rxn'][scatMFMT]['xsOut']
        return
    # Find union sparsity pattern
    for (mf,mt) in mfmtFullList:
        rxnDict = data['rxn'][(mf,mt)]
        localRowStart = rxnDict['rowStart']
        localColSize = np.diff(rxnDict['indexPtr'])
        localRowEnd = localRowStart + localColSize
        mask = (localRowStart != -1)
        minRowStart[mask] = np.minimum(minRowStart[mask], localRowStart[mask])
        maxRowEnd = np.maximum(maxRowEnd, localRowEnd)
        maxColSize = np.maximum(maxColSize, maxRowEnd - minRowStart)

    # Fill in what would be holes if viewed from CSR perspective
    minRowStart, maxColSize = remove_sparse_holes(minRowStart, maxColSize)
    minRowStart[minRowStart == identityForMin] = -1
    unionIndexPtr = np.zeros(numGroups + 1, dtype=np.int)
    unionIndexPtr[1:] = np.cumsum(maxColSize)
    xsSize = unionIndexPtr[-1]
    mfmtXfer = (6,1)
    data['mfmts'].update([mfmtXfer])
    data['rxn'][mfmtXfer] = {}
    data['rxn'][mfmtXfer]['xsOut'] = np.zeros((numLegMoments, xsSize))
    data['rxn'][mfmtXfer]['rowStart'] = minRowStart
    data['rxn'][mfmtXfer]['indexPtr'] = unionIndexPtr
    data['rxn'][mfmtXfer]['numTherm'] = 1
    data['rxn'][mfmtXfer]['numSig0'] = 1
    data['rxn'][mfmtXfer]['numLegMoments'] = numLegMoments
    # Copy non-thermal and non-fission to union sparsity pattern
    xs = data['rxn'][mfmtXfer]['xsOut']
    for (mf,mt) in mfmtList:
        rxnDict = data['rxn'][(mf,mt)]
        localXS = rxnDict['xsOut']
        localRowStart = rxnDict['rowStart']
        localIndexPtr = rxnDict['indexPtr']
        localColSize = np.diff(localIndexPtr)
        for fromGroup in range(numGroups):
            localStrt, localEnd = localIndexPtr[fromGroup], localIndexPtr[fromGroup + 1]
            unionStrt = unionIndexPtr[fromGroup] + (localRowStart[fromGroup] - minRowStart[fromGroup])
            unionEnd = unionStrt + localColSize[fromGroup]
            xs[:, unionStrt:unionEnd] += localXS[:, localStrt:localEnd]
    # Overwrite thermal section (which may be null)
    thermalGroupStart = numGroups
    for (mf,mt) in thermalMFMTList:
        where = np.where(data['rxn'][(mf,mt)]['rowStart'] != -1)[0]
        localStart = numGroups
        if np.any(where):
            localStart = where[0]
        thermalGroupStart = min(thermalGroupStart, localStart)
    unionStrt = unionIndexPtr[thermalGroupStart]
    unionEnd = unionIndexPtr[-1]
    thermalScatTot = np.zeros(numGroups)
    if thermalMFMTList:
        xs[:, unionStrt:unionEnd] = 0.0
    for (mf, mt) in thermalMFMTList:
        rxnDict = data['rxn'][(mf,mt)]
        localXS = rxnDict['xsOut']
        localRowStart = rxnDict['rowStart']
        localIndexPtr = rxnDict['indexPtr']
        localColSize = np.diff(localIndexPtr)
        for fromGroup in range(thermalGroupStart, numGroups):
            localStrt, localEnd = localIndexPtr[fromGroup], localIndexPtr[fromGroup + 1]
            unionStrt = unionIndexPtr[fromGroup] + (localRowStart[fromGroup] - minRowStart[fromGroup])
            unionEnd = unionStrt + localColSize[fromGroup]
            xs[:, unionStrt:unionEnd] += localXS[:, localStrt:localEnd]
            thermalScatTot[fromGroup] += np.sum(localXS[0, localStrt:localEnd])
    # Add xsTherm to the total xs to be consistent. Requires (3,1) and (3,2) to be present.
    # Overwrite elastic scattering xs, (3,2), with thermal xs so MT1 - MT2 - MT4 = abs in PDT.
    totalMFMT = (3,1)
    scatMFMT = (3,2)
    data['rxn'][totalMFMT]['xsTherm'] = data['rxn'][totalMFMT]['xsOut'].copy()
    data['rxn'][scatMFMT]['xsTherm'] = data['rxn'][scatMFMT]['xsOut'].copy()
    totalXS = data['rxn'][totalMFMT]['xsTherm']
    scatXS = data['rxn'][scatMFMT]['xsTherm']
    if thermalMFMTList:
        totalXS[thermalGroupStart:] += (thermalScatTot[thermalGroupStart:] - scatXS[thermalGroupStart:])
        scatXS[thermalGroupStart:] = thermalScatTot[thermalGroupStart:]

def interpolate_T_sig0_xs(data, desiredT, desiredSig0Vec, outputDict, verbosity=False):
    '''Interpolate over temperature and sigma0 values and store in xsOut variable. The size of desiredSig0Vec should be 1 or numGroups. Return the total XS and flux (useful for Bondarenko iteration). Be careful about combining weights from different background cross sections.'''
    numGroups = data['numGroups']
    numSig0 = data['numSig0']
    numLegMoments = data['numLegMoments']
    if not(hasattr(desiredSig0Vec, '__iter__')):
        desiredSig0Vec = [desiredSig0Vec]
    if len(desiredSig0Vec) != numGroups:
        desiredSig0Vec = desiredSig0Vec[0] * np.ones(numGroups)
    thermList = np.array(sorted(data['thermList']))
    sqrtLambda = lambda x: np.sqrt(x)
    thermIndices, thermFractions = find_interpolation_indices_fractions(desiredT, thermList, sqrtLambda)
    sig0List = np.array(sorted(data['sig0List']))
    sig0FractionMat = find_interpolation_indices_fractions_array(desiredSig0Vec, sig0List, sqrtLambda, True)
    if verbosity:
        print 'Temperature interpolation:', desiredT, thermList[thermIndices], thermFractions
        print 'Sigma0 interpolation:'
        print sorted(sig0List, reverse=True)
        print desiredSig0Vec
    if verbosity > 3:
        for i in range(sig0FractionMat.shape[0]):
            for j in range(sig0FractionMat.shape[1]):
                print '{0:.3f}'.format(sig0FractionMat[i,j]),
            print ''
    # Bound group boundaries
    data['groupBdrs'][-1] = 0.0
    # Use one avg sig0 for all energies to avoid different normalizations at different sig0.
    #avgSig0 = np.exp(np.mean(np.log(desiredSig0Vec)))
    avgSig0 = np.median(desiredSig0Vec)
    avgSig0Index = np.argmin(np.abs(avgSig0 - np.asarray(data['sig0List'])))
    # Condense flux to fluxOut
    data['fluxOut'] = np.zeros(numGroups)
    for thermFrac, thermIndex in zip(thermFractions, thermIndices):
        data['fluxOut'] += thermFrac * data['flux'][thermIndex, avgSig0Index, :]
    outputDict['flux'] = data['fluxOut']
    # Condense xs for each reaction. Do MF 3 first so they're there for the fission spectrum calculation.
    mfmts = sorted(data['mfmts'])
    for (mf, mt) in mfmts:
        rxn = data['rxn'][(mf,mt)]
        if mf == 3:
            rxn['xsOut'] = np.zeros(numGroups)
            localNumTherm = rxn['numTherm']
            localNumSig0 = rxn['numSig0']
            if localNumTherm == 1 and localNumSig0 == 1:
                rxn['xsOut'] = rxn['xs'][0,0,:]
            elif localNumTherm == 1:
                for sig0Index in range(numSig0):
                    mask = sig0FractionMat[sig0Index, :] > 0
                    rxn['xsOut'][mask] += (sig0FractionMat[sig0Index, mask] *
                        rxn['xs'][0, sig0Index, mask])
            elif localNumSig0 == 1:
                for thermFrac, thermIndex in zip(thermFractions, thermIndices):
                    rxn['xsOut'] += thermFrac * rxn['xs'][thermIndex, 0, :]
            else:
                for thermFrac, thermIndex in zip(thermFractions, thermIndices):
                    for sig0Index in range(numSig0):
                        mask = sig0FractionMat[sig0Index, :] > 0
                        rxn['xsOut'][mask] += (thermFrac * sig0FractionMat[sig0Index, mask] *
                            rxn['xs'][thermIndex, sig0Index, mask])
            if mt == 1:
                outputDict['tot'] = rxn['xsOut']
        elif (mf, mt) == (6, 18):
            # Take a flux-weighted sum over all groups
            # Assumes you have (MF,MT) (3,18) (n,fission) and (3,452) (total nubar)
            flux = data['fluxOut']

            fissionXS = data['rxn'][(3,18)]['xsOut']
            nuBar = data['rxn'][(3,452)]['xsOut']
            cutoffGroup = rxn['highestHighEnergyGroup']
            lowEnergySpectrum = rxn['lowEnergySpectrum']
            highEnergyMatrix = rxn['highEnergyMatrix']
            effectiveSpectrum = np.zeros(numGroups)
            effectiveSpectrum += lowEnergySpectrum * np.sum(
                flux[cutoffGroup+1:] * fissionXS[cutoffGroup+1:] * nuBar[cutoffGroup+1:])
            for groupFrom in range(0,cutoffGroup+1):
                effectiveSpectrum += highEnergyMatrix[groupFrom, :] * flux[groupFrom]
            effectiveSpectrum /= np.sum(effectiveSpectrum)
            rxn['xsOut'] = effectiveSpectrum

            FissionMatrix = rxn['FissionMatrix']
            # Compute total fission source = F'*phi
            F_phi = np.dot(flux, FissionMatrix)
            # prompt Chi is normalized total fission source
            rxn['promptChi'] = F_phi/np.sum(F_phi)
            # prompt nu*sig_f (Prod) is row sum of F
            rxn['promptProd'] = np.sum(FissionMatrix, 1)

        elif mf == 6:
            #Convert internal CSC-like format to CSR format for PDT XS
            numVals = len(rxn['xs'][0,0,0,:])
            rxn['xsOut'] = np.zeros((numLegMoments, numVals))
            interpXS = rxn['xsOut']
            indptr = rxn['indexPtr']
            rowStart = rxn['rowStart']
            localNumTherm = rxn['numTherm']
            localNumSig0 = rxn['numSig0']
            rowsPerCol = np.diff(indptr)
            if localNumTherm == 1 and localNumSig0 == 1:
                interpXS[:,:] = rxn['xs'][0,0,:,:]
            elif localNumSig0 == 1:
                for thermFrac, thermIndex in zip(thermFractions, thermIndices):
                    interpXS += thermFrac * rxn['xs'][thermIndex, 0, :, :]
            else:
                for thermFrac, thermIndex in zip(thermFractions, thermIndices):
                    for sig0Index in range(numSig0):
                        for fromGroup in range(numGroups):
                            sig0Frac = sig0FractionMat[sig0Index, fromGroup]
                            if sig0Frac:
                                strt, end = indptr[fromGroup], indptr[fromGroup + 1]
                                interpXS[:, strt:end] += (thermFrac * sig0Frac *
                                    rxn['xs'][thermIndex, sig0Index, :, strt:end])
    if (6, 18) in mfmts and (5, 455) in mfmts:
        # steady-state fission nu and chi is not stored in mfmts
        data['rxn'][(3, 2452)] = {}
        data['rxn'][(3, 2018)] = {}

        flux = data['fluxOut']
        fission_xs = data['rxn'][(3 ,18)]['xsOut']
        fission_x_prompt = data['rxn'][(6, 18)]['FissionMatrix']
        promptProd = data['rxn'][(6, 18)]['promptProd']
        nu_prompt = promptProd/fission_xs
        nu_delayed = data['rxn'][(3 ,455)]['xsOut']
        chis_delayed = data['rxn'][(5 ,455)]['delayedChi']
        chi_delayed = np.sum(chis_delayed, axis=0)

        nu_ss = nu_prompt + nu_delayed
        n_per_gout = ( np.dot(flux, fission_x_prompt) + \
                   chi_delayed*np.sum(nu_delayed*fission_xs*flux) )
        chi_ss = n_per_gout/np.sum(n_per_gout)

        data['rxn'][(3, 2452)]['xsOut'] = nu_ss
        data['rxn'][(3, 2018)]['xsOut'] = chi_ss



def convert_to_scipy_sparse_scat_matrix(data, format='csr'):
    '''Convert from my internal CSC-like sparse format to an official Scipy format (either CSR or CSC). We have filled in holes so the conversion will be dense on each row.'''
    numGroups = data['numGroups']
    numLegMoments = data['numLegMoments']
    mfmts = [(mf,mt) for (mf,mt) in data['mfmts'] if (mf == 6 and mt != 18)]
    for (mf, mt) in mfmts:
        rxn = data['rxn'][(mf,mt)]
        rowStart = rxn['rowStart']
        indptr = rxn['indexPtr']
        rowsPerCol = np.diff(indptr)

        rowStart[rowStart == -1] = numGroups
        rowStart_holess = numGroups*np.ones(len(rowStart))
        rowsPerCol_holess = np.zeros(len(rowsPerCol), dtype=np.int)
        rowStart_holess, rowsPerCol_holess = remove_sparse_holes(rowStart, rowsPerCol)
        rowStart_holess[rowStart_holess == numGroups] = -1
        rowStart[rowStart == numGroups] = -1

        indptr_holess = np.zeros(numGroups+1, dtype=np.int)
        indptr_holess[1:] = np.cumsum(rowsPerCol_holess)
        # for i in range(len(indptr)):
        #     assert indptr[i] == indptr_holess[i], "indptr[{0}]={1}, indptr_holess[{0}]={2}".format(i, indptr[i], indptr_holess[i])

        # numVals = len(rxn['xsOut'][0,:])
        numVals = indptr_holess[-1]
        vals_holess = np.zeros((numLegMoments, numVals))
        indices_holess = np.zeros(numVals, dtype=np.int)
        for col in range(numGroups):
            assert (indptr_holess[col+1] - indptr_holess[col]) == rowsPerCol_holess[col]
            indices_holess[indptr_holess[col]:indptr_holess[col+1]] = np.arange(rowsPerCol_holess[col]) + rowStart_holess[col]
        for i in range(len(indices_holess)):
            assert indices_holess[i] >= -1 and indices_holess[i] < numGroups, "indices_holess[{0}] = {1}".format(i, indices_holess[i])


        for moment in range(numLegMoments):
            vals = rxn['xsOut'][moment, :]
            for fromGroup in range(numGroups):
                Strt, End = indptr[fromGroup], indptr[fromGroup + 1]
                holessStrt = indptr_holess[fromGroup] + (rowStart[fromGroup] - rowStart_holess[fromGroup])
                holessEnd = holessStrt + rowsPerCol[fromGroup]
                vals_holess[moment, holessStrt:holessEnd] = vals[Strt:End]



        xs = []
        for moment in range(numLegMoments):
#            vals = rxn['xsOut'][moment, :]
            if format.lower().strip() == 'csc':
                xs.append(sparse.csc_matrix(
                    (vals_holess[moment, :], indices_holess, indptr_holess), shape=(numGroups, numGroups)))
            elif format.lower().strip() == 'csr':
                xs.append(sparse.csc_matrix(
                    (vals_holess[moment, :], indices_holess, indptr_holess), shape=(numGroups, numGroups)).tocsr())
        rxn['xsOut'] = xs

def write_pdt_xs(filePath, data, temperature, format='csr', whichXS='usual', fromFile='barnfire'):
    '''Write a PDT XS from a xs dat. CSR/CSC format will print the matrices by row/column (incident/exiting group). whichXS determines what is written. "all" means all xs are written. "total" means only the flux and total cross section are written. Anything else means the normal PDT cross sections are written.'''
    timeStr = datetime.strftime(datetime.now(), '%c')
    mfmts = data['mfmts']
    numGroups = data['numGroups']
    numLegMoments = data['numLegMoments']
    xsType = 'multigroup' # Not always true
    groupBoundaries = data['groupBdrs']
    if whichXS.lower().strip() in ['separate', 'all']:
        numXS = np.sum([1 for (mf, mt) in mfmts if mf == 3])
        if whichXS.lower().strip() in 'separate':
            # separate - Include all XS except the combined scattering matrix
            numXfer = np.sum([1 for (mf, mt) in mfmts if (mf == 6 and mt != 1)])
        else:
            # all - Include all XS
            numXfer = np.sum([1 for (mf, mt) in mfmts if (mf == 6)])
        # The flux counts as a cross section
        numXS += 1
        if (6,18) in mfmts:
            # Add in cross sections for prompt fission spectrum, prod. Fission matrix is store as MT 2518
            numXS += 2
            numXfer += 1
        if (5,455) in mfmts:
            # Add decayConst and delayedChi for delay fission neutrons
            # Both are regarded as 1D xs (MF 3)
            numXS += 2
        if (6,18) in mfmts and (5,455) in mfmts:
            # If has both prompt and delayed neutron information, compute and add steady-state nu and chi
            numXS += 2
    elif whichXS.lower().strip() in 'total':
        #Include flux and total XS
        numXS = 2
        numXfer = 0
    else:
        # Include flux, total XS, and combined scattering matrix by default. Fission matrix is also included
        numXS = 2
        numXfer = 0
        if (6,1) in mfmts:
            numXfer += 1
        if (6,18) in mfmts:
            numXS += 2
            numXfer += 1
        if (3,18) in mfmts:
            # add MT 18
            numXS += 1
        if (3, 455) in mfmts:
            # add MT 455 (delayed nu)
            numXS += 1
        if (5,455) in mfmts:
            # Add decayConst and delayedChi for delay fission neutrons
            # Both are regarded as 1D xs (MF 3)
            numXS += 2
        if (6,18) in mfmts and (5,455) in mfmts:
            # If has both prompt and delayed neutron information, compute and add steady-state nu and chi
            numXS += 2
            numXfer += 1 # Fission matrix (ATT: already included?)
        if (3,259) in mfmts:
            # Include the inverse velocity XS, if available
            numXS += 1
        if whichXS.lower().strip() in 'abs':
            if (3,2) in mfmts:
                # Include the (1-D) elastic scattering XS, if available
                numXS += 1
            if (3,4) in mfmts:
                # Include the (1-D) inelastic scattering XS, if available
                numXS += 1
    # Get plurality of process[es]
    pluralMF3Str = ''
    if numXS != 1:
        pluralMF3Str = 'es'
    pluralMF6Str = ''
    if numXfer != 1:
        pluralMF6Str = 'es'
    with open(filePath, 'w') as fid:
        fid.write('PDT Format Material Data File created {0}\n'.format(timeStr))
        fid.write('\n')
        fid.write('This file is a {0} neutron library generated from {1}.\n'.format(xsType, fromFile))
        fid.write('1 temperatures, 1 densities, and {0} groups.\n'.format(numGroups))
        fid.write('\n')
        fid.write('{0:g} neutron process{1} and {2:g} transfer process{3}.\n'.format(numXS, pluralMF3Str, numXfer, pluralMF6Str))
        fid.write('Scattering order {0:g}\n'.format(numLegMoments-1))
        fid.write('\n')
        fid.write('Microscopic cross sections are in units of barns.\n')
        fid.write('\n')
        fid.write('Temperatures in Kelvin:\n')
        fid.write('{0:>15g}\n'.format(temperature))
        fid.write('\n')
        fid.write('Densities in g/cc:\n')
        fid.write('{0:>15}\n'.format(0))
        fid.write('\n')
        fid.write('Group boundaries in eV:\n')
        fid.write(multiline_string(groupBoundaries, 15, 5, 7))
        fid.write('\n')
        fid.write('T = {0:g} density = 0\n'.format(temperature))
        fid.write('---------------------------------------------------\n')
        # Write flux
        pdtMT = 1099
        xs = data['fluxOut']
        fid.write('MT {0}\n'.format(pdtMT))
        fid.write(multiline_string(xs, 20, 5, 12))
        # Write total
        mf, mt = (3,1)
        if whichXS.lower().strip() in 'separate':
            xs = data['rxn'][(mf,mt)]['xsOut']
        else:
            xs = data['rxn'][(mf,mt)]['xsTherm']
        pdtMT = get_pdt_mt(mf, mt)
        fid.write('MT {0}\n'.format(pdtMT))
        fid.write(multiline_string(xs, 20, 5, 12))
        if whichXS.lower().strip() in 'total':
            fid.write('\n')
            return
        # Write renormalized elastic scattering XS (vector) if desired
        if whichXS.lower().strip() in 'abs':
            mf, mt = (3,2)
            xs = data['rxn'][(mf,mt)]['xsTherm']
            pdtMT = get_pdt_mt(mf, mt)
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write(multiline_string(xs, 20, 5, 12))
        # Write XS (vectors)
        mtsForMF3 = []
        if (3,18) in mfmts:
            mtsForMF3 = [18]
        if (3,455) in mfmts:
            mtsForMF3.append(455)
        if whichXS.lower().strip() in ['separate', 'all']:
            mtsForMF3 = [mt for (mf,mt) in sorted(mfmts) if (mf == 3 and mt != 1)]
        if (3,259) in mfmts:
            mtsForMF3.append(259)
        if whichXS.lower().strip() in 'abs':
            # The elastic scattering XS has already been printed, so don't print twice
            if (3,4) in mfmts:
                mtsForMF3.append(4)
        mf = 3
        for mt in mtsForMF3:
            pdtMT = get_pdt_mt(mf, mt)
            xs = data['rxn'][(mf,mt)]['xsOut']
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write(multiline_string(xs, 20, 5, 12))
        # Write fission spectrum and nu * sig_f
        if (6,18) in mfmts:
            # write effective prompt fission spectrum (Chi)
            pdtMT = 1018
              #  xs = data['rxn'][(6,18)]['xsOut']
            xs = data['rxn'][(6,18)]['promptChi']
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write(multiline_string(xs, 20, 5, 12))
            # write effective prompt fission nu*sif_f (Prod)
            pdtMT = 1452
              # xs = data['rxn'][(3,18)]['xsOut'] * data['rxn'][(3,452)]['xsOut']
            xs = data['rxn'][(6,18)]['promptProd']
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write(multiline_string(xs, 20, 5, 12))
            # write total prompt fission transfer matrix (chi dot nusigf)
            pdtMT = 2518
            xs = data['rxn'][(6,18)]['FissionMatrix']
            fid.write('MT {0}, Moment {1}\n'.format(pdtMT, 0))
            for g in range(numGroups):
                fid.write('{0}, first, last: '.format('  Sink'))
                first = 0
                last = numGroups - 1
                vec = [g, first, last]
                fid.write(multiline_string(vec, 5, 3, 10))
                fid.write(multiline_string(xs[:, g], 20, 5, 12))
        # if (3,18) in mfmts and (3,452) in mfmts:
        #     pdtMT = 1452
        #     xs = data['rxn'][(3,18)]['xsOut'] * data['rxn'][(3,452)]['xsOut']
        #     fid.write('MT {0}\n'.format(pdtMT))
        #     fid.write(multiline_string(xs, 20, 5, 12))
        if (5,455) in mfmts:
            pdtMT = 1054
            numDNGs = data['rxn'][(5,455)]['numDNGs']
            decayConst = data['rxn'][(5,455)]['decayConst']
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write('  Number of delayed neutron groups: {0}\n'.format(numDNGs))
            fid.write(multiline_string(decayConst, 20, 5, 12))
            pdtMT = 2055
            delayedChi = data['rxn'][(5,455)]['delayedChi']
            fid.write('MT {0}\n'.format(pdtMT))
            fid.write('  Number of delayed neutron groups: {0}\n'.format(numDNGs))
            for iDNG in range(numDNGs):
                fid.write('  DNG {0}\n'.format(iDNG))
                fid.write(multiline_string(delayedChi[iDNG,:], 20, 5, 12))
        # Write steady-state nu and chi
        if (6, 18) in mfmts and (5, 455) in mfmts:
            pdtMT = 2452
            fid.write('MT {0}\n'.format(pdtMT))
            nu_ss = data['rxn'][(3,2452)]['xsOut']
            fid.write(multiline_string(nu_ss, 20, 5, 12))
            pdtMT = 2018
            fid.write('MT {0}\n'.format(pdtMT))
            chi_ss = data['rxn'][(3,2018)]['xsOut']
            fid.write(multiline_string(chi_ss, 20, 5, 12))
        # Write transfer matrices (sparse matrices)
        if whichXS.lower().strip() in 'separate':
            mtsForMF6 = [mt for (mf,mt) in sorted(mfmts) if (mf == 6 and mt != 1)]
        elif whichXS.lower().strip() in 'all': # ATT ???
            mtsForMF6 = [mt for (mf,mt) in sorted(mfmts) if (mf == 6 and mt != 18)]
        else:
            mtsForMF6 = [mt for (mf,mt) in sorted(mfmts) if (mf == 6 and mt == 1)]
        mf = 6
        for mt in mtsForMF6:
            pdtMT = get_pdt_mt(mf, mt)
            rxn = data['rxn'][(mf,mt)]
            if format.lower().strip() == 'csr':
                description = '  Sink'
            elif format.lower().strip() == 'csc':
                description = 'Source'
            for moment in range(numLegMoments):
                xs = rxn['xsOut'][moment].data
                indexPtr = rxn['xsOut'][moment].indptr
                indices = rxn['xsOut'][moment].indices
                fid.write('MT {0}, Moment {1}\n'.format(pdtMT, moment))
                for g in range(numGroups):
                    strt, end = indexPtr[g], indexPtr[g+1]
                    fid.write('{0}, first, last: '.format(description))
                    if strt == end:
                        vec = [g, -1, -2]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write('\n')
                    else:
                        first = indices[strt]
                        last = indices[end-1]
                        for i in range(strt+1,end):
                            index_prev = indices[i-1]
                            index = indices[i]
                            assert (index - index_prev) == 1, "index i-1: {0},  index i: {1}".format(index_prev, index)
                        assert (last - first + 1) == len(xs[strt:end]), "indices len {0} <--> data len {1}".format((last - first + 1), len(xs[strt:end]))
                        vec = [g, first, last]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write(multiline_string(xs[strt:end], 20, 5, 12))
        fid.write('\n')

def print_csc():
        for mt in mtsForMF6:
            pdtMT = get_pdt_mt(mf, mt)
            rxn = data['rxn'][(mf,mt)]
            xs = rxn['xsOut']
            indexPtr = rxn['indexPtr']
            rowStart = rxn['rowStart']
            for moment in range(numLegMoments):
                fid.write('MT {0}, Moment {1}\n'.format(pdtMT, moment))
                for gFrom in range(numGroups):
                    strt, end = indexPtr[gFrom], indexPtr[gFrom+1]
                    sz = end - strt
                    firstTo = rowStart[gFrom]
                    lastTo = firstTo + sz - 1
                    fid.write('Source, first, last: ')
                    if sz != 0:
                        vec = [gFrom, firstTo, lastTo]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write(multiline_string(xs[moment, strt:end], 20, 5, 12))
                    else:
                        vec = [gFrom, -1, -2]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write('\n')

####################################################################################
def finish_parsing(inputDict):
    '''Specify these in the calling program'''
    if 'energyMesh' not in inputDict:
        inputDict['energyMesh'] = None
    if 'flux' not in inputDict:
        inputDict['flux'] = None
    if 'sig0Vec' not in inputDict:
        inputDict['sig0Vec'] = None
    if 'zeroDmtDict' not in inputDict:
        inputDict['zeroDmtDict'] = None
    inputDict['finishedParsing'] = True

def define_input_parser():
    import argparse
    from directories import get_common_directories
    #
    dirDict = get_common_directories()
    thisDirr = os.path.abspath('.')
    pyDirr = dirDict['src']
    gendfDirr = dirDict['gendf']
    # If run from the Pyscripts directory, default xs directory is $gendf/pu-239. Otherwise, default is '.'
    if thisDirr == pyDirr:
        xsDefaultDirr = os.path.join(gendfDirr, 'pu-239')
    else:
        xsDefaultDirr = thisDirr
    xsInDirr = xsDefaultDirr
    xsOutDirr = xsDefaultDirr
    parser = argparse.ArgumentParser(description='GENDF to PDT XS file format converter.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=0, choices=[0,1,2,3,4], type=int)
    parser.add_argument('-w', '--workopt', help='What to do. gendf means read in GENDF file and write pdt xs file and pickle file. pickle means read in pickled xs file. rxn means read reaction list from GENDF file. scat means read scattering sizes from GENDF file. pot means read potential scattering from ENDF file. pen means read total cross section from pendf file. decay means read in decay data from endf file (mf8, mt457). fp means read in fission product yield from endf file (mf8, mt454/459). meta means read in metastable branching ratios or cross sections (mf8/9/10). qvalue means read in the Q value for the reactions given in mf3list.', choices=['gendf', 'pickle', 'rxn', 'scat', 'pen', 'pot', 'decay', 'fp', 'meta', 'qvalue'], default='gendf')
    #res means read resolved resonance region from ENDF file (not implemented).
    parser.add_argument('-T', '--temperature', help='Temperature in Kelvin.', type=float, default=0.)
    parser.add_argument('-Z', '--sig0', help='Background XS in barns.', type=float, default=1e10)
    parser.add_argument('-i', '--inputdir', help='Input directory with GENDF (or ENDF / PENDF / pickle) xs data. Default does not change based on other options.', default=xsInDirr)
    parser.add_argument('-o', '--outputdir', help='Output directory', default=xsOutDirr)
    parser.add_argument('-I', '--inname', help='Name of input GENDF (or ENDF / PENDF / pickle) file. Default does not change based on other options..', default='gendf')
    parser.add_argument('-P', '--picklename', help="Name for output/input pickle. Do not print if 'none'.", default='xs.p')
    parser.add_argument('-O', '--outname', help='Name for output PDT XS.', default='pdt_xs.data')
    parser.add_argument('-m', '--mf3list', help='List of MTs to keep for MF 3. MT of 0 means keep all available. Use ENDF numbering.', nargs='*', type=int, default=[0])
    parser.add_argument('-n', '--mf5list', help='List of MTs to keep for MF 5. MT of 0 means keep all available. Use ENDF numbering.', nargs='*', type=int, default=[0])
    parser.add_argument('-M', '--mf6list', help='List of MTs to keep for MF 6. MT of 0 means keep all available. Use ENDF numbering.', nargs='*', type=int, default=[0])
    parser.add_argument('-t', '--thermallist', help='List of thermal MTs to keep. Use unofficial ENDF numbering.', type=int, nargs='*', default=[])
    parser.add_argument('-f', '--format', help='Output format for scattering matrices column or row major).', choices=['csr', 'csc'], default='csr')
    parser.add_argument('-p', '--printopt', help='Which XS to print. Usual prints selected XS and the combined transfer matrix; abs prints the one-dimensional elastic and inelastic scattering cross sections so PDT can correctly calculate the net absorption cross section; total just prints the flux and total XS; separate prints all XS and transfer matrices except the combined transfer matrix; all prints everything; none prints no PDT XS.', choices=['none', 'total', 'usual', 'abs', 'separate', 'all'], default='usual')
    return parser

####################################################################################
if __name__=='__main__':
    parser = define_input_parser()
    inputDict = vars(parser.parse_args())
    finish_parsing(inputDict)
    if inputDict['verbosity'] > 1:
        print 'Summary of inputs:'
        print inputDict
    execute_reader(inputDict)
