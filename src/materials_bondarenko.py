'''
Andrew Till
Summer 2014

Bondarenko iteration utility for materials
'''

#STDLIB
import os
import shutil
#TPL
import numpy as np
#MINE
from materials_util import is_fissionable
import materials_util as util
from directories import get_common_directories
import Readgroupr as readgroupr
import PDTXS as pdtxs

def perform_bondarenko_iterations(inputDict, materials, verbosity):
    '''Driver for Bondarenko iterations'''
    maxIterations = inputDict['numberiterationsmax']
    maxError = inputDict['errormax']
    useSimpleRxn = inputDict['simplereactions']
    scatMatrixFormat = inputDict['format']
    rxnsToPrint = inputDict['printopt']
    energyMeshPath = inputDict['mesh']
    fluxBasePath = inputDict['fluxes']
    #
    dirDict = get_common_directories()
    rootDirr = dirDict['gendf']
    outDirr = dirDict['pdtxs']
    #
    energyMesh = None
    if energyMeshPath is not None:
        energyMesh = np.loadtxt(energyMeshPath, dtype=np.int, skiprows=2, usecols=[0])[:-1]
        numElements = len(np.unique(energyMesh))
        energyMeshFilenameOut = 'mesh_{0}.txt'.format(numElements)
        energyMeshPathOut = os.path.join(outDirr, energyMeshFilenameOut)
        shutil.copy2(energyMeshPath, energyMeshPathOut)
    #
    if verbosity:
        print '------- Bondarenko -------'
    fluxDict = read_fluxes(fluxBasePath, materials)
    for material in materials:
        backgroundXSDict = iterate_one_material(rootDirr, material, maxError, maxIterations, energyMesh, fluxDict, verbosity)
        if maxIterations < 0:
            unset_background_xs_dict(material, backgroundXSDict, verbosity)
        print_one_material(rootDirr, outDirr, material, backgroundXSDict, scatMatrixFormat, useSimpleRxn, rxnsToPrint, energyMesh, fluxDict, verbosity)
    if verbosity:
        key = backgroundXSDict.keys()[0]
        numGroups = len(backgroundXSDict[key])
        print 'Number of groups is', numGroups

def read_fluxes(fluxBasePath, materials):
    '''Read in flux files'''
    fluxDict = {}
    if fluxBasePath is None:
        return None
    for material in materials:
        shortName = material.shortName
        fluxPath = fluxBasePath.format(m=shortName)
        fluxDict[shortName] = np.loadtxt(fluxPath, skiprows=1, usecols=[1])
    return fluxDict

def print_one_material(rootDirr, outDirr, material, backgroundXSDict, scatMatrixFormat, useSimpleRxn, rxnsToPrint, energyMesh, fluxDict, verbosity):
    '''Print PDT XS for one material. Prints both the component-wise and combined material xs's. Requires unique shortName for each material to prevent over-writing.'''
    txs2mtDict = readgroupr.get_short2mt_dict(readgroupr.get_endf_mt_list())
    T = material.temperature
    ZAList = sorted(material.ZAList)
    for (Z,A) in ZAList:
        sig0Vec = backgroundXSDict[(Z,A)]
        numGroups = len(sig0Vec)
        Sab = material.SabDict[(Z,A)]
        sym = material.symDict[Z]
        shortName = material.shortName
        # Metastable isomeric states use the groundstate A + 400
        effA = A % 400
        metastableStr = ''
        if A // 400 > 0:
            metastableStr = 'm'
        leafDirr = util.get_nuclide_dirr(sym, effA, Sab, metastableStr) 
        inDirr = os.path.join(rootDirr, leafDirr)
        readerOpt = 'gendf'
        outName = 'xs_{0}_{1}-{2}_{3}.data'.format(shortName, sym.lower(), A, numGroups)
        pickleName = None
        thermalMTList = ['{0}'.format(txs2mtDict[txs]) for txs in material.thermalXSDict[(Z,A)]]
        thermalMTStr = ' '.join(thermalMTList)
        parser = readgroupr.define_input_parser()
        parseStr = '-i {i} -o {o} -O {O} -P {P} -w {w} -p {p} -t {t} -T {T} -f {f}'.format(
            i=inDirr, o=outDirr, O=outName, P=pickleName, p=rxnsToPrint, w=readerOpt, t=thermalMTStr, T=T, f=scatMatrixFormat)
        if useSimpleRxn:
            parseStr += ' -m 1 2 18 102 221 452 -M 2 18 221 -t 221'
        if verbosity > 2:
            print 'Calling ./Readgroupr', parseStr
        if verbosity:
            print 'Printing XS to {0}'.format(os.path.join(outDirr, outName))
        readerDict = vars(parser.parse_args(parseStr.split()))
        if fluxDict is not None:
            readerDict['flux'] = fluxDict[material.shortName]
        readerDict['energyMesh'] = energyMesh
        readerDict['sig0Vec'] = sig0Vec
        readgroupr.finish_parsing(readerDict)
        readgroupr.execute_reader(readerDict)
        if verbosity > 2:
            plot_bondarenko(rootDirr, backgroundXSDict)
    form_and_print_macroscopic_xs(outDirr, ZAList, material, numGroups, verbosity)

def form_and_print_macroscopic_xs(dirr, ZAList, material, numGroups, verbosity=False):
    '''Combine all microscopic component XS into one macroscopic material XS'''
    shortName = material.shortName
    MTinvel = 259
    MTdecay = 457
    MTfissEnergy = 458
    MTwgt = 1099
    MTchi = 1018
    MTnuSigF = 1452
    # Read in component cross sections
    xsDictIn = {}
    for (Z,A) in ZAList:
        sym = material.symDict[Z]
        inName = 'xs_{0}_{1}-{2}_{3}.data'.format(shortName, sym.lower(), A, numGroups)
        inPath = os.path.join(dirr, inName)
        xsDictIn[(Z,A)] = pdtxs.read_PDT_xs_generally(inPath)
    # Initialize material cross section dictionary (xsOut)
    key = (Z,A)
    t = xsDictIn[key]
    microStr = 'Macroscopic cross sections are in units of cm^-1.'
    xsDict = pdtxs.PDT_XS(t.G, t.M, t.T, t.typeStr, microStr, t.Eg, t.dE, {})
    xsOut = xsDict.xs
    # Keep a reaction if it appears in at least one component
    MTs = set()
    for (Z,A) in xsDictIn:
        MTs.update(xsDictIn[(Z,A)].xs.keys())
    # A material-average decay rate does not make sense, so do not compute one
    if MTdecay in MTs:
        MTs.remove(MTdecay)
    # Initialize 0D XS
    MTs0D = [MT for MT in MTs if MT in [MTdecay, MTfissEnergy]]
    for MT in MTs0D:
        xsOut[MT] = 0.
    # Initialize 1D XS
    MTs1D = [MT for MT in MTs if (MT < 2000 and MT not in MTs0D)]
    for MT in MTs1D:
        xsOut[MT] = np.zeros(xsDict.G)
    # Initialize transfer matrices
    MTsXfer = [MT for MT in MTs if MT >= 2000]
    for MT in MTsXfer:
        xsOut[MT] = np.zeros((xsDict.M, xsDict.G, xsDict.G))
    # Save denominators for XS that are averages instead of density-weighted sums
    MTsAvg = [MT for MT in MTs if MT in [MTwgt, MTchi, MTfissEnergy, MTinvel]]
    norms = {}
    for MT in MTsAvg:
        norms[MT] = 0.
    # Compute XS averages and sums by adding in the contribution of each component
    for (Z,A) in ZAList:
        xsIn = xsDictIn[(Z,A)].xs
        compDensity = material.atomDensity * material.elemAtomFracDict[Z] * material.abundanceDict[(Z,A)]

        # Compute flux weight and its sum for this component
        wgt = 0.
        wgtSum = 1.
        if MTwgt in xsIn:
            wgt = xsIn[MTwgt]
            wgtSum = np.sum(wgt)
            if not wgtSum:
                wgtSum = 1.
            xsOut[MTwgt] += compDensity * wgt / wgtSum
            norms[MTwgt] += compDensity
        # Compute fission rate sum for this component
        fissRate = 0.
        if MTnuSigF in xsIn:
            fissRate = np.sum(xsIn[MTnuSigF] * wgt) / wgtSum
        # Update numerator and denominator for energy per fission using fission-source weighting
        if MTfissEnergy in xsIn:
            xsOut[MTfissEnergy] += compDensity * fissRate * xsIn[MTfissEnergy]
            norms[MTfissEnergy] += compDensity * fissRate
        # Update numerator and denominator for chi using fission-source weighting
        if MTchi in xsIn:
            xsOut[MTchi] += compDensity * fissRate * xsIn[MTchi]
            norms[MTchi] += compDensity * fissRate
        # Update numerator and denominator for inverse velocity using density weighting
        if MTinvel in xsIn:
            xsOut[MTinvel] += compDensity * xsIn[MTinvel]
            norms[MTinvel] += compDensity

        # Compute cross sections that are density-weighted sums
        MTsSum = [MT for MT in MTs if MT not in MTsAvg]
        for MT in MTsSum:
            if MT in xsIn:
                xsOut[MT] += compDensity * xsIn[MT]

    # Normalize XS averages
    for MT in MTsAvg:
        if norms[MT]:
            xsOut[MT] /= norms[MT]

    # Print out material XS
    outName = 'xs_{0}_{1}.data'.format(shortName, numGroups)
    outPath = os.path.join(dirr, outName)
    pdtxs.write_PDT_xs_generally(outPath, xsDict)
    if verbosity:
        print 'Printing combined XS to {0}'.format(outPath)

def iterate_one_material(rootDirr, material, maxError, maxIterations, energyMesh=None, fluxDict=None, verbosity=False):
    '''Perform Bondarenko iteration on one material. Fine groups within an energy element share the same background cross section.'''
    sig0Vec = None
    if verbosity:
        print 'Performing Bondarenko iteration for material {0}'.format(material.longName)
    ZAList = sorted(material.ZAList)
    readerOpt = 'gendf'
    totalXSDict = {}
    backgroundXSDict = {}
    iterationCount = 0
    globalError = 1.0
    for (Z,A) in ZAList:
        read_one_total_xs(rootDirr, Z, A, material, totalXSDict, readerOpt, sig0Vec, energyMesh, fluxDict, verbosity)
    build_all_background_xs(material, totalXSDict, backgroundXSDict, verbosity)
    print_bondarenko(iterationCount, maxIterations, globalError, maxError, verbosity)
    readerOpt = 'pickle'
    while globalError > maxError and iterationCount < maxIterations:
        globalError = 0.0
        for (Z,A) in ZAList:
            sig0Vec = backgroundXSDict[(Z,A)]
            localError = read_one_total_xs(rootDirr, Z, A, material, totalXSDict, readerOpt, sig0Vec, energyMesh, fluxDict, verbosity)
            globalError = max(localError, globalError)
        build_all_background_xs(material, totalXSDict, backgroundXSDict, verbosity)
        iterationCount += 1
        print_bondarenko(iterationCount, maxIterations, globalError, maxError, verbosity)
    return backgroundXSDict

def unset_background_xs_dict(material, backgroundXSDict, verbosity):
    for key in backgroundXSDict:
        size = len(backgroundXSDict[key])
        sig0 = material.backgroundXSDict[key]
        if sig0 == np.inf:
            sig0 = 1e10
        backgroundXSDict[key] = sig0 * np.ones(size)

def read_one_total_xs(rootDirr, Z, A, material, totalXSDict, readerOpt='gendf', sig0Vec=None, energyMesh=None, fluxDict=None, verbosity=False):
    '''Read the total XS for one nuclide for one material'''
    T = material.temperature
    Sab = material.SabDict[(Z,A)]
    sym = material.symDict[Z]
    sig0 = material.backgroundXSDict[(Z,A)]
    if sig0 == np.inf:
        sig0 = 1e10
    # Metastable isomeric states use the groundstate A + 400
    effA = A % 400
    metastableStr = ''
    if A // 400 > 0:
        metastableStr = 'm'
    #
    leafDirr = util.get_nuclide_dirr(sym, effA, Sab, metastableStr) 
    fullDirr = os.path.join(rootDirr, leafDirr)
    parser = readgroupr.define_input_parser()
    parseStr = '-i {i} -o {o} -w {w} -p none -m 1 -M -t -T {T} -Z {Z}'.format(
        i=fullDirr, o=fullDirr, w=readerOpt, T=T, Z=sig0)
    if verbosity > 2:
        print 'Calling ./Readgroupr', parseStr
    readerDict = vars(parser.parse_args(parseStr.split()))
    if fluxDict is not None:
        readerDict['flux'] = fluxDict[material.shortName]
    readerDict['energyMesh'] = energyMesh
    readerDict['sig0Vec'] = sig0Vec
    readgroupr.finish_parsing(readerDict)
    xsDict = readgroupr.execute_reader(readerDict)
    totalXS = xsDict['tot']
    epsilon = 1E-11
    if (Z,A) in totalXSDict:
        err = np.linalg.norm((totalXSDict[(Z,A)] - totalXS) / (totalXS + epsilon), np.inf) / np.sqrt(len(totalXS))
        if verbosity > 1:
            print '>>> Error for ({0}, {1}) is {2}'.format(Z, A, err)
        if verbosity > 3:
            print (totalXSDict[(Z,A)] - totalXS) / (totalXS + epsilon)
    else:
        err = 1.0
    totalXSDict[(Z,A)] = totalXS
    return err

def build_all_background_xs(material, totalXSDict, backgroundXSDict, verbosity=False):
    '''Build the background XS for each nuclide in one material.'''
    ZAList = sorted(material.ZAList)
    numGroups = len(totalXSDict[ZAList[0]])
    numNuclides = len(totalXSDict)
    backgroundXSs = np.zeros((numNuclides, numGroups))
    # Build the background XS
    for i, (Z,A) in enumerate(ZAList):
        mask = np.ones(numNuclides, dtype=np.bool)
        mask[i] = False
        atomDensity = material.atomDensity * material.elemAtomFracDict[Z] * material.abundanceDict[(Z,A)]
        # Uses broadcasting
        backgroundXSs[mask, :] += atomDensity * totalXSDict[(Z,A)]
    backgroundXSs += material.chordLength
    for i, (Z,A) in enumerate(ZAList):
        atomDensity = material.atomDensity * material.elemAtomFracDict[Z] * material.abundanceDict[(Z,A)]
        backgroundXSs[i, :] /= atomDensity
        # Uses aliasing
        backgroundXSDict[(Z,A)] = backgroundXSs[i, :]

def print_bondarenko(iterationCount, maxIterations, error, maxError, verbosity):
    if verbosity:
        print 'Iteration {0:2g} (max {1})'.format(iterationCount, maxIterations),
        print 'Error {0:.2e} (max {1})'.format(error, maxError)

def plot_bondarenko(rootDirr, backgroundXSDict):
    for (Z,A) in backgroundXSDict.keys():
        from matplotlib import pyplot as plt
        plt.clf()
        y = backgroundXSDict[(Z,A)]
        x = np.arange(len(y))
        plt.semilogy(x,y)
        plt.xlabel('Group index')
        plt.ylabel('Background xs (b)')
        nom = 'sig0_{0}_{1}.pdf'.format(Z,A)
        path = os.path.join(rootDirr, nom)
        plt.savefig(path)



