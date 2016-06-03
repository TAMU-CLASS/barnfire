'''
Andrew Till
Summer 2014

Write NJOY script
'''

#STDLIB
import os
from datetime import datetime
#TPL
#import numpy as np
#MINE
from directories import try_mkdir, try_make_exe, get_common_directories

def create_thermal_ace_njoy_script(dat, tapes):
    '''Create NJOY file to generate a thermal ACE output. These are done for each thermal treatment and each temperature.'''


    #Directories
    dirDict = get_common_directories()
    pyDirr = dirDict['src']
    endfDirr = dirDict['endf']
    pendfRootDirr = dirDict['pendf']
    aceRootDirr = dirDict['ace']
    njoyPath = os.path.join(dirDict['njoyInstall'], 'xnjoy')
    #
    # TODO: CHANGE nuclideName
    pendfDirr = os.path.join(pendfRootDirr, str(dat.nuclideName))
    aceDirr = os.path.join(aceRootDirr, str(dat.thermalFileName))
    #
    pendfFile = 'pendf'
    endfPath = os.path.join(endfDirr, dat.endfFile)
    pendfPath = os.path.join(pendfDirr, pendfFile)
    #
    relPathNJOYAce = os.path.relpath(njoyPath, aceDirr)
    relPathEndfAce = os.path.relpath(endfPath, aceDirr)
    relPathPendfAce = os.path.relpath(pendfPath, aceDirr)
    #
    numThermals = len(dat.endfThermalFileList)
    numTemperatures = len(dat.thermList)
    pendfOutTape = get_pendf_out_tape(tapes, numThermals)
    endfThermalPathList = [os.path.join(endfDirr, thermalFile) for
        thermalFile in dat.endfThermalFileList if thermalFile]
    relPathThermalAceList = [os.path.relpath(thermalPath, aceDirr) for
        thermalPath in endfThermalPathList]
    endfThermalTapeList = map(lambda x: x + tapes.endfThermalStart, range(numThermals))
    endfNonFreeThermalTapeList = [tape for (mt,tape) in
        zip(dat.inelasticMTList, endfThermalTapeList) if mt != 221]
    #
    aceScriptFile = 'runNJOY.sh'
    aceScriptPath = os.path.join(aceDirr, aceScriptFile)
    # TODO: Generalize to use different extension root (instead of 2)
    aceFileTemplate = '{H}.{e}{s}t'

    # Script to run NJOY from PENDF to thermal ACE:
    deck = []
    deck.append(["#! /usr/bin/env bash"])
    deck.append(["echo 'NJOY Problem {0} {1} (PENDF to thermal ACE)'".format(dat.thermalNuclideName, dat.endfName)])
    deck.append(["echo 'Getting ENDF input tape, PENDF input tape and NJOY executable'"])
    deck.append(['ln -fs {0} xnjoy'.format(relPathNJOYAce)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathEndfAce, tapes.endf)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathPendfAce, abs(pendfOutTape))])
    for endfTape, relPath in zip(endfNonFreeThermalTapeList, relPathThermalAceList):
        deck.append(['ln -fs {0} tape{1:02g}'.format(relPath, endfTape)])
    deck.append(["echo 'Running NJOY'"])
    deck.append(["cat>input <<EOF"])
    create_moder_input(deck, dat, tapes.endf, tapes.bendf)
    deck.append(['stop'])
    deck.append(["EOF"])
    deck.append(["./xnjoy<input"])
    for iT in range(numTemperatures):
        # Do multiple calls to NJOY (bug in NJOY somewhere)
        deck.append(["cat>input <<EOF"])
        create_thermal_acer_input(deck, dat, tapes.bendf, pendfOutTape, tapes.aceStart,
            tapes.aceXSDirStart, iT)
        create_acer_visualization_input(deck, dat, tapes.aceStart, tapes.acerViewStart,
            tapes.viewrAceStart, iT)
        deck.append(['stop'])
        deck.append(["EOF"])
        deck.append(["./xnjoy<input"])
    deck.append(["echo 'Cleaning up and saving ACE files'"])
    aceFiles = []
    for i in range(numTemperatures):
        tapeAceOut, tapeXSDirOut, viewrOut = tapes.aceStart + i, tapes.aceXSDirStart + i, tapes.viewrAceStart + i
        aceFile = aceFileTemplate.format(H=dat.thermalFileName, s=i, e=dat.aceExt)
        aceFiles.append(aceFile)
        aceFileDir = 'xsdir.{}'.format(aceFile)
        plotName = 'p_xs_{}.pdf'.format(aceFile)
        deck.append(['cp -f tape{0:02g} {1}'.format(tapeAceOut, aceFile)])
        deck.append(['cp -f tape{0:02g} {1}'.format(tapeXSDirOut, aceFileDir)])
        deck.append(['ps2pdf tape{0:02g} {1}'.format(viewrOut, plotName)])
    deck.append(['rm -f xnjoy'])
    tapeNamesToRemove = ' '.join(get_temporary_tapes_ace(tapes, pendfOutTape, endfNonFreeThermalTapeList, numTemperatures, False))
    deck.append(['rm -f {0}'.format(tapeNamesToRemove)])
    try_mkdir(aceDirr)
    print_njoy_file(aceScriptPath, deck)

    return aceScriptPath, aceFiles

def create_njoy_script(dat, tapes):
    '''Create three NJOY files: one to generate a PENDF output, one to generate a GENDF output, and one to generate an ACE output'''

    #Directories
    dirDict = get_common_directories()
    pyDirr = dirDict['src']
    endfDirr = dirDict['endf']
    pendfRootDirr = dirDict['pendf']
    gendfRootDirr = dirDict['gendf']
    aceRootDirr = dirDict['ace']
    njoyPath = os.path.join(dirDict['njoyInstall'], 'xnjoy')
    #
    #TODO: Extend nuclideName to use thermalName
    pendfDirr = os.path.join(pendfRootDirr, str(dat.nuclideName))
    gendfDirr = os.path.join(gendfRootDirr, str(dat.nuclideName))
    aceDirr = os.path.join(aceRootDirr, str(dat.nuclideName))
    #
    relPathNJOYPendf = os.path.relpath(njoyPath, pendfDirr)
    relPathNJOYGendf = os.path.relpath(njoyPath, gendfDirr)
    relPathNJOYAce = os.path.relpath(njoyPath, aceDirr)
    #
    endfPath = os.path.join(endfDirr, dat.endfFile)
    relPathEndfPendf = os.path.relpath(endfPath, pendfDirr)
    relPathEndfGendf = os.path.relpath(endfPath, gendfDirr)
    relPathEndfAce = os.path.relpath(endfPath, aceDirr)
    #
    numThermals = len(dat.endfThermalFileList)
    numTemperatures = len(dat.thermList)
    pendfOutTape = get_pendf_out_tape(tapes, numThermals)
    endfThermalPathList = [os.path.join(endfDirr, thermalFile) for
        thermalFile in dat.endfThermalFileList if thermalFile]
    relPathThermalPendfList = [os.path.relpath(thermalPath, pendfDirr) for
        thermalPath in endfThermalPathList]
    endfThermalTapeList = map(lambda x: x + tapes.endfThermalStart, range(numThermals))
    endfNonFreeThermalTapeList = [tape for (mt,tape) in
        zip(dat.inelasticMTList, endfThermalTapeList) if mt != 221]
    #
    pendfScriptFile = 'runNJOY.sh'
    pendfScriptPath = os.path.join(pendfDirr, pendfScriptFile)
    pendfFile = 'pendf'
    pendfFileASCII = 'pendf_ascii.txt'
    pendfPath = os.path.join(pendfDirr, pendfFile)
    #
    gendfScriptFile = 'runNJOY.sh'
    gendfScriptPath = os.path.join(gendfDirr, gendfScriptFile)
    gendfFile = 'gendf'
    gendfPath = os.path.join(gendfDirr, gendfFile)
    relPathPendfGendf = os.path.relpath(pendfPath, gendfDirr)
    readgrouprPath = os.path.join(pyDirr, 'Readgroupr.py')
    #
    aceScriptFile = 'runNJOY.sh'
    aceScriptPath = os.path.join(aceDirr, aceScriptFile)
    # TODO: Generalize to use different extension root (instead of 9)
    aceFileTemplate = '{Z}{A:03}.{e}{s}c'
    relPathPendfAce = os.path.relpath(pendfPath, aceDirr)

    # Script to run NJOY from ENDF to PENDF:
    deck = []
    deck.append(['#! /usr/bin/env bash'])
    deck.append(["echo 'NJOY Problem {0} {1} (ENDF to PENDF)'".format(dat.nuclideName, dat.endfName)])
    deck.append(["echo 'Getting ENDF input tapes and NJOY executable'"])
    deck.append(['ln -fs {0} xnjoy'.format(relPathNJOYPendf)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathEndfPendf, tapes.endf)])
    for endfTape, relPath in zip(endfNonFreeThermalTapeList, relPathThermalPendfList):
        deck.append(['ln -fs {0} tape{1:02g}'.format(relPath, endfTape)])
    deck.append(["echo 'Running NJOY'"])
    deck.append(['cat>input <<EOF'])
    create_moder_input(deck, dat, tapes.endf, tapes.bendf)
    create_reconr_input(deck, dat, tapes.bendf, tapes.reconrOut)
    create_broadr_input(deck, dat, tapes.bendf, tapes.reconrOut, tapes.broadrOut)
    # Warning: 9547 (Am-242m) and (Np-238) have negative values here sometimes for low Sigma_0
    create_unresr_input(deck, dat, tapes.bendf, tapes.broadrOut, tapes.unresrOut)
    for i, endfThermal in enumerate(endfThermalTapeList):
        thermrInTape, thermrOutTape = get_thermr_tapes(tapes, i)
        if dat.inelasticMTList[i] == 221:
            endfThermal = 0
        create_thermr_input(deck, dat, i, endfThermal, thermrInTape, thermrOutTape)
    create_moder_input(deck, dat, pendfOutTape, tapes.pendfOutASCII)
    deck.append(['stop'])
    deck.append(['EOF'])
    deck.append(['./xnjoy<input'])
    deck.append(["echo 'Cleaning up and saving PENDF file'"])
    deck.append(['cp -f tape{0:02g} {1}'.format(abs(pendfOutTape), pendfFile)])
    deck.append(['cp -f tape{0:02g} {1}'.format(tapes.pendfOutASCII, pendfFileASCII)])
    deck.append(['rm -f xnjoy'])
    tapeNamesToRemove = ' '.join(get_temporary_tapes_pendf(tapes, endfNonFreeThermalTapeList, numThermals))
    deck.append(['rm -f {0}'.format(tapeNamesToRemove)])
    try_mkdir(pendfDirr)
    print_njoy_file(pendfScriptPath, deck)

    # Script to run NJOY from PENDF to GENDF:
    doPlotr = False
    deck = []
    deck.append(["#! /usr/bin/env bash"])
    deck.append(["echo 'NJOY Problem {0} {1} (PENDF to GENDF)'".format(dat.nuclideName, dat.endfName)])
    deck.append(["echo 'Getting ENDF input tape, PENDF input tape and NJOY executable'"])
    deck.append(['ln -fs {0} xnjoy'.format(relPathNJOYGendf)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathEndfGendf, tapes.endf)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathPendfGendf, abs(pendfOutTape))])
    deck.append(["echo 'Running NJOY'"])
    deck.append(["cat>input <<EOF"])
    create_moder_input(deck, dat, tapes.endf, tapes.bendf)
    create_groupr_input(deck, dat, tapes.bendf, pendfOutTape, tapes.grouprIn, tapes.grouprOut)
    create_moder_input(deck, dat, tapes.grouprOut, tapes.gendfOut)
    if False and dat.numGroups <= 200:
        doPlotr = True
        plotNames = create_plotr_inputs(deck, dat, tapes.gendfOut, tapes.plotrOut, tapes.viewrOut)
    deck.append(['stop'])
    deck.append(["EOF"])
    deck.append(["./xnjoy<input"])
    deck.append(["echo 'Cleaning up and saving GROUPR file'"])
    if doPlotr:
        for plotName, viewrOut in zip(plotNames, tapes.viewrOut):
            deck.append(['ps2pdf tape{0:02g} {1}'.format(viewrOut, plotName)])
    deck.append(['cp -f tape{0:02g} {1}'.format(tapes.gendfOut, gendfFile)])
    deck.append(['rm -f xnjoy'])
    tapeNamesToRemove = ' '.join(get_temporary_tapes_gendf(tapes, pendfOutTape, doPlotr))
    deck.append(['rm -f {0}'.format(tapeNamesToRemove)])
    optionalComment = ''
    if True and dat.numGroups > 200:
        optionalComment = '#'
    deck.append(["{0}echo 'Reading in XS file and pickling for future use'".format(optionalComment)])
    deck.append(['{0} {1} -I {2} -p all -f csc'.format(optionalComment, readgrouprPath, gendfFile)])
    try_mkdir(gendfDirr)
    print_njoy_file(gendfScriptPath, deck)

    # Script to run NJOY from PENDF to (non-thermal) ACE:
    deck = []
    deck.append(["#! /usr/bin/env bash"])
    deck.append(["echo 'NJOY Problem {0} {1} (PENDF to ACE)'".format(dat.nuclideName, dat.endfName)])
    deck.append(["echo 'Getting ENDF input tape, PENDF input tape and NJOY executable'"])
    deck.append(['ln -fs {0} xnjoy'.format(relPathNJOYAce)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathEndfAce, tapes.endf)])
    deck.append(['ln -fs {0} tape{1:02g}'.format(relPathPendfAce, abs(pendfOutTape))])
    deck.append(["echo 'Running NJOY'"])
    deck.append(["cat>input <<EOF"])
    create_moder_input(deck, dat, tapes.endf, tapes.bendf)
    if dat.usePURR:
        create_purr_input(deck, dat, tapes.bendf, pendfOutTape, tapes.purrOut)
    else:
        create_moder_input(deck, dat, pendfOutTape, tapes.purrOut)
    deck.append(['stop'])
    deck.append(["EOF"])
    deck.append(["./xnjoy<input"])
    for iT in range(numTemperatures):
        # Do multiple calls to NJOY (bug in NJOY somewhere)
        deck.append(["cat>input <<EOF"])
        create_acer_input(deck, dat, tapes.bendf, tapes.purrOut, tapes.aceStart,
            tapes.aceXSDirStart, iT)
        create_acer_visualization_input(deck, dat, tapes.aceStart, tapes.acerViewStart,
            tapes.viewrAceStart, iT)
        deck.append(['stop'])
        deck.append(["EOF"])
        deck.append(["./xnjoy<input"])
    deck.append(["echo 'Cleaning up and saving ACE files'"])
    aceFiles = []
    for i in range(numTemperatures):
        tapeAceOut, tapeXSDirOut, viewrOut = tapes.aceStart + i, tapes.aceXSDirStart + i, tapes.viewrAceStart + i
        aceFile = aceFileTemplate.format(A=dat.A, Z=dat.Z, s=i, e=dat.aceExt)
        aceFiles.append(aceFile)
        aceFileDir = 'xsdir.{}'.format(aceFile)
        plotName = 'p_xs_{}.pdf'.format(aceFile)
        deck.append(['cp -f tape{0:02g} {1}'.format(tapeAceOut, aceFile)])
        deck.append(['cp -f tape{0:02g} {1}'.format(tapeXSDirOut, aceFileDir)])
        deck.append(['ps2pdf tape{0:02g} {1}'.format(viewrOut, plotName)])
    deck.append(['rm -f xnjoy'])
    tapeNamesToRemove = ' '.join(get_temporary_tapes_ace(tapes, pendfOutTape, [], numTemperatures, True))
    deck.append(['rm -f {0}'.format(tapeNamesToRemove)])
    try_mkdir(aceDirr)
    print_njoy_file(aceScriptPath, deck)

    return pendfScriptPath, gendfScriptPath, aceScriptPath, aceFiles

#####################################################################################
def create_njoy_driver(pendfScriptPaths, gendfScriptPaths, aceScriptPaths):
    '''Create a driver script to run all the NJOY calculations in parallel'''
    dirDict = get_common_directories()
    runDirr = dirDict['xs']
    #
    pendfJobFile = 'runPendf.txt'
    pendfJobPath = os.path.join(runDirr, pendfJobFile)
    pendfDriverFile = 'RunPendf.sh'
    pendfDriverPath = os.path.join(runDirr, pendfDriverFile)
    pendfDirrs = [os.path.dirname(path) for path in pendfScriptPaths]
    relPathsPendfRun = [os.path.relpath(path, runDirr) for path in pendfDirrs]
    pendfScriptName = './{0}'.format(os.path.split(pendfScriptPaths[0])[-1])
    write_job_file(pendfJobPath, pendfScriptName, relPathsPendfRun)
    write_job_driver(pendfDriverPath, pendfJobFile)
    #
    gendfJobFile = 'runGendf.txt'
    gendfJobPath = os.path.join(runDirr, gendfJobFile)
    gendfDriverFile = 'RunGendf.sh'
    gendfDriverPath = os.path.join(runDirr, gendfDriverFile)
    gendfDirrs = [os.path.dirname(path) for path in gendfScriptPaths]
    relPathsPendfRun = [os.path.relpath(path, runDirr) for path in gendfDirrs]
    gendfScriptName = './{0}'.format(os.path.split(gendfScriptPaths[0])[-1])
    write_job_file(gendfJobPath, gendfScriptName, relPathsPendfRun)
    write_job_driver(gendfDriverPath, gendfJobFile)
    #
    aceJobFile = 'runAce.txt'
    aceJobPath = os.path.join(runDirr, aceJobFile)
    aceDriverFile = 'RunAce.sh'
    aceDriverPath = os.path.join(runDirr, aceDriverFile)
    aceDirrs = [os.path.dirname(path) for path in aceScriptPaths]
    relPathsPendfRun = [os.path.relpath(path, runDirr) for path in aceDirrs]
    aceScriptName = './{0}'.format(os.path.split(aceScriptPaths[0])[-1])
    write_job_file(aceJobPath, aceScriptName, relPathsPendfRun)
    write_job_driver(aceDriverPath, aceJobFile)
    #
    driverDriverFile = 'RunNJOY.sh'
    driverDriverPath = os.path.join(runDirr, driverDriverFile)
    driversToDrive = [pendfDriverFile, gendfDriverFile, aceDriverFile]
    driversToDrive = ['./{0}'.format(driver) for driver in driversToDrive]
    write_driver_driver(driverDriverPath, driversToDrive, 'Run ENDF -> PENDF -> GENDF / ACE in 2 steps.')

def create_ace_copier(aceFilesDict):
    '''Create a script that copies the ACE files to xdata'''
    dirDict = get_common_directories()
    scriptname = 'copyAce.sh'
    runDirr = dirDict['xdata']
    scriptPath = os.path.join(runDirr, scriptname)
    #
    deck = []
    deck.append(["#! /usr/bin/env bash"])
    deck.append([""])
    deck.append(["echo 'Copying ACE XS data'"])
    for nuclide in sorted(aceFilesDict):
        for aceFile in sorted(aceFilesDict[nuclide]):
            deck.append(["cp -f ../{}/{} .".format(nuclide, aceFile)])
    deck.append([""])
    deck.append(["echo 'Copying ACE directory data'"])
    deck.append(["rm -f xsdir_tail"])
    for nuclide in sorted(aceFilesDict):
        for aceFile in sorted(aceFilesDict[nuclide]):
            aceFileDir = 'xsdir.{}'.format(aceFile)
            deck.append(["awk '{{$3=\"{}\"; $4=\"0\"; print }}' ../{}/{} >> xsdir_tail".format(aceFile, nuclide, aceFileDir)])
    deck.append(["echo '{}' > xsdir_date".format(datetime.strftime(datetime.now(), '%m/%d/%Y'))])
    deck.append(["echo 'directory' >> xsdir_date"])
    deck.append(["cat xsdir_head xsdir_date xsdir_tail > xsdir"])
    deck.append([""])
    deck.append(["echo 'Copying ACE figures'"])
    for nuclide in sorted(aceFilesDict):
        for aceFile in sorted(aceFilesDict[nuclide]):
            plotName = 'p_xs_{}.pdf'.format(aceFile)
            deck.append(["cp -f ../{}/{} figures".format(nuclide, plotName)])
    print_njoy_file(scriptPath, deck)

#####################################################################################
def create_moder_input(deck, dat, tapeIn, tapeOut):
    # Convert from ascii to binary
    deck.append(['moder'])
    deck.append((tapeIn, tapeOut))

def create_reconr_input(deck, dat, tapeENDFIn, tapePENDFOut):
    # Reconstruct XS to within linear representation
    # See http://t2.lanl.gov/njoy/in-reconr.html
    ncards = 3
    deck.append(['reconr'])
    deck.append((tapeENDFIn, tapePENDFOut))
    deck.append(("'pendf tape for " + dat.nuclideName + " from " + dat.endfName + "'",'/'))
    deck.append((dat.mat, ncards, '/'))
    deck.append((dat.errTol, 0, dat.errTolMax, '/'))
    deck.append(("'" + dat.nuclideName + " from " + dat.endfName + " tape'",'/'))
    deck.append(("'processed by the njoy nuclear data processing system'",'/'))
    deck.append(("'see original " + dat.endfName + " tape for details of evaluation'",'/'))
    deck.append((0, '/'))

def create_broadr_input(deck, dat, tapeENDFIn, tapePENDFIn, tapePENDFOut):
    # Broaden XS
    # See http://t2.lanl.gov/njoy/in-broadr.html
    deck.append(['broadr'])
    deck.append((tapeENDFIn, tapePENDFIn, tapePENDFOut))
    deck.append((dat.mat, len(dat.thermList), 0, 0, 0, '/'))
    deck.append((dat.errTol, '/'))
    deck.append(dat.thermList)
    deck.append((0, '/'))

def create_unresr_input(deck, dat, tapeENDFIn, tapePENDFIn, tapePENDFOut):
    # Fix unresolved resonance region
    # http://t2.lanl.gov/njoy/in-unresr.html
    deck.append(['unresr'])
    deck.append((tapeENDFIn, tapePENDFIn, tapePENDFOut))
    deck.append((dat.mat, len(dat.thermList), len(dat.sig0List), 0))
    deck.append(dat.thermList)
    deck.append(dat.sig0List)
    deck.append((0, '/'))

def create_thermr_input(deck, dat, index, tapeENDFThermalIn, tapePENDFIn, tapePENDFOut):
    # Apply thermal cross sections
    # See http://t2.lanl.gov/njoy/in-thermr.html
    inelasticMT = dat.inelasticMTList[index]
    inelasticOpt = get_inelastic_option(inelasticMT)
    elasticOpt = get_elastic_option(inelasticMT)
    matThermal = dat.matThermalList[index]
    numAtom = get_num_atoms_in_molecule(inelasticMT)
    deck.append(['thermr'])
    deck.append((tapeENDFThermalIn, tapePENDFIn, tapePENDFOut))
    # For NJOY 2012, before numAtom, insert iform (probably should be 0).
    # NB NJOY2012 messes up the output of the thermal data
    # in the GENDF file, which messes up Readgroupr.py
    deck.append((matThermal, dat.mat, dat.numAngles, len(dat.thermList), inelasticOpt, elasticOpt, numAtom, inelasticMT, 1))
    deck.append(dat.thermList)
    deck.append((dat.thermrTol, dat.Emax))

def create_purr_input(deck, dat, tapeENDFThermalIn, tapePENDFIn, tapePENDFOut):
    # Create probability tables for the unresolved resonance range (for ACE files)
    # See http://t2.lanl.gov/njoy/in-purr.html
    deck.append(['purr'])
    deck.append((tapeENDFThermalIn, tapePENDFIn, tapePENDFOut))
    deck.append((dat.mat, len(dat.thermList), len(dat.sig0List), dat.numProbBins, dat.numLadders, 0, '/'))
    deck.append(dat.thermList)
    deck.append(dat.sig0List)
    deck.append((0, '/'))

def create_acer_input(deck, dat, tapeENDFIn, tapePENDFIn, tapeACEROutStart, tapeXSDIROutStart, tapeIndex):
    '''Create non-thermal ACE files'''
    tapeACEROut = tapeACEROutStart + tapeIndex
    tapeXSDIROut = tapeXSDIROutStart + tapeIndex
    thermStr = '{0:g}'.format(round(dat.thermList[tapeIndex],1))
    # May not match suffix in materials_materials.py
    suffixNonThermal = '.{}{}'.format(dat.aceExt, tapeIndex)
    #
    deck.append(['acer'])
    deck.append((tapeENDFIn, tapePENDFIn, 0, tapeACEROut, tapeXSDIROut))
    # first 1,0,1 = fast(non-thermal), no-print, ascii
    deck.append((1, 0, 1, suffixNonThermal, '/'))
    deck.append(("'{} ACE at {} K from {}'".format(dat.nuclideName, thermStr, dat.endfName), '/'))
    deck.append((dat.mat, thermStr, '/'))
    deck.append((dat.useNewAngularDist, 0, '/'))
    deck.append(( '/'))

def create_thermal_acer_input(deck, dat, tapeENDFIn, tapePENDFIn, tapeACEROutStart, tapeXSDIROutStart, tapeIndex):
    '''Create thermal ACE files'''
    tapeACEROut = tapeACEROutStart + tapeIndex
    tapeXSDIROut = tapeXSDIROutStart + tapeIndex
    thermStr = '{0:g}'.format(round(dat.thermList[tapeIndex],1))
    # May not match suffix in materials_materials.py
    suffixThermal = '.{}{}'.format(dat.aceExt, tapeIndex)
    #
    # A thermal ACE file is done for one thermal material only
    inelasticMT = dat.inelasticMTList[0]
    elasticMT = get_elasticMT(inelasticMT)
    elasticCoh= get_elastic_coherence(inelasticMT)
    #
    inelasticOpt = get_inelastic_option(inelasticMT)
    elasticOpt = get_elastic_option(inelasticMT)
    # Thermal ACE is done for one bound thermal treatment only, so use the first item in the list
    matThermal = dat.matThermalList[0]
    numAtom = get_num_atoms_in_molecule(inelasticMT)
    #
    deck.append(['acer'])
    deck.append((tapeENDFIn, tapePENDFIn, 0, tapeACEROut, tapeXSDIROut))
    # first 2,0,1 = thermal, no-print, ascii
    deck.append((2, 0, 1, suffixThermal, len(dat.ZAIDs), '/'))
    deck.append(("'{} thermal ACE at {} K from {}'".format(dat.thermalNuclideName, thermStr, dat.endfName), '/'))
    zaidStr = ' '.join(['{} 0'.format(ZAID) for ZAID in dat.ZAIDs])
    deck.append([zaidStr])
    deck.append((dat.mat, thermStr, "'{}'".format(dat.thermalNuclideName), '/'))
    # ZAID of nuclides for which this thermal XS applies
    zaidStr = ' '.join(['{}'.format(ZAID) for ZAID in dat.ZAIDs[:3]]) + '/'
    deck.append([zaidStr])
    deck.append((inelasticMT, dat.numThermalAceBins, elasticMT, elasticCoh, dat.numMix, dat.emaxThermalAce, 2))

def create_acer_visualization_input(deck, dat, tapeACERInStart, tapeACEROutStart, tapeVIEWROutStart, tapeIndex):
    '''Create ps of the cross sections from the ACE files'''
    tapeACERIn = tapeACERInStart + tapeIndex
    tapeACEROut = tapeACEROutStart + tapeIndex
    tapeVIEWROut = tapeVIEWROutStart + tapeIndex
    acerOpt = 7
    #
    deck.append(['acer'])
    deck.append((0, tapeACERIn, tapeACEROut, 0, 0))
    deck.append((acerOpt, 0, 1, '/'))
    deck.append(['/'])
    deck.append(['viewr'])
    deck.append((tapeACEROut, tapeVIEWROut, '/'))

def create_plotr_inputs(deck, dat, tapePLOTRIn, tapesPLOTROut, tapesVIEWROut):
    '''Create ps of the flux and total cross section at extreme temperatures and 3 sigma0 values'''
    numTherm = len(dat.thermList)
    numSig0 = len(dat.sig0List)
    thermIndexList = [0, numTherm - 1]
    sig0IndexList = [0, min(numSig0 - 1, (1*numSig0) / 2), numSig0 - 1]
    thermList = [dat.thermList[i] for i in thermIndexList]
    sig0List = [dat.sig0List[i] for i in sig0IndexList]
    thermStrList = ['{0:g}'.format(round(val,1)) for val in thermList]
    sig0StrList = ['{0:.1e}'.format(val) for val in sig0List]
    lineStyles = [0, 1, 2]
    colors = [2, 3, 4]
    #
    whatToPlotList = [0, 1]
    axisStrList = ["'<f>flux (1/lethargy)'", '']
    titleStrList = ['<f>lux at', '<t>otal <xs> at']
    #
    tapeIndex = 0
    for thermIndex, thermVal, thermStr in zip(thermIndexList, thermList, thermStrList):
        for toPlotNum, axisStr, titleStr in zip(whatToPlotList, axisStrList, titleStrList):
            deck.append(['plotr'])
            deck.append((tapesPLOTROut[tapeIndex], '/'))
            deck.append(['/'])
            plotCount = 1
            for sig0Index, sig0Val, sig0Str in zip(sig0IndexList, sig0List, sig0StrList):
                deck.append((plotCount, '/'))
                if plotCount == 1:
                    deck.append(("'<endf/b-vii1 {0}>{1}'".format(dat.nuclideName[0], dat.nuclideName[1:]), '/'))
                    deck.append(("'{0} {1} K'".format(titleStr, thermStr), '/'))
                    deck.append((4, 0, 2, 1, '/'))
                    deck.append(['/']) #x-axis low/high
                    deck.append(['/'])
                    deck.append(['/']) #y-axis low/high
                    deck.append((axisStr, '/'))
                deck.append((1, tapePLOTRIn, dat.mat, 3, 1, thermVal, toPlotNum, sig0Index+1, 1, '/'))
                deck.append((0, 0, lineStyles[plotCount-1], colors[plotCount-1], 3, 0, '/'))
                deck.append(("'{0}'".format(sig0Str), '/'))
                plotCount += 1
            deck.append((99, '/'))
            deck.append(['viewr'])
            deck.append((tapesPLOTROut[tapeIndex], tapesVIEWROut[tapeIndex]))
            tapeIndex += 1
    names = []
    for thermStr in thermStrList:
        names.append('p_flux_{0}.pdf'.format(thermStr))
        names.append('p_xs_{0}.pdf'.format(thermStr))
    return names

def create_groupr_input(deck, dat, tapeENDFIn, tapePENDFIn, tapeGroupsIn, tapeGENDFOut):
    '''Create multigroup cross sections
    See http://t2.lanl.gov/njoy/in-groupr.html'''
    # Add in more general weighting options (iwt=1,-4)?

    # Pre-process groupr input
    # Add currently recognized names to rxts
    mtMap = dict(mt1="'total'", mt2="'elastic'", mt18="'fission'", mt102="'capture'",
                  mt221="'free gas thermal'")
    namedRxts = []
    #
    namedRxts.append((3, '/'))
    namedRxts.append((3, 259, "'invel'", '/'))
    for mt in dat.thermalMTList:
        namedRxts.append((3,mt,"'(thermal)'",'/'))
    if dat.isFissionable:
        namedRxts.append((3, 452, "'nu'", '/'))
    if dat.includeMF6:
        namedRxts.append((6, '/'))
        for mt in dat.thermalMTList:
            namedRxts.append((6,mt,"'(thermal)'",'/'))
    elif dat.isFissionable:
        # To get chi, currently need (6,18), so use this even if not including other MF6
        namedRxts.append((6, 18,"'fission'",'/'))

    namedRxts = sorted(namedRxts)
    # repeat total, elastic, fission, capture, thermal for each temperature
    toRepeat = [(3,1,"'total'",'/'), (3,2,"'elastic scat'",'/'), (3,102,"'rad capture'",'/')]
    if dat.isFissionable:
        toRepeat.append((3,18,"'nu'",'/'))
    for mt in dat.thermalMTList:
        toRepeat.append((3,mt,"'(thermal)'",'/'))
    if dat.includeMF6:
        toRepeat.append((6,2,"'elastic scat'",'/'))
        for mt in dat.thermalMTList:
            toRepeat.append((6,mt,"'(thermal)'",'/'))
    elif dat.isFissionable:
        # To get chi, currently need (6,18), so use this even if not including other MF6
        toRepeat.append((6,18,"'fission'",'/'))
    toRepeat = sorted(toRepeat)
    groupMap = [-1, -1, 240, 30, 27, 50, 68, 100, 35, 69, 187, 70, 620, 80, 100, 640, 174, 175]
    dat.numGroups = groupMap[dat.groupOpt]
    if dat.groupOpt == 1:
        dat.numGroups = len(dat.groupBdrs) - 1
    #
    deck.append(['groupr'])
    deck.append((tapeENDFIn, tapePENDFIn, tapeGroupsIn, tapeGENDFOut))
    deck.append((dat.mat, dat.groupOpt, 0, dat.weightOpt, dat.legendreOrder, len(dat.thermList), len(dat.sig0List), 0, '/'))
    deck.append((dat.nuclideName,'/'))
    deck.append(dat.thermList)
    deck.append(dat.sig0List)
    if dat.groupOpt == 1:
        deck.append((dat.numGroups,'/'))
        groupFormat = []
        rowSize = 4
        rowPos = 1
        for group in dat.groupBdrs:
            # Single precision uses 23 bits in the mantissa (and one implicit bit).
            # Since it is base-2, this ~ log10(2**(23+1)) = 7.22 total decimals may be
            # specified. This translates into ~ 6 digits after the '.' in scientific notation.
            groupFormat.append("{0:.6e}".format(group))
            if rowPos > rowSize:
                deck.append(groupFormat)
                groupFormat = []
                rowPos = 1
            else:
                rowPos += 1
        groupFormat.append('/')
        if rowPos != 1:
            deck.append(groupFormat)
    for rxt in namedRxts:
        deck.append(rxt)
    deck.append((0, '/'))
    for T in dat.thermList[1:]:
        for rxt in toRepeat:
            deck.append(rxt)
        deck.append((0, '/'))
    deck.append((0, '/'))

###############################################################################
def print_njoy_file(filePath, deck):
    """Given a file path and deck, print out an NJOY file."""

    with open(filePath, 'w') as fid:
        for card in deck:
            sz = len(card)
            if sz == 1:
                fid.write(str(card[0])+'\n')
            else:
                for i,obj in enumerate(card):
                    val = str(obj)
                    if val == '/':
                        fid.write('/\n')
                    elif i == sz-1:
                        fid.write(' ' + val + '\n')
                    elif i == 0:
                        fid.write(val)
                    else:
                        fid.write(' ' + val)
        try_make_exe(filePath)

def write_job_file(filePath, fName, dirrs):
    with open(filePath, 'w') as fid:
        for dirr in dirrs:
            fid.write("echo 'Running {0}'\n".format(dirr))
            fid.write('cd {0}\n'.format(dirr))
            fid.write('{0}\n'.format(fName))
            fid.write('cd -\n')
            fid.write('\n')

def write_job_driver(driverPath, jobFile):
    '''Write a .sh file to execute a jobFile in parallel'''
    with open(driverPath, 'w') as fid:
        fid.write('#! /usr/bin/env bash\n')
        fid.write('\n')
        fid.write('# Run all of {0} in parallel, using a blank line as a delimiter\n'.format(jobFile))
        fid.write("cat {0} | parallel -d '\\n\\n' --eta -j-1\n".format(jobFile))
    try_make_exe(driverPath)

def write_driver_driver(driverDriverPath, driversToDrive, description='Run drivers'):
    '''Write a .sh file that calls each of driversToDrive, in order'''
    with open(driverDriverPath, 'w') as fid:
        fid.write('#! /usr/bin/env bash\n')
        fid.write('\n')
        fid.write('# {0}\n'.format(description))
        for driver in driversToDrive:
            fid.write('{0}\n'.format(driver))
    try_make_exe(driverDriverPath)

###############################################################################
class NJOYTape():
    def __init__(self):
        # The starts assume at most 10 temperatures are used
        self.endfThermalStart = 40
        self.aceStart = 50
        self.aceXSDirStart = 60
        self.acerViewStart = 70
        self.viewrAceStart = 80
        self.endf = 21
        self.bendf = -22
        self.reconrOut = -23
        self.broadrOut = -24
        self.unresrOut = -25
        self.thermrOutA = -26
        self.thermrOutB = -27
        self.grouprIn = 0
        self.grouprOut = -28
        self.gendfOut = 29
        self.pendfOutASCII = 30
        self.plotrOut = [31, 33, 35, 37]
        self.viewrOut = [32, 34, 36, 38]
        self.purrOut = -39


def get_temporary_tapes_pendf(tapes, endfThermalTapeList, numThermals):
    tapeNumsToRemove = [tapes.endf, tapes.bendf, tapes.reconrOut, tapes.broadrOut, tapes.unresrOut, tapes.pendfOutASCII]
    if numThermals:
        tapeNumsToRemove.append(tapes.thermrOutA)
    if numThermals > 1:
        tapeNumsToRemove.append(tapes.thermrOutB)
    for endfThermalTape in endfThermalTapeList:
        tapeNumsToRemove.append(endfThermalTape)
    tapeNamesToRemove = ['tape{0:02g}'.format(abs(tape)) for tape in tapeNumsToRemove]
    return tapeNamesToRemove

def get_temporary_tapes_gendf(tapes, pendfOutTape, doPlotr):
    tapeNumsToRemove = [tapes.endf, pendfOutTape, tapes.bendf, tapes.grouprOut, tapes.gendfOut]
    if doPlotr:
        tapeNumsToRemove += tapes.plotrOut
        tapeNumsToRemove += tapes.viewrOut
    tapeNamesToRemove = ['tape{0:02g}'.format(abs(tape)) for tape in tapeNumsToRemove]
    return tapeNamesToRemove

def get_pendf_out_tape(tapes, numThermals):
    if numThermals == 0:
        return tapes.unresrOut
    elif numThermals % 2 == 1:
        return tapes.thermrOutA
    else:
        return tapes.thermrOutB

def get_thermr_tapes(tapes, i):
    '''Shuffle between thermrOutA and thermrOutB for arbitrary number of thermr tapes'''
    if i == 0:
        thermrInTape = tapes.unresrOut
        thermrOutTape = tapes.thermrOutA
    elif i % 2 == 1:
        thermrInTape = tapes.thermrOutA
        thermrOutTape = tapes.thermrOutB
    else:
        thermrInTape = tapes.thermrOutB
        thermrOutTape = tapes.thermrOutA
    return  thermrInTape, thermrOutTape

def get_temporary_tapes_ace(tapes, pendfOutTape, endfThermalTapeList, numTemperatures, usePurr):
    tapeNumsToRemove = [tapes.endf, tapes.bendf, pendfOutTape]
    if usePurr:
        tapeNumsToRemove.append(tapes.purrOut)
    for endfThermalTape in endfThermalTapeList:
        tapeNumsToRemove.append(endfThermalTape)
    for i in range(numTemperatures):
        tapeNumsToRemove.append(tapes.aceStart + i)
    for i in range(numTemperatures):
        tapeNumsToRemove.append(tapes.aceXSDirStart + i)
    for i in range(numTemperatures):
        tapeNumsToRemove.append(tapes.acerViewStart + i)
    for i in range(numTemperatures):
        tapeNumsToRemove.append(tapes.viewrAceStart + i)
    tapeNamesToRemove = ['tape{0:02g}'.format(abs(tape)) for tape in tapeNumsToRemove]
    return tapeNamesToRemove

###############################################################################
class NJOYDat():
    def __init__(self, A=None, Z=None, nuclideName=None, thermalNuclideName=None, thermalFileName=None, endfName=None, endfFile=None, mat=None, aceExt=None, thermList=[0], sig0List=[1e10], groupBdrs=None, groupOpt=0, isFissionable=False, includeMF6=True, usePURR=True, legendreOrder=0, inelasticMTList=[], thermalMTList=[], matThermalList=[], endfThermalFileList=[]):
        # General
        # Nuclide information
        self.A = A
        self.Z = Z
        self.nuclideName = nuclideName
        self.thermalNuclideName = thermalNuclideName
        self.thermalFileName = thermalFileName
        # which ENDF library to use for the nuclide
        self.endfName = endfName
        self.endfFile = endfFile
        # which ENDF library to use for the thermal treatment
        self.endfThermalFileList = endfThermalFileList
        # List of material numbers used for the thermal reaction (not the MT number)
        self.matThermalList = matThermalList
        # Material number used for the nuclide
        self.mat = mat
        # Relative error tolerances in XS reconstruction
        self.errTol = 0.001
        self.errTolMax = 0.005
        # An array of the temperatures to use (in K)
        self.thermList = thermList
        # An array of the background XS to use (in b; use 1.e10 for infinite [no energy shielding])
        self.sig0List = sig0List
        #
        # THERMR
        # number of equi-probable angles
        self.numAngles = 16
        # List of MT numbers for (inelastic) reactions (use 221 for free, see get_endf_mt_list)
        self.inelasticMTList = inelasticMTList
        self.thermalMTList = thermalMTList
        # tolerance
        self.thermrTol = 0.001
        # maximum energy for thermal treatment (eV) (originally used 1.9; others use 4.0)
        # self.Emax = 2.9
        #self.Emax = 3.3 # works with wiggle room for SHEM-361 and edits-12
        #self.Emax = 3.9 # does not work for ""
        self.Emax = 3.6 # works for ""
        #
        # GROUPR
        self.includeMF6 = includeMF6
        self.isFissionable = isFissionable
        self.groupBdrs = groupBdrs
        # ign: 1 for custom structure (see manual)
        self.groupOpt = groupOpt
        # iwt: weighting spectrum (see manual)
        self.weightOpt = 5
        # How many Legendre orders of transfer matrices to output; n in Pn
        self.legendreOrder = legendreOrder
        #
        # PURR
        self.usePURR = usePURR 
        self.numProbBins = 20
        self.numLadders = 50
        #
        # ACE
        self.aceExt = aceExt
        self.useNewAngularDist = 1
        # Thermal ACE
        self.numMix = 1 #True except for BeO or C6H6 (Benzene)
        self.numThermalAceBins = 128
        self.emaxThermalAce = 1000.
        # ZAIDs for which the thermal ACE treatment applies
        self.ZAIDs = []

def get_elasticMT(inelasticMT):
    # From Table 25 in NJOY2012 manual (ref: njoy2012 in materials_utilities.py)
    if inelasticMT in [222, 228, 227]:
        return 0
    else:
        return inelasticMT + 1

def get_elastic_coherence(inelasticMT):
    # 0: coherent elastic (or no data); 1: incoherent elastic
    # From Table 4 in NJOY2012 manual
    if inelasticMT in [223, 225, 235]:
        return 1
    else:
        return 0

def get_num_atoms_in_molecule(inelasticMT):
    # It is unclear what n is for the NSC TRIGA reactor for the ZrH_n fuel (MT 225)
    numAtomsInMolecule = {221: 1, 222: 2, 223: 2, 225: 1, 229: 1, 235: 1, 239: 2, 241: 1, 243: 1, 245: 1}
    return numAtomsInMolecule[inelasticMT]

def get_elastic_option(inelasticMT):
    if inelasticMT == 221:
        # Free gas
        return 0
    else:
        # Read from thermal ENDF file
        return 1

def get_inelastic_option(inelasticMT):
    if inelasticMT == 221:
        # Free gas
        return 1
    else:
        # Read S(alpha, beta)
        return 4
    #else:
    #    return 0
