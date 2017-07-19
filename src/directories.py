'''
Andrew Till
Winter 2016

Functions to create and determine the locations of the directories used in barnfire.
For the scripts to work, the environment variable SCRATCH_BARN must be set and exported
'''

#STDLIB
import os
import sys
import shutil
import stat

#####################################################################################
def get_common_directories():
    dirDict = get_root_directories()
    populate_leaf_directories(dirDict)
    return dirDict

def get_root_directories():
    dirDict = {}
    headDirr = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # If this fails, you probably did not export the variable SCRATCH_BARN:
    scratchDirr = os.path.abspath(os.environ['SCRATCH_BARN'])
    # The cross-section data can be quite large. You may wish to locate this on a scratch drive:
    # If this fails, you probably did not export the variable ENDF:
    endfDirr = os.path.abspath(os.environ['ENDF'])
    #
    # This directory should contain 'njoy', the NJOY executable:
    # If this fails, you probably did not export the variable NJOY:
    njoyDirr = os.path.abspath(os.environ['NJOY'])
    if os.path.isfile(njoyDirr):
        njoyDirr = os.path.dirname(njoyDirr)
    #
    dirDict['head'] = headDirr
    dirDict['scratch'] = scratchDirr
    dirDict['endf'] = endfDirr
    dirDict['njoyBuild'] = njoyDirr
    return dirDict

def get_source_directories(headDirr):
    '''The head directory contains directories that contain the Python scripts, NJOY, and common data.'''
    dirDict = {}
    dirDict['src'] = os.path.join(headDirr, 'src')
    # There are two dat directories: the main one at ../dat, and the one in scratch
    mainDatDirr = os.path.join(headDirr, 'dat')
    dirDict['common_group_structures'] = os.path.join(mainDatDirr , 'energy_groups')
    # The initialize script should create and populate these
    dirDict['nuclideData'] = os.path.join(mainDatDirr, 'nuclide-data')
    dirDict['thermalXSSource'] = os.path.join(mainDatDirr, 'thermal_endf')
    dirDict['nuclideWeights'] = os.path.join(mainDatDirr, 'nuclide_weights')
    return dirDict

def get_scratch_directories(scratchDirr):
    '''Each problem writes its local output to scratch. These are the directories contained therein.'''
    scratchDict = {}
    # xs
    scratchDict['xs'] = os.path.join(scratchDirr, 'xs')
    scratchDict['pendf'] = os.path.join(scratchDirr, 'xs/pendf')
    scratchDict['gendf'] = os.path.join(scratchDirr, 'xs/gendf')
    scratchDict['ace'] = os.path.join(scratchDirr, 'xs/ace')
    scratchDict['xdata'] = os.path.join(scratchDirr, 'xs/ace/xdata')
    scratchDict['xdata/figures'] = os.path.join(scratchDirr, 'xs/ace/xdata/figures')
    scratchDict['pdtxs'] = os.path.join(scratchDirr, 'xs/pdtxs')
    scratchDict['njoyInstall'] = os.path.join(scratchDirr, 'xs/bin')
    # dat (output and temporary data files)
    scratchDatDirr = os.path.join(scratchDirr, 'dat')
    # energy_groups: clust files; indicators: flux, fluxe, wgt, aptn files
    scratchDict['dat/energy_groups'] = os.path.join(scratchDatDirr , 'energy_groups')
    scratchDict['dat/indicators'] = os.path.join(scratchDatDirr , 'indicators')
    # figures (output)
    figuresDirr = os.path.join(scratchDirr, 'figures')
    scratchDict['figures/indicators'] = os.path.join(figuresDirr , 'indicators')
    scratchDict['figures/clustering'] = os.path.join(figuresDirr , 'clustering')
    return scratchDict

def populate_leaf_directories(dirDict):
    dirDict.update(get_scratch_directories(dirDict['scratch']))
    dirDict.update(get_source_directories(dirDict['head']))

#####################################################################################
def make_scratch_directories():
    '''Create the scratch directories, if they are not already created.'''
    scratchDict = get_scratch_directories(get_root_directories()['scratch'])
    for dirr in scratchDict.values():
        try_mkdir(dirr)

def copy_xnjoy():
    '''Copy xnjoy to the scratch location'''
    dirDict = get_common_directories()
    njoyExecutableIn = 'njoy'
    njoyExecutableOut = 'xnjoy'
    inPath = os.path.join(dirDict['njoyBuild'], njoyExecutableIn)
    outPath = os.path.join(dirDict['njoyInstall'], njoyExecutableOut)
    try_cp(inPath, outPath)

def make_endf_directory():
    dirDict = get_common_directories()
    try_mkdir(dirDict['endf'])

def copy_thermal_endf_xs():
    dirDict = get_common_directories()
    sourceDirr = dirDict['thermalXSSource']
    destinationDirr = dirDict['endf']
    copy_contents(sourceDirr, destinationDirr)

def copy_common_group_structures():
    dirDict = get_common_directories()
    sourceDirr = dirDict['common_group_structures']
    destinationDirr = dirDict['dat/energy_groups']
    copy_contents(sourceDirr, destinationDirr)

def copy_xs_script(callingFile):
    '''If callingFile is a file, copy it to scratch to maintain a record of
    how the XS were created'''
    if callingFile is None:
        return
    callingFile = os.path.abspath(callingFile)
    dirDict = get_common_directories()
    outPath = os.path.join(dirDict['scratch'], 'script_sav_00.sh')
    i = 0
    while os.path.isfile(outPath):
        outPath = os.path.join(dirDict['scratch'], 'script_sav_{:02}.sh'.format(i))
        i += 1
    if os.path.isfile(callingFile):
        try_cp(callingFile, outPath)

def copy_nuclide_weights():
    dirDict = get_common_directories()
    inPath = os.path.join(dirDict['nuclideWeights'], 'atomic_wgt_ratios.txt')
    outPath = os.path.join(dirDict['xdata'], 'xsdir_head')
    try_cp(inPath, outPath)

#####################################################################################
def try_mkdir(dirr):
    '''Attempts to make directory dirr, if it does not already exist'''
    try:
        if not os.path.isdir(dirr):
            # Makes all intermediate directories as well
            os.makedirs(dirr)
    except OSError:
        pass

def try_cp(inPath, outPath):
    '''Attempts to copy file from inPath to outPath, if outPath does not already exist'''
    try:
        if os.path.isfile(inPath) and not os.path.isfile(outPath):
            shutil.copy2(inPath, outPath)
    except OSError:
        pass

def try_make_exe(filename):
    try:
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IEXEC)
    except OSError:
        pass

def copy_contents(sourceDirr, destinationDirr):
    sourceFiles = os.listdir(sourceDirr)
    for sourceFile in sourceFiles:
        inPath = os.path.join(sourceDirr, sourceFile)
        try_cp(inPath, destinationDirr)
