#! /usr/bin/env python

'''
Andrew Till
Summer 2014

Create energy grid based on an existing group structure, e.g., the 44-group structure.

Highly resolve specific portions of the resonance structure.
'''

#STDLIB
import os
import argparse
#TPL
import numpy as np
#MINE
from directories import get_common_directories
from readxs import read_group_file, read_xs_file

#####################################################################################
def resolve_group_structure(inputDat, inputGrid=None):
    if not hasattr(inputDat, 'finishedParsing'):
        raise ValueError('Inputs have not finished being parsed. Use finish_parsing_inputs().')
    verbosity = inputDat.verbosity
    inDirr = inputDat.inDirr
    inName = inputDat.inName
    outDirr = inputDat.outDirr
    outName = inputDat.outName
    noExtend = inputDat.noExtend
    resolvedRange = inputDat.resolvedRange
    fullRange = inputDat.fullRange
    desiredGroups = inputDat.groups
    fineLethargySpacing = inputDat.fineLethargySpacing
    boundaryLethargySpacing = inputDat.boundaryLethargySpacing
    #
    # Read inital group structure
    if inputGrid is None:
        filePath = os.path.join(inDirr, inName)
        eGrid = np.array(read_group_file(filePath))
    else:
        eGrid = np.array(inputGrid)
    eGridStart = eGrid.copy()
    # Compress or expand to full range
    if not noExtend:
        eGrid = crop_energy_range(fullRange, eGrid)
    # Determine spacing to use from groups and spacing
    existingGroups = len(eGrid) - 1
    desiredResGroups = desiredGroups - existingGroups - 2
    # STILL CAN USE SOME WORK
    if desiredResGroups > 0:
        # Insert new group boundaries for boundary regions
        startE = resolvedRange[0]
        stopE = resolvedRange[1]
        bdrRatio = np.exp(boundaryLethargySpacing)
        leftBdrStart = startE / bdrRatio
        leftBdrStop = startE
        rightBdrStart = stopE
        rightBdrStop = stopE * bdrRatio
        newBoundaries = [leftBdrStart, leftBdrStop, rightBdrStart, rightBdrStop]
        eGrid = add_group_boundaries(newBoundaries, eGrid)
        # Overwrite group boundaries for resolved regions
        eGrid = overwrite_group_boundaries(startE, stopE, desiredGroups, eGrid, verbosity)
    else:
        print 'Warning! Number of desired groups is too small compared to the existing group structure.'
    # Save new energy grid
    if outName is not None:
        filename = 'energy_out_{0}.txt'.format(outName)
        filePath = os.path.join(outDirr, filename)
        write_egrid(filePath, eGrid)
    return eGrid

#####################################################################################
def crop_energy_range(newRange, oldEGrid):
    '''Crop eGrid to newRange. oldEGrid is reverse sorted.'''
    lowE, highE = min(newRange), max(newRange)
    lowG = np.argmin(np.abs(oldEGrid - highE))
    highG = np.argmin(np.abs(oldEGrid - lowE))
    newEGrid = []
    if oldEGrid[lowG] != highE:
        newEGrid.append(highE)
    if oldEGrid[lowG] > highE:
        lowG += 1
    if oldEGrid[highG] != lowE:
        newEGrid.append(lowE)
    if oldEGrid[highG] < lowE:
        highG -= 1
    newEGrid += list(oldEGrid[lowG:highG+1])
    return np.array(sorted(newEGrid, reverse=True))

def add_group_boundaries(boundaryList, eGrid):
    '''Add boundaryList into eGrid. Do not worry about duplicates at this point.'''
    eGrid = np.unique(np.concatenate([eGrid, boundaryList]))[::-1]
    return eGrid

def overwrite_group_boundaries(lowE, highE, desiredGroups, eGrid, verbosity):
    '''Add the boundaries of the resolved regions, so min(abs(eGrid-low/highE))) == 0'''
    gStart = np.argmin(np.abs(eGrid - highE))
    gStop = np.argmin(np.abs(eGrid - lowE))
    gExisting = (len(eGrid) - gStop - 1) + (gStart + 1)
    desiredResGroups = desiredGroups - gExisting + 2
    uTop = np.log(highE / lowE)
    uInsert = np.linspace(0.0, uTop, desiredResGroups)
    eInsert = highE * np.exp(-uInsert)
    eInsert[-1] = lowE
    if verbosity > 2:
        i = 0
        print 'Unrefined energy grid (high-energy portion)'
        for e in sorted(set(eGrid[0:gStart+1]), reverse=True):
            print '{0:.4e},'.format(e),
            i += 1
        print ''
        j = 0
        print 'Unrefined energy grid (low-energy portion)'
        for e in sorted(set(eGrid[gStop:]), reverse=True):
            print '{0:.4e},'.format(e),
            j += 1
        print ''
        print 'Number of high-energy groups', i
        print 'Number of low-energy groups', j
    if verbosity > 3:
        print 'Added energy grid (resonance region)'
        print eInsert[1:-1]
    eGrid = np.concatenate([eGrid[0:gStart], eInsert, eGrid[gStop+1:]])
    return eGrid

#####################################################################################
def write_egrid(filePath, eGrid):
    numGroups = len(eGrid) - 1
    with open(filePath, 'w') as fid:
        fid.write('User-defined {0}-group library\n'.format(numGroups))
        fid.write('group upper bound region\n')
        for g, Eg in enumerate(eGrid):
            if Eg > 2.5E4:
                Estr = 'fast'
            elif Eg <= 3.0:
                Estr = 'thermal'
            else:
                Estr = 'resonance'
            fid.write('{0} {1:.10e} {2}\n'.format(g+1, Eg, Estr))
        fid.write('\n')

#####################################################################################
def set_input_dat_low(inputDat):
    inputDat.resolvedRange = [3.0, 55.6]
    inputDat.fineLethargySpacing = 1/500.0
    inputDat.outName = 'low'

def set_input_dat_med(inputDict):
    inputDat.resolvedRange = [20.5, 36.5]
    inputDat.fineLethargySpacing = 1/3000.0
    inputDat.outName = 'med'

def set_input_dat_high(inputDict):
    inputDat.resolvedRange = [1.4953E4, 2.0E4]
    inputDat.fineLethargySpacing = 1/5000.0
    inputDat.outName = 'high'

def define_input_parser():
    dirDict = get_common_directories()
    pyDirr = dirDict['src']
    groupDirr = dirDict['dat/energy_groups']
    parser = argparse.ArgumentParser(description='Group structure specification functions.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=0, choices=[0,1,2,3,4], type=int)
    parser.add_argument('-i', '--inputdir', dest='inDirr', help='Input directory with MG library specifications', default=groupDirr)
    parser.add_argument('-I', '--inputname', dest='inName', help='Name of MG library to start with.', default='energy_in_default.txt')
    parser.add_argument('-o', '--outputdir', dest='outDirr', help='Output directory.', default=groupDirr)
    parser.add_argument('-O', '--outname', dest='outName', help="Short name for output files. A value of 'none' means do not print", default='custom')
    parser.add_argument('-d', '--defaultuse', dest='useDefault', help="Use a pre-defined input scheme. If not 'none', overrides 'rangeresolved,' 'finelethargyspacing,' and 'outname.'", choices=['none', 'low', 'med', 'high'], default='none')
    parser.add_argument('-r', '--rangeresolved', dest='resolvedRange', help='The resolved resonance range to be refined (energies in eV). If multiple ranges are desired, this program may be run multiple times.', nargs=2, type=float, default=[3.0, 1000.0])
    parser.add_argument('-R', '--rangefull', dest='fullRange', help='The full energy range for the group structure.', nargs=2, type=float, default=[1e-5, 2e7])
    parser.add_argument('-e', '--noextend', dest='noExtend', help="If specified, do not crop to 'rangefull.'", action='store_true', default=False)
    energyMeshGroup = parser.add_argument_group('Number of groups specification', 'If both number of groups and lethargy spacing are provided, the finer of the two is used.')
    energyMeshGroup.add_argument('-g', '--groups', help='Total number of groups to use.', type=int, default=1)
    energyMeshGroup.add_argument('-f', '--finelethargyspacing', dest='fineLethargySpacing', help='Lethargy spacing to use inside the resolve resonance region.', type=float, default=1e4)
    energyMeshGroup.add_argument('-b', '--boundarylethargyspacing', dest='boundaryLethargySpacing', help='Lethargy spacing for the boundary between resolved and unresolved resonance regions. The boundary acts as a buffer for downscattering from medium-to-heavy nuclei so the two regions do not see each other directly. Use a value of 0 for no buffer.', type=float, default=0.3)
    return parser

def finish_parsing_inputs(inputDat):
    if inputDat.outName.lower().strip() == 'none':
        inputDat.outName = None
    inputDat.resolvedRange = sorted(inputDat.resolvedRange)
    if inputDat.useDefault == 'low':
        set_input_dat_low(inputDat)
    elif inputDat.useDefault == 'med':
        set_input_dat_med(inputDat)
    elif inputDat.useDefault == 'high':
        set_input_dat_high(inputDat)
    if inputDat.verbosity > 2:
        for key in sorted(inputDat.__dict__):
            print key, getattr(inputDat, key)
    elif inputDat.verbosity > 1:
        print inputDat
    inputDat.finishedParsing = True

if __name__ == '__main__':
    inputDat = define_input_parser().parse_args()
    finish_parsing_inputs(inputDat)
    resolve_group_structure(inputDat)
