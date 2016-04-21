#! /usr/bin/env python

'''
Andrew Till
Spring 2016

Merge two or more PDT cross sections of the same material at different temperatures
'''

#STDLIB
import os
import argparse
#TPL
import numpy as np
#MINE
from directories import get_common_directories
import PDTXS as pdtxs

def do_all(inputDict):
    verbosity = inputDict['verbosity']
    aggregationOpt = inputDict['aggregationopt']
    materialsList = inputDict['materialslist']
    numGroups = inputDict['groups']
    inDirr = inputDict['inputdir']
    outDirr = inputDict['outputdir']
    inFormat = inputDict['inputformat']
    outFormat = inputDict['outputformat']
    # Step 0: Echo input
    if verbosity:
        print 'Running XS Merger [o] [o] > [oo]'
    if verbosity > 1:
        print '-------------------------------'
        for k,v in sorted(inputDict.items()):
            print k, ":", v
    # Step 1: Aggregate materials
    materialsSetDict = {}
    if aggregationOpt == 'all':
        materialsSetDict['all'] = set(materialsList)
    elif aggregationOpt == 'auto':
        baseNames = set([name.split('_')[0] for name in materialsList])
        for baseName in baseNames:
            matchingNames = set([name for name in materialsList if name.split('_')[0] == baseName])
            materialsSetDict[baseName] = matchingNames
    if verbosity:
        print '-------------------------------'
        print 'basename : matching_materials'
        for k,v in sorted(materialsSetDict.items()):
            print k, ":", v
    # Step 2: For each aggregate, create merged cross section
    for baseName in materialsSetDict:
        if verbosity:
            print '-------------------------------'
        matchingNames = materialsSetDict[baseName]
        # Step 2a: Read input cross sections
        xsDict = {}
        for materialName in sorted(matchingNames):
            filename = inFormat.format(m=materialName, g=numGroups)
            inPath = os.path.join(inDirr, filename)
            if verbosity:
                print 'Reading {}'.format(inPath)
            xsDict[materialName] = pdtxs.read_PDT_xs_generally(inPath)
        # Step 2b: Make sure the cross sections have the same information
        cmpXS = xsDict[materialName]
        cmpNumGroups = cmpXS.G
        cmpNumMoments = cmpXS.M
        cmpGroupBounds = cmpXS.Eg
        cmpPDTMTs = set(cmpXS.xs.keys())
        for materialName in matchingNames:
            tstXS = xsDict[materialName]
            tstNumGroups = tstXS.G
            tstNumMoments = tstXS.M
            tstGroupBounds = tstXS.Eg
            tstPDTMTs = set(tstXS.xs.keys())
            if ((tstNumGroups != cmpNumGroups) or (tstNumMoments != cmpNumMoments) or 
                np.any(tstGroupBounds != cmpGroupBounds) or (tstPDTMTs != cmpPDTMTs)):
                print "ERROR! The XS of material {} does not match other XS to be merged. Quiting".format(materialName)
                exit(1)
        # Step 2c: Get dictionary of temperatures for materials in the aggregate
        # If two materials have the same temperature, the last one alphabetically is used
        temperatureToNameDict = {}
        for materialName in sorted(matchingNames):
            temperatureToNameDict[xsDict[materialName].T] = materialName
        temperatureList = sorted(temperatureToNameDict.keys())
        numMergedXS = len(temperatureList)
        # Step 2d: Write header of merged XS
        filename = outFormat.format(b=baseName, g=numGroups, n=numMergedXS)
        outPath = os.path.join(outDirr, filename)
        if verbosity:
            print 'Writing {}'.format(outPath)
        pdtxs.write_PDT_xs_header(outPath, cmpXS, temperatureList)
        # Step 2e: Write body of merged XS
        for T in temperatureList:
            materialName = temperatureToNameDict[T]
            xs = xsDict[materialName]
            pdtxs.write_PDT_xs_body(outPath, xs)
            if verbosity > 1:
                print '    Adding {} at {}K'.format(materialName, T)
        

###############################################################################
def define_input_parser():
    dirDict = get_common_directories()
    thisDirr = os.path.abspath('.')
    xsInDirr = dirDict['pdtxs']
    xsOutDirr = dirDict['pdtxs']
    #
    parser = argparse.ArgumentParser(description='Cross section merger', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=0, choices=[0,1,2,3,4], type=int)
    parser.add_argument('-m', '--materialslist', help='List of materials to use (for more advanced options, see materials_materials.py).', nargs='+', default=[])
    parser.add_argument('-a', '--aggregationopt', help="Which cross sections to merge. 'auto' means combine all cross sections with the same materialName but different temperatures, assuming materials are named like '{materialName}_{temperatureIndex}. 'all' means combine all cross sections in the materials list into one new cross section.", default='auto', choices=['auto', 'all'])
    parser.add_argument('-g', '--groups', help="Number of groups to use when determining the input name.", type=int, default=0)
    parser.add_argument('-i', '--inputdir', help='Input directory with PDT XS.', default=xsInDirr)
    parser.add_argument('-o', '--outputdir', help='Output directory', default=xsOutDirr)
    parser.add_argument('-I', '--inputformat', help='Format for determining the cross section. {m} is replaced with the material name. {g} is replaced by the number of groups', default='xs_{m}_{g}.data')
    parser.add_argument('-O', '--outputformat', help="Format for determining the output cross section. {b} is replaced with the base material name (no temperature index; for 'all' aggregation, base name is 'all'). {g} is replaced by the number of groups. {n} is replaced by the number of XS to be merged", default='xs_{b}_{n}T_{g}.data')
    return parser

if __name__=='__main__':
    inputDict = vars(define_input_parser().parse_args())
    do_all(inputDict)
