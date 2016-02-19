#! /usr/bin/env python

'''
Andrew Till
Summer 2014

Materials specification subroutines.

Uses nuclide_data python module by Joshua Hykes that has been modified to use python module uncertainties version 2. The data is based on NIST data.

Dependencies: numpy, uncertainties (ver 2), nuclide_data (https://github.com/jhykes/nuclide-data; updated to work with ver 2 of uncertainties), directories (dirDict), Readgroupr (to read potential xs and define some mt2name dicts), iapws (https://github.com/jjgomera/iapws/; for water densities)
'''

#TPL
import os
import argparse
#MINE
from directories import get_common_directories
import materials_njoy as njoy
import materials_materials as mat
import materials_util as util
import materials_global as glob
import materials_bondarenko as bon

def do_all(inputDict):
    if not 'finishedParsing' in inputDict:
        raise ValueError('Inputs have not finished being parsed. Use finish_parsing_inputs().')
    verbosity = inputDict['verbosity']
    useCommonBXS = inputDict['commonbackgroundxs']
    materials = []
    mat2funcDict = mat.get_materials_name2function_dict()
    for matl in inputDict['materialslist']:
        materials.append(mat2funcDict[matl]())
    #
    util.print_newline(verbosity)
    globalTDict, globalBXSDict, globalTXSDict = {}, {}, {}
    globalZAList, globalZList = glob.get_union_parameters(
        materials, globalTDict, globalBXSDict, globalTXSDict, not(useCommonBXS), verbosity)
    mat.print_materials(materials, verbosity)
    if inputDict['njoy']:
        njoyTDict, njoyBXSDict = {}, {}
        njoy.get_njoy_temperatures(globalZAList, njoyTDict, globalTDict, globalTXSDict)
        njoy.get_njoy_background_xs(globalZAList, njoyBXSDict, globalBXSDict, useCommonBXS)
        glob.print_globals(globalZAList, globalZList, globalTDict, globalBXSDict, globalTXSDict, verbosity)
        njoy.print_njoys(njoyTDict, njoyBXSDict, verbosity)
        njoy.create_njoy_decks(inputDict, globalZAList, njoyTDict, njoyBXSDict, globalTXSDict, verbosity)
    elif inputDict['bondarenko']:
        bon.perform_bondarenko_iterations(inputDict, materials, verbosity)
    util.print_newline(verbosity)

###############################################################################
def define_input_parser():
    parser = argparse.ArgumentParser(description='Materials specification functions.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=0, choices=[0,1,2,3,4], type=int)
    actionGroup = parser.add_mutually_exclusive_group()
    actionGroup.add_argument('-n', '--njoy', help="Create NJOY input decks. If neither '-n' nor '-b' is specified, do this.", action='store_true')
    actionGroup.add_argument('-b', '--bondarenko', help='Run Bondarenko iteration.', action='store_true')
    # Both options
    parser.add_argument('-m', '--materialslist', help='List of materials to use (for more advanced options, see materials_materials.py).', nargs='+', default=['hpu'], choices=mat.get_materials_name2function_dict().keys())
    # NJOY options
    njoyGroup = parser.add_argument_group('NJOY options', 'Options to use with --njoy.')
    njoyGroup.add_argument('-c', '--commonbackgroundxs', help='Do not use common background XS if specified.', action='store_false', default=True)
    njoyGroup.add_argument('-g', '--groups', help='Number of groups for NJOY calculation. Unless the groupopt is an NJOY groupset, the final number of groups will be the maximum of this argument and the native number of groups in the specified groupopt. Use 0 to specify no resolution.', type=int, default=0)
    njoyGroup.add_argument('-G', '--groupopt', help="Base coarse group structure for NJOY calculation. 'default' is read in default for makegroups and refine. 'equal' means do equal lethargy spacing with 'groups' groups. 'time' means do equal spacing of 1/v_g - 1/v_{g-1} for time-of-flight fidelity. For the NJOY options, no resolving is done. 'njoy-10' means use LANL 187-group structure. 'njoy-15' means use Sandia 640-group structure. Anything else is assumed to a filename. Some examples include 'c5g7', 'c5g7_alt', 'scale-44', 'helios-47' (just resonance groups), 'wims-69', 'wims-172', 'scale-238', 'low-1500' 'med-1500', 'high-1500', and 'res-N' where N=1,..9. Always looks in default directory of makegroups.", default='default')
    njoyGroup.add_argument('-E', '--energybounds', help='Low and high bounds for the cross sections (eV). The interpretation of this value depends on the value of groupopt.', nargs=2, type=float, default=[1E-5,2E7])
    njoyGroup.add_argument('-r', '--resolvedrange', help='The resolved resonance range to be refined (energies in eV). The interpretation of this value depends on the value of groupopt.', nargs=2, type=float, default=[3.0, 2.5E4])
    njoyGroup.add_argument('-L', '--legendreorder', help='Legendre order to use for the scattering kernel.', type=int, default=0)
    njoyGroup.add_argument('-S', '--smallscattering', help='Do not build transfer matrices when running GROUPR if specified.', action='store_true', default=False)
    # Bondarenko options
    bondarenkoGroup = parser.add_argument_group('Bondarenko options', 'Options to use with --bondarenko.')
    bondarenkoGroup.add_argument('-N', '--numberiterationsmax', help='Maximum number of Bondarenko iterations. If a value less than 0 is used, unshielded (infinite) background cross sections are used.', type=int, default=10)
    bondarenkoGroup.add_argument('-e', '--errormax', help='Maximum L^infty error for Bondarenko iterations.', type=float, default=1e-10)
    bondarenkoGroup.add_argument('-f', '--format', help='Output format for scattering matrices column or row major.', choices=['csr', 'csc'], default='csr')
    bondarenkoGroup.add_argument('-s', '--simplereactions', help='Use non-threshold reactions and free thermal scattering cross sections only.', action='store_true', default=False)
    bondarenkoGroup.add_argument('-p', '--printopt', help='Which XS to print. Usual prints selected XS and the combined transfer matrix, invel prints the usual in addition to the inverse velocity, abs prints the usual in addition to some cross sections necessary to determine the total absorption rate, total just prints the flux and total XS. All prints all XS and transfer matrices except the combined transfer matrix. None prints no PDT XS.', choices=['usual', 'invel', 'abs', 'total', 'all', 'none'], default='usual')
    bondarenkoGroup.add_argument('-M', '--mesh', help='File name for the mapping of fine groups to discontiguous energy elements. If no directory specified, looks in $SCRATCH_BARN/dat/energy_groups.  If no filetype specified, use .txt', default='none')
    bondarenkoGroup.add_argument('-F', '--fluxes', help='File name pattern for the material fluxes. {m} is replaced by the shortname of the material. If no directory specified, looks in $SCRATCH_BARN/dat/indicators. If no filetype specified, use .txt', default='wgt_{m}_9')
    return parser

def finish_parsing_inputs(inputDict):
    dirDict = get_common_directories()
    if inputDict['mesh'] in ['none', 'None']:
        inputDict['mesh'] = None
        inputDict['fluxes'] = None
    else:
        if not inputDict['mesh'].count('.'):
            inputDict['mesh'] += '.txt'
        if not inputDict['fluxes'].count('.'):
            inputDict['fluxes'] += '.txt'
        if not os.path.split(inputDict['mesh'])[0]:
            inputDict['mesh'] = os.path.join(dirDict['dat/energy_groups'], inputDict['mesh'])
        if not os.path.split(inputDict['fluxes'])[0]:
            inputDict['fluxes'] = os.path.join(dirDict['dat/indicators'], inputDict['fluxes'])
    if not(inputDict['njoy'] or inputDict['bondarenko']):
        inputDict['njoy'] = True
    inputDict['materialslist'] = set(inputDict['materialslist'])
    if inputDict['energybounds'][1] > 2E7:
        print 'Current ENDF data does not exist above 20 MeV. Changing energy upperbound to this value.'
        inputDict['energybounds'][1] = 2E7
    inputDict['energybounds'][0] = max(0, inputDict['energybounds'][0])
    if inputDict['verbosity'] > 1:
        print inputDict
    inputDict['finishedParsing'] = True

if __name__=='__main__':
    inputDict = vars(define_input_parser().parse_args())
    finish_parsing_inputs(inputDict)
    do_all(inputDict)
