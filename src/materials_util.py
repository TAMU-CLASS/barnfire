'''
Andrew Till
Summer 2014

Utility functions for materials.

References:
mcnp5: MCNP 5 Manual, Volume I, LA-UR-03-1987, Appendix G (MCNP Data Libraries, including S(alpha,beta))
mcnpdata: Listing of Available ACE Data Tables [for MCNP6], LA-UR-13-21822 (An updated version of the above)
njoy2016: The NJOY Nuclear Data Processing System, Version 2016, LA-UR-17-20093 (NJOY 2016 manual)
endf6: ENDF-6 Formats Manual, BNL-90365-2009 Rev.2 (Description of the ENDF-6 format)

Naming conventions:
'Thermal name' is specified in materials_materials.py and applies to a material (e.g., ZrH is 'zrh').
'Element thermal name' includes the (bound) thermal treatment of the material and the applicable
    element (e.g., H in Zr is 'hzrh'). Uses Hollerith strings specified in Table 21 of MATXSR chapter
    in NJOY manual (ref: njoy2016).
    Elsewhere, this is called Sab.
'Thermal XS name' includes the (bound) thermal treatment of the material, the applicable element,
    and the applicable type of thermal XS (e.g., the inelastic bound thermal XS for H in ZrH is 'hzrhinel')
'MCNP thermal name' is the string used by MCNP to refer to a bound thermal treatment of a material along
    with a consituent element (e.g., the bound thermal XS for H in ZrH is 'h/zr')
NB: 'free' is used for thermal name, element thermal name, and thermal xs name for free-gas thermal
    treatment
NB: ''/None/'none' is used for thermal name; 'none' is used for element thermal name; [] is used
    for the thermal xs list (for no thermal treatment)
'''

#TPL
import numpy as np

###############################################################################
def get_nearest_points(value, sortedArray):
    '''If value in sortedArray, return [value], else return two nearest points in sortedArray'''
    nearestIndex = np.argmin(np.abs(sortedArray - value))
    nearestPoints = [sortedArray[nearestIndex]]
    lastIndex = len(sortedArray) - 1
    if value < sortedArray[nearestIndex] and nearestIndex != 0:
        nearestPoints.append(sortedArray[nearestIndex-1])
    elif value > sortedArray[nearestIndex] and nearestIndex != lastIndex:
        nearestPoints.append(sortedArray[nearestIndex+1])
    return nearestPoints

def thin_list(inList, useLinSpacing=True):
    '''Thin a list by returning the value of the point which is nearest to its two neighbors.'''
    sortedList = np.array(sorted(inList))
    # Use [2:] and [:-2] to get sum of distances to neighbors. Add 1 for the indexing to be correct.
    if useLinSpacing:
        distance = sortedList[2:] - sortedList[:-2]
    else:
        distance = sortedList[2:] / sortedList[:-2]
    index = np.argmin(distance) + 1
    return sortedList[index]

###############################################################################
def calc_chord_length(fuelRadius):
    if fuelRadius == 'unshielded':
        return 1.e10
    elif fuelRadius:
        surfaceArea = 2 * np.pi * fuelRadius
        volume = np.pi * fuelRadius * fuelRadius
        chordLength = surfaceArea / (4 * volume)
        return chordLength
    else:
        return 0.0

def has_bondarenko_iteration(Z):
    if Z <= 5:
        return False
    else:
        return True

def is_fissionable((Z,A)):
    if Z >= 89:
        return True
    elif (Z,A) in [(88, 223), (88, 226)]:
        return True
    else:
        return False

def avogadros_number():
    '''In units of atoms / mole, but multiplied by 1E-24'''
    return 0.60221413

def get_inelastic_thermal_mt_list():
    '''Ref: Culled from Table 4 in the NJOY manual (ref: njoy2016)'''
    return set([221, 222, 223, 225, 227, 228, 229, 231, 233, 235, 237, 239, 241, 243, 245])

###############################################################################
def format_zaid():
    '''Use format_ZAID()(**dict) to apply to a dictionary'''
    return '{Z:d}{A:03d}'.format

def format_zaid_leading_zeros():
    return '{Z:02d}{A:03d}'.format

def format_thermal_filename():
    return 'endf_th_{Z:02d}_{thermalName}_vii1'.format

def get_nuclide_dirr(sym, A, elementThermalName, metastableStr=''):
    '''Returns the directory name for a nuclide'''
    nuclideName = '{0}-{1}{2}'.format(sym.lower(), A, metastableStr)
    et2t = get_element_thermal_name_to_thermal_name_dict()
    if elementThermalName in et2t:
        thermalName = et2t[elementThermalName]
        nuclideName = '{0}-{1}'.format(nuclideName, thermalName)
    return nuclideName

def get_ace_extension(elementThermalName):
    '''Returns the ACE extension given an element thermal name. The extension is 9 (as in .90c)
    if the *nuclide* does not have a bound thermal treatment.'''
    if elementThermalName in get_non_bound_names():
        return 9
    else:
        et2t = get_element_thermal_name_to_thermal_name_dict()
        t2ext = get_thermal_name_to_ace_ext_dict()
        return t2ext[et2t[elementThermalName]]

###############################################################################
def get_thermal_name_to_element_thermal_name_dict():
    '''Returns a dict that maps tuple (Z, thermal_name) to element thermal name (a Hollerith string).
    Free thermal treatment is not in dict'''
    return {(1,'h2o'): 'hh2o', (8,'uo2'): 'ouo2', (92,'uo2'): 'uuo2', (1,'zrh'): 'hzrh', (40,'zrh'): 'zrzrh', (6,'graphite'): 'graph', (13,'al'): 'al', (26,'fe'): 'fe'}

def get_element_thermal_name_to_nuclide_list_dict():
    '''Returns a dict that maps element thermal name to the (Z,A)'s of the nuclides to which it applies.
    Free thermal treatment is not in dict. Ref: mcnp5, mcnpdata'''
    return {'hh2o': [(1,1)], 'ouo2': [(8,16), (8,17), (8,18)], 'uuo2': [(92,238)], 'hzrh': [(1,1)], 'zrzrh': [(40,0), (40,90), (40,91), (40,92), (40,94), (40,96)], 'graph': [(6,0), (6,12)], 'al': [(13,27)], 'fe': [(26,56)]}

def get_element_thermal_name_to_inelastic_mt_number_dict():
    '''Returns a dict that maps element thermal name to the MT number corresponding to the
    inelastic thermal XS for that element thermal name. See Tables 4 and 25 in NJOY manual (ref: njoy2016).
    Should correspond to 'inel' endf numbers in Readgroupr.py's get_endf_mt_list() function.'''
    return {'free': 221, 'hh2o': 222, 'ouo2': 239, 'uuo2': 241, 'hzrh': 225, 'zrzrh': 235, 'graph': 229, 'al': 243, 'fe': 245}

def get_thermal_name_to_nuclide_list_dict():
    '''Returns a dict that maps thermal name to the (Z,A)'s of the nuclides to which it applies.
    Free thermal treatment is not in dict'''
    # Derive this dictionary from previous information
    Zthermal2elem = get_thermal_name_to_element_thermal_name_dict()
    elem2ZAs = get_element_thermal_name_to_nuclide_list_dict()
    # First, initialize each thermal name as an empty set of nuclides
    thermal2ZAs = {}
    for (Z,thermalName) in Zthermal2elem:
        thermal2ZAs[thermalName] = set()
    # Then, populate nuclides that correspond to that thermal name
    for (Z,thermalName) in Zthermal2elem:
        elem = Zthermal2elem[(Z,thermalName)]
        ZAList = elem2ZAs[elem]
        thermal2ZAs[thermalName].update(ZAList)
    return thermal2ZAs

def get_element_thermal_name_to_thermal_name_dict():
    '''Returns a dict that maps element thermal name (e.g., hh2o) to thermal name (h2o).
    Free thermal treatment is not in dict.'''
    # Derive this dictionary from previous information
    return {et: t for (Z,t), et in get_thermal_name_to_element_thermal_name_dict().items()}

###############################################################################
def get_non_bound_names():
    '''Returns a set of thermal names that do not use bound cross sections'''
    return set(['free', '', None, 'none'])

def get_element_thermal_name_to_thermal_xs_list_dict():
    '''Returns a dict that maps element thermal name to a list of thermal XS names.
    See Tables 4 and 25 in NJOY manual (ref: njoy2016)'''
    return {'free': ['free'], 'none': [], 'hh2o': ['hh2o'], 'ouo2': ['ouo2inel', 'ouo2elas'], 'uuo2': ['uuo2inel', 'uuo2elas'], 'hzrh': ['hzrhinel', 'hzrhelas'], 'zrzrh': ['zrzrhinel', 'zrzrhelas'], 'graph': ['graphinel', 'graphelas'], 'al': ['alinel', 'alelas'], 'fe': ['feinel', 'feelas'],}

def get_element_thermal_name_to_mat_number_dict():
    '''Returns a dict that maps element thermal name to thermal MAT number.
    See Table 4 in THERMR chapter of NJOY manual (ref: njoy2016).
    Change: Manuals say Al is material 45, but the data file has material 53'''
    return {'free': 0, 'none': 0, 'hh2o': 1, 'ouo2': 75, 'uuo2': 48, 'hzrh': 7, 'zrzrh': 58, 'graph': 31, 'al': 53, 'fe': 56}

def get_element_thermal_name_to_bnl_id_dict():
    '''Returns a dict that maps element thermal name to BNL's ID. Used in URL for automatic downloads.
    Free thermal treatment not in dict. Last verified: 05/31/2016.'''
    return {'hh2o': 15391, 'ouo2': 15395, 'uuo2': 15402, 'hzrh': 15392, 'zrzrh': 15403, 'graph': 15389, 'al': 15383, 'fe': 15384}

def get_element_thermal_name_to_endf_filename_dict():
    '''Returns a dict that maps element thermal name to bound thermal ENDF file name'''
    # Derive this dictionary from previous information
    ZName2elem = get_thermal_name_to_element_thermal_name_dict()
    ZName2filename = format_thermal_filename()
    elem2filename = {}
    for (Z,thermalName) in ZName2elem:
        elem = ZName2elem[(Z,thermalName)]
        elem2filename[elem] = ZName2filename(Z=Z, thermalName=thermalName)
    for key in ['free', None]:
        elem2filename[key] = None
    return elem2filename

###############################################################################
def get_thermal_name_to_ace_ext_dict():
    '''Returns a dict that maps thermal name to ACE extension (e.g., the 9 in 1001.92c).
    Reasonable extension options include poly (CH2) being 6, benzene (C6H6) being 7, and d2o (D2O) being 8, though these are currently unsupported'''
    return {'h2o': 0, 'uo2': 1, 'zrh': 2, 'graphite': 3, 'al': 4, 'fe': 5, 'free': 9}

def get_element_thermal_name_to_mcnp_thermal_name_dict():
    '''Returns a dict that maps element thermal name to MCNP thermal material name. Ref: mcnp5, mcnpdata'''
    return {'hh2o': 'lwtr', 'zrzrh': 'zr/h', 'hzrh': 'h/zr', 'uuo2': 'u/o2', 'ouo2': 'o2/u', 'graph': 'grph', 'fe': 'fe56', 'al': 'al27'}

def get_mcnp_thermal_name_to_main_za_dict():
    '''Returns a dict that maps MCNP bound thermal material name to the (Z,A) of the nuclide from which to read PENDF tape (for ACE only). Uses nuclide with highest atom fraction.'''
    return {'lwtr': (1,1), 'o2/u': (8,16), 'u/o2': (92,238), 'grph': (6,0), 'h/zr': (1,1), 'zr/h': (40,90), 'al27': (13,27), 'fe56': (26,56), }

def get_mcnp_thermal_name_to_zaid_list_dict():
    '''Returns a dict that maps MCNP bound thermal material name to a ZZAAA string of the nuclides that use this thermal option. Ref: mcnp5, mcnpdata'''
    # Derive this dictionary from previous information
    elem2mcnp = get_element_thermal_name_to_mcnp_thermal_name_dict()
    elem2ZAs = get_element_thermal_name_to_nuclide_list_dict()
    ZA2zaid = format_zaid()
    mcnp2zaid = {}
    for elem in elem2mcnp:
        mcnpName = elem2mcnp[elem]
        ZAList = elem2ZAs[elem]
        zaidList = [ZA2zaid(Z=Z, A=A) for (Z,A) in ZAList]
        mcnp2zaid[mcnpName] = zaidList
    return mcnp2zaid

###############################################################################
def print_newline(verbosity=False):
    if verbosity:
        print ''
