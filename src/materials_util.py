'''
Andrew Till
Summer 2014

Utility functions for materials.
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
    if fuelRadius:
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

def get_thermal_names():
    return ['endf_th_01_H2O_vii1', 'endf_th_01_ZrH_vii1', 'endf_th_06_C_vii1', 'endf_th_08_UO2_vii1', 'endf_th_13_Al_vii1', 'endf_th_40_ZrH_vii1', 'endf_th_92_UO2_vii1']

def get_thermal_name2filename_dict():
    return {'free': '', 'hh2o': 'endf_th_01_H2O_vii1', 'hzrhinel': 'endf_th_01_ZrH_vii1', 'graphinel': 'endf_th_06_C_vii1', 'ouo2inel': 'endf_th_08_UO2_vii1', 'alinel': 'endf_th_13_Al_vii1', 'feinel': 'endf_th_26_Fe_vii1', 'zrzrhinel': 'endf_th_40_ZrH_vii1', 'uuo2inel': 'endf_th_92_UO2_vii1'}

def get_inelastic_thermal_mt_list():
    return set([221, 222, 223, 225, 227, 228, 229, 231, 233, 235, 237, 239, 241, 243, 245])

def get_thermal_short2mat_dict():
    '''Return a dict that maps thermal XS shortname to thermal mat. See THERMR in NJOY manual.'''
    return {'free': 0, 'hh2o': 1, 'hzrhinel': 7, 'graphinel': 31, 'ouo2inel': 75, 'alinel': 45, 'zrzrhinel': 58, 'uuo2inel': 48}

###############################################################################
def print_newline(verbosity=False):
    if verbosity:
        print ''
