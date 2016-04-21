'''
Andrew Till
Summer 2014

Materials class and example materials.
'''

#MINE
from directories import get_common_directories

#STDLIB
import os
import sys
import copy
sys.path.insert(1, get_common_directories()['nuclideData'])

#TPL
import numpy as np
import nuclide_data as nd
from iapws import IAPWS97 as steam

#MINE
from materials_util import calc_chord_length
import materials_util as util

###############################################################################
def get_materials_name2function_dict():
    return {
        # Simple examples
        'hpu': get_hpu_slurry_material,
        'hheu': get_hheu_slurry_material,
        'hleu': get_hleu_slurry_material,
        'puMetal': get_pu_metal_material,
        'puMetalHot': get_hot_pu_metal_material,
        'uMetal': get_u_metal_material,
        # Simple multi-temperature examples
        'uo2ColdPin': get_cold_pin_uo2_material,
        'uo2Cold': get_cold_uo2_material,
        'uo2InnerHot': get_inner_hot_uo2_material,
        'uo2MiddleHot': get_middle_hot_uo2_material,
        'uo2OuterHot': get_outer_hot_uo2_material,
        'moxCold': get_cold_mox_material,
        'moxInnerHot': get_inner_hot_mox_material,
        'moxMiddleHot': get_middle_hot_mox_material,
        'moxOuterHot': get_outer_hot_mox_material,
        'h2oCold': get_cold_h2o_material,
        'h2oHot': get_hot_h2o_material,
        'graphite': get_graphite_material,
        # C5G7
        'clowMOX': get_c5g7_low_mox_material,
        'cmedMOX': get_c5g7_med_mox_material,
        'chighMOX': get_c5g7_high_mox_material,
        'cUO2': get_c5g7_uo2_material,
        'cMOD': get_c5g7_moderator_material,
        'cGUIDE': get_c5g7_guide_tube_material,
        'cFCHAMBER': get_c5g7_fission_chamber_material,
        'cCR': get_c5g7_control_rod_material,
        # TRIGA (BOL)
        'tFUEL': get_triga_fuel_material,
        'tZIRC': get_triga_zirconium_material,
        'tCLAD': get_triga_clad_material,
        'tMOD': get_triga_moderator_material,
        'tBORATEDGRAPHITE': get_triga_borated_graphite_material,
        'tB4C': get_triga_b4c_material,
        'tGRAPHITE': get_triga_graphite_material,
        'tAIRTUBE': get_triga_air_tube_material,
        'tIRRADIATIONTUBE': get_triga_irradiation_tube_material,
        'tGRIDPLATE': get_triga_grid_plate_material,
        # Iron (for time-dependent dissertation problem)
        'iron': get_iron_material,
        'thickiron': get_thick_iron_material,
        # Simplified PWR pincell (for Kord Smith)
        'kFUEL': get_kord_fuel_material,
        'kRFUEL': get_kord_rod_fuel_material,
        'kEFUEL': get_kord_enriched_fuel_material,
        'kCLAD': get_kord_clad_material,
        'kZR': get_kord_zirconium_material,
        'kMOD': get_kord_moderator_material,
        'kREFUEL': get_kord_enriched_rod_fuel_material,
        'kRMFUEL': get_kord_mox_rod_fuel_material,
        # UO2 (for Don Bruss)
        'debFUEL': get_bruss_enriched_rod_fuel_material,
        # Multi-temperature examples
        'mtH2O_0': get_multi_temperature_h2o_material_T0,
        'mtH2O_1': get_multi_temperature_h2o_material_T1,
        'mtH2O_2': get_multi_temperature_h2o_material_T2,
        'mtTFUEL_0': get_multi_temperature_triga_fuel_material_T0,
        'mtTFUEL_1': get_multi_temperature_triga_fuel_material_T1,
        'mtTFUEL_2': get_multi_temperature_triga_fuel_material_T2,
        'mtTFUEL_3': get_multi_temperature_triga_fuel_material_T3,
        'mtTFUEL_4': get_multi_temperature_triga_fuel_material_T4,
        'mtTFUEL_5': get_multi_temperature_triga_fuel_material_T5,
        'mtTFUEL_6': get_multi_temperature_triga_fuel_material_T6,
        'mtTFUEL_7': get_multi_temperature_triga_fuel_material_T7,
    }

###############################################################################

def get_hpu_slurry_material():
    shortName = 'hpu'
    longName = 'Pu-H slurry'
    massDensity = 5. #g/cc
    fuelRadius = 0. #cm
    temperature = 301. #K
    thermalOpt = 'free'
    puAtomFractionsDict = {239: 1.0}
    elemAtomFracDict = {'Pu': 0.07, 'H': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_hheu_slurry_material():
    shortName = 'hheu'
    longName = 'HEU-H slurry'
    massDensity = 5. #g/cc
    fuelRadius = 0. #cm
    temperature = 301. #K
    thermalOpt = 'free'
    uAtomFractionsDict = {235: 1.0}
    elemAtomFracDict = {'U': 0.07, 'H': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_hleu_slurry_material():
    shortName = 'hleu'
    longName = 'LEU-H slurry'
    massDensity = 5. #g/cc
    fuelRadius = 0. #cm
    temperature = 301. #K
    thermalOpt = 'free'
    uAtomFractionsDict = {235: 0.07, 238: 0.93}
    elemAtomFracDict = {'U': 1.0, 'H': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_pu_metal_material():
    shortName = 'Pu metal'
    longName = 'Pu metal'
    massDensity = 19.8 #g/cc
    fuelRadius = 0.0 #cm
    temperature = 301. #K
    thermalOpt = 'none'
    puAtomFractionsDict = {239: 1.0}
    elemAtomFracDict = {'Pu': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_hot_pu_metal_material(T=400.):
    material = get_pu_metal_material()
    material.update_temperature(T)
    material.update_names('puMetalHot', 'hot Pu metal')
    return material

def get_u_metal_material():
    shortName = 'U metal'
    longName = 'U metal'
    massDensity = 19.1 #g/cc
    fuelRadius = 0.0 #cm
    temperature = 301. #K
    thermalOpt = 'none'
    uAtomFractionsDict = {235: 0.05, 238: 0.95}
    elemAtomFracDict = {'U': 1}
    #
    chordLength = util.calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_uo2_material():
    shortName = 'UO2'
    longName = 'uranium dioxide'
    massDensity = 10.97 #g/cc
    fuelRadius = 0.47 #cm
    #fuelRadius = 0.0 #cm
    temperature = 301. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {235: 0.05, 238: 0.95}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_cold_pin_uo2_material(T=400.):
    material = get_uo2_material()
    material.update_temperature(T)
    #material.update_mass_density(material.massDensity / 1.1)
    material.update_names('uo2PinCold', 'cold uranium dioxide')
    return material

def get_cold_uo2_material(T=400.):
    material = get_uo2_material()
    material.update_temperature(T)
    material.update_names('uo2Cold', 'cold uranium dioxide')
    return material

def get_inner_hot_uo2_material(T=1000.):
    material = get_uo2_material()
    material.update_temperature(T)
    material.update_names('uo2InnerHot', 'inner hot uranium dioxide')
    return material

def get_middle_hot_uo2_material(T=800.):
    material = get_uo2_material()
    material.update_temperature(T)
    material.update_names('uo2MiddleHot', 'middle hot uranium dioxide')
    return material

def get_outer_hot_uo2_material(T=700.):
    material = get_uo2_material()
    material.update_temperature(T)
    material.update_names('uo2OuterHot', 'outer hot uranium dioxide')
    return material

def get_mox_material():
    shortName = 'MOX'
    longName = 'mixed oxide'
    massDensity = 10.97 #g/cc
    fuelRadius = 0.47 #cm
    #fuelRadius = 0.0 #cm
    temperature = 301. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {238: 1.0}
    puAtomFractionsDict = {239: 1.0}
    elemAtomFracDict = {'O': 2, 'U': 0.95, 'Pu': 0.05}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_cold_mox_material(T=400.):
    material = get_mox_material()
    material.update_temperature(T)
    material.update_names('moxCold', 'cold mixed oxide')
    return material

def get_inner_hot_mox_material(T=1000.):
    material = get_mox_material()
    material.update_temperature(T)
    material.update_names('moxInnerHot', 'inner hot mixed oxide')
    return material

def get_middle_hot_mox_material(T=800.):
    material = get_mox_material()
    material.update_temperature(T)
    material.update_names('moxMiddleHot', 'middle hot mixed oxide')
    return material

def get_outer_hot_mox_material(T=700.):
    material = get_mox_material()
    material.update_temperature(T)
    material.update_names('moxOuterHot', 'outer hot mixed oxide')
    return material

def get_graphite_material():
    shortName = 'graphite'
    longName = 'graphite'
    massDensity = 1.72 #g/cc
    fuelRadius = 0.0 #cm
    temperature = 296. #K
    thermalOpt = 'graphite'
    elemAtomFracDict = {'C': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances_as_elemental(ZAList, abundanceDict, 'C')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_thick_iron_material():
    shortName = 'thickiron'
    longName = 'thick iron'
    massDensity = 7.874 #g/cc
    fuelRadius = 1E-12 #cm (to give a large escape cross section)
    temperature = 296. #K
    thermalOpt = 'none'
    elemAtomFracDict = {'Fe': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_iron_material():
    shortName = 'iron'
    longName = 'iron'
    massDensity = 7.874 #g/cc
    # Unknown how I came about this:
    #fuelRadius = 13.41 #cm
    # Alternatively, we could use:
    #fuelRadius = 1E-12 #cm (make escape cross section large)
    # The iron reaction rate is proportional to (1-exp(-SigmaT*X)), where X is the iron thickness.
    # If we use SigmaEscape = 1/X, and our weighting spectrum is 1/(1+SigmaT*X),
    # then the energy shape of our iron reaction rate is (approximately):
    # SigmaT / (SigmaEscape + SigmaT) = SigmaT*X / (1 + SigmaT*X) ~= (1-exp(-SigmaT*X)).
    fuelRadius = 2.06 / 2 #cm (factor of 2 is to compensate for chord length being defined for cylinders)
    temperature = 296. #K
    thermalOpt = 'none'
    elemAtomFracDict = {'Fe': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material


def get_h2o_material():
    shortName = 'H2O'
    longName = 'light water'
    massDensity = 1.0 #g/cc
    fuelRadius = 0.47 #cm
    #fuelRadius = 0.0 #cm
    temperature = 293.6 #K
    thermalOpt = 'h2o'
    elemAtomFracDict = {'H': 2, 'O': 1}
    #
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    chordLength = calc_chord_length(fuelRadius)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material
    
def get_cold_h2o_material(T=400.):
    material = get_h2o_material()
    material.update_temperature(T)
    material.update_mass_density(0.749)
    material.update_names('h2oCold', 'cold light water')
    return material

def get_hot_h2o_material(T=550.):
    material = get_h2o_material()
    material.update_temperature(T)
    material.update_mass_density(0.749)
    material.update_names('h2oHot', 'hot light water')
    return material

def get_general_cold_h2o_material(T=293.6, P=0.1):
    #P in MPa, T in K, rho in g/cc
    rho = steam(P=P, T=T).rho / 1000
    material = get_h2o_material()
    material.update_temperature(T)
    material.update_mass_density(rho)
    material.update_chord_length(0.)
    material.update_names('h2oAtmCold', 'cold light water at atmospheric pressure')
    return material

def get_general_hot_h2o_material(T=570., P=15.5):
    #T: 275 - 315 C => 548 - 588 K
    #P in MPa, T in K, rho in g/cc
    rho = steam(P=P, T=T).rho / 1000
    material = get_h2o_material()
    material.update_temperature(T)
    material.update_mass_density(rho)
    material.update_chord_length(0.)
    material.update_names('h2oRxtHot', 'hot light water at high pressure')
    return material

###############################################################################
def get_c5g7_low_mox_material():
    shortName = 'clowMOX'
    longName = 'homogenized 4.3% MOX'
    atomDensity = 5.90443E-2
    fuelRadius = 0.54 #cm
    temperature = 296. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {235: 2.25734E-03, 238: 9.97743E-01}
    puAtomFractionsDict = {238: 1.51976E-02, 239: 5.87640E-01, 240: 2.43161E-01, 241: 9.92908E-02, 242: 5.47113E-02}
    amAtomFractionsDict = {241: 1.}
    elemAtomFracDict = {'U': 0.21573257, 'Pu': 0.009613, 'Am': 0.00012662, 'O': 0.45094438, 'Zr': 0.12712441, 'Al': 0.19645902}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    override_abundances(ZAList, abundanceDict, amAtomFractionsDict, 'Am')
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_med_mox_material():
    shortName = 'cmedMOX'
    longName = 'homogenized 7.0% MOX'
    atomDensity = 5.93894E-02
    fuelRadius = 0.54 #cm
    temperature = 296. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {235: 2.25734E-03, 238: 9.97743E-01}
    puAtomFractionsDict = {238: 1.51899E-02, 239: 5.88608E-01, 240: 2.46836E-01, 241: 9.62026E-02, 242: 5.31646E-02}
    amAtomFractionsDict = {241: 1.}
    elemAtomFracDict = {'U': 0.2144792, 'Pu': 0.01529919, 'Am': 0.00019366, 'O': 0.44832447, 'Zr': 0.12638584, 'Al': 0.19531763}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    override_abundances(ZAList, abundanceDict, amAtomFractionsDict, 'Am')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_high_mox_material():
    shortName = 'chighMOX'
    longName = 'homogenized 8.7% MOX'
    atomDensity = 5.96194E-02
    fuelRadius = 0.54 #cm
    temperature = 296. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {235: 2.25734E-03, 238: 9.97743E-01}
    puAtomFractionsDict = {238: 1.51899E-02, 239: 5.87342E-01, 240: 2.48101E-01, 241: 9.62025E-02, 242: 5.31645E-02}
    amAtomFractionsDict = {241: 1.}
    elemAtomFracDict = {'U': 0.21365168, 'Pu': 0.01905021, 'Am': 0.00024114, 'O': 0.44659472, 'Zr': 0.12589821, 'Al': 0.19456404}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    override_abundances(ZAList, abundanceDict, amAtomFractionsDict, 'Am')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_uo2_material():
    shortName = 'cUO2'
    longName = 'homogenized 3.7% UO2'
    atomDensity = 5.89782E-02
    fuelRadius = 0.54 #cm
    temperature = 296. #K
    thermalOpt = 'uo2'
    uAtomFractionsDict = {235: 3.74216E-02, 238: 9.62578E-01}
    elemAtomFracDict = {'U': 0.22538375, 'O': 0.45066999, 'Zr': 0.12726695, 'Al': 0.19667931}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_moderator_material():
    shortName = 'cMOD'
    longName = 'borated light water'
    atomDensity = 1.00528E-01
    fuelRadius = 0.54 #cm
    temperature = 293.6 #K
    thermalOpt = 'h2o'
    elemAtomFracDict = {'O': 0.33324115, 'H': 0.66648231, 'B': 2.765404E-04}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_guide_tube_material():
    shortName = 'cGUIDE'
    longName = 'homogenized guide tube'
    atomDensity = 7.60666E-02
    fuelRadius = 0.54 #cm
    temperature = 293.6 #K
    thermalOpt = 'h2o'
    elemAtomFracDict = {'O': 0.17459076, 'H':  0.34918152, 'B': 0.00014488, 'Al': 0.47608284}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_fission_chamber_material():
    shortName = 'cFCHAMBER'
    longName = 'homogenized fission chamber'
    atomDensity = 7.60666E-02
    fuelRadius = 0.54 #cm
    temperature = 293.6 #K
    thermalOpt = 'h2o'
    uAtomFractionsDict = {235: 1.0}
    elemAtomFracDict = {'U': 5.2117E-08, 'O': 0.174590749, 'H': 0.349181499, 'B': 0.000144884, 'Al': 0.476082816}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

def get_c5g7_control_rod_material():
    shortName = 'cCR'
    longName = 'homogenized control rod'
    atomDensity = 5.84389E-02
    fuelRadius = 0.54 #cm
    temperature = 293.6 #K
    thermalOpt = 'free'
    elemAtomFracDict = {'Al': 0.61969028, 'Ag': 0.30424777, 'In': 0.05704646, 'Cd': 0.01901549}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity)
    return material

###############################################################################
def get_triga_fuel_material():
    shortName = 'tFUEL'
    longName = 'U-ZrH fuel'
    atomDensity = 8.71115E-2
    fuelRadius = 1.7920 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'zrh'
    uAtomFractionsDict = {234: 0.000132187, 235: 0.017378221, 236: 0.000194346, 237: 1.0E-24, 238: 0.069406717}
    #npAtomFractionsDict = {237: 1.0E-24, 238: 1.0E-24, 239: 1.0E-24}
    #puAtomFractionsDict = {238: 1.0E-24, 239: 1.0E-24, 240: 1.0E-24, 241: 1.0E-24, 242: 1.0E-24}
    #amAtomFractionsDict = {241: 1.0E-24, 242: 1.0E-24, 243: 1.0E-24, 642: 1.0E-24}
    #cmAtomFractionsDict = {242: 1.0E-24, 243: 1.0E-24}
    #xeAtomFractionsDict = {135: 1.0E-24}
    #smAtomFractionsDict = {149: 1.0E-24}
    hAtomFractionsDict = {1: 0.999885, 2: 0.000115}
    #elemAtomFracDict = {'U': 0.062260112, 'Np': 1E-24, 'Pu': 1E-24, 'Am': 1E-24, 'Cm': 1E-24, 'Xe': 1E-24, 'Sm': 1E-24, 'Zr': 0.370555099, 'Hf': 2.22332E-05, 'Er': 0.002650567, 'C': 0.000205141, 'H': 0.564306848}
    elemAtomFracDict = {'U': 0.062260112, 'Zr': 0.370555099, 'Hf': 2.22332E-05, 'Er': 0.002650567, 'C': 0.000205141, 'H': 0.564306848}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #override_abundances(ZAList, abundanceDict, npAtomFractionsDict, 'Np')
    #override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    #override_abundances(ZAList, abundanceDict, amAtomFractionsDict, 'Am')
    #override_abundances(ZAList, abundanceDict, cmAtomFractionsDict, 'Cm')
    #override_abundances(ZAList, abundanceDict, xeAtomFractionsDict, 'Xe')
    #override_abundances(ZAList, abundanceDict, smAtomFractionsDict, 'Sm')
    override_abundances(ZAList, abundanceDict, hAtomFractionsDict, 'H')
    override_abundances_as_elemental(ZAList, abundanceDict, 'C')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_clad_material():
    shortName = 'tCLAD'
    longName = 'SS304 clad'
    atomDensity = 8.58765E-02
    fuelRadius = 1.7920 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'free'
    elemAtomFracDict = {'Cr': 0.202441879, 'Fe': 0.687731801, 'Ni': 0.089657823, 'Mn': 0.020168498}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_zirconium_material():
    shortName = 'tZIRC'
    longName = 'Zr fuel center'
    atomDensity = 4.29600E-2
    fuelRadius = 0.3175 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'free'
    elemAtomFracDict = {'Zr': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_graphite_material():
    shortName = 'tGRAPHITE'
    longName = 'graphite'
    atomDensity = 8.52100E-02
    fuelRadius = 2.19255 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'graphite'
    elemAtomFracDict = {'C': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances_as_elemental(ZAList, abundanceDict, 'C')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_borated_graphite_material():
    shortName = 'tBORATEDGRAPHITE'
    longName = 'borated graphite'
    atomDensity = 1.27794E-01
    fuelRadius = 1.7411 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'graphite'
    elemAtomFracDict = {'C': 0.788925928, 'B': 0.211074072}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances_as_elemental(ZAList, abundanceDict, 'C')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_b4c_material():
    shortName = 'tB4C'
    longName = 'B4C powder'
    atomDensity = 1.35143E-01
    fuelRadius = 1.5164 #cm
    temperature = 296. #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'free'
    elemAtomFracDict = {'C': 0.201616066, 'B': 0.798383934}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances_as_elemental(ZAList, abundanceDict, 'C')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_irradiation_tube_material():
    shortName = 'tIRRADIATIONTUBE'
    longName = 'homogenized air-filled, Al-lined tube in water'
    atomDensity = 4.90539E-02
    fuelRadius = 2.19255 #cm
    temperature = 293.6 #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'h2o'
    elemAtomFracDict = {'O': 0.29203693, 'H': 0.58387295, 'N': 0.00038432, 'Al': 0.12370580}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_moderator_material():
    shortName = 'tMOD'
    longName = 'light water'
    atomDensity = 1.00037E-1
    fuelRadius = 1.7920 #cm
    temperature = 293.6 #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'h2o'
    elemAtomFracDict = {'O': 1.0, 'H': 2.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_air_tube_material():
    shortName = 'tAIRTUBE'
    longName = 'homogenized air-filled tube in water'
    atomDensity = 2.73966E-02
    fuelRadius = 3.8 #cm
    temperature = 293.6 #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'h2o'
    elemAtomFracDict = {'O': 0.00017253, 'H': 0.00017253, 'N': 0.00017253}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material

def get_triga_grid_plate_material():
    shortName = 'tGRIDPLATE'
    longName = 'structural Al'
    #massDensity = 2.7 # g/cc
    atomDensity = 0.059195 #a/b-cm
    fuelRadius = 10. #cm (complete guess)
    temperature = 293.6 #K
    temperatureIndex = 0 # X in .9Xc
    thermalOpt = 'al'
    #thermalOpt = 'free'
    elemAtomFracDict = {'Al': 0.058693, 'Fe': 0.000502}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, atomDensity=atomDensity,
        temperatureIndex=temperatureIndex)
    return material


###############################################################################
def get_bruss_enriched_rod_fuel_material():
    shortName = 'debFUEL'
    longName = 'uranium dioxide (rod)'
    massDensity = 10.29769 #g/cc
    fuelRadius = 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    temperature = 500 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 0.95, 235: 0.05}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

###############################################################################
def get_kord_fuel_material():
    shortName = 'kFUEL'
    longName = 'depleted uranium dioxide'
    massDensity = 10.29769 #g/cc
    MCNPfactor = 2.01692089627977
    fuelRadius = MCNPfactor * 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    #temperature = 296. #K
    #thermalOpt = 'uo2'
    temperature = 293.6 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 1.0}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_rod_fuel_material():
    shortName = 'kRFUEL'
    longName = 'depleted uranium dioxide (rod)'
    massDensity = 10.29769 #g/cc
    # In a 2D pin cell, the fuel is a cylinder/rod, not a slab, and a different chord length is used
    fuelRadius = 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    #temperature = 296. #K
    #thermalOpt = 'uo2'
    temperature = 293.6 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 1.0}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_enriched_rod_fuel_material():
    shortName = 'kREFUEL'
    longName = 'uranium dioxide (rod)'
    massDensity = 10.29769 #g/cc
    # In a 2D pin cell, the fuel is a cylinder/rod, not a slab, and a different chord length is used
    fuelRadius = 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    #temperature = 296. #K
    #thermalOpt = 'uo2'
    temperature = 293.6 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 0.96, 235: 0.04}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_mox_rod_fuel_material():
    shortName = 'kRMFUEL'
    longName = 'MOX (rod)'
    massDensity = 10.29769 #g/cc
    # In a 2D pin cell, the fuel is a cylinder/rod, not a slab, and a different chord length is used
    fuelRadius = 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    #temperature = 296. #K
    #thermalOpt = 'uo2'
    temperature = 293.6 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 1.0}
    puAtomFractionsDict = {239: 1.0}
    elemAtomFracDict = {'O': 2, 'U': 0.95, 'Pu': 0.05}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    override_abundances(ZAList, abundanceDict, puAtomFractionsDict, 'Pu')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_enriched_fuel_material():
    shortName = 'kEFUEL'
    longName = 'enriched uranium dioxide'
    massDensity = 10.29769 #g/cc
    MCNPfactor = 2.01692089627977
    fuelRadius = MCNPfactor * 0.39218 #cm
    # Want a consistent thermal behavior with MCNP
    #temperature = 296. #K
    #thermalOpt = 'uo2'
    temperature = 293.6 #K
    thermalOpt = 'free'
    uAtomFractionsDict = {238: 0.96, 235: 0.04}
    elemAtomFracDict = {'O': 2, 'U': 1}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    override_abundances(ZAList, abundanceDict, uAtomFractionsDict, 'U')
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_zirconium_material():
    shortName = 'kZR'
    longName = 'zirconium cladding'
    # This density is a homogenization of clad and void from an OpenMC pincell
    massDensity = 5.8105 #g/cc
    fuelRadius = 0.45720 #cm
    temperature = 293.6 #K
    thermalOpt = 'free'
    elemAtomFracDict = {'Zr': 1.0}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_clad_material():
    shortName = 'kCLAD'
    longName = 'zircaloy cladding'
    # This density is a homogenization of clad and void from an OpenMC pincell
    massDensity = 5.8105 #g/cc
    fuelRadius = 0.45720 #cm
    temperature = 293.6 #K
    thermalOpt = 'free'
    elemAtomFracDict = {'O': 0.0070946, 'Cr': 0.0017464, 'Fe': 0.0034148, 'Zr': 0.9766522, 'Sn': 0.011092}
    #
    chordLength = calc_chord_length(fuelRadius)
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_kord_moderator_material():
    material = get_h2o_material()
    material.update_temperature(293.6) #K
    material.update_mass_density(0.740582)
    material.update_names('kMOD', 'cold light water')
    return material

###############################################################################
def get_multi_temperature_h2o_material_base():
    shortName = 'mtH2O'
    longName = 'light water'
    massDensity = 1.0 #g/cc
    fuelRadius = 1.7920 #cm
    temperature = 293.6 #K
    # HACK TO ALIGN THERMAL GRIDS!!!
    #thermalOpt = 'h2o'
    thermalOpt = 'free'
    elemAtomFracDict = {'H': 2, 'O': 1}
    #
    symDict, ZList, ZAList = get_all_isotopes(elemAtomFracDict)
    abundanceDict = lookup_natl_abundances(ZAList)
    chordLength = calc_chord_length(fuelRadius)
    #
    material = Material(
        shortName=shortName, longName=longName,
        temperature=temperature, thermalOpt=thermalOpt,
        symDict=symDict, ZList=ZList, ZAList=ZAList,
        abundanceDict=abundanceDict, chordLength=chordLength,
        elemAtomFracDict=elemAtomFracDict, massDensity=massDensity)
    return material

def get_multi_temperature_h2o_material_Tgrid():
    return [293.6, 350, 400, 450, 500, 550, 600, 650, 800] #K

def get_multi_temperature_h2o_material_Ti(iT):
    P = 0.101325 #MPa
    #P = 15.5 #MPa
    material = get_multi_temperature_h2o_material_base()
    Tgrid = get_multi_temperature_h2o_material_Tgrid()
    T = Tgrid[iT] 
    rho = steam(P=P, T=T).rho / 1000
    shortName = 'mtH2O_{}'.format(iT)
    longName = 'light water ({} K)'.format(T)
    material.update_temperature(T)
    material.update_temperature_index(iT) # X in .9Xc
    material.update_mass_density(rho)
    material.update_names(shortName, longName)
    return material
    
def get_multi_temperature_h2o_material_T0():
    return get_multi_temperature_h2o_material_Ti(0)

def get_multi_temperature_h2o_material_T1():
    return get_multi_temperature_h2o_material_Ti(1)

def get_multi_temperature_h2o_material_T2():
    return get_multi_temperature_h2o_material_Ti(2)

def get_multi_temperature_h2o_material_T3():
    return get_multi_temperature_h2o_material_Ti(3)

def get_multi_temperature_h2o_material_T4():
    return get_multi_temperature_h2o_material_Ti(4)

def get_multi_temperature_h2o_material_T5():
    return get_multi_temperature_h2o_material_Ti(5)

def get_multi_temperature_h2o_material_T6():
    return get_multi_temperature_h2o_material_Ti(6)

def get_multi_temperature_h2o_material_T7():
    return get_multi_temperature_h2o_material_Ti(7)

def get_multi_temperature_h2o_material_T8():
    return get_multi_temperature_h2o_material_Ti(8)

###############################################################################
def get_multi_temperature_triga_fuel_material_base():
    return get_triga_fuel_material()

def get_multi_temperature_triga_fuel_material_Tgrid():
    return [296, 400, 500, 600, 700, 800, 1000, 1200] #K

def get_multi_temperature_triga_fuel_material_Ti(iT):
    material = get_multi_temperature_triga_fuel_material_base()
    Tgrid = get_multi_temperature_triga_fuel_material_Tgrid()
    T = Tgrid[iT] 
    shortName = 'mtTFUEL_{}'.format(iT)
    longName = '{} ({} K)'.format(material.longName, T)
    material.update_temperature(T)
    material.update_temperature_index(iT) # X in .9Xc
    # MASS DENSITY SHOULD BE UPDATED TO ACCOUNT FOR THERMAL EXPANSION
    #material.update_mass_density(rho)
    material.update_names(shortName, longName)
    return material

def get_multi_temperature_triga_fuel_material_T0():
    return get_multi_temperature_triga_fuel_material_Ti(0)

def get_multi_temperature_triga_fuel_material_T1():
    return get_multi_temperature_triga_fuel_material_Ti(1)

def get_multi_temperature_triga_fuel_material_T2():
    return get_multi_temperature_triga_fuel_material_Ti(2)

def get_multi_temperature_triga_fuel_material_T3():
    return get_multi_temperature_triga_fuel_material_Ti(3)

def get_multi_temperature_triga_fuel_material_T4():
    return get_multi_temperature_triga_fuel_material_Ti(4)

def get_multi_temperature_triga_fuel_material_T5():
    return get_multi_temperature_triga_fuel_material_Ti(5)

def get_multi_temperature_triga_fuel_material_T6():
    return get_multi_temperature_triga_fuel_material_Ti(6)

def get_multi_temperature_triga_fuel_material_T7():
    return get_multi_temperature_triga_fuel_material_Ti(7)

def get_multi_temperature_triga_fuel_material_T8():
    return get_multi_temperature_triga_fuel_material_Ti(8)
    
###############################################################################
def get_all_isotopes(elemDict):
    symList = elemDict.keys()
    symDict = {}
    ZList = []
    ZAList = []
    for sym in symList:
        Z = nd.sym2z[sym.capitalize()]
        ZList.append(Z)
        symDict[Z] = sym
        if sym == 'Th':
            ZAList += get_Th_isotopes()
        elif sym == 'Pa':
            ZAList += get_Pa_isotopes()
        elif sym == 'U':
            ZAList += get_U_isotopes()
        elif sym == 'Np':
            ZAList += get_Np_isotopes()
        elif sym == 'Pu':
            ZAList += get_Pu_isotopes()
        elif sym == 'Am':
            ZAList += get_Am_isotopes()
        elif sym == 'Cm':
            ZAList += get_Cm_isotopes()
        elif sym == 'Bk':
            ZAList += get_Bk_isotopes()
        elif sym == 'Cf':
            ZAList += get_Cf_isotopes()
        elif sym == 'Es':
            ZAList += get_Es_isotopes()
        elif sym == 'Fm':
            ZAList += get_Fm_isotopes()
        else:
            ZAList += get_isotope_ZAs(Z)
    return symDict, ZList, ZAList

def get_isotope_ZAs(Z, cutoff=0.005):
    '''Get all isotopes with natural abundance at least cutoff for element Z.
    nd.isotopes does not return any metastable A's.'''
    AList = [A for A in nd.isotopes[Z] if nd.nuc(Z, A)['abundance'].nominal_value > cutoff]
    return [(Z, A) for A in AList]

def get_Th_isotopes():
    ZAList = []
    ZAList.append((90, 227))
    ZAList.append((90, 228))
    ZAList.append((90, 229))
    ZAList.append((90, 230))
    ZAList.append((90, 231))
    ZAList.append((90, 232))
    ZAList.append((90, 233))
    ZAList.append((90, 234))
    return ZAList

def get_Pa_isotopes():
    ZAList = []
    ZAList.append((91, 229))
    ZAList.append((91, 230))
    ZAList.append((91, 231))
    ZAList.append((91, 232))
    ZAList.append((91, 233))
    return ZAList

def get_U_isotopes():
    ZAList = []
    ZAList.append((92, 230))
    ZAList.append((92, 231))
    ZAList.append((92, 232))
    ZAList.append((92, 233))
    ZAList.append((92, 234))
    ZAList.append((92, 235))
    ZAList.append((92, 236))
    ZAList.append((92, 237))
    ZAList.append((92, 238))
    ZAList.append((92, 239))
    ZAList.append((92, 240))
    ZAList.append((92, 241))
    return ZAList

def get_Np_isotopes():
    ZAList = []
    ZAList.append((93, 234))
    ZAList.append((93, 235))
    ZAList.append((93, 236))
    ZAList.append((93, 237))
    ZAList.append((93, 238))
    ZAList.append((93, 239))
    return ZAList

def get_Pu_isotopes():
    ZAList = []
    ZAList.append((94, 236))
    ZAList.append((94, 237))
    ZAList.append((94, 238))
    ZAList.append((94, 239))
    ZAList.append((94, 240))
    ZAList.append((94, 241))
    ZAList.append((94, 242))
    ZAList.append((94, 243))
    ZAList.append((94, 244))
    ZAList.append((94, 246))
    return ZAList

def get_Am_isotopes():
    ZAList = []
    ZAList.append((95, 240))
    ZAList.append((95, 241))
    ZAList.append((95, 242))
    ZAList.append((95, 642)) # metastable
    ZAList.append((95, 243))
    ZAList.append((95, 244))
    ZAList.append((95, 644)) # metastable
    return ZAList

def get_Cm_isotopes():
    ZAList = []
    ZAList.append((96, 240))
    ZAList.append((96, 241))
    ZAList.append((96, 242))
    ZAList.append((96, 243))
    ZAList.append((96, 244))
    ZAList.append((96, 245))
    ZAList.append((96, 246))
    ZAList.append((96, 247))
    ZAList.append((96, 248))
    ZAList.append((96, 249))
    ZAList.append((96, 250))
    return ZAList

def get_Bk_isotopes():
    ZAList = []
    ZAList.append((97, 245))
    ZAList.append((97, 246))
    ZAList.append((97, 247))
    ZAList.append((97, 248))
    ZAList.append((97, 249))
    ZAList.append((97, 250))
    return ZAList

def get_Cf_isotopes():
    ZAList = []
    ZAList.append((98, 246))
    ZAList.append((98, 248))
    ZAList.append((98, 249))
    ZAList.append((98, 250))
    ZAList.append((98, 251))
    ZAList.append((98, 252))
    ZAList.append((98, 253))
    ZAList.append((98, 254))
    return ZAList

def get_Es_isotopes():
    ZAList = []
    ZAList.append((99, 251))
    ZAList.append((99, 252))
    ZAList.append((99, 253))
    ZAList.append((99, 254))
    ZAList.append((99, 654)) # metastable
    ZAList.append((99, 255))
    return ZAList

def get_Fm_isotopes():
    ZAList = []
    ZAList.append((100, 255))
    return ZAList

###############################################################################
def lookup_natl_abundances(ZAList):
    '''Following MCNP, metastable states have A increased by 400: e.g., Am-242m has an A of 642'''
    Zset = set([Z for (Z,A) in ZAList])
    abundanceDict = {}
    for Zelem in Zset:
        AList = [A for (Z,A) in ZAList if Z == Zelem]
        # Use ground-state abundances as default metastable abundances (may override later)
        abundanceList = [nd.nuc(Zelem,A%400)['abundance'].nominal_value for A in AList]
        maxIndex = np.argmax(abundanceList)
        norm = np.sum(abundanceList)
        if norm:
            abundanceList[maxIndex] *= (1.0 / norm)
        for A, abundance in zip(AList, abundanceList):
            abundanceDict[(Zelem,A)] = abundance
    return abundanceDict

def override_abundances_as_elemental(ZAList, abundanceDict, sym):
    '''Override isotopes and abundances for element sym to be elemental values (used for carbon)'''
    Zthis = nd.sym2z[sym]
    AList = [A for (Z,A) in ZAList if Z == Zthis]
    for A in AList:
        i = ZAList.index((Zthis, A))
        del(ZAList[i])
        del(abundanceDict[(Zthis, A)])
    ZAList.append((Zthis,0))
    abundanceDict[(Zthis, 0)] = 1.0

def override_abundances(ZAList, abundanceDict, atomFractionsDict, sym):
    '''Override isotopes and abundances for element given by sym'''
    Zthis = nd.sym2z[sym]
    AList = [A for (Z,A) in ZAList if Z == Zthis]
    #
    # Zero-out existing abundances
    for A in AList:
        abundanceDict[(Zthis, A)] = 0.0
    #
    # Replace with abundances in atomFractionsDict
    if not atomFractionsDict:
        raise ValueError('atomsFractionsDict for {0} must be non-empty.'.format(sym))
    # Python guarantees values() and keys() are congruent if no changes to the dict are made between calls
    atomFractionsList = np.array(atomFractionsDict.values())
    atomFractionsKeys = atomFractionsDict.keys()
    norm = np.sum(atomFractionsList)
    if norm == 0.:
        # If all atom fractions are zero, normalize by making each equal to 1/(# isotopes)
        numIsotopes = len(atomFractionsList)
        for A in AList:
            abundanceDict[(Zthis, A)] = 1. / numIsotopes
    elif norm > 1.1 or norm < 0.9:
        # If atom fractions are very unnormalized, multiplicatively normalize
        for A in atomFractionsKeys:
            abundanceDict[(Zthis, A)] = atomFractionsDict[A] / norm
    else:
        # Else, normalize by changing the abundance of the most abundant isotope
        majorA = atomFractionsKeys[np.argmax(atomFractionsList)]
        partialSum = np.sum([atomFractionsDict[A] for A in atomFractionsKeys if A != majorA])
        atomFractionsDict[majorA] = 1.0 - partialSum
        for A in atomFractionsKeys:
            abundanceDict[(Zthis, A)] = atomFractionsDict[A]
    #
    # Add in new isotopes to ZAList
    newAs = [A for A in atomFractionsKeys if (Zthis, A) not in ZAList]
    for A in newAs:
        ZAList.append((Zthis, A))
    #
    # Remove entries with zero abundances
    for A in AList:
        if abundanceDict[(Zthis, A)] == 0.0:
            i = ZAList.index((Zthis, A))
            del(ZAList[i])
            del(abundanceDict[(Zthis, A)])

###############################################################################
def print_materials(materials, verbosity=False):
    if verbosity:
        for i, material in enumerate(materials):
            print '------- Material {0} -------'.format(i)
            material.print_contents(verbosity)

###############################################################################
class Material():
    def __init__(self, shortName=None, longName=None, temperature=None, atomDensity=None, massDensity=None, chordLength=None, abundanceDict=None, elemAtomFracDict=None, elemMassFracDict=None, ZAList=None, ZList=None, symDict=None, backgroundXSDict=None, admixedModeratorDict=None, admixedModeratorList=None, thermalOpt=0, elemWeightDict=None, matlWeight=None, thermalXSDict=None, temperatureIndex=0):
        '''The following units are used:
        temperature in Kelvin is the temperature of the material
        temperatureIndex is the least-significant digit of an MCNP-like matl. name and refers to temperature
        atomDensity in atoms/barn-cm is the atom density of the material
        massDensity in g/cm^3 is the mass density of the material
        chordLength in cm is used to calculate the escape cross section; set to zero to use 0 escape XS
        abundanceDict is the atom fraction of each isotope in its element
        elemAtomFracDict is the atom fraction of each element in the material
        elemMassFracDict is the mass/weight fraction of each element in the material
        elemWeightDict is the elemental mass/weight in g/mole calculated using the correct abundances
        matlWeight in g/mole is the effective weight of the material
        ZAList is a list of all (Z,A) = (atomic number, atomic mass) pairs in the problem
        ZList is a list of all elements in the problem stored as atomic number
        symDict is the chemical symbol for each element in the material
        backgroundXSDict is an approximate background XS seen by each nuclide in the material
        admixedModeratorDict is currently not used
        admixedModeratorList is currently not used
        thermalOpt specifies which thermal treatment to use
        thermalXSDict is a list of thermal XS used for each element (list of shortName's in endfMTList)
        shortName is a succinct name without whitespaces
        longName is a longer name that may contain whitespace
        '''
        # locals() includes local names, including 'self'. Update Material to include all keywords above.
        varss = locals()
        self.__dict__.update(varss)
        del self.__dict__["self"]
        self.setup()

    def setup(self):
        '''Initialize attributes. Prefer atom fractions and mass densities over mass fractions and atom densities if both are given.'''
        self.make_lists_sets()
        self.strip_short_name()
        #
        self.calc_and_store_elem_weight()
        #
        if self.elemAtomFracDict:
            self.change_atom_frac_dict_keys_to_Z()
            self.calc_and_store_elem_mass_frac()
        elif self.elemMassFracDict:
            self.change_mass_frac_dict_keys_to_Z()
            self.calc_and_store_elem_atom_frac()
        else:
            raise ValueError('Must set either elemAtomFracDict or elemMassFracDict')
        #
        self.calc_and_store_matl_weight()
        #
        if self.massDensity:
            self.calc_and_store_atom_density()
        elif self.atomDensity:
            self.calc_and_store_mass_density()
        else:
            raise ValueError('Must set either atomDensity or massDensity')
        #
        self.set_thermal_xs_dict()
        self.check_internal_keys_consistency()

    def copy(self):
        return copy.deepcopy(self)

    def update_names(self, shortName, longName):
        self.shortName = shortName
        self.longName = longName

    def update_thermal_option(self, thermalOpt):
        self.thermalOpt = thermalOpt
        self.set_thermal_xs_dict()

    def update_temperature(self, temperature):
        self.temperature = temperature

    def update_temperature_index(self, temperatureIndex):
        self.temperatureIndex = temperatureIndex

    def update_chord_length(self, chordLength):
        self.chordLength = chordLength

    def update_mass_density(self, massDensity):
        self.massDensity = massDensity
        self.calc_and_store_atom_density()

    def update_atom_density(self, atomDensity):
        self.atomDensity = atomDensity
        self.calc_and_store_mass_density()

    def update_atom_fractions_and_atom_density(self, elemAtomFracDict, atomDensity):
        self.elemAtomFracDict = elemAtomFracDict
        self.calc_and_store_elem_mass_frac()
        self.calc_and_store_matl_weight()
        self.update_atom_density(atomDensity)

    def update_atom_fractions_and_mass_density(self, elemAtomFracDict, massDensity):
        self.elemAtomFracDict = elemAtomFracDict
        self.calc_and_store_elem_mass_frac()
        self.calc_and_store_matl_weight()
        self.update_mass_density(massDensity)

    def update_mass_fractions_and_atom_density(self, elemMassFracDict, atomDensity):
        self.elemMassFracDict = elemMassFracDict
        self.calc_and_store_elem_atom_frac()
        self.calc_and_store_matl_weight()
        self.update_atom_density(atomDensity)

    def update_mass_fractions_and_mass_density(self, elemMassFracDict, massDensity):
        self.elemMassFracDict = elemMassFracDict
        self.calc_and_store_elem_atom_frac()
        self.calc_and_store_matl_weight()
        self.update_mass_density(massDensity)

    ################################################################
    def print_contents(self, verbosity=0):
        if verbosity:
            strs = ['shortName', 'longName']
            for strr in strs:
                print strr, getattr(self, strr)
            #
            strr = 'thermalOpt'
            val = getattr(self, strr)
            print strr, val
            #
            strs = [('chordLength', 'cm'), ('temperature', 'K'), ('atomDensity', '1/b-cm'), ('massDensity', 'g/cc') ]
            for (strr, unit) in strs:
                print strr, getattr(self, strr), '({0})'.format(unit)
            #
        if verbosity > 1:
            strr, unit =  ('matlWeight', 'g/mole')
            print strr, getattr(self, strr), '({0})'.format(unit)
            #
        if verbosity:
            sortedSymList = [self.symDict[Z] for Z in sorted(self.ZList)]
            print 'symList', sortedSymList
            #
        if verbosity > 1:
            strs = ['ZList', 'ZAList']
            for strr in strs:
                print strr, sorted(getattr(self, strr))
            #
            strs = [('abundanceDict', 'atom fraction'), ('elemWeightDict', 'g/mole')]
            for (strr, unit) in strs:
                print strr, getattr(self, strr), '({0})'.format(unit)
            #
            strs = ['elemAtomFracDict', 'elemMassFracDict']
            for strr in strs:
                print strr, getattr(self, strr)
            #
            strs = ['thermalXSDict', 'backgroundXSDict']
            for strr in strs:
                print strr, getattr(self, strr)
            #
        if verbosity:
            a = getattr(self, 'abundanceDict')
            e = getattr(self, 'elemAtomFracDict')
            aD = getattr(self, 'atomDensity')
            T_Kelvin = getattr(self, 'temperature')
            T_MeV = 8.6173324E-11 * T_Kelvin
            iT = getattr(self, 'temperatureIndex')
            norm = sum([a[(Z,A)]*e[Z] for (Z,A) in a])
            print 'MCNP-style material input:'
            print '     atom_density={:.5f}'.format(aD)
            print '     tmp={:.3e}'.format(T_MeV)
            for Z,A in sorted(a.keys()):
                # NB: Some of the abundanceDict's may be unnormalized
                # because they neglect some of the isotopes
                atomFrac = a[(Z,A)] * e[Z]
                print '     {:02}{:03}.{}c    {:.8e}'.format(Z, A, 90+iT, atomFrac)
            t = getattr(self, 'thermalXSDict')
            therm2mcnpDict = util.get_thermal_short2mcnp_dict()
            for vList in t.values():
                for v in vList:
                    if v in therm2mcnpDict:
                        print '     {}.{}t'.format(therm2mcnpDict[v], 20+iT)

    ################################################################
    def make_lists_sets(self):
        self.ZAList = set(self.ZAList)
        self.ZList = set(self.ZList)

    def check_internal_keys_consistency(self):
        if self.shortName.count(' '):
            raise ValueError('shortName should not contain spaces')
        if set(self.abundanceDict.keys()) != self.ZAList:
            raise ValueError('ZAList should be the keys for abundanceDict')
        if set(self.elemAtomFracDict.keys()) != self.ZList:
            raise ValueError('ZList should be the keys for elemAtomFracDict')
        if set(self.elemMassFracDict.keys()) != self.ZList:
            raise ValueError('ZList should be the keys for elemMassFracDict')
        if set(self.elemWeightDict.keys()) != self.ZList:
            raise ValueError('ZList should be the keys for elemWeightDict')
        if set(self.thermalXSDict.keys()) != self.ZList:
            raise ValueError('ZList should be the keys for thermalXSDict')
        if set(self.symDict.keys()) != self.ZList:
            raise ValueError('ZList should be the keys for consistent with symDict')
        if set([Z for (Z, A) in self.ZAList]) != self.ZList:
            raise ValueError('ZList should be consistent with ZAList')

    def strip_short_name(self):
        '''Remove whitespace from shortName'''
        namePieces = self.shortName.split()
        namePieces = [piece.strip() for piece in namePieces]
        self.shortName = ''.join(namePieces)

    def change_atom_frac_dict_keys_to_Z(self):
        '''Change keys from using sym to using Z'''
        atomFracDictFromSym = self.elemAtomFracDict
        atomFracDictFromZ = {}
        for sym in atomFracDictFromSym.keys():
            atomFracDictFromZ[nd.sym2z[sym]] = atomFracDictFromSym[sym]
        self.elemAtomFracDict = atomFracDictFromZ

    def change_mass_frac_dict_keys_to_Z(self):
        '''Change keys from using sym to using Z'''
        massFracDictFromSym = self.elemMassFracDict
        massFracDictFromZ = {}
        for sym in massFracDictFromSym.keys():
            massFracDictFromZ[nd.sym2z[sym]] = massFracDictFromSym[sym]
        self.elemMassFracDict = massFracDictFromZ

    def calc_and_store_elem_weight(self):
        '''Calculate the molar mass of each element from its isotopic abundance.'''
        self.elemWeightDict = {}
        for Zelem in self.ZList:
            AList = [A for (Z,A) in self.ZAList if Z == Zelem]
            elemAtomFracList = [self.abundanceDict[(Zelem,A)] for A in AList]
            if A:
                # For metastable states (A increased by 400), use groundstate weights
                weightList = [nd.nuc(Zelem, A%400)['weight'].nominal_value for A in AList]
                self.elemWeightDict[Zelem] = np.sum(np.multiply(elemAtomFracList, weightList))
            else:
                # A of 0 indicates natural composition (used for carbon)
                self.elemWeightDict[Zelem] = nd.weight(Zelem)

    def normalize_elem_atom_frac(self):
        '''Make sum(elemAtomFracDcit) == 1'''
        norm = np.sum((float(value) for value in self.elemAtomFracDict.values()))
        for key in self.elemAtomFracDict:
            self.elemAtomFracDict[key] /= norm

    def normalize_elem_mass_frac(self):
        '''Make sum(elemMassFracDcit) == 1'''
        norm = np.sum((float(value) for value in self.elemMassFracDict.values()))
        for key in self.elemMassFracDict:
            self.elemMassFracDict[key] /= norm

    def calc_and_store_elem_atom_frac(self):
        '''Convert from weight fractions to atom fractions.'''
        self.normalize_elem_mass_frac()
        self.elemAtomFracDict = {}
        norm = np.sum((self.elemMassFracDict[Z] / self.elemWeightDict[Z] for Z in self.ZList))
        for Z in self.ZList:
            self.elemAtomFracDict[Z] = self.elemMassFracDict[Z] / (self.elemWeightDict[Z] * norm)
        self.normalize_elem_atom_frac()

    def calc_and_store_elem_mass_frac(self):
        '''Convert from atom fractions to weight fractions.'''
        self.normalize_elem_atom_frac()
        self.elemMassFracDict = {}
        norm = np.sum((self.elemAtomFracDict[Z] * self.elemWeightDict[Z] for Z in self.ZList))
        for Z in self.ZList:
            self.elemMassFracDict[Z] = self.elemAtomFracDict[Z] * self.elemWeightDict[Z] / norm
        self.normalize_elem_mass_frac()

    def calc_and_store_matl_weight(self):
        '''Calculate the molar mass of the material. Requires elemWeightDict and elemAtomFracDict be set.'''
        self.matlWeight = np.sum((self.elemWeightDict[Z] * self.elemAtomFracDict[Z] for Z in self.ZList))

    def calc_and_store_atom_density(self):
        '''Convert from g/cm^3 to 1/b-cm. Requires massDensity and matlWeight be set.'''
        self.atomDensity = self.massDensity * util.avogadros_number() / self.matlWeight

    def calc_and_store_mass_density(self):
        '''Convert from 1/b-cm to g/cm^3. Requires atomDensity and matlWeight be set.'''
        self.massDensity = self.matlWeight * self.atomDensity / util.avogadros_number()

    ################################################################
    def set_thermal_xs_dict(self):
        # WARNING: thermalXSDict should be a ZA list, not a Z list, because the thermal
        # option is often nuclide-specific. See util.get_thermal_mcnp2zaid_dict for
        # which nuclides should use which thermal treatments
        # Current cases where this script does not match MCNP:
        # H-2 and H-3 in 'h2o'
        # Non Fe-56 isotopes in 'fe'
        # U-235 (and anything other than U-238) in 'uo2'
        self.thermalXSDict = {}
        self.thermalOpt = self.thermalOpt.lower().strip()
        if self.thermalOpt == 'none':
            for Z in self.ZList:
                self.thermalXSDict[Z] = []
        else:
            for Z in self.ZList:
                self.thermalXSDict[Z] = ['free']
        if self.thermalOpt == 'h2o':
            for Z in self.ZList:
                if Z == 1:
                    self.thermalXSDict[Z] = ['hh2o']
        elif self.thermalOpt == 'al':
            for Z in self.ZList:
                if Z == 13:
                    self.thermalXSDict[Z] = ['alinel', 'alelas']
        elif self.thermalOpt == 'fe':
            for Z in self.ZList:
                if Z == 26:
                    # This XS only applies to Fe-56, but here we are using it for all Fe isotopes!
                    # This XS may only apply to Fe in a BCC configuration, notably not stainless steel
                    self.thermalXSDict[Z] = ['feinel', 'feelas']
        elif self.thermalOpt == 'uo2':
            for Z in self.ZList:
                if Z == 8:
                    self.thermalXSDict[Z] = ['ouo2inel', 'ouo2elas']
                if Z == 92:
                    # This XS only applies to U-238, but here we are using it for all U isotopes!
                    # (The UUO2 thermal treatment is very wrong for Pu due to low-lying resonances)
                    self.thermalXSDict[Z] = ['uuo2inel', 'uuo2elas']
        elif self.thermalOpt == 'graphite':
            for Z in self.ZList:
                if Z == 6:
                    self.thermalXSDict[Z] = ['graphinel', 'graphelas']
        elif self.thermalOpt == 'zrh':
            for Z in self.ZList:
                if Z == 1:
                    self.thermalXSDict[Z] = ['hzrhinel', 'hzrhelas']
                if Z == 40:
                    self.thermalXSDict[Z] = ['zrzrhinel', 'zrzrhelas']

    def check_background_xs_keys_consistency(self):
        if set(self.backgroundXSDict.keys()) != self.ZAList:
            raise ValueError('ZAList should be the keys for backgroundXSDict')
