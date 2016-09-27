#! /usr/bin/env python

'''
Andrew Till
Winter 2015

Python-based PDT cross section reader / writer
'''

#STDLIB
from datetime import datetime

#TPL
import numpy as np

#########################################################################################
def read_PDT_xs_generally(filePath):
    '''Read a PDT cross section and store in PDT_XS class. Must have 1 temperature and density'''
    with open(filePath, 'r') as fid:
        # >>>>>> HEADER <<<<<<
        # Skip initial lines
        for i in range(2):
            fid.readline()
        # Read xs type (MG/MB)
        t = fid.readline().split()
        xsType = t[4]
        # Read number of groups. Assume 1 temperature and 1 density
        t = fid.readline().split()
        numGroups = int(t[5])
        # Read number of 1D cross sections and transfer matrices
        fid.readline()
        t = fid.readline().split()
        num1D = int(t[0])
        numXfer = int(t[4])
        # Read number of Legendre moments
        t = fid.readline().split()
        numMoments = int(t[2]) + 1
        # Read xs units (micro/macro)
        fid.readline()
        xsUnits = fid.readline().strip()
        # Read temperature
        for i in range(2):
            fid.readline()
        t = fid.readline().split()
        temperature = float(t[0])
        for i in range(5):
            fid.readline()
        # Initialize number of delayed neutron groups as 0
        numDNGs = 0

        # >>>>>> GROUP BOUNDARIES <<<<<<
        # From now on, must deal with fixed-length rows (5 entries per row)
        loc = 0
        entriesPerLine = 5
        groupBdrs = np.zeros(numGroups+1)
        read_fixed_line(groupBdrs, numGroups+1, entriesPerLine, fid)
        groupWidths = - np.diff(groupBdrs)

        # >>>>>> XS <<<<<<
        # Get all 0D, 1D, and 3D xs
        for i in range(3):
            fid.readline()
        xsDict = {}
        read1D, readXfer = 0, 0
        while read1D < num1D or readXfer < numXfer:
            t = fid.readline().split()
            if not t:
                print 'File said it contained {0} cross sections and {1} transfer matrices, but only contained {2} and {3}, respectively.'.format(num1D, numXfer, read1D, readXfer)
                break
            MT = int(t[1].strip(','))
            if MT in [457, 458]:
                # Read a 0D value (halflife or energy per fission)
                read1D += 1
                t = fid.readline().split()
                xsDict[MT] = float(t[0])
            elif MT == 1054:
                # Read delay neutron decay constant. 1D value, size = numDNGs
                read1D += 1
                line = fid.readline().split()
                numDNGs = int(line[5])
                xsDict[MT] = np.zeros(numDNGs)
                read_fixed_line(xsDict[MT], numDNGs, entriesPerLine, fid)
            elif MT == 2055:
                # Read delay neutron spectra. numDNGs 1D value's
                read1D += 1
                line = fid.readline().split()
                numDNGs = int(line[5])
                xsDict[MT] = np.zeros((numDNGs, numGroups))
                for iDNGs in range(numDNGs):
                    line = fid.readline().split()
                    assert iDNGs == int(line[1])
                    sliceChi = xsDict[MT][iDNGs, :]
                    read_fixed_line(sliceChi, numGroups, entriesPerLine, fid)
            elif MT < 2500:
                # Read a 1D cross section
                read1D += 1
                xsDict[MT] = np.zeros(numGroups)
                read_fixed_line(xsDict[MT], numGroups, entriesPerLine, fid)
            elif MT == 2518:
                # Read total fission matrix which only has 0-th moment
                readXfer += 1
                xsDict[MT] = np.zeros((numGroups, numGroups))
                fissionXfer = xsDict[MT]
                for g2 in range(numGroups):
                    t = fid.readline().split()
                    sink, first, last = int(t[3]), int(t[4]), int(t[5])
                    if last < first:
                        fid.readline()
                    else:
                        sliceFission = fissionXfer[sink, first:last+1]
                        read_fixed_line(sliceFission, last-first+1, entriesPerLine, fid)
            else:
                # Read a 3D transfer matrix
                readXfer += 1
                xsDict[MT] = np.zeros((numMoments, numGroups, numGroups))
                # Index scatXfer by [moment, group to, group from]. Uses aliasing
                scatXfer = xsDict[MT]                
                for m in range(numMoments):
                    for g2 in range(numGroups):
                        t = fid.readline().split()
                        sink, first, last = int(t[3]), int(t[4]), int(t[5])
                        if last < first:
                            # No data for this row of the matrix
                            fid.readline()
                        else:
                            sliceScat = scatXfer[m, sink, first:last+1]
                            read_fixed_line(sliceScat, last-first+1, entriesPerLine, fid)
                    if m < (numMoments-1):
                        fid.readline()
    return PDT_XS(numGroups, numMoments, numDNGs, temperature, xsType, xsUnits, groupBdrs, groupWidths, xsDict)

def write_PDT_xs_generally(filePath, xsDat, fromStr='barnfire'):
    temperatureList = [xsDat.T]
    write_PDT_xs_header(filePath, xsDat, temperatureList, fromStr)
    write_PDT_xs_body(filePath, xsDat)
    
def write_PDT_xs_header(filePath, xsDat, temperatureList=[], fromStr='barnfire'):
    '''Write a PDT XS from a PDT_XS object'''
    # Get XS meta-information
    timeStr = datetime.strftime(datetime.now(), '%c')
    numGroups = xsDat.G
    numMoments = xsDat.M
    typeStr = xsDat.typeStr
    microStr = xsDat.microStr
    groupBoundaries = xsDat.Eg
    numTemperatures = len(temperatureList)

    # Print all reactions in xsDat, but print the weight first, if it's included
    mtWgt = 1099
    oneDMTOrder = sorted([key for key in xsDat.xs.keys() if (key != mtWgt and key < 2500)])
    xferMTOrder = sorted([key for key in xsDat.xs.keys() if key >= 2500])
    if mtWgt in xsDat.xs.keys():
        oneDMTOrder.insert(0, 1099)
    num1D = len(oneDMTOrder)
    numXfer = len(xferMTOrder)
    oneDStr = '{0} neutron process'.format(num1D)
    xferStr = '{0} transfer process'.format(numXfer)
    if num1D > 1:
        oneDStr += 'es'
    if numXfer > 1:
        xferStr += 'es'

    # Write XS in PDT format
    with open(filePath, 'w') as fid:
        fid.write('PDT Format Material Data File created {0}\n'.format(timeStr))
        fid.write('\n')
        fid.write('This file is a {0} neutron library generated from {1}.\n'.format(typeStr, fromStr))
        fid.write('{0} temperatures, 1 densities, and {1} groups.\n'.format(numTemperatures, numGroups))
        fid.write('\n')
        fid.write('{0} and {1}.\n'.format(oneDStr, xferStr))
        fid.write('Scattering order {0}\n'.format(numMoments-1))
        fid.write('\n')
        fid.write('{0}\n'.format(microStr))
        fid.write('\n')
        fid.write('Temperatures in Kelvin:\n')
        fid.write(multiline_string(temperatureList, 15, 5, 7))
        fid.write('\n')
        fid.write('Densities in g/cc:\n')
        fid.write('{0:>15}\n'.format(0))
        fid.write('\n')
        fid.write('Group boundaries in eV:\n')
        fid.write(multiline_string(groupBoundaries, 15, 5, 7))
        fid.write('\n') 

def write_PDT_xs_body(filePath, xsDat):
    '''Write a PDT XS from a PDT_XS object'''
    # Get XS meta-information
    timeStr = datetime.strftime(datetime.now(), '%c')
    numGroups = xsDat.G
    numDNGs = xsDat.D
    numMoments = xsDat.M
    temperature = xsDat.T

    # Define special reaction (MT) numbers; zeroDMT are MT numbers that have only 1 value
    zeroDMTList = [457, 458]

    # Print all reactions in xsDat, but print the weight first, if it's included
    mtWgt = 1099
    mtDecayConst = 1054
    mtDelayedChi = 2055
    mtFissionMatrix = 2518
    oneDMTOrder = sorted([key for key in xsDat.xs.keys() if (key not in [mtWgt, mtDecayConst, mtDelayedChi] and key < 2500)])
    xferMTOrder = sorted([key for key in xsDat.xs.keys() if key >= 2500])
    if mtWgt in xsDat.xs.keys():
        oneDMTOrder.insert(0, 1099)
    num1D = len(oneDMTOrder)
    numXfer = len(xferMTOrder)

    # Write XS in PDT format
    with open(filePath, 'a') as fid:
        fid.write('T = {0:g} density = 0\n'.format(temperature))
        fid.write('---------------------------------------------------\n')
        for MT in oneDMTOrder:
            fid.write('MT {0}\n'.format(MT))
            vectorAlias = xsDat.xs[MT]
            if not hasattr(vectorAlias, '__iter__'):
                print MT
                # If MT number corresponds to a 0D quantity, print it as a length-1 vector
                vectorAlias = np.array([vectorAlias])
                print vectorAlias
            fid.write(multiline_string(vectorAlias, 20, 5, 12))

        # write decay constants for delayed neutron groups
        if mtDecayConst in xsDat.xs.keys():
            MT = mtDecayConst
            fid.write('MT {0}\n'.format(MT))
            vectorAlias = xsDat.xs[MT]
            fid.write('  Number of delayed neutron groups: {0}\n'.format(numDNGs))
            fid.write(multiline_string(vectorAlias, 20, 5, 12))
        # write delayed neutron spectra
        if mtDelayedChi in xsDat.xs.keys():
            MT = mtDelayedChi
            fid.write('MT {0}\n'.format(MT))
            vectorAlias = xsDat.xs[MT]
            fid.write('  Number of delayed neutron groups: {0}\n'.format(numDNGs))
            for iDNG in range(numDNGs):
                    fid.write('  DNG {0}\n'.format(iDNG))
                    fid.write(multiline_string(vectorAlias[iDNG,:], 20, 5, 12))
        # write fission matrix
        if mtFissionMatrix in xferMTOrder:
            MT = mtFissionMatrix
            fissionMatrix = xsDat.xs[MT]
            fid.write('MT {0}, Moment {1}\n'.format(MT, 0))
            for g in range(numGroups):
                fid.write('  Sink, first, last: ')
                first = 0
                last = numGroups - 1
                vec = [g, first, last]
                fid.write(multiline_string(vec, 5, 3, 10))
                fid.write(multiline_string(fissionMatrix[:, g], 20, 5, 12))
        # write transfer matrices except for fission matrix
        for MT in [MT for MT in xferMTOrder if MT != mtFissionMatrix]:
            scatMatrix = xsDat.xs[MT]
            for m in range(numMoments):
                fid.write('MT {0}, Moment {1}\n'.format(MT, m))
                for gTo in range(numGroups):
                    scatSlice = scatMatrix[m,gTo,:]
                    nonzeroLeft = np.argmin(scatSlice==0)
                    nonzeroRight = numGroups - np.argmin(scatSlice[::-1]==0)
                    fid.write('  Sink, first, last: ')
                    if all(scatSlice==0):
                        vec = [gTo, -1, -2]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write('\n')
                    else:
                        vec = [gTo, nonzeroLeft, nonzeroRight-1]
                        fid.write(multiline_string(vec, 5, 3, 10))
                        fid.write(multiline_string(scatSlice[nonzeroLeft:nonzeroRight], 20, 5, 12))
        fid.write('\n')

#########################################################################################
def read_fixed_line(obj, objSize, numPerLine, fid):
    """
    Reads into obj from a file using readline(), where the file has
    at most numPerLine values per line. The final size of obj is objSize.

    obj is returned as a numpy array.
    """

    loc = 0
    requiredLines = objSize / numPerLine  # integer math!
    lastLineSize =  objSize % numPerLine
    for L in range(requiredLines):
        t = fid.readline().split()
        for i in range(5):
            obj[loc] = t[i]
            loc += 1
    if lastLineSize > 0:
        t = fid.readline().split()
#        print t
        for i in range(lastLineSize):
#            print i, t[i]
            obj[loc] = t[i]
            loc += 1

def multiline_string(vector, spacing, numberPerLine, decimals):
    outStr = ''
    N = int(np.ceil(len(vector)/float(numberPerLine)))*numberPerLine
    for i in range(numberPerLine, N+1, numberPerLine):
        strs = ['{0:>{1}.{2}g}'.format(vi, spacing, decimals) for vi in vector[i-numberPerLine:i]]
        outStr += ''.join(strs) + '\n'
    return outStr

#########################################################################################
class PDT_XS():
    def __init__(self, numGroups, numMoments,  numDNGs, temperature, typeStr, microStr, groupBdrs, groupWidths, xsDict):
        self.G = numGroups
        self.M = numMoments
        self.D = numDNGs
        self.T = temperature
        self.typeStr = typeStr
        self.microStr = microStr
        self.Eg = groupBdrs
        self.dE = groupWidths
        self.xs = xsDict

    def print_stats(self):
        print 'numGroups numMoments temperature type(MG/MB) type(micro/macro)'
        print self.G, self.M, self.D, self.T, self.typeStr.lower(), self.microStr.lower().split()[0]
        print 'MT list'
        print sorted(self.xs.keys())


#########################################################################################
def print_PDT_MT_enum():
    print '''
    These are the MT numbers used in PDT. The rule to move from PDT_MT to (MF,MT) used in ENDF is:
    if PDT_MT < 1000:
        MF = 3
        MT = PDT_MT
    elif PDT_MT >= 2502:
        MF = 6
        MT = PDT_MT - 2500
    else:
        PDT_MT is a derived quantity does not have an ENDF (MF,MT) classification

  // ===========================================================================
  //    Scalar neutron processes
  // ===========================================================================
  MT_half_life          , // MT = 457,    half life of the nuclide
  MT_E_per_fission      , // MT = 458,    total energy per fission (eV) minus neutrinos
  MT_N_SCALAR_COUNT     , // All neutron scalar values above here
  // ===========================================================================
  //    Single group (1D) neutron processes
  // ===========================================================================

  // common cross sections
  // =====================
  MT_total              , // MT =    1,   total cross section
  MT_elastic            , // MT =    2,   elastic scattering
  MT_nonelastic         , // MT =    3,   nonelastic, sig_t - sig_el
  MT_inelastic          , // MT =    4,   inelastic scattering
  MT_transfer           , // MT =    5,   transfer (sum over final group)
  MT_loss               , // MT =   15,   total - transfer
  MT_absorption         , // MT =   27,   absorption
  MT_n2n                , // MT =   16,   (n, 2n)
  MT_n3n                , // MT =   17,   (n, 3n)
  MT_n4n                , // MT =   37,   (n, 4n)
  MT_n_nalpha           , // MT =   22,   (n, n+alpha)
  MT_n_np               , // MT =   28,   (n, n+p)
  MT_n_gamma            , // MT =  102,   (n, gamma)
  MT_n_p                , // MT =  103,   (n, proton)
  MT_n_alpha            , // MT =  107,   (n, alpha)
  MT_n_disappear        , // MT =  101,   disappearance (no exit neutron)
  MT_inv_velocity       , // MT =  259,   flux-weighted inverse velocity

  MT_weight_func        , // MT = 1099,   group-averaged weight function

  // less-common cross sections
  // ==========================
  //MT_n_n3alpha          , // MT =   23,   (n, n+3alpha)
  //MT_n_2nalpha          , // MT =   24,   (n, 2n+alpha)
  //MT_n_3nalpha          , // MT =   25,   (n, 3n+alpha)
  //MT_n_n2alpha          , // MT =   23,   (n, n+2alpha)
  //MT_n_n2alpha          , // MT =   29,   (n, n+2alpha)
  //MT_n_2n2alpha         , // MT =   30,   (n, 2n+2alpha)
  //MT_n_ndeuteron        , // MT =   32,   (n, n+deuteron)
  //MT_n_ntriton          , // MT =   33,   (n, n+triton)
  //MT_n_nhe3             , // MT =   34,   (n, n+3He)
  //MT_n_ndeuteron2alpha  , // MT =   35,   (n, n+deuteron+2alpha)
  //MT_n_ntriton2alpha    , // MT =   36,   (n, n+triton+2alpha)
  MT_n_deuteron         , // MT =  104,   (n, deuteron)
  MT_n_triton           , // MT =  105,   (n, triton)
  //MT_n_he3              , // MT =  106,   (n, 3He)
  // To add more MT numbers, see
  // www.nndc.bnl.gov/exfor/help7.jsp

  // fission related cross sections
  // ==============================
  MT_nu_sig_f           , // MT = 1452,
  MT_fission            , // MT =   18,   total fission
  MT_nubar              , // MT =  452,   total nubar, average # n0 per fission
  MT_chi                , // MT = 1018,   total fission spectrum

  MT_lambda_del         , // MT = 1054,   decay constants of delayed neutron precursor
  MT_nubar_del          , // MT =  455,   nubar, delayed neutrons
  MT_chi_del            , // MT = 1055,   delayed neutron spectrum
  MT_chis_del           , // MT = 2055,   delayed neutron spectra for all delayed neutron groups
  MT_nubar_prompt       , // MT =  456,   nubar, prompt neutrons
  MT_chi_prompt         , // MT = 1056,   prompt neutron spectrum

  MT_f_first            , // MT =   19,   first chance fission (n, f)
  MT_nubar_p1           , // MT = 4561,   prompt nubar, first chance fission
  MT_chi_1              , // MT = 1019,   first chance fission spectrum

  MT_f_second           , // MT =   20,   second chance fission (n, n'f)
  MT_nubar_p2           , // MT = 4562,   prompt nubar, second chance fission
  MT_chi_2              , // MT = 1020,   second chance fission spectrum

  MT_f_third            , // MT =   21,   third chance fission (n, 2n'f)
  MT_nubar_p3           , // MT = 4563,   prompt nubar, third chance fission
  MT_chi_3              , // MT = 1021,   third chance fission spectrum

  MT_f_fourth           , // MT =   38,   fourth chance fission (n, 3n'f)
  MT_nubar_p4           , // MT = 4564,   prompt nubar, fourth chance fission
  MT_chi_4              , // MT = 1038,   fourth chance fission spectrum

  // inelastic scattering by discrete post-interaction nuclear state
  // ===============================================================
  MT_in_1               , // MT =   51,   inelastic, 1st level
  MT_in_2               , // MT =   52,   inelastic, 2nd level
  MT_in_3               , // MT =   53,   inelastic, 3rd ...
  MT_in_4               , // MT =   54,   inelastic, 4th ...
  MT_in_5               , // MT =   55,   inelastic, 5th ...

  MT_in_6 , MT_in_7 , MT_in_8 , MT_in_9 , MT_in_10,   // MT = level + 50
  MT_in_11, MT_in_12, MT_in_13, MT_in_14, MT_in_15,
  MT_in_16, MT_in_17, MT_in_18, MT_in_19, MT_in_20,
  MT_in_21, MT_in_22, MT_in_23, MT_in_24, MT_in_25,
  MT_in_26, MT_in_27, MT_in_28, MT_in_29, MT_in_30,
  MT_in_31, MT_in_32, MT_in_33, MT_in_34, MT_in_35,
  MT_in_36, MT_in_37, MT_in_38, MT_in_39, MT_in_40,
  MT_in_cont            , // MT =   91,   inelastic continuum not covered above

  // 1D thermal scattering xs
  // =====================================
  MT_th_free            , // MT = 221,  free gas model
  MT_th_h2o             , // MT = 222,  H in H2O
  MT_th_poly_incoh      , // MT = 223,  H in polyethylene (CH2) incoherent
  MT_th_poly_coh        , // MT = 224,  H in polyethylene (CH2) coherent
  MT_th_zrhyd_h_incoh   , // MT = 225,  H in ZrH incoherent
  MT_th_zrhyd_h_coh     , // MT = 226,  H in ZrH coherent
  MT_th_benz_incoh      , // MT = 227,  benzene incoherent
  MT_th_d2o             , // MT = 228,  D in D2O
  MT_th_graph_incoh     , // MT = 229,  C in graphite incoherent
  MT_th_graph_coh       , // MT = 230,  C in graphite coherent
  MT_th_be_incoh        , // MT = 231,  Be metal incoherent
  MT_th_be_coh          , // MT = 232,  Be metal coherent
  MT_th_beo_incoh       , // MT = 233,  BeO incoherent
  MT_th_beo_coh         , // MT = 234,  BeO coherent
  MT_th_zrhyd_zr_incoh  , // MT = 235,  Zr in ZrH incoherent
  MT_th_zrhyd_zr_coh    , // MT = 236,  Zr in ZrH coherent

  // ===========================================================================
  //    Transfer (group-to-group) processes - Neutron AND Gamma combined
  // ===========================================================================
  MT_x_transfer         , // MT = 2500,   total transfer (group to group)
  MT_x_scatter          , // MT = 2501,   total scattering transfer
  MT_x_not_fission      , // MT = 2519,   all transfer except fission

  // ===========================================================================
  //    Transfer (group-to-group) processes - Neutron
  // ===========================================================================
  MT_x_elastic          , // MT = 2502,   elastic scattering
  MT_x_inelastic        , // MT = 2504,   inelastic scattering

  MT_x_n2n              , // MT = 2516,   (n, 2n)
  MT_x_n3n              , // MT = 2517,   (n, 3n)

  MT_x_fission          , // MT = 2518,   total fission transfer matrix (chi and nusigf)

  // inelastic scattering by discrete post-interaction nuclear state
  // ===============================================================
  MT_x_1               , // MT = 2551,   inelastic, 1st level
  MT_x_2               , // MT = 2552,   inelastic, 2nd level
  MT_x_3               , // MT = 2553,   inelastic, 3rd ...
  MT_x_4               , // MT = 2554,   inelastic, 4th ...
  MT_x_5               , // MT = 2555,   inelastic, 5th ...

  MT_x_6 , MT_x_7 , MT_x_8 , MT_x_9 , MT_x_10,   // MT = level + 2550
  MT_x_11, MT_x_12, MT_x_13, MT_x_14, MT_x_15,
  MT_x_16, MT_x_17, MT_x_18, MT_x_19, MT_x_20,
  MT_x_21, MT_x_22, MT_x_23, MT_x_24, MT_x_25,
  MT_x_26, MT_x_27, MT_x_28, MT_x_29, MT_x_30,
  MT_x_31, MT_x_32, MT_x_33, MT_x_34, MT_x_35,
  MT_x_36, MT_x_37, MT_x_38, MT_x_39, MT_x_40,
  MT_x_cont            , // MT = 2591,   inelastic continuum not covered above

  // transfer thermal scattering processes
  // =====================================
  MT_x_th_free          , // MT = 2721,  free gas model
  MT_x_th_h2o           , // MT = 2722,  H in H2O
  MT_x_th_poly_incoh    , // MT = 2723,  H in polyethylene (CH2) incoherent
  MT_x_th_poly_coh      , // MT = 2724,  H in polyethylene (CH2) coherent
  MT_x_th_zrhyd_h_incoh , // MT = 2725,  H in ZrH incoherent
  MT_x_th_zrhyd_h_coh   , // MT = 2726,  H in ZrH coherent
  MT_x_th_benz_incoh    , // MT = 2727,  benzene incoherent
  MT_x_th_d2o           , // MT = 2728,  D in D2O
  MT_x_th_graph_incoh   , // MT = 2729,  C in graphite incoherent
  MT_x_th_graph_coh     , // MT = 2730,  C in graphite coherent
  MT_x_th_be_incoh      , // MT = 2731,  Be metal incoherent
  MT_x_th_be_coh        , // MT = 2732,  Be metal coherent
  MT_x_th_beo_incoh     , // MT = 2733,  BeO incoherent
  MT_x_th_beo_coh       , // MT = 2734,  BeO coherent
  MT_x_th_zrhyd_zr_incoh, // MT = 2735,  Zr in ZrH incoherent
  MT_x_th_zrhyd_zr_coh  , // MT = 2736,  Zr in ZrH coherent
'''

