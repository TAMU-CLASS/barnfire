'''
Remove this!
'''

#TODO: Have 90/20 be input parameter
#TODO: Use consistent iT? (Or warn about it?)
    def print_contents(self, verbosity=0):
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
            print '(Warning: x in .9xc and .2xt may not be accurate)'
    ################################################################
    # TODO: Populate material-local SabDict[(Z,A)]
    # Make thermalXSDict depend on Z and A
    # 
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
        # TODO: When making SabDict, use util's dicts to make this better. Allow for more 'none' options
        # Combine all 'none's into SabDict[(Z,A)] == 'none'
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




################################################################
#materials_util:


# Derivable:
#def get_thermal_name_to_Z_list_dict():
#    '''Returns a dict that maps thermal name to a list of atomic numbers corresponding to elements with bound thermal treatment in that material. Free thermal treatment is not in dict'''
#    return {'h2o': [1], 'uo2': [8,92], 'zrh': [1,40], 'graphite': [6], 'al': [13], 'fe': [26]}






################################################################
#materials_njoy:

    maxNumTemperatures = 10
    allowedTDict = {}
    get_allowed_njoy_temperatures(allowedTDict)
    for (Z, A) in globalZAList:
        allowedTList = set()
        njoyTDict[(Z, A)] = set()
        for thermalXS in globalTXSDict[Z]:
            if thermalXS != 'free':
                allowedTList |= allowedTDict[thermalXS]
        if allowedTList:
            # Use temperature grid where thermal XS's are defined
            sortedAllowedTList = np.array(sorted(allowedTList))
            for T in globalTDict[(Z, A)]:
                njoyTDict[(Z, A)] |= set(util.get_nearest_points(T, sortedAllowedTList))
        else:
            # No grid specified, free to use actual temperatures
            njoyTDict[(Z,A)] = copy.copy(globalTDict[(Z,A)])
        while len(njoyTDict[(Z,A)]) > maxNumTemperatures:
            # This may unintentially preferentially thin one thermal XS grid over another
            # in the case of multiple non-free thermal XS per nuclide
            itemToRemove = util.thin_list(njoyTDict[(Z,A)])
            njoyTDict[(Z,A)].remove(itemToRemove)





    #Using thermal XS names:
    allowedTDict['hh2o'] = Tmod
    allowedTDict['ouo2inel'] = Tfuel
    allowedTDict['ouo2elas'] = Tfuel
    allowedTDict['uuo2inel'] = Tfuel
    allowedTDict['uuo2elas'] = Tfuel
    allowedTDict['hzrhinel'] = Tfuel
    allowedTDict['hzrhelas'] = Tfuel
    allowedTDict['zrzrhinel'] = Tfuel
    allowedTDict['zrzrhelas'] = Tfuel
    allowedTDict['graphinel'] = Tgraphite
    allowedTDict['graphelas'] = Tgraphite
    allowedTDict['alinel'] = Tstructural
    allowedTDict['alelas'] = Tstructural
    allowedTDict['feinel'] = Tstructural
    allowedTDict['feelas'] = Tstructural




get_element_thermal_name_to_endf_filename_dict

def get_thermal_xs_name_to_mat_number_dict():
    '''Returns a dict that maps thermal XS name to thermal MAT number.
    Since MAT numbers are based on element thermal names and each element thermal name contains
    potentially multiple thermal XS names, the same MAT number will be pointed to by potentially
    multiple thermal XS names.'''
    # Derive this dictionary from previous information
    et2mat = get_element_thermal_name_to_mat_number_dict()
    et2txs = get_element_thermal_name_to_thermal_xs_list_dict()
    txs2mat = {}
    for et in et2mat:
        for txs in et2txs[et]:
            txs2mat[txs] = et2mat[et]
    return txs2mat
def get_thermal_xs_name_to_mat_number_dict():
    '''Returns a dict that maps thermal XS name to thermal MAT number.
    Since MAT numbers are based on element thermal names and each element thermal name contains
    potentially multiple thermal XS names, the same MAT number will be pointed to by potentially
    multiple thermal XS names.'''
    # Derive this dictionary from previous information
    et2mat = get_element_thermal_name_to_mat_number_dict()
    et2txs = get_element_thermal_name_to_thermal_xs_list_dict()
    txs2mat = {}
    for et in et2mat:
        for txs in et2txs[et]:
            txs2mat[txs] = et2mat[et]
    return txs2mat
def get_thermal_xs_name_to_mat_number_dict():
    '''Returns a dict that maps thermal XS name to thermal MAT number.
    Since MAT numbers are based on element thermal names and each element thermal name contains
    potentially multiple thermal XS names, the same MAT number will be pointed to by potentially
    multiple thermal XS names.'''
    # Derive this dictionary from previous information
    et2mat = get_element_thermal_name_to_mat_number_dict()
    et2txs = get_element_thermal_name_to_thermal_xs_list_dict()
    txs2mat = {}
    for et in et2mat:
        for txs in et2txs[et]:
            txs2mat[txs] = et2mat[et]
    return txs2mat
def get_thermal_xs_name_to_mat_number_dict():
    '''Returns a dict that maps thermal XS name to thermal MAT number.
    Since MAT numbers are based on element thermal names and each element thermal name contains
    potentially multiple thermal XS names, the same MAT number will be pointed to by potentially
    multiple thermal XS names.'''
    # Derive this dictionary from previous information
    et2mat = get_element_thermal_name_to_mat_number_dict()
    et2txs = get_element_thermal_name_to_thermal_xs_list_dict()
    txs2mat = {}
    for et in et2mat:
        for txs in et2txs[et]:
            txs2mat[txs] = et2mat[et]
    return txs2mat
