#! /usr/bin/env python
'''
Andrew Till

Read in 1-D XS and group boundary files

'''

import os
import string
import numpy as np

def read_xs_file(filename, Elower, Eupper):
    '''Reads a xs file and returns energy domain portion at least [Elower,Eupper]'''
    energy, sigma = np.loadtxt(filename, delimiter=',', skiprows=1, unpack=True)
    lowerIndex = np.argmin(np.abs(energy - Elower))
    upperIndex = np.argmin(np.abs(energy - Eupper))
    if energy[lowerIndex] > Elower:
        lowerIndex -= 1
    if energy[upperIndex] < Eupper:
        upperIndex += 1
    return (energy[lowerIndex:upperIndex+1], sigma[lowerIndex:upperIndex+1])

def read_group_file(filename):
    '''Reads a group boundary file'''
    groupBdr = np.loadtxt(filename, skiprows=2, usecols=[1])
    return groupBdr
    
class EnergyRange():
    def __init__(self, Elower, Eupper, EshortName, ElongName=''):
        self.Elower = Elower
        self.Eupper = Eupper
        self.EshortName = EshortName
        self.ElongName = ElongName

class XSType():
    def __init__(self, nuclide, rxn):
        self.nuclide = nuclide
        self.rxn = rxn
