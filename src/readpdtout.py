#! /usr/bin/env python

'''
Andrew Till
Spring 2014

Attempt to read the PDT .output_* format.
Steady state only.
'''

#TPL
import numpy as np

def call_example(filename):
    dat = read_pdt_output(filename)
    flatDat = make_flat_pdt_output(dat)
    print dat.numCells, dat.numGroups, dat.numTotElem, dat.eigval, dat.geoType, dat.probType

def make_flat_pdt_output(dat):
    '''Returns flattened data structure of PDTOutput'''
    xcoord = np.zeros(dat.numTotElem)
    ycoord = np.zeros(dat.numTotElem)
    zcoord = np.zeros(dat.numTotElem)
    vol = np.zeros(dat.numTotElem)
    flux = np.zeros((dat.numTotElem, dat.numGroups))
    ielem = 0
    for cell in dat.cells:
        for elem in cell.elems:
            xcoord[ielem] = elem.coords[0]
            ycoord[ielem] = elem.coords[1]
            zcoord[ielem] = elem.coords[2]
            vol[ielem] = elem.vol
            flux[ielem, :] = elem.flux
            ielem += 1
    flatDat = PDTOutputFlat(xcoord, ycoord, zcoord, vol, flux)
    return flatDat

def read_pdt_output(filename):
    '''Read .output_? file'''
    finished = False
    numCells = 0
    numGroups = 0
    numMoments = 0
    numTotElem = 0
    eigenvalue = 0
    probType = ''
    geoType = ''
    with open(filename, 'r') as fid:
        # Get number of groups, problem type, etc.
        line = fid.readline()
        t = fid.readline().strip().split()
        probType = t[0]
        geoType = t[2]
        numGroups = int(t[3])
        numMoments = int(t[4])
        if probType == 'K_EIGEN':
            line = fid.readline()
            line = fid.readline()
            eigenvalue = float(fid.readline())
        # Loop through all cells
        cells = []
        while not(finished):
            lineFull = fid.readline()
            lineStrip = lineFull.strip()
            if not(lineFull):
                finished = True
            elif lineStrip == 'BEGIN_RECORD TYPE = CELL':
                #Read a cell
                numCells += 1
                # Read cell geometric info
                line = fid.readline()
                t = fid.readline().strip().split()
                cellID = int(t[0])
                cellCenters = [float(t[1]), float(t[2]), float(t[3])]
                cellVol = float(t[4])
                numElem = int(t[5])
                numTotElem += numElem
                # Read elements' geometric info
                geoElements = []
                line = fid.readline()
                for ielem in range(numElem):
                    t = fid.readline().strip().split()
                    elemCoords = [float(t[1]), float(t[2]), float(t[3])]
                    elemVol = float(t[4])
                    geoElements.append(GeoElement(elemCoords, elemVol))
                # Read cells' flux info
                line = fid.readline()
                cellFlux = np.zeros(numGroups)
                for igroup in range(numGroups):
                    t = fid.readline().strip().split()
                    cellFlux[igroup] = float(t[2])
                for igroup in range(numGroups * (numMoments - 1)):
                    line = fid.readline()
                # Read elements' flux info
                fluxElements = []
                for ielem in range(numElem):
                    line = fid.readline()
                    elementFlux = np.zeros(numGroups)
                    for igroup in range(numGroups):
                        t = fid.readline().strip().split()
                        elementFlux[igroup] = float(t[2])
                    for igroup in range(numGroups * (numMoments - 1)):
                        line = fid.readline()
                    fluxElements.append(elementFlux)
                # Assemble elements
                elements = []
                for geoElem, fluxElem in zip(geoElements, fluxElements):
                    elements.append(Element(geoElem.coords, geoElem.vol, fluxElem))
                # Assemble cell
                cell = Cell(cellID, cellCenters, cellVol, numElem, elements, cellFlux)
                cells.append(cell)
    # Assemble output
    dat = PDTOutput(cells, numCells, numGroups, numTotElem, probType, geoType, eigenvalue)
    return dat


class PDTOutputFlat():
    def __init__(self, x, y, z, vol, flux):
        self.x = x
        self.y = y
        self.z = z
        self.vol = vol
        self.flux = flux

class PDTOutput():
    def __init__(self, cells, numCells, numGroups, numTotElem, probType, geoType, eigval):
        self.cells = cells
        self.numCells = numCells
        self.numGroups = numGroups
        self.numTotElem = numTotElem
        self.probType = probType
        self.geoType = geoType
        self.eigval = eigval

class Cell():
    def __init__(self, ID, centers, vol, numElem, elems, flux):
        self.ID = ID
        self.centers = centers
        self.vol = vol
        self.elems = elems
        self.numElem = numElem
        self.flux = flux

class GeoElement():
    def __init__(self, coords, vol):
        self.coords = coords
        self.vol = vol

class Element():
    def __init__(self, coords, vol, flux):
        self.coords = coords
        self.vol = vol
        self.flux = flux

