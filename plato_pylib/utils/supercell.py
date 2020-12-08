#!/usr/bin/env python3

import copy
import itertools as it
import math
import numpy as np
import plato_pylib.shared.ucell_class as UCell


def superCellFromUCell(unitCell, dims):

	xShifted = _getSuperCellOneDim(unitCell,0,dims[0])
	yShifted = _getSuperCellOneDim(xShifted,1,dims[1])
	zShifted = _getSuperCellOneDim(yShifted,2,dims[2])

	return zShifted



def getUnitCellSurroundedByNeighbourCells(unitCell, alongA=True, alongB=True, alongC=True, removeStartCoords=False):
	""" Gets the unit cell surrouned by images in each direction; meaning a total of 27 cells merged into a supercell. Note that the original cells cartesian co-ordinates will be unchanged; meaning that most atoms will have -ve x/y/z co-ordinates
	
	Args:
		unitCell: plato_pylib UnitCell object
		alongA (Optional, Bool): Whether to get the images along lattice vector a (default=True)
		alongB (Optional, Bool): Whether to get the images along lattice vector b (default=True)
		alongC (Optional, Bool): Whether to get the images along lattice vector c (default=True)
		removeStartCoords (Optional, Bool): Option to remove the co-ordinates from the original cell; useful if your looking for ONLY the neighbour co-ords
			 
	Returns
		outCell: a DIFFERENT plato_pylib UnitCell object containing the original surrounded by 26 image cells in total
 
	Raises:
	"""

	startCoords = copy.deepcopy(unitCell.cartCoords)
	nStartCoords = len(startCoords)
	xCell = _getCellWithImageAddedEachSide(unitCell,0) if alongA is True else unitCell
	yCell = _getCellWithImageAddedEachSide(xCell,1) if alongB is True else xCell
	outCell = _getCellWithImageAddedEachSide(yCell,2) if alongC is True else yCell

	if removeStartCoords:
		outCoords = list()
		tempCoords = copy.deepcopy(outCell.cartCoords)
		indicesToSkip = list()
		diffTol = 1e-4
		for idx,coord in enumerate(startCoords):
			diffVals = [abs(x-y) for x,y in zip(tempCoords[idx][:3],startCoords[idx][:3])]
			assert all([x<1e-4 for x in diffVals])
			indicesToSkip.append(idx)
		for idx,coord in enumerate(tempCoords):
			if idx not in indicesToSkip:
				outCoords.append(tempCoords[idx])
		outCell.cartCoords = outCoords

	return outCell


def _getCellWithImageAddedEachSide(startCell, dimIdx:"0,1,2 for x,y,z"):
	#Get translation vectors
	startCartCoords = copy.deepcopy( startCell.cartCoords )
	translationVector = startCell.lattVects[dimIdx]

	#get the new cartesian coordinates for both directions; +ve then negative
	positiveNewCarts, negativeNewCarts = list(), list()
	for coord in startCartCoords:
		currPositive = [a+b for a,b in it.zip_longest(coord[:3],translationVector[:3])] + [coord[-1]]
		currNegative = [a-b for a,b in it.zip_longest(coord[:3],translationVector[:3])] + [coord[-1]]
		positiveNewCarts.append(currPositive)
		negativeNewCarts.append(currNegative)

	newCartCoords = startCartCoords + positiveNewCarts + negativeNewCarts
 
	#Figure out the new lattice parameters
	lattParams = copy.deepcopy( startCell.getLattParamsList() )
	lattParams[dimIdx] *= 3
	lattAngles = copy.deepcopy( startCell.getLattAnglesList() )

	#Return a (new) unitCell with all these new values
	outCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = newCartCoords
	return outCell

def _getSuperCellOneDim(unitCell,dimIdx:"0,1,2 for x,y,z", multiple:"int, number of cells in this dimension (1 means do nothing)"):

	#step 1 = get the unit-vectors for the unit-cell
	startCellVects = copy.deepcopy(unitCell.lattVects)
	startCellParams = unitCell.getLattParamsList()
	unitVects = [_getUnitVectorFromVector(x) for x in startCellVects]
	startCartCoords = unitCell.cartCoords

	#Step 2 = figure out the translation vector required
	transVector = [startCellParams[dimIdx]*x for x in unitVects[dimIdx]]

	#Step 3 = Get the cartesian co-ordinates for new atoms (translate n-times from the translation vector)
	endCartCoords = [list(x) for x in startCartCoords]

	for atom in startCartCoords:
		for mult in range(1,multiple):
			newAtom = [x + (t*mult) for x,t in it.zip_longest(atom[:3],transVector)]
			endCartCoords.append(newAtom + [atom[3]])

	#Step 4 = Get new cell vectors
	newCellVects = list()
	for idx,(vect,uVect) in enumerate( it.zip_longest(startCellVects,unitVects) ):
		if idx == dimIdx:
			currVect = [x+(ux*(multiple-1)*startCellParams[dimIdx]) for x,ux in zip(vect,uVect)]
		else:
			currVect = list(vect)
		newCellVects.append(currVect)


	#Step 6 = Create the output object
	outCell = UCell.UnitCell.fromLattVects( newCellVects)
	outCell.cartCoords = endCartCoords
	outCell.putCAlongZ = unitCell.putCAlongZ

	return outCell



def _getUnitVectorFromVector(vector:"iter"):
	lenVect = math.sqrt( sum([x**2 for x in vector]) )
	return [x/lenVect for x in vector]



