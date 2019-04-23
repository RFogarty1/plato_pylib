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


def _getSuperCellOneDim(unitCell,dimIdx:"0,1,2 for x,y,z", multiple:"int, number of cells in this dimension (1 means do nothing)"):

	#step 1 = get the unit-vectors for the unit-cell
	startCellVects = copy.deepcopy(unitCell.lattVects)
	startCellParams = unitCell.getLattParamsList()
	unitVects = [_getUnitVectorFromVector(x) for x in startCellVects]
	startFractCoords = unitCell.fractCoords
	startCartCoords = UCell.getCartCoordsFromFractCoords(startCellVects, startFractCoords)

	#Step 2 = figure out the translation vector required
	transVector = [startCellParams[dimIdx]*x for x in unitVects[dimIdx]]

	#Step 3 = Get the cartesian co-ordinates for new atoms (translate n-times from the translation vector)
	endCartCoords = UCell.getCartCoordsFromFractCoords(startCellVects, startFractCoords)
	
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

	#Step 5 = Get new fractional co-ordinates
	fractCoords = getFractCoordsFromCartCoords(newCellVects,endCartCoords)

	#Step 6 = Create the output object
	return UCell.UnitCell.fromLattVects( newCellVects, fractCoords=fractCoords)



def _getUnitVectorFromVector(vector:"iter"):
	lenVect = math.sqrt( sum([x**2 for x in vector]) )
	return [x/lenVect for x in vector]



def getFractCoordsFromCartCoords(lattVects, cartCoords):
	outList = list()
	for currAtom in cartCoords:
		outList.append( _getFractCoordsFromCartOneAtom(lattVects, currAtom[:3]) + [currAtom[-1]] ) 
	return outList

def _getFractCoordsFromCartOneAtom(lattVects:"3x3 iterxiter", cartCoords:"3 item list"):
	vMatrix = np.array(lattVects).transpose()

	fractVals = ( np.linalg.inv(vMatrix) @ np.array(cartCoords) ).tolist()
	return fractVals

