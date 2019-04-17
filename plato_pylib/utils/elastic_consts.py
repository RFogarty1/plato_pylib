

import itertools as it
import numpy as np
import plato_pylib.shared.ucell_class as UCell

from scipy.optimize import minimize 
#Purpose of these functions are to help calculate elastic constants


_STRAIN_MATRIX_DICT = dict()


#---------->Functions to calculate the elastic constant from output calcs<----------

def calcElasticConstantFromStrainEnergies(strainVals, energyVals, refVol, n, m, nElastic=None, mElastic=None):
	#Step 1 = calculate 2nd deriv of strain vs energy curve
	objFunct = _createObjFunctFor2ndDerivFit(strainVals,energyVals)
	startVals = [0.0 for x in range(3)]
	fitRes = minimize(objFunct, startVals)
	print("Got a value of {} for the 2nd deriv".format(fitRes.x[0]))

	#Step 2 = calculate elastic constant based on that info
	secondDeriv = fitRes.x[0]
	leftSide = fitRes.x[0] * (1/refVol)
	if n==m:
		elastic = leftSide
	else:
		elastic = leftSide - mElastic - nElastic

	#Step 3 = create object with relevant info
	outObj = FitResultElastic( fitRes, elastic, refVol, (n,m) )
	return outObj

#For now limit to x**2 polynomial fit for the minimize interface
def _createObjFunctFor2ndDerivFit(strainVals,energies):
	def objFunct(params):
		predEnergies = list()
		for x in strainVals:
			result = (params[0]*(x**2)) + (params[1]*x) + (params[2])
			predEnergies.append(result)
		outVal = sum( [ (b-a)**2 for b,a in it.zip_longest(energies, predEnergies)] )
		return outVal
	return objFunct



class FitResultElastic:
	def __init__(self, fitObj, elasticConst, refVol, nmTuple):
		self.fitObj = fitObj
		self.elasticConst = elasticConst
		self.refVol = refVol
		self.nm = nmTuple


#---------->Functions related to generate the strained files<-------------------

def applyStrainToUCellForElasticConstant(unitCell, n, m, strainCoeff):
	''' Applies a strain to UnitCell class structure needed to get elastic constant C_nm '''
	currLattVects = unitCell.lattVects
	origFractCoords = [list(x) for x in unitCell.fractCoords]
	applyStrainForElasticConstant(currLattVects, n, m, strainCoeff)
	unitCell.lattVects = currLattVects
	unitCell.fractCoords = origFractCoords



def applyStrainForElasticConstant(cartVects:"3x3 iterxiter", n, m, strainCoeff):
	''' Applies strain needed to calculate elastic constant Cnm (works in place) '''
	strainMatrix = getStrainMatrix(n,m,strainCoeff=strainCoeff)
	vectShiftMatrix = (np.identity(3) + strainMatrix) @ np.array(cartVects)

	for idxA,iterA in enumerate(cartVects):
		for idxB,unused in enumerate(iterA):
			cartVects[idxA][idxB] = vectShiftMatrix[idxA][idxB]



def getStrainMatrix(n, m, strainCoeff=1):
	''' Get the relevant strain matrix to calculate C_{n,m} are Voight indices, see plato section on elastic constants '''
	if (n==m):
		return _STRAIN_MATRIX_DICT[(n,m)](strainCoeff)

	matN = _STRAIN_MATRIX_DICT[(n,n)](strainCoeff)
	matM = _STRAIN_MATRIX_DICT[(m,m)](strainCoeff)
	return matN + matM

	raise ValueError("C_nm with n={}, m={} is an invalid elastic constant".format(n,m))


def registerStrainMatrix(key:"2-tuple"):
	def decorate(funct):
		_STRAIN_MATRIX_DICT[key] = funct
		return funct
	return decorate


@registerStrainMatrix( (1,1) )
def _strainMatrix11(strainCoeff):
	outArray = np.array( ([strainCoeff,0.0,0.0],
	                      [0.0,0.0,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray


@registerStrainMatrix( (2,2) )
def _strainMatrix22(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,strainCoeff,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray

@registerStrainMatrix( (3,3) )
def _strainMatrix33(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,0.0,0.0],
	                      [0.0,0.0,strainCoeff]) )
	return outArray



@registerStrainMatrix( (4,4) )
def _strainMatrix44(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,0.0,0.5*strainCoeff],
	                      [0.0,0.5*strainCoeff,0.0]) )
	return outArray


@registerStrainMatrix( (5,5) )
def _strainMatrix55(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.5*strainCoeff],
	                      [0.0,0.0,0.0],
	                      [0.5*strainCoeff,0.0,0.0]) )
	return outArray



@registerStrainMatrix( (6,6) )
def _strainMatrix66(strainCoeff):
	outArray = np.array( ([0.0,0.5*strainCoeff,0.0],
	                      [0.5*strainCoeff,0.0,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray


