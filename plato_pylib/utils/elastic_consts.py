
import copy
import itertools as it
import numpy as np
import plato_pylib.shared.ucell_class as UCell

from scipy.optimize import minimize 
#Purpose of these functions are to help calculate elastic constants


_STRAIN_MATRIX_DICT = dict()




#-------->ATTEMPT TWO STARTS HERE <-------------------------------


#-------->Code to generate the strains<----------------
def getStrainedStructsForElasticConsts(uCell:"UnitCell obj", strainParams:"list", crystType=None):
	if crystType is None:
		raise ValueError("{} is an invalid crystal type (for now)".format(crystType))
	#Pass on to lower functions by just extracting the lattice vectors, which are all thats really needed. Grab a loada shifted lattice vectors at the end
	lattVects = uCell.lattVects
	strainedLattVects = getStrainedLattVects(lattVects, strainParams, crystType)
	outUCells = list()

	for idxA in range(len(strainedLattVects)):
		outUCells.append( list() )
		for idxB in range(len(strainedLattVects[idxA])):
			outUCells[idxA].append( UCell.UnitCell.fromLattVects( strainedLattVects[idxA][idxB], uCell.fractCoords) )

	return outUCells

def getStrainedLattVects(lattVects, strainParams, crystType):
	''' Returns: list of strained lattice vectors for each strain matrix. output[0] contains all strains for the 0-th strain matrix (e.g. strains along c),
	             while output[1] contains all lattice vectors for strains along 
	'''         
	unitStrainMatrices = getStrainMatrices(crystType,1) 

	strainedLattVects = list()
	for idx,sMatrix in enumerate(unitStrainMatrices):
		strainedLattVects.append(list())
		for sParam in strainParams:
			currLattVects = copy.deepcopy(lattVects)
			applyStrainToLattVects(currLattVects,sMatrix*sParam)
			strainedLattVects[idx].append(currLattVects)

	return strainedLattVects


def applyStrainToLattVects(lattVects:"3x3 mutable iter e.g [ [v1],[v2],[v3] ]", strainMatrix):
	''' Works in place '''
	vectShiftMatrix = (np.identity(3) + strainMatrix) @ np.array(lattVects)
	for idxA,iterA in enumerate(lattVects):
		for idxB,unused in enumerate(iterA):
			lattVects[idxA][idxB] = vectShiftMatrix[idxA][idxB]


def getStrainMatrices(structType,strainParam):
	typeToFunct = {"hexagonal":_hexStrainMatrices}
	try:
		return typeToFunct[structType.lower()](strainParam)
	except KeyError:
		raise ValueError("{} is an invalid crystal type, allowedVals are {}".format(structType,typeToFunct.keys()))


def _hexStrainMatrices(strainParam):
	outMatrices = list()
	outMatrices.append( _STRAIN_MATRIX_DICT[3](strainParam) ) 
	outMatrices.append( _STRAIN_MATRIX_DICT[1](strainParam) + _STRAIN_MATRIX_DICT[2](strainParam) )
	outMatrices.append( _STRAIN_MATRIX_DICT[1](strainParam) + _STRAIN_MATRIX_DICT[2](strainParam) + _STRAIN_MATRIX_DICT[3](strainParam) )
	outMatrices.append( _STRAIN_MATRIX_DICT[4](2*strainParam) + _STRAIN_MATRIX_DICT[5](2*strainParam) )
	outMatrices.append( 2* (_STRAIN_MATRIX_DICT[1](strainParam) + _STRAIN_MATRIX_DICT[2](strainParam) + _STRAIN_MATRIX_DICT[6](strainParam)) )

	return outMatrices



#--------->Code to get elastic constants from the correct calculations <----------

def calcElasticsFromStressStain(strainStress:"nx2 iter, e.g. [(1,1)], with elements (strain,stress)", crystType:str):
	polyFits = list()


	for x in strainStress:
		polyFits.append( _polyFitAndGetSecondDeriv(x)  )


	secondDerivs = [x.secondDeriv for x in polyFits]
	elasticKeys = getElasticKeysInOrder(crystType)
	coeffMatrix = getElasticConstCoeffMatrix(crystType)

	elasticConsts = np.linalg.inv(coeffMatrix) @ np.array(secondDerivs)

	outDict = {k:v for k,v in it.zip_longest(elasticKeys, elasticConsts)}

	return outDict





def _polyFitAndGetSecondDeriv(strainStress):
	inpVals = np.array(strainStress)

	normFactor = 1 / max(inpVals[:,1]) 

	objFunct = _createObjFunctFor2ndDerivFit(inpVals[:,0],inpVals[:,1]*normFactor)
	startVals = [0,0,0]
	fitRes = minimize(objFunct, startVals)

	outObj = QuadFitInfo([x/normFactor for x in fitRes.x], fitRes, inpVals[:,0], inpVals[:,1])

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


class QuadFitInfo:
	def __init__(self, params, fitObj, strainVals, stressVals):
		self.strainVals = strainVals
		self.stressVals = stressVals
		self.params = params #[a,b,c] for ax^2,bx+c
		self.fitObj = fitObj

	@property
	def secondDeriv(self):
		return 2*self.params[0]



def getElasticKeysInOrder(crystType):
	typeToFunct = {"hexagonal":_hexElasticKeys}
	try:
		return typeToFunct[crystType.lower()]()
	except KeyError:
		raise ValueError("{} is an invalid crystal type, allowedVals are {}".format(structType,typeToFunct.keys()))

def _hexElasticKeys():
	return [(1,1), (1,2), (1,3), (3,3), (4,4)]




def getElasticConstCoeffMatrix(crystType):
	typeToFunct = {"hexagonal":_hexElasticEquationMatrix}
	try:
		return typeToFunct[crystType.lower()]()
	except KeyError:
		raise ValueError("{} is an invalid crystal type, allowedVals are {}".format(structType,typeToFunct.keys()))



def _hexElasticEquationMatrix():
	return np.array( [ [0, 0, 0, 1, 0],
	                   [2, 2, 0, 0, 0],
	                   [2, 2, 4, 0, 0],
	                   [0, 0, 0, 0, 4],
	                   [6,-6, 0, 0, 0] ] )




#---------->Functions to calculate the elastic constant from output calcs<----------

def calcElasticConstantFromStrainEnergies(strainVals, energyVals, refVol, n, m, nElastic=None, mElastic=None):
	#Step 1 = calculate 2nd deriv of strain vs energy curve
	objFunct = _createObjFunctFor2ndDerivFit(strainVals,energyVals)
	startVals = [0.0 for x in range(3)]

	fitRes = minimize(objFunct, startVals)

	#Step 2 = calculate elastic constant based on that info
	secondDeriv = 2*fitRes.x[0]
	leftSide = secondDeriv * (1/refVol)
	if n==m:
		elastic = leftSide
	else:
		elastic = leftSide - mElastic - nElastic

	#Step 3 = create object with relevant info
	outObj = FitResultElastic( fitRes, elastic, refVol, (n,m), strainVals, energyVals )
	return outObj




class FitResultElastic:
	def __init__(self, fitObj, elasticConst, refVol, nmTuple, strainVals, energyVals):
		self.fitObj = fitObj
		self.elasticConst = elasticConst
		self.refVol = refVol
		self.nm = nmTuple
		self.strainVals = strainVals
		self.energyVals = energyVals

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


def registerStrainMatrix(key:"2-tuple"):
	def decorate(funct):
		_STRAIN_MATRIX_DICT[key] = funct
		return funct
	return decorate


@registerStrainMatrix( 1 )
def _strainMatrix11(strainCoeff):
	outArray = np.array( ([strainCoeff,0.0,0.0],
	                      [0.0,0.0,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray


@registerStrainMatrix( 2 )
def _strainMatrix22(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,strainCoeff,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray

@registerStrainMatrix( 3 )
def _strainMatrix33(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,0.0,0.0],
	                      [0.0,0.0,strainCoeff]) )
	return outArray



@registerStrainMatrix( 4 )
def _strainMatrix44(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.0],
	                      [0.0,0.0,0.5*strainCoeff],
	                      [0.0,0.5*strainCoeff,0.0]) )
	return outArray


@registerStrainMatrix( 5 )
def _strainMatrix55(strainCoeff):
	outArray = np.array( ([0.0,0.0,0.5*strainCoeff],
	                      [0.0,0.0,0.0],
	                      [0.5*strainCoeff,0.0,0.0]) )
	return outArray



@registerStrainMatrix( 6 )
def _strainMatrix66(strainCoeff):
	outArray = np.array( ([0.0,0.5*strainCoeff,0.0],
	                      [0.5*strainCoeff,0.0,0.0],
	                      [0.0,0.0,0.0]) )
	return outArray


