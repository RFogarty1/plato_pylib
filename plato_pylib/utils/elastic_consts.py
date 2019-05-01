
import copy
import itertools as it
import numpy as np
import plato_pylib.shared.ucell_class as UCell

from scipy.optimize import minimize 
#Purpose of these functions are to help calculate elastic constants


_STRAIN_MATRIX_DICT = dict()




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
	tranposedVects = np.array(lattVects).transpose()
	vectShiftMatrix = ( (np.identity(3) + strainMatrix) @ tranposedVects )
	vectShiftMatrix = vectShiftMatrix.transpose()

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
	outMatrices.append( _STRAIN_MATRIX_DICT[6](2*strainParam) )

	print("strainParam = {}".format(strainParam))
	print("First matrix is:")
	print(outMatrices[0])
	print("Second matrix is:")
	print(outMatrices[1])
	print("Third matrix is:")
	print(outMatrices[2])
	print("Fourth matrix is:")
	print(outMatrices[3])
	print("Fifth matrix is:")
	print(outMatrices[4])


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

	return ElasticOutputData(outDict,polyFits)



class ElasticOutputData:
	def __init__(self, elasticDict, polyFits):
		self.fits = polyFits
		self.elasticConstants = elasticDict

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





def getElasticConstCoeffMatrix(crystType):
	typeToFunct = {"hexagonal":_hexElasticEquationMatrix}
	try:
		return typeToFunct[crystType.lower()]()
	except KeyError:
		raise ValueError("{} is an invalid crystal type, allowedVals are {}".format(structType,typeToFunct.keys()))




def _hexElasticKeys():
	return [(1,1), (1,2), (1,3), (3,3), (4,4)]

def _hexElasticEquationMatrix():
	return np.array( [ [0, 0, 0, 1, 0],
	                   [2, 2, 0, 0, 0],
	                   [2, 2, 4, 1, 0],
	                   [0, 0, 0, 0, 8],
	                   [2,-2, 0, 0, 0] ] )




#------->Strain matrices<---------------

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


