#!/usr/bin/env python3

import itertools as it
import math

import numpy as np


ANG_TO_BOHR = 1.8897259886


def getTransformedFractCoords(origLattVects:"list of lists [[v1],[v2],[v3]]", finalLattVects, origFractCoords:"list of lists"):
	lattVectsOrig = np.array(origLattVects).transpose() #1 column = 1 vector
	lattVectsFinal = np.array(finalLattVects).transpose()

	finalFractCoords = list()
	for fracCoords in origFractCoords:
		currCoords = np.array(fracCoords)
		finalFractCoords.append( (currCoords@lattVectsOrig@np.linalg.inv(lattVectsFinal)).tolist() )

	return finalFractCoords

class UnitCell():
	def __init__(self,**kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.lattParams = self.listToLattParams( kwargs.get( "lattParams".lower(), None ) )
		self.lattAngles = self.listToLattAngles( kwargs.get( "lattAngles".lower(), None ) )
		self._fractCoords = kwargs.get("fractCoords".lower(), None)
		self._elementList = kwargs.get("elementList".lower(), None)

	@classmethod
	def fromLattVects(cls, lattVectors:"iterable, len=3"):
		lattParams, lattAngles = lattParamsAndAnglesFromLattVects(lattVectors)
		return cls(lattParams=lattParams, lattAngles = lattAngles)

	@property
	def fractCoords(self):
		if ( (self._fractCoords is None) or (self._elementList is None) ):
			return None
		outList = list()
		for fCoords, element in it.zip_longest(self._fractCoords, self._elementList):
			fCoords.append(element)
			outList.append(fCoords)
		return outList

	@fractCoords.setter
	def fractCoords(self, value: "list of [x,y,z,Element]"):
		fCoords = list()
		eList = list()
		for x in value:
			fCoords.append(x[:3])
			eList.append(x[3])
		self._elementList = eList
		self._fractCoords = fCoords


	#Non-property Setter Functions
	def setLattParams(self, lattList:list):
		oldLattVects = self.getLattVects()
		self.lattParams = listToLattParams(lattList)
		newLattVects = self.getLattVects()
		if self.fractCoords is not None:
			self.fractCoords = getTransformedFractCoords(oldLattVects, newLattVects, self.fractCoords)

	def setLattAngles(self, lattAngles:list):
		oldLattVects = self.getLattVects()
		self.lattAngles = listToLattAngles(lattAngles)
		newLattVects = self.getLattVects()
		if self.fractCoords is not None:
			self.fractCoords = getTransformedFractCoords(oldLattVects, newLattVects, self.fractCoords)

	#Non-property Getter Functions
	def getLattParamsList(self):
		return [self.lattParams["a"], self.lattParams["b"], self.lattParams["c"]]

	def getLattAnglesList(self):
		return [self.lattAngles["alpha"], self.lattAngles["beta"], self.lattAngles["gamma"]]

	def getLattVects(self, **kwargs):
		lattVects = lattParamsAndAnglesToLattVects(self.getLattParamsList() , self.getLattAnglesList(), **kwargs )
		return lattVects

	def listToLattParams(self, lattList):
		if lattList is None:
			return None
		assert len(lattList) == 3, "lattParams must be a list 3 elements long; {} is invalid".format(lattList)
		lattParams = dict()
		lattParams["a"], lattParams["b"], lattParams["c"] = [float(x) for x in lattList]
		return lattParams

	def listToLattAngles(self, lattAngleList):
		if lattAngleList is None:
			return None
		assert len(lattAngleList) == 3, "lattAngleList must be a list 3 elements long; {} is invalid".format(lattAngleList)
		lattAnglesDict = dict()
		lattAnglesDict["alpha"], lattAnglesDict["beta"], lattAnglesDict["gamma"] = [float(x) for x in lattAngleList]
		return lattAnglesDict


	def getVolume(self):
		if (self.lattParams is not None) and (self.lattAngles is not None):
			return self.calcVolumeFromLattParamsAngles(self.lattParams, self.lattAngles)
		else:
			raise ValueError("UnitCell.getVolume() failed since either self.lattParams or self.lattAngles"
			                 "are undefined for current object")

	def calcVolumeFromLattParamsAngles(self, lattParams, lattAngles):
		abcFactor = lattParams["a"]*lattParams["b"]*lattParams["c"]
		angles = [math.radians(lattAngles["alpha"]), math.radians(lattAngles["beta"]), math.radians(lattAngles["gamma"])]
		alpha, beta, gamma = angles
		sqrtTerm = ( 1 + (2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)) 
		             - (math.cos(alpha)**2) - (math.cos(beta)**2) - (math.cos(gamma)**2) )
		return math.sqrt(sqrtTerm)*abcFactor

	def convAngToBohr(self):
		for key in self.lattParams.keys():
			self.lattParams[key] *= ANG_TO_BOHR


	def getCastepCellStr(self, units="bohr"):
		lattVects = self.getLattVects()
		lattStrVects = list()
		for currVect in lattVects:
			 lattStrVects.append(" ".join([str(x) for x in currVect]))

		lattStrings =  "\n".join(lattStrVects)
		outStr = units + '\n' + lattStrings

		return outStr

def lattVectorsFromCellFile(cellFilePath):
	fileAsDict = tokenizeCastepCellFile(cellFilePath)
	lattVectorStr = fileAsDict["%block lattice_cart"]
	vectList = [x for x in lattVectorStr.splitlines()]
	lattVects = [x.strip().split() for x in vectList if len(x.strip().split())==3]

	for vectorIdx,unused in enumerate(lattVects):
		lattVects[vectorIdx] = [ float(x) for x in lattVects[vectorIdx] ]

	assert len(lattVects) == 3, "Only found {} lattice vectors in {}".format(len(lattVects), cellFilePath)

	return lattVects


def lattParamsAndAnglesFromLattVects(lattVectors):
	sqrSums = [ sum([x**2 for x in currVect]) for currVect in lattVectors ]
	lattParams = [math.sqrt(sqrSum) for sqrSum in sqrSums]

	#Calc lattice angles:
	bcDot = sum( [x*y for x,y in zip(lattVectors[1], lattVectors[2])] )
	acDot = sum( [x*y for x,y in zip(lattVectors[0], lattVectors[2])] )
	abDot = sum( [x*y for x,y in zip(lattVectors[0], lattVectors[1])] )

	bc = lattParams[1]*lattParams[2]
	ac = lattParams[0]*lattParams[2]
	ab = lattParams[0]*lattParams[1]

	alpha = math.degrees( math.acos(bcDot/bc) )
	beta  = math.degrees( math.acos(acDot/ac) )
	gamma = math.degrees( math.acos(abDot/ab) )
	lattAngles = [alpha, beta, gamma]

	return lattParams, lattAngles


def lattParamsAndAnglesToLattVects(lattParams:"[a,b,c]", lattAngles:"[alpha,beta,gamma]", **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	testInverse = kwargs.get("testInverse", True)
	testTol = kwargs.get("testTol".lower(), (0.001, 0.01) ) 
	uVectA = [1.0, 0.0, 0.0]
	uVectB, uVectC = [0.0,0.0,0.0], [0.0, 0.0, 0.0] 

	uVectB[0] = math.cos(math.radians( lattAngles[2] ))
	uVectB[1] = math.sqrt( 1 - (uVectB[0]**2) ) #fix z component of B to zero
	uVectC[0] = math.cos(math.radians( lattAngles[1] ))
	uVectC[1] = ( math.cos(math.radians( lattAngles[0] )) - (uVectB[0]*uVectC[0]) ) / uVectB[1]
	uVectC[2] = math.sqrt (1 - ( (uVectC[0]**2) + (uVectC[1]**2) ))

	lattVects = list()
	lattVects.append([x*lattParams[0] for x in uVectA])
	lattVects.append([x*lattParams[1] for x in uVectB])
	lattVects.append([x*lattParams[2] for x in uVectC])

	#Test that these cell vectors give the input lattice params/angles
	if testInverse:
		invLattParams, invLattAngles = lattParamsAndAnglesFromLattVects(lattVects)
		errorMsg = ("Problem with lattParamsAndAnglesToLattVects\n Input Params/Angles = {}\t{}\nOutput Cell Vectors = {}\n,"
                     "Output Params/Angles = {}\t{}".format(lattParams, lattAngles, lattVects, invLattParams, invLattAngles)) 
		if not all( [abs(x-y)<testTol[0] for x,y in zip(lattParams,invLattParams)]):
			raise ValueError(errorMsg)
		elif not all( [abs(x-y)<testTol[1] for x,y in zip(lattAngles,invLattAngles)]):
			raise ValueError(errorMsg)

	return lattVects


