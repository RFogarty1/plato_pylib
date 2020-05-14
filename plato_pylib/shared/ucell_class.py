#!/usr/bin/env python3

import itertools as it
import json
import math

import numpy as np


ANG_TO_BOHR = 1.8897259886




def _getTransformedFractCoordsWithElements(origLattVects:"list of lists [[v1],[v2],[v3]]", finalLattVects, origFractCoords:"list of lists"):
	tempFractCoords = list()
	for x in origFractCoords:
		tempFractCoords.append( list(x[:3]) )

	transCoords = getTransformedFractCoords(origLattVects, finalLattVects, tempFractCoords)

	finalFractCoords = list()
	for orig,trans in zip(origFractCoords,transCoords):
		currList = list(trans)
		currList.append( orig[-1] )
		finalFractCoords.append( currList )

	return finalFractCoords

def getTransformedFractCoords(origLattVects:"list of lists [[v1],[v2],[v3]]", finalLattVects, origFractCoords:"list of lists"):
	#I use transposes since i want to apply TC=C', where T is the transform matrix, C are orig. coords and C' are final coords
	lattVectsOrigT = np.array(origLattVects).transpose()
	lattVectsFinalT = np.array(finalLattVects).transpose()

	#Step 1 = get the matrix that transforms the lattice vectors
	transformMatrix = lattVectsFinalT @ np.linalg.inv(lattVectsOrigT)

	#Step 2 = Apply that transformation to the original cartesian co-ordinates
	origCartCoords = [ x for x in _getCartCoordsFromFract_NoElement(origLattVects, origFractCoords) ]
	origCartArrayT = np.array(origCartCoords).transpose()
	finalCartArray = (transformMatrix @ origCartArrayT).transpose()

	#Step 3 = Convert the transformed cartesian co-ordinates into fractional ones
	finalFractCoords = list()
	for currAtomCart in finalCartArray:
		finalFractCoords.append( _getFractCoordsFromCartOneAtom(finalLattVects, currAtomCart) )

	return finalFractCoords





class UnitCell():
	def __init__(self,**kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self._lattParams = self.listToLattParams( kwargs.get( "lattParams".lower(), None ) )
		self._lattAngles = self.listToLattAngles( kwargs.get( "lattAngles".lower(), None ) )
		self._fractCoords = kwargs.get("fractCoords".lower(), None)
		self._elementList = kwargs.get("elementList".lower(), None)
		self._eqTolPlaces = 5


	def __eq__(self,other):
		outVal = True
		if self._eqTolPlaces != other._eqTolPlaces:
			outVal = False

		allowedDiff = 10**(-1*self._eqTolPlaces)


		#Purely numberical arrays 
		relAttrs = ["_fractCoords"]
		for attr in relAttrs:
			attA,attB = getattr(self,attr), getattr(other,attr)
			if (attA is None) and (attB is None):
				pass
			elif (attA is None) or (attB is None):
				outVal = False
				break
			else:
				if len(attA) != len(attB):
					outVal = False
					break
				for x,y in it.zip_longest(attA,attB):
					if not all( [abs(a-b) < allowedDiff for a,b in it.zip_longest(x,y)] ):
						outVal = False
						break

		#Dictionaries of numeric vals
		relAttrs = ["lattParams","lattAngles"]
		for attr in relAttrs:
			attA,attB = getattr(self,attr), getattr(other,attr)
			if (attA is None) and (attB is None):
				pass
			elif (attA is None) or (attB is None):
				outVal = False
				break
			else:
				keysA, keysB = sorted(attA.keys()), sorted(attB.keys())
				for keyA,keyB in it.zip_longest(keysA,keysB):
					if keyA!=keyB:
						outVal = False
						break
				else:
					for key in keysA:
						if abs(attA[key]-attB[key]) > allowedDiff:
							outVal = False
							break

		#Other
		if (self._elementList is None) and (other._elementList is None):
			pass
		elif (self._elementList is None) or (other._elementList is None):
			outVal = False
		else:
			if self._elementList != other._elementList:
				outVal = False

		return outVal


	@classmethod
	def fromLattVects(cls, lattVectors:"iterable, len=3", fractCoords = None):
		lattParams, lattAngles = lattParamsAndAnglesFromLattVects(lattVectors)
		outObj = cls(lattParams=lattParams, lattAngles = lattAngles)
		if fractCoords is not None:
			finalLattVects = outObj.lattVects
			transFractCoords = _getTransformedFractCoordsWithElements(lattVectors, finalLattVects, fractCoords)
			outObj.fractCoords = transFractCoords
		return outObj

	@classmethod
	def fromFile(cls, inpPath):
		with open(inpPath) as f:
			inpDict = json.load(f)
		outObj = cls.fromLattVects( inpDict["lattvects"], fractCoords = inpDict["fractCoords".lower()] )
		return outObj

	def writeFile(self,outPath):
		outDict = dict()
		outDict["lattvects"] = self.lattVects
		outDict["fractCoords".lower()] = self.fractCoords
		with open(outPath,"wt") as f:
			f.write(json.dumps(outDict))

	@property
	def fractCoords(self):
		if ( (self._fractCoords is None) or (self._elementList is None) ):
			return None
		outList = list()
		for fCoords, element in it.zip_longest(self._fractCoords, self._elementList):
			tempList = list(fCoords)
			tempList.append(element)
			outList.append(tempList)
		return outList

	@fractCoords.setter
	def fractCoords(self, value: "list of [x,y,z,Element]"):
		fCoords = list()
		eList = list()
		for x in list(value):
			fCoords.append(list(x[:3]))
			eList.append(x[3])
		self._elementList = eList
		self._fractCoords = fCoords


	@property
	def cartCoords(self):
		return self._getCartCoords(sort=False)

	@cartCoords.setter
	def cartCoords(self, value: "list of [x,y,z,Element]"):
		fractCoords = getFractCoordsFromCartCoords(self.lattVects,value)
		self.fractCoords = fractCoords

	@property
	def volume(self):
		return self.getVolume()

	@volume.setter
	def volume(self, value:float):
		angularTerm = self._calcVolumeAngularTerm()
		currVolOverAngular = self.volume / angularTerm
		newVolOverAngular = value / angularTerm
		scaleFactor = (newVolOverAngular / currVolOverAngular)**(1/3)
		for key in self.lattParams.keys():
			self.lattParams[key] *= scaleFactor

	#Note that a getter is still exposed for use of keywords (e.g. disabling error checks)
	@property
	def lattVects(self, **kwargs):
		return self.getLattVects()

	@lattVects.setter
	def lattVects(self, value):
		lattParams, lattAngles = lattParamsAndAnglesFromLattVects(value)
		self.setLattParams(lattParams)
		self.setLattAngles(lattAngles)


	@property
	def lattParams(self):
		return self._lattParams

	@lattParams.setter
	def lattParams(self,value):
		if not isinstance(value,dict):
			raise ValueError("Can only set lattParams using a dict")
		oldLattVects = self.lattVects
		self._lattParams = dict(value)
		newLattVects = self.lattVects

		if self.fractCoords is not None:
			self._fractCoords = getTransformedFractCoords(oldLattVects, newLattVects, self._fractCoords)

	@property
	def lattAngles(self):
		return self._lattAngles


	@lattAngles.setter
	def lattAngles(self,value):
		if not isinstance(value,dict):
			raise ValueError("Can only set lattAngles using a dict")
		oldLattVects = self.lattVects
		self._lattAngles = dict(value)
		newLattVects = self.lattVects

		if self.fractCoords is not None:
			self._fractCoords = getTransformedFractCoords(oldLattVects, newLattVects, self._fractCoords)


	#Non-property Setter Functions
	def setLattParams(self, lattList:list):
		''' Set lattParams from a list [a,b,c]. Maintained for backwards compatability '''
		self.lattParams = self.listToLattParams(lattList)


	def setLattAngles(self, lattAngles:list):
		self.lattAngles = self.listToLattAngles(lattAngles)


	def _getCartCoords(self,sort=False):
		cartCoords = getCartCoordsFromFractCoords(self.lattVects, self.fractCoords)
		if sort is False:
			return cartCoords

		#Convert the numbers to an np array, sort that and convert back to list
		cCoords = np.array([x[:3] for x in cartCoords])
		sortOrder = np.lexsort( (cCoords[:,2],cCoords[:,1],cCoords[:,0]) )
		outputCoords = cCoords[sortOrder].tolist()
		outputElementList = [cartCoords[i][3] for i in sortOrder]
		
		outputCartCoords = list()
		for coords,element in it.zip_longest(outputCoords,outputElementList):
			currList = list(coords)
			currList.append(element)
			outputCartCoords.append(currList)

		return outputCartCoords


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
		angularTerm = self._calcVolumeAngularTerm()
		return angularTerm*abcFactor

	def _calcVolumeAngularTerm(self):
		angles = [math.radians(self.lattAngles["alpha"]), math.radians(self.lattAngles["beta"]), math.radians(self.lattAngles["gamma"])]
		alpha, beta, gamma = angles
		sqrtTerm = ( 1 + (2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)) 
		             - (math.cos(alpha)**2) - (math.cos(beta)**2) - (math.cos(gamma)**2) )
		return math.sqrt(sqrtTerm)

	def convAngToBohr(self):
		for key in self.lattParams.keys():
			self.lattParams[key] *= ANG_TO_BOHR

	def convBohrToAng(self):
		for key in self.lattParams.keys():
			self.lattParams[key] *= (1/ANG_TO_BOHR)


	#TODO: Eventually this needs to be removed. It can be accesed through a function in parseCastep at current
	def getCastepCellStr(self, units="bohr"):
		lattVects = self.getLattVects()
		lattStrVects = list()
		fmt = "{:.7g}"
		zeroTol = 1e-8 #If a component is below this, then we treat it as zero

		for currVect in lattVects:
			for idx,unused in enumerate(currVect):
				currVect[idx] = currVect[idx] if abs(currVect[idx])>zeroTol else 0

		for currVect in lattVects:
			 lattStrVects.append(" ".join([ fmt.format(x) for x in currVect]))

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

def getCartCoordsFromFractCoords(lattVects, fractCoords):
	coordsOnly = [x[:3] for x in fractCoords]
	cartCoords = _getCartCoordsFromFract_NoElement(lattVects, coordsOnly)
	for cart,fract in it.zip_longest(cartCoords,fractCoords):
		cart.append( fract[-1] )
	return cartCoords

def _getCartCoordsFromFract_NoElement(lattVects,fractCoords):
	outList = list()
	for currAtom in fractCoords:
		fromA  = [currAtom[0]*x for x in lattVects[0]]
		fromB  = [currAtom[1]*x for x in lattVects[1]]
		fromC  = [currAtom[2]*x for x in lattVects[2]]
		outVect = [a+b+c for a,b,c in zip(fromA,fromB,fromC)]
		outList.append(outVect)
	return outList

def getFractCoordsFromCartCoords(lattVects, cartCoords):
	outList = list()
	for currAtom in cartCoords:
		outList.append( _getFractCoordsFromCartOneAtom(lattVects, currAtom[:3]) + [currAtom[-1]] ) 
	return outList

def _getFractCoordsFromCartOneAtom(lattVects:"3x3 iterxiter", cartCoords:"3 item list"):
	vMatrix = np.array(lattVects).transpose()

	fractVals = ( np.linalg.inv(vMatrix) @ np.array(cartCoords) ).tolist()
	return fractVals



def getLattVectorsTransformedToAlignParamCWithZ(lattVectors):
	lattParams, lattAngles = lattParamsAndAnglesFromLattVects(lattVectors)
	lattParamA, lattParamB, lattParamC = lattParams
	alpha, beta, gamma = lattAngles

	#Each row is one lattice vector
	outVectors = np.zeros((3,3))
	outVectors[2,2] = lattParamC
	outVectors[0,2] = lattParamA*math.cos(math.radians(beta ))
	outVectors[1,2] = lattParamB*math.cos(math.radians(alpha))
	outVectors[0,0] = math.sqrt( (lattParamA**2) - (outVectors[0,2]**2) ) 
	outVectors[1,0] = ( (lattParamA*lattParamB*math.cos(math.radians(gamma))) - (outVectors[0,2]*outVectors[1,2]) )  / outVectors[0,0] 
	outVectors[1,1] = math.sqrt( (lattParamB**2) - (outVectors[1,0]**2) - (outVectors[1,2]**2) )


	#Convert to output
	output = list()
	for currRow in outVectors:
		output.append( list(currRow) )

	return output

