#!/usr/bin/env python3

import itertools as it
import json
import math

import numpy as np

from . import unit_convs as uConv

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
	finalFractCoords = [ x[:3] for x in getFractCoordsFromCartCoords(finalLattVects,finalCartArray) ]

	return finalFractCoords



class UnitCell():
	def __init__(self,**kwargs):
		""" Initializer
		
		Args:
			All args are optional keyword args; if none are passed the object is simply empty (and will throw various errors if you try to do certain things)
			lattParams: (len 3 float iter) [a,b,c] where a,b and c are the lattice parameters
			lattAngles: (len 3 float iter) [alpha,beta,gamma] where these are lattice angles as traditionally defined. Namely alpha,beta,gamma are the angles between bc,ac and ab respectivetly
			fractCoords: (iter of nx3 iter) Fractional coordinates for a list of atoms. e.g. [ [0.0,0.0,0.0], [0.5,0.5,0.5]]. Note that setting fractCoords OUTSIDE the initialiser also requires the atom symbols to be passed (e.g. [[0.0,0.0,0.0,"Mg"]]). Also note elementList SHOULD be set alongside fractCoords in the initialiser
			elementList: (str iter) Each entry contains a string representing the element of one atom in fractCoords. Therefore
			putCAlongZ: (Bool, default is False) If true then the 3rd lattice vector in self.lattVects will always be [0,0,c]. DEPRECATED/STUPID: PLEASE leave it as False			 

		"""
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self._lattParams = self.listToLattParams( kwargs.get( "lattParams".lower(), None ) )
		self._lattAngles = self.listToLattAngles( kwargs.get( "lattAngles".lower(), None ) )
		self._fractCoords = kwargs.get("fractCoords".lower(), None)
		self._elementList = kwargs.get("elementList".lower(), None)
		self.putCAlongZ = kwargs.get("putCAlongZ".lower(), False)
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


	def toDict(self):
		outDict = dict()
		outDict["lattParams"] = self.getLattParamsList()
		outDict["lattAngles"] = self.getLattAnglesList()
		outDict["fractCoords"] = self.fractCoords
		return outDict

	@classmethod
	def fromDict(self, inpDict):
		outObj = UnitCell(lattParams=inpDict["lattParams"], lattAngles=inpDict["lattAngles"])
		if inpDict.get("fractCoords",None) is not None:
			outObj.fractCoords = inpDict["fractCoords"]
		return outObj

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
		""" Deprecated; no idea if it works properly """
		with open(inpPath) as f:
			inpDict = json.load(f)
		outObj = cls.fromLattVects( inpDict["lattvects"], fractCoords = inpDict["fractCoords".lower()] )
		return outObj

	def writeFile(self,outPath):
		""" Deprecated; no idea if it works propert """
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
	def lattVects(self):
		return self.getLattVects(putCAlongZ=self.putCAlongZ)

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

		if (self.fractCoords is not None) and (self.fractCoords!=list()):
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
		lattVects = self.lattVects
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
	testInverse = kwargs.get("testInverse".lower(), True)
	putCAlongZ = kwargs.get("putCAlongZ".lower(), False)
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

	if putCAlongZ:
		lattVects = getLattVectorsTransformedToAlignParamCWithZ(lattVects)

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
	lVectArray = np.array(lattVects).transpose()
	fCoordArray = np.array(fractCoords).transpose()
	outArray = np.dot(lVectArray, fCoordArray)
	outList = outArray.transpose().tolist()
	return outList

def getFractCoordsFromCartCoords(lattVects, cartCoords):
	if len(cartCoords)==0:
		return list()

	outList = list()
	eleList = [x[-1] for x in cartCoords]
	coordList = [x[:3] for x in cartCoords]

	#Convert x,y,z into fract coords
	cartCoordsT = np.array(coordList).transpose()
	vMatrix = np.array(lattVects).transpose()
	fractVals = ( np.linalg.inv(vMatrix) @ np.array(cartCoordsT) ).transpose().tolist()

	#Construct the output list
	for fVals,eleKey in it.zip_longest(fractVals,eleList):
		outList.append( fVals+[eleKey] )

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

def foldAtomicPositionsIntoCell(cellObj, tolerance=1e-2):
	#This part is now the slow step (>50% of the runtime)
	startFractCoords = cellObj.fractCoords
	if startFractCoords is None:
		return None

	useFractCoords = np.array( cellObj._fractCoords ) #Avoid splitting up eles/coords

	#Transform
	foldFractCoordArrayToValsBetweenZeroAndOne(useFractCoords, tolerance=tolerance)

	cellObj._fractCoords = useFractCoords

	#The old, simple implementation. This is slower + less flexible than the newer numpy implementation
#	for idxA,currFract in enumerate(startFractCoords):
#		for idxB,x in enumerate(currFract[:3]):
#			if (x<-1*tolerance):
#				shiftVal = math.ceil(abs(x))
#				startFractCoords[idxA][idxB] += shiftVal
#			elif (x>1+tolerance):
#				shiftVal = -1* math.floor(x)
#				startFractCoords[idxA][idxB] += shiftVal

def foldFractCoordArrayToValsBetweenZeroAndOne(fractCoords, tolerance=1e-5):
	""" Function meant for optimisation purposees really; so foldAtomicPositionsIntoCell (and various inefficiencies associated with UnitCell class) can be bypassed
	
	Args:
		fractCoords: (nx3 numpy array) Fraction Coordinates
		tolerance: (float) We only apply a shift if the fraction coord x is -tolerance > x > 1+tolerance
			 
	Returns
		Nothing; works in place	
 
	"""
	#When a value is >1 we subtract math.floor(x); else leave it
	outArray = np.where(fractCoords>1+tolerance, fractCoords+(-1*np.floor(fractCoords)), fractCoords)

	#When a value is <1 we need to add math.ceil(x); else leave it
	outArray = np.where(outArray<-1*tolerance, outArray+np.ceil(np.abs(outArray)), outArray)

	np.copyto(fractCoords,outArray)



def applyTranslationVectorToFractionalCoords(inpCell, tVect, foldInAfter=False):
	""" Applies a translation vector (given in terms of fractional co-ordinates) to all atoms in the cell
	
	Args:
		inpCell: (UnitCell object) The cell we want to apply the translation to
		tVect: (len-3 iter) Translation (in fractional co-ordinates) for x,y,z
		foldInAfter: (Bool) If true fold all atomic positions back into cell after applying the transformation (so all output fractional co-ords are between 0 and 1)

	Returns
		Nothing. Acts in place.
 
	"""
#	startFractCoords = inpCell.fractCoords

	#Separate into np array(for fast translate) and indices (to map back)
#	atomSymbols = [x[-1] for x in startFractCoords]
#	useFractCoords = np.array([x[:3] for x in startFractCoords])
#	outCoordArray = np.add(useFractCoords,np.array(tVect)) #Will add to all elements
#
#	outFractCoords = outCoordArray.tolist()
#	for row, ele in it.zip_longest(outFractCoords, atomSymbols):
#		row.append(ele)



#	outFractCoords = list()
#	for fCoord in startFractCoords:
#		currCoord = [x+t for x,t in it.zip_longest(fCoord[:3],tVect)] + [fCoord[-1]]
#		outFractCoords.append(currCoord)

#	inpCell.fractCoords = outFractCoords	


	#OPTIMISED VERSION: We directly access _fractCoords and work with that, since access through the property is slow (and means we have to deal with the element symbols). This was about 4x as fast
	useFractCoords = np.array(inpCell._fractCoords)
	useFractCoords = np.add(useFractCoords, tVect)
	inpCell._fractCoords = useFractCoords.tolist()

	if foldInAfter:
		foldAtomicPositionsIntoCell(inpCell)

def applyTranslationVectorToCartCoords(inpCell, tVect, foldInAfter=False):
	""" Applies a translation vector (given in terms of cartesian co-ordinates) to all atoms in the cell
	
	Args:
		inpCell: (UnitCell object) The cell we want to apply the translation to
		tVect: (len-3 iter) Translation (in cartesian co-ordinates) for x,y,z
		foldInAfter: (Bool) If true fold all atomic positions back into cell after applying the transformation (so all output fractional co-ords are between 0 and 1)

	Returns
		Nothing. Acts in place.
 
	"""
	startCartCoords = inpCell.cartCoords
	outCartCoords = list()
	for cCoord in startCartCoords:
		currCoord = [x+t for x,t in it.zip_longest(cCoord[:3],tVect)] + [cCoord[-1]]
		outCartCoords.append( currCoord )
	
	inpCell.cartCoords = outCartCoords

	if foldInAfter:
		foldAtomicPositionsIntoCell(inpCell)

def getDensityFromUCellObj(uCellObj, massDict=None, lenConvFactor=1):
	""" Calculate the density in units of g/length**3 for a UnitCell object. If you want different mass units, then passing your own massDict is a way to do it
	
	Args:
		uCellObj: (UnitCell) Default length units will be those used here
		massDict: (dict) Keys are capitalized element keys ("mg".capitalize()="Mg") while values are floats in units of mass per mole (usually g/mol). Default is a dictionary of g/mol
		lenConvFactor: (float) Multiply the length units by this factor
 
	Returns
		 density: (float) Usually in units of g/(len**3) where len are the units in the UnitCell
 
	"""

	massDict = massDict if massDict is not None else getEleKeyToMassDictStandard()
	totalMass = 0
	volume = uCellObj.volume*(lenConvFactor**3)
	for cartCoord in uCellObj.cartCoords:
		eleKey = cartCoord[-1]
		totalMass += massDict[eleKey.capitalize()]/uConv.AVOGADRO_NUMBER

	outDensity = totalMass/volume
	return outDensity


def moveIndicesToTopOfGeomForUnitCell(inpCell, inpIndices):
	""" Re-arranges co-ordinates to put certain indices on top. Original use-case is putting all indices involving collective vars at the top of the geometry
	
	Args:
		inpCell: (UnitCell object)
		inpIndices: (iter of ints) Indices of atoms (base-zero) we want to move to the top of the geom coords
			 
	Returns
		Nothing; works in place
 
	"""
	fCoords = inpCell.fractCoords

	topCoords = [fCoords[idx] for idx in inpIndices]
	otherCoords = [fCoords[idx] for idx in range(len(fCoords)) if idx not in inpIndices]
	outFractCoords = topCoords + otherCoords
	inpCell.fractCoords = outFractCoords


#Taken from https://gist.github.com/lukasrichters14/c862644d4cbcf2d67252a484b7c6049c
def getEleKeyToMassDictStandard():
	elements_dict = {'H' : 1.008,'D' : 2.014,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
	             'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
	             'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
	             'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
	             'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
	             'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
	             'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
	             'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
	             'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
	             'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
	             'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
	             'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
	             'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
	             'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
	             'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
	             'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
	             'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
	             'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
	             'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
	             'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
	             'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
	             'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
	             'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
	             'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
	             'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
	             'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
	             'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
	             'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
	             'OG' : 294}
	outDict = dict()
	for key,val in elements_dict.items():
		outDict[key.capitalize()] = val
	return outDict


