#!/usr/bin/python3

import math
import numpy as np
import itertools

ANG_TO_BOHR = 1.8897259886


def tokenizeCastepCellFile(cellFilePath:str)->"lower case dict":
	with open(cellFilePath, "rt") as f:
		fileAsList = f.readlines()
	removeCommentsFromFileAsList(fileAsList)
	fileAsList = [x for x in fileAsList if x.strip()!=''] #remove blank lines

	lineIdx = 0
	outputDict = dict()
	while lineIdx < len(fileAsList):
		line = fileAsList[lineIdx]
		if line.lower().find("%block") != -1:
			currArg = " ".join( line.strip().lower().split()[:2] )
			lineIdx += 1
			outputDict[currArg], lineIdx = parseCellBlock(fileAsList, lineIdx)
		else:
			currKey = line.strip().lower().split()[0]
			currVal = " ".join( line.strip().lower().split()[1:] )
			outputDict[currKey] = currVal
			lineIdx += 1

	return outputDict


def removeCommentsFromFileAsList(fileAsList:"list; each entry=1 line of a file"):

	for lineIdx,unused in enumerate(fileAsList):
		commentPosStart = fileAsList[lineIdx].find('#')
		if commentPosStart != -1:
			fileAsList[lineIdx] = fileAsList[lineIdx][:commentPosStart] + "\n"

def parseCellBlock(fileAsList, lineIdx):
	parsedStr = ""
	while lineIdx < len(fileAsList):
		line = fileAsList[lineIdx]
		if line.lower().find("endblock") == -1:
			parsedStr += fileAsList[lineIdx]
		else:
			lineIdx += 1
			break
		lineIdx += 1
	return parsedStr.strip(), lineIdx


def modCastepCellFile(inpFilePath:str, fieldValDict: "dict {keyword:value}"):
	fieldValDict = {k.lower():v for k,v in fieldValDict.items()}
	tokenizedOutFile = tokenizeCastepCellFile(inpFilePath)

	for key,val in fieldValDict.items():
		outKey = key
		testStrB, testStrC = "%block "+key, key.replace("%block","")
		if (testStrB in tokenizedOutFile.keys()):
			outKey = testStrB
		elif (testStrC in tokenizedOutFile.keys()):
			outKey = testStrC
		tokenizedOutFile[outKey] = val

	writeCastepCellFileFromTokens(inpFilePath, tokenizedOutFile)

def writeCastepCellFileFromTokens(outFilePath, tokens: "dict {keyword:value}"):
	fileStr = "\n"
	for key, val in tokens.items():
		endVal = "\n\n"
		if val.strip().count('\n') > 0 or key.find("%block") != -1:
			key = "%block " + key.replace("%block ","") + '\n'
			endVal = "\n%endblock " + key.replace("%block ","") + "\n\n"
		else: 
			key += ' '
		fileStr +=  key + val + endVal

	with open(outFilePath,"wt") as f:
		f.write(fileStr)	


###########################################################################
# Function takes a *.castep file from castep and extracts k points, 
# eigenvalues and the fermi energy in dictionary format
# Output keys:
#	'energy' : The final energy in the file without basis set corr. or 
#		   free energy correction (just grabs line with "Final energy"
#		   str present
#	'numbAtoms': Number of atoms in the primitive cell 
###########################################################################
def parseCastepOutfile(inpFile: str) -> dict:
	
	# Read in file
	with open(inpFile,'rt') as f:
		fileList = f.readlines()

	# Loop over looking for info we want
	for lineIdx, line in enumerate(fileList):
		if line.find('Total number of ions in cell') != -1:
			numbAtoms = int(line.strip().split()[-1])
		elif line.find('Final energy') != -1:
			energy = float(line.strip().split()[-2])
		elif line.find("Unit Cell") != -1:
			unitCellObj, unusedIdx = parseCastepUnitCellSection(fileList, lineIdx)
		elif line.find("MP grid size for SCF calculation") != -1:
			scf_k_grid = [ int(x) for x in line.split()[-3:] ]
		elif line.find("Number of kpoints used") != -1:
			scf_numb_k = int( line.split()[-1] )
		elif line.find("plane wave basis set cut-off") != -1:
			kin_cut = float( line.split()[-2] )

	# Create the dictionary and return it
	outDict = { 'numbAtoms' : numbAtoms,
		        'energy' : energy,
	            'unitCell': unitCellObj,
	            'scf_numb_k': scf_numb_k,
 	            'scf_kgrid': scf_k_grid,
	            'kin_cut' : kin_cut}
	return outDict



def parseCastepUnitCellSection(fileList:list, startPos:int):
	lineIdx = startPos
	
	lattParams = list()
	lattAngles = list()
	fractCoords = list()
	lattVectsFromFile = list()
	atomsInOrder = list()
	while fileList[lineIdx].find("Details") == -1:
		if fileList[lineIdx].find("SCF loop") != -1:
			break
		elif fileList[lineIdx].find("BFGS") != -1:
			break

		if fileList[lineIdx].find("Real Lattice") != -1:
			lineIdx += 1
			for currLineIdx in range(lineIdx, lineIdx+3):
				currLine = fileList[currLineIdx].strip().split()
				lattVectsFromFile.append( [float(x) for x in currLine[:3]] )

		if fileList[lineIdx].find("Lattice parameters") != -1:
			lineIdx += 1
			for currLineIdx in range(lineIdx, lineIdx+3):
				currLine = fileList[currLineIdx].strip().split()
				lattParams.append( currLine[2]  )
				lattAngles.append( currLine[-1] )

		if fileList[lineIdx].find("Fractional coordinates") != -1:
			lineIdx += 3
			while fileList[lineIdx].find("xx") == -1:
				currLine = fileList[lineIdx].strip().split()
				fractCoords.append([ float(currLine[3]), float(currLine[4]), float(currLine[5])])
				atomsInOrder.append( currLine[1] )
				lineIdx += 1

		else:
			lineIdx+=1

	#Transform fract co-ords the same way as the input cell vectors (output cell vectors are directly 
	#from the angles and lengths
	unitCellObj = UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	lattVectsFromCell = unitCellObj.getLattVects()
	transFractCoords = getTransformedFractCoords(lattVectsFromFile, lattVectsFromCell, fractCoords)
	for idx, unused in enumerate(transFractCoords):
		transFractCoords[idx].append( atomsInOrder[idx] )
	unitCellObj.fractCoords = transFractCoords

	return unitCellObj, lineIdx


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

	@classmethod
	def fromCellFile(cls, cellFilePath):
		lattVects = lattVectorsFromCellFile(cellFilePath)
		return cls.fromLattVects(lattVects)

	@property
	def fractCoords(self):
		if ( (self._fractCoords is None) or (self._elementList is None) ):
			return None
		outList = list()
		for fCoords, element in itertools.zip_longest(self._fractCoords, self._elementList):
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


