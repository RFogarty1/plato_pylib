
import math
import itertools as it
import re
from ..plato import parse_gau_files as parseGau

def parseCP2KBasisFile(inpPath):
	fileAsList = _readInpFileIntoIter(inpPath)
	allBasisFileAsList = _getAllBasisSetSectionsFromCP2KCleanedFileAsList(fileAsList)
	allBasisSets = list()
	for x in allBasisFileAsList:
		allBasisSets.append( _parseSingleBasisFileAsList(x) )

	return ParsedBasisFileCP2K(inpPath,allBasisSets)


#Here purely so i can mock it
def _readInpFileIntoIter(inpPath):
	with open(inpPath,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList

#Currently a super inefficient implementation. File needs traversing at least twice and the whole file always needs to be read, even if i only want one basis thats at the top. Generator expression with sensible terminating conditions would be best way, but annoyingly complicated
def _getAllBasisSetSectionsFromCP2KCleanedFileAsList(fileAsList):
	""" Takes CP2K basis file in special list format(see below) and returns an iter of the subsections for each basis set present in the file
	
	Args:
		fileAsList: str iter, each element contains one line of the basis input file
			 
	Returns
		basisSects: iter of str iters, each element contains the lines corresponding to one basis set
 
	"""
	cleanedFileAsList = [x for x in fileAsList if ( (not x.startswith("#")) and (x.strip()!="") )]

	anyLetterRegExp = re.compile("[A-Z,a-z]+") #Letters only appear in the first line of the basis definition

	allBasisStartLineIndices = [0] #At this point we have no separations between the basis sets
	for idx,line in enumerate(cleanedFileAsList[1:],1):
		if _regExpInStr(anyLetterRegExp,line):
 			allBasisStartLineIndices.append(idx)

	allLists = list()
	for idx,x in enumerate(allBasisStartLineIndices):
		startIdx = allBasisStartLineIndices[idx]
		if idx+1 < len(allBasisStartLineIndices):
			endIdx = allBasisStartLineIndices[idx+1]
		else:
			endIdx = len(cleanedFileAsList)
		allLists.append( cleanedFileAsList[startIdx:endIdx] )

	return allLists



def _regExpInStr(regExpr,inpStr):
	if len( re.findall(regExpr,inpStr) ) == 0:
		return False
	else:
		return True


def _parseSingleBasisFileAsList(fileAsList):
	""" Takes a list containing lines defining a single CP2K basis set and parses it into a convenient object
	
	Args:
		fileAsList: Each element in a line in the CP2K basis set format
			 
	Returns
		 basisSetCP2K: (basisSetCP2K  object) Full representation of the CP2K basis set

	"""
	element, *basisName = fileAsList[0].strip().split()
	basisExpansions = list()
	nSets = int(fileAsList[1].strip())

	startLine = 2
	for idx in range(nSets):
		currSet,nLines =  _parseNextBasisExponentSet(fileAsList[startLine:])
		startLine += nLines
		basisExpansions.append(currSet)

	return BasisSetCP2K(element,basisName,basisExpansions)

#Needs to start on the header line; doesnt matter if list contains MORE than needed
def _parseNextBasisExponentSet(fileAsList):
	nVal, minL, maxL, nExponents, *nEachLVal = [int(x) for x in fileAsList[0].strip().split()]
	nLines  = nExponents + 1

	#Figure out the l-values corresponding to each coefficient
	lVals = list()
	for lVal,nVals in enumerate(nEachLVal,minL):
		currLVals = [lVal for x in range(nVals)]
		lVals.extend( currLVals )

	#Figure out exponents and coefficients
	exponents, coeffs = list(), list()
	
	for idx,line in enumerate(fileAsList[1:nExponents+1]):
		exponents.append( float(line.strip().split()[0]) )
		coeffs.append( [float(x) for x in line.strip().split()[1:]] )

	#Convert coeffs such that its one set per basis function rather than one per exponent
	outCoeffs = list()
	numbBasisFuncts = len(coeffs[0])
	for idx in range(numbBasisFuncts):
		currCoeffs = [x[idx] for x in coeffs]
		outCoeffs.append(currCoeffs)

	return ExponentSetCP2K(exponents,outCoeffs,lVals,nVal), nLines


class ParsedBasisFileCP2K():
	""" Represents a parsed basis set file for CP2K
	"""

	def __init__(self,inpPath, basisSets):
		""" Initializer
		
		Args:
			inpPath: (str) Full path to the input file
			basisSets: (iter of BasisSetCP2K) These each represent one full basis set. The iter should be ordered in the same way they appear in the basis-set file
				 
		"""
		self.inpPath = inpPath
		self.basisSets = list(basisSets)


	def getUniqueBasisSet(self, element, basisName):
		""" Gets basis set corresponding to input element and basis-set name
		
		Args:
			element: (str, case insensitive) The element symbol for the basis set required (e.g. "Mg")
			basisName: (str, case insensitive) The name of the basis set

		Returns
			outBasis: (BasisSetCP2K object) Representation of the basis set
 
		Raises:
			AssertionError: If more than one basis set matches the crietria. This shouldnt really ever happen
			KeyError: If the basis set is not present
		"""

		searchEle, searchBasis = element.lower(), basisName.lower()

		outBasisSets = list()
		for bSet in self.basisSets:
			currEle, currBasisNames = bSet.element.lower(), [x.lower() for x in bSet.basisNames]
			if (searchEle==currEle) and (searchBasis in currBasisNames):
				outBasisSets.append( bSet )

		assert len(outBasisSets) < 2, "Searched basis sets should be unique, but {} found for ele={}, name={}".format( len(outBasisSets), element, basisName)
		if len(outBasisSets) == 0:
			raise KeyError("Basis set {} {} not found".format(element, basisName))

		return outBasisSets[0]



	def __eq__(self,other):
		attrs = ["inpPath","basisSets"]
		for currAttr in attrs:
			if getattr(self,currAttr) != getattr(other,currAttr):
				return False
		return True



class BasisSetCP2K():
	"""Class representing a single basis set in a format convenient for CP2K

	"""
	def __init__(self, element, basisNames, exponentSets):
		""" Initializer
		
		Args:
			element: (str), Elemental symbol (e.g. Mg, H, He)
			basisNames: (str list) The names given to the basis set in the CP2K file (e.g. DZV-GTH-q1). One basis set can have more than one name, hence why its a list
			exponentSets: (iter of ExponentSetCP2K objects) Each represents all basis functions associated with one set of exponents

		"""
		self.element = element
		self.basisNames = list(basisNames)
		self.exponentSets = exponentSets

	def __eq__(self, other):
		attrs = ["element", "exponentSets"]
		for currAttr in attrs:
			if getattr(self, currAttr) != getattr(other, currAttr):
				return False

		if sorted(self.basisNames) != sorted(other.basisNames):
			return False

		return True


class ExponentSetCP2K():
	"""Class representing all basis functions which share a set of exponents
	"""

	def __init__(self, exponents, coeffs, lVals, nVal):
		""" Initializer
		
		Args:
			exponents: (iter of floats) The list of exponents shared between these basis functions
			coeffs: (iter of iters) Each element is the list of coefficients for one basis function
			lVals: (iter of ints) Each element is the orbital angular momentum (i.e. l quantum number; s=0, p=1, d=2 etc.)

		"""
		self._eqTol = 1e-7
		self.nVal = nVal
		self.exponents = exponents
		self.coeffs = coeffs
		self.lVals = lVals
		self._checkAllCoeffListsSameLengthAsExponents()

	def _checkAllCoeffListsSameLengthAsExponents(self):
		assert all([len(x)==len(self.exponents) for x in self.coeffs]), "All coeff lengths should be {} (the number of exponents)".format(len(self.exponents))

	def __eq__(self,other):
		tolerance = min(self._eqTol,other._eqTol)

		if self.nVal != other.nVal:
			return False

		#Check exponents
		if not self._listsSameWithinErrorTol(self.exponents,other.exponents, tolerance):
			return False

		#Check coeffs and lVals in the same manner ish
		if len(self.coeffs) != len(other.coeffs):
			return False

		for coeffsA, coeffsB in zip(self.coeffs,other.coeffs):
			if not self._listsSameWithinErrorTol(coeffsA,coeffsB,tolerance):
				return False

		if not self._listsSameWithinErrorTol(self.lVals,other.lVals,tolerance):
			return False

		return True


	def _listsSameWithinErrorTol(self, listA, listB, tolerance):
		if (len(listA) != len(listB)):
			return False

		diffVals = [abs(x-y) for x,y in zip(listA,listB)]
		if any([x>tolerance for x in diffVals]):
			return False

		return True


def writeBasisFileFromParsedBasisFileObj(filePath, parsedBasisFile):
	""" Writes a file with all basis sets in parsedBasisFile object
	
	Args:
		filePath: (str) Path to the file you want to write to
		parsedBasisFile: (ParsedBasisFileCP2K) Contains information on all the basis sets to write to file
 
	"""
	fileAsList = _getFileAsListFromParsedBasisFile(parsedBasisFile)
	_writeFileAsListToPath(filePath, fileAsList)

def _writeFileAsListToPath(filePath, fileAsList):
	with open(filePath,"wt") as f:
		f.write("\n".join(fileAsList))

def _getFileAsListFromParsedBasisFile(parsedBasisFile):
	""" Gets a list of strings, 1 per line, representing the output file for the basis sets in input argument 
	
	Args:
		parsedBasisFile: (ParsedBasisFileCP2K object) contains all info from a CP2K basis file
			 
	Returns
		 fileAsList: (iter of strs) Each element is 1 line of the output file in the format used in CP2K
 
	"""
	basisSetsAsLists = list()
	for x in parsedBasisFile.basisSets:
		basisSetsAsLists.extend( _getFileAsListForCP2KBasis(x) )
		basisSetsAsLists.append('')

	return basisSetsAsLists

def _getFileAsListForCP2KBasis(inpCP2KBasis):
	""" Gets the string (used in CP2K file) for a single CP2K basis set (in the format of a list, see below)
	
	Args:
		inpCP2KBasis: (BasisSetCP2K object)
			 
	Returns
		 outStr: (iter of strs) String used to represent this basis set in a CP2K file, each element is one line
 
	"""
	basisAsList = list()
	firstLine = inpCP2KBasis.element + " " + " ".join([x for x in inpCP2KBasis.basisNames])
	basisAsList.append(firstLine)
	basisAsList.append( " {}".format(len(inpCP2KBasis.exponentSets)) )

	for x in inpCP2KBasis.exponentSets:
		basisAsList.extend( _getFileAsListForExponentSet(x) )


	return basisAsList

def _getFileAsListForExponentSet(inpExponentSet):
	expSet = _getExponentSetWithAngMomValsInCorrectOrder(inpExponentSet)
	outList = list()

	#Figure out the first line
	nVal, minL, maxL, numbExp = inpExponentSet.nVal, min(inpExponentSet.lVals), max(inpExponentSet.lVals), len(inpExponentSet.exponents)

	numberL = list()
	for currL in range(minL,maxL+1):
		currNumbL = len([x for x in inpExponentSet.lVals if x==currL])
		numberL.append(currNumbL)

	firstLine = " ".join( [str(x) for x in [nVal, minL, maxL, numbExp] + numberL] )
	outList.append(firstLine)

	lineFmt = "\t" + " ".join(["{:.8f}" for x in range(len(inpExponentSet.coeffs)+1) ])
	for idx,exp in enumerate(inpExponentSet.exponents):
		currCoeffs = [x[idx] for x in inpExponentSet.coeffs]
		outLine = lineFmt.format( *([exp] + currCoeffs) )
		outList.append(outLine)

	return outList


def _getExponentSetWithAngMomValsInCorrectOrder(exponentSet):
	uniqueLVals = sorted( list(set([x for x in exponentSet.lVals])) )
	outCoeffs, outLVals = list(), list()
	for uniqueL in uniqueLVals:
		for l,coeffs in zip(exponentSet.lVals, exponentSet.coeffs):
			if l==uniqueL:
				outLVals.append(l)
				outCoeffs.append(coeffs)	

	outExponents = ExponentSetCP2K(exponentSet.exponents, outCoeffs, outLVals, exponentSet.nVal)

	return outExponents

def getCP2KBasisFromPlatoOrbitalGauPolyBasisExpansion(gauPolyBasisObjs, angMomVals, eleName, basisNames=None, shareExp=True, nVals=None):
	""" Gets a BasisSetCP2K object, with coefficients normalised, from an iter of GauPolyBasis objects in plato format
	
	Args:
		gauPolyBasisObjs: (iter plato_pylib GauPolyBasis object) Each element is the Plato representation of a basis function
		angMomVals: (iter of int) Angular momentum values for each orbital
		eleName: (str) Label for the element used in this basis set
		basisNames: (iter of str) Names used to specify this basis set in the CP2K input file (more than one allowed)
		shareExp: (Bool, Optional) If True, will try to exploit shared exponents when generating the basis set

	Returns
		 outBasis: (BasisSetCP2K Object) Convenient representation of a basis set for CP2K; this is the object that would be parsed from a CP2K basis file
		 
	"""
	if basisNames is None:
		basisNames = ["basis_set_a"]

	if nVals is None:
		nVals = [1 for x in gauPolyBasisObjs]

	splitExponentSets = [getCP2KExponentSetFromGauPolyBasis(gaus, angMom, nVal) for gaus,angMom,nVal in it.zip_longest(gauPolyBasisObjs, angMomVals, nVals)]

	if shareExp:
		outExponentSets = _getExponentSetsWithSharedExponentPartsMerged(splitExponentSets)
	else:
		outExponentSets = splitExponentSets
	
	outObj = BasisSetCP2K(eleName, basisNames, outExponentSets)
	return outObj


def _getExponentSetsWithSharedExponentPartsMerged(exponentSets, expTol=1e-4):
	""" Takes an iter of ExponentSetCP2K objects and returns iter (shorter or same length) where any objects with shared exponents are merged
	
	Args:
		exponentSets: (iter of ExponentSetCP2K objects) Effectively contain (at least) a single basis function each
			 
	Returns
		outExponentSets: (iter of ExponentSetCP2K objects) Same as input except each should now have a unique set of exponents (shared exponents being merged into the same object)
 
	"""
	
	mergedIndices = list()
	newOutput = list()
	for idxA, expSet in enumerate(exponentSets):
		if idxA in mergedIndices:
			pass
		elif idxA-1==(len(exponentSets)): #Only need to append if we havnt already merged
			newOutput.append( exponentSets[idxA] )
		else:
			expToMatch = expSet.exponents
			for idxB in range(idxA+1, len(exponentSets)):
				currExpSet = exponentSets[idxB]
				if len(expToMatch) == len(currExpSet.exponents):
					diffs = [abs(x-y) for x,y in zip(expToMatch, currExpSet.exponents)]
					if all([x < expTol for x in diffs]):
						_appendExponentSetBToA( exponentSets[idxA], currExpSet) 
						mergedIndices.append(idxB)
			newOutput.append(exponentSets[idxA])

	return newOutput

def _appendExponentSetBToA( expSetA, expSetB ):
	expSetA.coeffs.extend(expSetB.coeffs)
	expSetA.lVals.extend( expSetB.lVals )


def getCP2KExponentSetFromGauPolyBasis(gauPolyBasis, angMom, nVal=1):
	""" Returns ExponentSetCP2K object, with normalisation applied from gauPolyBasis. Meant as part of transforming basis sets from plato to cp2k format
	
	Args:
		gauPolyBasis: (plato_pylib GauPolyBasis object) Plato representation of a basis set
		angMom: (integer) 0,1,2 for s,p,d; Needed to carry out the normalisation
			 
	Returns
		outExponentSet: (ExponentSetCP2K Object) Representation of this basis function in a useful format for CP2K. Note this includes the normalisation of coefficients
 
	Raises:
		 AssertionError: If gauPolyBasis has more than r^0 terms in it
	"""
	assert gauPolyBasis.nPoly==0, "nPoly=0 required, but nPoly={} found".format(gauPolyBasis.nPoly)
	exponents = [x for x in gauPolyBasis.exponents]
	normFactors = [1/calcNormConstantForCP2KOnePrimitive(x,angMom) for x in exponents]
	coeffs = [x*normFactor for x,normFactor in it.zip_longest(gauPolyBasis.r0Coeffs,normFactors)]

#	outCoeffs = [[x] for x in coeffs]
	outObj = ExponentSetCP2K(exponents, [coeffs], [angMom], nVal)
	return outObj

#Need to divide by this when going Plato->CP2K
def calcNormConstantForCP2KOnePrimitive(exponent:float, angMom:int):
	expZet = 0.25 * ((2*angMom)+3)
	preFac = (2**angMom) * ((2/math.pi)**0.75) 
	normFactor = preFac*(exponent**expZet)
	return normFactor


def getGauPolyBasisFunctionsFromCP2KBasisSet(cp2kBasisSet):
	""" Gets a list of GauPolyBasis format objects from a CP2K Basis set (includes reversing the normalisation needed for CP2K basis sets
	
	Args:
		cp2kBasisSet: (BasisSetCP2K Object) Represents a full basis set in CP2K format
			 
	Returns
		gauPolyExpansions: (iter of GauPolyBasis objects) Each represents ONE basis function in the standard format for plato
 
	"""
	outObjs = list()
	for x in cp2kBasisSet.exponentSets:
		currObjs = getGauPolyBasisFunctionsFromCP2KExponentSet(x)
		outObjs.extend(currObjs)
	return outObjs

def getGauPolyBasisFunctionsFromCP2KExponentSet(cp2kExponentSet):
	""" Gets basis set in the GauPolyBasis format from CP2K format (including reversing the normalisation needed for a CP2K basis
	
	Args:
		cp2kExponentSet: (ExponentSetCP2K object) Representation of basis functions sharing exponents in CP2K
			 
	Returns
		gauPolyExpansions: (iter of GauPolyBasis objects) Each represents ONE basis function in the standard format for plato
 
	"""
	exponents = cp2kExponentSet.exponents

	outObjs = list()

#	for idx, angMom in enumerate(cp2kExponentSet.lVals):
#		coeffs = [x[idx] for x in cp2kExponentSet.coeffs] #For the CURRENT orbital
#		currNormFactors = [calcNormConstantForCP2KOnePrimitive(exp,angMom) for exp in exponents]
#		currCoeffs = [c*normFactor for c,normFactor in it.zip_longest(coeffs,currNormFactors)]
#		currObj = parseGau.GauPolyBasis(exponents, [currCoeffs], label=angMom)
#		outObjs.append(currObj)
#	return outObjs

	for coeffs,angMom in it.zip_longest(cp2kExponentSet.coeffs, cp2kExponentSet.lVals):
		currNormFactors = [calcNormConstantForCP2KOnePrimitive(exp,angMom) for exp in exponents]
		currCoeffs = [c*normFactor for c,normFactor in it.zip_longest(coeffs,currNormFactors)]
		currObj = parseGau.GauPolyBasis(exponents, [currCoeffs], label=angMom)
		outObjs.append(currObj)
	return outObjs

