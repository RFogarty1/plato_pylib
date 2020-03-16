
import re


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

	return ExponentSetCP2K(exponents,coeffs,lVals,nVal), nLines


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


