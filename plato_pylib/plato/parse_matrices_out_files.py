
import math
import re
import types
import numpy as np


class InfoNotFoundError(RuntimeError):
	pass


def _getMatricesFromDft2Outfile(outPath):
	""" Parses all hamiltonian and overlap matrices from a dft2 file
	
	Args:
		outPath: Full path to the output file
			
	Returns
		allMatrices: types.SimpleNamespace with two fields, hamil and overlap. Both are lists of square numpy arrays (dtype=complex).
		             The indices refer to the individual k-points (in order)
	
	Raises:
		InfoNotFoundError(RuntimeError) - raised if matrices not present in the outfile
	"""

	with open(outPath,"rt") as f:
		fileAsList = f.readlines()

	outLists = None
	startStr = "Hamiltonian (Ry) and overlap"
	for lineIdx,currStr in enumerate(fileAsList):
		if startStr in currStr:
			outLists = _parseHamilOverlapSectionFromDft2(fileAsList, lineIdx)

	if outLists is None:
		raise InfoNotFoundError("Could not find hamil/overlap matrices in file {}".format(outPath))

	return outLists



def _parseHamilOverlapSectionFromDft2(fileAsList, startLine):
	currLineIdx = startLine
	outLists = list()

	while currLineIdx < len(fileAsList):
		if "State" in fileAsList[currLineIdx]:
			currHamilOverlap, currLineIdx = _parseSingleStateHamilOverlap(fileAsList, currLineIdx)
			outLists.append(currHamilOverlap)
		elif "Eigenvalues" in fileAsList[currLineIdx]:
			break
		else:
			currLineIdx += 1 

	return outLists

#Termination condition bugged
def _parseSingleStateHamilOverlap(fileAsList, startLine):
	currLineIdx = startLine+1
	startedMatricesSection = False

	allSubMatrices = list()
	while currLineIdx < len(fileAsList):
		currStr = fileAsList[currLineIdx]
		if "atom" in currStr:
			currAtomPairMatrices, currLineIdx = parseSingleAtomPairStateHamilOverlap(fileAsList, currLineIdx)
			strAsList = currStr.strip().split()
			atomA, atomB = int(strAsList[2]), int(strAsList[-1])
			currAtomPairMatrices.atomA, currAtomPairMatrices.atomB = atomA, atomB
			allSubMatrices.append(currAtomPairMatrices)
		if ("Eigenvalues" in currStr) or ("State" in currStr):
			break
		else:
			currLineIdx += 1

	#Combine matrices + Create new minimal namespace
	hamil = _getCombinedArrayFromNamespaceList(allSubMatrices,"hamil")
	overlap = _getCombinedArrayFromNamespaceList(allSubMatrices,"overlap")
	outVals = types.SimpleNamespace(hamil=hamil,overlap=overlap)

	return outVals, currLineIdx 

def _getCombinedArrayFromNamespaceList(matricesNamespace, attrib:"hamil or overlap"):
	numbAtoms = max( [max(x.atomA,x.atomB) for x in matricesNamespace] )

	#Step 1 = handle all the horizontal stacking
	allRowMatrices = list()
	for rIdx in range(1,numbAtoms+1): #Need this to ALWAYS run once at least, hence the +1 at the end
		currArray = getattr( _getTwoAtomsSubMatrixObj(matricesNamespace,rIdx,1), attrib )
		for cIdx in range(1,numbAtoms):
			otherArray = getattr( _getTwoAtomsSubMatrixObj(matricesNamespace,rIdx,cIdx+1), attrib ) #+1 due to atom numbering starting at 1
			currArray = np.hstack([currArray, otherArray])
		allRowMatrices.append(currArray)

	#Step 2 = handle the vertical stacking
	outArray = np.vstack(allRowMatrices)

	return outArray


def _getTwoAtomsSubMatrixObj(matricesNamespace, atomAIdx, atomBIdx):
	for x in matricesNamespace:
		if (x.atomA==atomAIdx) and (x.atomB==atomBIdx):
			return x
	return None

def parseSingleAtomPairStateHamilOverlap(fileAsList, startLine):
	currLineIdx = startLine + 1
	startedSection = False

	#Step 1 is to figure out the start/end indices, + figure out the number of rows/cols from that
	while currLineIdx < len(fileAsList):
		currStr = fileAsList[currLineIdx]
		if ("[" in currStr) and (not startedSection):
			startedSection = True
			startLineIdx = currLineIdx
		elif ("[" not in currStr) and (startedSection):
			endLineIdx = currLineIdx-1
			break			
		else:
			currLineIdx += 1	


	#Step 2 is to actually parse it
	numbElements = _findIntegerSquareRoot(currLineIdx-startLineIdx) 
	hamilMatrix = np.zeros(   (numbElements,numbElements), dtype=complex )
	overlapMatrix = np.zeros( (numbElements,numbElements), dtype=complex )

	for currLine in fileAsList[startLineIdx:endLineIdx+1]:
		parsedLine = _parseSingleLineHamilOverlapDft2(currLine)
		rIdx,cIdx = parsedLine.idx
		hamilMatrix[rIdx,cIdx] = parsedLine.hVal
		overlapMatrix[rIdx,cIdx] = parsedLine.sVal

	outVal = types.SimpleNamespace(hamil=hamilMatrix, overlap=overlapMatrix)

	return outVal, endLineIdx+1


def _findIntegerSquareRoot(inpVal):
	outVal = inpVal
	y = (outVal + 1) // 2
	while (y < outVal):
		outVal = y
		y = (outVal + inpVal // outVal) // 2
	
	assert outVal*outVal == inpVal
	return outVal



#(premature) Optimisation attempt. We compile the regexp patterns once and keep fixed across funct calls
def _fixedRegExpDecorator(funct):
	regExpFindMatrixPos = re.compile("\[ \s*[0-9]+\s*,\s*[0-9]+\s*\]")
	regExpFindNumber = re.compile("[0-9]+")
	def outFunct(*args,**kwargs):
		setattr(outFunct, "regExpFindMatrixPos", regExpFindMatrixPos)
		setattr(outFunct, "regExpFindNumber", regExpFindNumber)
		outVal = funct(*args,**kwargs)
		return outVal
	return outFunct

@_fixedRegExpDecorator
def _parseSingleLineHamilOverlapDft2(inpLine):
	#Get the index in the matrix
	findMatrixPosRegExp = _parseSingleLineHamilOverlapDft2.regExpFindMatrixPos
	findNumberRegExp = _parseSingleLineHamilOverlapDft2.regExpFindNumber
	idxPattern = re.findall(findMatrixPosRegExp,inpLine)[0]
	idx = [int(x) for x in re.findall(findNumberRegExp,idxPattern)]

	#Get the values
	subbedLine = inpLine.replace(idxPattern,"").replace("+ i","").strip().split()
	hVal = complex( float(subbedLine[0]), float(subbedLine[1]) )
	sVal = complex( float(subbedLine[2]), float(subbedLine[3]) )

	#Combine all into output arg
	parsedLine = types.SimpleNamespace(hVal=hVal, sVal=sVal, idx=idx)

	return parsedLine






#------------------------>Parsing .ham files below here <---------------------------------------

def getOnSiteDiagHamilTermsFromHamFile(filePath):
	diagTerms = list()

	atomHamils = parseHamFile(filePath)["atomHamils".lower()]

	selfAtomHamils = getOnSiteHamils(atomHamils)
	for x in selfAtomHamils:
		diagTerms.append( x.getDiagTerms() )

	return diagTerms



def parseHamFile(filePath):
	outDict = dict()
	
	with open(filePath,"rt") as f:
		fileAsList = f.readlines()

	startLine, unused = parseHamHeader(fileAsList)
	atomHamils = parseAllAtomHamilBlocks(fileAsList,startLine)

	outDict["atomHamils"] = atomHamils

	outDict = {k.lower():v for k,v in outDict.items()}
	return outDict

#eventually actually return the info, but for now just skip it
def parseHamHeader(fileAsList:list):
	listPos = 0
	inHeader = True
	while (listPos<len(fileAsList)) and inHeader:
		currLine = fileAsList[listPos].strip()
		if len(currLine.split()) == 3: #First atom info; [atomIdx, numbOrbs, numbNebs]
			break
		listPos += 1

	return listPos, None




#Needs to start at line with atomIdx, numbOrbs, numbNebs
def parseAllAtomHamilBlocks(fileAsList:list, listPos:int):
	allAtomHamilBlocks = list() #list of lists. [atomAIdx][NebIdx]
	while listPos < len(fileAsList):
		listPos, currHamilBlocks = parseSingleAtomHamilBlocks(fileAsList,listPos)
		allAtomHamilBlocks.append( currHamilBlocks )

	return allAtomHamilBlocks

#Needs to start at line with atomIdx, numbOrbs, numbNebs
def parseSingleAtomHamilBlocks(fileAsList:list, listPos:int):
	atomInfoLine = fileAsList[listPos].strip().split()
	atomIdx, atomOrbs, numbNebs= [int(x) for x in atomInfoLine[:3]]

	listPos += 1
	hamilBlocks = list()

	sameAtom = True
	while sameAtom and ( listPos<len(fileAsList) ):
		currLine = fileAsList[listPos].strip()
		if len( currLine.split() ) == 3:
			break
		else:
			listPos, currBlock = parseSingleAtomPairHamil(fileAsList, listPos, atomIdx, atomOrbs)
			hamilBlocks.append( currBlock )	

	return listPos, hamilBlocks

#Start at displacement vector line
def parseSingleAtomPairHamil(fileAsList:list, listPos:int, atomAIdx:int, atomAOrbs:int):
	currLine = fileAsList[listPos].strip().split()
	atomBIdx, dispVector = currLine[0],currLine[1:]
	sameBlock = True
	listPos += 1

	hValsByCol = list() 

	while sameBlock and ( listPos<len(fileAsList) ):
		currLine = fileAsList[listPos].strip().split()
		if int(currLine[0]) == 1: #Always start at first row, so if first numb isnt 1 it must be the next neb idx
			hValsByCol.append( list() )
			for idx in range(atomAOrbs):
				currLine = fileAsList[listPos].strip().split()
				hValsByCol[-1].append( currLine[2] )
				listPos+=1
		else:
			break

	numbOrbsB = len(hValsByCol)
	hamilData = np.zeros( (atomAOrbs,numbOrbsB) )
	for row in range(atomAOrbs):
		for col in range(numbOrbsB):
			hamilData[row][col] = hValsByCol[col][row]

	atomPairHamil = AtomPairHamilSubBlock(atomAIdx, atomBIdx, dispVector, hamilData)

	return listPos, atomPairHamil


def getOnSiteHamils(atomHamils:"list of AtomPairHamilSubBlock objs"):
	onSiteHamils = list()
	for atomIdx,unused in enumerate(atomHamils):
		for hamil in atomHamils[atomIdx]:
			if (hamil.distance==0) and (hamil.atomA == hamil.atomB):
				onSiteHamils.append(hamil)
	return onSiteHamils


class AtomPairHamilSubBlock():
	def __init__(self, atomAType, atomBType, dispVector, hVals:"np array"):
		self.atomA = int(atomAType)
		self.atomB = int(atomBType)
		self.dispVector = [float(x) for x in dispVector]
		self.h = hVals

	@property
	def distance(self):
		return ( sum( [(x)**2 for x in self.dispVector] ) )**0.5

	def getDiagTerms(self):
		return list(self.h.diagonal())






