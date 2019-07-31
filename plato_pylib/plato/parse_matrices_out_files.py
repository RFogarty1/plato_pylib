
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
