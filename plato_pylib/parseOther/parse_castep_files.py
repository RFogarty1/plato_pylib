#!/usr/bin/python3

import math
import itertools

from plato_pylib.shared.ucell_class import UnitCell, getTransformedFractCoords


def unitCellObjFromCastepCellFile(cellFilePath:str):
	tokFile = tokenizeCastepCellFileAndRemoveBlockFromKeys(cellFilePath)
	return _getUnitCellObjFromTokenizedCastepCellFile(tokFile)


def _getUnitCellObjFromTokenizedCastepCellFile(tokFile:dict):

	#Step 1 = Get lattice vectors
	lattVectStr = [x for x in tokFile["lattice_cart"].split("\n")[1:]]
	lattVects = list()
	for vectStr in lattVectStr:
		currVect = [float(x) for x in vectStr.replace("\t"," ").split()]
		lattVects.append( currVect )

	#Step 2 = Get the fractional co-ords
	fractCoords = list()
	fractCoordStrList = [x for x in tokFile["positions_frac"].split("\n")]

	for fractStr in fractCoordStrList:
		element, x, y, z = fractStr.replace("\t"," ").strip().split()
		currFract = [float(x), float(y), float(z), element]
		fractCoords.append(currFract)

	return UnitCell.fromLattVects(lattVects,fractCoords=fractCoords)

def tokenizeCastepCellFileAndRemoveBlockFromKeys(cellFilePath:str)->"lower case dict":
	startDict = tokenizeCastepCellFile(cellFilePath)
	outDict = dict()
	for key in startDict.keys():
		outDict[key.lower().strip().replace("%block ","")] = startDict[key]
	return outDict

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

def writeCastepCellFileFromTokens(outFilePath, inpTokens: "dict {keyword:value}"):
	tokens = {k.lower():v for k,v in inpTokens.items()}
	fileStr = "\n"
	reqBlockKeys = ["species_pot"]
	for key, val in tokens.items():
		endVal = "\n\n"
		if val.strip().count('\n') > 0 or key.find("%block") != -1 or (key in reqBlockKeys):
			key = "%block " + key.replace("%block ","") + '\n'
			endVal = "\n%endblock " + key.replace("%block ","") + "\n\n"
		else: 
			key += ' '
		fileStr +=  key + val + endVal

	with open(outFilePath,"wt") as f:
		f.write(fileStr)	



def getCellGeomDictSectionFromUCell(unitCell, units='bohr'):
	outDict = dict()
	outDict["lattice_cart"] = _getCastepLattCartFromUCellClass(unitCell, units=units)
	outDict["positions_frac"] = _getCastepFractPosStrFromUCellClass(unitCell)
	return outDict

def _getCastepLattCartFromUCellClass(unitCell, units='bohr'):
	return unitCell.getCastepCellStr(units=units)

def _getCastepFractPosStrFromUCellClass(unitCell):
	defFormat = unitCell.fractCoords
	outStr = ""
	fmt = "{} {:.6f} {:.6f} {:.6f}"
	for x in defFormat:
		currElement = x[-1]
		outStr += fmt.format(currElement,*x[:3]) + '\n'
	outStr = outStr.strip()
	return outStr
		

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
#	unitCellObj = UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	for idx, unused in enumerate(fractCoords):
		fractCoords[idx].append( atomsInOrder[idx] )

	unitCellObj = UnitCell.fromLattVects(lattVectsFromFile, fractCoords=fractCoords)
#	lattVectsFromCell = unitCellObj.getLattVects()



#	transFractCoords = getTransformedFractCoords(lattVectsFromFile, lattVectsFromCell, fractCoords)
#	for idx, unused in enumerate(transFractCoords):
#		transFractCoords[idx].append( atomsInOrder[idx] )
#	unitCellObj.fractCoords = transFractCoords

	return unitCellObj, lineIdx



