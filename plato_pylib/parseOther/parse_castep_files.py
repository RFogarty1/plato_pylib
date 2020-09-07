#!/usr/bin/python3

import math
import itertools as it

import plato_pylib.shared.ucell_class as UCell

from ..shared.energies_class import EnergyVals

#------------->Functions for parsing the *.geom files<------------------

def parseCastepGeomFile(castepGeomFile:str)->"list of dicts":
	with open(castepGeomFile,"rt") as f:
		fileAsList = f.readlines()

	allIters = list()	
	lIdx, maxIdx = 0, len(fileAsList)
	while lIdx < maxIdx:
		line = fileAsList[lIdx]
		if line.find("<--") != -1:
			lIdx, parsedSection = _parseSingleIterGeomFile(fileAsList,lIdx)
			allIters.append(parsedSection)
		else:
			lIdx += 1
	return allIters


def _parseSingleIterGeomFile(fileAsList,filePos)->"int, newFilePos":

	maxIdx =  len(fileAsList)
	#Get cell vectors (later can make it get energies too)
	filePos += 1
	cellVecs = list()
	for x in range(3):
		currLine = fileAsList[filePos].strip()
		currVect = [float(x) for x in currLine.split()[:3]]
		cellVecs.append(currVect)
		filePos += 1

	#Get the Cartesian co-ordinates + skip forces for now
	elementList = list()
	cartCoordList = list()
	while filePos < maxIdx:
		line = fileAsList[filePos]
		if line.find("<--") == -1:
			break
		elif line.find("<-- R") != -1:
			currCart = [float(x) for x in line.strip().split()[2:5]]
			cartCoordList.append(currCart)
			elementList.append( line.strip().split()[0] )
			filePos += 1
		else:
			filePos += 1

	#Combine info into a unitCell object + dict
	fractCoords = list()
	fractXYZ = UCell.getFractCoordsFromCartCoords(cellVecs,cartCoordList)

	for fCoord,element in it.zip_longest(fractXYZ,elementList):
		currFracts = [x for x in fCoord] + [element]
		fractCoords.append(currFracts)

	outDict = dict()
	outDict["unitCell"] = UCell.UnitCell.fromLattVects(cellVecs,fractCoords=fractCoords)

	return filePos, outDict


#------------->Functions for parsing the *.cell files<------------------


def unitCellObjFromCastepCellFile(cellFilePath:str):
	tokFile = tokenizeCastepCellFileAndRemoveBlockFromKeys(cellFilePath)
	return _getUnitCellObjFromTokenizedCastepCellFile(tokFile)


def _getUnitCellObjFromTokenizedCastepCellFile(tokFile:dict):

	#Step 1 = Get lattice vectors
	lattVectStr = [x for x in tokFile["lattice_cart"].split("\n")]
	if len(lattVectStr)==4:
		lattVectStr = lattVectStr[1:] #this should only happen if units are specified on the first line
	lattVects = list()
	for vectStr in lattVectStr:
		currVect = [float(x) for x in vectStr.replace("\t"," ").split()]
		lattVects.append( currVect )

	#Step 2 = Get the fractional co-ords
	try:
		fractCoords = _getFractCoordsFromTokenizedCellFile(tokFile)
		outCell = UCell.UnitCell.fromLattVects(lattVects, fractCoords=fractCoords)
	except KeyError:
		cartCoords =_getCartCoordsFromTokenizedCellFile(tokFile)
		outCell = UCell.UnitCell.fromLattVects(lattVects)
		outCell.cartCoords = cartCoords

	return outCell

#Return None if not found
def _getFractCoordsFromTokenizedCellFile(tokFile):
	fractCoordStrList = [x for x in tokFile["positions_frac"].split("\n")]
	return _getCoordsFromStrList(fractCoordStrList)

def _getCartCoordsFromTokenizedCellFile(tokFile):
	cartCoordStrList = [x for x in tokFile["positions_abs"].split("\n")]
	return _getCoordsFromStrList(cartCoordStrList)

def _getCoordsFromStrList(strList):
	coords = list()
	for coordStr in strList:
		element, x, y, z = coordStr.replace("\t"," ").strip().split()
		currFract = [float(x), float(y), float(z), element]
		coords.append(currFract)

	return coords

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
	for key in sorted(tokens.keys()):
		val = tokens[key]
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
		
#------------->Functions for parsing the *.param files<------------------

def tokenizeCastepParamFile(inpParamPath:str):
	with open(inpParamPath,"rt") as f:
		fileAsList = f.readlines()

	outDict = dict()
	for line in fileAsList:
		if ":" in line:
			key = line.strip().split(":")[0]
			val = line.strip().split(":")[1]
			outDict[key.strip().lower()] = val.strip().lower()

	return outDict

def modCastepParamFile(filePath:str, modKeyVals:dict):
	tokenisedFile = tokenizeCastepParamFile(filePath)
	for key in modKeyVals.keys():
		tokenisedFile[key.lower()] = str( modKeyVals[key] ).lower()
	_writeCastepParamFileFromDict(filePath,tokenisedFile)


def writeCastepParamFileFromDict(filePath, inpDict):
	_writeCastepParamFileFromDict(filePath, inpDict)

def _writeCastepParamFileFromDict(filePath, inpDict):
	outFormat = "{} : {}\n"
	with open(filePath,"wt") as f:
		for key in sorted(inpDict.keys()):
			f.write( outFormat.format(key, inpDict[key]) )

#------------->Functions for parsing the *.castep files<------------------



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
	numbAtoms = energy = energies = unitCellObj = scf_numb_k = scf_k_grid = kin_cut = None

	for lineIdx, line in enumerate(fileList):
		if line.find('Total number of ions in cell') != -1:
			numbAtoms = int(line.strip().split()[-1])
		elif line.find('Final energy') != -1:
			energy = float(line.strip().split()[-2])
			energies = EnergyVals(castepTotalElectronic=energy)
		elif line.find("Unit Cell") != -1:
			unitCellObj, unusedIdx = parseCastepUnitCellSection(fileList, lineIdx)
			if unitCellObj.fractCoords is not None:
				lastFractCoords = unitCellObj.fractCoords
		elif line.find("Cell Contents") != -1: #Sometimes (e.g constant vol opts) this comes without the Unit Cell section
			fractCoords, unused = _parseCellContentsSection(fileList, lineIdx)
			lastFractCoords = fractCoords 
		elif line.find("MP grid size for SCF calculation") != -1:
			scf_k_grid = [ int(x) for x in line.split()[-3:] ]
		elif line.find("Number of kpoints used") != -1:
			scf_numb_k = int( line.split()[-1] )
		elif line.find("plane wave basis set cut-off") != -1:
			kin_cut = float( line.split()[-2] )

	#Deal with case of certain restrained geometry optimisations (where fract coords are only printed at the start)
	unitCellObj.fractCoords = lastFractCoords

	# Create the dictionary and return it
	outDict = { 'numbAtoms' : numbAtoms,
		        'energy' : energy,
	            'energies' : energies,
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

		if fileList[lineIdx].find("Cell Contents") != -1:
			fractCoords, lineIdx = _parseCellContentsSection(fileList, lineIdx) 

		else:
			lineIdx+=1


	#Create the output unitCell
	if len(fractCoords) == 0:
		fractCoords = None

	unitCellObj = UCell.UnitCell.fromLattVects(lattVectsFromFile, fractCoords=fractCoords)

	return unitCellObj, lineIdx


def _parseCellContentsSection(fileAsList, lineIdx):
	""" returns fractCoords from Cell Contents section of castep
		
	  Args:
	    fileAsList(str list): Each entry is 1 line of the castep input file
	    lineIdx(int): The index containing the line "cell contents"

	  Returns
	    fractCoords: nx4 iter with each containing [x,y,z,symbol]. Used to init UnitCell objects
	 
	"""
	finished = False
	while not finished:
		currLine = fileAsList[lineIdx].strip()
		if "Fractional coord" in fileAsList[lineIdx]:
			lineIdx = lineIdx + 3
			fractCoords = list()
			while "xx" not in currLine:
				currLine = fileAsList[lineIdx].strip()
				splitLine = currLine.split()
				if len(splitLine) == 1:
					break
				currCoords = [float(x) for x in splitLine[3:6]] + [splitLine[1]]
				fractCoords.append(currCoords)
				lineIdx = lineIdx + 1
			break
		else:
			lineIdx = lineIdx+1


	return fractCoords, lineIdx




