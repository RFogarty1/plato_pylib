#!/usr/bin/python3

import itertools as it
import re
import types
from plato_pylib.shared.ucell_class import UnitCell
from plato_pylib.shared.energies_class import EnergyVals 
from . import parse_xyz_files as parseXyzHelp
from ..shared import custom_errors as errorHelp

import pycp2k

RYD_TO_EV = 13.6056980659
HART_TO_EV = 2*RYD_TO_EV


def parseGeomFromCpInputFile(inpFile):
	""" Gets a plato_pylib UnitCell object when passed a path to cp2k input file
	
	Args:
		inpFile: (str) Path to cp2k input file
			 
	Returns
		 outCell: (plato_pylib) UnitCell object containing the geometry in the input file
 
	IMPORTANT:
		Current implementation only works if the unit cell is represented by cell vectors in the input file
		Also currently no attempt to deal with units; you get whatever units you express the geometry in
	"""
	parser = pycp2k.inputparser.CP2KInputParser()
	inpObj = pycp2k.CP2K()
	pyCp2kObj = parser.parse(inpObj, inpFile)
	cellVects = _getCellVectsFromPyCp2kObj(pyCp2kObj)
	useFractCoords, coords = _getCoordsFromPyCp2kObj(pyCp2kObj)
	outCell = UnitCell.fromLattVects(cellVects)
	if useFractCoords:
		outCell.fractCoords = coords
	else:
		outCell.cartCoords = coords

	return outCell

def _getCellVectsFromPyCp2kObj(pyCp2kObj):
	cellSection = pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.CELL
	rawVects = [getattr(cellSection,x) for x in ["A","B","C"]]
	outVects = list()
	for rawVect in rawVects:
		outVect = [float(x) for x in re.sub( '\[[A-Z,a-z]*\]','',rawVect).strip().split()]
		outVects.append(outVect)
	return outVects

def _getCoordsFromPyCp2kObj(pyCp2kObj):
	coordSection = pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COORD
	scaled = coordSection.Scaled
	outCoords = list()
	for currAtom in coordSection.Default_keyword:
		currCoord = [float(x) for x in currAtom.strip().split()[1:]]
		currCoord.append( currAtom.strip().split()[0] )
		outCoords.append(currCoord)
	return scaled,outCoords

def parseCpout(outFile):
	fileAsList = _getFileAsListFromInpFile(outFile)
	parser = _getStandardCpoutParser()
	try:
		outDict = parser.getOutDictFromFileAsList(fileAsList)
	except Exception as e:
	    raise errorHelp.PlatoPylibParseFileError("Something went wrong when parsing the current CP2K output file {}".format(outFile)) from e
	return outDict

def _getStandardCpoutParser():
	outParser = CpoutFileParser()
	condNumbDeco = getDecoToAttachSectionParserToCpoutParser("OVERLAP MATRIX CONDITION NUMBER AT GAMMA POINT", _parseOverlapCondSection)
	condNumbDeco(outParser)
	return outParser

def _getFileAsListFromInpFile(inpFile):
	with open(inpFile,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList

#Want to introduce a way to add a new section to parse without directly modifying the parse source code
#(justified by open-closed principle)
def getDecoToAttachSectionParserToCpoutParser(pattern, parseFunction):
	""" Attaches a function to inpCls (which should be CpoutFileParser) for parsing a section of the output file
	
	Args:
		pattern: (str) The pattern to search for in a single line of a cpout file. Finding the pattern should trigger the parse function
		parseFunction: Function with interface parsedDict, lineIdx = parseFunction(fileAsList, lineIdx). lineIdx is the index in the file where the initial arg is passed (when inp-arg) and where the section is over (when it appears as outArg). ParsedDict is simply a dictionary containing key:val pairs for this section; this is used to update the main dictionary the parser outputs.
 
	Returns
		parseSectionDeco: Decorator for attaching this section parser to the overall CpoutFileParser. After calling parseSectionDeco(CpoutFileParser) any parser instances should be apply "parseFunction" upon finding "pattern" in any file its passed. Thus, this is essentially acting as a hook function for the parser behaviour.

	"""
	def decoFunct(inpCls):
		inpCls.extraSingleLinePatterns.append(pattern)
		inpCls.extraFunctsToParseFromSingleLine.append(parseFunction)
	return decoFunct


#Parse from single line has the interface outDict, lineIdx = parseSectionStartFromLine(fileAsList, lineIdx)
class CpoutFileParser():
	"""Class used to parse CP2K files; NOT meant to be called directly in code; At time of writing _getStandardCpoutParser() is the most sensible way to create this object while the parseCpout function is the best way to parse a CP2K output file

	"""
	extraSingleLinePatterns = list()
	extraFunctsToParseFromSingleLine = list()

	def getOutDictFromFileAsList(self, fileAsList):
		try:
			outDict = self._getOutDictFromFileAsList(fileAsList)
		except Exception as e:
			raise errorHelp.PlatoPylibParseFileError("Something went wrong when parsing the current CP2K output file") from e
		return outDict

	def _getOutDictFromFileAsList(self, fileAsList):
		outDict = self._getInitCp2kOutDict()
		lineIdx=0

		while lineIdx < len(fileAsList):
			currLine = fileAsList[lineIdx].strip()
			if currLine.find("CELL|") != -1:
				outDict["unitCell"], lineIdx = parseCellSectionCpout(fileAsList,lineIdx)
			elif currLine.find("Total energy:") != -1:
				totalE = float( currLine.split()[-1] ) * HART_TO_EV
				outDict["energy"] = totalE
				outDict["energies"] = EnergyVals(dftTotalElectronic=totalE)
				lineIdx += 1
			elif currLine.find("Number of atoms:") != -1:
				outDict["numbAtoms"] += int( currLine.split()[-1]  ) 
				lineIdx += 1
			elif currLine.find("PROGRAM STARTED AT") !=-1: #Reset certain counters every time we find a new start of file
				outDict = self._getInitCp2kOutDict()
				lineIdx += 1
			elif currLine.find("OPTIMIZATION STEP") != -1:
				outDict["multiple_geom_present"] = True
				lineIdx += 1
			elif currLine.find("PROGRAM ENDED") != -1:
				outDict["terminate_flag_found"] = True
				lineIdx += 1
			elif self._patternInExtraSingleLinePatterns(currLine):
				lineIdx = self._updateDictBasedOnFindingSingleLinePatterns(fileAsList, lineIdx, outDict)
			else:
				lineIdx +=1

		if outDict["terminate_flag_found"] is False:
			raise ValueError("Termination flag not found in current cp2k output file")
	
		return outDict

	def _patternInExtraSingleLinePatterns(self, currLine):
		for x in self.extraSingleLinePatterns:
			if currLine.find(x) != -1:
				return True
		return False

	#Should work with multiple parse-functions on the same input pattern; though unlikely that would ever be a good idea (and returned lineIdx will just be that of the LAST matching pattern)
	def _updateDictBasedOnFindingSingleLinePatterns(self, fileAsList, lineIdx, inpDict):
		outLineIdx = lineIdx
		for funct,pattern in it.zip_longest(self.extraFunctsToParseFromSingleLine,self.extraSingleLinePatterns):
			if fileAsList[lineIdx].find(pattern) != -1:
				updateDict, outLineIdx = funct(fileAsList, lineIdx)
				inpDict.update(updateDict)
		return outLineIdx


	def _getInitCp2kOutDict(self):
		outDict = dict()
		outDict["numbAtoms"] = 0
		outDict["multiple_geom_present"] = False #Probably actually a useless output
		outDict["terminate_flag_found"] = False
		return outDict

def parseCellSectionCpout(fileAsList, lineIdx):
	lattParams = list()
	lattAngles = list()

	while lineIdx < len(fileAsList):
		currLine = fileAsList[lineIdx].strip()
		if currLine.find("CELL|") == -1:
			break
		elif currLine.find("Vector")!=-1 and currLine.find("angstrom")!=-1:	
			lattParams.append( currLine.split()[-1] )
			lineIdx += 1
		elif currLine.find("Angle")!=-1 and currLine.find("degree")!=-1:
			lattAngles.append( currLine.split()[-1] )
			lineIdx += 1
		else:
			lineIdx += 1


	unitCell = UnitCell(lattParams=lattParams,lattAngles=lattAngles)

	return unitCell,lineIdx



def parseMOInfo(inpFile:"normal cpout but print MO keyword must be used"):
	with open(inpFile,"rt") as f:
		fileAsList = f.readlines()

	lineIdx = 0
	startSectStr = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
	notInStartSectStr = "MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER SCF STEP 0"
	outDict = _getInitDictForParseMO()

	while lineIdx < len(fileAsList):
		currLine = fileAsList[lineIdx]
		if (startSectStr in currLine) and (notInStartSectStr not in currLine):
			lineIdx = _parseSingleMoKpointSection(fileAsList, lineIdx, outDict)
		else:
			lineIdx += 1


	#Change lists to None if empty

	return outDict

def _getInitDictForParseMO():
	outDict = {"eigenvals":list(),
	           "occvals": list(),
	           "efermi":None}

	return outDict

def _parseSingleMoKpointSection(fileAsList, lineIdx, inpDict):
	allEigs, allOccs = list(), list()
	sectStartStr = "MO index"
	sectEndStr = "Sum"
	sectStart = False
	while lineIdx < len(fileAsList):
		currLine = fileAsList[lineIdx]
		if sectStartStr in currLine:
			sectStart = True
		elif sectEndStr in currLine:
			sectStart = False
		elif sectStart:
			splitLine = currLine.strip().split()
			allEigs.append( float(splitLine[1])*HART_TO_EV ), allOccs.append( float(splitLine[2]) )
		elif "Fermi energy" in currLine:
			eFermi = float( currLine.strip().split()[-1] ) * HART_TO_EV
			break

		lineIdx += 1

	inpDict["eigenvals"].append(allEigs)
	inpDict["occvals"].append(allOccs)
	inpDict["efermi"] = eFermi
	return lineIdx


def _parseOverlapCondSection(fileAsList, lineIdx):
	outDict = dict()
	outDict["overlap_condition_number"] = None
	retOutObj = False
	outObj = types.SimpleNamespace( estimate=types.SimpleNamespace(oneNorm=None),
	                                diag=types.SimpleNamespace(oneNorm=None, twoNorm=None) )  
 
	endStr = "Number of electrons"
	while endStr not in fileAsList[lineIdx]:
		if "1-Norm Condition Number (Estimate)" in fileAsList[lineIdx]:
			lineIdx += 1
			outObj.estimate.oneNorm = float( fileAsList[lineIdx].split("=")[1].strip().split()[0] )
			retOutObj = True
		elif "Condition Numbers using Diagonalization" in fileAsList[lineIdx]:
			lineIdx += 1
			outObj.diag.oneNorm = float( fileAsList[lineIdx].split("=")[-1].strip().split()[0] )
			lineIdx += 1
			outObj.diag.twoNorm = float( fileAsList[lineIdx].split("=")[-1].strip().split()[0] )
			retOutObj = True
		else:
			lineIdx += 1
	
	if retOutObj:
		outDict["overlap_condition_number"] = outObj

	return outDict,lineIdx-1


def parseXyzFromGeomOpt(inpFile):
	outFileStr = _getFileStrFromInpFile(inpFile)
	fileAsList = [x for x in outFileStr.split("\n") if x.strip()!='']

	#Step 1 is to split the file up into individual strings for an xyz parser
	lineIdx = 0
	while lineIdx < len(fileAsList):
		nAtomsLine, commentLine = fileAsList[lineIdx], fileAsList[lineIdx+1]
		commentLine = commentLine.replace(","," ")
		nAtoms = int( nAtomsLine.strip() )
		geomIdx = int(commentLine.strip().split()[2])
		if (geomIdx==1): #This means we only take geometries from the current job; not previous jobs of the same name
			startLines, endLines = list(), list()

		startLines.append(lineIdx)
		endLines.append(lineIdx+nAtoms+1)
		lineIdx += nAtoms + 2 #1 line per atom, 1 for comment line and we then need 1 more to move onto the next step

	#Get the parsed xyz dicts for each
	parsedGeoms = list()
	for start,end in it.zip_longest(startLines,endLines):
		currXyzAsList = fileAsList[start:end+1]
		parsedGeoms.append( parseXyzHelp._parseStandardXyzFile(currXyzAsList) )

	#Convert into the outDict format we want
	outDict = dict()
	outDict["all_geoms"] = parsedGeoms

	return outDict

def _getFileStrFromInpFile(inpFile):
	with open(inpFile,"rt") as f:
		outStr = f.read()
	return outStr
