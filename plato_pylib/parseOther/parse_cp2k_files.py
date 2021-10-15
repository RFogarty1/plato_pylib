#!/usr/bin/python3

import itertools as it
import re
import types
from plato_pylib.shared.ucell_class import UnitCell
from plato_pylib.shared.energies_class import EnergyVals


from . import parse_xyz_files as parseXyzHelp
from ..shared import custom_errors as errorHelp
from ..shared import unit_convs as uConvHelp
import pycp2k

RYD_TO_EV = uConvHelp.RYD_TO_EV
HART_TO_EV = uConvHelp.HA_TO_EV


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

def parseCpout(outFile, ThrowIfTerminateFlagMissing=True):
	fileAsList = _getFileAsListFromInpFile(outFile)
	parser = _getStandardCpoutParser()

	#TODO: Some way to maintain the ACTUAL terminate flag may be nice
	if ThrowIfTerminateFlagMissing is False:
		def _finalSetTerminateFlagToTrue(instance):
			instance.outDict["terminate_flag_found"] = True
		parser.finalStepsFunctions.append(_finalSetTerminateFlagToTrue)

	try:
		outDict = parser.getOutDictFromFileAsList(fileAsList)
	except Exception as e:
	    raise errorHelp.PlatoPylibParseFileError("Something went wrong when parsing the current CP2K output file {}".format(outFile)) from e
	return outDict

def _getStandardCpoutParser():
	outParser = CpoutFileParser()
	_addSearchWordAndFunctToParserObj("OVERLAP MATRIX CONDITION NUMBER AT GAMMA POINT", _parseOverlapCondSection, outParser)
	_addSearchWordAndFunctToParserObj("BSSE RESULTS", _parseBSSESection, outParser)
	_addSearchWordAndFunctToParserObj("Core Hamiltonian energy", _parseEnergiesSection, outParser)
	_addSearchWordAndFunctToParserObj("T I M I N G", _parseTimingSection, outParser)
	_addSearchWordAndFunctToParserObj("Total number of message passing", _parseNumbProcsSection, outParser)
	_addSearchWordAndFunctToParserObj("CP2K| version string", _parseCompileInfoSection, outParser)
	_addSearchWordAndFunctToParserObj("BSSE CALCULATION", _parseBSSEFragmentsInfo, outParser, handleParsedDictFunct=_handleParsedBSSEFragsInfo)
	_addSearchWordAndFunctToParserObj("Hirshfeld Charges", _parseHirshfeldChargesSection, outParser, handleParsedDictFunct=_handleHirshfeldChargesInfo)
	_addSearchWordAndFunctToParserObj("Mulliken Population Analysis", _parseHirshfeldChargesSection, outParser, handleParsedDictFunct=_handleMullikenChargesInfo)
	_addSearchWordAndFunctToParserObj("ATOMIC FORCES in [a.u.]", _parseAtomicForcesSection, outParser, handleParsedDictFunct=_handleAtomicForcesSection)
	outParser.finalStepsFunctions.append(_parseBSSEFragmentsFinalStepFunct)
	return outParser

def _getFileAsListFromInpFile(inpFile):
	with open(inpFile,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList


def _addSearchWordAndFunctToParserObj(searchWord, funct, parserObj, handleParsedDictFunct=None):
	decoObj = getDecoToAttachSectionParserToCpoutParser(searchWord, funct, handleParsedDictFunct=handleParsedDictFunct)
	decoObj(parserObj)

#Want to introduce a way to add a new section to parse without directly modifying the parse source code
#(justified by open-closed principle)
def getDecoToAttachSectionParserToCpoutParser(pattern, parseFunction, handleParsedDictFunct=None):
	""" Attaches a function to inpCls (which should be CpoutFileParser INSTANCE) for parsing a section of the output file
	
	Args:
		pattern: (str) The pattern to search for in a single line of a cpout file. Finding the pattern should trigger the parse function
		parseFunction: Function with interface parsedDict, lineIdx = parseFunction(fileAsList, lineIdx). lineIdx is the index in the file where the initial arg is passed (when inp-arg) and where the section is over (when it appears as outArg). ParsedDict is simply a dictionary containing key:val pairs for this section; this is used to update the main dictionary the parser outputs.
		handleParsedDictFunct: f(instance, parsedDict) Default of None means we simply update the output dict with the dict parsed from this section (usual desired behaviour). But setting this function explicitly allows for things such as parsing a series of values (e.g. temperature at each MD step) and saving ALL of them into the outptu dict (instance.outDict)
 
	Returns
		parseSectionDeco: Decorator for attaching this section parser to the overall CpoutFileParser. After calling parseSectionDeco(CpoutFileParser) any parser instances should be apply "parseFunction" upon finding "pattern" in any file its passed. Thus, this is essentially acting as a hook function for the parser behaviour.

	"""
	def decoFunct(inpCls):
		inpCls.extraSingleLinePatterns.append(pattern)
		inpCls.extraFunctsToParseFromSingleLine.append(parseFunction)
		inpCls.extraHandleParsedOutputFuncts.append(handleParsedDictFunct)
	return decoFunct


#Parse from single line has the interface outDict, lineIdx = parseSectionStartFromLine(fileAsList, lineIdx)
class CpoutFileParser():
	"""Class used to parse CP2K files; NOT meant to be called directly in code; At time of writing _getStandardCpoutParser() is the most sensible way to create this object while the parseCpout function is the best way to parse a CP2K output file

	"""
	def __init__(self):
		self.extraSingleLinePatterns = list() #Search strings that trigger us to parse a section
		self.extraFunctsToParseFromSingleLine = list() #Functions to parse the relevant sections and return a dictionary AND lineIdx (so we dont re-read lines in this section) 
		self.extraHandleParsedOutputFuncts = list() #These functions map the parsed-dicts to the "global" self.outDict. If set to None then we simply do self.outDict.update(parsedDict) for each section.
		self.finalStepsFunctions = list()

	def getOutDictFromFileAsList(self, fileAsList):
		try:
			outDict = self._getOutDictFromFileAsList(fileAsList)
		except Exception as e:
			raise errorHelp.PlatoPylibParseFileError("Something went wrong when parsing the current CP2K output file") from e
		return outDict

	def _getOutDictFromFileAsList(self, fileAsList):
		self.outDict = self._getInitCp2kOutDict() #Attach to class so we can access it with hook functions
		lineIdx=0

		while lineIdx < len(fileAsList):
			currLine = fileAsList[lineIdx].strip()
			if currLine.find("CELL|") != -1:
				self.outDict["unitCell"], lineIdx = parseCellSectionCpout(fileAsList,lineIdx)
			elif currLine.find("Number of atoms:") != -1:
				self.outDict["numbAtoms"] += int( currLine.split()[-1]  ) 
				lineIdx += 1
			elif currLine.find("PROGRAM STARTED AT") !=-1: #Reset certain counters every time we find a new start of file
				self.outDict = self._getInitCp2kOutDict()
				lineIdx += 1
			elif currLine.find("OPTIMIZATION STEP") != -1:
				self.outDict["multiple_geom_present"] = True
				lineIdx += 1
			elif currLine.find("PROGRAM ENDED") != -1:
				self.outDict["terminate_flag_found"] = True
				lineIdx += 1
			elif self._patternInExtraSingleLinePatterns(currLine):
				lineIdx = self._updateDictBasedOnFindingSingleLinePatterns(fileAsList, lineIdx, self.outDict)
			else:
				lineIdx +=1


		self._applyFinalStepsFunctions()

		if self.outDict["terminate_flag_found"] is False:
			raise ValueError("Termination flag not found in current cp2k output file")

	
		return self.outDict

	def _applyFinalStepsFunctions(self):
		for funct in self.finalStepsFunctions:
			funct(self)

	def _patternInExtraSingleLinePatterns(self, currLine):
		for x in self.extraSingleLinePatterns:
			if currLine.find(x) != -1:
				return True
		return False

	#TODO: Add the ability to change the update function from outside (needed for getting lists)
	#Should work with multiple parse-functions on the same input pattern; though unlikely that would ever be a good idea (and returned lineIdx will just be that of the LAST matching pattern)
	def _updateDictBasedOnFindingSingleLinePatterns(self, fileAsList, lineIdx, inpDict):
		outLineIdx = lineIdx
		for funct,pattern,handleFunct in it.zip_longest(self.extraFunctsToParseFromSingleLine,self.extraSingleLinePatterns, self.extraHandleParsedOutputFuncts):
			if fileAsList[lineIdx].find(pattern) != -1:
				updateDict, outLineIdx = funct(fileAsList, lineIdx)
				if handleFunct is None:
					inpDict.update(updateDict)
				else:
					handleFunct(self, updateDict)
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


#This is the RESULTS section, which comes last
def _parseBSSESection(fileAsList, lineIdx):
	outDict = dict()
	outDict["bsse"] = None
	retOutObj = False
	outObj = types.SimpleNamespace(cpCorrectedTotalEnergy=None)

	endStr = "BSSE-free interaction energy"
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		if "CP-corrected Total energy" in fileAsList[lineIdx]:
			corrE = float( fileAsList[lineIdx].strip().split()[-2] ) * HART_TO_EV
			outObj.cpCorrectedTotalEnergy = corrE
			retOutObj = True
		lineIdx += 1

	if retOutObj:
		outDict["bsse"] = outObj

	return outDict, lineIdx-1

def _handleParsedBSSEFragsInfo(parserInstance, outDict):
	if parserInstance.outDict.get("bsse_fragments",None) is None:
		parserInstance.outDict["bsse_fragments"] = list()

	parserInstance.outDict["bsse_fragments"].append(outDict)
	if parserInstance.outDict.get("energies",None) is not None:
		parserInstance.outDict["bsse_fragments"][-2]["energies"] = parserInstance.outDict["energies"]


def _parseBSSEFragmentsInfo(fileAsList, lineIdx):
	outDict = dict()
	endStr = "-----------------------------"
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		currLine = fileAsList[lineIdx]
		if "FRAGMENT CONF:" in currLine:
			outDict["conf"] = currLine.strip().split()[5]
			outDict["frag_sub_conf"] = currLine.strip().split()[8]
		elif "CHARGE" in currLine:
			outDict["charge"] = int( currLine.strip().split()[3] )
			outDict["multiplicity"] = int( currLine.strip().split()[6] )
		elif "ATOM INDEX" in currLine:
			lineIdx += 2
			atomIndices, atomKinds = list(), list()
			while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
				currLine = fileAsList[lineIdx]
				atomIndices.append( int(currLine.strip().split()[1]) )
				atomKinds.append( currLine.strip().split()[-2] )
				lineIdx+=1
			break

		lineIdx += 1

	outDict["indices"], outDict["kinds"] = atomIndices, atomKinds

	return outDict, lineIdx

def _parseBSSEFragmentsFinalStepFunct(parserInstance):
	if parserInstance.outDict.get("bsse_fragments",None) is not None:
		parserInstance.outDict["bsse_fragments"][-1]["energies"] = parserInstance.outDict["energies"]


def _parseEnergiesSection(fileAsList, lineIdx):
	outDict = dict()
	dftTotalElectronic, dispVal, entropy, fermiE = None, None, None, None
	endStr = "Total energy:"
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		if "Electronic entropic energy" in fileAsList[lineIdx]:
			entropy = float( fileAsList[lineIdx].split()[-1] ) * HART_TO_EV
		if "Dispersion energy" in fileAsList[lineIdx]:
			dispVal = float( fileAsList[lineIdx].split()[-1] ) * HART_TO_EV
		if "Fermi energy:" in fileAsList[lineIdx]:
			pass
			fermiE = float( fileAsList[lineIdx].split()[-1] ) * HART_TO_EV
		lineIdx += 1

	dftTotalElectronic = float( fileAsList[lineIdx].split()[-1] ) * HART_TO_EV
	lineIdx += 1

	outDict["energies"] = EnergyVals(dispersion=dispVal, entropy=entropy, dftTotalElectronic=dftTotalElectronic)
	outDict["energy"] = dftTotalElectronic
	if fermiE is not None:
		outDict["fermi_energy"] = fermiE

	return outDict,lineIdx


def _parseTimingSection(fileAsList, lineIdx):
	outDict = dict()
	endStr = "The number of warnings"
	timingDict = dict()
	subroutineTotals = dict()
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		if "CP2K  " in fileAsList[lineIdx]:
			timingDict["CP2K_total"] = float(fileAsList[lineIdx].strip().split()[-1])
		if "-" not in fileAsList[lineIdx]:
			line = fileAsList[lineIdx]
			if ("SUBROUTINE" not in line) and ("MAXIMUM" not in line) and (line.strip()!=""):
				currKey = line.strip().split()[0]
				currVal = float(line.strip().split()[-1])
				subroutineTotals[currKey] = currVal 
		lineIdx+=1

	timingDict["subroutineTotals"] = subroutineTotals
	outDict["timings"] = types.SimpleNamespace(**timingDict)

	return outDict, lineIdx-1


def _parseNumbProcsSection(fileAsList, lineIdx):
	outDict = dict()
	endStr = "This output is from"
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		currLine = fileAsList[lineIdx]
		if "Total number of message passing processes" in currLine:
			outDict["nMPI"] = int(currLine.strip().split()[-1])
		elif "Number of threads" in currLine:
			outDict["nThreads"] = int(currLine.strip().split()[-1])

		lineIdx +=1

	return outDict, lineIdx


#NOTE: This probably works for ALL charges
def _parseHirshfeldChargesSection(fileAsList, lineIdx):
	outDict = dict()
	endStr = "!-----"
	lineIdx += 1

	outCharges = list()
	parseCharges = True #Flag invented to deal with annoying case of Mulliken charges being mixed with orbital population
	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		currLine = fileAsList[lineIdx]
		if currLine.strip() == "":
			pass
		elif "Orbital" in currLine:
			parseCharges = False #Dont try to parse anything now; but cant break since i want to get to the endStr symbol first
		elif "Atom" in currLine:
			pass
		elif ("Total Charge".lower() in currLine.lower()) and parseCharges:
			outDict["total"] = float( currLine.strip().split()[-1] )
		elif parseCharges:
			currCharge = float( currLine.strip().split()[-1] )
			outCharges.append(currCharge)

		lineIdx += 1

	if parseCharges:
		outDict["charges"] = outCharges

	return outDict, lineIdx


def _handleHirshfeldChargesInfo(parserInstance, outDict):
	parserInstance.outDict["hirshfeld_charges_final"] = outDict

def _handleMullikenChargesInfo(parserInstance, outDict):
	parserInstance.outDict["mulliken_charges_final"] = outDict


def _parseAtomicForcesSection(fileAsList, lineIdx):
	outDict = dict()
	endStr = "SUM OF ATOMIC FORCES"
	lineIdx+=3

	outForces = list()
	while (lineIdx<len(fileAsList)) and (endStr not in fileAsList[lineIdx]):
		currLine = fileAsList[lineIdx]
		splitLine = currLine.strip().split()
		currVals = [float(x) for x in splitLine[-3:]]
		outForces.append(currVals)
		lineIdx+=1

	outDict["forces"] = outForces

	return outDict, lineIdx


def _handleAtomicForcesSection(parserInstance, outDict):
	parserInstance.outDict["forces_final"] = outDict["forces"]


def _parseCompileInfoSection(fileAsList, lineIdx):
	outDict = dict()
	endStr = "is freely available from"

	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		currLine = fileAsList[lineIdx]
		if "CP2K| version string:" in currLine:
			outDict["version_string"] = currLine.replace("CP2K| version string:","").strip()
		if "CP2K| source code revision number:" in currLine:
			outDict["source_code_number"] = currLine.replace("CP2K| source code revision number:","").strip()
		if "CP2K| cp2kflags:" in currLine:
			tempDict, lineIdx = _parseCp2kCompileFlags(fileAsList, lineIdx)
			outDict.update(tempDict)
		lineIdx += 1

	return {"cp2k_compile_info":outDict}, lineIdx

def _parseCp2kCompileFlags(fileAsList, lineIdx):
	outStr = ""
	endStr = "is freely available from"

	startIdx = fileAsList[lineIdx].find("CP2K| cp2kflags: ") + len("CP2K| cp2kflags: ")

	while (endStr not in fileAsList[lineIdx]) and (lineIdx<len(fileAsList)):
		outStr += fileAsList[lineIdx][startIdx:].strip("\n")
		lineIdx += 1
		

	return {"cp2kflags":outStr},lineIdx-1

def parseXyzFromGeomOpt(inpFile, startGeomIdx=1):
	""" Description of function
	
	Args:
		inpFile: Path to xyz file
		startGeomIdx: (int, optional) The index which denotes the FIRST step in a geometry optimisation. This is 1 for geo_opt but 0 in nudged elastic band calculations. The list of out geoms resets when this index is found (such that we ONLY parse results from the most recent optimisation contained in the file)
			 
	Returns
		outDict: Contains "all_geoms" key which contains a list of geometries
 
	Raises:
		 Errors
	"""
	outFileStr = _getFileStrFromInpFile(inpFile)
	fileAsList = [x for x in outFileStr.split("\n") if x.strip()!='']

	#Step 1 is to split the file up into individual strings for an xyz parser
	lineIdx = 0
	startLines, endLines = list(), list()
	while lineIdx < len(fileAsList):
		nAtomsLine, commentLine = fileAsList[lineIdx], fileAsList[lineIdx+1]
		commentLine = commentLine.replace(","," ")
		nAtoms = int( nAtomsLine.strip() )
		geomIdx = int(commentLine.strip().split()[2])
		if (geomIdx==startGeomIdx): #This means we only take geometries from the current job; not previous jobs of the same name
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
