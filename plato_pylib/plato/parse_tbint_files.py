#/usr/bin/python3

'''Purpose of these functions are to help parse the *.adt/*.bdt etc. that come from tbint'''

import itertools
import os

import numpy as np

from .parse_bas_files import parseBasFile

#Format4 file headings and how they map to the python3 parser output
BDT_KEYS_TO_PARSER_KEYS = dict()
BDT_KEYS_TO_PARSER_KEYS["hopping"] = "hopping"
BDT_KEYS_TO_PARSER_KEYS["overlap"] = "overlap"
BDT_KEYS_TO_PARSER_KEYS["xtalFieldXc".lower()] = "crystalFieldXc" #This isnt actually present in the bdt files
BDT_KEYS_TO_PARSER_KEYS["xtalFieldNoXc".lower()] = "crystalFieldNonXc" 
BDT_KEYS_TO_PARSER_KEYS["xtalFieldTotal".lower()] = "crystalFieldTotal"
BDT_KEYS_TO_PARSER_KEYS["kinetic"] = "kinetic"
BDT_KEYS_TO_PARSER_KEYS["snIntegrals".lower()] = "snInts"
BDT_KEYS_TO_PARSER_KEYS["pairPot".lower()] = "pairPot"
BDT_KEYS_TO_PARSER_KEYS["nijIntegrals".lower()] = "nijInts"
BDT_KEYS_TO_PARSER_KEYS["hop3B2C".lower()] = "hop3B2C"
BDT_KEYS_TO_PARSER_KEYS["pairFunct".lower()] = "pairFunct"
BDT_KEYS_TO_PARSER_KEYS["nonLocPP".lower()] = "nonLocPP"
BDT_KEYS_TO_PARSER_KEYS["HopXcContrib".lower()] = "HopXcContrib"

PARSER_TO_BDT_KEYS = {v:k for k,v in BDT_KEYS_TO_PARSER_KEYS.items()}

#Maintain a dictionary linkng keywords(in bdt file) to the type of integral (i.e. atom based or orbital based)
BDT_FORM4_INT_TYPES = dict()
BDT_FORM4_INT_TYPES["hopping"] = "orb"
BDT_FORM4_INT_TYPES["overlap"] = "orb"
BDT_FORM4_INT_TYPES["xtalFieldXc".lower()] = "orb"
BDT_FORM4_INT_TYPES["xtalFieldNoXc".lower()] = "orb"
BDT_FORM4_INT_TYPES["xtalFieldTotal".lower()] = "orb"
BDT_FORM4_INT_TYPES["kinetic"] = "orb"
BDT_FORM4_INT_TYPES["snIntegrals".lower()] = "orb"
BDT_FORM4_INT_TYPES["pairPot".lower()] = "atom"
BDT_FORM4_INT_TYPES["nijIntegrals".lower()] = "atom"
BDT_FORM4_INT_TYPES["hop3B2C".lower()] = "orb"
BDT_FORM4_INT_TYPES["pairFunct".lower()] = "atom"
BDT_FORM4_INT_TYPES["nonLocPP".lower()] = "orb"
BDT_FORM4_INT_TYPES["HopXcContrib".lower()] = "orb"


#----------------------------------These functions parse the relevant files------------------#

#Generally the parsing function that should be called from the outside
def getIntegralsFromBdt(bdtFilePath, inclPP=True):
	formType = getFormatTypeBdtFile(bdtFilePath)
	if formType == 3:
		return getIntegralsFromBdt_format3(bdtFilePath, inclPP=inclPP)
	elif formType ==4:
		return parseBdtForm4(bdtFilePath)
	else:
		raise ValueError("formType={} is invalid. This error should never happen".format(formType))
	return None

def getFormatTypeBdtFile(inpBdtPath):
	with open(inpBdtPath,"r") as f:
		fileAsStr = f.read()
	if fileAsStr.find("format_3") != -1:
		return 3
	elif fileAsStr.find("format_4")!= -1:
		return 4
	else:
		raise ValueError("{} is neither format_3 or format_4 type of *.bdt file".format(inpBdtPath))


def parseBdtForm4(inpBdtFile: str):

	with open(inpBdtFile) as f:
		fileAsStr = f.read()

	with open(inpBdtFile) as f:
		fileAsList = f.readlines()

	if fileAsStr.find('format_4') == -1:
		raise ValueError("File {} is not a format4 bdt file".format(inpBdtFile))

	#Remove comment lines
	fileAsList = [x for x in fileAsList if not x.startswith('#')]

	outDict = dict()
	for key in BDT_FORM4_INT_TYPES:
		outDict[BDT_KEYS_TO_PARSER_KEYS[key]] = parseIntegralSetInBdtFormat4(key, fileAsList, BDT_FORM4_INT_TYPES[key])

	#Add the pathname and atom names to each
	atomA, atomB = getAtomNamesFromInpBdtFile(inpBdtFile)
	for key in outDict:
		if outDict[key] is not None:
			for val in outDict[key]:
				val.atomAName = atomA
				val.atomBName = atomB
				val.inpFilePath = inpBdtFile

	#Subtract non-xc xtal field from total to get the xc contribution
	if (outDict["crystalFieldNonXc"] is not None) and (outDict["crystalFieldTotal"] is not None):
		xcXtal = list()
		for totXtal, vnaXtal in itertools.zip_longest( outDict["crystalFieldTotal"], outDict["crystalFieldNonXc"]):
			xcXtal.append( comboSimilarIntegrals(totXtal, vnaXtal, "sub") )
		outDict["crystalFieldXc"] = xcXtal

	return outDict



def parseIntegralSetInBdtFormat4(key:str, fileAsList:list, orbType:str):
	#Figure out starting Position for these integrals
	startIdx = [idx for idx,currLine in enumerate(fileAsList) if key.lower() in currLine.lower()]
	if len(startIdx) == 0:
		return None #Integrals not found in this case
	elif len(startIdx) > 1:
		raise ValueError("keyword {} found more than once in bdt file", key)
	print("The value of key = {}".format(key))
	#Figure out how mnay tables we need
	currPos = int(startIdx[0])
	nTables = int(fileAsList[currPos + 1].strip().split()[0])
	if nTables == 0:
		return None
	currPos = currPos + 2 #Gets us to the start of the first table
	#Parse all the integrals
	allInts = list()
	for idx in range(nTables):
		currInts, currPos = parseSingleTableBdtFormat4(fileAsList, currPos, orbType)
		allInts.append( currInts )

	return allInts

def parseSingleTableBdtFormat4(fileAsList, currPos:int, intType):
	if (intType!="orb") and (intType!="atom"):
		raise ValueError("intType = {} is invalid for function parseSingleTableBdtFormat4".format(orbType))


	#Need to parse the orbital info if relevant
	if intType=="orb":
		orbIdxA, orbIdxB = int(fileAsList[currPos].strip().split()[0]), int(fileAsList[currPos].strip().split()[1])
		currPos += 1
		orbAngMomA, orbAngMomB = int(fileAsList[currPos].strip().split()[0]), int(fileAsList[currPos].strip().split()[1])
		currPos += 1
		orbSubIdx = int(fileAsList[currPos].strip().split()[0]) + 1
		currPos += 1

	numbPoints = int( fileAsList[currPos].strip().split()[0] )

	currPos += 1
	intVals, currPos = parseNextIntegralSet(fileAsList, currPos, numbPoints)

	if intType == "atom":
		outObj = TbintIntegrals(integrals=intVals) 
	else:
		outObj = TbintIntegrals(shellA=orbIdxA, shellB=orbIdxB, angMomA=orbAngMomA, angMomB=orbAngMomB, orbSubIdx=orbSubIdx,
		                        integrals=intVals)

	return outObj, currPos


def getIntegralsFromBdt_format3(bdtFilePath, inclPP=False):
	modelFilePath = os.path.join( os.path.split(bdtFilePath)[0], "model.dat" )
	modelParams = parseModelFile(modelFilePath)
	adtFilePathA, adtFilePathB = getAdtFilePathsFromBdt(bdtFilePath)
	if inclPP:
		nlPPInfo = parseBasFile( adtFilePathB.replace(".adt",".bas") ) ["nlPP"]
	else:
		nlPPInfo = None
	shellToAngMomA = parseAdtFile(adtFilePathA)["shellToAngMom"]
	shellToAngMomB = parseAdtFile(adtFilePathB)["shellToAngMom"]
	allBdtInts = parseBdtFile(bdtFilePath, shellToAngMomA, shellToAngMomB, modelParams, nlPPInfo=nlPPInfo)
	return allBdtInts



def parseModelFile(inpModelFile):
	allModelParams = dict()

	#Keys in same order as in input file
	keyOrder = ["model", "pairFunct", "overlap", "crystalField",
	            "threeCentre", "nijInts", "snInts","hop3B2C"]

	with open(inpModelFile,"rt") as f:
		fileList = f.readlines()

	keyCounter = 0
	for line in fileList:
		if line.find('#')==-1 and line.strip()!='':
			allModelParams[ keyOrder[keyCounter] ] = int( line.strip() )
			keyCounter += 1

	if keyCounter < len(keyOrder):
		for key in keyOrder:
			if key not in allModelParams.keys():
				allModelParams[key] = 0

	return allModelParams


def parseAdtFile(inpAdtFile):
	with open(inpAdtFile,"rt") as f:
		inpFileList = f.readlines()

	allOutVals = dict()
	allOutVals["symbol"] = inpFileList[0].strip().split()[0]
	allOutVals["coreCharge"] = float(inpFileList[1].strip().split()[0])
	allOutVals["numbShells"] = int(inpFileList[2].strip().split()[0])
	allOutVals["numbOrbitals"] = int(inpFileList[3].strip().split()[0])
	allOutVals["orbRadius"] = float(inpFileList[4].strip().split()[0])
	allOutVals["atomEnergy"] = float(inpFileList[5].strip().split()[0])

	orbInfoStartIdx = 6
	orbInfoEndIdx = orbInfoStartIdx + allOutVals["numbOrbitals"]

	shellIdxToAngMom = dict()
	shellCounter = 0
	currNL = (-1,-1) #n/l quantum numbers; should never be able to take these vals
	for currLineIdx in range(orbInfoStartIdx, orbInfoEndIdx):
		currLine = inpFileList[currLineIdx].strip().split()
		newNL = (int(currLine[0]), int(currLine[1]))
		if newNL!=currNL:
			currNL = newNL
			shellIdxToAngMom[shellCounter] = newNL[1]
			shellCounter += 1

	if shellCounter != allOutVals["numbShells"]:
		raise ValueError(" Problem parsing angular momentum of individual shells.File says {} shell, i count {}".format(allOutVals["numbShells"],shellCounter))

	allOutVals["shellToAngMom"] = shellIdxToAngMom

	return allOutVals

#This is the main parsing function for format 3 files. This shouldnt really be called from outside
def parseBdtFile(inpBdtFile, shellIdxToAngMomAtomA:dict, shellIdxToAngMomAtomB, modelParams, nlPPInfo=None):
	allBondIntegrals = dict()

	with open(inpBdtFile,"rt") as f:
		inpFileList = f.readlines()
	inpFileList.append("\n") #Blank line at end or i get index errors later if NL-PP not included

	currLineIdx = 1 #Start at 2nd line (first just says format_3)
	currInts, currLineIdx = parseOrbitalBasedIntegrals(inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB)
	if modelParams["overlap"] == 0:
		overlapInts = None
		hopInts = currInts[0]
		kineticInts = None
	else:
		overlapInts, hopInts, kineticInts = currInts[0], currInts[1], currInts[2] 
		if len(kineticInts) == 0:
			kineticInts = None

	if modelParams["crystalField"] == 1:
		crystalFieldTotal = list()
		currInts, currLineIdx = parseOrbitalBasedIntegrals (inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomA) 
		crystalFieldNonXcInts, crystalFieldXcInts = currInts[0], currInts[1]
		for intIdx in range(len(crystalFieldNonXcInts)):
			crystalFieldTotal.append( comboSimilarIntegrals(crystalFieldNonXcInts[intIdx], crystalFieldXcInts[intIdx]) )
	else:
		crystalFieldNonXcInts = None
		crystalFieldXcInts = None
		crystalFieldTotal = None

	pairPot, currLineIdx = parseAtomBasedIntegrals(inpFileList, currLineIdx)
	allBondIntegrals["pairPot"] = [pairPot]
	if pairPot.integrals.shape[0] == 0:
		print("PairPot empty table found")
		allBondIntegrals["pairPot"] = None


	if modelParams["pairFunct"] == 1:
		pairFunct, currLineIdx = parseAtomBasedIntegrals(inpFileList, currLineIdx)
		allBondIntegrals["pairFunct"] = [pairFunct]
	else:
		allBondIntegrals["pairFunct"] = None

	if modelParams["nijInts"] == 1:
		nijInts, currLineIdx = parseAtomBasedIntegrals(inpFileList, currLineIdx)
		allBondIntegrals["nijInts"] = [nijInts]
	else:
		allBondIntegrals["nijInts"] = None

	if modelParams["snInts"] == 1:
		currInts, currLineIdx = parseOrbitalBasedIntegrals(inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomA)
		snInts = currInts[0]
	else:
		snInts = None

	if modelParams["hop3B2C"] == 1:
		currInts, currLineIdx = parseOrbitalBasedIntegrals(inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB)
		hop3B2C = currInts[0]
	else:
		hop3B2C = None

	if nlPPInfo is not None:
		shellIdxToAngMomNlPP = {k:v for k,v in enumerate(nlPPInfo["lVals"]) }
		currInts, currLineIdx = parseOrbitalBasedIntegrals(inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomNlPP)
		nlPPInts = currInts[0]
	else:
		nlPPInts = None

	allBondIntegrals["hopping"] = hopInts 
	allBondIntegrals["overlap"] = overlapInts
	allBondIntegrals["crystalFieldNonXc"] = crystalFieldNonXcInts
	allBondIntegrals["crystalFieldXc"] = crystalFieldXcInts
	allBondIntegrals["crystalFieldTotal"] = crystalFieldTotal
	allBondIntegrals["snInts"] = snInts
	allBondIntegrals["kinetic"] = kineticInts #often an empty list
	allBondIntegrals["hop3B2C"] = hop3B2C
	allBondIntegrals["nonLocPP"] = nlPPInts

	atomAName, atomBName = getAtomNamesFromInpBdtFile(inpBdtFile)
	for key in allBondIntegrals.keys():
		currIntSet = allBondIntegrals[key]
		if currIntSet is not None:
			for currInt in currIntSet:
				currInt.inpFilePath = inpBdtFile
				currInt.atomAName = atomAName
				currInt.atomBName = atomBName
				currInt.intType = key

	return allBondIntegrals	

#currLineIdx must poit to a line with info on orbital indices
def parseOrbitalBasedIntegrals(inpFileList, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB):
	sectionEnd = False
	setCInts = list()
	setBInts = list()
	setAInts = list()

	loopCounter = 0
	while not sectionEnd:
		currLineLength = len( inpFileList[currLineIdx].strip().split() )
		if currLineLength == 2:
			shellA, shellB = inpFileList[currLineIdx].strip().split()
			shellA, shellB = int(shellA), int(shellB)
			orbSubIdx = 1
			currLineIdx += 1
		elif currLineLength == 1:
			maxOrbSubIdx = min( shellIdxToAngMomAtomA[shellA], shellIdxToAngMomAtomB[shellB] ) + 1 #+1 due to me numbering base 1 vs plato base 0
			if orbSubIdx < maxOrbSubIdx:
				orbSubIdx += 1
			else:
				break
		elif currLineLength == 0:
			break
		else:
			raise ValueError("In parseHopOverlapIntegrals len of line {} = {}".format(inpFileList[currLineIdx], currLineLength) )

		minNumberLoops = min( shellIdxToAngMomAtomA[0], shellIdxToAngMomAtomB[0] ) + 1	
		if (loopCounter>=minNumberLoops) and (shellA==0) and (shellB==0):
			sectionEnd=True
			currLineIdx -= 1
			break
		loopCounter += 1

		numbPoints = int(inpFileList[currLineIdx].strip().split()[0])
		currLineIdx += 1 #Gets to the r=0 integral
		currInts, currLineIdx = parseNextIntegralSet(inpFileList, currLineIdx, numbPoints)

		angMomA = shellIdxToAngMomAtomA[shellA]
		angMomB = shellIdxToAngMomAtomB[shellB]
		currSetAInts = currInts[:,:2]
		setAInts.append( TbintIntegrals(shellA=shellA, shellB=shellB, angMomA=angMomA, angMomB=angMomB, orbSubIdx=orbSubIdx, integrals = currSetAInts)  )
		if currInts.shape[1] > 2:
			currSetBInts = np.array( (currInts[:,0], currInts[:,2]) ).transpose()
			setBInts.append( TbintIntegrals(shellA=shellA, shellB=shellB, angMomA=angMomA, angMomB=angMomB, orbSubIdx=orbSubIdx, integrals = currSetBInts)  )
		if currInts.shape[1] > 3:
			currSetCInts = np.array( (currInts[:,0], currInts[:,3]) ).transpose()
			setCInts.append( TbintIntegrals(shellA=shellA, shellB=shellB, angMomA=angMomA, angMomB=angMomB, orbSubIdx=orbSubIdx, integrals = currSetCInts)  )

	allIntegrals = (setAInts, setBInts, setCInts)
	return allIntegrals , currLineIdx

def parseAtomBasedIntegrals(inpFileList, currLineIdx):
	numbPoints = int(inpFileList[currLineIdx].strip().split()[0])
	currLineIdx += 1 #Gets to the r=0 integral
	currInts, currLineIdx = parseNextIntegralSet(inpFileList, currLineIdx, numbPoints)
	atomBasedIntegrals = TbintIntegrals(integrals = currInts)

	return atomBasedIntegrals, currLineIdx

def parseNextIntegralSet(inpFileList, currLineIdx, numbPoints):
	numbCols = len( inpFileList[currLineIdx].strip().split() )
	intSet = np.ones(( numbPoints,numbCols )) * np.nan

	for currRow in range(numbPoints):
		currIntVals = inpFileList[currLineIdx].strip().split()
		for idx, currVal in enumerate(currIntVals):
			intSet[currRow, idx] = float(currVal)
		currLineIdx += 1

	return intSet, currLineIdx

def getBdtRcut(tbintFilePath):
	adtFilePathA, adtFilePathB = getAdtFilePathsFromBdt(tbintFilePath)
	rcutA = parseAdtFile(adtFilePathA)["orbRadius"]
	rcutB = parseAdtFile(adtFilePathB)["orbRadius"]
	rcutAB = rcutA + rcutB
	return rcutAB


#----------------------------------These functions write out relevant files------------------#

#allInts is in the format output by the parsers
def writeBdtFileFormat4(allInts:dict, outPath:str):
	outStr = "format_4\n"
	for key in allInts:
		if allInts[key] is not None:
			outStr += getBdtStr_singleIntSet_format4( allInts[key] , key)
	with open(outPath,"w") as f:
		f.write(outStr)


def getBdtStr_singleIntSet_format4(intsList, intType):
	#Get header info
	intTypeStr = PARSER_TO_BDT_KEYS[intType]
	atomOrOrbType = BDT_FORM4_INT_TYPES[intTypeStr] 
	nTables = len(intsList)
	outList = list()
	outStr = ""

	#Add header to str
	outList.append( intTypeStr )
	outList.append( "#Number of Tables" )
	outList.append( str(nTables) )

	#Loop over integrals to write each individual table
	for idx,currInt in enumerate(intsList):
		if idx != 0:
			outList.append("#{}".format(intTypeStr))
		outList.append( getBdtStrOneTable(currInt, atomOrOrbType) )
	
	return "\n".join(outList) + "\n"

def getBdtStrOneTable(integrals, atomOrOrb):
	outList = list()
	if atomOrOrb == "orb":
		outList.append( "#Orbital Shell Indices" )
		outList.append( "{} {}".format(integrals.shellA, integrals.shellB) )
		outList.append( "#Orbital Angular Momentum" )
		outList.append( "{} {}".format(integrals.angMomA, integrals.angMomB) )
		outList.append( "#Axial Angular Momentum" )
		outList.append( "{}".format(integrals.orbSubIdx - 1) )

	nPoints = integrals.integrals.shape[0]
	outList.append( "#Number of points in table" )
	outList.append( "{}".format(nPoints) )
	formStr = "{:17.10g} {:17.10g}"
	for idx, currInt in enumerate(integrals.integrals):
		outList.append( formStr.format( currInt[0], currInt[1] ) )

	return "\n".join(outList)


#----------------------------------These functions deal with file paths------------------#

def getAdtFilePathsFromBdt(inpBdtFilePath:str):
	baseFolder, fileName = os.path.split(inpBdtFilePath)
	baseFileName = fileName.replace(".bdt","")

	baseFileNameA, baseFileNameB = baseFileName.split('_')

	filePathA = os.path.join(baseFolder, baseFileNameA + ".adt")
	filePathB = os.path.join(baseFolder, baseFileNameB + ".adt")

	return filePathA, filePathB

def getModelFilePathFromBdt(inpBdtFilePath:str):

	baseFolder = os.path.split(inpBdtFilePath)[0]
	modelFilePath = os.path.join(baseFolder, "model.dat")

	return modelFilePath


def getAtomNamesFromInpBdtFile(inpBdtFilePath:str):
	baseFolder, baseFile = os.path.split(inpBdtFilePath)
	atomAName, atomBName = baseFile.replace(".bdt","").split('_')
	return atomAName, atomBName



#Purpose of this function is to create a single TbintIntegrals object
# by summing or subtracting integrals in two objects iff those objects have the same properties
# (except the actual integrals)
def comboSimilarIntegrals(objA:"TbintIntegrals object", objB:"TbintIntegrals object", operator=None):
	if operator is None:
		operator = "add"

	#Check the objects really are "similar" (same except actual integral values)
	errorFormString = "{} does not equal {}"
	numberMisMatches = 0

	if objA.inpFilePath!=objB.inpFilePath:
		numberMisMatches += 1
	elif objA.shellA!=objB.shellA:
		numberMisMatches += 1
	elif objA.shellB!=objB.shellB:
		numberMisMatches += 1
	elif objA.orbSubIdx!=objB.orbSubIdx:
		numberMisMatches += 1
	elif objA.atomAName!=objB.atomAName:
		numberMisMatches += 1
	elif objA.atomBName!=objB.atomBName:
		numberMisMatches += 1
	elif objA.angMomA != objB.angMomA:
		numberMisMatches += 1
	elif objA.angMomB != objB.angMomB:
		numberMisMatches += 1
	elif objA.integrals.shape!=objB.integrals.shape:
		numberMisMatches += 1

	if numberMisMatches != 0:
		print("Problem adding following two objects")
		objA.printProperties()
		objB.printProperties()
		raise ValueError("The two TbintIntegrals objects displayed above are not similar, but an attempt was made to add them")

	#Add the integrals
	integralsA = objA.integrals
	integralsB = objB.integrals
	sumIntegrals = np.array(objA.integrals)
	if operator == "add": 
		sumIntegrals[:,1] = sumIntegrals[:,1] + objB.integrals[:,1] 
	elif operator == "sub":
		sumIntegrals[:,1] = sumIntegrals[:,1] - objB.integrals[:,1]
	else:
		raise ValueError("{} is an invalid value for operator".format(operator))

	#Create new object
	newObj = TbintIntegrals(inpFilePath=objA.inpFilePath, shellA=objA.shellA, shellB=objA.shellB, angMomA=objA.angMomA, 
                            angMomB = objA.angMomB, orbSubIdx=objA.orbSubIdx,atomAName=objA.atomAName,
	                        atomBName=objB.atomBName, integrals=sumIntegrals)
	return newObj


#Not properly tested yet; also not even sure on the interface. Hence consider this a temporary function
def _findIntObjInList(shellA,shellB,orbSubIdx, integralList):
	for currInt in integralList:
		if (currInt.shellA==shellA) and (currInt.shellB==shellB) and (currInt.orbSubIdx==orbSubIdx):
			return currInt
	print("Could not find integrals with shellA={},shellB={},orbSubIdx={}",shellA,shellB,orbSubIdx)
	return None


#----------------------------------Storage class for individual integral tables------------------#

class TbintIntegrals:
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.inpFilePath = kwargs.get("inpFilePath".lower(), None) #str
		self.shellA = kwargs.get("shellA".lower(), None) #int
		self.shellB = kwargs.get("shellB".lower(), None)
		self.angMomA = kwargs.get("angMomA".lower(), None)
		self.angMomB = kwargs.get("angMomB".lower(), None)
		self.orbSubIdx =  kwargs.get("orbSubIdx".lower(), 1) #pp_sigma=1, pp_pi=2; similar for d etc.
		self.integrals = kwargs.get ("integrals".lower(), None) #nx2 np array
		self.atomAName = kwargs.get("atomAName".lower(), None)
		self.atomBName = kwargs.get("atomBName".lower(), None)
		self.intType  = kwargs.get("intType".lower(), None)

	def printProperties(self):
		print("Printing TbintIntegrals Properties")
		print("self.inpFilePath = {}".format(self.inpFilePath))
		print("self.shellA = {}".format(self.shellA))
		print("self.shellB = {}".format(self.shellB))
		print("self.orbSubIdx {}".format(self.orbSubIdx))
		print("self.atomAName = {}".format(self.atomAName))
		print("self.atomBName = {}".format(self.atomBName))
		print("self.intType = {}".format(self.intType))
		print("self.integrals = ")
		print(self.integrals) 

	def __eq__(self,other):
		errorTol = 1e-8
		truthVals = list()
		truthVals.append( self.shellA == other.shellA)
		truthVals.append( self.shellB == other.shellB)
		truthVals.append( self.orbSubIdx == other.orbSubIdx )
		if (self.inpFilePath is not None) and (other.inpFilePath is not None):
			truthVals.append( os.path.abspath(self.inpFilePath) == os.path.abspath(other.inpFilePath) )
		else:
			truthVals.append( self.inpFilePath == other.inpFilePath )
		truthVals.append( self.atomAName == other.atomAName )
		truthVals.append( self.atomBName == other.atomBName )
		truthVals.append( self.angMomA == other.angMomA )
		truthVals.append( self.angMomB == other.angMomB )

		if isinstance(self.integrals,np.ndarray) and isinstance(other.integrals,np.ndarray):
			truthVals.append( np.allclose(self.integrals, other.integrals, atol=errorTol) )
		elif ~isinstance(self.integrals,np.ndarray) and ~isinstance(other.integrals,np.ndarray):
			truthVals.append( self.integrals == other.integrals )
		else:
			return False

#		print("truthVals = {}".format(truthVals))
		return all(truthVals)


	# FORMAT 3 FILES ONLY AT THE MOMENT
	def replaceIntsTbintFile(self, tbIntFile=None):
		if tbIntFile is None:
			tbIntFile = self.inpFilePath
#		if self.intType.lower() != "pairPot".lower():
#			raise ValueError("replaceIntsTbintFile currently only works for pairPot integrals not {} integrals".format(self.intType))

		if self.intType == "crystalFieldTotal":
			print("WARNING: NO INTEGRALS REPLACED BY replaceIntsTbintFile for intType = crystalFieldTotal")
			return 0

		modelFilePath = os.path.join( os.path.split(tbIntFile)[0], "model.dat" )
		modelParams = parseModelFile(modelFilePath)
		adtFilePathA, adtFilePathB = getAdtFilePathsFromBdt(tbIntFile)
		shellIdxToAngMomAtomA = parseAdtFile(adtFilePathA)["shellToAngMom"]
		shellIdxToAngMomAtomB = parseAdtFile(adtFilePathB)["shellToAngMom"]
	
		with open(tbIntFile,"rt") as f:
			inpFileListForm = f.readlines()

		currLineIdx = 1

		#Hopping or overlap
		if self.intType.lower() == "hopping":
			orbCol = 1 #column containing these integrals (numbering starts from zero)
			replaceOrbitalBasedInts(inpFileListForm, currLineIdx, self.integrals, tbIntFile, self.shellA, self.shellB, self.orbSubIdx, orbCol)
			return 0
		elif self.intType.lower() == "overlap":
			orbCol = 0 #column containing these integrals (numbering starts from zero)
			replaceOrbitalBasedInts(inpFileListForm, currLineIdx, self.integrals, tbIntFile, self.shellA, self.shellB, self.orbSubIdx, orbCol)
			return 0
		else:
			unused ,currLineIdx = parseOrbitalBasedIntegrals(inpFileListForm, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB)

		#Xtal field
		if modelParams["crystalField"] == 1:
			if self.intType.lower() == "crystalFieldNonXcInts":
				orbCol = 0 #column containing these integrals (numbering starts from zero)
				replaceOrbitalBasedInts(inpFileListForm, currLineIdx, self.integrals, tbIntFile, self.shellA, self.shellB, self.orbSubIdx, orbCol)
				return 0
			elif self.intType.lower() == "crystalFieldXcInts":
				orbCol = 1 #column containing these integrals (numbering starts from zero)
				replaceOrbitalBasedInts(inpFileListForm, currLineIdx, self.integrals, tbIntFile, self.shellA, self.shellB, self.orbSubIdx, orbCol)
				return 0
			else:
				unused, currLineIdx = parseOrbitalBasedIntegrals(inpFileListForm, currLineIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB)

		#Pair pot
		if self.intType.lower() == "pairpot":
			replaceAtomBasedInts(inpFileListForm, currLineIdx, self.integrals, tbIntFile)
			return 0
		else:
			unused, currLineIdx = parseAtomBasedIntegrals(inpFileListForm, currLineIdx)

		#Nij Ints
		if modelParams["nijInts"] == 1:
			if self.intType.lower() == "nijInts".lower():
				raise ValueError("Havnt implemented replacing nijInts yet")
				return 0
			else:
				unused, currLineIdx = parseAtomBasedIntegrals(inpFileListForm, currLineIdx)
		#SN Ints
		if modelParams["snInts"] == 1:
			if self.intType.lower() == "snInts".lower():
				raise ValueError("Havnt implemented replacing snInts yet")
				return 0
		
		self.printProperties()
		raise ValueError("Function replaceIntsTbintFile couldnt find correct integral to replace for above object")



#-----------Non-Interface functions to replace integrals in FORMAT 3 FILES ONLY ------------------#

def replaceAtomBasedInts(inpFileListForm, startLineIdx, replacementIntegrals:"np array", inpBdtFilePath):
	numbPointsOrig = int(inpFileListForm[startLineIdx].strip().split()[0])

	currInts, endLineIdx = parseNextIntegralSet(inpFileListForm, startLineIdx+1, numbPointsOrig)

	startFileStr = "".join(inpFileListForm[:startLineIdx]) #Part of file before the integrals
	endFileStr = "".join(inpFileListForm[endLineIdx:]) #Part of the file after the integrals

	newNumbPoints  = replacementIntegrals.shape[0]
	newIntegralsList = [ "{:17.10g} {:17.10g}".format(x,y) for x,y in replacementIntegrals ]
	newIntegralsFileStr = str(newNumbPoints) + "\n" + "\n".join( newIntegralsList ) + "\n"

	fullFileStr = startFileStr + newIntegralsFileStr + endFileStr
	
	with open(inpBdtFilePath,"wt") as f:
		f.write(fullFileStr)


#TODO: Refactor as a direct class function, since i need to access shells/sub orbs etc.
def replaceOrbitalBasedInts(inpFileListForm, startLineIdx, replacementIntegrals:"np array", inpBdtFilePath, shellAIdx, shellBIdx, orbSubIdx, orbColumn):
	#Function starts at the start of a large section (e.g. all hopping+overlap integrals). But we need the start line for 
	#this one specific orbital (e.g. pp\sigma{})
	adtFilePathA, adtFilePathB = getAdtFilePathsFromBdt(inpBdtFilePath)
	shellIdxToAngMomAtomA, shellIdxToAngMomAtomB = parseAdtFile(adtFilePathA)["shellToAngMom"], parseAdtFile(adtFilePathB)["shellToAngMom"]
	currIntStartLineIdx = findStartLineOrbBasedInts(inpFileListForm, startLineIdx, shellAIdx, shellBIdx, orbSubIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB)

	#Once we've found the integrals (Mostly lifted from the replace AtomBasedInts)
	numbPointsOrig = int(inpFileListForm[currIntStartLineIdx].strip().split()[0])
	currInts, endLineIdx = parseNextIntegralSet(inpFileListForm, currIntStartLineIdx+1, numbPointsOrig)

	startFileStr = "".join(inpFileListForm[:currIntStartLineIdx]) #Part of file before the integrals
	endFileStr = "".join(inpFileListForm[endLineIdx:]) #Part of the file after the integrals

	newNumbPoints  = replacementIntegrals.shape[0]

	#Replace 1 column of currInts
	assert np.allclose( replacementIntegrals[:,0], currInts[:,0] ), "replacementIntegrals[:,0] = {}\ncurrInts[:,0]={}".format(replacementIntegrals[:,0],currInts[:,0])
	for rowIdx in range(currInts.shape[0]):
		currInts[rowIdx, orbColumn+1] = replacementIntegrals[rowIdx, 1]
	
	newIntegralsList = [ "{:17.10g} {:17.10g} {:17.10g}".format(x,y,z) for x,y,z in currInts ]
	newIntegralsFileStr = str(newNumbPoints) + "\n" + "\n".join( newIntegralsList ) + "\n"

	fullFileStr = startFileStr + newIntegralsFileStr + endFileStr
	
	with open(inpBdtFilePath,"wt") as f:
		f.write(fullFileStr)


#similar to parseOrbitalBasedInts; but just gets the start line Idx for integrals of a certain
#shell and orbital sub-idx.
def findStartLineOrbBasedInts(inpFileListForm, startLineIdx, shellAIdx, shellBIdx, targetOrbSubIdx, shellIdxToAngMomAtomA, shellIdxToAngMomAtomB):
	currLineIdx = startLineIdx
	loopCounter = 0
	sectionEnd = False
	while not sectionEnd:
		currLineLength = len( inpFileListForm[currLineIdx].strip().split() )
		if currLineLength == 2: #This is the case where a new shell is found (so only called once-should refactor out of this nasty loop)
			currShellA, currShellB = ( int(x) for x in inpFileListForm[currLineIdx].strip().split() )
			orbSubIdx = 1
			maxOrbSubIdx = min(shellIdxToAngMomAtomA[currShellA],shellIdxToAngMomAtomB[currShellB]) + 1 #+1 is because im using base 1 vs platos base 0
			currLineIdx += 1 #gets us to the number of points line
		elif currLineLength == 1: #Only true for line stating the number of data points for next integral
			
			if orbSubIdx < maxOrbSubIdx:
				orbSubIdx += 1
			else:
				break

		else:
			raise ValueError("In findStartLineOrbBasedInts len of line {} = {}".format(inpFileListForm[currLineIdx], currLineLength) )
		#Now check if we've hit the next set of orbitals (quit if so)
		if (loopCounter!=0) and (currShellA==0) and (int(currShellB)==0):
			raise ValueError("Could not find integrals for shellA={},shellB={},subIdx={}".format(shellAIdx,shellBIdx,orbSubIdx))
		loopCounter += 1

		if (currShellA==shellAIdx) and (currShellB==shellBIdx) and (orbSubIdx==targetOrbSubIdx):
			return currLineIdx
		else:
			numbPoints = int( inpFileListForm[currLineIdx].strip() )
			currLineIdx += 1 #Gets to the r=0 integral
			unusedInts, currLineIdx = parseNextIntegralSet(inpFileListForm, currLineIdx, numbPoints)

	raise ValueError("Could not find integrals for shellA={},shellB={},subIdx={}".format(shellAIdx,shellBIdx,orbSubIdx))


