#!/usr/bin/python3

import math
from collections import namedtuple

def parseBasFile(basFile:str):
	with open(basFile,"rt") as f:
		fileAsList = f.readlines()	

	#Initialise the output dictionary(NOT ALL KEYS HERE)
	outDict = {"nlpp":None, "density":None, "orbitals":None, "potential":None, "path":basFile}

	firstLine = fileAsList[0].strip().split()	
	if int(firstLine[2]) != 0:
		raise ValueError("pcc != 0 for file {}. Cant correctly parse this file currently".format(basFile))

	lIdx = 0
	while lIdx<len(fileAsList):
		if lIdx==0:
			lineOneDict = parseBasLineOne(fileAsList[lIdx])
			outDict.update(lineOneDict)
			lIdx+=1
		elif lIdx==1:
			orbObjs, lIdx =  parseOrbGenInfo(fileAsList, lIdx,  outDict["numbShells"])  
			nlInfo, lIdx = parseNLPPInfo(fileAsList, lIdx)
		else:#Grid data starts here
			gridInfoDict = parseBasGridVals(fileAsList, lIdx,  orbObjs, nlInfo)
			nlInfo["grids"] = gridInfoDict["nlPP"]
			break
	outDict.update(gridInfoDict)
	outDict["nlPP"] = nlInfo
	outDict["orbinfo"] = orbObjs

	return outDict


def parseBasLineOne(line:str):
	firstLine = line.strip().split()
	outDict = {"numbShells": int(firstLine[0]),
	           "ngridpoints": int(firstLine[1]),
	           "pcc":int(firstLine[2]),
	           "energy": float(firstLine[3]),
	           "cutoff": float(firstLine[4]),
	           "zCore": int(firstLine[5]),
	           "dNu": float(firstLine[6]),
	           "hubbard": float(firstLine[7]),
	           "stoner": float(firstLine[8]),
	           "xcfunctional":int(firstLine[9])}
	return outDict


def parseOrbGenInfo(fileAsList, lIdx, numbShells):
	OrbInfoStruct = namedtuple("OrbInfoStruct",['n', 'l', 'zetaIdx', 'occ', 'energy', 'cutoff'])
	orbGenInfo = list()

	for orbIdx in range(numbShells):
		currLine = fileAsList[lIdx].strip().split()
		n, l, occ = int(currLine[0]), int(currLine[1]), float(currLine[2])
		energy, cutoff = float(currLine[3]), float(currLine[4])
		zetaIdx = countNumbSameOrbNL(orbGenInfo,n,l) + 1
		orbGenInfo.append( OrbInfoStruct(n,l,zetaIdx,occ,energy,cutoff) )
		lIdx+=1

	return orbGenInfo, lIdx

#check number of orbitals already found with quantum numbers n and l
def countNumbSameOrbNL(orbGenInfo, n, l):
	if len(orbGenInfo)==0:
		return 0
	else:
		counter = 0
		for orb in orbGenInfo:
			if (orb.n == n) and (orb.l == l):
				counter += 1
		return counter


def parseNLPPInfo(fileAsList, lIdx):
	outDict = dict()
	firstLineList = fileAsList[lIdx].strip().split()
	lIdx+=1
	secondLineList = fileAsList[lIdx].strip().split()
	lIdx+=1

	outDict["lVals"] = [int(x) for x in firstLineList[1:]]
	outDict["signVals"] = [int(x) for x in secondLineList]
	return outDict, lIdx


def parseBasGridVals(fileAsList, lIdx, orbObjs, nlInfo):
	# Get info into list-of-lists format (all same length)
	allRows = list()
	for line in fileAsList[lIdx:]:
		allRows.append([x for x in line.strip().split()])

	outDict = dict()
	#parse the orbitals
	orbitals = list()
	for orbIdx in range( len(orbObjs) ):
		currGridVals = [ (float(allRows[x][0]), float(allRows[x][orbIdx+1]) ) for x in range(len(allRows)) ]
		currOrbs = basGridInfo.fromOrbInfoTuple(orbObjs[orbIdx], currGridVals)
		orbitals.append(currOrbs)
	outDict["orbitals"] = orbitals

	#parse the density
	densIdx = len(orbObjs)+1
	currGridVals = [ (float(allRows[x][0]), float(allRows[x][densIdx]) ) for x in range(len(allRows)) ]
	outDict["density"] = basGridInfo(currGridVals,"density")

	#Parse the local potential
	locPot = densIdx + 1
	currGridVals = [ (float(allRows[x][0]), float(allRows[x][locPot]) ) for x in range(len(allRows)) ]
	outDict["vloc"] = basGridInfo(currGridVals,"vloc")

	#Parse the nl-pseudopot grid data
	numbNlPP = len(nlInfo["lVals"])
	nlPPGrids = list()
	startIdx = len(orbObjs) + 3 #+3 for density and vloc skipping
	for nlIdx in range( numbNlPP ):
		currGridVals = [ (float(allRows[x][0]), float(allRows[x][startIdx + nlIdx]) ) for x in range(len(allRows)) ]
		nlPPGrids.append( basGridInfo(currGridVals , "nlPP", l=nlInfo["lVals"][nlIdx])  )
	outDict["nlPP"] = nlPPGrids

	#Parse the neutral atom potential
	vnaIdx = len(orbObjs) + numbNlPP + 3 #+3 is to skip over the density and vloc
	currGridVals = [ (float(allRows[x][0]), float(allRows[x][vnaIdx])) for x in range(len(allRows)) ]
	outDict["vna"] = basGridInfo(currGridVals,"vna")

	return outDict





#--------->Functions to rewrite a basis file from a parsed dictionary<---------------

def writeBasFileFromParsedDict(outPath, parsedDict):
	headerLine = _getFirstLineBasFileFromParsedDict(parsedDict)
	shellInfoLines = _shellInfoStringsFromParsedDict(parsedDict)
	nlPPLines = _ppInfoStringsFromParsedDict(parsedDict)
	allGridLines = _allGridStringsFromParsedDict(parsedDict)

	with open(outPath,"wt") as f:
		f.write(headerLine + "\n")
		f.write(shellInfoLines + "\n")
		f.write(nlPPLines + "\n")
		f.write(allGridLines)

def _getFirstLineBasFileFromParsedDict(parsedDict):
	attrsFormStr = [ ("numbShells","{:}"),
	                 ("ngridpoints","{:}"),
	                 ("pcc","{:}"),
	                 ("energy","{:21.14g}"),
	                 ("cutoff","{:21.14g}"),
	                 ("zCore", "{:21.14g}"),
	                 ("dNu", "{:21.14g}"),
	                 ("hubbard", "{:21.14g}"),
	                 ("stoner",  "{:21.14g}"),
	                 ("xcfunctional",  "{:}") ]

	outStr = ""
	for attr,formStr in attrsFormStr:
		outStr += formStr.format( parsedDict[attr] ) + " "
	return outStr.strip() #need to get rid of final space

def _shellInfoStringsFromParsedDict(parsedDict):
	outStr = ""
	orbitals = parsedDict["orbitals"]
	outStrFormat = "{} {} {:21.14g} {:21.14g} {:21.14g}"
	for currOrb in orbitals:
		outStr += outStrFormat.format(currOrb.n, currOrb.l, currOrb.occ, currOrb.energy, currOrb.cutoff) + "\n"

	outStr = outStr.strip() #To remove the final \n

	return outStr

def _ppInfoStringsFromParsedDict(parsedDict):
	outStr = ""
	lVals, signVals = parsedDict["nlPP"]["lVals"], parsedDict["nlPP"]["signVals"]

	formStrA = " ".join(["{}" for unused in range(len(lVals)+1)])	
	formStrB = " ".join(["{}" for unused in range(len(lVals))])


	outStr += formStrA.format(len(lVals),*lVals) + "\n" + " " + formStrB.format(*signVals)
	return outStr

def _allGridStringsFromParsedDict(parsedDict):
	outStr = ""
	nOrbs = len(parsedDict["orbitals"])
	nProj = len(parsedDict["nlPP"]["lVals"])
	nOther = 4 #distance, density, vna, local pot.
	nRows = parsedDict["ngridpoints"]

	formStr = " ".join(["{:21.14g}" for unused in range(nOrbs+nProj+nOther)])
	for rIdx in range(nRows): 
		distVal = parsedDict["density"].gridVals[rIdx][0]
		denVal = parsedDict["density"].gridVals[rIdx][1]
		orbVals = [parsedDict["orbitals"][x].gridVals[rIdx][1] for x in range(nOrbs)]
		nlPPVals = [parsedDict["nlPP"]["grids"][x].gridVals[rIdx][1] for x in range(nProj)]
		vLocVal = parsedDict["vloc"].gridVals[rIdx][1]
		vNaVal = parsedDict["vna"].gridVals[rIdx][1]
		outStr += formStr.format( distVal, *orbVals, denVal, vLocVal, *nlPPVals, vNaVal ) + "\n"

	return outStr

#1 object of this class should hold 1 column of the *.bas info in effect (except the column
#containing the grid co-ords, which is present in all objects)
class basGridInfo():
	def __init__(self, gridVals, valType, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.gridVals = gridVals
		self.valType = valType
		self.n = kwargs.get("n",None)
		self.l = kwargs.get("l",None)
		self.zetaIdx = kwargs.get("zetaIdx".lower(),None)
		self.occ = kwargs.get("occ", None)
		self.energy = kwargs.get("energy", None)
		self.cutoff = kwargs.get("cutoff",None)
	@classmethod
	def fromOrbInfoTuple(cls,orbInfo,gridVals):
		return cls(gridVals, "orb", n=orbInfo.n, l=orbInfo.l, zetaIdx=orbInfo.zetaIdx,
		                            occ=orbInfo.occ, energy=orbInfo.energy, cutoff=orbInfo.cutoff)



	def getGridValsPlatoAOrep(self):
		''' Plato needs the radial contrib to have a prefactor ( sqrt(2l+1)/(r^{l}) ) to
	        properly work with the cubic harmonics '''
		prefactor = ( ((self.l*2) + 1)**0.5 )
		outGrid = list()
		for row in self.gridVals:
			yVal = ( prefactor/(row[0]**self.l) ) * row[1]
			outGrid.append( (row[0], yVal) )
		return outGrid

	def getGridValsPlatoNlPPRep(self):
		''' For non-local pseudopots plato needs to multiply gridpoints by ((2l+1)/4pi)^0.5 '''
		prefactor = math.sqrt( ((self.l*2) + 1)/ (4*math.pi) )
		outGrid = list()
		for row in self.gridVals:
			yVal = ( prefactor/(row[0]**self.l) ) * row[1]
			outGrid.append( (row[0], yVal) )
		return outGrid




