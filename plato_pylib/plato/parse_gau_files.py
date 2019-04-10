#!/usr/bin/python3

import itertools
import math
import os

import numpy as np

def parseGauFile(gauFile:str):
	with open(gauFile,"rt") as f:
		fileAsList = f.readlines()

	lIdx=0
	orbitalFits = list()
	nlPPFits = list()
	weightFits = list()
	outDict = dict()

	outDict["orbitals"] = None
	outDict["nlpp"] = None
	outDict["weightfuncts"]=None

	while lIdx < len(fileAsList):
		currLine = fileAsList[lIdx]
		if lIdx==0:
			outDict["element"] = currLine.split()[0]
			lIdx += 1
		elif lIdx==2:
			orbInfo, lIdx = parseOrbInfo(fileAsList,lIdx)
			outDict.update( orbInfo )
		#ORBITALS
		elif currLine.find("Fitting parameters - wavefunction") != -1:
			currFit, lIdx = parseFit(fileAsList, lIdx+3)
			try:
				orbitalFits.append(  GauPolyBasis.fromIterable(currFit)  )
				outDict["orbitals"] = orbitalFits
			except IndexError:
				pass
		#DENSITY
		elif currLine.find("Fitting parameters - density") != -1:
			currFit, lIdx = parseFit(fileAsList, lIdx+2)
			try:
				densityFit = GauPolyBasis.fromIterable(currFit)
				outDict["density"] = densityFit
			except IndexError:
				print("Could not parse density from file {}".format(gauFile))
		#NEUTRAL ATOM POTENTIAL
		elif currLine.find("neutral atom potential") != -1:
			currFit,lIdx = parseFit(fileAsList, lIdx+2)
			try:
				neutAtomFit = GauPolyBasis.fromIterable(currFit)
				outDict["NeutAtom".lower()] = neutAtomFit
			except IndexError:
				print("Could not parse Neutral Atom Potential from file {}".format(gauFile))
		#NON-LOCAL PSEUDOPOTENTIAL
		elif currLine.find("pseudopotential") != -1:
			currFit,lIdx = parseFit(fileAsList, lIdx+3)
			try:
				nlFit = GauPolyBasis.fromIterable(currFit)
				nlPPFits.append(nlFit)
				outDict["nlpp"] = nlPPFits
			except IndexError:
				pass

		#WEIGHT FUNCTIONS
		elif currLine.find("McWeda Weight Functions") != -1:
			currFit, lIdx = parseFit(fileAsList, lIdx+3)
			try:
				weightFits.append(  GauPolyBasis.fromIterable(currFit)  )
				outDict["weightfuncts"] = weightFits
			except IndexError:
				pass

		else:
			lIdx += 1

	return outDict



def parseFit(fileAsList:list, lineIdx:int):
	fit = list()
	numbGauss = int(fileAsList[lineIdx].strip().split()[0])
	lineIdx+=1

	startIdx = lineIdx
	while lineIdx < startIdx+numbGauss:
		currFit = [float(x) for x in fileAsList[lineIdx].strip().split()]
		fit.append(currFit)
		lineIdx += 1

	return fit, lineIdx


def parseOrbInfo(fileAsList,lineIdx):
	shellLVals, shellNVals = list(), list()
	outDict = dict()

	while (len(fileAsList[lineIdx].strip().split()) != 1) and ( lineIdx < len(fileAsList) ):
		currLine = fileAsList[lineIdx]
		if int(currLine.strip().split()[2]) == 0: #m quantum number; exactly one of these are zero per shell
			shellNVals.append( int(currLine.strip().split()[0])  )
			shellLVals.append( int(currLine.strip().split()[1])  )
		lineIdx += 1
	outDict["shellNVals".lower()], outDict["shellAngMoms".lower()] = shellNVals, shellLVals

	return outDict, lineIdx


def getHeaderStrFromGauFile(filePath):
	''' Used to get section before Fitting parameters listed; this is so i can write '''
	''' a *.gau file from just the fits as long as i have a template file to use '''

	with open(filePath,'rt') as f:
		fileAsList = f.readlines()
	lIdx = 0
	outStr = ""
	while lIdx < len(fileAsList):
		if fileAsList[lIdx].find('#Fitting') == -1:
			outStr += fileAsList[lIdx]
			lIdx += 1
		else:
			break
	return outStr



def getHeaderStrFromParsedBasFile(parsedBasFile:dict):
	#Get all line one info
	element = os.path.split(parsedBasFile["path"])[1].replace('.bas','')
	if len(element) > 2:
		element = element[:2]

	nValenceOrbs = 0
	orbInfoObjs = parsedBasFile["orbinfo"]
	ppInfo = parsedBasFile["nlPP"]
	for orbObj in orbInfoObjs:
		nValenceOrbs += (2*orbObj.l)+1

	outStr = ""
	outStr += "{} {:f} {} {}\n".format(element, parsedBasFile["zCore"], parsedBasFile["numbShells"], nValenceOrbs)
	outStr += "{:f} {:f}\n".format(parsedBasFile["cutoff"], parsedBasFile["energy"])
	for orbObj in orbInfoObjs:
		for idx in range( 2*orbObj.l + 1 ):
			currOcc = orbObj.occ / (2*orbObj.l + 1)
			outStr += "{} {} {} {:17.10g} {:27.20g}\n".format(orbObj.n, orbObj.l, idx, orbObj.energy , currOcc)
	outStr += "{}\n".format( len(ppInfo) )
	
	for idx in range( len(ppInfo) ):
		outStr += "{} {}\n".format( ppInfo["lVals"][idx], ppInfo["signVals"][idx] )

	return outStr

def writeGauFile(filePath, gauData:"dict, format same as the parser", fileHeader):
	keyToHeader = {"orbitals": "#Fitting parameters - wavefunction\n"
	                           "#n = ? l = ? occ = ?.?????? eigenvalue =     -?.??????????\n"
	                           "#               a                           c\n",
	               "density":  "#Fitting parameters - density\n"
	                           "#               a                           c\n",
	               "neutAtom".lower(): "#Fitting parameters - neutral atom potential\n"
	                           "#               a                           c\n",
	               "nlpp":     "#Fitting parameters - non-local pseudopotential\n"
	                           "# l = ?\n"
	                           "#               a                           c\n",
	               "weightfuncts":  "#Fitting parameters - McWeda Weight Functions\n"
	                                "#n(orig orbital) = 3 l(orig orbital) = 0\n"
	                                "#               a                           c\n"
	              }
	outStr = fileHeader
	#Orbitals first to be written
	for currOrb in gauData["orbitals"]:
		outStr += keyToHeader["orbitals"] + currOrb.toGauStr(orb=True)
	#Density/Pot
	outStr += keyToHeader["density"] + gauData["density"].toGauStr()
	outStr += keyToHeader["neutAtom".lower()] + gauData["neutAtom".lower()].toGauStr()
	#Nl-PP
	for currNl in gauData["nlpp"]:
		outStr += keyToHeader["nlpp"] + currNl.toGauStr()
	#dexc string, always constant
	dexcStr = ("#Fitting parameters - dExc\n" +
	           "#               a                           c\n" +
	           "0 0\n")
	outStr += dexcStr
	#weight functs
	if gauData["weightfuncts"] is None:
		pass
	else:
		for currWeight in gauData["weightfuncts"]:
			outStr += keyToHeader["weightfuncts"] + currWeight.toGauStr(orb=True)

	with open(filePath,'wt') as f:
		f.write(outStr)

def parseGauCsv(filePath):
	outDict = {"density":None, "neutatom":None, "orbitals":None, "nlpp":None, "weightfuncts":None}
	with open(filePath,"rt") as f:
		fileAsList = f.readlines()
	
	lIdx = 0
	orbList, nlList, weightList = list(), list(), list()
	while lIdx < len(fileAsList):
		if fileAsList[lIdx].find("Wavefunction")!=-1:
			lIdx, currOrbFit = parseSectionGauCsv(fileAsList, lIdx+2)
			orbList.append(currOrbFit)
			outDict["orbitals"] = orbList
		elif fileAsList[lIdx].find("Density")!=-1:
			lIdx, currFit = parseSectionGauCsv(fileAsList,lIdx+2)
			outDict["density"] = currFit
		elif fileAsList[lIdx].find("Neutral atom") != -1:
			lIdx, currFit = parseSectionGauCsv(fileAsList,lIdx+2)
			outDict["neutatom"] = currFit
		elif fileAsList[lIdx].find("Non-local pseudopotential") != -1:
			lIdx, currFit = parseSectionGauCsv(fileAsList,lIdx+2)
			nlList.append(currFit)
			outDict["nlpp"] = nlList
		elif fileAsList[lIdx].find("McWeda") != -1:
			lIdx, currFit = parseSectionGauCsv(fileAsList,lIdx+2)
			weightList.append(currFit)
			outDict["weightfuncts"] = weightList
		else:
			lIdx += 1

	return outDict

def parseSectionGauCsv(fileAsList,lIdx):
	maxIdx = len(fileAsList)

	fitData = list()
	actData = list()
	while (lIdx<maxIdx):
		splitLine = fileAsList[lIdx].strip().split(',')
		if len(splitLine) != 3:
			break
	
		splitLine = [float(x) for x in splitLine]
		actData.append( (splitLine[0], splitLine[1]) )
		fitData.append( (splitLine[0], splitLine[2]) )
			
		lIdx += 1

	outData = GauCsvGridInfo( actData, fitData )

	return lIdx,outData

def writeGauCsvFile(outPath, gridData:"dict of GauCsvGridInfo objects"):
	outStr = ""

	#Orbitals
	if gridData["orbitals"] is not None:
		for currOrb in gridData["orbitals"]:
			outStr += _getGauCsvSectionStr("orbital", currOrb)

	#Density
	outStr += _getGauCsvSectionStr("density", gridData["density"])

	#Neutal atom pot
	outStr += _getGauCsvSectionStr("neutatom", gridData["neutatom"])

	#Non-loc pseudopot
	if gridData["nlpp"] is not None:
		for nlpp in gridData["nlpp"]:
			outStr += _getGauCsvSectionStr("nlpp", nlpp)

	#Weightfuncts
	if gridData["weightfuncts"] is not None:
		for wfunct in gridData["weightfuncts"]:
			outStr += _getGauCsvSectionStr("weightfunct", wfunct )

	with open(outPath,"wt") as f:
		f.write(outStr)	


def _getGauCsvSectionStr(keyVal:str, gauGridObj:"GauCsvGridInfo obj"):
	keyToHeader = {"orbital": "Wavefunction, n=?, l=?, occ=?, eigenvalue=    ?, ngauss=?",
	               "density": "Density, ngauss=?, npoly=?",
	               "neutatom": "Neutral atom potential, ngauss=?, npoly=?",
	               "nlpp": "Non-local pseudopotential, ngauss=?, npoly=?",
	               "weightfunct": "McWeda, n=?, l=?, , ngauss=?"}
	if gauGridObj is None:
		return ""

	outStr = ""
	outStr += keyToHeader[keyVal] + '\n'
	outStr += "R, Original, Fit" + '\n'
	outStr += gauGridObj.toStr() + '\n'
	return outStr



class GauCsvGridInfo():
	def __init__(self, actVals:iter, fitVals:iter):
		self.actVals = [(x,y) for x,y in actVals]
		self.fitVals = [(x,y) for x,y in fitVals]
		self._eqTol = 1e-8 #Allowance for float-errors when comparing objs for equality

	def toStr(self):
		actVals = np.array(self.actVals)
		fitVals = np.array(self.fitVals)
		outStr = ""
		for x,act,fit in zip( actVals[:,0], actVals[:,1], fitVals[:,1] ):
			outStr += "{:12.5g}, {:12.5g}, {:12.5g}\n".format(x,act,fit)
		return outStr


	def __eq__(self,other):
		retVal = True

		#Only sensible to compare if both objects use the same tolerance. 
		if not abs(self._eqTol - other._eqTol)<1e-9:
			retVal = False

		#Compare actual grid vals
		for valPairActA, valPairActB in itertools.zip_longest(self.actVals,other.actVals):
			if abs(valPairActA[0]-valPairActB[0])>self._eqTol:
				retVal=False
				break
			if abs(valPairActA[1]-valPairActB[1])>self._eqTol:
				retVal=False
				break

		#Compared fitted grid vals
		for valPairFitA, valPairFitB in itertools.zip_longest(self.fitVals, other.fitVals):
			if abs(valPairFitA[0]-valPairFitB[0])>self._eqTol:
				retVal = False
				break
			if abs(valPairFitA[1]-valPairFitB[1])>self._eqTol:
				retVal=False
				break

		return retVal

class GauPolyBasis():
	def __init__(self, exponents:"list", coeffs:"list of lists", label=None):
		self.label = label
		self.exponents = exponents
		self.parseCoeffs(coeffs)
		self.nGauss = len(self.exponents)
		self._eqTol = 1e-8 #Allowance for float-errors when comparing objs for equality

	@classmethod
	def fromIterable(cls,inpIter):
		exponents = list()
		numbCoeffs = len( inpIter[0] )-1
		coeffs = [list() for x in range(numbCoeffs)]

		for row in inpIter:
			exponents.append(row[0])
			for idx in range(len(coeffs)):
				coeffs[idx].append(row[idx+1])
		return cls(exponents,coeffs)

	def parseCoeffs(self,coeffs):
		self.nPoly = 0
		self.r0Coeffs = coeffs[0]
		try:
			self.r1Coeffs = coeffs[1]
		except IndexError:
			self.r1Coeffs = None
			self.r2Coeffs = None
			return 0
		self.nPoly = 1

		try:
			self.r2Coeffs = coeffs[2]
		except IndexError:
			self.r2Coeffs = None
			pass
		self.nPoly = 2


	def removeSmallCoeffsAndExponents(self, minCoeff, minExp):
		self.removeSmallCoeffs(minCoeff)
		self.removeSmallExponents(minExp)

	def removeSmallCoeffs(self, minCoeff):
		''' Removes coeffs/exponents where abs(coeff) < minCoeff for ALL lists (r0,r1,r2) '''
		removeIndices = list()
		for idx,val in enumerate(self.r0Coeffs):
			if self._isCoeffIdxBelowCutoff(idx, minCoeff): 
				removeIndices.append(idx)
		for idx in reversed(removeIndices): #Need to pop large indices before small indices
			self._popIndex(idx)

	def removeSmallExponents(self,minExp):
		removeIndices = list()
		for idx,val in enumerate(self.exponents):
			if abs(val) < minExp:
				removeIndices.append(idx)
		for idx in reversed(removeIndices): #Need to pop large indices before small indices
			self._popIndex(idx)

	def _isCoeffIdxBelowCutoff(self, idx, minCoeff):
		if abs(self.r0Coeffs[idx]) < minCoeff:
			return True
		if self.r1Coeffs is not None:
			if abs(self.r1Coeffs[idx]) < minCoeff:
				return True
		if self.r2Coeffs is not None:
			if abs(self.r2Coeffs[idx]) < minCoeff:
				return True
		return False

	def _popIndex(self,idx):
		self.exponents.pop(idx)
		self.r0Coeffs.pop(idx)
		if self.r1Coeffs is not None:
			self.r1Coeffs.pop(idx)
		if self.r2Coeffs is not None:
			self.r2Coeffs.pop(idx)

	def toGauStr(self,orb=False):
		outStr = ""

		if orb:
			outStr += "{}\n".format(len(self.exponents))
		else:
			outStr += "{} {}\n".format(len(self.exponents), self._getNumbPowersR()+1) 

		lineForm = "{:17.10g} " + " ".join(["{:27.20g}" for x in range(self._getNumbPowersR()+1)])

		allCoeffsZipped = self._zipCoeffs()
		allCoeffs = list()
		for line in allCoeffsZipped:
			allCoeffs.append([x for x in line])

		for idx in range(len(self.exponents)):
			outStr += lineForm.format( self.exponents[idx], *([x for x in allCoeffs[idx]]) ) + "\n"

		return outStr

	def _zipCoeffs(self):
		nPowers = self._getNumbPowersR()
		if nPowers == 0:
			return [[x] for x in self.r0Coeffs]
		elif nPowers == 1:
			return zip( self.r0Coeffs, self.r1Coeffs )
		elif nPowers == 2:
			return zip( self.r0Coeffs, self.r1Coeffs, self.r2Coeffs )
		return None 

	def _getNumbPowersR(self):
		if self.r0Coeffs is None:
			return -1
		if self.r1Coeffs is None:
			return 0
		if self.r2Coeffs is None:
			return 1
		return 2

	#TODO: This function needs to be optimised (combine all the coeffs/exponents in 1 funct).
	def toGaussFunct(self):
		#Firstly need to make GaussPrimitives
		allCoeffs = 1.0
		gFuncts = [GaussPrim(exp, 1).funct for exp in self.exponents]

		def totGaussFunct(x):
			primVals = [f(x) for f in gFuncts]
			rSqrTerm = x**2
			rQuadTerm = x**4
			total = 0.0

			#r^0 contributions
			for idx,c in enumerate(self.r0Coeffs):
				total += c*primVals[idx]

			#Higher r power contribs (2 and 4)
			if self.r1Coeffs is not None:
				for idx,c in enumerate(self.r1Coeffs):
					total += c*primVals[idx]*rSqrTerm
			if self.r2Coeffs is not None:
				for idx,c in enumerate(self.r2Coeffs):
					total += c*primVals[idx]*rQuadTerm
			return total
		return totGaussFunct

	def __eq__(self,other):
		retVal = True
		if self.nPoly != other.nPoly:
			retVal = False

		#Only sensible to compare if both objects use the same tolerance. 
		if not abs(self._eqTol - other._eqTol)<1e-9:
			retVal = False
		
		expDiffs = [abs(a-b) for a,b in itertools.zip_longest(self.exponents, other.exponents)]
		if not all([x<self._eqTol for x in expDiffs]):
			retVal = False

		#Check r0 coeffs
		r0Diffs = [abs(a-b) for a,b in itertools.zip_longest(self.r0Coeffs, other.r0Coeffs)]
		if not all([x<self._eqTol for x in r0Diffs]):
			retVal = False

		#Check r1 Coeffs if present
		if (self.r1Coeffs is None) and (other.r1Coeffs is None):
			if (self.r2Coeffs is None) and (other.r2Coeffs is None):
				pass
			else:
				retVal = False
		elif (self.r1Coeffs is None) or (other.r1Coeffs is None):
			retVal = False
		else: #Final case means r1Coeffs present in both objects
			r1Diffs = [abs(a-b) for a,b in itertools.zip_longest(self.r1Coeffs, other.r1Coeffs)]
			if not all([x < self._eqTol for x in r1Diffs]):
				retVal = False

		#Check r2 Coeffs if present
		if (self.r2Coeffs is None) and (other.r2Coeffs is None):
			pass
		elif (self.r2Coeffs is None) or (other.r2Coeffs is None):
			retVal = False
		else:
			r2Diffs = [abs(a-b) for a,b in itertools.zip_longest(self.r2Coeffs, other.r2Coeffs)]
			if not all( [x<self._eqTol for x in r2Diffs] ):
				retVal = False

		return retVal



#Holds parameters for each prim. gaussian plus a function f(R) that gives val
#of the primitive at position R
class GaussPrim:
	def __init__(self,exp:"float, the exponent",coeff:"float, the coeff"):
		self._exp = exp #cant use self.exp 1st time since updateFunct needs both exp and coeff
		self._coeff = coeff
		self.updateFunct()

	def updateFunct(self):
		def gaussFunct(dist):
#			print("self.exp = {}".format(self.exp))
			val = self.coeff * math.exp( -1*self.exp*(dist**2) )
			return val
		self.funct = gaussFunct

	@property
	def exp(self):
		return self._exp
	@exp.setter
	def exp(self,value):
		self._exp = value
		self.updateFunct()

	@property
	def coeff(self):
		return self._coeff
	@coeff.setter
	def coeff(self,value):
		self._coeff = value
		self.updateFunct()



#	def getPolyParams(self,polyIdx:int)->"list of lists [ [exp1,c1], [exp2,c2] ]":
#		currParams = list()
#		for exp,coeff in itertools.zip_longest(self.exponents[polyIdx], self.coeffs[polyIdx]):
#			currParams.append( [exp,coeff] )
#		return currParams

