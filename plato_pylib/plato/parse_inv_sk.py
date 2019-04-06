#!/usr/bin/python3

import math
import itertools

FLIPSHELLS=False #Plato used to label shells the wrong way round; this corrects for that if True

def parseInvSK(inpFile):
	with open(inpFile,"rt") as f:
		fileAsList = f.readlines()

	allParsedLines = list()
	screenFunctPresent = _lineHasScreenFunctField(fileAsList[0])
	for line in fileAsList[1:]:
		if line.split(",")[0]=='x1':
			allParsedLines = list() #Another run was appended to this one, so we reset the parsed data
			screenFunctPresent = _lineHasScreenFunctField(line)
		else:
			allParsedLines.append( InvSKField.fromOutputFileStr(line, screenFunctPresent) )

	return InvSKAllData(allParsedLines)


def _lineHasScreenFunctField(line:str):
	splitLine = line.strip().split(",")
	if all([x.find("Screen Funct")==-1 for x in splitLine]):
		return False
	else:
		return True


class InvSKAllData:
	def __init__(self, invSKFieldObjs:iter):
		self.invSKObjs = list(invSKFieldObjs)


	def __eq__(self,other):
		if type(other) is type(self):
			return self.__dict__ == other.__dict__
		else:
			return NotImplemented
		

	def getShellIndices(self,atomIdx=0):
		atomIdxToShellAttr = {0:"shellA",1:"shellB"}
		shellAttr = atomIdxToShellAttr[atomIdx]
		shellIndices = list()
		for obj in self.invSKObjs:
			currIdx = getattr(obj,shellAttr)
			if currIdx not in shellIndices:
				shellIndices.append(currIdx)

		return shellIndices

	def getAllValsOrbPair(self,valType,shellA,shellB, bondType=None):
		rVsValType = list()

		if bondType is not None:
			valType = self._modValTypeBasedOnBondType(valType,bondType)

		orbIdStrFormat = "{} {}" #used to find instances with same shells/ang mom involved
		orbIdStr = orbIdStrFormat.format(int(shellA),int(shellB))
		for obj in self.invSKObjs:
			idStr = obj.getOrbRepStr(fmt=orbIdStrFormat)
			if idStr==orbIdStr:
				rVsValType.append(   (obj.dist,obj.__getattribute__(valType))   )
		return rVsValType

	def _modValTypeBasedOnBondType(self,valType, bondType):
		if valType.lower() == "hVal".lower():
			if bondType.lower() == "sigma":
				return "hValSigma"
			elif bondType.lower() == "pi":
				return "hValPi"
			elif bondType.lower() == "delta":
				return "hValDelta"
		elif valType.lower() == "sVal".lower():
			if bondType.lower() == "sigma":
				return "sValSigma"
			elif bondType.lower() == "pi":
				return "sValPi"
			elif bondType.lower() == "delta":
				return "sValDelta"
		elif valType.lower() == "screenFunctAngDep".lower():
			if bondType.lower() == "sigma":
				return "screenFunctSigma"
			elif bondType.lower() == "pi":
				return "screenFunctPi"
			elif bondType.lower() == "delta":
				return "screenFunctDelta"

		print("Warning: bondType {} passed with valType {}, which doesnt seem to depend on bondType".format(bondType,valType))
		return valType

	def appendInvSKField(self, invSKField):
		self.invSKObjs.append(invSKField)

	def addInvSKParsedFileData(self, parsedFile:"obj of same class"):	
		self.invSKObjs.extend ( parsedFile.invSKObjs )

	def removeXtalFieldTerms(self, distTol=1e-9):
		newList = list()
		for x in self.invSKObjs:
			if x.dist > distTol:
				newList.append(x)
		self.invSKObjs = newList


class InvSKField:
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.posA = kwargs.get("posA".lower())
		self.posB = kwargs.get("posB".lower())
		self.shellA = kwargs.get("shellA".lower())
		self.shellB = kwargs.get("shellB".lower())
		self.lA = kwargs.get("lA".lower())
		self.lB = kwargs.get("lB".lower())
		self.dist = kwargs.get("dist")
		self.sError = kwargs.get("sError".lower())
		self.hError = kwargs.get("hError".lower())
		self.hValSigma = kwargs.get("hValSigma".lower())
		self.hValPi = kwargs.get("hValPi".lower(),None)
		self.hValDelta = kwargs.get("hValDelta".lower(),None)
		self.sValSigma = kwargs.get("sValSigma".lower(),None)
		self.sValPi = kwargs.get("sValPi".lower(),None)
		self.sValDelta = kwargs.get("sValDelta".lower(),None)
		self.screenFunct = kwargs.get("screenFunct".lower(),None)
		self._eqTol = 1e-7 #Two numbers considered equal if < this

		#Set the ang-dep screening function
		screenFunctAngDep = kwargs.get("screenFunctAngDep".lower(),None) #Should be a size 3 iter (e.g list most sensible)
		angDepScreenAttrs = ["screenFunctSigma", "screenFunctPi", "screenFunctDelta"]
		if screenFunctAngDep is not None:
			[setattr(self,att,val) for att,val in itertools.zip_longest( angDepScreenAttrs,screenFunctAngDep )]
		else:
			[setattr(self,att,None) for att in angDepScreenAttrs]


	def __eq__(self,other):
		#Check comparing the same type
		if type(other) is type(self):
			pass
		else:
			return NotImplemented

		#Check the iterables individually [TODO: these can be refactored with getattr]
		if self.posA is not None:
			posADiff = [abs(a-b) for a,b in itertools.zip_longest(self.posA, other.posA)]
			if not all([x<self._eqTol for x in posADiff]):
				return False
		elif other.posA is not None:
			return False #Means posA defined on other but not self, hence fail

		if self.posB is not None:
			posBDiff = [abs(a-b) for a,b in itertools.zip_longest(self.posB, other.posB)]
			if not all([x<self._eqTol for x in posBDiff]):
				return False
		elif self.posB is not None:
			return False

		#Check all the other attributes; which are either floats or ints
		if len(self.__dict__.keys())!=len(other.__dict__.keys()):
			return False
		for key in self.__dict__.keys():
			if key not in other.__dict__.keys():
				return False

		excludedKeys = ["posA","posB","_screenFunctAngDep"]
		for key in self.__dict__.keys():
			if key in excludedKeys:
				pass
			else:
				if (getattr(self,key) is None) and (getattr(self,key) is None):
					pass
				elif (getattr(self,key) is None) ^ (getattr(self,key) is None):
					return False
				elif ( abs(getattr(self,key)-getattr(other,key)) > self._eqTol ):
					return False

		return True

	@classmethod
	def fromOutputFileStr(cls,outStr, screenFunctPresent=False):
		asList = outStr.strip().split(",")
		posA = [ float(x) for x in asList[:3] ]
		posB = [ float(x) for x in asList[3:6] ]

		if FLIPSHELLS:
			shellB, shellA = int(asList[6]), int(asList[7])
			lB,lA = int(asList[8]), int(asList[9])
		else:
			shellA, shellB = int(asList[6]), int(asList[7])
			lA,lB = int(asList[8]), int(asList[9])
		dist = float(asList[10])

		#Certain fields get shifted up if a screen function is present in the output.
		if screenFunctPresent:
			screenFunct = float(asList[11])
			screenFunctShift = 1
		else:
			screenFunct = None
			screenFunctShift = 0

		#Parse all the matrix elements
		sError, hError = float(asList[11+screenFunctShift]), float(asList[12+screenFunctShift])
		sSigma, hSigma = float(asList[13+screenFunctShift]), float(asList[14+screenFunctShift])	
		try:
			sPi, hPi = float(asList[15+screenFunctShift]), float(asList[16+screenFunctShift])
			if math.isnan(sPi):
				sPi,hPi = None,None
		except IndexError:
			sPi,hPi = None,None

		try:
			sDelta,hDelta = float(asList[17+screenFunctShift]), float(asList[18+screenFunctShift])
			if math.isnan(sDelta):
				sDelta,hDelta = None,None
		except IndexError:
			sDelta,hDelta = None,None

		#Parse the ang-dep part of the screen function
		try:
			angDepScreenFunct = [x for x in [float(y) for y in asList[19+screenFunctShift:22+screenFunctShift]] ]
			for idx,val in enumerate(angDepScreenFunct):
				if math.isnan(val):
					angDepScreenFunct[idx] = None
		except IndexError:
			angDepScreenFunct = None

		outObj = cls(posA=posA, posB=posB, shellA=shellA, shellB=shellB, lA=lA, lB=lB, dist=dist,
		             sError=sError, hError=hError, sValSigma=sSigma, hValSigma=hSigma, sValPi=sPi, hValPi=hPi,
		             sValDelta=sDelta, hValDelta=hDelta, screenFunct=screenFunct, screenFunctAngDep=angDepScreenFunct )
		return outObj

	def getOrbRepStr(self,fmt=None):
		''' String which helps identify the orbitals involved in the interaction '''
		if fmt is None:
			outFormat = "{} {} {} {}"
		else:
			outFormat = fmt

		return outFormat.format(self.shellA,self.shellB)


