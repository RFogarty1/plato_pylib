#!/usr/bin/python3

''' Purpose of these functions is to aid in parsing files from tb1 and tb2 '''

from ..shared.ucell_class import UnitCell

def parsePlatoOutFile(inpFilePath):
	with open(inpFilePath,"rt") as f:
		inpFileList = f.readlines()

	if _isDftFile(inpFileList):
		return parseDftFile(inpFileList)

	outDict = dict()
	lineIdx = 0
	while lineIdx < len(inpFileList):
		if (inpFileList[lineIdx].find("E0") != -1) and (inpFileList[lineIdx].find("Eatom") == -1):
			lineIdx, outDict["energies"] = parseEnergiesPlatoOutFile(inpFileList, lineIdx) #tb1-specific
		elif inpFileList[lineIdx].find("Energy\n") != -1:
			lineIdx, outDict["energies"] = parseEnergiesPlatoTb2(inpFileList, lineIdx)
		elif inpFileList[lineIdx].find("Primitive translation vectors") != -1:
			lineIdx, outDict["unitCell"] = parseUnitCellTb1(inpFileList, lineIdx)
		elif inpFileList[lineIdx].find("Cell vectors in") != -1:
			lineIdx, outDict["unitCell"] = parseUnitCellTb2(inpFileList, lineIdx)
		elif inpFileList[lineIdx].find("Number of atoms") != -1:
			outDict["numbAtoms"] = int( inpFileList[lineIdx].strip().split()[-1] )
			lineIdx += 1
		elif inpFileList[lineIdx].find("Geometry (xyz format in Angstroms)") != -1:
			lineIdx += 2
			outDict["numbAtoms"] = int( inpFileList[lineIdx].strip().split()[-1] )
			lineIdx += 1
		else:
			lineIdx += 1

	return outDict


def parseEnergiesPlatoOutFile(inpFileList:list, lineIdx:int):
	objDict = dict() #stores all keys/values needed to create EnergyVals object

	while inpFileList[lineIdx].strip() != '':
		if inpFileList[lineIdx].find("E0") != -1:
			objDict["e0"] = float( inpFileList[lineIdx].strip().split()[2] )
		elif inpFileList[lineIdx].find("E1") != -1:
			objDict["e1"] = float( inpFileList[lineIdx].strip().split()[2] )
		lineIdx += 1

	energies = EnergyVals(**objDict)
	return lineIdx, energies

def parseEnergiesPlatoTb2(inpFileList:list, lineIdx:int):
	objDict = dict()

	passedCohesive = False
	passedSection = False
	while not passedSection:
		line = inpFileList[lineIdx]
		if passedCohesive and line.strip() == '':
			passedSection = True
		elif line.find("Zeroth order energy") != -1:
			objDict["e0"] = float( line.strip().split()[-2] )
		elif line.find("First order energy") != -1:
			objDict["e1"] = float( line.strip().split()[-2] )
		elif line.find("Electron entropy") != -1:
			objDict["entropy"] = float( line.strip().split()[-2] )
		elif line.find('Cohesive') != -1:
			objDict["tb2CohesiveFree"] = float( line.strip().split()[3] )
			passedCohesive = True
		elif line.find("Total energy") != -1:
			objDict["tb2TotalElectronic"] = float( line.strip().split()[3] )
		elif line.find("Atom energy") != -1:
			objDict["atomEnergy"] = float( line.strip().split()[3] )
		lineIdx += 1

	energies = EnergyVals(**objDict)

	return lineIdx, energies

def parseUnitCellTb1(inpFileList:list, lineIdx:int):
	lattVects = list()
	lineIdx+=2
	for x in range(3):
		currLineList = inpFileList[lineIdx].strip().split()
		lattVects.append([float(x) for x in currLineList])
		lineIdx += 1

	unitCell = UnitCell.fromLattVects(lattVects)

	return lineIdx, unitCell

def parseUnitCellTb2(inpFileList:list, lineIdx:int):
	lattVects = list()
	lineIdx+=2
	for x in range(3):
		currLineStr = inpFileList[lineIdx].replace('(','').replace(')','')
		currLineList = currLineStr.replace(',','').split()
		lattVects.append( [ float(y) for y in currLineList[-3:] ] )
		lineIdx += 1

	unitCell = UnitCell.fromLattVects(lattVects)

	return lineIdx, unitCell




def _isDftFile(fileAsList):
	for line in fileAsList:
		if line.find("LCAO density functional theory") != -1:
			return True
	return False

def parseDftFile(fileAsList):
	outDict = dict()
	for idx,line in enumerate(fileAsList):
		if line.find("Sum") != -1:
			totalElectronicE  = float(line.strip().split()[-2])
		elif line.find("Eatom") != -1:
			atomEnergy = float( line.strip().split()[-2] ) #NOT REALLY THE ATOMIC ENERGY(somethine weird instead)
		elif line.find("Number of atoms") != -1:
			outDict["numbAtoms"] = int(line.strip().split()[-1])
		elif line.find("Primitive translation vectors") != -1:
			unusedIdx, outDict["unitCell"] = parseUnitCellTb1(fileAsList, idx)

	outDict["energies"] = EnergyVals(dftTotalElectronic=totalElectronicE)

	return outDict



class EnergyVals():
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.e0 = kwargs.get("e0", None)
		self.e1 = kwargs.get("e1", None)
		self.e2 = kwargs.get("e2", None)
		self.entropy = kwargs.get("entropy", None)
		self.tb2TotalElectronic = kwargs.get("tb2TotalElectronic".lower(), None)
		self.tb2CohesiveFree = kwargs.get("tb2CohesiveFree".lower(), None) 
		self.dftTotalElectronic = kwargs.get("dftTotalElectronic".lower(),None)
		self.atomEnergy = kwargs.get("atomEnergy".lower(), None)


	@property
	def electronicCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree + self.entropy
		elif (self.e0 is not None) and (self.e1 is not None): #In this case the file is Tb1. Meaning these values are relative to the atomic ones
			return self.e0 + self.e1
		else:
			raise ValueError("No information on electronic Cohesive Energy appears to be "
			                 "held in current EnergyVals object")

	@property
	def electronicTotalE(self):
		if self.dftTotalElectronic is not None:
			return self.dftTotalElectronic
		elif self.tb2CohesiveFree is not None:
			return self.e0 + self.e1
		else:
			raise ValueError("No information on total electronic energy is present in current "
			                 "EnergyVals object")

	@property
	def freeCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree
		else:
			raise ValueError("No information on free Cohesive Energy appears to be "
			                 "held in current EnergyVals object")
		


