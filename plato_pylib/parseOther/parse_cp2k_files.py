#!/usr/bin/python3


from plato_pylib.shared.ucell_class import UnitCell
from plato_pylib.shared.energies_class import EnergyVals 

RYD_TO_EV = 13.6056980659
HART_TO_EV = 2*RYD_TO_EV


def parseCpout(outFile):
	with open(outFile,"rt") as f:
		fileAsList = f.readlines()

	outDict = dict()
	outDict["numbAtoms"] = 0

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
			outDict["numbAtoms"] = 0
			lineIdx += 1
		else:
			lineIdx +=1

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

