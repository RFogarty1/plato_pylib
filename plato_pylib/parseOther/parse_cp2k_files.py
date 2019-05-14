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


