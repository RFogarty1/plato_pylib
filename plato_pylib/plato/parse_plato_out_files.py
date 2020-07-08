#!/usr/bin/python3

''' Purpose of these functions is to aid in parsing files from tb1 and tb2 '''

import re

from ..shared.ucell_class import UnitCell
from ..shared.energies_class import EnergyVals

import numpy as np


def parsePlatoOutFile_energiesInEv(inpFilePath):
	outDict = parsePlatoOutFile(inpFilePath)
	outDict["energies"].convRydToEv()
	return outDict

def parsePlatoOutFile(inpFilePath):
	with open(inpFilePath,"rt") as f:
		inpFileList = f.readlines()

	#Initialize the dictionary we write to
	outDict = dict()
	outDict["scf_is_converged"] = None

	if _isDftFile(inpFileList):
		dftDict = parseDftFile(inpFileList)
		outDict.update(dftDict)
		return outDict

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
		elif inpFileList[lineIdx].find("Self Consistency Cycles") != -1:
			lineIdx += 1
			lineIdx, scfSectionDict = _parseScfSection(inpFileList, lineIdx)
			outDict.update( scfSectionDict )
		else:
			lineIdx += 1

	return outDict


def _parseScfSection(inpFileList:list, lineIdx:int):

	outDict = dict()

	while "SCF Converged" not in inpFileList[lineIdx]:
		lineIdx += 1

	#Get True/False on whether scf was converged
	pattern = r"SCF Converged: \b([a-zA-Z]+)\b"
	isConvStrList = re.findall(pattern, inpFileList[lineIdx])
	assert len(isConvStrList) == 1, "found {} matches for whether scf was converged".format( len(isConvStrList) )
	outDict["scf_is_converged"] = True if isConvStrList[0].lower()=="true" else False

	return lineIdx, outDict

def parseEnergiesPlatoOutFile(inpFileList:list, lineIdx:int):
	objDict = dict() #stores all keys/values needed to create EnergyVals object

	while inpFileList[lineIdx].strip() != '':
		if inpFileList[lineIdx].find("E0") != -1:
			objDict["e0coh"] = float( inpFileList[lineIdx].strip().split()[2] )
		elif inpFileList[lineIdx].find("E1") != -1:
			objDict["e1"] = float( inpFileList[lineIdx].strip().split()[2] )
		elif inpFileList[lineIdx].find("Entropy") != -1:
			objDict["entropy"] = float( inpFileList[lineIdx].strip().split()[2] )
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
			objDict["e0tot"] = float( line.strip().split()[-2] )
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

	objDict["e0coh"] = objDict["e0tot"] - objDict["atomEnergy"]

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



def parseOccFile(inpPath):
	with open(inpPath) as f:
		file_list = f.readlines()

	# Initialise array for k-path and eigen_vals
	numb_k_points = int(file_list[0].strip().split()[0])
	numb_eigens = int(file_list[0].strip().split()[1])

	#TODO: shorten slightly with refactor
	if numb_k_points > 0:
		eigen_vals = np.nan * np.ones(( numb_k_points, numb_eigens  )) # each row is all eigenvals for 1 k point
		k_path = np.nan * np.ones(( numb_k_points, 3  )) # [kpoint index, x, y, z] format
		k_weights = np.nan * np.ones (( numb_k_points, 1 ))
		occs = np.nan * np.ones(( numb_k_points, numb_eigens ))
	else:
		eigen_vals = np.nan * np.ones(( 1 , numb_eigens ))
		k_path = np.array([0.0,0.0,0.0])
		k_weights = np.array((1.0))
		occs = np.nan * np.ones(( 1, numb_eigens))

	
	eigen_count = 0
	for line in file_list[1:]:
		if line.find('K-point') != -1:
			curr_k_point = int(line.strip().split()[1])
			eigen_count = 0
			k_path[curr_k_point - 1,:] = line.strip().split()[2:5]
			k_weights[curr_k_point - 1,:] = line.strip().split()[-1]
		elif line.strip()!='':
			if numb_k_points>0:
				eigen_vals[curr_k_point - 1, eigen_count] = line.strip().split()[0]
				occs[curr_k_point - 1, eigen_count] = line.strip().split()[1]
			else:
				eigen_vals[0, eigen_count] = line.strip().split()[0]
				occs[0, eigen_count] = line.strip().split()[1]
			eigen_count +=1	

	# Check for inconsitencies (e.g. remaining np.nan entries)
	if np.any ( np.isnan(eigen_vals) ):
		raise ValueError("""NaN value was detected while parsing eigenvalues in file {:}.
				  Check this is a Plato *.occ file""".format(inpPath))
	if np.any ( np.isnan(k_path) ):
		raise ValueError("""NaN value was deteced while parsing the k-path in file {:} 
				 (but eigenvalues appear consistent). Check this is a 
				 Plato *.occ file""".format(inpPath))

	# Return dictionary of values
	out_dict = {'k_path' : k_path,
		    'eigen_vals' : eigen_vals,
	        'occs':occs,
	        'k_weights':k_weights}
	return out_dict



