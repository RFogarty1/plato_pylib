
import itertools as it

import numpy as np

def parseCubeFile(inpPath):
	""" Parses the input cube file. Specification for a cube file format is given here:

	https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

	Note: This will only work for a subset of cube files.

	
	Args:
		inpPath: (str) Path to the input *.cube file
			 
	Returns
		outDict: (dict) All info in a dict format.
 
	"""
	fileAsList = _readFileIntoList(inpPath)
	outDict = dict()

	#Parse the header
	outDict["header_a"], outDict["header_b"] = fileAsList[0].strip(), fileAsList[1].strip()
	outDict["n_atoms"] = int( fileAsList[2].strip().split()[0] )
	outDict["origin"] = [ float(x) for x in fileAsList[2].strip().split()[1:] ] 

	outDict["n_x"] = int( fileAsList[3].strip().split()[0] )
	outDict["n_y"] = int( fileAsList[4].strip().split()[0] )
	outDict["n_z"] = int( fileAsList[5].strip().split()[0] )

	outDict["step_x"] = [float(x) for x in fileAsList[3].strip().split()[1:]]
	outDict["step_y"] = [float(x) for x in fileAsList[4].strip().split()[1:]]
	outDict["step_z"] = [float(x) for x in fileAsList[5].strip().split()[1:]]

	outDict["atomic_numbers"], outDict["atomic_coords"], outDict["atomic_charges"] = list(), list(), list()
	startGeomLine, endGeomLine = 6, 6+outDict["n_atoms"]
	geomLines = fileAsList[startGeomLine:endGeomLine]
	for line in geomLines:
		splitLine = line.strip().split()
		outDict["atomic_numbers"].append( int(splitLine[0]) )
		outDict["atomic_charges"].append( float(splitLine[1]) )
		outDict["atomic_coords"].append( [float(x) for x in splitLine[2:]] ) 

	#Parse the data grid part
	#Step 1 = get all the numbers in order
	finalLineIdx = len(fileAsList)
	gridNumbersOneDim = list()

	for line in fileAsList[endGeomLine:finalLineIdx]:
		currVals = [float(x) for x in line.strip().split()]
		gridNumbersOneDim.extend(currVals)

	#Step 2 = put them in a grid
	outGrid = np.zeros( (outDict["n_x"],outDict["n_y"],outDict["n_z"]) )

	totalIdx = 0
	for xIdx in range(outDict["n_x"]):
		for yIdx in range(outDict["n_y"]):
			for zIdx in range(outDict["n_z"]):
				outGrid[xIdx][yIdx][zIdx] = gridNumbersOneDim[totalIdx]
				totalIdx += 1
	outDict["data_grid"] = outGrid.tolist() 

	return outDict


def _readFileIntoList(inpPath):
	with open(inpPath,"rt") as f:
		outList = f.readlines()
	return outList







