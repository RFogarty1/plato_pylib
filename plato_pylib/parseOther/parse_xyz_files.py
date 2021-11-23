
import os
import pathlib
import re

from ..shared import ucell_class as uCell

def parseXyzFile(xyzFile):
	""" Extract geometry from xyz file into a unit cell
	
	Args:
		xyzFile: (str) Path to the .xyz file
			 
	Returns
		 outCell: (plato_pylib UnitCell object) This contains the cartesian co-ordinates in outCell.cartCoords. Also contains a unit-cell, which will be cubic
 
	"""
	fileAsList = _loadFileIntoList(xyzFile)

	#We need to remove empty lines at the end, but possibly preserve an empty comment line
	fileAsList = fileAsList[:2] + [x for x in fileAsList[2:] if x.strip()!=""]


	return _parseStandardXyzFile(fileAsList)


def _parseStandardXyzFile(fileAsList):

	expNumbAtoms = int( fileAsList[0].strip() )

	allCoords = list()
	for line in fileAsList[2:]:
		atom, *coords = line.strip().split()
		currCoords = [float(x) for x in coords] + [atom]
		allCoords.append( currCoords )

	outCell =   uCell.UnitCell.fromLattVects( [ [50.0,  0.0 ,  0.0],
	                                            [ 0.0, 50.0 ,  0.0],
	                                            [ 0.0,  0.0 , 50.0] ] )
	outCell.cartCoords = allCoords

	actNumbAtoms = len(allCoords)

	assert expNumbAtoms==actNumbAtoms, "xyz file contains {} atoms but expected {}".format(actNumbAtoms, expNumbAtoms)

	return outCell



def parseExtendedXyzFile_singleGeom(inpPath):
	""" Parses the first geometry from an extended xyz file
	
	Args:
		inpPath: (str) Path to input *.exyz file
			 
	Returns
		outCell: (plato_pylib UnitCell object)

	Raises:
		AssertionError: If there are more lines in file than expected from number of atoms in 1st line (e.g. due to multiple geometries)
 
	"""
	fileAsList = _loadFileIntoList(inpPath)
	outCartCoords = _parseStandardXyzFile(fileAsList).cartCoords
	outCell = _parseExtendedXyzCommentLine(fileAsList[1])
	outCell.cartCoords = outCartCoords
	return outCell
	


def dumpExtendedXyzFile_singleGeom(outPath, outGeom):
	""" Writes an extended xyz file (which includes cell vectors) for a single geometry
	
	Args:
		outPath: (str) Path to write the exyz file to
		outGeom: (plato_pylib UnitCell object)
 
	"""
	#Initialize
	outCart = outGeom.cartCoords
	outList = [ str(len(outCart)) ]

	#Get the comment line
	lattVectStr = ""
	for lVect in outGeom.lattVects:
		lattVectStr += " ".join([str(x) for x in lVect]) + " "

	commentLine = "Lattice=\"{}\"".format(lattVectStr)
	commentLine += " Properties=species:S:1:pos:R:3"
	outList.append(commentLine)

	#Get the cart coords
	for coord in outCart:
		currStr = "{} {:.8f} {:.8f} {:.8f}".format( coord[-1], *coord[0:3] )
		outList.append(currStr)

	#Merge all into a str
	outStr = "\n".join(outList)

	#Dump
	outFolder = os.path.split(outPath)[0]
	pathlib.Path(outFolder).mkdir(parents=True, exist_ok=True)
	with open(outPath,"wt") as f:
		f.write(outStr)


#Taken from gen_basis_helpers traj_io
def _parseExtendedXyzCommentLine(commentLine):
	#1) Get the unit cell parameters
	pattern = "Lattice=\".*\""
	lattStr = re.findall(pattern, commentLine)[0].replace("Lattice=","").replace("\"","")
	lattVals = [float(x) for x in lattStr.split()]
	assert len(lattVals)==9
	lattVects = [ lattVals[:3], lattVals[3:6], lattVals[6:9] ]
	unitCell = uCell.UnitCell.fromLattVects(lattVects)
	return unitCell


def _loadFileIntoList(inpPath):
	with open(inpPath,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList

