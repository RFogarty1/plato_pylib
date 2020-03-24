
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




def _loadFileIntoList(inpPath):
	with open(inpPath,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList

