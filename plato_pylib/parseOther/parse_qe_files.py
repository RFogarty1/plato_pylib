
import io
import ase.io.espresso as aseQeParser

from ..utils import ase_conversions as aseConv
from ..shared import energies_class as energiesHelp

def parseQuantumEspressoOutfile(filePath):
	""" Parse contents of a quantum espresso output file
	
	Args:
		filePath: (str) Fulle path to the output file
			 
	Returns
		 parsedFileDict: (dict) Contains all information parsed from the file. Example keys below
 
	Example output keys:
		unitCell: (plato_pylib UnitCell object) The final geometry from the file. Units=bohr
		energies: (plato_pylib energies object) Contains the output energy of the final structure. Each attr is a different type of energy, e.g. "electronicTotalE". Units= eV
		numbAtoms: (int) Number of atoms in the calculation

	"""
	#Define the dicionary we return
	outKeys = ["unitCell", "energies", "numbAtoms"]
	outDict = {k:None for k in outKeys}

	#Get representations of the file
	fileStr = _getFileStrFromFilePath(filePath)
	fileAsList = [x.strip() for x in fileStr.split("\n") if x.strip()!=""]

	outDict["unitCell"] = _getUCellFromFileStr(fileStr)
	outDict.update( _parseQuantumEspressoFileAsList(fileAsList) )
	return outDict


def _getFileStrFromFilePath(inpFilePath):
	with open(inpFilePath,"rt") as f:
		outStr = f.read()
	return outStr


def _parseQuantumEspressoFileAsList(fileAsList):
	outKeys = ["energies"]
	outDict = {k:None for k in outKeys}

	for line in fileAsList:
		if "total energy              =" in line:
			totalEnergy = float( line.strip().split()[-2] )

	outDict["energies"] = energiesHelp.EnergyVals(dftTotalElectronic=totalEnergy)
	outDict["energies"].convRydToEv()

	return outDict




def _getUCellFromFileStr(fileStr):
	mockedFileObj = io.StringIO(fileStr)
	lastItemSlice = slice(-1,None,None)
	aseGeomObj = list(aseQeParser.read_espresso_out(mockedFileObj, index=lastItemSlice)) [-1]
	uCellObj = aseConv.getUnitCellObjectFromASEAtomsObject(aseGeomObj)
	uCellObj.convAngToBohr() #ASE uses angstroms
	return uCellObj
