
import itertools as it
from ..shared import ucell_class as uCell


def getUnitCellObjectFromASEAtomsObject(aseAtomsObj):
	""" Get a plato_pylib UnitCell object from ASE Atoms object (essentially the equivalent ASE object). Note at time of writing im only planning to test thes on objects with pbcs, so be careful if using for anything else
	
	Args:
		aseAtomsObj: (Atoms object from ASE) 
			 
	Returns
		 outCell: (plato_pylib UnitCell object)
 
	"""

	fractCoords = aseAtomsObj.get_scaled_positions()
	symbols = aseAtomsObj.get_chemical_symbols()
	lattLengths = aseAtomsObj.get_cell_lengths_and_angles()[:3]
	lattAngles = aseAtomsObj.get_cell_lengths_and_angles()[3:]

	outCell = uCell.UnitCell( lattParams=lattLengths, lattAngles=lattAngles )
	outCell.fractCoords = [list(x)+[y] for x,y in it.zip_longest(fractCoords,symbols)]
	return outCell


