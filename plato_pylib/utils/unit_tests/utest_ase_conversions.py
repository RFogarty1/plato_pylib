
import itertools as it
import unittest
import unittest.mock as mock

import ase.atoms as aseAtoms
import plato_pylib.shared.ucell_class as uCell

import plato_pylib.utils.ase_conversions as tCode

class TestUCellFromAtomsObj(unittest.TestCase):

	def setUp(self):
		self.lattVectA = [-5.0, 0.0, 5.0]
		self.lattVectB = [ 0.0, 5.0, 5.0]
		self.lattVectC = [-5.0, 5.0, 0.0]
		self.fractCoordsA = [ [0.0 , 0.0 , 0.0 ],
		                      [0.25, 0.25, 0.25] ]
		self.elesA = ["Si","Si"]

		self.createTestObjs()

	def createTestObjs(self):
		lattVects = [self.lattVectA, self.lattVectB, self.lattVectC]
		aseKwargDict = {"symbols":self.elesA, "scaled_positions":self.fractCoordsA, "cell":lattVects, "pbc":True}
		uCellFractCoords = [ x+[y] for x,y in it.zip_longest(self.fractCoordsA,self.elesA)]		
		self.aseObjA = aseAtoms.Atoms( **aseKwargDict )
		self.uCellA = uCell.UnitCell.fromLattVects(lattVects, fractCoords=uCellFractCoords)

	def testExpCellFromAseObjA(self):
		expUCell = self.uCellA
		actUCell = tCode.getUnitCellObjectFromASEAtomsObject(self.aseObjA)
		self.assertEqual(expUCell,actUCell)


