#!/usr/bin/env python3

import itertools as it
import unittest

import plato_pylib.shared.ucell_class as UCell
#import plato_pylib.utils.supercell as supCell
import plato_pylib.utils.defects as tCode

class testCreateVacancy(unittest.TestCase):

	def setUp(self):
		lattVectsA = [ [3.0,0.0,0.0], [0.0,3.0,0.0], [0.0,0.0,3.0] ]
		fractCoordsA = [[0.3,0.3,0.3,"Mg"], [0.4,0.4,0.4,"Mg"], [0.57,0.57,0.8,"Mg"]]
		self.uCellA = UCell.UnitCell.fromLattVects( lattVectsA, fractCoordsA )

	def testCentralMethod_UCellInterface(self):

		expFractCoords = list(self.uCellA.fractCoords)
		expFractCoords.pop(1) #Expect the middle fract co-ord to disapear

		tCode.makeVacancyUnitCell(self.uCellA ,method="central")
		actFractCoords = self.uCellA.fractCoords

		for exp,act in it.zip_longest(expFractCoords,actFractCoords):
			[self.assertAlmostEqual(x,y) for x,y in it.zip_longest(exp[:3],act[:3])]
			self.assertEqual(exp[-1],act[-1])


if __name__ == '__main__':
	unittest.main()

