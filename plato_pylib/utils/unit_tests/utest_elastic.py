#!/usr/bin/python3

import itertools as it
import numpy as np
import unittest

import plato_pylib.utils.elastic_consts as tCode
import plato_pylib.shared.ucell_class as UCell


class TestCalcElastic(unittest.TestCase):

	def setUp(self):
		self.crystTypeA = "hexagonal"
		strainParams = [-0.01,0.0,0.01]
		stressVals = [0.00013, 0.00021, 0.00024, 0.0003, 0.00033]
		self.strainStresssA = [ [[x,y] for x,y in zip(strainParams,[s,0,s])] for s in stressVals ] 
		self.expElasticA = {(1,1):2.7, (1,2):-0.6, (1,3):-0.5, (3,3):2.6, (4,4):0.75}

	def testCalcElasticHexagonalA(self):
		actElastic = tCode.calcElasticsFromStressStain(self.strainStresssA, self.crystTypeA).elasticConstants
		for key in self.expElasticA.keys():
			self.assertAlmostEqual( self.expElasticA[key], actElastic[key],places=5 )


class TestApplyStrain(unittest.TestCase):
	def setUp(self):
		self.lattVectsA = [ [5, 0, 0],
		                    [-3, 2, 0],
		                    [0, 0, 9] ]
		self.strainParamA = 0.5
		self.expectedShearStrainedA = [ [5, 2.5, 0],
		                                [-2, 0.5, 0],
		                                [0, 0, 9] ]

		self.startCartCoordsA = [ [0.0,0.0,0.0,"Mg"],
		                          [2.0, 2.0, 4.0, "Mg"] ]

		self.fractCoordsUCellB = [ [0.0,0.0,0.0,"Mg"],
		                           [0.33,0.33,0.5,"Mg"] ]

		self.uCellB = UCell.UnitCell.fromLattVects(self.lattVectsA, fractCoords = self.fractCoordsUCellB)

		
	def testXyShearStrainA_fromList_lattVects(self):
		strainMatrix = tCode._STRAIN_MATRIX_DICT[6](2*self.strainParamA)
		expLattVects = np.array(self.expectedShearStrainedA) 
		actLattVects = list(self.lattVectsA)
		tCode.applyStrainToLattVects(actLattVects, strainMatrix)
		self.assertTrue( np.allclose ( expLattVects , np.array(actLattVects) ) ) 

	def testXXStrain_ucellInterface_withAtomCoords(self):
		#Note: I dont expect the fractional co-ords to change with xx-strain. Meaning doing nothing to atomic-coords would pass this test.
		strainMatrix = tCode._STRAIN_MATRIX_DICT[1](1)
		expLattVects = [ [7.5, 0, 0],
		                 [-4.5, 2, 0],
		                 [0, 0, 9] ]
		expUCell = UCell.UnitCell.fromLattVects(expLattVects, fractCoords = self.fractCoordsUCellB) #The fractional co-ordinates shouldnt change for a non-shear strain
		actUCell = tCode.getStrainedUnitCellStructsForUnitStrainVects(self.uCellB, [self.strainParamA], [strainMatrix])[0][0]

		self.assertTrue(expUCell == actUCell)


	def testXYStrain_ucellInterface_withAtomCoords(self):
		strainMatrix = tCode._STRAIN_MATRIX_DICT[6](2) #xy strain; each ,matrix element is 0.5 by convention hence i set strain param to 2
		expLattVects = np.array(self.expectedShearStrainedA)
		expCartCoords = [[0.0,0.0,0.0,"Mg"],[0.99,0.99,4.5,"Mg"]]

		expUCell = UCell.UnitCell.fromLattVects(expLattVects)
		expUCell.cartCoords = expCartCoords

		actUCell = tCode.getStrainedUnitCellStructsForUnitStrainVects(self.uCellB, [self.strainParamA], [strainMatrix])[0][0]

		self.assertTrue(expUCell == actUCell)


if __name__ == '__main__':
	unittest.main()



