#!/usr/bin/python3

import itertools as it
import numpy as np
import unittest

import plato_pylib.utils.elastic_consts as tCode
import plato_pylib.shared.ucell_class as UCell

class TestGetStrainMatrix(unittest.TestCase):

	def test33(self):
		strainCoeff = 2
		inpVals = (3,3)
		expMatrix = np.zeros( (3,3) )
		expMatrix[2,2] = strainCoeff
		actMatrix = tCode.getStrainMatrix(*inpVals, strainCoeff=strainCoeff)
		self.assertTrue(np.allclose(expMatrix,actMatrix))

	def test56(self):
		strainCoeff = 2
		inpVals = (5,6)
		expMatrix = np.array( [[0.0            , 0.5*strainCoeff, 0.5*strainCoeff],
		                       [0.5*strainCoeff, 0.0            , 0.0            ],
		                       [0.5*strainCoeff, 0.0            , 0.0            ]] )
		actMatrix = tCode.getStrainMatrix(*inpVals, strainCoeff=strainCoeff)
		self.assertTrue(np.allclose(expMatrix,actMatrix))


class TestApplyStrain(unittest.TestCase):
	def setUp(self):
		self.testVectsA = [ [1.0, 0.0, 0.0],
		                    [0.0, 1.4, 0.0],
		                    [0.0, 2.4, 4.6] ]

		self.outputVectsA_strain3_36 = [ [1.0, 2.1, 0.0],
		                                 [1.5, 1.4, 0.0],
		                                 [0.0, 9.6, 18.4] ]

		self.outputVectsA_strain3_11 = [ [4.0, 0.0, 0.0],
		                                 [0.0, 1.4, 0.0],
		                                 [0.0, 2.4, 4.6] ]

		fakeFractCoords = [[0.5,0.5,0.5,"Mg"], [0.75,0.75,0.75,"Mg"]]
		self.testUCellA = UCell.UnitCell.fromLattVects(self.testVectsA)
		self.testUCellA.fractCoords = fakeFractCoords


	def testFromListA_11(self):
		strainCoeff = 3
		n,m = 1,1
		tCode.applyStrainForElasticConstant(self.testVectsA, n, m, strainCoeff)
		self.assertTrue(np.allclose(self.outputVectsA_strain3_11, self.testVectsA))

	def testFromListA_36(self):
		strainCoeff = 3
		n,m = 3,6
		tCode.applyStrainForElasticConstant(self.testVectsA, n, m, strainCoeff)

		self.assertTrue(np.allclose(self.outputVectsA_strain3_36, self.testVectsA))

	def testFromNpArrayA_36(self):
		strainCoeff = 3
		n,m = 3,6
		inpVects = np.array(self.testVectsA)
		tCode.applyStrainForElasticConstant(inpVects, n, m, strainCoeff)
		self.assertTrue(np.allclose(self.outputVectsA_strain3_36, inpVects))

	def testInterfaceToUCellClass_11(self):
		strainCoeff = 3
		n,m = 1,1

		#We compare parameters/angles due to the semi-arbitrary nature of the lattice vectors
		expFractCoords = [list(x) for x in self.testUCellA.fractCoords]
		expAngles , expParams = self.testUCellA.getLattAnglesList(), self.testUCellA.getLattParamsList()
		expParams[0] = (strainCoeff*expParams[0]) + expParams[0] #ONLY this should change for the xx strain
		tCode.applyStrainToUCellForElasticConstant(self.testUCellA, n, m, strainCoeff)
		actAngles, actParams = self.testUCellA.getLattAnglesList(), self.testUCellA.getLattParamsList()
		actFractCoords = self.testUCellA.fractCoords

		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expParams,actParams)]
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expAngles,actAngles)]
		for expFract, actFract in it.zip_longest(expFractCoords,actFractCoords):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expFract,actFract)]
	
	def testInterfaceToUCellClass_36(self):
		''' Test correct elastic strain is applied to unit-cell class for calculating C36 constant '''
		strainCoeff = 3
		n,m = 3,6


		expFractCoords = [list(x) for x in self.testUCellA.fractCoords]
		expUCell = UCell.UnitCell.fromLattVects(self.outputVectsA_strain3_36)
		expAngles, expParams = expUCell.getLattAnglesList(), expUCell.getLattParamsList()

		tCode.applyStrainToUCellForElasticConstant(self.testUCellA, n, m, strainCoeff)
		actAngles, actParams = 	self.testUCellA.getLattAnglesList(), self.testUCellA.getLattParamsList()
		actFractCoords = self.testUCellA.fractCoords

		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expParams,actParams)]
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expAngles,actAngles)]
		for expFract, actFract in it.zip_longest(expFractCoords,actFractCoords):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expFract,actFract)]


if __name__ == '__main__':
	unittest.main()



