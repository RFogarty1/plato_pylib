#!/usr/bin/python3

import itertools as it
import numpy as np
import unittest

import plato_pylib.utils.elastic_consts as tCode
import plato_pylib.shared.ucell_class as UCell


@unittest.skip("")
class TestCalcElastic(unittest.TestCase):

	def setUp(self):
		self.crystTypeA = "hexagonal"
		strainParams = [-0.01,0.0,0.01]
		stressVals = [0.00013, 0.00021, 0.00024, 0.0003, 0.00033]
		self.strainStresssA = [ [[x,y] for x,y in zip(strainParams,[s,0,s])] for s in stressVals ] 
		self.expElasticA = {(1,1):1.6, (1,2):0.5, (1,3):0.15, (3,3):2.6, (4,4):1.5}

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
		
	def testXyShearStrainA_fromList(self):
		strainMatrix = tCode._STRAIN_MATRIX_DICT[6](2*self.strainParamA)
		expLattVects = np.array(self.expectedShearStrainedA) 
		actLattVects = list(self.lattVectsA)
		tCode.applyStrainToLattVects(actLattVects, strainMatrix)
		self.assertTrue( np.allclose ( expLattVects , np.array(actLattVects) ) ) 




if __name__ == '__main__':
	unittest.main()



