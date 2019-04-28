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
		self.expElasticA = {(1,1):1.6, (1,2):0.5, (1,3):0.15, (3,3):2.6, (4,4):1.5}

	def testCalcElasticHexagonalA(self):
		actElastic = tCode.calcElasticsFromStressStain(self.strainStresssA, self.crystTypeA)
		for key in self.expElasticA.keys():
			self.assertAlmostEqual( self.expElasticA[key], actElastic[key],places=5 )

if __name__ == '__main__':
	unittest.main()



