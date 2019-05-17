#!/usr/bin/python3


import itertools as it
import numpy as np
import unittest

import plato_pylib.utils.dosplot as tCode

class TestGenDosData(unittest.TestCase):

	def setUp(self):
		self.testEigValsA = [ [0.5,0.3] ]
		self.testOccValsA = [ [0.8,1.3] ]
		self.minX = 0.0
		self.maxX = 0.2
		self.stepX = 0.1 - 1e-8 #Make sure we get 3 points exactly
		self.smearWidth = 0.1
		self.expXVals = [self.minX,self.minX+self.stepX,self.maxX]

	def testForTwoGauFuncts(self):
		expYVals = [0.0576259231, 0.7029532065, 3.181074206]
		expOutVals = [(x,y) for x,y in it.zip_longest(self.expXVals,expYVals)]
		actOutVals = tCode.genDosData(self.testEigValsA, self.testOccValsA, minX=self.minX, maxX=self.maxX,
		                              step=self.stepX, smearWidth=self.smearWidth)

		expArray, actArray = np.array(expOutVals), np.array(actOutVals)

		self.assertEqual(expArray.shape,actArray.shape)
		self.assertTrue( np.allclose(expArray,actArray) )

	def testForTwoGauFunctsWithKPtWeights(self):
		self.testEigValsA.append( ([0.7,0.1]) )
		self.testOccValsA.append( ([1.4,0.3]) )
		testKWeights = [0.4,0.6]
		expYVals = [0.4585976735, 0.9992774384, 1.707989475]

		expOutVals = [(x,y) for x,y in it.zip_longest(self.expXVals,expYVals)]
		actOutVals = tCode.genDosData(self.testEigValsA, self.testOccValsA, minX=self.minX, maxX=self.maxX,
		                              step=self.stepX, smearWidth=self.smearWidth, kptweights=testKWeights)

		expArray, actArray = np.array(expOutVals), np.array(actOutVals)

		self.assertEqual(expArray.shape,actArray.shape)
		self.assertTrue( np.allclose(expArray,actArray) )



if __name__ == '__main__':
	unittest.main()

