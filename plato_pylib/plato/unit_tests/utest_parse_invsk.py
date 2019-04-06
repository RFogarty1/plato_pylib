#!/usr/bin/python3

import itertools
import os
import sys
import unittest


import plato_pylib.plato.parse_inv_sk as tCode


class testInvSKParse(unittest.TestCase):
	def setUp(self):
		self.partialFileA = createPartialInvSKFileA()

	def testVsKnownHVals_ssOrbs(self):
		shellA, shellB = 0,0
		bondType = "sigma"
		expectedHVals = [ (5.437502,-0.189830),
		                  (12.196023,-0.000494) ]
		invSkParsedObj = tCode.parseInvSK(self.partialFileA)

		actualHVals = invSkParsedObj.getAllValsOrbPair("hVal", shellA, shellB, bondType=bondType)

		for expRow,actRow in itertools.zip_longest(expectedHVals,actualHVals):
			[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expRow,actRow)]

	def testVsKnowSVals_ppOrbs(self):
		shellA, shellB = 1,1
		bondType = "pi"
		expectedSVals = [ (5.437502,0.171492),
		                  (12.196023,0.000043) ]
		invSKParsedObj = tCode.parseInvSK(self.partialFileA)
		actualSVals = invSKParsedObj.getAllValsOrbPair("sVal",shellA,shellB,bondType=bondType)

		for expRow,actRow in itertools.zip_longest(expectedSVals,actualSVals):
			[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expRow,actRow)]

	def testVsKnownHVals_spOrbs(self):
		shellA, shellB = 0,1 
		bondType = "sigma"
		expectedHVals = [ (5.437502, 0.214655),
		                  (12.196023, 0.001242) ]

		invSKParsedObj = tCode.parseInvSK(self.partialFileA)
		actualHVals = invSKParsedObj.getAllValsOrbPair("hVal",shellA,shellB,bondType=bondType)

		for expRow,actRow in itertools.zip_longest(expectedHVals,actualHVals):
			[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expRow,actRow)]

class testInvSKParseWithScreenFunct(unittest.TestCase):
	def setUp(self):
		self.partialFileB = createPartialInvSKFileB()

	def testHValsSS(self):
		shellA, shellB = 0,0
		bondType = "sigma"
		expVal = [(12.184001, -0.000076)]

		parsedObj = tCode.parseInvSK(self.partialFileB)
		actVal = parsedObj.getAllValsOrbPair("hVal", shellA, shellB, bondType=bondType)

		for expRow, actRow in itertools.zip_longest(expVal,actVal):
			[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(expRow,actRow)]


class testInvSKParseNewerNaNPadded(unittest.TestCase):
	def setUp(self):
		self.testFile=createPartialInvSkFilePaddedNoScreen()
	def testVsKnownHVals_ppPiOrbs(self):
		shellA, shellB = 1,1
		bondType = "pi"
		expVal = [(4.0, -0.232661)] #distance vs value, since thats the way the getter works most easily
	
		invSkParsedObj = tCode.parseInvSK(self.testFile)
		actVal = invSkParsedObj.getAllValsOrbPair("hVal",shellA,shellB,bondType=bondType)

		for expRow, actRow in itertools.zip_longest(expVal,actVal):
			[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(expRow,actRow)]
		

class testInvSKParseAngDepScreenFunct(unittest.TestCase):
	def setUp(self):
		self.testFile = createPartialInvSkFileAngDepScreen()
	
	def testVsKnownHVals_ppPiOrbs(self):
		shellA, shellB = 1,1
		bondType = "pi"
		expVal = [(4.0, -0.232661)] #distance vs value, since thats the way the getter works most easily
	
		invSkParsedObj = tCode.parseInvSK(self.testFile)
		actVal = invSkParsedObj.getAllValsOrbPair("hVal",shellA,shellB,bondType=bondType)

		for expRow, actRow in itertools.zip_longest(expVal,actVal):
			[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(expRow,actRow)]

	def testParsedAngDepScreenPP(self):
		shellA, shellB = 1,1
		expVal = [(4.0,-0.123241), (4.0,0.200541)] # [sigma,pi]

		invSkParsedObj = tCode.parseInvSK(self.testFile)

		actVal = list()
		actVal.append( invSkParsedObj.getAllValsOrbPair("screenFunctAngDep", shellA, shellB, bondType="sigma")[0] )
		actVal.append( invSkParsedObj.getAllValsOrbPair("screenFunctAngDep", shellA, shellB, bondType="pi")[0] )

		for expRow, actRow in itertools.zip_longest(expVal,actVal):
			[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(expRow,actRow)]



class testRemoveXtalTermsFromParsedObj(unittest.TestCase):
	def setUp(self):
		fakeDistsA = [0.0, 1.0, 4.0, 0.0, 3.3]
		fakeFieldObjs = [tCode.InvSKField(dist=x,posA=[0.0,0.0,0.0], posB=[0.0,0.0,x]) for x in fakeDistsA]
		fakeFieldObjsNoZeroDists = [tCode.InvSKField(dist=x,posA=[0.0,0.0,0.0], posB=[0.0,0.0,x]) for x in fakeDistsA if x>1e-9]
		self.fakeParsedFileA = tCode.InvSKAllData(fakeFieldObjs)
		self.fakeParsedFileAWithoutXtal = tCode.InvSKAllData(fakeFieldObjsNoZeroDists)

	def testForObjA(self):
		self.fakeParsedFileA.removeXtalFieldTerms()
		self.assertTrue( self.fakeParsedFileA==self.fakeParsedFileAWithoutXtal )

#from a compressed hcp Mg
def createPartialInvSKFileA():
	fileName = "partialInvSK.csv"
	filePath = os.path.join(os.getcwd(), fileName)
	fileStr = "x1, y1, z1, x2, y2, z2, shell 1, shell 2, l1, l2, r, S error, H error, S sigma, H sigma, S pi, H pi, S delta, H delta\n0.000000, 0.000000, 0.000000, -0.000000, 3.151409, 4.431145, 0, 0, 0, 0, 5.437502, 0.000000, 0.000000, 0.238795, -0.189830\n0.000000, 0.000000, 0.000000, -0.000000, 3.151409, 4.431145, 0, 1, 0, 1, 5.437502, 0.000000, 0.000009, -0.355209, 0.214655\n0.000000, 0.000000, 0.000000, -0.000000, 3.151409, 4.431145, 1, 0, 1, 0, 5.437502, 0.000000, 0.000007, 0.355209, -0.214655\n0.000000, 0.000000, 0.000000, -0.000000, 3.151409, 4.431145, 1, 1, 1, 1, 5.437502, 0.000000, 0.000205, -0.445739, 0.155412, 0.171492, -0.109021\n0.000000, 0.000000, 0.000000, -2.729200, -11.029931, -4.431145, 0, 0, 0, 0, 12.196023, 0.000000, 0.000000, 0.000134, -0.000494\n0.000000, 0.000000, 0.000000, -2.729200, -11.029931, -4.431145, 0, 1, 0, 1, 12.196023, 0.000000, 0.000004, -0.000367, 0.001242\n0.000000, 0.000000, 0.000000, -2.729200, -11.029931, -4.431145, 1, 0, 1, 0, 12.196023, 0.000000, 0.000004, 0.000367, -0.001242\n0.000000, 0.000000, 0.000000, -2.729200, -11.029931, -4.431145, 1, 1, 1, 1, 12.196023, 0.000000, 0.000018, -0.000993, 0.003084, 0.000043, -0.000172\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createPartialInvSKFileB():
	fileName = "partialInvSkB.csv"
	filePath = os.path.join(os.getcwd(),fileName)
	fileStr = "x1, y1, z1, x2, y2, z2, shell 1, shell 2, l1, l2, r, Screen Funct, S error, H error, S sigma, H sigma, S pi, H pi, S delta, H delta\n0.000000, 0.000000, 0.000000, -4.061333, -5.743593, -9.948195, 0, 0, 0, 0, 12.184001, 0.999491, 0.000000, 0.000000, 0.000004, -0.000076\n0.000000, 0.000000, 0.000000, -4.061333, -5.743593, -9.948195, 0, 1, 0, 1, 12.184001, 0.999491, 0.000000, 0.000000, 0.000011, -0.000172\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createPartialInvSkFilePaddedNoScreen():
	fileName = "partialInvSk_NaNpaddedNoScreen"
	filePath = os.path.join(os.getcwd(),fileName)
	fileStr = "x1, y1, z1, x2, y2, z2, shell 1, shell 2, l1, l2, r, Screen Funct, S error, H error, S sigma, H sigma, S pi, H pi, S delta, H delta\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 0, 0, 0, 0, 4.000000, 0.237697, 0.000000, 0.000000, 0.277741, -0.451181, NaN, NaN, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 0, 1, 0, 1, 4.000000, 0.237697, 0.000000, 0.029459, 0.402767, -0.508958, NaN, NaN, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 1, 0, 1, 0, 4.000000, 0.237697, 0.000000, 0.031368, -0.402767, 0.567173, NaN, NaN, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 1, 1, 1, 1, 4.000000, 0.237697, 0.000000, 0.040934, -0.460012, 0.442668, 0.216593, -0.232661, NaN, NaN\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createPartialInvSkFileAngDepScreen():
	fileName = "partialInvSk_angDepScreenFunct.csv"
	filePath = os.path.join(os.getcwd(),fileName)
	fileStr = "x1, y1, z1, x2, y2, z2, shell 1, shell 2, l1, l2, r, Screen Funct, S error, H error, S sigma, H sigma, S pi, H pi, S delta, H delta, Screen Sigma, Screen Pi, Screen Delta\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 0, 0, 0, 0, 4.000000, 0.237697, 0.000000, 0.000000, 0.277741, -0.451181, NaN, NaN, NaN, NaN, 0.237697, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 0, 1, 0, 1, 4.000000, 0.237697, 0.000000, 0.029459, 0.402767, -0.508958, NaN, NaN, NaN, NaN, 0.061981, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 1, 0, 1, 0, 4.000000, 0.237697, 0.000000, 0.031368, -0.402767, 0.567173, NaN, NaN, NaN, NaN, 0.052475, NaN, NaN\n0.000000, 0.000000, 0.000000, 4.000000, 0.000000, 0.000000, 1, 1, 1, 1, 4.000000, 0.237697, 0.000000, 0.040934, -0.460012, 0.442668, 0.216593, -0.232661, NaN, NaN, -0.123241, 0.200541, NaN\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

if __name__ == '__main__':
	unittest.main()

