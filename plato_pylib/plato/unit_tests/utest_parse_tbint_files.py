#!/usr/bin/python3

import itertools
import os
import sys
import unittest
import numpy as np


sys.path.append('../..')
import plato_pylib.plato.parse_tbint_files as tCode
from utest_parse_bas_files import createPartialBasFileA

import plato_pylib.plato.private.tbint_test_data as tData

DOUBLEABSTOL = 1e-7


class TestParseTbintFiles(unittest.TestCase):
	AllCreatedFilePaths = list()

	def setUp(self):
		self.AllCreatedFilePaths = list()
		self.AllCreatedFilePaths.append( tData.createTestModelFileA() )
		self.AllCreatedFilePaths.append( tData.createTestBdtFileA() )
		self.AllCreatedFilePaths.append( tData.createTestAdtFileA() )
		basFile = createPartialBasFileA()
		newBasName = os.path.join(os.getcwd(),"Mg.bas")
		os.rename(basFile, os.path.join(os.getcwd(),"Mg.bas"))
		self.AllCreatedFilePaths.append(newBasName)

	def tearDown(self):
		for currFilePath in self.AllCreatedFilePaths:
			os.remove(currFilePath)

	def testParseModelFile(self):
		expectedVals = {"model":1,
		                "pairFunct":0,
		                "overlap":1,
		                "crystalField":1,
		                "threeCentre":0,
		                "nijInts":1,
		                "snInts":1}

		inpModelFile = "model.dat"
		actualVals = tCode.parseModelFile("model.dat")
		for key in expectedVals.keys():
			self.assertEqual(expectedVals[key], actualVals[key])


	def testParseAdtFile(self):
		expectedshellIdxToAngMom = {0:0, 1:1}

		expectedVals = {"symbol":"Mg",
		                "coreCharge":2,
		                "numbShells":2,
		                "numbOrbitals":4,
		                "orbRadius":7.5,
		                "atomEnergy":-1.090379,
						"shellToAngMom":expectedshellIdxToAngMom}

		inpAdtFile = "Mg.adt"
		actualVals = tCode.parseAdtFile(inpAdtFile)

		for key in expectedVals:
			#NOTE: This tests for actual equality first; hence will only throw
			#      an error related to type (e.g. not a float) if that initial test fails
			self.assertAlmostEqual(expectedVals[key],actualVals[key]) 

	def testParseBdtFile(self):

		expectedVals = tData.loadTestBdtFileAExpectedVals()
		inpBdtFile = "Mg_Mg.bdt"
		shellIdxToAngMom = {0:0, 1:1}
		modelInfo = {"model":1,
		             "pairFunct":0,
		             "overlap":1,
		             "crystalField":1,
		             "threeCentre":0,
		             "nijInts":1,
		             "snInts":1,
		             "hop3B2C":0}

		actualVals = tCode.parseBdtFile(inpBdtFile, shellIdxToAngMom, shellIdxToAngMom, modelInfo)
		#Note: Assumes order is the same between expected/actual results
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

	def testParseWithOverallInterfaceFunct(self):
		expectedVals = tData.loadTestBdtFileAExpectedVals()
		expectedVals["nonLocPP"] = tData.loadTestBdtFileAExpectedVals_nlPPOnly()
		actualVals = tCode.getIntegralsFromBdt("Mg_Mg.bdt")

		#Note: Assumes order in each list (NOT DICT) is the same between expected/actual results
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

	def testGetAdtFilePathsFromBdt(self):
		fakeBaseFolder = os.path.join("fake","dir")
		testFilePath = os.path.join(fakeBaseFolder,"Mg_S.bdt")
		expectedFilePaths = [os.path.join(fakeBaseFolder,"Mg.adt"), os.path.join(fakeBaseFolder,"S.adt")]

		actualFilePaths = ["",""]
		actualFilePaths[0], actualFilePaths[1] = tCode.getAdtFilePathsFromBdt(testFilePath)

		self.assertEqual(expectedFilePaths,actualFilePaths)

	def testGetModelFilePathFromBdt(self):
		fakeBaseFolder = os.path.join("fake","dir")
		expectedFilePath = os.path.join(fakeBaseFolder, "model.dat")
		
		testFilePath = os.path.join(fakeBaseFolder,"Mg_S.bdt")
		actualFilePath = tCode.getModelFilePathFromBdt(testFilePath)

		self.assertEqual(expectedFilePath, actualFilePath)

	def testGetAtomNamesFromInpBdtFile(self):
		fakeBaseFolder = os.path.join("fake","dir")
		expectedAtomNames = ["Mg", "S"]
		
		actualAtomNames = ["",""]
		testFilePath = os.path.join(fakeBaseFolder,"Mg_S.bdt")
		actualAtomNames[0], actualAtomNames[1] = tCode.getAtomNamesFromInpBdtFile(testFilePath)

		self.assertEqual(expectedAtomNames, actualAtomNames)

	
	def checkNumpyArraysEqual(self, arrayA, arrayB): #DEPECRATED.
		self.assertEqual(arrayA.shape, arrayB.shape)
		self.assertEqual(arrayB.shape, arrayB.shape)
		if arrayA.shape == arrayB.shape:
			self.assertTrue( np.allclose(arrayA, arrayB, atol=DOUBLEABSTOL) )



class TestParseWithPairFunctNoOverlap(unittest.TestCase):
	def setUp(self):
		self.modFile = tData.createTestModelPairFunct_NoOverlapA()
		self.adt = tData.createTestAdtFile_PairFunctNoOverlapA()
		self.bdt = tData.createTestBdtFile_PairFunctNoOverlapA()
		self.basFile = createPartialBasFileA()

	def tearDown(self):
		os.remove(self.modFile)
		os.remove(self.adt)
		os.remove(self.bdt)
		os.remove(self.basFile)

	def testVsExpected(self):
		expInts = self.loadExpectedBdtResults()
		actualInts = tCode.getIntegralsFromBdt(self.bdt, inclPP=False)
		for key in expInts:
			self.assertTrue( actualInts[key] == expInts[key] )


	def loadExpectedBdtResults(self):
		outDict = dict()
		expPath = self.bdt
		shellA, shellB = 0,0
		angMomA, angMomB = 0,0
		atomAName, atomBName = "Xc", "Xc"

		pairPotInts = np.array(( [0,11962.26124], [0.01889725989,11430.03085], [0.03779451977,10921.29474] ))
		pairFunctInts = np.array(( [0,95.4944988], [0.01889725989,94.26066582], [0.03779451977,93.04277452] ))
		hopInts = np.array(( [0,1],[0.01889725989,2],[0.03779451977,3] ))

		outDict["pairFunct"] = [tCode.TbintIntegrals(inpFilePath=expPath, atomAName=atomAName, atomBName=atomBName, integrals=pairFunctInts)]
		outDict["pairPot"] = [tCode.TbintIntegrals(inpFilePath=expPath, atomAName=atomAName, atomBName=atomBName, integrals=pairPotInts)]
		outDict["hopping"] = [tCode.TbintIntegrals(inpFilePath=expPath, atomAName=atomAName, atomBName=atomBName, shellA=shellA, shellB=shellB, angMomA=angMomA, angMomB=angMomB, orbSubIdx=1, integrals= hopInts)]
		return outDict



class TestParseFormat4BdtFile(unittest.TestCase):
	def setUp(self):
		self.bdtFile = tData.createFormat4BdtFile_setAData()
	def tearDown(self):
		os.remove(self.bdtFile)

	def testParseBdtForm4_nonInterfaceFunct(self):
		expectedVals = tData.loadTestBdtFileAExpectedVals_format4()
		actualVals = tCode.parseBdtForm4(self.bdtFile)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

	def testParseBdtForm4_usingInterfaceFunct(self):
		expectedVals = tData.loadTestBdtFileAExpectedVals_format4()
		actualVals = tCode.getIntegralsFromBdt(self.bdtFile)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

class TestWriteOutputBdtFiles(unittest.TestCase):
	def setUp(self):
		self.bdtFile = tData.createFormat4BdtFile_setAData() #I'm only doing this to get the filepath really
		os.remove(self.bdtFile) #Need to delete the file, otherwise some tests will pass as long as its not overwritten 
	def tearDown(self):
		os.remove(self.bdtFile)
	def testWriteFormat4File(self):
		outputPath = self.bdtFile #Needs to be the same or parsed atom identities are wrong
		expectedVals = tData.loadTestBdtFileAExpectedVals_format4()
		tCode.writeBdtFileFormat4(expectedVals, outputPath)	
		actualVals = tCode.parseBdtForm4(outputPath)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )



class TestParseTBIntWithKineticAndHop3B2C(unittest.TestCase):

	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["modelFileA"] = tData.createModelFileB()
		self.filePathDict["adtFileA"] = tData.createTestAdtFileWithHop3B2CEigVals()
		self.filePathDict["bdtFileA"] = tData.createTestBdtFile_KineticHop3B2C_A()
	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

#	def testParseAdtEigVals(self):
#		expectedEigVals = [-3.917606593,22.49435857,22.49435857,22.49435857]
#		parsedDict = tCode.parseAdtFile(self.filePathDict["adtFileA"])
#		for exp,val in itertools.zip_longest(expectedEigVals, parsedDict["eigVals"]):
#			self.assertAlmostEqual(exp,val)

	def testParseKineticHoppingOverlap(self):
		testedInts = ["kinetic","overlap","hopping","hop3B2C"]
		expectedInts = tData.loadExpIntsBdtFile_KineticHop3B2C_A()
		actualInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"],inclPP=False)

		for key in testedInts:
			self.assertTrue( expectedInts[key] == actualInts[key] )


class testModIntegralsBdt(unittest.TestCase):
	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["modFileA"] = tData.createTestModelFileA()
		self.filePathDict["adtFileA"] = tData.createTestAdtFileA()
		self.filePathDict["bdtFileA"] = tData.createTestBdtFileA()
		basFile = createPartialBasFileA()
		newBasName = os.path.join(os.getcwd(),"Mg.bas")
		os.rename(basFile, os.path.join(os.getcwd(),"Mg.bas"))
		self.filePathDict["basFileA"] = newBasName


	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]
	
	def testVsKnownInput_PairPot(self):
		''' Test bdt file correctly modified for pair pot '''
		allInitInts = tData.loadTestBdtFileAExpectedVals()
		fakePairPotInts = np.array(( [0.0,4.0], [1.0,3.0], [2.0,2.0] ))
		allInitInts["pairPot"][0].integrals = fakePairPotInts
		allInitInts["pairPot"][0].intType = "pairPot"

		allInitInts["pairPot"][0].replaceIntsTbintFile( self.filePathDict["bdtFileA"] )
		newInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"], inclPP=False)

		for key in allInitInts.keys():
			self.assertTrue( allInitInts[key] == newInts[key] )


	def testVsKnownInput_HopPPpi(self):
		testShellA = 1
		testShellB = 1
		testSubIdx = 2

		allInitInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"])
		fakeInts = np.array(( [0.0,2.0], [6.0,3.0], [12.0,5.0] ))
		for intList in allInitInts.values():
			if intList is not None:
				for intSet in intList:
					if intSet.intType=="hopping" and  intSet.shellA == testShellA and intSet.shellB==testShellB and intSet.orbSubIdx==testSubIdx:
						intSet.integrals = fakeInts
						intSet.replaceIntsTbintFile( self.filePathDict["bdtFileA"] )

		newInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"])

		for key in newInts.keys():
			self.assertTrue( allInitInts[key] == newInts[key] )


	def testVsKnownInput_spSigma(self):
		testShellA = 0
		testShellB = 1
		testSubIdx = 1

		allInitInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"])
		fakeInts = np.array(( [0.0,3.0], [6.0,4.0], [12.0,9.2] ))

		for intList in allInitInts.values():
			if intList is not None:
				for intSet in intList:
					if intSet.intType=="hopping" and intSet.shellA==testShellA and intSet.shellB==testShellB and intSet.orbSubIdx==testSubIdx:
						intSet.integrals = fakeInts
						intSet.replaceIntsTbintFile( self.filePathDict["bdtFileA"] )

		newInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"])

		for key in newInts.keys():
			self.assertTrue( allInitInts[key] == newInts[key] )


if __name__ == '__main__':
	unittest.main()
