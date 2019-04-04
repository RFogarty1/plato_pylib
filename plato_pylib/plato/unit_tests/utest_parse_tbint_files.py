#!/usr/bin/python3

import itertools
import os
import sys
import unittest
import numpy as np


sys.path.append('../..')
import plato_pylib.plato.parse_tbint_files as tCode
from utest_parse_bas_files import createPartialBasFileA

DOUBLEABSTOL = 1e-7


class TestParseTbintFiles(unittest.TestCase):
	AllCreatedFilePaths = list()

	def setUp(self):
		self.AllCreatedFilePaths = list()
		self.AllCreatedFilePaths.append( createTestModelFileA() )
		self.AllCreatedFilePaths.append( createTestBdtFileA() )
		self.AllCreatedFilePaths.append( createTestAdtFileA() )
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

		expectedVals = loadTestBdtFileAExpectedVals()
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
		expectedVals = loadTestBdtFileAExpectedVals()
		expectedVals["nonLocPP"] = loadTestBdtFileAExpectedVals_nlPPOnly()
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
		self.modFile = createTestModelPairFunct_NoOverlapA()
		self.adt = createTestAdtFile_PairFunctNoOverlapA()
		self.bdt = createTestBdtFile_PairFunctNoOverlapA()
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
		self.bdtFile = createFormat4BdtFile_setAData()
	def tearDown(self):
		os.remove(self.bdtFile)

	def testParseBdtForm4_nonInterfaceFunct(self):
		expectedVals = loadTestBdtFileAExpectedVals_format4()
		actualVals = tCode.parseBdtForm4(self.bdtFile)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

	def testParseBdtForm4_usingInterfaceFunct(self):
		expectedVals = loadTestBdtFileAExpectedVals_format4()
		actualVals = tCode.getIntegralsFromBdt(self.bdtFile)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )

class TestWriteOutputBdtFiles(unittest.TestCase):
	def setUp(self):
		self.bdtFile = createFormat4BdtFile_setAData() #I'm only doing this to get the filepath really
		os.remove(self.bdtFile) #Need to delete the file, otherwise some tests will pass as long as its not overwritten 
	def tearDown(self):
		os.remove(self.bdtFile)
	def testWriteFormat4File(self):
		outputPath = self.bdtFile #Needs to be the same or parsed atom identities are wrong
		expectedVals = loadTestBdtFileAExpectedVals_format4()
		tCode.writeBdtFileFormat4(expectedVals, outputPath)	
		actualVals = tCode.parseBdtForm4(outputPath)
		for key in expectedVals:
			self.assertTrue( actualVals[key] == expectedVals[key] )



class TestParseTBIntWithKineticAndHop3B2C(unittest.TestCase):

	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["modelFileA"] = createModelFileB()
		self.filePathDict["adtFileA"] = createTestAdtFileWithHop3B2CEigVals()
		self.filePathDict["bdtFileA"] = createTestBdtFile_KineticHop3B2C_A()
	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

#	def testParseAdtEigVals(self):
#		expectedEigVals = [-3.917606593,22.49435857,22.49435857,22.49435857]
#		parsedDict = tCode.parseAdtFile(self.filePathDict["adtFileA"])
#		for exp,val in itertools.zip_longest(expectedEigVals, parsedDict["eigVals"]):
#			self.assertAlmostEqual(exp,val)

	def testParseKineticHoppingOverlap(self):
		testedInts = ["kinetic","overlap","hopping","hop3B2C"]
		expectedInts = loadExpIntsBdtFile_KineticHop3B2C_A()
		actualInts = tCode.getIntegralsFromBdt(self.filePathDict["bdtFileA"],inclPP=False)

		for key in testedInts:
			self.assertTrue( expectedInts[key] == actualInts[key] )


class testModIntegralsBdt(unittest.TestCase):
	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["modFileA"] = createTestModelFileA()
		self.filePathDict["adtFileA"] = createTestAdtFileA()
		self.filePathDict["bdtFileA"] = createTestBdtFileA()
		basFile = createPartialBasFileA()
		newBasName = os.path.join(os.getcwd(),"Mg.bas")
		os.rename(basFile, os.path.join(os.getcwd(),"Mg.bas"))
		self.filePathDict["basFileA"] = newBasName


	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]
	
	def testVsKnownInput_PairPot(self):
		''' Test bdt file correctly modified for pair pot '''
		allInitInts = loadTestBdtFileAExpectedVals()
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



def loadTestBdtFileAExpectedVals_format4():
	expectedVals = loadTestBdtFileAExpectedVals()
	for key in expectedVals:
		if expectedVals[key] is not None:
			for vals in expectedVals[key]:
				vals.atomAName = "Xa"
				vals.atomBName = "Xb"
				vals.inpFilePath = "Xa_Xb.bdt"

	return expectedVals
	

def loadTestBdtFileAExpectedVals():
	allBondIntegrals = dict()
	inpFilePath = "Mg_Mg.bdt"
	atomAName = "Mg"
	atomBName = "Mg"
	#Set overlap integrals
	overlapInts = list()

	overlapOrb00 = np.array(( [0.0, 0.9999999998], [6.0, 0.1716366773], [12.0, 0.0001915229276] ))
	overlapOrb01 = np.array(( [0.0, 0.0], [6.0, 0.2773564949], [12.0, 0.0005180942813] ))
	overlapOrb10 = np.array(( [0.0, 0.0], [6.0, -0.2773564949], [12.0, -0.0005180942813] ))
	overlapOrb11a = np.array(( [0.0, 1.0], [6.0, -0.3966863916], [12.0, -0.001386782855] ))
	overlapOrb11b = np.array(( [0.0, 1.0], [6.0, 0.1156877715], [12.0, 6.308189738e-05] ))

	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= overlapOrb00) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= overlapOrb01) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= overlapOrb10) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= overlapOrb11a) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= overlapOrb11b) )
	
	#Set Hopping Integrals
	hoppingInts = list()

	hoppingOrb00 = np.array(( [0.0, -0.4991995575], [6.0, -0.1101223541], [12.0, -0.0005213464141] ))
	hoppingOrb01 = np.array(( [0.0, -5.187901168e-16], [6.0, -0.126228632], [12.0, -0.001302054971] ))
	hoppingOrb10 = np.array(( [0.0, -5.159296592e-16], [6.0, 0.1262286317], [12.0, 0.001302054971] ))
	hoppingOrb11a = np.array(( [0.0, -0.1687854488], [6.0, 0.08381266327], [12.0, 0.003191656523] ))
	hoppingOrb11b = np.array(( [0.0, -0.1687854488], [6.0, -0.04564922265], [12.0, -0.0001780407443] ))

	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= hoppingOrb00) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= hoppingOrb01) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= hoppingOrb10) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= hoppingOrb11a) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= hoppingOrb11b) )

	#Set crystal field
	crystalNonXcInts = list()
	crystalXcInts = list()
	crystalFieldTotal = list()

	crystalNonXcOrb00 = np.array(( [0.0, -0.1365756498], [6.0, -0.00470605341], [12.0, -2.904952645e-09] ))
	crystalNonXcOrb01 = np.array(( [0.0, 0.0], [6.0, 0.008588245847], [12.0, 7.966734036e-09] ))
	crystalNonXcOrb10 = np.array(( [0.0, 0.0], [6.0, 0.008588245847], [12.0, 7.966734036e-09] ))
	crystalNonXcOrb11a = np.array(( [0.0, -0.1153456975], [6.0, -0.01611813044], [12.0, -2.197077218e-08] ))
	crystalNonXcOrb11b = np.array(( [0.0, -0.1153456975], [6.0, -0.001578543428], [12.0, -4.998967891e-10] ))

	crystalXcOrb00 = np.array(( [0.0, -0.1068463789], [6.0, -0.00815782228], [12.0, 3.38597624e-05] ))
	crystalXcOrb01 = np.array(( [0.0, -9.421579236e-17], [6.0, 0.01506649147], [12.0, -9.662080981e-05] ))
	crystalXcOrb10 = np.array(( [0.0, -9.416973494e-17], [6.0, 0.01506649147], [12.0, -9.662080981e-05] ))
	crystalXcOrb11a = np.array(( [0.0, -0.09604535093], [6.0, -0.02865339315], [12.0, 0.0002768477647] ))
	crystalXcOrb11b = np.array(( [0.0, -0.09604535093], [6.0, -0.003671575451], [12.0, 3.351807523e-05] ))

	crystalTotalOrb00 = np.array(crystalNonXcOrb00)
	crystalTotalOrb00[:,1] = crystalTotalOrb00[:,1] + crystalXcOrb00[:,1]

	crystalTotalOrb01 = np.array(crystalNonXcOrb01)
	crystalTotalOrb01[:,1] = crystalTotalOrb01[:,1] + crystalXcOrb01[:,1]

	crystalTotalOrb10 = np.array(crystalNonXcOrb10)
	crystalTotalOrb10[:,1] = crystalTotalOrb10[:,1] + crystalXcOrb10[:,1]

	crystalTotalOrb11a = np.array(crystalNonXcOrb11a)
	crystalTotalOrb11a[:,1] = crystalTotalOrb11a[:,1] + crystalXcOrb11a[:,1]

	crystalTotalOrb11b = np.array(crystalNonXcOrb11b)
	crystalTotalOrb11b[:,1] = crystalTotalOrb11b[:,1] + crystalXcOrb11b[:,1]

	crystalNonXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= crystalNonXcOrb00) )
	crystalNonXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= crystalNonXcOrb01) )
	crystalNonXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= crystalNonXcOrb10) )
	crystalNonXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= crystalNonXcOrb11a) )
	crystalNonXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= crystalNonXcOrb11b) )

	crystalXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= crystalXcOrb00) )
	crystalXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= crystalXcOrb01) )
	crystalXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= crystalXcOrb10) )
	crystalXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= crystalXcOrb11a) )
	crystalXcInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= crystalXcOrb11b) )

	crystalFieldTotal.append(  tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= crystalTotalOrb00) )
	crystalFieldTotal.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= crystalTotalOrb01) )
	crystalFieldTotal.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= crystalTotalOrb10) )
	crystalFieldTotal.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= crystalTotalOrb11a) )
	crystalFieldTotal.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= crystalTotalOrb11b) )

	#Set pair pot
	pairPotIntVals = np.array(( [1.0, 4.112339165], [7.0, -0.007517321281], [13.0, -4.787878002e-07] ))
	pairPot = [ tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, integrals=pairPotIntVals) ]

	#Set Nij Ints
	nijIntVals = np.array(( [0.0, 0.008730606763], [6.0, 0.0004110792408], [12.0, 1.363735317e-09] ))
	nijInts = [ tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, integrals=nijIntVals) ]

	#Set SN Ints
	snInts = list()

	snOrb00 = np.array(( [0.0, 0.008730996469], [6.0, 0.0004110655087], [12.0, 1.227666375e-09] ))
	snOrb01 = np.array(( [0.0, 0.0], [6.0, -0.0006938033489], [12.0, -3.324516136e-09] ))
	snOrb10 = np.array(( [0.0, 0.0], [6.0, -0.000693803349], [12.0, -3.324516136e-09] ))
	snOrb11a = np.array(( [0.0, 0.006955233076], [6.0, 0.001233863203], [12.0, 9.035170431e-09] ))
	snOrb11b = np.array(( [0.0, 0.006955233076], [6.0, 0.0001635580284], [12.0, 2.203649604e-10] ))

	snInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= snOrb00) )
	snInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= snOrb01) )
	snInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= snOrb10) )
	snInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= snOrb11a) )
	snInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= snOrb11b) )

	allBondIntegrals["overlap"] = overlapInts
	allBondIntegrals["hopping"] = hoppingInts
	allBondIntegrals["crystalFieldNonXc"] = crystalNonXcInts
	allBondIntegrals["crystalFieldXc"] = crystalXcInts
	allBondIntegrals["crystalFieldTotal"] = crystalFieldTotal
	allBondIntegrals["pairPot"] = pairPot
	allBondIntegrals["nijInts"] = nijInts
	allBondIntegrals["snInts"] = snInts
	allBondIntegrals["kinetic"] = None
	allBondIntegrals["hop3B2C"] = None
	allBondIntegrals["nonLocPP"] = None

	return allBondIntegrals


def loadTestBdtFileAExpectedVals_nlPPOnly():
	nlPPInts = list()
	inpFilePath = "Mg_Mg.bdt"
	atomAName = "Mg"
	atomBName = "Mg"
	
	nlPP_ssA = np.array((  [0.0, -0.443383718], [6.0, -0.02508443685], [12.0, -1.611075295e-09]  ))
	nlPP_ssB = np.array((  [0.0, -0.2069554578], [6.0, -0.01545985107 ], [12.0, 1.506580082e-10]  ))
	nlPP_sp =  np.array((  [0.0, 0], [6.0, 0.02542208488], [12.0, 6.128735821e-09] ))

	nlPP_psA = np.array(( [0.0, 0], [6.0, 0.06353056931], [12.0, 6.724517099e-09]  ))
	nlPP_psB = np.array(( [0.0, 0], [6.0, 0.04123458029], [12.0, -5.463310771e-10]  )) 
	nlPP_ppSigma = np.array(( [0.0, 0.3666159822], [6.0, -0.06026926636], [12.0, -2.455293202e-08]  )) 
	nlPP_ppPi = np.array(( [0.0, 0.3666159822], [6.0, 0.00816017747], [12.0, 5.185124571e-10]  )) 

	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= nlPP_ssA) )
	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=0, orbSubIdx=1, integrals= nlPP_ssB) )
	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=2, angMomA=0, angMomB=1, orbSubIdx=1, integrals= nlPP_sp) )

	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals=nlPP_psA) )
	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=0, orbSubIdx=1, integrals=nlPP_psB) )
	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=2, angMomA=1, angMomB=1, orbSubIdx=1, integrals=nlPP_ppSigma) )
	nlPPInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=2, angMomA=1, angMomB=1, orbSubIdx=2, integrals=nlPP_ppPi) )
	
	return nlPPInts


def loadExpIntsBdtFile_KineticHop3B2C_A():
	allBondIntegrals = dict()
	inpFilePath = "Mg_Mg.bdt"
	atomAName = "Mg"
	atomBName = "Mg"

	allBondIntegrals = dict()

	#Set overlap integrals
	overlapInts = list()

	overlapOrb00 = np.array(( [0.0, 0.999999992], [5.0, 0.298705153], [10.0, 0.003694447329], [15.0, 8.236051154e-08] ))
	overlapOrb01 = np.array(( [0.0, 0], [5.0, 0.4143005739], [10.0, 0.008969061077], [15.0, 2.876602768e-07] ))
	overlapOrb10 = np.array(( [0.0, 0], [5.0, -0.4143005739], [10.0, -0.008969061077], [15.0, -2.876602768e-07] ))
	overlapOrb11a = np.array(( [0.0, 1.000000003], [5.0, -0.4592864963], [10.0, -0.02129329409], [15.0, -9.910150161e-07] ))
	overlapOrb11b = np.array(( [0.0, 1.000000003], [5.0, 0.2249463523], [10.0, 0.001553751521], [15.0, 2.214349792e-08] ))

	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= overlapOrb00) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= overlapOrb01) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= overlapOrb10) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= overlapOrb11a) )
	overlapInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= overlapOrb11b) )

	#Set Hopping Integrals
	hoppingInts = list()

	hoppingOrb00 = np.array(( [0.0, -0.4988381778], [5.0, -0.1730445788], [10.0, -0.005594261718], [15.0, -6.86586574e-07] ))
	hoppingOrb01 = np.array(( [0.0, -5.23850093e-16], [5.0, -0.1627657523], [10.0, -0.01167103048], [15.0, -2.280812488e-06] ))
	hoppingOrb10 = np.array(( [0.0, -5.238311183e-16], [5.0, 0.1627657522], [10.0, 0.01167103048], [15.0, 2.280812488e-06] ))
	hoppingOrb11a = np.array(( [0.0, -0.1680618227], [5.0, 0.05890914432], [10.0, 0.02292371474], [15.0, 7.46059019e-06] ))
	hoppingOrb11b = np.array(( [0.0, -0.1680618227], [5.0, -0.06966393588], [10.0, -0.002334175974], [15.0, -1.809824886e-07] ))

	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= hoppingOrb00) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= hoppingOrb01) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= hoppingOrb10) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= hoppingOrb11a) )
	hoppingInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= hoppingOrb11b) )

	#Set Kinetic Integrals
	kineticInts = list()

	kineticOrb00 = np.array(( [0.0, 0.355720564], [5.0, 0.02246507682], [10.0, -0.004819182291], [15.0, -6.865849725e-07] ))
	kineticOrb01 = np.array(( [0.0, 0], [5.0, 0.1142679659], [10.0, -0.009780990376], [15.0, -2.280806806e-06] ))
	kineticOrb10 = np.array(( [0.0, 0], [5.0, -0.1142679658], [10.0, 0.009780990376], [15.0, 2.280806806e-06] ))
	kineticOrb11a = np.array(( [0.0, 0.5877085516], [5.0, -0.2768894815], [10.0, 0.01844809263], [15.0, 7.460571671e-06] ))
	kineticOrb11b = np.array(( [0.0, 0.5877085516], [5.0, 0.05599558099], [10.0, -0.00206491015], [15.0, -1.80982194e-07] ))

	kineticInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= kineticOrb00) )
	kineticInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= kineticOrb01) )
	kineticInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= kineticOrb10) )
	kineticInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= kineticOrb11a) )
	kineticInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= kineticOrb11b) )

	#Set Hop3B2C Integrals
	hop3B2CInts = list()

	hop3B2COrb00 = np.array(( [0.0, -0.6106894707], [5.0, -0.1377182742], [10.0, -0.001823848825], [15.0, -8.007372204e-13] ))
	hop3B2COrb01 = np.array(( [0.0, -4.305480713e-16], [5.0, -0.2150811139], [10.0, -0.003534695771], [15.0, -3.600359355e-12] ))
	hop3B2COrb10 = np.array(( [0.0, -4.326097098e-16], [5.0, 0.1744404505], [10.0, 0.005377167039], [15.0, 2.081136735e-12] ))
	hop3B2COrb11a = np.array(( [0.0, -0.5431737511], [5.0, 0.2275954553], [10.0, 0.009963562597], [15.0, 9.259682063e-12] ))
	hop3B2COrb11b = np.array(( [0.0, -0.5431737511], [5.0, -0.08964552501], [10.0, -0.001115323552], [15.0, -1.472760153e-13] ))

	hop3B2CInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=0, angMomA=0, angMomB=0, orbSubIdx=1, integrals= hop3B2COrb00) )
	hop3B2CInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=0, shellB=1, angMomA=0, angMomB=1, orbSubIdx=1, integrals= hop3B2COrb01) )
	hop3B2CInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=0, angMomA=1, angMomB=0, orbSubIdx=1, integrals= hop3B2COrb10) )
	hop3B2CInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=1, integrals= hop3B2COrb11a) )
	hop3B2CInts.append( tCode.TbintIntegrals(inpFilePath=inpFilePath, atomAName=atomAName, atomBName=atomBName, shellA=1, shellB=1, angMomA=1, angMomB=1, orbSubIdx=2, integrals= hop3B2COrb11b) )


	allBondIntegrals["overlap"] = overlapInts
	allBondIntegrals["hopping"] = hoppingInts
	allBondIntegrals["kinetic"] = kineticInts
	allBondIntegrals["hop3B2C"] = hop3B2CInts
	allBondIntegrals["nonLocPP"] = None


	return allBondIntegrals







def createTestModelFileA():
	filePath = os.path.join(os.getcwd(), "model.dat")
	with open(filePath,"wt") as f:
		f.write("# Model: Empirical potential = 0, Tabulated Tight Binding = 1, Gaussian Tight Binding = 2, Gaussian Moments = 3\n1\n# Pair functional\n0\n# Overlap\n1\n# Crystal field\n1\n# Three centre\n0\n# n_[ij} integrals (for many body xchange-correlation correction): 0=not present, 1=present\n1\n# Sankey (SN) integrals for treating many body xchange-correlation correction: 0=not present, 1=present \n1\n")
	return filePath

def createTestAdtFileA():
	filePath = os.path.join(os.getcwd(), "Mg.adt")
	with open(filePath,"wt") as f:
		f.write("Mg\n2.000000\n2\n4\n7.500000\n-1.090379\n3 0 0     -0.2102929769                 2               7.5\n3 1 0       0.121544206                 0               7.5\n3 1 1       0.121544206                 0               7.5\n3 1 2       0.121544206                 0               7.5\n\n     0.999998475917781     0.000000000000000     0.000000000000000     0.000000000000000\n     0.000000000000000     0.999994654553022     0.000000000000000     0.000000000000000\n     0.000000000000000     0.000000000000000     0.999994654553022     0.000000000000000\n     0.000000000000000     0.000000000000000     0.000000000000000     0.999994654553022\n\n    -0.255975428172679     0.000000000000000     0.000000000000000     0.000000000000000\n     0.000000000000000     0.042820105438336     0.000000000000000     0.000000000000000\n     0.000000000000000     0.000000000000000     0.042820105438336     0.000000000000000\n     0.000000000000000     0.000000000000000     0.000000000000000     0.042820105438336\n\n     0.500000000000000      0.050000000000000\n\n")
	return filePath



def createTestBdtFileA():
	filePath = os.path.join(os.getcwd(), "Mg_Mg.bdt")
	with open(filePath,"wt") as f:
		f.write("format_3\n0 0\n3\n                0      0.9999999998     -0.4991995575\n                6      0.1716366773     -0.1101223541\n               12   0.0001915229276  -0.0005213464141\n0 1\n3\n                0                 0  -5.187901168e-16\n                6      0.2773564949      -0.126228632\n               12   0.0005180942813   -0.001302054971\n1 0\n3\n                0                 0  -5.159296592e-16\n                6     -0.2773564949      0.1262286317\n               12  -0.0005180942813    0.001302054971\n1 1\n3\n                0                 1     -0.1687854488\n                6     -0.3966863916     0.08381266327\n               12   -0.001386782855    0.003191656523\n3\n                0                 1     -0.1687854488\n                6      0.1156877715    -0.04564922265\n               12   6.308189738e-05  -0.0001780407443\n0 0\n3\n                0     -0.1365756498     -0.1068463789\n                6    -0.00470605341    -0.00815782228\n               12  -2.904952645e-09    3.38597624e-05\n0 1\n3\n                0                 0  -9.421579236e-17\n                6    0.008588245847     0.01506649147\n               12   7.966734036e-09  -9.662080981e-05\n1 0\n3\n                0                 0  -9.416973494e-17\n                6    0.008588245847     0.01506649147\n               12   7.966734036e-09  -9.662080981e-05\n1 1\n3\n                0     -0.1153456975    -0.09604535093\n                6    -0.01611813044    -0.02865339315\n               12  -2.197077218e-08   0.0002768477647\n3\n                0     -0.1153456975    -0.09604535093\n                6   -0.001578543428   -0.003671575451\n               12  -4.998967891e-10   3.351807523e-05\n3\n                1       4.112339165\n                7   -0.007517321281\n               13  -4.787878002e-07\n3\n                0    0.008730606763\n                6   0.0004110792408\n               12   1.363735317e-09\n0 0\n3\n                0    0.008730996469\n                6   0.0004110655087\n               12   1.227666375e-09\n0 1\n3\n                0                 0\n                6  -0.0006938033489\n               12  -3.324516136e-09\n1 0\n3\n                0                 0\n                6   -0.000693803349\n               12  -3.324516136e-09\n1 1\n3\n                0    0.006955233076\n                6    0.001233863203\n               12   9.035170431e-09\n3\n                0    0.006955233076\n                6   0.0001635580284\n               12   2.203649604e-10\n0 0\n3\n                0      -0.443383718\n                6    -0.02508443685\n               12  -1.611075295e-09\n0 1\n3\n                0     -0.2069554578\n                6    -0.01545985107\n               12   1.506580082e-10\n0 2\n3\n                0                 0\n                6     0.02542208488\n               12   6.128735821e-09\n1 0\n3\n                0                 0\n                6     0.06353056931\n               12   6.724517099e-09\n1 1\n3\n                0                 0\n                6     0.04123458029\n               12  -5.463310771e-10\n1 2\n3\n                0      0.3666159822\n                6    -0.06026926636\n               12  -2.455293202e-08\n3\n                0      0.3666159822\n                6     0.00816017747\n               12   5.185124571e-10\n")

	return filePath




def createTestBdtFile_KineticHop3B2C_A():
	filePath = os.path.join(os.getcwd(), "Mg_Mg.bdt")
	fileStr = "format_3\n0 0\n4\n                0       0.999999992     -0.4988381778       0.355720564\n                5       0.298705153     -0.1730445788     0.02246507682\n               10    0.003694447329   -0.005594261718   -0.004819182291\n               15   8.236051154e-08   -6.86586574e-07  -6.865849725e-07\n0 1\n4\n                0                 0   -5.23850093e-16                 0\n                5      0.4143005739     -0.1627657523      0.1142679659\n               10    0.008969061077    -0.01167103048   -0.009780990376\n               15   2.876602768e-07  -2.280812488e-06  -2.280806806e-06\n1 0\n4\n                0                 0  -5.238311183e-16                 0\n                5     -0.4143005739      0.1627657522     -0.1142679658\n               10   -0.008969061077     0.01167103048    0.009780990376\n               15  -2.876602768e-07   2.280812488e-06   2.280806806e-06\n1 1\n4\n                0       1.000000003     -0.1680618227      0.5877085516\n                5     -0.4592864963     0.05890914432     -0.2768894815\n               10    -0.02129329409     0.02292371474     0.01844809263\n               15  -9.910150161e-07    7.46059019e-06   7.460571671e-06\n4\n                0       1.000000003     -0.1680618227      0.5877085516\n                5      0.2249463523    -0.06966393588     0.05599558099\n               10    0.001553751521   -0.002334175974    -0.00206491015\n               15   2.214349792e-08  -1.809824886e-07   -1.80982194e-07\n0 0\n4\n                0     -0.1366817755     -0.1071874956\n                5    -0.01428833251    -0.01910484966\n               10  -1.678389766e-06  -2.577920574e-05\n               15  -1.662308546e-13                 0\n0 1\n4\n                0                 0  -9.591238612e-17\n                5     0.02128862238     0.03044061564\n               10   4.105577806e-06   5.448751559e-05\n               15   4.747232819e-13                 0\n1 0\n4\n                0                 0  -9.626178721e-17\n                5     0.02128862238     0.03044061564\n               10   4.105577806e-06   5.448751559e-05\n               15   4.747232819e-13                 0\n1 1\n4\n                0     -0.1154622835    -0.09713433976\n                5    -0.03238608161    -0.05175981284\n               10  -1.020826742e-05  -0.0001158626389\n               15  -1.370809868e-12                 0\n4\n                0     -0.1154622835    -0.09713433976\n                5   -0.006531423293    -0.01085244035\n               10  -3.940497012e-07   2.166939183e-06\n               15  -3.353452506e-14                 0\n3\n                1       4.121963454\n                6    -0.01233929885\n               11    0.000267243982\n4\n                0    0.008729574959\n                5    0.001157711652\n               10   4.206289384e-07\n               15    1.43859688e-13\n0 0\n4\n                0    0.008731161449\n                5    0.001157952064\n               10   3.784674529e-07\n               15   2.356842078e-14\n0 1\n4\n                0                 0\n                5   -0.001703227799\n               10  -8.994691384e-07\n               15  -6.967370115e-14\n1 0\n4\n                0                 0\n                5   -0.001703227799\n               10  -8.994691384e-07\n               15  -6.967370115e-14\n1 1\n4\n                0    0.006956117659\n                5    0.002716976243\n               10   2.170399827e-06\n               15   2.081041112e-13\n4\n                0    0.006956117659\n                5   0.0005424796711\n               10   9.635563572e-08\n               15   4.370379516e-15\n0 0\n4\n                0     -0.6106894707\n                5     -0.1377182742\n               10   -0.001823848825\n               15  -8.007372204e-13\n0 1\n4\n                0  -4.305480713e-16\n                5     -0.2150811139\n               10   -0.003534695771\n               15  -3.600359355e-12\n1 0\n4\n                0  -4.326097098e-16\n                5      0.1744404505\n               10    0.005377167039\n               15   2.081136735e-12\n1 1\n4\n                0     -0.5431737511\n                5      0.2275954553\n               10    0.009963562597\n               15   9.259682063e-12\n4\n                0     -0.5431737511\n                5    -0.08964552501\n               10   -0.001115323552\n               15  -1.472760153e-13\n0 0\n4\n                0     -0.4433974117\n                5    -0.07766735975\n               10  -1.277958369e-06\n               15  -8.040769796e-17\n0 1\n4\n                0     -0.2069664966\n                5    -0.05410482182\n               10   4.631271198e-09\n               15   2.200698479e-17\n0 2\n4\n                0                 0\n                5     0.06025546141\n               10   3.682238894e-06\n               15    7.97559683e-16\n1 0\n4\n                0                 0\n                5      0.1762438117\n               10   5.015107154e-06\n               15   4.821171444e-16\n1 1\n4\n                0                 0\n                5      0.1296349285\n               10   7.979411692e-08\n               15  -1.292144709e-16\n1 2\n4\n                0       0.366651209\n                5     -0.1203356194\n               10  -1.378853524e-05\n               15  -4.549841863e-15\n4\n                0       0.366651209\n                5     0.02769635098\n               10    4.06097591e-07\n               15   5.566006014e-17\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createTestAdtFileWithHop3B2CEigVals():
	filePath = os.path.join(os.getcwd(),"Mg.adt")
	fileStr = "Mg\n2.000000\n2\n4\n7.500000\n-1.090324\n3 0 0     -0.2102929769                 2               7.5\n3 1 0       0.121544206                 0               7.5\n3 1 1       0.121544206                 0               7.5\n3 1 2       0.121544206                 0               7.5\n\n     0.999999396114689     0.000000000000000     0.000000000000000     0.000000000000000\n     0.000000000000000     0.999997503861023     0.000000000000000     0.000000000000000\n     0.000000000000000     0.000000000000000     0.999997503861023     0.000000000000000\n     0.000000000000000     0.000000000000000     0.000000000000000     0.999997503861023\n\n    -0.255257891819027     0.000000000000000     0.000000000000000     0.000000000000000\n     0.000000000000000     0.044455590806116     0.000000000000000     0.000000000000000\n     0.000000000000000     0.000000000000000     0.044455590806116     0.000000000000000\n     0.000000000000000     0.000000000000000     0.000000000000000     0.044455590806116\n\n     0.500000000000000      0.050000000000000\n\n     -3.917606593\n      22.49435857\n      22.49435857\n      22.49435857\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createModelFileB():
	filePath = os.path.join(os.getcwd(),"model.dat")
	fileStr = "# Model: Empirical potential = 0, Tabulated Tight Binding = 1, Gaussian Tight Binding = 2, Gaussian Moments = 3\n1\n# Pair functional\n0\n# Overlap\n1\n# Crystal field\n1\n# Three centre\n0\n# n_[ij} integrals (for many body xchange-correlation correction): 0=not present, 1=present\n1\n# Sankey (SN) integrals for treating many body xchange-correlation correction: 0=not present, 1=present\n1\n# Hop3B2C Flag (2 center integrals used for a 3-body correction of hopping terms)\n1\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createTestModelPairFunct_NoOverlapA():
	filePath = os.path.join(os.getcwd(),"model.dat")
	fileStr = "# Model: Empirical potential = 0, Tabulated Tight Binding = 1, Gaussian Tight Binding = 2, Gaussian Moments = 3\n1\n# Pair functional\n1\n# Overlap\n0\n# Crystal field\n0\n# Three centre\n0\n# n_[ij} integrals (for many body xchange-correlation correction): 0=not present, 1=present\n0\n# Sankey (SN) integrals for treating many body xchange-correlation correction: 0=not present, 1=present\n0\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createTestAdtFile_PairFunctNoOverlapA():
	filePath = os.path.join(os.getcwd(),"Xc.adt")
	fileStr = "V\n0.0\n1\n1\n6.141609\n0.0\n1 0 0 0.0 0.0 6.141609\n\n1.0\n\n0.0\n\n0.0 0.0\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createTestBdtFile_PairFunctNoOverlapA():
	filePath = os.path.join(os.getcwd(),"Xc_Xc.bdt")
	fileStr = "format_3\n0 0\n3\n                0                 1\n    0.01889725989                 2\n    0.03779451977                 3\n3\n                0       11962.26124\n    0.01889725989       11430.03085\n    0.03779451977       10921.29474\n3\n                0        95.4944988\n    0.01889725989       94.26066582\n    0.03779451977       93.04277452\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createFormat4BdtFile_setAData():
	filePath = os.path.join(os.getcwd(), "Xa_Xb.bdt")
	fileStr = "format_4\nOverlap\n#Number of Tables\n5\n#Orbital Shell Indices\n0 0\n#Orbital Angular Momentum\n0 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0.9999999998\n6	0.1716366773\n12	0.0001915229\n#Overlap\n#Orbital Shell Indices\n0 1\n#Orbital Angular Momentum\n0 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	0.2773564949\n12	0.0005180943\n#Overlap\n#Orbital Shell Indices\n1 0\n#Orbital Angular Momentum\n1 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	-0.2773564949\n12	-0.0005180943\n#Overlap\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	1\n6	-0.3966863916\n12	-0.0013867829\n#Overlap\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n1\n#Number of points in table\n3\n0	1\n6	0.1156877715\n12	6.308189738E-05\nHopping\n#Number of Tables\n5\n#Orbital Shell Indices\n0 0\n#Orbital Angular Momentum\n0 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.4991995575\n6	-0.1101223541\n12	-0.0005213464\n#Hopping\n#Orbital Shell Indices\n0 1\n#Orbital Angular Momentum\n0 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-5.19E-16\n6	-1.26228632E-01\n12	-1.30205497E-03\n#Hopping\n#Orbital Shell Indices\n1 0\n#Orbital Angular Momentum\n1 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-5.16E-16\n6.0	1.26228632E-01\n12.0	1.30205497E-03\n#Hopping\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.1687854488\n6	0.0838126633\n12	0.0031916565\n#Hopping\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n1\n#Number of points in table\n3\n0	-0.1687854488\n6	-0.04564922265\n12	-0.0001780407443\nxtalFieldTotal\n#Number of Tables\n5\n#Orbital Shell Indices\n0 0\n#Orbital Angular Momentum\n0 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.2434220287\n6	-0.0128638757\n12	3.3856857447355E-05\n#xtalFieldTotal\n#Orbital Shell Indices\n0 1\n#Orbital Angular Momentum\n0 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-9.421579236E-17\n6	0.0236547373\n12	-9.6612843075964E-05\n#xtalFieldTotal\n#Orbital Shell Indices\n1 0\n#Orbital Angular Momentum\n1 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-9.416973494E-17\n6	0.0236547373\n12	-9.6612843075964E-05\n#xtalFieldTotal\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.2113910484\n6	-0.0447715236\n12	0.0002768258\n#xtalFieldTotal\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n1\n#Number of points in table\n3\n0	-0.2113910484\n6	-0.0052501189\n12	3.35175753332109E-05\nxtalFieldNoXc\n#Number of Tables\n5\n#Orbital Shell Indices\n0 0\n#Orbital Angular Momentum\n0 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.1365756498\n6	-0.0047060534\n12	-2.90E-09\n#xtalFieldNoXc\n#Orbital Shell Indices\n0 1\n#Orbital Angular Momentum\n0 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	0.0085882458\n12	7.97E-09\n#xtalFieldNoXc\n#Orbital Shell Indices\n1 0\n#Orbital Angular Momentum\n1 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	0.0085882458\n12	7.97E-09\n#xtalFieldNoXc\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	-0.1153456975\n6	-0.0161181304\n12	-2.20E-08\n#xtalFieldNoXc\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n1\n#Number of points in table\n3\n0	-0.1153456975\n6	-0.0015785434\n12	-5.00E-10\nsnIntegrals\n#Number of Tables\n5\n#Orbital Shell Indices\n0 0\n#Orbital Angular Momentum\n0 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0.0087309965\n6	0.0004110655\n12	1.23E-09\n#snIntegrals\n#Orbital Shell Indices\n0 1\n#Orbital Angular Momentum\n0 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	-0.0006938033\n12	-3.32E-09\n#snIntegrals\n#Orbital Shell Indices\n1 0\n#Orbital Angular Momentum\n1 0\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0\n6	-0.0006938033\n12	-3.32E-09\n#snIntegrals\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n0\n#Number of points in table\n3\n0	0.0069552331\n6	0.0012338632\n12	9.04E-09\n#snIntegrals\n#Orbital Shell Indices\n1 1\n#Orbital Angular Momentum\n1 1\n#Axial Angular Momentum\n1\n#Number of points in table\n3\n0	0.0069552331\n6	0.000163558\n12	2.20E-10\nPairPot\n#Number of Tables\n1\n#Number of points in table\n3\n                1       4.112339165\n                7   -0.007517321281\n               13  -4.787878002e-07\nnijIntegrals\n#Number of Tables\n1\n#Number of points in table\n3\n                0    0.008730606763\n                6   0.0004110792408\n               12   1.363735317e-09\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

if __name__ == '__main__':
	unittest.main()
