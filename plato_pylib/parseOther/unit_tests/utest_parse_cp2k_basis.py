
import copy
import unittest
import unittest.mock as mock

import plato_pylib.parseOther.parse_cp2k_basis as tCode



class TestSplitAllBasisSetsUp(unittest.TestCase):

	def setUp(self):
		self.fileAsListA = _createTestFileAsListA()
		self.expBasisListA = _createExpBasisAsListA()
		self.expBasisListB = _createExpBasisAsListB()


	def testSimpleFileSplitsCorrectly(self):
		expResult = [self.expBasisListA, self.expBasisListB]
		actResult = tCode._getAllBasisSetSectionsFromCP2KCleanedFileAsList(self.fileAsListA)
		self.assertEqual(expResult,actResult)


class TestParseCP2kSingleBasisSet(unittest.TestCase):

	def setUp(self):
		self.basisListA = _createExpBasisAsListA()
		self.expA = [0.18835110]
		self.expB = [0.25684646,0.41223507]
		self.coeffsA = [ [8.44145278, 4.38527172, 1.46989022] ]
		self.coeffsB =  [ [6.29193240, 10.71646800, 4.67136283],
		                  [-33.73264167, -80.11784971, -35.36689601], ]
		self.lValsA = [0,1,2]
		self.lValsB = [0,1,2]
		self.nValA = 3

		self.eleName = "Mg"
		self.basisName = "spd-2z-rc7pt0-r05pt5-1"

		self.createTestObjs()

	def createTestObjs(self):
		self.expExponentSetA = tCode.ExponentSetCP2K(self.expA, self.coeffsA, self.lValsA, self.nValA)
		self.expExponentSetB = tCode.ExponentSetCP2K(self.expB, self.coeffsB, self.lValsB, self.nValA)
		self.expExponentSets = [self.expExponentSetA, self.expExponentSetB]
		self.expBasisSetObj = tCode.BasisSetCP2K(self.eleName, self.basisName, self.expExponentSets)

	def runTestFunct(self):
		return tCode._parseSingleBasisFileAsList(self.basisListA)

	def testNameAndElementParsedCorrectly(self):
		expElement = "Mg"
		expName = ["spd-2z-rc7pt0-r05pt5-1"]
		actObj = self.runTestFunct()
		self.assertEqual(expElement,actObj.element)
		self.assertEqual(expName,actObj.basisNames)

	def testBasisExpansionAsExpected(self):
		actObj = self.runTestFunct()
		actExponentSets = actObj.exponentSets
		self.assertEqual(self.expExponentSets,actExponentSets)


class TestExponentSetClass(unittest.TestCase):

	def setUp(self):
		self.exponentsA = [1,2]
		self.coeffsA = [[1,1],[2,3]]
		self.lValsA = [0,1]
		self.nValA = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ExponentSetCP2K(self.exponentsA,self.coeffsA,self.lValsA, self.nValA)

	def testEqualObjectsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertFalse(objA is objB)
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_singleCoeffChange(self):
		objA = copy.deepcopy(self.testObjA)
		self.coeffsA[0][0] += 2
		self.createTestObjs() #Probably not really needed
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_singleLValChange(self):
		objA = copy.deepcopy(self.testObjA)
		self.lValsA[0] += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffLengthExpLists(self):
		objA = copy.deepcopy(self.testObjA)
		self.exponentsA.append( 4 )
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)


class TestBasisSetClass(unittest.TestCase):

	def setUp(self):
		self.exponentSetsA = ["fake_exp_sets"]
		self.elementA = "Zr"
		self.basisNamesA = ["nameA","nameB"] #The order shouldnt matter
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.BasisSetCP2K(self.elementA, self.basisNamesA, self.exponentSetsA)

	def testEqualObjsCompareEqual_sameExceptCopied(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertTrue( objA is not objB )
		self.assertEqual(objA,objB)

	def testEqualObjsCompareEqual_basisNameOrderChanged(self):
		objA = copy.deepcopy(self.testObjA)
		newBasisNames = [x for x in reversed(self.basisNamesA)]
		self.assertNotEqual(newBasisNames, self.basisNamesA)
		self.basisNamesA = newBasisNames
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffExpSets(self):
		objA = copy.deepcopy(self.testObjA)
		newExpSets = ["fake_new_set"]
		self.assertNotEqual(newExpSets, self.exponentSetsA)
		self.exponentSetsA = newExpSets
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffBasisName(self):
		objA = copy.deepcopy(self.testObjA)
		newBasisNames = ["totally_new_basis"]
		self.assertNotEqual(self.basisNamesA, newBasisNames)
		self.basisNamesA = newBasisNames
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)


def _createExpBasisAsListA():
	 return "Mg spd-2z-rc7pt0-r05pt5-1\n2\n 3 0 2 1 1 1 1\n     0.18835110     8.44145278     4.38527172     1.46989022\n 3 0 2 2 1 1 1\n     0.25684646     6.29193240     10.71646800     4.67136283\n     0.41223507     -33.73264167     -80.11784971     -35.36689601".strip().split("\n")

def _createExpBasisAsListB():
	return "Mg fake-basis-a\n2\n 3 0 1 2 1 1\n     0.16410043     6.42183555     3.61896728\n     0.19845384     -7.80636538     -3.70059388\n 3 0 1 2 1 1\n     0.20037409     3.47476137     5.55019829\n     0.38555230     -15.08683172     -59.20857370".strip().split("\n")

def _createTestFileAsListA():
	fileStr = "#/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/create_cp2k_basis_mg/repo/src/create_basis/double_zeta_spd/vary_rc_r0/rc_7pt0\nMg spd-2z-rc7pt0-r05pt5-1\n2\n 3 0 2 1 1 1 1\n     0.18835110     8.44145278     4.38527172     1.46989022\n 3 0 2 2 1 1 1\n     0.25684646     6.29193240     10.71646800     4.67136283\n     0.41223507     -33.73264167     -80.11784971     -35.36689601\n\n\n#Just a comment\nMg fake-basis-a\n2\n 3 0 1 2 1 1\n     0.16410043     6.42183555     3.61896728\n     0.19845384     -7.80636538     -3.70059388\n 3 0 1 2 1 1\n     0.20037409     3.47476137     5.55019829\n     0.38555230     -15.08683172     -59.20857370\n\n"
	fileAsList = fileStr.split("\n")
	return fileAsList



