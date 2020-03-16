
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


class TestParsedBasisFileEquality(unittest.TestCase):

	def setUp(self):
		self.inpPath = "fake_path_a"
		self.basisSets = ["fake_obj_a","fake_obj_b"]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ParsedBasisFileCP2K(self.inpPath,self.basisSets)

	def testEqualObjsCompareEqual_copiedObj(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffBasisSets(self):
		objA = copy.deepcopy(self.testObjA)
		newBasisSets = ["this_is_new"]
		self.assertNotEqual(newBasisSets, self.basisSets)
		self.basisSets = newBasisSets
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)


class TestParseCP2KBasisSet(unittest.TestCase):

	def setUp(self):
		self.expInpPath = "fake_inp_path"
		self.expBasisA = _createExpectedBasisSetA()
		self.expBasisB = _createExpectedBasisSetB()
		self.createTestObjs()

	def createTestObjs(self):
		self.expObjA = tCode.ParsedBasisFileCP2K(self.expInpPath, [self.expBasisA,self.expBasisB]) 

	@mock.patch("plato_pylib.parseOther.parse_cp2k_basis._readInpFileIntoIter")
	def testExpectedOutputGiven(self, mockedFileReader):
		mockedFileReader.side_effect = lambda x: _createTestFileAsListA()
		actObj = tCode.parseCP2KBasisFile(self.expInpPath)
		mockedFileReader.assert_called_once_with(self.expInpPath)
		self.assertEqual(self.expObjA,actObj)


#Below is data for a fake file with 2 Mg basis sets
def _createExpectedBasisSetA():
	#Basic info
	element = "Mg"
	basisNames = ["spd-2z-rc7pt0-r05pt5-1"]

	#The exponent sets
	exponentsA = [0.18835110]
	coeffsA = [ [8.44145278, 4.38527172, 1.46989022] ]
	lValsA = [0,1,2]
	nValA = 3
	exponentSetA = tCode.ExponentSetCP2K(exponentsA, coeffsA, lValsA, nValA)

	exponentsB = [0.25684646, 0.41223507]
	coeffsB = [ [6.29193240, 10.71646800, 4.67136283], [-33.73264167, -80.11784971, -35.36689601] ]
	lValsB = [0,1,2]
	nValB = 3
	exponentSetB = tCode.ExponentSetCP2K(exponentsB, coeffsB, lValsB, nValB)

	return tCode.BasisSetCP2K(element, basisNames, [exponentSetA,exponentSetB])


def _createExpectedBasisSetB():
	#Basic info
	element = "Mg"
	basisNames = ["fake-basis-a"]

	#The exponent sets
	exponentsA = [0.16410043, 0.19845384]
	coeffsA = [ [6.42183555, 3.61896728], [-7.80636538, -3.70059388] ]
	lValsA = [0,1]
	nValA = 3
	exponentSetA = tCode.ExponentSetCP2K(exponentsA, coeffsA, lValsA, nValA)

	exponentsB = [0.20037409, 0.38555230]
	coeffsB = [ [3.47476137, 5.55019829], [-15.08683172, -59.20857370] ]
	lValsB = [0,1]
	nValB = 3
	exponentSetB = tCode.ExponentSetCP2K(exponentsB, coeffsB, lValsB, nValB)

	return tCode.BasisSetCP2K(element, basisNames, [exponentSetA,exponentSetB])


def _createExpBasisAsListA():
	 return "Mg spd-2z-rc7pt0-r05pt5-1\n2\n 3 0 2 1 1 1 1\n     0.18835110     8.44145278     4.38527172     1.46989022\n 3 0 2 2 1 1 1\n     0.25684646     6.29193240     10.71646800     4.67136283\n     0.41223507     -33.73264167     -80.11784971     -35.36689601".strip().split("\n")

def _createExpBasisAsListB():
	return "Mg fake-basis-a\n2\n 3 0 1 2 1 1\n     0.16410043     6.42183555     3.61896728\n     0.19845384     -7.80636538     -3.70059388\n 3 0 1 2 1 1\n     0.20037409     3.47476137     5.55019829\n     0.38555230     -15.08683172     -59.20857370".strip().split("\n")

def _createTestFileAsListA():
	fileStr = "#/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/create_cp2k_basis_mg/repo/src/create_basis/double_zeta_spd/vary_rc_r0/rc_7pt0\nMg spd-2z-rc7pt0-r05pt5-1\n2\n 3 0 2 1 1 1 1\n     0.18835110     8.44145278     4.38527172     1.46989022\n 3 0 2 2 1 1 1\n     0.25684646     6.29193240     10.71646800     4.67136283\n     0.41223507     -33.73264167     -80.11784971     -35.36689601\n\n\n#Just a comment\nMg fake-basis-a\n2\n 3 0 1 2 1 1\n     0.16410043     6.42183555     3.61896728\n     0.19845384     -7.80636538     -3.70059388\n 3 0 1 2 1 1\n     0.20037409     3.47476137     5.55019829\n     0.38555230     -15.08683172     -59.20857370\n\n"
	fileAsList = fileStr.split("\n")
	return fileAsList



