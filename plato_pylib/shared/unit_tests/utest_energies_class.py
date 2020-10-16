
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.energies_class as tCode

class TestEnergiesClassProps(unittest.TestCase):

	def setUp(self):
		self.dftTotalElectronic = 20
		self.entropy = 4
		self.dispersion = 3
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"dftTotalElectronic": self.dftTotalElectronic,
		          "entropy": self.entropy, "dispersion":self.dispersion}
		self.testObjA = tCode.EnergyVals(**kwargs)

	def testElectronicMinusEntropy_bothSet(self):
		expVal = self.dftTotalElectronic - self.entropy
		actVal = self.testObjA.electronicMinusEntropy
		self.assertEqual(expVal, actVal)

	def testElectronicMinusEntropy_entropyNone(self):
		self.entropy = None
		self.createTestObjs()
		expVal = None
		actVal = self.testObjA.electronicMinusEntropy
		self.assertEqual(expVal, actVal)

	def testEletronicMinusHalfEntrop(self):
		expVal = self.dftTotalElectronic - (0.5*self.entropy)
		actVal = self.testObjA.electronicMinusHalfEntropy
		self.assertEqual(expVal,actVal)

class TestEnergiesClassEqualities(unittest.TestCase):

	def setUp(self):
		self.e0tot, self.e0coh = 1, 2
		self.e1, self.e2 = 3, 4
		self.entropy, self.tb2CohesiveFree = 5, 6
		self.dftTotalElectronic, self.castepTotalElectronic = 7, 8
		self.dispersion = 9
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"e0tot":self.e0tot, "e0coh":self.e0coh, "e1":self.e1, "e2":self.e2,
		          "entropy":self.entropy, "tb2CohesiveFree":self.tb2CohesiveFree,
		          "dftTotalElectronic":self.dftTotalElectronic, "castepTotalElectronic":self.castepTotalElectronic,
		          "dispersion":self.dispersion}
		self.testObjA = tCode.EnergyVals(**kwargs)

	def testTwoEqualButDiffCompareEqual(self):
		objA = self.testObjA
		objB = copy.deepcopy(self.testObjA)
		self.assertEqual(objA, objB)

	def testTwoUnequalCompareUnequal_testA(self):
		objA = copy.deepcopy(self.testObjA)
		self.dispersion += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

class TestEnergiesClassToAndFromDict(unittest.TestCase):

	def setUp(self):
		self.e0tot, self.e0coh = 1, 2
		self.e1, self.e2 = 3, 4
		self.entropy, self.tb2CohesiveFree = 6, 6
		self.dftTotalElectronic, self.castepTotalElectronic = 7, 8
		self.dispersion = 9
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"e0tot":self.e0tot, "e0coh":self.e0coh, "e1":self.e1, "e2":self.e2,
		          "entropy":self.entropy, "tb2CohesiveFree":self.tb2CohesiveFree,
		          "dftTotalElectronic":self.dftTotalElectronic, 
		          "dispersion":self.dispersion}
		self.testObjA = tCode.EnergyVals(**kwargs)


	def testExpectedExtraKeysPresentInOutDict(self):
		expWithEntropyMinus, expWithHalfEntropyMinus = self.dftTotalElectronic-self.entropy, self.dftTotalElectronic-(0.5*self.entropy)
		expExtraDict = {"electronicTotalE":self.dftTotalElectronic, "electronicMinusEntropy": expWithEntropyMinus,
		                "electronicMinusHalfEntropy":expWithHalfEntropyMinus}
		actFullDict = self.testObjA.toDict()
		actExtraDict = {k:actFullDict[k] for k in expExtraDict.keys()}
		self.assertEqual(expExtraDict, actExtraDict)

	def testToAndFromDictCompatible(self):
		tempDict = self.testObjA.toDict()
		tempObj = tCode.EnergyVals(**tempDict)
		self.assertEqual(self.testObjA, tempObj)





