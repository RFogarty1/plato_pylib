
import unittest
import unittest.mock as mock

import plato_pylib.shared.energies_class as tCode

class TestEnergiesClass(unittest.TestCase):

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


