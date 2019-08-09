#!/usr/bin/python3

import copy
import os
import sys
import unittest

import plato_pylib.plato.plato_paths as tCode

class TestModelFoldersClass(unittest.TestCase):

	def setUp(self):
		self.testTightBindingBasePath = os.path.abspath( os.path.join("fake","path") )
		self.expPlatoPath = "faker"

		self.testAbsPathA = os.path.abspath( os.path.join( self.testTightBindingBasePath, self.expPlatoPath ) )
		self.testObjA = tCode.PlatoModelFolders(dft2Path=self.testAbsPathA)

		self.startModuleFunct = copy.deepcopy(tCode.getTightBindingDataPath) #copy probably overkill
		tCode.getTightBindingDataPath = lambda: self.testTightBindingBasePath

	def tearDown(self):
		tCode.getTightBindingDataPath = self.startModuleFunct

	def testExpectedDft2PlatoPathSimpleInput(self):
		actPlatoPath = self.testObjA.dft2PlatoPath
		self.testObjA.dft2PlatoPath
		self.assertEqual(self.expPlatoPath, actPlatoPath)



class TestGetAbsPathForTightBindingData(unittest.TestCase):

	def setUp(self):
		self.startFolder = os.path.abspath(os.path.join("fake","folder"))
		self.startRcPathFunction = copy.deepcopy(tCode.getPlatoRcPath)
		tCode.getPlatoRcPath = lambda: self.startFolder

		self.tightBindingPart = os.path.join("Data","TightBinding")
		self.dftPart = os.path.join("Data","SCF_LCAO")

		self.testDataSetStr = os.path.join("model","is","fake")
		self.dTypeArg = None

	def tearDown(self):
		tCode.getPlatoRcPath = self.startRcPathFunction

	def runFunction(self):
		if self.dTypeArg is None:
			return tCode.getAbsolutePathForPlatoTightBindingDataSet(self.testDataSetStr)
		else:
			return tCode.getAbsolutePathForPlatoTightBindingDataSet(self.testDataSetStr, dtype=self.dTypeArg)

	def testExpectedForSimpleInputTightBinding(self):
		expectedPath = os.path.join(self.startFolder,self.tightBindingPart,self.testDataSetStr)
		actualPath = self.runFunction()
		self.assertEqual(expectedPath,actualPath)

	def testExpectedForSimpleInputDft(self):
		self.dTypeArg = "dft"
		expectedPath = os.path.join(self.startFolder,self.dftPart,self.testDataSetStr)
		actualPath = self.runFunction()
		self.assertEqual(expectedPath,actualPath)

if __name__ == '__main__':
	unittest.main()

