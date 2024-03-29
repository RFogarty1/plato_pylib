
import os
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.parseOther.parse_xyz_files as tCode

class TestParseStandardXyzFile(unittest.TestCase):

	def setUp(self):
		self.testStrA = _loadTestFileStrA()
		self.fileAsListA =self.testStrA.split("\n")
		self.inpPathA = "fake_file_path.xyz"

	def _runTestFunct(self):
		return tCode.parseXyzFile(self.inpPathA)

	@mock.patch("plato_pylib.parseOther.parse_xyz_files.uCell.UnitCell")
	@mock.patch("plato_pylib.parseOther.parse_xyz_files._loadFileIntoList")
	def testCartCoordsFromSimpleFile(self, mockedFileOpener, mockedUnitCell):
		outUCell = mock.Mock()
		mockedFileOpener.side_effect = lambda *args,**kwargs : self.fileAsListA
		outUCell.fromLattVects.side_effect = lambda *args,**kwargs: outUCell
		
		expCartCoords = [ [0.47601,  1.17311, -0.01666, "Mg"],
		                  [2.49299,  1.29156, -0.08446, "O" ],
		                  [-1.52873, 1.05441,  0.04713, "O" ] ]

		parsedFile = self._runTestFunct()
		actCartCoords = parsedFile.cartCoords

		mockedFileOpener.assert_called_once_with(self.inpPathA)

		for expCoord,actCoord in it.zip_longest(expCartCoords, actCartCoords):
			self.assertTrue( len(expCoord)==len(actCoord) )
			for exp,act in zip( expCoord[:3], actCoord[:3] ):
				self.assertAlmostEqual(exp,act)
			self.assertEqual( expCoord[-1], actCoord[-1] )

	@mock.patch("plato_pylib.parseOther.parse_xyz_files._loadFileIntoList")
	def testRaisesIfNumbAtomsDoesntEqualFirstLine(self, mockedFileOpener):
		mockedFileOpener.side_effect = lambda *args,**kwargs : self.fileAsListA
		self.fileAsListA[0] = str( int(self.fileAsListA[0].strip()) + 1 )
		with self.assertRaises(AssertionError):
			self._runTestFunct()



def _loadTestFileStrA():
	outStr = "3\n	Energy:       1.1879058\nMg         0.47601        1.17311       -0.01666\nO          2.49299        1.29156       -0.08446\nO         -1.52873        1.05441        0.04713"
	return outStr

class TestDumpExyzForSingleUnitCell(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,11,12], [90,90,90]
		self.cartCoords = [ [0,1,1,"X"], [0,2,3,"Y"] ] 
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _dumpAndParse(self):
		tempFileName = "_temp_test_file.exyz"
		self._removeTempFile(tempFileName)
		tCode.dumpExtendedXyzFile_singleGeom(tempFileName, self.cellA)
		outGeom = tCode.parseExtendedXyzFile_singleGeom(tempFileName)
		self._removeTempFile(tempFileName)
		return outGeom

	def _removeTempFile(self, tempFileName):
		try:
			os.remove(tempFileName)
		except FileNotFoundError:
			pass

	def testDumpAndParseConsistentA(self):
		expCell = self.cellA
		actCell = self._dumpAndParse()
		self.assertEqual(expCell,actCell)




