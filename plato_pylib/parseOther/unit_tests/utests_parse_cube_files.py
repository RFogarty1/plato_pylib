

import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.parseOther.parse_cube_files as tCode

class TestParseCubeFiles(unittest.TestCase):

	def setUp(self):
		self.inpPath = "fake_path_a"
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = _loadFileAsListA()

	def _runTestFunct(self):
		return tCode.parseCubeFile(self.inpPath)

	def _loadExpDictA(self):
		outDict = dict()
		outDict["header_a"] = "-Quickstep-"
		outDict["header_b"] = "ELECTRON DENSITY"
		outDict["n_atoms"] = 3
		outDict["origin"] = [1.0,0.0,0.0]
		outDict["n_x"], outDict["n_y"], outDict["n_z"] = 2, 2, 3
		outDict["step_x"], outDict["step_y"], outDict["step_z"] = [0.25,0,0], [0,0.25,0], [0,0,0.3]
		outDict["atomic_numbers"] = [8,1,1]
		outDict["atomic_coords"] = [ [1,2,3], [4,5,6], [7,8,9] ]
		outDict["atomic_charges"] = [0.3, -2, 0.7]

		#Sort the 3-dimensional grid of data
		x0Grid = [ [7,8,9],
		           [4,6,2]  ]
		x1Grid = [ [4,1,3],
		           [2,3,8] ]

		outDict["data_grid"] =[ x0Grid, x1Grid ]
		return outDict


	@mock.patch("plato_pylib.parseOther.parse_cube_files._readFileIntoList")
	def testCubeFilesA(self, mockReadFileIntoList):
		mockReadFileIntoList.side_effect = lambda *args,**kwargs: self.fileAsListA

		expDict = self._loadExpDictA()
		actDict = self._runTestFunct()
		mockReadFileIntoList.assert_called_with(self.inpPath)
		self._checkExpAndActDictMatch(expDict, actDict)


	def _checkExpAndActDictMatch(self, expDict, actDict):
		directCmpAttrs = ["header_a", "header_b", "n_atoms", "n_x", "n_y", "n_z", "atomic_numbers"]
		oneDimFloatVectAttrs = ["step_x", "step_y", "step_z", "atomic_charges"]
		multiDimFloatVectAttrs = ["atomic_coords", "data_grid"]

		#Check things we can compare directly
		for attr in directCmpAttrs:
			expVal, actVal = expDict[attr], actDict[attr]
			self.assertEqual(expVal,actVal)

		#Vectors of floats
		for attr in oneDimFloatVectAttrs + multiDimFloatVectAttrs:
			expVal, actVal = np.array(expDict[attr]), np.array(actDict[attr])
			self.assertTrue( np.allclose(expVal,actVal) )



#NOTE: Im not 100% sure on the newlines in the grid data. But read-only implementations shouldnt care about that particularly
def _loadFileAsListA():
	outList = list()
	outList.append( "-Quickstep-")
	outList.append( " ELECTRON DENSITY" )
	outList.append( "3 1.0 0.0 0.0 " ) #NAtoms and origin co-ordinates
	outList.append( "2 0.25 0.0 0.0" ) #Number X, vector for one step along x [Bohr]
	outList.append( "2 0.00 0.25 0.0" ) #Number Y, vector for one step along y [Bohr]
	outList.append( "3 0.00 0.00 0.3" ) #Number Z, vector for one step along z [Bohr]
	outList.append( "8 0.3 1 2 3" ) #Atomic number, charge, [x,y,z] (Atom 1)
	outList.append( "1 -2 4 5 6") #Atomic number, charge, [x,y,z] (Atom 2)
	outList.append( "1 0.7 7 8 9") #Atomic number, charge, [x,y,z] (Atom 3)

	#The data
	outList.append( "7 8 9" ) #grid[0][0][0], grid[0][0][1], grid[0][0][2]
	outList.append( "4 6 2" ) #grid[0][1][0], grid[0][1][1], grid[0][1][2]
	outList.append( "4 1 3" ) #grid[1][0][0], grid[1][0][1], grid[1][0][2]
	outList.append( "2 3 8" ) #grid[1][1][0], grid[1][1][1], grid[1][1][2]

	return outList



