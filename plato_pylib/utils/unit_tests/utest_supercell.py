#!/usr/bin/env python3

import copy
import unittest
import plato_pylib.utils.supercell as tCode
import plato_pylib.shared.ucell_class as UCell


class testCreateSuperCells(unittest.TestCase):

	def setUp(self):
		fractCoordsA = [ [0.0,0.0,0.0,"Mg"], [0.33333333, 0.33333333, 0.5, "Mg"] ]
		origLattVectsA = [ [6.06, 0.00       , 0.00],
		                   [3.03, 5.248113941, 0.00],
		                   [0.00, 0.00       , 9.839043529] ]

		self.startUCellA = UCell.UnitCell.fromLattVects(origLattVectsA, fractCoordsA)



	def testCellA_2x2x1(self):
		#These are values from platos build (the fractCoords were reordered)
		lattVects_221 = [ [12.12, 0.00         , 0.00],
		                  [ 6.06, 10.4962278819, 0.00],
		                  [ 0.00, 0.00         , 9.839043529] ]


		fractCoords_221 = [[          0, 0.000000000 , 0.0 , "Mg"],
		                   [0.166666665, 0.166666665 , 0.5 , "Mg"],
		                   [0.5        , 0.0         , 0   , "Mg"],
		                   [0.666666665, 0.166666665 , 0.5 , "Mg"],
		                   [0.00       , 0.5         , 0.00, "Mg"],
		                   [0.166666665, 0.666666665 , 0.5 , "Mg"],
		                   [0.5        , 0.5         , 0.0 , "Mg"],
		                   [0.666666665, 0.666666665 , 0.5 , "Mg"]]

	

		expectedUCell = UCell.UnitCell.fromLattVects( lattVects_221, fractCoords_221 )
		supCell = tCode.superCellFromUCell( self.startUCellA , [2,2,1])

		self.assertTrue(expectedUCell==supCell)

	def testCellB_3x1x2(self):
		lattVects_312 = [ [18.18, 0.0        , 0],
		                  [ 3.0300000097, 5.248113941, 0],
		                  [ 0.00, 0.0        , 19.67808706] ]

		fractCoords_312 = [ [0, 0, 0, "Mg"],
		                    [0.11111111, 0.3333333, 0.25, "Mg"],
		                    [0.33333333, 0, 0, "Mg"],
		                    [0.66666667, 0, 0, "Mg"],
		                    [0.44444444, 0.3333333, 0.25, "Mg"],
		                    [0.77777778, 0.3333333, 0.25, "Mg"],
		                    [0, 0, 0.5, "Mg"],
		                    [0.11111111, 0.3333333, 0.75, "Mg"],
		                    [0.33333333, 0, 0.5, "Mg"],
		                    [0.66666667, 0, 0.5, "Mg"],
		                    [0.44444444, 0.3333333, 0.75, "Mg"],
		                    [0.77777778, 0.3333333, 0.75, "Mg"] ]
		
		expectedUCell = UCell.UnitCell.fromLattVects( lattVects_312, fractCoords_312 )
		supCell = tCode.superCellFromUCell( self.startUCellA , [3,1,2])


		self.assertTrue(expectedUCell==supCell)



class testSurroundCell(unittest.TestCase):

	def setUp(self):
		self.lattParams = [1,2,3]
		self.lattAngles = [90,90,90]
		self.fractPositions = [ [0.5,0.5,0.5] ]
		self.createTestObjs()

	def createTestObjs(self):
		eleList = ["X" for x in self.fractPositions]
		kwargsDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles, "fractCoords":self.fractPositions,
		             "elementList":eleList}
		self.testCellA = UCell.UnitCell(**kwargsDict)

	def testForSimpleOneAtomCubicCell(self):
		expCell = self._loadExpectedForSingleAtomCubicCell()
		actCell = tCode.getUnitCellSurroundedByNeighbourCells(self.testCellA)
		self.assertEqual(expCell,actCell)

	def testForSimpleOneAtomCubic_removedCentralAtom(self):
		expCell = self._loadExpectedForSingleAtomCubicCell()
		cartCoords = copy.deepcopy(expCell.cartCoords)
		expStartCart = [0.5,1,1.5]
		[self.assertAlmostEqual(e,a) for e,a in zip(cartCoords[0][:3],expStartCart)]
		cartCoords.pop(0)
		expCell.cartCoords = cartCoords
		actCell = tCode.getUnitCellSurroundedByNeighbourCells(self.testCellA,removeStartCoords=True)
		self.assertEqual(expCell,actCell)

	def _loadExpectedForSingleAtomCubicCell(self):
		a,b,c = self.lattParams
		outLattParams = [3*x for x in self.lattParams]
		outCell = UCell.UnitCell(lattParams=outLattParams, lattAngles=self.lattAngles)

		#x direction cells
		startCartCoords = [[0.5*a,0.5*b,0.5*c,"X"]]
		firstImageCoords  = [[0.5*a + a, 0.5*b, 0.5*c, "X"]]
		secondImageCoords = [[0.5*a - a, 0.5*b, 0.5*c, "X"]]
		xDirCartCoords = startCartCoords + firstImageCoords + secondImageCoords

		#y direction cells
		plusYCells = list()
		minusYCells = list()
		
		for coords in xDirCartCoords:
			newCoordsPlus  = [ coords[0], coords[1]+b, coords[2], coords[3] ]
			newCoordsMinus = [ coords[0], coords[1]-b, coords[2], coords[3] ]
			plusYCells.append ( newCoordsPlus )
			minusYCells.append( newCoordsMinus )
		allXyCells = xDirCartCoords + plusYCells + minusYCells

		#z direction cells
		plusZCells = list()
		minusZCells = list()

		for coords in allXyCells:
			newCoordsPlus  = [ coords[0], coords[1], coords[2]+c, coords[3] ]
			newCoordsMinus = [ coords[0], coords[1], coords[2]-c, coords[3] ]
			plusZCells.append(newCoordsPlus)
			minusZCells.append(newCoordsMinus)
		allCellsCartCoords = allXyCells + plusZCells + minusZCells

		outCell.cartCoords = allCellsCartCoords
		return outCell

if __name__ == '__main__':
	unittest.main()

