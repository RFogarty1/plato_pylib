#!/usr/bin/env python3

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

#	@unittest.skip("")
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


if __name__ == '__main__':
	unittest.main()

