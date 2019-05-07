#!/usr/bin/python3
import itertools as it
import os
import unittest

import plato_pylib.shared.ucell_class as tCode

BOHR_TO_ANG = 1.8897259886


class testUnitCellClass(unittest.TestCase):
	def setUp(self):
		lattParams = [3.746172, 3.746172, 3.746172]
		lattAngles = [60.0,60.0,60.0]
		self.testObjA = tCode.UnitCell(lattParams=lattParams, lattAngles=lattAngles)

		self.testLattVectsA = [ [2.0 ,0.0     , 0.0     ], #Not neccesarily relted to the above (though might be, im unsure)
		                        [0.0, 2.0, 0.0     ],
		                        [0.0     , 3.40000, 12.09628] ] 

		self.testLattVectsB = [ [4.8,0.0,0.0], #These are in the native format, hence can test the lattVects getter/setter
		                        [1.6,0.5,0.0],
		                        [3.5,1.7,12.2] ]


	def testCorrectLatticeParams_constructor(self):
		''' Test tCode.UnitCell class correctly converts lattice param list into a dict '''
		inpLatticeParams = [3.5, 2.04, 1.004]
		testObj = tCode.UnitCell( lattParams=inpLatticeParams )

		self.assertEqual(inpLatticeParams[0], testObj.lattParams["a"])
		self.assertEqual(inpLatticeParams[1], testObj.lattParams["b"])
		self.assertEqual(inpLatticeParams[2], testObj.lattParams["c"])

	def testCorrectLatticeAngles_constructor(self):
		''' Test tCode.UnitCell class correctly converts lattice angles list into a dict '''
		inpLatticeAngles = [20, 45, 90]
		testObj = tCode.UnitCell( lattAngles=inpLatticeAngles )

		self.assertEqual(inpLatticeAngles[0], testObj.lattAngles["alpha"])
		self.assertEqual(inpLatticeAngles[1], testObj.lattAngles["beta"])
		self.assertEqual(inpLatticeAngles[2], testObj.lattAngles["gamma"]) 


	def testCorrectLattParamsAnglesFromVectors(self):
		''' Test tCode.UnitCell class correctly converts vectors to lattice params and angles (hcp)'''
		
		latticeVectors = self.testLattVectsA

		expectedLattParams = [2, 2, 12.56503044]
		expectedLattAngles = [74.3004871146, 90.0, 90.0]
		
		testObj = tCode.UnitCell.fromLattVects(latticeVectors)
		actualLattParams = testObj.getLattParamsList()
		actualLattAngles = testObj.getLattAnglesList()

		[self.assertAlmostEqual(x,y) for x,y in zip(expectedLattParams, actualLattParams)]
		[self.assertAlmostEqual(x,y) for x,y in zip(expectedLattAngles, actualLattAngles)]

	def testLattVectsSetter(self):
		''' Test that setting lattice vectors works correctly (the getter returns the value set) '''
		self.testObjA.lattVects = self.testLattVectsB
		for exp,act in zip(self.testLattVectsB,self.testObjA.lattVects):
			[self.assertAlmostEqual(x,y) for x,y in zip(exp,act)]


	def testGetVolume_lattParamsAngles(self):
		''' Test getVolume works when the object has lattice parameters and angles defined '''
		lattParams = [ [3.746172, 3.746172, 3.746172],
		               [3.209401, 3.209331, 5.210802] ]
		lattAngles = [ [60.0,60.0,60.0],
		               [90.0, 90.0, 120.000728] ]
		allTestObjs = [tCode.UnitCell( lattParams=x, lattAngles=y) for x,y in zip(lattParams, lattAngles)]

		expectedVolumes = [37.174744, 46.480470] #Come from castep calcs
		actualVolumes = [currObj.getVolume() for currObj in allTestObjs]

		[self.assertAlmostEqual(expected, actual, places=4) for expected, actual in zip(expectedVolumes, actualVolumes)]
		[print("Expected Volume:" + str(expected) + "\tActual Volume:" + str(actual)) for expected, actual in zip(expectedVolumes, actualVolumes)]

	def testVolumeSetter(self):
		''' Test that the volume setter works (set a volume, then get it and check its the set value) '''
		lattParams = [1.0,2.1,3.2]
		lattAngles = [90,90,120]
		testVolume = 2000.5
		testObj = tCode.UnitCell(lattParams=lattParams, lattAngles=lattAngles)

		testObj.volume = testVolume*2 #so we definitely dont start near the expected final volume
		testObj.volume = testVolume
		actVolume = testObj.volume
		self.assertAlmostEqual(testVolume,actVolume)

	def testConvAngstromToBohr(self):
		''' Test tCode.UnitCell.convAngstromToBohr correctly works on lattice parameters and resultant volumes (volumes'''
		''' depend on the tCode.UnitCell.getVolume() implementation) ''' 
		lattParams = [ [3.746172, 3.746172, 3.746172],
		               [3.209401, 3.209331, 5.210802] ]
		lattAngles = [ [60.0,60.0,60.0],
		               [90.0, 90.0, 120.000728] ]
		allTestObjs = [tCode.UnitCell( lattParams=x, lattAngles=y) for x,y in zip(lattParams, lattAngles)]

		expectedLattParams = [ [7.07923858616564, 7.07923858616564, 7.07923858616564],
		                       [6.06488847753883, 6.06475619671963, 9.84698796084886] ]
		expectedVolumes = [250.867643296291, 313.665691831337] #in angstron 37.1747573672726 and 46.4804700799773; 

		[currCell.convAngToBohr() for currCell in allTestObjs]
		actualLattParams = [currCell.getLattParamsList() for currCell in allTestObjs]
		actualVolumes = [currCell.getVolume() for currCell in allTestObjs]

		for expected,actual in zip(expectedLattParams, actualLattParams):
			self.assertAlmostEqual(expected[0], actual[0])
			self.assertAlmostEqual(expected[1], actual[1])
			self.assertAlmostEqual(expected[2], actual[2])
		[self.assertAlmostEqual(expected, actual) for expected, actual in zip(expectedVolumes, actualVolumes)]

	def testGetLattVectors(self):
		testLattParams = [x*BOHR_TO_ANG for x in [3.209407, 3.209440, 5.210808] ]
		testLattAngles = [90.0, 90.0, 120.0]

		testObj = tCode.UnitCell( lattParams=testLattParams, lattAngles=testLattAngles )
		expectedLattVects = [ [6.0649000000,  0.0000000000,  0.0000000000],
		                           [-3.0325000000, 5.2524000000,  0.0000000000],
		                           [0.0000000000,  0.0000000000,  9.8470000000] ]

		actualLattVects = testObj.getLattVects()

		for vectIdx in range(len(expectedLattVects)):
			currExp, currAct = expectedLattVects[vectIdx], actualLattVects[vectIdx]
			[self.assertAlmostEqual(exp,act,3) for exp,act in it.zip_longest(currExp, currAct)]



class testCartCoords(unittest.TestCase):

	def setUp(self):
		self.lattVects_312A = [ [18.18, 0.0        , 0],
		                        [ 3.03, 5.248113941, 0],
		                        [ 0.00, 0.0        , 19.67808706] ]

		self.fractCoords_312A = [ [0, 0, 0, "Mg"],
		                          [0, 0, 0.5, "Mg"],
		                          [0.11111111, 0.33333333, 0.25, "Mg"],
		                          [0.11111111, 0.33333333, 0.75, "Mg"],
		                          [0.3333333333, 0, 0, "Mg"],
		                          [0.3333333333, 0, 0.5, "Mg"],
		                          [0.4444444433, 0.33333333, 0.25, "Mg"],
		                          [0.4444444433, 0.33333333, 0.75, "Mg"],
		                          [0.6666666667, 0, 0, "Mg"],
		                          [0.6666666667, 0, 0.5, "Mg"],
		                          [0.7777777767, 0.33333333, 0.25, "X"],
		                          [0.7777777767, 0.33333333, 0.75, "X"] ]

		self.expCartCoords_312A = [ [0, 0, 0, "Mg"],
		                            [0, 0, 9.839059786, "Mg"],
		                            [3.0300049797, 1.7493741865, 4.919529893, "Mg"],
		                            [3.0300049797, 1.7493741865, 14.7585896771, "Mg"],
		                            [6.0600100124, 0, 0, "Mg"],
		                            [6.0600100124, 0, 9.839059786, "Mg"],
		                            [9.0900149921, 1.7493741865, 4.919529893, "Mg"],
		                            [9.0900149921, 1.7493741865, 14.7585896771, "Mg"],
		                            [12.1200200248, 0, 0, "Mg"],
		                            [12.1200200248, 0, 9.839059786, "Mg"],
		                            [15.1500250045, 1.7493741865, 4.919529893, "X"],
		                            [15.1500250045, 1.7493741865, 14.7585896771, "X"] ]


	def testGetCartCoordsNoSortNeeded(self):
		testUCell = tCode.UnitCell.fromLattVects( self.lattVects_312A, self.fractCoords_312A )
		actCartCoords = testUCell._getCartCoords(sort=True)
		actCartCoords = testUCell.cartCoords

		for expCart,actCart in it.zip_longest(self.expCartCoords_312A,actCartCoords):
			[self.assertAlmostEqual(exp,act,places=4) for exp,act in it.zip_longest(expCart[:3],actCart[:3])]
			self.assertEqual(expCart[3],actCart[3])


	def testSetCartCoords(self):
		testUCell = tCode.UnitCell.fromLattVects( self.lattVects_312A )
		testUCell.cartCoords = self.expCartCoords_312A
		actCartCoords = testUCell.cartCoords

		for expCart,actCart in it.zip_longest(self.expCartCoords_312A,actCartCoords):
			[self.assertAlmostEqual(exp,act,places=4) for exp,act in it.zip_longest(expCart[:3],actCart[:3])]
			self.assertEqual(expCart[3],actCart[3])


class testReadWriteFiles(unittest.TestCase):

	def setUp(self):
		self.fractCoordsA = [ [0.0,0.0,0.0], [0.3,0.5,0.2], [0.4,0.7,0.8] ]
		self.lattVectsA = [ [3.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,2.0] ]
		self.unitCellA = tCode.UnitCell.fromLattVects( self.lattVectsA, fractCoords = self.fractCoordsA )
		self.testPathA = "ucellTestA.ucell"
		self.unitCellA.writeFile(self.testPathA)

		self.unitCellA_noFract = tCode.UnitCell.fromLattVects(self.lattVectsA)
		self.testPathA_noFract = "ucellTestA_noFract.ucell"
		self.unitCellA_noFract.writeFile(self.testPathA_noFract)

	def tearDown(self):
		os.remove(self.testPathA)
		os.remove(self.testPathA_noFract)

	def testReadWriteConsistent_fractCoords_setA(self):
		''' Test that the .readFile method creates an equal object to that which was used to make the file with .writeFile '''
		expObj = self.unitCellA
		actObj = tCode.UnitCell.fromFile(self.testPathA)
		self.assertEqual(expObj,actObj)

	def testReadWriteConsistent_noFractCoords_setA(self):
		expObj = self.unitCellA_noFract
		actObj = tCode.UnitCell.fromFile(self.testPathA_noFract)
		self.assertEqual(expObj,actObj)
		

if __name__ == '__main__':
	unittest.main()

