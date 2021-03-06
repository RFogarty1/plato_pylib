#!/usr/bin/python3

import copy
import itertools as it
import os
import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as tCode

BOHR_TO_ANG = 1.8897259886


class TestUnitCellClass(unittest.TestCase):
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

		self.fractCoordsA = [[0.0,0.0,0.0,"Mg"],[0.333,0.666,0.5,"Mg"]]
		self.testObjAPlusFractCoords = tCode.UnitCell.fromLattVects( self.testObjA.lattVects, fractCoords=self.fractCoordsA )

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


	def testLattParamSetter_fractCoords(self):
		''' Test that changing a lattice parameter does not alter the fractional co-ordinates '''
		#Note this WOULD be trivially easy to pass, but at time of writing i'm using the same function to
		#transform fractional co-ordinates whether lattice angles or params are changed. I'm doing this so
		#i can easily test the case where fract-coords SHOULDNT change
		expFractCoords = copy.deepcopy(self.testObjAPlusFractCoords.fractCoords)
		origLattParams = self.testObjA.getLattParamsList()
		self.testObjAPlusFractCoords.lattParams = {key:2*x for key,x in it.zip_longest( ["a","b","c"], origLattParams)}
		actFractCoords = self.testObjAPlusFractCoords.fractCoords
		for exp,act in it.zip_longest(expFractCoords,actFractCoords):
			self.assertEqual(exp[-1],act[-1]) # The element string
			[self.assertAlmostEqual(e,a) for e,a in zip(exp[:3],act[:3])]

	def testLattParamSetterGetterConsistent(self):
		''' simple test to see that lattParam getter returns newly set lattParams '''
		testParams = {key:2*val for key,val in self.testObjA.lattParams.items()}
		self.testObjA.lattParams = testParams
		actParams = self.testObjA.lattParams

		for key in self.testObjA.lattParams.keys():
			self.assertAlmostEqual(testParams[key], actParams[key])

	def testLattAngleSetterGetterConsistent(self):
		testAngles = {key:2+val for key,val in self.testObjA.lattAngles.items()}
		self.testObjA.lattAngles = testAngles
		expAngles = dict(testAngles)
		testAngles['alpha'] *= 2 #I DONT want this to affect the value in the object. 
		actAngles = self.testObjA.lattAngles
		for key in expAngles.keys():
			self.assertAlmostEqual( expAngles[key], actAngles[key] )

	def testToAndFromDictAreConsistent(self):
		lattParams, lattAngles = [2,2,3], [90,90,120]
		expObj = tCode.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
		expObj.fractCoords = self.fractCoordsA
		tempDict = expObj.toDict()
		actObj = tCode.UnitCell.fromDict(tempDict)
		self.assertEqual(expObj, actObj)

	def testToAndFromDictAreConsistent_noCoords(self):
		lattParams, lattAngles = [2,2,3], [90,90,120]
		expObj = tCode.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
		tempDict = expObj.toDict()
		actObj = tCode.UnitCell.fromDict(tempDict)
		self.assertEqual(expObj, actObj)

	def testNoErrorWhenChangingLattParamsWithEmptyCoordList(self):
		expLattParams = [10,10,10]
		self.assertFalse( all([abs(x-y)<1e-1 for x,y in it.zip_longest(expLattParams,self.testObjA.getLattParamsList())]) )
		self.testObjA.cartCoords = list()
		self.testObjA.setLattParams(expLattParams)
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expLattParams, self.testObjA.getLattParamsList())]


class TestCartCoords(unittest.TestCase):

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


class TestReadWriteFiles(unittest.TestCase):

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
		


class TestTransformLattVectsToAlignCWithZ(unittest.TestCase):

	def setUp(self):
		self.lattVectA = [2.0, 0.0, 0.0]
		self.lattVectB = [0.0, 3.0, 0.0]
		self.lattVectC = [0.0, 0.0, 5.0]

		#Note i generated this as:
		# H = np.random.rand(3,3)
		# u,*unused = np.linalg.svd(H)
		self.unitaryTransformMatrix = np.array([[-0.60826506,  0.6932105 ,  0.38661715],
		                                        [-0.72740725, -0.29191765, -0.62101754],
		                                        [-0.31763551, -0.65897138,  0.68180965]])

		self.createTestObjs()

	def createTestObjs(self):
		self.startLattVects = [self.lattVectA, self.lattVectB, self.lattVectC]

	def testExpectedOutput(self):
		transformedLattVectors = self.unitaryTransformMatrix.dot( np.array(self.startLattVects) )
		actLattVects = tCode.getLattVectorsTransformedToAlignParamCWithZ(transformedLattVectors)

		#The output latt vectors SHOULD lead to the same angles and parameters as the input. Hence we compare unit-cells made from these
		expUCell = tCode.UnitCell.fromLattVects(transformedLattVectors)
		actUCell = tCode.UnitCell.fromLattVects(actLattVects)

		self.assertAlmostEqual( actLattVects[-1][0], 0 )
		self.assertAlmostEqual( actLattVects[-1][1], 0 )
		self.assertEqual( expUCell, actUCell )

	def testVerySimpleCaseExpectedOutput(self):
		expLattVects = self.startLattVects
		actLattVects = tCode.getLattVectorsTransformedToAlignParamCWithZ(self.startLattVects)

		expUCell = tCode.UnitCell.fromLattVects(expLattVects)
		actUCell = tCode.UnitCell.fromLattVects(actLattVects)

		self.assertAlmostEqual( actLattVects[-1][0], 0 )
		self.assertAlmostEqual( actLattVects[-1][1], 0 )
		self.assertEqual( expUCell, actUCell )


class TestUcellWithLattVectsWithCAlongZ(unittest.TestCase):

	def setUp(self):
		self.lattParams = [2,2,2]
		self.lattAngles = [90,90,90]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = tCode.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)

	@mock.patch("plato_pylib.shared.ucell_class.getLattVectorsTransformedToAlignParamCWithZ")
	def testLattParamsToVectsFunctionReturnsExpected(self, mockedLattVectTransform):
		expVects = mock.Mock()
		mockedLattVectTransform.side_effect = [expVects]
		actVects = tCode.lattParamsAndAnglesToLattVects(self.lattParams, self.lattAngles, testInverse=False, putCAlongZ=True)
		self.assertTrue(mockedLattVectTransform.called)
		self.assertEqual(expVects, actVects)

	@mock.patch("plato_pylib.shared.ucell_class.lattParamsAndAnglesToLattVects")
	def testUCellReturnsExpLattVects(self, mockedGetLattVects):
		self.testCellA.putCAlongZ = True
		expVects = mock.Mock()
		mockedGetLattVects.side_effect = [ expVects ] 
		actVects = self.testCellA.lattVects
		self.assertTrue(mockedGetLattVects.called)
		self.assertEqual(expVects, actVects)


class TestFoldInAtomCoords(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [1,2,3]
		self.lattAnglesA = [90,90,90]
		self.fractCoordsA = [ [0.5,0.5,0.5] ]
		self.eleListA = [ "Mg" ]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":list(self.lattParamsA), "lattAngles": list(self.lattAnglesA),
		             "fractCoords":list(self.fractCoordsA), "elementList": self.eleListA}
		self.testCellA = tCode.UnitCell(**kwargDict)

	def testDoesntAlterWhenNotNeeded(self):
		expCell = copy.deepcopy(self.testCellA)
		self.createTestObjs()
		actCell = self.testCellA
		tCode.foldAtomicPositionsIntoCell(actCell)
		self.assertEqual(expCell, actCell)

	def testShiftsWhenCoordsAreNegative(self):
		self.fractCoordsA[0] = [-1*x for x in self.fractCoordsA[0]]
		self._runStandardTest()

	def testShiftsWhenCoordsMagnitudeGreaterThanTwo(self):
		self.fractCoordsA[0] = [x+2 for x in self.fractCoordsA[0]]
		self._runStandardTest()

	def _runStandardTest(self):
		expCell = copy.deepcopy(self.testCellA)
		self.createTestObjs()
		actCell = self.testCellA
		self.assertNotEqual(expCell, actCell)
		tCode.foldAtomicPositionsIntoCell(actCell)
		self.assertEqual(expCell,actCell)

class TestGetDensityFromUCellObj(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.fractCoordsA = [ [0.5,0.5,0.5,"O"],
		                      [0.6,0.6,0.6,"C"] ]
		self.massDict = {"O":16,"C":12}
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = tCode.UnitCell(**currKwargs)
		self.testCellA.fractCoords = self.fractCoordsA

	@mock.patch("plato_pylib.shared.ucell_class.uConv")
	def testWithNoConvFactorsA_mockedAvogadroNumber(self, mockedUnitConvs):
		mockedUnitConvs.AVOGADRO_NUMBER = 2
		totalMass = (16/mockedUnitConvs.AVOGADRO_NUMBER) + (12/mockedUnitConvs.AVOGADRO_NUMBER)
		expVal = totalMass/(10**3)
		actVal = tCode.getDensityFromUCellObj(self.testCellA, massDict=self.massDict, lenConvFactor=1)
		self.assertAlmostEqual(expVal,actVal)


class TestApplyTranslationVector(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.fractCoords = [ [0.9,0.8,0.7,"X"],
		                     [0.5,0.4,0.3,"Y"],
		                     [0.3,0.2,0.1,"Z"] ]
		self.foldInAfter = False
		self.tVect = [-0.5,-0.1,-0.2]
		self.createTestObjs()

	def createTestObjs(self):
		self.geomA = tCode.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.geomA.fractCoords = self.fractCoords

	def _runTestFunct(self):
		tCode.applyTranslationVectorToFractionalCoords(self.geomA, self.tVect, foldInAfter=self.foldInAfter)

	def testForSimpleCaseA(self):
		expFractCoords = [ [ 0.4,  0.7,  0.5, "X"],
		                   [ 0.0,  0.3,  0.1, "Y"],
		                   [-0.2,  0.1, -0.1, "Z"] ]
		self._runTestFunct()
		actFractCoords = self.geomA.fractCoords
		self.checkExpAndActFractCoordsMatch(expFractCoords,actFractCoords)

	def testFoldInAfterCaseA(self):
		self.foldInAfter = True
		expFractCoords = [ [ 0.4,  0.7,  0.5, "X"],
		                   [ 0.0,  0.3,  0.1, "Y"],
		                   [ 0.8,  0.1,  0.9, "Z"] ]

		self.foldInAfter = True
		self._runTestFunct()
		actFractCoords = self.geomA.fractCoords
		self.checkExpAndActFractCoordsMatch(expFractCoords, actFractCoords)

	def checkExpAndActFractCoordsMatch(self, expCoords, actCoords):
		expGeom, actGeom = copy.deepcopy(self.geomA), copy.deepcopy(self.geomA)
		expGeom.fractCoords = expCoords
		actGeom.fractCoords = actCoords
		self.assertEqual(expGeom,actGeom)

class TestApplyTransVectorToCartCoords(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords  = [ [1, 2, 7, "X"],
		                     [3, 3, 3, "Y"],
		                     [7, 6, 1, "Z"] ]
		self.foldInAfter = False
		self.tVect = [2,3,-4]
		self.createTestObjs()

	def createTestObjs(self):
		self.geomA = tCode.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.geomA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		tCode.applyTranslationVectorToCartCoords(self.geomA, self.tVect, foldInAfter=self.foldInAfter)

	def testForSimpleCaseA(self):
		expCartCoords = [ [3,  5,  3, "X"],
		                  [5,  6, -1, "Y"],
		                  [9,  9, -3, "Z"] ]
		self._runTestFunct()
		self.checkExpAndActCoordsMatch(expCartCoords, self.geomA.cartCoords)

	def testForFoldInAfterCaseA(self):
		expCartCoords = [ [3,  5,  3, "X"],
		                  [5,  6,  9, "Y"],
		                  [9,  9,  7, "Z"] ]
		self.foldInAfter = True
		self._runTestFunct()
		self.checkExpAndActCoordsMatch(expCartCoords, self.geomA.cartCoords)

	def checkExpAndActCoordsMatch(self, expCoords, actCoords):
		expGeom, actGeom = copy.deepcopy(self.geomA), copy.deepcopy(self.geomA)
		expGeom.cartCoords, actGeom.cartCoords = expCoords, actCoords
		self.assertEqual(expGeom, actGeom)


class TestMoveIndicesToTopOfGeom(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coordsA = [ [2,3,4,"A"],
		                 [5,6,7,"B"],
		                 [8,9,1,"C"],
		                 [3,4,5,"D"] ]
		self.indices = [2,1]
		self.createTestObjs()

	def createTestObjs(self):
		self.geomA = tCode.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.geomA.cartCoords = self.coordsA

	def _runTestFunct(self):
		tCode.moveIndicesToTopOfGeomForUnitCell(self.geomA, self.indices)

	def testMoveMiddleTwoToTop(self):
		expCoords = [ [8,9,1,"C"],
		              [5,6,7,"B"],
		              [2,3,4,"A"],
		              [3,4,5,"D"] ]
		self._runTestFunct()
		actCoords = self.geomA.cartCoords
		self.checkExpAndActCoordsMatch(expCoords, self.geomA.cartCoords)

	def checkExpAndActCoordsMatch(self, expCoords, actCoords):
		expGeom, actGeom = copy.deepcopy(self.geomA), copy.deepcopy(self.geomA)
		expGeom.cartCoords, actGeom.cartCoords = expCoords, actCoords
		self.assertEqual(expGeom, actGeom)


if __name__ == '__main__':
	unittest.main()







