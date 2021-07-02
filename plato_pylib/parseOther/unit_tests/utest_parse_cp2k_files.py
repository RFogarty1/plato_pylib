#!/usr/bin/python3

import itertools
import os
import sys
import types
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.parseOther.parse_cp2k_files as tCode
import plato_pylib.shared.custom_errors as errorHelp

import numpy as np


#This used to be neccesary when i was attaching to the CLASS rather than instances of the class
def createStubClassForTestingAttachingFunctions():
	class StubClassForTestingAttachingFuncts():
		def __init__(self):
			self.extraSingleLinePatterns = list()
			self.extraFunctsToParseFromSingleLine = list()
			self.extraHandleParsedOutputFuncts = list()
	return StubClassForTestingAttachingFuncts()

#NOTE: I can actually also test the attaching to instances rather than the full class maybe. May or may not require
#that i set an initializer to (likely shallow) copy the initial lists
class testAttachExtraCommandToParser(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.testClsA = createStubClassForTestingAttachingFunctions()

	def testPatternGetsAppended(self):
		testPattern, testFunct = mock.Mock(), mock.Mock()
		testHandleDict = mock.Mock()
		decoObj = tCode.getDecoToAttachSectionParserToCpoutParser(testPattern, testFunct, handleParsedDictFunct=testHandleDict)
		decoObj(self.testClsA)
		self.assertEqual( testPattern,self.testClsA.extraSingleLinePatterns[-1] )
		self.assertEqual( testFunct, self.testClsA.extraFunctsToParseFromSingleLine[-1] )
		self.assertEqual( testHandleDict, self.testClsA.extraHandleParsedOutputFuncts[-1] )

#Tests related to parsing of the output file (e.g. outfile.cpout in cp2k.sopt -o outfile.cpout *.inp)
class testCPoutParsing(unittest.TestCase):

	def setUp(self):
		self.fullFilePathA = createCP2KFullFileA() #hcp Mg with PBCs on
		self.fullFilePathB = createCP2KFullFile_fccOpt()

	def tearDown(self):
		os.remove( self.fullFilePathA )
		os.remove( self.fullFilePathB )

	def testUnitCellParsing(self):
		expectedVolume = 46.369 #in angstroms
		expectedAngles = [90,90,120] #alpha,beta,gamma in degrees
		expectedLengths = [3.207, 3.207, 5.207] #a,b,c in angstroms

		actualUnitCell = tCode.parseCpout(self.fullFilePathA)["unitCell"]
		actualVolume = actualUnitCell.getVolume()
		actualLengths, actualAngles = actualUnitCell.getLattParamsList(), actualUnitCell.getLattAnglesList()

		self.assertAlmostEqual(expectedVolume,actualVolume, places=1)
		[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expectedAngles,actualAngles)]
		[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expectedLengths,actualLengths)]

	def testEnergyParsing(self):
		expectedEnergy = -1.57065744839809 * tCode.HART_TO_EV
		actualEnergy = tCode.parseCpout(self.fullFilePathA)["energy"] #Here for backwards compat really.
		otherTotalElectron = tCode.parseCpout(self.fullFilePathA)["energies"].electronicTotalE #Proper way to get electronic energy
		self.assertAlmostEqual(expectedEnergy, actualEnergy)
		self.assertAlmostEqual(expectedEnergy, otherTotalElectron)

	def testAtomsParsing(self):
		expectedNumbAtoms = 2
		actualNumbAtoms = tCode.parseCpout(self.fullFilePathA)["numbAtoms"]
		self.assertEqual(expectedNumbAtoms,actualNumbAtoms)

	def testDetectedSingleGeomCase(self):
		expFlag = False
		actFlag = tCode.parseCpout(self.fullFilePathA)["multiple_geom_present"]
		self.assertEqual(expFlag,actFlag)

	def testDetectsMultipleGeomCase(self):
		expectedMultipleGeomFlag = True
		actFlag = tCode.parseCpout(self.fullFilePathB)["multiple_geom_present"]
		self.assertEqual(expectedMultipleGeomFlag, actFlag)

	def testCorrectLattParamsInCellOpt(self):
		expectedLattParams = [3.216, 3.216, 3.216]
		actLattParams = tCode.parseCpout(self.fullFilePathB)["unitCell"].getLattParamsList()
		[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expectedLattParams,actLattParams)]


class testMOInfoParsing(unittest.TestCase):

	def setUp(self):
		self.fileWithKPoints = createCP2KFullOutFileB_kPts_moPrint()
		self.fileNoKPoints = createCP2KFullOutFileC_noKPts_moPrint()
	def tearDown(self):
		os.remove(self.fileWithKPoints)
		os.remove(self.fileNoKPoints)

	def testWithKPoints(self):
		#1 row = 1 k-point
		actDict = tCode.parseMOInfo(self.fileWithKPoints)
		expDict = {"eigenVals".lower(): [ [-2.1723129646, 11.1090796822, 11.328675649, 11.3287844946, 11.795922532],
		                                  [-0.376850625 , 6.0465354889 , 7.3657983961, 12.9275084401, 14.5907233945] ],
		           "occVals".lower():  [ [2.0, 0.0, 0.0, 0.0, 0.0],
		                                 [2.0, 0.0, 0.0, 0.0, 0.0] ],
		           "eFermi".lower(): 0.7123943507 }

		self.checkExpAndActDictEqual(expDict,actDict)

	def testWithNoKPts(self):
		actDict = tCode.parseMOInfo(self.fileNoKPoints)
		expDict = {"eigenvals": [ [-4.0214633887, 15.4634745026, 15.4635289254, 15.4635289254, 17.0066872] ],
		           "occvals":   [ [2.0, 0.0, 0.0, 0.0, 0.0] ],
		           "efermi": 14.9486893106}

		self.checkExpAndActDictEqual(expDict,actDict)

	def checkExpAndActDictEqual(self,expDict,actDict,places=7):
		singleValKeys = ["efermi"]
		twoDimArrayKeys = ["eigenvals", "occvals"]
		for key in expDict.keys():
			self.assertTrue(key in actDict.keys())
			if key in singleValKeys:
				self.assertAlmostEqual( expDict[key], actDict[key], places=places)
			elif key in twoDimArrayKeys:
				expArray, actArray = np.array( expDict[key] ), np.array( actDict[key] )
				self.assertTrue( expArray.shape == actArray.shape )
				self.assertTrue( np.allclose(expArray,actArray) )
			else:
				raise ValueError("{} is not an expected key".format(key))


class TestCP2KOverlapConditionParsing(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sectionA = self._loadSectionNoDiag()
		self.fileAsListA = self.sectionA.split("\n")
		self.startIdxA = 6 #zeroth energy is a blank line

		self.sectionB = self._loadSectionWithDiag()
		self.fileAsListB = self.sectionB.split("\n")
		self.startIdxB = 6

	#NOTE: I need to load the CP2K file termination flag (PROGRAM ENDED) to avoid an error from it missing
	def _loadSectionNoDiag(self):
		return """
 RS_GRID| Information for grid number                                          4
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3            -18      17                Points:          36

 OVERLAP MATRIX CONDITION NUMBER AT GAMMA POINT
 1-Norm Condition Number (Estimate)
   CN : |A|*|A^-1|:  3.448E+001 *  2.507E+003   = 8.645E+004  Log(1-CN):  4.9368

 Number of electrons:                                                         16
 Number of occupied orbitals:                                                  8
 Number of molecular orbitals:                                                16

PROGRAM ENDED
"""

	def _loadSectionWithDiag(self):
		return """
 RS_GRID| Information for grid number                                          4
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3            -18      17                Points:          36

 OVERLAP MATRIX CONDITION NUMBER AT GAMMA POINT
 1-Norm Condition Number (Estimate)
   CN : |A|*|A^-1|:  1.840E+001 *  1.468E+003   = 2.701E+004  Log(1-CN):  4.4315
 1-Norm and 2-Norm Condition Numbers using Diagonalization
   CN : |A|*|A^-1|:  1.840E+001 *  1.468E+003   = 2.701E+004  Log(1-CN):  4.4315
   CN : max/min ev:  1.133E+001 /  1.036E-003   = 1.093E+004  Log(2-CN):  4.0386

 Number of electrons:                                                         16

"""

	def testExpOutputForNoDiagCase(self):
		expEndIdx = 9
		expObj = types.SimpleNamespace( estimate=types.SimpleNamespace(oneNorm=8.645E+004) )
		actDict, actEndIdx = tCode._parseOverlapCondSection(self.fileAsListA, self.startIdxA)
		actObj = actDict["overlap_condition_number"]
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertAlmostEqual(expObj.estimate.oneNorm, actObj.estimate.oneNorm)

	def testExpOutputForDiagCase(self):
		expEndIdx = 12
		expOneNorm, expTwoNorm = 2.701E+004, 1.093E+004
		actDict, actEndIdx = tCode._parseOverlapCondSection(self.fileAsListB, self.startIdxB)
		actObj = actDict["overlap_condition_number"]
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertAlmostEqual(expOneNorm, actObj.diag.oneNorm)
		self.assertAlmostEqual(expTwoNorm, actObj.diag.twoNorm)

	@mock.patch("plato_pylib.parseOther.parse_cp2k_files._getFileAsListFromInpFile")
	def testImportedClsObjParses(self, mockedReader):
		mockedReader.side_effect = lambda *args: self._loadSectionNoDiag().split("\n")
		expOneNorm = 8.645E+004
		actObj = tCode.parseCpout( mock.Mock() )["overlap_condition_number"]
		self.assertAlmostEqual( expOneNorm, actObj.estimate.oneNorm )


class TestParseBSSE(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sectionA = self._loadSectionBSSE_a()
		self.startIdxA = 2
		self.fileAsListA = self.sectionA.split("\n")

	def _loadSectionBSSE_a(self):
		return """
 -------------------------------------------------------------------------------
 -                                                                             -
 -                                 BSSE RESULTS                                -
 -                                                                             -
 -                 CP-corrected Total energy:       33.604063                  -
 -                                                                             -
 -                       1-body contribution:       69.221187                  -
 -                       1-body contribution:      -30.604558                  -
 -                       1-body contribution:       -3.741475                  -
 -                                                                             -
 -                       2-body contribution:       -2.753861                  -
 -                       2-body contribution:       -1.832959                  -
 -                       2-body contribution:        0.406903                  -
 -                       3-body contribution:        2.908826                  -
 -                 BSSE-free interaction energy:       -1.271091               -
 -------------------------------------------------------------------------------
"""

	def testExpectedCorrectedEnergyA(self):
		expCorrEnergy = 33.604063*tCode.HART_TO_EV
		expEndIdx = 14
		actDict, actEndIdx = tCode._parseBSSESection(self.fileAsListA, self.startIdxA)
		actCorrEnergy = actDict["bsse"].cpCorrectedTotalEnergy
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertAlmostEqual(expCorrEnergy, actCorrEnergy)


class TestParseTimingSection(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sectionA = self._loadTestSectionA()
		self.startIdxA = 3
		self.fileAsListA = self.sectionA.split("\n")

	def _loadTestSectionA(self):
		return """
 -------------------------------------------------------------------------------
 -                                                                             -
 -                                T I M I N G                                  -
 -                                                                             -
 -------------------------------------------------------------------------------
 SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    0.008    0.008  101.536  101.536
 cp_cell_opt                          1  2.0    0.000    0.000  101.482  101.482
 geoopt_bfgs                          1  3.0    0.001    0.001  101.481  101.481
 cp_eval_at                           6  4.0    0.001    0.001  101.476  101.476
 qs_energies                          6  5.8    0.009    0.009   96.070   96.070
 scf_env_do_scf                       6  6.8    0.000    0.000   93.674   93.674
 scf_env_do_scf_inner_loop           48  7.8    0.016    0.016   93.674   93.674
 qs_forces                            5  5.0    0.002    0.002   86.009   86.009
 qs_scf_new_mos_kp                   48  8.8    0.000    0.000   41.202   41.202
 do_general_diag_kp                  48  9.8    0.472    0.472   41.201   41.201
 rebuild_ks_matrix                   53  9.6    0.000    0.000   30.080   30.080
 qs_ks_build_kohn_sham_matrix        53 10.6    0.034    0.034   30.080   30.080
 sum_up_and_integrate                53 11.6    0.000    0.000   29.876   29.876
 integrate_v_rspace                  53 12.6   28.935   28.935   29.876   29.876
 qs_rho_update_rho                   54  8.8    0.000    0.000   26.456   26.456
 calculate_rho_elec                  54  9.8   26.441   26.441   26.456   26.456
 qs_ks_update_qs_env                 48  8.8    0.000    0.000   25.682   25.682
 rskp_transform                   48000 10.8   15.279   15.279   15.279   15.279
 kpoint_density_transform            53 10.6    0.163    0.163   11.302   11.302
 copy_dbcsr_to_fm                 98898 10.8    0.645    0.645    9.582    9.582
 transform_dmat                   26500 11.6    8.052    8.052    8.052    8.052
 dbcsr_desymmetrize_deep         194898 11.3    2.045    2.045    6.474    6.474
 dbcsr_complete_redistribute     151898 12.1    2.440    2.440    6.062    6.062
 qs_ks_update_qs_env_forces           5  6.0    0.000    0.000    4.507    4.507
 dbcsr_finalize                  515784 12.4    1.003    1.003    4.057    4.057
 copy_fm_to_dbcsr                 53000 11.6    0.126    0.126    2.948    2.948
 -------------------------------------------------------------------------------

 The number of warnings for this run is : 1
 
 -------------------------------------------------------------------------------
"""

	def testExpectedTotalTimeCaseA(self):
		expEndIdx = 35
		expTotal = 101.536
		actDict, actEndIdx = tCode._parseTimingSection(self.fileAsListA, self.startIdxA)
		actTotal = actDict["timings"].CP2K_total
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertAlmostEqual(expTotal, actTotal)

	def testSomeSubroutineTotalTimingsGiven(self):
		expDict = { "do_general_diag_kp": 41.201,
		            "sum_up_and_integrate": 29.876 }
		fullDict, actEndIdx = tCode._parseTimingSection(self.fileAsListA, self.startIdxA)
		actDict = fullDict["timings"].subroutineTotals
		for key in expDict.keys():
			exp, act = expDict[key], actDict[key]
			self.assertAlmostEqual(exp,act)



class TestParseNumbProcsSection(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sectionA = self._loadTestSectionA()
		self.startIdxA = 3
		self.fileAsListA = self.sectionA.split("\n")

	def _loadTestSectionA(self):
		return """
 GLOBAL| FFTs using library dependent lengths                                  F
 GLOBAL| Global print level                                               MEDIUM
 GLOBAL| Total number of message passing processes                            43
 GLOBAL| Number of threads for this process                                    1
 GLOBAL| This output is from process                                           0
 GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-1620 v4 @ 3.50GHz
		"""

	def testExpectedFromSimpleCaseA(self):
		expEndIdx = 5
		expDict = {"nMPI": 43, "nThreads": 1}
		actDict, actEndIdx = tCode._parseNumbProcsSection(self.fileAsListA, self.startIdxA)
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertEqual(expDict, actDict)


class TestParseHirshfeldOrMullikenChargesSection(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sectionA = self._loadTestSectionA()
		self.startIdxA = 3
		self.fileAsListA = self.sectionA.split("\n")

		self.sectionMulliken = self._loadTestSectionMullikenA()
		self.startIdxMulliken = 3
		self.fileAsListMulliken = self.sectionMulliken.split("\n")

	def _loadTestSectionA(self):
		return """
 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       O      1       6.000          7.115                        -1.115
      2       H      2       1.000          0.439                         0.561
      3       H      2       1.000          0.439                         0.561

  Total Charge                                                            0.007
 !-----------------------------------------------------------------------------!
 """

	def _loadTestSectionMullikenA(self):
		return """
 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     O        1          6.443303                             -0.443303
       2     H        2          0.778348                              0.221652
       3     H        2          0.778348                              0.221652
 # Total charge                              8.000000                 -0.000000

 !-----------------------------------------------------------------------------!
"""

	def testExpectedFromSimpleCaseA(self):
		expEndIdx = 10
		expChargeDict = {"charges":[-1.115, 0.561, 0.561], "total":0.007}
		actDict, actEndIdx = tCode._parseHirshfeldChargesSection(self.fileAsListA, self.startIdxA)
		self.assertEqual(expEndIdx, actEndIdx)
		self._checkExpAndActChargeDictsMatch(expChargeDict, actDict)

	def _checkExpAndActChargeDictsMatch(self, expDict, actDict):
		numbListAttrs = ["charges"]
		numbAttrs = ["total"]

		for key in numbListAttrs:
			for exp,act in itertools.zip_longest(expDict[key], actDict[key]):
				self.assertAlmostEqual(exp,act)

		for key in numbAttrs:
			self.assertAlmostEqual( expDict[key], actDict[key] )

	def testHandleChargesInfo(self):
		parserInstance = mock.Mock()
		parserInstance.outDict = dict()
		expChargeDict =  {"charges":[-1.115, 0.561, 0.561], "total":0.007}
		tCode._handleHirshfeldChargesInfo(parserInstance, expChargeDict)
		actOutDict = parserInstance.outDict
		self._checkExpAndActChargeDictsMatch( expChargeDict, actOutDict["hirshfeld_charges_final"] )

	def testParserAlsoWorksForMulliken(self):
		expEndIdx = 10
		expChargeDict = {"charges": [-0.443303, 0.221652, 0.221652], "total":0}
		actDict, actEndIdx = tCode._parseHirshfeldChargesSection(self.fileAsListMulliken, self.startIdxMulliken)
		self.assertEqual(expEndIdx, actEndIdx)
		self._checkExpAndActChargeDictsMatch(expChargeDict, actDict)


class TestParseCP2kGeomOutputXyzFiles(unittest.TestCase):

	def setUp(self):
		self.testFileStrA = getGeomOptXyzFileStrA()
		#Used as tools to compare cartesian co-ordinates
		self.testCellA = uCellHelp.UnitCell(lattParams=[60,60,60], lattAngles=[90,90,90])
		self.testCellB = uCellHelp.UnitCell(lattParams=[60,60,60], lattAngles=[90,90,90])

	@mock.patch("plato_pylib.parseOther.parse_cp2k_files._getFileStrFromInpFile")
	def testExpectedGeomsTestFileA(self, mockedReader):
		mockedReader.side_effect = lambda *args: self.testFileStrA
		expCartCoordsA = [ [0.0,-0.3, 0.0,"Mg"] ]
		expCartCoordsB = [ [0.0, 0.4, 0.0,"Mg"] ]
		expCartCoordsTotal = [expCartCoordsA,expCartCoordsB]
		actCartCoordsTotal = [x.cartCoords for x in tCode.parseXyzFromGeomOpt(None)["all_geoms"]]

		#Compare cart coords
		for expCoords,actCoords in itertools.zip_longest(expCartCoordsTotal, actCartCoordsTotal):
			expCell,actCell = self.testCellA, self.testCellB
			expCell.cartCoords = expCoords
			actCell.cartCoords = actCoords
			self.assertEqual(expCell,actCell)	
	
class TestParseCpoutRaisesParseFileError(unittest.TestCase):

	def setUp(self):
		self.fakeFileStrA = "Not really\n A CP2K file"

	@mock.patch("plato_pylib.parseOther.parse_cp2k_files._getFileAsListFromInpFile")
	def testRaisesWhenTerminateFlagMissing(self, mockedReader):
		mockedReader.side_effect = lambda *args: self.fakeFileStrA
		with self.assertRaises(errorHelp.PlatoPylibParseFileError):
			tCode.parseCpout(None)


class TestParseExtraEnergiesSection(unittest.TestCase):

	def setUp(self):
		self.sectionA = self._loadEnergySectionSmearOnA()
		self.sectionB = self._loadEnergySectionDispersionOn()
		self.startIdxA, self.startIdxB = 2, 2
		self.fileAsListA = self.sectionA.split("\n")
		self.fileAsListB = self.sectionB.split("\n")

	def _loadEnergySectionSmearOnA(self):
		return """

  Overlap energy of the core charge distribution:               0.00000000000077
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.48619741993392
  Hartree energy:                                               1.02122004110676
  Exchange-correlation energy:                                 -0.43144128983853
  Electronic entropic energy:                                  -0.00001257113930
  Fermi energy:                                                 0.11183551444274

  Total energy:                                                -0.87976787980010

		"""

	def _loadEnergySectionDispersionOn(self):
		return """

  Overlap energy of the core charge distribution:               0.00000249522169
  Self energy of the core charge distribution:                -37.95380752448555
  Core Hamiltonian energy:                                     10.92685808005664
  Hartree energy:                                              17.85044956271290
  Exchange-correlation energy:                                 -5.94473660189153
  Dispersion energy:                                           -0.00355123783846

  Total energy:                                               -15.12478522622432

		"""

	def testExpectedWhenSmearingOn(self):
		expEntropy = -0.00001257113930*tCode.HART_TO_EV
		expEnergy = -0.87976787980010*tCode.HART_TO_EV
		expFermiEnergy = 0.11183551444274*tCode.HART_TO_EV
		expEndIdx = 11
		actDict, actEndIdx = tCode._parseEnergiesSection(self.fileAsListA, self.startIdxA)
		actEntropy, actEnergy, actFermiEnergy = actDict["energies"].entropy, actDict["energies"].electronicTotalE, actDict["fermi_energy"]
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertAlmostEqual(expEntropy, actEntropy)
		self.assertAlmostEqual(expEnergy, actEnergy)

	def testExpectedWhenDispersionCorrOn(self):
		expEntropy = None
		expEnergy = -15.12478522622432*tCode.HART_TO_EV
		expDispersion = -0.00355123783846*tCode.HART_TO_EV
		expEndIdx = 10
		actDict, actEndIdx = tCode._parseEnergiesSection(self.fileAsListB, self.startIdxB)
		actEntropy, actEnergy = actDict["energies"].entropy, actDict["energies"].electronicTotalE
		actDispersion = actDict["energies"].dispersion
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertEqual(expEntropy,actEntropy)
		self.assertAlmostEqual(expEnergy, actEnergy)
		self.assertAlmostEqual(expDispersion,actDispersion)


class TestGetVersionAndCompilationInfo(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = self._loadSectionAStr().split("\n")
		self.startIdxA = 0

	def _loadSectionAStr(self):
		outStr = """ CP2K| version string:                                          CP2K version 6.1
 CP2K| source code revision number:                                    svn:18464
 CP2K| cp2kflags: libint fftw3 parallel mpi3 scalapack has_no_shared_glibc max_c
 CP2K|            ontr=4 mkl
 CP2K| is freely available from                            https://www.cp2k.org/
 CP2K| Program compiled at                          Thu 12 Jul 13:52:12 BST 2018
 CP2K| Program compiled on                                  build-2.hpc.ic.ac.uk
 CP2K| Program compiled for                                   Linux-x86-64-intel
 CP2K| Data directory path                                   /tmp/cp2k/cp2k/data
 CP2K| Input file name                                    Mg_spd_1z_template.inp
 
		"""
		return outStr

	def testExpectedDictA(self):
		expDict = {"version_string": "CP2K version 6.1",
		           "cp2kflags": "libint fftw3 parallel mpi3 scalapack has_no_shared_glibc max_contr=4 mkl",
		           "source_code_number":"svn:18464"}
		expEndIdx = 4
		outDict, actEndIdx = tCode._parseCompileInfoSection(self.fileAsListA, self.startIdxA)
		actDict = outDict["cp2k_compile_info"]
		self.assertEqual(expEndIdx,actEndIdx)
		self.assertEqual(expDict, actDict)  

class TestParseForcesSection(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = self._loadSectionAStr().split("\n")
		self.startIdxA = 3

	def _loadSectionAStr(self):
		outStr = """ ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.988512016628202


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg         -2.40000000     1.20000000    -4.20000000
      2      1      Mg          3.00000000    -0.00000000     0.00000000
 SUM OF ATOMIC FORCES           1.00000000    -0.00000000     0.00000000     1.00000000

 MD_ENERGIES| Initialization proceeding"""
		return outStr

	def testExpectedDictA(self):
		expEndIdx = 8
		expForces = [ [-2.4,1.2,-4.2], [3.0,0,0] ]
		actDict, actEndIdx = tCode._parseAtomicForcesSection(self.fileAsListA, self.startIdxA)
		actForces = actDict["forces"]
		self.assertEqual(expEndIdx,actEndIdx)
		self._checkExpAndActForcesEqual(expForces, actForces)

	def testHandleAtomicForcesOutDict(self):
		mockParser = mock.Mock()
		mockParser.outDict = dict()
		expForces = [ [-2.4,1.2,-4.2], [3.0,0,0] ]
		mockOutDict = {"forces":expForces}
		tCode._handleAtomicForcesSection(mockParser, mockOutDict)
		actOutDict = mockParser.outDict
		actForces = mockParser.outDict["forces_final"]
		self._checkExpAndActForcesEqual(expForces, actForces)

	def _checkExpAndActForcesEqual(self, expForces, actForces):
		for exp,act in itertools.zip_longest(expForces,actForces):
			[self.assertAlmostEqual(e,a) for e,a in itertools.zip_longest(exp,act)]


#File contains mock results of two optimisations; of which only the second are of interest really
def getGeomOptXyzFileStrA():
	fileStr = """
       1
 i =        1, E =        -0.6797723573
 Mg         0.0000000000       -0.1000000000        0.0000000000
       1
 i =        2, E =        -0.7799071863
 Mg        -0.0000000000        0.2000000000       -0.0000000000
       1
 i =        1, E =        -0.8797723573
 Mg         0.0000000000       -0.3000000000        0.0000000000
       1
 i =        2, E =        -0.9799071863
 Mg        -0.0000000000        0.4000000000       -0.0000000000

"""
	return fileStr

def createCP2KFullFile_fccOpt():
	filePath = os.path.join( os.getcwd(), "full_file_fcc_mg_opt.cpout" )
	fileStr = _getCP2KFullFileStr_fccOpt()
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath





def createCP2KFullFileA():
	filePath = os.path.join( os.getcwd(), "full_file_hcp_mg.cpout" )
	fileStr = " DBCSR| Multiplication driver                                               BLAS\n DBCSR| Multrec recursion limit                                              512\n DBCSR| Multiplication stack size                                           1000\n DBCSR| Maximum elements for images                                    UNLIMITED\n DBCSR| Multiplicative factor virtual images                                   1\n DBCSR| Multiplication size stacks                                             3\n DBCSR| Number of 3D layers                                               SINGLE\n DBCSR| Use MPI memory allocation                                              T\n DBCSR| Use RMA algorithm                                                      F\n DBCSR| Use Communication thread                                               T\n DBCSR| Communication thread load                                             87\n\n\n  **** **** ******  **  PROGRAM STARTED AT               2018-09-13 10:45:39.061\n ***** ** ***  *** **   PROGRAM STARTED ON         cx1-138-14-1.cx1.hpc.ic.ac.uk\n **    ****   ******    PROGRAM STARTED BY                                 rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 14076\n  **** **  *******  **  PROGRAM STARTED IN /tmp/pbs.1998026.cx1/workDir/cutoff_6\n                                           000\n\n CP2K| version string:                                          CP2K version 6.1\n CP2K| source code revision number:                                    svn:18464\n CP2K| cp2kflags: libint fftw3 parallel mpi3 scalapack has_no_shared_glibc max_c\n CP2K|            ontr=4 mkl\n CP2K| is freely available from                            https://www.cp2k.org/\n CP2K| Program compiled at                          Thu 12 Jul 13:52:12 BST 2018\n CP2K| Program compiled on                                  build-2.hpc.ic.ac.uk\n CP2K| Program compiled for                                   Linux-x86-64-intel\n CP2K| Data directory path                                   /tmp/cp2k/cp2k/data\n CP2K| Input file name                                    Mg_spd_1z_template.inp\n \n GLOBAL| Force Environment number                                              1\n GLOBAL| Basis set file name                                 TIGHT_BINDING_BASIS\n GLOBAL| Potential file name                                      GTH_POTENTIALS\n GLOBAL| MM Potential file name                                     MM_POTENTIAL\n GLOBAL| Coordinate file name                                      __STD_INPUT__\n GLOBAL| Method name                                                        CP2K\n GLOBAL| Project name                                                         Mg\n GLOBAL| Preferred FFT library                                             FFTW3\n GLOBAL| Preferred diagonalization lib.                                       SL\n GLOBAL| Run type                                                         ENERGY\n GLOBAL| All-to-all communication in single precision                          F\n GLOBAL| FFTs using library dependent lengths                                  F\n GLOBAL| Global print level                                               MEDIUM\n GLOBAL| Total number of message passing processes                            12\n GLOBAL| Number of threads for this process                                    1\n GLOBAL| This output is from process                                           0\n GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-2650 0 @ 2.00GHz\n\n MEMORY| system memory details [Kb]\n MEMORY|                        rank 0           min           max       average\n MEMORY| MemTotal            132001080     132001080     132001080     132001080\n MEMORY| MemFree              74114120      74114120      74114120      74114120\n MEMORY| Buffers                     0             0             0             0\n MEMORY| Cached               50686260      50686260      50686260      50686260\n MEMORY| Slab                  1372172       1372172       1372172       1372172\n MEMORY| SReclaimable          1244172       1244172       1244172       1244172\n MEMORY| MemLikelyFree       126044552     126044552     126044552     126044552\n\n\n *** Fundamental physical constants (SI units) ***\n\n *** Literature: B. J. Mohr and B. N. Taylor,\n ***             CODATA recommended values of the fundamental physical\n ***             constants: 2006, Web Version 5.1\n ***             http://physics.nist.gov/constants\n\n Speed of light in vacuum [m/s]                             2.99792458000000E+08\n Magnetic constant or permeability of vacuum [N/A**2]       1.25663706143592E-06\n Electric constant or permittivity of vacuum [F/m]          8.85418781762039E-12\n Planck constant (h) [J*s]                                  6.62606896000000E-34\n Planck constant (h-bar) [J*s]                              1.05457162825177E-34\n Elementary charge [C]                                      1.60217648700000E-19\n Electron mass [kg]                                         9.10938215000000E-31\n Electron g factor [ ]                                     -2.00231930436220E+00\n Proton mass [kg]                                           1.67262163700000E-27\n Fine-structure constant                                    7.29735253760000E-03\n Rydberg constant [1/m]                                     1.09737315685270E+07\n Avogadro constant [1/mol]                                  6.02214179000000E+23\n Boltzmann constant [J/K]                                   1.38065040000000E-23\n Atomic mass unit [kg]                                      1.66053878200000E-27\n Bohr radius [m]                                            5.29177208590000E-11\n\n *** Conversion factors ***\n\n [u] -> [a.u.]                                              1.82288848426455E+03\n [Angstrom] -> [Bohr] = [a.u.]                              1.88972613288564E+00\n [a.u.] = [Bohr] -> [Angstrom]                              5.29177208590000E-01\n [a.u.] -> [s]                                              2.41888432650478E-17\n [a.u.] -> [fs]                                             2.41888432650478E-02\n [a.u.] -> [J]                                              4.35974393937059E-18\n [a.u.] -> [N]                                              8.23872205491840E-08\n [a.u.] -> [K]                                              3.15774647902944E+05\n [a.u.] -> [kJ/mol]                                         2.62549961709828E+03\n [a.u.] -> [kcal/mol]                                       6.27509468713739E+02\n [a.u.] -> [Pa]                                             2.94210107994716E+13\n [a.u.] -> [bar]                                            2.94210107994716E+08\n [a.u.] -> [atm]                                            2.90362800883016E+08\n [a.u.] -> [eV]                                             2.72113838565563E+01\n [a.u.] -> [Hz]                                             6.57968392072181E+15\n [a.u.] -> [1/cm] (wave numbers)                            2.19474631370540E+05\n [a.u./Bohr**2] -> [1/cm]                                   5.14048714338585E+03\n \n\n CELL_TOP| Volume [angstrom^3]:                                           46.369\n CELL_TOP| Vector a [angstrom     3.207     0.000     0.000    |a| =       3.207\n CELL_TOP| Vector b [angstrom    -1.603     2.777     0.000    |b| =       3.207\n CELL_TOP| Vector c [angstrom     0.000     0.000     5.207    |c| =       5.207\n CELL_TOP| Angle (b,c), alpha [degree]:                                   90.000\n CELL_TOP| Angle (a,c), beta  [degree]:                                   90.000\n CELL_TOP| Angle (a,b), gamma [degree]:                                  120.000\n CELL_TOP| Numerically orthorhombic:                                          NO\n \n GENERATE|  Preliminary Number of Bonds generated:                             0\n GENERATE|  Achieved consistency in connectivity generation.\n\n CELL| Volume [angstrom^3]:                                               46.369\n CELL| Vector a [angstrom]:       3.207     0.000     0.000    |a| =       3.207\n CELL| Vector b [angstrom]:      -1.603     2.777     0.000    |b| =       3.207\n CELL| Vector c [angstrom]:       0.000     0.000     5.207    |c| =       5.207\n CELL| Angle (b,c), alpha [degree]:                                       90.000\n CELL| Angle (a,c), beta  [degree]:                                       90.000\n CELL| Angle (a,b), gamma [degree]:                                      120.000\n CELL| Numerically orthorhombic:                                              NO\n\n CELL_REF| Volume [angstrom^3]:                                           46.369\n CELL_REF| Vector a [angstrom     3.207     0.000     0.000    |a| =       3.207\n CELL_REF| Vector b [angstrom    -1.603     2.777     0.000    |b| =       3.207\n CELL_REF| Vector c [angstrom     0.000     0.000     5.207    |c| =       5.207\n CELL_REF| Angle (b,c), alpha [degree]:                                   90.000\n CELL_REF| Angle (a,c), beta  [degree]:                                   90.000\n CELL_REF| Angle (a,b), gamma [degree]:                                  120.000\n CELL_REF| Numerically orthorhombic:                                          NO\n\n *** WARNING in cryssym.F:163 :: Symmetry library SPGLIB not available ***\n\n\n *******************************************************************************\n                                    Kpoints\n *******************************************************************************\n BRILLOUIN| K-point scheme                                        Monkhorst-Pack\n BRILLOUIN| K-Point grid                                            10   10    6\n BRILLOUIN| Accuracy in Symmetry determination                      0.100000E-05\n BRILLOUIN| K-Point point group symmetrization                               OFF\n BRILLOUIN| Wavefunction type                                            COMPLEX\n BRILLOUIN| List of Kpoints [2 Pi/Bohr]                                      300\n BRILLOUIN| Number           Weight            X              Y              Z\n BRILLOUIN|     1           0.00333       -0.45000       -0.45000       -0.41667\n BRILLOUIN|     2           0.00333       -0.45000       -0.45000       -0.25000\n BRILLOUIN|     3           0.00333       -0.45000       -0.45000       -0.08333\n BRILLOUIN|     4           0.00333       -0.45000       -0.45000        0.08333\n BRILLOUIN|     5           0.00333       -0.45000       -0.45000        0.25000\n BRILLOUIN|     6           0.00333       -0.45000       -0.45000        0.41667\n BRILLOUIN|     7           0.00333       -0.45000       -0.35000       -0.41667\n BRILLOUIN|     8           0.00333       -0.45000       -0.35000       -0.25000\n BRILLOUIN|     9           0.00333       -0.45000       -0.35000       -0.08333\n BRILLOUIN|    10           0.00333       -0.45000       -0.35000        0.08333\n BRILLOUIN|    11           0.00333       -0.45000       -0.35000        0.25000\n BRILLOUIN|    12           0.00333       -0.45000       -0.35000        0.41667\n BRILLOUIN|    13           0.00333       -0.45000       -0.25000       -0.41667\n BRILLOUIN|    14           0.00333       -0.45000       -0.25000       -0.25000\n BRILLOUIN|    15           0.00333       -0.45000       -0.25000       -0.08333\n BRILLOUIN|    16           0.00333       -0.45000       -0.25000        0.08333\n BRILLOUIN|    17           0.00333       -0.45000       -0.25000        0.25000\n BRILLOUIN|    18           0.00333       -0.45000       -0.25000        0.41667\n BRILLOUIN|    19           0.00333       -0.45000       -0.15000       -0.41667\n BRILLOUIN|    20           0.00333       -0.45000       -0.15000       -0.25000\n BRILLOUIN|    21           0.00333       -0.45000       -0.15000       -0.08333\n BRILLOUIN|    22           0.00333       -0.45000       -0.15000        0.08333\n BRILLOUIN|    23           0.00333       -0.45000       -0.15000        0.25000\n BRILLOUIN|    24           0.00333       -0.45000       -0.15000        0.41667\n BRILLOUIN|    25           0.00333       -0.45000       -0.05000       -0.41667\n BRILLOUIN|    26           0.00333       -0.45000       -0.05000       -0.25000\n BRILLOUIN|    27           0.00333       -0.45000       -0.05000       -0.08333\n BRILLOUIN|    28           0.00333       -0.45000       -0.05000        0.08333\n BRILLOUIN|    29           0.00333       -0.45000       -0.05000        0.25000\n BRILLOUIN|    30           0.00333       -0.45000       -0.05000        0.41667\n BRILLOUIN|    31           0.00333       -0.45000        0.05000       -0.41667\n BRILLOUIN|    32           0.00333       -0.45000        0.05000       -0.25000\n BRILLOUIN|    33           0.00333       -0.45000        0.05000       -0.08333\n BRILLOUIN|    34           0.00333       -0.45000        0.05000        0.08333\n BRILLOUIN|    35           0.00333       -0.45000        0.05000        0.25000\n BRILLOUIN|    36           0.00333       -0.45000        0.05000        0.41667\n BRILLOUIN|    37           0.00333       -0.45000        0.15000       -0.41667\n BRILLOUIN|    38           0.00333       -0.45000        0.15000       -0.25000\n BRILLOUIN|    39           0.00333       -0.45000        0.15000       -0.08333\n BRILLOUIN|    40           0.00333       -0.45000        0.15000        0.08333\n BRILLOUIN|    41           0.00333       -0.45000        0.15000        0.25000\n BRILLOUIN|    42           0.00333       -0.45000        0.15000        0.41667\n BRILLOUIN|    43           0.00333       -0.45000        0.25000       -0.41667\n BRILLOUIN|    44           0.00333       -0.45000        0.25000       -0.25000\n BRILLOUIN|    45           0.00333       -0.45000        0.25000       -0.08333\n BRILLOUIN|    46           0.00333       -0.45000        0.25000        0.08333\n BRILLOUIN|    47           0.00333       -0.45000        0.25000        0.25000\n BRILLOUIN|    48           0.00333       -0.45000        0.25000        0.41667\n BRILLOUIN|    49           0.00333       -0.45000        0.35000       -0.41667\n BRILLOUIN|    50           0.00333       -0.45000        0.35000       -0.25000\n BRILLOUIN|    51           0.00333       -0.45000        0.35000       -0.08333\n BRILLOUIN|    52           0.00333       -0.45000        0.35000        0.08333\n BRILLOUIN|    53           0.00333       -0.45000        0.35000        0.25000\n BRILLOUIN|    54           0.00333       -0.45000        0.35000        0.41667\n BRILLOUIN|    55           0.00333       -0.45000        0.45000       -0.41667\n BRILLOUIN|    56           0.00333       -0.45000        0.45000       -0.25000\n BRILLOUIN|    57           0.00333       -0.45000        0.45000       -0.08333\n BRILLOUIN|    58           0.00333       -0.45000        0.45000        0.08333\n BRILLOUIN|    59           0.00333       -0.45000        0.45000        0.25000\n BRILLOUIN|    60           0.00333       -0.45000        0.45000        0.41667\n BRILLOUIN|    61           0.00333       -0.35000       -0.45000       -0.41667\n BRILLOUIN|    62           0.00333       -0.35000       -0.45000       -0.25000\n BRILLOUIN|    63           0.00333       -0.35000       -0.45000       -0.08333\n BRILLOUIN|    64           0.00333       -0.35000       -0.45000        0.08333\n BRILLOUIN|    65           0.00333       -0.35000       -0.45000        0.25000\n BRILLOUIN|    66           0.00333       -0.35000       -0.45000        0.41667\n BRILLOUIN|    67           0.00333       -0.35000       -0.35000       -0.41667\n BRILLOUIN|    68           0.00333       -0.35000       -0.35000       -0.25000\n BRILLOUIN|    69           0.00333       -0.35000       -0.35000       -0.08333\n BRILLOUIN|    70           0.00333       -0.35000       -0.35000        0.08333\n BRILLOUIN|    71           0.00333       -0.35000       -0.35000        0.25000\n BRILLOUIN|    72           0.00333       -0.35000       -0.35000        0.41667\n BRILLOUIN|    73           0.00333       -0.35000       -0.25000       -0.41667\n BRILLOUIN|    74           0.00333       -0.35000       -0.25000       -0.25000\n BRILLOUIN|    75           0.00333       -0.35000       -0.25000       -0.08333\n BRILLOUIN|    76           0.00333       -0.35000       -0.25000        0.08333\n BRILLOUIN|    77           0.00333       -0.35000       -0.25000        0.25000\n BRILLOUIN|    78           0.00333       -0.35000       -0.25000        0.41667\n BRILLOUIN|    79           0.00333       -0.35000       -0.15000       -0.41667\n BRILLOUIN|    80           0.00333       -0.35000       -0.15000       -0.25000\n BRILLOUIN|    81           0.00333       -0.35000       -0.15000       -0.08333\n BRILLOUIN|    82           0.00333       -0.35000       -0.15000        0.08333\n BRILLOUIN|    83           0.00333       -0.35000       -0.15000        0.25000\n BRILLOUIN|    84           0.00333       -0.35000       -0.15000        0.41667\n BRILLOUIN|    85           0.00333       -0.35000       -0.05000       -0.41667\n BRILLOUIN|    86           0.00333       -0.35000       -0.05000       -0.25000\n BRILLOUIN|    87           0.00333       -0.35000       -0.05000       -0.08333\n BRILLOUIN|    88           0.00333       -0.35000       -0.05000        0.08333\n BRILLOUIN|    89           0.00333       -0.35000       -0.05000        0.25000\n BRILLOUIN|    90           0.00333       -0.35000       -0.05000        0.41667\n BRILLOUIN|    91           0.00333       -0.35000        0.05000       -0.41667\n BRILLOUIN|    92           0.00333       -0.35000        0.05000       -0.25000\n BRILLOUIN|    93           0.00333       -0.35000        0.05000       -0.08333\n BRILLOUIN|    94           0.00333       -0.35000        0.05000        0.08333\n BRILLOUIN|    95           0.00333       -0.35000        0.05000        0.25000\n BRILLOUIN|    96           0.00333       -0.35000        0.05000        0.41667\n BRILLOUIN|    97           0.00333       -0.35000        0.15000       -0.41667\n BRILLOUIN|    98           0.00333       -0.35000        0.15000       -0.25000\n BRILLOUIN|    99           0.00333       -0.35000        0.15000       -0.08333\n BRILLOUIN|   100           0.00333       -0.35000        0.15000        0.08333\n BRILLOUIN|   101           0.00333       -0.35000        0.15000        0.25000\n BRILLOUIN|   102           0.00333       -0.35000        0.15000        0.41667\n BRILLOUIN|   103           0.00333       -0.35000        0.25000       -0.41667\n BRILLOUIN|   104           0.00333       -0.35000        0.25000       -0.25000\n BRILLOUIN|   105           0.00333       -0.35000        0.25000       -0.08333\n BRILLOUIN|   106           0.00333       -0.35000        0.25000        0.08333\n BRILLOUIN|   107           0.00333       -0.35000        0.25000        0.25000\n BRILLOUIN|   108           0.00333       -0.35000        0.25000        0.41667\n BRILLOUIN|   109           0.00333       -0.35000        0.35000       -0.41667\n BRILLOUIN|   110           0.00333       -0.35000        0.35000       -0.25000\n BRILLOUIN|   111           0.00333       -0.35000        0.35000       -0.08333\n BRILLOUIN|   112           0.00333       -0.35000        0.35000        0.08333\n BRILLOUIN|   113           0.00333       -0.35000        0.35000        0.25000\n BRILLOUIN|   114           0.00333       -0.35000        0.35000        0.41667\n BRILLOUIN|   115           0.00333       -0.35000        0.45000       -0.41667\n BRILLOUIN|   116           0.00333       -0.35000        0.45000       -0.25000\n BRILLOUIN|   117           0.00333       -0.35000        0.45000       -0.08333\n BRILLOUIN|   118           0.00333       -0.35000        0.45000        0.08333\n BRILLOUIN|   119           0.00333       -0.35000        0.45000        0.25000\n BRILLOUIN|   120           0.00333       -0.35000        0.45000        0.41667\n BRILLOUIN|   121           0.00333       -0.25000       -0.45000       -0.41667\n BRILLOUIN|   122           0.00333       -0.25000       -0.45000       -0.25000\n BRILLOUIN|   123           0.00333       -0.25000       -0.45000       -0.08333\n BRILLOUIN|   124           0.00333       -0.25000       -0.45000        0.08333\n BRILLOUIN|   125           0.00333       -0.25000       -0.45000        0.25000\n BRILLOUIN|   126           0.00333       -0.25000       -0.45000        0.41667\n BRILLOUIN|   127           0.00333       -0.25000       -0.35000       -0.41667\n BRILLOUIN|   128           0.00333       -0.25000       -0.35000       -0.25000\n BRILLOUIN|   129           0.00333       -0.25000       -0.35000       -0.08333\n BRILLOUIN|   130           0.00333       -0.25000       -0.35000        0.08333\n BRILLOUIN|   131           0.00333       -0.25000       -0.35000        0.25000\n BRILLOUIN|   132           0.00333       -0.25000       -0.35000        0.41667\n BRILLOUIN|   133           0.00333       -0.25000       -0.25000       -0.41667\n BRILLOUIN|   134           0.00333       -0.25000       -0.25000       -0.25000\n BRILLOUIN|   135           0.00333       -0.25000       -0.25000       -0.08333\n BRILLOUIN|   136           0.00333       -0.25000       -0.25000        0.08333\n BRILLOUIN|   137           0.00333       -0.25000       -0.25000        0.25000\n BRILLOUIN|   138           0.00333       -0.25000       -0.25000        0.41667\n BRILLOUIN|   139           0.00333       -0.25000       -0.15000       -0.41667\n BRILLOUIN|   140           0.00333       -0.25000       -0.15000       -0.25000\n BRILLOUIN|   141           0.00333       -0.25000       -0.15000       -0.08333\n BRILLOUIN|   142           0.00333       -0.25000       -0.15000        0.08333\n BRILLOUIN|   143           0.00333       -0.25000       -0.15000        0.25000\n BRILLOUIN|   144           0.00333       -0.25000       -0.15000        0.41667\n BRILLOUIN|   145           0.00333       -0.25000       -0.05000       -0.41667\n BRILLOUIN|   146           0.00333       -0.25000       -0.05000       -0.25000\n BRILLOUIN|   147           0.00333       -0.25000       -0.05000       -0.08333\n BRILLOUIN|   148           0.00333       -0.25000       -0.05000        0.08333\n BRILLOUIN|   149           0.00333       -0.25000       -0.05000        0.25000\n BRILLOUIN|   150           0.00333       -0.25000       -0.05000        0.41667\n BRILLOUIN|   151           0.00333       -0.25000        0.05000       -0.41667\n BRILLOUIN|   152           0.00333       -0.25000        0.05000       -0.25000\n BRILLOUIN|   153           0.00333       -0.25000        0.05000       -0.08333\n BRILLOUIN|   154           0.00333       -0.25000        0.05000        0.08333\n BRILLOUIN|   155           0.00333       -0.25000        0.05000        0.25000\n BRILLOUIN|   156           0.00333       -0.25000        0.05000        0.41667\n BRILLOUIN|   157           0.00333       -0.25000        0.15000       -0.41667\n BRILLOUIN|   158           0.00333       -0.25000        0.15000       -0.25000\n BRILLOUIN|   159           0.00333       -0.25000        0.15000       -0.08333\n BRILLOUIN|   160           0.00333       -0.25000        0.15000        0.08333\n BRILLOUIN|   161           0.00333       -0.25000        0.15000        0.25000\n BRILLOUIN|   162           0.00333       -0.25000        0.15000        0.41667\n BRILLOUIN|   163           0.00333       -0.25000        0.25000       -0.41667\n BRILLOUIN|   164           0.00333       -0.25000        0.25000       -0.25000\n BRILLOUIN|   165           0.00333       -0.25000        0.25000       -0.08333\n BRILLOUIN|   166           0.00333       -0.25000        0.25000        0.08333\n BRILLOUIN|   167           0.00333       -0.25000        0.25000        0.25000\n BRILLOUIN|   168           0.00333       -0.25000        0.25000        0.41667\n BRILLOUIN|   169           0.00333       -0.25000        0.35000       -0.41667\n BRILLOUIN|   170           0.00333       -0.25000        0.35000       -0.25000\n BRILLOUIN|   171           0.00333       -0.25000        0.35000       -0.08333\n BRILLOUIN|   172           0.00333       -0.25000        0.35000        0.08333\n BRILLOUIN|   173           0.00333       -0.25000        0.35000        0.25000\n BRILLOUIN|   174           0.00333       -0.25000        0.35000        0.41667\n BRILLOUIN|   175           0.00333       -0.25000        0.45000       -0.41667\n BRILLOUIN|   176           0.00333       -0.25000        0.45000       -0.25000\n BRILLOUIN|   177           0.00333       -0.25000        0.45000       -0.08333\n BRILLOUIN|   178           0.00333       -0.25000        0.45000        0.08333\n BRILLOUIN|   179           0.00333       -0.25000        0.45000        0.25000\n BRILLOUIN|   180           0.00333       -0.25000        0.45000        0.41667\n BRILLOUIN|   181           0.00333       -0.15000       -0.45000       -0.41667\n BRILLOUIN|   182           0.00333       -0.15000       -0.45000       -0.25000\n BRILLOUIN|   183           0.00333       -0.15000       -0.45000       -0.08333\n BRILLOUIN|   184           0.00333       -0.15000       -0.45000        0.08333\n BRILLOUIN|   185           0.00333       -0.15000       -0.45000        0.25000\n BRILLOUIN|   186           0.00333       -0.15000       -0.45000        0.41667\n BRILLOUIN|   187           0.00333       -0.15000       -0.35000       -0.41667\n BRILLOUIN|   188           0.00333       -0.15000       -0.35000       -0.25000\n BRILLOUIN|   189           0.00333       -0.15000       -0.35000       -0.08333\n BRILLOUIN|   190           0.00333       -0.15000       -0.35000        0.08333\n BRILLOUIN|   191           0.00333       -0.15000       -0.35000        0.25000\n BRILLOUIN|   192           0.00333       -0.15000       -0.35000        0.41667\n BRILLOUIN|   193           0.00333       -0.15000       -0.25000       -0.41667\n BRILLOUIN|   194           0.00333       -0.15000       -0.25000       -0.25000\n BRILLOUIN|   195           0.00333       -0.15000       -0.25000       -0.08333\n BRILLOUIN|   196           0.00333       -0.15000       -0.25000        0.08333\n BRILLOUIN|   197           0.00333       -0.15000       -0.25000        0.25000\n BRILLOUIN|   198           0.00333       -0.15000       -0.25000        0.41667\n BRILLOUIN|   199           0.00333       -0.15000       -0.15000       -0.41667\n BRILLOUIN|   200           0.00333       -0.15000       -0.15000       -0.25000\n BRILLOUIN|   201           0.00333       -0.15000       -0.15000       -0.08333\n BRILLOUIN|   202           0.00333       -0.15000       -0.15000        0.08333\n BRILLOUIN|   203           0.00333       -0.15000       -0.15000        0.25000\n BRILLOUIN|   204           0.00333       -0.15000       -0.15000        0.41667\n BRILLOUIN|   205           0.00333       -0.15000       -0.05000       -0.41667\n BRILLOUIN|   206           0.00333       -0.15000       -0.05000       -0.25000\n BRILLOUIN|   207           0.00333       -0.15000       -0.05000       -0.08333\n BRILLOUIN|   208           0.00333       -0.15000       -0.05000        0.08333\n BRILLOUIN|   209           0.00333       -0.15000       -0.05000        0.25000\n BRILLOUIN|   210           0.00333       -0.15000       -0.05000        0.41667\n BRILLOUIN|   211           0.00333       -0.15000        0.05000       -0.41667\n BRILLOUIN|   212           0.00333       -0.15000        0.05000       -0.25000\n BRILLOUIN|   213           0.00333       -0.15000        0.05000       -0.08333\n BRILLOUIN|   214           0.00333       -0.15000        0.05000        0.08333\n BRILLOUIN|   215           0.00333       -0.15000        0.05000        0.25000\n BRILLOUIN|   216           0.00333       -0.15000        0.05000        0.41667\n BRILLOUIN|   217           0.00333       -0.15000        0.15000       -0.41667\n BRILLOUIN|   218           0.00333       -0.15000        0.15000       -0.25000\n BRILLOUIN|   219           0.00333       -0.15000        0.15000       -0.08333\n BRILLOUIN|   220           0.00333       -0.15000        0.15000        0.08333\n BRILLOUIN|   221           0.00333       -0.15000        0.15000        0.25000\n BRILLOUIN|   222           0.00333       -0.15000        0.15000        0.41667\n BRILLOUIN|   223           0.00333       -0.15000        0.25000       -0.41667\n BRILLOUIN|   224           0.00333       -0.15000        0.25000       -0.25000\n BRILLOUIN|   225           0.00333       -0.15000        0.25000       -0.08333\n BRILLOUIN|   226           0.00333       -0.15000        0.25000        0.08333\n BRILLOUIN|   227           0.00333       -0.15000        0.25000        0.25000\n BRILLOUIN|   228           0.00333       -0.15000        0.25000        0.41667\n BRILLOUIN|   229           0.00333       -0.15000        0.35000       -0.41667\n BRILLOUIN|   230           0.00333       -0.15000        0.35000       -0.25000\n BRILLOUIN|   231           0.00333       -0.15000        0.35000       -0.08333\n BRILLOUIN|   232           0.00333       -0.15000        0.35000        0.08333\n BRILLOUIN|   233           0.00333       -0.15000        0.35000        0.25000\n BRILLOUIN|   234           0.00333       -0.15000        0.35000        0.41667\n BRILLOUIN|   235           0.00333       -0.15000        0.45000       -0.41667\n BRILLOUIN|   236           0.00333       -0.15000        0.45000       -0.25000\n BRILLOUIN|   237           0.00333       -0.15000        0.45000       -0.08333\n BRILLOUIN|   238           0.00333       -0.15000        0.45000        0.08333\n BRILLOUIN|   239           0.00333       -0.15000        0.45000        0.25000\n BRILLOUIN|   240           0.00333       -0.15000        0.45000        0.41667\n BRILLOUIN|   241           0.00333       -0.05000       -0.45000       -0.41667\n BRILLOUIN|   242           0.00333       -0.05000       -0.45000       -0.25000\n BRILLOUIN|   243           0.00333       -0.05000       -0.45000       -0.08333\n BRILLOUIN|   244           0.00333       -0.05000       -0.45000        0.08333\n BRILLOUIN|   245           0.00333       -0.05000       -0.45000        0.25000\n BRILLOUIN|   246           0.00333       -0.05000       -0.45000        0.41667\n BRILLOUIN|   247           0.00333       -0.05000       -0.35000       -0.41667\n BRILLOUIN|   248           0.00333       -0.05000       -0.35000       -0.25000\n BRILLOUIN|   249           0.00333       -0.05000       -0.35000       -0.08333\n BRILLOUIN|   250           0.00333       -0.05000       -0.35000        0.08333\n BRILLOUIN|   251           0.00333       -0.05000       -0.35000        0.25000\n BRILLOUIN|   252           0.00333       -0.05000       -0.35000        0.41667\n BRILLOUIN|   253           0.00333       -0.05000       -0.25000       -0.41667\n BRILLOUIN|   254           0.00333       -0.05000       -0.25000       -0.25000\n BRILLOUIN|   255           0.00333       -0.05000       -0.25000       -0.08333\n BRILLOUIN|   256           0.00333       -0.05000       -0.25000        0.08333\n BRILLOUIN|   257           0.00333       -0.05000       -0.25000        0.25000\n BRILLOUIN|   258           0.00333       -0.05000       -0.25000        0.41667\n BRILLOUIN|   259           0.00333       -0.05000       -0.15000       -0.41667\n BRILLOUIN|   260           0.00333       -0.05000       -0.15000       -0.25000\n BRILLOUIN|   261           0.00333       -0.05000       -0.15000       -0.08333\n BRILLOUIN|   262           0.00333       -0.05000       -0.15000        0.08333\n BRILLOUIN|   263           0.00333       -0.05000       -0.15000        0.25000\n BRILLOUIN|   264           0.00333       -0.05000       -0.15000        0.41667\n BRILLOUIN|   265           0.00333       -0.05000       -0.05000       -0.41667\n BRILLOUIN|   266           0.00333       -0.05000       -0.05000       -0.25000\n BRILLOUIN|   267           0.00333       -0.05000       -0.05000       -0.08333\n BRILLOUIN|   268           0.00333       -0.05000       -0.05000        0.08333\n BRILLOUIN|   269           0.00333       -0.05000       -0.05000        0.25000\n BRILLOUIN|   270           0.00333       -0.05000       -0.05000        0.41667\n BRILLOUIN|   271           0.00333       -0.05000        0.05000       -0.41667\n BRILLOUIN|   272           0.00333       -0.05000        0.05000       -0.25000\n BRILLOUIN|   273           0.00333       -0.05000        0.05000       -0.08333\n BRILLOUIN|   274           0.00333       -0.05000        0.05000        0.08333\n BRILLOUIN|   275           0.00333       -0.05000        0.05000        0.25000\n BRILLOUIN|   276           0.00333       -0.05000        0.05000        0.41667\n BRILLOUIN|   277           0.00333       -0.05000        0.15000       -0.41667\n BRILLOUIN|   278           0.00333       -0.05000        0.15000       -0.25000\n BRILLOUIN|   279           0.00333       -0.05000        0.15000       -0.08333\n BRILLOUIN|   280           0.00333       -0.05000        0.15000        0.08333\n BRILLOUIN|   281           0.00333       -0.05000        0.15000        0.25000\n BRILLOUIN|   282           0.00333       -0.05000        0.15000        0.41667\n BRILLOUIN|   283           0.00333       -0.05000        0.25000       -0.41667\n BRILLOUIN|   284           0.00333       -0.05000        0.25000       -0.25000\n BRILLOUIN|   285           0.00333       -0.05000        0.25000       -0.08333\n BRILLOUIN|   286           0.00333       -0.05000        0.25000        0.08333\n BRILLOUIN|   287           0.00333       -0.05000        0.25000        0.25000\n BRILLOUIN|   288           0.00333       -0.05000        0.25000        0.41667\n BRILLOUIN|   289           0.00333       -0.05000        0.35000       -0.41667\n BRILLOUIN|   290           0.00333       -0.05000        0.35000       -0.25000\n BRILLOUIN|   291           0.00333       -0.05000        0.35000       -0.08333\n BRILLOUIN|   292           0.00333       -0.05000        0.35000        0.08333\n BRILLOUIN|   293           0.00333       -0.05000        0.35000        0.25000\n BRILLOUIN|   294           0.00333       -0.05000        0.35000        0.41667\n BRILLOUIN|   295           0.00333       -0.05000        0.45000       -0.41667\n BRILLOUIN|   296           0.00333       -0.05000        0.45000       -0.25000\n BRILLOUIN|   297           0.00333       -0.05000        0.45000       -0.08333\n BRILLOUIN|   298           0.00333       -0.05000        0.45000        0.08333\n BRILLOUIN|   299           0.00333       -0.05000        0.45000        0.25000\n BRILLOUIN|   300           0.00333       -0.05000        0.45000        0.41667\n *******************************************************************************\n\n *******************************************************************************\n *******************************************************************************\n **                                                                           **\n **     #####                         ##              ##                      **\n **    ##   ##            ##          ##              ##                      **\n **   ##     ##                       ##            ######                    **\n **   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **\n **   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **\n **   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **\n **    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **\n **     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **\n **           ##                                                    ##        **\n **                                                                           **\n **                                                ... make the atoms dance   **\n **                                                                           **\n **            Copyright (C) by CP2K developers group (2000 - 2018)           **\n **                                                                           **\n *******************************************************************************\n\n DFT| Spin restricted Kohn-Sham (RKS) calculation                            RKS\n DFT| Multiplicity                                                             1\n DFT| Number of spin states                                                    1\n DFT| Charge                                                                   0\n DFT| Self-interaction correction (SIC)                                       NO\n DFT| Cutoffs: density                                              1.000000E-10\n DFT|          gradient                                             1.000000E-10\n DFT|          tau                                                  1.000000E-10\n DFT|          cutoff_smoothing_range                               0.000000E+00\n DFT| XC density smoothing                                                  NONE\n DFT| XC derivatives                                                          PW\n FUNCTIONAL| ROUTINE=NEW\n FUNCTIONAL| PBE:\n FUNCTIONAL| J.P.Perdew, K.Burke, M.Ernzerhof, Phys. Rev. Letter, vol. 77, n 18,\n FUNCTIONAL|  pp. 3865-3868, (1996){spin unpolarized}                           \n\n QS| Method:                                                                 GPW\n QS| Density plane wave grid type                        NON-SPHERICAL FULLSPACE\n QS| Number of grid levels:                                                    4\n QS| Density cutoff [a.u.]:                                                220.5\n QS| Multi grid cutoff [a.u.]: 1) grid level                               220.5\n QS|                           2) grid level                                73.5\n QS|                           3) grid level                                24.5\n QS|                           4) grid level                                 8.2\n QS| Grid level progression factor:                                          3.0\n QS| Relative density cutoff [a.u.]:                                     18374.7\n QS| Consistent realspace mapping and integration \n QS| Interaction thresholds: eps_pgf_orb:                                1.0E-05\n QS|                         eps_filter_matrix:                          0.0E+00\n QS|                         eps_core_charge:                            1.0E-12\n QS|                         eps_rho_gspace:                             1.0E-10\n QS|                         eps_rho_rspace:                             1.0E-10\n QS|                         eps_gvg_rspace:                             1.0E-05\n QS|                         eps_ppl:                                    1.0E-02\n QS|                         eps_ppnl:                                   1.0E-07\n\n\n ATOMIC KIND INFORMATION\n\n  1. Atomic kind: Mg                                    Number of atoms:       2\n\n     Orbital Basis Set                               SPD-1Z-14-14-14-constrained\n\n       Number of orbital shell sets:                                           3\n       Number of orbital shells:                                               3\n       Number of primitive Cartesian functions:                               42\n       Number of Cartesian basis functions:                                   10\n       Number of spherical basis functions:                                    9\n       Norm type:                                                              2\n\n       Normalised Cartesian orbitals:\n\n                        Set   Shell   Orbital            Exponent    Coefficient\n\n                          1       1    3s                0.211334       1.477450\n                                                         0.263508      -2.000027\n                                                         0.328563      -0.343831\n                                                         0.409678       2.000027\n                                                         0.510819      -0.479560\n                                                         0.636930      -1.363869\n                                                         0.794175       0.458932\n                                                         0.990240       0.811588\n                                                         1.234710      -0.330687\n                                                         1.539530      -0.599919\n                                                         1.919610       0.255446\n                                                         2.393530       0.389734\n                                                         2.984440      -0.353964\n                                                         3.721230       0.081129\n\n                          2       1    3px               0.228352       1.086738\n                                                         0.289625      -1.999848\n                                                         0.401779       1.952662\n                                                         0.568835      -1.204299\n                                                         0.569196      -0.285228\n                                                         0.859259       0.965738\n                                                         1.263320      -0.701610\n                                                         1.705140       0.378680\n                                                         2.237760       0.038726\n                                                         2.960640      -0.316935\n                                                         4.023710       0.428984\n                                                         5.407310      -0.329659\n                                                         7.196990       0.133456\n                                                         9.587470      -0.022141\n                          2       1    3py               0.228352       1.086738\n                                                         0.289625      -1.999848\n                                                         0.401779       1.952662\n                                                         0.568835      -1.204299\n                                                         0.569196      -0.285228\n                                                         0.859259       0.965738\n                                                         1.263320      -0.701610\n                                                         1.705140       0.378680\n                                                         2.237760       0.038726\n                                                         2.960640      -0.316935\n                                                         4.023710       0.428984\n                                                         5.407310      -0.329659\n                                                         7.196990       0.133456\n                                                         9.587470      -0.022141\n                          2       1    3pz               0.228352       1.086738\n                                                         0.289625      -1.999848\n                                                         0.401779       1.952662\n                                                         0.568835      -1.204299\n                                                         0.569196      -0.285228\n                                                         0.859259       0.965738\n                                                         1.263320      -0.701610\n                                                         1.705140       0.378680\n                                                         2.237760       0.038726\n                                                         2.960640      -0.316935\n                                                         4.023710       0.428984\n                                                         5.407310      -0.329659\n                                                         7.196990       0.133456\n                                                         9.587470      -0.022141\n\n                          3       1    3dx2              0.238188       0.888463\n                                                         0.268929      -1.492528\n                                                         0.307376       0.706661\n                                                         1.729730       0.293067\n                                                         2.119890      -0.414806\n                                                         2.793100       0.204577\n                                                         6.883380      -0.056976\n                                                         9.807220       0.057477\n                                                        13.505700      -0.030453\n                                                        18.244300       0.007665\n                                                        39.342200      -0.000264\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n                          3       1    3dxy              0.238188       1.538864\n                                                         0.268929      -2.585134\n                                                         0.307376       1.223972\n                                                         1.729730       0.507606\n                                                         2.119890      -0.718466\n                                                         2.793100       0.354338\n                                                         6.883380      -0.098686\n                                                         9.807220       0.099553\n                                                        13.505700      -0.052746\n                                                        18.244300       0.013276\n                                                        39.342200      -0.000458\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n                          3       1    3dxz              0.238188       1.538864\n                                                         0.268929      -2.585134\n                                                         0.307376       1.223972\n                                                         1.729730       0.507606\n                                                         2.119890      -0.718466\n                                                         2.793100       0.354338\n                                                         6.883380      -0.098686\n                                                         9.807220       0.099553\n                                                        13.505700      -0.052746\n                                                        18.244300       0.013276\n                                                        39.342200      -0.000458\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n                          3       1    3dy2              0.238188       0.888463\n                                                         0.268929      -1.492528\n                                                         0.307376       0.706661\n                                                         1.729730       0.293067\n                                                         2.119890      -0.414806\n                                                         2.793100       0.204577\n                                                         6.883380      -0.056976\n                                                         9.807220       0.057477\n                                                        13.505700      -0.030453\n                                                        18.244300       0.007665\n                                                        39.342200      -0.000264\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n                          3       1    3dyz              0.238188       1.538864\n                                                         0.268929      -2.585134\n                                                         0.307376       1.223972\n                                                         1.729730       0.507606\n                                                         2.119890      -0.718466\n                                                         2.793100       0.354338\n                                                         6.883380      -0.098686\n                                                         9.807220       0.099553\n                                                        13.505700      -0.052746\n                                                        18.244300       0.013276\n                                                        39.342200      -0.000458\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n                          3       1    3dz2              0.238188       0.888463\n                                                         0.268929      -1.492528\n                                                         0.307376       0.706661\n                                                         1.729730       0.293067\n                                                         2.119890      -0.414806\n                                                         2.793100       0.204577\n                                                         6.883380      -0.056976\n                                                         9.807220       0.057477\n                                                        13.505700      -0.030453\n                                                        18.244300       0.007665\n                                                        39.342200      -0.000264\n                                                       122.198000       0.000000\n                                                       310.401000       0.000000\n                                                       526.124000       0.000000\n\n     GTH Potential information for                                    GTH-PBE-q2\n\n       Description:                       Goedecker-Teter-Hutter pseudopotential\n                                           Goedecker et al., PRB 54, 1703 (1996)\n                                          Hartwigsen et al., PRB 58, 3641 (1998)\n                                                      Krack, TCA 114, 145 (2005)\n\n       Gaussian exponent of the core charge distribution:               1.502029\n       Electronic configuration (s p d ...):                                   2\n\n       Parameters of the local part of the GTH pseudopotential:\n\n                          rloc        C1          C2          C3          C4\n                        0.576960   -2.690407\n\n       Parameters of the non-local part of the GTH pseudopotential:\n\n                   l      r(l)      h(i,j,l)\n\n                   0    0.593924    3.503211   -0.716772\n                                   -0.716772    0.925348\n                   1    0.707157    0.831158\n\n\n MOLECULE KIND INFORMATION\n\n\n All atoms are their own molecule, skipping detailed information\n\n\n TOTAL NUMBERS AND MAXIMUM NUMBERS\n\n  Total number of            - Atomic kinds:                                   1\n                             - Atoms:                                          2\n                             - Shell sets:                                     6\n                             - Shells:                                         6\n                             - Primitive Cartesian functions:                 84\n                             - Cartesian basis functions:                     20\n                             - Spherical basis functions:                     18\n\n  Maximum angular momentum of- Orbital basis functions:                        2\n                             - Local part of the GTH pseudopotential:          0\n                             - Non-local part of the GTH pseudopotential:      2\n\n\n MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n\n  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n\n       1     1 Mg  12    0.000000    0.000000    0.000000      2.00      24.3050\n       2     1 Mg  12    0.000018    1.851434    2.603299      2.00      24.3050\n\n\n\n\n SCF PARAMETERS         Density guess:                                    ATOMIC\n                        --------------------------------------------------------\n                        max_scf:                                               1\n                        max_scf_history:                                       0\n                        max_diis:                                              4\n                        --------------------------------------------------------\n                        eps_scf:                                        1.00E-07\n                        eps_scf_history:                                0.00E+00\n                        eps_diis:                                       1.00E-01\n                        eps_eigval:                                     1.00E-05\n                        --------------------------------------------------------\n                        level_shift [a.u.]:                                 0.00\n                        added MOs                                         4    0\n                        --------------------------------------------------------\n                        Mixing method:                            BROYDEN_MIXING\n                                                charge density mixing in g-space\n                        --------------------------------------------------------\n                        Smear method:                                FERMI_DIRAC\n                        Electronic temperature [K]:                        157.9\n                        Electronic temperature [a.u.]:                  5.00E-04\n                        Accuracy threshold:                             1.00E-10\n                        --------------------------------------------------------\n                        No outer SCF\n\n PW_GRID| Information for grid number                                          1\n PW_GRID| Grid distributed over                                    12 processors\n PW_GRID| Real space group dimensions                                    12    1\n PW_GRID| the grid is blocked:                                                NO\n PW_GRID| Cutoff [a.u.]                                                    220.5\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1            -22      22                Points:          45\n PW_GRID|   Bounds   2            -22      22                Points:          45\n PW_GRID|   Bounds   3            -36      35                Points:          72\n PW_GRID| Volume element (a.u.^3)  0.2146E-02     Volume (a.u.^3)       312.9133\n PW_GRID| Grid span                                                    FULLSPACE\n PW_GRID|   Distribution                         Average         Max         Min\n PW_GRID|   G-Vectors                            12150.0       12150       12150\n PW_GRID|   G-Rays                                 270.0         270         270\n PW_GRID|   Real Space Points                    12150.0       12960        9720\n\n PW_GRID| Information for grid number                                          2\n PW_GRID| Grid distributed over                                    12 processors\n PW_GRID| Real space group dimensions                                    12    1\n PW_GRID| the grid is blocked:                                                NO\n PW_GRID| Cutoff [a.u.]                                                     73.5\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1            -12      11                Points:          24\n PW_GRID|   Bounds   2            -12      11                Points:          24\n PW_GRID|   Bounds   3            -20      19                Points:          40\n PW_GRID| Volume element (a.u.^3)  0.1358E-01     Volume (a.u.^3)       312.9133\n PW_GRID| Grid span                                                    FULLSPACE\n PW_GRID|   Distribution                         Average         Max         Min\n PW_GRID|   G-Vectors                             1920.0        1944        1896\n PW_GRID|   G-Rays                                  80.0          81          79\n PW_GRID|   Real Space Points                     1920.0        1920        1920\n\n PW_GRID| Information for grid number                                          3\n PW_GRID| Grid distributed over                                    12 processors\n PW_GRID| Real space group dimensions                                    12    1\n PW_GRID| the grid is blocked:                                                NO\n PW_GRID| Cutoff [a.u.]                                                     24.5\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -7       7                Points:          15\n PW_GRID|   Bounds   2             -7       7                Points:          15\n PW_GRID|   Bounds   3            -12      11                Points:          24\n PW_GRID| Volume element (a.u.^3)  0.5795E-01     Volume (a.u.^3)       312.9133\n PW_GRID| Grid span                                                    FULLSPACE\n PW_GRID|   Distribution                         Average         Max         Min\n PW_GRID|   G-Vectors                              450.0         495         405\n PW_GRID|   G-Rays                                  30.0          33          27\n PW_GRID|   Real Space Points                      450.0         720         360\n\n PW_GRID| Information for grid number                                          4\n PW_GRID| Grid distributed over                                    12 processors\n PW_GRID| Real space group dimensions                                     3    4\n PW_GRID| the grid is blocked:                                                NO\n PW_GRID| Cutoff [a.u.]                                                      8.2\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -4       3                Points:           8\n PW_GRID|   Bounds   2             -4       3                Points:           8\n PW_GRID|   Bounds   3             -7       7                Points:          15\n PW_GRID| Volume element (a.u.^3)  0.3260         Volume (a.u.^3)       312.9133\n PW_GRID| Grid span                                                    FULLSPACE\n PW_GRID|   Distribution                         Average         Max         Min\n PW_GRID|   G-Vectors                               80.0          96          72\n PW_GRID|   G-Rays                                  10.0          12           9\n PW_GRID|   Real Space Points                       80.0          90          60\n\n POISSON| Solver                                                        PERIODIC\n POISSON| Periodicity                                                        XYZ\n\n RS_GRID| Information for grid number                                          1\n RS_GRID|   Bounds   1            -22      22                Points:          45\n RS_GRID|   Bounds   2            -22      22                Points:          45\n RS_GRID|   Bounds   3            -36      35                Points:          72\n RS_GRID| Real space fully replicated\n RS_GRID| Group size                                                           1\n\n RS_GRID| Information for grid number                                          2\n RS_GRID|   Bounds   1            -12      11                Points:          24\n RS_GRID|   Bounds   2            -12      11                Points:          24\n RS_GRID|   Bounds   3            -20      19                Points:          40\n RS_GRID| Real space fully replicated\n RS_GRID| Group size                                                           1\n\n RS_GRID| Information for grid number                                          3\n RS_GRID|   Bounds   1             -7       7                Points:          15\n RS_GRID|   Bounds   2             -7       7                Points:          15\n RS_GRID|   Bounds   3            -12      11                Points:          24\n RS_GRID| Real space fully replicated\n RS_GRID| Group size                                                           1\n\n RS_GRID| Information for grid number                                          4\n RS_GRID|   Bounds   1             -4       3                Points:           8\n RS_GRID|   Bounds   2             -4       3                Points:           8\n RS_GRID|   Bounds   3             -7       7                Points:          15\n RS_GRID| Real space fully replicated\n RS_GRID| Group size                                                           1\n \n KPOINTS| Number of kpoint groups                                             12\n KPOINTS| Size of each kpoint group                                            1\n KPOINTS| Number of kpoints per group                                         25\n\n Number of electrons:                                                          4\n Number of occupied orbitals:                                                  2\n Number of molecular orbitals:                                                 6\n\n Number of orbital functions:                                                 18\n Number of independent orbital functions:                                     18\n\n Extrapolation method: initial_guess\n\n Atomic guess: The first density matrix is obtained in terms of atomic orbitals\n               and electronic configurations assigned to each atomic kind\n\n Guess for atomic kind: Mg\n\n Electronic structure\n    Total number of core electrons                                         10.00\n    Total number of valence electrons                                       2.00\n    Total number of electrons                                              12.00\n    Multiplicity                                                   not specified\n    S   [  2.00  2.00] 2.00\n    P   [  6.00]\n \n\n *******************************************************************************\n                  Iteration          Convergence                     Energy [au]\n *******************************************************************************\n                          1         0.00000                      -0.700690874452\n\n Energy components [Hartree]           Total Energy ::           -0.700690874452\n                                        Band Energy ::           -0.578994953267\n                                     Kinetic Energy ::            0.460313113317\n                                   Potential Energy ::           -1.161003987768\n                                      Virial (-V/T) ::            2.522204895279\n                                        Core Energy ::           -1.026242406922\n                                          XC Energy ::           -0.392219170285\n                                     Coulomb Energy ::            0.717770702756\n                       Total Pseudopotential Energy ::           -1.518980860266\n                       Local Pseudopotential Energy ::           -1.829807554801\n                    Nonlocal Pseudopotential Energy ::            0.310826694535\n                                        Confinement ::            0.324253400279\n\n Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]\n\n                       1     0          2.000      -0.289497           -7.877627\n \n\n Total Electron Density at R=0:                                         0.000012\n Re-scaling the density matrix to get the right number of electrons\n                  # Electrons              Trace(P)               Scaling factor\n                            4                 4.000                        1.000\n\n\n SCF WAVEFUNCTION OPTIMIZATION\n\n  Step     Update method      Time    Convergence         Total energy    Change\n  ------------------------------------------------------------------------------\n     1 NoMix/Diag. 0.40E+00   25.3     1.42568183        -1.5706574484 -1.57E+00\n\n  Leaving inner SCF loop after reaching     1 steps.\n\n\n  Electronic density on regular grids:         -4.0000000107       -0.0000000107\n  Core density on regular grids:                3.9999999968       -0.0000000032\n  Total charge density on r-space grids:       -0.0000000140\n  Total charge density g-space grids:          -0.0000000140\n\n  Overlap energy of the core charge distribution:               0.00000000000099\n  Self energy of the core charge distribution:                 -3.91146295972394\n  Core Hamiltonian energy:                                      1.40535228093216\n  Hartree energy:                                               1.80852671945613\n  Exchange-correlation energy:                                 -0.87307348906342\n  Electronic entropic energy:                                  -0.00001885713412\n  Fermi energy:                                                 0.10974707283205\n\n  Total energy:                                                -1.57065744839809\n\n *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***\n\n\n !-----------------------------------------------------------------------------!\n                     Mulliken Population Analysis\n\n #  Atom  Element  Kind  Atomic population                           Net charge\n       1     Mg       1          2.000000                              0.000000\n       2     Mg       1          2.000000                             -0.000000\n # Total charge                              4.000000                 -0.000000\n\n !-----------------------------------------------------------------------------!\n\n !-----------------------------------------------------------------------------!\n                           Hirshfeld Charges\n\n  #Atom  Element  Kind  Ref Charge     Population                    Net charge\n      1       Mg     1       2.000          2.000                         0.000\n      2       Mg     1       2.000          2.000                        -0.000\n\n  Total Charge                                                           -0.000\n !-----------------------------------------------------------------------------!\n\n ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.570657448398086\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                DBCSR STATISTICS                             -\n -                                                                             -\n -------------------------------------------------------------------------------\n COUNTER                                    TOTAL       BLAS       SMM       ACC\n flops inhomo. stacks                           0       0.0%      0.0%      0.0%\n flops total                         0.000000E+00       0.0%      0.0%      0.0%\n flops max/rank                      0.000000E+00       0.0%      0.0%      0.0%\n matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%\n matmuls total                                  0       0.0%      0.0%      0.0%\n number of processed stacks                     0       0.0%      0.0%      0.0%\n average stack size                                     0.0       0.0       0.0\n marketing flops                     0.000000E+00\n -------------------------------------------------------------------------------\n \n MEMORY| Estimated peak process memory [MiB]                                  83\n\n -------------------------------------------------------------------------------\n ----                             MULTIGRID INFO                            ----\n -------------------------------------------------------------------------------\n count for grid        1:         122610          cutoff [a.u.]          220.50\n count for grid        2:              0          cutoff [a.u.]           73.50\n count for grid        3:              0          cutoff [a.u.]           24.50\n count for grid        4:              0          cutoff [a.u.]            8.17\n total gridlevel count  :         122610\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                         MESSAGE PASSING PERFORMANCE                         -\n -                                                                             -\n -------------------------------------------------------------------------------\n\n ROUTINE             CALLS      AVE VOLUME [Bytes]\n MP_Group             1288\n MP_Bcast            10830                  10561.\n MP_Allreduce          357                    212.\n MP_Sync                 4\n MP_Alltoall         25611\n MP_ISendRecv          165                  45840.\n MP_Wait              7515\n MP_ISend             7650                     66.\n MP_IRecv            12600                     43.\n MP_Recv               410                    432.\n -------------------------------------------------------------------------------\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                           R E F E R E N C E S                               -\n -                                                                             -\n -------------------------------------------------------------------------------\n \n CP2K version 6.1, the CP2K developers group (2018).\n CP2K is freely available from https://www.cp2k.org/ .\n\n Schuett, Ole; Messmer, Peter; Hutter, Juerg; VandeVondele, Joost. \n Electronic Structure Calculations on Graphics Processing Units, John Wiley & Sons, Ltd, 173-190 (2016). \n GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory.\n http://dx.doi.org/10.1002/9781118670712.ch8\n\n\n Borstnik, U; VandeVondele, J; Weber, V; Hutter, J. \n PARALLEL COMPUTING, 40 (5-6), 47-58 (2014). \n Sparse matrix multiplication: The distributed block-compressed sparse\n row library.\n http://dx.doi.org/10.1016/j.parco.2014.03.012\n\n\n Hutter, J; Iannuzzi, M; Schiffmann, F; VandeVondele, J. \n WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE, 4 (1), 15-25 (2014). \n CP2K: atomistic simulations of condensed matter systems.\n http://dx.doi.org/10.1002/wcms.1159\n\n\n Krack, M. \n THEORETICAL CHEMISTRY ACCOUNTS, 114 (1-3), 145-152 (2005). \n Pseudopotentials for H to Kr optimized for gradient-corrected\n exchange-correlation functionals.\n http://dx.doi.org/10.1007/s00214-005-0655-y\n\n\n VandeVondele, J; Krack, M; Mohamed, F; Parrinello, M; Chassaing, T;\n Hutter, J. COMPUTER PHYSICS COMMUNICATIONS, 167 (2), 103-128 (2005). \n QUICKSTEP: Fast and accurate density functional calculations using a\n mixed Gaussian and plane waves approach.\n http://dx.doi.org/10.1016/j.cpc.2004.12.014\n\n\n Frigo, M; Johnson, SG. \n PROCEEDINGS OF THE IEEE, 93 (2), 216-231 (2005). \n The design and implementation of FFTW3.\n http://dx.doi.org/10.1109/JPROC.2004.840301\n\n\n Hartwigsen, C; Goedecker, S; Hutter, J. \n PHYSICAL REVIEW B, 58 (7), 3641-3662 (1998). \n Relativistic separable dual-space Gaussian pseudopotentials from H to Rn.\n http://dx.doi.org/10.1103/PhysRevB.58.3641\n\n\n Lippert, G; Hutter, J; Parrinello, M. \n MOLECULAR PHYSICS, 92 (3), 477-487 (1997). \n A hybrid Gaussian and plane wave density functional scheme.\n http://dx.doi.org/10.1080/002689797170220\n\n\n Perdew, JP; Burke, K; Ernzerhof, M. \n PHYSICAL REVIEW LETTERS, 77 (18), 3865-3868 (1996). \n Generalized gradient approximation made simple.\n http://dx.doi.org/10.1103/PhysRevLett.77.3865\n\n\n Goedecker, S; Teter, M; Hutter, J. \n PHYSICAL REVIEW B, 54 (3), 1703-1710 (1996). \n Separable dual-space Gaussian pseudopotentials.\n http://dx.doi.org/10.1103/PhysRevB.54.1703\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                T I M I N G                                  -\n -                                                                             -\n -------------------------------------------------------------------------------\n SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME\n                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM\n CP2K                                 1  1.0    0.005    0.011   58.640   58.641\n qs_energies                          1  2.0    0.001    0.001   58.547   58.547\n scf_env_do_scf                       1  3.0    0.000    0.000   57.776   57.777\n scf_env_do_scf_inner_loop            1  4.0    0.000    0.000   57.776   57.776\n rs_pw_transfer                      15  8.5    0.000    0.000   28.842   32.983\n mp_waitall_1                      7515  8.4   28.840   32.979   28.840   32.979\n rs_pw_transfer_RS2PW_230             4  8.2    0.005    0.005   28.837   32.977\n qs_rho_update_rho                    2  5.0    0.000    0.000   32.822   32.822\n calculate_rho_elec                   2  6.0    4.116   32.802   32.822   32.822\n density_rs2pw                        2  7.0    0.000    0.000   28.704   32.820\n qs_scf_new_mos_kp                    1  5.0    0.000    0.000   22.296   25.300\n do_general_diag_kp                   1  6.0    0.011    0.016   22.296   25.300\n dbcsr_desymmetrize_deep           2482  7.5    0.062    0.068   21.307   24.328\n mp_alltoall_i22                   4364  8.8   21.223   24.272   21.223   24.272\n qs_ks_update_qs_env                  1  5.0    0.000    0.000    3.048   23.932\n rebuild_ks_matrix                    1  6.0    0.000    0.000    3.046   23.929\n qs_ks_build_kohn_sham_matrix         1  7.0    0.000    0.001    3.046   23.929\n sum_up_and_integrate                 1  8.0    0.000    0.000    3.016   23.899\n integrate_v_rspace                   1  9.0    3.005   23.888    3.016   23.899\n -------------------------------------------------------------------------------\n\n The number of warnings for this run is : 2\n \n -------------------------------------------------------------------------------\n  **** **** ******  **  PROGRAM ENDED AT                 2018-09-13 10:46:37.882\n ***** ** ***  *** **   PROGRAM RAN ON             cx1-138-14-1.cx1.hpc.ic.ac.uk\n **    ****   ******    PROGRAM RAN BY                                     rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 14076\n  **** **  *******  **  PROGRAM STOPPED IN /tmp/pbs.1998026.cx1/workDir/cutoff_6\n                                           000\n"

	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createCP2KFullOutFileB_kPts_moPrint():
	filePath = os.path.abspath( os.path.join( os.getcwd(), "full_file_fcc_mg_moprint_kpts.cpout" ) )
	fileStr = " DBCSR| Multiplication driver                                               BLAS\n DBCSR| Multrec recursion limit                                              512\n DBCSR| Multiplication stack size                                           1000\n DBCSR| Maximum elements for images                                    UNLIMITED\n DBCSR| Multiplicative factor virtual images                                   1\n DBCSR| Multiplication size stacks                                             3\n\n\n  **** **** ******  **  PROGRAM STARTED AT               2019-05-17 10:21:32.420\n ***** ** ***  *** **   PROGRAM STARTED ON                              mt-rf614\n **    ****   ******    PROGRAM STARTED BY                                 rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 36923\n  **** **  *******  **  PROGRAM STARTED IN /media/ssd1/rf614/Work/Documents/jobs\n                                           /Corrosion_Work/random/learning_cp2k/\n                                           get_dos\n\n CP2K| version string:                                          CP2K version 6.1\n CP2K| source code revision number:                                    svn:18464\n CP2K| cp2kflags: libint fftw3 libxc libderiv_max_am1=5 libint_max_am=7 max_cont\n CP2K|            r=4\n CP2K| is freely available from                            https://www.cp2k.org/\n CP2K| Program compiled at                          Thu 13 Sep 09:07:18 BST 2018\n CP2K| Program compiled on                                              mt-rf614\n CP2K| Program compiled for                                                local\n CP2K| Data directory path             /media/ssd1/rf614/Work/CP2K/cp2k-6.1/data\n CP2K| Input file name                                               rel_600.inp\n\n GLOBAL| Force Environment number                                              1\n GLOBAL| Basis set file name                                 TIGHT_BINDING_BASIS\n GLOBAL| Potential file name                                      GTH_POTENTIALS\n GLOBAL| MM Potential file name                                     MM_POTENTIAL\n GLOBAL| Coordinate file name                                      __STD_INPUT__\n GLOBAL| Method name                                                        CP2K\n GLOBAL| Project name                                                    rel_600\n GLOBAL| Preferred FFT library                                             FFTW3\n GLOBAL| Preferred diagonalization lib.                                       SL\n GLOBAL| Run type                                                         ENERGY\n GLOBAL| All-to-all communication in single precision                          F\n GLOBAL| FFTs using library dependent lengths                                  F\n GLOBAL| Global print level                                               MEDIUM\n GLOBAL| Total number of message passing processes                             1\n GLOBAL| Number of threads for this process                                    1\n GLOBAL| This output is from process                                           0\n GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-1620 v4 @ 3.50GHz\n\n MEMORY| system memory details [Kb]\n MEMORY|                        rank 0           min           max       average\n MEMORY| MemTotal             65865396      65865396      65865396      65865396\n MEMORY| MemFree              48782684      48782684      48782684      48782684\n MEMORY| Buffers                679500        679500        679500        679500\n MEMORY| Cached                9912856       9912856       9912856       9912856\n MEMORY| Slab                   776260        776260        776260        776260\n MEMORY| SReclaimable           643112        643112        643112        643112\n MEMORY| MemLikelyFree        60018152      60018152      60018152      60018152\n\n\n *** Fundamental physical constants (SI units) ***\n\n *** Literature: B. J. Mohr and B. N. Taylor,\n ***             CODATA recommended values of the fundamental physical\n ***             constants: 2006, Web Version 5.1\n ***             http://physics.nist.gov/constants\n\n Speed of light in vacuum [m/s]                             2.99792458000000E+08\n Magnetic constant or permeability of vacuum [N/A**2]       1.25663706143592E-06\n Electric constant or permittivity of vacuum [F/m]          8.85418781762039E-12\n Planck constant (h) [J*s]                                  6.62606896000000E-34\n Planck constant (h-bar) [J*s]                              1.05457162825177E-34\n Elementary charge [C]                                      1.60217648700000E-19\n Electron mass [kg]                                         9.10938215000000E-31\n Electron g factor [ ]                                     -2.00231930436220E+00\n Proton mass [kg]                                           1.67262163700000E-27\n Fine-structure constant                                    7.29735253760000E-03\n Rydberg constant [1/m]                                     1.09737315685270E+07\n Avogadro constant [1/mol]                                  6.02214179000000E+23\n Boltzmann constant [J/K]                                   1.38065040000000E-23\n Atomic mass unit [kg]                                      1.66053878200000E-27\n Bohr radius [m]                                            5.29177208590000E-11\n\n *** Conversion factors ***\n\n [u] -> [a.u.]                                              1.82288848426455E+03\n [Angstrom] -> [Bohr] = [a.u.]                              1.88972613288564E+00\n [a.u.] = [Bohr] -> [Angstrom]                              5.29177208590000E-01\n [a.u.] -> [s]                                              2.41888432650478E-17\n [a.u.] -> [fs]                                             2.41888432650478E-02\n [a.u.] -> [J]                                              4.35974393937059E-18\n [a.u.] -> [N]                                              8.23872205491840E-08\n [a.u.] -> [K]                                              3.15774647902944E+05\n [a.u.] -> [kJ/mol]                                         2.62549961709828E+03\n [a.u.] -> [kcal/mol]                                       6.27509468713739E+02\n [a.u.] -> [Pa]                                             2.94210107994716E+13\n [a.u.] -> [bar]                                            2.94210107994716E+08\n [a.u.] -> [atm]                                            2.90362800883016E+08\n [a.u.] -> [eV]                                             2.72113838565563E+01\n [a.u.] -> [Hz]                                             6.57968392072181E+15\n [a.u.] -> [1/cm] (wave numbers)                            2.19474631370540E+05\n [a.u./Bohr**2] -> [1/cm]                                   5.14048714338585E+03\n \n\n CELL_TOP| Volume [angstrom^3]:                                           23.193\n CELL_TOP| Vector a [angstrom     3.201     0.000     0.000    |a| =       3.201\n CELL_TOP| Vector b [angstrom     1.601     2.772     0.000    |b| =       3.201\n CELL_TOP| Vector c [angstrom     1.601     0.924     2.614    |c| =       3.201\n CELL_TOP| Angle (b,c), alpha [degree]:                                   60.000\n CELL_TOP| Angle (a,c), beta  [degree]:                                   60.000\n CELL_TOP| Angle (a,b), gamma [degree]:                                   60.000\n CELL_TOP| Numerically orthorhombic:                                          NO\n\n GENERATE|  Preliminary Number of Bonds generated:                             0\n GENERATE|  Achieved consistency in connectivity generation.\n\n CELL| Volume [angstrom^3]:                                               23.193\n CELL| Vector a [angstrom]:       3.201     0.000     0.000    |a| =       3.201\n CELL| Vector b [angstrom]:       1.601     2.772     0.000    |b| =       3.201\n CELL| Vector c [angstrom]:       1.601     0.924     2.614    |c| =       3.201\n CELL| Angle (b,c), alpha [degree]:                                       60.000\n CELL| Angle (a,c), beta  [degree]:                                       60.000\n CELL| Angle (a,b), gamma [degree]:                                       60.000\n CELL| Numerically orthorhombic:                                              NO\n\n CELL_REF| Volume [angstrom^3]:                                           23.193\n CELL_REF| Vector a [angstrom     3.201     0.000     0.000    |a| =       3.201\n CELL_REF| Vector b [angstrom     1.601     2.772     0.000    |b| =       3.201\n CELL_REF| Vector c [angstrom     1.601     0.924     2.614    |c| =       3.201\n CELL_REF| Angle (b,c), alpha [degree]:                                   60.000\n CELL_REF| Angle (a,c), beta  [degree]:                                   60.000\n CELL_REF| Angle (a,b), gamma [degree]:                                   60.000\n CELL_REF| Numerically orthorhombic:                                          NO\n\n *** WARNING in cryssym.F:163 :: Symmetry library SPGLIB not available ***\n\n\n *******************************************************************************\n                                    Kpoints\n *******************************************************************************\n BRILLOUIN| K-point scheme                                        Monkhorst-Pack\n BRILLOUIN| K-Point grid                                             2    2    1\n BRILLOUIN| Accuracy in Symmetry determination                      0.100000E-05\n BRILLOUIN| K-Point point group symmetrization                               OFF\n BRILLOUIN| Wavefunction type                                            COMPLEX\n BRILLOUIN| List of Kpoints [2 Pi/Bohr]                                        2\n BRILLOUIN| Number           Weight            X              Y              Z\n BRILLOUIN|     1           0.50000       -0.25000       -0.25000        0.00000\n BRILLOUIN|     2           0.50000       -0.25000        0.25000        0.00000\n *******************************************************************************\n\n *******************************************************************************\n *******************************************************************************\n **                                                                           **\n **     #####                         ##              ##                      **\n **    ##   ##            ##          ##              ##                      **\n **   ##     ##                       ##            ######                    **\n **   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **\n **   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **\n **   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **\n **    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **\n **     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **\n **           ##                                                    ##        **\n **                                                                           **\n **                                                ... make the atoms dance   **\n **                                                                           **\n **            Copyright (C) by CP2K developers group (2000 - 2018)           **\n **                                                                           **\n *******************************************************************************\n DFT| Spin restricted Kohn-Sham (RKS) calculation                            RKS\n DFT| Multiplicity                                                             1\n DFT| Number of spin states                                                    1\n DFT| Charge                                                                   0\n DFT| Self-interaction correction (SIC)                                       NO\n DFT| Cutoffs: density                                              1.000000E-10\n DFT|          gradient                                             1.000000E-10\n DFT|          tau                                                  1.000000E-10\n DFT|          cutoff_smoothing_range                               0.000000E+00\n DFT| XC density smoothing                                                  NONE\n DFT| XC derivatives                                                          PW\n FUNCTIONAL| ROUTINE=NEW\n FUNCTIONAL| PBE:\n FUNCTIONAL| J.P.Perdew, K.Burke, M.Ernzerhof, Phys. Rev. Letter, vol. 77, n 18,\n FUNCTIONAL|  pp. 3865-3868, (1996){spin unpolarized}                           \n\n QS| Method:                                                                 GPW\n QS| Density plane wave grid type                        NON-SPHERICAL FULLSPACE\n QS| Number of grid levels:                                                    4\n QS| Density cutoff [a.u.]:                                                 18.4\n QS| Multi grid cutoff [a.u.]: 1) grid level                                18.4\n QS|                           2) grid level                                 6.1\n QS|                           3) grid level                                 2.0\n QS|                           4) grid level                                 0.7\n QS| Grid level progression factor:                                          3.0\n QS| Relative density cutoff [a.u.]:                                        22.0\n QS| Consistent realspace mapping and integration \n QS| Interaction thresholds: eps_pgf_orb:                                1.0E-05\n QS|                         eps_filter_matrix:                          0.0E+00\n QS|                         eps_core_charge:                            1.0E-12\n QS|                         eps_rho_gspace:                             1.0E-10\n QS|                         eps_rho_rspace:                             1.0E-10\n QS|                         eps_gvg_rspace:                             1.0E-05\n QS|                         eps_ppl:                                    1.0E-02\n QS|                         eps_ppnl:                                   1.0E-07\n\n\n ATOMIC KIND INFORMATION\n\n  1. Atomic kind: Mg                                    Number of atoms:       1\n\n     Orbital Basis Set                                 spd-2z-2setsExp-5exp-6exp\n\n       Number of orbital shell sets:                                           2\n       Number of orbital shells:                                               6\n       Number of primitive Cartesian functions:                               11\n       Number of Cartesian basis functions:                                   20\n       Number of spherical basis functions:                                   18\n       Norm type:                                                              2\n\n       Normalised Cartesian orbitals:\n\n                        Set   Shell   Orbital            Exponent    Coefficient\n\n                          1       1    3s                0.172856       1.423308\n                                                         0.208729      -2.344040\n                                                         0.255247       1.168454\n                                                         0.527695      -0.219714\n                                                         1.699878      -0.026291\n\n                          1       2    4px               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n                          1       2    4py               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n                          1       2    4pz               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n\n                          1       3    5dx2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n                          1       3    5dxy              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dxz              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dy2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n                          1       3    5dyz              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dz2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n\n                          2       1    3s                0.220648       0.765795\n                                                         0.396852      -3.213702\n                                                         0.480094     -11.290279\n                                                         0.578055      33.333097\n                                                         0.703605     -27.251732\n                                                         0.862658       8.390367\n\n                          2       2    4px               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n                          2       2    4py               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n                          2       2    4pz               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n\n                          2       3    5dx2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n                          2       3    5dxy              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dxz              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dy2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n                          2       3    5dyz              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dz2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n\n     GTH Potential information for                                    GTH-PBE-q2\n\n       Description:                       Goedecker-Teter-Hutter pseudopotential\n                                           Goedecker et al., PRB 54, 1703 (1996)\n                                          Hartwigsen et al., PRB 58, 3641 (1998)\n                                                      Krack, TCA 114, 145 (2005)\n\n       Gaussian exponent of the core charge distribution:               1.502029\n       Electronic configuration (s p d ...):                                   2\n\n       Parameters of the local part of the GTH pseudopotential:\n\n                          rloc        C1          C2          C3          C4\n                        0.576960   -2.690407\n\n       Parameters of the non-local part of the GTH pseudopotential:\n\n                   l      r(l)      h(i,j,l)\n\n                   0    0.593924    3.503211   -0.716772\n                                   -0.716772    0.925348\n                   1    0.707157    0.831158\n\n\n MOLECULE KIND INFORMATION\n\n\n All atoms are their own molecule, skipping detailed information\n\n\n TOTAL NUMBERS AND MAXIMUM NUMBERS\n\n  Total number of            - Atomic kinds:                                   1\n                             - Atoms:                                          1\n                             - Shell sets:                                     2\n                             - Shells:                                         6\n                             - Primitive Cartesian functions:                 11\n                             - Cartesian basis functions:                     20\n                             - Spherical basis functions:                     18\n\n  Maximum angular momentum of- Orbital basis functions:                        2\n                             - Local part of the GTH pseudopotential:          0\n                             - Non-local part of the GTH pseudopotential:      2\n\n\n MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n\n  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n\n       1     1 Mg  12    0.000000    0.000000    0.000000      2.00      24.3050\n\n\n\n\n SCF PARAMETERS         Density guess:                                    ATOMIC\n                        --------------------------------------------------------\n                        max_scf:                                               1\n                        max_scf_history:                                       0\n                        max_diis:                                              4\n                        --------------------------------------------------------\n                        eps_scf:                                        1.00E-07\n                        eps_scf_history:                                0.00E+00\n                        eps_diis:                                       1.00E-01\n                        eps_eigval:                                     1.00E-05\n                        --------------------------------------------------------\n                        level_shift [a.u.]:                                 0.00\n                        added MOs                                         4    0\n                        --------------------------------------------------------\n                        Mixing method:                            BROYDEN_MIXING\n                                                charge density mixing in g-space\n                        --------------------------------------------------------\n                        Smear method:                                FERMI_DIRAC\n                        Electronic temperature [K]:                        157.9\n                        Electronic temperature [a.u.]:                  5.00E-04\n                        Accuracy threshold:                             1.00E-10\n                        --------------------------------------------------------\n                        No outer SCF\n\n PW_GRID| Information for grid number                                          1\n PW_GRID| Cutoff [a.u.]                                                     18.4\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -6       5                Points:          12\n PW_GRID|   Bounds   2             -6       5                Points:          12\n PW_GRID|   Bounds   3             -6       5                Points:          12\n PW_GRID| Volume element (a.u.^3)  0.9057E-01     Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          2\n PW_GRID| Cutoff [a.u.]                                                      6.1\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -4       3                Points:           8\n PW_GRID|   Bounds   2             -4       3                Points:           8\n PW_GRID|   Bounds   3             -4       3                Points:           8\n PW_GRID| Volume element (a.u.^3)  0.3057         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          3\n PW_GRID| Cutoff [a.u.]                                                      2.0\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -2       1                Points:           4\n PW_GRID|   Bounds   2             -2       1                Points:           4\n PW_GRID|   Bounds   3             -2       1                Points:           4\n PW_GRID| Volume element (a.u.^3)   2.446         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          4\n PW_GRID| Cutoff [a.u.]                                                      0.7\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -2       1                Points:           4\n PW_GRID|   Bounds   2             -2       1                Points:           4\n PW_GRID|   Bounds   3             -2       1                Points:           4\n PW_GRID| Volume element (a.u.^3)   2.446         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n POISSON| Solver                                                        PERIODIC\n POISSON| Periodicity                                                        XYZ\n\n RS_GRID| Information for grid number                                          1\n RS_GRID|   Bounds   1             -6       5                Points:          12\n RS_GRID|   Bounds   2             -6       5                Points:          12\n RS_GRID|   Bounds   3             -6       5                Points:          12\n\n RS_GRID| Information for grid number                                          2\n RS_GRID|   Bounds   1             -4       3                Points:           8\n RS_GRID|   Bounds   2             -4       3                Points:           8\n RS_GRID|   Bounds   3             -4       3                Points:           8\n\n RS_GRID| Information for grid number                                          3\n RS_GRID|   Bounds   1             -2       1                Points:           4\n RS_GRID|   Bounds   2             -2       1                Points:           4\n RS_GRID|   Bounds   3             -2       1                Points:           4\n\n RS_GRID| Information for grid number                                          4\n RS_GRID|   Bounds   1             -2       1                Points:           4\n RS_GRID|   Bounds   2             -2       1                Points:           4\n RS_GRID|   Bounds   3             -2       1                Points:           4\n\n Number of electrons:                                                          2\n Number of occupied orbitals:                                                  1\n Number of molecular orbitals:                                                 5\n\n Number of orbital functions:                                                 18\n Number of independent orbital functions:                                     18\n\n Extrapolation method: initial_guess\n\n Atomic guess: The first density matrix is obtained in terms of atomic orbitals\n               and electronic configurations assigned to each atomic kind\n\n Guess for atomic kind: Mg\n\n Electronic structure\n    Total number of core electrons                                         10.00\n    Total number of valence electrons                                       2.00\n    Total number of electrons                                              12.00\n    Multiplicity                                                   not specified\n    S   [  2.00  2.00] 2.00\n    P   [  6.00]\n\n\n *******************************************************************************\n                  Iteration          Convergence                     Energy [au]\n *******************************************************************************\n                          1        0.450866E-02                  -0.735938424170\n                          2        0.173646E-03                  -0.735943231943\n                          3        0.224554E-07                  -0.735943239082\n\n Energy components [Hartree]           Total Energy ::           -0.735943239082\n                                        Band Energy ::           -0.170223083884\n                                     Kinetic Energy ::            0.383696410659\n                                   Potential Energy ::           -1.119639649741\n                                      Virial (-V/T) ::            2.918035245154\n                                        Core Energy ::           -1.043063546624\n                                          XC Energy ::           -0.369655695440\n                                     Coulomb Energy ::            0.676776002981\n                       Total Pseudopotential Energy ::           -1.463503440618\n                       Local Pseudopotential Energy ::           -1.742048485333\n                    Nonlocal Pseudopotential Energy ::            0.278545044716\n                                        Confinement ::            0.367434833351\n\n Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]\n\n                       1     0          2.000      -0.085112           -2.316003\n\n\n Total Electron Density at R=0:                                         0.000119\n Re-scaling the density matrix to get the right number of electrons\n                  # Electrons              Trace(P)               Scaling factor\n                            2                 2.000                        1.000\n\n\n SCF WAVEFUNCTION OPTIMIZATION\n\n  Step     Update method      Time    Convergence         Total energy    Change\n  ------------------------------------------------------------------------------\n\n\n MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER SCF STEP 0\n\n# MO index          MO eigenvalue [a.u.]            MO occupation\n         1                      0.000000                 2.000000\n         2                      0.000000                 0.000000\n         3                      0.000000                 0.000000\n         4                      0.000000                 0.000000\n         5                      0.000000                 0.000000\n# Sum                                                    2.000000\n\n  Fermi energy:                 0.000000\n\n  HOMO-LUMO gap:                0.000000 =   0.00 eV\n\n     1 NoMix/Diag. 0.40E+00    0.4     1.37565204        -0.8487613780 -8.49E-01\n\n  Leaving inner SCF loop after reaching     1 steps.\n\n\n  Electronic density on regular grids:         -2.0000000004       -0.0000000004\n  Core density on regular grids:                1.9999999986       -0.0000000014\n  Total charge density on r-space grids:       -0.0000000018\n  Total charge density g-space grids:          -0.0000000018\n\n  Overlap energy of the core charge distribution:               0.00000000000049\n  Self energy of the core charge distribution:                 -1.95573147986197\n  Core Hamiltonian energy:                                      0.60581667193016\n  Hartree energy:                                               0.93214075526534\n  Exchange-correlation energy:                                 -0.43098732537217\n  Electronic entropic energy:                                  -0.00000000000000\n  Fermi energy:                                                 0.02618039273851\n\n  Total energy:                                                -0.84876137803815\n\n *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***\n\n\n\n MO EIGENVALUES AND MO OCCUPATION NUMBERS FOR K-POINT:     1\n\n# MO index          MO eigenvalue [a.u.]            MO occupation\n         1                     -0.079831                 2.000000\n         2                      0.408251                 0.000000\n         3                      0.416321                 0.000000\n         4                      0.416325                 0.000000\n         5                      0.433492                 0.000000\n# Sum                                                    2.000000\n\n  Fermi energy:                 0.026180\n\n  HOMO-LUMO gap:                0.488081 =  13.28 eV\n\n\n\n MO EIGENVALUES AND MO OCCUPATION NUMBERS FOR K-POINT:     2\n\n# MO index          MO eigenvalue [a.u.]            MO occupation\n         1                     -0.013849                 2.000000\n         2                      0.222206                 0.000000\n         3                      0.270688                 0.000000\n         4                      0.475077                 0.000000\n         5                      0.536199                 0.000000\n# Sum                                                    2.000000\n\n  Fermi energy:                 0.026180\n\n  HOMO-LUMO gap:                0.236054 =   6.42 eV\n\n\n *** WARNING in qs_scf_post_gpw.F:2123 :: Projected density of states not ***\n *** implemented for k-points.                                            ***\n\n\n !-----------------------------------------------------------------------------!\n                     Mulliken Population Analysis\n\n #  Atom  Element  Kind  Atomic population                           Net charge\n       1     Mg       1          2.000000                              0.000000\n # Total charge                              2.000000                  0.000000\n\n !-----------------------------------------------------------------------------!\n\n !-----------------------------------------------------------------------------!\n                           Hirshfeld Charges\n\n  #Atom  Element  Kind  Ref Charge     Population                    Net charge\n      1       Mg     1       2.000          2.000                        -0.000\n\n  Total Charge                                                           -0.000\n !-----------------------------------------------------------------------------!\n\n ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.848761378038147\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                DBCSR STATISTICS                             -\n -                                                                             -\n -------------------------------------------------------------------------------\n COUNTER                                    TOTAL       BLAS       SMM       ACC\n flops inhomo. stacks                           0       0.0%      0.0%      0.0%\n flops total                         0.000000E+00       0.0%      0.0%      0.0%\n flops max/rank                      0.000000E+00       0.0%      0.0%      0.0%\n matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%\n matmuls total                                  0       0.0%      0.0%      0.0%\n number of processed stacks                     0       0.0%      0.0%      0.0%\n average stack size                                     0.0       0.0       0.0\n marketing flops                     0.000000E+00\n -------------------------------------------------------------------------------\n\n MEMORY| Estimated peak process memory [MiB]                                  48\n\n -------------------------------------------------------------------------------\n ----                             MULTIGRID INFO                            ----\n -------------------------------------------------------------------------------\n count for grid        1:          11459          cutoff [a.u.]           18.37\n count for grid        2:             48          cutoff [a.u.]            6.12\n count for grid        3:              0          cutoff [a.u.]            2.04\n count for grid        4:              0          cutoff [a.u.]            0.68\n total gridlevel count  :          11507\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                           R E F E R E N C E S                               -\n -                                                                             -\n -------------------------------------------------------------------------------\n \n CP2K version 6.1, the CP2K developers group (2018).\n CP2K is freely available from https://www.cp2k.org/ .\n\n Schuett, Ole; Messmer, Peter; Hutter, Juerg; VandeVondele, Joost. \n Electronic Structure Calculations on Graphics Processing Units, John Wiley & Sons, Ltd, 173-190 (2016). \n GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory.\n http://dx.doi.org/10.1002/9781118670712.ch8\n\n\n Borstnik, U; VandeVondele, J; Weber, V; Hutter, J. \n PARALLEL COMPUTING, 40 (5-6), 47-58 (2014). \n Sparse matrix multiplication: The distributed block-compressed sparse\n row library.\n http://dx.doi.org/10.1016/j.parco.2014.03.012\n\n\n Hutter, J; Iannuzzi, M; Schiffmann, F; VandeVondele, J. \n WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE, 4 (1), 15-25 (2014). \n CP2K: atomistic simulations of condensed matter systems.\n http://dx.doi.org/10.1002/wcms.1159\n\n\n Krack, M. \n THEORETICAL CHEMISTRY ACCOUNTS, 114 (1-3), 145-152 (2005). \n Pseudopotentials for H to Kr optimized for gradient-corrected\n exchange-correlation functionals.\n http://dx.doi.org/10.1007/s00214-005-0655-y\n\n\n VandeVondele, J; Krack, M; Mohamed, F; Parrinello, M; Chassaing, T;\n Hutter, J. COMPUTER PHYSICS COMMUNICATIONS, 167 (2), 103-128 (2005). \n QUICKSTEP: Fast and accurate density functional calculations using a\n mixed Gaussian and plane waves approach.\n http://dx.doi.org/10.1016/j.cpc.2004.12.014\n\n\n Frigo, M; Johnson, SG. \n PROCEEDINGS OF THE IEEE, 93 (2), 216-231 (2005). \n The design and implementation of FFTW3.\n http://dx.doi.org/10.1109/JPROC.2004.840301\n\n\n Hartwigsen, C; Goedecker, S; Hutter, J. \n PHYSICAL REVIEW B, 58 (7), 3641-3662 (1998). \n Relativistic separable dual-space Gaussian pseudopotentials from H to Rn.\n http://dx.doi.org/10.1103/PhysRevB.58.3641\n\n\n Lippert, G; Hutter, J; Parrinello, M. \n MOLECULAR PHYSICS, 92 (3), 477-487 (1997). \n A hybrid Gaussian and plane wave density functional scheme.\n http://dx.doi.org/10.1080/002689797170220\n\n\n Perdew, JP; Burke, K; Ernzerhof, M. \n PHYSICAL REVIEW LETTERS, 77 (18), 3865-3868 (1996). \n Generalized gradient approximation made simple.\n http://dx.doi.org/10.1103/PhysRevLett.77.3865\n\n\n Goedecker, S; Teter, M; Hutter, J. \n PHYSICAL REVIEW B, 54 (3), 1703-1710 (1996). \n Separable dual-space Gaussian pseudopotentials.\n http://dx.doi.org/10.1103/PhysRevB.54.1703\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                T I M I N G                                  -\n -                                                                             -\n -------------------------------------------------------------------------------\n SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME\n                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM\n CP2K                                 1  1.0    0.003    0.003    1.111    1.111\n qs_energies                          1  2.0    0.001    0.001    1.066    1.066\n scf_env_do_scf                       1  3.0    0.000    0.000    0.860    0.860\n scf_env_do_scf_inner_loop            1  4.0    0.000    0.000    0.860    0.860\n qs_rho_update_rho                    2  5.0    0.000    0.000    0.480    0.480\n calculate_rho_elec                   2  6.0    0.479    0.479    0.480    0.480\n qs_ks_update_qs_env                  1  5.0    0.000    0.000    0.411    0.411\n rebuild_ks_matrix                    1  6.0    0.000    0.000    0.411    0.411\n qs_ks_build_kohn_sham_matrix         1  7.0    0.000    0.000    0.411    0.411\n sum_up_and_integrate                 1  8.0    0.000    0.000    0.409    0.409\n integrate_v_rspace                   1  9.0    0.402    0.402    0.409    0.409\n qs_energies_init_hamiltonians        1  3.0    0.000    0.000    0.135    0.135\n build_core_hamiltonian_matrix        1  4.0    0.001    0.001    0.112    0.112\n init_scf_run                         1  3.0    0.000    0.000    0.060    0.060\n scf_env_initial_rho_setup            1  4.0    0.000    0.000    0.059    0.059\n build_core_ppl                       1  5.0    0.041    0.041    0.041    0.041\n build_overlap_matrix                 1  5.0    0.029    0.029    0.033    0.033\n build_kinetic_matrix                 1  5.0    0.026    0.026    0.029    0.029\n qs_env_update_s_mstruct              1  4.0    0.000    0.000    0.022    0.022\n -------------------------------------------------------------------------------\n\n The number of warnings for this run is : 3\n \n -------------------------------------------------------------------------------\n  **** **** ******  **  PROGRAM ENDED AT                 2019-05-17 10:21:33.592\n ***** ** ***  *** **   PROGRAM RAN ON                                  mt-rf614\n **    ****   ******    PROGRAM RAN BY                                     rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 36923\n  **** **  *******  **  PROGRAM STOPPED IN /media/ssd1/rf614/Work/Documents/jobs\n                                           /Corrosion_Work/random/learning_cp2k/\n                                           get_dos\n"

	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath



def createCP2KFullOutFileC_noKPts_moPrint():
	filePath = os.path.abspath( os.path.join( os.getcwd(), "full_file_fcc_mg_moprint_nokpts.cpout" ) )
	fileStr = " DBCSR| Multiplication driver                                               BLAS\n DBCSR| Multrec recursion limit                                              512\n DBCSR| Multiplication stack size                                           1000\n DBCSR| Maximum elements for images                                    UNLIMITED\n DBCSR| Multiplicative factor virtual images                                   1\n DBCSR| Multiplication size stacks                                             3\n\n\n  **** **** ******  **  PROGRAM STARTED AT               2019-05-17 09:50:12.639\n ***** ** ***  *** **   PROGRAM STARTED ON                              mt-rf614\n **    ****   ******    PROGRAM STARTED BY                                 rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 33982\n  **** **  *******  **  PROGRAM STARTED IN /media/ssd1/rf614/Work/Documents/jobs\n                                           /Corrosion_Work/random/learning_cp2k/\n                                           get_dos/gamma_point\n\n CP2K| version string:                                          CP2K version 6.1\n CP2K| source code revision number:                                    svn:18464\n CP2K| cp2kflags: libint fftw3 libxc libderiv_max_am1=5 libint_max_am=7 max_cont\n CP2K|            r=4\n CP2K| is freely available from                            https://www.cp2k.org/\n CP2K| Program compiled at                          Thu 13 Sep 09:07:18 BST 2018\n CP2K| Program compiled on                                              mt-rf614\n CP2K| Program compiled for                                                local\n CP2K| Data directory path             /media/ssd1/rf614/Work/CP2K/cp2k-6.1/data\n CP2K| Input file name                                               rel_600.inp\n\n GLOBAL| Force Environment number                                              1\n GLOBAL| Basis set file name                                 TIGHT_BINDING_BASIS\n GLOBAL| Potential file name                                      GTH_POTENTIALS\n GLOBAL| MM Potential file name                                     MM_POTENTIAL\n GLOBAL| Coordinate file name                                      __STD_INPUT__\n GLOBAL| Method name                                                        CP2K\n GLOBAL| Project name                                                    rel_600\n GLOBAL| Preferred FFT library                                             FFTW3\n GLOBAL| Preferred diagonalization lib.                                       SL\n GLOBAL| Run type                                                         ENERGY\n GLOBAL| All-to-all communication in single precision                          F\n GLOBAL| FFTs using library dependent lengths                                  F\n GLOBAL| Global print level                                               MEDIUM\n GLOBAL| Total number of message passing processes                             1\n GLOBAL| Number of threads for this process                                    1\n GLOBAL| This output is from process                                           0\n GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-1620 v4 @ 3.50GHz\n\n MEMORY| system memory details [Kb]\n MEMORY|                        rank 0           min           max       average\n MEMORY| MemTotal             65865396      65865396      65865396      65865396\n MEMORY| MemFree              48970396      48970396      48970396      48970396\n MEMORY| Buffers                679444        679444        679444        679444\n MEMORY| Cached                9842604       9842604       9842604       9842604\n MEMORY| Slab                   776080        776080        776080        776080\n MEMORY| SReclaimable           642852        642852        642852        642852\n MEMORY| MemLikelyFree        60135296      60135296      60135296      60135296\n\n\n *** Fundamental physical constants (SI units) ***\n\n *** Literature: B. J. Mohr and B. N. Taylor,\n ***             CODATA recommended values of the fundamental physical\n ***             constants: 2006, Web Version 5.1\n ***             http://physics.nist.gov/constants\n\n Speed of light in vacuum [m/s]                             2.99792458000000E+08\n Magnetic constant or permeability of vacuum [N/A**2]       1.25663706143592E-06\n Electric constant or permittivity of vacuum [F/m]          8.85418781762039E-12\n Planck constant (h) [J*s]                                  6.62606896000000E-34\n Planck constant (h-bar) [J*s]                              1.05457162825177E-34\n Elementary charge [C]                                      1.60217648700000E-19\n Electron mass [kg]                                         9.10938215000000E-31\n Electron g factor [ ]                                     -2.00231930436220E+00\n Proton mass [kg]                                           1.67262163700000E-27\n Fine-structure constant                                    7.29735253760000E-03\n Rydberg constant [1/m]                                     1.09737315685270E+07\n Avogadro constant [1/mol]                                  6.02214179000000E+23\n Boltzmann constant [J/K]                                   1.38065040000000E-23\n Atomic mass unit [kg]                                      1.66053878200000E-27\n Bohr radius [m]                                            5.29177208590000E-11\n\n *** Conversion factors ***\n\n [u] -> [a.u.]                                              1.82288848426455E+03\n [Angstrom] -> [Bohr] = [a.u.]                              1.88972613288564E+00\n [a.u.] = [Bohr] -> [Angstrom]                              5.29177208590000E-01\n [a.u.] -> [s]                                              2.41888432650478E-17\n [a.u.] -> [fs]                                             2.41888432650478E-02\n [a.u.] -> [J]                                              4.35974393937059E-18\n [a.u.] -> [N]                                              8.23872205491840E-08\n [a.u.] -> [K]                                              3.15774647902944E+05\n [a.u.] -> [kJ/mol]                                         2.62549961709828E+03\n [a.u.] -> [kcal/mol]                                       6.27509468713739E+02\n [a.u.] -> [Pa]                                             2.94210107994716E+13\n [a.u.] -> [bar]                                            2.94210107994716E+08\n [a.u.] -> [atm]                                            2.90362800883016E+08\n [a.u.] -> [eV]                                             2.72113838565563E+01\n [a.u.] -> [Hz]                                             6.57968392072181E+15\n [a.u.] -> [1/cm] (wave numbers)                            2.19474631370540E+05\n [a.u./Bohr**2] -> [1/cm]                                   5.14048714338585E+03\n \n\n CELL_TOP| Volume [angstrom^3]:                                           23.193\n CELL_TOP| Vector a [angstrom     3.201     0.000     0.000    |a| =       3.201\n CELL_TOP| Vector b [angstrom     1.601     2.772     0.000    |b| =       3.201\n CELL_TOP| Vector c [angstrom     1.601     0.924     2.614    |c| =       3.201\n CELL_TOP| Angle (b,c), alpha [degree]:                                   60.000\n CELL_TOP| Angle (a,c), beta  [degree]:                                   60.000\n CELL_TOP| Angle (a,b), gamma [degree]:                                   60.000\n CELL_TOP| Numerically orthorhombic:                                          NO\n\n GENERATE|  Preliminary Number of Bonds generated:                             0\n GENERATE|  Achieved consistency in connectivity generation.\n\n CELL| Volume [angstrom^3]:                                               23.193\n CELL| Vector a [angstrom]:       3.201     0.000     0.000    |a| =       3.201\n CELL| Vector b [angstrom]:       1.601     2.772     0.000    |b| =       3.201\n CELL| Vector c [angstrom]:       1.601     0.924     2.614    |c| =       3.201\n CELL| Angle (b,c), alpha [degree]:                                       60.000\n CELL| Angle (a,c), beta  [degree]:                                       60.000\n CELL| Angle (a,b), gamma [degree]:                                       60.000\n CELL| Numerically orthorhombic:                                              NO\n\n CELL_REF| Volume [angstrom^3]:                                           23.193\n CELL_REF| Vector a [angstrom     3.201     0.000     0.000    |a| =       3.201\n CELL_REF| Vector b [angstrom     1.601     2.772     0.000    |b| =       3.201\n CELL_REF| Vector c [angstrom     1.601     0.924     2.614    |c| =       3.201\n CELL_REF| Angle (b,c), alpha [degree]:                                   60.000\n CELL_REF| Angle (a,c), beta  [degree]:                                   60.000\n CELL_REF| Angle (a,b), gamma [degree]:                                   60.000\n CELL_REF| Numerically orthorhombic:                                          NO\n\n *******************************************************************************\n *******************************************************************************\n **                                                                           **\n **     #####                         ##              ##                      **\n **    ##   ##            ##          ##              ##                      **\n **   ##     ##                       ##            ######                    **\n **   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **\n **   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **\n **   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **\n **    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **\n **     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **\n **           ##                                                    ##        **\n **                                                                           **\n **                                                ... make the atoms dance   **\n **                                                                           **\n **            Copyright (C) by CP2K developers group (2000 - 2018)           **\n **                                                                           **\n *******************************************************************************\n DFT| Spin restricted Kohn-Sham (RKS) calculation                            RKS\n DFT| Multiplicity                                                             1\n DFT| Number of spin states                                                    1\n DFT| Charge                                                                   0\n DFT| Self-interaction correction (SIC)                                       NO\n DFT| Cutoffs: density                                              1.000000E-10\n DFT|          gradient                                             1.000000E-10\n DFT|          tau                                                  1.000000E-10\n DFT|          cutoff_smoothing_range                               0.000000E+00\n DFT| XC density smoothing                                                  NONE\n DFT| XC derivatives                                                          PW\n FUNCTIONAL| ROUTINE=NEW\n FUNCTIONAL| PBE:\n FUNCTIONAL| J.P.Perdew, K.Burke, M.Ernzerhof, Phys. Rev. Letter, vol. 77, n 18,\n FUNCTIONAL|  pp. 3865-3868, (1996){spin unpolarized}                           \n\n QS| Method:                                                                 GPW\n QS| Density plane wave grid type                        NON-SPHERICAL FULLSPACE\n QS| Number of grid levels:                                                    4\n QS| Density cutoff [a.u.]:                                                 18.4\n QS| Multi grid cutoff [a.u.]: 1) grid level                                18.4\n QS|                           2) grid level                                 6.1\n QS|                           3) grid level                                 2.0\n QS|                           4) grid level                                 0.7\n QS| Grid level progression factor:                                          3.0\n QS| Relative density cutoff [a.u.]:                                        22.0\n QS| Consistent realspace mapping and integration \n QS| Interaction thresholds: eps_pgf_orb:                                1.0E-05\n QS|                         eps_filter_matrix:                          0.0E+00\n QS|                         eps_core_charge:                            1.0E-12\n QS|                         eps_rho_gspace:                             1.0E-10\n QS|                         eps_rho_rspace:                             1.0E-10\n QS|                         eps_gvg_rspace:                             1.0E-05\n QS|                         eps_ppl:                                    1.0E-02\n QS|                         eps_ppnl:                                   1.0E-07\n\n\n ATOMIC KIND INFORMATION\n\n  1. Atomic kind: Mg                                    Number of atoms:       1\n\n     Orbital Basis Set                                 spd-2z-2setsExp-5exp-6exp\n\n       Number of orbital shell sets:                                           2\n       Number of orbital shells:                                               6\n       Number of primitive Cartesian functions:                               11\n       Number of Cartesian basis functions:                                   20\n       Number of spherical basis functions:                                   18\n       Norm type:                                                              2\n\n       Normalised Cartesian orbitals:\n\n                        Set   Shell   Orbital            Exponent    Coefficient\n\n                          1       1    3s                0.172856       1.423308\n                                                         0.208729      -2.344040\n                                                         0.255247       1.168454\n                                                         0.527695      -0.219714\n                                                         1.699878      -0.026291\n\n                          1       2    4px               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n                          1       2    4py               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n                          1       2    4pz               0.172856       0.638775\n                                                         0.208729      -0.966873\n                                                         0.255247       0.482746\n                                                         0.527695      -0.067933\n                                                         1.699878       0.017809\n\n                          1       3    5dx2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n                          1       3    5dxy              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dxz              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dy2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n                          1       3    5dyz              0.172856       0.321859\n                                                         0.208729      -0.369439\n                                                         0.255247       0.136209\n                                                         0.527695       0.068982\n                                                         1.699878       0.063572\n                          1       3    5dz2              0.172856       0.185825\n                                                         0.208729      -0.213296\n                                                         0.255247       0.078640\n                                                         0.527695       0.039827\n                                                         1.699878       0.036703\n\n                          2       1    3s                0.220648       0.765795\n                                                         0.396852      -3.213702\n                                                         0.480094     -11.290279\n                                                         0.578055      33.333097\n                                                         0.703605     -27.251732\n                                                         0.862658       8.390367\n\n                          2       2    4px               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n                          2       2    4py               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n                          2       2    4pz               0.220648       1.497214\n                                                         0.396852     -29.770249\n                                                         0.480094      74.918767\n                                                         0.578055     -77.081909\n                                                         0.703605      38.422885\n                                                         0.862658      -8.196574\n\n                          2       3    5dx2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n                          2       3    5dxy              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dxz              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dy2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n                          2       3    5dyz              0.220648       1.070877\n                                                         0.396852     -29.002058\n                                                         0.480094      82.272792\n                                                         0.578055     -95.160048\n                                                         0.703605      52.871905\n                                                         0.862658     -12.596418\n                          2       3    5dz2              0.220648       0.618271\n                                                         0.396852     -16.744346\n                                                         0.480094      47.500218\n                                                         0.578055     -54.940679\n                                                         0.703605      30.525609\n                                                         0.862658      -7.272546\n\n     GTH Potential information for                                    GTH-PBE-q2\n\n       Description:                       Goedecker-Teter-Hutter pseudopotential\n                                           Goedecker et al., PRB 54, 1703 (1996)\n                                          Hartwigsen et al., PRB 58, 3641 (1998)\n                                                      Krack, TCA 114, 145 (2005)\n\n       Gaussian exponent of the core charge distribution:               1.502029\n       Electronic configuration (s p d ...):                                   2\n\n       Parameters of the local part of the GTH pseudopotential:\n\n                          rloc        C1          C2          C3          C4\n                        0.576960   -2.690407\n\n       Parameters of the non-local part of the GTH pseudopotential:\n\n                   l      r(l)      h(i,j,l)\n\n                   0    0.593924    3.503211   -0.716772\n                                   -0.716772    0.925348\n                   1    0.707157    0.831158\n\n\n MOLECULE KIND INFORMATION\n\n\n All atoms are their own molecule, skipping detailed information\n\n\n TOTAL NUMBERS AND MAXIMUM NUMBERS\n\n  Total number of            - Atomic kinds:                                   1\n                             - Atoms:                                          1\n                             - Shell sets:                                     2\n                             - Shells:                                         6\n                             - Primitive Cartesian functions:                 11\n                             - Cartesian basis functions:                     20\n                             - Spherical basis functions:                     18\n\n  Maximum angular momentum of- Orbital basis functions:                        2\n                             - Local part of the GTH pseudopotential:          0\n                             - Non-local part of the GTH pseudopotential:      2\n\n\n MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n\n  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n\n       1     1 Mg  12    0.000000    0.000000    0.000000      2.00      24.3050\n\n\n\n\n SCF PARAMETERS         Density guess:                                    ATOMIC\n                        --------------------------------------------------------\n                        max_scf:                                               1\n                        max_scf_history:                                       0\n                        max_diis:                                              4\n                        --------------------------------------------------------\n                        eps_scf:                                        1.00E-07\n                        eps_scf_history:                                0.00E+00\n                        eps_diis:                                       1.00E-01\n                        eps_eigval:                                     1.00E-05\n                        --------------------------------------------------------\n                        level_shift [a.u.]:                                 0.00\n                        added MOs                                         4    0\n                        --------------------------------------------------------\n                        Mixing method:                            BROYDEN_MIXING\n                                                charge density mixing in g-space\n                        --------------------------------------------------------\n                        Smear method:                                FERMI_DIRAC\n                        Electronic temperature [K]:                        157.9\n                        Electronic temperature [a.u.]:                  5.00E-04\n                        Accuracy threshold:                             1.00E-10\n                        --------------------------------------------------------\n                        No outer SCF\n\n PW_GRID| Information for grid number                                          1\n PW_GRID| Cutoff [a.u.]                                                     18.4\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -6       5                Points:          12\n PW_GRID|   Bounds   2             -6       5                Points:          12\n PW_GRID|   Bounds   3             -6       5                Points:          12\n PW_GRID| Volume element (a.u.^3)  0.9057E-01     Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          2\n PW_GRID| Cutoff [a.u.]                                                      6.1\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -4       3                Points:           8\n PW_GRID|   Bounds   2             -4       3                Points:           8\n PW_GRID|   Bounds   3             -4       3                Points:           8\n PW_GRID| Volume element (a.u.^3)  0.3057         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          3\n PW_GRID| Cutoff [a.u.]                                                      2.0\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -2       1                Points:           4\n PW_GRID|   Bounds   2             -2       1                Points:           4\n PW_GRID|   Bounds   3             -2       1                Points:           4\n PW_GRID| Volume element (a.u.^3)   2.446         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n PW_GRID| Information for grid number                                          4\n PW_GRID| Cutoff [a.u.]                                                      0.7\n PW_GRID| spherical cutoff:                                                   NO\n PW_GRID|   Bounds   1             -2       1                Points:           4\n PW_GRID|   Bounds   2             -2       1                Points:           4\n PW_GRID|   Bounds   3             -2       1                Points:           4\n PW_GRID| Volume element (a.u.^3)   2.446         Volume (a.u.^3)       156.5131\n PW_GRID| Grid span                                                    FULLSPACE\n\n POISSON| Solver                                                        PERIODIC\n POISSON| Periodicity                                                        XYZ\n\n RS_GRID| Information for grid number                                          1\n RS_GRID|   Bounds   1             -6       5                Points:          12\n RS_GRID|   Bounds   2             -6       5                Points:          12\n RS_GRID|   Bounds   3             -6       5                Points:          12\n\n RS_GRID| Information for grid number                                          2\n RS_GRID|   Bounds   1             -4       3                Points:           8\n RS_GRID|   Bounds   2             -4       3                Points:           8\n RS_GRID|   Bounds   3             -4       3                Points:           8\n\n RS_GRID| Information for grid number                                          3\n RS_GRID|   Bounds   1             -2       1                Points:           4\n RS_GRID|   Bounds   2             -2       1                Points:           4\n RS_GRID|   Bounds   3             -2       1                Points:           4\n\n RS_GRID| Information for grid number                                          4\n RS_GRID|   Bounds   1             -2       1                Points:           4\n RS_GRID|   Bounds   2             -2       1                Points:           4\n RS_GRID|   Bounds   3             -2       1                Points:           4\n\n Number of electrons:                                                          2\n Number of occupied orbitals:                                                  1\n Number of molecular orbitals:                                                 5\n\n Number of orbital functions:                                                 18\n Number of independent orbital functions:                                     18\n\n Extrapolation method: initial_guess\n\n Atomic guess: The first density matrix is obtained in terms of atomic orbitals\n               and electronic configurations assigned to each atomic kind\n\n Guess for atomic kind: Mg\n\n Electronic structure\n    Total number of core electrons                                         10.00\n    Total number of valence electrons                                       2.00\n    Total number of electrons                                              12.00\n    Multiplicity                                                   not specified\n    S   [  2.00  2.00] 2.00\n    P   [  6.00]\n\n\n *******************************************************************************\n                  Iteration          Convergence                     Energy [au]\n *******************************************************************************\n                          1        0.450866E-02                  -0.735938424170\n                          2        0.173646E-03                  -0.735943231943\n                          3        0.224554E-07                  -0.735943239082\n\n Energy components [Hartree]           Total Energy ::           -0.735943239082\n                                        Band Energy ::           -0.170223083884\n                                     Kinetic Energy ::            0.383696410659\n                                   Potential Energy ::           -1.119639649741\n                                      Virial (-V/T) ::            2.918035245154\n                                        Core Energy ::           -1.043063546624\n                                          XC Energy ::           -0.369655695440\n                                     Coulomb Energy ::            0.676776002981\n                       Total Pseudopotential Energy ::           -1.463503440618\n                       Local Pseudopotential Energy ::           -1.742048485333\n                    Nonlocal Pseudopotential Energy ::            0.278545044716\n                                        Confinement ::            0.367434833351\n\n Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]\n\n                       1     0          2.000      -0.085112           -2.316003\n\n\n Total Electron Density at R=0:                                         0.000119\n Re-scaling the density matrix to get the right number of electrons\n                  # Electrons              Trace(P)               Scaling factor\n                            2                 5.492                        0.364\n\n\n SCF WAVEFUNCTION OPTIMIZATION\n\n  Step     Update method      Time    Convergence         Total energy    Change\n  ------------------------------------------------------------------------------\n\n\n MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER SCF STEP 0\n\n# MO index          MO eigenvalue [a.u.]            MO occupation\n         1                     -0.147786                 2.000000\n         2                      0.568272                 0.000000\n         3                      0.568274                 0.000000\n         4                      0.568274                 0.000000\n         5                      0.624984                 0.000000\n# Sum                                                    2.000000\n\n  Fermi energy:                 0.549354\n\n  HOMO-LUMO gap:                0.716058 =  19.48 eV\n\n     1 NoMix/Diag. 0.40E+00    0.4     0.02871189        -1.1844755581 -1.18E+00\n\n  Leaving inner SCF loop after reaching     1 steps.\n\n\n  Electronic density on regular grids:         -2.0000000000       -0.0000000000\n  Core density on regular grids:                1.9999999986       -0.0000000014\n  Total charge density on r-space grids:       -0.0000000015\n  Total charge density g-space grids:          -0.0000000015\n\n  Overlap energy of the core charge distribution:               0.00000000000049\n  Self energy of the core charge distribution:                 -1.95573147986197\n  Core Hamiltonian energy:                                      0.19076838188720\n  Hartree energy:                                               1.00742354502905\n  Exchange-correlation energy:                                 -0.42693600510872\n  Electronic entropic energy:                                  -0.00000000000000\n  Fermi energy:                                                 0.54935437851948\n\n  Total energy:                                                -1.18447555805395\n\n *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***\n\n\n\n MO EIGENVALUES AND MO OCCUPATION NUMBERS\n\n# MO index          MO eigenvalue [a.u.]            MO occupation\n         1                     -0.147786                 2.000000\n         2                      0.568272                 0.000000\n         3                      0.568274                 0.000000\n         4                      0.568274                 0.000000\n         5                      0.624984                 0.000000\n# Sum                                                    2.000000\n\n  Fermi energy:                 0.549354\n\n  HOMO-LUMO gap:                0.716058 =  19.48 eV\n\n\n   Calculate PDOS at iteration step                                  0\n  Reached convergence in            1  iterations \n\n   Compute           13    additional unoccupied KS orbitals\n\n              ---- PDOS: start iteration on the KS states --- \n\n !-----------------------------------------------------------------------------!\n                     Mulliken Population Analysis\n\n #  Atom  Element  Kind  Atomic population                           Net charge\n       1     Mg       1          2.000000                              0.000000\n # Total charge                              2.000000                  0.000000\n\n !-----------------------------------------------------------------------------!\n\n !-----------------------------------------------------------------------------!\n                           Hirshfeld Charges\n\n  #Atom  Element  Kind  Ref Charge     Population                    Net charge\n      1       Mg     1       2.000          2.000                        -0.000\n\n  Total Charge                                                           -0.000\n !-----------------------------------------------------------------------------!\n\n ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.184475558053947\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                DBCSR STATISTICS                             -\n -                                                                             -\n -------------------------------------------------------------------------------\n COUNTER                                    TOTAL       BLAS       SMM       ACC\n flops     5 x    13 x     5                  650     100.0%      0.0%      0.0%\n flops     5 x     5 x    18                  900     100.0%      0.0%      0.0%\n flops     5 x    13 x    18                 2340     100.0%      0.0%      0.0%\n flops    18 x    13 x     5                 2340     100.0%      0.0%      0.0%\n flops    18 x    18 x     5                 3240     100.0%      0.0%      0.0%\n flops    18 x     5 x    18                 3240     100.0%      0.0%      0.0%\n flops    13 x    13 x    18                 6084     100.0%      0.0%      0.0%\n flops    18 x    18 x    18                11664     100.0%      0.0%      0.0%\n flops    18 x    13 x    13                24336     100.0%      0.0%      0.0%\n flops    18 x    13 x    18                58968     100.0%      0.0%      0.0%\n flops inhomo. stacks                           0       0.0%      0.0%      0.0%\n flops total                       113.762000E+03     100.0%      0.0%      0.0%\n flops max/rank                    113.762000E+03     100.0%      0.0%      0.0%\n matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%\n matmuls total                                 19     100.0%      0.0%      0.0%\n number of processed stacks                    19     100.0%      0.0%      0.0%\n average stack size                                     1.0       0.0       0.0\n marketing flops                   123.338000E+03\n -------------------------------------------------------------------------------\n # multiplications                             21\n max memory usage/rank              47.071232E+06\n # max total images/rank                        1\n # max 3D layers                                1\n # MPI messages exchanged                       0\n MPI messages size (bytes):\n  total size                         0.000000E+00\n  min size                           0.000000E+00\n  max size                           0.000000E+00\n  average size                       0.000000E+00\n MPI breakdown and total messages size (bytes):\n             size <=      128                   0                        0\n       128 < size <=     8192                   0                        0\n      8192 < size <=    32768                   0                        0\n     32768 < size <=   131072                   0                        0\n    131072 < size <=  4194304                   0                        0\n   4194304 < size <= 16777216                   0                        0\n  16777216 < size                               0                        0\n -------------------------------------------------------------------------------\n\n MEMORY| Estimated peak process memory [MiB]                                  45\n\n -------------------------------------------------------------------------------\n ----                             MULTIGRID INFO                            ----\n -------------------------------------------------------------------------------\n count for grid        1:          11459          cutoff [a.u.]           18.37\n count for grid        2:             48          cutoff [a.u.]            6.12\n count for grid        3:              0          cutoff [a.u.]            2.04\n count for grid        4:              0          cutoff [a.u.]            0.68\n total gridlevel count  :          11507\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                           R E F E R E N C E S                               -\n -                                                                             -\n -------------------------------------------------------------------------------\n \n CP2K version 6.1, the CP2K developers group (2018).\n CP2K is freely available from https://www.cp2k.org/ .\n\n Schuett, Ole; Messmer, Peter; Hutter, Juerg; VandeVondele, Joost. \n Electronic Structure Calculations on Graphics Processing Units, John Wiley & Sons, Ltd, 173-190 (2016). \n GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory.\n http://dx.doi.org/10.1002/9781118670712.ch8\n\n\n Borstnik, U; VandeVondele, J; Weber, V; Hutter, J. \n PARALLEL COMPUTING, 40 (5-6), 47-58 (2014). \n Sparse matrix multiplication: The distributed block-compressed sparse\n row library.\n http://dx.doi.org/10.1016/j.parco.2014.03.012\n\n\n Hutter, J; Iannuzzi, M; Schiffmann, F; VandeVondele, J. \n WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE, 4 (1), 15-25 (2014). \n CP2K: atomistic simulations of condensed matter systems.\n http://dx.doi.org/10.1002/wcms.1159\n\n\n Krack, M. \n THEORETICAL CHEMISTRY ACCOUNTS, 114 (1-3), 145-152 (2005). \n Pseudopotentials for H to Kr optimized for gradient-corrected\n exchange-correlation functionals.\n http://dx.doi.org/10.1007/s00214-005-0655-y\n\n\n VandeVondele, J; Krack, M; Mohamed, F; Parrinello, M; Chassaing, T;\n Hutter, J. COMPUTER PHYSICS COMMUNICATIONS, 167 (2), 103-128 (2005). \n QUICKSTEP: Fast and accurate density functional calculations using a\n mixed Gaussian and plane waves approach.\n http://dx.doi.org/10.1016/j.cpc.2004.12.014\n\n\n Frigo, M; Johnson, SG. \n PROCEEDINGS OF THE IEEE, 93 (2), 216-231 (2005). \n The design and implementation of FFTW3.\n http://dx.doi.org/10.1109/JPROC.2004.840301\n\n\n VandeVondele, J; Hutter, J. \n JOURNAL OF CHEMICAL PHYSICS, 118 (10), 4365-4369 (2003). \n An efficient orbital transformation method for electronic structure\n calculations.\n http://dx.doi.org/10.1063/1.1543154\n\n\n Hartwigsen, C; Goedecker, S; Hutter, J. \n PHYSICAL REVIEW B, 58 (7), 3641-3662 (1998). \n Relativistic separable dual-space Gaussian pseudopotentials from H to Rn.\n http://dx.doi.org/10.1103/PhysRevB.58.3641\n\n\n Lippert, G; Hutter, J; Parrinello, M. \n MOLECULAR PHYSICS, 92 (3), 477-487 (1997). \n A hybrid Gaussian and plane wave density functional scheme.\n http://dx.doi.org/10.1080/002689797170220\n\n\n Perdew, JP; Burke, K; Ernzerhof, M. \n PHYSICAL REVIEW LETTERS, 77 (18), 3865-3868 (1996). \n Generalized gradient approximation made simple.\n http://dx.doi.org/10.1103/PhysRevLett.77.3865\n\n\n Goedecker, S; Teter, M; Hutter, J. \n PHYSICAL REVIEW B, 54 (3), 1703-1710 (1996). \n Separable dual-space Gaussian pseudopotentials.\n http://dx.doi.org/10.1103/PhysRevB.54.1703\n\n\n -------------------------------------------------------------------------------\n -                                                                             -\n -                                T I M I N G                                  -\n -                                                                             -\n -------------------------------------------------------------------------------\n SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME\n                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM\n CP2K                                 1  1.0    0.003    0.003    1.342    1.342\n qs_energies                          1  2.0    0.000    0.000    1.299    1.299\n qs_rho_update_rho                    2  5.0    0.000    0.000    0.805    0.805\n calculate_rho_elec                   2  6.0    0.804    0.804    0.805    0.805\n scf_env_do_scf                       1  3.0    0.000    0.000    0.752    0.752\n scf_env_do_scf_inner_loop            1  4.0    0.000    0.000    0.752    0.752\n init_scf_run                         1  3.0    0.000    0.000    0.420    0.420\n scf_env_initial_rho_setup            1  4.0    0.000    0.000    0.420    0.420\n qs_ks_update_qs_env                  1  5.0    0.000    0.000    0.364    0.364\n rebuild_ks_matrix                    1  6.0    0.000    0.000    0.364    0.364\n qs_ks_build_kohn_sham_matrix         1  7.0    0.000    0.000    0.364    0.364\n sum_up_and_integrate                 1  8.0    0.000    0.000    0.363    0.363\n integrate_v_rspace                   1  9.0    0.363    0.363    0.363    0.363\n qs_energies_init_hamiltonians        1  3.0    0.000    0.000    0.114    0.114\n build_core_hamiltonian_matrix        1  4.0    0.000    0.000    0.096    0.096\n build_core_ppl                       1  5.0    0.042    0.042    0.042    0.042\n build_overlap_matrix                 1  5.0    0.027    0.027    0.027    0.027\n -------------------------------------------------------------------------------\n\n The number of warnings for this run is : 1\n \n -------------------------------------------------------------------------------\n  **** **** ******  **  PROGRAM ENDED AT                 2019-05-17 09:50:14.054\n ***** ** ***  *** **   PROGRAM RAN ON                                  mt-rf614\n **    ****   ******    PROGRAM RAN BY                                     rf614\n ***** **    ** ** **   PROGRAM PROCESS ID                                 33982\n  **** **  *******  **  PROGRAM STOPPED IN /media/ssd1/rf614/Work/Documents/jobs\n                                           /Corrosion_Work/random/learning_cp2k/\n                                           get_dos/gamma_point\n"

	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath


def _getCP2KFullFileStr_fccOpt():
	return """
 DBCSR| Multiplication driver                                               BLAS
 DBCSR| Multrec recursion limit                                              512
 DBCSR| Multiplication stack size                                           1000
 DBCSR| Maximum elements for images                                    UNLIMITED
 DBCSR| Multiplicative factor virtual images                                   1
 DBCSR| Multiplication size stacks                                             3


  **** **** ******  **  PROGRAM STARTED AT               2020-07-15 11:33:14.052
 ***** ** ***  *** **   PROGRAM STARTED ON                              mt-rf614
 **    ****   ******    PROGRAM STARTED BY                                 rf614
 ***** **    ** ** **   PROGRAM PROCESS ID                                140639
  **** **  *******  **  PROGRAM STARTED IN /media/ssd1/rf614/Work/Documents/jobs
                                           /Corrosion_Work/random/learning_cp2k/
                                           mg_opt/fcc_opt/constrained_angles/tem
                                           p

 CP2K| version string:                                          CP2K version 6.1
 CP2K| source code revision number:                                    svn:18464
 CP2K| cp2kflags: fftw3 max_contr=4                                             
 CP2K| is freely available from                            https://www.cp2k.org/
 CP2K| Program compiled at                          Tue 30 Jul 13:04:07 BST 2019
 CP2K| Program compiled on                                              mt-rf614
 CP2K| Program compiled for                                Linux-x86-64-gfortran
 CP2K| Data directory path             /media/ssd1/rf614/Work/CP2K/cp2k-6.1/data
 CP2K| Input file name                                          vol_152pt735.inp

 GLOBAL| Force Environment number                                              1
 GLOBAL| Basis set file name                                         PLATO_BASIS
 GLOBAL| Potential file name                                      GTH_POTENTIALS
 GLOBAL| MM Potential file name                                     MM_POTENTIAL
 GLOBAL| Coordinate file name                                      __STD_INPUT__
 GLOBAL| Method name                                                        CP2K
 GLOBAL| Project name                                               vol_152pt735
 GLOBAL| Preferred FFT library                                             FFTW3
 GLOBAL| Preferred diagonalization lib.                                       SL
 GLOBAL| Run type                                                       CELL_OPT
 GLOBAL| All-to-all communication in single precision                          F
 GLOBAL| FFTs using library dependent lengths                                  F
 GLOBAL| Global print level                                               MEDIUM
 GLOBAL| Total number of message passing processes                             1
 GLOBAL| Number of threads for this process                                    1
 GLOBAL| This output is from process                                           0
 GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-1620 v4 @ 3.50GHz

 MEMORY| system memory details [Kb]
 MEMORY|                        rank 0           min           max       average
 MEMORY| MemTotal             65866104      65866104      65866104      65866104
 MEMORY| MemFree              59290372      59290372      59290372      59290372
 MEMORY| Buffers                313416        313416        313416        313416
 MEMORY| Cached                3836760       3836760       3836760       3836760
 MEMORY| Slab                   322868        322868        322868        322868
 MEMORY| SReclaimable           239368        239368        239368        239368
 MEMORY| MemLikelyFree        63679916      63679916      63679916      63679916


 *** Fundamental physical constants (SI units) ***

 *** Literature: B. J. Mohr and B. N. Taylor,
 ***             CODATA recommended values of the fundamental physical
 ***             constants: 2006, Web Version 5.1
 ***             http://physics.nist.gov/constants

 Speed of light in vacuum [m/s]                             2.99792458000000E+08
 Magnetic constant or permeability of vacuum [N/A**2]       1.25663706143592E-06
 Electric constant or permittivity of vacuum [F/m]          8.85418781762039E-12
 Planck constant (h) [J*s]                                  6.62606896000000E-34
 Planck constant (h-bar) [J*s]                              1.05457162825177E-34
 Elementary charge [C]                                      1.60217648700000E-19
 Electron mass [kg]                                         9.10938215000000E-31
 Electron g factor [ ]                                     -2.00231930436220E+00
 Proton mass [kg]                                           1.67262163700000E-27
 Fine-structure constant                                    7.29735253760000E-03
 Rydberg constant [1/m]                                     1.09737315685270E+07
 Avogadro constant [1/mol]                                  6.02214179000000E+23
 Boltzmann constant [J/K]                                   1.38065040000000E-23
 Atomic mass unit [kg]                                      1.66053878200000E-27
 Bohr radius [m]                                            5.29177208590000E-11

 *** Conversion factors ***

 [u] -> [a.u.]                                              1.82288848426455E+03
 [Angstrom] -> [Bohr] = [a.u.]                              1.88972613288564E+00
 [a.u.] = [Bohr] -> [Angstrom]                              5.29177208590000E-01
 [a.u.] -> [s]                                              2.41888432650478E-17
 [a.u.] -> [fs]                                             2.41888432650478E-02
 [a.u.] -> [J]                                              4.35974393937059E-18
 [a.u.] -> [N]                                              8.23872205491840E-08
 [a.u.] -> [K]                                              3.15774647902944E+05
 [a.u.] -> [kJ/mol]                                         2.62549961709828E+03
 [a.u.] -> [kcal/mol]                                       6.27509468713739E+02
 [a.u.] -> [Pa]                                             2.94210107994716E+13
 [a.u.] -> [bar]                                            2.94210107994716E+08
 [a.u.] -> [atm]                                            2.90362800883016E+08
 [a.u.] -> [eV]                                             2.72113838565563E+01
 [a.u.] -> [Hz]                                             6.57968392072181E+15
 [a.u.] -> [1/cm] (wave numbers)                            2.19474631370540E+05
 [a.u./Bohr**2] -> [1/cm]                                   5.14048714338585E+03
 

 CELL_TOP| Volume [angstrom^3]:                                           22.633
 CELL_TOP| Vector a [angstrom     3.175     0.000     0.000    |a| =       3.175
 CELL_TOP| Vector b [angstrom     1.588     2.750     0.000    |b| =       3.175
 CELL_TOP| Vector c [angstrom     1.588     0.917     2.592    |c| =       3.175
 CELL_TOP| Angle (b,c), alpha [degree]:                                   60.000
 CELL_TOP| Angle (a,c), beta  [degree]:                                   60.000
 CELL_TOP| Angle (a,b), gamma [degree]:                                   60.000
 CELL_TOP| Numerically orthorhombic:                                          NO

 GENERATE|  Preliminary Number of Bonds generated:                             0
 GENERATE|  Achieved consistency in connectivity generation.

 CELL| Volume [angstrom^3]:                                               22.633
 CELL| Vector a [angstrom]:       3.175     0.000     0.000    |a| =       3.175
 CELL| Vector b [angstrom]:       1.588     2.750     0.000    |b| =       3.175
 CELL| Vector c [angstrom]:       1.588     0.917     2.592    |c| =       3.175
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 CELL_REF| Volume [angstrom^3]:                                           22.633
 CELL_REF| Vector a [angstrom     3.175     0.000     0.000    |a| =       3.175
 CELL_REF| Vector b [angstrom     1.588     2.750     0.000    |b| =       3.175
 CELL_REF| Vector c [angstrom     1.588     0.917     2.592    |c| =       3.175
 CELL_REF| Angle (b,c), alpha [degree]:                                   60.000
 CELL_REF| Angle (a,c), beta  [degree]:                                   60.000
 CELL_REF| Angle (a,b), gamma [degree]:                                   60.000
 CELL_REF| Numerically orthorhombic:                                          NO

 *** WARNING in cryssym.F:163 :: Symmetry library SPGLIB not available ***


 *******************************************************************************
                                    Kpoints
 *******************************************************************************
 BRILLOUIN| K-point scheme                                        Monkhorst-Pack
 BRILLOUIN| K-Point grid                                            1   1   1
 BRILLOUIN| Accuracy in Symmetry determination                      0.100000E-05
 BRILLOUIN| K-Point point group symmetrization                               OFF
 BRILLOUIN| Wavefunction type                                            COMPLEX
 BRILLOUIN| List of Kpoints [2 Pi/Bohr]                                      500
 BRILLOUIN| Number           Weight            X              Y              Z
 BRILLOUIN|     1           0.00200       -0.45000       -0.45000       -0.45000
 *******************************************************************************

 *******************************************************************************
 *******************************************************************************
 **                                                                           **
 **     #####                         ##              ##                      **
 **    ##   ##            ##          ##              ##                      **
 **   ##     ##                       ##            ######                    **
 **   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **
 **   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **
 **   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **
 **    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **
 **     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **
 **           ##                                                    ##        **
 **                                                                           **
 **                                                ... make the atoms dance   **
 **                                                                           **
 **            Copyright (C) by CP2K developers group (2000 - 2018)           **
 **                                                                           **
 *******************************************************************************

 DFT| Spin restricted Kohn-Sham (RKS) calculation                            RKS
 DFT| Multiplicity                                                             1
 DFT| Number of spin states                                                    1
 DFT| Charge                                                                   0
 DFT| Self-interaction correction (SIC)                                       NO
 DFT| Cutoffs: density                                              1.000000E-10
 DFT|          gradient                                             1.000000E-10
 DFT|          tau                                                  1.000000E-10
 DFT|          cutoff_smoothing_range                               0.000000E+00
 DFT| XC density smoothing                                                  NONE
 DFT| XC derivatives                                                          PW
 FUNCTIONAL| ROUTINE=NEW
 FUNCTIONAL| PBE:
 FUNCTIONAL| J.P.Perdew, K.Burke, M.Ernzerhof, Phys. Rev. Letter, vol. 77, n 18,
 FUNCTIONAL|  pp. 3865-3868, (1996){spin unpolarized}                           

 QS| Method:                                                                 GPW
 QS| Density plane wave grid type                        NON-SPHERICAL FULLSPACE
 QS| Number of grid levels:                                                    4
 QS| Density cutoff [a.u.]:                                                 29.4
 QS| Multi grid cutoff [a.u.]: 1) grid level                                29.4
 QS|                           2) grid level                                 9.8
 QS|                           3) grid level                                 3.3
 QS|                           4) grid level                                 1.1
 QS| Grid level progression factor:                                          3.0
 QS| Relative density cutoff [a.u.]:                                        18.4
 QS| Consistent realspace mapping and integration 
 QS| Interaction thresholds: eps_pgf_orb:                                1.0E-05
 QS|                         eps_filter_matrix:                          0.0E+00
 QS|                         eps_core_charge:                            1.0E-12
 QS|                         eps_rho_gspace:                             1.0E-10
 QS|                         eps_rho_rspace:                             1.0E-10
 QS|                         eps_gvg_rspace:                             1.0E-05
 QS|                         eps_ppl:                                    1.0E-02
 QS|                         eps_ppnl:                                   1.0E-07


 ATOMIC KIND INFORMATION

  1. Atomic kind: Mg                                    Number of atoms:       1

     Orbital Basis Set                  spd-1z-s-2z-rc7pt5-r06pt0-plus-diffuse-s

       Number of orbital shell sets:                                           3
       Number of orbital shells:                                               5
       Number of primitive Cartesian functions:                               12
       Number of Cartesian basis functions:                                   12
       Number of spherical basis functions:                                   11
       Norm type:                                                              2

       Normalised Cartesian orbitals:

                        Set   Shell   Orbital            Exponent    Coefficient

                          1       1    3s                0.164100       1.180053
                                                         0.198454      -1.654259
                                                         0.264453       0.732902
                                                         0.510115      -0.228584
                                                         1.615631      -0.028470

                          1       2    4px               0.164100       0.538781
                                                         0.198454      -0.698692
                                                         0.264453       0.319622
                                                         0.510115      -0.076876
                                                         1.615631       0.017972
                          1       2    4py               0.164100       0.538781
                                                         0.198454      -0.698692
                                                         0.264453       0.319622
                                                         0.510115      -0.076876
                                                         1.615631       0.017972
                          1       2    4pz               0.164100       0.538781
                                                         0.198454      -0.698692
                                                         0.264453       0.319622
                                                         0.510115      -0.076876
                                                         1.615631       0.017972

                          1       3    5dx2              0.164100       0.165032
                                                         0.198454      -0.173803
                                                         0.264453       0.061115
                                                         0.510115       0.031075
                                                         1.615631       0.037517
                          1       3    5dxy              0.164100       0.285843
                                                         0.198454      -0.301036
                                                         0.264453       0.105855
                                                         0.510115       0.053823
                                                         1.615631       0.064981
                          1       3    5dxz              0.164100       0.285843
                                                         0.198454      -0.301036
                                                         0.264453       0.105855
                                                         0.510115       0.053823
                                                         1.615631       0.064981
                          1       3    5dy2              0.164100       0.165032
                                                         0.198454      -0.173803
                                                         0.264453       0.061115
                                                         0.510115       0.031075
                                                         1.615631       0.037517
                          1       3    5dyz              0.164100       0.285843
                                                         0.198454      -0.301036
                                                         0.264453       0.105855
                                                         0.510115       0.053823
                                                         1.615631       0.064981
                          1       3    5dz2              0.164100       0.165032
                                                         0.198454      -0.173803
                                                         0.264453       0.061115
                                                         0.510115       0.031075
                                                         1.615631       0.037517

                          2       1    3s                0.200374       0.741680
                                                         0.385552      -5.261029
                                                         0.464015      -4.514705
                                                         0.564438      25.414730
                                                         0.689075     -23.536951
                                                         0.843766       7.869628

                          3       1    3s                0.051676       0.077246

     GTH Potential information for                                    GTH-PBE-q2

       Description:                       Goedecker-Teter-Hutter pseudopotential
                                           Goedecker et al., PRB 54, 1703 (1996)
                                          Hartwigsen et al., PRB 58, 3641 (1998)
                                                      Krack, TCA 114, 145 (2005)

       Gaussian exponent of the core charge distribution:               1.502029
       Electronic configuration (s p d ...):                                   2

       Parameters of the local part of the GTH pseudopotential:

                          rloc        C1          C2          C3          C4
                        0.576960   -2.690407

       Parameters of the non-local part of the GTH pseudopotential:

                   l      r(l)      h(i,j,l)

                   0    0.593924    3.503211   -0.716772
                                   -0.716772    0.925348
                   1    0.707157    0.831158


 MOLECULE KIND INFORMATION


 All atoms are their own molecule, skipping detailed information


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   1
                             - Atoms:                                          1
                             - Shell sets:                                     3
                             - Shells:                                         5
                             - Primitive Cartesian functions:                 12
                             - Cartesian basis functions:                     12
                             - Spherical basis functions:                     11

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          0
                             - Non-local part of the GTH pseudopotential:      2


 MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom

  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass

       1     1 Mg  12    0.000000    0.000000    0.000000      2.00      24.3050




 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                             300
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-07
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        added MOs                                         4    0
                        --------------------------------------------------------
                        Mixing method:                            BROYDEN_MIXING
                                                charge density mixing in g-space
                        --------------------------------------------------------
                        Smear method:                                FERMI_DIRAC
                        Electronic temperature [K]:                        157.9
                        Electronic temperature [a.u.]:                  5.00E-04
                        Accuracy threshold:                             1.00E-10
                        --------------------------------------------------------
                        No outer SCF

 PW_GRID| Information for grid number                                          1
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4525E-01     Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          2
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2095         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          3
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7071         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          4
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.386         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                          1
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                          2
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                          3
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                          4
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 CELL_OPT| Pressure tolerance [bar]:                                       100.0
 CELL_OPT| Keep angles between the cell vectors:                             YES
 CELL_OPT| Keep cell symmetry:                                                NO
 CELL_OPT| Constraint:                                                      NONE

 BFGS| Use rational function optimization for step estimation:                NO
 BFGS| Use model Hessian for initial guess:                                   NO
 BFGS| Restart Hessian:                                                       NO
 BFGS| Trust radius:                                                       0.472

 *******************************************************************************
 ***                     STARTING   CELL   OPTIMIZATION                      ***
 ***                                   BFGS                                  ***
 *******************************************************************************

 CELL| Volume [angstrom^3]:                                               22.633
 CELL| Vector a [angstrom]:       3.175     0.000     0.000    |a| =       3.175
 CELL| Vector b [angstrom]:       1.588     2.750     0.000    |b| =       3.175
 CELL| Vector c [angstrom]:       1.588     0.917     2.592    |c| =       3.175
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:           1460          cutoff [a.u.]           29.40
 count for grid        2:            544          cutoff [a.u.]            9.80
 count for grid        3:             16          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :           2024

 PW_GRID| Information for grid number                                          5
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4525E-01     Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          6
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2095         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          7
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7071         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                          8
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.386         Volume (a.u.^3)       152.7350
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                          5
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                          6
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                          7
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                          8
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 NoMix/Diag. 0.40E+00    1.5     0.49648194        -0.9448081364 -9.45E-01
     2 Broy./Diag. 0.40E+00    2.0     0.00890232        -0.9277480659  1.71E-02
     3 Broy./Diag. 0.40E+00    2.0     0.02520229        -0.8802131651  4.75E-02
     4 Broy./Diag. 0.40E+00    2.0     0.00014835        -0.8802406772 -2.75E-05
     5 Broy./Diag. 0.40E+00    2.0     0.00010866        -0.8797684327  4.72E-04
     6 Broy./Diag. 0.40E+00    2.0     0.00000114        -0.8797680580  3.75E-07
     7 Broy./Diag. 0.40E+00    2.0     0.00000747        -0.8797677555  3.02E-07
     8 Broy./Diag. 0.40E+00    2.0     0.00000004        -0.8797678798 -1.24E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -2.0000000007       -0.0000000007
  Core density on regular grids:                1.9999999985       -0.0000000015
  Total charge density on r-space grids:       -0.0000000021
  Total charge density g-space grids:          -0.0000000021

  Overlap energy of the core charge distribution:               0.00000000000077
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.48619741993392
  Hartree energy:                                               1.02122004110676
  Exchange-correlation energy:                                 -0.43144128983853
  Electronic entropic energy:                                  -0.00001257113930
  Fermi energy:                                                 0.11183551444274

  Total energy:                                                -0.87976787980010

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                             -0.000000
 # Total charge                              2.000000                 -0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879767892541844


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00000000    -0.00000000     0.00000000
 SUM OF ATOMIC FORCES           0.00000000    -0.00000000     0.00000000     0.00000000

 --------  Informations at step =     0 ------------
  Optimization Method        =                 BFGS
  Total Energy               =        -0.8797678925
  Internal Pressure [bar]    =     14208.9926713290
  Used time                  =               17.995
 ---------------------------------------------------

 --------------------------
 OPTIMIZATION STEP:      1
 --------------------------

 CELL| Volume [angstrom^3]:                                               22.647
 CELL| Vector a [angstrom]:       3.176     0.000     0.000    |a| =       3.176
 CELL| Vector b [angstrom]:       1.588     2.750     0.000    |b| =       3.176
 CELL| Vector c [angstrom]:       1.588     0.917     2.593    |c| =       3.176
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:           7149          cutoff [a.u.]           29.40
 count for grid        2:           5242          cutoff [a.u.]            9.80
 count for grid        3:            547          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :          12942

 PW_GRID| Information for grid number                                          9
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4528E-01     Volume (a.u.^3)       152.8283
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         10
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2096         Volume (a.u.^3)       152.8283
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         11
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7075         Volume (a.u.^3)       152.8283
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         12
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.388         Volume (a.u.^3)       152.8283
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                          9
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                         10
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                         11
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                         12
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.5     0.49639799        -0.9447329952 -9.45E-01
     2 Broy./Diag. 0.40E+00    2.0     0.00890410        -0.9277691308  1.70E-02
     3 Broy./Diag. 0.40E+00    2.0     0.02519389        -0.8802178088  4.76E-02
     4 Broy./Diag. 0.40E+00    2.0     0.00014784        -0.8802456196 -2.78E-05
     5 Broy./Diag. 0.40E+00    2.0     0.00010853        -0.8797729025  4.73E-04
     6 Broy./Diag. 0.40E+00    2.0     0.00000115        -0.8797725174  3.85E-07
     7 Broy./Diag. 0.40E+00    2.0     0.00000747        -0.8797722197  2.98E-07
     8 Broy./Diag. 0.40E+00    2.0     0.00000004        -0.8797723443 -1.25E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -2.0000000007       -0.0000000007
  Core density on regular grids:                1.9999999986       -0.0000000014
  Total charge density on r-space grids:       -0.0000000021
  Total charge density g-space grids:          -0.0000000021

  Overlap energy of the core charge distribution:               0.00000000000076
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.48600412320470
  Hartree energy:                                               1.02132630340221
  Exchange-correlation energy:                                 -0.43135868731541
  Electronic entropic energy:                                  -0.00001260372588
  Fermi energy:                                                 0.11171438776458

  Total energy:                                                -0.87977234429731

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                              0.000000
 # Total charge                              2.000000                  0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879772357321208


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg         -0.00000000     0.00000000    -0.00000000
 SUM OF ATOMIC FORCES          -0.00000000     0.00000000    -0.00000000     0.00000000

 --------  Informations at step =     1 ------------
  Optimization Method        =                 BFGS
  Total Energy               =        -0.8797723573
  Internal Pressure [bar]    =     13966.2825496846
  Real energy change         =        -0.0000044648
  Predicted change in energy =        -0.0000022353
  Scaling factor             =         0.0000000000
  Step size                  =         0.0012207431
  Trust radius               =         0.4724315332
  Decrease in energy         =                  YES
  Used time                  =               17.795

  Convergence check :
  Max. step size             =         0.0012207431
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                  YES
  RMS step size              =         0.0007047989
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                  YES
  Max. gradient              =         0.0012002314
  Conv. limit for gradients  =         0.0004500000
  Conv. for gradients        =                   NO
  RMS gradient               =         0.0006929565
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. for gradients        =                   NO
  Pressure Deviation [bar]   =     13866.2825496846
  Pressure Tolerance [bar]   =       100.0000000000
  Conv. for  PRESSURE        =                   NO
 ---------------------------------------------------

 --------------------------
 OPTIMIZATION STEP:      2
 --------------------------

 CELL| Volume [angstrom^3]:                                               23.465
 CELL| Vector a [angstrom]:       3.214     0.000     0.000    |a| =       3.214
 CELL| Vector b [angstrom]:       1.607     2.783     0.000    |b| =       3.214
 CELL| Vector c [angstrom]:       1.607     0.928     2.624    |c| =       3.214
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:           7149          cutoff [a.u.]           29.40
 count for grid        2:           5242          cutoff [a.u.]            9.80
 count for grid        3:            547          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :          12942

 PW_GRID| Information for grid number                                         13
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4692E-01     Volume (a.u.^3)       158.3508
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         14
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2172         Volume (a.u.^3)       158.3508
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         15
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7331         Volume (a.u.^3)       158.3508
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         16
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.474         Volume (a.u.^3)       158.3508
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                         13
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                         14
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                         15
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                         16
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.4     0.49129089        -0.9396218008 -9.40E-01
     2 Broy./Diag. 0.40E+00    1.9     0.00902073        -0.9288669082  1.08E-02
     3 Broy./Diag. 0.40E+00    1.9     0.02473319        -0.8803626505  4.85E-02
     4 Broy./Diag. 0.40E+00    1.9     0.00011783        -0.8804090373 -4.64E-05
     5 Broy./Diag. 0.40E+00    1.9     0.00009983        -0.8799079908  5.01E-04
     6 Broy./Diag. 0.40E+00    1.9     0.00000178        -0.8799070196  9.71E-07
     7 Broy./Diag. 0.40E+00    1.9     0.00000762        -0.8799070148  4.80E-09
     8 Broy./Diag. 0.40E+00    1.9     0.00000003        -0.8799071584 -1.44E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -1.9999999998        0.0000000002
  Core density on regular grids:                1.9999999984       -0.0000000016
  Total charge density on r-space grids:       -0.0000000014
  Total charge density g-space grids:          -0.0000000014

  Overlap energy of the core charge distribution:               0.00000000000039
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.47508636724696
  Hartree energy:                                               1.02734079612643
  Exchange-correlation energy:                                 -0.42658784297437
  Electronic entropic energy:                                  -0.00001499891312
  Fermi energy:                                                 0.10472292938713

  Total energy:                                                -0.87990715837542

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                              0.000000
 # Total charge                              2.000000                  0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879907186254556


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg         -0.00000000    -0.00000000     0.00000000
 SUM OF ATOMIC FORCES          -0.00000000    -0.00000000     0.00000000     0.00000000

 --------  Informations at step =     2 ------------
  Optimization Method        =                 BFGS
  Total Energy               =        -0.8799071863
  Internal Pressure [bar]    =       748.7289634804
  Real energy change         =        -0.0001348289
  Predicted change in energy =        -0.0001286024
  Scaling factor             =         0.0000000000
  Step size                  =         0.0714314811
  Trust radius               =         0.4724315332
  Decrease in energy         =                  YES
  Used time                  =               16.757

  Convergence check :
  Max. step size             =         0.0714314811
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                   NO
  RMS step size              =         0.0412411364
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                   NO
  Max. gradient              =         0.0000574856
  Conv. limit for gradients  =         0.0004500000
  Conv. in gradients         =                  YES
  RMS gradient               =         0.0000331961
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. in RMS gradients     =                  YES
  Pressure Deviation [bar]   =       648.7289634804
  Pressure Tolerance [bar]   =       100.0000000000
  Conv. for  PRESSURE        =                   NO
 ---------------------------------------------------

 --------------------------
 OPTIMIZATION STEP:      3
 --------------------------

 CELL| Volume [angstrom^3]:                                               23.507
 CELL| Vector a [angstrom]:       3.215     0.000     0.000    |a| =       3.215
 CELL| Vector b [angstrom]:       1.608     2.785     0.000    |b| =       3.215
 CELL| Vector c [angstrom]:       1.608     0.928     2.625    |c| =       3.215
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:           6897          cutoff [a.u.]           29.40
 count for grid        2:           5050          cutoff [a.u.]            9.80
 count for grid        3:            475          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :          12426

 PW_GRID| Information for grid number                                         17
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4700E-01     Volume (a.u.^3)       158.6321
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         18
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2176         Volume (a.u.^3)       158.6321
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         19
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7344         Volume (a.u.^3)       158.6321
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         20
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.479         Volume (a.u.^3)       158.6321
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                         17
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                         18
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                         19
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                         20
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.4     0.49105468        -0.9393682450 -9.39E-01
     2 Broy./Diag. 0.40E+00    1.9     0.00902669        -0.9289154120  1.05E-02
     3 Broy./Diag. 0.40E+00    1.9     0.02470980        -0.8803635802  4.86E-02
     4 Broy./Diag. 0.40E+00    1.9     0.00011625        -0.8804109295 -4.73E-05
     5 Broy./Diag. 0.40E+00    1.9     0.00009936        -0.8799084236  5.03E-04
     6 Broy./Diag. 0.40E+00    1.9     0.00000181        -0.8799074235  1.00E-06
     7 Broy./Diag. 0.40E+00    1.9     0.00000762        -0.8799074343 -1.08E-08
     8 Broy./Diag. 0.40E+00    1.9     0.00000003        -0.8799075788 -1.45E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -1.9999999999        0.0000000001
  Core density on regular grids:                1.9999999985       -0.0000000015
  Total charge density on r-space grids:       -0.0000000014
  Total charge density g-space grids:          -0.0000000014

  Overlap energy of the core charge distribution:               0.00000000000038
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.47455700181866
  Hartree energy:                                               1.02763308495014
  Exchange-correlation energy:                                 -0.42635105074495
  Electronic entropic energy:                                  -0.00001513499395
  Fermi energy:                                                 0.10437639272444

  Total energy:                                                -0.87990757883134

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                             -0.000000
 # Total charge                              2.000000                 -0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879907607331958


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg         -0.00000000    -0.00000000    -0.00000000
 SUM OF ATOMIC FORCES          -0.00000000    -0.00000000    -0.00000000     0.00000000

 --------  Informations at step =     3 ------------
  Optimization Method        =                 BFGS
  Total Energy               =        -0.8799076073
  Internal Pressure [bar]    =       132.7550848565
  Real energy change         =        -0.0000004211
  Predicted change in energy =        -0.0000003100
  Scaling factor             =         0.0000000000
  Step size                  =         0.0035940888
  Trust radius               =         0.4724315332
  Decrease in energy         =                  YES
  Used time                  =               16.735

  Convergence check :
  Max. step size             =         0.0035940888
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                   NO
  RMS step size              =         0.0020750628
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                   NO
  Max. gradient              =         0.0000028947
  Conv. limit for gradients  =         0.0004500000
  Conv. in gradients         =                  YES
  RMS gradient               =         0.0000016782
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. in RMS gradients     =                  YES
  Pressure Deviation [bar]   =        32.7550848565
  Pressure Tolerance [bar]   =       100.0000000000
  Conv. for  PRESSURE        =                  YES
 ---------------------------------------------------

 --------------------------
 OPTIMIZATION STEP:      4
 --------------------------

 CELL| Volume [angstrom^3]:                                               23.509
 CELL| Vector a [angstrom]:       3.216     0.000     0.000    |a| =       3.216
 CELL| Vector b [angstrom]:       1.608     2.785     0.000    |b| =       3.216
 CELL| Vector c [angstrom]:       1.608     0.928     2.625    |c| =       3.216
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:           6873          cutoff [a.u.]           29.40
 count for grid        2:           5050          cutoff [a.u.]            9.80
 count for grid        3:            475          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :          12402

 PW_GRID| Information for grid number                                         21
 PW_GRID| Cutoff [a.u.]                                                     29.4
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -7       7                Points:          15
 PW_GRID|   Bounds   2             -7       7                Points:          15
 PW_GRID|   Bounds   3             -7       7                Points:          15
 PW_GRID| Volume element (a.u.^3)  0.4701E-01     Volume (a.u.^3)       158.6471
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         22
 PW_GRID| Cutoff [a.u.]                                                      9.8
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -4       4                Points:           9
 PW_GRID|   Bounds   2             -4       4                Points:           9
 PW_GRID|   Bounds   3             -4       4                Points:           9
 PW_GRID| Volume element (a.u.^3)  0.2176         Volume (a.u.^3)       158.6471
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         23
 PW_GRID| Cutoff [a.u.]                                                      3.3
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -3       2                Points:           6
 PW_GRID|   Bounds   2             -3       2                Points:           6
 PW_GRID|   Bounds   3             -3       2                Points:           6
 PW_GRID| Volume element (a.u.^3)  0.7345         Volume (a.u.^3)       158.6471
 PW_GRID| Grid span                                                    FULLSPACE

 PW_GRID| Information for grid number                                         24
 PW_GRID| Cutoff [a.u.]                                                      1.1
 PW_GRID| spherical cutoff:                                                   NO
 PW_GRID|   Bounds   1             -2       1                Points:           4
 PW_GRID|   Bounds   2             -2       1                Points:           4
 PW_GRID|   Bounds   3             -2       1                Points:           4
 PW_GRID| Volume element (a.u.^3)   2.479         Volume (a.u.^3)       158.6471
 PW_GRID| Grid span                                                    FULLSPACE

 POISSON| Solver                                                        PERIODIC
 POISSON| Periodicity                                                        XYZ

 RS_GRID| Information for grid number                                         21
 RS_GRID|   Bounds   1             -7       7                Points:          15
 RS_GRID|   Bounds   2             -7       7                Points:          15
 RS_GRID|   Bounds   3             -7       7                Points:          15

 RS_GRID| Information for grid number                                         22
 RS_GRID|   Bounds   1             -4       4                Points:           9
 RS_GRID|   Bounds   2             -4       4                Points:           9
 RS_GRID|   Bounds   3             -4       4                Points:           9

 RS_GRID| Information for grid number                                         23
 RS_GRID|   Bounds   1             -3       2                Points:           6
 RS_GRID|   Bounds   2             -3       2                Points:           6
 RS_GRID|   Bounds   3             -3       2                Points:           6

 RS_GRID| Information for grid number                                         24
 RS_GRID|   Bounds   1             -2       1                Points:           4
 RS_GRID|   Bounds   2             -2       1                Points:           4
 RS_GRID|   Bounds   3             -2       1                Points:           4

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.4     0.49104212        -0.9393547576 -9.39E-01
     2 Broy./Diag. 0.40E+00    1.9     0.00902701        -0.9289179763  1.04E-02
     3 Broy./Diag. 0.40E+00    1.9     0.02470855        -0.8803636132  4.86E-02
     4 Broy./Diag. 0.40E+00    1.9     0.00011617        -0.8804110138 -4.74E-05
     5 Broy./Diag. 0.40E+00    1.9     0.00009934        -0.8799084302  5.03E-04
     6 Broy./Diag. 0.40E+00    1.9     0.00000181        -0.8799074285  1.00E-06
     7 Broy./Diag. 0.40E+00    1.9     0.00000763        -0.8799074401 -1.16E-08
     8 Broy./Diag. 0.40E+00    1.9     0.00000003        -0.8799075847 -1.45E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -1.9999999999        0.0000000001
  Core density on regular grids:                1.9999999985       -0.0000000015
  Total charge density on r-space grids:       -0.0000000014
  Total charge density g-space grids:          -0.0000000014

  Overlap energy of the core charge distribution:               0.00000000000038
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.47452887072397
  Hartree energy:                                               1.02764861868715
  Exchange-correlation energy:                                 -0.42633845203270
  Electronic entropic energy:                                  -0.00001514225165
  Fermi energy:                                                 0.10435795693020

  Total energy:                                                -0.87990758473446

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                             -0.000000
 # Total charge                              2.000000                 -0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879907613267063


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00000000     0.00000000    -0.00000000
 SUM OF ATOMIC FORCES           0.00000000     0.00000000    -0.00000000     0.00000000

 --------  Informations at step =     4 ------------
  Optimization Method        =                 BFGS
  Total Energy               =        -0.8799076133
  Internal Pressure [bar]    =       100.0845768870
  Real energy change         =        -0.0000000059
  Predicted change in energy =        -0.0000000008
  Scaling factor             =         0.0000000000
  Step size                  =         0.0001913462
  Trust radius               =         0.4724315332
  Decrease in energy         =                  YES
  Used time                  =               16.731

  Convergence check :
  Max. step size             =         0.0001913462
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                  YES
  RMS step size              =         0.0001104814
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                  YES
  Max. gradient              =         0.0000000399
  Conv. limit for gradients  =         0.0004500000
  Conv. in gradients         =                  YES
  RMS gradient               =         0.0000000179
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. in RMS gradients     =                  YES
  Pressure Deviation [bar]   =         0.0845768870
  Pressure Tolerance [bar]   =       100.0000000000
  Conv. for  PRESSURE        =                  YES
 ---------------------------------------------------

 *******************************************************************************
 ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
 *******************************************************************************

                    Reevaluating energy at the minimum

 CELL| Volume [angstrom^3]:                                               23.509
 CELL| Vector a [angstrom]:       3.216     0.000     0.000    |a| =       3.216
 CELL| Vector b [angstrom]:       1.608     2.785     0.000    |b| =       3.216
 CELL| Vector c [angstrom]:       1.608     0.928     2.625    |c| =       3.216
 CELL| Angle (b,c), alpha [degree]:                                       60.000
 CELL| Angle (a,c), beta  [degree]:                                       60.000
 CELL| Angle (a,b), gamma [degree]:                                       60.000
 CELL| Numerically orthorhombic:                                              NO

 Number of electrons:                                                          2
 Number of occupied orbitals:                                                  1
 Number of molecular orbitals:                                                 5

 Number of orbital functions:                                                 11
 Number of independent orbital functions:                                     11

 Extrapolation method: initial_guess

 Atomic guess: The first density matrix is obtained in terms of atomic orbitals
               and electronic configurations assigned to each atomic kind

 Guess for atomic kind: Mg

 Electronic structure
    Total number of core electrons                                         10.00
    Total number of valence electrons                                       2.00
    Total number of electrons                                              12.00
    Multiplicity                                                   not specified
    S   [  2.00  2.00] 2.00
    P   [  6.00]


 *******************************************************************************
                  Iteration          Convergence                     Energy [au]
 *******************************************************************************
                          1        0.431953E-01                  -0.762803637831
                          2        0.215465E-02                  -0.765747311497
                          3        0.339278E-04                  -0.765754574588
                          4        0.110408E-05                  -0.765754576316
                          5        0.545122E-06                  -0.765754576318

 Energy components [Hartree]           Total Energy ::           -0.765754576318
                                        Band Energy ::           -0.262283723627
                                     Kinetic Energy ::            0.281098603880
                                   Potential Energy ::           -1.046853180198
                                      Virial (-V/T) ::            3.724149340299
                                        Core Energy ::           -1.039168153629
                                          XC Energy ::           -0.328640767724
                                     Coulomb Energy ::            0.602054345036
                       Total Pseudopotential Energy ::           -1.369004284385
                       Local Pseudopotential Energy ::           -1.590718742515
                    Nonlocal Pseudopotential Energy ::            0.221714458130
                                        Confinement ::            0.487375268766

 Orbital energies  State     L     Occupation   Energy[a.u.]          Energy[eV]

                       1     0          2.000      -0.131142           -3.568552


 Total Electron Density at R=0:                                         0.000044
 Re-scaling the density matrix to get the right number of electrons
                  # Electrons              Trace(P)               Scaling factor
                            2                 2.000                        1.000


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.4     0.49104212        -0.9393547648 -9.39E-01
     2 Broy./Diag. 0.40E+00    1.9     0.00902701        -0.9289179763  1.04E-02
     3 Broy./Diag. 0.40E+00    1.9     0.02470855        -0.8803636132  4.86E-02
     4 Broy./Diag. 0.40E+00    1.9     0.00011617        -0.8804110138 -4.74E-05
     5 Broy./Diag. 0.40E+00    1.9     0.00009934        -0.8799084302  5.03E-04
     6 Broy./Diag. 0.40E+00    1.9     0.00000181        -0.8799074285  1.00E-06
     7 Broy./Diag. 0.40E+00    1.9     0.00000763        -0.8799074401 -1.16E-08
     8 Broy./Diag. 0.40E+00    1.9     0.00000003        -0.8799075847 -1.45E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -1.9999999999        0.0000000001
  Core density on regular grids:                1.9999999985       -0.0000000015
  Total charge density on r-space grids:       -0.0000000014
  Total charge density g-space grids:          -0.0000000014

  Overlap energy of the core charge distribution:               0.00000000000038
  Self energy of the core charge distribution:                 -1.95573147986197
  Core Hamiltonian energy:                                      0.47452887072397
  Hartree energy:                                               1.02764861868715
  Exchange-correlation energy:                                 -0.42633845203270
  Electronic entropic energy:                                  -0.00001514225165
  Fermi energy:                                                 0.10435795693020

  Total energy:                                                -0.87990758473446

 !-----------------------------------------------------------------------------!
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     Mg       1          2.000000                             -0.000000
 # Total charge                              2.000000                 -0.000000

 !-----------------------------------------------------------------------------!

 !-----------------------------------------------------------------------------!
                           Hirshfeld Charges

  #Atom  Element  Kind  Ref Charge     Population                    Net charge
      1       Mg     1       2.000          2.000                        -0.000

  Total Charge                                                           -0.000
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -0.879907584734461


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                DBCSR STATISTICS                             -
 -                                                                             -
 -------------------------------------------------------------------------------
 COUNTER                                    TOTAL       BLAS       SMM       ACC
 flops inhomo. stacks                           0       0.0%      0.0%      0.0%
 flops total                         0.000000E+00       0.0%      0.0%      0.0%
 flops max/rank                      0.000000E+00       0.0%      0.0%      0.0%
 matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%
 matmuls total                                  0       0.0%      0.0%      0.0%
 number of processed stacks                     0       0.0%      0.0%      0.0%
 average stack size                                     0.0       0.0       0.0
 marketing flops                     0.000000E+00
 -------------------------------------------------------------------------------

 MEMORY| Estimated peak process memory [MiB]                                  62

 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:          12286          cutoff [a.u.]           29.40
 count for grid        2:           9556          cutoff [a.u.]            9.80
 count for grid        3:            934          cutoff [a.u.]            3.27
 count for grid        4:              4          cutoff [a.u.]            1.09
 total gridlevel count  :          22780

 -------------------------------------------------------------------------------
 -                                                                             -
 -                           R E F E R E N C E S                               -
 -                                                                             -
 -------------------------------------------------------------------------------
 
 CP2K version 6.1, the CP2K developers group (2018).
 CP2K is freely available from https://www.cp2k.org/ .

 Schuett, Ole; Messmer, Peter; Hutter, Juerg; VandeVondele, Joost. 
 Electronic Structure Calculations on Graphics Processing Units, John Wiley & Sons, Ltd, 173-190 (2016). 
 GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory.
 http://dx.doi.org/10.1002/9781118670712.ch8


 Borstnik, U; VandeVondele, J; Weber, V; Hutter, J. 
 PARALLEL COMPUTING, 40 (5-6), 47-58 (2014). 
 Sparse matrix multiplication: The distributed block-compressed sparse
 row library.
 http://dx.doi.org/10.1016/j.parco.2014.03.012


 Hutter, J; Iannuzzi, M; Schiffmann, F; VandeVondele, J. 
 WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE, 4 (1), 15-25 (2014). 
 CP2K: atomistic simulations of condensed matter systems.
 http://dx.doi.org/10.1002/wcms.1159


 Krack, M. 
 THEORETICAL CHEMISTRY ACCOUNTS, 114 (1-3), 145-152 (2005). 
 Pseudopotentials for H to Kr optimized for gradient-corrected
 exchange-correlation functionals.
 http://dx.doi.org/10.1007/s00214-005-0655-y


 VandeVondele, J; Krack, M; Mohamed, F; Parrinello, M; Chassaing, T;
 Hutter, J. COMPUTER PHYSICS COMMUNICATIONS, 167 (2), 103-128 (2005). 
 QUICKSTEP: Fast and accurate density functional calculations using a
 mixed Gaussian and plane waves approach.
 http://dx.doi.org/10.1016/j.cpc.2004.12.014


 Frigo, M; Johnson, SG. 
 PROCEEDINGS OF THE IEEE, 93 (2), 216-231 (2005). 
 The design and implementation of FFTW3.
 http://dx.doi.org/10.1109/JPROC.2004.840301


 Hartwigsen, C; Goedecker, S; Hutter, J. 
 PHYSICAL REVIEW B, 58 (7), 3641-3662 (1998). 
 Relativistic separable dual-space Gaussian pseudopotentials from H to Rn.
 http://dx.doi.org/10.1103/PhysRevB.58.3641


 Lippert, G; Hutter, J; Parrinello, M. 
 MOLECULAR PHYSICS, 92 (3), 477-487 (1997). 
 A hybrid Gaussian and plane wave density functional scheme.
 http://dx.doi.org/10.1080/002689797170220


 Perdew, JP; Burke, K; Ernzerhof, M. 
 PHYSICAL REVIEW LETTERS, 77 (18), 3865-3868 (1996). 
 Generalized gradient approximation made simple.
 http://dx.doi.org/10.1103/PhysRevLett.77.3865


 Goedecker, S; Teter, M; Hutter, J. 
 PHYSICAL REVIEW B, 54 (3), 1703-1710 (1996). 
 Separable dual-space Gaussian pseudopotentials.
 http://dx.doi.org/10.1103/PhysRevB.54.1703


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                T I M I N G                                  -
 -                                                                             -
 -------------------------------------------------------------------------------
 SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    0.008    0.008  101.536  101.536
 cp_cell_opt                          1  2.0    0.000    0.000  101.482  101.482
 geoopt_bfgs                          1  3.0    0.001    0.001  101.481  101.481
 cp_eval_at                           6  4.0    0.001    0.001  101.476  101.476
 qs_energies                          6  5.8    0.009    0.009   96.070   96.070
 scf_env_do_scf                       6  6.8    0.000    0.000   93.674   93.674
 scf_env_do_scf_inner_loop           48  7.8    0.016    0.016   93.674   93.674
 qs_forces                            5  5.0    0.002    0.002   86.009   86.009
 qs_scf_new_mos_kp                   48  8.8    0.000    0.000   41.202   41.202
 do_general_diag_kp                  48  9.8    0.472    0.472   41.201   41.201
 rebuild_ks_matrix                   53  9.6    0.000    0.000   30.080   30.080
 qs_ks_build_kohn_sham_matrix        53 10.6    0.034    0.034   30.080   30.080
 sum_up_and_integrate                53 11.6    0.000    0.000   29.876   29.876
 integrate_v_rspace                  53 12.6   28.935   28.935   29.876   29.876
 qs_rho_update_rho                   54  8.8    0.000    0.000   26.456   26.456
 calculate_rho_elec                  54  9.8   26.441   26.441   26.456   26.456
 qs_ks_update_qs_env                 48  8.8    0.000    0.000   25.682   25.682
 rskp_transform                   48000 10.8   15.279   15.279   15.279   15.279
 kpoint_density_transform            53 10.6    0.163    0.163   11.302   11.302
 copy_dbcsr_to_fm                 98898 10.8    0.645    0.645    9.582    9.582
 transform_dmat                   26500 11.6    8.052    8.052    8.052    8.052
 dbcsr_desymmetrize_deep         194898 11.3    2.045    2.045    6.474    6.474
 dbcsr_complete_redistribute     151898 12.1    2.440    2.440    6.062    6.062
 qs_ks_update_qs_env_forces           5  6.0    0.000    0.000    4.507    4.507
 dbcsr_finalize                  515784 12.4    1.003    1.003    4.057    4.057
 copy_fm_to_dbcsr                 53000 11.6    0.126    0.126    2.948    2.948
 -------------------------------------------------------------------------------

 The number of warnings for this run is : 1
 
 -------------------------------------------------------------------------------
  **** **** ******  **  PROGRAM ENDED AT                 2020-07-15 11:34:55.649
 ***** ** ***  *** **   PROGRAM RAN ON                                  mt-rf614
 **    ****   ******    PROGRAM RAN BY                                     rf614
 ***** **    ** ** **   PROGRAM PROCESS ID                                140639
  **** **  *******  **  PROGRAM STOPPED IN /media/ssd1/rf614/Work/Documents/jobs
                                           /Corrosion_Work/random/learning_cp2k/
                                           mg_opt/fcc_opt/constrained_angles/tem
                                           p
"""




if __name__ == '__main__':
	unittest.main()


