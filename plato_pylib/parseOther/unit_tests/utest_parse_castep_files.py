#!/usr/bin/python3

import itertools
import unittest
import os
import sys
import shutil

sys.path.append('..')


import plato_pylib.parseOther.parse_castep_files as tCode
import plato_pylib.shared.ucell_class as UCell

class testParseCastepOutFile(unittest.TestCase):

	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["fileA"], unusedStr = createCastepOutfileNaCl()
		self.filePathDict["fileB"], unusedStr = createCastepOutfileHcpMgPartial()

	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

	def testParseEnergies(self):
		''' test tCode.parseCastepOutfile correctly reads in total energy '''
		filePathKeys = ["fileA", "fileB"]
		expectedEnergiesDict = {"fileA":-1713.025929482, "fileB":-3192.239986240}

		actualEnergiesDict = dict()
		for fileKey in filePathKeys:
			currFilePath = self.filePathDict[fileKey]
			actualEnergiesDict[fileKey] = tCode.parseCastepOutfile(currFilePath)["energy"]

		for fileKey in filePathKeys:
			self.assertAlmostEqual( expectedEnergiesDict[fileKey], actualEnergiesDict[fileKey] )


	def testParseNumberAtoms(self):
		''' test tCode.parseCastepOutfile correctly reads number of atoms '''
		
		filePathKeys = ["fileA", "fileB"]
		expectedNumberAtomsDict = {"fileA":2, "fileB":2}

		actualNumberAtomsDict = dict()
		for fileKey in filePathKeys:
			currFilePath = self.filePathDict[fileKey]
			actualNumberAtomsDict[fileKey] = tCode.parseCastepOutfile(currFilePath)["numbAtoms"]

		for fileKey in filePathKeys:
			self.assertAlmostEqual( expectedNumberAtomsDict[fileKey], actualNumberAtomsDict[fileKey] )


	def testParseScfKgrid(self):
		''' test tCode.parseCastepOutfile correctly reads info for grid used in scf-calc '''

		filePaths = [self.filePathDict["fileA"], self.filePathDict["fileB"]] #[NaCl, hcp]
		expectedKgrids = [[4,4,4], [4,4,2]]
		expectedNumbK = [10,6]

		parsedFiles = [tCode.parseCastepOutfile(x) for x in filePaths]
		actualKgrids = [x["scf_kgrid"] for x in parsedFiles]
		actualNumbK = [x["scf_numb_k"] for x in parsedFiles]

		for expGrid, actGrid in itertools.zip_longest(expectedKgrids, actualKgrids):
			[self.assertEqual(exp,act) for exp, act in itertools.zip_longest(expGrid, actGrid) ]

		[self.assertEqual(exp,act) for exp,act in itertools.zip_longest(expectedNumbK, actualNumbK)]

	def testParseUnitCellParams(self):
		''' test parseCastepOutFile correctly reads unit cell parameters '''

		filePathKeys = ["fileA", "fileB"]

		expectedLattParams = { "fileA": {"a":3.746172, "b":3.746172, "c":3.746172} ,
		                       "fileB": {"a":3.209401, "b":3.209331, "c":5.210802} }

		expectedLattAngles = {"fileA": {"alpha":60.0, "beta":60.0, "gamma":60.0},
		                      "fileB": {"alpha":90.0, "beta":90.0, "gamma":120.000728}}

		expectedLatticeVects = {"fileA": [ [3.746172,0,0], [1.873086,3.2442801189,0], [1.873086,1.0814267063,3.0587366296] ],
		                        "fileB": [ [3.209401, 0.0, 0.0], [-1.6047008144378012, 2.779341786053608, 0.0], [0.0, 0.0, 5.210802] ]}

		expectedFractCoords = {"fileA": [ [0.0,0.0,0.0,"Na"], [0.70710678,0.40824829,0.28867513,"Cl"] ],
		                       "fileB": [ [0.0,0.0,0.0,"Mg"], [0.333333,0.666667,0.50,"Mg"] ]}

		actualLattParams, actualLattAngles, actualLatticeVects, actualFractCoords = dict(), dict(), dict(), dict()
		for fileKey in filePathKeys:
			currFilePath = self.filePathDict[fileKey]
			currObj = tCode.parseCastepOutfile(currFilePath)["unitCell"]
			actualLattParams[fileKey], actualLattAngles[fileKey] = currObj.lattParams, currObj.lattAngles
			actualFractCoords[fileKey] = currObj.fractCoords	
			actualLatticeVects[fileKey] = currObj.getLattVects()
	
		for filekey in filePathKeys:
			for pA,pB in itertools.zip_longest(expectedLattParams[filekey], actualLattParams[filekey]):
				self.assertAlmostEqual(pA,pB)
			for angA,angB in itertools.zip_longest(expectedLattAngles[filekey],actualLattAngles[filekey]):
				self.assertAlmostEqual(angA,angB)
			for vectA,vectB in itertools.zip_longest(expectedLatticeVects[filekey], actualLatticeVects[filekey]):
				[self.assertAlmostEqual(a,b,places=5) for a,b in itertools.zip_longest(vectA,vectB)]
			for coA,coB in itertools.zip_longest(expectedFractCoords[filekey], actualFractCoords[filekey]):
				[self.assertAlmostEqual(a,b, delta=1e-7) for a,b in itertools.zip_longest(coA,coB)]

	def testParseKinCutoff(self):
		filePaths = [ self.filePathDict["fileA"], self.filePathDict["fileB"] ]
		expectedCutoffs = [330.0, 950.0]
		actualCutoffs = [ tCode.parseCastepOutfile(x)["kin_cut"] for x in filePaths ]
		[self.assertAlmostEqual(exp,act) for exp,act in itertools.zip_longest(expectedCutoffs,actualCutoffs)]



class testParseCellFile(unittest.TestCase):

	def setUp(self):
		self.testCellFileA = createCastepCellFileA()
	def tearDown(self):
		os.remove(self.testCellFileA)

	def testVsKnownVals(self):
		expectedVals = {"%block lattice_cart": "bohr\n 6.0649000000	0.0000000000	0.0000000000\n-3.0325000000	5.2524000000	0.0000000000\n 0.0000000000	0.0000000000	9.8470000000",
		                "%block cell_constraints":"1   2   3\n  0   0   0",
		                "%block positions_frac":"Mg 0.0    0.0    0.0\nMg 0.33333333    0.66666667    0.5",
		                "%block species_pot":"Mg Mg_OTF_PBE_mine.usp",
		                "symmetry_generate": "",
		                "kpoint_mp_grid": "10 10 6",
		                "%block bs_kpoint_list":"0.000 0.000 0.000 1.000"}

		actualVals = tCode.tokenizeCastepCellFile( self.testCellFileA )

		for key in expectedVals.keys():
			self.assertEqual(expectedVals[key], actualVals[key])

	def testWriteCastepFileFromTokens(self):
		testOutFile = "fakeCellFile.txt"
		testTokens = tCode.tokenizeCastepCellFile(self.testCellFileA)
		expectedVals = testTokens

		newFile = tCode.writeCastepCellFileFromTokens(testOutFile, testTokens)
		actualVals = tCode.tokenizeCastepCellFile(testOutFile)

		for key in expectedVals.keys():
			self.assertEqual(expectedVals[key], actualVals[key])		

		os.remove(testOutFile)

	def testModifyCastepCellFile_blockmod(self):
		testField = "cell_constraints"
		testVal = "bohr\n 0.00 0.00 0.00\n 0.00 0.00 0.00\n 0.00 0.00 0.00\n"
		testFieldValDict = {testField:testVal}
		expectedVal = testVal.strip()
	
		tCode.modCastepCellFile(self.testCellFileA, testFieldValDict)
		tokFile = tCode.tokenizeCastepCellFile(self.testCellFileA)
		actualVal =  tokFile["%block " + testField.lower()]

		self.assertEqual(expectedVal, actualVal)

	def testModifyCastepCellFile_singleLineMod(self):
		testField = "symmetry_generate"
		testVal = "Hello".lower()
		testFieldValDict = {testField:testVal}
		expectedVal = testVal.strip()

		tCode.modCastepCellFile(self.testCellFileA, testFieldValDict)
		actualVal = tCode.tokenizeCastepCellFile(self.testCellFileA)[testField.lower()]

		self.assertEqual(expectedVal, actualVal)

	def testModifyCellFile_addNewField(self):
		''' Test correctly modifies a castep file with a field not originally present '''
		testField = "fake_field"
		testVal = "this is fake"
		testFieldValDict = {testField:testVal}
		expectedVal = testVal.strip()

		tCode.modCastepCellFile(self.testCellFileA, testFieldValDict)
		actualVal = tCode.tokenizeCastepCellFile(self.testCellFileA)[testField.lower()]

		self.assertEqual(expectedVal, actualVal)

	def testModifyCellFile_multipleFields(self):
		testFieldValDict = {"fakeA":"fake", "fakeB":"fake"}
		expectedVals = testFieldValDict

		tCode.modCastepCellFile(self.testCellFileA, testFieldValDict)
		actualVals = tCode.tokenizeCastepCellFile(self.testCellFileA)

		print("actualVals = {}".format(actualVals))

		for key in expectedVals.keys():
			self.assertEqual(expectedVals[key], actualVals[key.lower()])

	def testGetUnitCellObjFromCellFile(self):
		lattVects = [ [ 6.0649000000, 0.0000000000, 0.0000000000],
		              [-3.0325000000, 5.2524000000, 0.0000000000],
		              [ 0.0000000000, 0.0000000000, 9.8470000000] ]

		expFractCoords = [ [0.0, 0.0, 0.0, "Mg"],
		                   [0.33333333, 0.66666667, 0.5, "Mg"] ]

		expObj = UCell.UnitCell.fromLattVects(lattVects, fractCoords = expFractCoords )
		actObj = tCode.unitCellObjFromCastepCellFile(self.testCellFileA)
		self.assertEqual(expObj,actObj)


class testGeomStringsFromUnitCell(unittest.TestCase):
	''' Test we can correctly get strings for castep cell file from the UnitCell Class '''
	def setUp(self):
		testLattParamsA = [3.209401, 3.209331, 5.210802] #Actually angstrom, but code assumes bohr and no need to convert in this test code
		testLattAnglesA = [90.0, 90.0, 120.000728]
		testFractCoordsA = [ [0.0,0.0,0.0,"Mg"], [0.333333,0.666667,0.50,"Mg"] ]
		self.testCellA = UCell.UnitCell(lattParams=testLattParamsA , lattAngles=testLattAnglesA)
		self.testCellA.fractCoords = testFractCoordsA
		self.expLattCartStrA = "bohr\n3.209401 0 0\n-1.604701 2.779342 0\n0 0 5.210802"
		self.expFractCoordsStrA = "Mg 0.000000 0.000000 0.000000\nMg 0.333333 0.666667 0.500000"
		self.expDictSectionA = {"lattice_cart": self.expLattCartStrA,
		                        "positions_frac": self.expFractCoordsStrA}

	def testGetCellGeomDictSectionFromUCell(self):
		actDict = tCode.getCellGeomDictSectionFromUCell(self.testCellA)
		for key in self.expDictSectionA.keys():
			self.assertEqual( self.expDictSectionA[key], actDict[key] )

	def testGetFractPosBlockStr(self):
		actStr = tCode._getCastepFractPosStrFromUCellClass(self.testCellA)
		self.assertEqual(self.expFractCoordsStrA,actStr)

	def testGetLattCartBlock(self):
		actStr = tCode._getCastepLattCartFromUCellClass(self.testCellA)
		self.assertEqual(self.expLattCartStrA,actStr)



class testParseCastepGeomFile(unittest.TestCase):

	def setUp(self):
		self.testFileA = createCastepGeomFileA()

	def tearDown(self):
		os.remove(self.testFileA)


	def testFirstIterParsingUCell(self):
		expLattVects = [ [6.06, 0.0, 0.0],
		                 [-3.0299497916016769, 5.2480560910267373, 0.0],
		                 [0.0, 0.0, 9.8390435291853642] ]
		
		expFractCoords =  [ [0.0, 0.0, 0.0, "Mg"],
		                    [0.33333333, 0.66666667, 0.5, "Mg"] ]

		expUCell = UCell.UnitCell.fromLattVects(expLattVects, fractCoords=expFractCoords)
		actUCell = tCode.parseCastepGeomFile(self.testFileA)[0]["unitCell"]

		self.assertEqual(expUCell,actUCell)


class testParseCastepParamFile(unittest.TestCase):

	def setUp(self):
		self.testFileA = createCastepParamFileA()
		self.partialDictFileA = {"task":"geometryoptimisation",
		                          "xc_functional": "PBE".lower(),
		                          "cut_off_energy":"600"}
	def tearDown(self):
		os.remove(self.testFileA)

	@unittest.skip("")
	def testTokenizeCastepParamFile(self):
		actFullDict = tCode.tokenizeCastepParamFile(self.testFileA)
		for key in self.partialDictFileA.keys():
			self.assertEqual( self.partialDictFileA[key], actFullDict[key])

	def testModValInParamFile_cutoffEnergy(self):
		testCutoff = 500
		self.partialDictFileA["cut_off_energy"] = testCutoff
		tCode.modCastepParamFile(self.testFileA,self.partialDictFileA) #Testing this function
		self.partialDictFileA["cut_off_energy"] = str(testCutoff) #Wanted it as int to more fully test funct, but expected result is for it to be a str
		tokenizedFile = tCode.tokenizeCastepParamFile(self.testFileA) #Using this function only to help test the other
		for key in self.partialDictFileA.keys():
			self.assertEqual( self.partialDictFileA[key], tokenizedFile[key] )



def createCastepCellFileA():
	filePath = os.path.join(os.getcwd(), "cellFileA.cell")
	fileStr = "# Change lattice parameters here\n%BLOCK LATTICE_CART\nbohr\n 6.0649000000	0.0000000000	0.0000000000\n-3.0325000000	5.2524000000	0.0000000000\n 0.0000000000	0.0000000000	9.8470000000\n%ENDBLOCK LATTICE_CART\n\n%BLOCK CELL_CONSTRAINTS\n  1   2   3\n  0   0   0\n%ENDBLOCK CELL_CONSTRAINTS\n\n# Change elements here\n%BLOCK POSITIONS_FRAC\nMg 0.0    0.0    0.0\nMg 0.33333333    0.66666667    0.5\n%ENDBLOCK POSITIONS_FRAC\n\n# You will need to get a potential file for magnesium and have it in the same folder (and refernce to it here)\n# I am pretty sure you can get them from here: http://cmt.dur.ac.uk/Pseudopotentials/\n%BLOCK SPECIES_POT\nMg Mg_OTF_PBE_mine.usp\n%ENDBLOCK SPECIES_POT\n\nsymmetry_generate\n\nkpoint_mp_grid 10 10 6\n\n#Your path though k-space from Plato (without the weight)\n%BLOCK BS_KPOINT_LIST\n0.000 0.000 0.000 1.000\n%ENDBLOCK BS_KPOINT_LIST\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath


def createCastepParamFileA():
	filePath = os.path.abspath( os.path.join(os.getcwd(), "paramFileA.param" ) )
	fileStr = "task :				geometryoptimisation\nxc_functional :			PBE\nspin_polarized :		false\nspin :				0\ncut_off_energy :		600\ngrid_scale :			2.000000000000000\nfine_grid_scale :		2.300000000000000\nfinite_basis_corr :		2\nelec_energy_tol :		1.000000000000000e-006\nmax_scf_cycles :		150\nfix_occupancy :			false\nmetals_method :			dm\nmixing_scheme :			Pulay\nsmearing_scheme: 		FermiDirac\nsmearing_width :        	0.0136 eV\nmix_charge_amp :		0.500000000000000\nmix_spin_amp :			2.000000000000000\nmix_charge_gmax :		1.500000000000000\nmix_spin_gmax :			1.500000000000000\nmix_history_length :		20\nperc_extra_bands :		40\nspin_fix :			6\ngeom_energy_tol :		1.000000000000000e-005\ngeom_force_tol :		0.010000000000000\ngeom_stress_tol :		0.050000000000000\ngeom_disp_tol :			5.000000000000000e-004\ngeom_max_iter :			200\ngeom_method :			BFGS\nfixed_npw :			false\ngeom_modulus_est :      	100.000000000000000  GPa\ncalculate_ELF :			false\ncalculate_stress :		true\npopn_calculate :		true\ncalculate_hirshfeld :		true\ncalculate_densdiff :		false\npdos_calculate_weights : 	false\nnum_dump_cycles :		0\nnum_backup_iter :		2\nopt_strategy :			speed\npage_wvfns :			0\nrun_time :			259100\ncharge :			0\n\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath

def createCastepGeomFileA():
	filePath = os.path.join(os.getcwd(),"geomFileA.geom")
	fileStr = " BEGIN header\n  \n END header\n  \n                                      0\n                     -1.2423646407988753E+002   -1.2423646407988753E+002                             <-- E\n                      6.0599999999999996E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- h\n                     -3.0299497916016769E+000    5.2480560910267373E+000    0.0000000000000000E+000  <-- h\n                      6.0246765823441672E-016    1.0435017001568253E-015    9.8390435291853642E+000  <-- h\n                      1.2536669721016249E-005    1.1993626171350869E-010    1.3003737642768126E-021  <-- S\n                      1.1993626171350869E-010    1.2536393309559805E-005    2.2523071824947772E-021  <-- S\n                      1.3003737642768126E-021    2.2523071824947772E-021    2.1236715193020276E-005  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3441965716053763E-005    3.4987040781780125E+000    4.9195217645926821E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      1\n                     -1.2423648361329813E+002   -1.2423648361329813E+002                             <-- E\n                      6.0525494010867629E+000   -7.1278657017122572E-008   -7.7281793962142931E-019  <-- h\n                     -3.0262246156033141E+000    5.2416039321280552E+000   -7.7280834840665398E-019  <-- h\n                      6.0047217826623655E-016    1.0400454795107369E-015    9.8185518948892216E+000  <-- h\n                      6.2257935287971063E-006    5.9487827644999277E-011    8.2949135219161241E-022  <-- S\n                      5.9487827644999277E-011    6.2256564300010843E-006    1.4371126374083785E-021  <-- S\n                      8.2949135219161241E-022    1.4371126374083785E-021    1.3563332628283939E-005  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3359697465174452E-005    3.4944026151311651E+000    4.9092759474446108E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      2\n                     -1.2423649139311006E+002   -1.2423649139311006E+002                             <-- E\n                      6.0488616295277797E+000   -1.0655896228802952E-007   -1.2651835524717242E-018  <-- h\n                     -3.0243807909311040E+000    5.2384103515370590E+000   -1.2651678506530170E-018  <-- h\n                      5.9930613757194439E-016    1.0380258644012049E-015    9.8054965870848179E+000  <-- h\n                      1.2006137161184351E-006    1.1464924393711919E-011    4.5041431156489958E-022  <-- S\n                      1.1464924393711919E-011    1.2005872934468897E-006    7.8114485655588434E-022  <-- S\n                      4.5041431156489958E-022    7.8114485655588434E-022    7.3694155923667525E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3318977715823709E-005    3.4922735496330870E+000    4.9027482935424089E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      3\n                     -1.2423649401879969E+002   -1.2423649401879969E+002                             <-- E\n                      6.0481518391308562E+000   -1.1334941052519963E-007   -1.5323764873051865E-018  <-- h\n                     -3.0240259074940430E+000    5.2377956787464504E+000   -1.5323574694340257E-018  <-- h\n                      5.9880175543322601E-016    1.0371522553781538E-015    9.7984118396606021E+000  <-- h\n                      2.8493477471534414E-006    2.7205843346097138E-011    4.9308317826987391E-022  <-- S\n                      2.7205843346097138E-011    2.8492850471269788E-006    8.5447239829073007E-022  <-- S\n                      4.9308317826987391E-022    8.5447239829073007E-022    8.0685001472674851E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3311140330860220E-005    3.4918637655071496E+000    4.8992059198303011E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      4\n                     -1.2423649664374317E+002   -1.2423649664374317E+002                             <-- E\n                      6.0464687490037869E+000   -1.2945125754357399E-007   -1.8248474839828538E-018  <-- h\n                     -3.0231843903197255E+000    5.2363381361962311E+000   -1.8248248363406947E-018  <-- h\n                      5.9815956747330960E-016    1.0360399658797643E-015    9.7906568324346015E+000  <-- h\n                      7.9122364917404268E-007    7.5525675479553870E-012    3.2066121442602407E-022  <-- S\n                      7.5525675479553870E-012    7.9120624312714741E-007    5.5623547306507955E-022  <-- S\n                      3.2066121442602407E-022    5.5623547306507955E-022    5.2485725893818598E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3292555935059847E-005    3.4908920651015296E+000    4.8953284162173007E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      5\n                     -1.2423649781008835E+002   -1.2423649781008835E+002                             <-- E\n                      6.0460018783944767E+000   -1.3391773204044196E-007   -2.0149944754942537E-018  <-- h\n                     -3.0229509627512314E+000    5.2359338299932707E+000   -2.0149694679940464E-018  <-- h\n                      5.9780442826045685E-016    1.0354248510857768E-015    9.7856149947483910E+000  <-- h\n                      2.3120746185618186E-006    2.2068027558973438E-011    3.7319692460071080E-022  <-- S\n                      2.2068027558973438E-011    2.3120237594175427E-006    6.4656440760634029E-022  <-- S\n                      3.7319692460071080E-022    6.4656440760634029E-022    6.1089567904230697E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3287400828939913E-005    3.4906225261427171E+000    4.8928074973741955E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      6\n                     -1.2423649899035129E+002   -1.2423649899035129E+002                             <-- E\n                      6.0446384202607133E+000   -1.4696171223657286E-007   -2.2362775693493142E-018  <-- h\n                     -3.0222692562771862E+000    5.2347530861855320E+000   -2.2362498155701058E-018  <-- h\n                      5.9730959871704475E-016    1.0345677913223748E-015    9.7797475685374806E+000  <-- h\n                     -3.7803991756630201E-007   -3.6074548648313439E-012    1.7744046298961861E-022  <-- S\n                     -3.6074548648313439E-012   -3.7803160363489407E-007    3.0836309320724544E-022  <-- S\n                      1.7744046298961861E-022    3.0836309320724544E-022    2.9052319604609918E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3272345754863382E-005    3.4898353592522953E+000    4.8898737842687403E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      7\n                     -1.2423649932432403E+002   -1.2423649932432403E+002                             <-- E\n                      6.0448611708534683E+000   -1.4483069387331832E-007   -2.3414655446787083E-018  <-- h\n                     -3.0223806278825318E+000    5.2349459864122077E+000   -2.3414364854426625E-018  <-- h\n                      5.9716096062854568E-016    1.0343103420456952E-015    9.7769584594701406E+000  <-- h\n                      2.0671328685530037E-007    1.9726394915202643E-012    1.6866155432187839E-022  <-- S\n                      1.9726394915202643E-012    2.0670874060592522E-007    2.9272717150145604E-022  <-- S\n                      1.6866155432187839E-022    2.9272717150145604E-022    2.7613945301766011E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3274805329185297E-005    3.4899639601143946E+000    4.8884792297350703E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n                                      8\n                     -1.2423649962350710E+002   -1.2423649962350710E+002                             <-- E\n                      6.0447394002326131E+000   -1.4599565341815054E-007   -2.4414530539991855E-018  <-- h\n                     -3.0223197445898733E+000    5.2348405340312514E+000   -2.4414227538477481E-018  <-- h\n                      5.9698651417095413E-016    1.0340081934498078E-015    9.7743072432261435E+000  <-- h\n                     -5.8679611248907208E-007   -5.5996098148374308E-012    9.3741286500518324E-023  <-- S\n                     -5.5996098148374308E-012   -5.8678320733192977E-007    1.6318735266550080E-022  <-- S\n                      9.3741286500518324E-023    1.6318735266550080E-022    1.5348020665153149E-006  <-- S\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- R\n Mg              2    3.3273460758508189E-005    3.4898936581384192E+000    4.8871536216130718E+000  <-- R\n Mg              1    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n Mg              2    0.0000000000000000E+000    0.0000000000000000E+000    0.0000000000000000E+000  <-- F\n  \n"

	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createCastepOutfileNaCl():
	fileName = "NaCl.castep"
	filePath = os.path.join( os.getcwd(), fileName )
	fileStr = "\n \n Compiled for linux_x86_64_ifort15 on Tue, 24 Feb 2015 12:57:19 +0000\n from code version 7a66ba17f5c9+  Mon, 08 Dec 2014 20:23:19 +0000\n Compiler: Intel Fortran 15.0.1.133; Optimisation: intermediate\n MATHLIBS: Intel MKL(11.2.1) (LAPACK version 3.5.0)\n FFT Lib : mkl\n Fundamental constants values: CODATA 2010\n \n Run started: Wed, 07 Mar 2018 16:13:52 +0000\n \n Pseudo atomic calculation performed for Na 2s2 2p6 3s1\n \n Converged in 22 iterations to a total energy of -1299.0842 eV\n \n \n Pseudo atomic calculation performed for Cl 3s2 3p5\n \n Converged in 19 iterations to a total energy of -405.9090 eV\n \n Calculation parallelised over 5 processes.\n Data is distributed by k-point(5-way)\n\n ************************************ Title ************************************\n \n\n ***************************** General Parameters ******************************\n  \n output verbosity                               : normal  (1)\n write checkpoint data to                       : /var/tmp/pbs.1430725.cx1/NaCl.check\n type of calculation                            : single point energy\n stress calculation                             : off\n density difference calculation                 : on\n electron localisation func (ELF) calculation   : off\n Hirshfeld analysis                             : off\n unlimited duration calculation\n timing information                             : on\n memory usage estimate                          : on\n write final potential to formatted file        : off\n write final density to formatted file          : on\n write BibTeX reference list                    : on\n checkpoint writing                             : both castep_bin and check files\n  \n output         length unit                     : A\n output           mass unit                     : amu\n output           time unit                     : ps\n output         charge unit                     : e\n output           spin unit                     : hbar/2\n output         energy unit                     : eV\n output          force unit                     : eV/A\n output       velocity unit                     : A/ps\n output       pressure unit                     : GPa\n output     inv_length unit                     : 1/A\n output      frequency unit                     : cm-1\n output force constant unit                     : eV/A**2\n output         volume unit                     : A**3\n output   IR intensity unit                     : (D/A)**2/amu\n output         dipole unit                     : D\n output         efield unit                     : eV/A/e\n output        entropy unit                     : J/mol/K\n  \n wavefunctions paging                           : none\n random number generator seed                   : randomised (161352628)\n data distribution                              : optimal for this architecture\n optimization strategy                          : maximize speed(+++)\n\n *********************** Exchange-Correlation Parameters ***********************\n  \n using functional                               : Local Density Approximation\n Divergence correction                          : off\n relativistic treatment                         : Koelling-Harmon\n DFT+D: Semi-empirical dispersion correction    : off\n\n ************************* Pseudopotential Parameters **************************\n  \n pseudopotential representation                 : reciprocal space\n <beta|phi> representation                      : reciprocal space\n\n **************************** Basis Set Parameters *****************************\n  \n basis set accuracy                             : MEDIUM\n plane wave basis set cut-off                   :   330.0000   eV\n size of standard grid                          :     1.7500\n size of   fine   gmax                          :    16.2867   1/A\n finite basis set correction                    : none\n\n **************************** Electronic Parameters ****************************\n  \n number of  electrons                           :  16.00    \n net charge of system                           :  0.000    \n net spin   of system                           :  0.000    \n number of  up  spins                           :  8.000    \n number of down spins                           :  8.000    \n treating system as non-spin-polarized\n number of bands                                :          8\n\n ********************* Electronic Minimization Parameters **********************\n  \n Method: Treating system as non-metallic,\n         and number of  SD  steps               :          1\n         and number of  CG  steps               :          4\n  \n total energy / atom convergence tol.           : 0.1000E-04   eV\n max force / atom convergence tol.              : ignored\n convergence tolerance window                   :          3   cycles\n max. number of SCF cycles                      :         30\n periodic dipole correction                     : NONE\n\n *********************** Population Analysis Parameters ************************\n  \n Population analysis with cutoff                :  3.000       A\n\n *******************************************************************************\n  \n \n                           -------------------------------\n                                      Unit Cell\n                           -------------------------------\n        Real Lattice(A)                      Reciprocal Lattice(1/A)\n   2.6489433   2.6489433   0.0000000        1.1859796   1.1859796  -1.1859796\n   2.6489433   0.0000000   2.6489433        1.1859796  -1.1859796   1.1859796\n   0.0000000   2.6489433   2.6489433       -1.1859796   1.1859796   1.1859796\n \n                       Lattice parameters(A)       Cell Angles\n                    a =    3.746172          alpha =   60.000000\n                    b =    3.746172          beta  =   60.000000\n                    c =    3.746172          gamma =   60.000000\n \n                       Current cell volume =   37.174744       A**3\n \n                           -------------------------------\n                                     Cell Contents\n                           -------------------------------\n                         Total number of ions in cell =    2\n                      Total number of species in cell =    2\n                        Max number of any one species =    1\n \n            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n            x  Element    Atom        Fractional coordinates of atoms  x\n            x            Number           u          v          w      x\n            x----------------------------------------------------------x\n            x  Na           1         0.000000   0.000000   0.000000   x \n            x  Cl           1         0.500000   0.500000   0.500000   x \n            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n \n \n                         No user defined ionic velocities\n \n                           -------------------------------\n                                   Details of Species\n                           -------------------------------\n \n                               Mass of species in AMU\n                                    Na   22.9897700\n                                    Cl   35.4527000\n \n                          Electric Quadrupole Moment (Barn)\n                                    Na    0.1040000 Isotope 23\n                                    Cl   -0.0816500 Isotope 35\n \n                          Files used for pseudopotentials:\n                                    Na Na_00.usp\n                                    Cl Cl_00.usp\n \n                           -------------------------------\n                              k-Points For BZ Sampling\n                           -------------------------------\n                       MP grid size for SCF calculation is  4  4  4\n                            with an offset of   0.000  0.000  0.000\n                       Number of kpoints used =            10\n \n                           -------------------------------\n                               Symmetry and Constraints\n                           -------------------------------\n \n                      Maximum deviation from symmetry =  0.00000         ANG\n \n                      Number of symmetry operations   =          48\n                      Number of ionic constraints     =           3\n                      Point group of crystal =    32: Oh, m-3m, 4/m -3 2/m\n                      Space group of crystal =   225: Fm-3m, -F 4 2 3\n \n             Set iprint > 1 for details on symmetry rotations/translations\n \n                         Centre of mass is constrained\n             Set iprint > 1 for details of linear ionic constraints\n \n                         Number of cell constraints= 5\n                         Cell constraints are:  1 1 1 0 0 0\n \n                         External pressure/stress (GPa)\n                          0.00000   0.00000   0.00000\n                                    0.00000   0.00000\n                                              0.00000\n  \n+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER PROCESS --------------+\n|                                                     Memory          Disk    |\n| Model and support data                               20.3 MB         0.0 MB |\n| Electronic energy minimisation requirements           1.2 MB         0.0 MB |\n|                                               ----------------------------- |\n| Approx. total storage required per process           21.6 MB         0.0 MB |\n|                                                                             |\n| Requirements will fluctuate during execution and may exceed these estimates |\n+-----------------------------------------------------------------------------+\n------------------------------------------------------------------------ <-- SCF\nSCF loop      Energy                           Energy gain       Timer   <-- SCF\n                                               per atom          (sec)   <-- SCF\n------------------------------------------------------------------------ <-- SCF\nInitial  -2.76157985E+002                                          0.65  <-- SCF\n      1  -1.69343189E+003                    7.08636953E+002       0.80  <-- SCF\n      2  -1.71277231E+003                    9.67021035E+000       0.95  <-- SCF\n      3  -1.71302453E+003                    1.26111679E-001       1.09  <-- SCF\n      4  -1.71302592E+003                    6.91259192E-004       1.23  <-- SCF\n      5  -1.71302593E+003                    6.17634417E-006       1.38  <-- SCF\n      6  -1.71302593E+003                    8.70912307E-008       1.52  <-- SCF\n------------------------------------------------------------------------ <-- SCF\n \nFinal energy =  -1713.025929482     eV\n(energy not corrected for finite basis set)\n \n\nWriting analysis data to /var/tmp/pbs.1430725.cx1/NaCl.castep_bin\n\nWriting model to /var/tmp/pbs.1430725.cx1/NaCl.check\n \n ***************** Symmetrised Forces *****************\n *                                                    *\n *            Cartesian components (eV/A)             *\n * -------------------------------------------------- *\n *                   x            y            z      *\n *                                                    *\n * Na        1      0.00000      0.00000      0.00000 *\n * Cl        1      0.00000      0.00000      0.00000 *\n *                                                    *\n ******************************************************\n \n Pseudo atomic calculation performed for Na 2s2 2p6 3s1\n \n Converged in 22 iterations to a total energy of -1299.0842 eV\n \n \n Pseudo atomic calculation performed for Cl 3s2 3p5\n \n Converged in 19 iterations to a total energy of -405.9090 eV\n \nCharge spilling parameter for spin component 1 = 0.07%\n \n     Atomic Populations (Mulliken)\n     -----------------------------\nSpecies   Ion     s      p      d      f     Total  Charge (e)\n==============================================================\n  Na       1     2.21   6.32   0.00   0.00   8.53     0.47\n  Cl       1     1.92   5.54   0.00   0.00   7.47    -0.47\n==============================================================\n \n            Bond              Population      Length (A)\n============================================================\n         Na 1 -- Cl 1              0.77        2.64894\n============================================================\n \n\nWriting analysis data to /var/tmp/pbs.1430725.cx1/NaCl.castep_bin\n\nWriting model to /var/tmp/pbs.1430725.cx1/NaCl.check\n \n A BibTeX formatted list of references used in this run has been written to \n /var/tmp/pbs.1430725.cx1/NaCl.bib\n \nInitialisation time =      0.39 s\nCalculation time    =      1.42 s\nFinalisation time   =      0.05 s\nTotal time          =      1.87 s\nPeak Memory Use     = 126084 kB\n  \nOverall parallel efficiency rating: Very good (80%)                             \n  \nData was distributed by:-\nk-point (5-way); efficiency rating: Very good (81%)                             \n  \nParallel notes:\n1) Calculation only took 1.8 s, so efficiency estimates may be inaccurate.      \n')\n\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath, fileStr



def createCastepOutfileHcpMgPartial():
	fileName = "hcpMg.castep"
	filePath = os.path.join( os.getcwd(), fileName )
	fileStr = " Compiled for linux_x86_64_ifort15 on Tue, 24 Feb 2015 12:57:19 +0000\n from code version 7a66ba17f5c9+  Mon, 08 Dec 2014 20:23:19 +0000\n Compiler: Intel Fortran 15.0.1.133; Optimisation: intermediate\n MATHLIBS: Intel MKL(11.2.1) (LAPACK version 3.5.0)\n FFT Lib : mkl\n Fundamental constants values: CODATA 2010\n \n Run started: Thu, 15 Mar 2018 14:34:58 +0000\n Warning in parameters_read: it appears you may not have specified enough \n                           - extra bands so you may have trouble with metals\n                           - convergence as a consequence. Suggest you \n                           - consider setting nextra_bands=4 at least.\n \n Pseudo atomic calculation performed for Mg 2s2 2p6 3s2\n \n Converged in 23 iterations to a total energy of -1594.7271 eV\n \n Calculation parallelised over 12 processes.\n Data is distributed by G-vector(2-way) and k-point(6-way)\n\n ************************************ Title ************************************\n \n\n ***************************** General Parameters ******************************\n  \n output verbosity                               : normal  (1)\n write checkpoint data to                       : /var/tmp/pbs.1525038.cx1/Mg_hcp_bands_10el_recpot.check\n type of calculation                            : band structure\n stress calculation                             : on\n density difference calculation                 : off\n electron localisation func (ELF) calculation   : off\n Hirshfeld analysis                             : on\n calculation limited to maximum                 :     259100   seconds\n timing information                             : on\n memory usage estimate                          : on\n write final potential to formatted file        : off\n write final density to formatted file          : off\n write BibTeX reference list                    : on\n checkpoint writing                             : both castep_bin and check files\n  \n output         length unit                     : A\n output           mass unit                     : amu\n output           time unit                     : ps\n output         charge unit                     : e\n output           spin unit                     : hbar/2\n output         energy unit                     : eV\n output          force unit                     : eV/A\n output       velocity unit                     : A/ps\n output       pressure unit                     : GPa\n output     inv_length unit                     : 1/A\n output      frequency unit                     : cm-1\n output force constant unit                     : eV/A**2\n output         volume unit                     : A**3\n output   IR intensity unit                     : (D/A)**2/amu\n output         dipole unit                     : D\n output         efield unit                     : eV/A/e\n output        entropy unit                     : J/mol/K\n  \n wavefunctions paging                           : none\n random number generator seed                   : randomised (143458565)\n data distribution                              : optimal for this architecture\n optimization strategy                          : maximize speed(+++)\n\n *********************** Exchange-Correlation Parameters ***********************\n  \n using functional                               : Perdew Burke Ernzerhof\n Divergence correction                          : off\n relativistic treatment                         : Koelling-Harmon\n DFT+D: Semi-empirical dispersion correction    : off\n\n ************************* Pseudopotential Parameters **************************\n  \n pseudopotential representation                 : reciprocal space\n <beta|phi> representation                      : reciprocal space\n\n **************************** Basis Set Parameters *****************************\n  \n plane wave basis set cut-off                   :   950.0000   eV\n size of standard grid                          :     2.0000\n size of   fine   grid                          :     2.3000\n size of   fine   gmax                          :    36.3185   1/A\n largest prime factor in FFT                    :          5\n finite basis set correction                    : automatic\n number of sample energies                      :          3\n           sample  spacing                      :     5.0000   eV\n\n **************************** Electronic Parameters ****************************\n  \n number of  electrons                           :  20.00    \n net charge of system                           :  0.000    \n net spin   of system                           :  0.000    \n number of  up  spins                           :  10.00    \n number of down spins                           :  10.00    \n treating system as non-spin-polarized\n number of bands                                :         12\n\n ********************* Electronic Minimization Parameters **********************\n  \n Method: Treating system as metallic with density mixing treatment of electrons,\n         and number of  SD  steps               :          1\n         and number of  CG  steps               :          4\n  \n total energy / atom convergence tol.           : 0.1000E-05   eV\n eigen-energy convergence tolerance             : 0.1667E-06   eV\n max force / atom convergence tol.              : ignored\n convergence tolerance window                   :          3   cycles\n max. number of SCF cycles                      :        150\n number of fixed-spin iterations                :          6\n smearing scheme                                : cold smearing\n smearing width                                 : 0.1000       eV\n Fermi energy convergence tolerance             : 0.2721E-13   eV\n periodic dipole correction                     : NONE\n\n ************************** Density Mixing Parameters **************************\n  \n density-mixing scheme                          : Pulay\n max. length of mixing history                  :         20\n charge density mixing amplitude                : 0.5000    \n cut-off energy for mixing                      :  950.0       eV\n charge density mixing g-vector                 :  1.500       1/A\n\n *********************** Population Analysis Parameters ************************\n  \n Population analysis with cutoff                :  3.000       A\n\n ************************** Band Structure Parameters **************************\n  \n max. number of iterations                      :         60\n max. CG steps in BS calc                       :         25\n number of bands / k-point                      :         30\n band convergence tolerance                     : 0.1000E-05   eV\n write orbitals file                            : on\n\n *******************************************************************************\n  \n \n                           -------------------------------\n                                      Unit Cell\n                           -------------------------------\n        Real Lattice(A)                      Reciprocal Lattice(1/A)\n   3.2094013   0.0000000   0.0000000        1.9577438   1.1303370   0.0000000\n  -1.6047006   2.7793415   0.0000000        0.0000000   2.2606741   0.0000000\n   0.0000000   0.0000000   5.2108021        0.0000000   0.0000000   1.2058000\n \n                       Lattice parameters(A)       Cell Angles\n                    a =    3.209401          alpha =   90.000000\n                    b =    3.209331          beta  =   90.000000\n                    c =    5.210802          gamma =  120.000728\n \n                       Current cell volume =   46.480470       A**3\n \n                           -------------------------------\n                                     Cell Contents\n                           -------------------------------\n                         Total number of ions in cell =    2\n                      Total number of species in cell =    1\n                        Max number of any one species =    2\n \n            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n            x  Element    Atom        Fractional coordinates of atoms  x\n            x            Number           u          v          w      x\n            x----------------------------------------------------------x\n            x  Mg           1         0.000000   0.000000   0.000000   x \n            x  Mg           2         0.333333   0.666667   0.500000   x \n            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n \n \n                         No user defined ionic velocities\n \n                           -------------------------------\n                                   Details of Species\n                           -------------------------------\n \n                               Mass of species in AMU\n                                    Mg   24.3050000\n \n                          Electric Quadrupole Moment (Barn)\n                                    Mg    0.1994000 Isotope 25\n \n                          Files used for pseudopotentials:\n                                    Mg Mg_00PBE_OP.recpot\n \n                           -------------------------------\n                              k-Points For BZ Sampling\n                           -------------------------------\n                       MP grid size for SCF calculation is  4  4  2\n                            with an offset of   0.000  0.000  0.000\n                       Number of kpoints used =             6\n \n                           -------------------------------\n                               Symmetry and Constraints\n                           -------------------------------\n \n                      Maximum deviation from symmetry = 0.320940E-07     ANG\n \n !-----------------------------------------------------------------------------!\n !                                                                             !\n !                     * * *      W A R N I N G     * * *                      !\n !                                                                             !\n !   Atomic positions are slightly inconsistent with symmetry operations,      !\n !   (although the discrepancy is lower than value of SYMMETRY_TOL).           !\n  !   In severe cases this could lead to seriously inaccurate results.         !\n !                                                                             !\n !   It is recommended that symmetries are satisfied to a tolerance of 1e-8au  !\n !   or better. Consider one or more of the following:                         !\n !     - Adjust input so that atomic positions and symmetry are consistent.    !\n !     - Use SNAP_TO_SYMMETRY keyword in the cell file to adjust co-ordinates. !\n !     - Reduce the value of SYMMETRY_TOL in the cell file.                    !\n !     - Remove incorrect symmetry operations from .cell file.                 !\n !                                                                             !\n !-----------------------------------------------------------------------------!\n \n \n                      Number of symmetry operations   =          24\n                      Number of ionic constraints     =           3\n                      Point group of crystal =    27: D6h, 6/mmm, 6/m 2/m 2/m\n                      Space group of crystal =   194: P6_3/mmc, -P 6c 2c\n \n             Set iprint > 1 for details on symmetry rotations/translations\n \n                         Centre of mass is constrained\n             Set iprint > 1 for details of linear ionic constraints\n \n                         Number of cell constraints= 3\n                         Cell constraints are:  1 2 3 0 0 0\n \n                         External pressure/stress (GPa)\n                          0.00000   0.00000   0.00000\n                                    0.00000   0.00000\n                                              0.00000\n  \n+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER PROCESS --------------+\n|                                                     Memory          Disk    |\n| Model and support data                               21.9 MB         0.0 MB |\n| Electronic energy minimisation requirements           5.3 MB         0.0 MB |\n|                                               ----------------------------- |\n| Approx. total storage required per process           27.2 MB         0.0 MB |\n|                                                                             |\n| Requirements will fluctuate during execution and may exceed these estimates |\n+-----------------------------------------------------------------------------+\nCalculating finite basis set correction with  3 cut-off energies.\nCalculating total energy with cut-off of  940.000eV.\n------------------------------------------------------------------------ <-- SCF\nSCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF\n                               energy          per atom          (sec)   <-- SCF\n------------------------------------------------------------------------ <-- SCF\nInitial  -3.23117569E+003  0.00000000E+000                         0.84  <-- SCF\n      1  -3.15232116E+003  7.50066098E+000  -3.94272621E+001       1.00  <-- SCF\n      2  -3.19180736E+003  4.30800366E+000   1.97430974E+001       1.09  <-- SCF\n      3  -3.19225744E+003  4.16183547E+000   2.25042827E-001       1.18  <-- SCF\n      4  -3.19224523E+003  4.12812261E+000  -6.10870907E-003       1.44  <-- SCF\n      5  -3.19223715E+003  4.10378219E+000  -4.03602880E-003       1.68  <-- SCF\n      6  -3.19223672E+003  4.09935587E+000  -2.15274160E-004       1.93  <-- SCF\n      7  -3.19223672E+003  4.09932287E+000  -1.44811178E-006       2.19  <-- SCF\n      8  -3.19223672E+003  4.09931145E+000  -8.31207028E-007       2.41  <-- SCF\n      9  -3.19223672E+003  4.09930820E+000  -2.87002676E-008       2.63  <-- SCF\n------------------------------------------------------------------------ <-- SCF\n \nFinal energy, E             =  -3192.239833908     eV\nFinal free energy (E-TS)    =  -3192.236718311     eV\n(energies not corrected for finite basis set)\n \nNB est. 0K energy (E-2TS/3)      =  -3192.237756843     eV\n \nCalculating total energy with cut-off of  945.000eV.\n------------------------------------------------------------------------ <-- SCF\nSCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF\n                               energy          per atom          (sec)   <-- SCF\n------------------------------------------------------------------------ <-- SCF\nInitial  -3.19223986E+003  0.00000000E+000                         2.91  <-- SCF\n      1  -3.19223679E+003  4.09930594E+000  -1.53341658E-003       3.01  <-- SCF\n      2  -3.19223679E+003  4.09930594E+000   2.66434261E-009       3.23  <-- SCF\n      3  -3.19223679E+003  4.09930194E+000   1.27145987E-009       3.44  <-- SCF\n      4  -3.19223679E+003  4.09930554E+000   2.57482224E-009       3.66  <-- SCF\n------------------------------------------------------------------------ <-- SCF\n \nFinal energy, E             =  -3192.239909208     eV\nFinal free energy (E-TS)    =  -3192.236793610     eV\n(energies not corrected for finite basis set)\n \nNB est. 0K energy (E-2TS/3)      =  -3192.237832143     eV\n \nCalculating total energy with cut-off of  950.000eV.\n------------------------------------------------------------------------ <-- SCF\nSCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF\n                               energy          per atom          (sec)   <-- SCF\n------------------------------------------------------------------------ <-- SCF\nInitial  -3.19223994E+003  0.00000000E+000                         3.93  <-- SCF\n      1  -3.19223687E+003  4.09929524E+000  -1.53254318E-003       4.03  <-- SCF\n      2  -3.19223687E+003  4.09929524E+000   3.13224603E-010       4.25  <-- SCF\n      3  -3.19223687E+003  4.09928758E+000  -1.83990452E-009       4.46  <-- SCF\n      4  -3.19223687E+003  4.09929560E+000   4.47253797E-009       4.68  <-- SCF\n------------------------------------------------------------------------ <-- SCF\n \nFinal energy, E             =  -3192.239986240     eV\nFinal free energy (E-TS)    =  -3192.236870642     eV\n(energies not corrected for finite basis set)\n \nNB est. 0K energy (E-2TS/3)      =  -3192.237909175     eV\n \n For future reference: finite basis dEtot/dlog(Ecut) =      -0.015141eV\n Total energy corrected for finite basis set =   -3192.236871 eV\n  \n+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER PROCESS --------------+\n|                                                     Memory          Disk    |\n| Model and support data                               71.8 MB         0.0 MB |\n| Band structure calculation requirements             239.8 MB         0.0 MB |\n|                                               ----------------------------- |\n| Approx. total storage required per process          311.6 MB         0.0 MB |\n|                                                                             |\n| Requirements will fluctuate during execution and may exceed these estimates |\n+-----------------------------------------------------------------------------+\n\nWriting model to /var/tmp/pbs.1525038.cx1/Mg_hcp_bands_10el_recpot.orbitals\n  + Calculation re-parallelised over 12 processes.                    +\n  + Data distributed by G-vector(2-way), k-point(6-way)               +\n \n ***************** Symmetrised Forces *****************\n *                                                    *\n *            Cartesian components (eV/A)             *\n * -------------------------------------------------- *\n *                   x            y            z      *\n *                                                    *\n * Mg        1      0.00000     -0.00000      0.00000 *\n * Mg        2     -0.00000      0.00000      0.00000 *\n *                                                    *\n ******************************************************\n \n *********** Symmetrised Stress Tensor ***********\n *                                               *\n *          Cartesian components (GPa)           *\n * --------------------------------------------- *\n *             x             y             z     *\n *                                               *\n *  x     -2.675983      0.000000      0.000000  *\n *  y      0.000000     -2.675826      0.000000  *\n *  z      0.000000      0.000000      2.558953  *\n *                                               *\n *  Pressure:    0.9310                          *\n *                                               *\n *************************************************\n \n Pseudo atomic calculation performed for Mg 2s2 2p6 3s2\n \n Converged in 23 iterations to a total energy of -1594.7271 eV\n \nCharge spilling parameter for spin component 1 = 0.06%\n \n     Atomic Populations (Mulliken)\n     -----------------------------\nSpecies   Ion     s      p      d      f     Total  Charge (e)\n==============================================================\n  Mg       1     2.77   7.23   0.00   0.00  10.00     0.00\n  Mg       2     2.77   7.23   0.00   0.00  10.00     0.00\n==============================================================\n \n            Bond              Population      Length (A)\n============================================================\n============================================================\n \n \n     Hirshfeld Analysis\n     ------------------\nSpecies   Ion     Hirshfeld Charge (e)  Spin (hbar/2)\n===================================================\n  Mg       1                 0.00        0.00\n  Mg       2                -0.00        0.00\n===================================================\n \n\nWriting analysis data to /var/tmp/pbs.1525038.cx1/Mg_hcp_bands_10el_recpot.castep_bin\n\nWriting model to /var/tmp/pbs.1525038.cx1/Mg_hcp_bands_10el_recpot.check\n \n A BibTeX formatted list of references used in this run has been written to \n /var/tmp/pbs.1525038.cx1/Mg_hcp_bands_10el_recpot.bib\n \nInitialisation time =      0.65 s\nCalculation time    =    679.56 s\nFinalisation time   =      0.30 s\nTotal time          =    680.51 s\nPeak Memory Use     = 455564 kB\n  \nOverall parallel efficiency rating: Poor (43%)                                  \n  \nData was distributed by:-\nG-vector (2-way); efficiency rating: Poor (43%)                                 \nk-point (6-way); efficiency rating: Ideal (100%)                                \n)\n\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath, fileStr

if __name__ == '__main__':
	unittest.main()


