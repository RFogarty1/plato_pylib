#!/usr/bin/python3


import os
import unittest

import numpy as np

import plato_pylib.plato.parse_plato_out_files as tCode

class testExtractEnergiesFromOutFile(unittest.TestCase):
	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["tb1fileA"], unusedStr = createTb1OutFilePartialA()
		self.filePathDict["tb2fileA"], unusedStr = createTb2OutFilePartialA()
	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

	def testCohesiveElectronicEnergy_tb1(self):
		''' Test the electronic energy from tb1 is parsed correctly; not this is a cohesive energy measure '''
		testInpFile = self.filePathDict["tb1fileA"]
		expectedElectronic = 1.0089780 + -1.8915282 
		actualElectronic = tCode.parsePlatoOutFile(testInpFile)["energies"].electronicCohesiveE
		self.assertAlmostEqual(expectedElectronic,actualElectronic)

	def testCohesiveElectronicEnergy_tb2(self):
		testInpFile = self.filePathDict["tb2fileA"]
#		expectedCohesive = -0.03136604
		expectedCohesive = -0.031328838
		actualCohesive = tCode.parsePlatoOutFile(testInpFile)["energies"].electronicCohesiveE
		self.assertAlmostEqual(expectedCohesive, actualCohesive)


	def testCohesiveFreeEnergy_tb2(self):
		testInpFile = self.filePathDict["tb2fileA"]
		expectedCohesive = -0.03136604
		actualCohesive = tCode.parsePlatoOutFile(testInpFile)["energies"].freeCohesiveE
		self.assertAlmostEqual(expectedCohesive, actualCohesive, places=5)

	def testE0CohesiveFreeEnergy_tb2(self):
		testInpFile = self.filePathDict["tb2fileA"]
		expectedE0Cohesive = -0.167288
		actualE0Cohesive = tCode.parsePlatoOutFile(testInpFile)["energies"].e0Coh
		print("actual cohesive = {}".format(actualE0Cohesive))
		self.assertAlmostEqual(expectedE0Cohesive, actualE0Cohesive, places=5)

	def testCohesiveFreeEnergy_tb1(self):
		testInpFile = self.filePathDict["tb1fileA"]
		expectedCohesive = 1.0089780 + -1.8915282 + 0.0010000 
		actualCohesive = tCode.parsePlatoOutFile(testInpFile)["energies"].freeCohesiveE
		self.assertAlmostEqual(expectedCohesive, actualCohesive)

class testExtractLatticeFromOutFile(unittest.TestCase):
	def setUp(self):
		self.filePathDict = dict()
		self.filePathDict["tb1FileOptA"], unusedStr = createTb1OutFilePartialB_lattparamsOpt()
		self.filePathDict["tb2FileCellA"], unusedStr = createTb2OutFilePartialB_lattparams()
	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

	def testLattParamsAngles_opt(self):
		''' Test tCode.parsePlatoOutFile produces expected lattice parameters/angles for
		    an optimisation (tb1) '''
		testInpFile = self.filePathDict["tb1FileOptA"]
		expectedLattParams = [7.45026, 7.45025616801194, 12.09628]

		expectedLattAngles = [90.0, 90.0, 120.0]
		expectedVolume = 581.466851 #Comes from plato output file
		allFileInfo = tCode.parsePlatoOutFile(testInpFile)

		actualLattParams = allFileInfo["unitCell"].getLattParamsList()
		actualLattAngles = allFileInfo["unitCell"].getLattAnglesList()
		actualVolume = allFileInfo["unitCell"].getVolume()
		[self.assertAlmostEqual(x,y, places=4) for x,y in zip(expectedLattParams, actualLattParams)]
		[self.assertAlmostEqual(x,y, places=4) for x,y in zip(expectedLattAngles, actualLattAngles)]
		self.assertAlmostEqual(expectedVolume, actualVolume, places=3)

	def testLattParamsAngles_tb2(self):
		''' Test tCode.parsePlatoOutFile gets correct lattice params/angles from tb2 '''
		testInpFile = self.filePathDict["tb2FileCellA"]
		expectedLattParams = [6.064889, 6.064889, 9.846989]
		expectedLattAngles = [90.0,90.0,120.0]
		allFileInfo = tCode.parsePlatoOutFile(testInpFile)

		actualLattParams = allFileInfo["unitCell"].getLattParamsList()
		actualLattAngles = allFileInfo["unitCell"].getLattAnglesList()
		[self.assertAlmostEqual(x,y, places=4) for x,y in zip(expectedLattParams, actualLattParams)]
		[self.assertAlmostEqual(x,y, places=4) for x,y in zip(expectedLattAngles, actualLattAngles)]


class testParseDftFile(unittest.TestCase):

	def setUp(self):
		self.filePaths = dict()
		self.filePaths["dftA"] = createDftOutFileA()

	def tearDown(self):
		pass

	def testParseEnergy(self):
		expectedEnergy = -1.7340183
		actualEnergy = tCode.parsePlatoOutFile(self.filePaths["dftA"])["energies"].electronicTotalE
		self.assertAlmostEqual(expectedEnergy, actualEnergy)

	def testDoesntReturnCohesive(self):
		with self.assertRaises(ValueError):
			tCode.parsePlatoOutFile(self.filePaths["dftA"])["energies"].electronicCohesiveE

	def testParseNumbAtoms(self):
		expectedNumbAtoms = 1
		actualNumbAtoms = tCode.parsePlatoOutFile(self.filePaths["dftA"])["numbAtoms"]
		self.assertEqual(expectedNumbAtoms, actualNumbAtoms)


class testParseOccFile(unittest.TestCase):

	def setUp(self):
		self.pbcsFile = createPartialOccFileWithPBCs()
		self.noPbcsFile = createPartialOccFileNoPBCs()

	def tearDown(self):
		os.remove(self.pbcsFile)
		os.remove(self.noPbcsFile)

	def testParseOccNoPBCs(self):
		''' parseOccFile should read values correctly when no periodic boundary conditions are applied '''
		expDict = dict()
		expDict["k_path"] = np.array([0.0, 0.0, 0.0])
		expDict["eigen_vals"] = np.array([1.88352, 4.92246, 4.94033, 4.9407, 8.68629, 8.68861, 10.4299, 32.6306])
		expDict["occs"] = np.array(( [2.0000000000, 1.0677221221, 0.4710287010, 0.4612491767, 0.0, 0.0, 0.0, 0.0] ))
		expDict["k_weights"] = np.array(( [1] ))

		parsedOccDict = tCode.parseOccFile(self.noPbcsFile)

		for key in expDict.keys():
			self.assertTrue( np.allclose(expDict[key], parsedOccDict[key]) )

	def testParseOccPBCs(self):
		''' parseOccFile should read values correctly when periodic boundary conditions are applied '''
		expDict = dict()
		expDict["k_path"] = np.array(([0.0,0.0,0.0],[0.0,0.0,0.00820]))
		expDict["eigen_vals"] = np.array(([-5.86279, -1.08205, -0.895486, 8.63832, 8.63832, 10.4732, 10.4732, 10.611],
		                             [-5.86251, -1.11491, -0.861833, 8.63859, 8.63859, 10.4729, 10.4729, 10.6091]))
		expDict["occs"] = np.array(([1.0,1.0,1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
		                         [1.0,1.0,1.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
		expDict["k_weights"] = np.array([0.0016000000,0.0016000000])

		parsedOccDict = tCode.parseOccFile(self.pbcsFile)

		for key in expDict.keys():
			self.assertTrue( np.allclose(expDict[key], parsedOccDict[key]) )

class testEnergyValsClass(unittest.TestCase):
	def testTotalElectronicTb1(self):
		testObj = tCode.EnergyVals(e0coh=2.9,e1=1.1)
		expectedElectronic = 4.0
		actualElectronic = testObj.electronicCohesiveE
		self.assertAlmostEqual(expectedElectronic, actualElectronic)

class testParseDft2FileWithScfInfo(unittest.TestCase):

	def setUp(self):
		self.testPathNoInfo, unused = createTb2OutFilePartialA()
		self.testPathNotConverged, unused = createDft2OutFileWithScfNotConvergedA() #NOTE: This actualy WRITES the file out
		self.testPathConverged, unused = createDft2OutFileWithScfConvergedA()

	def tearDown(self):
		os.remove(self.testPathNoInfo)
		os.remove(self.testPathNotConverged)
		os.remove(self.testPathConverged)

	def testCorrectlyFiguresOutNotConveged(self):
		parsedFile = tCode.parsePlatoOutFile(self.testPathNotConverged)
		expConvVal = False
		actConvVal = parsedFile["scf_is_converged"]
		self.assertEqual(expConvVal,actConvVal)

	def testWeGetNoneForAnOldStyleFile(self):
		parsedFile = tCode.parsePlatoOutFile(self.testPathNoInfo)
		expConvVal = None
		actConvVal = parsedFile["scf_is_converged"]
		self.assertEqual(expConvVal,actConvVal)

	def testCorrectlyFiguresOutConvergence(self):
		parsedFile = tCode.parsePlatoOutFile(self.testPathConverged)
		expConvVal = True
		actConvVal = parsedFile["scf_is_converged"]
		self.assertEqual(expConvVal,actConvVal)


def createPartialOccFileNoPBCs():
	filePath = os.path.join( os.getcwd(), "partialOccNoPBC.occ" )
	with open(filePath,"wt") as f:
		f.write("    0     8              4  0\n      1.88352  2.0000000000\n      4.92246  1.0677221221\n      4.94033  0.4710287010\n       4.9407  0.4612491767\n      8.68629  0.0000000000\n      8.68861  0.0000000000\n      10.4299  0.0000000000\n      32.6306  0.0000000000\n")
	return filePath

def createPartialOccFileWithPBCs():
	filePath = os.path.join( os.getcwd(), "partialOccWithPBCs.occ")
	with open(filePath,"wt") as f:
		f.write("  2     8              4  0\nK-point 1   0.00000   0.00000   0.00000 0.0016000000\n     -5.86279  1.0000000000\n     -1.08205  1.0000000000\n    -0.895486  1.0000000000\n      8.63832  0.0000000000\n      8.63832  0.0000000000\n      10.4732  0.0000000000\n      10.4732  0.0000000000\n       10.611  0.0000000000\nK-point 2   0.00000   0.00000   0.00820 0.0016000000\n     -5.86251  1.0000000000\n     -1.11491  1.0000000000\n    -0.861833  1.0000000000\n      8.63859  0.0000000000\n      8.63859  0.0000000000\n      10.4729  0.0000000000\n      10.4729  0.0000000000\n      10.6091  0.0000000000\n")
	return filePath



#Note i modified the entropy (from exactly 0) without modifying other bits
def createTb1OutFilePartialA():
	filePath = os.path.join( os.getcwd(), "tb1outPartialA.out")
	fileStr ="Model characteristics\n---------------------\nTabulated tight binding\nOverlap      : 0\nCrystal field: 0\nThree centre : 0\nN_{ij} Integrals : 0\nSN Integrals: 0\nSelf-consistency parameters:\n----------------------------\nNumber of loops         = 30\nEnergy tolerance        = 0.000000\nDensity tolerance       = 0.001000\nMixing scheme           = 1\nPreconditioning type    = 1\nThomas-Fermi K          = 1.000000\nMix factor              = 0.200000\nNumber of mixing levels = 5\n\nOther parameters:\n-----------------\nExcess number of electrons = 0.000000\n\nAfter relaxation:\n=================\n\nPrimitive translation vectors (a0):\n-----------------------------------\n     5.13147      5.13147      0.00000\n     5.13147      0.00000      5.13147\n     0.00000      5.13147      5.13147\nVolume       = 270.244193 a0^3\nCell repeats =    4    3    3\n\nNumber of atoms = 2\nNeighbor cutoff =      7.00000 a0\n\nAtomic positions (a0):\n----------------------\nSi     10.25745     10.25745     10.25745\nSi      2.57124      2.57124      2.57124\n\nTotal forces (Ry/a0):\n---------------------\n -1.417724538e-05  -1.417724538e-05  -1.417724538e-05\n  1.417724538e-05   1.417724538e-05   1.417724538e-05\n\nStress Tensor (MBars):\n----------------------\n   0.004974960316    0.009437011637    0.009437011637\n   0.009437011637    0.004974960316    0.009437011637\n   0.009437011637    0.009437011637    0.004974960316\n\nEFermi =     -0.4109111 Ry\n(kT    =      0.0010000 Ry)\n\nE0               =      1.0089780 Ry\nE1               =     -1.8915282 Ry\nEntropy          =      0.0010000 Ry\n------------------------------------\nTotal (Free)     =     -0.8825502 Ry\n                 =    -12.0078015 eV\n                 =   -277.6061662 kCal/mol\n                 =  -1161.4360667 kJ/mol\n------------------------------------\nTotal (Internal) =     -0.8825502 Ry\n                 =    -12.0078015 eV\n                 =   -277.6061662 kCal/mol\n                 =  -1161.4360667 kJ/mol\n------------------------------------\nTotal (kT->0)    =     -0.8825502 Ry\n                 =    -12.0078015 eV\n                 =   -277.6061662 kCal/mol\n                 =  -1161.4360667 kJ/mol\n\nNet atomic charges (Mulliken):\n------------------------------\nSi(   1) =     -0.00000\nSi(   2) =     -0.00000\n\nOrbital populations (Mulliken):\n-------------------------------\nSi(   1):  1.47  0.84  0.84  0.84 \nSi(   2):  1.47  0.84  0.84  0.84 \n\nTotal multipole moments:\n------------------------\n\nN.B. See the manual for the ordering of the moments.\n\nDipole  [e*a0]:  -0.00000000  -0.00000000  -0.00000000\nDipole [Debye]:  -0.00000000  -0.00000000  -0.00000000\n\nQuadrupole  [e*a0^2]:  -0.00000000  -0.00000000  -0.00000000  -0.00000000  -0.00000000  -0.00000000\nQuadrupole [Debye*A]:  -0.00000000  -0.00000000  -0.00000000  -0.00000000  -0.00000000  -0.00000000\n\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr

def createTb2OutFilePartialA():
	filePath = os.path.join( os.getcwd(), "tb2outPartialA.out" )
	fileStr = "Cell vectors in Bohr radii\n==========================\na1 = (      3.372899,       3.372899,      -3.372899)\na2 = (     -3.372899,       3.372899,       3.372899)\na3 = (      3.372899,      -3.372899,       3.372899)\n\nb1 = (      0.931422,       0.931422,              0)\nb2 = (             0,       0.931422,       0.931422)\nb3 = (      0.931422,              0,       0.931422)\n\nGeometry (xyz format in Angstroms)\n==================================\n1\n\nMg              0              0              0\n\nExchange and correlation functional\n===================================\nPBE\n\nIntegral mesh\n=============\nType: regular mesh\ndx  : 0.180000\n\nElectron chemical potential (Ry)\n================================\nMu =     -0.1279633\n\nEnergy\n======\nBreakdown of zeroth order energy:\nKinetic energy           =      0.7101726 Ry\nNeutral atom energy      =     -0.9010047 Ry\nElectrostatic energy     =      -1.231837 Ry\nNon-local energy         =      0.5164655 Ry\nExchange and correlation =     -0.8636517 Ry\n--------------------------------------------\nZeroth order energy      =      -1.769856 Ry\n\nBreakdown of total energy:\nZeroth order energy      =      -1.769856 Ry\nFirst order energy       =      0.1359589 Ry\nElectron entropy         =   3.720202e-05 Ry\n--------------------------------------------\nTotal energy             =      -1.633934 Ry\n\nBreakdown of cohesive energy:\nTotal energy             =      -1.633934 Ry\nAtom energy              =      -1.602568 Ry\n--------------------------------------------\nCohesive energy          =    -0.03136604 Ry\n\nCharge Multipoles\n=================\n             s  \n[   0]  -0.00000\n\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr


def createTb1OutFilePartialB_lattparamsOpt():
	filePath = os.path.join( os.getcwd(), "tb1PartialBLattParams.out" )
	fileStr = "Electrostatic potential:\n------------------------\nMg(   1):              0\nMg(   2):              0\n\nDebugging Information:\n======================\nmmicom: 2\nTestStress: -0.0000043197\nUserExit: 0\nNSteps: 23\n\n\n\nAfter relaxation:\n=================\n\nPrimitive translation vectors (a0):\n-----------------------------------\n     7.45026      0.00000      0.00000\n    -3.72513      6.45211      0.00000\n     0.00000      0.00000     12.09628\nVolume       = 581.466851 a0^3\nCell repeats =    5    4    3\n\nNumber of atoms = 2\nNeighbor cutoff =     15.00000 a0\n\nAtomic positions (a0):\n----------------------\nMg      0.00000      0.00000      0.00000\nMg      3.20399      4.76903      6.04814\n\nTotal forces (Ry/a0):\n---------------------\n                0                 0                 0\n                0                 0                 0\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr

def createTb2OutFilePartialB_lattparams():
	filePath = os.path.join( os.getcwd(), "tb2PartialBLattParams.out" )
	fileStr = "# l =  1 ng =   1 np =   1\n#               a                           c\n     0.99985718149445101943       1.8374666607438268073\n\nCell vectors in Bohr radii\n==========================\na1 = (      6.064889,              0,              0)\na2 = (     -3.032445,       5.252348,              0)\na3 = (             0,              0,       9.846989)\n\nb1 = (      1.035993,      0.5981311,             -0)\nb2 = (             0,       1.196262,              0)\nb3 = (             0,             -0,      0.6380819)\n\nGeometry (xyz format in Angstroms)\n==================================\n2\n\nMg              0              0              0\nMg  -1.955318e-08       1.852948         2.6054\n\nExchange and correlation functional\n===================================\nPBE\n\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr


def createDft2OutFileWithScfNotConvergedA():
	filePath = os.path.join( os.getcwd(), "Dft2OutFileWithScfNotConvergedA.out" )
	fileStr = """
    -----------------------------
    | Density Functional Theory |
    -----------------------------


Model characteristics
---------------------
Tabulated tight binding
Overlap                 : 1
Crystal field           : 1
Crystal field correction: 0
Three centre            : 0

Self Consistency Cycles
===================
SCF: [  0] E =      -22.10673 Ry  Residue =       1.6
SCF: [  1] E =      -24.74322 Ry  Residue =      0.45

SCF Converged: False (Self-Consistency used)

Geometry (xyz format in Angstroms)
==================================
4

Zr              0              0              0
Zr              0              0       1.322942
Zr              0              0       2.645885
Zr              0              0       3.968828

Exchange and correlation functional
===================================
PBE
Exact functional

Integral mesh
=============
Type                   : atom centered spherical mesh
Number of r points     : 50
Number of theta points : 20
Number of phi points   : 20

Electron chemical potential (Ry)
================================
Mu =      -0.802983

Energy
======

Contribution from each order to the total energy:
Zeroth order energy      =      -13.39443 Ry
First order energy       =      -11.63092 Ry
Electron entropy         =              0 Ry
Second order energy      =      0.2821341 Ry
--------------------------------------------
Total energy             =      -24.74322 Ry

Calculation of cohesive energy:
Total energy             =      -24.74322 Ry
Atom energy              =      -21.99314 Ry
--------------------------------------------
Cohesive energy          =      -2.750075 Ry

Charge Multipoles
=================
             s  
[   0]  -0.90608
[   1]   0.61294
[   2]   0.10020
[   3]   0.19294

Time distribution (s):
======================
Total                       0.87
Build Hamiltonian           0.74
Diagonalisation             0.00
Build density matrix        0.00
Numerical integrals         0.60
Gaussian integrals          0.25
"""
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr

def createDft2OutFileWithScfConvergedA():
	filePath = os.path.join( os.getcwd(), "Dft2OutFileWithScfConvergedA.out" )
	fileStr = """
    -----------------------------
    | Density Functional Theory |
    -----------------------------


Model characteristics
---------------------
Tabulated tight binding
Overlap                 : 1
Crystal field           : 1
Crystal field correction: 0
Three centre            : 0

Self Consistency Cycles
===================
SCF: [  0] E =      -22.10673 Ry  Residue =       1.6
SCF: [  1] E =      -24.74322 Ry  Residue =      0.45
SCF: [  2] E =      -21.64941 Ry  Residue =      0.24
SCF: [  3] E =      -24.21376 Ry  Residue =   4.3e-10

SCF Converged: True (Self-Consistency used)

Geometry (xyz format in Angstroms)
==================================
4

Zr              0              0              0
Zr              0              0       1.322942
Zr              0              0       2.645885
Zr              0              0       3.968828

Exchange and correlation functional
===================================
PBE
Exact functional

Integral mesh
=============
Type                   : atom centered spherical mesh
Number of r points     : 50
Number of theta points : 20
Number of phi points   : 20

Electron chemical potential (Ry)
================================
Mu =     -0.7755432

Energy
======

Contribution from each order to the total energy:
Zeroth order energy      =      -13.39443 Ry
First order energy       =      -10.91655 Ry
Electron entropy         =              0 Ry
Second order energy      =     0.09722013 Ry
--------------------------------------------
Total energy             =      -24.21376 Ry

Calculation of cohesive energy:
Total energy             =      -24.21376 Ry
Atom energy              =      -21.99314 Ry
--------------------------------------------
Cohesive energy          =      -2.220618 Ry

Charge Multipoles
=================
             s  
[   0]  -0.37506
[   1]   0.41837
[   2]   0.16522
[   3]  -0.20853

Time distribution (s):
======================
Total                       0.95
Build Hamiltonian           0.72
Diagonalisation             0.00
Build density matrix        0.00
Numerical integrals         0.66
Gaussian integrals          0.28
"""
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath, fileStr


def createDftOutFileA():
	filePath = os.path.join( os.getcwd(), "dftOutFile.out" )
	fileStr = "    ----------------------------------\n    | LCAO density functional theory |\n    ----------------------------------\n\nSelf-consistency parameters:\n----------------------------\nNumber of loops         = 30\nEnergy tolerance        = 0.000000\nDensity tolerance       = 0.001000\nMixing scheme           = 1\nParallelise FFT flag    = 0\nSplit density flag      = 0\nPreconditioning type    = 1\nThomas-Fermi K          = 1.000000\nMix metric type         = 0\nMix factor              = 0.200000\nNumber of mixing levels = 5\n\nEnergy functional\n-----------------\nGeneralized gradient approximation\nNo spin polarization\n\nOther parameters:\n-----------------\nDiameter of non-local pot. = 0.000000\nBasis flag                 = 0\nExcess number of electrons = 0.000000\nExternal field type        = 0\n\nK points:\n---------\n-0.4500 -0.4500 -0.4500             0.002\n 0.0500 -0.4500 -0.4500             0.002\n 0.0500 -0.4500 -0.3500             0.002\n 0.0500 -0.4500 -0.2500             0.002\n 0.0500 -0.4500 -0.1500             0.002\n 0.0500 -0.4500 -0.0500             0.002\n 0.0500 -0.4500  0.0500             0.002\n 0.0500 -0.4500  0.1500             0.002\n 0.0500 -0.4500  0.2500             0.002\n 0.0500 -0.4500  0.3500             0.002\n 0.0500 -0.4500  0.4500             0.002\n 0.0500 -0.3500 -0.4500             0.002\n 0.0500 -0.3500 -0.3500             0.002\n 0.0500 -0.3500 -0.2500             0.002\n 0.0500 -0.3500 -0.1500             0.002\n 0.0500 -0.3500 -0.0500             0.002\n 0.0500 -0.3500  0.0500             0.002\n 0.0500 -0.3500  0.1500             0.002\n 0.0500 -0.3500  0.2500             0.002\n 0.0500 -0.3500  0.3500             0.002\n 0.0500 -0.3500  0.4500             0.002\n 0.0500 -0.2500 -0.4500             0.002\n 0.0500 -0.2500 -0.3500             0.002\n 0.0500 -0.2500 -0.2500             0.002\n 0.0500 -0.2500 -0.1500             0.002\n 0.0500 -0.2500 -0.0500             0.002\n 0.0500 -0.2500  0.0500             0.002\n 0.0500 -0.2500  0.1500             0.002\n 0.0500 -0.2500  0.2500             0.002\n 0.0500 -0.2500  0.3500             0.002\n 0.0500 -0.2500  0.4500             0.002\n 0.0500 -0.1500 -0.4500             0.002\n 0.0500 -0.1500 -0.3500             0.002\n 0.0500 -0.1500 -0.2500             0.002\n 0.0500 -0.1500 -0.1500             0.002\n 0.0500 -0.1500 -0.0500             0.002\n 0.0500 -0.1500  0.0500             0.002\n 0.0500 -0.1500  0.1500             0.002\n 0.0500 -0.1500  0.2500             0.002\n 0.0500 -0.1500  0.3500             0.002\n 0.0500 -0.1500  0.4500             0.002\n 0.0500 -0.0500 -0.4500             0.002\n 0.0500 -0.0500 -0.3500             0.002\n 0.0500 -0.0500 -0.2500             0.002\n 0.0500 -0.0500 -0.1500             0.002\n 0.0500 -0.0500 -0.0500             0.002\n 0.0500 -0.0500  0.0500             0.002\n 0.0500 -0.0500  0.1500             0.002\n 0.0500 -0.0500  0.2500             0.002\n 0.0500 -0.0500  0.3500             0.002\n 0.0500 -0.0500  0.4500             0.002\n 0.0500  0.0500 -0.4500             0.002\n 0.0500  0.0500 -0.3500             0.002\n 0.0500  0.0500 -0.2500             0.002\n 0.0500  0.0500 -0.1500             0.002\n 0.0500  0.0500 -0.0500             0.002\n 0.0500  0.0500  0.0500             0.002\n 0.0500  0.0500  0.1500             0.002\n 0.0500  0.0500  0.2500             0.002\n 0.0500  0.0500  0.3500             0.002\n 0.0500  0.0500  0.4500             0.002\n 0.0500  0.1500 -0.4500             0.002\n 0.0500  0.1500 -0.3500             0.002\n 0.0500  0.1500 -0.2500             0.002\n 0.0500  0.1500 -0.1500             0.002\n 0.0500  0.1500 -0.0500             0.002\n 0.0500  0.1500  0.0500             0.002\n 0.0500  0.1500  0.1500             0.002\n 0.0500  0.1500  0.2500             0.002\n 0.0500  0.1500  0.3500             0.002\n 0.0500  0.1500  0.4500             0.002\n 0.0500  0.2500 -0.4500             0.002\n 0.0500  0.2500 -0.3500             0.002\n 0.0500  0.2500 -0.2500             0.002\n 0.0500  0.2500 -0.1500             0.002\n 0.0500  0.2500 -0.0500             0.002\n 0.0500  0.2500  0.0500             0.002\n 0.0500  0.2500  0.1500             0.002\n 0.0500  0.2500  0.2500             0.002\n 0.0500  0.2500  0.3500             0.002\n 0.0500  0.2500  0.4500             0.002\n 0.0500  0.3500 -0.4500             0.002\n 0.0500  0.3500 -0.3500             0.002\n 0.0500  0.3500 -0.2500             0.002\n 0.0500  0.3500 -0.1500             0.002\n 0.0500  0.3500 -0.0500             0.002\n 0.0500  0.3500  0.0500             0.002\n 0.0500  0.3500  0.1500             0.002\n 0.0500  0.3500  0.2500             0.002\n 0.0500  0.3500  0.3500             0.002\n 0.0500  0.3500  0.4500             0.002\n 0.0500  0.4500 -0.4500             0.002\n 0.0500  0.4500 -0.3500             0.002\n 0.0500  0.4500 -0.2500             0.002\n 0.0500  0.4500 -0.1500             0.002\n 0.0500  0.4500 -0.0500             0.002\n 0.0500  0.4500  0.0500             0.002\n 0.0500  0.4500  0.1500             0.002\n 0.0500  0.4500  0.2500             0.002\n 0.0500  0.4500  0.3500             0.002\n 0.0500  0.4500  0.4500             0.002\n 0.1500 -0.4500 -0.4500             0.002\n 0.1500 -0.4500 -0.3500             0.002\n 0.1500 -0.4500 -0.2500             0.002\n 0.1500 -0.4500 -0.1500             0.002\n 0.1500 -0.4500 -0.0500             0.002\n 0.1500 -0.4500  0.0500             0.002\n 0.1500 -0.4500  0.1500             0.002\n 0.1500 -0.4500  0.2500             0.002\n 0.1500 -0.4500  0.3500             0.002\n 0.1500 -0.4500  0.4500             0.002\n 0.1500 -0.3500 -0.4500             0.002\n 0.1500 -0.3500 -0.3500             0.002\n 0.1500 -0.3500 -0.2500             0.002\n 0.1500 -0.3500 -0.1500             0.002\n 0.1500 -0.3500 -0.0500             0.002\n 0.1500 -0.3500  0.0500             0.002\n 0.1500 -0.3500  0.1500             0.002\n 0.1500 -0.3500  0.2500             0.002\n 0.1500 -0.3500  0.3500             0.002\n 0.1500 -0.3500  0.4500             0.002\n 0.1500 -0.2500 -0.4500             0.002\n 0.1500 -0.2500 -0.3500             0.002\n 0.1500 -0.2500 -0.2500             0.002\n 0.1500 -0.2500 -0.1500             0.002\n 0.1500 -0.2500 -0.0500             0.002\n 0.1500 -0.2500  0.0500             0.002\n 0.1500 -0.2500  0.1500             0.002\n 0.1500 -0.2500  0.2500             0.002\n 0.1500 -0.2500  0.3500             0.002\n 0.1500 -0.2500  0.4500             0.002\n 0.1500 -0.1500 -0.4500             0.002\n 0.1500 -0.1500 -0.3500             0.002\n 0.1500 -0.1500 -0.2500             0.002\n 0.1500 -0.1500 -0.1500             0.002\n 0.1500 -0.1500 -0.0500             0.002\n 0.1500 -0.1500  0.0500             0.002\n 0.1500 -0.1500  0.1500             0.002\n 0.1500 -0.1500  0.2500             0.002\n 0.1500 -0.1500  0.3500             0.002\n 0.1500 -0.1500  0.4500             0.002\n 0.1500 -0.0500 -0.4500             0.002\n 0.1500 -0.0500 -0.3500             0.002\n 0.1500 -0.0500 -0.2500             0.002\n 0.1500 -0.0500 -0.1500             0.002\n 0.1500 -0.0500 -0.0500             0.002\n 0.1500 -0.0500  0.0500             0.002\n 0.1500 -0.0500  0.1500             0.002\n 0.1500 -0.0500  0.2500             0.002\n 0.1500 -0.0500  0.3500             0.002\n 0.1500 -0.0500  0.4500             0.002\n 0.1500  0.0500 -0.4500             0.002\n 0.1500  0.0500 -0.3500             0.002\n 0.1500  0.0500 -0.2500             0.002\n 0.1500  0.0500 -0.1500             0.002\n 0.1500  0.0500 -0.0500             0.002\n 0.1500  0.0500  0.0500             0.002\n 0.1500  0.0500  0.1500             0.002\n 0.1500  0.0500  0.2500             0.002\n 0.1500  0.0500  0.3500             0.002\n 0.1500  0.0500  0.4500             0.002\n 0.1500  0.1500 -0.4500             0.002\n 0.1500  0.1500 -0.3500             0.002\n 0.1500  0.1500 -0.2500             0.002\n 0.1500  0.1500 -0.1500             0.002\n 0.1500  0.1500 -0.0500             0.002\n 0.1500  0.1500  0.0500             0.002\n 0.1500  0.1500  0.1500             0.002\n 0.1500  0.1500  0.2500             0.002\n 0.1500  0.1500  0.3500             0.002\n 0.1500  0.1500  0.4500             0.002\n 0.1500  0.2500 -0.4500             0.002\n 0.1500  0.2500 -0.3500             0.002\n 0.1500  0.2500 -0.2500             0.002\n 0.1500  0.2500 -0.1500             0.002\n 0.1500  0.2500 -0.0500             0.002\n 0.1500  0.2500  0.0500             0.002\n 0.1500  0.2500  0.1500             0.002\n 0.1500  0.2500  0.2500             0.002\n 0.1500  0.2500  0.3500             0.002\n 0.1500  0.2500  0.4500             0.002\n 0.1500  0.3500 -0.4500             0.002\n 0.1500  0.3500 -0.3500             0.002\n 0.1500  0.3500 -0.2500             0.002\n 0.1500  0.3500 -0.1500             0.002\n 0.1500  0.3500 -0.0500             0.002\n 0.1500  0.3500  0.0500             0.002\n 0.1500  0.3500  0.1500             0.002\n 0.1500  0.3500  0.2500             0.002\n 0.1500  0.3500  0.3500             0.002\n 0.1500  0.3500  0.4500             0.002\n 0.1500  0.4500 -0.4500             0.002\n 0.1500  0.4500 -0.3500             0.002\n 0.1500  0.4500 -0.2500             0.002\n 0.1500  0.4500 -0.1500             0.002\n 0.1500  0.4500 -0.0500             0.002\n 0.1500  0.4500  0.0500             0.002\n 0.1500  0.4500  0.1500             0.002\n 0.1500  0.4500  0.2500             0.002\n 0.1500  0.4500  0.3500             0.002\n 0.1500  0.4500  0.4500             0.002\n 0.2500 -0.4500 -0.4500             0.002\n 0.2500 -0.4500 -0.3500             0.002\n 0.2500 -0.4500 -0.2500             0.002\n 0.2500 -0.4500 -0.1500             0.002\n 0.2500 -0.4500 -0.0500             0.002\n 0.2500 -0.4500  0.0500             0.002\n 0.2500 -0.4500  0.1500             0.002\n 0.2500 -0.4500  0.2500             0.002\n 0.2500 -0.4500  0.3500             0.002\n 0.2500 -0.4500  0.4500             0.002\n 0.2500 -0.3500 -0.4500             0.002\n 0.2500 -0.3500 -0.3500             0.002\n 0.2500 -0.3500 -0.2500             0.002\n 0.2500 -0.3500 -0.1500             0.002\n 0.2500 -0.3500 -0.0500             0.002\n 0.2500 -0.3500  0.0500             0.002\n 0.2500 -0.3500  0.1500             0.002\n 0.2500 -0.3500  0.2500             0.002\n 0.2500 -0.3500  0.3500             0.002\n 0.2500 -0.3500  0.4500             0.002\n 0.2500 -0.2500 -0.4500             0.002\n 0.2500 -0.2500 -0.3500             0.002\n 0.2500 -0.2500 -0.2500             0.002\n 0.2500 -0.2500 -0.1500             0.002\n 0.2500 -0.2500 -0.0500             0.002\n 0.2500 -0.2500  0.0500             0.002\n 0.2500 -0.2500  0.1500             0.002\n 0.2500 -0.2500  0.2500             0.002\n 0.2500 -0.2500  0.3500             0.002\n 0.2500 -0.2500  0.4500             0.002\n 0.2500 -0.1500 -0.4500             0.002\n 0.2500 -0.1500 -0.3500             0.002\n 0.2500 -0.1500 -0.2500             0.002\n 0.2500 -0.1500 -0.1500             0.002\n 0.2500 -0.1500 -0.0500             0.002\n 0.2500 -0.1500  0.0500             0.002\n 0.2500 -0.1500  0.1500             0.002\n 0.2500 -0.1500  0.2500             0.002\n 0.2500 -0.1500  0.3500             0.002\n 0.2500 -0.1500  0.4500             0.002\n 0.2500 -0.0500 -0.4500             0.002\n 0.2500 -0.0500 -0.3500             0.002\n 0.2500 -0.0500 -0.2500             0.002\n 0.2500 -0.0500 -0.1500             0.002\n 0.2500 -0.0500 -0.0500             0.002\n 0.2500 -0.0500  0.0500             0.002\n 0.2500 -0.0500  0.1500             0.002\n 0.2500 -0.0500  0.2500             0.002\n 0.2500 -0.0500  0.3500             0.002\n 0.2500 -0.0500  0.4500             0.002\n 0.2500  0.0500 -0.4500             0.002\n 0.2500  0.0500 -0.3500             0.002\n 0.2500  0.0500 -0.2500             0.002\n 0.2500  0.0500 -0.1500             0.002\n 0.2500  0.0500 -0.0500             0.002\n 0.2500  0.0500  0.0500             0.002\n 0.2500  0.0500  0.1500             0.002\n 0.2500  0.0500  0.2500             0.002\n 0.2500  0.0500  0.3500             0.002\n 0.2500  0.0500  0.4500             0.002\n 0.2500  0.1500 -0.4500             0.002\n 0.2500  0.1500 -0.3500             0.002\n 0.2500  0.1500 -0.2500             0.002\n 0.2500  0.1500 -0.1500             0.002\n 0.2500  0.1500 -0.0500             0.002\n 0.2500  0.1500  0.0500             0.002\n 0.2500  0.1500  0.1500             0.002\n 0.2500  0.1500  0.2500             0.002\n 0.2500  0.1500  0.3500             0.002\n 0.2500  0.1500  0.4500             0.002\n 0.2500  0.2500 -0.4500             0.002\n 0.2500  0.2500 -0.3500             0.002\n 0.2500  0.2500 -0.2500             0.002\n 0.2500  0.2500 -0.1500             0.002\n 0.2500  0.2500 -0.0500             0.002\n 0.2500  0.2500  0.0500             0.002\n 0.2500  0.2500  0.1500             0.002\n 0.2500  0.2500  0.2500             0.002\n 0.2500  0.2500  0.3500             0.002\n 0.2500  0.2500  0.4500             0.002\n 0.2500  0.3500 -0.4500             0.002\n 0.2500  0.3500 -0.3500             0.002\n 0.2500  0.3500 -0.2500             0.002\n 0.2500  0.3500 -0.1500             0.002\n 0.2500  0.3500 -0.0500             0.002\n 0.2500  0.3500  0.0500             0.002\n 0.2500  0.3500  0.1500             0.002\n 0.2500  0.3500  0.2500             0.002\n 0.2500  0.3500  0.3500             0.002\n 0.2500  0.3500  0.4500             0.002\n 0.2500  0.4500 -0.4500             0.002\n 0.2500  0.4500 -0.3500             0.002\n 0.2500  0.4500 -0.2500             0.002\n 0.2500  0.4500 -0.1500             0.002\n 0.2500  0.4500 -0.0500             0.002\n 0.2500  0.4500  0.0500             0.002\n 0.2500  0.4500  0.1500             0.002\n 0.2500  0.4500  0.2500             0.002\n 0.2500  0.4500  0.3500             0.002\n 0.2500  0.4500  0.4500             0.002\n 0.3500 -0.4500 -0.4500             0.002\n 0.3500 -0.4500 -0.3500             0.002\n 0.3500 -0.4500 -0.2500             0.002\n 0.3500 -0.4500 -0.1500             0.002\n 0.3500 -0.4500 -0.0500             0.002\n 0.3500 -0.4500  0.0500             0.002\n 0.3500 -0.4500  0.1500             0.002\n 0.3500 -0.4500  0.2500             0.002\n 0.3500 -0.4500  0.3500             0.002\n 0.3500 -0.4500  0.4500             0.002\n 0.3500 -0.3500 -0.4500             0.002\n 0.3500 -0.3500 -0.3500             0.002\n 0.3500 -0.3500 -0.2500             0.002\n 0.3500 -0.3500 -0.1500             0.002\n 0.3500 -0.3500 -0.0500             0.002\n 0.3500 -0.3500  0.0500             0.002\n 0.3500 -0.3500  0.1500             0.002\n 0.3500 -0.3500  0.2500             0.002\n 0.3500 -0.3500  0.3500             0.002\n 0.3500 -0.3500  0.4500             0.002\n 0.3500 -0.2500 -0.4500             0.002\n 0.3500 -0.2500 -0.3500             0.002\n 0.3500 -0.2500 -0.2500             0.002\n 0.3500 -0.2500 -0.1500             0.002\n 0.3500 -0.2500 -0.0500             0.002\n 0.3500 -0.2500  0.0500             0.002\n 0.3500 -0.2500  0.1500             0.002\n 0.3500 -0.2500  0.2500             0.002\n 0.3500 -0.2500  0.3500             0.002\n 0.3500 -0.2500  0.4500             0.002\n 0.3500 -0.1500 -0.4500             0.002\n 0.3500 -0.1500 -0.3500             0.002\n 0.3500 -0.1500 -0.2500             0.002\n 0.3500 -0.1500 -0.1500             0.002\n 0.3500 -0.1500 -0.0500             0.002\n 0.3500 -0.1500  0.0500             0.002\n 0.3500 -0.1500  0.1500             0.002\n 0.3500 -0.1500  0.2500             0.002\n 0.3500 -0.1500  0.3500             0.002\n 0.3500 -0.1500  0.4500             0.002\n 0.3500 -0.0500 -0.4500             0.002\n 0.3500 -0.0500 -0.3500             0.002\n 0.3500 -0.0500 -0.2500             0.002\n 0.3500 -0.0500 -0.1500             0.002\n 0.3500 -0.0500 -0.0500             0.002\n 0.3500 -0.0500  0.0500             0.002\n 0.3500 -0.0500  0.1500             0.002\n 0.3500 -0.0500  0.2500             0.002\n 0.3500 -0.0500  0.3500             0.002\n 0.3500 -0.0500  0.4500             0.002\n 0.3500  0.0500 -0.4500             0.002\n 0.3500  0.0500 -0.3500             0.002\n 0.3500  0.0500 -0.2500             0.002\n 0.3500  0.0500 -0.1500             0.002\n 0.3500  0.0500 -0.0500             0.002\n 0.3500  0.0500  0.0500             0.002\n 0.3500  0.0500  0.1500             0.002\n 0.3500  0.0500  0.2500             0.002\n 0.3500  0.0500  0.3500             0.002\n 0.3500  0.0500  0.4500             0.002\n 0.3500  0.1500 -0.4500             0.002\n 0.3500  0.1500 -0.3500             0.002\n 0.3500  0.1500 -0.2500             0.002\n 0.3500  0.1500 -0.1500             0.002\n 0.3500  0.1500 -0.0500             0.002\n 0.3500  0.1500  0.0500             0.002\n 0.3500  0.1500  0.1500             0.002\n 0.3500  0.1500  0.2500             0.002\n 0.3500  0.1500  0.3500             0.002\n 0.3500  0.1500  0.4500             0.002\n 0.3500  0.2500 -0.4500             0.002\n 0.3500  0.2500 -0.3500             0.002\n 0.3500  0.2500 -0.2500             0.002\n 0.3500  0.2500 -0.1500             0.002\n 0.3500  0.2500 -0.0500             0.002\n 0.3500  0.2500  0.0500             0.002\n 0.3500  0.2500  0.1500             0.002\n 0.3500  0.2500  0.2500             0.002\n 0.3500  0.2500  0.3500             0.002\n 0.3500  0.2500  0.4500             0.002\n 0.3500  0.3500 -0.4500             0.002\n 0.3500  0.3500 -0.3500             0.002\n 0.3500  0.3500 -0.2500             0.002\n 0.3500  0.3500 -0.1500             0.002\n 0.3500  0.3500 -0.0500             0.002\n 0.3500  0.3500  0.0500             0.002\n 0.3500  0.3500  0.1500             0.002\n 0.3500  0.3500  0.2500             0.002\n 0.3500  0.3500  0.3500             0.002\n 0.3500  0.3500  0.4500             0.002\n 0.3500  0.4500 -0.4500             0.002\n 0.3500  0.4500 -0.3500             0.002\n 0.3500  0.4500 -0.2500             0.002\n 0.3500  0.4500 -0.1500             0.002\n 0.3500  0.4500 -0.0500             0.002\n 0.3500  0.4500  0.0500             0.002\n 0.3500  0.4500  0.1500             0.002\n 0.3500  0.4500  0.2500             0.002\n 0.3500  0.4500  0.3500             0.002\n 0.3500  0.4500  0.4500             0.002\n 0.4500 -0.4500 -0.4500             0.002\n 0.4500 -0.4500 -0.3500             0.002\n 0.4500 -0.4500 -0.2500             0.002\n 0.4500 -0.4500 -0.1500             0.002\n 0.4500 -0.4500 -0.0500             0.002\n 0.4500 -0.4500  0.0500             0.002\n 0.4500 -0.4500  0.1500             0.002\n 0.4500 -0.4500  0.2500             0.002\n 0.4500 -0.4500  0.3500             0.002\n 0.4500 -0.4500  0.4500             0.002\n 0.4500 -0.3500 -0.4500             0.002\n 0.4500 -0.3500 -0.3500             0.002\n 0.4500 -0.3500 -0.2500             0.002\n 0.4500 -0.3500 -0.1500             0.002\n 0.4500 -0.3500 -0.0500             0.002\n 0.4500 -0.3500  0.0500             0.002\n 0.4500 -0.3500  0.1500             0.002\n 0.4500 -0.3500  0.2500             0.002\n 0.4500 -0.3500  0.3500             0.002\n 0.4500 -0.3500  0.4500             0.002\n 0.4500 -0.2500 -0.4500             0.002\n 0.4500 -0.2500 -0.3500             0.002\n 0.4500 -0.2500 -0.2500             0.002\n 0.4500 -0.2500 -0.1500             0.002\n 0.4500 -0.2500 -0.0500             0.002\n 0.4500 -0.2500  0.0500             0.002\n 0.4500 -0.2500  0.1500             0.002\n 0.4500 -0.2500  0.2500             0.002\n 0.4500 -0.2500  0.3500             0.002\n 0.4500 -0.2500  0.4500             0.002\n 0.4500 -0.1500 -0.4500             0.002\n 0.4500 -0.1500 -0.3500             0.002\n 0.4500 -0.1500 -0.2500             0.002\n 0.4500 -0.1500 -0.1500             0.002\n 0.4500 -0.1500 -0.0500             0.002\n 0.4500 -0.1500  0.0500             0.002\n 0.4500 -0.1500  0.1500             0.002\n 0.4500 -0.1500  0.2500             0.002\n 0.4500 -0.1500  0.3500             0.002\n 0.4500 -0.1500  0.4500             0.002\n 0.4500 -0.0500 -0.4500             0.002\n 0.4500 -0.0500 -0.3500             0.002\n 0.4500 -0.0500 -0.2500             0.002\n 0.4500 -0.0500 -0.1500             0.002\n 0.4500 -0.0500 -0.0500             0.002\n 0.4500 -0.0500  0.0500             0.002\n 0.4500 -0.0500  0.1500             0.002\n 0.4500 -0.0500  0.2500             0.002\n 0.4500 -0.0500  0.3500             0.002\n 0.4500 -0.0500  0.4500             0.002\n 0.4500  0.0500 -0.4500             0.002\n 0.4500  0.0500 -0.3500             0.002\n 0.4500  0.0500 -0.2500             0.002\n 0.4500  0.0500 -0.1500             0.002\n 0.4500  0.0500 -0.0500             0.002\n 0.4500  0.0500  0.0500             0.002\n 0.4500  0.0500  0.1500             0.002\n 0.4500  0.0500  0.2500             0.002\n 0.4500  0.0500  0.3500             0.002\n 0.4500  0.0500  0.4500             0.002\n 0.4500  0.1500 -0.4500             0.002\n 0.4500  0.1500 -0.3500             0.002\n 0.4500  0.1500 -0.2500             0.002\n 0.4500  0.1500 -0.1500             0.002\n 0.4500  0.1500 -0.0500             0.002\n 0.4500  0.1500  0.0500             0.002\n 0.4500  0.1500  0.1500             0.002\n 0.4500  0.1500  0.2500             0.002\n 0.4500  0.1500  0.3500             0.002\n 0.4500  0.1500  0.4500             0.002\n 0.4500  0.2500 -0.4500             0.002\n 0.4500  0.2500 -0.3500             0.002\n 0.4500  0.2500 -0.2500             0.002\n 0.4500  0.2500 -0.1500             0.002\n 0.4500  0.2500 -0.0500             0.002\n 0.4500  0.2500  0.0500             0.002\n 0.4500  0.2500  0.1500             0.002\n 0.4500  0.2500  0.2500             0.002\n 0.4500  0.2500  0.3500             0.002\n 0.4500  0.2500  0.4500             0.002\n 0.4500  0.3500 -0.4500             0.002\n 0.4500  0.3500 -0.3500             0.002\n 0.4500  0.3500 -0.2500             0.002\n 0.4500  0.3500 -0.1500             0.002\n 0.4500  0.3500 -0.0500             0.002\n 0.4500  0.3500  0.0500             0.002\n 0.4500  0.3500  0.1500             0.002\n 0.4500  0.3500  0.2500             0.002\n 0.4500  0.3500  0.3500             0.002\n 0.4500  0.3500  0.4500             0.002\n 0.4500  0.4500 -0.4500             0.002\n 0.4500  0.4500 -0.3500             0.002\n 0.4500  0.4500 -0.2500             0.002\n 0.4500  0.4500 -0.1500             0.002\n 0.4500  0.4500 -0.0500             0.002\n 0.4500  0.4500  0.0500             0.002\n 0.4500  0.4500  0.1500             0.002\n 0.4500  0.4500  0.2500             0.002\n 0.4500  0.4500  0.3500             0.002\n\nAtomic relaxation parameters:\n-----------------------------\nRelaxation method       = 3\nForce tolerance         = 0.000100Ry/a0\nMaximum number of steps = 1000\n\nAtomic data:\n------------\nAtom Type #    1\n----------------\nChemical symbol           : Mg\nCore charge               : 2.000000\nNumber of angular momenta : 3\nNumber of orbitals        : 9\nAtomic radius             : 6.000000\nAtomic energy correction  : -1.435150\n  n l m     energy     occupancy    radius\n  3 0 0   -0.043279     2.000000     6.000000\n  3 1 0   0.349847     0.000000     6.000000\n  3 1 1   0.349847     0.000000     6.000000\n  3 1 2   0.349847     0.000000     6.000000\n  3 2 0   0.878371     0.000000     6.000000\n  3 2 1   0.878371     0.000000     6.000000\n  3 2 2   0.878371     0.000000     6.000000\n  3 2 3   0.878371     0.000000     6.000000\n  3 2 4   0.878371     0.000000     6.000000\n\nOrbital and electron count:\n---------------------------\nNumber of orbitals  = 9\nNumber of electrons = 2.000000\n\nDiagonalizer information:\n-------------------------\nUsing Lapack\n\nFFT integral grid being used.\n-----------------------------\nFFT grid spacing = 0.548509\nFFT grid size    =   12   12   12\n\nscf :  0) E =      -1.735770599 Ry  dRho = 28.0277\nscf :  1) E =      -1.737147192 Ry  dRho = 24.1157\nscf :  2) E =      -1.738884976 Ry  dRho = 2.95901\nscf :  3) E =       -1.73931942 Ry  dRho = 0.869798\nscf :  4) E =      -1.739322087 Ry  dRho = 0.179336\nscf :  5) E =      -1.739322264 Ry  dRho = 0.0386397\nscf :  6) E =      -1.739322269 Ry  dRho = 0.00844228\nscf :  7) E =       -1.73932227 Ry  dRho = 0.000491891\n\nPrimitive translation vectors (a0):\n-----------------------------------\n     5.59200      0.00000      0.00000\n    -1.86400      5.27219      0.00000\n    -1.86400     -2.63609      4.56585\nVolume       = 134.610679 a0^3\nCell repeats =    6    5    4\n\nNumber of atoms = 1\nNeighbor cutoff =     12.00000 a0\n\nAtomic positions (a0):\n----------------------\nMg      0.00000      0.00000      0.00000\n\nTotal forces (Ry/a0):\n---------------------\n                0                 0                 0\n\nStress Tensor (MBars):\n----------------------\n    -0.1017685675   3.274190267e-06   5.671358449e-06\n  3.274190267e-06     -0.1017662523   8.019479171e-06\n  5.671358425e-06    8.01947913e-06      -0.101756992\n\nEFermi =     -0.0767951 Ry\n(kT    =      0.0010000 Ry)\n\nEatom            =     -1.4351500 Ry\nE0               =      0.0384321 Ry\nBand             =     -0.6107085 Ry\nEntropy          =      0.0000486 Ry\nDouble Count     =      0.2681527 Ry\n------------------------------------\nTotal (Free)     =     -1.7393223 Ry\n                 =    -23.6648709 eV\n                 =   -547.1038200 kCal/mol\n                 =  -2288.9481073 kJ/mol\n------------------------------------\nTotal (Internal) =     -1.7392737 Ry\n                 =    -23.6642099 eV\n                 =   -547.0885368 kCal/mol\n                 =  -2288.8841662 kJ/mol\n------------------------------------\nTotal (kT->0)    =     -1.7392980 Ry\n                 =    -23.6645404 eV\n                 =   -547.0961784 kCal/mol\n                 =  -2288.9161368 kJ/mol\n\nTotal energy breakdown\n----------------------\nKinetic                    =      0.7712240 Ry\nNeutral atom               =     -0.6391693 Ry\nElectrostatic (atomic)     =     -1.3967179 Ry\nElectrostatic (correction) =      0.0127755 Ry\nExchange and correlation   =     -0.8886980 Ry\nNon-local pseudopotential  =      0.4065674 Ry\n---------------------------------------------\nSum                        =     -1.7340183 Ry\n\nNet atomic charges (Mulliken):\n------------------------------\nMg(   1) =      0.00000\n\nOrbital populations (Mulliken):\n-------------------------------\nMg(   1):  0.87  0.31  0.31  0.31  0.04  0.04  0.04  0.03  0.03 \n\nMemory allocation (KBytes):\n===========================\nFFTGridR                41\nFFTGridG                54\nMixStack               135\nFFTRho0                 14\nFFTRho                  27\nFFTV                    27\nFFTVna                  14\nhh                       0\nss                       0\nh_k                      1\ns_k                      1\nEigVec                   0\nEigVecK                633\nH                       75\nFock                    37\nRho                     37\nS                       37\nEMatrix                 37\nIntTab                3140\nDiag work arr            0\nACGGridR                 0\nACGGridWeight            0\nACGGridVna               0\nACGGridVtot              0\nACGGridRho               0\nACGGridZeta              0\nACGGridRho0              0\nACGGridRhoFit            0\nACGGridRhoPW             0\nACGGridStack             0\nACGGridMoments           0\nFFTVxc                  14\nFFTVHa                  14\nFFTZeta                  0\nFFTSpinMix               0\n--------------------------\nTotal                 4338\n--------------------------\n\nTime distribution (s):\n======================\nBuildFock                   0.00\nDiagFock                    1.00\nBuildForce                  0.00\nBuild Initial Fock          0.00\nBuild Rho                   0.00\nBuild Density               1.00\nBuild Potential             0.00\n\n\nWall time =       2.0s        0.0m        0.0h\nCPU time  =       2.1s        0.0m        0.0h\n\n"
	with open(filePath, "wt") as f:
		f.write(fileStr)
	return filePath






if __name__=='__main__':
	unittest.main()

