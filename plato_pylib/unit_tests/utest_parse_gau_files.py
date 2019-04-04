#!/usr/bin/python3

import itertools
import os
import sys
import unittest
sys.path.append('../..')
import plato_pylib.parse_gau_files as tCode
import plato_pylib.parse_bas_files as parseBas #To get info for a header file

from utest_parse_bas_files import createPartialBasFileA

class testOrbitalParsing(unittest.TestCase):
	def setUp(self):
		self.filePathGauR0_A = createGauFileR0_A()
		self.filePathGauR1_A = createGauFileR1_A()
	def tearDown(self):
		os.remove(self.filePathGauR0_A)
		os.remove(self.filePathGauR1_A)

	def testR0FitParams_OrbsA(self):
		expectedOrbAVals = [(0.2256533391134165345  , 37.734285273902244739), #(exponent,coefficient)
		                    (0.27166909389199239699 , -225.58343293227176218),
		                    (0.32706848862095144748 ,  654.62927000596084781),
		                    (0.39376505702677788712 , -1200.7690751643528984),
		                    (0.47406254509279360798 ,  1526.929211213746612)] 	

		expectedOrbBVals = [(0.2256533391134165345 ,  17.814973037417860979), #(exponent,coefficient)
		                    (0.270830897196620346  , -104.93193089850436195),
		                    (0.32505335469226059875,  301.52446638262307488),
		                    (0.39013157099274675677, -548.60645017933700274),
		                    (0.46823895366151246922,  691.88850762815638973)] 	

		parsedFile = tCode.parseGauFile(self.filePathGauR0_A)
		actualOrbAVals = parsedFile["orbitals"][0]
		actualOrbBVals = parsedFile["orbitals"][1]

		allExpectedOrbs = [expectedOrbAVals, expectedOrbBVals]
		allActualOrbs = parsedFile["orbitals"]

		for expectedOrbVals,actualOrbVals in itertools.zip_longest(allExpectedOrbs,allActualOrbs):
			for primIdx,expPrimVals in enumerate(expectedOrbVals):
				expExponent, expCoeff = expPrimVals[0], expPrimVals[1]		
				actExponent, actCoeff = actualOrbVals.exponents[primIdx],actualOrbVals.r0Coeffs[primIdx]

				self.assertAlmostEqual( expExponent, actExponent )
				self.assertAlmostEqual( expCoeff, actCoeff )

	def testR0Density(self):
		expectedDensity = [ (0.24977446411892198497, 0.18906107877553832153),
		                    (0.29973495506030567448, -0.27126779867644273958),
		                    (0.38482187439714704569, 0.17122695197285703328),
		                    (0.65934295531661524237, -0.18788457858635970732),
		                    (0.88972219525501938797, 0.16118101005396684444) ]
    
		actualDensity = tCode.parseGauFile(self.filePathGauR0_A)["density"]

		for idx,exp in enumerate(expectedDensity):
			self.assertAlmostEqual( exp[0], actualDensity.exponents[idx] )
			self.assertAlmostEqual( exp[1], actualDensity.r0Coeffs[idx] )


	def testR1FitParams_DenAndNeutAtom(self):
		expDenExp = (0.25977874512795540163, 1.189557405762461384)
		expDenCoeffs_R0  = (0.02679475556046136131, -0.026859074525597253363)
		expDenCoeffs_R1  = (0.0046827657367609657954, -0.027478265238950479177)

		expNeutAtomExp = (0.34263537115218700713, 1.3551837469820726678)
		expNeutAtomCoeffs_R0 = (-2.0638099246154166799, -7.4099430193128172917)
		expNeutAtomCoeffs_R1 = (0.029208685956966910086, 0.20618555433380203468)

		expDenObj = tCode.GauPolyBasis( expDenExp, [expDenCoeffs_R0, expDenCoeffs_R1] )
		expNaObj = tCode.GauPolyBasis( expNeutAtomExp, [expNeutAtomCoeffs_R0, expNeutAtomCoeffs_R1] )

		parsedFile = tCode.parseGauFile(self.filePathGauR1_A)
		actualDensity = parsedFile["density"]
		actualNeutAtom = parsedFile["neutAtom".lower()]

		self.assertTrue( actualDensity==expDenObj )
		self.assertTrue( actualNeutAtom==expNaObj )

	def testParseNlPP(self):
		expNlPP = list()
		expExponentsA = [1.4174540372986133008]
		expCoeffsA = [ [-0.2826563852303170532], [-1.595515250025759002] ]
		expExponentsB = [1.4174540055315427711]
		expCoeffsB = [ [-2.4343440773146363121], [0.92424384058112241114] ]
		expExponentsC = [0.99985718149534019705]
		expCoeffsC = [[1.8374666607470948598]]

		expNlA = tCode.GauPolyBasis( expExponentsA, expCoeffsA )
		expNlB = tCode.GauPolyBasis ( expExponentsB, expCoeffsB )
		expNlC = tCode.GauPolyBasis ( expExponentsC, expCoeffsC )

		expNlAll = [expNlA, expNlB, expNlC]
		actNlAll = tCode.parseGauFile(self.filePathGauR1_A)["nlpp"]

		for exp,act in itertools.zip_longest( expNlAll, actNlAll ):
			self.assertTrue( exp==act )

	@unittest.skip("")
	def testWhenSomeFieldsAreEmpty(self):
		''' Need a test for when (for example) orbitals wernt fit (ngauss set to 0 in the *.gin file) '''
		self.assertTrue(False)


class testGeneralInfoParsing(unittest.TestCase):
	def setUp(self):
		self.filePathGauR0_A = createGauFileR0_A()
	def tearDown(self):
		os.remove(self.filePathGauR0_A)

	def testVsKnownValsA(self):
		expectedElement = "Mg"
		expectedShellLVals = [0,1]
		expectedShellNVals = [3,3]

		parsedFile = tCode.parseGauFile(self.filePathGauR0_A)
		actualElement  = parsedFile["element"]
		actualShellLVals = parsedFile["shellAngMoms".lower()]
		actualShellNVals = parsedFile["shellNVals".lower()]


		self.assertEqual( expectedElement, actualElement )
		[self.assertEqual(exp,act) for exp,act in itertools.zip_longest(expectedShellLVals, actualShellLVals)]
		[self.assertEqual(exp,act) for exp,act in itertools.zip_longest(expectedShellNVals, actualShellNVals)]

class testWriteOutputFile(unittest.TestCase):
	def setUp(self):
		self.filePathGauR1_A = createGauFileR1_A()
		self.filePathGauR0_A = createGauFileR0_A()
		self.outFilePath = "testWriteOutFile.gau"
	def tearDown(self):
		os.remove(self.filePathGauR1_A)
		os.remove(self.filePathGauR0_A)
		os.remove(self.outFilePath)

	def testWriteGauFile_gauR1(self):
		fileData = tCode.parseGauFile(self.filePathGauR1_A)
		fileHeader = tCode.getHeaderStrFromGauFile(self.filePathGauR1_A)
		outFilePath = self.outFilePath

		tCode.writeGauFile( outFilePath, fileData, fileHeader )
		writtenData = tCode.parseGauFile(outFilePath)
		
		for key in fileData:
			self.assertEqual( fileData[key], writtenData[key] )

	def testWriteGauFile_gauR0(self):
		fileData = tCode.parseGauFile(self.filePathGauR0_A)
		fileHeader = tCode.getHeaderStrFromGauFile(self.filePathGauR0_A)
		outFilePath = self.outFilePath

		tCode.writeGauFile( outFilePath, fileData, fileHeader )
		writtenData = tCode.parseGauFile(outFilePath)
		
		for key in fileData:
			self.assertEqual( fileData[key], writtenData[key] )


class testGauPolyBasToFunction(unittest.TestCase):

	def testForDensityR1Powers(self):
		expResults = [0.0011548514, 0.0077382712, 0.0147058341]

		testDistances = [0.5,1.0,1.5]
		testExponents = [0.259778745127955,1.18955740576246]
		testCoeffsR0 = [0.0267947555604614, -0.0268590745255972]
		testCoeffsR1 = [0.00468276573676097, -0.0274782652389505]
		polyBas = tCode.GauPolyBasis(testExponents,[testCoeffsR0,testCoeffsR1])
		
		createdFunct = polyBas.toGaussFunct()
		actResults = [createdFunct(x) for x in testDistances]

		print("actResults = {}".format(actResults))
		print("expResults = {}".format(expResults))
		for exp,act in itertools.zip_longest(expResults, actResults):
			self.assertAlmostEqual(exp, act)

class testRemoveSmallCoeffsAndExponents(unittest.TestCase):

	def setUp(self):
		testExpA =      [1e-5, 2.0, 3.0 , 1e-6, 3.0]
		testCoeffR0_A = [2.0 , 3.0, 0.1 , 5.0, 1e-5]
		testCoeffR1_A = [1.0 , 6.0, 0.2 , 4.0, 2.0]
		testCoeffR2_A = [0.0 , 2.0, 1e-6, 2.0, 3.0]
		self.testStructA = tCode.GauPolyBasis( testExpA, [testCoeffR0_A, testCoeffR1_A, testCoeffR2_A] )
		self.expOutputStructA = tCode.GauPolyBasis( [2.0], [ [3.0], [6.0], [2.0] ] )

	def testRemoveSmallCoeffsAndExponents(self):
		minCoeff = 1e-3
		minExp = 1e-4
		self.testStructA.removeSmallCoeffsAndExponents( minCoeff, minExp )
		self.assertEqual( self.testStructA, self.expOutputStructA )

class testParseGauCsv(unittest.TestCase):
	def setUp(self):
		self.csvFileA = createGauCsvFileA()
	def tearDown(self):
		os.remove(self.csvFileA)

	def testAllFieldsFileA(self):
		expVals = loadExpDataCsvFileA()
		actVals = tCode.parseGauCsv(self.csvFileA)
		for key in expVals.keys():
			self.assertEqual( expVals[key], actVals[key] )

class testGetGauHeaderStr(unittest.TestCase):
	def setUp(self):
		basFile = createPartialBasFileA()
		newBasPath = os.path.join( os.getcwd(), "Mg.bas" )
		os.rename( basFile, newBasPath)
		self.parsedFile = parseBas.parseBasFile(newBasPath)
		os.remove(newBasPath)

	def testForBasFileA(self):
		expectedStr = self._loadExpectedGauStr()
		actualStr = tCode.getHeaderStrFromParsedBasFile(self.parsedFile)
		self.assertEqual(expectedStr, actualStr)


	def _loadExpectedGauStr(self):
		fileStr = "Mg 2.000000 6 18\n6.000000 -1.287902\n3 0 0    -0.04327913674                           2\n3 1 0      0.3498469233                           0\n3 1 1      0.3498469233                           0\n3 1 2      0.3498469233                           0\n3 2 0      0.8783708167                           0\n3 2 1      0.8783708167                           0\n3 2 2      0.8783708167                           0\n3 2 3      0.8783708167                           0\n3 2 4      0.8783708167                           0\n3 0 0                 0                           0\n3 1 0                 0                           0\n3 1 1                 0                           0\n3 1 2                 0                           0\n3 2 0                 0                           0\n3 2 1                 0                           0\n3 2 2                 0                           0\n3 2 3                 0                           0\n3 2 4                 0                           0\n3\n0 1\n0 1\n1 1\n"

		return fileStr

def loadExpDataCsvFileA():

	denAct = [(1e-13, 2.4461e-05), (1.9983, 0.020204), (3.112, 0.0050693)]
	denFit = [(1e-13, 0.084333), (1.9983, 0.018447), (3.112,  0.0021148)]
 
	neutAtomAct = [(1e-13, -9.062), (1.9983, -0.35785), (3.112, -0.024727)]
	neutAtomFit = [(1e-13, -8.5913), (1.9983, -0.28392), (3.112, -0.0022019)]

	orbActA =  [(1e-13, 0.0034972), (1.9983, 0.10051), (3.112, 0.050346)]
	orbFitA = [(1e-13, 0.34538), (1.9983, 0.075536), (3.112, 0.0086574)]

	orbActB = [(1e-13, 0.14562), (1.9983, 0.079541), (3.112, 0.031147)]
	orbFitB = [(1e-13, 0.42622), (1.9983, 0.093202), (3.112, 0.010679)]

	orbActC = [(1e-13, 0.19476), (1.9983, 0.045897), (3.112, 0.014045)]
	orbFitC = [(1e-13, 0.3036), (1.9983, 0.066399), (3.112, 0.00761)]

	nlActA =[(1e-13, -0.28266), (1.9983, -0.023163), (3.112, -1.7187e-05)]
	nlFitA = [(1e-13, -0.49161), (1.9983, -0.052189), (3.112,-0.0021353)]

	nlActB = [(1e-13, -2.4343), (1.9983, 0.0043739), (3.112, 7.1182e-06)]
	nlFitB = [(1e-13, -2.4739), (1.9983, -0.00078248), (3.112, -8.0493e-09)]

	outDict = dict()
	densityVals = tCode.GauCsvGridInfo(denAct,denFit) 
	neutAtomVals = tCode.GauCsvGridInfo(neutAtomAct,neutAtomFit)
	orbitalVals = [ tCode.GauCsvGridInfo(orbActA,orbFitA),
	                tCode.GauCsvGridInfo(orbActB,orbFitB),
	                tCode.GauCsvGridInfo(orbActC,orbFitC) ]
	nlVals = [ tCode.GauCsvGridInfo(nlActA,nlFitA),
	           tCode.GauCsvGridInfo(nlActB,nlFitB) ]

	outDict["density"] = densityVals
	outDict["neutAtom".lower()] = neutAtomVals
	outDict["orbitals"] = orbitalVals
	outDict["nlpp"] = nlVals

	return outDict

def createGauFileR0_A():
	filePath = os.path.join( os.getcwd(), "gaufileR0_Only.gau" )
	fileStr = "Mg 2.000000 2 4\n7.500000 -1.512379\n3 0 0     -0.2102929769                           2\n3 1 0       0.121544206                           0\n3 1 1       0.121544206                           0\n3 1 2       0.121544206                           0\n3\n0 1\n0 1\n1 1\n#Fitting parameters - wavefunction\n#n = 3 l = 0 occ = 2.000000 eigenvalue =     -0.2102929769\n#               a                           c\n5\n      0.2256533391134165345       37.734285273902244739\n     0.27166909389199239699      -225.58343293227176218\n     0.32706848862095144748       654.62927000596084781\n     0.39376505702677788712      -1200.7690751643528984\n     0.47406254509279360798        1526.929211213746612\n#Fitting parameters - wavefunction\n#n = 3 l = 1 occ = 0.000000 eigenvalue =       0.121544206\n#               a                           c\n5\n      0.2256533391134165345       17.814973037417860979\n       0.270830897196620346      -104.93193089850436195\n     0.32505335469226059875       301.52446638262307488\n     0.39013157099274675677      -548.60645017933700274\n     0.46823895366151246922       691.88850762815638973\n#Fitting parameters - density\n#               a                           c\n5 1\n     0.24977446411892198497      0.18906107877553832153\n     0.29973495506030567448     -0.27126779867644273958\n     0.38482187439714704569      0.17122695197285703328\n     0.65934295531661524237     -0.18788457858635970732\n     0.88972219525501938797      0.16118101005396684444\n#Fitting parameters - neutral atom potential\n#               a                           c\n5 1\n     0.22031261386712974737     -0.12219239318591164356\n     0.30806498025338135971      -1.1844526542350564124\n     0.43822264944404770715      0.77419438726540679152\n     0.59723068430107395521      -3.4482242270734384526\n     0.77866229002160136652       3.9000745852127369773\n#Fitting parameters - non-local pseudopotential\n# l = 0\n#               a                           c\n1 2\n      1.4174540372999995252     -0.28265638523037667218       -1.595515250029995169\n#Fitting parameters - non-local pseudopotential\n# l = 0\n#               a                           c\n1 2\n      1.4174540055336761757       -2.434344077312082355      0.92424384057699382478\n#Fitting parameters - non-local pseudopotential\n# l = 1\n#               a                           c\n1 1\n     0.99985718149187763348       1.8374666607343681513\n#Fitting parameters - dExc\n#               a                           c\n0 0\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath


def createGauFileR1_A():
	filePath = os.path.join( os.getcwd(), "gaufileR1_A.gau" )
	fileStr = "Mg 2.000000 2 4\n7.500000 -1.512379\n3 0 0     -0.2102929769                           2\n3 1 0       0.121544206                           0\n3 1 1       0.121544206                           0\n3 1 2       0.121544206                           0\n3\n0 1\n0 1\n1 1\n#Fitting parameters - wavefunction\n#n = 3 l = 0 occ = 2.000000 eigenvalue =     -0.2102929769\n#               a                           c\n1\n     0.24563635546590129044      0.24867353955701998469\n#Fitting parameters - wavefunction\n#n = 3 l = 1 occ = 0.000000 eigenvalue =       0.121544206\n#               a                           c\n1\n      0.2456647894110553787      0.24652940695742253663\n#Fitting parameters - density\n#               a                           c\n2 2\n     0.25977874512795540163      0.02679475556046136131    0.0046827657367609657954\n       1.189557405762461384    -0.026859074525597253363    -0.027478265238950479177\n#Fitting parameters - neutral atom potential\n#               a                           c\n2 2\n     0.34263537115218700713      -2.0638099246154166799     0.029208685956966910086\n      1.3551837469820726678      -7.4099430193128172917      0.20618555433380203468\n#Fitting parameters - non-local pseudopotential\n# l = 0\n#               a                           c\n1 2\n      1.4174540372986133008      -0.2826563852303170532       -1.595515250025759002\n#Fitting parameters - non-local pseudopotential\n# l = 0\n#               a                           c\n1 2\n      1.4174540055315427711      -2.4343440773146363121      0.92424384058112241114\n#Fitting parameters - non-local pseudopotential\n# l = 1\n#               a                           c\n1 1\n     0.99985718149534019705       1.8374666607470948598\n#Fitting parameters - dExc\n#               a                           c\n0 0\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

def createGauCsvFileA():
	filePath = os.path.join( os.getcwd(), "gauFileTestA_gau.csv" )
	fileStr = "Wavefunction, n=3, l=0, occ=2.000000, eigenvalue=    0.06155133048, ngauss=1\nR, Original, Fit\n       1e-13,    0.0034972,      0.34538\n      1.9983,      0.10051,     0.075536\n       3.112,     0.050346,    0.0086574\n\nWavefunction, n=3, l=1, occ=0.000000, eigenvalue=     0.4858438401, ngauss=1\nR, Original, Fit\n       1e-13,      0.14562,      0.42622\n      1.9983,     0.079541,     0.093202\n       3.112,     0.031147,     0.010679\n\nWavefunction, n=3, l=2, occ=0.000000, eigenvalue=      1.060697807, ngauss=1\nR, Original, Fit\n       1e-13,      0.19476,       0.3036\n      1.9983,     0.045897,     0.066399\n       3.112,     0.014045,      0.00761\n\nDensity, ngauss=1, npoly=1\nR, Original, Fit\n       1e-13,   2.4461e-05,     0.084333\n      1.9983,     0.020204,     0.018447\n       3.112,    0.0050693,    0.0021148\n\nNeutral atom potential, ngauss=1, npoly=1\nR, Original, Fit\n       1e-13,       -9.062,      -8.5913\n      1.9983,     -0.35785,     -0.28392\n       3.112,    -0.024727,   -0.0022019\n\nNon-local pseudopotential, ngauss=1, npoly=1\nR, Original, Fit\n       1e-13,     -0.28266,     -0.49161\n      1.9983,    -0.023163,    -0.052189\n       3.112,  -1.7187e-05,   -0.0021353\n\nNon-local pseudopotential, ngauss=1, npoly=1\nR, Original, Fit\n       1e-13,      -2.4343,      -2.4739\n      1.9983,    0.0043739,  -0.00078248\n       3.112,   7.1182e-06,  -8.0493e-09\n\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)
	return filePath

if __name__ == '__main__':
	unittest.main()

