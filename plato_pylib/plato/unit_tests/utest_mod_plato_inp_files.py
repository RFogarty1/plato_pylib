#!/usr/bin/python3

import os
import sys
import unittest

sys.path.append('../..')
import plato_pylib.plato.mod_plato_inp_files as tCode
import plato_pylib.shared.ucell_class as UCell


class testModPlatoInpOption(unittest.TestCase):
	def setUp(self):
		self.filePathDict  = dict()
		self.fileStrDict = dict()
		self.filePathDict["testInpFragA"], self.fileStrDict["testInpFragA"] = createPartialPlatoInpFileA()

	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]

	def testOneLineMod(self):
		''' Test modPlatoInpOption works for a simple 1 line case (RFlag). NOTE: dependent '''
		''' on tCode.tokenizePlatoInpFile '''
		testInpFile = self.filePathDict["testInpFragA"]
		testOptVal = {'RFlag':5}
		key = 'RFlag'.lower()
		expectedVal = "5"
		tCode.modPlatoInpOptions(testInpFile, testOptVal)

		parsedOptions = tCode.tokenizePlatoInpFile(testInpFile)
		self.assertEqual( expectedVal, parsedOptions[key] )

	def testMultiLineMod(self):
		''' Test modPlatoInpOption works for replacing option that runs over multiple lines'''
		testInpFile = self.filePathDict["testInpFragA"]
		expectedVal = "2\n1\n2"
		testOptVal = {'Pinned':expectedVal}
		key = 'Pinned'.lower()
		tCode.modPlatoInpOptions(testInpFile, testOptVal)

		parsedOptions = tCode.tokenizePlatoInpFile(testInpFile)
		self.assertEqual( expectedVal, parsedOptions[key] )

	def testMutliMods(self):
		''' Test modPlatoInpOption works when replacing >1 option '''
		testInpFile = self.filePathDict["testInpFragA"]
		expectedVals = ["5", "2\n1\n2"]
		keys = ["pinned", "rflag"]
		testOptVals = {k:v for k,v in zip(keys, expectedVals)}
		tCode.modPlatoInpOptions(testInpFile, testOptVals)

		parsedOptions = tCode.tokenizePlatoInpFile(testInpFile)
		[self.assertEqual(expectedVals[x], parsedOptions[keys[x]]) for x in range(len(expectedVals))]


		
class testTokenizePlatoInpFile(unittest.TestCase):
	def setUp(self):
		self.filePathDict  = dict()
		self.fileStrDict = dict()
		self.filePathDict["testInpFragA"], self.fileStrDict["testInpFragA"] = createPartialPlatoInpFileA()

	def tearDown(self):
		[os.remove(x) for x in self.filePathDict.values()]
	
	def testTokenizeInpFragmentA(self):
		testInpFile = self.filePathDict["testInpFragA"]
		expectedVals = { "rflag" : "3", 
		                 "cellrelaxmode" : "0",
		                 "FTol".lower() : "0.0001", 
		                 "MaxDisplacement".lower(): "0.1",
		                 "StressTol".lower(): "0.001",
		                 "Pinned".lower(): "0\n1",
		                 "Constrain".lower(): "0\n1 1.0 1.0 0.0",
		                 "ReactionCoordinate".lower():"0"}
		actualVals = tCode.tokenizePlatoInpFile(testInpFile)
		self.assertEqual(expectedVals, actualVals)


class testGetGeomStrFromUnitCell(unittest.TestCase):

	def setUp(self):
		fractCoordsA = [ [0.0000000, 0.0000000, .00000000,  "Si"], [0.25, 0.25, 0.25, "Si"] ]
		lattVectsA = [ [2.0, 0.0, 0.0], [0.0,2.0,0.0], [0.0,0.0,2.0] ]
		self.testUCellA = UCell.UnitCell.fromLattVects(lattVectsA,fractCoords=fractCoordsA)

		cellVecForm = "{:.8f} {:.8f} {:.8f}\n"
		fractForm = "{:.8f} {:.8f} {:.8f} {}\n"
		cellSizeForm = "{:.8f} {:.8f} {:.8f}"
		self.expDictA = {"cellvec": cellVecForm.format(1.0,0.0,0.0) +
		                            cellVecForm.format(0.0,1.0,0.0) +
		                            cellVecForm.format(0.0,0.0,1.0),
		                 "cellsize":cellSizeForm.format(2.0,2.0,2.0),
		                 "natom":"2",
		                 "format":"0",
		                 "atoms": fractForm.format(0.0, 0.0, 0.0, "Si") + 
		                          fractForm.format(0.25,0.25,0.25,"Si") }

	def testForGeomA(self):
		actDict = tCode.getPlatoGeomDictFromUnitCell(self.testUCellA)
		self.assertEqual(self.expDictA, actDict)



class testGetDefaultDicts(unittest.TestCase):

	def setUp(self):
		self.expTb1StrDict = getDefaultTb1StrDict()
		self.expTb2StrDict = getDefaultTb2StrDict() #Subset of the tb1 dict
		self.expDftStrDict = getDefaultDftStrDict()

	def testDefaultTb1StrDict(self):
		optDict = tCode.getDefOptDict("tb1")
		actStrDict = tCode.getStrDictFromOptDict(optDict,"tb1")
		self.assertEqual(actStrDict,self.expTb1StrDict)

	def testDefaultTb2StrDict(self):
		optDict = tCode.getDefOptDict("tb2")
		actStrDict = tCode.getStrDictFromOptDict(optDict,"tb2")
		self.assertEqual(actStrDict, self.expTb2StrDict)

	def testDefaultDftStrDict(self):
		optDict = tCode.getDefOptDict("dft")
		actStrDict = tCode.getStrDictFromOptDict(optDict,"dft")
		self.assertEqual(actStrDict, self.expDftStrDict)


def createPartialPlatoInpFileA():
	filePath = os.path.join( os.getcwd(), "platoInpFileBasic" )
	fileStr = "################################################################################\n# ------- Parameters for relaxation simulations -------\n#\n# RFlag is the relaxation method flag.\n#  1  ==>  Steepest descent\n#  2  ==>  Conjugate gradient\n#  3  ==>  Variable metric\n#  4  ==>  Find transition state using mode following\n#  5  ==>  Relax cell and atomic coordinates\n#  6  ==>  BFGS\n# Making RFlag negative means numeric rather than analytic forces are used.    \n#\n# CellRelaxMode defines the way the cell vectors are relaxed:\n#  0  ==>  Free relaxation\n#  1  ==>  Volume relaxed at fixed shape\n#  2  ==>  Relax along axis 1 only\n#  3  ==>  Relax along axis 2 only\n#  4  ==>  Relax along axis 3 only\n#\n# FTol is the highest force allowed in a relaxed structure. Units are\n# rydbergs/bohr radius.\n#\n# MaxDisplacement is the maximum distance any atom can move in any one update.\n# Currently only used by variable metric minimization.\n#\n# StressTol is the highest trace of the stress allowed in a relaxed structure.\n# Units are MBars.\n#\n# Pinned determines which atoms are pinned. Variables are:\n# Number of pinned atoms and the index for each pinned atom.\n#\n# Constrain determines constraints to be applied to atomic motion. Variables\n# are:\n# Number of constrained atoms \n# The index of atom and normal to plane of constraint\n#\n# ReactionCoordinate defines the constraint direction for a drag calculation of\n# a minimum energy pathway. Variables are:\n# Type\n# Coordinates of atoms\n#\n# MaxNSteps: see molecular dynamics section\n################################################################################\nRFlag\n3\n\nCellRelaxMode\n0\n\nFTol\n0.0001\n\nMaxDisplacement\n0.1\n\nStressTol\n0.001\n\nPinned\n0\n1\n\nConstrain\n0\n1 1.0 1.0 0.0\n\nReactionCoordinate\n0\n"
	with open(filePath, "wt") as f: 
		f.write(fileStr)
	return filePath, fileStr


def getDefaultTb1StrDict():
	outDict = {'atomtemperature': '300.0',
	 'autosave': '10',
	 'basisflag': '0',
	 'blochstates': '0',
	 'cellrelaxmode': '0',
	 'cellrepeat': '1 1 1',
	 'constrain': '0',
	 'crystalfieldxcweight': '1.0',
	 'electronexcess': '0.0',
	 'electrontemperature': '0.001',
	 'etol': '0.0',
	 'excmbcorr': '0',
	 'externalfield': '0',
	 'ftol': '0.0001',
	 'hamiltonianflag': '1',
	 'integralmesh': '1\n100 30 30',
	 'job': '1',
	 'kickatom': '0',
	 'lsdapuflag': '0',
	 'magneticfield': '0.0',
	 'maxdisplacement': '0.1',
	 'maximumtime': '1000000.0',
	 'maxnsteps': '1000',
	 'mdtimestep': '1.0',
	 'mixclipping': '0.01 0.1',
	 'mixfactor': '0.2',
	 'mixlevels': '5',
	 'mixscheme': '1',
	 'movierate': '0',
	 'nbands': '0',
	 'nloops': '30',
	 'nonetforceflag': '1',
	 'nosehoover': '-1.0',
	 'occupyflag': '1',
	 'occupytolerance': '1e-12',
	 'onsiteelectrostatics': '0',
	 'partialpopulationflag': '0',
	 'pbcflag': '1',
	 'pinned': '0',
	 'preconditioningflag': '1 1',
	 'reactioncoordinate': '0',
	 'residuetol': '0.001',
	 'restartflag': '0',
	 'rflag': '3',
	 'scfflag': '0',
	 'seed': '12345',
	 'spinflag': '0',
	 'spinmixfactor': '0.2',
	 'stresstol': '0.001',
	 'temperaturetolerance': '0.1',
	 'testoverlap': '1',
	 'vanderwaals': '0',
	 'verbosity': '1',
	 'vxcmbcorrxtal': '0 1.00000000',
	 'writeeigenvects': '0',
	 'writeenergyflag': '0',
	 'xcfunctional': '1'}
	
	return outDict



def getDefaultTb2StrDict():
	outDict = {'atomtemperature': '300.0',
	 'autosave': '10',
	 'basisflag': '0',
	 'blochstates': '0',
	 'cellrelaxmode': '0',
	 'cellrepeat': '1 1 1',
	 'constrain': '0',
	 'e0method': '0',
	 'electronexcess': '0.0',
	 'electrontemperature': '0.001',
	 'etol': '0.0',
	 'excmbcorr': '0',
	 'externalfield': '0',
	 'ftol': '0.0001',
	 'hamiltonianflag': '1',
	 'hopxcmethod' : '0',
	 'integralmesh': '1\n100 30 30',
	 'inversesk': '1',
	 'job': '1',
	 'kickatom': '0',
	 'lsdapuflag': '0',
	 'magneticfield': '0.0',
	 'maxdisplacement': '0.1',
	 'maximumtime': '1000000.0',
	 'maxnsteps': '1000',
	 'mcwedaxcflag': '0',
	 'mdtimestep': '1.0',
	 'mixclipping': '0.01 0.1',
	 'mixfactor': '0.2',
	 'mixlevels': '5',
	 'mixscheme': '1',
	 'movierate': '0',
	 'nbands': '0',
	 'nloops': '30',
	 'nonetforceflag': '1',
	 'nosehoover': '-1.0',
	 'occupyflag': '1',
	 'occupytolerance': '1e-12',
	 'onsiteelectrostatics': '0',
	 'partialpopulationflag': '0',
	 'pbcflag': '1',
	 'pinned': '0',
	 'preconditioningflag': '1 1',
	 'reactioncoordinate': '0',
	 'residuetol': '0.001',
	 'restartflag': '0',
	 'rflag': '3',
	 'scfflag': '0',
	 'seed': '12345',
	 'spinflag': '0',
	 'spinmixfactor': '0.2',
	 'stresstol': '0.001',
	 'temperaturetolerance': '0.1',
	 'testoverlap': '1',
	 'verbosity': '1',
	 'writeeigenvects': '0',
	 'writeenergyflag': '0',
	 'xcfunctional': '1',
	 'xtalxcmethod': '0'}
	
	return outDict


def getDefaultDftStrDict():
	outDict = {'acgangularmeshtype': '1',
	 'acgewaldparameters': '0.4 10.0',
	 'acggridspacing': '30 29 3',
	 'acgmeshreduction': '6 5.0',
	 'acgminpartitionwt': '0.0',
	 'atomtemperature': '300.0',
	 'autosave': '10',
	 'basisflag': '0',
	 'blochstates': '0',
	 'cellrelaxmode': '0',
	 'cellrepeat': '1 1 1',
	 'constrain': '0',
	 'densityfit': '1 0',
	 'densityrwflag': '0',
	 'diagonalisationmethod': '1',
	 'diagonalisationworksize': '500',
	 'diameternlv': '0.0',
	 'electronexcess': '0.0',
	 'electrontemperature': '0.001',
	 'etol': '0.0',
	 'excmbcorr': '0',
	 'externalfield': '0',
	 'fftgridspacing': '0.7',
	 'ftol': '0.0001',
	 'hamiltonianflag': '1',
	 'hyperfineflag': '0',
	 'integralmesh': '0',
	 'integralpartitionflag': '0',
	 'job': '1',
	 'kickatom': '0',
	 'lsdapuflag': '0',
	 'magneticfield': '0.0',
	 'maxbond': '6.0',
	 'maxdisplacement': '0.1',
	 'maximumtime': '1000000.0',
	 'maxnsteps': '1000',
	 'mdtimestep': '1.0',
	 'mixclipping': '0.01 0.1',
	 'mixfactor': '0.2',
	 'mixlevels': '5',
	 'mixmetric': '0 20.0',
	 'mixscheme': '1',
	 'mixthreshold': '100',
	 'model': '2',
	 'movierate': '0',
	 'nbands': '0',
	 'nloops': '30',
	 'nonetforceflag': '1',
	 'nosehoover': '-1.0',
	 'numericvna': '0',
	 'occupyflag': '1',
	 'occupytolerance': '1e-12',
	 'onsiteelectrostatics': '0',
	 'optimizemesh': '1',
	 'parallelisefft': '0',
	 'partialpopulationflag': '0',
	 'pbcflag': '1',
	 'pinned': '0',
	 'preconditioningflag': '1 1',
	 'reactioncoordinate': '0',
	 'residuetol': '0.001',
	 'restartflag': '0',
	 'rflag': '3',
	 'scfflag': '0',
	 'seed': '12345',
	 'spinflag': '0',
	 'spinmixfactor': '0.2',
	 'spinmixlevels': '5',
	 'spinmixscheme': '1',
	 'splitdensityflag': '0',
	 'stresstol': '0.001',
	 'temperaturetolerance': '0.1',
	 'testoverlap': '1',
	 'verbosity': '1',
	 'wavefunctionflag': '0',
	 'writeatomdensityflag': '0',
	 'writeeigenvects': '0',
	 'writeenergyflag': '0',
	 'xcflag': '1',
	 'xcfunctional': '1'}
	
	return outDict


if __name__ == '__main__':
	unittest.main()




