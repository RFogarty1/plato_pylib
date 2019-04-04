#!/usr/bin/python3

from collections import namedtuple
import os
import itertools
import sys
import unittest

import numpy as np
sys.path.append('../..')
import plato_pylib.plato.parse_bas_files as tCode

from collections import namedtuple


class testWriteBasFile(unittest.TestCase):
	def setUp(self):
		self.partBasFile = createPartialBasFileA()

	def tearDown(self):
		os.remove(self.partBasFile)

	def testWrittenFileSameAsParsed(self):
		''' get str-rep of partial file, parse it, write new file based on parsed file, get str-rep of newly written file and see if strings match '''
		tempPath = "temp_file.bas"
		with open(self.partBasFile,"rt") as f:
			expectedStr = f.read()

		parsedFile = tCode.parseBasFile(self.partBasFile)
		tCode.writeBasFileFromParsedDict(tempPath,parsedFile)
		
		with open(tempPath,"rt") as f:
			actualStr = f.read()

		self.assertEqual(expectedStr,actualStr)
	
		os.remove(tempPath)



class testParseBasFiles(unittest.TestCase):

	def setUp(self):
		self.partBasAPath = createPartialBasFileA()
		self.partBasBPath = createPartialBasFileB()

	def tearDown(self):
		os.remove(self.partBasAPath)
		os.remove(self.partBasBPath)

	def testGenInfoParse(self):
		expectedInfo = {"numbShells":6,
		                "ngridpoints":5,
		                "pcc":0,
		                "energy":-1.287901814,
		                "cutoff":6,
		                "zCore": 2,
		                "dNu": 0.03172536568,
		                "hubbard": 0.5,
		                "stoner": 0.05,
		                "xcfunctional":1}
		actualInfo = tCode.parseBasFile(self.partBasAPath)

		for key in expectedInfo.keys():
			self.assertAlmostEqual(expectedInfo[key], actualInfo[key])

	def testParseOrbitals(self):
		# Need to have a load of expected n/l/idx values (3 vals that uniquely define an orbital)
		#Then for each need to attach the correct grid to a dict
		OrbInfoStruct = namedtuple('OrbInfoStruct',['n', 'l', 'zetaIdx', 'occ', 'energy', 'cutoff', 'gridVals'])
		
		gridVals = [ 1.741018362, 1.797138313, 1.855067232, 1.914863431, 1.9765871]
		orbS1GridVals = [ (x,y) for x,y in zip(gridVals, [0.0989413852, 0.0986739085, 0.0981262476, 0.0972919597, 0.0961667723]) ]
		orbP1GridVals = [ (x,y) for x,y in zip(gridVals, [0.0835960866, 0.0841920978, 0.0846214765, 0.0848668996, 0.0849112177]) ]
		orbD1GridVals = [ (x,y) for x,y in zip(gridVals, [0.0682501546, 0.0693816424, 0.0704222457, 0.0713585989, 0.0721760637]) ]
		orbS2GridVals = [ (x,y) for x,y in zip(gridVals, [-0.0325282329, -0.038983002, -0.0440527676, -0.0477287834, -0.0500195748]) ]
		orbP2GridVals = [ (x,y) for x,y in zip(gridVals, [-0.1104856274, -0.1070121274, -0.102854406, -0.0979695372, -0.0923217034]) ]
		orbD2GridVals = [ (x,y) for x,y in zip(gridVals, [-0.1083726475, -0.106048856, -0.1032277781, -0.0998785009, -0.0959702762]) ]

		orbS1 = OrbInfoStruct(n=3, l=0, zetaIdx=1, occ=2, energy=-0.04327913674, cutoff=6, gridVals=orbS1GridVals)
		orbP1 = OrbInfoStruct(n=3, l=1, zetaIdx=1, occ=0, energy=0.3498469233,   cutoff=6, gridVals=orbP1GridVals)
		orbD1 = OrbInfoStruct(n=3, l=2, zetaIdx=1, occ=0, energy=0.8783708167,   cutoff=6, gridVals=orbD1GridVals)
		orbS2 = OrbInfoStruct(n=3, l=0, zetaIdx=2, occ=0, energy=0.0, cutoff=6, gridVals=orbS2GridVals )
		orbP2 = OrbInfoStruct(n=3, l=1, zetaIdx=2, occ=0, energy=0.0, cutoff=6, gridVals=orbP2GridVals )
		orbD2 = OrbInfoStruct(n=3, l=2, zetaIdx=2, occ=0, energy=0.0, cutoff=6, gridVals=orbD2GridVals )

		expectedInfo = [orbS1, orbP1, orbD1, orbS2, orbP2, orbD2]
		actualInfo = tCode.parseBasFile(self.partBasAPath)["orbitals"]

		for exp,act in itertools.zip_longest(expectedInfo,actualInfo):
			self.assertEqual( exp.n, act.n )
			self.assertEqual( exp.l, act.l )
			self.assertEqual( exp.zetaIdx, act.zetaIdx )
			self.assertAlmostEqual( exp.occ, act.occ )
			self.assertAlmostEqual( exp.energy, act.energy )
			self.assertAlmostEqual( exp.cutoff, act.cutoff )
			gridValsExp, gridValsAct = exp.gridVals, act.gridVals
			for rowExp, rowAct in itertools.zip_longest(gridValsExp, gridValsAct):
				[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(rowExp,rowAct)]

	def testGetPlatoRadialRepr(self):
		gridXVals = [5.001975975, 5.162760327, 5.328712959, 5.5]

		orbSGridVals = [ (x,y) for x,y in zip(gridXVals, [0.0023404,0.001394,0.00063806,0.0]) ]
		orbPGridVals = [ (x,y) for x,y in zip(gridXVals, [0.0010947,0.0006362,0.00028348,0.0]) ]
		orbDGridVals = [ (x,y) for x,y in zip(gridXVals, [0.00038715,0.00021987,9.55E-05,0.0]) ]

		expectedInfo = [orbSGridVals, orbPGridVals, orbDGridVals]
		actualOrbs = tCode.parseBasFile(self.partBasBPath)["orbitals"]
		actualInfo = [x.getGridValsPlatoAOrep() for x in actualOrbs]

		for exp,act in itertools.zip_longest(expectedInfo,actualInfo):
			for expRow,actRow in zip(exp,act):
				self.assertAlmostEqual(expRow[0],actRow[0])
				self.assertAlmostEqual(expRow[1],actRow[1])


	def testParseDensity(self):
		gridXVals = [ 1.741018362, 1.797138313, 1.855067232, 1.914863431, 1.9765871]
		expectedDensityVals = [(x,y) for x,y in zip(gridXVals,[0.019578795395808, 0.019473080431912, 0.019257520949328, 0.018931450826999, 0.018496096177719])]

		actualDensities = tCode.parseBasFile(self.partBasAPath)["density"].gridVals

		for exp,act in zip(expectedDensityVals,actualDensities):
			self.assertAlmostEqual(exp[0],act[0])
			self.assertAlmostEqual(exp[1],act[1])

	def testParseVNA(self):
		gridXVals = [ 1.741018362, 1.797138313, 1.855067232, 1.914863431, 1.9765871]
		expectedVNAVals = [-0.5938349748, -0.526176039, -0.4645687067, -0.4085946799, -0.3578461053]
		expectedParsedVNA = [(x,y) for x,y in zip(gridXVals, expectedVNAVals)]
		actualParsedVNA = tCode.parseBasFile(self.partBasAPath)["vna"].gridVals

		for exp,act in itertools.zip_longest(expectedParsedVNA, actualParsedVNA):
			self.assertAlmostEqual(exp[0],act[0])
			self.assertAlmostEqual(exp[1],act[1])

	def testParseVNlGrids(self):
		# Load expected results
		gridXVals = [ 1.741018362, 1.797138313, 1.855067232, 1.914863431, 1.9765871]
		expectedNlPPVals = [ [-0.2470763092, -0.1980057457, -0.1558129653, -0.1202483716, -0.0908952528],
		                     [0.0177225219, 0.0200599937, 0.0201398618, 0.0187164856, 0.0164123448],
		                     [0.316115369, 0.2675482131, 0.2235169044, 0.1841622352, 0.149512821] ]

		expectedParsedNlPP = list()
		for nlGrid in expectedNlPPVals:
			currList = [(x,y) for x,y in zip(gridXVals,nlGrid)]
			expectedParsedNlPP.append( currList ) 
		#Run code
		parsedBas = tCode.parseBasFile(self.partBasAPath)
		actualParsedNlPP = [x.gridVals for x in parsedBas["nlPP"]["grids"]]

		for expList, actList in itertools.zip_longest(expectedParsedNlPP, actualParsedNlPP):
			for exp,act in itertools.zip_longest(expList, actList):
				[self.assertAlmostEqual(x,y) for x,y in itertools.zip_longest(exp,act)]

	def testParsePPInfo(self):
		expLVals, expSigns = [0,0,1], [1,1,1]
		expDict = {"lVals":expLVals, "signVals":expSigns}
		actVals = tCode.parseBasFile(self.partBasAPath)["nlPP"]
		print("actVals = {}".format(actVals))
		print("expVals = {}".format(expDict))
		self.assertTrue( expDict["lVals"] == actVals["lVals"] )
		self.assertTrue( expDict["signVals"] == actVals["signVals"] )

def createPartialBasFileA():
	filePath = os.path.join(os.getcwd(),"partialBasFile.bas")
	fileStr = "6 5 0          -1.287901814                     6                     2         0.03172536568                   0.5                  0.05 1\n3 0                     2        -0.04327913674                     6\n3 1                     0          0.3498469233                     6\n3 2                     0          0.8783708167                     6\n3 0                     0                     0                     6\n3 1                     0                     0                     6\n3 2                     0                     0                     6\n3 0 0 1\n 1 1 1\n          1.741018362      0.09894138516265     0.083596086636076     0.068250154575109    -0.032528232875049     -0.11048562735379     -0.10837264754094     0.019578795395808          -2.348349006         -0.2470763092         0.01772252193           0.316115369         -0.5938349748\n          1.797138313     0.098673908486266     0.084192097775887     0.069381642418538    -0.038983001961287      -0.1070121273882     -0.10604885602777     0.019473080431912          -2.263740608         -0.1980057457         0.02005999374          0.2675482131          -0.526176039\n          1.855067232     0.098126247633668     0.084621476512955     0.070422245692053    -0.044052767633384     -0.10285440595618     -0.10322777805343     0.019257520949328          -2.184067198         -0.1558129653         0.02013986175          0.2235169044         -0.4645687067\n          1.914863431     0.097291959654945     0.084866899576995     0.071358598854195    -0.047728783449601     -0.09796953723218    -0.099878500940624     0.018931450826999          -2.108859893         -0.1202483716         0.01871648558          0.1841622352         -0.4085946799\n            1.9765871     0.096166772270154     0.084911217657724     0.072176063723369      -0.0500195748216    -0.092321703397074    -0.095970276184012     0.018496096177719          -2.037665641        -0.09089525275         0.01641234478           0.149512821         -0.3578461053\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)

	return filePath


def createPartialBasFileB():
	filePath = os.path.join(os.getcwd(), "partialBasFileB.bas")
	fileStr = "3 4 0          -1.134685843                   5.5                     2          0.0316383543                   0.5                  0.05 1\n3 0                     2          0.0615513302                   5.5\n3 1                     0          0.4858438401                   5.5\n3 2                     0           1.060697807                   5.5\n3 0 0 1\n 1 1 1\n          5.001975975    0.0023404211031394    0.0031614977717075    0.0043318895846903   1.0955141880041e-05         -0.7996839689                     0                     0                     0                     0\n          5.162760327    0.0013939629542743       0.0018963228318    0.0026209257169819    3.886265435778e-06         -0.7747793324                     0                     0                     0                     0\n          5.328712959   0.00063805692978326   0.00087212758356027    0.0012124216032395   8.1423329128889e-07         -0.7506503035                     0                     0                     0                     0\n                  5.5                     0                     0                     0                     0     -0.72727272727273                     0                     0                     0                     0\n"
	with open(filePath,"wt") as f:
		f.write(fileStr)

	return filePath

if __name__=='__main__':
	unittest.main()

