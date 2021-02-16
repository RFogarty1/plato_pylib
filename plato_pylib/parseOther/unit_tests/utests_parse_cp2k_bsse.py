
import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.parseOther.parse_cp2k_files as tCode
import plato_pylib.shared.energies_class as energiesHelp
import plato_pylib.shared.unit_convs as uConvHelp


class TestParseBSSEFragmentsFullFile(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = _loadTestFileA().split("\n")

	@mock.patch("plato_pylib.parseOther.parse_cp2k_files._getFileAsListFromInpFile")
	def testExpectedOutDictA(self, getFileAsListMock):
		getFileAsListMock.side_effect = lambda *args,**kwargs: self.fileAsListA
		outDict = tCode.parseCpout(mock.Mock())
		expFragList = self._getExpFragmentsDicts()
		actFragList = outDict["bsse_fragments"]
		self._checkExpAndActDictsMatch(expFragList, actFragList)


	def _getExpFragmentsDicts(self):
		#First get all the total energies in anoying energies objects
		allEnergies = [-14.44625651841509, -17.91750786640011, -14.44633305634415, -17.91768064987254, -32.36423479187624]
		energyObjs = [energiesHelp.EnergyVals(dftTotalElectronic=x*uConvHelp.RYD_TO_EV*2) for x in allEnergies]

		#Now get all the fragment info
		fragA = { "conf": "10", "frag_sub_conf": "10", "charge": 0, "multiplicity": 0, 
		          "kinds": ["O","H","H"], "indices":[1,2,3] }
		fragB = { "conf": "01", "frag_sub_conf": "01", "charge": 0, "multiplicity": 0,
		          "kinds": ["O","H","H"], "indices":[4,5,6] }
		fragC = { "conf": "11", "frag_sub_conf": "10", "charge": 0, "multiplicity": 0,
		          "kinds": ["O", "H", "H", "O_ghost", "H_ghost", "H_ghost"], "indices": [1,2,3,4,5,6] }
		fragD = { "conf": "11", "frag_sub_conf": "01", "charge": 0, "multiplicity": 0,
		          "kinds": ["O_ghost", "H_ghost", "H_ghost", "O", "H", "H"], "indices":[1,2,3,4,5,6] }
		fragE = { "conf": "11", "frag_sub_conf": "11", "charge": 0, "multiplicity": 0,
		          "kinds": ["O","H","H","O","H","H"], "indices":[1,2,3,4,5,6] }
		allFrags = [fragA, fragB, fragC, fragD, fragE]

		outDicts = list()
		for eObj, frag in it.zip_longest(energyObjs, allFrags):
			frag["energies"] = eObj
			outDicts.append(frag)
		return outDicts

	def _checkExpAndActDictsMatch(self, exp, act):
		self.assertEqual( len(exp), len(act) )

		cmpKeys = ["conf", "frag_sub_conf", "charge", "multiplicity", "kinds", "indices"]

		for e,a in zip(exp,act):
			for key in cmpKeys:
				self.assertEqual( e[key], a[key] )
			self.assertAlmostEqual( e["energies"].dftTotalElectronic, a["energies"].dftTotalElectronic )


class TestParseSingleBSSEFragment(unittest.TestCase):

	def setUp(self):
		self.fragA = """ -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 01        FRAGMENT SUBCONF: 01     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      4                                   O                  -
 -                      5                                   H                  -
 -                      6                                   H                  -
 -------------------------------------------------------------------------------

"""
		self.lineIdxA = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsList = self.fragA.split("\n")

	def testExpectedA(self):
		expOutDict = {"kinds":["O","H","H"], "indices":[4,5,6], "conf":"01", "charge":0,
		              "frag_sub_conf":"01", "multiplicity":0}
		expEndLineIdx = 10
		actOutDict, actEndLineIdx = tCode._parseBSSEFragmentsInfo(self.fileAsList, self.lineIdxA)
		self.assertEqual(expEndLineIdx, actEndLineIdx)
		self.assertEqual(expOutDict, actOutDict)


class TestHandleBSSEFragParsedDict(unittest.TestCase):

	def setUp(self):
		self.currOutDict = {}
		self.parsedOutDict = {"kinds":["O","H","H"]}
		self.createTestObjs()

	def createTestObjs(self):
		self.parser = mock.Mock()
		self.parser.outDict = self.currOutDict	

	def _runTestFunct(self):
		tCode._handleParsedBSSEFragsInfo(self.parser, self.parsedOutDict)

	def testFirstFrag(self):
		expDict = {"bsse_fragments":[self.parsedOutDict]}
		self._runTestFunct()
		actDict = self.parser.outDict
		self.assertEqual(expDict,actDict)

	def testExpectedA_energiesPresent(self):
		expEnergies = mock.Mock()
		self.currOutDict["energies"] = expEnergies
		self.currOutDict["bsse_fragments"] = [ {} ]
#		expSecondDict = copy.deepcopy(self.parsedOutDict)
#		expSecondDict["energies"] = expEnergies

		expDicts = [ {"energies":expEnergies}, self.parsedOutDict ]
		self._runTestFunct()
		actDicts = self.parser.outDict["bsse_fragments"]
		self.assertEqual(expDicts,actDicts)
	

def _loadTestFileA():
	return """
 DBCSR| Multiplication driver                                               BLAS
 DBCSR| Multrec recursion limit                                              512
 DBCSR| Multiplication stack size                                           1000
 DBCSR| Maximum elements for images                                    UNLIMITED
 DBCSR| Multiplicative factor virtual images                                   1
 DBCSR| Multiplication size stacks                                             3


  **** **** ******  **  PROGRAM STARTED AT               2020-10-10 10:56:35.508
 ***** ** ***  *** **   PROGRAM STARTED ON                              mt-rf614
 **    ****   ******    PROGRAM STARTED BY                                 rf614
 ***** **    ** ** **   PROGRAM PROCESS ID                                186225
  **** **  *******  **  PROGRAM STARTED IN /media/ssd1/rf614/Work/temp/jul_2020/
                                           cp2k_bsse

 CP2K| version string:                                          CP2K version 6.1
 CP2K| source code revision number:                                    svn:18464
 CP2K| cp2kflags: fftw3 max_contr=4                                             
 CP2K| is freely available from                            https://www.cp2k.org/
 CP2K| Program compiled at                          Wed 23 Sep 11:44:25 BST 2020
 CP2K| Program compiled on                                              mt-rf614
 CP2K| Program compiled for                                Linux-x86-64-gfortran
 CP2K| Data directory path             /media/ssd1/rf614/Work/CP2K/cp2k-6.1/data
 CP2K| Input file name                                             2H2O_bsse.inp

 GLOBAL| Force Environment number                                              1
 GLOBAL| Basis set file name                                           BASIS_SET
 GLOBAL| Potential file name                                           POTENTIAL
 GLOBAL| MM Potential file name                                     MM_POTENTIAL
 GLOBAL| Coordinate file name                                      __STD_INPUT__
 GLOBAL| Method name                                                        CP2K
 GLOBAL| Project name                                                   H2O_BSSE
 GLOBAL| Preferred FFT library                                             FFTW3
 GLOBAL| Preferred diagonalization lib.                                       SL
 GLOBAL| Run type                                                           BSSE
 GLOBAL| All-to-all communication in single precision                          F
 GLOBAL| FFTs using library dependent lengths                                  F
 GLOBAL| Global print level                                                  LOW
 GLOBAL| Total number of message passing processes                             1
 GLOBAL| Number of threads for this process                                    1
 GLOBAL| This output is from process                                           0
 GLOBAL| CPU model name :  Intel(R) Xeon(R) CPU E5-1620 v4 @ 3.50GHz

 MEMORY| system memory details [Kb]
 MEMORY|                        rank 0           min           max       average
 MEMORY| MemTotal             65866104      65866104      65866104      65866104
 MEMORY| MemFree              57879552      57879552      57879552      57879552
 MEMORY| Buffers                357856        357856        357856        357856
 MEMORY| Cached                1322268       1322268       1322268       1322268
 MEMORY| Slab                   315048        315048        315048        315048
 MEMORY| SReclaimable           216280        216280        216280        216280
 MEMORY| MemLikelyFree        59775956      59775956      59775956      59775956


 GENERATE|  Preliminary Number of Bonds generated:                             4
 GENERATE|  Achieved consistency in connectivity generation.
 GENERATE|  Number of Bonds generated:                                         4
 GENERATE|  Preliminary Number of Bends generated:                             2
 GENERATE|  Number of Bends generated:                                         2
 GENERATE|  Number of UB generated:                                            2
 GENERATE|  Preliminary Number of Torsions generated:                          0

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   2
                             - Atoms:                                          6
                             - Shell sets:                                    12
                             - Shells:                                        22
                             - Primitive Cartesian functions:                 30
                             - Cartesian basis functions:                     48
                             - Spherical basis functions:                     46

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 10        FRAGMENT SUBCONF: 10     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      1                                   O                  -
 -                      2                                   H                  -
 -                      3                                   H                  -
 -------------------------------------------------------------------------------

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   2
                             - Atoms:                                          3
                             - Shell sets:                                     6
                             - Shells:                                        11
                             - Primitive Cartesian functions:                 15
                             - Cartesian basis functions:                     24
                             - Spherical basis functions:                     23

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 Number of electrons:                                                          8
 Number of occupied orbitals:                                                  4
 Number of molecular orbitals:                                                 4

 Number of orbital functions:                                                 23
 Number of independent orbital functions:                                     23

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 P_Mix/Diag. 0.40E+00    0.0     1.25505937       -14.3100140078 -1.43E+01
     2 P_Mix/Diag. 0.40E+00    0.0     0.72506679       -14.3963683093 -8.64E-02
     3 P_Mix/Diag. 0.40E+00    0.0     0.43171406       -14.4462565184 -4.99E-02

  Leaving inner SCF loop after reaching     3 steps.


  Electronic density on regular grids:         -7.9810272194        0.0189727806
  Core density on regular grids:                8.0604493642        0.0604493642
  Total charge density on r-space grids:        0.0794221448
  Total charge density g-space grids:           0.0794221448

  Overlap energy of the core charge distribution:               0.00000001844485
  Self energy of the core charge distribution:                -43.83289054591484
  Core Hamiltonian energy:                                     12.37678636914670
  Hartree energy:                                              21.01357980336000
  Exchange-correlation energy:                                 -4.00373216345181

  Total energy:                                               -14.44625651841509

 *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***


 -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 01        FRAGMENT SUBCONF: 01     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      4                                   O                  -
 -                      5                                   H                  -
 -                      6                                   H                  -
 -------------------------------------------------------------------------------

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   2
                             - Atoms:                                          3
                             - Shell sets:                                     6
                             - Shells:                                        11
                             - Primitive Cartesian functions:                 15
                             - Cartesian basis functions:                     24
                             - Spherical basis functions:                     23

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 Number of electrons:                                                          8
 Number of occupied orbitals:                                                  4
 Number of molecular orbitals:                                                 4

 Number of orbital functions:                                                 23
 Number of independent orbital functions:                                     23

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 P_Mix/Diag. 0.40E+00    0.0     1.14232298       -17.7900734306 -1.78E+01
     2 P_Mix/Diag. 0.40E+00    0.0     0.84069166       -17.8687264272 -7.87E-02
     3 P_Mix/Diag. 0.40E+00    0.0     0.49328190       -17.9175078664 -4.88E-02

  Leaving inner SCF loop after reaching     3 steps.


  Electronic density on regular grids:         -8.0135999514       -0.0135999514
  Core density on regular grids:                7.9991039523       -0.0008960477
  Total charge density on r-space grids:       -0.0144959991
  Total charge density g-space grids:          -0.0144959991

  Overlap energy of the core charge distribution:               0.00000001844463
  Self energy of the core charge distribution:                -43.83289054591484
  Core Hamiltonian energy:                                     13.08197907885014
  Hartree energy:                                              16.94243485590102
  Exchange-correlation energy:                                 -4.10903127368106

  Total energy:                                               -17.91750786640011

 *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***


 -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 11        FRAGMENT SUBCONF: 10     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      1                                   O                  -
 -                      2                                   H                  -
 -                      3                                   H                  -
 -                      4                                   O_ghost            -
 -                      5                                   H_ghost            -
 -                      6                                   H_ghost            -
 -------------------------------------------------------------------------------

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   4
                             - Atoms:                                          6
                             - Shell sets:                                    12
                             - Shells:                                        22
                             - Primitive Cartesian functions:                 30
                             - Cartesian basis functions:                     48
                             - Spherical basis functions:                     46

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 Number of electrons:                                                          8
 Number of occupied orbitals:                                                  4
 Number of molecular orbitals:                                                 4

 Number of orbital functions:                                                 46
 Number of independent orbital functions:                                     46

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 P_Mix/Diag. 0.40E+00    0.0     1.25452389       -14.3100140078 -1.43E+01
     2 P_Mix/Diag. 0.40E+00    0.0     0.72504403       -14.3964077769 -8.64E-02
     3 P_Mix/Diag. 0.40E+00    0.0     0.43185434       -14.4463330563 -4.99E-02

  Leaving inner SCF loop after reaching     3 steps.


  Electronic density on regular grids:         -7.9810289552        0.0189710448
  Core density on regular grids:                8.0604493642        0.0604493642
  Total charge density on r-space grids:        0.0794204090
  Total charge density g-space grids:           0.0794204090

  Overlap energy of the core charge distribution:               0.00000001844485
  Self energy of the core charge distribution:                -43.83289054591484
  Core Hamiltonian energy:                                     12.37552147261268
  Hartree energy:                                              21.01443216690922
  Exchange-correlation energy:                                 -4.00339616839607

  Total energy:                                               -14.44633305634415

 *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***


 -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 11        FRAGMENT SUBCONF: 01     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      1                                   O_ghost            -
 -                      2                                   H_ghost            -
 -                      3                                   H_ghost            -
 -                      4                                   O                  -
 -                      5                                   H                  -
 -                      6                                   H                  -
 -------------------------------------------------------------------------------

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   4
                             - Atoms:                                          6
                             - Shell sets:                                    12
                             - Shells:                                        22
                             - Primitive Cartesian functions:                 30
                             - Cartesian basis functions:                     48
                             - Spherical basis functions:                     46

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 Number of electrons:                                                          8
 Number of occupied orbitals:                                                  4
 Number of molecular orbitals:                                                 4

 Number of orbital functions:                                                 46
 Number of independent orbital functions:                                     46

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 P_Mix/Diag. 0.40E+00    0.0     1.14395715       -17.7900734306 -1.78E+01
     2 P_Mix/Diag. 0.40E+00    0.0     0.84171158       -17.8688061495 -7.87E-02
     3 P_Mix/Diag. 0.40E+00    0.0     0.49387470       -17.9176806499 -4.89E-02

  Leaving inner SCF loop after reaching     3 steps.


  Electronic density on regular grids:         -8.0136045107       -0.0136045107
  Core density on regular grids:                7.9991039523       -0.0008960477
  Total charge density on r-space grids:       -0.0145005584
  Total charge density g-space grids:          -0.0145005584

  Overlap energy of the core charge distribution:               0.00000001844463
  Self energy of the core charge distribution:                -43.83289054591484
  Core Hamiltonian energy:                                     13.07931650005626
  Hartree energy:                                              16.94402834918780
  Exchange-correlation energy:                                 -4.10813497164639

  Total energy:                                               -17.91768064987254

 *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***


 -------------------------------------------------------------------------------
 -                                                                             -
 -  BSSE CALCULATION         FRAGMENT CONF: 11        FRAGMENT SUBCONF: 11     -
 -                           CHARGE =     0           MULTIPLICITY =     0     -
 -                                                                             -
 -                 ATOM INDEX                              ATOM NAME           -
 -                 ----------                              ---------           -
 -                      1                                   O                  -
 -                      2                                   H                  -
 -                      3                                   H                  -
 -                      4                                   O                  -
 -                      5                                   H                  -
 -                      6                                   H                  -
 -------------------------------------------------------------------------------

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


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   2
                             - Atoms:                                          6
                             - Shell sets:                                    12
                             - Shells:                                        22
                             - Primitive Cartesian functions:                 30
                             - Cartesian basis functions:                     48
                             - Spherical basis functions:                     46

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          2
                             - Non-local part of the GTH pseudopotential:      0


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                               3
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-04
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        --------------------------------------------------------
                        Mixing method:                           DIRECT_P_MIXING
                        --------------------------------------------------------
                        No outer SCF

 Number of electrons:                                                         16
 Number of occupied orbitals:                                                  8
 Number of molecular orbitals:                                                 8

 Number of orbital functions:                                                 46
 Number of independent orbital functions:                                     46

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 P_Mix/Diag. 0.40E+00    0.0     1.25324632       -32.1031446416 -3.21E+01
     2 P_Mix/Diag. 0.40E+00    0.0     0.84543275       -32.2667100744 -1.64E-01
     3 P_Mix/Diag. 0.40E+00    0.0     0.49687917       -32.3642347919 -9.75E-02

  Leaving inner SCF loop after reaching     3 steps.


  Electronic density on regular grids:        -15.9946919530        0.0053080470
  Core density on regular grids:               16.0595533165        0.0595533165
  Total charge density on r-space grids:        0.0648613635
  Total charge density g-space grids:           0.0648613635

  Overlap energy of the core charge distribution:               0.00000003688949
  Self energy of the core charge distribution:                -87.66578109182967
  Core Hamiltonian energy:                                     25.44280098284285
  Hartree energy:                                              37.97079241838002
  Exchange-correlation energy:                                 -8.11204713815893

  Total energy:                                               -32.36423479187624

 *** WARNING in qs_scf.F:542 :: SCF run NOT converged ***


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                 BSSE RESULTS                                -
 -                                                                             -
 -                 CP-corrected Total energy:      -32.363985                  -
 -                                                                             -
 -                       1-body contribution:      -14.446257                  -
 -                       1-body contribution:      -17.917508                  -
 -                                                                             -
 -                       2-body contribution:       -0.000221                  -
 -                 BSSE-free interaction energy:       -0.000221               -
 -------------------------------------------------------------------------------

 -------------------------------------------------------------------------------
 -                                                                             -
 -                                DBCSR STATISTICS                             -
 -                                                                             -
 -------------------------------------------------------------------------------
 COUNTER                                    TOTAL       BLAS       SMM       ACC
 flops     5 x     5 x     8                12000     100.0%      0.0%      0.0%
 flops    13 x     5 x     8                12480     100.0%      0.0%      0.0%
 flops     5 x    13 x     8                12480     100.0%      0.0%      0.0%
 flops    13 x     5 x     4                15600     100.0%      0.0%      0.0%
 flops     5 x    13 x     4                15600     100.0%      0.0%      0.0%
 flops     5 x     5 x     4                15600     100.0%      0.0%      0.0%
 flops    13 x    13 x     8                24336     100.0%      0.0%      0.0%
 flops    13 x    13 x     4                32448     100.0%      0.0%      0.0%
 flops inhomo. stacks                           0       0.0%      0.0%      0.0%
 flops total                       140.544000E+03     100.0%      0.0%      0.0%
 flops max/rank                    140.544000E+03     100.0%      0.0%      0.0%
 matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%
 matmuls total                                225     100.0%      0.0%      0.0%
 number of processed stacks                    60     100.0%      0.0%      0.0%
 average stack size                                     3.8       0.0       0.0
 marketing flops                   228.528000E+03
 -------------------------------------------------------------------------------
 # multiplications                             15
 max memory usage/rank              72.966144E+06
 # max total images/rank                        1
 # max 3D layers                                1
 # MPI messages exchanged                       0
 MPI messages size (bytes):
  total size                         0.000000E+00
  min size                           0.000000E+00
  max size                           0.000000E+00
  average size                       0.000000E+00
 MPI breakdown and total messages size (bytes):
             size <=      128                   0                        0
       128 < size <=     8192                   0                        0
      8192 < size <=    32768                   0                        0
     32768 < size <=   131072                   0                        0
    131072 < size <=  4194304                   0                        0
   4194304 < size <= 16777216                   0                        0
  16777216 < size                               0                        0
 -------------------------------------------------------------------------------

 MEMORY| Estimated peak process memory [MiB]                                  70

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
 CP2K                                 1  1.0    0.007    0.007    0.766    0.766
 qs_energies                          5  2.0    0.000    0.000    0.426    0.426
 scf_env_do_scf                       5  3.0    0.000    0.000    0.348    0.348
 scf_env_do_scf_inner_loop           15  4.0    0.001    0.001    0.348    0.348
 qs_ks_update_qs_env                 15  5.0    0.000    0.000    0.208    0.208
 rebuild_ks_matrix                   15  6.0    0.000    0.000    0.207    0.207
 qs_ks_build_kohn_sham_matrix        15  7.0    0.001    0.001    0.207    0.207
 pw_transfer                        210  9.0    0.004    0.004    0.204    0.204
 fft_wrap_pw1pw2                    180  9.9    0.001    0.001    0.194    0.194
 fft_wrap_pw1pw2_30                  75 10.5    0.014    0.014    0.183    0.183
 create_qs_kind_set                   6  2.0    0.000    0.000    0.165    0.165
 read_qs_kind                        16  3.0    0.096    0.096    0.165    0.165
 qs_init_subsys                       6  2.0    0.002    0.002    0.165    0.165
 qs_rho_update_rho                   20  5.0    0.000    0.000    0.161    0.161
 calculate_rho_elec                  20  6.0    0.048    0.048    0.161    0.161
 qs_env_setup                         6  3.0    0.000    0.000    0.150    0.150
 qs_env_rebuild_pw_env               11  3.5    0.000    0.000    0.150    0.150
 pw_env_rebuild                       6  5.0    0.001    0.001    0.150    0.150
 fft3d_s                            181 11.9    0.125    0.125    0.130    0.130
 density_rs2pw                       20  7.0    0.001    0.001    0.112    0.112
 pw_grid_setup                       24  6.0    0.000    0.000    0.105    0.105
 pw_grid_setup_internal              24  7.0    0.007    0.007    0.105    0.105
 sum_up_and_integrate                15  8.0    0.001    0.001    0.100    0.100
 integrate_v_rspace                  15  9.0    0.051    0.051    0.099    0.099
 pw_grid_sort                        24  8.0    0.058    0.058    0.084    0.084
 parser_read_line                 37224  4.0    0.018    0.018    0.068    0.068
 init_scf_run                         5  3.0    0.000    0.000    0.052    0.052
 scf_env_initial_rho_setup            5  4.0    0.000    0.000    0.051    0.051
 parser_read_line_low                51  5.0    0.050    0.050    0.050    0.050
 potential_pw2rs                     15 10.0    0.000    0.000    0.047    0.047
 compute_max_radius                   6  6.0    0.042    0.042    0.042    0.042
 qs_vxc_create                       15  8.0    0.000    0.000    0.039    0.039
 xc_vxc_pw_create                    15  9.0    0.005    0.005    0.039    0.039
 xc_rho_set_and_dset_create          15 10.0    0.000    0.000    0.034    0.034
 xc_functional_eval                  15 11.0    0.032    0.032    0.032    0.032
 pw_gather_s                        100 11.3    0.030    0.030    0.030    0.030
 pw_poisson_solve                    15  8.0    0.015    0.015    0.026    0.026
 sort_shells                         24  9.0    0.026    0.026    0.026    0.026
 qs_energies_init_hamiltonians        5  3.0    0.000    0.000    0.025    0.025
 qs_env_update_s_mstruct              5  4.0    0.000    0.000    0.021    0.021
 qs_scf_new_mos                      15  5.0    0.000    0.000    0.021    0.021
 pw_scatter_s                        80 12.7    0.020    0.020    0.020    0.020
 calculate_rho_core                   5  5.0    0.002    0.002    0.019    0.019
 -------------------------------------------------------------------------------

 The number of warnings for this run is : 5
 
 -------------------------------------------------------------------------------
  **** **** ******  **  PROGRAM ENDED AT                 2020-10-10 10:56:36.337
 ***** ** ***  *** **   PROGRAM RAN ON                                  mt-rf614
 **    ****   ******    PROGRAM RAN BY                                     rf614
 ***** **    ** ** **   PROGRAM PROCESS ID                                186225
  **** **  *******  **  PROGRAM STOPPED IN /media/ssd1/rf614/Work/temp/jul_2020/
                                           cp2k_bsse
"""




