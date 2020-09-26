import os
import tempfile
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.parseOther.parse_cp2k_files as tCode


class testParseInputFile(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		expLattVectsA = [ [5.80999959 , 0.00000000, 0.00000000],
		                 [-2.90495161, 5.03155172, 0.00000000],
		                 [0.00000000 , 0.00000000, 9.43314175] ]
		expFractCoordsA = [ [0.0, 0.0, 0.0, "Mg"],
		                    [1/3, 2/3, 0.5, "Mg"] ]
		self.expCellA = uCellHelp.UnitCell.fromLattVects(expLattVectsA, fractCoords=expFractCoordsA)
		self.fileStrA = _loadInputFileA()


	def testExpectedGeometry(self):

		with tempfile.NamedTemporaryFile('w+t', suffix='.csv', prefix=os.path.basename(__file__)) as f:
			f.write(self.fileStrA)
			f.seek(0)
			currDir = os.path.dirname(f.name)
			currPath = os.path.join( currDir, f.name )
			actGeom = tCode.parseGeomFromCpInputFile(currPath)
			f.close()
		
		self.assertEqual(self.expCellA, actGeom)


def _loadInputFileA():
	return """
&GLOBAL
  PROJECT_NAME vol_275pt762
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &COORD
      Mg 0.0 0.0 0.0
      Mg 0.333333 0.666667 0.5
      SCALED TRUE
    &END COORD
    &CELL
      B [bohr] -2.90495161 5.03155172 0.00000000
      A [bohr] 5.80999959 0.00000000 0.00000000
      C [bohr] 0.00000000 0.00000000 9.43314175
    &END CELL
    &KIND Mg
      POTENTIAL GTH-PBE-q2
      ELEMENT Mg
      BASIS_SET TZV2P-MOLOPT-SR-GTH-q2
    &END KIND
  &END SUBSYS
  &DFT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UCL
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &MGRID
      REL_CUTOFF [eV] 650
      CUTOFF [eV] 800
      NGRIDS 4
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      EPS_SCF 1.0E-7
      MAX_SCF 300
      ADDED_MOS 4
      SCF_GUESS ATOMIC
      &MIXING T
        ALPHA 0.4
        METHOD BROYDEN_MIXING
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 157.9
      &END SMEAR
      &DIAGONALIZATION ON
        ALGORITHM Standard
      &END DIAGONALIZATION
    &END SCF
    &KPOINTS
      SCHEME MONKHORST-PACK 20 20 12
    &END KPOINTS
  &END DFT
  &PRINT
    &FORCES On
    &END FORCES
  &END PRINT
&END FORCE_EVAL
"""

