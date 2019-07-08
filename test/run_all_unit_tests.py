#!/usr/bin/python3

import os
import sys
import unittest

def main():
	#Need to add certain paths to pythons search. This is purely becuase the u-tests sometimes import from each other
	startDir = os.path.split( os.path.abspath(os.getcwd(),) )[0]
	uTestDirs = [os.path.join(startDir,"plato_pylib","plato","unit_tests"),
	             os.path.join(startDir,"plato_pylib","parseOther","unit_tests")]
	for tDir in uTestDirs:
		sys.path.append(tDir)

	#Find all unit-tests and run them
	loader = unittest.TestLoader()
	suite = loader.discover(startDir, pattern='*utest*.py')
	runner = unittest.TextTestRunner(verbosity=2)
	result = runner.run(suite)

	sys.exit( int(not result.wasSuccessful()) ) #return 0 on success, else 1

if __name__ == '__main__':
	main()
