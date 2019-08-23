#!/usr/bin/python3

import os
import unittest

import plato_pylib.utils.job_running_functs as tCode

class TestInvSkBashComms(unittest.TestCase):

	def setUp(self):
		self.fakeInpFolder = os.path.abspath( os.getcwd() )#Code will take the absPath + i'd prefer not to mock it
		self.fakeInpPath = os.path.abspath( os.path.join(self.fakeInpFolder,"fake_inp_file.in") ) 

	def runTestFunction(self):
		outVals = tCode.invSkInputPathsToBashComms([self.fakeInpPath])
		return outVals

	def testExpCommandsForSingleFile(self):
		actCommands = self.runTestFunction()
		firstExpCommand = "cd {}".format(self.fakeInpFolder)
		otherComms = "mkdir tempdir_fake_inp_file;cp fake_inp_file.in tempdir_fake_inp_file;cd tempdir_fake_inp_file;dft2 fake_inp_file;for file in *.csv;do mv $file fake_inp_file_$file;done;cp *.csv ..;cp *.out ..;rm *;cd ..;rmdir tempdir_fake_inp_file"
		expCommands = [firstExpCommand + ";" + otherComms]
		self.assertEqual(expCommands,actCommands)


if __name__ == '__main__':
	unittest.main()

