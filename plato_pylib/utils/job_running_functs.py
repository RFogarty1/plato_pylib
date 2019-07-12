#!/usr/bin/python3

import itertools
import pathlib
import os
import shutil
import subprocess
from concurrent import futures

import plato_pylib.plato.mod_plato_inp_files as platoInp

def executeRunCommsParralel(allRunComms:list, maxThreads:int, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	
	quietMode = kwargs.get("quiet",None)

	jobsLeft = len(allRunComms)
	if quietMode is None:
		print("A total of {} jobs will be run".format(jobsLeft))

	numbThreads = min(jobsLeft, maxThreads)
	with futures.ThreadPoolExecutor(numbThreads) as executor:
		for i in executor.map(_runOneComm, allRunComms):
			jobsLeft -= 1
			if quietMode is None:
				print("{} Jobs remaining".format(jobsLeft))
	


def _runOneComm(runComm:str):
	subprocess.check_call(runComm, shell=True)


def removeTb1OutputFilesFromFolder(folder):
	outExtensions = [".occ", ".out", ".ham"]
	allFilePaths = [os.path.join(folder,x) for x in os.listdir(folder)]
	outFilePaths = list()
	for currPath in allFilePaths:
		if any([currPath.endswith(y) for y in outExtensions]):
			outFilePaths.append(currPath)
	[os.remove(x) for x in outFilePaths]


def pathListToPlatoRunComms(pathList: list, platoComm: str)->list:
	runComms = list()

	for path in pathList:
		currFolder, currFile = os.path.split(path)
		currRunComm = "cd {:};{:} {:} > outFile".format(currFolder, platoComm, currFile.replace('.in',''))
		runComms.append(currRunComm)

	return runComms



#---------------->Inverse SK specific <------------------

def runInvSkParralel(inpFilePaths,nCores):
	actStartDir = os.getcwd()
	startDirs = list()
	absFilePaths = [os.path.abspath(x) for x in inpFilePaths]
	for inpPath in absFilePaths:
		currStart = os.path.split(inpPath)[0]
		startDirs.append(currStart)

	jobsLeft = len(absFilePaths)
	print("A total of {} jobs will be run".format(jobsLeft))
	numbThreads = min( jobsLeft, nCores )

	with futures.ThreadPoolExecutor(numbThreads) as executor:
		for i in executor.map(runInvSk, absFilePaths, startDirs):
			jobsLeft -= 1
			print("jobsLeft = {}".format(jobsLeft))

	os.chdir(actStartDir)


def runInvSk(inpPath, startDir=None):
	if startDir is None:
		startDir = os.getcwd()

	#Step 1 = get base file name
	folder, fullFName = os.path.split(inpPath)
	baseFName, unused = os.path.splitext(fullFName)
	if folder=="":
		folder = startDir

	#Step 2 = create a temporary folder + copy file over
	tempDir = os.path.join(folder, baseFName)
	runFilePath = os.path.join(tempDir,fullFName)
	pathlib.Path(tempDir).mkdir(exist_ok=True)
	shutil.copy2(inpPath, runFilePath)

	#Step 3 = modify inversesk flag in file [to make sure its set] 
	_modInvSkFlagToOn(runFilePath)

	#Step 4 = Figure out what the output inv-sk filepaths will be
	outputInvSkFileNames = _getOutputInvSkFileNames(runFilePath)

	#Step 5 = run the actual inv-sk calculation [Mod the file to use inverse-SK if not already present]
	with ChDir(startDir, tempDir):
		subprocess.check_call("cd {};dft2 {}".format(tempDir,baseFName),shell=True) #Often breaks when the cd {} is missing

	#Step 6 = rename the relevant inv-sk output files to something unique
	for x in outputInvSkFileNames:
		filePath = os.path.join(tempDir,x)
		outPath = os.path.join(tempDir, baseFName + "_" + x)
		os.rename( filePath, outPath )

	#Step 7 = move back all files except the input file (NOT copy)
	for fName in os.listdir(tempDir):
		if fName.endswith(".in"):
			os.remove( os.path.join(tempDir,fName) )
		else:
			fPath = os.path.join(tempDir,fName)
			os.rename(fPath, os.path.join(folder,fName))

	#Step 8 = delete temporary directory
	os.rmdir(tempDir)


def _modInvSkFlagToOn(inpFilePath):
	tokenizedFile = platoInp.tokenizePlatoInpFile(inpFilePath)
	outOptDict = {k.lower():v for k,v in tokenizedFile.items()}
	outOptDict["inversesk"] = "1"
	platoInp.writePlatoOutFileFromDict(inpFilePath,outOptDict)

def _getOutputInvSkFileNames(inpFilePath):
	tokenizedFile = {k.lower():v for k,v in platoInp.tokenizePlatoInpFile(inpFilePath).items()}
	if int(tokenizedFile["format"]) != 0:
		raise ValueError("Geometry format needs to be {}; but format={}".format(0,tokenizedFile["format"]))	

	#Get the elements present
	nAtoms = tokenizedFile["natom"]
	atoms = tokenizedFile["atoms"].split("\n")
	elementsPresent = list()
	for x in atoms:
		currElement = x.split()[-1]
		if currElement not in (elementsPresent):
			elementsPresent.append(currElement)

	#Get the file names from them
	outFileNames = [ "{}_{}_SK.csv".format(a,b) for a,b in itertools.product(elementsPresent,elementsPresent)   ]

	return outFileNames

class ChDir:
	"""
	Step into a directory temporarily. Propagates exceptions upwards (since __exit__ always returns None)
	"""
	def __init__(self, startDir, workDir):
		self.old_dir = startDir
		self.new_dir = workDir
 
	def __enter__(self):
		os.chdir(self.new_dir)
 
	def __exit__(self, errorType, value, traceback):
		os.chdir(self.old_dir)
		if errorType is not None: #An error has occured
			return False
