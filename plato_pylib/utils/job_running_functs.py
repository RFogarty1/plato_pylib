#!/usr/bin/python3

import os
import subprocess
from concurrent import futures

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

