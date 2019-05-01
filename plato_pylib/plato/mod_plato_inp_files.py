#!/usr/bin/python3

import itertools as it
import os
import math

def replaceStrInFile(inpFilePath: str, inpStr: str, replaceStr: str):
	with open(inpFilePath, "rt") as f:
		fileStr = f.read()
	fileStr = fileStr.replace(inpStr, replaceStr)
	with open(inpFilePath, "wt") as f:
		f.write(fileStr)


def tokenizePlatoInpFile(inpFilePath:str)->dict:
	with open(inpFilePath, "rt") as f:
		fileList = f.readlines()
	platoInpOptions = dict()

	lineIdx = 0
	while lineIdx < len(fileList):
		currLine = fileList[lineIdx]
		if currLine.find('#') != -1:
			pass
		elif currLine.strip()=='':
			pass
		else:
			currKey = currLine.strip().lower()
			platoInpOptions[currKey], lineIdx = parseSingleOption(fileList, lineIdx+1)
		lineIdx += 1

	return platoInpOptions


def parseSingleOption(inpFileList:list ,lineIdx: int)->"str,int":
	optStr = ""
	currLine = inpFileList[lineIdx]
	while (currLine.strip()!= "") and (lineIdx<len(inpFileList)):
		optStr += inpFileList[lineIdx].strip()+'\n'
		if lineIdx+1 < len(inpFileList):
			lineIdx += 1
			currLine = inpFileList[lineIdx]
		else:
			break

	optStr = optStr.strip()

	return optStr, lineIdx


def modPlatoInpOptions(inpFilePath:str, optValDict: dict):
	with open(inpFilePath, "rt") as f:
		fileList = f.readlines()

	lineIdx = 0
	optValDictLower = {k.lower():v for k,v in optValDict.items()}


	while lineIdx < len(fileList):
		currLine = fileList[lineIdx].strip().lower()
		if currLine in optValDictLower.keys():
			lineIdx += 1
			modSingleInpOption(fileList, optValDictLower[currLine], lineIdx)
		lineIdx += 1

	with open(inpFilePath,"wt") as f:
		f.writelines(fileList)
	


def modSingleInpOption(inpFileList:"Str from single file in list format", modVal, lineIdx):
	currLine  = inpFileList[lineIdx]
	while (currLine.strip()!= "") and (lineIdx<len(inpFileList)):
		inpFileList.pop(lineIdx)
		if lineIdx+1 < len(inpFileList):
			currLine = inpFileList[lineIdx]
		else:
			break
	inpFileList[lineIdx] = str(modVal) +"\n\n"

	return lineIdx

def writePlatoOutFileFromDict(filePath, optValDict:"dict, all keys/vals should be str"):
	with open(filePath,"wt") as f:
		for key in optValDict.keys():
			f.write(key + "\n" + optValDict[key] + "\n\n")



def getPlatoGeomDictFromUnitCell(uCell:"UnitCell object"):

	#Get all info from unit-cell
	nAtoms = len(uCell.fractCoords)
	lattVects = uCell.lattVects
	fractCoords = uCell.fractCoords

	print("fractCoords = {}".format(fractCoords))
	#Convert info into nice format for plato
	cellSizes = [_getMagVector(x) for x in lattVects]
	unitVects = list()
	uVectsStr = ""
	fractCoordStr = ""

	for mag,vect in it.zip_longest(cellSizes,lattVects):
		currVect = [x/mag for x in vect]
		unitVects.append(currVect)
		uVectsStr += " ".join(["{:.8f}".format(x) for x in currVect]) + '\n'

	for currCoords in fractCoords:
		currStr = "{:.8f} {:.8f} {:.8f} {}".format(*currCoords)
		fractCoordStr += currStr + '\n'

	#Create dict with all relevant kwargs/strings
	outDict = dict()
	outDict["natom"] = str(nAtoms)
	outDict["cellsize"] = " ".join(["{:.8f}".format(x) for x in cellSizes])
	outDict["cellvec"] = uVectsStr
	outDict["atoms"] = fractCoordStr
	outDict["format"] = "0"

	return outDict

def _getMagVector(inpVect:iter):
	return math.sqrt( sum( [x**2 for x in inpVect] ) )


#Functions dealing with paths to Plato data files
def getDataFolderPathForPlatoInpFile():
    fullPath = getTightBindingDataPath()
    basePath = getPlatoRcPath()
    return os.path.relpath(fullPath,basePath)

def getTightBindingDataPath():
	rcPath = getPlatoRcPath()
	tbExt = os.path.join("Data","TightBinding")
	return os.path.join(rcPath, tbExt)

def getDftDataPath():
	rcPath = getPlatoRcPath()
	dftExt = os.path.join("Data","SCF_LCAO")
	return os.path.join(rcPath,dftExt)


def getPlatoRcPath():
	rcFile = os.path.join( os.path.expanduser('~'), ".platorc" )
	with open(rcFile,"rt") as f:
		rcPath = f.read()
	return rcPath.strip()

