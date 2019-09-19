# File to hold default options for input to plato. This SHOULD NOT be imported

#from plato_pylib.plato.mod_plato_inp_files import registerPlatoStrToOptDictFuntction

import copy

PLATO_PROG_STRS = set() #Some defaults are loaded at the bottom of this file
PROG_STR_TO_DEF_DICT_OBJ = dict() #Objects constructed at bottom on the file


def registerPlatoStrToDefDictObj(platoStr):
	def decorate(factoryFunct):
		PLATO_PROG_STRS.add(platoStr)
		PROG_STR_TO_DEF_DICT_OBJ[platoStr] = factoryFunct()
		return factoryFunct
	return decorate


def loadDefaultTb1OptDict():
	outDict = _loadDefaultDictAllPlatoProgs()
	
	outDict["CrystalFieldXCWeight"] = 1.0 #Full xtal field
	outDict["VxcMBCorrXtal"] = (0, 1.0) #first=flag, 2nd=weight attached to correction term
	outDict["VanDerWaals"] = 0

	outDict["addCorrectingPPFromBdt"] = 0
	outDict["addCorrectingHopFromBdt"] = 0

	#Want everything lower-case
	outDict = {k.lower():v for k,v in outDict.items()}
	return outDict


def loadDefaultTb2OptDict():
	outDict = _loadDefaultDictAllPlatoProgs()
	outDict["InverseSK"] = 1
	outDict["McWedaXcFlag"] = 0
	outDict["e0Method"] = 0
	outDict["xtalxcmethod"] = 0
	outDict["hopxcmethod"] = 0
	outDict["xtalVnaMethod"] = 0
	outDict["hopVnaMethod"] = 0
	outDict["xtalVnlMethod"] = 0
	outDict["hopVNLMethod"] = 0

	outDict["addCorrectingPPFromBdt"] = 0
	outDict["E0NonLocalXcCorr"] = 0
	outDict["e1XtalNonLocalXcCorr"] = 0
	outDict["e1HopNonLocalXcCorr"] = 0

	outDict = {k.lower():v for k,v in outDict.items()}
	return outDict



def getDefaultOptDictDft():
	outDict = _loadDefaultDictAllPlatoProgs()

	#DFT specific flags
	outDict["integralmesh"] = 0
	outDict["model"] = "lcao" #will correspond to 2
	outDict["optimizemesh"] = 0
	outDict["fftGridSpacing"] = 0.7
	outDict["ACGGridSpacing"] = [30,29,3]
	outDict["paralleliseFFT"] = 0
	outDict["SpinMixLevels"] = 5
	outDict["diagonalisationMethod"] = 1
	outDict["DiagonalisationWorkSize"] = 500
	outDict["MaxBond"] = 6.0
	outDict["DensityRWFlag"] = 0
	outDict["writeAtomDensityFlag"] = 0
	outDict["WavefunctionFlag"] = 0
	outDict["HyperfineFlag"] = 0

	#Parameters that should seldom be chnaged
	outDict["ACGEwaldParameters"] = [0.4,10.0]
	outDict["ACGMeshReduction"] = [6, 5.0]
	outDict["ACGMinPartitionWt"] = 0.0
	outDict["ACGAngularMeshType"] = 1
	outDict["DiameterNLV"] = 0.0
	outDict["IntegralPartitionFlag"] = 0
	outDict["DensityFit"] = [1,0]
	outDict["NumericVna"] = 0
	outDict["SpinMixScheme"] = 1
	outDict["MixThresHold"] = 100
	outDict["MixMetric"] = [0, 20.0]
	outDict["SplitDensityFlag"] = 0
	outDict["scfFlag"] = 1

	outDict = {k.lower():v for k,v in outDict.items()}

	return outDict



def _loadDefaultDictAllPlatoProgs():
	''' Any keyword used in >1 program should go here '''
	outDict = dict()

	#Electornic structure method parameters
	outDict["spinflag"] = 0
	outDict["LSDApUFlag"] = 0
	outDict["OccupyFlag"] = "fermi"
	outDict["OccupyTolerance"] = 1.0e-12
	outDict["electronTemperature"] = 0.001
	outDict["BlochStates"] = None
	outDict["ElectronExcess"] = 0.0
	outDict["MagneticField"] = 0.0
	outDict["ExternalField"] = 0

	#Scf parameters
	outDict["scfFlag"] = 0
	outDict["Nloops"] = 30
	outDict["ETol"] = 0.0
	outDict["ResidueTol"] = 0.001
	outDict["MixScheme"] = 1
	outDict["MixLevels"] = 5
	outDict["MixFactor"] = 0.2
	outDict["SpinMixFactor"] = 0.2
	outDict["PreconditioningFlag"] = [1,1]
	outDict["mixClipping"] = [0.01,0.1]

	#Hairy probes stuff - Currently not on by default since its not obvious to me (yet) what strings to write out by default
#	outDict["OpenBoundaryTerminals"] = None
#	outDict["OpenBoundaryCurrent"] = None
#	outDict["OpenBoundaryTransmission"] = None
	outDict["NBands"] = 0

	#"Characteristics of the job" parameters
	outDict["job"] = "spe"
	outDict["restartFlag"] = 0
	outDict["MaximumTime"] = 1000000.0 #Number of hours it can run
	outDict["NoNetForceFlag"] = 1

	#Geometry specific parameters
	outDict["cellRepeat"] = [1,1,1]
	outDict["pbcflag"] = 1

	#For relaxation simulations
	outDict["RFlag"] = 3 #just the algorithm for relaxation
	outDict["cellRelaxMode"] = "free"
	outDict["FTol"] = 0.0001
	outDict["MaxDisplacement"] = 0.1
	outDict["StressTol"] = 0.001
	outDict["Pinned"] = None
	outDict["Constrain"] = None
	outDict["ReactionCoordinate"] = None

	#For molecular dynamics
	outDict["atomTemperature"] = 300.0
	outDict["temperatureTolerance"] = 0.1
	outDict["NoseHoover"] = -1.0
	outDict["KickAtom"] = None
	outDict["MaxNSteps"] = 1000
	outDict["MDTimeStep"] = 1.0
	outDict["AutoSave"] = 10

	#Tb1 input/output stuff
	outDict["verbosity"] = 1
	outDict["movierate"] = 0
	outDict["writeEigenVects"] = 0 #Needs translating later
	outDict["PartialPopulationFlag"] = 0
	outDict["writeEnergyFlag"] = 0
	outDict["hamiltonianFlag"] = 1

	#"Parameters that should seldom be changed" 
	outDict["basisFlag"] = 0
	outDict["seed"] = 12345
	outDict["OnSiteElectrostatics"] = 0

	#TB2 specific parameters
	outDict["XCFunctional"] = "PBE"
	outDict["IntegralMeshType"] = "atom"
	outDict["IntegralMeshSpacing"] = [100, 30, 30]
	outDict["testOverlap"] = 1
	outDict["excMbCorr"] = 0

	outDict = {k.lower():v for k,v in outDict.items()}

	return outDict


#Functions to convert opt dicts to string dicts
def getStrDictFromOptDict_tb1OrTb2(inpDict):
	optDict = copy.deepcopy(inpDict)
	_doPartialOptToStrDictAnyPlatoProg(optDict)
	_partialOptToStrDictTb2OnlyKeywords(optDict)

	#Convert any remaining non-strings. 
	for k in optDict.keys():
		optDict[k] = str(optDict[k])

	return optDict

def _partialOptToStrDictTb2OnlyKeywords(optDict):
	#flags -1 to 4 just take 1 value. flag=5 takes 3 additional lines
	try:
		optDict["xtalxcmethod"] = " ".join([str(x) for x in optDict["xtalxcmethod"]])
	except TypeError:
		optDict["xtalxcmethod"] = str(optDict["xtalxcmethod"])
	except KeyError:
		pass

def getStrDictFromOptDict_dft(inpDict):
	optDict = copy.deepcopy(inpDict)
	_doPartialOptToStrDictAnyPlatoProg(optDict)

	modelToFlag = {"lcao":"2"}

	if "XCFunctional".lower() in inpDict.keys():
		optDict["XCFlag".lower()] = optDict["XCFunctional".lower()]

	if "model" in inpDict.keys():
		optDict["model"] = modelToFlag[ optDict["model"] ]

	if "ACGGridSpacing".lower() in inpDict.keys():
		optDict["ACGGridSpacing".lower()] = " ".join([str(x) for x in optDict["ACGGridSpacing".lower()] ])

	if "ACGEwaldParameters".lower() in inpDict.keys():
		optDict["ACGEwaldParameters".lower()] = " ".join([str(x) for x in optDict["ACGEwaldParameters".lower()] ])

	if "ACGMeshReduction".lower() in inpDict.keys():
		optDict["ACGMeshReduction".lower()] = " ".join([str(x)for x in optDict["ACGMeshReduction".lower()] ])

	if "DensityFit".lower() in inpDict.keys():
		optDict["DensityFit".lower()] = " ".join([str(x) for x in optDict["DensityFit".lower() ]])

	if "MixMetric".lower() in inpDict.keys():
		optDict["MixMetric".lower()] = " ".join([str(x) for x in optDict["MixMetric".lower() ]])

	#Convert any remaining non-strings
	for k in optDict.keys():
		optDict[k] = str(optDict[k])


	optDict["integralmesh"] = str(inpDict["integralmesh"])
	#"XCFlag" keyword used for the functional in dft
	#IntegralMesh is same keyword but different fields
	return optDict


def _doPartialOptToStrDictAnyPlatoProg(optDict):
	#Translate to plato Strings
	functToFlag = {"none":-1,"lda":0,"pbe":1}
	occupyToFlag = {"fermi":1}
	jobTypeToFlag = {"spe":1}
	cellRelaxToFlag = {"free":0}
	meshTypeToFlag = {"uniform":0,"atom":1}

	if "BlochStates".lower() in optDict.keys():
		k = "BlochStates".lower()
		if optDict[k] is None:
			optDict[k] = str(0)
		elif len(optDict[k]) == 3:
			optDict[k] = "-1\n" + " ".join([str(x) for x in optDict[k]])
		elif len(optDict[k]) == 2:
			blochStr = str( optDict[k][0] ) + "\n"
			for kPoint in optDict[k][1:]:
				blochStr += " ".join ([str(x) for x in kPoint]) + "\n"
			optDict[k] = blochStr
		else:
			raise ValueError("{} is an invalid value for blochstates".format( optDict[k] ))

	if "occupyflag" in optDict.keys():
		optDict["OccupyFlag".lower()] = occupyToFlag[ optDict["OccupyFlag".lower()] ]

	if "PreconditioningFlag".lower() in optDict.keys():
		optDict["PreconditioningFlag".lower()] = " ".join( [str(x) for x in optDict["PreconditioningFlag".lower()]] )

	if "mixClipping".lower() in optDict.keys():
		optDict["mixClipping".lower()] = " ".join( [str(x) for x in optDict["mixClipping".lower()]] )

	if "job" in optDict.keys():
		optDict["job"] = jobTypeToFlag[ optDict["job"] ]

	if "cellrepeat" in optDict.keys():
		optDict["cellrepeat"] = " ".join([str(x) for x in optDict["cellrepeat"]])

	if "cellRelaxMode".lower() in optDict.keys():
		optDict["cellRelaxMode".lower()] = cellRelaxToFlag[ optDict["cellRelaxMode".lower()] ]

	if "Pinned".lower() in optDict.keys():
		k = "Pinned".lower()
		if optDict[k] is None:
			optDict[k] = 0

	if "Constrain".lower() in optDict.keys():
		k = "Constrain".lower()
		if optDict[k] is None:
			optDict[k] = 0

	if "ReactionCoordinate".lower() in optDict.keys():
		k = "ReactionCoordinate".lower()
		if optDict[k] is None:
			optDict[k] = 0

	if "KickAtom".lower() in optDict.keys():
		k = "KickAtom".lower()
		if optDict[k] is None:
			optDict[k] = 0

	if "XCFunctional".lower() in optDict.keys():
		optDict["XCFunctional".lower()] = functToFlag[optDict["XCFunctional".lower()].lower()]

	if "IntegralMeshType".lower() in optDict.keys():
		meshSpacingStr = " ".join( [str(x) for x in optDict["IntegralMeshSpacing".lower()]] )
		meshTypeStr =  str( meshTypeToFlag[ optDict["IntegralMeshType".lower()] ] )
		optDict.pop("IntegralMeshSpacing".lower())
		optDict.pop("IntegralMeshType".lower())
		optDict["IntegralMesh".lower()] = meshTypeStr + "\n" + meshSpacingStr

	if "cellsize" in optDict.keys():
		optDict["cellsize"] = " ".join( ["{:.8f}".format(x) for x in optDict["cellsize"]] )

	if "atoms" in optDict.keys():
		optDict["atoms"] = " ".join( [str(x) for x in optDict["atoms"]] )

	if "cellvec" in optDict.keys():
		cellVecStr = ""
		for vect in optDict["cellvec"]:
			cellVecStr += " ".join( ["{:.8f}".format(x) for x in vect] ) + "\n"
		optDict["cellvec"] = cellVecStr

	if "VxcMBCorrXtal".lower() in optDict.keys():
		flag, fract = optDict["VxcMBCorrXtal".lower()]
		optDict["VxcMBCorrXtal".lower()] = "{} {:.8f}".format( int(flag), float(fract) ) 

	if "excmbcorr".lower() in optDict.keys():
		try:
			flag, fract = optDict["excMbCorr".lower()]
			optDict["excmbcorr"] = "{} {:.8f}".format( int(flag), float(fract) ) 
		except TypeError:
			pass #For dft2 we dont need to specify the fraction of xc at time of writing



class DefaultPlatoDicts():
	def __init__(self, optDictFunct, optDictToStrDictFunct):
		self.optDictFunct = optDictFunct
		self.optDictToStrFunct = optDictToStrDictFunct


@registerPlatoStrToDefDictObj("tb1")
def createTb1Obj():
	optDictFunct = loadDefaultTb1OptDict
	optDictToStrFunct = getStrDictFromOptDict_tb1OrTb2
	return DefaultPlatoDicts(optDictFunct, optDictToStrFunct)


@registerPlatoStrToDefDictObj("tb2")
@registerPlatoStrToDefDictObj("dft2")
def createTb2Obj():
	optDictFunct = loadDefaultTb2OptDict
	optDictToStrFunct = getStrDictFromOptDict_tb1OrTb2
	return DefaultPlatoDicts(optDictFunct, optDictToStrFunct)


@registerPlatoStrToDefDictObj("dft")
def createDftObj():
	optDictFunct = getDefaultOptDictDft
	optDictToStrFunct = getStrDictFromOptDict_dft
	return DefaultPlatoDicts(optDictFunct, optDictToStrFunct)

