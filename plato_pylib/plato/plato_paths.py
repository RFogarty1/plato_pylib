
import os
#import plato_pylib.plato.mod_plato_inp_files as modInp


""" Module for manipulating paths related to plato; main original use is finding model data paths """



class PlatoPathDescriptor():
	
	def __init__(self,attrName, dataType="tightbinding"):
		self.attrName = attrName
		self.dType = dataType

	#Separate function for this (vs fix on initialisation) lets me mock out the calls
	#(and therefore test more easily)
	def _getStartPath(self):
		if self.dType == "tightbinding":
			startPath = getTightBindingDataPath()
		elif self.dType == "dft":
			startPath = getDftDataPath()
		else:
			raise ValueError("{} is an invalid option for dataType".format(self.dType))
		return startPath

	#TODO: Check this throws sensible error in the case of a bad rc path(check for .. in outPat)
	def __get__(self,instance,owner):
		startPath = self._getStartPath()
		absPath = os.path.abspath(getattr(instance,self.attrName))
		outPath = os.path.relpath(absPath,startPath)
		return outPath

	def __set__(self,val):
		raise NotImplementedError("Cannot set attribute {}".format(self.name))



class PlatoModelFolders():

	tb1PlatoPath = PlatoPathDescriptor("tb1Model")
	dft2PlatoPath = PlatoPathDescriptor("dft2Model")
	dftPlatoPath = PlatoPathDescriptor("dftModel", dataType="dft")

	def __init__(self, tb1Path=None, dft2Path=None, dftPath=None):
		""" 
		Args:
			tb1Path : Absolute path to folder containing model for tb1
			dft2Path: Absolute path to folder containing model for dft2
			dftPath : Absolute path to folder containing model for dft
		"""
		self._tb1Path = tb1Path
		self._dft2Path = dft2Path
		self._dftPath = dftPath

	@property
	def tb1Model(self):
		""" Full path to the FOLDER containing the model files to use for tb1 """
		if self._tb1Path is not None:
			return os.path.abspath(self._tb1Path)
		else:
			return None

	@tb1Model.setter
	def tb1Model(self,value):
		self._tb1Path = os.path.abspath(value)

	@property
	def dft2Model(self):
		""" Full path to the FOLDER containing the model files to use for dft2 """
		if self._dft2Path is not None:
			return os.path.abspath(self._dft2Path)
		else:
			return None

	@dft2Model.setter
	def dft2Model(self,value):
		self._dft2Path = os.path.abspath(value)

	@property
	def dftModel(self):
		""" Full path to the FOLDER containing the model files to use for dft """
		if self._dftPath is not None:
			return os.path.abspath(self._dftPath)
		else:
			return None

	@dftModel.setter
	def dftModel(self,value):
		self._dftPath = os.path.abspath(value)



#Low level functions dealing with model folder paths
def getAbsolutePathForPlatoTightBindingDataSet(dSetPath,dtype="tightbinding"):
	""" Takes the relative path in the dataset of a plato input file and converts to an absolute path
	to the relevant model data folder
	
	Args:
		dsetPath: The path you pass to dataset field of a plato input file
		dtype(optional): Whether tightbinding or dft model. The dft models are located in a different folder (SCF_LCAO)
	Returns
		dfolderAbsPath: Absolute path to the dataset that plato will use when given dsetPath as dataset
	
	"""
	if dtype=="tightbinding":
		basePath = getTightBindingDataPath()
	elif dtype=="dft":
		basePath = getDftDataPath()
	else:
		raise ValueError("{} is an invalid argument for dtype".format(dtype))
	return os.path.abspath(os.path.join(basePath,dSetPath))

def getTightBindingDataPath():
	""" Return absolute path to the folder plato uses to store tight binding models
	"""
	rcPath = getPlatoRcPath()
	tbExt = os.path.join("Data","TightBinding")
	return os.path.join(rcPath, tbExt)

def getDftDataPath():
	""" Return absolute path to the folder plato uses to store models for the dft code
	"""
	rcPath = getPlatoRcPath()
	dftExt = os.path.join("Data","SCF_LCAO")
	return os.path.join(rcPath,dftExt)



def getPlatoRcPath():
	""" Returns absolute path found in .platorc file (which links to where plato is located)	
	"""
	rcFile = os.path.join( os.path.expanduser('~'), ".platorc" )
	with open(rcFile,"rt") as f:
		rcPath = f.read()
	return rcPath.strip()

