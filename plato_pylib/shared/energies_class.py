#!/usr/bin/python3


RYD_TO_EV = 13.6056980659

class EnergyVals():
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self._eqTol = 1e-5
		self._e0Tot = kwargs.get("e0tot", None)
		self._e0Coh = kwargs.get("e0coh", None)
		self.e1 = kwargs.get("e1", None)
		self.e2 = kwargs.get("e2", None)
		self.entropy = kwargs.get("entropy", None)
		self.tb2CohesiveFree = kwargs.get("tb2CohesiveFree".lower(), None) 
		self.dftTotalElectronic = kwargs.get("dftTotalElectronic".lower(),None)
		self.castepTotalElectronic = kwargs.get("castepTotalElectronic".lower(), None)
		self.dispersion = kwargs.get("dispersion",None)

		#These are numerical-attributes which form the base-data for the class. Keep them
		#here since they are used both to convert the object into a dictionary AND
		#to check equality between objects
		self.numAttrs = ["e0Tot", "e0Coh", "e1", "e2", "entropy", "tb2CohesiveFree",
		                 "dftTotalElectronic", "castepTotalElectronic", "dispersion"]

	def convRydToEv(self):
		for key in self.numAttrs:
			if getattr(self,key) is not None:
				setattr(self, key, getattr(self, key)*RYD_TO_EV)

	def toDict(self):
		extraAttrs = ["electronicTotalE", "electronicMinusEntropy", "electronicMinusHalfEntropy"]
		basicDict = {attr:getattr(self,attr) for attr in self.numAttrs if getattr(self,attr) is not None}
		extraAttrDict = {attr:getattr(self,attr) for attr in extraAttrs if getattr(self,attr) is not None}
		outDict = dict()
		outDict.update(basicDict)
		outDict.update(extraAttrDict)
		return outDict

	@property
	def e0Tot(self):
		return self._e0Tot

	@property
	def e0Coh(self):
		return self._e0Coh

	@property
	def electronicCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree + self.entropy
		elif (self.e0Coh is not None) and (self.e1 is not None): #In this case the file is Tb1. Meaning these values are relative to the atomic ones
			return self.e0Coh + self.e1
		else:
			raise ValueError("No information on electronic Cohesive Energy appears to be "
			                 "held in current EnergyVals object")

	@property
	def electronicTotalE(self):
		""" Note: Castep/CP2K (at least) include the entropy contribution here """
		if self.castepTotalElectronic is not None:
			return self.castepTotalElectronic
		if self.dftTotalElectronic is not None:
			return self.dftTotalElectronic
		elif self.tb2CohesiveFree is not None:
			return self.e0Tot + self.e1
		else:
			raise ValueError("No information on total electronic energy is present in current "
			                 "EnergyVals object")

	@property
	def electronicMinusEntropy(self):
		if self.entropy is None:
			outVal = None
		else:
			outVal = self.electronicTotalE - self.entropy
		return outVal

	@property
	def electronicMinusHalfEntropy(self):
		""" This is what castep uses to get 0K estimates """
		if self.entropy is None:
			outVal = None
		else:
			outVal = self.electronicTotalE - (0.5*self.entropy)
		return outVal

	@property
	def freeCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree
		else:
			return self.e0Coh + self.e1 + self.entropy
#		else:
#			raise ValueError("No information on free Cohesive Energy appears to be "
#			                 "held in current EnergyVals object")


	def __eq__(self, other):
		eqTol = min([abs(x._eqTol) for x in [self,other]])

		for currAttr in self.numAttrs:
			if not self._attrsAreTheSameOnOtherObj(currAttr, other, eqTol):
				return False

		return True

	def _attrsAreTheSameOnOtherObj(self, attr, other, eqTol):
		valA, valB = getattr(self, attr), getattr(other,attr)
		if (valA is None) and (valB is None):
			return True
		elif (valA is None) or (valB is None):
			return False

		if abs(valA-valB)>eqTol:
			return False
		return True
