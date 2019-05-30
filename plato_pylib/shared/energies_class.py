#!/usr/bin/python3


RYD_TO_EV = 13.6056980659

class EnergyVals():
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self._e0Tot = kwargs.get("e0tot", None)
		self._e0Coh = kwargs.get("e0coh", None)
		self.e1 = kwargs.get("e1", None)
		self.e2 = kwargs.get("e2", None)
		self.entropy = kwargs.get("entropy", None)
		self.tb2CohesiveFree = kwargs.get("tb2CohesiveFree".lower(), None) 
		self.dftTotalElectronic = kwargs.get("dftTotalElectronic".lower(),None)
		self.castepTotalElectronic = kwargs.get("castepTotalElectronic".lower(), None)

	def convRydToEv(self):
		for key in self.__dict__:
			if getattr(self,key) is not None:
				setattr(self, key, getattr(self, key)*RYD_TO_EV)


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
	def freeCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree
		else:
			raise ValueError("No information on free Cohesive Energy appears to be "
			                 "held in current EnergyVals object")

