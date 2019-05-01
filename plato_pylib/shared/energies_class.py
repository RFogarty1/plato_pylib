#!/usr/bin/python3


class EnergyVals():
	def __init__(self, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.e0 = kwargs.get("e0", None)
		self.e1 = kwargs.get("e1", None)
		self.e2 = kwargs.get("e2", None)
		self.entropy = kwargs.get("entropy", None)
		self.tb2TotalElectronic = kwargs.get("tb2TotalElectronic".lower(), None)
		self.tb2CohesiveFree = kwargs.get("tb2CohesiveFree".lower(), None) 
		self.dftTotalElectronic = kwargs.get("dftTotalElectronic".lower(),None)
		self.atomEnergy = kwargs.get("atomEnergy".lower(), None)
		self.castepTotalElectronic = kwargs.get("castepTotalElectronic".lower(), None)

	@property
	def electronicCohesiveE(self):
		if self.tb2CohesiveFree is not None:
			return self.tb2CohesiveFree + self.entropy
		elif (self.e0 is not None) and (self.e1 is not None): #In this case the file is Tb1. Meaning these values are relative to the atomic ones
			return self.e0 + self.e1
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
			return self.e0 + self.e1
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

