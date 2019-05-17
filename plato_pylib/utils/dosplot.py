
import itertools as it
import math
import numpy as np


def genDosData(eigVals, occVals, **kwargs):
	''' For eigVals and occVals, each row is 1 k-pt, each column is one orbital '''
	dataOpts = _createDosDataOpts(eigVals, occVals, **kwargs)
	return _createDosData(dataOpts)



def _createDosDataOpts(eigVals, occVals, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	minX = kwargs.get("minx", None)
	maxX = kwargs.get("maxx", None)
	step = kwargs.get("step", 0.01)
	smearWidth = kwargs.get("smearwidth", 0.2)
	kPtWeights = kwargs.get("kptweights",None)

	minEig, maxEig = min(min(eigVals)), max(max(eigVals))

	if minX is None:
		minX = minEig - (5*smearWidth)
	if maxX is None:
		maxX = maxEig + (5*smearWidth)
	if kPtWeights is None:
		kPtWeights = [1/len(eigVals) for x in range(len(eigVals))]

	return DosDataOptions(eigVals, occVals, minX, maxX, smearWidth, step, kPtWeights)


def _createDosData(dataOpts):
	assert len(dataOpts.eigVals) == len(dataOpts.occVals)
	assert len(dataOpts.eigVals[0]) == len(dataOpts.occVals[0])

	eigVals = np.array( dataOpts.eigVals )
	occVals = np.array( dataOpts.occVals )
	assert eigVals.shape == occVals.shape

	#Generate the Gaussian functions
	allGauFuncts = list()
	for kPtIdx, unused in enumerate(dataOpts.eigVals):
		allGauFuncts.append( list() )
		for eigIdx, unused in enumerate(dataOpts.eigVals[kPtIdx]):
			currGauFunct = _createNormalisedGaussianFunct( eigVals[kPtIdx,eigIdx], dataOpts.smearWidth, coeff=occVals[kPtIdx,eigIdx]*dataOpts.kPtWeights[kPtIdx] ) 
			allGauFuncts[-1].append( currGauFunct )

	comboFunct = _combineAllGauFunctsIntoOne(allGauFuncts)

	#Generate the output values
	xVals = np.arange(dataOpts.minX, dataOpts.maxX, dataOpts.step)
	yVals = np.zeros( (xVals.shape) )
	for idx, unused in enumerate(yVals):
#		yVals[idx] = _getSumOfTwoDimGauFunctArray(xVals[idx],allGauFuncts)
		yVals[idx] = comboFunct(xVals[idx])

	return [(x,y) for x,y in it.zip_longest(xVals,yVals)]


def _createNormalisedGaussianFunct(centre,smearWidth,coeff=1.0):
	alpha = 1 / (2* (smearWidth**2) )
	normFactor =  math.sqrt(alpha / math.pi)
	cVal = normFactor*coeff
	gauFunct = _GauFunctOneDim(alpha, cVal, centre)
	return gauFunct


def _getSumOfTwoDimGauFunctArray(xVal, gFuncts):
	outVal = 0.0
	for idxA, unused in enumerate(gFuncts):
		for idxB, unused in enumerate(gFuncts[idxA]):
			outVal += gFuncts[idxA][idxB](xVal)
	return outVal


def _combineAllGauFunctsIntoOne( twoDimGaus ):
	#Flatten array + extract useful params
	centreVals = list()
	expVals = list()
	cVals = list() 
	for idxA, unused in enumerate(twoDimGaus):
		for idxB, unused in enumerate(twoDimGaus[idxA]):
			currGau = twoDimGaus[idxA][idxB]
			centreVals.append(currGau.centre)
			cVals.append(currGau.coeff)
			expVals.append(currGau.exponent)
	
	def combinedFunct(xVal):
		outVal = 0.0
		for r,a,c in zip(centreVals, expVals, cVals):
			sqrDiff = (r-xVal)**2
			outVal += c*math.exp( (-1*a *sqrDiff) )
		return outVal
	return combinedFunct

class DosDataOptions:
	def __init__(self, eigVals, occVals, minX, maxX, smearWidth, step, kPtWeights):
		self.eigVals = eigVals
		self.occVals = occVals
		self.minX = minX
		self.maxX = maxX
		self.smearWidth = smearWidth
		self.step = step
		self.kPtWeights = kPtWeights


class _GauFunctOneDim():
	def __init__(self, exponent, coeff, centre):
		self.exponent = exponent
		self.coeff = coeff
		self.centre = centre

	def __call__(self,xVal):
		expTerm = math.exp( (-1*self.exponent * ((xVal-self.centre)**2)) )
		return expTerm*self.coeff


