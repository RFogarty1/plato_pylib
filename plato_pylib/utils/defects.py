
import itertools as it
import math


def calcVacancyE(energyNoVac, energyVac, nAtomsOrig, nVacancies=1):
	return energyVac - ( ((nAtomsOrig-nVacancies)/nAtomsOrig) * energyNoVac )


def makeVacancyUnitCell(inpUCell:"UnitCell class object", method="central"):
	''' Makes vacancy on UnitCell object in place. method currently limited to central,
	    meaning the atom nearest the center of the unit cell is removed '''
	if method != "central":
		raise ValueError("method = {} is an invalid argument".format(method))

	outFractCoords = _getStructWithVacancyFractCoords(inpUCell.fractCoords, method)
	inpUCell.fractCoords = outFractCoords


def _getStructWithVacancyFractCoords(coordList, method):
	if method == "central":
		outCoords = _getStructWithVacancyFractCoords_central(coordList)
	else:
		raise ValueError("{} is an invalid option".format(method) )
	return outCoords

def _getStructWithVacancyFractCoords_central(coordList):
	onlyCoords = [x[:3] for x in coordList]
	centIdx = _getCentralIdxFractCoordList(onlyCoords)
	outCoords = [x for idx,x in enumerate(coordList) if idx!=centIdx]
	return outCoords


def _getCentralIdxFractCoordList(coordList:"nx3 iter, e.g [[x1,y1,z1],[x2,y2,z2]]"):
	center = [0.5,0.5,0.5]
	allDists = [_distanceTwoVects(x,center) for x in coordList]
	return allDists.index(min(allDists))


def _distanceTwoVects(vectA:"1 dim iter", vectB:"1 dim iter"):
	return math.sqrt(sum( [(a-b)**2 for a,b in it.zip_longest(vectA,vectB)] ))

