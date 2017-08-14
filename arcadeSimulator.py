#!/home/charizard/anaconda3/bin/python3
import numpy as np
from scipy import sparse
from arcadeUtils import displayer
import arcadePopulation as aP
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# NOTE : This is an script. Not the real object that will be used later
sns.set_style('whitegrid')
sns.set_palette('Set2')
np.random.seed(22)

# initialization

size_x = int(sys.argv[1])
size_y = int(sys.argv[2])
time   = int(sys.argv[3])
crops  = 1
disp = displayer(crops=crops, cropDays=time)
randomSeed = np.random.randint(0, high=size_x*size_y, size=2)

#parametersDict = dict(pathotypes  = 'P0|P12', transmision = dict(P0 = 1.0, P12 = 0.75))
parametersDict = dict(
	pathotypes = ['P0', 'P12'],
	beta = dict(
		P0 = 1.0,
		P12 = 0.75
	),
	n = dict(
		P0 = 15.0,
		P12 = 19.0
	),
	P = dict(
		P0 = 1000,
		P12 = 800
	),
	alpha = dict(
		P0 = 1.0 - 3.50e-2,
		P12 = 1.0 - 1.25e-2
	),
	C = np.array([
		[1.0, -0.45],
		[-0.19, 0.20]
	])
)

population = aP.arcadePopulation(size_x, size_y, parametersDict)
population.setSeed(randomSeed[0], patho='P0')
population.setSeed(randomSeed[1], patho='P12')


statistics = dict(
	P0 = np.empty(0, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('alive',      'float32'),
			('inoculum',  'float32')
		]),
	P12 = np.empty(0, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('alive',      'float32'),
			('inoculum',  'float32')
		])
)

# building interaction matrix

rx, ry = population.getCoordinates()
I = sparse.lil_matrix((size_x*size_y, size_x*size_y))
S = sparse.lil_matrix((size_x*size_y, size_x*size_y))
S_precalc = np.zeros(size_x*size_y)
iZ = 0

for ix in np.arange(size_x):
	for iy in np.arange(size_y):



		x_coordinates = np.arange(size_x, dtype='float32')
		x_coordinates -= ix
		y_coordinates = np.arange(size_y, dtype='float32')
		y_coordinates -= iy

		xx, yy = np.meshgrid(x_coordinates, y_coordinates)
		dist = np.sqrt((xx ** 2) + (yy ** 2))
		## dc holds for direct contact
		## sc holds for stochastic contact
		dc = np.exp(-1.25*dist)
		sc = np.exp(-0.15*dist)

		dc[dc < 1e-2] = 0.0
		sc[sc < 4e-1] = 0.0

		I[iZ,:] = dc.T.reshape((size_x * size_y))
		I[iZ,iZ] = 0

		S[iZ, :] = 1e-2 * sc.T.reshape((size_x * size_y))
		S[iZ, iZ] = 0

		iZ += 1


# simulation body
disp.start()
for crop in range(crops):
	for i in range(time):
		disp.update()
		population.updateAlive()
		population.primaryInfection()
		if i % 1 == 0:
			S_precalc = S.dot(population.getTranmission())
		population.updateExposition(I.dot(population.getTranmission()), S_precalc)
		population.updateInfective()
		population.updateInoculum()
		statistics['P0'] = np.concatenate((statistics['P0'],population.getStatistics('P0')))
		statistics['P12'] = np.concatenate((statistics['P12'],population.getStatistics('P12')))
		if i % 30 == 0:
			population.plotPopulation('P0', time=i, show=False)
			#population.plotPopulation('P12', time=i, show=False)

			#population.plotInfectiveAreas('P0', 'P12')
		#if i % 50 == 0:
			#population.plotInfectiveAreas('P0', 'P12')
			pass
	for j in range(365 - time):
		population.updateInoculum()
	population.nextCrop()

plt.plot(statistics['P0']['infective'], label="P0")
plt.plot(statistics['P12']['infective'], label="P12")

plt.xlabel('Days')
plt.ylabel('Value')
plt.legend(loc=0)
plt.show()

