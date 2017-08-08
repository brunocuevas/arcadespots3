#!/home/charizard/anaconda3/bin/ipython3
import numpy as np
from scipy import sparse
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
iZ = 0

for ix in np.arange(size_x):
	for iy in np.arange(size_y):



		x_coordinates = np.arange(size_x, dtype='float32')
		x_coordinates -= ix
		y_coordinates = np.arange(size_y, dtype='float32')
		y_coordinates -= iy

		xx, yy = np.meshgrid(x_coordinates, y_coordinates)
		dist = np.sqrt((xx ** 2) + (yy ** 2))
		dist = np.exp(-1.50*dist)
		dist[dist < 1e-2] = 0.0

		I[iZ,:] = dist.T.reshape((size_x * size_y))
		I[iZ,iZ] = 0
		iZ += 1


# simulation body
for crop in range(1):
	for i in range(time):
		print("day %d" % i)
		population.updateAlive()
		population.primaryInfection()
		population.updateExposition(I.dot(population.getTranmission()))
		population.updateInfective()
		population.updateInoculum()
		statistics['P0'] = np.concatenate((statistics['P0'],population.getStatistics('P0')))
		statistics['P12'] = np.concatenate((statistics['P12'],population.getStatistics('P12')))
		#if i == 156:
		#	population.plotPopulation('P12')

	for i in range(365 - time):
		population.updateInoculum()
	population.nextCrop()

plt.plot(statistics['P0']['infective'], label="P0")
plt.plot(statistics['P12']['infective'], label="P12")

plt.xlabel('Days')
plt.ylabel('Value')
plt.legend(loc=0)
plt.show()

