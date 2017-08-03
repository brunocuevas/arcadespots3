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
np.random.seed(54)

# initialization

size_x = int(sys.argv[1])
size_y = int(sys.argv[2])
time   = int(sys.argv[3])
randomSeed = np.random.randint(0, high=size_x*size_y, size=1)

population = aP.arcadePopulation(size_x, size_y)
population.setSeed(randomSeed)

statistics = np.empty(0, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('alive',      'float32'),
			('inoculum',  'float32')
		])

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
		dist = np.exp(-2.50*dist)
		dist[dist < 1e-2] = 0.0

		I[iZ,:] = dist.T.reshape((size_x * size_y))
		I[iZ,iZ] = 0
		iZ += 1


# simulation body
for crop in range(15):
	for i in range(time):
		print("day %d" % (i))
		population.updateAlive()
		population.updateExposition(I.dot(population.getExposition()))
		population.primaryInfection()
		population.updateInfective()
		population.updateInoculum()
		statistics = np.concatenate((statistics,population.getStatistics()))
	for i in range(365 - time):
		population.updateInoculum()
	population.nextCrop()

plt.plot(statistics['exposition'], label="exposition")
plt.plot(statistics['infective'], label="infective")
plt.plot(statistics['alive'], label="alive")
plt.plot(statistics['inoculum'], label="inoculum")
plt.xlabel('Days')
plt.ylabel('Value')
plt.legend(loc=0)
plt.show()

