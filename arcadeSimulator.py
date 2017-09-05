#!/home/charizard/anaconda3/bin/python3
import numpy as np
from scipy import sparse
from arcadeUtils import displayer
import arcadePopulation as aP
import matplotlib.pyplot as plt
import seaborn as sns


# NOTE : This is an script. Not the real object that will be used later
sns.set_style('whitegrid')
sns.set_palette('Set2')


# initialization
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
		]),
		size_x = 50,
		size_y = 50,
		crop_time = 180,
		crops  = 2
	)
class arcadeSimulator:
	def __init__(self, parameter_dictionary, seed = 42):
		np.random.seed(seed)
		self.__parameter_dictionary = parameter_dictionary
		try :

			size_x = 	parameter_dictionary['size_x']
			size_y = parameter_dictionary['size_y']
			self.__sx = size_x
			self.__sy = size_y

			self.__cropTime = parameter_dictionary['crop_time']
			self.__numberCrops = parameter_dictionary['crops']
		except KeyError :
			raise IOError("some of the parameters is not well specified: size_x, size_y, time, crops")
		self.__disp = displayer(crops=self.__numberCrops, cropDays=self.__cropTime)
		randomSeed = np.random.randint(0, high=size_x*size_y,
									   size=len(parameter_dictionary['pathotypes']))

		self.__population = aP.arcadePopulation(size_x, size_y, parameter_dictionary)
		self.__population.setSeed(randomSeed[0], patho='P0')
		self.__population.setSeed(randomSeed[1], patho='P12')

	def __setStatistics(self):
		statistics = dict()
		for patho in self.__parameter_dictionary['pathotypes'] :
			statistics[patho] = np.empty(0, dtype=[
				('exposition', 'float32'),
				('infective', 'float32'),
				('alive', 'float32'),
				('inoculum', 'float32')
			])
		self.__stats = statistics

	def __updateStatistics(self):
		for patho in self.__parameter_dictionary['pathotypes'] :
			self.__stats[patho] = np.concatenate((self.__stats[patho], self.__population.getStatistics(patho)))

	def __buildMatrix(self):
		size_xy = self.__sx*self.__sy
		self.__rx, self.__ry = self.__population.getCoordinates()
		self.__I = sparse.lil_matrix((size_xy, size_xy))
		self.__S = sparse.lil_matrix((size_xy, size_xy))

		for i in np.arange(self.__sx * self.__sy):
			x_coordinates, y_coordinates = self.__population.getCoordinates()
			x_coordinates -= x_coordinates[i]
			y_coordinates -= y_coordinates[i]

			# xx, yy = np.meshgrid(x_coordinates, y_coordinates)
			dist = np.sqrt((x_coordinates ** 2) + (y_coordinates ** 2))
			## dc holds for direct contact
			## sc holds for stochastic contact
			dc = np.exp(-1.25 * dist)
			sc = np.exp(-0.1 * dist)

			dc[dc < 1e-2] = 0.0
			sc[sc < 0.25] = 0.0

			self.__I[i, :] = dc  # .T.reshape((size_x * size_y))
			self.__I[i, i] = 0

			self.__S[i, :] = 1e-5 * sc  # .T.reshape((size_x * size_y))
			self.__S[i, i] = 0

	def simulate(self, refresh = 10):
		self.__setStatistics()
		self.__buildMatrix()
		S_precalc = np.zeros(self.__sx * self.__sy)
		self.__disp.start()
		for crop in range(self.__numberCrops):
			for i in range(self.__cropTime):
				self.__disp.update()

				self.__population.updateAlive()
				self.__population.primaryInfection()
				if i % refresh == 0:
					S_precalc = self.__S.dot(self.__population.getTranmission())
				self.__population.updateExposition(self.__I.dot(self.__population.getTranmission()))
				self.__population.stochasticSecondaryInfections(S_precalc, method='tau_leap')
				self.__population.updateInfective()
				self.__population.updateInoculum()
				self.__updateStatistics()
			for j in range(365 - self.__cropTime):
				self.__population.updateInoculum()
			self.__population.nextCrop()
	def plotStatistics(self, stat='infective'):
		for item in self.__parameter_dictionary['pathotypes'] :
			try :
				plt.plot(self.__stats[item][stat], label=item)
			except KeyError:
				raise IOError("the stat was not specified")
		plt.xlabel('Days')
		plt.ylabel(stat)
		plt.legend(loc=0)
		plt.show()

