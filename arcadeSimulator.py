#!/home/charizard/anaconda3/bin/python3
import numpy as np
from scipy import sparse
from arcadeUtils import displayer
#from memory_profiler import profile
import arcadePopulation as aP
import matplotlib.pyplot as plt
import seaborn as sns


# NOTE : This is an script. Not the real object that will be used later
sns.set_style('whitegrid')
sns.set_palette('Set2')


# initialization
parametersDict = dict(
		global_parameters = dict(
			pathotypes = ['P0', 'P12'],
			size_x = 50,
			size_y = 50,
			crop_time = 180,
			crops  = 2
		),
		specific_parameters = dict(
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
			k = dict(
				P0 = 0.32,
				P12 = 0.32
			),
			s = dict(
				P0 = 4.158,
				P12 = 4.888
			)
		),
		interspecific_parameters = dict(
			C = dict(
				P0 = dict(
					P0  =  1.00,
					P12 = -0.05
				),
				P12 = dict(
					P0  = -0.35,
					P12 =  0.60
				)
			),
			genotype_probability = list([0.5, 0.3, 0.2])
		),
		metaparameters = dict(
			mortality = 'lognormal_model',
			simulations = 5,
			seeds = [1,2,3,4,5],
			outname = 'trial',
			verbose = False
		)

	)
def parameters2JSON (dict2store, name):
	import json
	f = open(name, 'w')
	f.write(json.dumps(dict2store, sort_keys=True, indent=4))
	f.close()

class arcadeSimulator:
	def __init__(self, parameter_dictionary, sim=0):
		self.__sim = sim
		try :
			seed = parameter_dictionary['metaparameters']['seeds'][sim]
		except IndexError:
			seed = sim
		verbose = parameter_dictionary['metaparameters']['verbose']
		self.__v = verbose
		if verbose  :
			("\tworking with seed {0}".format(seed))
		np.random.seed(seed)
		self.__parameter_dictionary = parameter_dictionary
		try :
			global_parameters   = parameter_dictionary['global_parameters']
			size_x = global_parameters['size_x']
			size_y = global_parameters['size_y']
			self.__sx = size_x
			self.__sy = size_y
			self.__cropTime    = global_parameters['crop_time']
			self.__numberCrops = global_parameters['crops']
		except KeyError :
			raise IOError("some of the parameters is not well specified: size_x, size_y, time, crops")
		self.__disp = displayer(crops=self.__numberCrops, cropDays=self.__cropTime)
		randomSeed = np.random.randint(0, high=size_x*size_y,
									   size=len(global_parameters['pathotypes']))

		self.__population = aP.arcadePopulation(size_x, size_y, parameter_dictionary)
		self.__population.setSeed(randomSeed[0], patho='P0')
		self.__population.setSeed(randomSeed[1], patho='P12')

	def __setStatistics(self):
		statistics = dict()
		for patho in self.__parameter_dictionary['global_parameters']['pathotypes'] :
			statistics[patho] = np.empty(0, dtype=[
				('crop', 'int8'),
				('time', 'int64'),
				('exposition', 'float32'),
				('infective', 'float32'),
				('alive', 'float32'),
				('inoculum', 'float32')
			])
		self.__stats = statistics

	def __updateStatistics(self, crop, time):
		for patho in self.__parameter_dictionary['global_parameters']['pathotypes'] :
			st = self.__population.getStatistics(patho)
			st_plus = np.zeros(1,
							   dtype=[
								   ('crop', 'int8'),('time', 'int64'),
								   ('exposition', 'float32'),
								   ('infective', 'float32'),
								   ('alive', 'float32'),('inoculum', 'float32')
							   ]
							   )
			st_plus[0]['crop'] = np.array(crop)
			st_plus[0]['time'] = np.array(time)
			st_plus[0]['exposition'] = st['exposition']
			st_plus[0]['infective'] = st['infective']
			st_plus[0]['alive'] = st['alive']
			st_plus[0]['inoculum'] = st['inoculum']
			self.__stats[patho] = np.concatenate((self.__stats[patho], st_plus))

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
	#@profile()
	def simulate(self, refresh = 10):
		self.__setStatistics()
		self.__buildMatrix()
		S_precalc = np.zeros(self.__sx * self.__sy)
		if self.__v :
			self.__disp.start()
		for crop in range(self.__numberCrops):
			for i in range(self.__cropTime):
				if self.__v :
					self.__disp.update()

				self.__population.updateAlive()
				self.__population.primaryInfection()
				if i % refresh == 0:
					S_precalc = self.__S.dot(self.__population.getTranmission())
				self.__population.updateExposition(self.__I.dot(self.__population.getTranmission()))
				self.__population.stochasticSecondaryInfections(S_precalc, method='tau_leap')
				self.__population.updateInfective()
				self.__population.updateInoculum()
				self.__updateStatistics(crop, i)
				if i % 50 == 0 :
					self.__population.spatialMap('P0', i, crop)
			self.__population.nextCrop()
			for j in range(365 - self.__cropTime):
				self.__population.updateInoculum()

	def plotStatistics(self, stat='infective'):
		for item in self.__parameter_dictionary['global_parameters']['pathotypes'] :
			try :
				plt.plot(self.__stats[item][stat], label=item)
			except KeyError:
				raise IOError("the stat was not specified")
		plt.xlabel('Days')
		plt.ylabel(stat)
		plt.legend(loc=0)
		plt.show()
	def getPopulation(self):
		return self.__population
	def saveReport(self, folderName, compress = False, header=True):
		"""

		:param folderName:
		:param compress
		:param header
		:return:
		"""
		import os
		import pandas as pd
		try :
			os.makedirs(folderName)
		except FileExistsError:
			pass
		for patho in self.__parameter_dictionary['global_parameters']['pathotypes'] :
			statistics = pd.DataFrame(self.__stats[patho])
			statistics['sim '] = self.__sim
			f = open('./%s/timeSeriesStatistics_%s.csv' % (folderName, patho), 'a+')
			if self.__sim == 0 :
				f.write(statistics.to_csv(sep=",", index=False, float_format='%8.4f',
										  header=header))
			else:
				f.write(statistics.to_csv(sep=",", index=False, header=False,
										  float_format='%8.4f'))
			f.close()
		if compress:
			pass
	def printHeader(self, folderName):
		import os
		for patho in self.__parameter_dictionary['global_parameters']['pathotypes']:
			try:
				os.makedirs(folderName)
			except FileExistsError:
				pass
			f = open('./%s/timeSeriesStatistics_%s.csv' % (folderName, patho), 'w')
			f.write('crop,time,exposition,infective,alive,inoculum,sim\n')
			f.close()


