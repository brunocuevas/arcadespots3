#!/home/charizard/anaconda3/bin/python3
import numpy as np
from scipy import sparse
#from arcadeUtils import displayer
#from memory_profiler import profile
import arcadePopulation as aP
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import time

# NOTE : This is an script. Not the real object that will be used later
sns.set_style('whitegrid')
sns.set_palette('Set2')


parametersDict = list([
	dict(
		global_parameters = dict(
			pathotypes = ['P0', 'P12'],
			size_x = 50,
			size_y = 50,
			crop_time = 180,
			crops  = 2,
			k_beta = 1.25,
			k_alpha = 0.1,
			c = 1e-5,
			epsilon_beta  = 0.01,
			epsilon_alpha = 0.25,
			dpi = 7
		),
		specific_parameters = dict(
			omega = dict(
				P0 = 3e13,
				P12 = 2.4e12,
				P123 = 1.9e12
			),
			k_omega = dict(
				P0 = 0.906,
				P12 = 0.920,
				P123 = 0.917
			),
			sigma = dict(
				P0 = 0.32,
				P12 = 0.32,
				P123 = 0.32
			),
			mu = dict(
				P0 = 3.46,
				P12 = 4.19,
				P123 = 3.74
			)
		),
		interspecific_parameters = dict(
			C = dict(
				P0 = dict(P0  =  1.00,P12 =  0.0,P123 = 0.0),
				P12 = dict(P0  = -0.409,P12 =  0.489,P123 =-0.275),
				P123 = dict(P0 = -0.251,P12 = -0.155,P123 = 0.253)),
			genotype_probability = list([0.5, 0.3, 0.2])
		),
		metaparameters = dict(
			mortality = 'lognormal_model',
			map = False,
			simulations = 5,
			seeds = [1,2,3,4,5],
			outfile = 'trial',
			verbose = False,
			placement='regular_model'
		)
	)
])
# initialization
def parameters2JSON (dict2store, name):
	import json
	f = open(name, 'w')
	f.write(json.dumps(dict2store, sort_keys=True, indent=4))
	f.close()

class arcadeSimulator:
	def __init__(self, parameter_dictionary, sim=0, verbose = False):
		self.__verbose = verbose
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

		global_parameters   = parameter_dictionary['global_parameters']



		self.__k_stochastic   = parameter_dictionary['global_parameters']['k_alpha']
		self.__k_determinist  = parameter_dictionary['global_parameters']['k_beta']
		self.__e_stochastic   = parameter_dictionary['global_parameters']['epsilon_alpha']
		self.__e_determinist  = parameter_dictionary['global_parameters']['epsilon_beta']
		self.__stochastic_factor = parameter_dictionary['global_parameters']['c']

		self.__cropTime    = parameter_dictionary['global_parameters']['crop_time']
		self.__numberCrops = parameter_dictionary['global_parameters']['crops']

		self.__population = aP.arcadePopulation(parameter_dictionary)
		size_x, size_y    = self.__population.getShape()
		self.__sx = size_x
		self.__sy = size_y

		for patho in global_parameters['pathotypes']:
			self.__population.setSeed(patho=patho)

	def __log(self, message):

		if self.__verbose :
			ltime = time.localtime()
			print("[arcadeSpots3] %s %02d:%02d:%02d" % (message, ltime.tm_hour, ltime.tm_min, ltime.tm_sec))
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

	def __buildMatrix(self, ):
		size_xy = self.__sx*self.__sy
		self.__rx, self.__ry = self.__population.getCoordinates()
		#self.__I = sparse.csc_matrix((size_xy, size_xy))
		#self.__S = sparse.csc_matrix((size_xy, size_xy))
		self.__I = np.zeros((size_xy, size_xy))
		self.__S = np.zeros((size_xy, size_xy))

		for i in np.arange(self.__sx * self.__sy):
			x_coordinates, y_coordinates = self.__population.getCoordinates()
			x_coordinates -= x_coordinates[i]
			y_coordinates -= y_coordinates[i]


			dist = np.sqrt((x_coordinates ** 2) + (y_coordinates ** 2))
			## dc holds for direct contact
			## sc holds for stochastic contact
			dc = np.exp(-self.__k_determinist * dist)
			sc = np.exp(-self.__k_stochastic * dist)

			dc[dc < self.__e_determinist] = 0.0
			sc[sc < self.__e_stochastic] = 0.0

			self.__I[i, :] = dc
			self.__I[i, i] = 0

			self.__S[i, :] = self.__stochastic_factor * sc
			self.__S[i, i] = 0
		self.__I = sparse.csc_matrix(self.__I)
		self.__S = sparse.csc_matrix(self.__S)

	def getMatrix(self):
		try:
			return self.__I, self.__S
		except AttributeError:
			self.__buildMatrix()
			return self.__I, self.__S
	def simulate(self, refresh = 10):
		import  os
		self.__setStatistics()
		self.__buildMatrix()
		S_precalc = np.zeros(self.__sx * self.__sy)
		self.__log("starting simulation. refresh = %d" % refresh)
		for crop in range(self.__numberCrops):
			for i in range(self.__cropTime):
				self.__log("simulation day %d cycle %d" % (i, crop))
				self.__population.updateAlive()
				self.__population.primaryInfection()
				if i % refresh == 0:
					S_precalc = self.__S.dot(self.__population.getTranmission())
				tmp_transmission = self.__population.getTranmission()
				self.__population.updateExposition(self.__I.dot(tmp_transmission))
				self.__population.stochasticSecondaryInfections(S_precalc, method='tau_leap')
				self.__population.updateInfective()
				self.__population.updateInoculum()
				self.__updateStatistics(crop, i)
				if i % refresh == 0 and self.__parameter_dictionary['metaparameters']['map']:
					for patho in self.__parameter_dictionary['global_parameters']['pathotypes']:
						folder = self.__parameter_dictionary['metaparameters']['outfile']
						if os.path.exists(folder):
							pass
						else:
							os.makedirs(folder)
						name = './%s/spatial_map_sim%s_patho%s_crop%d_time%d.png'
						name = name % (folder, self.__sim, patho, crop, i)
						self.__population.spatialMap(patho, time=i, crop=crop, filename=name)
			self.__log("next crop")
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
	def saveReport(self):
		"""
		:param
		:return:
		"""
		import os
		import pandas as pd
		folderName = self.__parameter_dictionary['metaparameters']['outfile']
		try :
			os.makedirs(folderName)
		except FileExistsError:
			pass
		for patho in self.__parameter_dictionary['global_parameters']['pathotypes'] :
			statistics = pd.DataFrame(self.__stats[patho])
			statistics['sim '] = self.__sim
			filename = './%s/timeSeriesStatistics_%s.csv' % (folderName, patho)
			if os.path.exists(filename):
				f = open('./%s/timeSeriesStatistics_%s.csv' % (folderName, patho), 'a+')
				f.write(statistics.to_csv(sep=",", index=False, float_format='%12.4f',
										  header=False))
			else:
				f = open('./%s/timeSeriesStatistics_%s.csv' % (folderName, patho), 'w')
				f.write(statistics.to_csv(sep=",", index=False, header=True,
										  float_format='%12.4f'))
			f.close()


# I think that by now this wil be placed outside the object
def arcadeOutput(params, population):
	"""

	:param params:
	:param population:
	:return:
	"""
	import pandas as pd
	import json
	folderName = params['metaparameters']['outfile']
	excel = pd.ExcelWriter(folderName + '/timeSeries.xlsx')
	for patho in params['global_parameters']['pathotypes'] :
		temporal_table = pd.read_csv('./%s/timeSeriesStatistics_%s.csv' % (folderName, patho))
		temporal_table.to_excel(excel_writer=excel, sheet_name=patho, index=False)
	excel.save()
	f = open('./%s/simulation_parameters.json' % folderName, 'w')
	f.write(json.dumps(params, indent=4, sort_keys=True))
	f.close()
	coord_x, coord_y = population.getCoordinates()
	genotypes        = population.getGenotypes()
	coordinates_dataframe = pd.DataFrame.from_items([
		('x', coord_x),
		('y', coord_y),
		('genotype', genotypes)
	])
	f = open('./%s/host_coordinates_file.csv' % folderName, 'w')
	f.write(coordinates_dataframe.to_csv(
		index=False,
		sep=","
	))




