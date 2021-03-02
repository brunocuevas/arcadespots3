import numpy as np
from scipy import sparse
import arcadePopulation as aP
from arcadeUtils import save_statistics
import json
import os

parametersDict = list([
	dict(
		global_parameters=dict(
			pathotypes=['P0', 'P12'],
			size_x=50,
			size_y=50,
			crop_time=180,
			crops=2,
			k_beta=1.25,
			k_alpha=0.1,
			c=1e-5,
			epsilon_beta=0.01,
			epsilon_alpha=0.25,
			dpi=7
		),
		specific_parameters=dict(
			omega=dict(
				P0=3e13,
				P12=2.4e12,
				P123=1.9e12
			),
			k_omega=dict(
				P0=0.906,
				P12=0.920,
				P123=0.917
			),
			sigma=dict(
				P0=0.32,
				P12=0.32,
				P123=0.32
			),
			mu=dict(
				P0=3.46,
				P12=4.19,
				P123=3.74
			)
		),
		interspecific_parameters=dict(
			C=dict(
				P0=dict(P0=1.00, P12=0.0, P123=0.0),
				P12=dict(P0=-0.409, P12=0.489, P123=-0.275),
				P123=dict(P0=-0.251, P12=-0.155, P123=0.253)),
			genotype_probability=list([0.5, 0.3, 0.2])
		),
		metaparameters=dict(
			mortality='lognormal_model',
			map=False,
			simulations=5,
			seeds=[1, 2, 3, 4, 5],
			outfile='trial',
			verbose=False,
			placement='regular_model'
		)
	)
])


# initialization
def parameters_to_json(dict2store, name):
	f = open(name, 'w')
	f.write(json.dumps(dict2store, sort_keys=True, indent=4))
	f.close()


class ArcadeSimulator:
	def __init__(self, parameter_dictionary, sim=0, verbose=False):
		self.verbose = verbose
		self.sim = sim
		try:
			seed = parameter_dictionary['metaparameters']['seeds'][sim]
		except IndexError:
			seed = sim
		verbose = parameter_dictionary['metaparameters']['verbose']
		self.__v = verbose
		if verbose:
			("\tworking with seed {0}".format(seed))
		np.random.seed(seed)
		self.param_dict = parameter_dictionary

		global_parameters = parameter_dictionary['global_parameters']
		self.k_stochastic = parameter_dictionary['global_parameters']['k_alpha']
		self.k_determinist = parameter_dictionary['global_parameters']['k_beta']
		self.e_stochastic = parameter_dictionary['global_parameters']['epsilon_alpha']
		self.e_determinist = parameter_dictionary['global_parameters']['epsilon_beta']
		self.stochastic_factor = parameter_dictionary['global_parameters']['c']

		self.cropTime = parameter_dictionary['global_parameters']['crop_time']
		self.numberCrops = parameter_dictionary['global_parameters']['crops']

		self.population = aP.ArcadePopulation(parameter_dictionary)
		size_x, size_y = self.population.get_shape()
		self.sx = size_x
		self.sy = size_y
		self.stats = None
		self.rx = None
		self.ry = None
		self.infected = None
		self.susceptible = None

		for patho in global_parameters['pathotypes']:
			self.population.set_seed(patho=patho)

	def __log(self, message):

		if self.verbose:
			print("[arcadeSpots3] %s " % message)

	def set_statistics(self):
		statistics = dict()
		for patho in self.param_dict['global_parameters']['pathotypes']:
			statistics[patho] = np.empty(0, dtype=[
				('crop', 'int8'),
				('time', 'int64'),
				('exposition', 'float32'),
				('infective', 'float32'),
				('infective_ac', 'float32'),
				('alive', 'float32'),
				('inoculum', 'float32'),
				('coinfected', 'float32')
			])
		self.stats = statistics

	def update_statistics(self, crop, time):
		for patho in self.param_dict['global_parameters']['pathotypes']:
			st = self.population.get_statistics(patho)
			st_plus = np.zeros(
				1, dtype=[
					('crop', 'int8'),
					('time', 'int64'),
					('exposition', 'float32'),
					('infective', 'float32'),
					('infective_ac', 'float32'),
					('alive', 'float32'),
					('inoculum', 'float32'),
					('coinfected', 'float32')
				]
			)
			st_plus[0]['crop'] = np.array(crop)
			st_plus[0]['time'] = np.array(time)
			st_plus[0]['exposition'] = st['exposition']
			st_plus[0]['infective'] = st['infective']
			st_plus[0]['infective_ac'] = st['infective_ac']
			st_plus[0]['alive'] = st['alive']
			st_plus[0]['inoculum'] = st['inoculum']
			st_plus[0]['coinfected'] = st['coinfected']
			self.stats[patho] = np.concatenate((self.stats[patho], st_plus))

	def build_matrix(self):

		size_xy = self.sx * self.sy
		self.rx, self.ry = self.population.get_coordinates()
		self.infected = np.zeros((size_xy, size_xy))
		self.susceptible = np.zeros((size_xy, size_xy))

		for i in np.arange(self.sx * self.sy):
			x_coordinates, y_coordinates = self.population.get_coordinates()
			x_coordinates -= x_coordinates[i]
			y_coordinates -= y_coordinates[i]

			dist = np.sqrt((x_coordinates ** 2) + (y_coordinates ** 2))
			dc = np.exp(-self.k_determinist * dist)
			sc = np.exp(-self.k_stochastic * dist)

			dc[dc < self.e_determinist] = 0.0
			sc[sc < self.e_stochastic] = 0.0

			self.infected[i, :] = dc
			self.infected[i, i] = 0

			self.susceptible[i, :] = self.stochastic_factor * sc
			self.susceptible[i, i] = 0
		self.infected = sparse.csc_matrix(self.infected)
		self.susceptible = sparse.csc_matrix(self.susceptible)

	def get_matrix(self):
		try:
			return self.infected, self.susceptible
		except AttributeError:
			self.build_matrix()
			return self.infected, self.susceptible

	def simulate(self, refresh=10):
		self.set_statistics()
		self.build_matrix()
		s_precalc = np.zeros(self.sx * self.sy)
		self.__log("starting simulation. refresh = %d" % refresh)
		for crop in range(self.numberCrops):
			for i in range(self.cropTime):
				self.__log("simulation day %d cycle %d" % (i, crop))
				self.population.update_alive()
				self.population.primary_infection()
				if i % refresh == 0:
					s_precalc = self.susceptible.dot(self.population.get_transmision())
				tmp_transmission = self.population.get_transmision()
				self.population.update_exposition(self.infected.dot(tmp_transmission))
				self.population.stochastic_secondary_infections(s_precalc, method='tau_leap')
				self.population.update_infective()
				self.population.update_inoculum()
				self.update_statistics(crop, i)
				if i % refresh == 0 and self.param_dict['metaparameters']['map']:
					for patho in self.param_dict['global_parameters']['pathotypes']:
						folder = self.param_dict['metaparameters']['outfile']
						if os.path.exists(folder):
							pass
						else:
							os.makedirs(folder)
						# name = './%s/spatial_map_sim%s_patho%s_crop%d_time%d.png'
						# name = name % (folder, self.sim, patho, crop, i)
						# self.population.spatialMap(patho, time=i, crop=crop, filename=name)
						self.population.dump_spatial_data(
							filename='./%s/%s' % (folder, patho), patho=patho,
							time=i, crop=crop
						)
				if 'overlap' in self.param_dict['metaparameters']:
					if i % refresh == 0 and self.param_dict['metaparameters']['overlap']:
						folder = self.param_dict['metaparameters']['outfile']
						try:
							os.makedirs(folder)
						except FileExistsError:
							pass
						self.population.dump_coinfection_data(
							filename='./{:s}/coinfection.csv'.format(folder, ),
							time=i, crop=crop
						)
			self.__log("next crop")
			self.population.next_crop()
			for j in range(365 - self.cropTime):
				self.population.update_inoculum()

	def get_population(self):
		return self.population

	def save_report(self):
		"""
		:param
		:return:
		"""
		import os
		folder_name = self.param_dict['metaparameters']['outfile']
		try:
			os.makedirs(folder_name)
		except FileExistsError:
			pass
		for patho in self.param_dict['global_parameters']['pathotypes']:
			statistics = self.stats[patho]
			filename = './%s/timeSeriesStatistics_%s.csv' % (folder_name, patho)
			save_statistics(
				filename=filename,
				statistics=statistics, sim=self.sim
			)
