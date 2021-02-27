import numpy as np
import pandas as pd
import scipy.stats as stats


class ArcadePopulation:

	def __init__(self, parameters_dict):

		primary_inoculum = list()
		k_omega = list()
		sigma = list()
		mu = list()

		global_parameters = parameters_dict['global_parameters']
		specific_parameters = parameters_dict['specific_parameters']
		interspecific_parameters = parameters_dict['interspecific_parameters']
		verbose = parameters_dict['metaparameters']['verbose']

		pathotypes_list = sorted(global_parameters['pathotypes'])
		self.genotype_probability = np.array(interspecific_parameters['genotype_probability'][:len(pathotypes_list) + 1])
		self.pathotypes = pathotypes_list
		self.loc = dict()

		i = 0

		for item in pathotypes_list:

			primary_inoculum.append(specific_parameters['omega'][item])
			k_omega.append(specific_parameters['k_omega'][item])
			sigma.append(specific_parameters['sigma'][item])
			mu.append(specific_parameters['mu'][item])
			self.loc[item] = i
			i += 1

		self.primary_inoculum_release = np.array(primary_inoculum)
		self.primary_inoculum_decay = np.array(k_omega)
		self.virulence_variance = np.array(sigma)
		self.virulence_mean = np.array(mu)
		self.dpi = parameters_dict['global_parameters']['dpi']

		titre = np.zeros((len(pathotypes_list), len(pathotypes_list)))
		for patho1_index in range(len(pathotypes_list)):
			for patho2_index in range(len(pathotypes_list)):
				patho1 = pathotypes_list[patho1_index]
				patho2 = pathotypes_list[patho2_index]
				titre[patho1_index, patho2_index] = interspecific_parameters['C'][patho1][patho2]
		self.titre = titre

		genotypes = np.zeros((len(pathotypes_list) + 1, len(pathotypes_list)))
		for i in range(len(pathotypes_list)):
			genotypes[i, i:] = 1

		self.__genotypes = genotypes

		if parameters_dict['metaparameters']['placement'] == 'regular_model':
			size_x = global_parameters['size_x']
			size_y = global_parameters['size_y']
			self.x = size_x
			self.y = size_y
			size = size_x * size_y
			xx_coords, yy_coords = self.set_grid_values(size_x, size_y)
		else:
			placement_file = parameters_dict['metaparameters']['placement']
			xx_coords, yy_coords = self.set_grid_values_from_file(placement_file)
			size_x = xx_coords.size
			size_y = 1
			size = xx_coords.size
			self.x = size_x
			self.y = size_y

		if verbose:
			print("arcadePopulation instance defined")
			print("setting a population of %d hosts" % (int(size_x*size_y)))

		self.__genotypeList = None
		common_fields = [
			('index', 'uint32'), ('rx', 'float32'), ('ry', 'float32'),
			('A', 'b1'), ('G', 'int8', len(pathotypes_list)),
			('I', 'b1'), ('NI', 'int8')
		]
		self.hosts = self.set_population_list(size_x, size_y, common_fields)
		self.hosts['rx'] = xx_coords
		self.hosts['ry'] = yy_coords
		self.hosts['A'][:] = True

		mortality_functions = dict(
			linear_model=self.linear_mortality_model,
			gaussian_model=self.gaussian_mortality_model,
			lognormal_model=self.lognormal_mortality_model,
			sanity_model=self.sanity_model
		)
		self.exposition = np.zeros((size, len(pathotypes_list)))
		self.infectious = np.zeros((size, len(pathotypes_list)))
		self.infectious_acumulative = np.zeros((size, len(pathotypes_list)))
		self.days_post_infection = np.zeros((size, len(pathotypes_list)))
		self.primary_inoculum = np.zeros((size, len(pathotypes_list)))
		self.life_span = mortality_functions[parameters_dict['metaparameters']['mortality']]()
		self.life_span = self.life_span.astype(int)

		try:
			self.coinfection_limit = parameters_dict['global_parameters']['coinfection']
		except KeyError:
			self.coinfection_limit = 2

	# STATIC METHODS
	# setPopulationList
	# setGridValues
	# createAttributeList

	def random_genotyping(self, population_list):
		if np.sum(self.genotype_probability) < 1:
			self.genotype_probability /= np.sum(self.genotype_probability)
			# raise UserWarning("genotypes of the host were not properly defined")
		size = self.x * self.y
		indexes = np.arange(size)
		random_genotypes = np.zeros(size, dtype='int8')
		indexes = np.random.permutation(indexes)
		split_vector = (self.genotype_probability * size).astype(int)
		split_vector = np.cumsum(split_vector)
		split_indexes = np.split(indexes, split_vector)
		for i in range(len(split_indexes)-1):
			random_genotypes[split_indexes[i]] = i
		self.__genotypeList = random_genotypes
		try:
			population_list['G'] = self.__genotypes[random_genotypes, :]
		except ValueError:
			temp = self.__genotypes[random_genotypes, :]
			population_list['G'] = temp.reshape(temp.size)

	def set_population_list(self, size_x, size_y, fields):
		total_size = size_y * size_x
		population = np.zeros(total_size, dtype=fields)
		population['index'] = np.arange(size_x*size_y)
		self.random_genotyping(population_list=population)
		return population

	@staticmethod
	def set_grid_values_from_file(filename):

		coords = pd.read_csv(filename)
		return coords['x'].values, coords['y'].values

	@staticmethod
	def set_grid_values(size_x, size_y, asimetric_y=2.0, separation_x=0.5):
		x_coordinates = np.arange(size_x, dtype='float32')
		y_coordinates = np.arange(size_y, dtype='float32')
		x_coordinates *= separation_x
		y_coordinates *= separation_x * asimetric_y
		xx, yy = np.meshgrid(x_coordinates,y_coordinates)
		xx = xx.astype('float32')
		yy = yy.astype('float32')
		return xx.reshape(size_x*size_y), yy.reshape(size_x*size_y)

	@staticmethod
	def create_attribute_list(list_attributes, dtype='float32'):
		attributes = list()
		for item in list_attributes:
			attributes.append((item, dtype))
		return attributes

	def set_seed(self, patho, method ='random', number=1):
		"""
		arcadeSpots3 - arcadePopulation - setSeed()
		:param method:
		:param patho:
		:param number:
		:return:

		It sets the exposition value of a given position to 1.0,
		setting the beginning of an epidemic. It must be specified
		for each pathotype of the epidemic
		"""
		try:
			patho_loc = self.loc[patho]
		except KeyError:
			raise KeyError("There is no {0} pathotype specified".format(patho))
		if method == 'random' :
			try:
				possible_hosts = np.where(self.hosts['G'][:, patho_loc] == 1)[0]
			except IndexError:
				possible_hosts = np.where(self.hosts['G'][:] == 1)[0]
			try:
				seed_value = np.random.choice(possible_hosts, number)
				self.exposition[seed_value, patho_loc] = 1.0
			except ValueError:
				return

	def get_pathotypes(self):
		"""
		arcadeSpots3 - arcadePopulation - getPathotypes()
		:return:

		It returns a list with the pathototypes that are included
		within the population.
		"""
		return self.pathotypes

	def get_coordinates(self):
		"""
		arcadeSpots3 - arcadePopulation - getCoordinates()
		:return: numpy array, numpy array

		It returns two arrays (x,y) with the coordinates of each host
		"""
		return np.copy(self.hosts['rx']), np.copy(self.hosts['ry'])

	def get_genotypes(self):
		"""

		:return:
		"""
		try:
			return np.copy(self.__genotypeList)
		except AttributeError:
			raise IOError("There was some issue dealing with the genotypes")

	def get_shape(self):
		return self.x, self.y

	def get_alive(self):
		"""
		arcadeSpots3 - arcadePopulation - getAlive()
		:return: numpy array

		It returns an array with the alive state of the hosts
		"""
		return self.hosts['A']

	def get_exposition(self):
		"""
		arcadeSpots3 - arcadePopulation - getExposition()
		:return: numpy array

		It returns an array with the exposition states of the hosts
		"""
		return self.exposition

	def get_infective(self):
		"""
		arcadeSpots3 - arcadePopulation - getInfective()
		:return: numpy array

		It returns an array with the infective states of the hosts
		"""
		return self.infectious

	def get_transmision(self):

		transmission = self.infectious.dot(self.titre.T)
		transmission[transmission < 0] = 0.0
		return transmission

	def get_days(self):
		"""
		arcadeSpots3 - arcadePopulation - getDays()
		:return: numpy array

		It returns an array with the exposition states of the hosts
		"""
		return self.days_post_infection

	def update_infective(self):
		"""
		arcadeSpots3 - arcadePopulation - getDays()
		:return:

		It updates the infective states depending on the days post infection,
		the alive state, and the number of coinfections
		"""
		# This is necessary to avoid double or triple infections

		control_infective = np.sum(self.infectious, axis=1)
		try:
			self.days_post_infection += (self.exposition >= 1.0) * self.hosts['G']
		except ValueError:
			self.days_post_infection += (self.exposition >= 1.0) * self.hosts['G'].reshape((self.hosts['G'].size, 1))
		try:
			delta_infectious = ((self.days_post_infection == self.dpi) * (control_infective < self.coinfection_limit))
			self.infectious += delta_infectious
			self.infectious_acumulative += delta_infectious
		except ValueError:
			delta_infectious = (
				(self.days_post_infection == self.dpi) * (control_infective < self.coinfection_limit).reshape(
					(control_infective.size, 1)
				)
			)
			self.infectious += delta_infectious
			self.infectious_acumulative += delta_infectious
		self.infectious = (self.hosts['A'] * self.infectious.T).T
		self.hosts['NI'][:] = np.sum(self.infectious, axis=1)

	def update_exposition(self, exposition_delta):
		"""
		arcadeSpots3 - arcadePopulation - updateExposition()
		:param exposition_delta: numpy array
		:return:

		Once provided a vector that specifies the contacts between infective hosts and
		other hosts, it updates the exposition value by multiplying the contact by
		the rate of transmission of each pathotype.
		"""
		self.exposition += self.hosts['A'].reshape((self.hosts['A'].size, 1)) * exposition_delta
		# sSIP : secondary Stochastic Infection Probability
		# sSIE : secondary Stochastic Infection Event
		self.exposition = np.clip(self.exposition, a_max=1.0, a_min=0)

	def stochastic_secondary_infections(self, S, method='tau_leap'):
		"""
		arcadeSpots3 - arcadePopulation - stochasticSecondaryInfections()
		:param S:
		:param method:
		:return:
		"""
		t = 0.0
		if method == 'gillespie':
			while t < 1.0:
				a0 = np.sum(S, axis=0)
				if np.all(a0 == 0):
					break
				tau_reactions = (1/a0) * np.log(1/np.random.rand(1))
				tau = np.min(tau_reactions)
				reaction = np.argmin(tau_reactions)
				t += tau
				au = S / np.sum(S, axis=0)
				k = np.random.choice(np.arange(len(au[:, reaction])), 1, replace=False, p=au[:, reaction])

				self.exposition[k] = 1.0
		elif method == 'tau_leap':
			for patho in range(S.shape[1]) :
				events = self.tau_leap_ssi(S[:, patho], 1.0)
				self.exposition[events, patho] = 1.0

	@staticmethod
	def tau_leap_ssi(probability_vector, time_step):
		"""

		:param probability_vector:
		:param time_step:
		:return:
		"""
		if len(probability_vector.shape) == 1:
			oA = np.sum(probability_vector)
			if oA == 0.0:
				return np.array([], dtype='int32')
			events = np.arange(probability_vector.size)
			n = np.random.poisson(oA*time_step)
			chosen_events = np.random.choice(events, n, replace=False, p=probability_vector / oA)
			return chosen_events
		else:
			raise ValueError("probability Vector must be unidimensional")

	def update_alive(self):
		"""
		arcadeSpots3 - arcadePopulation - updateAlive()
		:return:

		Updates the primary inoculum, which gets increased when hosts die, and it modifies
		the values of the alive state
		"""

		ref_vals = np.array([self.titre[i, i] for i in range(len(self.pathotypes))])
		ref_vals = ref_vals.reshape(1, len(self.pathotypes))
		self.primary_inoculum += (self.days_post_infection == self.life_span) * (self.get_transmision() / ref_vals) * self.primary_inoculum_release
		self.hosts['A'][np.sum(self.days_post_infection >= self.life_span, axis=1) >= 1.0] = False

	def update_inoculum(self):
		"""
		arcadeSpots3 - arcadePopulation - updateAlive()
		:return:

		Updates the primary inoculum upon degradation
		"""
		self.primary_inoculum *= self.primary_inoculum_decay

	def primary_infection(self):
		"""
		arcadeSpots3 - arcadePopulation - primaryInfection()
		:return:

		Updates the exposition depending on the linear relationship between
		the concentration of primary inoculum in the soil.
		"""
		base_prob = 5.4e-4
		random_values = (
				np.random.rand(self.x * self.y, len(self.pathotypes)) < (self.primary_inoculum / 1000000) * base_prob)
		primary_infections = (self.hosts['A'] * random_values.T).T
		self.exposition += primary_infections
		self.exposition[self.exposition > 1.0] = 1.0

	def next_crop(self):
		"""
		arcadeSpots3 - arcadePopulation - nextCrop()
		:return:

		Resets all the values but the primary inoculum
		"""
		self.primary_inoculum[:, :] += self.primary_inoculum_release * (
				(self.hosts['A'] == True) * (self.infectious == True).T).T
		self.exposition[:, :] = 0.0
		self.infectious[:, :] = False
		self.infectious_acumulative[:, :] = False
		self.days_post_infection[:, :] = 0
		self.hosts['A'][:] = True

	def modify_population(self, ind, parameter, value):
		try:
			self.hosts[parameter]
		except KeyError:
			raise IOError("there is not such an attribute %s" % parameter)
		try:
			self.hosts[parameter][ind]
		except IndexError:
			raise IOError("there is not such an individual %d" % ind)
		self.hosts[parameter][ind] = value

	def dump_spatial_data(self, filename, patho, time, crop):

		try:
			patho_loc = self.loc[patho]
		except KeyError:
			raise KeyError("There is no such a {0} pathotype".format(patho))

		rx = self.hosts['rx']
		ry = self.hosts['ry']
		exposition = self.infectious[:, patho_loc]
		alive = self.hosts['A']
		inoculum = self.primary_inoculum[:, patho_loc]

		filename = filename + '_{:03d}_{:03d}.spatial.csv'.format(crop, time)
		with open(filename, 'w') as f:
			for i in range(len(rx)):
				f.write('{:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f}\n'.format(
					rx[i], ry[i], exposition[i], alive[i], inoculum[i]
				))

	def dump_coinfection_data(self, filename, time, crop):

		coinf = []
		for i, p in enumerate(self.pathotypes):
			for j, q in enumerate(self.pathotypes):
				if i >= j:
					continue
				else:
					inf_i = self.infectious[:, i]
					inf_j = self.infectious[:, j]
					coinf.append(float((inf_i * inf_j).sum() / (self.y * self.x)))
		with open(filename, 'a') as f:
			f.write('{:d} {:d} '.format(crop, time))
			f.write(''.join(['{:8.4f} '.format(x) for x in coinf]))
			f.write('\n')

	def get_statistics(self, patho):
		statistics = np.zeros(1, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('infective_ac', 'float32'),
			('alive',      'float32'),
			('inoculum',   'float32'),
			('coinfected', 'float32')
		])
		total_pop = self.y * self.x
		try:
			patho_loc = self.loc[patho]
		except KeyError:
			raise KeyError("There is no such a {0} pathotype".format(patho))
		statistics['exposition'] = np.sum(self.exposition[:, patho_loc]) / total_pop
		statistics['infective'] = np.sum(self.infectious[:, patho_loc]) / total_pop
		statistics['alive'] = np.sum(self.hosts['A']) / total_pop
		statistics['inoculum'] = np.sum(self.primary_inoculum[:, patho_loc]) / (total_pop * 1000)
		statistics['coinfected'] = np.sum(self.hosts['NI'] > 1) / total_pop
		statistics['infective_ac'] = np.sum(self.infectious_acumulative[:, patho_loc]) / total_pop
		return statistics

	def get_index(self, k):
		i = int(k / self.x)
		j = k - i*self.x
		return i, j

	def linear_mortality_model(self):
		x = np.random.rand(self.x * self.y, len(self.pathotypes))
		tau = np.log(x) / (-self.virulence_variance)
		tau = np.round(tau, 0)
		tau += self.virulence_mean
		tau += self.dpi
		return tau

	def gaussian_mortality_model(self):

		x = np.random.rand(self.x * self.y, len(self.pathotypes))
		u = stats.norm.ppf(x)
		return u * self.virulence_variance + self.virulence_mean + self.dpi

	def lognormal_mortality_model(self):
		x = np.random.rand(self.x * self.y, len(self.pathotypes))
		u = stats.norm.ppf(x)
		return np.exp(u * self.virulence_variance + self.virulence_mean) + self.dpi

	def sanity_model(self):
		hazard = self.virulence_mean * np.exp(self.virulence_variance)
		x = np.random.rand(self.x * self.y, len(self.pathotypes))
		u = np.log(x)*(1 / -hazard)
		return u.astype(int) + self.dpi
