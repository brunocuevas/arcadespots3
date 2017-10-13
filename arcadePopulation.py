# arcadePopulation
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


class arcadePopulation:

	def __init__(self, parametersDict):


		primaryInoc = list()
		k_omega = list()
		sigma = list()
		mu = list()

		global_parameters   = parametersDict['global_parameters']
		specific_parameters = parametersDict['specific_parameters']
		interspecific_parameters = parametersDict['interspecific_parameters']
		verbose = parametersDict['metaparameters']['verbose']



		pathotypesList = sorted(global_parameters['pathotypes'])
		self.__gP = np.array(interspecific_parameters['genotype_probability'][:len(pathotypesList)+1])
		self.__patho = pathotypesList
		self.__loc   = dict()

		i = 0

		for item in pathotypesList:

			primaryInoc.append(specific_parameters['omega'][item])
			k_omega.append(specific_parameters['k_omega'][item])
			sigma.append(specific_parameters['sigma'][item])
			mu.append(specific_parameters['mu'][item])
			self.__loc[item] = i
			i += 1


		self.__omega   = np.array(primaryInoc)
		self.__k_omega = np.array(k_omega)
		self.__sigma = np.array(sigma)
		self.__mu = np.array(mu)
		self.__dpi = parametersDict['global_parameters']['dpi']


		C = np.zeros((len(pathotypesList), len(pathotypesList)))
		for patho1_index in range(len(pathotypesList)) :
			for patho2_index in range(len(pathotypesList)) :
				patho1 = pathotypesList[patho1_index]
				patho2 = pathotypesList[patho2_index]
				C[patho1_index, patho2_index] = interspecific_parameters['C'][patho1][patho2]
		self.__C    = C


		genotypes = np.zeros((len(pathotypesList) + 1, len(pathotypesList)))
		for i in range(len(pathotypesList)):
			genotypes[i,i:] = 1

		self.__genotypes = genotypes

		if parametersDict['metaparameters']['placement'] == 'regular_model' :
			size_x = global_parameters['size_x']
			size_y = global_parameters['size_y']
			self.__sx = size_x
			self.__sy = size_y
			size = size_x * size_y
			xx_coords, yy_coords = self.__setGridValues(size_x, size_y)
		else:
			placement_file = parametersDict['metaparameters']['placement']
			xx_coords, yy_coords = self.__setGridValuesFromFile(placement_file)
			size_x = xx_coords.size
			size_y = 1
			size   = xx_coords.size
			self.__sx = size_x
			self.__sy = size_y

		if verbose :
			print("arcadePopulation instance defined")
			print("setting a population of %d hosts" %(int(size_x*size_y)))

		self.__genotypeList = None
		commonFields = [('index', 'uint32'),
						('rx', 'float32'),
						('ry', 'float32'),
						('A', 'b1'),
						('G', 'int8', len(pathotypesList)	),
						('I', 'b1')]
		self.__P = self.setPopulationList(size_x, size_y, commonFields)

		self.__P['rx'] = xx_coords
		self.__P['ry'] = yy_coords
		self.__P['A'][:] = True

		if verbose : print("common attributes set")
		mortalityFunctions = dict(
			linear_model   = self.__linearMortalityModel,
			gaussian_model = self.__gaussianMortalityModel,
			lognormal_model = self.__logNormalMortalityModel,
			sanity_model   = self.__sanityModel
		)
		self.__x = np.zeros((size, len(pathotypesList)))
		self.__y = np.zeros((size, len(pathotypesList)))
		self.__d = np.zeros((size, len(pathotypesList)))
		self.__w = np.zeros((size, len(pathotypesList)))
		self.__mortality = mortalityFunctions[parametersDict['metaparameters']['mortality']]()
		self.__mortality = self.__mortality.astype(int)


		print("population set")

	# STATIC METHODS
	## setPopulationList
	## setGridValues
	## createAttributeList

	def randomGenotyping(self, populationList):
		if np.sum(self.__gP) < 1 :
			self.__gP /= np.sum(self.__gP)
			#raise UserWarning("genotypes of the host were not properly defined")
		S = self.__sx*self.__sy
		randomGenotypes = np.random.choice(np.arange(len(self.__gP)), size=S, p=self.__gP)
		self.__genotypeList = randomGenotypes
		try:
			populationList['G'] = self.__genotypes[randomGenotypes, :]
		except ValueError:
			temp = self.__genotypes[randomGenotypes, :]
			populationList['G'] = temp.reshape(temp.size)




	def setPopulationList(self, size_x, size_y, fields):
		total_size = size_y * size_x
		population = np.zeros(total_size, dtype=fields)
		population['index'] = np.arange(size_x*size_y)
		self.randomGenotyping(populationList=population)
		print(population)
		return population

	@staticmethod
	def __setGridValuesFromFile(filename):

		coords = pd.read_csv(filename)
		return coords['x'].values, coords['y'].values

	@staticmethod
	def __setGridValues(size_x, size_y, asimetric_y = 2.0, separation_x = 0.5):
		x_coordinates = np.arange(size_x, dtype='float32')
		y_coordinates = np.arange(size_y, dtype='float32')
		x_coordinates *= separation_x
		y_coordinates *= separation_x * asimetric_y
		xx, yy = np.meshgrid(x_coordinates,y_coordinates)
		xx = xx.astype('float32')
		yy = yy.astype('float32')
		return xx.reshape(size_x*size_y), yy.reshape(size_x*size_y)

	@staticmethod
	def __createAttributeList(listAttributes, numpyType='float32'):
		attributes = list()
		for item in listAttributes:
			attributes.append((item, numpyType))
		return attributes

	def setSeed(self, patho, method = 'random', number = 1):
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
			pathoLoc = self.__loc[patho]
		except KeyError:
			raise KeyError("There is no {0} pathotype specified".format(patho))
		if method == 'random' :
			try:
				possible_hosts = np.where(self.__P['G'][:,pathoLoc] == 1)[0]
			except IndexError:
				possible_hosts = np.where(self.__P['G'][:] == 1)[0]
			try:
				seedValue = np.random.choice(possible_hosts, number)
				print("patho = {0}, seed = {1}".format(patho, seedValue))
				self.__x[seedValue, pathoLoc] = 1.0
			except ValueError:
				return


	def getPathotypes(self):
		"""
		arcadeSpots3 - arcadePopulation - getPathotypes()
		:return:

		It returns a list with the pathototypes that are included
		within the population.
		"""
		return self.__patho

	def getCoordinates(self):
		"""
		arcadeSpots3 - arcadePopulation - getCoordinates()
		:return: numpy array, numpy array

		 It returns two arrays (x,y) with the coordinates of each host
		"""
		return np.copy(self.__P['rx']), np.copy(self.__P['ry'])

	def getGenotypes(self):
		"""

		:return:
		"""
		try:
			return np.copy(self.__genotypeList)
		except AttributeError:
			raise IOError("There was some issue dealing with the genotypes")

	def getShape(self):
		return self.__sx, self.__sy

	def getAlive(self):
		"""
		arcadeSpots3 - arcadePopulation - getAlive()
		:return: numpy array

		It returns an array with the alive state of the hosts
		"""
		return self.__P['A']

	def getExposition(self):
		"""
		arcadeSpots3 - arcadePopulation - getExposition()
		:return: numpy array

		It returns an array with the exposition states of the hosts
		"""
		return self.__x

	def getInfective(self):
		"""
		arcadeSpots3 - arcadePopulation - getInfective()
		:return: numpy array

		It returns an array with the infective states of the hosts
		"""
		return self.__y
	#
	# In development
	def getTranmission(self):

		transmission = self.__y.dot(self.__C.T)
		transmission[transmission < 0] = 0.0
		return transmission

	def getDays(self):
		"""
		arcadeSpots3 - arcadePopulation - getDays()
		:return: numpy array

		 It returns an array with the exposition states of the hosts
		"""
		return self.__d

	def updateInfective(self):
		"""
		arcadeSpots3 - arcadePopulation - getDays()
		:return:

		It updates the infective states depending on the days post infection,
		the alive state, and the number of coinfections
		"""
		# This is necessary to avoid double or triple infections

		controlInfective = np.sum(self.__y, axis=1)

		self.__d += self.__x >= 1.0
		self.__y +=  ((self.__d == self.__dpi).T * (controlInfective < 2).T * self.__P['G'].T).T
		self.__y =  (self.__P['A'] * self.__y.T).T

	def updateExposition(self, I):
		"""
		arcadeSpots3 - arcadePopulation - updateExposition()
		:param I: numpy array
		:return:

		Once provided a vector that specifies the contacts between infective hosts and
		other hosts, it updates the exposition value by multiplying the contact by
		the rate of transmission of each pathotype.
		"""

		#self.__x += (self.__P['A'] *
		#			 (self.__beta * I).T).T.reshape((self.__sx*self.__sy, len(self.__patho)))

		self.__x += self.__P['A'].reshape((self.__P['A'].size,1)) * I
		# sSIP : secondary Stochastic Infection Probability
		# sSIE : secondary Stochastic Infection Event
		#sSIP = (self.__P['A']*S.T).T
		#sSIE = sSIP > np.random.rand(sSIP.shape[0], sSIP.shape[1])
		#self.__x[sSIE == True]  += 1.0
		self.__x[self.__x > 1.0] = 1.0

	def stochasticSecondaryInfections(self, S, method = 'tau_leap'):
		"""
		arcadeSpots3 - arcadePopulation - stochasticSecondaryInfections()
		:param S:
		:param method:
		:return:
		"""
		t = 0.0
		if method == 'gillespie' :
			while t < 1.0 :
				a0 =  np.sum(S, axis = 0)
				if np.all(a0 == 0) :
					break
				tau_reactions = (1/a0) * np.log(1/np.random.rand(1))
				tau = np.min(tau_reactions)
				reaction = np.argmin(tau_reactions)
				t +=  tau
				au = S / np.sum(S, axis = 0)
				k  = np.random.choice(np.arange(len(au[:,reaction])),
									  1, replace=False, p=au[:, reaction])

				self.__x[k] = 1.0
		elif method == 'tau_leap' :
			for patho in range(S.shape[1]) :
				events = self.__tau_leap_sSI(S[:, patho], 1.0)
				self.__x[events, patho] = 1.0

	@staticmethod
	def __tau_leap_sSI(probabilityVector, time_step):
		"""

		:param probabilityVector:
		:param time_step:
		:return:
		"""
		if len(probabilityVector.shape) == 1 :
			oA = np.sum(probabilityVector)
			if oA == 0.0 :
				return np.array([], dtype='int32')
			events             = np.arange(probabilityVector.size)
			n                  = np.random.poisson(oA*time_step)
			chosen_events      = np.random.choice(events, n, replace = False,
												  p = probabilityVector/oA)
			return chosen_events
		else:
			raise ValueError("probability Vector must be unidimensional")
	def updateAlive(self):
		"""
		arcadeSpots3 - arcadePopulation - updateAlive()
		:return:

		Updates the primary inoculum, which gets increased when hosts die, and it modifies
		the values of the alive state
		"""

		ref_vals = np.array([self.__C[i,i] for i in range(len(self.__patho))])
		ref_vals = ref_vals.reshape(1, len(self.__patho))
		self.__w += (self.__d == self.__mortality)*(self.getTranmission() / ref_vals)*self.__omega
		self.__P['A'][np.sum(self.__d >= self.__mortality, axis=1) >= 1.0] = False

	def updateInoculum(self):
		"""
		arcadeSpots3 - arcadePopulation - updateAlive()
		:return:

		Updates the primary inoculum upon degradation
		"""
		self.__w *= self.__k_omega

	def primaryInfection(self):
		"""
		arcadeSpots3 - arcadePopulation - primaryInfection()
		:return:

		Updates the exposition depending on the linear relationship between
		the concentration of primary inoculum in the soil.
		"""
		base_prob = 5.4e-4
		randomValues = (np.random.rand(self.__sx*self.__sy, len(self.__patho)) < (self.__w/1000000)*base_prob)
		primaryInfections = (self.__P['A'] * randomValues.T).T
		self.__x         += primaryInfections
		self.__x[self.__x > 1.0] = 1.0

	def nextCrop(self):
		"""
		arcadeSpots3 - arcadePopulation - nextCrop()
		:return:

		Resets all the values but the primary inoculum
		"""
		self.__w[:,:] += self.__omega * ((self.__P['A'] == True) * (self.__y == True).T).T
		self.__x[:,:] = 0.0
		self.__y[:,:] = False
		self.__d[:,:] = 0
		self.__P['A'][:] = True

	def modifyPopulation(self, ind, parameter, value):
		try :
			self.__P[parameter]
		except KeyError:
			raise IOError("there is not such an attribute %s" % parameter)
		try :
			self.__P[parameter][ind]
		except IndexError:
			raise IOError("there is not such an individual %d" % ind)
		self.__P[parameter][ind] = value

	def plotPopulation(self, patho, time=None, show=True):
		try:
			pathoLoc = self.__loc[patho]
		except KeyError:
			raise KeyError("There is no such a {0} pathotype".format(patho))
		exposition = self.__x[:,pathoLoc].reshape((self.__sx, self.__sy))
		infective  = self.__y[:,pathoLoc].reshape((self.__sx, self.__sy))
		alive      = self.__P['A'].reshape((self.__sx, self.__sy))
		primaryInoc= self.__w[:,pathoLoc].reshape((self.__sx, self.__sy))
		#plt.figure(figsize=(10,5))
		plt.subplot(221)
		sns.heatmap(exposition, cmap='viridis', xticklabels=False, yticklabels=False)
		plt.title('Exposition')
		plt.subplot(222)
		sns.heatmap(infective, cmap='viridis', xticklabels=False, yticklabels=False)
		plt.title('infective')
		plt.subplot(223)
		sns.heatmap(alive, cmap='viridis', xticklabels=False, yticklabels=False)
		plt.title('alive')
		plt.subplot(224)
		sns.heatmap(primaryInoc, cmap='viridis', xticklabels=False, yticklabels=False)
		plt.title('primary Inoculum')
		if show:
			plt.show()
		else:
			plt.savefig("patho_{0}_day_{1}.png".format(patho, time), dpi=300)
			plt.close()
	def spatialMap(self, patho, time, crop, show=False, filename = False):
		try :
			pathoLoc = self.__loc[patho]
		except KeyError:
			raise KeyError("There is no such a {0} pathotype".format(patho))
		sns.set_style('darkgrid')
		axis = [
			self.__P['rx'].min() - 1,
			self.__P['rx'].max() + 1,
			self.__P['ry'].min() - 1,
			self.__P['ry'].max() + 1

		]
		plt.suptitle('time = %d' % time)
		##------------------------------------------------------------#
		## exposition
		##------------------------------------------------------------#
		plt.subplot(221)
		plt.title('exposition')
		ax=plt.scatter(self.__P['rx'], self.__P['ry'], c=self.__x[:,pathoLoc],
					cmap='viridis')
		plt.ylim(self.__P['ry'].min() - 1, self.__P['ry'].max() + 1)
		plt.xlim(self.__P['rx'].min() - 1, self.__P['rx'].max() + 1)
		plt.colorbar(ax)
		##------------------------------------------------------------#
		## exposition
		##------------------------------------------------------------#
		plt.subplot(222)
		plt.title('infective')
		ax =plt.hexbin(self.__P['rx'][self.__y[:,pathoLoc] == True],
						self.__P['ry'][self.__y[:, pathoLoc] == True],
						cmap='viridis', extent=axis, gridsize=20)
		plt.colorbar(ax)
		plt.ylim(self.__P['ry'].min() - 1, self.__P['ry'].max() + 1)
		plt.xlim(self.__P['rx'].min() - 1, self.__P['rx'].max() + 1)
		##------------------------------------------------------------#
		## alive
		##------------------------------------------------------------#
		plt.subplot(223)
		plt.title('mortality')
		ax = plt.hexbin(self.__P['rx'][self.__P['A'] == False],
						 self.__P['ry'][self.__P['A']== False],
						 cmap='viridis', extent=axis, gridsize=20)
		plt.colorbar(ax)
		plt.ylim(self.__P['ry'].min() - 1, self.__P['ry'].max() + 1)
		plt.xlim(self.__P['rx'].min() - 1, self.__P['rx'].max() + 1)
		##------------------------------------------------------------#
		## primary inoc
		##------------------------------------------------------------#
		plt.subplot(224)
		plt.title('primary inoculum')
		ax = plt.scatter(self.__P['rx'],self.__P['ry'], c=np.log10(self.__w[:,pathoLoc]),
					cmap='viridis')
		plt.ylim(self.__P['ry'].min() - 1, self.__P['ry'].max() + 1)
		plt.xlim(self.__P['rx'].min() - 1, self.__P['rx'].max() + 1)
		plt.colorbar(ax)
		if filename :
			plt.savefig(filename, dpi=300)
		else:
			plt.savefig('space_map_patho_%s_crop_%d_time_%d.png' % (patho, crop, time), dpi=300)
		if show:
			plt.show()
		else:
			plt.close()
	def plotInfectiveAreas(self, patho1, patho2):
		sns.set_style('white')
		infectiveAreas = np.zeros((self.__sx, self.__sy,3))
		patho1_loc = self.__loc[patho1]
		patho2_loc = self.__loc[patho2]
		infectiveAreas[:, :, 0] = self.__y[:, patho1_loc].reshape((self.__sx, self.__sy))
		infectiveAreas[:, :, 1] = self.__y[:, patho2_loc].reshape((self.__sx, self.__sy))
		plt.imshow(infectiveAreas, interpolation=None)
		plt.grid(b=None)
		plt.show()

	def getStatistics(self, patho):
		statistics = np.zeros(1, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('alive',      'float32'),
			('inoculum',   'float32')
		])
		total_pop = self.__sy * self.__sx
		try:
			pathoLoc = self.__loc[patho]
		except KeyError:
			raise KeyError("There is no such a {0} pathotype".format(patho))
		statistics['exposition'] = np.sum(self.__x[:,pathoLoc]) / total_pop
		statistics['infective']  = np.sum(self.__y[:,pathoLoc]) / total_pop
		statistics['alive']      = np.sum(self.__P['A']) / total_pop
		statistics['inoculum']   = np.sum(self.__w[:,pathoLoc]) / (total_pop*1000)
		return statistics



	def __getIndex(self, k):
		i = int(k/self.__sx)
		j = k - i*self.__sx
		return i,j

	def __linearMortalityModel(self):
		x = np.random.rand(self.__sx * self.__sy, len(self.__patho))
		tau = np.log(x) / (-self.__sigma)
		tau = np.round(tau, 0)
		tau += self.__mu
		tau += self.__dpi
		return tau

	def __gaussianMortalityModel(self):

		x = np.random.rand(self.__sx * self.__sy, len(self.__patho))
		u = stats.norm.ppf(x)
		return u*self.__sigma + self.__mu + self.__dpi


	def __logNormalMortalityModel(self):
		x = np.random.rand(self.__sx * self.__sy, len(self.__patho))
		u = stats.norm.ppf(x)
		return np.exp(u * self.__sigma + self.__mu) + self.__dpi

	def __sanityModel(self):
		hazard = self.__mu * np.exp(self.__sigma)
		x = np.random.rand(self.__sx * self.__sy, len(self.__patho))
		u = np.log(x)*(1 / -hazard)
		return u.astype(int) + self.__dpi