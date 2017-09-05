# arcadePopulation
import numpy as np
import cProfile
def logistic(t, k, a):
	return np.exp(k*t - a)/(1+np.exp(k*t - a))

class arcadePopulation:
	commonFields = [('index', 'uint32'), ('rx', 'float32'), ('ry', 'float32'), ('A', 'b1'), ('G', 'int8', 2), ('I', 'b1')]

	dpi = 5
	mortality = 35
	genotypes = np.array([[1,1],[0,1],[0,0]])
	def __init__(self, size_x, size_y, parametersDict):

		beta = list()
		primaryInoc = list()
		alpha = list()
		pathotypesList = sorted(parametersDict['pathotypes'])
		self.__patho = pathotypesList
		self.__loc   = dict()
		i = 0
		for item in pathotypesList:
			beta.append(parametersDict['beta'][item])
			primaryInoc.append(parametersDict['P'][item])
			alpha.append(parametersDict['alpha'][item])
			self.__loc[item] = i
			i += 1
		self.__beta = np.array(beta)
		self.__pI   = np.array(primaryInoc)
		self.__alpha = np.array(alpha)
		self.__C    = parametersDict['C']
		print("arcadePopulation instance defined")
		print("setting a population of %d hosts" %(int(size_x*size_y)))

		self.__sx = size_x
		self.__sy = size_y
		size = size_x * size_y

		self.__P = self.setPopulationList(size_x, size_y, self.commonFields)
		xx_coords, yy_coords = self.__setGridValues(size_x, size_y)

		self.__P['rx'] = xx_coords
		self.__P['ry'] = yy_coords
		self.__P['A'][:] = True

		print("common attributes set")

		self.__x = np.zeros((size, len(pathotypesList)))
		self.__y = np.zeros((size, len(pathotypesList)))
		self.__d = np.zeros((size, len(pathotypesList)))
		self.__w = np.zeros((size, len(pathotypesList)))

		print("population set")

	# STATIC METHODS
	## setPopulationList
	## setGridValues
	## createAttributeList

	def randomGenotyping(self, populationList):
		S = self.__sx*self.__sy
		randomGenotypes = np.random.choice((0,1,2), size=S, p=(0.8, 0.1, 0.1))
		populationList['G'] = self.genotypes[randomGenotypes,:]


	def setPopulationList(self, size_x, size_y, fields):
		total_size = size_y * size_x
		population = np.zeros(total_size, dtype=fields)
		population['index'] = np.arange(size_x*size_y)
		self.randomGenotyping(populationList=population)
		return population

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

	def setSeed(self, seedValue, patho):
		"""
		arcadeSpots3 - arcadePopulation - setSeed()
		:param seedValue:
		:param patho:
		:return:

		It sets the exposition value of a given position to 1.0,
		 setting the beginning of an epidemic. It must be specified
		 for each pathotype of the epidemic
		"""
		try:
			pathoLoc = self.__loc[patho]
		except KeyError:
			raise KeyError("There is no {0} pathotype specified".format(patho))
		self.__x[seedValue, pathoLoc] = 1.0

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
		self.__y +=  ((self.__d == self.dpi).T *(controlInfective < 1)).T * self.__P['G']
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

		self.__x += (self.__P['A'] *
					 (self.__beta * I).T).T.reshape((self.__sx*self.__sy, len(self.__patho)))

		# sSIP : secondary Stochastic Infection Probability
		# sSIE : secondary Stochastic Infection Event
		#sSIP = (self.__P['A']*S.T).T
		#sSIE = sSIP > np.random.rand(sSIP.shape[0], sSIP.shape[1])
		#self.__x[sSIE == True]  += 1.0
		self.__x[self.__x > 1.0] = 1.0

	def stochasticSecondaryInfections(self, S, method = 'gillespie'):
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
		self.__w += (self.__d == self.mortality)*(self.__y == True)*self.__pI
		self.__P['A'][np.sum(self.__d >= self.mortality, axis=1) >= 1.0] = False

	def updateInoculum(self):
		"""
		arcadeSpots3 - arcadePopulation - updateAlive()
		:return:

		Updates the primary inoculum upon degradation
		"""
		self.__w *= self.__alpha

	def primaryInfection(self):
		"""
		arcadeSpots3 - arcadePopulation - primaryInfection()
		:return:

		Updates the exposition depending on the linear relationship between
		the concentration of primary inoculum in the soil.
		"""
		randomValues = (np.random.rand(self.__sx*self.__sy, len(self.__patho)) < (self.__w/10000))
		primaryInfections = (self.__P['A'] * randomValues.T).T
		self.__x         += primaryInfections
		self.__x[self.__x > 1.0] = 1.0

	def nextCrop(self):
		"""
		arcadeSpots3 - arcadePopulation - nextCrop()
		:return:

		Resets all the values but the primary inoculum
		"""
		self.__w[:,:] += self.__pI*((self.__P['A'] == True)*(self.__y == True).T).T
		self.__x[:,:] = 0.0
		self.__y[:,:] = False
		self.__d[:,:] = 0
		self.__P['A'][:] = True

	def plotPopulation(self, patho, time=None, show=True):
		import matplotlib.pyplot as plt
		import seaborn as sns
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
	def plotInfectiveAreas(self, patho1, patho2):
		import matplotlib.pyplot as plt
		import seaborn as sns
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



