# arcadePopulation
import numpy as np

class arcadePopulation:
	fields = [
		('rx', 'float32'),
		('ry', 'float32'),
		('x', 'float32'),
		('y', 'b1'),
		('d', 'uint64'),
		('w', 'float32'),
		('A', 'b1'),
		('G', 'int8')

	]
	dpi = 15
	mortality = 30
	def __init__(self, size_x, size_y, genotypes = 'P0', scattering='Regular'):
		print("arcadePopulation instance defined")
		print("setting a population of %d hosts" %(int(size_x*size_y)))
		self.__sx = size_x
		self.__sy = size_y
		self.__P = self.setPopulationList(size_x, size_y)
		xx_coords, yy_coords = self.__setGridValues(size_x, size_y, res=0.5, eliptical=1.0)
		self.__P['rx'] = xx_coords
		self.__P['ry'] = yy_coords
		self.__P['A']  = np.ones(size_y*size_x)
		print("population set")

	def setPopulationList(self, size_x, size_y):
		total_size = size_y * size_x
		population = np.zeros(total_size, dtype=self.fields)
		return population

	@staticmethod
	def __setGridValues(size_x, size_y, res = 1.0, eliptical = 1.0):
		x_coordinates = np.arange(size_x)
		y_coordinates = np.arange(size_y)
		xx, yy = np.meshgrid(x_coordinates,y_coordinates)
		xx = xx.astype('float32')
		yy = yy.astype('float32')
		#xx *= res
		#yy *= res
		#yy *= eliptical
		return xx.reshape(size_x*size_y), yy.reshape(size_x*size_y)

	def getCoordinates(self):
		return self.__P['rx'], self.__P['ry']

	def getAlive(self):
		return self.__P['A']

	def getExposition(self):
		return self.__P['x']

	def getInfective(self):
		return self.__P['y']

	def getDays(self):
		return self.__P['d']

	def setSeed(self, seedValue):
		self.__P['x'][seedValue] = 1.0

	def updateInfective(self):
		days = np.copy(self.__P['d'])
		self.__P['d'][self.__P['x'] >= 1.0] += 1
		self.__P['y'][self.__P['d'] >= self.dpi] = True
		self.__P['y'] = self.__P['A'] * self.__P['y']
		postDays = np.copy(self.__P['d'])
		if any(postDays < days):
			raise RuntimeError("no way man!")
	def updateExposition(self, deltaX):
		self.__P['x'] += deltaX*self.__P['A']
		self.__P['x'][self.__P['x'] > 1.0] = 1.0

	def updateAlive(self):

		self.__P['A'][self.__P['d'] > self.mortality] = False
		self.__P['w'][self.__P['d'] == self.mortality + 1] = 1000


	def updateInoculum(self):
		self.__P['w'] *= (1 - 34e-3)

	def primaryInfection(self):
		primaryInfections = self.__P['A']*(np.random.rand(self.__sx*self.__sy) < (self.__P['w']/10000))
		self.__P['x'] += primaryInfections
		self.__P['x'][self.__P['x'] > 1.0] = 1.0

	def nextCrop(self):
		self.__P['x'][:] = 0.0
		self.__P['y'][:] = False
		self.__P['A'][:] = True
		self.__P['d'][:] = 0



	def plotPopulation(self):
		import matplotlib.pyplot as plt
		import seaborn as sns
		exposition = self.__P['x'].reshape((self.__sx, self.__sy))
		infective  = self.__P['y'].reshape((self.__sx, self.__sy))
		alive      = self.__P['A'].reshape((self.__sx, self.__sy))
		primaryInoc= self.__P['w'].reshape((self.__sx, self.__sy))
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
		plt.show()

	def getStatistics(self):
		statistics = np.zeros(1, dtype=[
			('exposition', 'float32'),
			('infective',  'float32'),
			('alive',      'float32'),
			('inoculum',  'float32')
		])
		total_pop = self.__sy * self.__sx
		statistics['exposition'] = np.sum(self.__P['x']) / total_pop
		statistics['infective'] = np.sum(self.__P['y']) / total_pop
		statistics['alive'] = np.sum(self.__P['A']) / total_pop
		statistics['inoculum'] = np.sum(self.__P['w']) / (total_pop*1000)
		return statistics





