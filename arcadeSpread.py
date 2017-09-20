#!/home/bruno/anaconda3/bin/python3 -W ignore
import numpy as np

example_parameters = dict(
	size = 500,
	time = 100,
	diffusion_constant = 0.5,
	time_step = 1/(4*0.5) - 0.1,
	reaction_constant  = 0.1,
	maximum_density = 4.0,
	initial_seeds = 4,
	seasons = 5,
	outfile = 'wild_population_coordinates.csv'

)

class arcadeSpread :
	def __init__(self, parameters):
		try :
			self.__s = parameters['size']
			self.__time = parameters['time']
			self.__time_step = parameters['time_step']
			self.__diffusion = parameters['diffusion_constant']
			self.__reaction  = parameters['reaction_constant']
			self.__maximum_density = parameters['maximum_density']
			self.__initial_seeds   = parameters['initial_seeds']
			self.__seasons = parameters['seasons']
			self.__outfile = parameters['outfile']
			coordinates_list = np.random.rand(parameters['initial_seeds'], 2)
			coordinates_list *= parameters['size']
			self.__coords = coordinates_list

		except KeyError :
			raise IOError("parameters where not correctly specified")
	@staticmethod
	def laplacian(gg):

		x_left  = gg[1:-1, :-2]
		x_right = gg[1:-1, 2:]
		x_up    = gg[:-2, 1:-1]
		x_down  = gg[2:, 1:-1]
		x_center= gg[1:-1, 1:-1]

		return (x_left + x_right + x_up + x_down) - 4*x_center


	def events(self, gg):

		probability_sum = np.sum(gg)
		events_number = np.random.poisson(lam=probability_sum, size=1)
		flattened_probability = gg.reshape(gg.size)
		flattened_probability /= np.sum(flattened_probability)
		chosen_elements = np.random.choice(np.arange(flattened_probability.size),
										   size= events_number,
										   p=flattened_probability
										   )

		new_coordinates = np.zeros((int(events_number), 2))


		new_coordinates[:, 0] = np.array([int(elem / self.__s)  for elem in chosen_elements])
		new_coordinates[:, 1] = np.array([elem % self.__s for elem in chosen_elements])

		new_coordinates += np.random.rand(int(events_number), 2)
		self.__coords = np.concatenate((self.__coords, new_coordinates))


	def saturation(self):
		coords_histomap, axis_x, axis_y = np.histogram2d(self.__coords[:,0],
														 self.__coords[:,1],
														 bins = self.__s, range=[[0, self.__s]]*2)
		del axis_x, axis_y
		saturation_matrix = (self.__maximum_density - coords_histomap) / self.__maximum_density
		saturation_matrix[saturation_matrix < 0] = 0
		return saturation_matrix


	def simulate(self):


		histo_map, axis_x, axis_y = np.histogram2d(self.__coords[:,0],
												   self.__coords[:,1],
												   bins=self.__s, range=[[0, self.__s]]*2)
		del axis_x, axis_y
		grid = np.zeros((self.__s, self.__s))

		for season in range(self.__seasons) :
			for i in range(self.__time) :
				difusion_term = self.laplacian(grid) * self.__diffusion
				reaction_term =  histo_map[1:-1, 1:-1] * self.__reaction
				grid[1:-1, 1:-1] += self.__time_step*(difusion_term+reaction_term)
			self.events(grid*self.saturation())
			histo_map, axis_x, axis_y = np.histogram2d(self.__coords[:, 0],
													   self.__coords[:, 1],
													   bins=self.__s, range=[[0, self.__s]] * 2)
			grid = np.zeros((self.__s, self.__s))

	def write(self):

		f = open(self.__outfile, 'w')
		f.write("x,y\n")
		for line in range(self.__coords.shape[0]):
			f.write("%8.3f,%8.3f\n" % (self.__coords[line,0], self.__coords[line, 1]))
		f.close()







