import numpy as np

class displayer:
	def __init__(self, crops, cropDays, smallSpace=5, largeSpace = 30):
		self.__crops = crops
		self.__cD    = cropDays
		self.__currentCrop = 1
		self.__currentDay  = 1
		self.__sS = smallSpace
		self.__lS = largeSpace

	def start(self):
		print("START\n %4d ." % 0, end="")
	def update(self):
		self.__currentDay += 1
		if self.__currentDay < self.__cD:
			if self.__currentDay % self.__lS == 0:
				eol = "\n %4d " % self.__currentDay
			else:
				if self.__currentDay % self.__sS == 0 :
					eol = " "
				else:
					eol = ""
			print(".", end=eol)
		elif self.__currentDay == self.__cD:
			print(".", end="\n")
		else:
			if self.__currentCrop < self.__crops:
				print("next crop")
				self.__currentDay = 1
				self.__currentCrop += 1
				print(" %4d ." % 0, end="")
			else:
				print("\nEND")

class writer :

	def __init__(self, folder_name):
		self.__folder_name = folder_name
	def getFolderName(self):
		return self.__folder_name
	def write(self):
		import os
		os.makedirs(self.__folder_name)


def save_statistics(filename, statistics, sim=0):

	with open(filename, 'w') as f:

		f.write('crop, time, exposition, infective, alive, inoculm, sim\n')
		for line in statistics:
			f.write(
				'{:4d}, {:4d}, {:12.4e}, {:12.4e}, {:12.4e}, {:12.4e}, {:4d}\n'.format(
					line['crop'], line['time'], line['exposition'], line['infective'],
					line['alive'], line['inoculum'], sim
				)
			)
