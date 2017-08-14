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

