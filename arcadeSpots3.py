#!/home/charizard/anaconda3/bin/ipython3

import json
import arcadeSimulator as aS
import sys


def parseArgumentsFile(fileName):
	f = open(fileName)
	json_structure = json.load(f)
	terms = ['global_parameters', 'specific_parameters', 'interspecific_parameters', 'metaparameters']
	for batch in json_structure:
		if not all([k in batch.keys() for k in terms]):
			raise IOError("parameters were not correctly specified. global, specific and interspecific")
	return json_structure

def printReport(parameters):
	print("_" * 25)
	print("\tarcadeSpots3, by Bruno Cuevas")
	print("_" * 25)
	print("\tOutfile = {0}".format(parameters['metaparameters']['outfile']))
	print("\tNumber of simulations = {0}".format(parameters['metaparameters']['simulations']))
	print("\tNumber of pathotypes in the simulation = {0}".format(len(parameters['global_parameters']['pathotypes'])))
	print("\tMortality model = {0}".format(parameters['metaparameters']['mortality']))
	print("\tPlacement model = {0}".format(parameters['metaparameters']['placement']))
	print("\tCrops = {0}".format(parameters['global_parameters']['crops']))

parameters_file = None
try :
	parameters_file = sys.argv[1]
except IndexError:
	print("missing parameters file")
	quit()
params_batch = parseArgumentsFile(parameters_file)

for params in params_batch :
	if params['metaparameters']['verbose'] :
		printReport(params)
	aaa = None
	for i in range(params['metaparameters']['simulations']):
		aaa = aS.arcadeSimulator(params, sim=i, verbose=False)
		aaa.simulate()
		aaa.saveReport()
	aS.arcadeOutput(params, aaa.getPopulation())

