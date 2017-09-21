#!/home/charizard/anaconda3/bin/ipython3

import json
import arcadeSimulator as aS
import sys


def parseArgumentsFile(fileName):
	f = open(fileName)
	json_structure = json.load(f)
	try :
		assert json_structure['global_parameters']
		assert json_structure['specific_parameters']
		assert json_structure['interspecific_parameters']
	except KeyError :
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

try :
	parameters_file = sys.argv[1]
except IndexError:
	print("missing parameters file")
	quit()
params = parseArgumentsFile('params.json')
if params['metaparameters']['verbose'] :
	printReport(params)
aaa = None
for i in range(params['metaparameters']['simulations']):
	aaa = aS.arcadeSimulator(params, sim=i)
	aaa.simulate()
	if i == 0:
		aaa.printHeader()
	aaa.saveReport(header=True)
aS.arcadeOutput(params, aaa.getPopulation())
