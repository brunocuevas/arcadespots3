import json
import arcadeSimulator as aS
import sys


def parse_arguments_file(file_name):
	f = open(file_name)
	json_structure = json.load(f)
	terms = ['global_parameters', 'specific_parameters', 'interspecific_parameters', 'metaparameters']
	for batch in json_structure:
		if not all([k in batch.keys() for k in terms]):
			raise IOError("parameters were not correctly specified. global, specific and interspecific")
	return json_structure


def print_report(parameters):
	print("_" * 25)
	print("\tarcadeSpots3, by Bruno Cuevas")
	print("_" * 25)
	print("\tOutfile = {0}".format(parameters['metaparameters']['outfile']))
	print("\tNumber of simulations = {0}".format(parameters['metaparameters']['simulations']))
	print("\tNumber of pathotypes in the simulation = {0}".format(len(parameters['global_parameters']['pathotypes'])))
	print("\tMortality model = {0}".format(parameters['metaparameters']['mortality']))
	print("\tPlacement model = {0}".format(parameters['metaparameters']['placement']))
	print("\tCrops = {0}".format(parameters['global_parameters']['crops']))


if __name__ == '__main__':
	parameters_file = None
	try:
		parameters_file = sys.argv[1]
	except IndexError:
		print("missing parameters file")
		quit()
	params_batch = parse_arguments_file(parameters_file)

	for params in params_batch:
		if params['metaparameters']['verbose']:
			print_report(params)
		aaa = None
		for i in range(params['metaparameters']['simulations']):
			aaa = aS.ArcadeSimulator(params, sim=i, verbose=True)
			aaa.simulate()
			aaa.save_report()
