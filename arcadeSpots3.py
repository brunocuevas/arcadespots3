#!/home/charizard/anaconda3/bin/ipython3

import json
import arcadeSimulator as aS



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

params = parseArgumentsFile('params.txt')
#for sim in range(params['metaparameters']['simulations']) :
sim = 1
aaa = aS.arcadeSimulator(params, sim=sim)
aaa.simulate()
#aaa.plotStatistics()
aaa.printHeader('trials_1')
aaa.saveReport('trials_1', header=True)
