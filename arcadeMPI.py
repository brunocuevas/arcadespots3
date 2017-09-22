from mpi4py import MPI
import json
import arcadeSimulator as aS
import sys
import numpy as np

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

def distribute(task_to_split, dist_size):
	tasks_vector = np.arange(task_to_split)
	rest  = task_to_split % (dist_size - 1)
	distributed_tasks = np.split(tasks_vector[:-rest], dist_size-1)
	distributed_tasks.append(tasks_vector[-rest:])
	return distributed_tasks

input_arguments = sys.argv[1]
arguments = None
try:
	arguments = parseArgumentsFile(input_arguments)
except FileExistsError :
	print("file not found")
	quit()


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0 :
	simulation_batch = arguments['metaparameters']['simulations']
	tasks = distribute(simulation_batch, size)
elif rank != 0 :
	tasks = None

tasks = comm.scatter(tasks, root=0)
control_header = True
for sim in tasks :
	aaa = aS.arcadeSimulator(arguments, sim=sim)
	if rank == 1 and control_header :
		aaa.printHeader()
		control_header = False
	aaa.simulate()
	aaa.saveReport(header=False)



