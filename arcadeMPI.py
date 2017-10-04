from mpi4py import MPI
import json
import arcadeSimulator as aS
import sys
import numpy as np

def parseArgumentsFile(fileName):
	f = open(fileName)
	json_structure = json.load(f)
	terms = ['global_parameters', 'specific_parameters', 'interspecific_parameters', 'metaparameters']
	for batch in json_structure:
		if not all([k in batch.keys() for k in terms]):
			raise IOError("parameters were not correctly specified. global, specific and interspecific")
	return json_structure

def distribute(task_to_split, dist_size):
	max_div = dist_size * int(task_to_split / dist_size)
	task_order  =  np.arange(task_to_split)
	tasks_list = np.split(task_order[:max_div], dist_size)
	j = 0
	for i in range(max_div, task_to_split) :
		tasks_list[j] = np.append(tasks_list[j], task_order[i])
		j += 1
	return tasks_list

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

#for batch in arguments:
batch = arguments[0]
tasks = None
if rank == 0 :
	simulation_batch = batch['metaparameters']['simulations']
	tasks = distribute(simulation_batch, size)

tasks = comm.scatter(tasks, root=0)
for sim in tasks :
	aaa = aS.arcadeSimulator(batch, sim=sim)
	aaa.simulate()
	aaa.saveReport()



