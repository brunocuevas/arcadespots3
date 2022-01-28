# arcadeSpots 3
![arcadeSpots3_gif](./doc/epidemic.gif)

Providing epidemiological models for plant-virus pathosystems

## Summary

![arcadeSpots3_cover](./doc/arcadeSpots3_cover.png)

Plant pathogens are important due to their socioeconomical impact.
Since experimentation with pathogens at environmental level is
technically difficult and expensive, epidemiologial 
models are a reasonable approach that allow to study the main
features of the pathogen-host system.

arcadeSpots3 allows to perform agent based simulations where
host and pathogen populations are genetically diverse, and it
includes features as spatial location of individuals, virulence,
long and short kernel spread terms, and seasonality.

## Requirements

arcadeSpots3 is implemented in python3 and it requires
some standard libraries that can be found in the anaconda
distribution, as Numpy, Scipy, Matplotlib or Seaborn.

To run the parallel implementation of the program, both
mpi4py python library and mpi library are needed. Consult
your system administrator.

## Usage

arcadeSpots3 requires specific parameters for each simulation,
that must be stored in the form of a JSON dictionary. The parameters employed
in the simulations of the published article are available at:

	parameters/_2022/
 
![arcadeSpots3_parameters](./doc/arcadeSpots3_inputParameters.png)

Once parameters are specified, it is possible to
run simulations using just

    python arcadeSpots3.py sample_parameters.json
    
If more than one processor is available (for instance,
let's consider tweenty processors), then a parallel
version of the program can be executed, using

    mpirun -np 20 arcadeMPI.py sample_parameters.json
    
Both programs produce an output that can be analyzed in R,
Python or Excel.

## FAQ

Bruno Cuevas Zuviría - Centro de Biotecnología y Genómica de Plantas (UPM-INIA), Madrid, Spain
bruno.czuviria@upm.es
