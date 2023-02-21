# MetaSqueeze model

MetaSqueeze is a mechanistic, process-based, spatially-explicit metapopulation model used in the following article: 

**Souto-Veiga, R., Groeneveld, J., Enright, N. J., Fontaine, J. B., & Jeltsch, 
F. (20XX). Climate change may shift metapopulations towards unstable 
source-sink dynamics in a fire-killed, serotinous shrub**

Corresponding author:

Rodrigo Souto-Veiga 

rsoutoveiga@uni-potsdam.de 

ORCID: https://orcid.org/0000-0001-8639-620X

The TRACE document provides a detailed description and supporting evidence of the model.

## 

## Requirements

The model was developed in C++ 17. The gcc compiler must be 9.2 or above.

If the gcc version is older, between 8.2 and 9.1, then the random seed generator must be taken by clock time or another way, not with std::random_device. Random_device shows deterministic value in older versions (older than 9.2). Also, if the gcc is below 8.2, the std library "filesystem" is unavailable. This library is used mainly in the class Output (probably there are also other places).

The model was developed principally using Qt for open source development under GNU General Public License v3.0.

[Open Source Development | Open Source License | Qt](https://www.qt.io/download-open-source)

## Setup

1. Install CMake and Make in Linux

2. Call "make" in the terminal to build the program (i.e., MetaSqueeze model) in debug and/or release mode.



## Project directory

* build
  
  * debug
  
  * release

* data
  
  * in (input simulation files)
  
  * out (output simulation results)

* src (header and source cpp files)

* tools
  
  * data_analysis_visualization
  
  * experiment_generator
  
  * parallel
  
  * README.md



## Model run

- The model can be run using the sh file "functions.sh"

- The sh file "job_metasqueeze.sh" was used to submit jobs in a cluster and run multiple simulations simultaneously (parallel).

- The main input simulation folder is in **data/in/sim**. The rest of the folders in data/in are a group of related parameters. 



### Data input files



Description of input simulation files in **data/in**

- /sim: in this folder are located the main input simulation files.

- /study_rep: Files contain a list of seeds for generating random numbers in order to create the metapopulation (i.e., random patches). This file is read from the main input file, i.e., /sim.

- /scenarios: The main parameters that define the simulation. This file is read from the main input file, i.e., /sim.

- /study_size: the parameter that define the study area.

- /metapopulation: File containes the (sub)populations characteristics. This file is read from a scenario file.

- /climate: Climate scenario (i.e., baseline or current). This file is read from a scenario file.

- /fire: Fire parameters. This file is read from a scenario file.

- /dispersal: Dispersal parameters. This file is read from a scenario file.

- /species: Demographic parameters. This file is read from a scenario file.

- /habitat_quality: Habitat quality of the patches. This file is read from a scenario file.

- /plant_types: Intraspecific variability of plants (plant performance class). This file is read from a scenario file.

- /fuzzy_sets: Definitions of fuzzy sets and their membership functions for flower count data. This file is read from a scenario file.

- /flower_distributions: Probability density functions for each membership function and plant performance class. This file is read from a scenario file.



For a more detailed description of input parameters, please refer to the cpp source code in Parameter_reader class, or to the TRACE document.
