# How to create simulation input files using the generator script?


In the folder "sims", you can find all the input text files used to create the 
simulation input files for all the experiments used in our study.


- "# #" is used at the beginning of the file to name the generated input file, 
e.g.: #01_experiment#.
- "@ @" sets the parameter configuration in the main input simulation file,
 e.g.: @"sim_rep", (30)@.
- "{ }" is used to include a change to a specific parameter, e.g.: 
{"follicles", (9.97)}
- "< >" is used to create a list of parameter combinations. Each parameter can 
have a list of values using "( )" or a list of values by specifying a range 
with "[ ]" (e.g., [min max step], [5 35 1])

Note 1: parameter names are within quotation marks
Note 2: elements from lists are separated by a whitespace.

The input file for the sensitivity generator differs slightly from the 
experiment generator:

- "$ $" is used to set the percentage changes (e.g. $5 10$).
- "< >" to include here the names of the parameters to simulate the percentage 
changes (e.g., <carrying_capacity init_cones init_seeds long_term_rain>)

Note 3: names of the parameters between "< >" are without quotation marks in 
the sensitivity generator.

