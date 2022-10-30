This repository contains the scripts for the model presented in "Using a Plant Hydraulic Model to Design More Resilient Rehabilitated Landscapes in Arid Ecosystems". The packages needed for running the scripts are listed in requirements.txt. The two main scripts are run_parameterization.py and run_traectory.py. The parameterization script runs the model over the supplied parameter sets and saves the performance metrics of the model run as compared to the obervation data from Lamoureux et al. (2018), which is contained in observation_data.py. The traectory script runs the model using the supplied top performing parameter sets under historical precipitation trajectories saved in Ninghan_Station.csv (originally downloaded from Station ID 7068, http://www.bom.gov.au/climate/data/).

Both run_parameterization.py and run_trajectory.py are set up to run multiple parameter sets in parallel, and use ipyparallel for starting and managing tasks across the engines.

The model is run using a variety of functions defined in utility_functions.py, which references the model framework and functions in soil_plant_model.py.

The files params_constants.py, params_soil.py, and WA_hydraulic_traits_SebData.xls contain various parameters and constants that are fed into the model. Please note that some of these values are overwritten by those defined in the individual parameter sets.

Questions or comments please contact Jeannie Wilkening (jvwilkening@berkeley.edu). Last Updated: October 2022