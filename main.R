
# models are in ./Models_functions

#load multi-model wrapper functions
source("./Models_functions/multi_model.R")

#TODO: Q has now a transient temperature sensitivity but linear, maybe it can be updated
#TODO: Yasso 07 and 15  are initialized via spinup


#For testing, run
source("./Models_functions/synthetic_data_template.R")
#it contains both a function to generate synthetic data to use as a data requirement template
# and a test run 