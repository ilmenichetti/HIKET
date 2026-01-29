
# models are in ./Models_functions

#load multi-model wrapper functions
source("./Models_functions/multi_model.R")

# TODO: models still give a bit too different results, find ways for testing
# TODO: Q has now a transient temperature sensitivity but linear, maybe it can be updated
# TODO: Yasso 07 and 15  are initialized via spinup
# TODO documentation is still inconsistent, needs more worl

#For testing, run
source("./Models_functions/synthetic_data_template.R")
#it contains both a function to generate synthetic data to use as a data requirement template
# and a test run 