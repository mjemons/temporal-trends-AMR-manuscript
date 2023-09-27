## This is a masterfile for the temporal trends project. 

set.seed(1234)

#### 0. Load Utilities #####

source("code/Utils.R")

#### 1. Data Formating #####

## Included only for reference. You can skip running this: the formated files have been saved.

# Generate AMC and AMR summary files from raw data files.
#source("code/preprocessing/0_2_european_data_cleaning.R")

# Minor formatting changes to variable names

source("code/data_cleaning.R")

#### 2. Data Cleaning ####

# eliminate outliers from AMC data
source("code/remove_amc_outliers.R")

##### 3. Assigning temporal trends #####
source("code/logistic_fitting_functions.R")
source("code/fit_temporal_trend.R")

# function that takes three arguments: patient type, min number of years for inclusion, min data points per year, min number of R+I isolates
logistic_fits_wrapper("INPAT", 5, 30)
logistic_fits_wrapper("OUTPAT", 5, 30)

##### 4. Analysis of temporal trends #####

# plotting what fits look like and goodness of fits
source("analysis/goodness_of_fit.R")

# summary stats for trajectories
source("analysis/trajectory_stats.R")

# Plotting number of pathogen-drug-country in each category, using predictors
source("analysis/Analysis_temporal_trends.R")

# Plot the coefficient of the binomial predictors as a map of Europe
source("analysis/Europe_map.R")

# Analysis of the rate of change in increasing vs stabilising trajectories
source("analysis/speed_of_increase.R")

##### 5. Correlate with AMC, analyse and visualise the correlations

source("code/correlate_amc_amr.R")

source("code/analyse_temporal_spatial_correlations.R")

#### 6. Trend in AMC vs trend in resistance

source("analysis/consumption_temporal_trend.R")

##### 7. Supporting Information #####

#add the session info to have all packages that were used for the analysis
sessionInfo()
