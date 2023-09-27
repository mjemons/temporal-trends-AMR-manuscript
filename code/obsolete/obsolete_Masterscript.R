## This is a masterfile for the temporal trends project. 

library(minpack.lm)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(svglite)
library(nnet)
library(lmtest)
# library(forestplot)
library(gridExtra)
library(DescTools)

#setwd("/Users/francois.blanquart/ownCloud/AMR_Sonja_Martin/temporal_trends_AMR")


#### 1. Data Formating #####

## Included only for reference. You can skip running this: the firmated files have been saved.

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
# note: still need to work out best way to convey this info.
source("analysis/goodness_of_fit.R")

# Plotting temporal trends

source("analysis/plot_logistic_slopes.R")

# Plotting categories, analysing predictors
# Note: this is work in progress.
source("analysis/Analysis_temporal_trends.R")

##### 5. Correlate with AMC

source("code/correlate_amc.R")

source("analysis/analysis_correlations.R")

##### 7. Supporting Information #####
