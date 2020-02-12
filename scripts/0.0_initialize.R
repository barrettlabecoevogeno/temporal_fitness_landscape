##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: 
# Requires 
# NOTES: 
##########################################################################################

# Load libraries ----------------------------------------------------------

library(data.table)
library(dplyr)
library(lubridate)
library(MASS)
library(mgcv)
library(tidyr)
library(vegan)


# Source scripts ----------------------------------------------------------
source('scripts/0.1_misc.R')
# Load functions 
source('./scripts/PCA_custom.R')
source("./scripts/logit.r")
# Prepare data script 
source('scripts/prep.data.R')