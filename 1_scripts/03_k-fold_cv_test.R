# R code to run k-fold cross validation for modeling 

# Load required package
library(dplyr)
library(tidyr)
library(data.table)
library(jagsUI)


load("0_data/processed/osfl_data_package.Rdata")


### set MCMC specifications ###

nc = 3  # number of chains (5)
nb = 500  # number of adaptation iterations (1000)
nt = 1
ni = 1000 # chain length (1000)

#Check that the order of BS factor values is FY, P, S
print(levels(as.factor(daily_sc$standardized_bs)))

traindat <- 



