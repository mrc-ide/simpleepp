#' This script will implement a simple EPP model, with CD4 structure and a varying r paramter through time 

################################################################################################################################
## Implement simple EPP structure in R to simulate data, then fit STAN Model to this simulated data ############################
################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)
require(deSolve)
require(dplyr)
require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################################################################################################################
## Now lets set our initial conditions for the model ###########################################################################
################################################################################################################################

