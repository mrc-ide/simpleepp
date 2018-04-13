library(rstan)
library(splines)
library(ggplot2)
library(reshape2)
library(ggpubr)

expose_stan_functions("stan_files/simpleepp_expose.stan")


####

#' Fixed model parameters 
mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
mu_a <- c(0.002, 0.006, 0.006, 0.03)     # Mortality by stage, on ART
omega <- 0.7                             # Average effect of ART on reducing transmission

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total

#' Simulation of EPP model with constant kappa = 0.3
#' This assumes constant population size and intiailly a constant r that can be changed. 
#' This returns a matrix with r in first column, incidence in the second column and prevalence in the third column
#' Iota refers to the initial proportion of the population infected. 
mod <- simpleepp(kappa=rep(0.3, length(xout)), iota=0.005, alpha=rep(0, nsteps), mu, sigma, mu_i, mu_a, omega, dt)

#' Simulation of EPP classic model

mod_classic <- simpleepp_classic(r=1.5, f0=0.4, iota=0.0001, phi=0, mu, sigma, mu_i, dt, nsteps)

############################################################################################################################
## Now we are inputting prevalence and ART coverage data from which to create simulated data to fit the two models to ######
############################################################################################################################
#' # True prevalence and ART coverage
#' 

prev <- structure(c(0, 0, 1e-04, 2e-04, 4e-04, 7e-04, 0.0012, 0.002, 
                    0.0033, 0.0052, 0.0079, 0.0116, 0.0163, 0.0222, 0.0294, 0.0374, 
                    0.0464, 0.0564, 0.0672, 0.079, 0.0912, 0.1032, 0.1146, 0.1249, 
                    0.1339, 0.1409, 0.1464, 0.1504, 0.1526, 0.1533, 0.1522, 0.1495, 
                    0.1454, 0.1403, 0.1345, 0.1288, 0.1234, 0.1188, 0.1149, 0.1115, 
                    0.1085, 0.1056, 0.1026, 0.0998, 0.097, 0.0942, 0.0913, 0.0884, 
                    0.0856, 0.0827, 0.0799),
                  .Names = 1970:2020)

artcov <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0126, 0.0344, 
                      0.0713, 0.119, 0.1749, 0.2284, 0.283, 0.3595, 0.4411, 0.4974, 
                      0.5639, 0.6071, 0.687, 0.6818, 0.7223, 0.7621, 0.8009),
                    .Names=1970:2020)

alpha <- approx(names(artcov), artcov, step_vector)$y

#' # Simulate some data from 1985 to 2015

set.seed(12631569)
t_obs <- as.character(1980:2015)                    ## This creates time of samples, realistically this would occur once a year
n_obs <- rep(1500, length(t_obs))                   ## This samples 1500 people from the population every year         
x_obs <- rbinom(length(t_obs), n_obs, prev[t_obs])  ## This creates our binomally sampled data given our prevalence data
idx_obs <- match(t_obs, xout)                       ## This links back to our intial time series, giving us the row number of the model data that corresponds with the simulated data



#' # Fit models

#' Negative log likelihood
#' nll is the function that we input initial paramter values into and for which it then runs the model and calcualtes the
#' likelihood of those paramaters given our predicted data. The transmission paramter is allowed to vary with time. 
#' @param theta Vector of parameter values.
#' @param simmod Function for simulating model given parameters
#' @param ... Additional arguments to \code{simmod}
nll <- function(theta, simmod, ...){
  mod <- simmod(theta, ...)
  -sum(dbinom(x_obs, n_obs, mod[idx_obs, 3], log=TRUE))
}
  
#' ## Fit the EPP classic model
#' in this case we are creating the simmod function with which to call the above nll function


sim_classic <- function(theta){
  r <- exp(theta[1])
  f0 <- plogis(theta[2])
  iota <- exp(theta[3])
  phi <- theta[4]
  simpleepp_classic(r, f0, iota, phi, mu, sigma, mu_i, dt, nsteps)
}

## initial parameters
theta0 <- c(r=log(1.5), f0=qlogis(0.4), iota=log(0.001), phi=0)
nll(theta0, sim_classic)

## fit model via maximum likelihood
fit <- optim(theta0, nll, simmod=sim_classic, method="BFGS")
fit$simmod <- sim_classic
fit$mod <- fit$simmod(fit$par)

## Creating a df of the maximum likelihood model output for easier plotting 

plot_df<-fit$mod
plot_df<-data.frame(plot_df)
names(plot_df)<-c("Transmission_rate","Incidence","Prevalence")
plot_df$prev_percent<-plot_df$Prevalence * 100
plot_df$Time<-xout

## creating df of the sampled data 
sample_df<-data.frame(cbind(as.numeric(t_obs),as.numeric(x_obs)))
names(sample_df)<-c("Time","Infected")
sample_df$Prevalence<-(sample_df$Infected / n_obs) * 100

prevalence_plot<- ggplot(data = plot_df) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "ML fitting of Classic EPP to simulated data")

plot(prevalence_plot)

### Plotting all the values output on the same graph, in this case prevalence is rescaled to a fraction

melt_df<-plot_df[,-c(4)]
melted_mod<-melt(melt_df,id="Time")

incidence_and_rho_plot<-ggplot(data = melted_mod) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="ML fitting of Classic EPP")

plot(incidence_and_rho_plot)

############################################################################################################################
## So we have updated the plotting of the ML fitting methods, now we'll implement the r logisitic curve model for the ######
## transmission rate paramter ##############################################################################################
############################################################################################################################


#' ## Fit the rlogistic model
#' ## The r_logistic function creates the logistic curve for the infection transmission parameter, has more of a biological 
#' basis and can be more easily interpreted. 



rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

sim_rlogistic <- function(theta){
  kappa <- exp(rlogistic(step_vector, theta[1:4]))
  iota <- exp(theta[5])
  simpleepp(kappa, iota, alpha, mu, sigma, mu_i, mu_a, omega, dt)
}

theta0 <- c(log(0.5), log(0.1), 0.3, 1995, log(0.001))   ## here the first four values are used in forming the logistic curve for the r and the last character is the new start infection number
sim_rlogistic(theta0)
nll(theta0, sim_rlogistic)

fit <- optim(theta0, nll, simmod=sim_rlogistic, method="BFGS")
fit$simmod <- sim_rlogistic
fit$mod <- fit$simmod(fit$par)

############################################################################################################################
## Now we have fitted the data using this method we can plot the results below #############################################
############################################################################################################################

log_plot_df<-data.frame(fit$mod)
names(log_plot_df)<-c("Transmission rate", "Incidence", "Prevalence")
log_plot_df$prev_percent<-log_plot_df$Prevalence * 100
log_plot_df$Time<-xout

prevalence_plot_r_log<- ggplot(data = log_plot_df) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "R Logistic fitting of Classic EPP to simulated data")

plot(prevalence_plot_r_log)

melt_r_log_df<- log_plot_df[,-c(4)]
melted_mod_log<-melt(melt_r_log_df,id="Time")

incidence_and_rho_plot_r_log<-ggplot(data = melted_mod_log) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="Logistic curve fitting of R fitting of Classic EPP")

plot(incidence_and_rho_plot_r_log)


fit_rlogistic <- fit

############################################################################################################################
## Now moving on to fit the random walk model for the R paramter ###########################################################
############################################################################################################################

#' ## Fit random walk modeL

#' Generic function for simulating semi-parametric $\kappa(t)$
#'
#' This can be used for any specification of $\log \kappa(t)$
#' which can be summarized as a linear combination of basis
#' functions $\log \kappa(t) = X \beta$ and generic penaltization
#' expressed as $D \beta ~ Normal(0, \sigma)$
#' 
#' @param theta Vector of parameters
#' @param X A design matrix 
sim_semipar <- function(theta, X){
  beta <- theta[1:ncol(X)]
  kappa <- c(exp(X %*% beta))
  iota <- exp(theta[ncol(X)+2])
  simpleepp(kappa, iota, alpha, mu, sigma, mu_i, mu_a, omega, dt)
}

#' Negative log posterior for penalized model
nlp_semipar <- function(theta, X, D){
  beta <- theta[1:ncol(X)]
  sigma_pen <- exp(theta[ncol(X)+1])

  log_penalty <- sum(dnorm(D %*% beta, 0, sigma_pen, log=TRUE))
  nll(theta, sim_semipar, X) - log_penalty
}

#' ### First-order random walk (peicewise-linear)
Xrw <- splineDesign(1969:2021, step_vector, ord=2)
Drw1 <- diff(diag(ncol(Xrw)), diff=1)

beta0 <- log(fit$mod[0:50*10+1, 1]) # use log kappa values from previous fit
theta0 <- c(beta = beta0, sigma_pen = log(0.1), iota = log(0.004))
sim_semipar(theta0, Xrw)
nll(theta0, sim_semipar, Xrw)
nlp_semipar(theta0, Xrw, Drw1)
            
fit <- optim(theta0, nlp_semipar, X=Xrw, D=Drw1, method="BFGS")
fit$simmod <- sim_semipar
fit$X <- Xrw
fit$D <- Drw1
fit$mod <- fit$simmod(fit$par, Xrw)

log_plot_rw_first<-data.frame(fit$mod)
names(log_plot_rw_first)<-c("Transmission rate", "Incidence", "Prevalence")
log_plot_rw_first$prev_percent<-log_plot_rw_first$Prevalence * 100
log_plot_rw_first$Time<-xout


prevalence_plot_r_rw_first<- ggplot(data = log_plot_rw_first) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "R RW first order fitting of Classic EPP to simulated data")

plot(prevalence_plot_r_rw_first)

melt_r_rw_first_df<- log_plot_rw_first[,-c(4)]
melted_mod_rw_first<-melt(melt_r_rw_first_df,id="Time")

incidence_and_rho_plot_r_RW_plot<-ggplot(data = melted_mod_rw_first) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="RW First order fitting of R fitting of Classic EPP")

plot(incidence_and_rho_plot_r_RW_plot)

############################################################################################################################
## Now moving on to fitting a RW based on second order differences, which will maintain previous curve value in absence of #
## Data, so when predicting future incidence ###############################################################################
############################################################################################################################

#' ### Second order random walk

Drw2 <- diff(diag(ncol(Xrw)), diff=2)

beta0 <- log(fit$mod[0:50*10+1, 1]) # use log kappa values from previous fit
theta0 <- c(beta = beta0, sigma_pen = log(0.1), iota = log(0.004))
sim_semipar(theta0, Xrw)
nll(theta0, sim_semipar, Xrw)
nlp_semipar(theta0, Xrw, Drw2)
            
fit <- optim(theta0, nlp_semipar, X=Xrw, D=Drw2, method="BFGS")
fit$simmod <- sim_semipar
fit$X <- Xrw
fit$D <- Drw2
fit$mod <- fit$simmod(fit$par, Xrw)

log_plot_rw_second<-data.frame(fit$mod)
names(log_plot_rw_second)<-c("Transmission rate", "Incidence", "Prevalence")
log_plot_rw_second$prev_percent<-log_plot_rw_second$Prevalence * 100
log_plot_rw_second$Time<-xout


prevalence_plot_r_rw_second<- ggplot(data = log_plot_rw_second) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "R RW 2nd order fitting of Classic EPP to simulated data")

plot(prevalence_plot_r_rw_second)

melt_r_rw_second_df<- log_plot_rw_second[,-c(4)]
melted_mod_rw_second<-melt(melt_r_rw_second_df,id="Time")

incidence_and_rho_plot_r_RW_second_plot<-ggplot(data = melted_mod_rw_second) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="RW Second order fitting of R fitting of Classic EPP")

plot(incidence_and_rho_plot_r_RW_second_plot)

############################################################################################################################
## Now we will model R with a spline function, firstly a first order penalized spline ######################################
############################################################################################################################


#' ## Fit spline model
#'
#' ### First order penalty
#'
nk <- 7 # number of splines
dk <- diff(range(xout))/(nk-3)
knots <- xstart + -3:nk*dk

Xsp <- splineDesign(knots, step_vector, ord=4)
Dsp1 <- diff(diag(nk), diff=1)

beta0 <- lm(log(fit_rlogistic$mod[-1,1]) ~ -1+Xsp)$coef
theta0 <- c(beta = beta0, sigma_pen = log(0.1), iota = log(0.004))
sim_semipar(theta0, Xsp)
nll(theta0, sim_semipar, Xsp)
nlp_semipar(theta0, Xsp, Dsp1)
            
fit <- optim(theta0, nlp_semipar, X=Xsp, D=Dsp1, method="BFGS")
fit$simmod <- sim_semipar
fit$X <- Xsp
fit$D <- Dsp1
fit$mod <- fit$simmod(fit$par, Xsp)

log_plot_spline_first_spline<-data.frame(fit$mod)
names(log_plot_spline_first_spline)<-c("Transmission rate", "Incidence", "Prevalence")
log_plot_spline_first_spline$prev_percent<-log_plot_spline_first_spline$Prevalence * 100
log_plot_spline_first_spline$Time<-xout


prevalence_plot_spline_first_10<- ggplot(data = log_plot_spline_first_spline) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "R First Order Spline . Classic EPP to simulated data")

plot(prevalence_plot_spline_first_10)

melt_r_spline_first_10_df<- log_plot_spline_first_spline[,-c(4)]
melted_mod_spline_first_10<-melt(melt_r_spline_first_10_df,id="Time")

incidence_and_rho_plot_r_spline_first_10<-ggplot(data = melted_mod_spline_first_10) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="R spline 10 First Order fitting of R fitting of Classic EPP")

plot(incidence_and_rho_plot_r_spline_first_10)

############################################################################################################################
## Now we will model with a second order penalized spline ##################################################################
############################################################################################################################


#' ### Second order penalty
#'
Dsp2 <- diff(diag(nk), diff=2)

beta0 <- lm(log(fit_rlogistic$mod[-1,1]) ~ -1+Xsp)$coef
theta0 <- c(beta = beta0, sigma_pen = log(0.1), iota = log(0.004))
sim_semipar(theta0, Xsp)
nll(theta0, sim_semipar, Xsp)
nlp_semipar(theta0, Xsp, Dsp2)
            
fit <- optim(theta0, nlp_semipar, X=Xsp, D=Dsp2, method="BFGS")
fit$simmod <- sim_semipar
fit$X <- Xsp
fit$D <- Dsp2
fit$mod <- fit$simmod(fit$par, Xsp)

log_plot_spline_second<-data.frame(fit$mod)
names(log_plot_spline_second)<-c("Transmission rate", "Incidence", "Prevalence")
log_plot_spline_second$prev_percent<-log_plot_spline_second$Prevalence * 100
log_plot_spline_second$Time<-xout


prevalence_plot_spline_second<- ggplot(data = log_plot_spline_second) + geom_line(aes(x=Time, y=prev_percent),colour="midnightblue",size=1.2)+
  geom_point(data = sample_df,aes(x=Time,y=Prevalence),colour="red",fill="red")+
  labs(x="Year",y="Prevalence (%)",title= "R Second Order Spline fitting of Classic EPP to simulated data")

plot(prevalence_plot_spline_second)

melt_r_spline_second_df<- log_plot_spline_second[,-c(4)]
melted_mod_spline_second<-melt(melt_r_spline_second_df,id="Time")

incidence_and_rho_plot_r_spline_second<-ggplot(data = melted_mod_spline_second) + geom_line(aes(x=Time, y= value, colour=variable),size=1.2)+
  labs(x="Year",title="R spline Second Order fitting of R fitting of Classic EPP")

plot(incidence_and_rho_plot_r_spline_second)

############################################################################################################################
## Now we will plot all these together in the same output for ease of comparison, we'll start by plotting the incidence, ###
## prevalence and rho plots for each of the different methods together into one plot #######################################
############################################################################################################################

model_output_plot<-ggarrange(incidence_and_rho_plot,incidence_and_rho_plot_r_log,incidence_and_rho_plot_r_RW_plot,
                             incidence_and_rho_plot_r_RW_second_plot,incidence_and_rho_plot_r_spline_first_10,
                             incidence_and_rho_plot_r_spline_second,ncol = 3,nrow = 2)

plot(model_output_plot)

prevalence_output_plots<-ggarrange(prevalence_plot,prevalence_plot_r_log,prevalence_plot_r_rw_first,
                                   prevalence_plot_r_rw_second,prevalence_plot_spline_first_10,prevalence_plot_spline_second)

plot(prevalence_output_plots)

## Now we will plot the transmission curves produced via the different methods together 

ML_time_varying<-plot_df$Transmission_rate
ml_varying<-cbind.data.frame(ML_time_varying,rep("ml_varying",nrow(plot_df)))
names(ml_varying)<-c("rate","method")

logistic_r_transmiss<-cbind.data.frame(log_plot_df$`Transmission rate`,rep("logistic_r_model",nrow(log_plot_df)))
names(logistic_r_transmiss)<-c("rate","method")

rw_first<-cbind.data.frame(log_plot_rw_first$`Transmission rate`,rep("RW first order",nrow(log_plot_rw_first)))
names(rw_first)<-c("rate","method")

rw_second<-cbind.data.frame(log_plot_rw_second$`Transmission rate`,rep("RW second order",nrow(log_plot_rw_second)))
names(rw_second)<-c("rate","method")

spline_first<-cbind.data.frame(log_plot_spline_first$`Transmission rate`,
                               rep("Spline First Order",nrow(log_plot_spline_first)))
names(spline_first)<-c("rate","method")

spline_second<-cbind.data.frame(log_plot_spline_second$`Transmission rate`,
                                rep("Spline Second Order",nrow(log_plot_spline_second)))
names(spline_second)<-c("rate","method")

transmission_rate_data<-rbind.data.frame(ml_varying,logistic_r_transmiss,rw_first,rw_second,spline_first,spline_second)
time_for_transmission<-rep(xout,6)
transmission_rate_data$time<-time_for_transmission

transmission_rate_comparison_plot<-ggplot(transmission_rate_data)+geom_line(aes(x=time,y=rate,colour=method),size=1.05)+
  labs(x="Year",y="Rate",title="Comparison of methods for estimating transmission rate")
plot(transmission_rate_comparison_plot)

################################################################################################################################
## NOw we will plot the incidence curves produced via each method ##############################################################
################################################################################################################################

ml_incidence<-cbind.data.frame(plot_df$Incidence,rep("ML_varying",nrow(plot_df)))
names(ml_incidence)<-c("incidence","method")

logistic_incidence<-cbind.data.frame(log_plot_df$Incidence,rep("R Logistic",nrow(log_plot_df)))
names(logistic_incidence)<-c("incidence","method")

rw_first_incidence<-cbind.data.frame(log_plot_rw_first$Incidence,rep("RW first order",nrow(log_plot_rw_first)))
names(rw_first_incidence)<-c("incidence","method")

rw_second_inference<-cbind.data.frame(log_plot_rw_second$Incidence,rep("RW second order",nrow(log_plot_rw_second)))
names(rw_second_inference)<-c("incidence","method")

spline_first_incidence<-cbind.data.frame(log_plot_spline_first_spline$Incidence,
                                         rep("Spline first order",nrow(log_plot_spline_first_spline)))
names(spline_first_incidence)<-c("incidence","method")

spline_second_incidence<-cbind.data.frame(log_plot_spline_second$Incidence,
                                          rep("spline second order",nrow(log_plot_spline_second)))
names(spline_second_incidence)<-c("incidence","method")

incidence_rate_comparison_df<-rbind.data.frame(ml_incidence,logistic_incidence,rw_first_incidence,rw_second_inference,
                                               spline_first_incidence,spline_second_incidence)
incidence_time<-rep(xout,6)
incidence_rate_comparison_df$time<-incidence_time

incidence_comparion_plot<-ggplot(data = incidence_rate_comparison_df)+geom_line(aes(x=time,y=incidence,colour=method),size=1.05)+
  labs(x="Year",y="incidence rate")
plot(incidence_comparion_plot)

#' ## r-trend model
#'
#' ## r-spline with equilibrium prior
#' 
#' ## rlogistic-RW model
#'
#' ## rlogistic-AR(1) model
#' 

