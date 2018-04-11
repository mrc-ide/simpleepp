library(rstan)
library(splines)

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
xstart <- 1970
xx <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, xx)

#' Simulation of EPP model with constant kappa = 0.3
mod <- simpleepp(kappa=rep(0.3, 500), iota=0.005, alpha=rep(0, nsteps), mu, sigma, mu_i, mu_a, omega, dt)

#' Simulation of EPP classic model

mod_classic <- simpleepp_classic(r=1.5, f0=0.4, iota=0.0001, phi=0, mu, sigma, mu_i, dt, nsteps)

#' # True prevalence and ART coverage

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

alpha <- approx(names(artcov), artcov, xx)$y

#' # Simulate some data from 1985 to 2015

set.seed(12631569)
t_obs <- as.character(1980:2015)
n_obs <- rep(1500, length(t_obs))  
x_obs <- rbinom(length(t_obs), n_obs, prev[t_obs])
idx_obs <- match(t_obs, xout)


#' # Fit models

#' Negative log likelihood
#' @param theta Vector of parameter values.
#' @param simmod Function for simulating model given parameters
#' @param ... Additional arguments to \code{simmod}
nll <- function(theta, simmod, ...){
  mod <- simmod(theta, ...)
  -sum(dbinom(x_obs, n_obs, mod[idx_obs, 3], log=TRUE))
}
  
#' ## Fit the EPP classic model

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

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate


#' ## Fit the rlogistic model

rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

sim_rlogistic <- function(theta){
  kappa <- exp(rlogistic(xx, theta[1:4]))
  iota <- exp(theta[5])
  simpleepp(kappa, iota, alpha, mu, sigma, mu_i, mu_a, omega, dt)
}

theta0 <- c(log(0.5), log(0.1), 0.3, 1995, log(0.001))
sim_rlogistic(theta0)
nll(theta0, sim_rlogistic)

fit <- optim(theta0, nll, simmod=sim_rlogistic, method="BFGS")
fit$simmod <- sim_rlogistic
fit$mod <- fit$simmod(fit$par)

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate

fit_rlogistic <- fit


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
Xrw <- splineDesign(1969:2021, xx, ord=2)
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

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate


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

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate


#' ## Fit spline model
#'
#' ### First order penalty
#'
nk <- 7 # number of splines
dk <- diff(range(xout))/(nk-3)
knots <- xstart + -3:nk*dk

Xsp <- splineDesign(knots, xx, ord=4)
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

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate


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

plot(xout, fit$mod[,3], type="l", lwd=3, col="blue") # prevlance
points(t_obs, x_obs / n_obs, pch=19)

plot(xout, fit$mod[,2], type="l", lwd=3, col="blue") # incidence
plot(xout, fit$mod[,1], type="l", lwd=3, col="blue") # transmission rate


#' ## r-trend model
#'
#' ## r-spline with equilibrium prior
#' 
#' ## rlogistic-RW model
#'
#' ## rlogistic-AR(1) model
#' 

