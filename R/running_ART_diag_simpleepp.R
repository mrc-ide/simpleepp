###################################################################################################################################
## Testing out how the new art diseased model work ################################################################################
###################################################################################################################################

require(ggplot2)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=T)

expose_stan_functions("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model.STAN")

rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

step_vector<-seq(1970.1,by = 0.1,length.out = 500)
kappa_params<-c(log(0.5),log(0.1),0.3,1995)

kappa<- exp(rlogistic(step_vector,kappa_params))                       ## This is our transmission parameter for the output

iota_val<-0.0001                                                       ## This is our initial proportion infected

alpha<- c(0.7,0.75,0.8,0.85)                                          ## This is our diagnosed to ART treatment rate

mu <- 1/35                                                             ## Population level death rate
            
sigma<- c(1/3.16,1/2.13,1/3.2)                                         ## Progression along Cd4 stages among untreated

mu_i <- c(0.003, 0.008, 0.035, 0.27)                                    # Mortality by stage, no ART

mu_d <- c(0.003,0.008,0.035,0.27)                                      ## Mortality be stage on diagnosed

mu_a <- c(0.002, 0.006, 0.006, 0.03)                                   ## Mortality by stage, on ART

omega <- 0.7                                                           ## Reduction in trnasmissability on art

theta <- 0.2                                                          ## Reduction when know diagnosed

dt <- 0.1

start<- 1970

diag_start<- 1981

art_start<-2019

diag<-c(0.3,0.4,0.7,0.85)                                             ## Proportion to be diagnosed at each infected cd4 stage

art_prog<-c(1/10,1/10,1/10)

test_run<-simpleepp_art_diag(kappa = kappa,iota = iota,alpha = alpha,mu = mu,sigma = sigma, mu_i = mu_i, mu_d = mu_d,mu_a = mu_a,
                             omega = omega,theta = theta,dt = dt,start = start,diag_start = diag_start,art_start = art_start,
                             diag = diag,art_prog = art_prog)
plot(test_run[,2])
