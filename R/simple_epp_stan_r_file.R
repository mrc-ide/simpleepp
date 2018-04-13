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
require(gridExtra)
require(grid)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################################################################################################################
## Now lets set our initial conditions for the model ###########################################################################
################################################################################################################################

cd4_1<-100                                  ## Initial number of the population infected with 500>cd4>350
cd4_2<-0                                   ## Initial number of the population infected with 350>cd4>200
cd4_3<-0                                   ## Initial number of the population infected with 200>cd4>50
cd4_4<-0                                   ## Initial number of the population infected with cd4<50
s0<-1000000-(cd4_1+cd4_2+cd4_3+cd4_4)      ## Initial susceptible proportion of the population
n0<-s0+cd4_1+cd4_2+cd4_3+cd4_4             ## Initial population size
i0<-(cd4_1+cd4_2+cd4_3+cd4_4)              ## Initial infected size in total
prev0<-i0/n0                               ## Initial Prevalence
incide0<-i0/n0                             ## Initial Incidence 

inits_hiv_cd4<-c(s0,cd4_1,cd4_2,cd4_3,cd4_4,i0,n0,prev0,incide0)

################################################################################################################################
## Now we will set up the intial paramter values for the model #################################################################
################################################################################################################################

eta<- 20000                                ## The birth rate into the population
kappa<- 0.3                                ## The transmission rate paramter
sigma<-c(1/3.16,1/2.13,1/3.2)              ## Progression rates from each Cd4 class
mu_i<-c(0.003,0.008,0.035,0.27)            ## Vector of death rates in each cd4 class         
mu_s<-1/35                                 ## Normal population death rate 

params_cd4_hiv<-list(et=et,kappa=kappa,sigma=sigma,mu_i=mu_i,mu_s=mu_s)

################################################################################################################################
## Now we will create the time series over which to integrate the functions ####################################################
################################################################################################################################

t_min<-0
t_max<-100
times<-t_min:t_max

################################################################################################################################
## Now lets create the system of ODEs to model the change in the population over time ##########################################
################################################################################################################################

epp_cd4<-function(t, y, parms,...) {
  with(as.list(c(params_cd4_hiv, y)), {
    
    dS = eta - ((kappa * y[6] * y[1]) / y[7]) - (mu_s * y[1])
    
    dcd4_1 = ((kappa * y[6] * y[1]) / y[7]) - (mu_i[1] * y[2]) - (sigma[1] * y[2])
    
    dcd4_2 = (sigma[1]*y[2]) - (mu_i[2] * y[3]) - (sigma[2] * y[3])
    
    dcd4_3 = (sigma[2] * y[3]) - (mu_i[3] * y[4]) - (sigma[3] * y[4])
    
    dcd4_4 = (sigma[3] * y[4]) - (mu_i[4] * y[5])
    
    dInfec = dcd4_1 + dcd4_2 +dcd4_3 + dcd4_4
    
    dN = dS + dInfec
    
    dprev = ((dInfec + y[6]) / (dN + y[7])) - y[8]
    
    dincidence = -y[9] + (((kappa * y[6] * y[1]) / y[7])/y[7])
    
    res <- c(dS,dcd4_1,dcd4_2,dcd4_3,dcd4_4,dInfec,dN,dprev,dincidence)
    list(res)
    
    
    })

}


################################################################################################################################
## Now we have written the function we will try to use the ode solver in deSolve to iterate through ############################
################################################################################################################################

out_epp_cd4_hiv <- ode(inits_hiv_cd4, times, epp_cd4, params_cd4_hiv, method="ode45")

out_epp_cd4_hiv<-data.frame(out_epp_cd4_hiv)

names(out_epp_cd4_hiv)<-c("time","susceptible","cd4 > 500","500 > cd4 > 350",
                          "350 > cd4 > 200", "200 > cd4", "infected", "n", "prevalence","incidence")
out_epp_cd4_hiv$prev_percent<-out_epp_cd4_hiv$prevalence * 100

out_epp_cd4_hiv

################################################################################################################################
## Now we have some data we can plot our prevalence curve and incidence curve ##################################################
################################################################################################################################

prev_plot<-ggplot(data = out_epp_cd4_hiv, aes(x=time,y=prev_percent))+geom_line(colour="blue")+
  labs(x="Time",y="Prevalence (%)")

total_infected_plot<-ggplot(data = out_epp_cd4_hiv,aes(x=time,y=incidence))+geom_line(colour="red")+
  labs(x="Time",y="Incidence rate")

total_pop_plot<-ggplot(data = out_epp_cd4_hiv,aes(x=time,y=n))+geom_line(colour="forest green")+
  labs(x="Time",y="Total Population size")

grid::grid.draw(rbind(ggplotGrob(prev_plot), ggplotGrob(total_infected_plot), ggplotGrob(total_pop_plot),  size = "last"))

melt_incidence_prev<-out_epp_cd4_hiv[,-c(2:8,11)]
melted_incidence_prevalence<-melt(melt_incidence_prev,id="time")

incidence_prevalence_plot<-ggplot(data = melted_incidence_prevalence)+geom_line(aes(x=time,y=value,colour=variable),size=1.1)+
  labs(x="Time")
plot(incidence_prevalence_plot)

################################################################################################################################
## Now we will draw samples from our simulated epidemic to then fit our stan model to ##########################################
################################################################################################################################

sample_function<-function(number_of_years_to_sample,people_t0_sample,simulated_df,prevalence_column_id,t_max){
  sample_years_hiv <- number_of_years_to_sample # number of days sampled throughout the epidemic
  sample_n <- people_t0_sample # number of host individuals sampled per day
  
  # Choose which days the samples were taken. 
  # Ideally this would be daily, but we all know that is difficult.
  sample_time_hiv = sort(sample(1:t_max, sample_years_hiv, replace=F))
  
  # Extract the "true" fraction of the population that is infected on each of the sampled days:
  sample_propinf_hiv = simulated_df[simulated_df$time %in% sample_time_hiv, prevalence_column_id]
  
  ## this just samples our prevalence, to get a probability that the sample we take is HIV infected then we need to divide
  ## by 100
  
  #sample_propinf_hiv<-sample_propinf_hiv/100
  
  # Generate binomially distributed data.
  # So, on each day we sample a given number of people (sample_n), and measure how many are infected.
  # We expect binomially distributed error in this estimate, hence the random number generation.
  sample_y_hiv_prev = rbinom(sample_years_hiv, sample_n, sample_propinf_hiv)
  sample_prev_hiv<-(sample_y_hiv_prev/sample_n)*100
  
  ## lets have a ggplot of the y (infected) and out sample of Y over time 
  sample_df_100<-data.frame(cbind(sample_time_hiv,sample_prev_hiv))
  return(sample_df_100)  
}
sample_df_100<-sample_function(100,25,simulated_df = out_epp_cd4_hiv,prevalence_column_id = 9,t_max = 100)

plot_sample<-function(sample_df,simulated_df){
a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
  geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv),colour="red",size=1)

return(plot(a))
}

plot_sample(simulated_df = out_epp_cd4_hiv,sample_df = sample_df_100)


ggplot(data = sample_df_100,aes(x=sample_time_hiv,y=sample_prev_hiv))+geom_point(colour="red",size=1.5)

################################################################################################################################
## Now we have our sample data we can the data to a model in STAN ##############################################################
################################################################################################################################

stan_d_hiv_prev = list(n_obs = sample_years_hiv,
                       n_params = 2,
                       n_difeq = length(inits_hiv_cd4),
                       n_sample = sample_n,
                       n_fake = length(1:t_max),
                       y = sample_y_hiv_prev,
                       t0 = 0,
                       ts = sample_time_hiv,
                       fake_ts = c(1:t_max))

# Which parameters to monitor in the model:
params_monitor_hiv = c("y_hat", "y0", "params", "fake_I")

# Test / debug the model:
test_hiv_100_year = stan("hiv_project/simpleepp/stan_files/chunks/cd4_stan.stan",
                                                             data = stan_d_hiv_prev,
                                                             pars = params_monitor_hiv,
                                                             chains = 1, iter = 10)

## Run the model 

mod_hiv_prev = stan(fit = test_hiv_100_year,data = stan_d_hiv_prev,pars = params_monitor,chains = 3,warmup = 500,iter = 1500,
                                               control = list(adapt_delta = 0.85))
# Extract the posterior samples to a structured list:
posts_hiv <- extract(mod_hiv_prev)

apply(posts_hiv$params, 2, median)


apply(posts_hiv$y0, 2, median)[1:9]

# These should match well. 

#################
# Plot model fit:

# Proportion infected from the synthetic data:

#sample_prop = sample_y / sample_n

# Model predictions across the sampling time period.
# These were generated with the "fake" data and time series.
mod_median = apply(posts_hiv$fake_I[,,2], 2, median)
mod_low = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.025))
mod_high = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.975))
mod_time = stan_d_hiv$fake_ts

prev_median<-(apply(posts_hiv$fake_I[,,8],2,median))*100
prev_low<-(apply(posts_hiv$fake_I[,,8],2,quantile,probs=c(0.025)))*100
prev_high<-(apply(posts_hiv$fake_I[,,8],2,quantile,probs=c(0.975)))*100

# Combine into two data frames for plotting
#df_sample = data.frame(sample_prop, sample_time)
df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, stan_d_hiv_prev$ts)
names(df_fit_prevalence)<-c("median","low","high","time")

# Plot the synthetic data with the model predictions
# Median and 95% Credible Interval

n_25_plot<-ggplot(sample_df_100, aes(x=sample_df_100$sample_time_hiv, y=sample_df_100$sample_prev_hiv)) +
  geom_point(col="red", shape = 19, size = 1.5) +
  geom_line(data = df_fit_prevalence, aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),
              colour="midnightblue",alpha=0.2,fill="midnightblue")+
  coord_cartesian(ylim = c(0,100),xlim=c(0,100))+labs(x="Time",y="Prevalence (%)", title="N = 25 plot")
plot(n_25_plot)

grid.arrange(rbind(ggplotGrob(n_25_plot),ggplotGrob(n_100_plot),ggplotGrob(n_200_plot), size="last"))





