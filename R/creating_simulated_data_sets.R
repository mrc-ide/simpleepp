##############################################################################################################################
## Creating datasets for fitting to ##########################################################################################
##############################################################################################################################
require(ggplot2)
require(ggpubr)
require(reshape2)
require(rstan)
require(devtools)
require(splines)

find_rtools()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/cd4_spline_model.stan")

########################################################################################################################################
## Run the simulated model to get the output ###########################################################################################
########################################################################################################################################


run_simulated_model<-function(params,times){
  
  mu <- params$mu                               # Non HIV mortality / exit from population
  sigma <- params$sigma           # Progression from stages of infection
  mu_i <- params$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params$iota
  
  
  dt <- times$dt                                     # time step
  nsteps <- as.integer(times$years/dt)                   # number of steps
  xstart <- times$start                                # the start of the epidemic
  step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  rlogistic <- function(t, p) {
    p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
  }
  
  kappa_params<-c(log(params$kappa[1]),log(params$kappa[2]),params$kappa[3],params$kappa[4])
  kappa<-exp(rlogistic(step_vector,kappa_params))
  
  mod<-simpleepp_no_art(kappa = kappa,iota,mu,sigma,mu_i,dt)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa))
  
  
}

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_output<-run_simulated_model(params_sim,times_sim)

sim_plot<-function(sim_df){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=lambda),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  melt_df<-sim_df[,-c(4)]
  melted_df_funky<-melt(melt_df,id="time")
  
  three_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot,incidence=incidence_plot,kappa=kappa_plot,whole=three_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)

#######################################################################################################################################
## Now lets form our loop to sample through and fit to stan data ######################################################################
#######################################################################################################################################

sample_range<-2000:2015
## Need to change this to length of time_points_to_sample if sporadic
sample_n<-5000

exponential_decay_function<-function(N0,t,lambda){
  
  nt<-N0*exp(-t*lambda)
  
  return(nt)
}

decay<- 0  #exponential_decay_function(N0 = 90,t=(seq(0,50,0.1)),lambda = 0.05) # Must be put to equal 0 for no underreporting

penalty_order<-1

time_points_to_sample<- 0 #seq(1970,2015,3)                                   ## Must also be equal to 0 for complete reporting 


sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

sample_years<-length(rows_to_evaluate)

sample_df_tot<-NULL

total_timo<-NULL

iterations<-100

for(i in 1:iterations){
  
  timo<-Sys.time()
  
  sample_function<-function(year_range,number_of_years_to_sample,people_t0_sample,simulated_df,prevalence_column_id,decay=0,time_points_to_sample=0){
    sample_years_hiv <- number_of_years_to_sample # number of days sampled throughout the epidemic
    sample_n <- people_t0_sample # number of host individuals sampled per day
    
    # Choose which days the samples were taken. 
    # Ideally this would be daily, but we all know that is difficult.
    
    sample_time_hiv<-time_points_to_sample
    
    if(time_points_to_sample[1] == 0){
      sample_time_hiv = sort(sample(year_range, sample_years_hiv, replace=F))
    }
    # Extract the "true" fraction of the population that is infected on each of the sampled days:
    sample_propinf_hiv = simulated_df[simulated_df$time %in% sample_time_hiv, prevalence_column_id]
    
    ## this just samples our prevalence, to get a probability that the sample we take is HIV infected then we need to divide
    ## by 100
    
    #sample_propinf_hiv<-sample_propinf_hiv/100
    
    # Generate binomially distributed data.
    # So, on each day we sample a given number of people (sample_n), and measure how many are infected.
    # We expect binomially distributed error in this estimate, hence the random number generation.
    sample_y_hiv_prev = rbinom(length(sample_time_hiv), sample_n, sample_propinf_hiv)
    
    if(decay[1] != 0){
      proportion_reported= 1 - (decay/100)
      times_to_sample<-(sample_time_hiv - 1970) * 10 + 1
      proportion_reported = proportion_reported[times_to_sample]
      
      sample_y_hiv_prev = round(proportion_reported * sample_y_hiv_prev) 
    }
    
    
    
    sample_prev_hiv_percentage<-(sample_y_hiv_prev/sample_n)*100
    
    ## lets have a ggplot of the y (infected) and out sample of Y over time 
    sample_df_100<-data.frame(cbind(sample_time_hiv,sample_prev_hiv_percentage,sample_y_hiv_prev,sample_propinf_hiv))
    return(sample_df_100)  
  }
  
  
  
  sample_df_100_random_second<-sample_function(sample_range,sample_years,sample_n,
                                               simulated_df = sim_model_output$sim_df,prevalence_column_id = 3,decay = decay,
                                               time_points_to_sample = time_points_to_sample)
  
  sample_df_100_random_second$iteration<-rep(i,nrow(sample_df_100_random_second))
  
  sample_df_tot<-rbind.data.frame(sample_df_tot,sample_df_100_random_second)
  
  print(i/iterations * 100)
  
}

plot(sample_df_tot[sample_df_tot$iteration==64,2])
for(i in 1:100){
  colour="midnightblue"
  if( ((i+4)/5) == round((i+4)/5)){
    colour="red"
  }
  if(((i+3)/5) == round((i+3)/5)){
    colour="forestgreen"
  }
  if(((i+2)/5) == round((i+2)/5)){
    colour="yellow"
  }
  if(((i+1)/5) == round((i+1)/5)){
    colour="orange"
  }
  
  lines(sample_df_tot[sample_df_tot$iteration==i,2],col=colour)
  Sys.sleep(0.5)
  print(i)
}

  
sampled_n_5000_00_data<-sample_df_tot

save(sampled_n_5000_00_data,file = "C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_2000_runs/N_5000_sampled_data")  
save(sim_model_output,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/true_epidemic")
