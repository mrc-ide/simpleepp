###########################################################################################################################################
## Script for looping through multiple fittings of the random walk model ##################################################################
###########################################################################################################################################

require(ggplot2)
require(ggpubr)
require(reshape2)
require(rstan)
require(devtools)
require(splines)

find_rtools()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

expose_stan_functions("stan_files/chunks/cd4_matrix_random_walk.stan")

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

sample_range<-1970:2015
                                                         ## Need to change this to length of time_points_to_sample if sporadic
sample_n<-100

exponential_decay_function<-function(N0,t,lambda){
  
  nt<-N0*exp(-t*lambda)
  
  return(nt)
}

decay<- 0  #exponential_decay_function(N0 = 90,t=(seq(0,50,0.1)),lambda = 0.05) # Must be put to equal 0 for no underreporting

penalty_order<-1

time_points_to_sample<- 0 #seq(1970,2015,3)                                   ## Must also be equal to 0 for complete reporting 

rows_to_evaluate<- 0:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

sample_years<-length(rows_to_evaluate)

prev_df_tot<-NULL

incidence_df_tot<-NULL

kappa_df_tot<-NULL

iota_values_tot<-NULL

sample_df_tot<-NULL

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



xout<-seq(1970,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 




stan_data_discrete<-list(
  n_obs = sample_years,
  n_sample = sample_n,
  y = as.array(sample_df_100_random_second$sample_y_hiv_prev),
  time_steps_euler = 501,
  penalty_order = penalty_order,
  estimate_period = 5,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  


mod_hiv_prev <- stan("stan_files/chunks/cd4_matrix_random_walk.stan", data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.85))



plot_stan_model_fit<-function(model_output,sim_sample,sim_output,plot_name,xout){
  
  posts_hiv <- rstan::extract(model_output)
  
  
  iota_dist<-posts_hiv$iota
  params<-median(posts_hiv$iota)
  params_low<-quantile(posts_hiv$iota,c(0.025))
  params_high<-quantile(posts_hiv$iota,c(0.975))
  
  params_df<-rbind.data.frame(params_low,params,params_high)
  names(params_df)<-c("iota")
  
  sigma_pen_dist<-posts_hiv$sigma_pen
  sigma_values<-median(posts_hiv$sigma_pen)
  sigma_low<-quantile(posts_hiv$sigma_pen,c(0.025))
  sigma_high<-quantile(posts_hiv$sigma_pen,probs=c(0.975))
  sigma_df<-rbind.data.frame(sigma_low,sigma_values,sigma_high)
  names(sigma_df)<-c("sigma_pen")
  
  
  
  # These should match well. 
  
  #################
  # Plot model fit:
  
  # Proportion infected from the synthetic data:
  
  #sample_prop = sample_y / sample_n
  
  # Model predictions across the sampling time period.
  # These were generated with the "fake" data and time series.
  #mod_median = apply(posts_hiv$fake_I[,,2], 2, median)
  #mod_low = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.025))
  #mod_high = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.975))
  mod_time = xout
  
  prev_median<-(apply(posts_hiv$fitted_output[,,3],2,median))*100
  prev_low<-(apply(posts_hiv$fitted_output[,,3],2,quantile,probs=c(0.025)))*100
  prev_high<-(apply(posts_hiv$fitted_output[,,3],2,quantile,probs=c(0.975)))*100
  
  
  incidence_median<-apply(posts_hiv$fitted_output[,,2],2,median)
  incidence_low<-apply(posts_hiv$fitted_output[,,2],2,quantile,probs=c(0.025))
  incidence_high<-apply(posts_hiv$fitted_output[,,2],2,quantile,probs=c(0.975))
  
  r_median<-apply(posts_hiv$fitted_output[,,1],2,median)
  r_low<-apply(posts_hiv$fitted_output[,,1],2,quantile,probs=c(0.025))
  r_high<-apply(posts_hiv$fitted_output[,,1],2,quantile,probs=c(0.975))
  
  # Combine into two data frames for plotting
  #df_sample = data.frame(sample_prop, sample_time)
  df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, xout )
  names(df_fit_prevalence)<-c("median","low","high","time")
  df_fit_prevalence$credible_skew<-(df_fit_prevalence$high - df_fit_prevalence$median) - (df_fit_prevalence$median - df_fit_prevalence$low)
  
  df_fit_incidence<-data.frame(incidence_low,incidence_median,incidence_high,xout)
  names(df_fit_incidence)<-c("low","median","high","time")
  
  r_fit<-data.frame(r_low,r_median,r_high,xout)
  names(r_fit)<-c("low","median","high","time")
  # Plot the synthetic data with the model predictions
  # Median and 95% Credible Interval
  
  
 
  return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
              r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist))
  
  
}

xout<-seq(1970,2020.1,0.1)


stan_output_random_walk_second_n_100<-plot_stan_model_fit(model_output = mod_hiv_prev,sim_sample = sample_df_100_random_second,
                                                          plot_name = "Random walk second order, n = 100",xout = xout,
                                                          sim_output = sim_model_output$sim_df)

prev_df<-stan_output_random_walk_second_n_100$df_output
prev_df$iteration<-rep(i,nrow(prev_df))

incidence_df<-stan_output_random_walk_second_n_100$incidence_df
incidence_df$iteration<-rep(i,nrow(incidence_df))


kappa_df<-stan_output_random_walk_second_n_100$r_fit_df
kappa_df$iteration<-rep(i,nrow(kappa_df))


iota_value<-stan_output_random_walk_second_n_100$iota_value
iota_value$iteration<-rep(i,nrow(iota_value))

prev_df_tot<-rbind(prev_df_tot,prev_df)

incidence_df_tot<-rbind(incidence_df_tot,incidence_df)

kappa_df_tot<-rbind(kappa_df_tot,kappa_df)

iota_values_tot<-rbind(iota_values_tot,iota_value)

sample_df_100_random_second$iteration<-rep(i,nrow(sample_df_100_random_second))

sample_df_tot<-rbind(sample_df_tot,sample_df_100_random_second)

plot(sample_df_100_random_second$sample_prev_hiv_percentage,colour="red")
lines(prev_df$median,colour="midnightblue")

second_timo<-Sys.time()

print( )

print((i/iterations)*100)

print("length of this iteration:",quote = F)

print(second_timo-timo)

print("How Long to go:",quote = F)

timo_diff<-second_timo - timo

total_timo<-c(total_timo,timo_diff)

if(i <= 20){

print((second_timo-timo) * (iterations - i))
}

if(i > 20){
  
  print(mean(total_timo) * (iterations - i) )
}


}

random_walk_first_order_n_100_complete<-list(prev=prev_df_tot,incidence=incidence_df_tot,
                                             kappa=kappa_df_tot,sampled=sample_df_tot,iota=iota_values_tot)


save(random_walk_first_order_n_100_complete,file = "../stan_objects_from_simpleepp_R/loop_RW_first_n_100_complete")

prev_df_tot$time[503]

mean_value_function<-function(iterations,nrow_per_iteration,data_frame){
data_prev<-NULL
iter_value<-iterations - 1
nrow_value<-nrow_per_iteration
for(i in 1:nrow_value){
  values<-data_frame$median[0:iter_value*nrow_value+i]
  low<-data_frame$low[0:iter_value*nrow_value+i]
  high<-data_frame$high[0:iter_value*nrow_value+i]
  time_of_values<-data_frame$time[0:iter_value*nrow_value+i]
  
  row<-cbind(mean(low),mean(values),mean(high),mean(time_of_values))
  
  data_prev<-rbind.data.frame(data_prev,row)
  
}

names(data_prev)<-c("low","median","high","time")

return(data_prev)

}

a<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = prev_df_tot)
plot(a$median)


prev_df_tot$median[0:99*502+502]
prev_df_tot[0:99*502+502,]

plot(data_prev[,2])
lines(sim_model_output$sim_df$prev_percent,col="red")
