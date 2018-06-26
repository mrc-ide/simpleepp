##################################################################################################################################
## Modelling incidence as a spline/ RW now #######################################################################################
##################################################################################################################################

require(splines)
require(rstan)
require(ggplot2)
require(reshape2)

rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores())

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/simpleepp_foi.stan")

bell_curve_func<-function(a,b,c,x){
  f_t<-a*exp(-((x-b)^2)/(2*c^2))
  
}
bell_curve_func(0.1,0.2,0.3,seq(0,0.998,0.002))

sim_model_output$sim_df$lambda
f_t<-bell_curve_func(0.1,0.2,0.3,seq(0,0.998,0.002))

f_t<-sim_model_output$sim_df$lambda
f_t[1]<-f_t[2]
iota<-0.005
mu<-1/35
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
dt<-0.1
foi_flag<-1

test<-simpleepp_foi(f_t = f_t,iota = iota,mu = mu,mu_i=mu_i,sigma = sigma,dt = dt,foi_flag = foi_flag,kappa_init = 0.5)
test
plot(test[,3])
plot(test[,2])
plot(test[,1])
lines(test[,1])

run_simulated_model<-function(params,times,inc_vector){
  
  mu <- params$mu                               # Non HIV mortality / exit from population
  sigma <- params$sigma           # Progression from stages of infection
  mu_i <- params$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params$iota
  foi_flag<-params$foi_flag
  
  
  dt <- times$dt                                     # time step
  nsteps <- as.integer(times$years/dt)                   # number of steps
  xstart <- times$start                                # the start of the epidemic
  step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  foi_flag<-1
  f_t<-inc_vector
  f_t[1]<-f_t[2]
  f_t<-f_t[-501]
  mod<-simpleepp_foi(f_t = f_t,iota = iota,mu = mu,mu_i=mu_i,sigma = sigma,dt = dt,foi_flag = foi_flag,kappa_init = params$kappa_init)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa))
  
  
}

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/no_art_simpleepp/original_data_run",verbose = T)
f_t<-sim_model_output$sim_df$lambda
f_t[1]<-f_t[2]
iota<-0.005
mu<-1/35
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
dt<-0.1
foi_flag<-1

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa_init=0.5,iota=iota,foi_flag=1)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_foi<-run_simulated_model(params_sim,times_sim,sim_model_output$sim_df$lambda)

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

plotted_sim<-sim_plot(sim_model_foi$sim_df)
plot(plotted_sim$whole)

# save(sim_model_foi,file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/true_epidemic")


###################################################################################################################################
## Now we have simulated through our model we can extract some samples form it ####################################################
###################################################################################################################################

sample_range<-1970:2015
sample_years<-46
sample_n<-1000


sample_function<-function(year_range,number_of_years_to_sample,people_t0_sample,simulated_df,prevalence_column_id){
  sample_years_hiv <- number_of_years_to_sample # number of days sampled throughout the epidemic
  sample_n <- people_t0_sample # number of host individuals sampled per day
  
  # Choose which days the samples were taken. 
  # Ideally this would be daily, but we all know that is difficult.
  sample_time_hiv = sort(sample(year_range, sample_years_hiv, replace=F))
  
  # Extract the "true" fraction of the population that is infected on each of the sampled days:
  sample_propinf_hiv = simulated_df[simulated_df$time %in% sample_time_hiv, prevalence_column_id]
  
  ## this just samples our prevalence, to get a probability that the sample we take is HIV infected then we need to divide
  ## by 100
  
  #sample_propinf_hiv<-sample_propinf_hiv/100
  
  # Generate binomially distributed data.
  # So, on each day we sample a given number of people (sample_n), and measure how many are infected.
  # We expect binomially distributed error in this estimate, hence the random number generation.
  sample_y_hiv_prev = rbinom(sample_years_hiv, sample_n, sample_propinf_hiv)
  sample_prev_hiv_percentage<-(sample_y_hiv_prev/sample_n)*100
  
  ## lets have a ggplot of the y (infected) and out sample of Y over time 
  sample_df_100<-data.frame(cbind(sample_time_hiv,sample_prev_hiv_percentage,sample_y_hiv_prev,sample_time_hiv,sample_propinf_hiv))
  return(sample_df_100)  
}



sample_df_1000_second<-sample_function(sample_range,sample_years,sample_n,
                                       simulated_df = sim_model_foi$sim_df,prevalence_column_id = 3)




plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = sim_model_foi$sim_df,sample_df = sample_df_1000_second)


ggplot(data = sample_df_1000_second,aes(x=sample_time_hiv,y=sample_prev_hiv_percentage))+geom_point(colour="red",size=1.5)

###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################

splines_creator_rw_o_splino<-function(knot_number,penalty_order,type,step_vector){
  
  if(type == "spline"){
    mat_ord<-4
    nk <- knot_number # number of knots
    dk <- diff(range(xout))/(nk-3)
    knots <- xstart + -3:nk*dk
    spline_matrix<- splineDesign(knots, step_vector, ord=mat_ord)
    penalty_matrix <- diff(diag(nk), diff=penalty_order)
  }  else{
    mat_ord<-2
    spline_matrix<-splineDesign(1969:2021,step_vector,ord = mat_ord)
    penalty_matrix <- diff(diag(ncol(spline_matrix)), diff=penalty_order)
    
  }
  
  
  
  return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
  
}

knot_number= 51
penalty_order= 2
typo<-"rw"

splines_matrices<-splines_creator_rw_o_splino(knot_number,penalty_order,type = typo,step_vector = seq(1970.1,2020,0.1))
rows_to_evaluate<-0:45*10+1


stan_data_discrete_foi<-list(
  n_obs = sample_years,
  n_sample = sample_n,
  y = as.array(sample_df_1000_second$sample_y_hiv_prev),
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  knot_number = knot_number,
  estimate_years = 5,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  foi_flag = foi_flag,
  kappa_init = 0.5
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

test_stan_hiv<- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/simpleepp_foi.stan",
                     data = stan_data_discrete_foi,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/simpleepp_foi.stan",
                     data = stan_data_discrete_foi, pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.85))




plot_stan_model_fit<-function(model_output,sim_sample,plot_name,xout,sim_output){
  
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
  
  beta_median<-apply(posts_hiv$beta,2,median)
  beta_low<-apply(posts_hiv$beta,2,quantile,probs=c(0.025))
  beta_high<-apply(posts_hiv$beta,2,quantile,probs=c(0.975))
  beta_df<-rbind.data.frame(beta_low,beta_median,beta_high)
  names(beta_df)<-c("1","2","3","4")
  
  
  
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
  
  
  plotter<-ggplot(sim_sample) +
    geom_point(aes(x=sample_time_hiv.1, y=sample_prev_hiv_percentage),col="red", shape = 19, size = 1.5) +
    geom_line(data = df_fit_prevalence, aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),
                colour="midnightblue",alpha=0.2,fill="midnightblue")+
    geom_line(data = sim_output,aes(x=time,y=prev_percent),colour="yellow",size=1)+
    coord_cartesian(xlim=c(1965,2025))+labs(x="Time",y="Prevalence (%)", title=plot_name)
  
  incidence_plot<-ggplot(data=df_fit_incidence)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",alpha=0.2,colour="midnightblue")+
    geom_line(data = sim_output,aes(x=time,y=lambda),colour="yellow",size=1)+
    labs(x="time",y="incidence",title="incidence_plot")
  
  r_plot<-ggplot(data = r_fit)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",colour="midnightblue",alpha=0.2)+
    geom_line(data = sim_output,aes(x=time,y=kappa),colour="yellow",size=1)
  labs(x="Time",y="r value through time",title="R logistic through time")
  
  return(list(prevalence_plot=(plotter),df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
              r_fit_df=r_fit,incidence_plot=incidence_plot,r_plot=r_plot,sigma_pen_values=sigma_df,iota_value=params_df,
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,beta_values=beta_df))
  
  
}

xout<-seq(1970,2020,0.1)


stan_output_second_order_n_1000<-plot_stan_model_fit(model_output = mod_hiv_prev,
                                                     sim_sample = sample_df_1000_second,plot_name = "Spline Second order, n = 1000",xout = xout,
                                                     sim_output = sim_model_foi$sim_df)
stan_output_second_order_n_1000$prevalence_plot
stan_output_second_order_n_1000$incidence_plot
stan_output_second_order_n_1000$r_plot

spline_first_order_plots<-list(stan_output_second_order_n_1000$prevalence_plot,
                               stan_output_second_order_n_1000$incidence_plot,
                               stan_output_second_order_n_1000$r_plot)
spline_second_order_plots<-list(stan_output_second_order_n_1000$prevalence_plot,
                                stan_output_second_order_n_1000$incidence_plot,
                                stan_output_second_order_n_1000$r_plot)
rw_first_order_plots<-list(stan_output_second_order_n_1000$prevalence_plot,
                           stan_output_second_order_n_1000$incidence_plot,
                           stan_output_second_order_n_1000$r_plot)
rw_second_order_plots<-list(stan_output_second_order_n_1000$prevalence_plot,
                            stan_output_second_order_n_1000$incidence_plot,
                            stan_output_second_order_n_1000$r_plot)

spline_first_order_plots[[1]]$labels$title<-"Spline first order n = 500"
spline_second_order_plots[[1]]$labels$title<-"Spline second order n = 500"
rw_first_order_plots[[1]]$labels$title<-" RW first order n = 500"
rw_second_order_plots[[1]]$labels$title<-"RW second order n = 500"

require(ggpubr)

prev<-ggarrange(spline_first_order_plots[[1]],spline_second_order_plots[[1]],
                rw_first_order_plots[[1]],rw_second_order_plots[[1]],ncol = 2,nrow = 2)
inc<-ggarrange(spline_first_order_plots[[2]],spline_second_order_plots[[2]],
               rw_first_order_plots[[2]],rw_second_order_plots[[2]],ncol = 2,nrow = 2)
kappa<-ggarrange(spline_first_order_plots[[3]],spline_second_order_plots[[3]],
                 rw_first_order_plots[[3]],rw_second_order_plots[[3]],ncol = 2,nrow = 2)

prev
inc
kappa

########################################################################################################################
## Creating the simulated datasets from which we will fit our model too ###############################################
#######################################################################################################################
require(ggplot2)
require(ggpubr)
require(reshape2)
require(rstan)
require(devtools)
require(splines)

find_rtools()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

expose_stan_functions("stan_files/chunks/cd4_spline_model.stan")

########################################################################################################################################
## Run the simulated model to get the output ###########################################################################################
########################################################################################################################################


run_simulated_model_foi<-function(params,times,inc_vector){
  
  mu <- params$mu                               # Non HIV mortality / exit from population
  sigma <- params$sigma           # Progression from stages of infection
  mu_i <- params$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params$iota
  foi_flag<-params$foi_flag
  
  
  dt <- times$dt                                     # time step
  nsteps <- as.integer(times$years/dt)                   # number of steps
  xstart <- times$start                                # the start of the epidemic
  step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  foi_flag<-1
  f_t<-inc_vector
  f_t[1]<-f_t[2]
  f_t<-f_t[-501]
  mod<-simpleepp_foi(f_t = f_t,iota = iota,mu = mu,mu_i=mu_i,sigma = sigma,dt = dt,foi_flag = foi_flag,kappa_init = params$kappa_init)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa))
  
  
}

f_t<-sim_model_output$sim_df$lambda
f_t[1]<-f_t[2]
iota<-0.005
mu<-1/35
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
dt<-0.1
foi_flag<-1

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa_init=0.5,iota=iota,foi_flag=1)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_foi<-run_simulated_model_foi(params_sim,times_sim,sim_model_output$sim_df$lambda)

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

plotted_sim<-sim_plot(sim_model_foi$sim_df)
plot(plotted_sim$whole)

#######################################################################################################################################
## Now lets form our loop to sample through and fit to stan data ######################################################################
#######################################################################################################################################

sample_range<-1970:2015
## Need to change this to length of time_points_to_sample if sporadic
sample_n<5000

exponential_decay_function<-function(N0,t,lambda){
  
  nt<-N0*exp(-t*lambda)
  
  return(nt)
}

decay<- 0  #exponential_decay_function(N0 = 90,t=(seq(0,50,0.1)),lambda = 0.05) # Must be put to equal 0 for no underreporting

penalty_order<-1

time_points_to_sample<- 0 #seq(1970,2015,3)                                   ## Must also be equal to 0 for complete reporting 

rows_to_evaluate<- 0:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

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
                                               simulated_df = sim_model_foi$sim_df,
                                               prevalence_column_id = 3,decay = decay,
                                               time_points_to_sample = time_points_to_sample)
  
  sample_df_100_random_second$iteration<-rep(i,nrow(sample_df_100_random_second))
  
  sample_df_tot<-rbind.data.frame(sample_df_tot,sample_df_100_random_second)
  
  print(i/iterations * 100)
  
}

plot(sample_df_tot[sample_df_tot$iteration==19,2])
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
  Sys.sleep(0.33)
  print(i)
}



sampled_n_5000_complete_data_foi_model<-sample_df_tot

save(sampled_n_5000_complete_data_foi_model,
     file = "C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_5000_sampled_data")  
save(sim_model_output_changed_to_bell_curve,
     file = "C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/true_epidemic_data_double_peak")




