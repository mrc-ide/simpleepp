#################################################################################################################################
## Playing around with some functions for double peaks in the transmission parameter ############################################
#################################################################################################################################

require(splines)
require(rstan)
require(ggplot2)
require(reshape2)
require(grid)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/cd4_spline_model.stan")

#######################################################################################################################################
## Now I'll create the R script which to call the stan script for sampling with a first order random walk value #######################
#######################################################################################################################################

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
  
  splines_double_peak<-function(knot_number,step_vector,beta){
    
    mat_ord<-4
    nk <- knot_number # number of knots
    dk <- diff(range(xout))/(nk-3)
    knots <- xstart + -3:nk*dk
    spline_matrix<- splineDesign(knots, step_vector, ord=mat_ord)
    
    curve<- spline_matrix %*% beta
    plot_curve<-plot(curve)
    
    
    return(list(spline_matrix=spline_matrix,curve=curve,plot=plot_curve))
    
  }
  
  
 
  run_through<-splines_double_peak(params$knot_number,step_vector,params$beta)
  kappa<-run_through$curve
  
  #kappa<-bell_curve_func(0.5,0.2,0.2,seq(0,0.998,0.002))
  
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
iota<-0.0001
knot_number<-7
beta<-c(0.70,0.7,0.6,0.3,0.45,0.45,0.45)

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,knot_number=knot_number,iota=iota,beta=beta)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_output_changed_to_bell_curve<-run_simulated_model(params_sim,times_sim)

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

plotted_sim<-sim_plot(sim_model_output_changed_to_bell_curve$sim_df)
plot(plotted_sim$whole)

plotted_sim$incidence

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
                                       simulated_df = sim_model_output_changed_to_bell_curve$sim_df,
                                       prevalence_column_id = 3)




plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = sim_model_output_changed_to_bell_curve$sim_df,
            sample_df = sample_df_1000_second)


ggplot(data = sample_df_1000_second,aes(x=sample_time_hiv,y=sample_prev_hiv_percentage))+geom_point(colour="red",size=1.5)

###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################

splines_creator<-function(knot_number,penalty_order){
  
  nk <- knot_number # number of knots
  dk <- diff(range(xout))/(nk-3)
  knots <- xstart + -3:nk*dk
  spline_matrix<- splineDesign(knots, step_vector, ord=4)
  penalty_matrix <- diff(diag(nk), diff=penalty_order)
  
  return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
  
}

knot_number= 7
penalty_order= 2

splines_matrices<-splines_creator(knot_number,penalty_order)
rows_to_evaluate<-0:45*10+1


stan_data_discrete<-list(
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
  rows_to_interpret = as.array(rows_to_evaluate)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

test_stan_hiv<- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/cd4_spline_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/cd4_spline_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
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
  names(beta_df)<-c("1","2","3","4","5","6","7")
  
  
  
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
  
  return(list(prevalence_plot=(plotter),inits=inits,df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
              r_fit_df=r_fit,incidence_plot=incidence_plot,r_plot=r_plot,sigma_pen_values=sigma_df,iota_value=params_df,
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,beta_values=beta_df))
  
  
}

xout<-seq(1970,2020,0.1)


stan_output_second_order_n_1000<-plot_stan_model_fit(model_output = mod_hiv_prev,
                                                     sim_sample = sample_df_1000_second,plot_name = "Spline Second order, n = 1000",xout = xout,
                                                     sim_output = sim_model_output_changed_to_bell_curve$sim_df)
stan_output_second_order_n_1000$prevalence_plot
stan_output_second_order_n_1000$incidence_plot
stan_output_second_order_n_1000$r_plot
