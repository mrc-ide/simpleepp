############################################################################################
## Altering the R scripts into a series of functions to run on the cluster #################
############################################################################################

## This function creates the true epidemic curve

make_true_epidemic <- function(){

rstan::expose_stan_functions("simpleepp/stan_files/chunks/cd4_matrix_random_walk.stan")
 
  mu <- 1/35                               # Non HIV mortality / exit from population
  sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
  mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  kappa<-c(0.5,0.1,0.3,1995)
  iota<-0.0001
  
  dt <- 0.1                                     # time step
  nsteps <- base::as.integer(50/dt)                   # number of steps
  xstart <- 1970                                # the start of the epidemic
  step_vector <- base::seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)                #the vector of steps in total
  
  
  
  params_sim <- base::list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota)
  times_sim <- base::list(dt=dt,years=50,start=xstart)
  
  
  mu <- params_sim$mu                               # Non HIV mortality / exit from population
  sigma <- params_sim$sigma           # Progression from stages of infection
  mu_i <- params_sim$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params_sim$iota
  
  
  dt <- times_sim$dt                                     # time step
  nsteps <- base::as.integer(times_sim$years/dt)                   # number of steps
  xstart <- times_sim$start                                # the start of the epidemic
  step_vector <- base::seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  rlogistic <- function(t, p) {
    p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
  }
  
  kappa_params <- c(log(params_sim$kappa[1]),log(params_sim$kappa[2]),params_sim$kappa[3],params_sim$kappa[4])
  kappa <- exp(rlogistic(step_vector,kappa_params))
  
  mod <- simpleepp_no_art(kappa = kappa,iota,mu,sigma,mu_i,dt)
  sim_mod <- data.frame(mod)
  names(sim_mod) <- c("kappa","lambda","prevalence")
  sim_mod$prev_percent <- sim_mod$prevalence * 100
  sim_mod$time <- c(xstart,step_vector)
  
  return(list(sim_df = sim_mod,kappa_values = kappa))
  
  
}

## This function then fits the data
## sample_data_frame has to be the full length 100 set sampled from true epidemic
## iteration number is an integer between 1:100
## data_about_sampling is a list containg: penalty_order, sample_years (this is 46), sample_n for the dataset and Rows_to_evaluate
## params is another list containing mu, sigma and mu_i from the simulated true epidemic 
## simulated true_df is the output from the above function 


fitting_data_function<-function(samples_data_frame,iteration_number,data_about_sampling,params,simulated_true_df){

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sample_df_100_random_second<-samples_data_frame[samples_data_frame$iteration==iteration_number,]



xout<-seq(1970,2020,0.1)
spline_matrix<-splines::splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=data_about_sampling$penalty_order)        ## This matrix creates the differences between your kappa values 




stan_data_discrete<-list(
  n_obs = data_about_sampling$sample_years,
  n_sample = data_about_sampling$sample_n,
  y = as.array(sample_df_100_random_second$sample_y_hiv_prev),
  time_steps_euler = 501,
  penalty_order = data_about_sampling$penalty_order,
  estimate_period = 5,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = params$mu,
  sigma = params$sigma,
  mu_i = params$mu_i,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(data_about_sampling$rows_to_evaluate)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  


mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/cd4_matrix_random_walk.stan", data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.85))



plot_stan_model_fit<-function(model_output,sim_sample,sim_output,plot_name,xout){
  
  posts_hiv <- rstan::extract(model_output)
  
  
  iota_dist<-posts_hiv$iota
  params<-stats::median(posts_hiv$iota)
  params_low<-stats::quantile(posts_hiv$iota,c(0.025))
  params_high<-stats::quantile(posts_hiv$iota,c(0.975))
  
  params_df<-rbind.data.frame(params_low,params,params_high)
  names(params_df)<-c("iota")
  
  sigma_pen_dist<-posts_hiv$sigma_pen
  sigma_values<-stats::median(posts_hiv$sigma_pen)
  sigma_low<-stats::quantile(posts_hiv$sigma_pen,c(0.025))
  sigma_high<-stats::quantile(posts_hiv$sigma_pen,probs=c(0.975))
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
  
  prev_median<-(apply(posts_hiv$fitted_output[,,3],2,stats::median))*100
  prev_low<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.025)))*100
  prev_high<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.975)))*100
  
  
  incidence_median<-apply(posts_hiv$fitted_output[,,2],2,stats::median)
  incidence_low<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.025))
  incidence_high<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.975))
  
  r_median<-apply(posts_hiv$fitted_output[,,1],2,stats::median)
  r_low<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.025))
  r_high<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.975))
  
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
                                                          sim_output = simulated_true_df$sim_df)

prev_df<-stan_output_random_walk_second_n_100$df_output
prev_df$iteration<-rep(iteration_number,nrow(prev_df))

incidence_df<-stan_output_random_walk_second_n_100$incidence_df
incidence_df$iteration<-rep(iteration_number,nrow(incidence_df))


kappa_df<-stan_output_random_walk_second_n_100$r_fit_df
kappa_df$iteration<-rep(iteration_number,nrow(kappa_df))


iota_value<-stan_output_random_walk_second_n_100$iota_value
iota_value$iteration<-rep(iteration_number,nrow(iota_value))

return(list(prevalence_df = prev_df, incidence_df = incidence_df, transmission_parameter_df = kappa_df, starting_fraction = iota_value))

}

## now lets create a function that will iteritively store the above data 

fitting_data_function_loop<-function(samples_data_frame,iteration_number,data_about_sampling,params,simulated_true_df){
  
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
 
 
  xout<-seq(1970,2020,0.1)
  spline_matrix<-splines::splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
  penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=data_about_sampling$penalty_order)        ## This matrix creates the differences between your kappa values 
  
  plot_stan_model_fit<-function(model_output,sim_sample,sim_output,plot_name,xout){
    
    posts_hiv <- rstan::extract(model_output)
    
    
    iota_dist<-posts_hiv$iota
    params<-stats::median(posts_hiv$iota)
    params_low<-stats::quantile(posts_hiv$iota,c(0.025))
    params_high<-stats::quantile(posts_hiv$iota,c(0.975))
    
    params_df<-rbind.data.frame(params_low,params,params_high)
    names(params_df)<-c("iota")
    
    sigma_pen_dist<-posts_hiv$sigma_pen
    sigma_values<-stats::median(posts_hiv$sigma_pen)
    sigma_low<-stats::quantile(posts_hiv$sigma_pen,c(0.025))
    sigma_high<-stats::quantile(posts_hiv$sigma_pen,probs=c(0.975))
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
    
    prev_median<-(apply(posts_hiv$fitted_output[,,3],2,stats::median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,2],2,stats::median)
    incidence_low<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.975))
    
    r_median<-apply(posts_hiv$fitted_output[,,1],2,stats::median)
    r_low<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.025))
    r_high<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.975))
    
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
  
  prev_df_tot <- NULL
  incidence_df_tot <- NULL
  kappa_df_tot<-NULL
  iota_values_tot<-NULL
  
  
  
  for( i in 1:iteration_number){
  
  
  sample_df_100_random_second<-samples_data_frame[samples_data_frame$iteration==i,]
  
  
  stan_data_discrete<-list(
    n_obs = data_about_sampling$sample_years,
    n_sample = data_about_sampling$sample_n,
    y = as.array(sample_df_100_random_second$sample_y_hiv_prev),
    time_steps_euler = 501,
    penalty_order = data_about_sampling$penalty_order,
    estimate_period = 5,
    time_steps_year = 51,
    X_design = spline_matrix,
    D_penalty = penalty_matrix,
    mu = params$mu,
    sigma = params$sigma,
    mu_i = params$mu_i,
    dt = 1,
    dt_2 = 0.1,
    rows_to_interpret = as.array(data_about_sampling$rows_to_evaluate)
  )
  
  params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  
  
  
  mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/cd4_matrix_random_walk.stan", data = stan_data_discrete,
                              pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
  
  
  

  
  xout<-seq(1970,2020.1,0.1)
  
  
  stan_output_random_walk_second_n_100<-plot_stan_model_fit(model_output = mod_hiv_prev,sim_sample = sample_df_100_random_second,
                                                            plot_name = "Random walk second order, n = 100",xout = xout,
                                                            sim_output = simulated_true_df)
  
  prev_df<-stan_output_random_walk_second_n_100$df_output
  prev_df$iteration<-rep(i,nrow(prev_df))
  
  incidence_df<-stan_output_random_walk_second_n_100$incidence_df
  incidence_df$iteration<-rep(i,nrow(incidence_df))
  
  
  kappa_df<-stan_output_random_walk_second_n_100$r_fit_df
  kappa_df$iteration<-rep(i,nrow(kappa_df))
  
  
  iota_value<-stan_output_random_walk_second_n_100$iota_value
  iota_value$iteration<-rep(i,nrow(iota_value))
  
  
  prev_df_tot<-rbind(prev_df_tot,prev_df)
  incidence_df_tot <- rbind(incidence_df_tot,incidence_df)
  kappa_df_tot <- rbind(kappa_df_tot,kappa_df)
  iota_values_tot <- rbind(iota_values_tot,iota_value)
  }
  
  data_about_sampling$simul_type<-"Random_Walk"
  
  return(list(prev=prev_df_tot,incidence=incidence_df_tot,kappa=kappa_df_tot,iota_values=iota_values_tot,
              data_about_run = data_about_sampling))
  
}

fitting_data_function_spline_loop<-function(samples_data_frame,iteration_number,data_about_sampling,params,simulated_true_df){
  
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  steppy<-seq(1970.1,by=0.1,length.out = 500)
  
  splines_creator<-function(knot_number,penalty_order,step_vector){
    
    nk <- knot_number # number of knots
    dk <- diff(range(xout))/(nk-3)
    knots <- samples_data_frame$sample_time_hiv[1] + -3:nk*dk
    spline_matrix<- splineDesign(knots, step_vector, ord=4)
    penalty_matrix <- diff(diag(nk), diff=penalty_order)
    
    return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
    
  }                                   ## This creates our spline matrix
  
  knot_number= data_about_sampling$knot_number                                                ## Our number of knots
  penalty_order= data_about_sampling$penalty_order                                            ## Our penalty order defined on input
  
  xout<-seq(1970,2020,0.1)
  
  splines_matrices<-splines_creator(knot_number,penalty_order,steppy)                                ## Creates the penalty and the spline matrix
  
  rows_to_evaluate<-data_about_sampling$rows_to_evaluate                                      ## Rows for y_hat to fit to
  
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
    

    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,beta_values=beta_df))
    
    
  }       ## Extracts the results from stan fit
  
  
  
  prev_df_tot <- NULL                                                                         ## INitialize the data frames to take loop results into
  incidence_df_tot <- NULL                                                                   
  kappa_df_tot<-NULL
  iota_values_tot<-NULL
  beta_values_tot<-NULL
  
  for(i in 1:iteration_number){
    
    sample_df_1000_second<-samples_data_frame[samples_data_frame$iteration==i,]               ## Change the dataset each iteration
  
  stan_data_discrete<-list(
    n_obs = data_about_sampling$sample_years,
    n_sample = data_about_sampling$sample_n ,
    y = as.array(sample_df_1000_second$sample_y_hiv_prev),
    time_steps_euler = length(xout),
    penalty_order = data_about_sampling$penalty_order ,
    knot_number = knot_number,
    estimate_years = 5,
    time_steps_year = 51,
    X_design = splines_matrices$spline_matrix,
    D_penalty = splines_matrices$penalty_matrix,
    mu = params$mu,
    sigma = params$sigma,
    mu_i = params$mu_i,
    dt = 1,
    dt_2 = 0.1,
    rows_to_interpret = as.array(rows_to_evaluate)
  )
  
  params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  
  

  mod_hiv_prev <- stan("simpleepp/stan_files/chunks/cd4_spline_model.stan", data = stan_data_discrete,
                       pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                       control = list(adapt_delta = 0.85))
  
  
  stan_output_second_order_n_1000<-plot_stan_model_fit(model_output = mod_hiv_prev,
                                                       sim_sample = sample_df_1000_second,plot_name = "Spline Second order, n = 1000",xout = xout,
                                                       sim_output = simulated_true_df)
  
  

  prev_df<-stan_output_second_order_n_1000$df_output
  prev_df$iteration<-rep(i,nrow(prev_df))
  
  incidence_df<-stan_output_second_order_n_1000$incidence_df
  incidence_df$iteration<-rep(i,nrow(incidence_df))
  
  
  kappa_df<-stan_output_second_order_n_1000$r_fit_df
  kappa_df$iteration<-rep(i,nrow(kappa_df))
  
  
  iota_value<-stan_output_second_order_n_1000$iota_value
  iota_value$iteration<-rep(i,nrow(iota_value))
  
  beta_values<-stan_output_second_order_n_1000$beta_values
  beta_values$iteration<-rep(i,nrow(beta_values))
  
  
  
  prev_df_tot<-rbind(prev_df_tot,prev_df)
  incidence_df_tot <- rbind(incidence_df_tot,incidence_df)
  kappa_df_tot <- rbind(kappa_df_tot,kappa_df)
  iota_values_tot <- rbind(iota_values_tot,iota_value)
  beta_values_tot <- rbind(beta_values_tot, beta_values)
  }
  
  data_about_sampling$simul_type<-"Spline"
  
  return(list(prev=prev_df_tot,incidence=incidence_df_tot,kappa=kappa_df_tot,iota_values=iota_values_tot,beta_vals=beta_values_tot,
              data_about_run = data_about_sampling))
  
}

