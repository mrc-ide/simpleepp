fitting_data_function_loop_foi<-function(samples_data_frame,iteration_number,data_about_sampling,params){
  
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  xout<-seq(1970,2020.1,0.1)
  
  
  splines_creator_rw_o_splino<-function(knot_number,penalty_order,type,step_vector){
    
    if(type == "spline"){
      mat_ord<-4
      nk <- knot_number # number of knots
      dk <- diff(range(xout))/(nk-3)
      knots <- 1970 + -3:nk*dk
      spline_matrix<- splineDesign(knots, step_vector, ord=mat_ord)
      penalty_matrix <- diff(diag(nk), diff=penalty_order)
    }  else{
      mat_ord<-2
      spline_matrix<-splineDesign(1969:2021,step_vector,ord = mat_ord)
      penalty_matrix <- diff(diag(ncol(spline_matrix)), diff=penalty_order)
      
    }
    
    
    
    return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
    
  }
  
  knot_number= data_about_sampling$knot_number
  penalty_order= data_about_sampling$penalty_order
  typo<- data_about_sampling$model_type
  
  splines_matrices<-splines_creator_rw_o_splino(knot_number,penalty_order,type = typo,
                                                step_vector = seq(1970.1,2020,0.1))
  rows_to_evaluate<-0:45*10+1
  
 
 
  
  plot_stan_model_fit<-function(model_output,sim_sample,plot_name,xout){
    
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
    xout<-seq(1970,2020,0.1)
    
    stan_data_discrete_foi<-list(
      n_obs = data_about_sampling$sample_years,
      n_sample = data_about_sampling$sample_n,
      y = as.array(sample_df_100_random_second$sample_y_hiv_prev),
      time_steps_euler = length(xout),
      penalty_order = data_about_sampling$penalty_order,
      knot_number = data_about_sampling$knot_number,
      estimate_years = 5,
      time_steps_year = 51,
      X_design = splines_matrices$spline_matrix,
      D_penalty = splines_matrices$penalty_matrix,
      mu = params$mu,
      sigma = params$sigma,
      mu_i = params$mu_i,
      dt = 1,
      dt_2 = 0.1,
      rows_to_interpret = as.array(rows_to_evaluate),
      foi_flag = params$foi_flag,
      kappa_init = 0.5
    )
    params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  
    
    
    mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/simpleepp_foi.stan", 
                                data = stan_data_discrete_foi,
                                pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
    
    
    stan_output_random_walk_second_n_100<-plot_stan_model_fit(model_output = mod_hiv_prev,
                                                              sim_sample = sample_df_100_random_second,
                                                              plot_name = "Random walk second order, n = 100",
                                                              xout = xout)
    
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
  
  data_about_sampling$simul_type<-data_about_sampling$type
  
  return(list(prev=prev_df_tot,incidence=incidence_df_tot,kappa=kappa_df_tot,iota_values=iota_values_tot,
              data_about_run = data_about_sampling))
  
}

