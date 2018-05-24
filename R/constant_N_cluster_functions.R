count_data_fitting_spline<-function(simulated_dataset,stan_data,params_used_for_sim){
  
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  
  stan_data_discrete<-stan_data
  
  params_monitor_stan<-c("y_hat","iota","fitted_output","iota","beta","sigma_pen")
  
  mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/ART_mod_spline_COUNT_constant_N.stan", data = stan_data_discrete,
                              pars = params_monitor_stan,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
  
  #############################################################################################################################
  ## SO that has produced our model fit, we will now collect some summary statistics from the HMC fitting process #############
  #############################################################################################################################
  
  test_summary_fitted<-rstan::summary(mod_hiv_prev,pars=c("fitted_output"))$summary
  test_summary_beta<-rstan::summary(mod_hiv_prev,pars=c("beta"))$summary
  
  test_summary_fitted<-data.frame(test_summary_fitted)
  test_summary_beta<-data.frame(test_summary_beta)
  
  r_hat_fitted_summary<-base::summary(test_summary_fitted$Rhat)
  r_hat_beta_summary<-base::summary(test_summary_beta$Rhat)
  
  ### Now we will plot our fits against the sample data 
  
  plot_stan_model_fit<-function(model_output,sim_data){
    
    posts_hiv <- rstan::extract(model_output)
    
    xout<- seq(1970,2020,0.1)                         
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
    
    
    mod_time = seq(1970,2020,0.1)
    
    prev_median<-(apply(posts_hiv$fitted_output[,,7],2,median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,6],2,median)
    incidence_low<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.975))
    
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
    prev_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=prev_percent),colour="red",size=1.05)+
      geom_line(data = df_fit_prevalence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="prevalence %",title="fitted and true prevalence, fitting to count data of ART numbers")  
    
    inc_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=incidence),colour="red",size=1.05)+
      geom_line(data = df_fit_incidence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_incidence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="incidence",title="fitted and true incidence, fitting to count data of ART numbers")  
    
    
    kappa_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=kappa),colour="red",size=1.05)+
      geom_line(data = r_fit,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = r_fit,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="Kappa",title="fitted and true kappa, fitting to count data of ART numbers")  
    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,plot_prev=prev_plot,
                plot_inc=inc_plot,plot_kappa=kappa_plot))
  }
  
  list_of_results<-plot_stan_model_fit(mod_hiv_prev,simulated_dataset)
  
  return(list(fitting_results=list_of_results,simulated_data = simulated_dataset,sim_params=params_used_for_sim,
              rstan_summary_fitted=test_summary_fitted,rstan_summary_beta=test_summary_beta,
              rhat_fitted=r_hat_fitted_summary,rhat_beta=r_hat_beta_summary))
  
  
}

#################################################################################################################################
## NOW FOR THE PREVALENCE SPLINE MODEL *&*^%^(*&^%$£"£$%^&*((*&^%$£$%^&*(((((((((((((((((((((((((()))))))))))))))))))))))))))))##
#################################################################################################################################

prev_data_fitting_spline<-function(simulated_dataset,stan_data,params_used_for_sim){
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  
  stan_data_discrete<-stan_data
  
  params_monitor_stan<-c("y_hat","iota","fitted_output","iota","beta","sigma_pen")
  
  mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/ART_mod_spline_prevalence_constant_N.stan", data = stan_data_discrete,
                              pars = params_monitor_stan,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
  
  test_summary_fitted<-rstan::summary(mod_hiv_prev,pars=c("fitted_output"))$summary
  test_summary_beta<-rstan::summary(mod_hiv_prev,pars=c("beta"))$summary
  
  test_summary_fitted<-data.frame(test_summary_fitted)
  test_summary_beta<-data.frame(test_summary_beta)
  
  r_hat_fitted_summary<-base::summary(test_summary_fitted$Rhat)
  r_hat_beta_summary<-base::summary(test_summary_beta$Rhat)
  
  
  plot_stan_model_fit<-function(model_output,sim_data){#,sim_sample,sim_output,plot_name,xout){
    
    posts_hiv <- rstan::extract(model_output)
    
    xout<- seq(1970,2020,0.1)                         
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
    
    
    mod_time = seq(1970,2020,0.1)
    
    prev_median<-(apply(posts_hiv$fitted_output[,,7],2,median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,6],2,median)
    incidence_low<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.975))
    
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
    prev_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=prev_percent),colour="red",size=1.05)+
      geom_line(data = df_fit_prevalence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="prevalence %",title="fitted and true prevalence, fitting to count data of ART numbers")  
    
    inc_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=incidence),colour="red",size=1.05)+
      geom_line(data = df_fit_incidence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_incidence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="incidence",title="fitted and true incidence, fitting to count data of ART numbers")  
    
    
    kappa_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=kappa),colour="red",size=1.05)+
      geom_line(data = r_fit,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = r_fit,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="Kappa",title="fitted and true kappa, fitting to count data of ART numbers")  
    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,plot_prev=prev_plot,
                plot_inc=inc_plot,plot_kappa=kappa_plot))
    
  }
  
  list_of_results<-plot_stan_model_fit(mod_hiv_prev,simulated_dataset)
  
  
  return(list(fitting_results=list_of_results,simulated_data = simulated_dataset,sim_params=params_used_for_sim,
              rstan_summary_fitted=test_summary_fitted,rstan_summary_beta=test_summary_beta,
              rhat_fitted=r_hat_fitted_summary,rhat_beta=r_hat_beta_summary))
}

#################################################################################################################################
## NOW FOR THE RW MODEL STARTINF WITH THE PREVALENCE FITTING FUNCTION ###########################################################
#################################################################################################################################

prev_data_fitting_RW<-function(simulated_dataset,stan_data,params_used_for_sim){
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  
  stan_data_discrete<-stan_data
  
  params_monitor_stan<-c("y_hat","iota","fitted_output","iota","beta","sigma_pen")
  
  mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/ART_mod_RW_Prevalence_constant_N.stan", data = stan_data_discrete,
                              pars = params_monitor_stan,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
  
  test_summary_fitted<-rstan::summary(mod_hiv_prev,pars=c("fitted_output"))$summary
  test_summary_beta<-rstan::summary(mod_hiv_prev,pars=c("beta"))$summary
  
  test_summary_fitted<-data.frame(test_summary_fitted)
  test_summary_beta<-data.frame(test_summary_beta)
  
  r_hat_fitted_summary<-base::summary(test_summary_fitted$Rhat)
  r_hat_beta_summary<-base::summary(test_summary_beta$Rhat)
  
  
  plot_stan_model_fit<-function(model_output,sim_data){
    
    posts_hiv <- rstan::extract(model_output)
    
    xout<- seq(1970,2020,0.1)                         
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
    
    
    mod_time = seq(1970,2020,0.1)
    
    prev_median<-(apply(posts_hiv$fitted_output[,,7],2,median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,6],2,median)
    incidence_low<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.975))
    
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
    
    prev_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=prev_percent),colour="red",size=1.05)+
      geom_line(data = df_fit_prevalence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="prevalence %",title="fitted and true prevalence, fitting to count data of ART numbers")  
    
    inc_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=incidence),colour="red",size=1.05)+
      geom_line(data = df_fit_incidence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_incidence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="incidence",title="fitted and true incidence, fitting to count data of ART numbers")  
    
    
    kappa_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=kappa),colour="red",size=1.05)+
      geom_line(data = r_fit,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = r_fit,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="Kappa",title="fitted and true kappa, fitting to count data of ART numbers")  
    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,plot_prev=prev_plot,
                plot_inc=inc_plot,plot_kappa=kappa_plot))
    
  }
  
  list_of_results<-plot_stan_model_fit(mod_hiv_prev,simulated_dataset)
  
  return(list(fitting_results=list_of_results,simulated_data = simulated_dataset,sim_params=params_used_for_sim,
              rstan_summary_fitted=test_summary_fitted,rstan_summary_beta=test_summary_beta,
              rhat_fitted=r_hat_fitted_summary,rhat_beta=r_hat_beta_summary))
}

################################################################################################################################
## NOW FINALLY FOR THE COUNT FITTING FUNCTION FOR THE DATA #####################################################################
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&##


count_data_fitting_RW<-function(simulated_dataset,stan_data,params_used_for_sim){
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  
  stan_data_discrete<-stan_data        
  
  params_monitor_stan<-c("y_hat","iota","fitted_output","iota","beta","sigma_pen")
  
  mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/ART_mod_RW_constant_N.stan",
                              data = stan_data_discrete,
                              pars = params_monitor_stan,chains = 3,warmup = 500,iter = 1500,
                              control = list(adapt_delta = 0.85))
  
  test_summary_fitted<-rstan::summary(mod_hiv_prev,pars=c("fitted_output"))$summary
  test_summary_beta<-rstan::summary(mod_hiv_prev,pars=c("beta"))$summary
  
  test_summary_fitted<-data.frame(test_summary_fitted)
  test_summary_beta<-data.frame(test_summary_beta)
  
  r_hat_fitted_summary<-base::summary(test_summary_fitted$Rhat)
  r_hat_beta_summary<-base::summary(test_summary_beta$Rhat)
  
  
  plot_stan_model_fit<-function(model_output,sim_data){
    
    posts_hiv <- rstan::extract(model_output)
    
    xout<- seq(1970,2020,0.1)                         
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
    
    
    mod_time = seq(1970,2020,0.1)
    
    prev_median<-(apply(posts_hiv$fitted_output[,,7],2,median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,7],2,quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,6],2,median)
    incidence_low<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,6],2,quantile,probs=c(0.975))
    
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
    
    prev_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=prev_percent),colour="red",size=1.05)+
      geom_line(data = df_fit_prevalence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="prevalence %",title="fitted and true prevalence, fitting to count data of ART numbers")  
    
    inc_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=incidence),colour="red",size=1.05)+
      geom_line(data = df_fit_incidence,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = df_fit_incidence,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="incidence",title="fitted and true incidence, fitting to count data of ART numbers")  
    
    
    kappa_plot<-ggplot(data=sim_data)+geom_line(aes(x=time,y=kappa),colour="red",size=1.05)+
      geom_line(data = r_fit,aes(x=time,y=median),colour="midnightblue",size=1.02)+
      geom_ribbon(data = r_fit,aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.2)+
      labs(x="time",y="Kappa",title="fitted and true kappa, fitting to count data of ART numbers")  
    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist,plot_prev=prev_plot,
                plot_inc=inc_plot,plot_kappa=kappa_plot))
  }
  
  list_of_results<-plot_stan_model_fit(mod_hiv_prev,simulated_dataset)
  
  return(list(fitting_results=list_of_results,simulated_data = simulated_dataset,sim_params=params_used_for_sim,
              rstan_summary_fitted=test_summary_fitted,rstan_summary_beta=test_summary_beta,
              rhat_fitted=r_hat_fitted_summary,rhat_beta=r_hat_beta_summary))
}
