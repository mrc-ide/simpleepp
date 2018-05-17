###################################################################################################################################
## Testing out how the new art diseased model work ################################################################################
###################################################################################################################################

require(ggplot2)
require(rstan)
require(reshape2)
require(splines)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=T)

expose_stan_functions("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model_random_walk.stan")

rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

step_vector<-seq(1970.1,by = 0.1,length.out = 500)
kappa_params<-c(log(0.5),log(0.1),0.3,1995)

kappa<- exp(rlogistic(step_vector,kappa_params))                       ## This is our transmission parameter for the output

iota_val<-0.0001                                                       ## This is our initial proportion infected

mu <- 1/35                                                             ## Population level death rate
            
sigma<- c(1/3.16,1/2.13,1/3.2)                                         ## Progression along Cd4 stages among untreated

mu_i <- c(0.003, 0.008, 0.035, 0.27)                                   ## Mortality by stage, no ART

mu_d <- c(0.003,0.008,0.035,0.27)                                      ## Mortality be stage on diagnosed

mu_a <- c(0.002, 0.006, 0.006, 0.03)                                   ## Mortality by stage, on ART

omega <- 0.1                                                          ## Reduction in transmissability on art

theta <- 0.2                                                           ## Reduction when know diagnosed

dt <- 0.1

start<- 1970

diag_start<- 1975

art_start<-1976

sigmoid_curve_function<-function(x,lp){                               ## This is our function to produce the change in 
  (1 / (1 + exp(-lp[1]*x))) * lp[2]                                   ## diag rates and art uptake rates through time 
}

zero_time_art<-rep(0,((art_start-start)*10))                          ## How long initially art is not available       

zero_time_diag<-rep(0,((diag_start-start)*10))                        ## How long initially diagnoses not available

art_length<-length(kappa)+1 - length(zero_time_art)                   ## how long art is available for 

diag_time<-length(kappa)+1 - length(zero_time_diag)                   ## how long diag is available for

elll<-c(0.05,0.7)                                                     ## first number represents the curve on the uptake,
ellm<-c(0.05,0.7)                                                     ## second number represents the peak fraction on ART per time step
elln<-c(0.05,0.8)
ello<-c(0.1,0.9)

art_col_1<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),l = elll)) ## These are our columns for the Alpha matirx
art_col_2<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),ellm))     ## 1 is cd4>500, 4 is cd4<200
art_col_3<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),elln))
art_col_4<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),ello))

diag_params<-c(0.1,0.9)
diag_params_aids<-c(0.5,0.9)

diag_col_1<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))
diag_col_2<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))   ## Proportion to be diagnosed at each infected cd4 stage
diag_col_3<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))
diag_col_4<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params_aids))

alpha<-cbind(art_col_1,art_col_2,art_col_3,art_col_4)
diag<-cbind(diag_col_1,diag_col_2,diag_col_3,diag_col_4)

art_prog<-c(1/10,1/10,1/10)

obem<-(60*24*365)/(60*10^6)*5

test_run<-simpleepp_art_diag(kappa = kappa,iota = iota,alpha = alpha,mu = mu,sigma = sigma, mu_i = mu_i, mu_d = mu_d,mu_a = mu_a,
                             omega = omega,theta = theta,dt = dt,start = start,diag_start = diag_start,art_start = art_start,
                             diag = diag,art_prog = art_prog,onem = obem)
sim_data_art<-data.frame(test_run)
sim_data_art$time<-seq(start,by = dt,length.out = nrow(alpha))
names(sim_data_art)[1:7]<-c("kappa","art","diag","N","ART_inc", "incidence","prevalence")
sim_data_art$prev_percent<-sim_data_art$prevalence * 100

sim_plot<-function(sim_df,col_ids_to_get_rid){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=incidence),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  melt_df<-sim_df[,-c(col_ids_to_get_rid)]
  melted_df_funky<-melt(melt_df,id="time")
  
  three_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot,incidence=incidence_plot,kappa=kappa_plot,whole=three_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_data_art,c(2:5,9))
plot(plotted_sim$whole)
plot(plotted_sim$incidence)
plot(plotted_sim$prevalence)

###############################################################################################################################
## Now we will perform our sampling from the population #######################################################################
###############################################################################################################################


sample_range<-1970:2015
sample_years<-46
sample_n<-10000


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



sample_df_prev<-sample_function(sample_range,sample_years,sample_n,
                                       simulated_df = sim_data_art,prevalence_column_id = 7)


plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = sim_data_art,sample_df = sample_df_prev)

diagnosed_cd4_sample<-function(year_range_to_sample_from,integer_number_of_years_we_sampled_from,simulated_data,
                              col_id_time_and_inc){
  data_output<-NULL
  
  for(i in year_range_to_sample_from[1]:year_range_to_sample_from[integer_number_of_years_we_sampled_from]){
    overall_year<- i
    tot_for_indep_year<-0
    for(j in 1:10){
      add_on<-(j-1)/10
      overall_year_dec<-overall_year + add_on
      art_inc<-simulated_data[simulated_data$time %in% overall_year_dec, col_id_time_and_inc$inc]
      tot_for_indep_year<-tot_for_indep_year + art_inc
    }
    
    year_data<-cbind(overall_year,tot_for_indep_year)
    data_output<-rbind.data.frame(data_output,year_data)
  }
  
  poisson_samples<-rpois(nrow(data_output),data_output$tot_for_indep_year)
  
  data_output$samples<-poisson_samples
  names(data_output)<-c("year","true","sampled")
  
  sample_plot<-ggplot(data = data_output)+geom_line(aes(x=year,y=true),colour="midnightblue")+
    geom_point(aes(x=year,y=sampled),colour="red")+labs(x="year",y="Number of ART starting with CD4 > 500",
                                                        title="Sampled numbers of new ART starters with CD4 > 500")
  
  return(list(sample_data=data_output,sample_plot=sample_plot))
  
}

year_seq<- art_start:2015
number_year_obs<-length(year_seq)
col_ids<-list(time=8,inc=5)

poisson_sampled_data<-diagnosed_cd4_sample(year_seq,number_year_obs,simulated_data = sim_data_art,col_id_time_and_inc = col_ids)

poisson_sampled_data$sample_plot

sneaky_sample<-(art_start-start):(2015-start)*10+1
sneaky_poiss<-rpois(length(sneaky_sample),sim_data_art[sneaky_sample,5])

###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################
penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-0:45*10+1

stan_data_discrete<-list(
  n_obs = sample_years,
  y = as.array(sample_df_prev$sample_y_hiv_prev),
  n_sample = sample_n,
  time_steps_euler = length(xout)+1,
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

prev_optim<-stan_model("hiv_project/simpleepp/stan_files/chunks/ART_diag_prev_fitting.stan")
prev_optim_results<-optimizing(prev_optim,stan_data_discrete,as_vector=F)

## Lets plot the optim results

optim_plotter<-function(simulated_df,optimed_df){
 data_from_optim<-data.frame(optimed_df$par$fitted_output)
 data_from_optim$time<-seq(1970,2020,0.1)
 names(data_from_optim)[1:7]<-c("kappa","art","diag","N","ART_inc", "incidence","prevalence")
 
 prev_plot<-ggplot(data = simulated_df)+geom_line(aes(x=time,y=prevalence),colour="red",size=1)+
   geom_line(data = data_from_optim, aes(x=time,y=prevalence),colour="midnightblue",size=1)+
   labs(x="time",y="prevalence as fraction",title="Prevalence prediction from optimizer in STAN, Red is TRUE data")
 inc_plot<-ggplot(data = simulated_df)+geom_line(aes(x=time,y=incidence),colour="red",size=1)+
   geom_line(data = data_from_optim,aes(x=time,y=incidence),colour="midnightblue",size=1)+
   labs(x="time",y="incidence",title="Incidence prediction from optimizer in Stan Red is TRUE data")
 kappa_plot<-ggplot(data = simulated_df)+geom_line(aes(x=time,y=kappa),colour="red",size=1)+
   geom_line(data = data_from_optim,aes(x=time,y=kappa),colour="midnightblue",size=1)+
   labs(x="time",y="kappa",title="Kappa prediction from optimizer in Stan Red is TRUE data")
 ART_inc_plot<-ggplot(data = simulated_df)+geom_line(aes(x=time,y=ART_inc),colour="red",size=1)+
   geom_line(data = data_from_optim,aes(x=time,y=ART_inc),colour="midnightblue",size=1)+
   labs(x="time",y="incidence of ART uptake in cd4>500",title="ART incidence with CD4 > 500, RED is true Data")
 
 prev_rmse<- sqrt(mean((simulated_df$prevalence - data_from_optim$prevalence)^2))
 inc_rmse<- sqrt(mean((simulated_df$incidence - data_from_optim$incidence)^2))
 kappa_rmse<- sqrt(mean((simulated_df$kappa - data_from_optim$kappa)^2))
 art_inc_rmse<- sqrt(mean((simulated_df$ART_inc - data_from_optim$ART_inc)^2))
 
 rmse_df<-cbind.data.frame(prev_rmse,inc_rmse,kappa_rmse,art_inc_rmse)
 
 return(list(prev_plot=prev_plot,inc_plot=inc_plot,kappa_plot=kappa_plot,art_plot=ART_inc_plot, rmse_df=rmse_df))
 
}

optim_results<-optim_plotter(sim_data_art,prev_optim_results)
optim_results$prev_plot
optim_results$inc_plot
optim_results$kappa_plot
optim_results$art_plot
optim_results$rmse_df
require(ggpubr)

total_plots_prev<-ggarrange(optim_results$prev_plot,
                            optim_results$inc_plot,
                            optim_results$kappa_plot,
                            optim_results$art_plot,
                            ncol = 2, nrow = 2)


test_stan_hiv<- stan("hiv_project/simpleepp/stan_files/chunks/ART_diag_prev_fitting.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("hiv_project/simpleepp/stan_files/chunks/ART_diag_prev_fitting.stan", data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.8))

###############################################################################################################################
## This is our poisson sampling ###############################################################################################
###############################################################################################################################

penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-seq((art_start-start)*10+1,(year_seq[length(year_seq)]-start+1)*10)
sneaky_rows_evaluate<-sneaky_sample

stan_data_discrete<-list(
  n_obs = number_year_obs,
  y = as.array(sneaky_poiss),
  time_steps_euler = length(xout)+1,
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(sneaky_rows_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

stan_mod_count<-stan_model("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model_random_walk_ART_count_fit.stan")
optim_poiss<-optimizing(stan_mod_count,stan_data_discrete,as_vector=F)
optim_results_poisson<-optim_plotter(sim_data_art,optim_poiss)

optim_results_poisson$prev_plot
optim_results_poisson$inc_plot
optim_results_poisson$kappa_plot
optim_results_poisson$art_plot
optim_results_poisson$rmse_df
require(ggpubr)
total_plots<-ggarrange(optim_results_poisson$prev_plot,
                       optim_results_poisson$inc_plot,
                       optim_results_poisson$kappa_plot,
                       optim_results_poisson$art_plot,
                       ncol = 2, nrow = 2)

test_stan_hiv<- stan("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model_random_walk_ART_count_fit.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model_random_walk_ART_count_fit.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.8))


###############################################################################################################################
## Now we'll fit it with a spline model of kappa '#############################################################################
###############################################################################################################################
xout<-seq(1970,2020,0.1)
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
  y = as.array(sample_df_prev$sample_y_hiv_prev),
  n_sample = sample_n,
  knot_number = knot_number,
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)
params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

stan_mod_spline_prev_fit<-stan_model("hiv_project/simpleepp/stan_files/chunks/ART_diag_spline_prev_fitting.stan")
optim_fit_spline_mod_prev<-optimizing(stan_mod_spline_prev_fit,stan_data_discrete,as_vector=F)
optim_results_spline_prev<-optim_plotter(sim_data_art,optim_fit_spline_mod_prev)

optim_results_spline_prev$prev_plot
optim_results_spline_prev$inc_plot
optim_results_spline_prev$kappa_plot
optim_results_spline_prev$art_plot


test_stan_hiv<- stan("hiv_project/simpleepp/stan_files/chunks/ART_diag_spline_prev_fitting.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("hiv_project/simpleepp/stan_files/chunks/ART_diag_spline_prev_fitting.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.85))


stan_data_discrete<-list(
  n_obs = length(sneaky_sample),
  y = as.array(sneaky_sample),
  n_sample = sample_n,
  knot_number = knot_number,
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(sneaky_rows_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

optim_count_spline<-stan_model("hiv_project/simpleepp/stan_files/chunks/art_diag_model_count_fit_spline.stan")
optim_fit_spline_count<-optimizing(optim_count_spline,stan_data_discrete,as_vector=F)
optim_spline_results_count<-optim_plotter(sim_data_art,optim_fit_spline_count)

optim_spline_results_count$prev_plot
optim_spline_results_count$inc_plot
optim_spline_results_count$kappa_plot
optim_spline_results_count$art_plot 


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
  
  
  plotter<-ggplot(sim_sample) +
    geom_point(aes(x=sample_time_hiv, y=sample_prev_hiv_percentage),col="red", shape = 19, size = 1.5) +
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
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist))
  
  
}













