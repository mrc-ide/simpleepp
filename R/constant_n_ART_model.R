##################################################################################################################################
## COnstant pop size model with ART ##############################################################################################
###################################################################################################################################

require(rstan)
require(ggplot2)
require(reshape2)
require(ggpubr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=T)

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_RW_constant_N.stan")

sim_art_epidemic<-function(kappa_params,params){
  
  rlogistic <- function(t, p) {
    p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
  }
  
  step_vector<-seq(1970.1,by = 0.1,length.out = 500)
  kappa_params<-c(log(kappa_params$peak),log(kappa_params$low),kappa_params$grad,kappa_params$decline)
  
  kappa<- exp(rlogistic(step_vector,kappa_params))                       ## This is our transmission parameter for the output
  
  iota_val<-params$iota                                                       ## This is our initial proportion infected
  
  mu <- params$mu                                                             ## Population level death rate
  
  sigma<- params$sigma                                         ## Progression along Cd4 stages among untreated
  
  mu_i <- params$mu_i                                   ## Mortality by stage, no ART
  
  mu_d <- params$mu_d                                      ## Mortality be stage on diagnosed
  
  mu_a <- params$mu_a                                   ## Mortality by stage, on ART
  
  omega <- params$omega                                                          ## Reduction in transmissability on art
  
  theta <- params$theta                                                           ## Reduction when know diagnosed
  
  dt <- params$dt
  
  start<- params$start
  
  diag_start<- params$diag_start
  
  art_start<-params$art_start
  
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
  
  art_prog<-params$art_prog
  
  obem<-params$obem
  
  test_run<-simpleepp_art_diag_constant_N(kappa = kappa,iota = iota,alpha = alpha,mu = mu,sigma = sigma, mu_i = mu_i, mu_d = mu_d,mu_a = mu_a,
                               omega = omega,theta = theta,dt = dt,start = start,diag_start = diag_start,art_start = art_start,
                               diag = diag,art_prog = art_prog,onem = obem)
  sim_data_art_early_72_ART<-data.frame(test_run)
  sim_data_art_early_72_ART$time<-seq(start,by = dt,length.out = nrow(alpha))
  names(sim_data_art_early_72_ART)[1:7]<-c("kappa","art","diag","N","ART_inc", "incidence","prevalence")
  sim_data_art_early_72_ART$prev_percent<-sim_data_art_early_72_ART$prevalence * 100
  
  params$kappa<-kappa_params
  
  return(list(df=sim_data_art_early_72_ART,params_used=params,alpha_out=alpha,diag_out=diag))
  
}

kappa_params<-list(peak=0.5,low=0.1,grad=0.3,decline=1995)
iota<-0.0001
mu <- 1/35                                                             ## Population level death rate
sigma<- c(1/3.16,1/2.13,1/3.2)                                         ## Progression along Cd4 stages among untreated
mu_i <- c(0.003, 0.008, 0.035, 0.27)                                   ## Mortality by stage, no ART
mu_d <- c(0.003,0.008,0.035,0.27)                                      ## Mortality be stage on diagnosed
mu_a <- c(0.002, 0.006, 0.006, 0.03)                                   ## Mortality by stage, on ART
omega <- 0.8                                                          ## Reduction in transmissability on art
theta <- 0.2                                                           ## Reduction when know diagnosed
dt <- 0.1                                                              ## Time step difference
start<- 1970                                                           ## start of the epidemic
diag_start<- 1971                                                      ## start of diagnosis
art_start<-1972                                                       ## start of ART
art_prog<-c(1/10,1/10,1/10)                                            ## Progression rates through CD4 stages on ART
obem<-(60*24*365)/(60*10^6)*5                                          ## Birth rate in the population

params<-list(iota=iota,mu=mu,sigma=sigma,mu_i=mu_i,mu_d=mu_d,mu_a=mu_a,
             omega=omega,theta=theta,dt=dt,start=start,diag_start=diag_start,
             art_start=art_start,art_prog=art_prog,obem=obem)

art_epidemic<-sim_art_epidemic(kappa_params = kappa_params,params = params)


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

plotted_sim<-sim_plot(art_epidemic$df,c(2:5,9))
plot(plotted_sim$whole)

plot(art_epidemic$df$N)
plot(plotted_sim$incidence)
plot(plotted_sim$prevalence)

###############################################################################################################################
## Lets do some sampling !""!"!"!"!""!"!"?!>"?!>"@!":!?"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
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
                                simulated_df = art_epidemic$df,prevalence_column_id = 7)


plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = art_epidemic$df,sample_df = sample_df_prev)

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

poisson_sampled_data<-diagnosed_cd4_sample(year_seq,number_year_obs,
                                           simulated_data = art_epidemic$df,
                                           col_id_time_and_inc = col_ids)

poisson_sampled_data$sample_plot


###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################
penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-0:45*10+1

stan_data_discrete_prev_rw<-list(
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
  alpha = art_epidemic$alpha_out,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = art_epidemic$diag_out,
  art_prog = art_prog,
  onem = obem
  
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

prev_optim<-stan_model("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_RW_Prevalence_constant_N.stan")
prev_optim_results<-optimizing(prev_optim,stan_data_discrete_prev_rw,as_vector=F)

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

optim_results<-optim_plotter(art_epidemic$df,prev_optim_results)
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


test_stan_hiv<- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_RW_Prevalence_constant_N.stan",
                     data = stan_data_discrete_prev_rw,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

###############################################################################################################################
## NOw testing for the count datt fitting #####################################################################################
###############################################################################################################################

penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-seq((art_start-start)*10+1,(year_seq[length(year_seq)]-start+1)*10)

stan_data_discrete_count_rw<-list(
  n_obs = nrow(poisson_sampled_data$sample_data),
  y = as.array(poisson_sampled_data$sample_data$sampled),
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
  alpha = art_epidemic$alpha_out,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = art_epidemic$diag_out,
  art_prog = art_prog,
  onem = obem
  
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

stan_mod_count<-stan_model("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_RW_constant_N.stan")
optim_poiss<-optimizing(stan_mod_count,stan_data_discrete_count_rw,as_vector=F)
optim_results_poisson<-optim_plotter(art_epidemic$df,optim_poiss)

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

test_stan_hiv<- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_RW_constant_N.stan",
                     data = stan_data_discrete_count_rw,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
## Now we will fit with the splines for kappa parameter #######################################################################
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

stan_data_discrete_prev_spline<-list(
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
  alpha = art_epidemic$alpha_out,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = art_epidemic$diag_out,
  art_prog = art_prog,
  onem = obem
  
)
params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

stan_mod_spline_prev_fit<-stan_model("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_spline_prevalence_constant_N.stan")
optim_fit_spline_mod_prev<-optimizing(stan_mod_spline_prev_fit,stan_data_discrete_prev_spline,as_vector=F)
optim_results_spline_prev<-optim_plotter(art_epidemic$df,optim_fit_spline_mod_prev)

optim_results_spline_prev$prev_plot
optim_results_spline_prev$inc_plot
optim_results_spline_prev$kappa_plot
optim_results_spline_prev$art_plot

total_spline_prev<-ggarrange(optim_results_spline_prev$prev_plot,
                             optim_results_spline_prev$inc_plot,
                             optim_results_spline_prev$kappa_plot,
                             optim_results_spline_prev$art_plot,
                             ncol = 2, nrow = 2)

test_stan_hiv<- stan("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_spline_prevalence_constant_N.stan",
                     data = stan_data_discrete_prev_spline,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


#mod_hiv_prev <- stan("hiv_project/simpleepp/stan_files/chunks/ART_diag_spline_prev_fitting.stan",
#data = stan_data_discrete_prev_spline,
#pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
#control = list(adapt_delta = 0.85))


###############################################################################################################################
## Now we will fit with the count data ########################################################################################
###############################################################################################################################


rows_to_evaluate<-seq((art_start-start)*10+1,(year_seq[length(year_seq)]-start+1)*10)


stan_data_discrete_count_spline<-list(
  n_obs = nrow(poisson_sampled_data$sample_data),
  y = as.array(poisson_sampled_data$sample_data$sampled),
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
  alpha = art_epidemic$alpha_out,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = art_epidemic$diag_out,
  art_prog = art_prog,
  onem = obem
  
)

optim_count_spline<-stan_model("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/ART_mod_spline_COUNT_constant_N.stan")
optim_fit_spline_count<-optimizing(optim_count_spline,stan_data_discrete_count_spline,as_vector=F)
optim_spline_results_count<-optim_plotter(art_epidemic$df,optim_fit_spline_count)

optim_spline_results_count$prev_plot
optim_spline_results_count$inc_plot
optim_spline_results_count$kappa_plot
optim_spline_results_count$art_plot 

spline_count_total_plots<-ggarrange(optim_spline_results_count$prev_plot,
                                    optim_spline_results_count$inc_plot,
                                    optim_spline_results_count$kappa_plot,
                                    optim_spline_results_count$art_plot,
                                    ncol = 2,nrow = 2)


test_stan_hiv<- stan("hiv_project/simpleepp/stan_files/chunks/art_diag_model_count_fit_spline.stan",
                     data = stan_data_discrete_count_spline,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

 

