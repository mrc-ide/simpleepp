##################################################################################################################################
## Gettting the double peak simpleepp model running on the cluster ###############################################################
##################################################################################################################################

setwd("X:")
options(didehpc.username = "jd2117",didehpc.home = "X:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts_2"

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),
                            sources = "simpleepp/R/loop_functions_random_walk.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

##################################################################################################################################
## COnstant pop size model with ART ##############################################################################################
###################################################################################################################################
obj$task_status()
obj$task_times()
obj$cluster_load()
#

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

save(sim_model_output_changed_to_bell_curve,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/inital_simmed_data")

##############################################################################################################################
## Now we will get this epidemic fitting on the cluster ######################################################################
##############################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/true_epidemic_data_double_peak",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/sample_n_100_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/sample_n_500_data", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/sample_n_1000_data", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/double_peak_run/sample_n_5000_data", verbose = T)

plot(sampled_n_100_complete_data[sampled_n_100_complete_data$iteration==64,2])
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
  
  lines(sampled_n_100_complete_data[sampled_n_100_complete_data$iteration==i,2],col=colour)
  Sys.sleep(0.5)
  print(i)
}


params<-params_sim
params

sample_range<-1970:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
rows_to_evaluate<- 0:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                       name = "Double_peak_RW_1_100")

n_100_RW_first_order_loop$status()
n_100_RW_first_order_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_100_FIRST_11_49_MAY_29")


sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,
                                                                  params = params,
                                                                  simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                       name = "Double_RW_first_order_n_500")
n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_500_FIRST_11_50_MAY_29")

sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,
                                                                   params = params,
                                                                   simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                        name = "Double_RW_frist_1000")
n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_1000_FIRST_11_51_MAY_29")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,
                                                                   params = params,
                                                                   simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                        name = "Double_RW_first_5000")
n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_5000_FIRST_11_51_MAY_29")

##############################################################################################################################
## Now we will put on the second order spline data ###########################################################################
##############################################################################################################################

penalty_order <- 2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_second_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,
                                                                   params = params,
                                                                   simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                        name = "Double_RW_second_order_100")
n_100_RW_second_order_loop$status()
n_100_RW_second_order_loop_id<-n_100_RW_second_order_loop$id
save(n_100_RW_second_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_100_SECOND_11_52_MAY_29")

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_500_RW_second_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,
                                                                   params = params,
                                                                   simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                        name = "Double_RW_second_500")
n_500_RW_second_order_loop$status()
n_500_RW_second_order_loop_id<-n_500_RW_second_order_loop$id
save(n_500_RW_second_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_500_SECOND_11_53_MAY_29")

sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_1000_rw_second_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                    data_about_sampling = data_about_sampling,
                                                                    iteration_number = 100,
                                                                    params = params,
                                                                    simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                         name = "Double_RW_second_1000")
n_1000_rw_second_order_loop$status()
n_1000_rw_second_order_loop_id<-n_1000_rw_second_order_loop$id
save(n_1000_rw_second_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_1000_SECOND_11_53_MAY_29")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_5000_rw_second_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                    data_about_sampling = data_about_sampling,
                                                                    iteration_number = 100,
                                                                    params = params,
                                                                    simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                         name = "Double_RW_sec_5000")
n_5000_rw_second_order_loop$status()
n_5000_rw_second_order_loop_id<-n_5000_rw_second_order_loop$id
save(n_5000_rw_second_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/RW_5000_SECOND_11_54_MAY_29")

##############################################################################################################################
## Now we'll put the spline fitting on the cluster ###########################################################################
##############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                      name = "Double_Spline_first_100")
spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_100_FIRST_11_55_MAY_29")

sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,
                                                                        params = params,
                                                                        simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                      name = "Double_Spline_first_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_500_FIRST_11_56_MAY_29")

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,
                                                                         params = params,
                                                                         simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                       name = "Double_Spline_first_1000")
spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_1000_FIRST_11_57_MAY_29")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,
                                                                         params = params,
                                                                         simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                       name = "Double_Spline_first_5000")
spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_5000_FIRST_11_57_MAY_29")

##############################################################################################################################
## Now for the second order splines ##########################################################################################
##############################################################################################################################

penalty_order<-2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,
                                                                         params = params,
                                                                         simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                       name = "Double_Spline_second_order_n_100")
spline_second_order_n_100$status()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_100_SECOND_11_58_MAY_29")

sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                       data_about_sampling = data_about_sampling,
                                                                       iteration_number = 100,
                                                                       params = params,
                                                                       simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                     name = "Double_Spline_second_order_n_500")
spline_second_order_500$status()
spline_second_order_500_id<-spline_second_order_500$id
save(spline_second_order_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_500_SECOND_11_58_MAY_29")

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,
                                                                        params = params,
                                                                        simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                      name = "Double_Spline_second_order_n_1000")

spline_second_order_1000$status()
spline_second_order_1000_id<-spline_second_order_1000$id
save(spline_second_order_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_1000_SECOND_11_58_MAY_29")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,
                                                                        params = params,
                                                                        simulated_true_df = sim_model_output_changed_to_bell_curve$sim_df),
                                      name = "Double_Spline_second_order_n_5000")
spline_second_order_5000$status()
spline_second_order_5000_id<-spline_second_order_5000$id
save(spline_second_order_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/SPLINE_5000_SECOND_11_59_MAY_29")


###############################################################################################################################
## NOW WE WILL LOAD UP THE RESULTS FROM THE CLUSTER RUNS ######################################################################
###############################################################################################################################

##RW 1st and 2nd 

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/double_peak/cluster_ids/"
double_peak_ids <- list.files(path_name,full.names = T)
for(i in 1:length(double_peak_ids)){
  load(double_peak_ids[i],verbose = T)
}

rw_first_100<-obj$task_get(n_100_RW_first_order_id)
rw_first_500<-obj$task_get(n_500_RW_first_order_loop_id)
rw_first_1000<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_first_5000<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_sec_100<-obj$task_get(n_100_RW_second_order_loop_id)
rw_sec_500<-obj$task_get(n_500_RW_second_order_loop_id)
rw_sec_1000<-obj$task_get(n_1000_rw_second_order_loop_id)
rw_sec_5000<-obj$task_get(n_5000_rw_second_order_loop_id)


rw_first_100_res<-rw_first_100$result()
rw_first_500_res<-rw_first_500$result()
rw_first_1000_res<-rw_first_1000$result()
rw_first_5000_res<-rw_first_5000$result()

rw_sec_100_res<-rw_sec_100$result()
rw_sec_500_res<-rw_sec_500$result()
rw_sec_1000_res<-rw_sec_1000$result()
rw_sec_5000_res<-rw_sec_5000$result()

#### splino first and seccy


spline_first_100<-obj$task_get(spline_first_order_n_100_id)
spline_first_500<-obj$task_get(spline_first_order_n_500_id)
spline_first_1000<-obj$task_get(spline_first_order_n_1000_id)
spline_first_5000<-obj$task_get(spline_first_order_n_5000_id)

spline_sec_100<-obj$task_get(spline_second_order_n_100_id)
spline_sec_500<-obj$task_get(spline_second_order_500_id)
spline_sec_1000<-obj$task_get(spline_second_order_1000_id)
spline_sec_5000<-obj$task_get(spline_second_order_5000_id)


spline_first_100_res<-spline_first_100$result()
spline_first_500_res<-spline_first_500$result()
spline_first_1000_res<-spline_first_1000$result()
spline_first_5000_res<-spline_first_5000$result()

spline_sec_100_res<-spline_sec_100$result()
spline_sec_500_res<-spline_sec_500$result()
spline_sec_1000_res<-spline_sec_1000$result()
spline_sec_5000_res<-spline_sec_5000$result()
obj$cluster_load()


save(rw_first_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_100")
save(rw_first_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_500")
save(rw_first_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_1000")
save(rw_first_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_5000")

save(rw_sec_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_100")
save(rw_sec_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_500")
save(rw_sec_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_1000")
save(rw_sec_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_5000")

save(spline_first_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_100")
save(spline_first_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_500")
save(spline_first_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_1000")
save(spline_first_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_5000")

save(spline_sec_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_100")
save(spline_sec_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_500")
save(spline_sec_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_1000")
save(spline_sec_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_5000")




