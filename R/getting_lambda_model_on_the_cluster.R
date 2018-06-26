###########################################################################################################################
## Getting the FOI as modelled parameter equations running on the cluster #################################################
###########################################################################################################################
Sys.sleep(5)
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
                            sources = "simpleepp/R/lambda_cluster_functions.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

########################################################################################################################
## Now we will quicky test our access to the cluster ###################################################################
########################################################################################################################

t<-obj$enqueue(sessionInfo())
t$wait(10)

obj$task_status()

########################################################################################################################
## Now we will load up the 100 iterations of samples we have created ###################################################
########################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_100_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_500_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_1000_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_5000_sampled_data",verbose = T)

plot(sampled_n_5000_complete_data_foi_model[sampled_n_5000_complete_data_foi_model$iteration==64,2])
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
  
  lines(sampled_n_5000_complete_data_foi_model[sampled_n_5000_complete_data_foi_model$iteration==i,2],col=colour)
  Sys.sleep(0.3)
  print(i)
}


########################################################################################################################
## So now we will set up the first order rw model simulations for the lambda parameter #################################
########################################################################################################################

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)    
foi_flag<-as.integer(1)

knot_number<-51
penalty_order<-1
model_type<-"rw"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

rw_first_100_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "RW_first_100_foi")
rw_first_100_foi$wait(10)
rw_first_100_foi$status()
rw_first_100_foi_id<-rw_first_100_foi$id
save(rw_first_100_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_FIRST_100_14_46_MAY_29")

## 500

sample_n<-500
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_first_500_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "RW_first_500_foi")
rw_first_500_foi$status()
rw_first_500_foi_id<-rw_first_500_foi$id
save(rw_first_500_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_FRIST_500_14_49_MAY_29")

## 1K
sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_first_1000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "RW_first_1000_foi")
rw_first_1000_foi$status()
rw_first_1000_foi_id<-rw_first_1000_foi$id
save(rw_first_1000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_FRIST_1000_14_51_MAY_29")

## 5K

sample_n<-5000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_first_5000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "RW_first_5000_foi")
rw_first_5000_foi$status()
rw_first_5000_foi_id<-rw_first_5000_foi$id
save(rw_first_5000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_FRIST_5000_11_26_JUNE_4")
#######################################################################################################################
### now we will do the second order rw runs ###########################################################################
#######################################################################################################################

penalty_order<-2
sample_n<-100
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_first_100_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                              iteration_number = 100,
                                                              data_about_sampling = data_about_sampling,
                                                              params = params_cluster), name = "RW_sec_100_foi")
rw_first_100_foi$status()
rw_sec_100_foi_id<-rw_first_100_foi$id
save(rw_sec_100_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_100_SECOND_14_58_MAY_29")

##### 500 #####


sample_n<-500
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_sec_500_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "RW_sec_500_foi")
rw_sec_500_foi$status()
rw_sec_500_foi_id<-rw_sec_500_foi$id
save(rw_sec_500_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_500_SECOND_15_00_MAY_29")


##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_sec_1000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                           iteration_number = 100,
                                                           data_about_sampling = data_about_sampling,
                                                           params = params_cluster), name = "RW_sec_1000_foi")
rw_sec_1000_foi$status()
rw_sec_1000_foi_id<-rw_sec_1000_foi$id
save(rw_sec_1000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_500_SECOND_15_02_MAY_29")

##### 5k #####

sample_n<-5000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)
rw_sec_5000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                           iteration_number = 100,
                                                           data_about_sampling = data_about_sampling,
                                                           params = params_cluster), name = "RW_sec_5000_foi")
rw_sec_5000_foi$status()
rw_sec_5000_foi_id<-rw_sec_5000_foi$id
save(rw_sec_5000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/RW_5000_SECOND_11_28_JUNE_4")

######################################################################################################################
## NOw for the 7 knot splines ########################################################################################
######################################################################################################################
knot_number<-7
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_7_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                        iteration_number = 100,
                                                                        data_about_sampling = data_about_sampling,
                                                                        params = params_cluster), name = "SPLINE_first_100_foi_7_knot")
spline_first_100_foi_7_knot$status()
spline_first_100_foi_7_knot_id<-spline_first_100_foi_7_knot$id
save(spline_first_100_foi_7_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_first_100_7_KNOT_11_51_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_first_500_foi_7_knot")
spline_first_500_foi_7$status()
spline_first_500_foi_id_7<-spline_first_500_foi_7$id
save(spline_first_500_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_first_500_7_knot_11_52_JUNE_4")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_1000_foi_7_knot")
spline_first_1000_foi_id_7<-spline_first_1000_foi_7$id
save(spline_first_1000_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_first_1000_7_knot_11_54_JUNE_4")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_5000_foi_7_knot")
spline_first_5000_foi_7$status()
spline_first_5000_foi_id_7<-spline_first_5000_foi_7$id
save(spline_first_5000_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_first_5000_7_knot_11_55_JUNE_4")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_100_foi_7_knot")
spline_sec_100_foi_7$status()
spline_sec_100_foi_id_7<-spline_sec_100_foi_7$id
save(spline_sec_100_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_sec_100_7_knot_11_57_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_500_foi_7_knot")
spline_sec_500_foi_7$status()
spline_sec_500_foi_id_7<-spline_sec_500_foi_7$id
save(spline_sec_500_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_sec_500_7_knot_12_01_JUNE_4")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_1000_foi_7_knot")
spline_sec_1000_foi_7$status()
spline_sec_1000_foi_id_7<-spline_sec_1000_foi_7$id
save(spline_sec_1000_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_sec_1000_7_knot_12_03_JUNE_4")

##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_7<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_5000_foi_7_knot")
spline_sec_5000_foi_7$status()
spline_sec_5000_foi_id_7<-spline_sec_5000_foi_7$id
save(spline_sec_5000_foi_id_7,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/SPLINE_sec_5000_7_knot_12_04_JUNE_4")

#################################################################################
## Lets load up the results #####################################################
#################################################################################

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/RW/"
rw<-list.files(path_name,full.names = T)
for(i in 1:length(rw)){
  load(rw[i],verbose = T)
}

rw_f_foi_100<-obj$task_get(rw_first_100_foi_id)
rw_f_foi_500<-obj$task_get(rw_first_500_foi_id)
rw_f_foi_1000<-obj$task_get(rw_first_1000_foi_id)
rw_f_foi_5000<-obj$task_get(rw_first_5000_foi_id)

rw_f_foi_100_res <- rw_f_foi_100$result()
rw_f_foi_500_res <- rw_f_foi_500$result()
rw_f_foi_1000_res <- rw_f_foi_1000$result()
rw_f_foi_5000_res <- rw_f_foi_5000$result()

save(rw_f_foi_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_1st_100_foi_result")
save(rw_f_foi_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_1st_500_foi_result")
save(rw_f_foi_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_1st_1000_foi_result")
save(rw_f_foi_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_1st_5000_foi_result")

rw_s_foi_100<-obj$task_get(rw_sec_100_foi_id)
rw_s_foi_500<-obj$task_get(rw_sec_500_foi_id)
rw_s_foi_1000<-obj$task_get(rw_sec_1000_foi_id)
rw_s_foi_5000<-obj$task_get(rw_sec_5000_foi_id)

rw_s_foi_100_res <- rw_s_foi_100$result()
rw_s_foi_500_res <- rw_s_foi_500$result()
rw_s_foi_1000_res <- rw_s_foi_1000$result()
rw_s_foi_5000_res <- rw_s_foi_5000$result()

save(rw_s_foi_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_2nd_100_foi_result")
save(rw_s_foi_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_2nd_500_foi_result")
save(rw_s_foi_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_2nd_1000_foi_result")
save(rw_s_foi_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/RW/RW_2nd_5000_foi_result")


#######################################################################################################################
## Now for the splines ################################################################################################
#######################################################################################################################

knot_number<-8
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_8_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "SPLINE_first_100_foi_8_knot")
spline_first_100_foi_8_knot$status()
spline_first_100_foi_8_knot_id<-spline_first_100_foi_8_knot$id
save(spline_first_100_foi_8_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_100_8_KNOT_11_51_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_first_500_foi_8_knot")
spline_first_500_foi_8$status()
spline_first_500_foi_id_8<-spline_first_500_foi_8$id
save(spline_first_500_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_500_8_knot_11_52_JUNE_4")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_first_1000_foi_8_knot")
spline_first_1000_foi_8$status()
spline_first_1000_foi_id_8<-spline_first_1000_foi_8$id
save(spline_first_1000_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_1000_8_knot_11_54_JUNE_4")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                params = params_cluster), name = "SPLINE_first_5000_foi_8_knot")
spline_first_5000_foi_8$status()
spline_first_5000_foi_id_8<-spline_first_5000_foi_8$id
save(spline_first_5000_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_5000_8_knot_11_55_JUNE_4")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_100_foi_8_knot")
spline_sec_100_foi_8$status()
spline_sec_100_foi_id_8<-spline_sec_100_foi_8$id
save(spline_sec_100_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_100_8_knot_11_57_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_500_foi_8_knot")
spline_sec_500_foi_8$status()
spline_sec_500_foi_id_8<-spline_sec_500_foi_8$id
save(spline_sec_500_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_500_8_knot_12_01_JUNE_4")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_1000_foi_8_knot")
spline_sec_1000_foi_8$status()
spline_sec_1000_foi_id_8<-spline_sec_1000_foi_8$id
save(spline_sec_1000_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_1000_8_knot_12_03_JUNE_4")

##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_8<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_5000_foi_8_knot")
spline_sec_5000_foi_8$status()
spline_sec_5000_foi_id_8<-spline_sec_5000_foi_8$id
save(spline_sec_5000_foi_id_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_5000_8_knot_12_04_JUNE_4")


######################################################################################################################
## Now we will load up the results ###################################################################################
#####################################################################################################################

##### spline first #####
path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/7_knots/"
sev_knots <- list.files(path_name,full.names = T)
for(i in 1:length(sev_knots)){
  load(sev_knots[i],verbose = T)
}

spline_f_100_foi<-obj$task_get(spline_first_100_foi_7_knot_id)
spline_f_500_foi<-obj$task_get(spline_first_500_foi_id_7)
spline_f_1k_foi<-obj$task_get(spline_first_1000_foi_id_7)
spline_f_5k_foi<-obj$task_get(spline_first_5000_foi_id_7)

spline_f_100_foi_res<-spline_f_100_foi$result()
spline_f_500_foi_res<-spline_f_500_foi$result()
spline_f_1k_foi_res<-spline_f_1k_foi$result()
spline_f_5k_foi_res<-spline_f_5k_foi$result()

save(spline_f_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_1st_7_100_foi")
save(spline_f_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_1st_7_500_foi")
save(spline_f_1k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_1st_7_1000_foi")
save(spline_f_5k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_1st_7_5000_foi")


###### spline second order ######


spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_7)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_7)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_7)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_7)

spline_s_100_foi_res<-spline_s_100_foi$result()
spline_s_500_foi_res<-spline_s_500_foi$result()
spline_s_1k_foi_res<-spline_s_1k_foi$result()
spline_s_5k_foi_res<-spline_s_5k_foi$result()

save(spline_s_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_2nd_7_100_foi")
save(spline_s_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_2nd_7_500_foi")
save(spline_s_1000_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_2nd_7_1000_foi")
save(spline_s_5000_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/7_knots/sp_2nd_7_5000_foi")


##########################################################################################################################################
## NOw for the 8 knotters, wicked wicked junglist massive ################################################################################
##########################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_100_8_KNOT_11_51_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_500_8_knot_11_52_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_1000_8_knot_11_54_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_first_5000_8_knot_11_55_JUNE_4",
     verbose = T)

spline_f_100_foi_8<-obj$task_get(spline_first_100_foi_8_knot_id)
spline_f_500_foi_8<-obj$task_get(spline_first_500_foi_id_8)
spline_f_1k_foi_8<-obj$task_get(spline_first_1000_foi_id_8)
spline_f_5k_foi_8<-obj$task_get(spline_first_5000_foi_id_8)

spline_f_100_foi_res_8<-spline_f_100_foi_8$result()
spline_f_500_foi_res_8<-spline_f_500_foi_8$result()
spline_f_1k_foi_res_8<-spline_f_1k_foi_8$result()
spline_f_5k_foi_res_8<-spline_f_5k_foi_8$result()

save(spline_f_100_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/spline_100_1st_8_knot")
save(spline_f_500_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/spline_500_1st_8_knot")
save(spline_f_1k_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/spline_1000_1st_8_knot")
save(spline_f_5k_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/spline_5000_1st_8_knot")


###### spline second order ######

load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_100_8_knot_11_57_JUNE_4")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_500_8_knot_12_01_JUNE_4")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_1000_8_knot_12_03_JUNE_4")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/8_knots/SPLINE_sec_5000_8_knot_12_04_JUNE_4")

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_8)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_8)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_8)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_8)

spline_s_100_foi_res_8<-spline_s_100_foi$result()
spline_s_500_foi_res_8<-spline_s_500_foi$result()
spline_s_1k_foi_res_8<-spline_s_1k_foi$result()
spline_s_5k_foi_res_8<-spline_s_5k_foi$result()

save(spline_s_100_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/sp_sec_100_foi_8_knot")
save(spline_s_500_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/sp_sec_500_foi_8_knot")
save(spline_s_1k_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/sp_sec_1000_foi_8_knot")
save(spline_s_5k_foi_res_8,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/8_knots/sp_sec_5000_foi_8_knot")

#######################################################################################################################################
## Now for the 9 knotters, yo yo ######################################################################################################
#######################################################################################################################################

knot_number<-9
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_9_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                        iteration_number = 100,
                                                                        data_about_sampling = data_about_sampling,
                                                                        params = params_cluster), name = "SPLINE_first_100_foi_9_knot")
spline_first_100_foi_9_knot$status()
spline_first_100_foi_9_knot_id<-spline_first_100_foi_9_knot$id
save(spline_first_100_foi_9_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_first_100_9_KNOT_13_59_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_first_500_foi_9_knot")
spline_first_500_foi_9$status()
spline_first_500_foi_id_9<-spline_first_500_foi_9$id
save(spline_first_500_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_first_500_9_knot_14_00_JUNE_4")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_1000_foi_9_knot")
spline_first_1000_foi_9$status()
spline_first_1000_foi_id_9<-spline_first_1000_foi_9$id
save(spline_first_1000_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_first_1000_9_knot_14_00_JUNE_4")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_5000_foi_9_knot")
spline_first_5000_foi_9$status()
spline_first_5000_foi_id_9<-spline_first_5000_foi_9$id
save(spline_first_5000_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_first_5000_9_knot_14_01_JUNE_4")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_100_foi_9_knot")
spline_sec_100_foi_9$status()
spline_sec_100_foi_id_9<-spline_sec_100_foi_9$id
save(spline_sec_100_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_sec_100_9_knot_14_01_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_500_foi_9_knot")
spline_sec_500_foi_9$status()
spline_sec_500_foi_id_9<-spline_sec_500_foi_9$id
save(spline_sec_500_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_sec_500_9_knot_14_03_JUNE_4")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_1000_foi_9_knot")
spline_sec_1000_foi_9$status()
spline_sec_1000_foi_id_9<-spline_sec_1000_foi_9$id
save(spline_sec_1000_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_sec_1000_9_knot_14_04_JUNE_4")
##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_9<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_5000_foi_9_knot")
spline_sec_5000_foi_9$status()
spline_sec_5000_foi_id_9<-spline_sec_5000_foi_9$id
save(spline_sec_5000_foi_id_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/SPLINE_sec_5000_9_knot_14_04_JUNE_4")
#########################################################################################################################################
## Now to load up the results of the cluster runs #######################################################################################
#########################################################################################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/9_knots/"
nine_knotters <- list.files(path_name, full.names = T)
for(i in 1:length(nine_knotters)){
  load(nine_knotters[i], verbose = T)
}

spline_f_100_foi_9<-obj$task_get(spline_first_100_foi_9_knot_id)
spline_f_500_foi_9<-obj$task_get(spline_first_500_foi_id_9)
spline_f_1k_foi_9<-obj$task_get(spline_first_1000_foi_id_9)
spline_f_5k_foi_9<-obj$task_get(spline_first_5000_foi_id_9)

spline_f_100_foi_res_9<-spline_f_100_foi_9$result()
spline_f_500_foi_res_9<-spline_f_500_foi_9$result()
spline_f_1k_foi_res_9<-spline_f_1k_foi_9$result()
spline_f_5k_foi_res_9<-spline_f_5k_foi_9$result()

save(spline_f_100_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/spline_100_1st_9_knot")
save(spline_f_500_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/spline_500_1st_9_knot")
save(spline_f_1k_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/spline_1000_1st_9_knot")
save(spline_f_5k_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/spline_5000_1st_9_knot")


###### spline second order ######

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_9)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_9)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_9)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_9)

spline_s_100_foi_res_9<-spline_s_100_foi$result()
spline_s_500_foi_res_9<-spline_s_500_foi$result()
spline_s_1k_foi_res_9<-spline_s_1k_foi$result()
spline_s_5k_foi_res_9<-spline_s_5k_foi$result()

save(spline_s_100_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/sp_sec_100_foi_9_knot")
save(spline_s_500_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/sp_sec_500_foi_9_knot")
save(spline_s_1k_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/sp_sec_1000_foi_9_knot")
save(spline_s_5k_foi_res_9,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/9_knots/sp_sec_5000_foi_9_knot")


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$£££££££££
## Now for the 10 knotters ###########################################################################################################
#####())())()()()()()()()()()()()(!)(")(£)($)(%)(^)(&)(*)(1)(2)(3)(4)(5)(6)(7)(8)(9)()()()()()()(((((()))))))()()()()()()()()()()()###


knot_number<-10
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_ten_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                        iteration_number = 100,
                                                                        data_about_sampling = data_about_sampling,
                                                                        params = params_cluster), name = "SPLINE_first_100_foi_ten_knot")
spline_first_100_foi_ten_knot$status()
spline_first_100_foi_ten_knot_id<-spline_first_100_foi_ten_knot$id
save(spline_first_100_foi_ten_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_first_100_10_knots_14_16_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_first_500_foi_ten_knot")
spline_first_500_foi_ten$status()
spline_first_500_foi_id_ten<-spline_first_500_foi_ten$id
save(spline_first_500_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_first_500_10_knots_14_16_JUNE_4")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_1000_foi_ten_knot")
spline_first_1000_foi_ten$status()
spline_first_1000_foi_id_ten<-spline_first_1000_foi_ten$id
save(spline_first_1000_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_first_1000_10_knots_14_00_JUNE_4")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_5000_foi_ten_knot")
spline_first_5000_foi_ten$status()
spline_first_5000_foi_id_ten<-spline_first_5000_foi_ten$id
save(spline_first_5000_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_first_5000_10_knots_14_17_JUNE_4")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_100_foi_ten_knot")
spline_sec_100_foi_ten$status()
spline_sec_100_foi_id_ten<-spline_sec_100_foi_ten$id
save(spline_sec_100_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_sec_100_10_knots_14_17_JUNE_4")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_500_foi_ten_knot")
spline_sec_500_foi_ten$status()
spline_sec_500_foi_id_ten<-spline_sec_500_foi_ten$id
save(spline_sec_500_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_sec_500_10_knots_14_17_JUNE_4")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_1000_foi_ten_knot")
spline_sec_1000_foi_ten$status()
spline_sec_1000_foi_id_ten<-spline_sec_1000_foi_ten$id
save(spline_sec_1000_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_sec_1000_10_knots_14_18_JUNE_4")
##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_ten<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_5000_foi_ten_knot")
spline_sec_5000_foi_ten$status()
spline_sec_5000_foi_id_ten<-spline_sec_5000_foi_ten$id
save(spline_sec_5000_foi_id_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/SPLINE_sec_5000_10_knots_14_18_JUNE_4")
#########################################################################################################################################
## Now to load up the results of the cluster runs #######################################################################################
#########################################################################################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/10_knots/"
ten_knotters <- list.files(path_name, full.names = T)
for(i in 1:length(ten_knotters)){
  load(ten_knotters[i], verbose = T)
}

spline_f_100_foi_ten<-obj$task_get(spline_first_100_foi_ten_knot_id)
spline_f_500_foi_ten<-obj$task_get(spline_first_500_foi_id_ten)
spline_f_1k_foi_ten<-obj$task_get(spline_first_1000_foi_id_ten)
spline_f_5k_foi_ten<-obj$task_get(spline_first_5000_foi_id_ten)

spline_f_100_foi_res_ten<-spline_f_100_foi_ten$result()
spline_f_500_foi_res_ten<-spline_f_500_foi_ten$result()
spline_f_1k_foi_res_ten<-spline_f_1k_foi_ten$result()
spline_f_5k_foi_res_ten<-spline_f_5k_foi_ten$result()

save(spline_f_100_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/spline_100_1st_10_knots")
save(spline_f_500_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/spline_500_1st_10_knots")
save(spline_f_1k_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/spline_1000_1st_10_knots")
save(spline_f_5k_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/spline_5000_1st_10_knots")


###### spline second order ######

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_ten)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_ten)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_ten)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_ten)

spline_s_100_foi_res_ten<-spline_s_100_foi$result()
spline_s_500_foi_res_ten<-spline_s_500_foi$result()
spline_s_1k_foi_res_ten<-spline_s_1k_foi$result()
spline_s_5k_foi_res_ten<-spline_s_5k_foi$result()

save(spline_s_100_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/sp_sec_100_foi_10_knots")
save(spline_s_500_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/sp_sec_500_foi_10_knots")
save(spline_s_1k_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/sp_sec_1000_foi_10_knots")
save(spline_s_5k_foi_res_ten,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/10_knots/sp_sec_5000_foi_10_knots")

#######################################################################################################################################
## Now for the 11 knot splines ###############################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
##############################################################***********************************######################################

knot_number<-11
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_11_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                          iteration_number = 100,
                                                                          data_about_sampling = data_about_sampling,
                                                                          params = params_cluster), name = "SPLINE_first_100_foi_11_knot")
spline_first_100_foi_11_knot$status()
spline_first_100_foi_11_knot_id<-spline_first_100_foi_11_knot$id
save(spline_first_100_foi_11_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_first_100_11_knots_13_54_JUNE_5")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                     iteration_number = 100,
                                                                     data_about_sampling = data_about_sampling,
                                                                     params = params_cluster), name = "SPLINE_first_500_foi_11_knot")
spline_first_500_foi_11$status()
spline_first_500_foi_id_11<-spline_first_500_foi_11$id
save(spline_first_500_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_first_500_11_knots_13_55_JUNE_5")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                      iteration_number = 100,
                                                                      data_about_sampling = data_about_sampling,
                                                                      params = params_cluster), name = "SPLINE_first_1000_foi_11_knot")
spline_first_1000_foi_11$status()
spline_first_1000_foi_id_11<-spline_first_1000_foi_11$id
save(spline_first_1000_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_first_1000_11_knots_13_55_JUNE_5")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                      iteration_number = 100,
                                                                      data_about_sampling = data_about_sampling,
                                                                      params = params_cluster), name = "SPLINE_first_5000_foi_11_knot")
spline_first_5000_foi_11$status()
spline_first_5000_foi_id_11<-spline_first_5000_foi_11$id
save(spline_first_5000_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_first_5000_11_knots_13_56_JUNE_5")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_sec_100_foi_11_knot")
spline_sec_100_foi_11$status()
spline_sec_100_foi_id_11<-spline_sec_100_foi_11$id
save(spline_sec_100_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_sec_100_11_knots_13_56_JUNE_5")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_sec_500_foi_11_knot")
spline_sec_500_foi_11$status()
spline_sec_500_foi_id_11<-spline_sec_500_foi_11$id
save(spline_sec_500_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_sec_500_11_knots_13_56_JUNE_5")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_sec_1000_foi_11_knot")
spline_sec_1000_foi_11$status()
spline_sec_1000_foi_id_11<-spline_sec_1000_foi_11$id
save(spline_sec_1000_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_sec_1000_11_knots_13_57_JUNE_5")
##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_11<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_sec_5000_foi_11_knot")
spline_sec_5000_foi_11$status()
spline_sec_5000_foi_id_11<-spline_sec_5000_foi_11$id
save(spline_sec_5000_foi_id_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/SPLINE_sec_5000_11_knots_13_57_JUNE_5")
#########################################################################################################################################
## Now to load up the results of the cluster runs 11 11 11 11 11 11  ####################################################################
#########################################################################################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/11_knots/"
elev_knotters <- list.files(path_name,full.names = T)
for(i in 1:length(elev_knotters)){
  load(elev_knotters[i],verbose = T)
}

spline_f_100_foi_11<-obj$task_get(spline_first_100_foi_11_knot_id)
spline_f_500_foi_11<-obj$task_get(spline_first_500_foi_id_11)
spline_f_1k_foi_11<-obj$task_get(spline_first_1000_foi_id_11)
spline_f_5k_foi_11<-obj$task_get(spline_first_5000_foi_id_11)

spline_f_100_foi_res_11<-spline_f_100_foi_11$result()
spline_f_500_foi_res_11<-spline_f_500_foi_11$result()
spline_f_1k_foi_res_11<-spline_f_1k_foi_11$result()
spline_f_5k_foi_res_11<-spline_f_5k_foi_11$result()

save(spline_f_100_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/spline_100_1st_11_knots")
save(spline_f_500_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/spline_500_1st_11_knots")
save(spline_f_1k_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/spline_1000_1st_11_knots")
save(spline_f_5k_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/spline_5000_1st_11_knots")


###### spline second order ######

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_11)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_11)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_11)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_11)

spline_s_100_foi_res_11<-spline_s_100_foi$result()
spline_s_500_foi_res_11<-spline_s_500_foi$result()
spline_s_1k_foi_res_11<-spline_s_1k_foi$result()
spline_s_5k_foi_res_11<-spline_s_5k_foi$result()

save(spline_s_100_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/sp_sec_100_foi_11_knots")
save(spline_s_500_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/sp_sec_500_foi_11_knots")
save(spline_s_1k_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/sp_sec_1000_foi_11_knots")
save(spline_s_5k_foi_res_11,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/11_knots/sp_sec_5000_foi_11_knots")

########################################################################################################################################
##   12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12   ##
## Now we will run through the 12 knots#########12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12 #
##   12  12  12  12  12  12  12  12  12  12 12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12 #
## 12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  #
########################################################################################################################################
knot_number<-12
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi_12_knot<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                         iteration_number = 100,
                                                                         data_about_sampling = data_about_sampling,
                                                                         params = params_cluster), name = "SPLINE_first_100_foi_12_knot")
spline_first_100_foi_12_knot$status()
spline_first_100_foi_12_knot_id<-spline_first_100_foi_12_knot$id
save(spline_first_100_foi_12_knot_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_first_100_12_knots_14_06_JUNE_5")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                    iteration_number = 100,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params = params_cluster), name = "SPLINE_first_500_foi_12_knot")
spline_first_500_foi_12$status()
spline_first_500_foi_id_12<-spline_first_500_foi_12$id
save(spline_first_500_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_first_500_12_knots_14_07_JUNE_5")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                     iteration_number = 100,
                                                                     data_about_sampling = data_about_sampling,
                                                                     params = params_cluster), name = "SPLINE_first_1000_foi_12_knot")
spline_first_1000_foi_12$status()
spline_first_1000_foi_id_12<-spline_first_1000_foi_12$id
save(spline_first_1000_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_first_1000_12_knots_14_07_JUNE_5")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                     iteration_number = 100,
                                                                     data_about_sampling = data_about_sampling,
                                                                     params = params_cluster), name = "SPLINE_first_5000_foi_12_knot")
spline_first_5000_foi_12$status()
spline_first_5000_foi_id_12<-spline_first_5000_foi_12$id
save(spline_first_5000_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_first_5000_12_knots_14_08_JUNE_5")

#######################################################################################################################
## Now for the second order spline results ############################################################################
#######################################################################################################################

penalty_order<-2
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_100_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_100_foi_12_knot")
spline_sec_100_foi_12$status()
spline_sec_100_foi_id_12<-spline_sec_100_foi_12$id
save(spline_sec_100_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_sec_100_12_knots_14_08_JUNE_5")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                  iteration_number = 100,
                                                                  data_about_sampling = data_about_sampling,
                                                                  params = params_cluster), name = "SPLINE_sec_500_foi_12_knot")
spline_sec_500_foi_12$status()
spline_sec_500_foi_id_12<-spline_sec_500_foi_12$id
save(spline_sec_500_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_sec_500_12_knots_14_09_JUNE_5")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_sec_1000_foi_12_knot")
spline_sec_1000_foi_12$status()
spline_sec_1000_foi_id_12<-spline_sec_1000_foi_12$id
save(spline_sec_1000_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_sec_1000_12_knots_14_09_JUNE_5")
##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi_12<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                   iteration_number = 100,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params = params_cluster), name = "SPLINE_sec_5000_foi_12_knot")
spline_sec_5000_foi_12$status()
spline_sec_5000_foi_id_12<-spline_sec_5000_foi_12$id
save(spline_sec_5000_foi_id_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/SPLINE_sec_5000_12_knots_14_09_JUNE_5")
#########################################################################################################################################
## Now to load up the results of the cluster runs 12 12 12 12 12 12  ####################################################################
#########################################################################################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_ids/12_knots/"
twelve_knotters <- list.files(path_name,full.names = T)
for(i in 1:length(twelve_knotters)){
  load(twelve_knotters[i],verbose = T)
}

spline_f_100_foi_12<-obj$task_get(spline_first_100_foi_12_knot_id)
spline_f_500_foi_12<-obj$task_get(spline_first_500_foi_id_12)
spline_f_1k_foi_12<-obj$task_get(spline_first_1000_foi_id_12)
spline_f_5k_foi_12<-obj$task_get(spline_first_5000_foi_id_12)

spline_f_100_foi_res_12<-spline_f_100_foi_12$result()
spline_f_500_foi_res_12<-spline_f_500_foi_12$result()
spline_f_1k_foi_res_12<-spline_f_1k_foi_12$result()
spline_f_5k_foi_res_12<-spline_f_5k_foi_12$result()

save(spline_f_100_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/spline_100_1st_12_knots")
save(spline_f_500_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/spline_500_1st_12_knots")
save(spline_f_1k_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/spline_1000_1st_12_knots")
save(spline_f_5k_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/spline_5000_1st_12_knots")


###### spline second order ######

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id_12)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id_12)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id_12)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id_12)

spline_s_100_foi_res_12<-spline_s_100_foi$result()
spline_s_500_foi_res_12<-spline_s_500_foi$result()
spline_s_1k_foi_res_12<-spline_s_1k_foi$result()
spline_s_5k_foi_res_12<-spline_s_5k_foi$result()

save(spline_s_100_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/sp_sec_100_foi_12_knots")
save(spline_s_500_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/sp_sec_500_foi_12_knots")
save(spline_s_1k_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/sp_sec_1000_foi_12_knots")
save(spline_s_5k_foi_res_12,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/foi_as_modelled/cluster_results/12_knots/sp_sec_5000_foi_12_knots")








