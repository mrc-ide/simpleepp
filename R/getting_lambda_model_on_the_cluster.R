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

a<-load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_100_sampled_data")
b<-load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_500_sampled_data")
c<-load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_1000_sampled_data")
d<-load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/foi_as_modelled_parameter/n_5000_sampled_data")

plot(sampled_n_100_complete_data_foi_model[sampled_n_100_complete_data_foi_model$iteration==64,2])
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
  
  lines(sampled_n_100_complete_data_foi_model[sampled_n_100_complete_data_foi_model$iteration==i,2],col=colour)
  Sys.sleep(0.5)
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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FIRST_100_14_46_MAY_29")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_500_14_49_MAY_29")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_1000_14_51_MAY_29")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_5000_14_49_MAY_29")
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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_100_SECOND_14_58_MAY_29")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_500_SECOND_15_00_MAY_29")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_500_SECOND_15_02_MAY_29")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_5000_SECOND_15_03_MAY_29")


#######################################################################################################################
## Now for the splines ################################################################################################
#######################################################################################################################

knot_number<-7
penalty_order<-1
model_type<-"spline"
sample_years<-46
sample_n<-100

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_100_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                             iteration_number = 100,
                                                             data_about_sampling = data_about_sampling,
                                                             params = params_cluster), name = "SPLINE_first_100_foi")
spline_first_100_foi$wait(10)
spline_first_100_foi$status()
spline_first_100_foi_id<-spline_first_100_foi$id
save(spline_first_100_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_100_15_18_MAY_29")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_500_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_first_500_foi")
spline_first_100_foi$wait(10)
spline_first_500_foi$status()
spline_first_500_foi_id<-spline_first_500_foi$id
save(spline_first_500_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_500_15_18_MAY_29")

##### 1k #####

sample_n<-1000
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_1000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_first_1000_foi")
spline_first_100_foi$wait(10)
spline_first_1000_foi$status()
spline_first_1000_foi_id<-spline_first_1000_foi$id
save(spline_first_1000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_1000_15_18_MAY_29")

##### 5k #####

sample_n<-5000

data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_first_5000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_first_5000_foi")
spline_first_100_foi$wait(10)
spline_first_5000_foi$status()
spline_first_5000_foi_id<-spline_first_5000_foi$id
save(spline_first_5000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_5000_15_22_MAY_29")

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

spline_sec_100_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_100_complete_data_foi_model,
                                                                 iteration_number = 100,
                                                                 data_about_sampling = data_about_sampling,
                                                                 params = params_cluster), name = "SPLINE_sec_100_foi")
spline_sec_100_foi$wait(10)
spline_sec_100_foi$status()
spline_sec_100_foi_id<-spline_sec_100_foi$id
save(spline_sec_100_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_100_15_18_MAY_29")

##### 500 #####

sample_n<-500

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_500_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_500_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_500_foi")
spline_sec_500_foi$wait(10)
spline_sec_100_foi$status()
spline_sec_500_foi_id<-spline_sec_500_foi$id
save(spline_sec_500_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_500_15_35_MAY_29")

##### 1000 #####

sample_n<-1000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_1000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_1000_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_1000_foi")
spline_sec_1000_foi$wait(10)
spline_sec_100_foi$status()
spline_sec_1000_foi_id<-spline_sec_1000_foi$id
save(spline_sec_1000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_1000_15_36_MAY_29")

##### 5k #####

sample_n<-5000

params_cluster<-list(mu=mu,sigma=sigma,mu_i=mu_i,foi_flag=foi_flag)
data_about_sampling<-list(knot_number=knot_number,penalty_order=penalty_order,
                          model_type=model_type,sample_years=sample_years,
                          sample_n=sample_n)

spline_sec_5000_foi<-obj$enqueue(fitting_data_function_loop_foi(samples_data_frame = sampled_n_5000_complete_data_foi_model,
                                                               iteration_number = 100,
                                                               data_about_sampling = data_about_sampling,
                                                               params = params_cluster), name = "SPLINE_sec_5000_foi")
spline_sec_5000_foi$wait(12)
spline_sec_100_foi$status()
spline_sec_5000_foi_id<-spline_sec_5000_foi$id
save(spline_sec_5000_foi_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_5000_15_38_MAY_29")


######################################################################################################################
## Now we will load up the results ###################################################################################
#####################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FIRST_100_14_46_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_500_14_49_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_1000_14_51_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_FRIST_5000_14_49_MAY_29",
     verbose = T)

rw_f_100_foi<-obj$task_get(rw_first_100_foi_id)
rw_f_500_foi<-obj$task_get(rw_first_500_foi_id)
rw_f_1k_foi<-obj$task_get(rw_first_1000_foi_id)
rw_f_5k_foi<-obj$task_get(rw_first_5000_foi_id)

rw_f_100_foi_res<-rw_f_100_foi$result
rw_f_500_foi_res<-rw_f_500_foi$result
rw_f_1k_foi_res<-rw_f_1k_foi$result
rw_f_5k_foi_res<-rw_f_5k_foi$result


save(rw_f_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_100_foi")
save(rw_f_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_500_foi")
save(rw_f_1k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_1k_foi")
save(rw_f_5k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_5k_foi")

##### sec ord rw #####

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_100_SECOND_14_58_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_500_SECOND_15_00_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_500_SECOND_15_02_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/RW/RW_5000_SECOND_15_03_MAY_29",
     verbose = T)

rw_s_100_foi<-obj$task_get(rw_sec_100_foi_id)
rw_s_500_foi<-obj$task_get(rw_sec_500_foi_id)
rw_s_1k_foi<-obj$task_get(rw_sec_1000_foi_id)
rw_s_5k_foi<-obj$task_get(rw_sec_5000_foi_id)

rw_s_100_foi_res<-rw_s_100_foi$result
rw_s_500_foi_res<-rw_s_500_foi$result
rw_s_1k_foi_res<-rw_s_1k_foi$result
rw_s_5k_foi_res<-rw_s_5k_foi$result

save(rw_s_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_100_foi")
save(rw_s_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_500_foi")
save(rw_s_1k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_1k_foi")
save(rw_s_5k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_5k_foi")


##### spline first #####

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_100_15_18_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_500_15_18_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_1000_15_18_MAY_29",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_FIRST_5000_15_22_MAY_29",
     verbose = T)

spline_f_100_foi<-obj$task_get(spline_first_100_foi_id)
spline_f_500_foi<-obj$task_get(spline_first_500_foi_id)
spline_f_1k_foi<-obj$task_get(spline_first_1000_foi_id)
spline_f_5k_foi<-obj$task_get(spline_first_5000_foi_id)

spline_f_100_foi_res<-spline_f_100_foi$result
spline_f_500_foi_res<-spline_f_500_foi$result
spline_f_1k_foi_res<-spline_f_1k_foi$result
spline_f_5k_foi_res<-spline_f_5k_foi$result

save(spline_f_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_100_foi")
save(spline_f_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_500_foi")
save(spline_f_1k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_1k_foi")
save(spline_f_5k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_5k_foi")


###### spline second order ######

load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_100_15_18_MAY_29")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_500_15_35_MAY_29")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_1000_15_36_MAY_29")
load(verbose = T,
     "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline/SPLINE_SEC_5000_15_38_MAY_29")

spline_s_100_foi<-obj$task_get(spline_sec_100_foi_id)
spline_s_500_foi<-obj$task_get(spline_sec_500_foi_id)
spline_s_1k_foi<-obj$task_get(spline_sec_1000_foi_id)
spline_s_5k_foi<-obj$task_get(spline_sec_5000_foi_id)

spline_s_100_foi_res<-spline_s_100_foi$result
spline_s_500_foi_res<-spline_s_500_foi$result
spline_s_1k_foi_res<-spline_s_1k_foi$result
spline_s_5k_foi_res<-spline_s_5k_foi$result

save(spline_s_100_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_100_foi")
save(spline_s_500_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_500_foi")
save(spline_s_1k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_1k_foi")
save(spline_s_5k_foi_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_5k_foi")





