#######################################################
## Setting up work for the DIDE cluster ###############
#######################################################

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

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),sources = "simpleepp/R/loop_functions_random_walk.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

################################################################################################################################
## So the above lines of code give me access to the cluster, with the obj object giving me queue functionalaity, I've also #####
## asked for three cores for stan to use for each of the chains ################################################################
################################################################################################################################
obj$cluster_load()
obj$task_status()

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/true_epidemic_data",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_100_complete_data_no_art",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_500_complete_data_no_art", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_1000_complete_data_no_art",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_5000_complete_data_no_art", verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1984_simpleepp/N_100_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1984_simpleepp/N_500_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1984_simpleepp/N_1000_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1984_simpleepp/N_5000_samples",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1990_runs/N_100_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1990_runs/N_500_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1990_runs/N_1000_samples",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1990_runs/N_5000_samples",verbose = T)


mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

params<-list(mu=mu,mu_i=mu_i,sigma=sigma)
params

###################################################################################################################################
## Now we will run through the data from begininng ################################################################################
###################################################################################################################################

sample_range<-1970:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df), 
                                       name = "data_70_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_100_FIRST_12_16_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_70_rw_1st_n_500")

n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_500_FIRST_12_17_JUNE_11")

#### 1k ####
sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_70_rw_1st_n_1000")

n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_1000_FIRST_12_18_JUNE_11")

##### 5k ######
sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_70_rw_1st_n_5000")

n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_5000_FIRST_12_18_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_70_rw_2nd_n_100")

n_100_RW_sec_order_loop$status()
n_100_RW_sec_order_loop_id<-n_100_RW_sec_order_loop$id
save(n_100_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_100_second_12_19_JUNE_11")

#### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_70_rw_2nd_n_500")

n_500_RW_sec_order_loop$status()
n_500_RW_sec_order_loop_id<-n_500_RW_sec_order_loop$id
save(n_500_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_500_second_12_19_JUNE_11")


sample_n <- 1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_70_rw_2nd_n_1000")

n_1000_RW_sec_order_loop$status()
n_1000_RW_sec_order_loop_id<-n_1000_RW_sec_order_loop$id
save(n_1000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_1000_second_12_19_JUNE_11")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_70_rw_2nd_n_5000")

n_5000_RW_sec_order_loop$status()
n_5000_RW_sec_order_loop_id<-n_5000_RW_sec_order_loop$id
save(n_5000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/RW_5000_second_12_20_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_70_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_100_first_12_20_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_70_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_500_first_12_20_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_70_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_1000_first_12_21_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_70_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_5000_first_12_22_JUNE_11")



################################################################################################################################
## Now we'll go for the second order work ######################################################################################
################################################################################################################################

knot_number <- 7
penalty_order <- 2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_70_100")
spline_second_order_n_100$status()
spline_second_order_n_100$log()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_100_second_12_23_JUNE_11")


sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_70_500")
spline_second_order_n_500$status()
spline_second_order_n_500_id<-spline_second_order_n_500$id
save(spline_second_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_500_second_12_23_JUNE_11")

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                          data_about_sampling = data_about_sampling,
                                                                          iteration_number = 100,params = params,
                                                                          simulated_true_df = sim_model_output$sim_df),
                                        name = "Spline_2nd_70_1000")

spline_second_order_n_1000$status()
spline_second_order_n_1000_id<-spline_second_order_n_1000$id
save(spline_second_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_1000_second_12_24_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                          data_about_sampling = data_about_sampling,
                                                                          iteration_number = 100,params = params,
                                                                          simulated_true_df = sim_model_output$sim_df),
                                        name = "Spline_2nd_70_5000")
spline_second_order_n_5000$status()
spline_second_order_n_5000_id<-spline_second_order_n_5000$id
save(spline_second_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/spline_5000_second_12_24_JUNE_11")

####################################################################################################
## Load up the results #############################################################################
####################################################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_ids/"
seventies_data <- list.files(path_name,full.names = T)
for(i in 1:length(seventies_data)){
  load(seventies_data[i],verbose = T)
}


### RW first 

rw_70_1_100 <- obj$task_get(n_100_RW_first_order_loop_id)
rw_70_1_500 <- obj$task_get(n_500_RW_first_order_loop_id)
rw_70_1_1000 <- obj$task_get(n_1000_RW_first_order_loop_id)
rw_70_1_5000 <- obj$task_get(n_5000_RW_first_order_loop_id)

rw_70_2_100 <- obj$task_get(n_100_RW_sec_order_loop_id)
rw_70_2_500 <- obj$task_get(n_500_RW_sec_order_loop_id)
rw_70_2_1000 <- obj$task_get(n_1000_RW_sec_order_loop_id)
rw_70_2_5000 <- obj$task_get(n_5000_RW_sec_order_loop_id)

rw_70_1_100_res <- rw_70_1_100$result()
rw_70_1_500_res <- rw_70_1_500$result()
rw_70_1_1000_res <- rw_70_1_1000$result()
rw_70_1_5000_res <- rw_70_1_5000$result()

rw_70_2_100_res <- rw_70_2_100$result()
rw_70_2_500_res <- rw_70_2_500$result()
rw_70_2_1000_res <- rw_70_2_1000$result()
rw_70_2_5000_res <- rw_70_2_5000$result()

list_of_results <- list(rw_70_1_100_res,rw_70_1_500_res,rw_70_1_1000_res,rw_70_1_5000_res,
                        rw_70_2_100_res,rw_70_2_500_res,rw_70_2_1000_res,rw_70_2_5000_res)
for(i in 1:length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  
  Sys.sleep(1.5)
}

save(rw_70_1_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_1st_n_100")
save(rw_70_1_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_1st_n_500")
save(rw_70_1_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_1st_n_1000")
save(rw_70_1_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_1st_n_5000")
save(rw_70_2_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_2nd_n_100")
save(rw_70_2_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_2nd_n_500")
save(rw_70_2_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_2nd_n_1000")
save(rw_70_2_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/RW_70_2nd_n_5000")


## Now for the splines

sp_70_1_100 <- obj$task_get(spline_first_order_n_100_id)
sp_70_1_500 <- obj$task_get(spline_first_order_n_500_id)
sp_70_1_1000 <- obj$task_get(spline_first_order_n_1000_id)
sp_70_1_5000 <- obj$task_get(spline_first_order_n_5000_id)

sp_70_2_100 <- obj$task_get(spline_second_order_n_100_id)
sp_70_2_500 <- obj$task_get(spline_second_order_n_500_id)
sp_70_2_1000 <- obj$task_get(spline_second_order_n_1000_id)
sp_70_2_5000 <- obj$task_get(spline_second_order_n_5000_id)

sp_70_1_100_res <- sp_70_1_100$result()
sp_70_1_500_res <- sp_70_1_500$result()
sp_70_1_1000_res <- sp_70_1_1000$result()
sp_70_1_5000_res <- sp_70_1_5000$result()

sp_70_2_100_res <- sp_70_2_100$result()
sp_70_2_500_res <- sp_70_2_500$result()
sp_70_2_1000_res <- sp_70_2_1000$result()
sp_70_2_5000_res <- sp_70_2_5000$result()

list_of_results <- list(sp_70_1_100_res,sp_70_1_500_res,sp_70_1_5000_res,
                        sp_70_2_100_res,sp_70_2_500_res,sp_70_2_1000_res,sp_70_2_5000_res)
for(i in 1:length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  Sys.sleep(1.5)
}

save(sp_70_1_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_1st_70_n_100")
save(sp_70_1_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_1st_70_n_500")
save(sp_70_1_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_1st_70_n_1000")
save(sp_70_1_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_1st_70_n_5000")

save(sp_70_2_100_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_2nd_70_n_100")
save(sp_70_2_500_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_2nd_70_n_100")
save(sp_70_2_1000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_2nd_70_n_100")
save(sp_70_2_5000_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/simplepp_early_sampling/cluster_results/SP_2nd_70_n_100")

##################################################################################################################################
## Data form 1984 onwards ########################################################################################################
##################################################################################################################################

obj$cluster_load()

sample_range<-1984:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_84_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_84_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_100_FIRST_12_16_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_84_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_84_rw_1st_n_500")

n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_500_FIRST_12_17_JUNE_11")

#### 1k ####
sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_84_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_84_rw_1st_n_1000")

n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_1000_FIRST_12_18_JUNE_11")

##### 5k ######
sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_84_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_84_rw_1st_n_5000")

n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_5000_FIRST_12_18_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_84_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_84_rw_2nd_n_100")

n_100_RW_first_order_loop$status()
n_100_RW_sec_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_100_second_12_19_JUNE_11")

#### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_84_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_84_rw_2nd_n_500")

n_500_RW_sec_order_loop$status()
n_500_RW_sec_order_loop_id<-n_500_RW_sec_order_loop$id
save(n_500_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_500_second_12_19_JUNE_11")


sample_n <- 1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_84_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_84_rw_2nd_n_1000")

n_1000_RW_sec_order_loop$status()
n_1000_RW_sec_order_loop_id<-n_1000_RW_sec_order_loop$id
save(n_1000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_1000_second_12_19_JUNE_11")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_84_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_84_rw_2nd_n_5000")

n_5000_RW_sec_order_loop$status()
n_5000_RW_sec_order_loop_id<-n_5000_RW_sec_order_loop$id
save(n_5000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/RW_5000_second_12_20_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_84_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),name = "spline_84_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_100_first_12_20_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_84_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_84_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_500_first_12_20_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_84_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_84_n_1000")


spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_1000_first_12_21_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_84_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_84_n_5000")


spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_5000_first_12_22_JUNE_11")



################################################################################################################################
## Now we'll go for the second order work ######################################################################################
################################################################################################################################

knot_number <- 7
penalty_order <- 2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_84_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_84_100")
spline_second_order_n_100$status()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_100_second_12_23_JUNE_11")


sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_84_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_84_500")
spline_second_order_n_500$status()
spline_second_order_n_500_id<-spline_second_order_n_500$id
save(spline_second_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_500_second_12_23_JUNE_11")

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_84_data,
                                                                          data_about_sampling = data_about_sampling,
                                                                          iteration_number = 100,params = params,
                                                                          simulated_true_df = sim_model_output$sim_df),
                                        name = "Spline_2nd_84_1000")
spline_second_order_n_1000$status()
spline_second_order_n_1000_id<-spline_second_order_n_1000$id
save(spline_second_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_1000_second_12_24_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_84_data,
                                                                          data_about_sampling = data_about_sampling,
                                                                          iteration_number = 100,params = params,
                                                                          simulated_true_df = sim_model_output$sim_df),
                                        name = "Spline_2nd_84_5000")
spline_second_order_n_5000$status()
spline_second_order_n_5000_id<-spline_second_order_n_5000$id
save(spline_second_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/spline_5000_second_12_24_JUNE_11")

#############################################################################
## Load up the 1984 runs ####################################################
#############################################################################

path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1984/"
eighty_four <- list.files(path_name, full.names = T)
for(i in 1:length(eighty_four)){
  load(eighty_four[i], verbose = T)
}


rw_1_100_84 <- obj$task_get(n_100_RW_first_order_loop_id)
rw_1_500_84 <- obj$task_get(n_500_RW_first_order_loop_id)
rw_1_1000_84 <- obj$task_get(n_1000_RW_first_order_loop_id)
rw_1_5000_84 <- obj$task_get(n_5000_RW_first_order_loop_id)

rw_1_100_84_res <- rw_1_100_84$result()
rw_1_500_84_res <- rw_1_500_84$result()
rw_1_1000_84_res <- rw_1_1000_84$result()
rw_1_5000_84_res <- rw_1_5000_84$result()

rw_1_500_84$status()
rw_1_500_84$log()

rw_2_100_84 <- obj$task_get(n_100_RW_sec_order_loop_id)
rw_2_500_84 <- obj$task_get(n_500_RW_sec_order_loop_id)
rw_2_1000_84 <- obj$task_get(n_1000_RW_sec_order_loop_id)
rw_2_5000_84 <- obj$task_get(n_5000_RW_sec_order_loop_id)

rw_2_100_84_res <- rw_2_100_84$result()
rw_2_500_84_res <- rw_2_500_84$result()
rw_2_1000_84_res <- rw_2_1000_84$result()
rw_2_5000_84_res <- rw_2_5000_84$result()

list_of_results <- list(rw_1_100_84_res,rw_1_500_84_res,rw_1_1000_84_res,rw_1_5000_84_res,
                        rw_2_100_84_res,rw_2_1000_84_res,rw_2_5000_84_res)

for(i in 1: length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  Sys.sleep(1.4)
}

save(rw_1_100_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_1st_84_n_100")
save(rw_1_500_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_1st_84_n_500")
save(rw_1_1000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_1st_84_n_1000")
save(rw_1_5000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_1st_84_n_5000")

save(rw_2_100_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_2nd_84_n_100")
save(rw_2_500_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_2nd_84_n_500")
save(rw_2_1000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_2nd_84_n_1000")
save(rw_2_5000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/RW_2nd_84_n_5000")

## Now for splines 

sp_1_100_84 <- obj$task_get(spline_first_order_n_100_id)
sp_1_500_84 <- obj$task_get(spline_first_order_n_500_id)
sp_1_1000_84 <- obj$task_get(spline_first_order_n_1000_id)
sp_1_5000_84 <- obj$task_get(spline_first_order_n_5000_id)

sp_2_100_84 <- obj$task_get(spline_second_order_n_100_id)
sp_2_500_84 <- obj$task_get(spline_second_order_n_500_id)
sp_2_1000_84 <- obj$task_get(spline_second_order_n_1000_id)
sp_2_5000_84 <- obj$task_get(spline_second_order_n_5000_id)

sp_2_5000_84$status()
sp_2_5000_84$log()

sp_1_100_84_res <- sp_1_100_84$result()
sp_1_500_84_res <- sp_1_500_84$result()
sp_1_1000_84_res <- sp_1_1000_84$result()
sp_1_5000_84_res <- sp_1_5000_84$result()

sp_2_100_84_res <- sp_2_100_84$result()
sp_2_500_84_res <- sp_2_500_84$result()
sp_2_1000_84_res <- sp_2_1000_84$result()
sp_2_5000_84_res <- sp_2_5000_84$result()

list_of_results <- list(sp_1_100_84_res,sp_1_500_84_res,sp_1_1000_84_res,sp_1_5000_84_res,
                        sp_2_100_84_res,sp_2_500_84_res,sp_2_1000_84_res,sp_2_5000_84_res)
for(i in 1:length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  Sys.sleep(1.4)
}


save(sp_1_100_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_1st_84_n_100")
save(sp_1_500_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_1st_84_n_500")
save(sp_1_1000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_1st_84_n_1000")
save(sp_1_5000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_1st_84_n_5000")

save(sp_2_100_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_2nd_84_n_100")
save(sp_2_500_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_2nd_84_n_500")
save(sp_2_1000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_2nd_84_n_1000")
save(sp_2_5000_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1984/SP_2nd_84_n_5000")



##################################################################################################################################
## Data from 1990 onwards ########################################################################################################
##################################################################################################################################

obj$cluster_load()

sample_range<-1990:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_90_data,
                                            data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                            simulated_true_df = sim_model_output$sim_df), name = "data_90_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_100_FIRST_12_16_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_90_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                       name = "data_90_rw_1st_n_500")


n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_500_FIRST_12_17_JUNE_11")

#### 1k ####
sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_90_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_90_rw_1st_n_1000")

n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_1000_FIRST_12_18_JUNE_11")

##### 5k ######
sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_90_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_90_rw_1st_n_5000")

n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_5000_FIRST_12_18_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_90_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_90_rw_2nd_n_100")

n_100_RW_sec_order_loop$status()
n_100_RW_sec_order_loop_id<-n_100_RW_sec_order_loop$id
save(n_100_RW_sec_order_loop_id,
   file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_100_second_12_19_JUNE_11")

#### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_90_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_90_rw_2nd_n_500")

n_500_RW_sec_order_loop$status()
n_500_RW_sec_order_loop_id<-n_500_RW_sec_order_loop$id
save(n_500_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_500_second_12_19_JUNE_11")


sample_n <- 1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_90_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_90_rw_2nd_n_1000")

n_1000_RW_sec_order_loop$status()
n_1000_RW_sec_order_loop_id<-n_1000_RW_sec_order_loop$id
save(n_1000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_1000_second_12_19_JUNE_11")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_90_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_90_rw_2nd_n_5000")

n_5000_RW_sec_order_loop$status()
n_5000_RW_sec_order_loop_id<-n_5000_RW_sec_order_loop$id
save(n_5000_RW_sec_order_loop_id,
   file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/RW_5000_second_12_20_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_90_data,
                                  data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                  simulated_true_df = sim_model_output$sim_df),name = "spline_1st_90_n_100")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_100_first_12_20_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_90_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_90_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_500_first_12_20_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_90_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_90_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_1000_first_12_21_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_90_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_90_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_5000_first_12_22_JUNE_11")



################################################################################################################################
## Now we'll go for the second order work ######################################################################################
################################################################################################################################

knot_number <- 7
penalty_order <- 2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_90_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_90_100")
spline_second_order_n_100$status()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_100_second_12_23_JUNE_11")


sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_90_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_90_500")
spline_second_order_n_500$status()
spline_second_order_n_500_id<-spline_second_order_n_500$id
save(spline_second_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_500_second_12_23_JUNE_11")

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_90_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_90_1000")
spline_second_order_n_1000$status()
spline_second_order_n_1000_id<-spline_second_order_n_1000$id
save(spline_second_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_1000_second_12_24_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=length(rows_to_evaluate),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_90_data,
                                                                          data_about_sampling = data_about_sampling,
                                                                          iteration_number = 100,params = params,
                                                                          simulated_true_df = sim_model_output$sim_df),
                                        name = "Spline_2nd_90_5000")
spline_second_order_n_5000$status()
spline_second_order_n_5000_id<-spline_second_order_n_5000$id
save(spline_second_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/spline_5000_second_12_24_JUNE_11")


#######################################################################################################################################
## Now lets load up the results #######################################################################################################
#######################################################################################################################################
path_name <- "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1990/"
ninety_runs <- list.files(path_name,full.names = T)
for(i in 1:length(ninety_runs)){
  load(ninety_runs[i],verbose = T)
}


rw_100_1_90<-obj$task_get(n_100_RW_first_order_loop_id)
rw_500_1_90<-obj$task_get(n_500_RW_first_order_loop_id)
rw_1000_1_90<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_5000_1_90<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_100_2_90<-obj$task_get(n_100_RW_sec_order_loop_id)
rw_500_2_90<-obj$task_get(n_500_RW_sec_order_loop_id)
rw_1000_2_90<-obj$task_get(n_1000_RW_sec_order_loop_id)
rw_5000_2_90<-obj$task_get(n_5000_RW_sec_order_loop_id)

sp_100_1_90<-obj$task_get(spline_first_order_n_100_id)
sp_500_1_90<-obj$task_get(spline_first_order_n_500_id)
sp_1000_1_90<-obj$task_get(spline_first_order_n_1000_id)
sp_5000_1_90<-obj$task_get(spline_first_order_n_5000_id)

sp_100_2_90<-obj$task_get(spline_second_order_n_100_id)
sp_500_2_90<-obj$task_get(spline_second_order_n_500_id)
sp_1000_2_90<-obj$task_get(spline_second_order_n_1000_id)
sp_5000_2_90<-obj$task_get(spline_second_order_n_5000_id)


#### now for the results ####

rw_100_1_90_res<-rw_100_1_90$result()
rw_500_1_90_res<-rw_500_1_90$result()
rw_1000_1_90_res<-rw_1000_1_90$result()
rw_5000_1_90_res<-rw_5000_1_90$result()

rw_100_2_90_res<-rw_100_2_90$result()
rw_500_2_90_res<-rw_500_2_90$result()
rw_1000_2_90_res<-rw_1000_2_90$result()
rw_5000_2_90_res<-rw_5000_2_90$result()

sp_100_1_90_res<-sp_100_1_90$result()
sp_500_1_90_res<-sp_500_1_90$result()
sp_1000_1_90_res<-sp_1000_1_90$result()
sp_5000_1_90_res<-sp_5000_1_90$result()

sp_100_2_90_res<-sp_100_2_90$result()
sp_500_2_90_res<-sp_500_2_90$result()
sp_1000_2_90_res<-sp_1000_2_90$result()
sp_5000_2_90_res<-sp_5000_2_90$result()


list_of_results <- list(rw_100_1_90_res,rw_500_1_90_res,rw_1000_1_90_res,rw_5000_1_90_res,
                        rw_100_2_90_res,rw_500_2_90_res,rw_1000_2_90_res,rw_5000_2_90_res,
                        sp_100_1_90_res,sp_500_1_90_res,sp_1000_1_90_res,sp_5000_1_90_res,
                        sp_100_2_90_res,sp_500_2_90_res,sp_1000_2_90_res,sp_5000_2_90_res)

for(i in 1:length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  Sys.sleep(1.4)
}

save(rw_100_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_1st_90_100")
save(rw_500_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_1st_90_500")
save(rw_1000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_1st_90_1000")
save(rw_5000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_1st_90_5000")

save(rw_100_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_2nd_90_100")
save(rw_500_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_2nd_90_500")
save(rw_1000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_2nd_90_1000")
save(rw_5000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/RW_2nd_90_5000")

save(sp_100_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_1st_90_100")
save(sp_500_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_1st_90_500")
save(sp_1000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_1st_90_1000")
save(sp_5000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_1st_90_5000")

save(sp_100_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_2nd_90_100")
save(sp_500_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_2nd_90_500")
save(sp_1000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_2nd_90_1000")
save(sp_5000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1990/SP_2nd_90_5000")

###########################################################################################################################################
## Now we will run this from 1995 for sampling data, and then finally 2000 ####################111#########################################
###############################################################################################666#########################################

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1995_runs/N_100_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1995_runs/N_500_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1995_runs/N_1000_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_1995_runs/N_5000_sampled_data",verbose = T)

sample_range<-1995:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_95_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_95_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_100_FIRST_14_42_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_95_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_95_rw_1st_n_500")

n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_500_FIRST_14_43_JUNE_11")

#### 1k ####
sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_95_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_95_rw_1st_n_1000")

n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_1000_FIRST_14_43_JUNE_11")

##### 5k ######
sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_95_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_95_rw_1st_n_5000")

n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_5000_FIRST_14_44_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_95_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_95_rw_2nd_n_100")

n_100_RW_sec_order_loop$status()
n_100_RW_sec_order_loop_id<-n_100_RW_sec_order_loop$id
save(n_100_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_100_second_14_45_JUNE_11")

#### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_95_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_95_rw_2nd_n_500")

n_500_RW_sec_order_loop$status()
n_500_RW_sec_order_loop_id<-n_500_RW_sec_order_loop$id
save(n_500_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_500_second_14_45_JUNE_11")


sample_n <- 1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_95_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_95_rw_2nd_n_1000")

n_1000_RW_sec_order_loop$status()
n_1000_RW_sec_order_loop_id<-n_1000_RW_sec_order_loop$id
save(n_1000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_1000_second_14_45_JUNE_11")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_95_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_95_rw_2nd_n_5000")

n_5000_RW_sec_order_loop$status()
n_5000_RW_sec_order_loop_id<-n_5000_RW_sec_order_loop$id
save(n_5000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/RW_5000_second_14_46_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_95_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_95_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_100_first_14_47_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_95_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_95_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_500_first_14_47_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_95_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_95_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_1000_first_14_47_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_95_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_95_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_5000_first_14_48_JUNE_11")

########################## NOw for the second orderr splines ###########
penalty_order <- 2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_100_95<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_95_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params=params,
                                                                   iteration_number=100,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                 name = "spline_2nd_95_n_100")

spline_sec_n_100_95$status()
spline_sec_n_100_95_id<-spline_sec_n_100_95$id
save(spline_sec_n_100_95_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_100_second_14_56_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_500_95<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_95_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params=params,
                                                                   iteration_number=100,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                 name = "spline_2nd_95_n_500")

spline_sec_n_500_95$status()
spline_sec_n_500_95_id<-spline_sec_n_500_95$id
save(spline_sec_n_500_95_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_500_second_15_31_JUNE_11")

##### 1k #####
sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_1000_95<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_95_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params=params,
                                                                   iteration_number=100,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                 name = "spline_2nd_95_n_1000")

spline_sec_n_1000_95$status()
spline_sec_n_1000_95_id<-spline_sec_n_1000_95$id
save(spline_sec_n_1000_95_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_500_second_15_32_JUNE_11")

#### 5000 #####

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_5000_95<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_95_data,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params=params,
                                                                    iteration_number=100,
                                                                    simulated_true_df = sim_model_output$sim_df),
                                  name = "spline_2nd_95_n_5000")

spline_sec_n_5000_95$status()
spline_sec_n_5000_95_id<-spline_sec_n_5000_95$id
save(spline_sec_n_5000_95_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/spline_5000_second_15_33_JUNE_11")


###### lets load up the results for these cluster runs ###########

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/1995/"
ids_1995<-list.files(path_name,full.names = T)
for(i in 1:length(ids_1995)){
  load(ids_1995[i],verbose = T)
}

rw_100_1_95<-obj$task_get(n_100_RW_first_order_loop_id)
rw_500_1_95<-obj$task_get(n_500_RW_first_order_loop_id)
rw_1000_1_95<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_5000_1_95<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_100_2_95<-obj$task_get(n_100_RW_sec_order_loop_id)
rw_500_2_95<-obj$task_get(n_500_RW_sec_order_loop_id)
rw_1000_2_95<-obj$task_get(n_1000_RW_sec_order_loop_id)
rw_5000_2_95<-obj$task_get(n_5000_RW_sec_order_loop_id)

sp_100_1_95<-obj$task_get(spline_first_order_n_100_id)
sp_500_1_95<-obj$task_get(spline_first_order_n_500_id)
sp_1000_1_95<-obj$task_get(spline_first_order_n_1000_id)
sp_5000_1_95<-obj$task_get(spline_first_order_n_5000_id)

sp_100_2_95<-obj$task_get(spline_sec_n_100_95_id)
sp_500_2_95<-obj$task_get(spline_sec_n_500_95_id)
sp_1000_2_95<-obj$task_get(spline_sec_n_1000_95_id)
sp_5000_2_95<-obj$task_get(spline_sec_n_5000_95_id)


#### now for the results ####

rw_100_1_95_res<-rw_100_1_95$result()
rw_500_1_95_res<-rw_500_1_95$result()
rw_1000_1_95_res<-rw_1000_1_95$result()
rw_5000_1_95_res<-rw_5000_1_95$result()

rw_100_2_95_res<-rw_100_2_95$result()
rw_500_2_95_res<-rw_500_2_95$result()
rw_1000_2_95_res<-rw_1000_2_95$result()
rw_5000_2_95_res<-rw_5000_2_95$result()

sp_100_1_95_res<-sp_100_1_95$result()
sp_500_1_95_res<-sp_500_1_95$result()
sp_1000_1_95_res<-sp_1000_1_95$result()
sp_5000_1_95_res<-sp_5000_1_95$result()

sp_100_2_95_res<-sp_100_2_95$result()
sp_500_2_95_res<-sp_500_2_95$result()
sp_1000_2_95_res<-sp_1000_2_95$result()
sp_5000_2_95_res<-sp_5000_2_95$result()

list_of_results <- list(rw_100_1_95_res,rw_500_1_95_res,rw_1000_1_95_res,rw_5000_1_95_res,
                        rw_100_2_95_res,rw_500_2_95_res,rw_1000_2_95_res,rw_5000_2_95_res,
                        sp_100_1_95_res,sp_500_1_95_res,sp_1000_1_95_res,sp_5000_1_95_res,
                        sp_100_2_95_res,sp_500_2_95_res,sp_1000_2_95_res,sp_5000_2_95_res)

for(i in 1:length(list_of_results)){
  print(list_of_results[[i]]$data_about_run)
  Sys.sleep(1.4)
}


save(rw_100_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_1st_95_100")
save(rw_500_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_1st_95_500")
save(rw_1000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_1st_95_1000")
save(rw_5000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_1st_95_5000")

save(rw_100_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_2nd_95_100")
save(rw_500_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_2nd_95_500")
save(rw_1000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_2nd_95_1000")
save(rw_5000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/RW_2nd_95_5000")

save(sp_100_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_1st_95_100")
save(sp_500_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_1st_95_500")
save(sp_1000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_1st_95_1000")
save(sp_5000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_1st_95_5000")

save(sp_100_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_2nd_95_100")
save(sp_500_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_2nd_95_500")
save(sp_1000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_2nd_95_1000")
save(sp_5000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_results/1995/SP_2nd_95_5000")


############################################################################################################################################
## NOw we will do this for the 2000 data ###################################################################################################
############################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_2000_runs/N_100_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_2000_runs/N_500_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_2000_runs/N_1000_sampled_data",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/data_from_2000_runs/N_5000_sampled_data",verbose = T)

sample_range<-2000:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_00_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df), name = "data_00_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_100_FIRST_16_04_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_00_data,
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df),
                                       name = "data_00_rw_1st_n_500")

n_500_RW_first_order_loop$status()
n_500_RW_first_order_loop_id<-n_500_RW_first_order_loop$id
save(n_500_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_500_FIRST_16_05_JUNE_11")

#### 1k ####
sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_00_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_00_rw_1st_n_1000")

n_1000_RW_first_order_loop$status()
n_1000_RW_first_order_loop_id<-n_1000_RW_first_order_loop$id
save(n_1000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_1000_FIRST_16_05_JUNE_11")

##### 5k ######
sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_00_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_00_rw_1st_n_5000")

n_5000_RW_first_order_loop$status()
n_5000_RW_first_order_loop_id<-n_5000_RW_first_order_loop$id
save(n_5000_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_5000_FIRST_16_05_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_00_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_00_rw_2nd_n_100")

n_100_RW_sec_order_loop$status()
n_100_RW_sec_order_loop_id<-n_100_RW_sec_order_loop$id
save(n_100_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_100_second_16_06_JUNE_11")

#### 500 ####

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_00_data,
                                                                data_about_sampling = data_about_sampling,
                                                                iteration_number = 100,params = params,
                                                                simulated_true_df = sim_model_output$sim_df),
                                     name = "data_00_rw_2nd_n_500")

n_500_RW_sec_order_loop$status()
n_500_RW_sec_order_loop_id<-n_500_RW_sec_order_loop$id
save(n_500_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_500_second_16_06_JUNE_11")


sample_n <- 1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_1000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_00_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_00_rw_2nd_n_1000")

n_1000_RW_sec_order_loop$status()
n_1000_RW_sec_order_loop_id<-n_1000_RW_sec_order_loop$id
save(n_1000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_1000_second_16_06_JUNE_11")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_sec_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_00_data,
                                                                 data_about_sampling = data_about_sampling,
                                                                 iteration_number = 100,params = params,
                                                                 simulated_true_df = sim_model_output$sim_df),
                                      name = "data_00_rw_2nd_n_5000")

n_5000_RW_sec_order_loop$status()
n_5000_RW_sec_order_loop_id<-n_5000_RW_sec_order_loop$id
save(n_5000_RW_sec_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/RW_5000_second_16_06_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_00_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_00_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_100_first_16_07_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_00_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_00_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500_id<-spline_first_order_n_500$id
save(spline_first_order_n_500_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_500_first_16_07_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_00_data,
                                                                         data_about_sampling = data_about_sampling,
                                                                         iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_00_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000_id<-spline_first_order_n_1000$id
save(spline_first_order_n_1000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_1000_first_16_07_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_00_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_00_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000_id<-spline_first_order_n_5000$id
save(spline_first_order_n_5000_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_5000_first_16_08_JUNE_11")

########################## NOw for the second orderr splines ###########
penalty_order <- 2
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_100_00<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_00_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params=params,
                                                                   iteration_number=100,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                 name = "spline_2nd_00_n_100")

spline_sec_n_100_00$status()
spline_sec_n_100_00_id<-spline_sec_n_100_00$id
save(spline_sec_n_100_00_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_100_second_16_08_JUNE_11")

### 500 ####

sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_500_00<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_00_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   params=params,
                                                                   iteration_number=100,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                 name = "spline_2nd_00_n_500")

spline_sec_n_500_00$status()
spline_sec_n_500_00_id<-spline_sec_n_500_00$id
save(spline_sec_n_500_00_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_500_second_16_08_JUNE_11")

##### 1k #####
sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_1000_00<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_00_data,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params=params,
                                                                    iteration_number=100,
                                                                    simulated_true_df = sim_model_output$sim_df),
                                  name = "spline_2nd_00_n_1000")

spline_sec_n_1000_00$status()
spline_sec_n_1000_00_id<-spline_sec_n_1000_00$id
save(spline_sec_n_1000_00_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_500_second_16_09_JUNE_11")

#### 5000 #####

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_sec_n_5000_00<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_00_data,
                                                                    data_about_sampling = data_about_sampling,
                                                                    params=params,
                                                                    iteration_number=100,
                                                                    simulated_true_df = sim_model_output$sim_df),
                                  name = "spline_2nd_00_n_5000")

spline_sec_n_5000_00$status()
spline_sec_n_5000_00_id<-spline_sec_n_5000_00$id
save(spline_sec_n_5000_00_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/spline_5000_second_16_09_JUNE_11")

############################################################################################################################
###### lets load up the results for these cluster runs #####################################################################
############################################################################################################################

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/log_penalized/varying_data_start/cluster_ids/2000/"
ids_2000<-list.files(path_name,full.names = T)
for(i in 1:length(ids_2000)){
  load(ids_2000[i],verbose = T)
}

rw_100_1_00<-obj$task_get(n_100_RW_first_order_loop_id)
rw_500_1_00<-obj$task_get(n_500_RW_first_order_loop_id)
rw_1000_1_00<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_5000_1_00<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_100_2_00<-obj$task_get(n_100_RW_sec_order_loop_id)
rw_500_2_00<-obj$task_get(n_500_RW_sec_order_loop_id)
rw_1000_2_00<-obj$task_get(n_1000_RW_sec_order_loop_id)
rw_5000_2_00<-obj$task_get(n_5000_RW_sec_order_loop_id)

sp_100_1_00<-obj$task_get(spline_first_order_n_100_id)
sp_500_1_00<-obj$task_get(spline_first_order_n_500_id)
sp_1000_1_00<-obj$task_get(spline_first_order_n_1000_id)
sp_5000_1_00<-obj$task_get(spline_first_order_n_5000_id)

sp_100_2_00<-obj$task_get(spline_sec_n_100_00_id)
sp_500_2_00<-obj$task_get(spline_sec_n_500_00_id)
sp_1000_2_00<-obj$task_get(spline_sec_n_1000_00_id)
sp_5000_2_00<-obj$task_get(spline_sec_n_5000_00_id)


#### now for the results ####

rw_100_1_00_res<-rw_100_1_00$result()
rw_500_1_00_res<-rw_500_1_00$result()
rw_1000_1_00_res<-rw_1000_1_00$result()
rw_5000_1_00_res<-rw_5000_1_00$result()

rw_100_2_00_res<-rw_100_2_00$result()
rw_500_2_00_res<-rw_500_2_00$result()
rw_1000_2_00_res<-rw_1000_2_00$result()
rw_5000_2_00_res<-rw_5000_2_00$result()

sp_100_1_00_res<-sp_100_1_00$result()
sp_500_1_00_res<-sp_500_1_00$result()
sp_1000_1_00_res<-sp_1000_1_00$result()
sp_5000_1_00_res<-sp_5000_1_00$result()

sp_100_2_00_res<-sp_100_2_00$result()
sp_500_2_00_res<-sp_500_2_00$result()
sp_1000_2_00_res<-sp_1000_2_00$result()
sp_5000_2_00_res<-sp_5000_2_00$result()

sp_5000_2_00_res$data_about_run


save(rw_100_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_100_FIRSt_ORDER")
save(rw_500_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_500_FIRSt_ORDER")
save(rw_1000_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_1000_FIRSt_ORDER")
save(rw_5000_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_5000_FIRSt_ORDER")

save(rw_100_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_100_SECOND_ORDER")
save(rw_500_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_500_SECOND_ORDER")
save(rw_1000_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_1000_SECOND_ORDER")
save(rw_5000_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/RW_5000_SECOND_ORDER")

save(sp_100_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_100_FIRSt_ORDER")
save(sp_500_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_500_FIRSt_ORDER")
save(sp_1000_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_1000_FIRSt_ORDER")
save(sp_5000_1_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_5000_FIRSt_ORDER")

save(sp_100_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_100_SECOND_ORDER")
save(sp_500_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_500_SECOND_ORDER")
save(sp_1000_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_1000_SECOND_ORDER")
save(sp_5000_2_00_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/SPLINE_5000_SECOND_ORDER")



plot(sp_5000_1_00_res$prev$median[1:501])
plot(sp_5000_2_00_res$prev$median[1:501])
plot(rw_5000_1_00_res$prev$median[1:501])
plot(rw_5000_2_00_res$prev$median[1:501])


