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

t <- obj$enqueue(make_true_epidemic())

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_100_FIRST_12_16_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_500_FIRST_12_17_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_1000_FIRST_12_18_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_5000_FIRST_12_18_JUNE_11")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_90_data,
                                                                   data_about_sampling = data_about_sampling,
                                                                   iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df),
                                        name = "data_90_rw_2nd_n_100")

n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
   file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_100_second_12_19_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_500_second_12_19_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_1000_second_12_19_JUNE_11")

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
   file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_5000_second_12_20_JUNE_11")

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
                                  simulated_true_df = sim_model_output$sim_df),name = "spline_90_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100<-spline_first_order_n_100$id
save(spline_first_order_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_100_first_12_20_JUNE_11")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_90_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_90_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500<-spline_first_order_n_500$id
save(spline_first_order_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_500_first_12_20_JUNE_11")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_90_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_90_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000<-spline_first_order_n_1000$id
save(spline_first_order_n_1000,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_1000_first_12_21_JUNE_11")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_90_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_90_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000<-spline_first_order_n_5000$id
save(spline_first_order_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_5000_first_12_22_JUNE_11")



spline_first_order_n_100$status()
spline_first_order_n_500$status()
spline_first_order_n_1000$status()
spline_first_order_n_5000$status()

first_order_spline_n_100<-spline_first_order_n_100$result()
first_order_spline_n_500<-spline_first_order_n_500$result()
first_order_spline_n_1000<-spline_first_order_n_1000$result()
first_order_spline_n_5000<-spline_first_order_n_5000$result()

save(first_order_spline_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_100")
save(first_order_spline_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_500")
save(first_order_spline_n_1000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_1000")
save(first_order_spline_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_5000")


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
spline_second_order_n_100$log()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_100_second_12_23_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_500_second_12_23_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_1000_second_12_24_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_5000_second_12_24_JUNE_11")


#######################################################################################################################################
## Now lets load up the results #######################################################################################################
#######################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_100_FIRST_12_16_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_500_FIRST_12_17_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_1000_FIRST_12_18_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_5000_FIRST_12_18_JUNE_11",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_100_second_12_19_JUNE_11",
     verbose = T)
n_100_RW_sec_order_loop_id<-n_100_RW_first_order_loop_id
rm(n_100_RW_first_order_loop_id)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_100_FIRST_12_16_JUNE_11",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_500_second_12_19_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_1000_second_12_19_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/RW_5000_second_12_20_JUNE_11",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_100_first_12_20_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_500_first_12_20_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_1000_first_12_21_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_5000_first_12_22_JUNE_11",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_100_second_12_23_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_500_second_12_23_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_1000_second_12_24_JUNE_11",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/spline_5000_second_12_24_JUNE_11",
     verbose = T)

rw_100_1_90<-obj$task_get(n_100_RW_first_order_loop_id)
rw_500_1_90<-obj$task_get(n_500_RW_first_order_loop_id)
rw_1000_1_90<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_5000_1_90<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_100_2_90<-obj$task_get(n_100_RW_sec_order_loop_id)
rw_500_2_90<-obj$task_get(n_500_RW_sec_order_loop_id)
rw_1000_2_90<-obj$task_get(n_1000_RW_sec_order_loop_id)
rw_5000_2_90<-obj$task_get(n_5000_RW_sec_order_loop_id)

sp_100_1_90<-obj$task_get(spline_first_order_n_100)
sp_500_1_90<-obj$task_get(spline_first_order_n_500)
sp_1000_1_90<-obj$task_get(spline_first_order_n_1000)
sp_5000_1_90<-obj$task_get(spline_first_order_n_5000)

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


save(rw_100_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_100_FIRSt_ORDER")
save(rw_500_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_500_FIRSt_ORDER")
save(rw_1000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_1000_FIRSt_ORDER")
save(rw_5000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_5000_FIRSt_ORDER")

save(rw_100_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_100_SECOND_ORDER")
save(rw_500_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_500_SECOND_ORDER")
save(rw_1000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_1000_SECOND_ORDER")
save(rw_5000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/RW_5000_SECOND_ORDER")

save(sp_100_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_100_FIRSt_ORDER")
save(sp_500_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_500_FIRSt_ORDER")
save(sp_1000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_1000_FIRSt_ORDER")
save(sp_5000_1_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_5000_FIRSt_ORDER")

save(sp_100_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_100_SECOND_ORDER")
save(sp_500_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_500_SECOND_ORDER")
save(sp_1000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_1000_SECOND_ORDER")
save(sp_5000_2_90_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/SPLINE_5000_SECOND_ORDER")

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
                                                                  data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df), name = "data_95_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_100_FIRST_14_42_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_500_FIRST_14_43_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_1000_FIRST_14_43_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_5000_FIRST_14_44_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_100_second_14_45_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_500_second_14_45_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_1000_second_14_45_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/RW_5000_second_14_46_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_95_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),name = "spline_95_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_100_first_14_47_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_500_first_14_47_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_1000_first_14_47_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_5000_first_14_48_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_100_second_14_56_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_500_second_15_31_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_500_second_15_32_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/spline_5000_second_15_33_JUNE_11")


###### lets load up the results for these cluster runs ###########

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/"
ids_1995<-list.files(path_name,include.dirs = F,recursive = T)
ids_1995<-paste(path_name,ids_1995,sep = "")
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

sp_5000_2_95_res$data_about_sampling

save(rw_100_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_100_FIRSt_ORDER")
save(rw_500_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_500_FIRSt_ORDER")
save(rw_1000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_1000_FIRSt_ORDER")
save(rw_5000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_5000_FIRSt_ORDER")

save(rw_100_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_100_SECOND_ORDER")
save(rw_500_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_500_SECOND_ORDER")
save(rw_1000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_1000_SECOND_ORDER")
save(rw_5000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/RW_5000_SECOND_ORDER")

save(sp_100_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_100_FIRSt_ORDER")
save(sp_500_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_500_FIRSt_ORDER")
save(sp_1000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_1000_FIRSt_ORDER")
save(sp_5000_1_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_5000_FIRSt_ORDER")

save(sp_100_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_100_SECOND_ORDER")
save(sp_500_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_500_SECOND_ORDER")
save(sp_1000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_1000_SECOND_ORDER")
save(sp_5000_2_95_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/SPLINE_5000_SECOND_ORDER")


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
                                                                  data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df), name = "data_00_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_100_FIRST_16_04_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_500_FIRST_16_05_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_1000_FIRST_16_05_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_5000_FIRST_16_05_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_100_second_16_06_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_500_second_16_06_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_1000_second_16_06_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/RW_5000_second_16_06_JUNE_11")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_00_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),name = "spline_00_n_100_first")

spline_first_order_n_100$status()
spline_first_order_n_100_id<-spline_first_order_n_100$id
save(spline_first_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_100_first_16_07_JUNE_11")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_500_first_16_07_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_1000_first_16_07_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_5000_first_16_08_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_100_second_16_08_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_500_second_16_08_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_500_second_16_09_JUNE_11")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/spline_5000_second_16_09_JUNE_11")


###### lets load up the results for these cluster runs ###########

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/"
ids_2000<-list.files(path_name,include.dirs = F,recursive = T)
ids_2000<-paste(path_name,ids_2000,sep = "")
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


