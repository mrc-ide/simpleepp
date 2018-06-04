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


mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

params<-list(mu=mu,mu_i=mu_i,sigma=sigma)
params



sample_range<-1984:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
sample_start<-sample_range[1]-1970
rows_to_evaluate<- sample_start:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_84_data,
                                            data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                            simulated_true_df = sim_model_output$sim_df), name = "data_84_rw_1st_n_100")


n_100_RW_first_order_loop$status()
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_100_FIRST_09_24_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_500_FIRST_09_24_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_1000_FIRST_09_26_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_5000_FIRST_09_26_JUNE_4")

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
n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(n_100_RW_first_order_loop_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_100_second_09_27_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_500_second_09_28_JUNE_4")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_1000_second_09_28_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_5000_second_09_29_JUNE_4")

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
spline_first_order_n_100<-spline_first_order_n_100$id
save(spline_first_order_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_100_first_09_29_JUNE_4")


sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_84_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_84_n_500")

spline_first_order_n_500$status()
spline_first_order_n_500<-spline_first_order_n_500$id
save(spline_first_order_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_500_first_09_31_JUNE_4")

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_84_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                      name = "spline_1st_84_n_1000")

spline_first_order_n_1000$status()
spline_first_order_n_1000<-spline_first_order_n_1000$id
save(spline_first_order_n_1000,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_1000_first_09_32_JUNE_4")

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,
                          sample_years=length(sample_range),sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_84_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df),
                                       name = "spline_1st_84_n_5000")

spline_first_order_n_5000$status()
spline_first_order_n_5000<-spline_first_order_n_5000$id
save(spline_first_order_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_1000_first_09_33_JUNE_4")



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

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_84_data,
                                                                        data_about_sampling = data_about_sampling,
                                                                        iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df),
                                       name = "Spline_2nd_84_100")
spline_second_order_n_100$status()
spline_second_order_n_100$log()
spline_second_order_n_100_id<-spline_second_order_n_100$id
save(spline_second_order_n_100_id,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_100_second_09_38_JUNE_4")


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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_100_second_09_48_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_1000_second_09_49_JUNE_4")

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
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_5000_second_09_50_JUNE_4")


#######################################################################################################################################
## Now lets load up the results #######################################################################################################
#######################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_100_FIRST_09_24_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_500_FIRST_09_24_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_1000_FIRST_09_26_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_5000_FIRST_09_26_JUNE_4",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_100_second_09_27_JUNE_4",
     verbose = T)
n_100_RW_sec_order_loop_id<-n_100_RW_first_order_loop_id
rm(n_100_RW_first_order_loop_id)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_100_FIRST_09_24_JUNE_4",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_500_second_09_28_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_1000_second_09_28_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/RW_5000_second_09_29_JUNE_4",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_100_first_09_29_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_500_first_09_31_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_1000_first_09_32_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_5000_first_09_33_JUNE_4",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_100_second_09_38_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_500_second_09_48_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_1000_second_09_49_JUNE_4",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/spline_5000_second_09_50_JUNE_4",
     verbose = T)

rw_100_1_84<-obj$task_get(n_100_RW_first_order_loop_id)
rw_500_1_84<-obj$task_get(n_500_RW_first_order_loop_id)
rw_1000_1_84<-obj$task_get(n_1000_RW_first_order_loop_id)
rw_5000_1_84<-obj$task_get(n_5000_RW_first_order_loop_id)

rw_100_2_84<-obj$task_get(n_100_RW_sec_order_loop_id)
rw_500_2_84<-obj$task_get(n_500_RW_sec_order_loop_id)
rw_1000_2_84<-obj$task_get(n_1000_RW_sec_order_loop_id)
rw_5000_2_84<-obj$task_get(n_5000_RW_sec_order_loop_id)

sp_100_1_84<-obj$task_get(spline_first_order_n_100)
sp_500_1_84<-obj$task_get(spline_first_order_n_500)
sp_1000_1_84<-obj$task_get(spline_first_order_n_1000)
sp_5000_1_84<-obj$task_get(spline_first_order_n_5000)

#### now for the results ####

rw_100_1_84_res<-rw_100_1_84$result()
rw_500_1_84_res<-rw_500_1_84$result()
rw_1000_1_84_res<-rw_1000_1_84$result()
rw_5000_1_84_res<-rw_5000_1_84$result()

rw_100_2_84_res<-rw_100_2_84$result()
rw_500_2_84_res<-rw_500_2_84$result()
rw_1000_2_84_res<-rw_1000_2_84$result()
rw_5000_2_84_res<-rw_5000_2_84$result()

sp_100_1_84_res<-sp_100_1_84$result()
sp_500_1_84_res<-sp_500_1_84$result()
sp_1000_1_84_res<-sp_1000_1_84$result()
sp_5000_1_84_res<-sp_5000_1_84$result()

sp_100_2_84_res<-sp_100_2_84$result()
sp_500_2_84_res<-sp_500_2_84$result()
sp_1000_2_84_res<-sp_1000_2_84$result()
sp_5000_2_84_res<-sp_5000_2_84$result()


save(rw_100_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_100_FIRSt_ORDER")
save(rw_500_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_500_FIRSt_ORDER")
save(rw_1000_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_1000_FIRSt_ORDER")
save(rw_5000_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_5000_FIRSt_ORDER")

save(rw_100_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_100_SECOND_ORDER")
save(rw_500_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_500_SECOND_ORDER")
save(rw_1000_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_1000_SECOND_ORDER")
save(rw_5000_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_5000_SECOND_ORDER")

save(sp_100_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_100_FIRSt_ORDER")
save(sp_500_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_500_FIRSt_ORDER")
save(sp_1000_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_1000_FIRSt_ORDER")
save(sp_5000_1_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_5000_FIRSt_ORDER")

save(sp_100_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_100_SECOND_ORDER")
save(sp_500_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_500_SECOND_ORDER")
save(sp_1000_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_1000_SECOND_ORDER")
save(sp_5000_2_84_res,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_5000_SECOND_ORDER")


spline_second_order_n_100$status()
spline_second_order_n_500$status()
spline_second_order_n_1000$status()
spline_second_order_n_5000$status()

second_order_spline_n_100<-spline_second_order_n_100$result()
second_order_spline_n_500<-spline_second_order_n_500$result()
second_order_spline_n_1000<-spline_second_order_n_1000$result()
second_order_spline_n_5000<-spline_second_order_n_5000$result()

save(second_order_spline_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_100")
save(second_order_spline_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_500")
save(second_order_spline_n_1000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_1000")
save(second_order_spline_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_5000")





save.image(file = "C:/Users/josh/Dropbox/.RData")


blab<-fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                 data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                 simulated_true_df = sim_model_output$sim_df)

rw_seond_n_100$status()
obj$task_status()



obj$unsubmit(n_1000_RW_first_order_loop$id)
obj$unsubmit(n_5000_RW_first_order_loop$id)


load("C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_100")



obj$cluster_load()
obj$task_list()
test_fit_loop$id
didehpc::didehpc_config()

n_500_id<-obj$task_list()[2]
n_100_id<-obj$task_list()[10]

n_500_RW_first_order<-obj$task_get(n_500_id)
n_100_RW_first_order<-obj$task_get(n_100_id)

n_100_RW_first_order$status()
obj$task_list()
obj$task_status(n_500_id)

data_about_sampling$knot_number<-7
data_about_sampling$sample_n<-100

a<-fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                              data_about_sampling = data_about_sampling,iteration_number = 1,params = params,
                              simulated_true_df = sim_model_output$sim_df)



n_100_RW_first_order_results<-n_100_RW_first_order$result()
