#######################################################
## Setting up work for the DIDE cluster ###############
#######################################################

setwd("Y:")
options(didehpc.username = "jd2117",didehpc.home = "Y:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##


context::context_log_start()

root <- "contexts"

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),sources = "simpleepp/R/loop_functions_random_walk.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

################################################################################################################################
## So the above lines of code give me access to the cluster, with the obj object giving me queue functionalaity, I've also #####
## asked for three cores for stan to use for each of the chains ################################################################
################################################################################################################################



t <- obj$enqueue(make_true_epidemic())

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_1000_complete_data_no_art")

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

params<-list(mu=mu,mu_i=mu_i,sigma=sigma)
params

sample_range<-1970:2015
sample_n<-1000                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
rows_to_evaluate<- 0:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_1000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                            data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                            simulated_true_df = sim_model_output$sim_df))


load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_5000_complete_data_no_art")

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_5000_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                   data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df))


sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_100_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                   data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                   simulated_true_df = sim_model_output$sim_df))

sample_n<-500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                  data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df))

n_100_RW_first_order_loop$status()
n_500_RW_first_order_loop$status()
n_1000_RW_first_order_loop$status()
n_5000_RW_first_order_loop$status()


for (i in 1:100){
  RW_first_order_n_100$prev$iteration[1:502+((i-1)*502)]<-i
}


RW_first_order_n_100<-n_100_RW_first_order_loop$result()
RW_first_order_n_1000<-n_1000_RW_first_order_loop$result()
RW_first_order_n_500<-n_500_RW_first_order_loop$result()
RW_first_order_n_5000<-n_5000_RW_first_order_loop$result()

save(RW_first_order_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_RW_first_order_n_100")

save(RW_first_order_n_1000,
     file="C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_1000")



obj$unsubmit(n_1000_RW_first_order_loop$id)
obj$unsubmit(n_5000_RW_first_order_loop$id)






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

a<-fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                              data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                              simulated_true_df = sim_model_output$sim_df)



n_100_RW_first_order_results<-n_100_RW_first_order$result()
