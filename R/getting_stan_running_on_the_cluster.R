#######################################################
## Setting up work for the DIDE cluster ###############
#######################################################

didehpc::didehpc_config_global(didehpc.username = "DIDE\jd2117")
options(didehpc.username = "jd2117",didehpc.home = "Y:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config



context::context_log_start()

root <- "contexts"

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),sources = "simpleepp/R/loop_functions_random_walk.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

t <- obj$enqueue(make_true_epidemic())

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_100_complete_data_no_art")

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

params<-list(mu=mu,mu_i=mu_i,sigma=sigma)
params

sample_range<-1970:2015
sample_n<-500
penalty_order<-1
rows_to_evaluate<- 0:45*10+1   #(time_points_to_sample - 1970) * 10 + 1                 ## If using all data points must use 0:45*10+1

data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

n_500_RW_first_order_loop<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                            data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                            simulated_true_df = sim_model_output$sim_df))
obj$cluster_load()
obj$task_list()
test_fit_loop$id
didehpc::didehpc_config()
