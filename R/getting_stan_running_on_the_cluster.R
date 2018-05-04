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
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts"

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

load("C:/Users/josh/Dropbox/hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_1000_complete_data_no_art")

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001

params<-list(mu=mu,mu_i=mu_i,sigma=sigma)
params

sample_range<-1970:2015
sample_n<-100                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
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
n_500_RW_first_order_loop<-obj$task_get(obj$task_list()[7])
n_1000_RW_first_order_loop$status()
n_5000_RW_first_order_loop<-obj$task_get(obj$task_list()[5])

n_500_or_5k<-obj$task_get(obj$task_list()[6])
n_5k_o_500<-obj$task_get(obj$task_list()[4])

n_500_or_5k$status()
n_5k_o_500$status()

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

save(RW_first_order_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_5000")

save(RW_first_order_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_500")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## That is the sending off of the first order jobs to the cluster, now we will send off the second order RW jobs %%%%%%%%%%%%%##
##%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&%%%%%%%%%%%%%%&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%*************%%%%%%%%&&&&&&&^^^^^^££££££££££££££##

penalty_order<-2
sample_n<-100
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

rw_second_n_100<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_100_complete_data,
                                                       data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                       simulated_true_df = sim_model_output$sim_df))


sample_n <- 500
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

rw_second_n_500<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_500_complete_data,
                                                       data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                       simulated_true_df = sim_model_output$sim_df))

sample_n<-1000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

rw_second_n_1000<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                        simulated_true_df = sim_model_output$sim_df))

sample_n<-5000
data_about_sampling<-list(penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

rw_second_n_5000<-obj$enqueue(fitting_data_function_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                        simulated_true_df = sim_model_output$sim_df))



rw_second_n_100$status()
rw_second_n_500$status()
rw_second_n_1000$status()
rw_second_n_5000$status()

RW_second_order_n_100<-rw_second_n_100$result()
RW_second_order_n_1000<-rw_second_n_1000$result()
RW_second_order_n_500<-rw_second_n_500$result()
RW_second_order_n_5000<-rw_second_n_5000$result()

save(RW_second_order_n_5000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_5000")
save(RW_second_order_n_100,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_100")
save(RW_second_order_n_1000,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_1000")
save(RW_second_order_n_500,
     file = "C:/Users/josh/Dropbox/hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_500")

##***************************************************************************************************************************##
## Now we will run the spline functions on the cluster ########################################################################
###############################################################################################################################

knot_number <- 7
penalty_order <- 1
sample_n<-100
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)


spline_first_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                  data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                  simulated_true_df = sim_model_output$sim_df))

sample_n <- 500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_first_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df))

sample_n <- 1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_first_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df))

sample_n <- 5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)
spline_first_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df))



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
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_100<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_100_complete_data,
                                                                        data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                        simulated_true_df = sim_model_output$sim_df))

sample_n<-500
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_500<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_500_complete_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df))

sample_n<-1000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_1000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_1000_complete_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df))

sample_n<-5000
data_about_sampling<-list(knot_number = knot_number,penalty_order=penalty_order,sample_years=46,sample_n=sample_n,rows_to_evaluate=rows_to_evaluate)

spline_second_order_n_5000<-obj$enqueue(fitting_data_function_spline_loop(samples_data_frame = sampled_n_5000_complete_data,
                                                                         data_about_sampling = data_about_sampling,iteration_number = 100,params = params,
                                                                         simulated_true_df = sim_model_output$sim_df))


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
