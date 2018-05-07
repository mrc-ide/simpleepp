####################################################################################################################################
## Performing our bias test on the datasets to see how bias in the sample affects the outcome ######################################
####################################################################################################################################

require(ggplot2)
require(reshape2)

load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_100")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_500")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_1000")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/first_order/cluster_RW_first_order_n_5000")

load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_100")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_500")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_1000")
load("hiv_project/stan_objects_from_simpleepp_R/random_walk_loops/cluster_runs/second_order/rw_second_order_complete_data_n_5000")

load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_100")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_500")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_1000")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/first_order/first_order_complete_spline_n_5000")

load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_100")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_500")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_1000")
load("hiv_project/stan_objects_from_simpleepp_R/spline_runs/cluster_runs/second_order/spline_second_order_complete_n_5000")

load("hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_100_complete_data_no_art")
load("hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_500_complete_data_no_art")
load("hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_1000_complete_data_no_art")
load("hiv_project/simulated_data_sets/complete_data_simpleepp_no_art/n_5000_complete_data_no_art")

#################################################################################################################################
## Now we've loaded the data lets create our function to evaluate the mean absolute error at each point #########################
#################################################################################################################################

bias_function<-function(fitted_prev_df,sample_prev_df,true_df,time_to_evaluate=seq(1970,2015,0.1)){
  
  time_to_test<-seq((time_to_evaluate[1]-1970)*10+1,(time_to_evaluate[length(time_to_evaluate)]-2015)*10+501,1)
  
  time_to_sample<-seq((time_to_evaluate[1]-1970)+1,(time_to_evaluate[length(time_to_evaluate)]-2015)+46,1)
  
  
  for (i in 1:100){
    
    fit_df_test<-fitted_prev_df[fitted_prev_df$iteration == i,]
    sample_df<-sample_prev_df[sample]
    
    
  }
  
  
  
  
  
}




