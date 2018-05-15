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
  
  time_to_test<-seq((time_to_evaluate[1]-1970)*10+1,(time_to_evaluate[length(time_to_evaluate)]-2015)*10+451,1)
  
  time_to_sample<-seq((time_to_evaluate[1]-1970)+1,(time_to_evaluate[length(time_to_evaluate)]-2015)+46,1)
  
  sample_to_true_abs<- NULL
  sample_to_fitted_abs<- NULL
  
  sample_to_true_tot<- NULL
  sample_to_fitted_tot<- NULL
  
  rmse_sample_true<- NULL
  rmse_sample_fitted<- NULL
  
  for (i in 1:100){
    
    fit_df_test<-fitted_prev_df[fitted_prev_df$iteration == i,]
    sample_df<-sample_prev_df[sample_prev_df$iteration == i,]
    
    sample_to_true_diff<-sample_df$sample_prev_hiv_percentage - true_df$prev_percent[0:45*10+1]
    sample_to_fitted_diff<- sample_df$sample_prev_hiv_percentage - fit_df_test$median[0:45*10+1]
    
    sample_true_rmse<-sqrt(mean(sample_to_true_diff^2))
    sample_fitted_rmse<-sqrt(mean(sample_to_fitted_diff^2))
    
    tot_sample_to_true<-sum(sample_to_true_diff)
    tot_sample_to_fitted<-sum(sample_to_fitted_diff)
    
    abs_sample_to_true<-sum(abs(sample_to_true_diff))
    abs_sample_to_fitted<-sum(abs(sample_to_fitted_diff))
    
    sample_to_true_abs<-c(sample_to_true_abs,abs_sample_to_true)
    sample_to_fitted_abs<-c(sample_to_fitted_abs,abs_sample_to_fitted)
    
    sample_to_true_tot<-c(sample_to_true_tot,tot_sample_to_true)
    sample_to_fitted_tot<-c(sample_to_fitted_tot,tot_sample_to_fitted)
    
    rmse_sample_true<-c(rmse_sample_true,sample_true_rmse)
    rmse_sample_fitted<-c(rmse_sample_fitted,sample_fitted_rmse)
    
  }
  
  bias_df<-cbind.data.frame(sample_to_true_tot,sample_to_true_abs,rmse_sample_true,
                            sample_to_fitted_tot,sample_to_fitted_abs,rmse_sample_fitted)
  names(bias_df)<-c("total_sample_true","abs_sample_true", "rmse_sample_true",
                    "total_sample_fitted","abs_sample_fitted","rmse_sample_fitted")
  
  bias_df$diff_tot<-bias_df$total_sample_true - bias_df$total_sample_fitted
  bias_df$diff_abs<-bias_df$abs_sample_true - bias_df$abs_sample_fitted
  bias_df$diff_rmse<-bias_df$rmse_sample_true - bias_df$rmse_sample_fitted
  bias_df$iteration<-1:100
 
  melted_rmse_bias<-bias_df[,c(3,6,9,10)]
  
  melto_sample_bioas_df_rmse<-melt(melted_rmse_bias,id="iteration")
  
  rmse_bias_plot<-ggplot(data = melto_sample_bioas_df_rmse)+
    geom_line(aes(x=iteration,y=value,colour=variable),size=1)+
    labs(x="iteration",y="rmse between sample and dataset",title="RMSE differences to sample")
  
  melted_abs_bias<-bias_df[,c(2,5,8,10)]
  
  melto_abs_bias<-melt(melted_abs_bias , id="iteration")
  
  abs_bias_plot<- ggplot(data = melto_abs_bias)+geom_line(aes(x=iteration,y=value,colour=variable),size=1)+
    labs(x="iteration",y="abs_differences_to_sample_dataset",title="Abs differences to sample")
  
  melted_tot_bias<-bias_df[,c(1,4,7,10)]
  
  melto_tot_bias<-melt(melted_tot_bias,id="iteration")
  
  tot_bias_plot<-ggplot(data = melto_tot_bias)+geom_line(aes(x=iteration,y=value,colour=variable),size=1)+
    labs(x="iteration",y="tot_differences_to_sample_data",title="total differences from sample data ")
  
  mean_sample_true_tot<-mean(sample_to_true_tot)
  mean_sample_true_abs<-mean(sample_to_true_abs)
  mean_sample_true_rmse<-mean(rmse_sample_true)
  
  mean_sample_fitted_tot<-mean(sample_to_fitted_tot)
  mean_sample_fitted_abs<-mean(sample_to_fitted_abs)
  mean_sample_fitted_rmse<-mean(rmse_sample_fitted)
  
  mean_tot_diff<-mean(bias_df$diff_tot)
  mean_abs_diff<-mean(bias_df$diff_abs)
  mean_rmse_diff<-mean(bias_df$diff_rmse)
  
  mean_df<-cbind.data.frame(mean_sample_fitted_tot,mean_sample_fitted_abs,mean_sample_fitted_rmse,
                            mean_sample_true_tot,mean_sample_true_abs,mean_sample_true_rmse,
                            mean_tot_diff,mean_abs_diff,mean_rmse_diff)
  
  
  return(list(df=bias_df,rmse_plot=rmse_bias_plot,abs_plot=abs_bias_plot,tot_plot=tot_bias_plot,mean_df=mean_df))
  
  
}

bias_spline_f_100 <- bias_function(first_order_spline_n_100$prev,sample_prev_df = sampled_n_100_complete_data,
                                 true_df = sim_model_output$sim_df)
bias_spline_f_500 <- bias_function(first_order_spline_n_500$prev,sample_prev_df = sampled_n_500_complete_data,
                                 true_df = sim_model_output$sim_df)
bias_spline_f_1k <- bias_function(first_order_spline_n_1000$prev,sample_prev_df = sampled_n_1000_complete_data,
                                  true_df = sim_model_output$sim_df)
bias_spline_f_5k <- bias_function(first_order_spline_n_5000$prev,sample_prev_df = sampled_n_5000_complete_data,
                                  true_df = sim_model_output$sim_df)

######################################################################################################################################
## So thats the first order spline tested, now we will test the second order splines #################################################
######################################################################################################################################

bias_spline_s_100 <- bias_function(second_order_spline_n_100$prev,sample_prev_df = sampled_n_100_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_spline_s_500 <- bias_function(second_order_spline_n_500$prev,sample_prev_df = sampled_n_500_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_spline_s_1k <- bias_function(second_order_spline_n_1000$prev,sample_prev_df = sampled_n_1000_complete_data,
                                  true_df = sim_model_output$sim_df)
bias_spline_s_5k <- bias_function(second_order_spline_n_5000$prev,sample_prev_df = sampled_n_5000_complete_data,
                                  true_df = sim_model_output$sim_df)

######################################################################################################################################
## So thats's splines done, now we can look at the two randomw walks #################################################################
######################################################################################################################################
bias_rw_f_100 <- bias_function(RW_first_order_n_100$prev,sample_prev_df = sampled_n_100_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_rw_f_500 <- bias_function(RW_first_order_n_500$prev,sample_prev_df = sampled_n_500_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_rw_f_1k <- bias_function(RW_first_order_n_1000$prev,sample_prev_df = sampled_n_1000_complete_data,
                                  true_df = sim_model_output$sim_df)
bias_rw_f_5k <- bias_function(RW_first_order_n_5000$prev,sample_prev_df = sampled_n_5000_complete_data,
                                  true_df = sim_model_output$sim_df)

######################################################################################################################################
## So thats the first order rw tested, now we will test the second order rws #########################################################
######################################################################################################################################

bias_rw_s_100 <- bias_function(RW_second_order_n_100$prev,sample_prev_df = sampled_n_100_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_rw_s_500 <- bias_function(RW_second_order_n_500$prev,sample_prev_df = sampled_n_500_complete_data,
                                   true_df = sim_model_output$sim_df)
bias_rw_s_1k <- bias_function(RW_second_order_n_1000$prev,sample_prev_df = sampled_n_1000_complete_data,
                                  true_df = sim_model_output$sim_df)
bias_rw_s_5k <- bias_function(RW_second_order_n_5000$prev,sample_prev_df = sampled_n_5000_complete_data,
                                  true_df = sim_model_output$sim_df)

######################################################################################################################################
## Now we'll try implementing our measure of how overfit the fitted values are to the data sampled compared to the true values #######
######################################################################################################################################

overfit_function<-function(fit_df,sample_df,true_df){
  
  sample_true_rmse_tot<-NULL
  true_fitted_rmse_tot<-NULL
  overfit_stat_tot<-NULL
  overfit_ratio_tot<- NULL
  
  for (i in 1:100){
    sample_df_test<- sample_df[sample_df$iteration == i,]
    fitted_df_test<- fit_df[fit_df$iteration == i,]
    
    sample_true_error<- sample_df_test$sample_prev_hiv_percentage - true_df$prev_percent[0:45*10+1]
    true_fit_error<- true_df$prev_percent[0:45*10+1] - fitted_df_test$median[0:45*10+1]
    sample_fit_error<- sample_df_test$sample_prev_hiv_percentage - fitted_df_test$median[0:45*10+1]
    
    
    sample_true_rmse<-sqrt(mean(sample_true_error^2))
    true_fitted_rmse<-sqrt(mean(true_fit_error^2))
    sample_fitted_rmse<-sqrt(mean(sample_fit_error^2))
    
    overfit_ratio<- true_fitted_rmse / sample_fitted_rmse
    overfit_test<- (sample_true_rmse - true_fitted_rmse) / (sample_true_rmse)
    
    sample_true_rmse_tot<-c(sample_true_rmse_tot,sample_true_rmse)
    true_fitted_rmse_tot<-c(true_fitted_rmse_tot,true_fitted_rmse)
    overfit_stat_tot<-c(overfit_stat_tot,overfit_test)
    overfit_ratio_tot<-c(overfit_ratio_tot,overfit_ratio)
  }
  out_df<-cbind.data.frame(sample_true_rmse_tot,true_fitted_rmse_tot,overfit_stat_tot)
  names(out_df)<-c("sample true rmse","true fitted rmse", "overfit test")
  out_df$iteration<-1:100
  
  mean_overfit<-mean(overfit_stat_tot)
  mean_ratio<-mean(overfit_ratio_tot)
  
  return(list(df=out_df,overfit=mean_overfit,ovefit_ratio=mean_ratio))
  
}

overfit_spline_f_n_100<-overfit_function(first_order_spline_n_100$prev,true_df = sim_model_output$sim_df,
                                         sample_df = sampled_n_100_complete_data)
overfit_spline_f_n_500<-overfit_function(first_order_spline_n_500$prev,true_df = sim_model_output$sim_df,
                                         sample_df = sampled_n_500_complete_data)
overfit_spline_f_n_1K<-overfit_function(first_order_spline_n_1000$prev,true_df = sim_model_output$sim_df,
                                        sample_df = sampled_n_1000_complete_data)
overfit_spline_f_n_5K<-overfit_function(first_order_spline_n_5000$prev,true_df = sim_model_output$sim_df,
                                        sample_df = sampled_n_5000_complete_data)
overfit_rw_f_100<-overfit_function(RW_first_order_n_100$prev,true_df = sim_model_output$sim_df,
                                   sample_df = sampled_n_100_complete_data)
overfit_rw_f_500<-overfit_function(RW_first_order_n_500$prev,true_df = sim_model_output$sim_df,
                                   sample_df = sampled_n_500_complete_data)
overfit_rw_f_1k<-overfit_function(RW_first_order_n_1000$prev,true_df = sim_model_output$sim_df,
                                  sample_df = sampled_n_1000_complete_data)
overfit_rw_f_5k<-overfit_function(RW_first_order_n_5000$prev,true_df = sim_model_output$sim_df,
                                  sample_df = sampled_n_5000_complete_data)

overfit_spline_s_n_100<-overfit_function(second_order_spline_n_100$prev,true_df = sim_model_output$sim_df,
                                         sample_df = sampled_n_100_complete_data)
overfit_spline_s_n_500<-overfit_function(second_order_spline_n_500$prev,true_df = sim_model_output$sim_df,
                                         sample_df = sampled_n_500_complete_data)
overfit_spline_s_n_1K<-overfit_function(second_order_spline_n_1000$prev,true_df = sim_model_output$sim_df,
                                        sample_df = sampled_n_1000_complete_data)
overfit_spline_s_n_5K<-overfit_function(second_order_spline_n_5000$prev,true_df = sim_model_output$sim_df,
                                        sample_df = sampled_n_5000_complete_data)
overfit_rw_s_100<-overfit_function(RW_second_order_n_100$prev,true_df = sim_model_output$sim_df,
                                   sample_df = sampled_n_100_complete_data)
overfit_rw_s_500<-overfit_function(RW_second_order_n_500$prev,true_df = sim_model_output$sim_df,
                                   sample_df = sampled_n_500_complete_data)
overfit_rw_s_1k<-overfit_function(RW_second_order_n_1000$prev,true_df = sim_model_output$sim_df,
                                  sample_df = sampled_n_1000_complete_data)
overfit_rw_s_5k<-overfit_function(RW_second_order_n_5000$prev,true_df = sim_model_output$sim_df,
                                  sample_df = sampled_n_5000_complete_data)

spline_first_overfit<-c(overfit_spline_f_n_100$overfit,overfit_spline_f_n_500$overfit,
                        overfit_spline_f_n_1K$overfit,overfit_spline_f_n_5K$overfit)
spline_second_overfit<-c(overfit_spline_s_n_100$overfit,overfit_spline_f_n_500$overfit,
                         overfit_spline_s_n_1K$overfit,overfit_spline_s_n_5K$overfit)
rw_first_overfit<-c(overfit_rw_f_100$overfit,overfit_rw_f_500$overfit,
                    overfit_rw_f_1k$overfit,overfit_rw_f_5k$overfit)
rw_second_overfit<-c(overfit_rw_s_100$overfit,overfit_rw_s_500$overfit,
                     overfit_rw_s_1k$overfit,overfit_rw_s_5k$overfit)

overfit_satas<-cbind.data.frame(spline_first_overfit,spline_second_overfit,rw_first_overfit,rw_second_overfit)

spline_first_overfit_ratio<-c(overfit_spline_f_n_100$ovefit_ratio,overfit_spline_f_n_500$ovefit_ratio,
                        overfit_spline_f_n_1K$ovefit_ratio,overfit_spline_f_n_5K$ovefit_ratio)
spline_second_overfit_ratio<-c(overfit_spline_s_n_100$ovefit_ratio,overfit_spline_f_n_500$ovefit_ratio,
                         overfit_spline_s_n_1K$ovefit_ratio,overfit_spline_s_n_5K$ovefit_ratio)
rw_first_overfit_ratio<-c(overfit_rw_f_100$ovefit_ratio,overfit_rw_f_500$ovefit_ratio,
                    overfit_rw_f_1k$ovefit_ratio,overfit_rw_f_5k$ovefit_ratio)
rw_second_overfit_ratio<-c(overfit_rw_s_100$ovefit_ratio,overfit_rw_s_500$ovefit_ratio,
                     overfit_rw_s_1k$ovefit_ratio,overfit_rw_s_5k$ovefit_ratio)
overall_ratio<-cbind.data.frame(spline_first_overfit_ratio,spline_second_overfit_ratio,rw_first_overfit_ratio,rw_second_overfit_ratio)

overall_ratio

overfit_tests<-list(ratio=overall_ratio,initial_test_stat=overfit_satas)

save(overfit_tests,file = "hiv_project/analysis_of_cluster_run_datasets/no_art_simpleepp/overfit_analyses/overfit_list")
