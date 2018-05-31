######################################################################################################################################
## 95 % credible interval analysis ###################################################################################################
######################################################################################################################################

require(ggplot2)
require(reshape2)
require(plyr)

######################################################################################################################################
## Now lets load up the datasets for measuring whether they fall in the 95% confidence interval ######################################
######################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_100",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_500",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_1000",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_FIRST_ORDER_N_5000",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_100",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_500",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_1000",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/RW_SECOND_ORDER_N_5000",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_100",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_500",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_1000",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_FIRST_ORDER_N_5000",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_100",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_500",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_1000",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/results/SPLINE_SECOND_ORDER_N_5000",verbose = T)

######################################################################################################################################
## Now we will write the function to test what percentage of time the true value lies within the lower and upper bounds ##############
######################################################################################################################################


credible_interval_test<-function(true_df_vector,fitted_df,time_points=seq(1970,2020,0.1)){
  
  time_to_test<-seq((time_points[1]-1970)*10+1,(time_points[length(time_points)]-2020)*10+501,1)
  percent_each_iter<-NULL
  
  
  for(i in 1:100){
    
    iteration_fitted<-fitted_df[fitted_df$iteration==i,]
    tot_iter_presence<-NULL
    
    iter_presence<-ifelse(true_df_vector[time_to_test]>=iteration_fitted$low[time_to_test] & 
                          true_df_vector[time_to_test] <= iteration_fitted$high[time_to_test],"y","n")
  
    tot_y<-plyr::count(iter_presence)
    tot_x<-tot_y[tot_y$x=="y",]
    
    
    
    percent_y<-(tot_x[1,2]/sum(tot_y[,2]))*100
    
    percent_each_iter<-cbind(percent_each_iter,percent_y)
    
  }
  
  overall_percent<-mean(percent_each_iter)
  
  return(list(overall=overall_percent,each_iter=percent_each_iter))
}


spline_95_first_order_100<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_100_res$prev)
spline_95_first_order_500<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_500_res$prev)
spline_95_first_order_1000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_1000_res$prev)
spline_95_first_order_5000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_5000_res$prev)

spline_95_second_order_100<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_100_res$prev)
spline_95_second_order_500<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_500_res$prev)
spline_95_second_order_1000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_1000_res$prev)
spline_95_second_order_5000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,spline_first_5000_res$prev)

RW_95_first_order_100<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_first_100_res$prev)
RW_95_first_order_500<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_first_500_res$prev)
RW_95_first_order_1000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_first_1000_res$prev)
RW_95_first_order_5000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_first_5000_res$prev)

RW_95_second_order_100<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_sec_100_res$prev)
RW_95_second_order_500<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_sec_500_res$prev)
RW_95_second_order_1000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_sec_1000_res$prev)
RW_95_second_order_5000<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$prev_percent,rw_sec_5000_res$prev)


splino_firsto_95<-c(spline_95_first_order_100$overall,spline_95_first_order_500$overall,
                    spline_95_first_order_1000$overall,spline_95_first_order_5000$overall)
splino_secondo_95<-c(spline_95_second_order_100$overall,spline_95_second_order_500$overall,
                     spline_95_second_order_1000$overall,spline_95_second_order_5000$overall)
rw_firsto_95<-c(RW_95_first_order_100$overall,RW_95_first_order_500$overall,
                RW_95_first_order_1000$overall,RW_95_first_order_5000$overall)
rw_secondo_95<-c(RW_95_second_order_100$overall,RW_95_second_order_500$overall,
                 RW_95_second_order_1000$overall,RW_95_second_order_5000$overall)

overall_stats<-cbind.data.frame(splino_firsto_95,splino_secondo_95,rw_firsto_95,rw_secondo_95)
overall_stats
#########################################################################################################################################
## So that's our test of whether prevalence is within the 95% confidence intervals or not, now we will test incidence ###################
#########################################################################################################################################

spline_95_first_order_100_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_100_res$incidence)
spline_95_first_order_500_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_500_res$incidence)
spline_95_first_order_1000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_1000_res$incidence)
spline_95_first_order_5000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_5000_res$incidence)

spline_95_second_order_100_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_100_res$incidence)
spline_95_second_order_500_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_500_res$incidence)
spline_95_second_order_1000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_1000_res$incidence)
spline_95_second_order_5000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,spline_first_5000_res$incidence)

RW_95_first_order_100_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_first_100_res$incidence)
RW_95_first_order_500_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_first_500_res$incidence)
RW_95_first_order_1000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_first_1000_res$incidence)
RW_95_first_order_5000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_first_5000_res$incidence)

RW_95_second_order_100_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_sec_100_res$incidence)
RW_95_second_order_500_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_sec_500_res$incidence)
RW_95_second_order_1000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_sec_1000_res$incidence)
RW_95_second_order_5000_inc<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$lambda,rw_sec_5000_res$incidence)


splino_firsto_95_inc<-c(spline_95_first_order_100_inc$overall,spline_95_first_order_500_inc$overall,
                    spline_95_first_order_1000_inc$overall,spline_95_first_order_5000_inc$overall)
splino_secondo_95_inc<-c(spline_95_second_order_100_inc$overall,spline_95_second_order_500_inc$overall,
                     spline_95_second_order_1000_inc$overall,spline_95_second_order_5000_inc$overall)
rw_firsto_95_inc<-c(RW_95_first_order_100_inc$overall,RW_95_first_order_500_inc$overall,
                RW_95_first_order_1000_inc$overall,RW_95_first_order_5000_inc$overall)
rw_secondo_95_inc<-c(RW_95_second_order_100_inc$overall,RW_95_second_order_500_inc$overall,
                 RW_95_second_order_1000_inc$overall,RW_95_second_order_5000_inc$overall)

overall_stats_inc<-cbind.data.frame(splino_firsto_95_inc,splino_secondo_95_inc,rw_firsto_95_inc,rw_secondo_95_inc)

overall_stats_inc

##########################################################################################################################################
## Now we will test for the kappa term, how this is kept within the 95% confidence intervals #############################################
##########################################################################################################################################

spline_95_first_order_100_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_100_res$kappa)
spline_95_first_order_500_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_500_res$kappa)
spline_95_first_order_1000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_1000_res$kappa)
spline_95_first_order_5000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_5000_res$kappa)

spline_95_second_order_100_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_100_res$kappa)
spline_95_second_order_500_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_500_res$kappa)
spline_95_second_order_1000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_1000_res$kappa)
spline_95_second_order_5000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,spline_first_5000_res$kappa)

RW_95_first_order_100_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_first_100_res$kappa)
RW_95_first_order_500_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_first_500_res$kappa)
RW_95_first_order_1000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_first_1000_res$kappa)
RW_95_first_order_5000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_first_5000_res$kappa)

RW_95_second_order_100_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_sec_100_res$kappa)
RW_95_second_order_500_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_sec_500_res$kappa)
RW_95_second_order_1000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_sec_1000_res$kappa)
RW_95_second_order_5000_kappa<-credible_interval_test(sim_model_output_changed_to_bell_curve$sim_df$kappa,rw_sec_5000_res$kappa)


splino_firsto_95_kappa<-c(spline_95_first_order_100_kappa$overall,spline_95_first_order_500_kappa$overall,
                        spline_95_first_order_1000_kappa$overall,spline_95_first_order_5000_kappa$overall)
splino_secondo_95_kappa<-c(spline_95_second_order_100_kappa$overall,spline_95_second_order_500_kappa$overall,
                         spline_95_second_order_1000_kappa$overall,spline_95_second_order_5000_kappa$overall)
rw_firsto_95_kappa<-c(RW_95_first_order_100_kappa$overall,RW_95_first_order_500_kappa$overall,
                    RW_95_first_order_1000_kappa$overall,RW_95_first_order_5000_kappa$overall)
rw_secondo_95_kappa<-c(RW_95_second_order_100_kappa$overall,RW_95_second_order_500_kappa$overall,
                     RW_95_second_order_1000_kappa$overall,RW_95_second_order_5000_kappa$overall)

overall_stats_kappa<-cbind.data.frame(splino_firsto_95_kappa,splino_secondo_95_kappa,rw_firsto_95_kappa,rw_secondo_95_kappa)

overall_stats_kappa


overall_95_stats<-list(prev=overall_stats,inc=overall_stats_inc,kappa=overall_stats_kappa)
save(overall_95_stats,file = "hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/95_percent_conf_analysis/overall_95")




