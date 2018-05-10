######################################################################################################################################
## 95 % credible interval analysis ###################################################################################################
######################################################################################################################################

require(ggplot2)
require(reshape2)
require(plyr)

######################################################################################################################################
## Now lets load up the datasets for measuring whether they fall in the 95% confidence interval ######################################
######################################################################################################################################

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


spline_95_first_order_100<-credible_interval_test(sim_model_output$sim_df$prev_percent,first_order_spline_n_100$prev)
spline_95_first_order_500<-credible_interval_test(sim_model_output$sim_df$prev_percent,first_order_spline_n_500$prev)
spline_95_first_order_1000<-credible_interval_test(sim_model_output$sim_df$prev_percent,first_order_spline_n_1000$prev)
spline_95_first_order_5000<-credible_interval_test(sim_model_output$sim_df$prev_percent,first_order_spline_n_5000$prev)

spline_95_second_order_100<-credible_interval_test(sim_model_output$sim_df$prev_percent,second_order_spline_n_100$prev)
spline_95_second_order_500<-credible_interval_test(sim_model_output$sim_df$prev_percent,second_order_spline_n_500$prev)
spline_95_second_order_1000<-credible_interval_test(sim_model_output$sim_df$prev_percent,second_order_spline_n_1000$prev)
spline_95_second_order_5000<-credible_interval_test(sim_model_output$sim_df$prev_percent,second_order_spline_n_5000$prev)

RW_95_first_order_100<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_first_order_n_100$prev)
RW_95_first_order_500<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_first_order_n_500$prev)
RW_95_first_order_1000<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_first_order_n_1000$prev)
RW_95_first_order_5000<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_first_order_n_5000$prev)

RW_95_second_order_100<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_second_order_n_100$prev)
RW_95_second_order_500<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_second_order_n_500$prev)
RW_95_second_order_1000<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_second_order_n_1000$prev)
RW_95_second_order_5000<-credible_interval_test(sim_model_output$sim_df$prev_percent,RW_second_order_n_5000$prev)


splino_firsto_95<-c(spline_95_first_order_100$overall,spline_95_first_order_500$overall,
                    spline_95_first_order_1000$overall,spline_95_first_order_5000$overall)
splino_secondo_95<-c(spline_95_second_order_100$overall,spline_95_second_order_500$overall,
                     spline_95_second_order_1000$overall,spline_95_second_order_5000$overall)
rw_firsto_95<-c(RW_95_first_order_100$overall,RW_95_first_order_500$overall,
                RW_95_first_order_1000$overall,RW_95_first_order_5000$overall)
rw_secondo_95<-c(RW_95_second_order_100$overall,RW_95_second_order_500$overall,
                 RW_95_second_order_1000$overall,RW_95_second_order_5000$overall)

overall_stats<-cbind.data.frame(splino_firsto_95,splino_secondo_95,rw_firsto_95,rw_secondo_95)

#########################################################################################################################################
## So that's our test of whether prevalence is within the 95% confidence intervals or not, now we will test incidence ###################
#########################################################################################################################################

spline_95_first_order_100_inc<-credible_interval_test(sim_model_output$sim_df$lambda,first_order_spline_n_100$incidence)
spline_95_first_order_500_inc<-credible_interval_test(sim_model_output$sim_df$lambda,first_order_spline_n_500$incidence)
spline_95_first_order_1000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,first_order_spline_n_1000$incidence)
spline_95_first_order_5000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,first_order_spline_n_5000$incidence)

spline_95_second_order_100_inc<-credible_interval_test(sim_model_output$sim_df$lambda,second_order_spline_n_100$incidence)
spline_95_second_order_500_inc<-credible_interval_test(sim_model_output$sim_df$lambda,second_order_spline_n_500$incidence)
spline_95_second_order_1000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,second_order_spline_n_1000$incidence)
spline_95_second_order_5000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,second_order_spline_n_5000$incidence)

RW_95_first_order_100_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_first_order_n_100$incidence)
RW_95_first_order_500_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_first_order_n_500$incidence)
RW_95_first_order_1000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_first_order_n_1000$incidence)
RW_95_first_order_5000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_first_order_n_5000$incidence)

RW_95_second_order_100_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_second_order_n_100$incidence)
RW_95_second_order_500_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_second_order_n_500$incidence)
RW_95_second_order_1000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_second_order_n_1000$incidence)
RW_95_second_order_5000_inc<-credible_interval_test(sim_model_output$sim_df$lambda,RW_second_order_n_5000$incidence)


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

spline_95_first_order_100_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,first_order_spline_n_100$kappa)
spline_95_first_order_500_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,first_order_spline_n_500$kappa)
spline_95_first_order_1000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,first_order_spline_n_1000$kappa)
spline_95_first_order_5000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,first_order_spline_n_5000$kappa)

spline_95_second_order_100_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,second_order_spline_n_100$kappa)
spline_95_second_order_500_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,second_order_spline_n_500$kappa)
spline_95_second_order_1000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,second_order_spline_n_1000$kappa)
spline_95_second_order_5000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,second_order_spline_n_5000$kappa)

RW_95_first_order_100_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_first_order_n_100$kappa)
RW_95_first_order_500_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_first_order_n_500$kappa)
RW_95_first_order_1000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_first_order_n_1000$kappa)
RW_95_first_order_5000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_first_order_n_5000$kappa)

RW_95_second_order_100_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_second_order_n_100$kappa)
RW_95_second_order_500_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_second_order_n_500$kappa)
RW_95_second_order_1000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_second_order_n_1000$kappa)
RW_95_second_order_5000_kappa<-credible_interval_test(sim_model_output$sim_df$kappa,RW_second_order_n_5000$kappa)


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
