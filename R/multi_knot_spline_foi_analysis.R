###########################################################################################################################################
## Analysisng the multi-run knot values on the cluster values #############################################################################
###########################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)

############################################################################################################################################
## Now we will load up the datasets ########################################################################################################
############################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/true_epidemic",verbose = T)

## 7 knotters

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/"
seven_knots<-list.files(path_name,full.names = T)
lapply(seven_knots,load,verbose=T)

## 8 knotters 

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/"
eight_knots<-list.files(path_name,full.names = T)
lapply(eight_knots,load,verbose=T)

## 9 knotters 

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/"
nine_knots<-list.files(path_name,full.names = T)
lapply(nine_knots,load,verbose=T)

## 10 knotters 
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/"
ten_knots<-list.files(path_name,full.names = T)
lapply(ten_knots,load,verbose=T)

## 11 knotters

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/"
eleven_knots<-list.files(path_name,full.names = T)
lapply(eleven_knots,load,verbose=T)

## 12 knotters

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/"
twelve_knots<-list.files(path_name,full.names = T)
a<-lapply(twelve_knots,load,verbose=T)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


root_mean_error_function<-function(true_data,fitted_data,metric="prevalence",time_period=seq(1970,2020,0.1)){
  
  if(metric=="prevalence"){
    true_metric <- true_data$prevalence
    fitted_metric <- fitted_data$prev
    
  }
  
  if(metric=="incidence"){
    
    true_metric <- true_data$lambda
    fitted_metric <- fitted_data$incidence
  }
  
  if(metric=="kappa"){
    true_metric <- true_data$kappa
    fitted_metric <- fitted_data$kappa
  }
  
  
  time_to_test<-seq((time_period[1]-1970)*10+1,(time_period[length(time_period)]-2020)*10+501,1)
  
  mean_rmse_tot<-NULL
  error_tot<-NULL
  
  for (i in 1:100){
    
    
    fitted_metric_iter <- fitted_metric[fitted_metric$iteration == i,]
    
    if(metric=="prevalence"){
      fitted_metric_iter<-fitted_metric_iter/100
    }
    
    
    error <- (fitted_metric_iter$median[time_to_test]) - true_metric[time_to_test]
    
    rmse <- sqrt(mean(error^2))
    
    iter <- i
    
    mean_rmse <- cbind(rmse,iter)
    
    error_tot <- rbind(error,error_tot)
    
    mean_rmse_tot <- rbind(mean_rmse_tot,mean_rmse)
  }
  
  
  mean_overall_rmse <- mean(mean_rmse_tot[,1])
  
  
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot,error_df=error_tot))
  
  
  
}

###########################################################################################################################################
## Now we will run through our results to get the rmse values for the different knot values ###############################################
###########################################################################################################################################

## 7 knotters ##
seven_knot_results<-list(spline_f_100_foi_res,spline_f_500_foi_res,spline_f_1k_foi_res,spline_f_5k_foi_res,
                         spline_s_100_foi_res,spline_s_500_foi_res,spline_s_1k_foi_res,spline_s_5k_foi_res)
seven_knot_prev_rmse<-lapply(seven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
seven_knot_inc_rmse<-lapply(seven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
seven_knot_kappa_rmse<-lapply(seven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")

## 8 knotters ##
eight_knot_results<-list(spline_f_100_foi_res_8,spline_f_500_foi_res_8,spline_f_1k_foi_res_8,spline_f_5k_foi_res_8,
                         spline_s_100_foi_res_8,spline_s_500_foi_res_8,spline_s_1k_foi_res_8,spline_s_5k_foi_res_8)
eight_knot_prev_rmse<-lapply(eight_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
eight_knot_inc_rmse<-lapply(eight_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
eight_knot_kappa_rmse<-lapply(eight_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")


## 9 knotters ##
nine_knot_results<-list(spline_f_100_foi_res_9,spline_f_500_foi_res_9,spline_f_1k_foi_res_9,spline_f_5k_foi_res_9,
                        spline_s_100_foi_res_9,spline_s_500_foi_res_9,spline_s_1k_foi_res_9,spline_s_5k_foi_res_9)
nine_knot_prev_rmse<-lapply(nine_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
nine_knot_inc_rmse<-lapply(nine_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
nine_knot_kappa_rmse<-lapply(nine_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")

## 10 knotters ##
ten_knot_results<-list(spline_f_100_foi_res_ten,spline_f_500_foi_res_ten,spline_f_1k_foi_res_ten,spline_f_5k_foi_res_ten,
                       spline_s_100_foi_res_ten,spline_s_500_foi_res_ten,spline_s_1k_foi_res_ten,spline_s_5k_foi_res_ten)
ten_knot_prev_rmse<-lapply(ten_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
ten_knot_inc_rmse<-lapply(ten_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
ten_knot_kappa_rmse<-lapply(ten_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")

## 11 knotters ##
eleven_knot_results<-list(spline_f_100_foi_res_11,spline_f_500_foi_res_11,spline_f_1k_foi_res_11,spline_f_5k_foi_res_11,
                          spline_s_100_foi_res_11,spline_s_500_foi_res_11,spline_s_1k_foi_res_11,spline_s_5k_foi_res_11)
eleven_knot_prev_rmse<-lapply(eleven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
eleven_knot_inc_rmse<-lapply(eleven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
eleven_knot_kappa_rmse<-lapply(eleven_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")

## 12 knotters ##
twelve_knot_results<-list(spline_f_100_foi_res_12,spline_f_500_foi_res_12,spline_f_1k_foi_res_12,spline_f_5k_foi_res_12,
                          spline_s_100_foi_res_12,spline_s_500_foi_res_12,spline_s_1k_foi_res_12,spline_s_5k_foi_res_12)
twelve_knot_prev_rmse<-lapply(twelve_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df)
twelve_knot_inc_rmse<-lapply(twelve_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="incidence")
twelve_knot_kappa_rmse<-lapply(twelve_knot_results,root_mean_error_function,true_data=sim_model_foi$sim_df,metric="kappa")

### tot values ###

tot_prev_rmse<-list(seven_knot_prev_rmse,eight_knot_prev_rmse,nine_knot_prev_rmse,
                    ten_knot_prev_rmse,eleven_knot_prev_rmse,twelve_knot_prev_rmse)
tot_inc_rmse<-list(seven_knot_inc_rmse,eight_knot_inc_rmse,nine_knot_inc_rmse,
                   ten_knot_inc_rmse,eleven_knot_inc_rmse,twelve_knot_inc_rmse)
tot_kappa_rmse<-list(seven_knot_kappa_rmse,eight_knot_kappa_rmse,nine_knot_kappa_rmse,
                     ten_knot_kappa_rmse,eleven_knot_kappa_rmse,twelve_knot_kappa_rmse)

############################################################################################################################################
## Now lets create a function to take these outputs and give us a graph and a matrix/dataframe output ######################################
############################################################################################################################################

rmse_plotter_and_table_giver<-function(total_list_of_results,start_knot_number=7,plot_title){
  knot_values<-as.character(start_knot_number:(start_knot_number+length(total_list_of_results)-1))
  tot_vals<-NULL
  tot_plot_vals<-NULL
  for(i in 1:length(total_list_of_results)){
    new_line<-c(total_list_of_results[[i]][[1]][[1]],
                total_list_of_results[[i]][[2]][[1]],
                total_list_of_results[[i]][[3]][[1]],
                total_list_of_results[[i]][[4]][[1]],
                total_list_of_results[[i]][[5]][[1]],
                total_list_of_results[[i]][[6]][[1]],
                total_list_of_results[[i]][[7]][[1]],
                total_list_of_results[[i]][[8]][[1]])
    new_df<-cbind.data.frame(new_line,rep(knot_values[i],length(new_line)))
    tot_vals<-cbind(tot_vals,new_line)
    tot_plot_vals<-rbind(tot_plot_vals,new_df)
    
  }
  tot_vals<-data.frame(tot_vals)
  row.names(tot_vals)<-paste(c(rep("first",4),rep("second",4)),rep(c(100,500,1000,5000),2),sep = " ")
  
  tot_plot_vals<-data.frame(tot_plot_vals)
  names(tot_vals)<-paste("Knot_number",knot_values,sep = "_")
  names(tot_plot_vals)<-c("rmse","knot_number")
  tot_plot_vals$sample_size<-rep(c(100,500,1000,5000),nrow(tot_plot_vals)/4)
  tot_plot_vals$order<-rep(c(rep("first",4),rep("second",4)),nrow(tot_plot_vals)/8)
  tot_plot_vals$type<-paste(tot_plot_vals$knot_number,tot_plot_vals$order)
  
  rmse_plotto<-ggplot(data = tot_plot_vals,aes(x=sample_size,y=rmse,group=type))+
    geom_line(aes(colour=knot_number,linetype=order),size=1.05)+labs(x="Sample size",y="RMSE to true",title=plot_title)+
    geom_point(aes(shape=order,fill=knot_number,colour=knot_number),fill="white",size=3)+
    scale_shape_manual("Order",values = c("first"=21,"second"=23))+
    scale_linetype_manual("Order",
                          values = c("first"="solid","second"="dotdash"))
  
  
  return(list(df_rmse=tot_vals,plot_df=tot_plot_vals,plotto=rmse_plotto))
  
  
}

prevalence_vals<-rmse_plotter_and_table_giver(tot_prev_rmse,plot_title = "Prevalence RMSE across 6 different knot sizes")
prevalence_vals$plotto
prevalence_vals$df_rmse

incidence_vals<-rmse_plotter_and_table_giver(tot_inc_rmse,plot_title = "Incidence RMSE across 6 different knot sizes")
incidence_vals$plotto
incidence_vals$df_rmse

kappa_vals<-rmse_plotter_and_table_giver(tot_kappa_rmse,plot_title = "Kappa RMSE across 6 different knot sizes")
kappa_vals$plotto
kappa_vals$df_rmse

save(prevalence_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/rmse_analysis/complete_RMSE_for_multi_knot_splines")
save(incidence_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/rmse_analysis/incidence_RMSE_multi_knot_splines")
save(kappa_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/rmse_analysis/kappa_RMSE_multi_knot_splines")

###########################################################################################################################################
## Now we will prepare the 95 % confidence interval analysis tables and plots #############################################################
###########################################################################################################################################

credible_interval_test_multi_knotters<-function(true_df_vector,fitted_list,time_points=seq(1970,2020,0.1),metric="prevalence"){
  
  time_to_test<-seq((time_points[1]-1970)*10+1,(time_points[length(time_points)]-2020)*10+501,1)
  percent_each_iter<-NULL
  
  if(metric=="prevalence"){
    fitted_df<-fitted_list$prev
  }
  if(metric=="incidence"){
    fitted_df<-fitted_list$incidence
  }
  if(metric=="kappa"){
    fitted_df<-fitted_list$kappa
  }
  
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

plotter_and_table_giver<-function(total_list_of_results,start_knot_number=7,plot_title){
  knot_values<-as.character(start_knot_number:(start_knot_number+length(total_list_of_results)-1))
  tot_vals<-NULL
  tot_plot_vals<-NULL
  for(i in 1:length(total_list_of_results)){
    new_line<-c(total_list_of_results[[i]][[1]][[1]],
                total_list_of_results[[i]][[2]][[1]],
                total_list_of_results[[i]][[3]][[1]],
                total_list_of_results[[i]][[4]][[1]],
                total_list_of_results[[i]][[5]][[1]],
                total_list_of_results[[i]][[6]][[1]],
                total_list_of_results[[i]][[7]][[1]],
                total_list_of_results[[i]][[8]][[1]])
    new_df<-cbind.data.frame(new_line,rep(knot_values[i],length(new_line)))
    tot_vals<-cbind(tot_vals,new_line)
    tot_plot_vals<-rbind(tot_plot_vals,new_df)
    
  }
  tot_vals<-data.frame(tot_vals)
  row.names(tot_vals)<-paste(c(rep("first",4),rep("second",4)),rep(c(100,500,1000,5000),2),sep = " ")
  
  tot_plot_vals<-data.frame(tot_plot_vals)
  names(tot_vals)<-paste("Knot_number",knot_values,sep = "_")
  names(tot_plot_vals)<-c("percent","knot_number")
  tot_plot_vals$sample_size<-rep(c(100,500,1000,5000),nrow(tot_plot_vals)/4)
  tot_plot_vals$order<-rep(c(rep("first",4),rep("second",4)),nrow(tot_plot_vals)/8)
  tot_plot_vals$type<-paste(tot_plot_vals$knot_number,tot_plot_vals$order)
  
  rmse_plotto<-ggplot(data = tot_plot_vals,aes(x=sample_size,y=percent,group=type))+
    geom_line(aes(colour=knot_number,linetype=order),size=1.05)+labs(x="Sample size",y="Percent of times true value is within credible interval",title=plot_title)+
    geom_point(aes(shape=order,fill=knot_number,colour=knot_number),fill="white",size=3)+
    scale_shape_manual("Order",values = c("first"=21,"second"=23))+
    scale_linetype_manual("Order",
                          values = c("first"="solid","second"="dotdash"))
  
  
  return(list(df_rmse=tot_vals,plot_df=tot_plot_vals,plotto=rmse_plotto))
  
  
}

### 7 knotters ####
seven_knot_95_prev<-lapply(seven_knot_results,credible_interval_test_multi_knotters,
                           true_df_vector=sim_model_foi$sim_df$prev_percent)
seven_knot_95_inc<-lapply(seven_knot_results,credible_interval_test_multi_knotters,
                          true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
seven_knot_95_kappa<-lapply(seven_knot_results,credible_interval_test_multi_knotters,
                            true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")

#### 8 knotters #####
eight_knot_95_prev<-lapply(eight_knot_results,credible_interval_test_multi_knotters,
                           true_df_vector=sim_model_foi$sim_df$prev_percent)
eight_knot_95_inc<-lapply(eight_knot_results,credible_interval_test_multi_knotters,
                          true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
eight_knot_95_kappa<-lapply(eight_knot_results,credible_interval_test_multi_knotters,
                            true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")
##### 9 knotters ####
nine_knot_95_prev<-lapply(nine_knot_results,credible_interval_test_multi_knotters,
                           true_df_vector=sim_model_foi$sim_df$prev_percent)
nine_knot_95_inc<-lapply(nine_knot_results,credible_interval_test_multi_knotters,
                          true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
nine_knot_95_kappa<-lapply(nine_knot_results,credible_interval_test_multi_knotters,
                            true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")
##### 10 knotters ####
ten_knot_95_prev<-lapply(ten_knot_results,credible_interval_test_multi_knotters,
                           true_df_vector=sim_model_foi$sim_df$prev_percent)
ten_knot_95_inc<-lapply(ten_knot_results,credible_interval_test_multi_knotters,
                          true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
ten_knot_95_kappa<-lapply(ten_knot_results,credible_interval_test_multi_knotters,
                            true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")

##### 11 knotters #####
eleven_knot_95_prev<-lapply(eleven_knot_results,credible_interval_test_multi_knotters,
                         true_df_vector=sim_model_foi$sim_df$prev_percent)
eleven_knot_95_inc<-lapply(eleven_knot_results,credible_interval_test_multi_knotters,
                        true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
eleven_knot_95_kappa<-lapply(eleven_knot_results,credible_interval_test_multi_knotters,
                          true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")

##### 12 knotters #####
twelve_knot_95_prev<-lapply(twelve_knot_results,credible_interval_test_multi_knotters,
                            true_df_vector=sim_model_foi$sim_df$prev_percent)
twelve_knot_95_inc<-lapply(twelve_knot_results,credible_interval_test_multi_knotters,
                           true_df_vector=sim_model_foi$sim_df$lambda,metric="incidence")
twelve_knot_95_kappa<-lapply(twelve_knot_results,credible_interval_test_multi_knotters,
                             true_df_vector=sim_model_foi$sim_df$kappa,metric="kappa")


tot_95_prev<-list(seven_knot_95_prev,eight_knot_95_prev,nine_knot_95_prev,
                  ten_knot_95_prev,eleven_knot_95_prev,twelve_knot_95_prev)
tot_95_inc<-list(seven_knot_95_inc,eight_knot_95_inc,nine_knot_95_inc,
                 ten_knot_95_inc,eleven_knot_95_inc,twelve_knot_95_inc)
tot_95_kappa<-list(seven_knot_95_kappa,eight_knot_95_kappa,nine_knot_95_kappa,
                   ten_knot_95_kappa,eleven_knot_95_kappa,twelve_knot_95_kappa)


prev_95<-plotter_and_table_giver(tot_95_prev,plot_title = "Prevalence within 95%")
prev_95$plotto

inc_95<-plotter_and_table_giver(tot_95_inc,plot_title = "Incidence within 95%")
inc_95$plotto

kappa_95<-plotter_and_table_giver(tot_95_kappa,plot_title = "Kappa within 95%")
kappa_95$plotto

save(prev_95,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/95_conf_anlaysis/prevalence_95_multi")
save(inc_95,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/95_conf_anlaysis/incidence_95_multi")
save(kappa_95,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/multi_knot_analysis/95_conf_anlaysis/kappa_95_multi")
