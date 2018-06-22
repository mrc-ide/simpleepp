############################################################################################################################################
## Comparing the spline results for RMSE over time #########################################################################################
############################################################################################################################################
require(ggplot2)
require(reshape2)
require(ggpubr)

############################################################################################################################################
## lets load up the data ###################################################################################################################
############################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/true_epidemic",verbose = T)

## 7 knotters ####

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/"
seven_knots<-list.files(path_name)
seven_knots<-paste(path_name,seven_knots,sep = "")
lapply(seven_knots,load)
for(i in 1:length(seven_knots)){
  load(seven_knots[i],verbose = T)
}

## 8 knotters ###

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/"
eight_knots<-list.files("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/")
eight_knots<-sub("sp","C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/sp",eight_knots)
for(i in 1:length(eight_knots)){
  load(file = eight_knots[[i]],verbose = T)
}

## 9 knotters ##
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/"
nine_knots<-list.files(path_name)
nine_knots<-paste(path_name,nine_knots,sep = "")
for(i in 1:length(nine_knots)){
  load(nine_knots[i],verbose = T)
}

## 10 knotters
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/"
ten_knots<-list.files(path_name)
ten_knots<-paste(path_name,ten_knots,sep = "")
for(i in 1:length(ten_knots)){
  load(ten_knots[i],verbose = T)
}

## 11 knotters
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/"
eleven_knots<-list.files(path_name)
eleven_knots<-paste(path_name,eleven_knots,sep = "")
lapply(eleven_knots,load,verbose=T)
for(i in 1:length(eleven_knots)){
  load(eleven_knots[i],verbose = T)
}


##12 Knotters
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/"
twelve_knots<-list.files(path_name)
twelve_knots<-paste(path_name,twelve_knots,sep = "")
lapply(twelve_knots,load,verbose=T)
for(i in 1:length(twelve_knots)){
  load(twelve_knots[i],verbose = T)
}

############################################################################################################################################
## NOw we will write our function to get the rmse to the true for each time point for all the spline data ##################################
############################################################################################################################################

rmse_per_time_knotters<-function(true_data,fitted_data,metric="prevalence",year_time_series=seq(1970,2020,0.1),type_of_data){
  
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
  
  
  time_to_test<-seq((year_time_series[1]-1970)*10+1,(year_time_series[length(year_time_series)]-2020)*10+501,1)
  
  mean_rmse_tot<-NULL
  error_tot<-NULL
  
  for (i in time_to_test[1]:time_to_test[length(time_to_test)]){
    
    year<- (i-1)/10 + 1970
    
    fitted_metric_iter <- fitted_metric[fitted_metric$time == year,]
    true_metric_now <- true_metric[i]
    
    if(metric=="prevalence"){
      fitted_metric_iter<-fitted_metric_iter/100
    }
    
    
    error <- (fitted_metric_iter$median) - rep(true_metric_now,nrow(fitted_metric_iter))
    
    rmse <- sqrt(mean(error^2))
    
    lower_bound<-quantile(error,probs = c(0.025))
    
    upper_bound<-quantile(error,probs = c(0.975))
    
    median<-median(error)
    
    mean_rmse <- cbind(rmse,lower_bound,median,upper_bound, year)
    
    error_tot <- rbind(error,error_tot)
    
    mean_rmse_tot <- rbind(mean_rmse_tot,mean_rmse)
  }
  
  
  mean_overall_rmse <- mean(mean_rmse_tot[,1])
  mean_rmse_tot<-data.frame(mean_rmse_tot)
  mean_rmse_tot$type<-rep(type_of_data,nrow(mean_rmse_tot))
  order_type<-NULL
  if(grepl("first",mean_rmse_tot$type[1])== T){
    order_type<-"first"
  }
  if(grepl("second", mean_rmse_tot$type[1]) == T){
    order_type<-"second"
  }
  
  sample_size<-NULL
  if(grepl("100",mean_rmse_tot$type[1])==T){
    sample_size<-"100"
  }
  
  if(grepl("500",mean_rmse_tot$type[1])==T){
    sample_size<-"500"
  }
  if(grepl("1000",mean_rmse_tot$type[1])==T){
    sample_size<-"1000"
  }
  if(grepl("5000",mean_rmse_tot$type[1])==T){
    sample_size<-"5000"
  }
  
  knot_number<-NULL
  if(grepl("7 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"7"
  }
  if(grepl("8 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"8"
  }
  if(grepl("9 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"9"
  }
  if(grepl("10 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"10"
  }
  if(grepl("11 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"11"
  }
  if(grepl("12 knots",mean_rmse_tot$type[1])==T){
    knot_number<-"12"
  }
  
  
  mean_rmse_tot$sample_size<-sample_size
  mean_rmse_tot$order_type<-order_type
  mean_rmse_tot$knot_number<-knot_number
  names(mean_rmse_tot)<-c("RMSE","low","median","high","time","type","sample_size","order_type","knots")
  
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot,error_df=error_tot))
  
  
  
}

plotter_function_rmse<-function(list_of_rmse_results_dfs,plot_title,colour_by_sample_size=F,knot_type_plot=F){
  total_data<-list_of_rmse_results_dfs[[1]][[2]]
  
  for(i in 2:length(list_of_rmse_results_dfs)){
    total_data<-rbind.data.frame(total_data,list_of_rmse_results_dfs[[i]][[2]])
  }
  
  
  mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="RMSE",title=plot_title)
  
  if(colour_by_sample_size==T){
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=sample_size,linetype=order_type),size=1.2)+
      #scale_colour_manual("Fitting method",values = colour_vector)+
      labs(x="time",y="RMSE",title=plot_title)
  }
  
  if(knot_type_plot==T){
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=knots,linetype=sample_size),size=1.05)+
      labs(x="time",y="RMSE",title=plot_title)
  }
  error_plot<-ggplot(data = total_data,aes(x=time,y=median,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high,fill=type),alpha=0.36)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    #scale_fill_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="error",title=plot_title)
  
  return(list(mean_plot=mean_plot,error_plot=error_plot))
  
  
}

########################################################################################################################################
## Now we will apply the functions to our datasets #####################################################################################
########################################################################################################################################

## 7 knotters ####

seven_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res,
                                             type_of_data = "Spline first order 100 7 knots")
seven_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res,
                                             type_of_data = "Spline first order 500 7 knots")
seven_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res,
                                             type_of_data = "Spline first order 1000 7 knots")
seven_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res,
                                             type_of_data = "Spline first order 5000 7 knots")

seven_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res,
                                             type_of_data = "Spline second order 100 7 knots")
seven_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res,
                                             type_of_data = "Spline second order 500 7 knots")
seven_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res,
                                             type_of_data = "Spline second order 1000 7 knots")
seven_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res,
                                             type_of_data = "Spline second order 5000 7 knots")

### 8 knotters ###

eight_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_8,
                                             type_of_data = "Spline first order 100 8 knots")
eight_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_8,
                                             type_of_data = "Spline first order 500 8 knots")
eight_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_8,
                                              type_of_data = "Spline first order 1000 8 knots")
eight_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_8,
                                              type_of_data = "Spline first order 5000 8 knots")


eight_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_8,
                                             type_of_data = "Spline second order 100 8 knots")
eight_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_8,
                                             type_of_data = "Spline second order 500 8 knots")
eight_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_8,
                                              type_of_data = "Spline second order 1000 8 knots")
eight_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_8,
                                              type_of_data = "Spline second order 5000 8 knots")

### 9 knotters ####

nine_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_9,
                                             type_of_data = "Spline first order 100 9 knots")
nine_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_9,
                                             type_of_data = "Spline first order 500 9 knots")
nine_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_9,
                                              type_of_data = "Spline first order 1000 9 knots")
nine_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_9,
                                              type_of_data = "Spline first order 5000 9 knots")


nine_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_9,
                                             type_of_data = "Spline second order 100 9 knots")
nine_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_9,
                                             type_of_data = "Spline second order 500 9 knots")
nine_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_9,
                                              type_of_data = "Spline second order 1000 9 knots")
nine_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_9,
                                              type_of_data = "Spline second order 5000 9 knots")

### 10 knotters ####

ten_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_ten,
                                            type_of_data = "Spline first order 100 10 knots")
ten_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_ten,
                                            type_of_data = "Spline first order 500 10 knots")
ten_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_ten,
                                             type_of_data = "Spline first order 1000 10 knots")
ten_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_ten,
                                             type_of_data = "Spline first order 5000 10 knots")


ten_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_ten,
                                            type_of_data = "Spline second order 100 10 knots")
ten_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_ten,
                                            type_of_data = "Spline second order 500 10 knots")
ten_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_ten,
                                             type_of_data = "Spline second order 1000 10 knots")
ten_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_ten,
                                             type_of_data = "Spline second order 5000 10 knots")

### 11 knotters ###
eleven_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_11,
                                           type_of_data = "Spline first order 100 11 knots")
eleven_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_11,
                                           type_of_data = "Spline first order 500 11 knots")
eleven_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_11,
                                            type_of_data = "Spline first order 1000 11 knots")
eleven_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_11,
                                            type_of_data = "Spline first order 5000 11 knots")


eleven_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_11,
                                           type_of_data = "Spline second order 100 11 knots")
eleven_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_11,
                                           type_of_data = "Spline second order 500 11 knots")
eleven_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_11,
                                            type_of_data = "Spline second order 1000 11 knots")
eleven_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_11,
                                            type_of_data = "Spline second order 5000 11 knots")

### 12 knotters ####

twelve_knots_sp_1_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_12,
                                              type_of_data = "Spline first order 100 12 knots")
twelve_knots_sp_1_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_12,
                                              type_of_data = "Spline first order 500 12 knots")
twelve_knots_sp_1_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_12,
                                               type_of_data = "Spline first order 1000 12 knots")
twelve_knots_sp_1_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_12,
                                               type_of_data = "Spline first order 5000 12 knots")


twelve_knots_sp_2_100<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_12,
                                              type_of_data = "Spline second order 100 12 knots")
twelve_knots_sp_2_500<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_12,
                                              type_of_data = "Spline second order 500 12 knots")
twelve_knots_sp_2_1000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_12,
                                               type_of_data = "Spline second order 1000 12 knots")
twelve_knots_sp_2_5000<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_12,
                                               type_of_data = "Spline second order 5000 12 knots")
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########
### Now that we've loaded up the results and run them through our function to extract the RMSE for prevalence we can now plot the results #
##### !!!!! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ##

first_order_splines<-list(seven_knots_sp_1_100,seven_knots_sp_1_500,seven_knots_sp_1_1000,seven_knots_sp_1_5000,
                          eight_knots_sp_1_100,eight_knots_sp_1_500,eight_knots_sp_1_1000,eight_knots_sp_1_5000,
                          nine_knots_sp_1_100,nine_knots_sp_1_500,nine_knots_sp_1_1000,nine_knots_sp_1_5000,
                          ten_knots_sp_1_100,ten_knots_sp_1_500,ten_knots_sp_1_1000,ten_knots_sp_1_5000,
                          eleven_knots_sp_1_100,eleven_knots_sp_1_500,eleven_knots_sp_1_1000,eleven_knots_sp_1_5000,
                          twelve_knots_sp_1_100,twelve_knots_sp_1_500,twelve_knots_sp_1_1000,twelve_knots_sp_1_5000)
first_order_spline_plot<-plotter_function_rmse(list_of_rmse_results_dfs = first_order_splines,
                                               plot_title = "Comparison of different knot numbers",knot_type_plot = T)
first_order_spline_plot$mean_plot
first_order_spline_plot$error_plot

second_ord_splines<-list(seven_knots_sp_2_100,seven_knots_sp_2_500,seven_knots_sp_2_1000,seven_knots_sp_2_5000,
                         eight_knots_sp_2_100,eight_knots_sp_2_500,eight_knots_sp_2_1000,eight_knots_sp_2_5000,
                         nine_knots_sp_2_100,nine_knots_sp_2_500,nine_knots_sp_2_1000,nine_knots_sp_2_5000,
                         ten_knots_sp_2_100,ten_knots_sp_2_500,ten_knots_sp_2_1000,ten_knots_sp_2_5000,
                         eleven_knots_sp_2_100,eleven_knots_sp_2_500,eleven_knots_sp_2_1000,eleven_knots_sp_2_5000,
                         twelve_knots_sp_2_100,twelve_knots_sp_2_500,twelve_knots_sp_2_1000,twelve_knots_sp_2_5000)
sec_ord_spline_plot<-plotter_function_rmse(second_ord_splines,plot_title = "Second order splines comparison",
                                           knot_type_plot = T)
sec_ord_spline_plot$mean_plot + coord_cartesian(xlim = c(2015,2020),ylim = c(0,0.3))

tot_splines<-list(seven_knots_sp_1_100,seven_knots_sp_1_500,seven_knots_sp_1_1000,seven_knots_sp_1_5000,
                  eight_knots_sp_1_100,eight_knots_sp_1_500,eight_knots_sp_1_1000,eight_knots_sp_1_5000,
                  nine_knots_sp_1_100,nine_knots_sp_1_500,nine_knots_sp_1_1000,nine_knots_sp_1_5000,
                  ten_knots_sp_1_100,ten_knots_sp_1_500,ten_knots_sp_1_1000,ten_knots_sp_1_5000,
                  eleven_knots_sp_1_100,eleven_knots_sp_1_500,eleven_knots_sp_1_1000,eleven_knots_sp_1_5000,
                  twelve_knots_sp_1_100,twelve_knots_sp_1_500,twelve_knots_sp_1_1000,twelve_knots_sp_1_5000,
                  seven_knots_sp_2_100,seven_knots_sp_2_500,seven_knots_sp_2_1000,seven_knots_sp_2_5000,
                  eight_knots_sp_2_100,eight_knots_sp_2_500,eight_knots_sp_2_1000,eight_knots_sp_2_5000,
                  nine_knots_sp_2_100,nine_knots_sp_2_500,nine_knots_sp_2_1000,nine_knots_sp_2_5000,
                  ten_knots_sp_2_100,ten_knots_sp_2_500,ten_knots_sp_2_1000,ten_knots_sp_2_5000,
                  eleven_knots_sp_2_100,eleven_knots_sp_2_500,eleven_knots_sp_2_1000,eleven_knots_sp_2_5000,
                  twelve_knots_sp_2_100,twelve_knots_sp_2_500,twelve_knots_sp_2_1000,twelve_knots_sp_2_5000)
tot_splino_ploto<-plotter_function_rmse(tot_splines,plot_title = "COmparison of spline results",
                                        knot_type_plot = T)
tot_splino_ploto$mean_plot + coord_cartesian(xlim = c(2015,2020),ylim = c(0,0.3))


save(first_order_spline_plot,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/first_order_splines_compo_OBJECT") 

save(sec_ord_spline_plot,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/sec_order_splines_compo_OBJECT")

save(tot_splino_ploto,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/tot_splines_compo_OBJECT")

##### n = 100 #####

splines_n_100<-list(seven_knots_sp_1_100,seven_knots_sp_2_100,
                    eight_knots_sp_1_100,eight_knots_sp_2_100,
                    nine_knots_sp_1_100,nine_knots_sp_2_100,
                    ten_knots_sp_1_100,ten_knots_sp_2_100,
                    eleven_knots_sp_1_100,eleven_knots_sp_2_100,
                    twelve_knots_sp_1_100,twelve_knots_sp_2_100)
n_100_prev_plots<-plotter_function_rmse(splines_n_100,plot_title = "Prevalence fits at n 100",knot_type_plot = T)
a<-n_100_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.3))

##### n = 500 ######

splines_n_500<-list(seven_knots_sp_1_500,seven_knots_sp_2_500,
                    eight_knots_sp_1_500,eight_knots_sp_2_500,
                    nine_knots_sp_1_500,nine_knots_sp_2_500,
                    ten_knots_sp_1_500,ten_knots_sp_2_500,
                    eleven_knots_sp_1_500,eleven_knots_sp_2_500,
                    twelve_knots_sp_1_500,twelve_knots_sp_2_500)
n_500_prev_plots<-plotter_function_rmse(splines_n_500,plot_title = "Prevalence fits at n 500", knot_type_plot = T)
n_500_prev_plots$mean_plot
b<-n_500_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.3))

#### n = 1k #####

splines_n_1000<-list(seven_knots_sp_1_1000,seven_knots_sp_2_1000,
                     eight_knots_sp_1_1000,eight_knots_sp_2_1000,
                     nine_knots_sp_1_1000,nine_knots_sp_2_1000,
                     ten_knots_sp_1_1000,ten_knots_sp_2_1000,
                     eleven_knots_sp_1_1000,eleven_knots_sp_2_1000,
                     twelve_knots_sp_1_1000,twelve_knots_sp_2_1000)
n_1000_prev_plots<-plotter_function_rmse(splines_n_1000,plot_title = "Prevalence fits at n 1000",knot_type_plot = T)
n_1000_prev_plots$mean_plot
c<-n_1000_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.3))

##### n = 5k ######

splines_n_5000<-list(seven_knots_sp_1_5000,seven_knots_sp_2_5000,
                     eight_knots_sp_1_5000,eight_knots_sp_2_5000,
                     nine_knots_sp_1_5000,nine_knots_sp_2_5000,
                     ten_knots_sp_1_5000,ten_knots_sp_2_5000,
                     eleven_knots_sp_1_5000,eleven_knots_sp_2_5000,
                     twelve_knots_sp_1_5000,twelve_knots_sp_2_5000)
n_5000_prev_plots<-plotter_function_rmse(splines_n_5000,plot_title = "Prevalence fits at n 5000",knot_type_plot = T)
n_5000_prev_plots$mean_plot

d<-n_5000_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.3))

tot_plot_by_sample_size<-ggarrange(a,b,c,d,ncol = 1,nrow = 4,align = "hv")
tot_plot_by_sample_size

save(tot_plot_by_sample_size,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/by_sample_size_compo_plot_OBJECT")

######## 11111!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!##########
#### NOW WE WILL DO THE SAME FOR INCIDENCE, SEE HOW THAT IS AFFECTED BY THE KNOT NUMBER WE USED TO MODEL INCIDENCE WITH ####################
#########()()()()()()()()()()()()()()()()()()()()()()()()()(()!()!()!()()!(!)(!)(!)()(!)(!)(!)(!)(!)(!)(!)(!)(!)!()!()!(!)(!)(!)!()!()######

seven_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res,
                                             type_of_data = "Spline first order 100 7 knots",metric = "incidence")
seven_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res,
                                             type_of_data = "Spline first order 500 7 knots",metric = "incidence")
seven_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res,
                                              type_of_data = "Spline first order 1000 7 knots",metric = "incidence")
seven_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res,
                                              type_of_data = "Spline first order 5000 7 knots",metric = "incidence")

seven_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res,
                                             type_of_data = "Spline second order 100 7 knots",metric = "incidence")
seven_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res,
                                             type_of_data = "Spline second order 500 7 knots",metric = "incidence")
seven_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res,
                                              type_of_data = "Spline second order 1000 7 knots",metric = "incidence")
seven_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res,
                                              type_of_data = "Spline second order 5000 7 knots",metric = "incidence")

### 8 knotters ###

eight_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_8,
                                             type_of_data = "Spline first order 100 8 knots",metric = "incidence")
eight_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_8,
                                             type_of_data = "Spline first order 500 8 knots",metric = "incidence")
eight_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_8,
                                              type_of_data = "Spline first order 1000 8 knots",metric = "incidence")
eight_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_8,
                                              type_of_data = "Spline first order 5000 8 knots",metric = "incidence")


eight_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_8,
                                             type_of_data = "Spline second order 100 8 knots",metric = "incidence")
eight_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_8,
                                             type_of_data = "Spline second order 500 8 knots",metric = "incidence")
eight_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_8,
                                              type_of_data = "Spline second order 1000 8 knots",metric = "incidence")
eight_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_8,
                                              type_of_data = "Spline second order 5000 8 knots",metric = "incidence")

### 9 knotters ####

nine_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_9,
                                            type_of_data = "Spline first order 100 9 knots",metric = "incidence")
nine_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_9,
                                            type_of_data = "Spline first order 500 9 knots",metric = "incidence")
nine_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_9,
                                             type_of_data = "Spline first order 1000 9 knots",metric = "incidence")
nine_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_9,
                                             type_of_data = "Spline first order 5000 9 knots",metric = "incidence")


nine_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_9,
                                            type_of_data = "Spline second order 100 9 knots",metric = "incidence")
nine_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_9,
                                            type_of_data = "Spline second order 500 9 knots",metric = "incidence")
nine_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_9,
                                             type_of_data = "Spline second order 1000 9 knots",metric = "incidence")
nine_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_9,
                                             type_of_data = "Spline second order 5000 9 knots",metric = "incidence")

### 10 knotters ####

ten_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_ten,
                                           type_of_data = "Spline first order 100 10 knots",metric = "incidence")
ten_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_ten,
                                           type_of_data = "Spline first order 500 10 knots",metric = "incidence")
ten_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_ten,
                                            type_of_data = "Spline first order 1000 10 knots",metric = "incidence")
ten_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_ten,
                                            type_of_data = "Spline first order 5000 10 knots",metric = "incidence")


ten_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_ten,
                                           type_of_data = "Spline second order 100 10 knots",metric = "incidence")
ten_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_ten,
                                           type_of_data = "Spline second order 500 10 knots",metric = "incidence")
ten_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_ten,
                                            type_of_data = "Spline second order 1000 10 knots",metric = "incidence")
ten_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_ten,
                                            type_of_data = "Spline second order 5000 10 knots",metric = "incidence")

### 11 knotters ###
eleven_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_11,
                                              type_of_data = "Spline first order 100 11 knots",metric = "incidence")
eleven_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_11,
                                              type_of_data = "Spline first order 500 11 knots",metric = "incidence")
eleven_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_11,
                                               type_of_data = "Spline first order 1000 11 knots",metric = "incidence")
eleven_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_11,
                                               type_of_data = "Spline first order 5000 11 knots",metric = "incidence")


eleven_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_11,
                                              type_of_data = "Spline second order 100 11 knots",metric = "incidence")
eleven_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_11,
                                              type_of_data = "Spline second order 500 11 knots",metric = "incidence")
eleven_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_11,
                                               type_of_data = "Spline second order 1000 11 knots",metric = "incidence")
eleven_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_11,
                                               type_of_data = "Spline second order 5000 11 knots",metric = "incidence")

### 12 knotters ####

twelve_knots_sp_1_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_12,
                                              type_of_data = "Spline first order 100 12 knots",metric = "incidence")
twelve_knots_sp_1_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_12,
                                              type_of_data = "Spline first order 500 12 knots",metric = "incidence")
twelve_knots_sp_1_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_12,
                                               type_of_data = "Spline first order 1000 12 knots",metric = "incidence")
twelve_knots_sp_1_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_12,
                                               type_of_data = "Spline first order 5000 12 knots",metric = "incidence")


twelve_knots_sp_2_100_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_12,
                                              type_of_data = "Spline second order 100 12 knots",metric = "incidence")
twelve_knots_sp_2_500_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_12,
                                              type_of_data = "Spline second order 500 12 knots",metric = "incidence")
twelve_knots_sp_2_1000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_12,
                                               type_of_data = "Spline second order 1000 12 knots",metric = "incidence")
twelve_knots_sp_2_5000_inc<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_12,
                                               type_of_data = "Spline second order 5000 12 knots",metric = "incidence")

#### tot inc plots ######

first_ord_vals_inc<-list(seven_knots_sp_1_100_inc,seven_knots_sp_1_500_inc,seven_knots_sp_1_1000_inc,seven_knots_sp_1_5000_inc,
                         eight_knots_sp_1_100_inc,eight_knots_sp_1_500_inc,eight_knots_sp_1_1000_inc,eight_knots_sp_1_5000_inc,
                         nine_knots_sp_1_100_inc,nine_knots_sp_1_500_inc,nine_knots_sp_1_1000_inc,nine_knots_sp_1_5000_inc,
                         ten_knots_sp_1_100_inc,ten_knots_sp_1_500_inc,ten_knots_sp_1_1000_inc,ten_knots_sp_1_5000_inc,
                         eleven_knots_sp_1_100_inc,eleven_knots_sp_1_500_inc,eleven_knots_sp_1_1000_inc,eleven_knots_sp_1_5000_inc,
                         twelve_knots_sp_1_100_inc,twelve_knots_sp_1_500_inc,twelve_knots_sp_1_1000_inc,twelve_knots_sp_1_5000_inc)
first_order_inc_plots<-plotter_function_rmse(first_ord_vals_inc,plot_title = "Kappa comparison first order",knot_type_plot = T)
first_order_inc_plots$mean_plot 
a<-first_order_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.1))

second_ord_vals_inc<-list(seven_knots_sp_2_100_inc,seven_knots_sp_2_500_inc,seven_knots_sp_2_1000_inc,seven_knots_sp_2_5000_inc,
                          eight_knots_sp_2_100_inc,eight_knots_sp_2_500_inc,eight_knots_sp_2_1000_inc,eight_knots_sp_2_5000_inc,
                          nine_knots_sp_2_100_inc,nine_knots_sp_2_500_inc,nine_knots_sp_2_1000_inc,nine_knots_sp_2_5000_inc,
                          ten_knots_sp_2_100_inc,ten_knots_sp_2_500_inc,ten_knots_sp_2_1000_inc,ten_knots_sp_2_5000_inc,
                          eleven_knots_sp_2_100_inc,eleven_knots_sp_2_500_inc,eleven_knots_sp_2_1000_inc,eleven_knots_sp_2_5000_inc,
                          twelve_knots_sp_2_100_inc,twelve_knots_sp_2_500_inc,twelve_knots_sp_2_1000_inc,twelve_knots_sp_2_5000_inc)
second_ord_inc_plots<-plotter_function_rmse(second_ord_vals_inc,plot_title = "Inicidence comparison second order",knot_type_plot = T)
second_ord_inc_plots$mean_plot
b<-second_ord_inc_plots$mean_plot + coord_cartesian(ylim=c(0,0.1))

first_and_sec_inc<-ggarrange(a,b,ncol = 2,nrow = 1,align = "hv")

save(first_and_sec_inc,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/INCIDNECE,first_and_sec_orde_compo_OBJECT")

### 100 n ######

splines_n_100_inc<-list(seven_knots_sp_1_100_inc,seven_knots_sp_2_100_inc,
                    eight_knots_sp_1_100_inc,eight_knots_sp_2_100_inc,
                    nine_knots_sp_1_100_inc,nine_knots_sp_2_100_inc,
                    ten_knots_sp_1_100_inc,ten_knots_sp_2_100_inc,
                    eleven_knots_sp_1_100_inc,eleven_knots_sp_2_100_inc,
                    twelve_knots_sp_1_100_inc,twelve_knots_sp_2_100_inc)
n_100_inc_plots<-plotter_function_rmse(splines_n_100_inc,plot_title = "Incidence fits at n 100",knot_type_plot = T)
a<-n_100_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.1))

##### n = 500 ######

splines_n_500_inc<-list(seven_knots_sp_1_500_inc,seven_knots_sp_2_500_inc,
                    eight_knots_sp_1_500_inc,eight_knots_sp_2_500_inc,
                    nine_knots_sp_1_500_inc,nine_knots_sp_2_500_inc,
                    ten_knots_sp_1_500_inc,ten_knots_sp_2_500_inc,
                    eleven_knots_sp_1_500_inc,eleven_knots_sp_2_500_inc,
                    twelve_knots_sp_1_500_inc,twelve_knots_sp_2_500_inc)
n_500_inc_plots<-plotter_function_rmse(splines_n_500_inc,plot_title = "Incidence fits at n 500", knot_type_plot = T)
n_500_inc_plots$mean_plot
b<-n_500_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.1))

#### n = 1k #####

splines_n_1000_inc<-list(seven_knots_sp_1_1000_inc,seven_knots_sp_2_1000_inc,
                     eight_knots_sp_1_1000_inc,eight_knots_sp_2_1000_inc,
                     nine_knots_sp_1_1000_inc,nine_knots_sp_2_1000_inc,
                     ten_knots_sp_1_1000_inc,ten_knots_sp_2_1000_inc,
                     eleven_knots_sp_1_1000_inc,eleven_knots_sp_2_1000_inc,
                     twelve_knots_sp_1_1000_inc,twelve_knots_sp_2_1000_inc)
n_1000_inc_plots<-plotter_function_rmse(splines_n_1000_inc,plot_title = "Incidence fits at n 1000",knot_type_plot = T)
n_1000_inc_plots$mean_plot
c<-n_1000_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.1))

##### n = 5k ######

splines_n_5000_inc<-list(seven_knots_sp_1_5000_inc,seven_knots_sp_2_5000_inc,
                     eight_knots_sp_1_5000_inc,eight_knots_sp_2_5000_inc,
                     nine_knots_sp_1_5000_inc,nine_knots_sp_2_5000_inc,
                     ten_knots_sp_1_5000_inc,ten_knots_sp_2_5000_inc,
                     eleven_knots_sp_1_5000_inc,eleven_knots_sp_2_5000_inc,
                     twelve_knots_sp_1_5000_inc,twelve_knots_sp_2_5000_inc)
n_5000_inc_plots<-plotter_function_rmse(splines_n_5000_inc,plot_title = "Incidence fits at n 5000",knot_type_plot = T)
n_5000_inc_plots$mean_plot

d<-n_5000_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.1))

tot_plot_by_sample_size<-ggarrange(a,b,c,d,ncol = 1,nrow = 4,align = "hv")
tot_plot_by_sample_size

save(tot_plot_by_sample_size,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/INCIDENCE_by_sample_size_compo_OBJECT")

############################################################################################################################################
## Now we will plot out the results for our kappa parameter ################################################################################
############################################################################################################################################

seven_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res,
                                                 type_of_data = "Spline first order 100 7 knots",metric = "kappa")
seven_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res,
                                                 type_of_data = "Spline first order 500 7 knots",metric = "kappa")
seven_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res,
                                                  type_of_data = "Spline first order 1000 7 knots",metric = "kappa")
seven_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res,
                                                  type_of_data = "Spline first order 5000 7 knots",metric = "kappa")

seven_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res,
                                                 type_of_data = "Spline second order 100 7 knots",metric = "kappa")
seven_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res,
                                                 type_of_data = "Spline second order 500 7 knots",metric = "kappa")
seven_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res,
                                                  type_of_data = "Spline second order 1000 7 knots",metric = "kappa")
seven_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res,
                                                  type_of_data = "Spline second order 5000 7 knots",metric = "kappa")

### 8 knotters ###

eight_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_8,
                                                 type_of_data = "Spline first order 100 8 knots",metric = "kappa")
eight_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_8,
                                                 type_of_data = "Spline first order 500 8 knots",metric = "kappa")
eight_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_8,
                                                  type_of_data = "Spline first order 1000 8 knots",metric = "kappa")
eight_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_8,
                                                  type_of_data = "Spline first order 5000 8 knots",metric = "kappa")


eight_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_8,
                                                 type_of_data = "Spline second order 100 8 knots",metric = "kappa")
eight_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_8,
                                                 type_of_data = "Spline second order 500 8 knots",metric = "kappa")
eight_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_8,
                                                  type_of_data = "Spline second order 1000 8 knots",metric = "kappa")
eight_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_8,
                                                  type_of_data = "Spline second order 5000 8 knots",metric = "kappa")

### 9 knotters ####

nine_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_9,
                                                type_of_data = "Spline first order 100 9 knots",metric = "kappa")
nine_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_9,
                                                type_of_data = "Spline first order 500 9 knots",metric = "kappa")
nine_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_9,
                                                 type_of_data = "Spline first order 1000 9 knots",metric = "kappa")
nine_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_9,
                                                 type_of_data = "Spline first order 5000 9 knots",metric = "kappa")


nine_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_9,
                                                type_of_data = "Spline second order 100 9 knots",metric = "kappa")
nine_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_9,
                                                type_of_data = "Spline second order 500 9 knots",metric = "kappa")
nine_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_9,
                                                 type_of_data = "Spline second order 1000 9 knots",metric = "kappa")
nine_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_9,
                                                 type_of_data = "Spline second order 5000 9 knots",metric = "kappa")

### 10 knotters ####

ten_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_ten,
                                               type_of_data = "Spline first order 100 10 knots",metric = "kappa")
ten_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_ten,
                                               type_of_data = "Spline first order 500 10 knots",metric = "kappa")
ten_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_ten,
                                                type_of_data = "Spline first order 1000 10 knots",metric = "kappa")
ten_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_ten,
                                                type_of_data = "Spline first order 5000 10 knots",metric = "kappa")


ten_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_ten,
                                               type_of_data = "Spline second order 100 10 knots",metric = "kappa")
ten_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_ten,
                                               type_of_data = "Spline second order 500 10 knots",metric = "kappa")
ten_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_ten,
                                                type_of_data = "Spline second order 1000 10 knots",metric = "kappa")
ten_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_ten,
                                                type_of_data = "Spline second order 5000 10 knots",metric = "kappa")

### 11 knotters ###
eleven_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_11,
                                                  type_of_data = "Spline first order 100 11 knots",metric = "kappa")
eleven_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_11,
                                                  type_of_data = "Spline first order 500 11 knots",metric = "kappa")
eleven_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_11,
                                                   type_of_data = "Spline first order 1000 11 knots",metric = "kappa")
eleven_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_11,
                                                   type_of_data = "Spline first order 5000 11 knots",metric = "kappa")


eleven_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_11,
                                                  type_of_data = "Spline second order 100 11 knots",metric = "kappa")
eleven_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_11,
                                                  type_of_data = "Spline second order 500 11 knots",metric = "kappa")
eleven_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_11,
                                                   type_of_data = "Spline second order 1000 11 knots",metric = "kappa")
eleven_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_11,
                                                   type_of_data = "Spline second order 5000 11 knots",metric = "kappa")

### 12 knotters ####

twelve_knots_sp_1_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_100_foi_res_12,
                                                  type_of_data = "Spline first order 100 12 knots",metric = "kappa")
twelve_knots_sp_1_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_500_foi_res_12,
                                                  type_of_data = "Spline first order 500 12 knots",metric = "kappa")
twelve_knots_sp_1_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_1k_foi_res_12,
                                                   type_of_data = "Spline first order 1000 12 knots",metric = "kappa")
twelve_knots_sp_1_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_f_5k_foi_res_12,
                                                   type_of_data = "Spline first order 5000 12 knots",metric = "kappa")


twelve_knots_sp_2_100_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_100_foi_res_12,
                                                  type_of_data = "Spline second order 100 12 knots",metric = "kappa")
twelve_knots_sp_2_500_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_500_foi_res_12,
                                                  type_of_data = "Spline second order 500 12 knots",metric = "kappa")
twelve_knots_sp_2_1000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_1k_foi_res_12,
                                                   type_of_data = "Spline second order 1000 12 knots",metric = "kappa")
twelve_knots_sp_2_5000_kappa<-rmse_per_time_knotters(sim_model_foi$sim_df,fitted_data = spline_s_5k_foi_res_12,
                                                   type_of_data = "Spline second order 5000 12 knots",metric = "kappa")

#### tot inc plots ######

first_ord_vals_kappa<-list(seven_knots_sp_1_100_kappa,seven_knots_sp_1_500_kappa,seven_knots_sp_1_1000_kappa,seven_knots_sp_1_5000_kappa,
                         eight_knots_sp_1_100_kappa,eight_knots_sp_1_500_kappa,eight_knots_sp_1_1000_kappa,eight_knots_sp_1_5000_kappa,
                         nine_knots_sp_1_100_kappa,nine_knots_sp_1_500_kappa,nine_knots_sp_1_1000_kappa,nine_knots_sp_1_5000_kappa,
                         ten_knots_sp_1_100_kappa,ten_knots_sp_1_500_kappa,ten_knots_sp_1_1000_kappa,ten_knots_sp_1_5000_kappa,
                         eleven_knots_sp_1_100_kappa,eleven_knots_sp_1_500_kappa,eleven_knots_sp_1_1000_kappa,eleven_knots_sp_1_5000_kappa,
                         twelve_knots_sp_1_100_kappa,twelve_knots_sp_1_500_kappa,twelve_knots_sp_1_1000_kappa,twelve_knots_sp_1_5000_kappa)
first_order_kappa_plots<-plotter_function_rmse(first_ord_vals_kappa,plot_title = "Kappa comparison first order",knot_type_plot = T)
first_order_kappa_plots$mean_plot 
a<-first_order_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.4))

second_ord_vals_kappa<-list(seven_knots_sp_2_100_kappa,seven_knots_sp_2_500_kappa,seven_knots_sp_2_1000_kappa,seven_knots_sp_2_5000_kappa,
                          eight_knots_sp_2_100_kappa,eight_knots_sp_2_500_kappa,eight_knots_sp_2_1000_kappa,eight_knots_sp_2_5000_kappa,
                          nine_knots_sp_2_100_kappa,nine_knots_sp_2_500_kappa,nine_knots_sp_2_1000_kappa,nine_knots_sp_2_5000_kappa,
                          ten_knots_sp_2_100_kappa,ten_knots_sp_2_500_kappa,ten_knots_sp_2_1000_kappa,ten_knots_sp_2_5000_kappa,
                          eleven_knots_sp_2_100_kappa,eleven_knots_sp_2_500_kappa,eleven_knots_sp_2_1000_kappa,eleven_knots_sp_2_5000_kappa,
                          twelve_knots_sp_2_100_kappa,twelve_knots_sp_2_500_kappa,twelve_knots_sp_2_1000_kappa,twelve_knots_sp_2_5000_kappa)
second_ord_kappa_plots<-plotter_function_rmse(second_ord_vals_kappa,plot_title = "Inicidence comparison second order",knot_type_plot = T)
second_ord_kappa_plots$mean_plot
b<-second_ord_kappa_plots$mean_plot + coord_cartesian(ylim=c(0,0.4))

first_and_sec_kappa<-ggarrange(a,b,ncol = 2,nrow = 1,align = "hv")

save(first_and_sec_kappa,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/KAPPA_first_and_sec_orde_compo_OBJECT")

### 100 n ######

splines_n_100_kappa<-list(seven_knots_sp_1_100_kappa,seven_knots_sp_2_100_kappa,
                        eight_knots_sp_1_100_kappa,eight_knots_sp_2_100_kappa,
                        nine_knots_sp_1_100_kappa,nine_knots_sp_2_100_kappa,
                        ten_knots_sp_1_100_kappa,ten_knots_sp_2_100_kappa,
                        eleven_knots_sp_1_100_kappa,eleven_knots_sp_2_100_kappa,
                        twelve_knots_sp_1_100_kappa,twelve_knots_sp_2_100_kappa)
n_100_kappa_plots<-plotter_function_rmse(splines_n_100_kappa,plot_title = "Kappa fits at n 100",knot_type_plot = T)
a<-n_100_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.4))

##### n = 500 ######

splines_n_500_kappa<-list(seven_knots_sp_1_500_kappa,seven_knots_sp_2_500_kappa,
                        eight_knots_sp_1_500_kappa,eight_knots_sp_2_500_kappa,
                        nine_knots_sp_1_500_kappa,nine_knots_sp_2_500_kappa,
                        ten_knots_sp_1_500_kappa,ten_knots_sp_2_500_kappa,
                        eleven_knots_sp_1_500_kappa,eleven_knots_sp_2_500_kappa,
                        twelve_knots_sp_1_500_kappa,twelve_knots_sp_2_500)
n_500_kappa_plots<-plotter_function_rmse(splines_n_500_kappa,plot_title = "Kappa fits at n 500", knot_type_plot = T)
n_500_kappa_plots$mean_plot
b<-n_500_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.4))

#### n = 1k #####

splines_n_1000_kappa<-list(seven_knots_sp_1_1000_kappa,seven_knots_sp_2_1000_kappa,
                         eight_knots_sp_1_1000_kappa,eight_knots_sp_2_1000_kappa,
                         nine_knots_sp_1_1000_kappa,nine_knots_sp_2_1000_kappa,
                         ten_knots_sp_1_1000_kappa,ten_knots_sp_2_1000_kappa,
                         eleven_knots_sp_1_1000_kappa,eleven_knots_sp_2_1000_kappa,
                         twelve_knots_sp_1_1000_kappa,twelve_knots_sp_2_1000)
n_1000_kappa_plots<-plotter_function_rmse(splines_n_1000_kappa,plot_title = "Kappa fits at n 1000",knot_type_plot = T)
n_1000_kappa_plots$mean_plot
c<-n_1000_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.4))

##### n = 5k ######

splines_n_5000_kappa<-list(seven_knots_sp_1_5000_kappa,seven_knots_sp_2_5000_kappa,
                         eight_knots_sp_1_5000_kappa,eight_knots_sp_2_5000_kappa,
                         nine_knots_sp_1_5000_kappa,nine_knots_sp_2_5000_kappa,
                         ten_knots_sp_1_5000_kappa,ten_knots_sp_2_5000_kappa,
                         eleven_knots_sp_1_5000_kappa,eleven_knots_sp_2_5000_kappa,
                         twelve_knots_sp_1_5000_kappa,twelve_knots_sp_2_5000)
n_5000_kappa_plots<-plotter_function_rmse(splines_n_5000_kappa,plot_title = "Kappa fits at n 5000",knot_type_plot = T)
n_5000_kappa_plots$mean_plot

d<-n_5000_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.4))

tot_plot_by_sample_size<-ggarrange(a,b,c,d,ncol = 1,nrow = 4,align = "hv")
tot_plot_by_sample_size

save(tot_plot_by_sample_size,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/knot_value_splines/KAPPA_by_sample_size_compo_OBJECT")



