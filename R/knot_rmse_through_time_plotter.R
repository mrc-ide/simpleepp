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

##12 Knotters
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/"
twelve_knots<-list.files(path_name)
twelve_knots<-paste(path_name,twelve_knots,sep = "")
lapply(twelve_knots,load,verbose=T)

############################################################################################################################################
## NOw we will write our function to get the rmse to the true for each time point for all the spline data ##################################
############################################################################################################################################

rmse_per_time_knotters<-function(true_data,fitted_data,metric="prevalence",year_time_series=seq(1970,2020,0.36),type_of_data){
  
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
  if(grepl("8 knots"),mean_rmse_tot$type[1]==T){
    knot_number<-"8"
  }
  if(grepl("9 knots"),mean_rmse_tot$type[1]==T){
    knot_number<-"9"
  }
  if(grepl("10 knots"),mean_rmse_tot$type[1]==T){
    knot_number<-"10"
  }
  if(grepl("11 knots"),mean_rmse_tot$type[1]==T){
    knot_number<-"11"
  }
  if(grepl("12 knots"),mean_rmse_tot$type[1]==T){
    knot_number<-"12"
  }
  
  
  mean_rmse_tot$sample_size<-sample_size
  mean_rmse_tot$order_type<-order_type
  mean_rmse_tot$knot_number<-knot_number
  names(mean_rmse_tot)<-c("RMSE","low","median","high","time","type","sample_size","order_type","knots")
  
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot,error_df=error_tot))
  
  
  
}

plotter_function_rmse<-function(list_of_rmse_results_dfs,plot_title,colour_by_sample_size=F){
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
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=knots,linetype=sample_size,size=order_type))+
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



