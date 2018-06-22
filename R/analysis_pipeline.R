#######################################################################################################################################
## Pipeline for analysing cluster results #############################################################################################
#######################################################################################################################################

require(ggplot2)
require(reshape2)
require(gridExtra)
require(ggpubr)

#######################################################################################################################################
## Load the results and the true epidemic #############################################################################################
#######################################################################################################################################



#######################################################################################
## Now we'll get the mean value plots #################################################
#######################################################################################
mean_value_function<-function(iterations,nrow_per_iteration,data_list,metric="prevalence"){
  data_prev<-NULL
  iter_value<-iterations - 1
  nrow_value<-nrow_per_iteration
  if(metric == "prevalence"){
    data_frame <- data_list$prev
  }
  if(metric == "incidence"){
    data_frame <- data_list$incidence
  }
  if(metric == "kappa"){
    data_frame <- data_list$kappa
  }
  
  
  for(i in 1:nrow_value){
    values<-data_frame$median[0:iter_value*nrow_value+i]
    low<-data_frame$low[0:iter_value*nrow_value+i]
    high<-data_frame$high[0:iter_value*nrow_value+i]
    time_of_values<-data_frame$time[0:iter_value*nrow_value+i]
    
    row<-cbind(mean(low),mean(values),mean(high),mean(time_of_values))
    
    data_prev<-rbind.data.frame(data_prev,row)
    
  }
  
  names(data_prev)<-c("low","median","high","time")
  
  return(data_prev)
  
}

mean_value_plotter<-function(list_of_results,plot_titles,simul_datam,metric="prevalence",limits_to_y=c(0,50)){
  list_leng<-length(list_of_results)
  plot_list<-vector("list",list_leng)
  
  if(metric == "prevalence"){
    simmed_data<-simul_datam[,c(4,5)]
    
  }
  if(metric == "incidence"){
    simmed_data<-simul_datam[,c(2,5)]
    
  }
  if(metric == "kappa"){
    simmed_data <- simul_datam[,c(1,5)]
  }
  for(i in 1:length(list_of_results)){
    plotto<-ggplot(data = list_of_results[[i]])+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
      geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
      geom_line(data = simmed_data ,aes(x=time,y=simmed_data[,1]),colour="red")+
      labs(x="time",y="Metric",title=plot_titles[i]) + coord_cartesian(ylim = limits_to_y)
    plot_list[[i]]<-plotto
  }
  
  tot_plotto<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
                        plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                        ncol = 4, nrow = 2)
  
  return(list(plots = plot_list, arranged_plots = tot_plotto))
  
}

#######################################################################################
## Now we'll get the rmse through time plots ##########################################
#######################################################################################
rmse_per_time<-function(true_data,fitted_data,metric="prevalence",year_time_series=seq(1970,2020,0.36),type_of_data){
  
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
  
  
  mean_rmse_tot$sample_size<-sample_size
  mean_rmse_tot$order_type<-order_type
  names(mean_rmse_tot)<-c("RMSE","low","median","high","time","type","sample_size","order_type")
  
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot,error_df=error_tot))
  
  
  
}

plotter_function_rmse<-function(list_of_rmse_results_dfs,plot_title,colour_by_sample_size=F){
  total_data<-list_of_rmse_results_dfs[[1]][[2]]
  
  for(i in 2:length(list_of_rmse_results_dfs)){
    total_data<-rbind.data.frame(total_data,list_of_rmse_results_dfs[[i]][[2]])
  }
  
  
  names_vector<-c(list_of_rmse_results_dfs[[1]][[2]]$type[1],
                  list_of_rmse_results_dfs[[2]][[2]]$type[1],
                  list_of_rmse_results_dfs[[3]][[2]]$type[1],
                  list_of_rmse_results_dfs[[4]][[2]]$type[1])
  colour_vector<-c("dodgerblue","red","blueviolet","green2")
  names(colour_vector) <- names_vector
  mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="RMSE",title=plot_title)
  
  if(colour_by_sample_size==T){
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=sample_size,linetype=order_type),size=1.2)+
      #scale_colour_manual("Fitting method",values = colour_vector)+
      labs(x="time",y="RMSE",title=plot_title)
  }
  
  error_plot<-ggplot(data = total_data,aes(x=time,y=median,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high,fill=type),alpha=0.36)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    #scale_fill_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="error",title=plot_title)
  
  return(list(mean_plot=mean_plot,error_plot=error_plot))
  
  
}

#######################################################################################
## Now we'll get the rmse in total and during the prediction and peak periods #########
#######################################################################################
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

