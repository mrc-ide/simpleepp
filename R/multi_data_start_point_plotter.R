###############################################################################################################################################
## Multi-year analysis of the different fitting methods #######################################################################################
###############################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)
require(gridExtra)

###############################################################################################################################################
## lets load up all the different results from our dataset ####################################################################################
###############################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/true_epidemic",verbose = T)

### 1984 data ####

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/"
year_1984<-list.files(path_name,full.names = T)
for(i in 1:length(year_1984)){
  load(year_1984[i],verbose = T)
}

### 1990 data ####

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1990_runs/results/"
year_1990<-list.files(path_name,full.names = T)
for(i in 1:length(year_1990)){
  load(year_1990[i],verbose = T)
}

### 1995 data ####

path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1995_runs/results/"
year_1995<-list.files(path_name,full.names = T)
for(i in 1:length(year_1995)){
  load(year_1995[i],verbose = T)
}

### 2000 data ####
path_name<-"C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_2000_runs/results/"
year_2000<-list.files(path_name,full.names = T)
for(i in 1:length(year_2000)){
  load(year_2000[i],verbose = T)
}

#############################################################################################################################################
## Lets get our plots mean values ready #####################################################################################################
#############################################################################################################################################

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

#### splines prev #####
splines_1984 <- list(sp_100_1_84_res,sp_500_1_84_res,sp_1000_1_84_res,sp_5000_1_84_res,
                     sp_100_2_84_res,sp_500_2_84_res,sp_1000_2_84_res,sp_5000_2_84_res)
splines_84_mean_prev<-lapply(splines_1984,mean_value_function,iterations=100,nrow_per_iteration=501)
plot_titulos<-c("Splines 1st n 100 Prevalence","Splines 1st n 500 Prevalence",
                "Splines 1st n 1000 Prevalence","Splines 1st n 5000 Prevalence",
                "Splines 2nd n 100 Prevalence","Splines 2nd n 500 Prevalence",
                "Splines 2nd n 500 Prevalence","Splines 2nd n 5000 Prevalence")

splines_prev_84_plots<-mean_value_plotter(splines_84_mean_prev,plot_titles = plot_titulos,
                                          simul_datam = sim_model_output$sim_df,limits_to_y = c(0,80))

#### RW prev ########

rw_1984 <- list(rw_100_1_84_res,rw_500_1_84_res,rw_1000_1_84_res,rw_5000_1_84_res,
                rw_100_2_84_res,rw_500_2_84_res,rw_1000_2_84_res,rw_5000_2_84_res)

rw_84_mean_prev<-lapply(rw_1984, mean_value_function,iterations=100,nrow_per_iteration=502)
plot_titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
                "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
                "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
                "RW 2nd n 500 Prevalence","RW 2nd n 5000 Prevalence")
rw_prev_84_plots<-mean_value_plotter(rw_84_mean_prev,plot_titles = plot_titulos,
                                     simul_datam = sim_model_output$sim_df, limits_to_y = c(0,80))

tot_84_plots<-grid.arrange(splines_prev_84_plots$arranged_plots,rw_prev_84_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1984")


####### 1990 #########
splines_1990 <- list(sp_100_1_90_res,sp_500_1_90_res,sp_1000_1_90_res,sp_5000_1_90_res,
                     sp_100_2_90_res,sp_500_2_90_res,sp_1000_2_90_res,sp_5000_2_90_res)
splines_90_mean_prev<-lapply(splines_1990,mean_value_function,iterations=100,nrow_per_iteration=501)
plot_titulos<-c("Splines 1st n 100 Prevalence","Splines 1st n 500 Prevalence",
                "Splines 1st n 1000 Prevalence","Splines 1st n 5000 Prevalence",
                "Splines 2nd n 100 Prevalence","Splines 2nd n 500 Prevalence",
                "Splines 2nd n 500 Prevalence","Splines 2nd n 5000 Prevalence")

splines_prev_90_plots<-mean_value_plotter(splines_90_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,80))

#### RW prev ########

rw_1990 <- list(rw_100_1_90_res,rw_500_1_90_res,rw_1000_1_90_res,rw_5000_1_90_res,
                rw_100_2_90_res,rw_500_2_90_res,rw_1000_2_90_res,rw_5000_2_90_res)

rw_90_mean_prev<-lapply(rw_1990, mean_value_function,iterations=100,nrow_per_iteration=502)
plot_titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
                "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
                "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
                "RW 2nd n 500 Prevalence","RW 2nd n 5000 Prevalence")
rw_prev_90_plots<-mean_value_plotter(rw_90_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,80))

tot_90_plots<-grid.arrange(splines_prev_90_plots$arranged_plots,rw_prev_90_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1990")

##### 1995, 2005 seen it with my eyes, real dope lives #####
splines_1995 <- list(sp_100_1_95_res,sp_500_1_95_res,sp_1000_1_95_res,sp_5000_1_95_res,
                     sp_100_2_95_res,sp_500_2_95_res,sp_1000_2_95_res,sp_5000_2_95_res)
splines_95_mean_prev<-lapply(splines_1995,mean_value_function,iterations=100,nrow_per_iteration=501)
plot_titulos<-c("Splines 1st n 100 Prevalence","Splines 1st n 500 Prevalence",
                "Splines 1st n 1000 Prevalence","Splines 1st n 5000 Prevalence",
                "Splines 2nd n 100 Prevalence","Splines 2nd n 500 Prevalence",
                "Splines 2nd n 500 Prevalence","Splines 2nd n 5000 Prevalence")

splines_prev_95_plots<-mean_value_plotter(splines_95_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,80))

#### RW prev ########

rw_1995 <- list(rw_100_1_95_res,rw_500_1_95_res,rw_1000_1_95_res,rw_5000_1_95_res,
                rw_100_2_95_res,rw_500_2_95_res,rw_1000_2_95_res,rw_5000_2_95_res)

rw_95_mean_prev<-lapply(rw_1995, mean_value_function,iterations=100,nrow_per_iteration=502)
plot_titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
                "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
                "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
                "RW 2nd n 500 Prevalence","RW 2nd n 5000 Prevalence")
rw_prev_95_plots<-mean_value_plotter(rw_95_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,80))

tot_95_plots<-grid.arrange(splines_prev_95_plots$arranged_plots,rw_prev_95_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1995")


##### 2000 y2k ####
splines_2000 <- list(sp_100_1_00_res,sp_500_1_00_res,sp_1000_1_00_res,sp_5000_1_00_res,
                     sp_100_2_00_res,sp_500_2_00_res,sp_1000_2_00_res,sp_5000_2_00_res)
splines_00_mean_prev<-lapply(splines_2000,mean_value_function,iterations=100,nrow_per_iteration=501)
plot_titulos<-c("Splines 1st n 100 Prevalence","Splines 1st n 500 Prevalence",
                "Splines 1st n 1000 Prevalence","Splines 1st n 5000 Prevalence",
                "Splines 2nd n 100 Prevalence","Splines 2nd n 500 Prevalence",
                "Splines 2nd n 500 Prevalence","Splines 2nd n 5000 Prevalence")

splines_prev_00_plots<-mean_value_plotter(splines_00_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,80))

#### RW prev ########

rw_2000 <- list(rw_100_1_00_res,rw_500_1_00_res,rw_1000_1_00_res,rw_5000_1_00_res,
                rw_100_2_00_res,rw_500_2_00_res,rw_1000_2_00_res,rw_5000_2_00_res)

rw_00_mean_prev<-lapply(rw_2000, mean_value_function,iterations=100,nrow_per_iteration=502)
plot_titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
                "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
                "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
                "RW 2nd n 500 Prevalence","RW 2nd n 5000 Prevalence")
rw_prev_00_plots<-mean_value_plotter(rw_00_mean_prev,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,80))

tot_00_plots<-grid.arrange(splines_prev_00_plots$arranged_plots,rw_prev_00_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 2000")
###########################################################################################################################################
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!******************************!!!!!!!!!!!!!!!!!!!!!!!!***************!!!!!!!###########
#### Now we will do this all over again for incidence #####################################################################################
#### ()(")*&$^£()("*&££&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#######
###########################################################################################################################################
splines_1984_inc <- list(sp_100_1_84_res,sp_500_1_84_res,sp_1000_1_84_res,sp_5000_1_84_res,
                     sp_100_2_84_res,sp_500_2_84_res,sp_1000_2_84_res,sp_5000_2_84_res)
splines_84_mean_inc<-lapply(splines_1984,mean_value_function,iterations=100,nrow_per_iteration=501,metric="incidence")
plot_titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
                "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
                "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
                "Splines 2nd n 500 Incidence","Splines 2nd n 5000 Incidence")

splines_inc_84_plots<-mean_value_plotter(splines_84_mean_inc,plot_titles = plot_titulos,
                                          simul_datam = sim_model_output$sim_df,limits_to_y = c(0,0.8),metric = "incidence")

#### RW inc ########

rw_1984_inc <- list(rw_100_1_84_res,rw_500_1_84_res,rw_1000_1_84_res,rw_5000_1_84_res,
                rw_100_2_84_res,rw_500_2_84_res,rw_1000_2_84_res,rw_5000_2_84_res)

rw_84_mean_inc<-lapply(rw_1984, mean_value_function,iterations=100,nrow_per_iteration=502,metric="incidence")
plot_titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
                "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
                "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
                "RW 2nd n 500 Incidence","RW 2nd n 5000 Incidence")
rw_inc_84_plots<-mean_value_plotter(rw_84_mean_inc,plot_titles = plot_titulos,
                                     simul_datam = sim_model_output$sim_df, limits_to_y = c(0,0.8),metric="incidence")

tot_84_plots<-grid.arrange(splines_inc_84_plots$arranged_plots,rw_inc_84_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1984")


####### 1990 #########
splines_1990_inc <- list(sp_100_1_90_res,sp_500_1_90_res,sp_1000_1_90_res,sp_5000_1_90_res,
                     sp_100_2_90_res,sp_500_2_90_res,sp_1000_2_90_res,sp_5000_2_90_res)
splines_90_mean_inc<-lapply(splines_1990,mean_value_function,iterations=100,nrow_per_iteration=501,metric="incidence")
plot_titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
                "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
                "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
                "Splines 2nd n 500 Incidence","Splines 2nd n 5000 Incidence")

splines_inc_90_plots<-mean_value_plotter(splines_90_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,0.8),metric="incidence")

#### RW inc ########

rw_1990_inc <- list(rw_100_1_90_res,rw_500_1_90_res,rw_1000_1_90_res,rw_5000_1_90_res,
                rw_100_2_90_res,rw_500_2_90_res,rw_1000_2_90_res,rw_5000_2_90_res)

rw_90_mean_inc<-lapply(rw_1990, mean_value_function,iterations=100,nrow_per_iteration=502,metric="incidence")
plot_titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
                "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
                "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
                "RW 2nd n 500 Incidence","RW 2nd n 5000 Incidence")
rw_inc_90_plots<-mean_value_plotter(rw_90_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,0.8),metric="incidence")

tot_90_plots<-grid.arrange(splines_inc_90_plots$arranged_plots,rw_inc_90_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1990")

##### 1995, 2005 seen it with my eyes, real dope lives #####
splines_1995_inc <- list(sp_100_1_95_res,sp_500_1_95_res,sp_1000_1_95_res,sp_5000_1_95_res,
                     sp_100_2_95_res,sp_500_2_95_res,sp_1000_2_95_res,sp_5000_2_95_res)
splines_95_mean_inc<-lapply(splines_1995,mean_value_function,iterations=100,nrow_per_iteration=501,metric="incidence")
plot_titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
                "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
                "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
                "Splines 2nd n 500 Incidence","Splines 2nd n 5000 Incidence")

splines_inc_95_plots<-mean_value_plotter(splines_95_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,0.8),metric="incidence")

#### RW inc ########

rw_1995_inc <- list(rw_100_1_95_res,rw_500_1_95_res,rw_1000_1_95_res,rw_5000_1_95_res,
                rw_100_2_95_res,rw_500_2_95_res,rw_1000_2_95_res,rw_5000_2_95_res)

rw_95_mean_inc<-lapply(rw_1995, mean_value_function,iterations=100,nrow_per_iteration=502,metric="incidence")
plot_titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
                "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
                "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
                "RW 2nd n 500 Incidence","RW 2nd n 5000 Incidence")
rw_inc_95_plots<-mean_value_plotter(rw_95_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,0.8),metric="incidence")

tot_95_plots<-grid.arrange(splines_inc_95_plots$arranged_plots,rw_inc_95_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 1995")


##### 2000 y2k ####
splines_2000_inc <- list(sp_100_1_00_res,sp_500_1_00_res,sp_1000_1_00_res,sp_5000_1_00_res,
                     sp_100_2_00_res,sp_500_2_00_res,sp_1000_2_00_res,sp_5000_2_00_res)
splines_00_mean_inc<-lapply(splines_2000,mean_value_function,iterations=100,nrow_per_iteration=501,metric="incidence")
plot_titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
                "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
                "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
                "Splines 2nd n 500 Incidence","Splines 2nd n 5000 Incidence")

splines_inc_00_plots<-mean_value_plotter(splines_00_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                          limits_to_y = c(0,0.8),metric="incidence")

#### RW inc ########

rw_2000_inc <- list(rw_100_1_00_res,rw_500_1_00_res,rw_1000_1_00_res,rw_5000_1_00_res,
                rw_100_2_00_res,rw_500_2_00_res,rw_1000_2_00_res,rw_5000_2_00_res)

rw_00_mean_inc<-lapply(rw_2000, mean_value_function,iterations=100,nrow_per_iteration=502,metric="incidence")
plot_titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
                "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
                "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
                "RW 2nd n 500 Incidence","RW 2nd n 5000 Incidence")
rw_inc_00_plots<-mean_value_plotter(rw_00_mean_inc,plot_titles = plot_titulos,simul_datam = sim_model_output$sim_df,
                                     limits_to_y = c(0,0.8),metric="incidence")

tot_00_plots<-grid.arrange(splines_inc_00_plots$arranged_plots,rw_inc_00_plots$arranged_plots,nrow=2,
                           top="Fits produced with Data starting in 2000")

length(rw_00_mean_inc)

##################################################################################################################################################
## Now we'll build a function that displays the 4 different time starts for data on the same graph ###############################################
##################################################################################################################################################

combined_splines<-list(splines_84_mean_prev,splines_90_mean_prev,splines_95_mean_prev,splines_00_mean_prev)


combiner_mean_value_func<-function(list_of_lists,vector_of_years_sampling_starts,metric="prevalence",plot_titles,simul_datam,limits_to_y){
  if(metric == "prevalence"){
    simmed_data<-simul_datam[,c(4,5)]
    
  }
  if(metric == "incidence"){
    simmed_data<-simul_datam[,c(2,5)]
    
  }
  if(metric == "kappa"){
    simmed_data <- simul_datam[,c(1,5)]
  }
  
  simmed_data<-cbind.data.frame(simmed_data[,1],simmed_data[,1],simmed_data[,1],simmed_data[,2],rep("true",nrow(simmed_data)))
  names(simmed_data)<-c("low","median","high","time","date_first_sampling")
  for(i in 1:length(vector_of_years_sampling_starts)){
    for(j in 1:length(list_of_lists[[i]])){
      list_of_lists[[i]][[j]]$date_first_sampling<-rep(vector_of_years_sampling_starts[i],nrow(list_of_lists[[i]][[j]]))
    }
  }
  
  tot_data_list_combined<-vector("list",length = length(combined_splines[[1]]))
  
  for(i in 1:4){
    for(j in 1:length(list_of_lists)){
      tot_data_list_combined[[i]]<- rbind.data.frame(tot_data_list_combined[[i]],list_of_lists[[j]][[i]][1:501,])
      tot_data_list_combined[[i+4]] <- rbind.data.frame(tot_data_list_combined[[i+4]],list_of_lists[[j]][[i+4]][1:501,])
    }
    tot_data_list_combined[[i]]<-rbind.data.frame(tot_data_list_combined[[i]],simmed_data)
    tot_data_list_combined[[i+4]]<-rbind.data.frame(tot_data_list_combined[[i+4]],simmed_data)
  }
  
  tot_plot_list_combined<-vector("list",length(list_of_lists[[1]]))
  
  colour_values<-c("blueviolet","forestgreen","red","dodgerblue","yellow")
  names(colour_values)<-vector_of_years_sampling_starts
  line_values<-rep("solid",length(vector_of_years_sampling_starts))
  names(line_values)<-vector_of_years_sampling_starts
  
  for(i in 1:length(tot_plot_list_combined)){
   tot_plot_list_combined[[i]]<-ggplot(data = tot_data_list_combined[[i]],aes(x=time,y=median,group=date_first_sampling))+
     geom_line(aes(colour=date_first_sampling,linetype=date_first_sampling),size=1.05)+
     geom_ribbon(aes(ymin=low,ymax=high,fill=date_first_sampling),alpha=0.1)+
     # geom_line(data = simmed_data,aes(x=time,y=simmed_data[,1]),colour="red",size=1.04)+
      labs(x="time",y=metric,title=plot_titles[i])+
     scale_colour_manual("Sampling",
                         values = c("true"="black",colour_values))+
     scale_fill_manual("Sampling",
                       values = c("true"="black",colour_values))+
     scale_linetype_manual("Sampling",
                           values = c("true"="dotdash",line_values))+
     coord_cartesian(ylim = limits_to_y)
  
    
  }
  
  arranged_plots<-ggarrange(tot_plot_list_combined[[1]],tot_plot_list_combined[[2]],tot_plot_list_combined[[3]],tot_plot_list_combined[[4]],
                           tot_plot_list_combined[[5]],tot_plot_list_combined[[6]],tot_plot_list_combined[[7]],tot_plot_list_combined[[8]],
                           ncol = 4,nrow = 2)
  
  return(list(plot_list=tot_plot_list_combined,arranged_plot=arranged_plots,df=tot_data_list_combined,dd_df=list_of_lists))
}

titulos<-c("Splines 1st n 100 Prevalence","Splines 1st n 500 Prevalence",
           "Splines 1st n 1000 Prevalence","Splines 1st n 5000 Prevalence",
           "Splines 2nd n 100 Prevalence","Splines 2nd n 500 Prevalence",
           "Splines 2nd n 1000 Prevalence","Splines 2nd n 5000 Prevalence")
years_sampled<-c("1984","1990","1995","2000")
splines_combined<-combiner_mean_value_func(combined_splines,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                             simul_datam = sim_model_output$sim_df,limits_to_y = c(0,80))
splines_combined$arranged_plot


rw_combined<-list(rw_84_mean_prev,rw_90_mean_prev,rw_95_mean_prev,rw_00_mean_prev)
titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
           "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
           "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
           "RW 2nd n 1000 Prevalence","RW 2nd n 5000 Prevalence")
rw_combined_plots<-combiner_mean_value_func(rw_combined,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                            simul_datam = sim_model_output$sim_df,limits_to_y = c(0,80))
rw_combined_plots$arranged_plot


tot_tot_plots<-grid.arrange(splines_combined$arranged_plot,rw_combined_plots$arranged_plot,nrow=2,
                            top="Comparison of fits to prevalence when starting population sampling at different time points")

save(tot_tot_plots,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/different_start_year_analysis/mean_plots/prevalence")

#### now lets do this for incidence 

splines_inc_tot<-list(splines_84_mean_inc,splines_90_mean_inc,splines_95_mean_inc,splines_00_mean_inc)
rw_inc_tot<-list(rw_84_mean_inc,rw_90_mean_inc,rw_95_mean_inc,rw_00_mean_inc)

titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
           "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
           "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
           "Splines 2nd n 1000 Incidence","Splines 2nd n 5000 Incidence")
years_sampled<-c("1984","1990","1995","2000")
splines_combined_inc<-combiner_mean_value_func(splines_inc_tot,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                           simul_datam = sim_model_output$sim_df, metric = "incidence",limits_to_y = c(0,0.5))
splines_combined_inc$arranged_plot

titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
           "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
           "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
           "RW 2nd n 1000 Incidence","RW 2nd n 5000 Incidence")
rw_combined_inc<-combiner_mean_value_func(rw_inc_tot,years_sampled,titulos,sim_model_output$sim_df,metric = "incidence",
                                          limits_to_y = c(0, 0.5))

tot_tot_inc_plot<-grid.arrange(splines_combined_inc$arranged_plot,rw_combined_inc$arranged_plot,nrow=2,
                               top="Comparison of incidence fits when starting population sampling in different years")



####################################################################################################################################
## Let's do this for kappa now #####################################################################################################
####################################################################################################################################

splines_kappa_84<-lapply(splines_1984,mean_value_function,iterations=100,nrow_per_iteration=501, metric = "kappa")
splines_kappa_90 <- lapply(splines_1990, mean_value_function, iterations = 100, nrow_per_iteration = 501, metric = "kappa")
splines_kappa_95 <- lapply(splines_1995, mean_value_function, iterations =100, nrow_per_iteration = 501, metric="kappa")
splines_kappa_00 <- lapply(splines_2000, mean_value_function, iterations =100, nrow_per_iteration = 501, metric="kappa")

rw_kappa_84 <- lapply(rw_1984, mean_value_function, iterations= 100, nrow_per_iteration = 502, metric = "kappa")
rw_kappa_90 <- lapply(rw_1990, mean_value_function, iterations= 100, nrow_per_iteration = 502, metric = "kappa")
rw_kappa_95 <- lapply(rw_1995, mean_value_function, iterations= 100, nrow_per_iteration = 502, metric = "kappa")
rw_kappa_00 <- lapply(rw_2000, mean_value_function, iterations= 100, nrow_per_iteration = 502, metric = "kappa")

splines_kappa_tot <- list(splines_kappa_84, splines_kappa_90, splines_kappa_95, splines_kappa_00)
rw_kappa_tot <- list(rw_kappa_84, rw_kappa_90, rw_kappa_95, rw_kappa_00)

titulos<-c("Splines 1st n 100 kappa","Splines 1st n 500 kappa",
           "Splines 1st n 1000 kappa","Splines 1st n 5000 kappa",
           "Splines 2nd n 100 kappa","Splines 2nd n 500 kappa",
           "Splines 2nd n 1000 kappa","Splines 2nd n 5000 kappa")
years_sampled<-c("1984","1990","1995","2000")
splines_combined_kappa<-combiner_mean_value_func(splines_kappa_tot,vector_of_years_sampling_starts = years_sampled,
                                               plot_titles = titulos, simul_datam = sim_model_output$sim_df,
                                               metric = "kappa",limits_to_y = c(0.05,0.8))
splines_combined_kappa$arranged_plot

titulos<-c("RW 1st n 100 kappa","RW 1st n 500 kappa",
           "RW 1st n 1000 kappa","RW 1st n 5000 kappa",
           "RW 2nd n 100 kappa","RW 2nd n 500 kappa",
           "RW 2nd n 1000 kappa","RW 2nd n 5000 kappa")
rw_combined_kappa<-combiner_mean_value_func(rw_kappa_tot,years_sampled,titulos,sim_model_output$sim_df,metric = "kappa",
                                          limits_to_y = c(0.05, 0.8))

tot_tot_kappa_plot<-grid.arrange(splines_combined_kappa$arranged_plot,rw_combined_kappa$arranged_plot,nrow=2,
                               top="Comparison of kappa fits when starting population sampling in different years")
###################################################################################################################################
## Now lets analyse the different rmse through time plots, and the rmse tables ####################################################
###################################################################################################################################
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

rmse_plotter_and_table_giver<-function(total_list_of_results,year_vector,plot_title){
  year_values<-as.character(year_vector)                                                ## Cretaing a character vector for plot
  tot_vals<-NULL                                                                        ## initialzing dfs
  tot_plot_vals<-NULL
  for(i in 1:length(total_list_of_results)){
    new_line<-c(total_list_of_results[[i]][[1]][[1]],                                   ## extracting results from results list
                total_list_of_results[[i]][[2]][[1]],
                total_list_of_results[[i]][[3]][[1]],
                total_list_of_results[[i]][[4]][[1]],
                total_list_of_results[[i]][[5]][[1]],
                total_list_of_results[[i]][[6]][[1]],
                total_list_of_results[[i]][[7]][[1]],
                total_list_of_results[[i]][[8]][[1]],
                total_list_of_results[[i]][[9]][[1]],
                total_list_of_results[[i]][[10]][[1]],
                total_list_of_results[[i]][[11]][[1]],
                total_list_of_results[[i]][[12]][[1]],
                total_list_of_results[[i]][[13]][[1]],
                total_list_of_results[[i]][[14]][[1]],
                total_list_of_results[[i]][[15]][[1]],
                total_list_of_results[[i]][[16]][[1]])
    new_df<-cbind.data.frame(new_line,rep(year_values[i],length(new_line)))                ## adding in year for plot
    tot_vals<-cbind(tot_vals,new_line)
    tot_plot_vals<-rbind(tot_plot_vals,new_df)
    
  }
  tot_vals<-data.frame(tot_vals)
  row.names(tot_vals)<-paste(c(rep("Spline first",4),rep(" Splines second",4),rep("Rw first",4),rep("RW second",4)),
                               rep(c(100,500,1000,5000),4),sep = " ")
  
  tot_plot_vals<-data.frame(tot_plot_vals)
  names(tot_vals)<-paste("year_sampled",year_values,sep = "_")
  names(tot_plot_vals)<-c("rmse","year_samples")
  tot_plot_vals$sample_size<-rep(c(100,500,1000,5000),nrow(tot_plot_vals)/4)
  tot_plot_vals$model <- rep(c(rep("Spline",8),rep("RW",8)),nrow(tot_plot_vals)/16)
  tot_plot_vals$order<-rep(c(rep("first",4),rep("second",4)),nrow(tot_plot_vals)/8)
  tot_plot_vals$type<-paste(tot_plot_vals$year_samples,tot_plot_vals$model,tot_plot_vals$order)
  tot_plot_vals$order_mod <- paste(tot_plot_vals$model, tot_plot_vals$order,sep = "_")
  
  rmse_plotto<-ggplot(data = tot_plot_vals,aes(x=sample_size,y=rmse,group=type))+
    geom_line(aes(colour=year_samples,linetype=order_mod),size=1.05)+labs(x="Sample size",y="RMSE to true",title=plot_title)+
    geom_point(aes(shape=order,fill=year_samples,colour=year_samples),fill="white",size=3)+
    scale_shape_manual("Order",values = c("first"=21,"second"=23))+
    scale_linetype_manual("Model order",
                          values = c("Spline_first"="solid","Spline_second"="dotdash",
                                     "RW_first"="dotted","RW_second"="longdash"))
  
  
  return(list(df_rmse=tot_vals,plot_df=tot_plot_vals,plotto=rmse_plotto))
  
  
}

###################################################################################################################################
## So that's our functions sorted, now we will get the rmse tables ################################################################
###################################################################################################################################

## 1984
tot_1984_results<-list(sp_100_1_84_res, sp_500_1_84_res, sp_1000_1_84_res, sp_5000_1_84_res,
                       sp_100_2_84_res, sp_500_2_84_res, sp_1000_2_84_res, sp_5000_2_84_res,
                       rw_100_1_84_res, rw_500_1_84_res, rw_1000_1_84_res, rw_5000_1_84_res,
                       rw_100_2_84_res, rw_500_2_84_res, rw_1000_2_84_res, rw_5000_2_84_res)
rmse_1984_prev <- lapply(tot_1984_results, root_mean_error_function, true_data = sim_model_output$sim_df)
rmse_1984_inc <- lapply(tot_1984_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                        metric = "incidence")
rmse_1984_kappa <- lapply(tot_1984_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                          metric= "kappa")
## 1990 
tot_1990_results<-list(sp_100_1_90_res, sp_500_1_90_res, sp_1000_1_90_res, sp_5000_1_90_res,
                       sp_100_2_90_res, sp_500_2_90_res, sp_1000_2_90_res, sp_5000_2_90_res,
                       rw_100_1_90_res, rw_500_1_90_res, rw_1000_1_90_res, rw_5000_1_90_res,
                       rw_100_2_90_res, rw_500_2_90_res, rw_1000_2_90_res, rw_5000_2_90_res)
rmse_1990_prev <- lapply(tot_1990_results, root_mean_error_function, true_data = sim_model_output$sim_df)
rmse_1990_inc <- lapply(tot_1990_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                        metric = "incidence")
rmse_1990_kappa <- lapply(tot_1990_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                          metric= "kappa")

## 1995 
tot_1995_results<-list(sp_100_1_95_res, sp_500_1_95_res, sp_1000_1_95_res, sp_5000_1_95_res,
                       sp_100_2_95_res, sp_500_2_95_res, sp_1000_2_95_res, sp_5000_2_95_res,
                       rw_100_1_95_res, rw_500_1_95_res, rw_1000_1_95_res, rw_5000_1_95_res,
                       rw_100_2_95_res, rw_500_2_95_res, rw_1000_2_95_res, rw_5000_2_95_res)
rmse_1995_prev <- lapply(tot_1995_results, root_mean_error_function, true_data = sim_model_output$sim_df)
rmse_1995_inc <- lapply(tot_1995_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                        metric = "incidence")
rmse_1995_kappa <- lapply(tot_1995_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                          metric= "kappa")

## 2000
tot_2000_results<-list(sp_100_1_00_res, sp_500_1_00_res, sp_1000_1_00_res, sp_5000_1_00_res,
                       sp_100_2_00_res, sp_500_2_00_res, sp_1000_2_00_res, sp_5000_2_00_res,
                       rw_100_1_00_res, rw_500_1_00_res, rw_1000_1_00_res, rw_5000_1_00_res,
                       rw_100_2_00_res, rw_500_2_00_res, rw_1000_2_00_res, rw_5000_2_00_res)
rmse_2000_prev <- lapply(tot_2000_results, root_mean_error_function, true_data = sim_model_output$sim_df)
rmse_2000_inc <- lapply(tot_2000_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                        metric = "incidence")
rmse_2000_kappa <- lapply(tot_2000_results, root_mean_error_function, true_data = sim_model_output$sim_df,
                          metric= "kappa")
#### total results

tot_prev_rmse <- list(rmse_1984_prev, rmse_1990_prev, rmse_1995_prev, rmse_2000_prev)
tot_inc_rmse <- list(rmse_1984_inc, rmse_1990_inc, rmse_1995_inc, rmse_2000_inc)
tot_kappa_rmse <- list(rmse_1984_kappa, rmse_1990_kappa, rmse_1995_kappa, rmse_2000_kappa)

###### lets get our tables and plots out!! ####
year_vecto<-c(1984,1990,1995,2000)

prevalence_vals<-rmse_plotter_and_table_giver(tot_prev_rmse,year_vector = year_vecto,
                                              plot_title = "RMSE of fitted to true over different sampling start dates for Prevalence")
prevalence_vals$df_rmse
prevalence_vals$plotto

incidence_vals <- rmse_plotter_and_table_giver(tot_inc_rmse,year_vector = year_vecto,
                                               plot_title = "RMSE of fitted to true over different sampling start dates for Incidence")
incidence_vals$df_rmse
incidence_vals$plotto

kappa_vals <- rmse_plotter_and_table_giver(tot_kappa_rmse,year_vector = year_vecto,
                                           plot_title = "RMSE of fitted to true over different sampling start dates for Kappa")
kappa_vals$df_rmse
kappa_vals$plotto

save(prevalence_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/different_start_year_analysis/rmse/prevalence")
save(incidence_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/different_start_year_analysis/rmse/incidence")
save(kappa_vals,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/different_start_year_analysis/rmse/kappa")

###################################################################################################################################
## Now we will plot the rmse through time analysis ################################################################################
###################################################################################################################################

rmse_per_time_knotters<-function(true_data,fitted_data,metric="prevalence",year_time_series=seq(1970,2020,0.1)){
  
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
  names(mean_rmse_tot)<-c("RMSE","low","upper","time")

  
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot,error_df=error_tot))
  
  
  
}

plotter_function_rmse<-function(list_of_rmse_results_dfs,plot_title,indiv_year=T,indiv_model=F,
                                sample_size_vec, model_vec, order_vec, year_vec){
  total_data<-NULL
  if(indiv_year == T){
    year_vec <- rep(year_vec,length(list_of_rmse_results_dfs))
    sample_size_vec <- rep(sample_size_vec, length(list_of_rmse_results_dfs)/4)
    model_vec <- c(rep(model_vec[1],length(list_of_rmse_results_dfs)/2),rep(model_vec[2],length(list_of_rmse_results_dfs)/2))
    order_vec <- c(rep(c(rep(order_vec[1],4),rep(order_vec[2],4)),length(list_of_rmse_results_dfs/8)))
  }
  
  if(indiv_model == T){
    year_vec <- c(rep(year_vec[1],list_of_rmse_results_dfs/4),rep(year_vec[2],list_of_rmse_results_dfs/4),
                  rep(year_vec[3],list_of_rmse_results_dfs/4),rep(year_vec[4],list_of_rmse_results_dfs/4))
    sample_size_vec <- rep(sample_size_vec, length(list_of_rmse_results_dfs) / 4)
    model_vec <- 
  }
  
  type <- paste(rep(year_vec,tot_plot_vals$order))
  order_mod <- paste(tot_plot_vals$model, tot_plot_vals$order,sep = "_")
  
    
  for(i in 1:length(list_of_rmse_results_dfs)){
    list_of_rmse_results_dfs[[i]][[2]]$type <- type_vector[i]
    list_of_rmse_results_dfs[[i]][[2]]$year_sampled <- year_sampled_vec[i]
    list_of_rmse_results_dfs[[i]][[2]]$model_order <- model_order_vec[i]
    
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
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=year_sampled,linetype=model_order),size=1.05)+
      labs(x="time",y="RMSE",title=plot_title)
  }
  error_plot<-ggplot(data = total_data,aes(x=time,y=median,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high,fill=type),alpha=0.36)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    #scale_fill_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="error",title=plot_title)
  
  return(list(mean_plot=mean_plot,error_plot=error_plot))
  
  
}

################################################################################################################################
## Lets apply our functions to the list and see what we get for plotting the resulting rmse changes over time ##################
################################################################################################################################

rmse_through_time_1984_prev <- lapply(tot_1984_results, rmse_per_time_knotters, true_data = sim_model_output$sim_df)
