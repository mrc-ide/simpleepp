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
           "Splines 2nd n 500 Prevalence","Splines 2nd n 5000 Prevalence")
years_sampled<-c("1984","1990","1995","2000")
splines_combined<-combiner_mean_value_func(combined_splines,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                             simul_datam = sim_model_output$sim_df,limits_to_y = c(0,80))
splines_combined$arranged_plot


rw_combined<-list(rw_84_mean_prev,rw_90_mean_prev,rw_95_mean_prev,rw_00_mean_prev)
titulos<-c("RW 1st n 100 Prevalence","RW 1st n 500 Prevalence",
           "RW 1st n 1000 Prevalence","RW 1st n 5000 Prevalence",
           "RW 2nd n 100 Prevalence","RW 2nd n 500 Prevalence",
           "RW 2nd n 500 Prevalence","RW 2nd n 5000 Prevalence")
rw_combined_plots<-combiner_mean_value_func(rw_combined,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                            simul_datam = sim_model_output$sim_df,limits_to_y = c(0,80))
rw_combined_plots$arranged_plot


tot_tot_plots<-grid.arrange(splines_combined$arranged_plot,rw_combined_plots$arranged_plot,nrow=2,
                            top="Comparison of fits to prevalence when starting population sampling at different time points")


#### now lets do this for incidence 

splines_inc_tot<-list(splines_84_mean_inc,splines_90_mean_inc,splines_95_mean_inc,splines_00_mean_inc)
rw_inc_tot<-list(rw_84_mean_inc,rw_90_mean_inc,rw_95_mean_inc,rw_00_mean_inc)

titulos<-c("Splines 1st n 100 Incidence","Splines 1st n 500 Incidence",
           "Splines 1st n 1000 Incidence","Splines 1st n 5000 Incidence",
           "Splines 2nd n 100 Incidence","Splines 2nd n 500 Incidence",
           "Splines 2nd n 500 Incidence","Splines 2nd n 5000 Incidence")
years_sampled<-c("1984","1990","1995","2000")
splines_combined_inc<-combiner_mean_value_func(splines_inc_tot,vector_of_years_sampling_starts = years_sampled,plot_titles = titulos,
                                           simul_datam = sim_model_output$sim_df, metric = "incidence")
splines_combined_inc$arranged_plot

titulos<-c("RW 1st n 100 Incidence","RW 1st n 500 Incidence",
           "RW 1st n 1000 Incidence","RW 1st n 5000 Incidence",
           "RW 2nd n 100 Incidence","RW 2nd n 500 Incidence",
           "RW 2nd n 500 Incidence","RW 2nd n 5000 Incidence")
rw_combined_inc<-combiner_mean_value_func(rw_inc_tot,years_sampled,titulos,sim_model_output$sim_df,metric = "incidence")

tot_tot_inc_plot<-grid.arrange(splines_combined_inc$arranged_plot,rw_combined_inc$arranged_plot,nrow=2,
                               top="Comparison of incidence fits when starting population sampling in different years")





