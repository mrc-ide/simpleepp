## ################################# ##
## RMSE PLots through time ######### ##
## ################################# ##
require(ggpubr)
require(ggplot2)
require(reshape2)

########################################
## Now we will load up the datasets ####
########################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/true_epidemic",verbose = T)

## RW ##

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_100_FIRSt_ORDER", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_500_FIRSt_ORDER", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_1000_FIRSt_ORDER", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_5000_FIRSt_ORDER", verbose = T)


load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_100_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_500_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_1000_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/RW_5000_SECOND_ORDER",verbose = T)

## splines ##

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_100_FIRSt_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_500_FIRSt_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_1000_FIRSt_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_5000_FIRSt_ORDER",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_100_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_500_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_1000_SECOND_ORDER",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/results/SPLINE_5000_SECOND_ORDER",verbose = T)

## ################################# ##
## Now I'll write the function to ## ##
## plot the rmse at each point in ## ##
## time ############################ ##

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


#########################################################################################################################################
## Now lets run this on our dataset #####################################################################################################
#########################################################################################################################################

spline_fir_100_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_1_84_res,
                                   type_of_data = "Spline first 100")
spline_fir_500_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_1_84_res,
                                   type_of_data = "Spline first 500")
spline_fir_1k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_1_84_res,
                                  type_of_data = "Spline first 1000")
spline_fir_5k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_1_84_res,
                                  type_of_data = "Spline first 5000")

spline_sec_100_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_2_84_res,
                                   type_of_data = "Spline second 100")
spline_sec_500_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_2_84_res,
                                   type_of_data = "Spline second 500")
spline_sec_1k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_2_84_res,
                                  type_of_data = "Spline second 1000")
spline_sec_5k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_2_84_res,
                                  type_of_data = "Spline second 5000")

## RW 

RW_fir_100_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_1_84_res,
                               type_of_data = "RW first 100")
RW_fir_500_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_1_84_res,
                               type_of_data = "RW first 500")
RW_fir_1k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_1_84_res,
                              type_of_data = "RW first 1000")
RW_fir_5k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_1_84_res,
                              type_of_data = "RW first 5000")

RW_sec_100_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_2_84_res,
                               type_of_data = "RW second 100")
RW_sec_500_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_2_84_res,
                               type_of_data = "RW second 500")
RW_sec_1k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_2_84_res,
                              type_of_data = "RW second 1000")
RW_sec_5k_rmse<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_2_84_res,
                              type_of_data = "RW second 5000")

#########################################################################################################################################
## First lets look at all the splines together and all the w together then break it down by sample size for a multiplot window ##########
#########################################################################################################################################

splines_comp<-list(spline_fir_100_rmse,spline_fir_500_rmse,spline_fir_1k_rmse,spline_fir_5k_rmse,
                   spline_sec_100_rmse,spline_sec_500_rmse,spline_sec_1k_rmse,spline_sec_5k_rmse)
splies_total_plots<-plotter_function_rmse(splines_comp,plot_title = "Splines comparison",colour_by_sample_size =F)
splies_total_plots$mean_plot #+ coord_cartesian(ylim = c(0,0.36))
splies_total_plots$error_plot
a<-splies_total_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))


rw_comp<-list(RW_fir_100_rmse,RW_fir_500_rmse,RW_fir_1k_rmse,RW_fir_5k_rmse,
              RW_sec_100_rmse,RW_sec_500_rmse,RW_sec_1k_rmse,RW_sec_5k_rmse)
rw_comp_plots<-plotter_function_rmse(rw_comp,plot_title = "RW comparison",colour_by_sample_size = T)
rw_comp_plots$mean_plot #+ coord_cartesian(ylim = c(0,0.36))
rw_comp_plots$error_plot
b<-rw_comp_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))

rw_and_spline_tot_plots<-ggarrange(a, b ,ncol = 1, nrow = 2,align = "hv")
rw_and_spline_tot_plots

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 

save(rw_and_spline_tot_plots,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/rmse_over_course_of_epi/PREVALENCE_spline_and_rw_plot_by_sample_size_OBJECT")


tot_comp<-list(spline_fir_100_rmse,spline_fir_500_rmse,spline_fir_1k_rmse,spline_fir_5k_rmse,
               spline_sec_100_rmse,spline_sec_500_rmse,spline_sec_1k_rmse,spline_sec_5k_rmse,
               RW_fir_100_rmse,RW_fir_500_rmse,RW_fir_1k_rmse,RW_fir_5k_rmse,
               RW_sec_100_rmse,RW_sec_500_rmse,RW_sec_1k_rmse,RW_sec_5k_rmse)
tot_compo_plot<-plotter_function_rmse(tot_comp,plot_title = "Comparion of all methods for prevalence",colour_by_sample_size = T)
tot_compo_plot$mean_plot
tot_compo_plot$error_plot

##### All 100 prevs #####

all_100_dfs<-list(spline_fir_100_rmse,spline_sec_100_rmse,RW_fir_100_rmse,RW_sec_100_rmse)
sample_100_prev_plots<-plotter_function_rmse(all_100_dfs,plot_title = "Comparison of techniques fitting to 100 n for prevalence")
prev_100<-sample_100_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))
sample_100_prev_plots$error_plot

##### All 500 prevs #####

all_500_dfs<-list(spline_fir_500_rmse,spline_sec_500_rmse,RW_fir_500_rmse,RW_sec_500_rmse)
sample_500_prev_plots<-plotter_function_rmse(all_500_dfs,plot_title = "Comparison of techniques fitting to 500 n for prevalence")
prev_500<-sample_500_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))
sample_500_prev_plots$error_plot

##### All 1k prevs #####

all_1k_dfs<-list(spline_fir_1k_rmse,spline_sec_1k_rmse,RW_fir_1k_rmse,RW_sec_1k_rmse)
sample_1k_prev_plots<-plotter_function_rmse(all_1k_dfs, plot_title = "Comparison of techniques fitting to 1000 n for prevalence")
prev_1000<-sample_1k_prev_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))
prev_1000

###### All 5k prevs #####

all_5k_dfs<-list(spline_fir_5k_rmse,spline_sec_5k_rmse,RW_fir_5k_rmse,RW_sec_5k_rmse)
sample_5k_prev_plots<-plotter_function_rmse(all_5k_dfs,plot_title = "Comparison of techniques fitting to 5000 n for prevalence")
prev_5000<-sample_5k_prev_plots$mean_plot + coord_cartesian(ylim=c(0,0.36))

prev_plots_across_sample_sizes<-ggarrange(prev_100,prev_500,prev_1000,prev_5000,nrow = 4,ncol = 1,align = "hv")
prev_plots_across_sample_sizes

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 

save(prev_plots_across_sample_sizes,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/rmse_over_course_of_epi/PREVALENCE_by_sample_size_plot_OBJECT")
prev_plots_across_sample_sizes


########################################################################################################################################
## Now we'll look at results for incidence #############################################################################################
########################################################################################################################################

spline_fir_100_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_1_84_res, metric = "incidence",
                                   type_of_data = "Spline first 100")
spline_fir_500_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_1_84_res,metric = "incidence",
                                   type_of_data = "Spline first 500")
spline_fir_1k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_1_84_res,metric = "incidence",
                                  type_of_data = "Spline first 1000")
spline_fir_5k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_1_84_res,metric = "incidence",
                                  type_of_data = "Spline first 5000")

spline_sec_100_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_2_84_res,metric = "incidence",
                                   type_of_data = "Spline second 100")
spline_sec_500_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_2_84_res,metric = "incidence",
                                   type_of_data = "Spline second 500")
spline_sec_1k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_2_84_res,metric = "incidence",
                                  type_of_data = "Spline second 1000")
spline_sec_5k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_2_84_res,metric = "incidence",
                                  type_of_data = "Spline second 5000")

## RW 

RW_fir_100_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_1_84_res,metric = "incidence",
                               type_of_data = "RW first 100")
RW_fir_500_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_1_84_res,metric = "incidence",
                               type_of_data = "RW first 500")
RW_fir_1k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_1_84_res,metric = "incidence",
                              type_of_data = "RW first 1000")
RW_fir_5k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_1_84_res,metric = "incidence",
                              type_of_data = "RW first 5000")

RW_sec_100_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_2_84_res,metric = "incidence",
                               type_of_data = "RW second 100")
RW_sec_500_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_2_84_res,metric = "incidence",
                               type_of_data = "RW second 500")
RW_sec_1k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_2_84_res,metric = "incidence",
                              type_of_data = "RW second 1000")
RW_sec_5k_rmse_inc<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_2_84_res,metric = "incidence",
                              type_of_data = "RW second 5000")

#########################################################################################################################################
## Now lets get out plots out ###########################################################################################################
#####################################################################545454545454354353153452534534534tr4t43rwef454ty####################

splines_comp_inc<-list(spline_fir_100_rmse_inc,spline_fir_500_rmse_inc,spline_fir_1k_rmse_inc,spline_fir_5k_rmse_inc,
                   spline_sec_100_rmse_inc,spline_sec_500_rmse_inc,spline_sec_1k_rmse_inc,spline_sec_5k_rmse_inc)
splies_total_plots_inc<-plotter_function_rmse(splines_comp_inc,plot_title = "Splines comparison Incidence",colour_by_sample_size = T)
splies_total_plots_inc$mean_plot #+ coord_cartesian(ylim = c(0,0.05))
splies_total_plots_inc$error_plot
a_inc<-splies_total_plots_inc$mean_plot + coord_cartesian(ylim = c(0,0.05))


rw_comp_inc<-list(RW_fir_100_rmse_inc,RW_fir_500_rmse_inc,RW_fir_1k_rmse_inc,RW_fir_5k_rmse_inc,
              RW_sec_100_rmse_inc,RW_sec_500_rmse_inc,RW_sec_1k_rmse_inc,RW_sec_5k_rmse_inc)
rw_comp_plots_inc<-plotter_function_rmse(rw_comp_inc,plot_title = "RW comparison Incidence",colour_by_sample_size = T)
rw_comp_plots_inc$mean_plot #+ coord_cartesian(ylim = c(0,0.05))
rw_comp_plots_inc$error_plot
b_inc<-rw_comp_plots_inc$mean_plot + coord_cartesian(ylim = c(0,0.05))

rw_and_spline_tot_plots_inc<-ggarrange(a_inc, b_inc ,ncol = 1, nrow = 2,align = "hv")
rw_and_spline_tot_plots_inc

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 

save(rw_and_spline_tot_plots_inc,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/rmse_over_course_of_epi/INCIDENCE_rw_and_splines_plot_by__SAMPLE_SIZE_OBJECT")

##### All 100 incs #####

all_100_dfs<-list(spline_fir_100_rmse_inc,spline_sec_100_rmse_inc,RW_fir_100_rmse_inc,RW_sec_100_rmse_inc)
sample_100_inc_plots<-plotter_function_rmse(all_100_dfs,plot_title = "Comparison of techniques fitting to 100 n for incidence")
inc_100<-sample_100_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.05))
sample_100_inc_plots$error_plot

##### All 500 incs #####

all_500_dfs<-list(spline_fir_500_rmse_inc,spline_sec_500_rmse_inc,RW_fir_500_rmse_inc,RW_sec_500_rmse_inc)
sample_500_inc_plots<-plotter_function_rmse(all_500_dfs,plot_title = "Comparison of techniques fitting to 500 n for incidence")
inc_500<-sample_500_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.05))
sample_500_inc_plots$error_plot

##### All 1k incs #####

all_1k_dfs<-list(spline_fir_1k_rmse_inc,spline_sec_1k_rmse_inc,RW_fir_1k_rmse_inc,RW_sec_1k_rmse_inc)
sample_1k_inc_plots<-plotter_function_rmse(all_1k_dfs, plot_title = "Comparison of techniques fitting to 1000 n for incidence")
inc_1000<-sample_1k_inc_plots$mean_plot + coord_cartesian(ylim = c(0,0.05))

###### All 5k incs #####

all_5k_dfs<-list(spline_fir_5k_rmse_inc,spline_sec_5k_rmse_inc,RW_fir_5k_rmse_inc,RW_sec_5k_rmse_inc)
sample_5k_inc_plots<-plotter_function_rmse(all_5k_dfs,plot_title = "Comparison of techniques fitting to 5000 n for incidence")
inc_5000<-sample_5k_inc_plots$mean_plot + coord_cartesian(ylim=c(0,0.05))

inc_plots_across_sample_sizes<-ggarrange(inc_100,inc_500,inc_1000,inc_5000,nrow = 4,ncol = 1,align = "hv")
inc_plots_across_sample_sizes

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 

save(inc_plots_across_sample_sizes,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rmse_over_course_of_epi/INCIDENCE_by_sample_size_plot_object")

#########################################################################################################################################
## NOw we will do the same for the kappa parameter ######################################################################################
#########################################################################################################################################

spline_fir_100_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_1_84_res, metric = "kappa",
                                       type_of_data = "Spline first 100")
spline_fir_500_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_1_84_res,metric = "kappa",
                                       type_of_data = "Spline first 500")
spline_fir_1k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_1_84_res,metric = "kappa",
                                      type_of_data = "Spline first 1000")
spline_fir_5k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_1_84_res,metric = "kappa",
                                      type_of_data = "Spline first 5000")

spline_sec_100_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_100_2_84_res,metric = "kappa",
                                       type_of_data = "Spline second 100")
spline_sec_500_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_500_2_84_res,metric = "kappa",
                                       type_of_data = "Spline second 500")
spline_sec_1k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_1000_2_84_res,metric = "kappa",
                                      type_of_data = "Spline second 1000")
spline_sec_5k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = sp_5000_2_84_res,metric = "kappa",
                                      type_of_data = "Spline second 5000")

## RW 

RW_fir_100_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_1_84_res,metric = "kappa",
                                   type_of_data = "RW first 100")
RW_fir_500_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_1_84_res,metric = "kappa",
                                   type_of_data = "RW first 500")
RW_fir_1k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_1_84_res,metric = "kappa",
                                  type_of_data = "RW first 1000")
RW_fir_5k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_1_84_res,metric = "kappa",
                                  type_of_data = "RW first 5000")

RW_sec_100_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_100_2_84_res,metric = "kappa",
                                   type_of_data = "RW second 100")
RW_sec_500_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_500_2_84_res,metric = "kappa",
                                   type_of_data = "RW second 500")
RW_sec_1k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_1000_2_84_res,metric = "kappa",
                                  type_of_data = "RW second 1000")
RW_sec_5k_rmse_kappa<-rmse_per_time(true_data = sim_model_foi$sim_df,fitted_data = rw_5000_2_84_res,metric = "kappa",
                                  type_of_data = "RW second 5000")

#########################################################################################################################################
## Now lets get out plots out ###########################################################################################################
#####################################################################545454545454354353153452534534534tr4t43rwef454ty####################

splines_comp_kappa<-list(spline_fir_100_rmse_kappa,spline_fir_500_rmse_kappa,spline_fir_1k_rmse_kappa,spline_fir_5k_rmse_kappa,
                       spline_sec_100_rmse_kappa,spline_sec_500_rmse_kappa,spline_sec_1k_rmse_kappa,spline_sec_5k_rmse_kappa)
splies_total_plots_kappa<-plotter_function_rmse(splines_comp_kappa,plot_title = "Splines comparison kappa",colour_by_sample_size = T)
splies_total_plots_kappa$mean_plot #+ coord_cartesian(ylim = c(0,0.36))
splies_total_plots_kappa$error_plot
a_kappa<-splies_total_plots_kappa$mean_plot  + coord_cartesian(ylim = c(0,0.36))


rw_comp_kappa<-list(RW_fir_100_rmse_kappa,RW_fir_500_rmse_kappa,RW_fir_1k_rmse_kappa,RW_fir_5k_rmse_kappa,
                  RW_sec_100_rmse_kappa,RW_sec_500_rmse_kappa,RW_sec_1k_rmse_kappa,RW_sec_5k_rmse_kappa)
rw_comp_plots_kappa<-plotter_function_rmse(rw_comp_kappa,plot_title = "RW comparison kappa",colour_by_sample_size = T)
rw_comp_plots_kappa$mean_plot #+ coord_cartesian(ylim = c(0,0.36))
rw_comp_plots_kappa$error_plot
b_kappa<-rw_comp_plots_kappa$mean_plot + coord_cartesian(ylim = c(0,0.36))

rw_and_spline_tot_plots_kappa<-ggarrange(a_kappa, b_kappa ,ncol = 1, nrow = 2,align = "hv")
rw_and_spline_tot_plots_kappa

tot_plots<-list(spline_fir_100_rmse_kappa,spline_fir_500_rmse_kappa,spline_fir_1k_rmse_kappa,spline_fir_5k_rmse_kappa,
                spline_sec_100_rmse_kappa,spline_sec_500_rmse_kappa,spline_sec_1k_rmse_kappa,spline_sec_5k_rmse_kappa,
                RW_fir_100_rmse_kappa,RW_fir_500_rmse_kappa,RW_fir_1k_rmse_kappa,RW_fir_5k_rmse_kappa,
                RW_sec_100_rmse_kappa,RW_sec_500_rmse_kappa,RW_sec_1k_rmse_kappa,RW_sec_5k_rmse_kappa)
tot_compo_plot_kappa<-plotter_function_rmse(tot_plots,plot_title = "Total kappa comparison",colour_by_sample_size = T)
tot_compo_plot_kappa$error_plot

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 

save(rw_and_spline_tot_plots_kappa,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/rmse_over_course_of_epi/KAPPA_rw_and_spline_plot__colour_by_SAMPLE_SIZE_object")

##### All 100 kappas #####

all_100_dfs<-list(spline_fir_100_rmse_kappa,spline_sec_100_rmse_kappa,RW_fir_100_rmse_kappa,RW_sec_100_rmse_kappa)
sample_100_kappa_plots<-plotter_function_rmse(all_100_dfs,plot_title = "Comparison of techniques fitting to 100 n for kappa")
kappa_100<-sample_100_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))
sample_100_kappa_plots$error_plot

##### All 500 kappas #####

all_500_dfs<-list(spline_fir_500_rmse_kappa,spline_sec_500_rmse_kappa,RW_fir_500_rmse_kappa,RW_sec_500_rmse_kappa)
sample_500_kappa_plots<-plotter_function_rmse(all_500_dfs,plot_title = "Comparison of techniques fitting to 500 n for kappa")
kappa_500<-sample_500_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))
sample_500_kappa_plots$error_plot

##### All 1k kappas #####

all_1k_dfs<-list(spline_fir_1k_rmse_kappa,spline_sec_1k_rmse_kappa,RW_fir_1k_rmse_kappa,RW_sec_1k_rmse_kappa)
sample_1k_kappa_plots<-plotter_function_rmse(all_1k_dfs, plot_title = "Comparison of techniques fitting to 1000 n for kappa")
kappa_1000<-sample_1k_kappa_plots$mean_plot + coord_cartesian(ylim = c(0,0.36))

###### All 5k kappas #####

all_5k_dfs<-list(spline_fir_5k_rmse_kappa,spline_sec_5k_rmse_kappa,RW_fir_5k_rmse_kappa,RW_sec_5k_rmse_kappa)
sample_5k_kappa_plots<-plotter_function_rmse(all_5k_dfs,plot_title = "Comparison of techniques fitting to 5000 n for kappa")
kappa_5000<-sample_5k_kappa_plots$mean_plot + coord_cartesian(ylim=c(0,0.36))

kappa_plots_across_sample_sizes<-ggarrange(kappa_100,kappa_500,kappa_1000,kappa_5000,nrow = 4,ncol = 1,align = "hv")
kappa_plots_across_sample_sizes

###### STOP! IN THE NAME OF LOVE, BEFORE YOU BREAK MY PREVIOUS SAVE LOCATION 


save(kappa_plots_across_sample_sizes,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/data_from_1984_simpleepp/rmse_over_course_of_epi/KAPPA_by_sample_size_plot_OBJECT")






