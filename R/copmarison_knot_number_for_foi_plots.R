#############################################################################################################################################
## Analyzing the multiple knot techniques used when modelling foi as a parameter ############################################################
#############################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)

#############################################################################################################################################
## loading up the required datasets #########################################################################################################
#############################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/true_epidemic",verbose = T)

## 7 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_5k_foi",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_5k_foi",verbose = T)

## 8 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/spline_100_1st_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/spline_500_1st_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/spline_1000_1st_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/spline_5000_1st_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/sp_sec_100_foi_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/sp_sec_500_foi_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/sp_sec_1000_foi_8_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/8_knot_splines/results/sp_sec_5000_foi_8_knot",
     verbose = T)

## 9 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/spline_100_1st_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/spline_500_1st_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/spline_1000_1st_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/spline_5000_1st_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/sp_sec_100_foi_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/sp_sec_500_foi_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/sp_sec_1000_foi_9_knot",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/9_knot_splines/results/sp_sec_5000_foi_9_knot",
     verbose = T)

## 10 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/spline_100_1st_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/spline_500_1st_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/spline_1000_1st_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/spline_5000_1st_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/sp_sec_100_foi_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/sp_sec_500_foi_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/sp_sec_1000_foi_10_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/10_knots/results/sp_sec_5000_foi_10_knots",
     verbose = T)

### 11 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/spline_100_1st_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/spline_500_1st_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/spline_1000_1st_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/spline_5000_1st_11_knots",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/sp_sec_100_foi_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/sp_sec_500_foi_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/sp_sec_1000_foi_11_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/11_knots/results/sp_sec_5000_foi_11_knots",
     verbose = T)

### 12 knotters

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/spline_100_1st_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/spline_500_1st_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/spline_1000_1st_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/spline_5000_1st_12_knots",
     verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/sp_sec_100_foi_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/sp_sec_500_foi_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/sp_sec_1000_foi_12_knots",
     verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/12_knots/results/sp_sec_5000_foi_12_knots",
     verbose = T)

#######################################################################################################################################
## Now lets write a function that outputs a graph comparing these different outcomes over differernt scales ###########################
#######################################################################################################################################

knotter_compo_after_2_hits<-function(list_of_total_data_sets,graph_titles,true_df,true_col_id,include_credible = T){
  
  mean_value_function<-function(iterations,nrow_per_iteration,data_frame){
    data_prev<-NULL
    iter_value<-iterations - 1
    nrow_value<-nrow_per_iteration
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
  
  tot_df_values<-NULL
  knot_vals<-NULL
  for (i in 1:length(list_of_total_data_sets)){
    mean_df<-mean_value_function(100,501,data_frame = list_of_total_data_sets[[i]])
    
    mean_df$knot<-rep(as.character(i+6),nrow(mean_df))
    
    tot_df_values<-rbind.data.frame(tot_df_values,mean_df)
    
    knot_vals<-c(knot_vals,mean_df$knot[1])
    
  }
  
  plot_colours<-c("dodgerblue","red","blueviolet","forestgreen","yellow",
                  "indianred4","springgreen2","chocolate4","azure2","midnightblue")
  
  random_colour<-round(runif(1,min = 1,max = length(plot_colours)))
  
  true_df_binding<-cbind.data.frame(true_col_id,true_col_id,true_col_id,true_df$time,rep("true",nrow(true_df)))
  names(true_df_binding)<-c("low","median","high","time","knot")
  tot_df_values<-rbind.data.frame(tot_df_values,true_df_binding)
  
  colour_values<-NULL
  
  for(i in 1:length(list_of_total_data_sets)){
    colour_values<-c(colour_values,plot_colours[i])
  }
  
  names(colour_values)<-knot_vals
  size_values<-rep(1.05,length(list_of_total_data_sets))
  names(size_values)<-knot_vals
  
  if(include_credible == T){
 total_plot<-ggplot(data = tot_df_values, aes(x=tot_df_values$time,y=tot_df_values$median,group=tot_df_values$knot))+
   geom_line(aes(colour=tot_df_values$knot,size=tot_df_values$knot))+
   geom_ribbon(aes(x=tot_df_values$time,
                   ymin=tot_df_values$low,
                   ymax=tot_df_values$high,
                   fill=tot_df_values$knot),
                   #colour=tot_df_values$knot),
               alpha= 0.2)+
   #geom_line(data = true_df,aes(x=true_df$time,y=true_col_id),size=1,colour=plot_colours[random_colour])+
   labs(x="time",y=graph_titles[1],title=graph_titles[2])+
   scale_colour_manual("Knot Number",
                       values = c("true"="black",colour_values))+
   scale_fill_manual("Knot Number",
                     values = c("true"="black",colour_values))+
   scale_size_manual("Knot Number",
                     values = c("true"=1.5,size_values))

  }else{
    total_plot<-ggplot(data = tot_df_values, aes(x=tot_df_values$time,y=tot_df_values$median,group=tot_df_values$knot))+
      geom_line(aes(colour=tot_df_values$knot,size=tot_df_values$knot)) + labs(x="time",y=graph_titles[1],title=graph_titles[2])+
      scale_colour_manual("Knot Number",
                          values = c("true"="black",colour_values))+
      scale_size_manual("Knot Number",
                        values = c("true"=1.5,size_values))
    
      
      
    }
  
  
  
  return(list(df=tot_df_values,compo_plot=total_plot,col_vals=colour_values))
  
}

#######################################################################################################################################
## We'll compare the results for prevalence now #######################################################################################
#######################################################################################################################################

graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 100")
df_tot<-list(spline_f_100_foi_res$prev,spline_f_100_foi_res_8$prev,spline_f_100_foi_res_9$prev,spline_f_100_foi_res_ten$prev,
             spline_f_100_foi_res_11$prev,spline_f_100_foi_res_12$prev)

n_100_knot_compo<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                              include_credible = T)

## 500
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 500")
df_tot<-list(spline_f_500_foi_res$prev,spline_f_500_foi_res_8$prev,spline_f_500_foi_res_9$prev,spline_f_500_foi_res_ten$prev,
             spline_f_500_foi_res_11$prev,spline_f_500_foi_res_12$prev)

n_500_knot_compo<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                             include_credible = T)

## 1000
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 1000")
df_tot<-list(spline_f_1k_foi_res$prev,spline_f_1k_foi_res_8$prev,spline_f_1k_foi_res_9$prev,spline_f_1k_foi_res_ten$prev,
             spline_f_1k_foi_res_11$prev,spline_f_1k_foi_res_12$prev)

n_1k_knot_compo<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                             include_credible = T)

## 5000
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 5000")
df_tot<-list(spline_f_5k_foi_res$prev,spline_f_5k_foi_res_8$prev,spline_f_5k_foi_res_9$prev,spline_f_5k_foi_res_ten$prev,
             spline_f_5k_foi_res_11$prev,spline_f_5k_foi_res_12$prev)

n_5k_knot_compo<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                            include_credible = T)


prev_compo_plots<-ggarrange(n_100_knot_compo$compo_plot,n_500_knot_compo$compo_plot,n_1k_knot_compo$compo_plot,
                            n_5k_knot_compo$compo_plot,ncol = 2,nrow = 2)
prev_compo_plots

##' 2nd order ####
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 100 2nd Order")
df_tot<-list(spline_s_100_foi_res$prev,spline_s_100_foi_res_8$prev,spline_s_100_foi_res_9$prev,spline_s_100_foi_res_ten$prev,
             spline_s_100_foi_res_11$prev,spline_s_100_foi_res_12$prev)

n_100_knot_compo_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                             include_credible = T)

## 500
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 500 2nd Order")
df_tot<-list(spline_s_500_foi_res$prev,spline_s_500_foi_res_8$prev,spline_s_500_foi_res_9$prev,spline_s_500_foi_res_ten$prev,
             spline_s_500_foi_res_11$prev,spline_s_500_foi_res_12$prev)

n_500_knot_compo_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                             include_credible = T)

## 1000
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 1000 2nd Order")
df_tot<-list(spline_s_1k_foi_res$prev,spline_s_1k_foi_res_8$prev,spline_s_1k_foi_res_9$prev,spline_s_1k_foi_res_ten$prev,
             spline_s_1k_foi_res_11$prev,spline_s_1k_foi_res_12$prev)

n_1k_knot_compo_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                            include_credible = T)

## 5000
graph_title<-c("prev","Comparison of number of knots predicting prevalence n = 5000 2nd Order")
df_tot<-list(spline_s_5k_foi_res$prev,spline_s_5k_foi_res_8$prev,spline_s_5k_foi_res_9$prev,spline_s_5k_foi_res_ten$prev,
             spline_s_5k_foi_res_11$prev,spline_s_5k_foi_res_12$prev)

n_5k_knot_compo_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$prev_percent,
                                            include_credible = T)


prev_compo_plots_sec<-ggarrange(n_100_knot_compo_sec$compo_plot,n_500_knot_compo_sec$compo_plot,
                                n_1k_knot_compo_sec$compo_plot,n_5k_knot_compo_sec$compo_plot,ncol = 2,nrow = 2)
prev_compo_plots_sec

prev_1_and_2_prev_compo<-ggarrange(n_100_knot_compo$compo_plot,n_100_knot_compo_sec$compo_plot,
                                   n_500_knot_compo$compo_plot,n_500_knot_compo_sec$compo_plot,
                                   n_1k_knot_compo$compo_plot,n_1k_knot_compo_sec$compo_plot,
                                   n_5k_knot_compo$compo_plot,n_5k_knot_compo_sec$compo_plot,
                                   ncol = 2,nrow = 4)

########################################################################################################################################
## INC now #############################################################################################################################
########################################################################################################################################
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 100")
df_tot<-list(spline_f_100_foi_res$incidence,spline_f_100_foi_res_8$incidence,spline_f_100_foi_res_9$incidence,
             spline_f_100_foi_res_ten$incidence,spline_f_100_foi_res_11$incidence,spline_f_100_foi_res_12$incidence)

n_100_knot_compo_inc<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                             include_credible = T)

## 500
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 500")
df_tot<-list(spline_f_500_foi_res$incidence,spline_f_500_foi_res_8$incidence,spline_f_500_foi_res_9$incidence,
             spline_f_500_foi_res_ten$incidence,spline_f_500_foi_res_11$incidence,spline_f_500_foi_res_12$incidence)

n_500_knot_compo_inc<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                             include_credible = T)

## 1000
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 1000")
df_tot<-list(spline_f_1k_foi_res$incidence,spline_f_1k_foi_res_8$incidence,spline_f_1k_foi_res_9$incidence,
             spline_f_1k_foi_res_ten$incidence,spline_f_1k_foi_res_11$incidence,spline_f_1k_foi_res_12$incidence)

n_1k_knot_compo_inc<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                            include_credible = T)

## 5000
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 5000")
df_tot<-list(spline_f_5k_foi_res$incidence,spline_f_5k_foi_res_8$incidence,spline_f_5k_foi_res_9$incidence,
             spline_f_5k_foi_res_ten$incidence,spline_f_5k_foi_res_11$incidence,spline_f_5k_foi_res_12$incidence)

n_5k_knot_compo_inc<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                            include_credible = T)


inc_compo_inc_plots<-ggarrange(n_100_knot_compo_inc$compo_plot,n_500_knot_compo_inc$compo_plot,n_1k_knot_compo_inc$compo_plot,
                            n_5k_knot_compo_inc$compo_plot,ncol = 2,nrow = 2)
inc_compo_inc_plots

##' 2nd order ####
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 100 2nd Order")
df_tot<-list(spline_s_100_foi_res$incidence,spline_s_100_foi_res_8$incidence,spline_s_100_foi_res_9$incidence,spline_s_100_foi_res_ten$incidence,
             spline_s_100_foi_res_11$incidence,spline_s_100_foi_res_12$incidence)

n_100_knot_compo_inc_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                                 include_credible = T)

## 500
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 500 2nd Order")
df_tot<-list(spline_s_500_foi_res$incidence,spline_s_500_foi_res_8$incidence,spline_s_500_foi_res_9$incidence,spline_s_500_foi_res_ten$incidence,
             spline_s_500_foi_res_11$incidence,spline_s_500_foi_res_12$incidence)

n_500_knot_compo_inc_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                                 include_credible = T)

## 1000
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 1000 2nd Order")
df_tot<-list(spline_s_1k_foi_res$incidence,spline_s_1k_foi_res_8$incidence,spline_s_1k_foi_res_9$incidence,spline_s_1k_foi_res_ten$incidence,
             spline_s_1k_foi_res_11$incidence,spline_s_1k_foi_res_12$incidence)

n_1k_knot_compo_inc_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                                include_credible = T)

## 5000
graph_title<-c("inc","Comparison of number of knots predicting incidence n = 5000 2nd Order")
df_tot<-list(spline_s_5k_foi_res$incidence,spline_s_5k_foi_res_8$incidence,spline_s_5k_foi_res_9$incidence,spline_s_5k_foi_res_ten$incidence,
             spline_s_5k_foi_res_11$incidence,spline_s_5k_foi_res_12$incidence)

n_5k_knot_compo_inc_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$lambda,
                                                include_credible = T)


inc_compo_plots_sec<-ggarrange(n_100_knot_compo_inc_sec$compo_plot,n_500_knot_compo_inc_sec$compo_plot,
                                n_1k_knot_compo_inc_sec$compo_plot,n_5k_knot_compo_inc_sec$compo_plot,ncol = 2,nrow = 2)
inc_compo_plots_sec

inc_1_and_2_inc_compo<-ggarrange(n_100_knot_compo_inc$compo_plot,n_100_knot_compo_inc_sec$compo_plot,
                                   n_500_knot_compo_inc$compo_plot,n_500_knot_compo_inc_sec$compo_plot,
                                   n_1k_knot_compo_inc$compo_plot,n_1k_knot_compo_inc_sec$compo_plot,
                                   n_5k_knot_compo_inc$compo_plot,n_5k_knot_compo_inc_sec$compo_plot,
                                   ncol = 2,nrow = 4)
inc_1_and_2_inc_compo

#########################################################################################################################################
## Now for the kappa terms ##############################################################################################################
###### kappa kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  #####
###### kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa kappa ####
###### kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kapp kappa ####
##### kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa  kappa ######
#########################################################################################################################################

graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 100")
df_tot<-list(spline_f_100_foi_res$kappa,spline_f_100_foi_res_8$kappa,spline_f_100_foi_res_9$kappa,spline_f_100_foi_res_ten$kappa,
             spline_f_100_foi_res_11$kappa,spline_f_100_foi_res_12$kappa)

n_100_knot_compo_kappa<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                 include_credible = T)

## 500
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 500")
df_tot<-list(spline_f_500_foi_res$kappa,spline_f_500_foi_res_8$kappa,spline_f_500_foi_res_9$kappa,spline_f_500_foi_res_ten$kappa,
             spline_f_500_foi_res_11$kappa,spline_f_500_foi_res_12$kappa)

n_500_knot_compo_kappa<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                 include_credible = T)

## 1000
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 1000")
df_tot<-list(spline_f_1k_foi_res$kappa,spline_f_1k_foi_res_8$kappa,spline_f_1k_foi_res_9$kappa,spline_f_1k_foi_res_ten$kappa,
             spline_f_1k_foi_res_11$kappa,spline_f_1k_foi_res_12$kappa)

n_1k_knot_compo_kappa<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                include_credible = T)

## 5000
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 5000")
df_tot<-list(spline_f_5k_foi_res$kappa,spline_f_5k_foi_res_8$kappa,spline_f_5k_foi_res_9$kappa,spline_f_5k_foi_res_ten$kappa,
             spline_f_5k_foi_res_11$kappa,spline_f_5k_foi_res_12$kappa)

n_5k_knot_compo_kappa<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                include_credible = T)


inc_compo_kappa_plots<-ggarrange(n_100_knot_compo_kappa$compo_plot,n_500_knot_compo_kappa$compo_plot,n_1k_knot_compo_kappa$compo_plot,
                               n_5k_knot_compo_kappa$compo_plot,ncol = 2,nrow = 2)
inc_compo_kappa_plots

##' 2nd order ####
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 100 2nd Order")
df_tot<-list(spline_s_100_foi_res$kappa,spline_s_100_foi_res_8$kappa,spline_s_100_foi_res_9$kappa,spline_s_100_foi_res_ten$kappa,
             spline_s_100_foi_res_11$kappa,spline_s_100_foi_res_12$kappa)

n_100_knot_compo_kappa_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                     include_credible = T)

## 500
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 500 2nd Order")
df_tot<-list(spline_s_500_foi_res$kappa,spline_s_500_foi_res_8$kappa,spline_s_500_foi_res_9$kappa,spline_s_500_foi_res_ten$kappa,
             spline_s_500_foi_res_11$kappa,spline_s_500_foi_res_12$kappa)

n_500_knot_compo_kappa_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                     include_credible = T)

## 1000
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 1000 2nd Order")
df_tot<-list(spline_s_1k_foi_res$kappa,spline_s_1k_foi_res_8$kappa,spline_s_1k_foi_res_9$kappa,spline_s_1k_foi_res_ten$kappa,
             spline_s_1k_foi_res_11$kappa,spline_s_1k_foi_res_12$kappa)

n_1k_knot_compo_kappa_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                    include_credible = T)

## 5000
graph_title<-c("kappa","Comparison of number of knots predicting incidence n = 5000 2nd Order")
df_tot<-list(spline_s_5k_foi_res$kappa,spline_s_5k_foi_res_8$kappa,spline_s_5k_foi_res_9$kappa,spline_s_5k_foi_res_ten$kappa,
             spline_s_5k_foi_res_11$kappa,spline_s_5k_foi_res_12$kappa)

n_5k_knot_compo_kappa_sec<-knotter_compo_after_2_hits(df_tot,graph_title,true_df = sim_model_foi$sim_df,true_col_id = sim_model_foi$sim_df$kappa,
                                                    include_credible = T)


kappa_compo_plots_sec<-ggarrange(n_100_knot_compo_kappa_sec$compo_plot,n_500_knot_compo_kappa_sec$compo_plot,
                               n_1k_knot_compo_kappa_sec$compo_plot,n_5k_knot_compo_kappa_sec$compo_plot,ncol = 2,nrow = 2)
kappa_compo_plots_sec

kappa_1_and_2_kappa_compo<-ggarrange(n_100_knot_compo_kappa$compo_plot,n_100_knot_compo_kappa_sec$compo_plot,
                                 n_500_knot_compo_kappa$compo_plot,n_500_knot_compo_kappa_sec$compo_plot,
                                 n_1k_knot_compo_kappa$compo_plot,n_1k_knot_compo_kappa_sec$compo_plot,
                                 n_5k_knot_compo_kappa$compo_plot,n_5k_knot_compo_kappa_sec$compo_plot,
                                 ncol = 2,nrow = 4)
kappa_1_and_2_kappa_compo


