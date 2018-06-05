###################################################################################################################################
## Analysis of our runs from the cluster ##########################################################################################
###################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)

###################################################################################################################################
## Lets do some plotting of the mean results over the range of 100 different datasets #############################################
###################################################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/true_epidemic",verbose = T)

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

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_first_5k_foi",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/rw_results/rw_sec_5k_foi",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_first_5k_foi",verbose = T)

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_100_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_500_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_1k_foi",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/foi_as_modelled/spline_results/sp_sec_5k_foi",verbose = T)





spline_firsty_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_100_foi_res$prev)

spline_first_n_100<-ggplot(data=spline_firsty_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 100")

spline_firsty_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_500_foi_res$prev)

spline_first_n_500<-ggplot(data=spline_firsty_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 500")

spline_firsty_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_1k_foi_res()$prev)

spline_first_n_1000<-ggplot(data=spline_firsty_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 1000")

spline_firsty_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_5k_foi_res()$prev)

spline_first_n_5000<-ggplot(data=spline_firsty_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 5000")

first_order_splines<-ggarrange(spline_first_n_100,spline_first_n_500,spline_first_n_1000,spline_first_n_5000,ncol = 2,nrow = 2)
plot(first_order_splines)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_100_foi_res()$prev)

spline_second_n_100<-ggplot(data=spline_secondy_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 100")

spline_secondy_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_500_foi_res()$prev)

spline_second_n_500<-ggplot(data=spline_secondy_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 500")

spline_secondy_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_1k_foi_res()$prev)

spline_second_n_1000<-ggplot(data=spline_secondy_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 1000")

spline_secondy_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_5k_foi_res()$prev)

spline_second_n_5000<-ggplot(data=spline_secondy_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 5000")

second_order_splines<-ggarrange(spline_second_n_100,spline_second_n_500,
                                spline_second_n_1000,spline_second_n_5000,ncol = 2,nrow = 2)
plot(second_order_splines)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = rw_f_100_foi_res()$prev)

RW_first_n_100<-ggplot(data=RW_firsty_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 100")

RW_firsty_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_f_500_foi_res()$prev)

RW_first_n_500<-ggplot(data=RW_firsty_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 500")

RW_firsty_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_f_1k_foi_res()$prev)

RW_first_n_1000<-ggplot(data=RW_firsty_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 1000")

RW_firsty_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_f_5k_foi_res()$prev)

RW_first_n_5000<-ggplot(data=RW_firsty_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 5000")

first_order_RW<-ggarrange(RW_first_n_100,RW_first_n_500,RW_first_n_1000,RW_first_n_5000,ncol = 2,nrow = 2)
plot(first_order_RW)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = rw_s_100_foi_res()$prev)

RW_second_n_100<-ggplot(data=RW_secondy_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 100")

RW_secondy_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_s_500_foi_res()$prev)

RW_second_n_500<-ggplot(data=RW_secondy_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 500")

RW_secondy_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_s_1k_foi_res()$prev)

RW_second_n_1000<-ggplot(data=RW_secondy_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 1000")

RW_secondy_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,rw_s_5k_foi_res()$prev)

RW_second_n_5000<-ggplot(data=RW_secondy_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=0.95)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 5000")

second_order_RW<-ggarrange(RW_second_n_100,RW_second_n_500,RW_second_n_1000,RW_second_n_5000,ncol = 2,nrow = 2)
plot(second_order_RW)


total_plots<-ggarrange(first_order_splines,second_order_splines,first_order_RW,second_order_RW,ncol = 2,nrow = 2)
plot(total_plots)

total_plots_by_n<-ggarrange(spline_first_n_100,spline_first_n_500,spline_first_n_1000,spline_first_n_5000,
                            spline_second_n_100,spline_second_n_500,spline_second_n_1000,spline_second_n_5000,
                            RW_first_n_100,RW_first_n_500,RW_first_n_1000,RW_first_n_5000,
                            RW_second_n_100,RW_second_n_500,RW_second_n_1000,RW_second_n_5000,
                            nrow = 4,ncol = 4)
plot(total_plots_by_n)

################################################################################################################################
## So thats our eyeball test of the fit by each method we will now do the same for incidence and transmission parameter ########
################################################################################################################################

spline_firsty_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_100_foi_res()$incidence)

spline_first_n_100_inc<-ggplot(data=spline_firsty_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Incidence",title="Spline First Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))


spline_firsty_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_500_foi_res()$incidence)

spline_first_n_500_inc<-ggplot(data=spline_firsty_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Incidence",title="Spline First Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

spline_firsty_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_1k_foi_res()$incidence)

spline_first_n_1000_inc<-ggplot(data=spline_firsty_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Incidence",title="Spline First Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

spline_firsty_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_5k_foi_res() $incidence)

spline_first_n_5000_inc<-ggplot(data=spline_firsty_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Incidence",title="Spline First Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

first_order_splines_inc<-ggarrange(spline_first_n_100_inc,spline_first_n_500_inc,spline_first_n_1000_inc,spline_first_n_5000_inc,
                                   ncol = 2,nrow = 2)
plot(first_order_splines_inc)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_100_foi_res()$incidence)

spline_second_n_100_inc<-ggplot(data=spline_secondy_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

spline_secondy_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_500_foi_res()$incidence)

spline_second_n_500_inc<-ggplot(data=spline_secondy_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

spline_secondy_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_1k_foi_res()$incidence)

spline_second_n_1000_inc<-ggplot(data=spline_secondy_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

spline_secondy_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_5k_foi_res()$incidence)

spline_second_n_5000_inc<-ggplot(data=spline_secondy_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

second_order_splines_inc<-ggarrange(spline_second_n_100_inc,spline_second_n_500_inc,
                                    spline_second_n_1000_inc,spline_second_n_5000_inc,ncol = 2,nrow = 2)
plot(second_order_splines_inc)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = rw_f_100_foi_res()$incidence)

RW_first_n_100_inc<-ggplot(data=RW_firsty_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))


RW_firsty_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_500_foi_res()$incidence)

RW_first_n_500_inc<-ggplot(data=RW_firsty_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

RW_firsty_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_1k_foi_res()$incidence)

RW_first_n_1000_inc<-ggplot(data=RW_firsty_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

RW_firsty_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_5k_foi_res()$incidence)

RW_first_n_5000_inc<-ggplot(data=RW_firsty_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

first_order_RW_inc<-ggarrange(RW_first_n_100_inc,RW_first_n_500_inc,RW_first_n_1000_inc,RW_first_n_5000_inc,ncol = 2,nrow = 2)
plot(first_order_RW_inc)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = rw_s_100_foi_res()$incidence)

RW_second_n_100_inc<-ggplot(data=RW_secondy_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

RW_secondy_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_500_foi_res()$incidence)

RW_second_n_500_inc<-ggplot(data=RW_secondy_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

RW_secondy_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_1k_foi_res()$incidence)

RW_second_n_1000_inc<-ggplot(data=RW_secondy_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

RW_secondy_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_5k_foi_res()$incidence)

RW_second_n_5000_inc<-ggplot(data=RW_secondy_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0,0.7))

second_order_RW_inc<-ggarrange(RW_second_n_100_inc,RW_second_n_500_inc,RW_second_n_1000_inc,RW_second_n_5000_inc,
                               ncol = 2,nrow = 2)
plot(second_order_RW_inc)


total_plots<-ggarrange(first_order_splines,second_order_splines,first_order_RW,second_order_RW,ncol = 2,nrow = 2)
plot(total_plots)

total_plots_by_n_inc<-ggarrange(spline_first_n_100_inc,spline_first_n_500_inc,spline_first_n_1000_inc,spline_first_n_5000_inc,
                                spline_second_n_100_inc,spline_second_n_500_inc,spline_second_n_1000_inc,spline_second_n_5000_inc,
                                RW_first_n_100_inc,RW_first_n_500_inc,RW_first_n_1000_inc,RW_first_n_5000_inc,
                                RW_second_n_100_inc,RW_second_n_500_inc,RW_second_n_1000_inc,RW_second_n_5000_inc,
                                nrow = 4,ncol = 4)
plot(total_plots_by_n_inc)


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~~~!~~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#|-|-|-|-|-|-|-|-|-|-|/-\|

################################################################################################################################
## Now lets go for the kappa parameter and see what fit that has to the data ###################################################
################################################################################################################################



spline_firsty_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_100_foi_res()$kappa)

spline_first_n_100_kappa<-ggplot(data=spline_firsty_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline First Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_firsty_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_500_foi_res()$kappa)

spline_first_n_500_kappa<-ggplot(data=spline_firsty_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline First Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_firsty_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_1k_foi_res()$kappa)

spline_first_n_1000_kappa<-ggplot(data=spline_firsty_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline First Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_firsty_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_f_5k_foi_res() $kappa)

spline_first_n_5000_kappa<-ggplot(data=spline_firsty_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline First Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

first_order_splines_kappa<-ggarrange(spline_first_n_100_kappa,spline_first_n_500_kappa,spline_first_n_1000_kappa,spline_first_n_5000_kappa,
                                     ncol = 2,nrow = 2)
plot(first_order_splines_kappa)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_100_foi_res()$kappa)

spline_second_n_100_kappa<-ggplot(data=spline_secondy_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline Second Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_secondy_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_500_foi_res()$kappa)

spline_second_n_500_kappa<-ggplot(data=spline_secondy_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline Second Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_secondy_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_1k_foi_res()$kappa)

spline_second_n_1000_kappa<-ggplot(data=spline_secondy_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline Second Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

spline_secondy_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = spline_s_5k_foi_res()$kappa)

spline_second_n_5000_kappa<-ggplot(data=spline_secondy_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="Spline Second Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

second_order_splines_kappa<-ggarrange(spline_second_n_100_kappa,spline_second_n_500_kappa,
                                      spline_second_n_1000_kappa,spline_second_n_5000_kappa,ncol = 2,nrow = 2)
plot(second_order_splines_kappa)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = rw_f_100_foi_res()$kappa)

RW_first_n_100_kappa<-ggplot(data=RW_firsty_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW First Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

RW_firsty_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_500_foi_res()$kappa)

RW_first_n_500_kappa<-ggplot(data=RW_firsty_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW First Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

RW_firsty_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_1k_foi_res()$kappa)

RW_first_n_1000_kappa<-ggplot(data=RW_firsty_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW First Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

RW_firsty_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_f_5k_foi_res()$kappa)

RW_first_n_5000_kappa<-ggplot(data=RW_firsty_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW First Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.0,0.7))

first_order_RW_kappa<-ggarrange(RW_first_n_100_kappa,RW_first_n_500_kappa,RW_first_n_1000_kappa,RW_first_n_5000_kappa,ncol = 2,nrow = 2)
plot(first_order_RW_kappa)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = rw_s_100_foi_res()$kappa)

RW_second_n_100_kappa<-ggplot(data=RW_secondy_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW Second Order n = 100")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.3,1))

RW_secondy_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_500_foi_res()$kappa)

RW_second_n_500_kappa<-ggplot(data=RW_secondy_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW Second Order n = 500")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.3,1))

RW_secondy_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_1k_foi_res()$kappa)

RW_second_n_1000_kappa<-ggplot(data=RW_secondy_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW Second Order n = 1000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.3,1))

RW_secondy_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,rw_s_5k_foi_res()$kappa)

RW_second_n_5000_kappa<-ggplot(data=RW_secondy_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_foi$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Kappa",title="RW Second Order n = 5000")+
  coord_cartesian(xlim = c(1970,2020),ylim = c(0.3,1))

second_order_RW_kappa<-ggarrange(RW_second_n_100_kappa,RW_second_n_500_kappa,RW_second_n_1000_kappa,RW_second_n_5000_kappa,
                                 ncol = 2,nrow = 2)
plot(second_order_RW_kappa)


total_plots<-ggarrange(first_order_splines,second_order_splines,first_order_RW,second_order_RW,ncol = 2,nrow = 2)
plot(total_plots)

total_plots_by_n_kappa<-ggarrange(spline_first_n_100_kappa,spline_first_n_500_kappa,spline_first_n_1000_kappa,spline_first_n_5000_kappa,
                                  spline_second_n_100_kappa,spline_second_n_500_kappa,spline_second_n_1000_kappa,spline_second_n_5000_kappa,
                                  RW_first_n_100_kappa,RW_first_n_500_kappa,RW_first_n_1000_kappa,RW_first_n_5000_kappa,
                                  RW_second_n_100_kappa,RW_second_n_500_kappa,RW_second_n_1000_kappa,RW_second_n_5000_kappa,
                                  nrow = 4,ncol = 4)
plot(total_plots_by_n_kappa)


################################################################################################################################
## So that's our eyeball tests performed for each metric we're interested in, now we need to do some more formal testing of ####
## the association between the produced lines and the actual lines #############################################################
################################################################################################################################

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


##############################################################################################################################
## So that's our function for evaluating the RMSE of the fitts, we will now go through each fitted data set and get the RMSE #
## for each of the fitts to our data #########################################################################################
##############################################################################################################################

RMSE_full_data_spline_first_order_n_100_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_f_100_foi_res())
RMSE_full_data_spline_first_order_n_100_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      fitted_data = spline_f_100_foi_res())
RMSE_full_data_spline_first_order_n_100_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        fitted_data = spline_f_100_foi_res())
spline_first_order_analysis_n_100<-list(rmse_prev=RMSE_full_data_spline_first_order_n_100_prev,
                                  rmse_inc=RMSE_full_data_spline_first_order_n_100_inc,
                                  rmse_kappa=RMSE_full_data_spline_first_order_n_100_kappa,
                                  prev_mean_plot=spline_first_n_100,inc_mean_plot=spline_first_n_100_inc,
                                  kappa_mean_plot=spline_first_n_100_kappa)

RMSE_full_data_spline_first_order_n_500_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_f_500_foi_res())
RMSE_full_data_spline_first_order_n_500_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      spline_f_500_foi_res())
RMSE_full_data_spline_first_order_n_500_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        spline_f_500_foi_res())

spline_first_order_analysis_n_500<-list(rmse_prev=RMSE_full_data_spline_first_order_n_500_prev,
                                        rmse_inc=RMSE_full_data_spline_first_order_n_500_inc,
                                        rmse_kappa=RMSE_full_data_spline_first_order_n_500_kappa,
                                        prev_mean_plot=spline_first_n_500,inc_mean_plot=spline_first_n_500_inc,
                                        kappa_mean_plot=spline_first_n_500_kappa)



RMSE_full_data_spline_first_order_n_1000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_f_1k_foi_res())
RMSE_full_data_spline_first_order_n_1000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      fitted_data = spline_f_1k_foi_res())
RMSE_full_data_spline_first_order_n_1000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        fitted_data = spline_f_1k_foi_res())
spline_first_order_analysis_n_1000<-list(rmse_prev=RMSE_full_data_spline_first_order_n_1000_prev,
                                        rmse_inc=RMSE_full_data_spline_first_order_n_1000_inc,
                                        rmse_kappa=RMSE_full_data_spline_first_order_n_1000_kappa,
                                        prev_mean_plot=spline_first_n_1000,inc_mean_plot=spline_first_n_1000_inc,
                                        kappa_mean_plot=spline_first_n_1000_kappa)


RMSE_full_data_spline_first_order_n_5000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_f_5k_foi_res() )
RMSE_full_data_spline_first_order_n_5000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      fitted_data = spline_f_5k_foi_res() )
RMSE_full_data_spline_first_order_n_5000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        fitted_data = spline_f_5k_foi_res() )
spline_first_order_analysis_n_5000<-list(rmse_prev=RMSE_full_data_spline_first_order_n_5000_prev,
                                        rmse_inc=RMSE_full_data_spline_first_order_n_5000_inc,
                                        rmse_kappa=RMSE_full_data_spline_first_order_n_5000_kappa,
                                        prev_mean_plot=spline_first_n_5000,inc_mean_plot=spline_first_n_5000_inc,
                                        kappa_mean_plot=spline_first_n_5000_kappa)

################################################################################################################################
## So thats the first order splines done, we'll now move on to the second order splines complete data ##########################
################################################################################################################################

RMSE_full_data_spline_second_order_n_100_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_s_100_foi_res())
RMSE_full_data_spline_second_order_n_100_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      fitted_data = spline_s_100_foi_res())
RMSE_full_data_spline_second_order_n_100_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        fitted_data = spline_s_100_foi_res())
spline_second_order_analysis_n_100<-list(rmse_prev=RMSE_full_data_spline_second_order_n_100_prev,
                                        rmse_inc=RMSE_full_data_spline_second_order_n_100_inc,
                                        rmse_kappa=RMSE_full_data_spline_second_order_n_100_kappa,
                                        prev_mean_plot=spline_second_n_100,inc_mean_plot=spline_second_n_100_inc,
                                        kappa_mean_plot=spline_second_n_100_kappa)

RMSE_full_data_spline_second_order_n_500_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                       spline_s_500_foi_res())
RMSE_full_data_spline_second_order_n_500_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                      spline_s_500_foi_res())
RMSE_full_data_spline_second_order_n_500_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                        spline_s_500_foi_res())

spline_second_order_analysis_n_500<-list(rmse_prev=RMSE_full_data_spline_second_order_n_500_prev,
                                        rmse_inc=RMSE_full_data_spline_second_order_n_500_inc,
                                        rmse_kappa=RMSE_full_data_spline_second_order_n_500_kappa,
                                        prev_mean_plot=spline_second_n_500,inc_mean_plot=spline_second_n_500_inc,
                                        kappa_mean_plot=spline_second_n_500_kappa)



RMSE_full_data_spline_second_order_n_1000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                        spline_s_1k_foi_res())
RMSE_full_data_spline_second_order_n_1000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                       fitted_data = spline_s_1k_foi_res())
RMSE_full_data_spline_second_order_n_1000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                         fitted_data = spline_s_1k_foi_res())
spline_second_order_analysis_n_1000<-list(rmse_prev=RMSE_full_data_spline_second_order_n_1000_prev,
                                         rmse_inc=RMSE_full_data_spline_second_order_n_1000_inc,
                                         rmse_kappa=RMSE_full_data_spline_second_order_n_1000_kappa,
                                         prev_mean_plot=spline_second_n_1000,inc_mean_plot=spline_second_n_1000_inc,
                                         kappa_mean_plot=spline_second_n_1000_kappa)


RMSE_full_data_spline_second_order_n_5000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                        spline_s_5k_foi_res())
RMSE_full_data_spline_second_order_n_5000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                       fitted_data = spline_s_5k_foi_res())
RMSE_full_data_spline_second_order_n_5000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                         fitted_data = spline_s_5k_foi_res())
spline_second_order_analysis_n_5000<-list(rmse_prev=RMSE_full_data_spline_second_order_n_5000_prev,
                                         rmse_inc=RMSE_full_data_spline_second_order_n_5000_inc,
                                         rmse_kappa=RMSE_full_data_spline_second_order_n_5000_kappa,
                                         prev_mean_plot=spline_second_n_5000,inc_mean_plot=spline_second_n_5000_inc,
                                         kappa_mean_plot=spline_second_n_5000_kappa)


first_order_col<-c(spline_first_order_analysis_n_100$rmse_prev$mean_rmse,
                   spline_first_order_analysis_n_500$rmse_prev$mean_rmse,
                   spline_first_order_analysis_n_1000$rmse_prev$mean_rmse,
                   spline_first_order_analysis_n_5000$rmse_prev$mean_rmse)

second_order_col<-c(spline_second_order_analysis_n_100$rmse_prev$mean_rmse,
                    spline_second_order_analysis_n_500$rmse_prev$mean_rmse,
                    spline_second_order_analysis_n_1000$rmse_prev$mean_rmse,
                    spline_second_order_analysis_n_5000$rmse_prev$mean_rmse)

mean_rmse_data_frame<-cbind.data.frame(first_order_col,second_order_col)
names(mean_rmse_data_frame)<-c("Spline First Order","Spline Second Order")


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::~@~@~@~~~?@@?@{}{@:>:<:@~@{>{}]]]]}}

################################################################################################################################
## So that's splines taken care of we can now look at random walks #############################################################
################################################################################################################################

RMSE_full_data_RW_first_order_n_100_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                          rw_f_100_foi_res())
RMSE_full_data_RW_first_order_n_100_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                       fitted_data = rw_f_100_foi_res())
RMSE_full_data_RW_first_order_n_100_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                         fitted_data = rw_f_100_foi_res())
RW_first_order_analysis_n_100<-list(rmse_prev=RMSE_full_data_RW_first_order_n_100_prev,
                                         rmse_inc=RMSE_full_data_RW_first_order_n_100_inc,
                                         rmse_kappa=RMSE_full_data_RW_first_order_n_100_kappa,
                                         prev_mean_plot=RW_first_n_100,inc_mean_plot=RW_first_n_100_inc,
                                         kappa_mean_plot=RW_first_n_100_kappa)

RMSE_full_data_RW_first_order_n_500_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                        rw_f_500_foi_res())
RMSE_full_data_RW_first_order_n_500_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                       rw_f_500_foi_res())
RMSE_full_data_rw_f_500_foi_res_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                         rw_f_500_foi_res())

RW_first_order_analysis_n_500<-list(rmse_prev=RMSE_full_data_RW_first_order_n_500_prev,
                                         rmse_inc=RMSE_full_data_RW_first_order_n_500_inc,
                                         rmse_kappa=RMSE_full_data_rw_f_500_foi_res_kappa,
                                         prev_mean_plot=RW_first_n_500,inc_mean_plot=RW_first_n_500_inc,
                                         kappa_mean_plot=RW_first_n_500_kappa)



RMSE_full_data_RW_first_order_n_1000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                         rw_f_1k_foi_res())
RMSE_full_data_RW_first_order_n_1000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                        fitted_data = rw_f_1k_foi_res())
RMSE_full_data_RW_first_order_n_1000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                          fitted_data = rw_f_1k_foi_res())
RW_first_order_analysis_n_1000<-list(rmse_prev=RMSE_full_data_RW_first_order_n_1000_prev,
                                          rmse_inc=RMSE_full_data_RW_first_order_n_1000_inc,
                                          rmse_kappa=RMSE_full_data_RW_first_order_n_1000_kappa,
                                          prev_mean_plot=RW_first_n_1000,inc_mean_plot=RW_first_n_1000_inc,
                                          kappa_mean_plot=RW_first_n_1000_kappa)


RMSE_full_data_RW_first_order_n_5000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                         rw_f_5k_foi_res())
RMSE_full_data_RW_first_order_n_5000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                        fitted_data = rw_f_5k_foi_res())
RMSE_full_data_RW_first_order_n_5000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                          fitted_data = rw_f_5k_foi_res())
RW_first_order_analysis_n_5000<-list(rmse_prev=RMSE_full_data_RW_first_order_n_5000_prev,
                                          rmse_inc=RMSE_full_data_RW_first_order_n_5000_inc,
                                          rmse_kappa=RMSE_full_data_RW_first_order_n_5000_kappa,
                                          prev_mean_plot=RW_first_n_5000,inc_mean_plot=RW_first_n_5000_inc,
                                          kappa_mean_plot=RW_first_n_5000_kappa)


## Second order RW now 

RMSE_full_data_RW_second_order_n_100_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                   rw_s_100_foi_res())
RMSE_full_data_RW_second_order_n_100_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                  fitted_data = rw_s_100_foi_res())
RMSE_full_data_RW_second_order_n_100_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                    fitted_data = rw_s_100_foi_res())
RW_second_order_analysis_n_100<-list(rmse_prev=RMSE_full_data_RW_second_order_n_100_prev,
                                    rmse_inc=RMSE_full_data_RW_second_order_n_100_inc,
                                    rmse_kappa=RMSE_full_data_RW_second_order_n_100_kappa,
                                    prev_mean_plot=RW_second_n_100,inc_mean_plot=RW_second_n_100_inc,
                                    kappa_mean_plot=RW_second_n_100_kappa)

RMSE_full_data_RW_second_order_n_500_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                   rw_s_500_foi_res())
RMSE_full_data_RW_second_order_n_500_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                  rw_s_500_foi_res())
RMSE_full_data_RW_second_order_n_500_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                    rw_s_500_foi_res())

RW_second_order_analysis_n_500<-list(rmse_prev=RMSE_full_data_RW_second_order_n_500_prev,
                                    rmse_inc=RMSE_full_data_RW_second_order_n_500_inc,
                                    rmse_kappa=RMSE_full_data_RW_second_order_n_500_kappa,
                                    prev_mean_plot=RW_second_n_500,inc_mean_plot=RW_second_n_500_inc,
                                    kappa_mean_plot=RW_second_n_500_kappa)



RMSE_full_data_RW_second_order_n_1000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                    rw_s_1k_foi_res())
RMSE_full_data_RW_second_order_n_1000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                   fitted_data = rw_s_1k_foi_res())
RMSE_full_data_RW_second_order_n_1000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                     fitted_data = rw_s_1k_foi_res())
RW_second_order_analysis_n_1000<-list(rmse_prev=RMSE_full_data_RW_second_order_n_1000_prev,
                                     rmse_inc=RMSE_full_data_RW_second_order_n_1000_inc,
                                     rmse_kappa=RMSE_full_data_RW_second_order_n_1000_kappa,
                                     prev_mean_plot=RW_second_n_1000,inc_mean_plot=RW_second_n_1000_inc,
                                     kappa_mean_plot=RW_second_n_1000_kappa)


RMSE_full_data_RW_second_order_n_5000_prev<-root_mean_error_function(sim_model_foi$sim_df,
                                                                    rw_s_5k_foi_res())
RMSE_full_data_RW_second_order_n_5000_inc<-root_mean_error_function(sim_model_foi$sim_df,metric = "incidence",
                                                                   fitted_data = rw_s_5k_foi_res())
RMSE_full_data_RW_second_order_n_5000_kappa<-root_mean_error_function(sim_model_foi$sim_df,metric = "kappa",
                                                                     fitted_data = rw_s_5k_foi_res())
RW_second_order_analysis_n_5000<-list(rmse_prev=RMSE_full_data_RW_second_order_n_5000_prev,
                                     rmse_inc=RMSE_full_data_RW_second_order_n_5000_inc,
                                     rmse_kappa=RMSE_full_data_RW_second_order_n_5000_kappa,
                                     prev_mean_plot=RW_second_n_5000,inc_mean_plot=RW_second_n_5000_inc,
                                     kappa_mean_plot=RW_second_n_5000_kappa)





first_order_col<-c(RW_first_order_analysis_n_100$rmse_prev$mean_rmse,
                   RW_first_order_analysis_n_500$rmse_prev$mean_rmse,
                   RW_first_order_analysis_n_1000$rmse_prev$mean_rmse,
                   RW_first_order_analysis_n_5000$rmse_prev$mean_rmse)

second_order_col<-c(RW_second_order_analysis_n_100$rmse_prev$mean_rmse,
                    RW_second_order_analysis_n_500$rmse_prev$mean_rmse,
                    RW_second_order_analysis_n_1000$rmse_prev$mean_rmse,
                    RW_second_order_analysis_n_5000$rmse_prev$mean_rmse)

mean_rmse_data_frame<-cbind.data.frame(mean_rmse_data_frame,first_order_col,second_order_col)
names(mean_rmse_data_frame)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")




overall_fitting_analysis<-list(sp_first_100=spline_first_order_analysis_n_100,sp_first_500=spline_first_order_analysis_n_500,
                               sp_first_1000=spline_first_order_analysis_n_1000,sp_first_5k=spline_first_order_analysis_n_5000,
                               sp_sec_100=spline_second_order_analysis_n_100,sp_sec_500=spline_second_order_analysis_n_500,
                               sp_sec_1k=spline_second_order_analysis_n_1000,sp_sec_5k=spline_second_order_analysis_n_5000,
                               rw_first_100=RW_first_order_analysis_n_100,rw_first_500=RW_first_order_analysis_n_500,
                               rw_first_1k=RW_first_order_analysis_n_1000,rw_first_5k=RW_first_order_analysis_n_5000,
                               rw_sec_100=RW_second_order_analysis_n_100,rw_sec_500=RW_second_order_analysis_n_500,
                               rw_sec_1k=RW_second_order_analysis_n_1000,rw_sec_5k=RW_second_order_analysis_n_5000,
                               mean_rmse_tot=mean_rmse_data_frame)


################################################################################################################################
## Going to ease this into creating a single function for all four datasets in one #############################################
################################################################################################################################


RMSE_dataset_extraction<-function(overall_analysis_list=overall_fitting_analysis,metric="prevalence"){

if(metric == "prevalence"){
  first_col<-c(overall_analysis_list[[1]][[1]][[1]],
               overall_analysis_list[[2]][[1]][[1]],
               overall_analysis_list[[3]][[1]][[1]],
               overall_analysis_list[[4]][[1]][[1]])
  
  sec_col<-c(overall_analysis_list[[5]][[1]][[1]],
             overall_analysis_list[[6]][[1]][[1]],
             overall_analysis_list[[7]][[1]][[1]],
             overall_analysis_list[[8]][[1]][[1]])
  
  third_col<-c(overall_analysis_list[[9]][[1]][[1]],
               overall_analysis_list[[10]][[1]][[1]],
               overall_analysis_list[[11]][[1]][[1]],
               overall_analysis_list[[12]][[1]][[1]])
  
  fourth_col<-c(overall_analysis_list[[13]][[1]][[1]],
                overall_analysis_list[[14]][[1]][[1]],
                overall_analysis_list[[15]][[1]][[1]],
                overall_analysis_list[[16]][[1]][[1]])
  
  df_analysis<-cbind.data.frame(first_col,sec_col,third_col,fourth_col)
  names(df_analysis)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
}

  if(metric == "incidence"){
    first_col<-c(overall_analysis_list[[1]][[2]][[1]],
                 overall_analysis_list[[2]][[2]][[1]],
                 overall_analysis_list[[3]][[2]][[1]],
                 overall_analysis_list[[4]][[2]][[1]])
    
    sec_col<-c(overall_analysis_list[[5]][[2]][[1]],
               overall_analysis_list[[6]][[2]][[1]],
               overall_analysis_list[[7]][[2]][[1]],
               overall_analysis_list[[8]][[2]][[1]])
    
    third_col<-c(overall_analysis_list[[9]][[2]][[1]],
                 overall_analysis_list[[10]][[2]][[1]],
                 overall_analysis_list[[11]][[2]][[1]],
                 overall_analysis_list[[12]][[2]][[1]])
    
    fourth_col<-c(overall_analysis_list[[13]][[2]][[1]],
                  overall_analysis_list[[14]][[2]][[1]],
                  overall_analysis_list[[15]][[2]][[1]],
                  overall_analysis_list[[16]][[2]][[1]])
    
    df_analysis<-cbind.data.frame(first_col,sec_col,third_col,fourth_col)
    names(df_analysis)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
  }
    
  if(metric == "kappa"){
    
    first_col<-c(overall_analysis_list[[1]][[3]][[1]],
                 overall_analysis_list[[2]][[3]][[1]],
                 overall_analysis_list[[3]][[3]][[1]],
                 overall_analysis_list[[4]][[3]][[1]])
    
    sec_col<-c(overall_analysis_list[[5]][[3]][[1]],
               overall_analysis_list[[6]][[3]][[1]],
               overall_analysis_list[[7]][[3]][[1]],
               overall_analysis_list[[8]][[3]][[1]])
    
    third_col<-c(overall_analysis_list[[9]][[3]][[1]],
                 overall_analysis_list[[10]][[3]][[1]],
                 overall_analysis_list[[11]][[3]][[1]],
                 overall_analysis_list[[12]][[3]][[1]])
    
    fourth_col<-c(overall_analysis_list[[13]][[3]][[1]],
                  overall_analysis_list[[14]][[3]][[1]],
                  overall_analysis_list[[15]][[3]][[1]],
                  overall_analysis_list[[16]][[3]][[1]])
    
    df_analysis<-cbind.data.frame(first_col,sec_col,third_col,fourth_col)
    names(df_analysis)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
  }
  
  return(df_analysis)
}

prev_mean_complete_df<-RMSE_dataset_extraction()

inc_mean_complete_df<-RMSE_dataset_extraction(metric = "incidence")

kappa_mean_complete_df<-RMSE_dataset_extraction(metric = "kappa")

complete_rmse_data<-list(prev=prev_mean_complete_df,inc=inc_mean_complete_df,kappa=kappa_mean_complete_df)

save(complete_rmse_data,file = "hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/RMSE_analysis/prev_inc_kappa_average_RMSE")

overall_fitting_analysis$mean_rmse_tot
complete_rmse_data$prev
complete_rmse_data$inc
complete_rmse_data$kappa


cool_plot_function<-function(base_vector,simmed_df,simmed_vector_col_id,simmed_iteration){
plot(base_vector)
for(i in 1:100){
  colour="yellow"
  if( ((i+4)/5) == round((i+4)/5)){
    colour="red"
  }
  if(((i+3)/5) == round((i+3)/5)){
    colour="forestgreen"
  }
  if(((i+2)/5) == round((i+2)/5)){
    colour="dodgerblue"
  }
  if(((i+1)/5) == round((i+1)/5)){
    colour="orange"
  }
  
  lines((simmed_df[simmed_iteration==i,simmed_vector_col_id] ),col=colour)
  Sys.sleep(0.3)
  print(i)
}
}
cool_plot_function(sim_model_foi$sim_df$kappa,simmed_df = spline_f_100_foi_res()$kappa,
                   simmed_vector_col_id = 2,simmed_iteration = spline_f_100_foi_res()$kappa$iteration)
cool_plot_function(sim_model_foi$sim_df$lambda,simmed_df = rw_f_1k_foi_res()$incidence,
                   simmed_vector_col_id = 2, simmed_iteration = rw_f_1k_foi_res()$incidence$iteration)
cool_plot_function(sim_model_foi$sim_df$prev_percent,simmed_df = spline_f_100_foi_res$prev,
                   simmed_vector_col_id = 1, simmed_iteration = spline_f_100_foi_res$prev$iteration)

################################################################################################################################
## Now lets come up with a function that can do all the above rmse stuff, just for the prediction period #######################
################################################################################################################################

overall_fitted<-list(spline_first_100=spline_f_100_foi_res(),spline_first_500=spline_f_500_foi_res(),
                     spline_first_1k=spline_f_1k_foi_res(),spline_first_5k=spline_f_5k_foi_res() ,
                     spline_second_100=spline_s_100_foi_res(),spline_second_500=spline_s_500_foi_res(),
                     spline_second_1000=spline_s_1k_foi_res(),spline_second_5000=spline_s_5k_foi_res(),
                     RW_first_100=rw_f_100_foi_res(),RW_first_500=rw_f_500_foi_res(),
                     RW_first_1k=rw_f_1k_foi_res(),RW_first_5k=rw_f_5k_foi_res(),
                     RW_second_100=rw_s_100_foi_res(),RW_second_500=rw_s_500_foi_res(),
                     RW_second_1k=rw_s_1k_foi_res(),RW_second_5k=rw_s_5k_foi_res())


getting_overall_rmse_for_data<-function(list_of_fitted_outputs,true_df,time_period_to_test_over=seq(1970,2020,0.1)){
  
  prev_col<-NULL
  inc_col <- NULL
  kappa_col <- NULL
  
  for(i in 1:length(list_of_fitted_outputs)){
  
  sp_1_100<-root_mean_error_function(true_data = true_df,fitted_data = list_of_fitted_outputs[[i]],
                                     time_period = time_period_to_test_over)
  sp_1_100_inc<-root_mean_error_function(true_data = true_df,fitted_data = list_of_fitted_outputs[[i]],
                                         metric = "incidence",time_period = time_period_to_test_over)
  sp_1_100_kappa<-root_mean_error_function(true_data = true_df,fitted_data = list_of_fitted_outputs[[i]],
                                           metric = "kappa",time_period = time_period_to_test_over)
  
  prev_col<-c(prev_col,sp_1_100[[1]])
  inc_col<-c(inc_col,sp_1_100_inc[[1]])
  kappa_col<-c(kappa_col,sp_1_100_kappa[[1]])
  test_stat<-ncol(sp_1_100[[3]])
  }

  spline_first_prev<-prev_col[1:4]
  spline_second_prev<-prev_col[5:8]
  rw_first_prev<-prev_col[9:12]
  rw_sec_prev<-prev_col[13:16]
  
  spline_first_inc<-inc_col[1:4]
  spline_sec_inc<-inc_col[5:8]
  rw_first_inc<-inc_col[9:12]
  rw_sec_inc<-inc_col[13:16]
  
  spline_first_kappa<-kappa_col[1:4]
  spline_sec_kappa<-kappa_col[5:8]
  rw_first_kappa<-kappa_col[9:12]
  rw_sec_kappa<-kappa_col[13:16]
  
  prev_df<-cbind.data.frame(spline_first_prev,spline_second_prev,rw_first_prev,rw_sec_prev)
  names(prev_df)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
  prev_df$n<-c(100,500,1000,5000)
  
  
  inc_df<-cbind.data.frame(spline_first_inc,spline_sec_inc,rw_first_inc,rw_sec_inc)
  names(inc_df)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
  inc_df$n<-c(100,500,1000,5000)
  
  kappa_df<-cbind.data.frame(spline_first_kappa,spline_sec_kappa,rw_first_kappa,rw_sec_kappa)
  names(kappa_df)<-c("Spline First Order","Spline Second Order","RW first order","RW second order")
  kappa_df$n<-c(100,500,1000,5000)
  
  return(list(prevalence=prev_df,incidence=inc_df,kappa=kappa_df,test_stat=test_stat))
}

prediction_period_rmse<-getting_overall_rmse_for_data(list_of_fitted_outputs = overall_fitted,true_df = sim_model_foi$sim_df,
                                                      time_period_to_test_over = seq(2015.1,2020,0.1))
prediction_period_rmse$prevalence
prediction_period_rmse$incidence
prediction_period_rmse$kappa

save(prediction_period_rmse,file = "hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/RMSE_analysis/prediction_period_rmse_dp")

plot(overall_fitting_analysis$sp_first_100$inc_mean_plot)

peak_epidemic_period<-getting_overall_rmse_for_data(list_of_fitted_outputs = overall_fitted,true_df = sim_model_foi$sim_df,
                                                    time_period_to_test_over = seq(1990,2000,0.1))
peak_epidemic_period$prevalence
peak_epidemic_period$incidence
peak_epidemic_period$kappa

save(peak_epidemic_period,file = "hiv_project/analysis_of_cluster_run_datasets/double_peak_simple_epp/RMSE_analysis/peak_epidemic_rmse_dp")

###############################################################################################################################
## Now we will plot a graph of the mean error for each iteration for each of the 100 datasets #################################
###############################################################################################################################

error_df_prev<-cbind.data.frame(overall_fitting_analysis$sp_first_100$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_first_500$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_first_1000$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_first_5k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_sec_100$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_sec_500$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_sec_1k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$sp_sec_5k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_first_100$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_first_500$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_first_1k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_first_5k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_sec_100$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_sec_500$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_sec_1k$rmse_prev$rmse_df[,1],
                                overall_fitting_analysis$rw_sec_5k$rmse_prev$rmse_df[,1])
names(error_df_prev)<-c("spline_first_100","spline_first_500","spline_first_1k","spline_first_5k",
                        "spline_sec_100","spline_sec_500","spline_sec_1k","spline_sec_5k",
                        "rw_first_100","rw_first_500","rw_first_1k","rw_first_5k",
                        "rw_sec_100","rw_sec_500","rw_sec_1k","rw_sec_5k")
error_df_prev$iteration<-seq(1,100,1)

melted_data_prev_error<-melt(error_df_prev,id="iteration")

error_plot<-ggplot(data = melted_data_prev_error)+geom_line(aes(x=iteration,y=value,colour=variable),size=1)+
  labs(x="iteration",y="RMSE value",title="RMSE over 100 different datasets for prevalence")+
  coord_cartesian(ylim = c())

error_df_prev_100<-cbind.data.frame(overall_fitting_analysis$sp_first_100$rmse_prev$rmse_df[,1],
                                    overall_fitting_analysis$sp_sec_100$rmse_prev$rmse_df[,1],
                                    overall_fitting_analysis$rw_first_100$rmse_prev$rmse_df[,1],
                                    overall_fitting_analysis$rw_sec_100$rmse_prev$rmse_df[,1])
names(error_df_prev_100)<-c("spline_first_100",
                        "spline_sec_100",
                        "rw_first_100",
                        "rw_sec_100")
error_df_prev_100$iteration<-seq(1,100,1)

melted_prev_100<-melt(error_df_prev_100,id="iteration")

error_plot_prev_100<-ggplot(data = melted_prev_100)+geom_line(aes(x=iteration,y=value,colour=variable),size=1)+
  labs(x="iteration",y="RMSE value",title="RMSE over 100 different datasets for prevalence")+
  coord_cartesian(ylim = c())



mean_error<-sqrt(mean((spline_firsty_n_100_inc$median[200:301] - sim_model_foi$sim_df$lambda[200:301])^2))  
mean_ezza<-sqrt(mean((RW_firsty_n_100_inc$median[200:301] - sim_model_foi$sim_df$lambda[200:301])^2))


