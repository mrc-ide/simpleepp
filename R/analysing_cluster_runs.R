###################################################################################################################################
## Analysis of our runs from the cluster ##########################################################################################
###################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)

###################################################################################################################################
## Lets do some plotting of the mean results over the range of 100 different datasets #############################################
###################################################################################################################################

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

spline_firsty_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_100$prev)

spline_first_n_100<-ggplot(data=spline_firsty_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 100")

spline_firsty_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_500$prev)

spline_first_n_500<-ggplot(data=spline_firsty_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 500")

spline_firsty_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_1000$prev)

spline_first_n_1000<-ggplot(data=spline_firsty_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 1000")

spline_firsty_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_5000$prev)

spline_first_n_5000<-ggplot(data=spline_firsty_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 5000")

first_order_splines<-ggarrange(spline_first_n_100,spline_first_n_500,spline_first_n_1000,spline_first_n_5000,ncol = 2,nrow = 2)
plot(first_order_splines)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_100$prev)

spline_second_n_100<-ggplot(data=spline_secondy_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 100")

spline_secondy_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_500$prev)

spline_second_n_500<-ggplot(data=spline_secondy_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 500")

spline_secondy_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_1000$prev)

spline_second_n_1000<-ggplot(data=spline_secondy_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 1000")

spline_secondy_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_5000$prev)

spline_second_n_5000<-ggplot(data=spline_secondy_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline Second Order n = 5000")

second_order_splines<-ggarrange(spline_second_n_100,spline_second_n_500,
                                spline_second_n_1000,spline_second_n_5000,ncol = 2,nrow = 2)
plot(second_order_splines)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_first_order_n_100$prev)

RW_first_n_100<-ggplot(data=RW_firsty_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 100")

RW_firsty_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_500$prev)

RW_first_n_500<-ggplot(data=RW_firsty_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 500")

RW_firsty_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_1000$prev)

RW_first_n_1000<-ggplot(data=RW_firsty_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 1000")

RW_firsty_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_5000$prev)

RW_first_n_5000<-ggplot(data=RW_firsty_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW First Order n = 5000")

first_order_RW<-ggarrange(RW_first_n_100,RW_first_n_500,RW_first_n_1000,RW_first_n_5000,ncol = 2,nrow = 2)
plot(first_order_RW)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_second_order_n_100$prev)

RW_second_n_100<-ggplot(data=RW_secondy_n_100)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 100")

RW_secondy_n_500<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_500$prev)

RW_second_n_500<-ggplot(data=RW_secondy_n_500)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 500")

RW_secondy_n_1000<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_1000$prev)

RW_second_n_1000<-ggplot(data=RW_secondy_n_1000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
  labs(x="Time",y="Prevalence",title="RW Second Order n = 1000")

RW_secondy_n_5000<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_5000$prev)

RW_second_n_5000<-ggplot(data=RW_secondy_n_5000)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=prev_percent),colour="red")+
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

spline_firsty_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_100$incidence)

spline_first_n_100_inc<-ggplot(data=spline_firsty_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 100 incidence")

spline_firsty_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_500$incidence)

spline_first_n_500_inc<-ggplot(data=spline_firsty_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 500 incidence")

spline_firsty_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_1000$incidence)

spline_first_n_1000_inc<-ggplot(data=spline_firsty_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 1000 incidence")

spline_firsty_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_5000$incidence)

spline_first_n_5000_inc<-ggplot(data=spline_firsty_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 5000 incidence")

first_order_splines_inc<-ggarrange(spline_first_n_100_inc,spline_first_n_500_inc,spline_first_n_1000_inc,spline_first_n_5000_inc,
                                   ncol = 2,nrow = 2)
plot(first_order_splines_inc)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_100$incidence)

spline_second_n_100_inc<-ggplot(data=spline_secondy_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 100 incidence")

spline_secondy_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_500$incidence)

spline_second_n_500_inc<-ggplot(data=spline_secondy_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 500 incidence")

spline_secondy_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_1000$incidence)

spline_second_n_1000_inc<-ggplot(data=spline_secondy_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 1000 incidence")

spline_secondy_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_5000$incidence)

spline_second_n_5000_inc<-ggplot(data=spline_secondy_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="Spline Second Order n = 5000 incidence")

second_order_splines_inc<-ggarrange(spline_second_n_100_inc,spline_second_n_500_inc,
                                spline_second_n_1000_inc,spline_second_n_5000_inc,ncol = 2,nrow = 2)
plot(second_order_splines_inc)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_first_order_n_100$incidence)

RW_first_n_100_inc<-ggplot(data=RW_firsty_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 100 incidence")

RW_firsty_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_500$incidence)

RW_first_n_500_inc<-ggplot(data=RW_firsty_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 500 incidence")

RW_firsty_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_1000$incidence)

RW_first_n_1000_inc<-ggplot(data=RW_firsty_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 1000 incidence")

RW_firsty_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_5000$incidence)

RW_first_n_5000_inc<-ggplot(data=RW_firsty_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW First Order n = 5000 incidence")

first_order_RW_inc<-ggarrange(RW_first_n_100_inc,RW_first_n_500_inc,RW_first_n_1000_inc,RW_first_n_5000_inc,ncol = 2,nrow = 2)
plot(first_order_RW_inc)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_second_order_n_100$incidence)

RW_second_n_100_inc<-ggplot(data=RW_secondy_n_100_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 100 incidence")

RW_secondy_n_500_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_500$incidence)

RW_second_n_500_inc<-ggplot(data=RW_secondy_n_500_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 500 incidence")

RW_secondy_n_1000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_1000$incidence)

RW_second_n_1000_inc<-ggplot(data=RW_secondy_n_1000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 1000 incidence")

RW_secondy_n_5000_inc<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_5000$incidence)

RW_second_n_5000_inc<-ggplot(data=RW_secondy_n_5000_inc)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=lambda),colour="red")+
  labs(x="Time",y="incidence",title="RW Second Order n = 5000 incidence")

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



spline_firsty_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_100$kappa)

spline_first_n_100_kappa<-ggplot(data=spline_firsty_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 100 kappa")

spline_firsty_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_500$kappa)

spline_first_n_500_kappa<-ggplot(data=spline_firsty_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 500 kappa")

spline_firsty_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_1000$kappa)

spline_first_n_1000_kappa<-ggplot(data=spline_firsty_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 1000 kappa")

spline_firsty_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = first_order_spline_n_5000$kappa)

spline_first_n_5000_kappa<-ggplot(data=spline_firsty_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="Prevalence",title="Spline First Order n = 5000 kappa")

first_order_splines_kappa<-ggarrange(spline_first_n_100_kappa,spline_first_n_500_kappa,spline_first_n_1000_kappa,spline_first_n_5000_kappa,
                                   ncol = 2,nrow = 2)
plot(first_order_splines_kappa)

#################################################################################################################################
## Now lets plot the second order splines average fits to the data ##############################################################
#################################################################################################################################

spline_secondy_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_100$kappa)

spline_second_n_100_kappa<-ggplot(data=spline_secondy_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="Spline Second Order n = 100 kappa")

spline_secondy_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_500$kappa)

spline_second_n_500_kappa<-ggplot(data=spline_secondy_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="Spline Second Order n = 500 kappa")

spline_secondy_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_1000$kappa)

spline_second_n_1000_kappa<-ggplot(data=spline_secondy_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="Spline Second Order n = 1000 kappa")

spline_secondy_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 501,data_frame = second_order_spline_n_5000$kappa)

spline_second_n_5000_kappa<-ggplot(data=spline_secondy_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="Spline Second Order n = 5000 kappa")

second_order_splines_kappa<-ggarrange(spline_second_n_100_kappa,spline_second_n_500_kappa,
                                    spline_second_n_1000_kappa,spline_second_n_5000_kappa,ncol = 2,nrow = 2)
plot(second_order_splines_kappa)

################################################################################################################################
## Now lets do random walk #####################################################################################################
################################################################################################################################

RW_firsty_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_first_order_n_100$kappa)

RW_first_n_100_kappa<-ggplot(data=RW_firsty_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW First Order n = 100 kappa")

RW_firsty_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_500$kappa)

RW_first_n_500_kappa<-ggplot(data=RW_firsty_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW First Order n = 500 kappa")

RW_firsty_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_1000$kappa)

RW_first_n_1000_kappa<-ggplot(data=RW_firsty_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW First Order n = 1000 kappa")

RW_firsty_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_first_order_n_5000$kappa)

RW_first_n_5000_kappa<-ggplot(data=RW_firsty_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW First Order n = 5000 kappa")

first_order_RW_kappa<-ggarrange(RW_first_n_100_kappa,RW_first_n_500_kappa,RW_first_n_1000_kappa,RW_first_n_5000_kappa,ncol = 2,nrow = 2)
plot(first_order_RW_kappa)

################################################################################################################################
## Now lets do random walk Second order ########################################################################################
################################################################################################################################
RW_secondy_n_100_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,data_frame = RW_second_order_n_100$kappa)

RW_second_n_100_kappa<-ggplot(data=RW_secondy_n_100_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW Second Order n = 100 kappa")

RW_secondy_n_500_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_500$kappa)

RW_second_n_500_kappa<-ggplot(data=RW_secondy_n_500_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW Second Order n = 500 kappa")

RW_secondy_n_1000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_1000$kappa)

RW_second_n_1000_kappa<-ggplot(data=RW_secondy_n_1000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW Second Order n = 1000 kappa")

RW_secondy_n_5000_kappa<-mean_value_function(iterations = 100,nrow_per_iteration = 502,RW_second_order_n_5000$kappa)

RW_second_n_5000_kappa<-ggplot(data=RW_secondy_n_5000_kappa)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="midnightblue",fill="midnightblue",alpha=0.25)+
  geom_line(data=sim_model_output$sim_df,aes(x=time,y=kappa),colour="red")+
  labs(x="Time",y="kappa",title="RW Second Order n = 5000 kappa")

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
  
  if(metric == "kappa"){
    true_metric <- true_data$kapp
    fitted_metric <- fitted_data$kappa
  }
  
  
 time_to_test<-seq((time_period[1]-1970)*10+1,(time_period[length(time_period)]-2020)*10+501,1)
  
 mean_rmse_tot<-NULL
 
  for (i in 1:100){
    
    
    fitted_metric_iter <- fitted_metric[fitted_metric$iteration == i,]
    
    error <- (fitted_metric_iter$median[time_to_test] /100) - true_metric[time_to_test]
    
    rmse <- sqrt(error^2)
    
    mean_rmse <- mean(rmse)
    iter <- i
    
    mean_rmse <- cbind(mean_rmse,iter)
    
    mean_rmse_tot <- rbind(mean_rmse_tot,mean_rmse)
  }
  
 mean_overall_rmse <- mean(mean_rmse_tot[,1])
 
 
  return(list(mean_rmse=mean_overall_rmse,rmse_df=mean_rmse_tot))
  
  
  
}

rw_first_order_rmse_n_100<-root_mean_error_function(sim_model_output$sim_df,fitted_data = RW_first_order_n_100)
rw_first_order_rmse_n_100_inc<-root_mean_error_function(sim_model_output$sim_df,
                                                        metric = "incidence",fitted_data = RW_first_order_n_100)
rw_first_order_rmse_n_100_inc
rw_first_order_rmse_n_100$mean_rmse
rw_first_order_rmse_n_5000<-root_mean_error_function(sim_model_output$sim_df,fitted_data = RW_first_order_n_5000)
rw_first_order_rmse_n_5000$mean_rmse







