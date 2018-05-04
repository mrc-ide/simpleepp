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




