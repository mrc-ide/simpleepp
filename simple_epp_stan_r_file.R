#' This script will implement a simple EPP model, with CD4 structure and a varying r paramter through time 

################################################################################################################################
## Implement simple EPP structure in R to simulate data, then fit STAN Model to this simulated data ############################
################################################################################################################################

require(ggplot2)
require(reshape2)
require(ggpubr)
require(deSolve)
require(dplyr)
require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################################################################################################################
## Now lets set our initial conditions for the model ###########################################################################
################################################################################################################################

cd4_1<-100                                  ## Initial number of the population infected with 500>cd4>350
cd4_2<-0                                   ## Initial number of the population infected with 350>cd4>200
cd4_3<-0                                   ## Initial number of the population infected with 200>cd4>50
cd4_4<-0                                   ## Initial number of the population infected with cd4<50
s0<-1000000-(cd4_1+cd4_2+cd4_3+cd4_4)      ## Initial susceptible proportion of the population
n0<-s0+cd4_1+cd4_2+cd4_3+cd4_4             ## Initial population size
i0<-(cd4_1+cd4_2+cd4_3+cd4_4)              ## Initial infected size in total
prev0<-i0/n0                               ## Initial Prevalence
incide0<-i0/n0                             ## Initial Incidence 

inits_hiv_cd4<-c(s0,cd4_1,cd4_2,cd4_3,cd4_4,i0,n0,prev0,incide0)

################################################################################################################################
## Now we will set up the intial paramter values for the model #################################################################
################################################################################################################################

eta<- 20000                                ## The birth rate into the population
kappa<- 0.3                                ## The transmission rate paramter
sigma<-c(1/3.16,1/2.13,1/3.2)              ## Progression rates from each Cd4 class
mu_i<-c(0.003,0.008,0.035,0.27)            ## Vector of death rates in each cd4 class         
mu_s<-1/35                                 ## Normal population death rate 

params_cd4_hiv<-list(et=et,kappa=kappa,sigma=sigma,mu_i=mu_i,mu_s=mu_s)

################################################################################################################################
## Now we will create the time series over which to integrate the functions ####################################################
################################################################################################################################

t_min<-0
t_max<-100
times<-t_min:t_max

################################################################################################################################
## Now lets create the system of ODEs to model the change in the population over time ##########################################
################################################################################################################################

epp_cd4<-function(t, y, parms,...) {
  with(as.list(c(params_cd4_hiv, y)), {
    
    dS = eta - ((kappa * y[6] * y[1]) / y[7]) - (mu_s * y[1])
    
    dcd4_1 = ((kappa * y[6] * y[1]) / y[7]) - (mu_i[1] * y[2]) - (sigma[1] * y[2])
    
    dcd4_2 = (sigma[1]*y[2]) - (mu_i[2] * y[3]) - (sigma[2] * y[3])
    
    dcd4_3 = (sigma[2] * y[3]) - (mu_i[3] * y[4]) - (sigma[3] * y[4])
    
    dcd4_4 = (sigma[3] * y[4]) - (mu_i[4] * y[5])
    
    dInfec = dcd4_1 + dcd4_2 +dcd4_3 + dcd4_4
    
    dN = dS + dInfec
    
    dprev = ((dInfec + y[6]) / (dN + y[7])) - y[8]
    
    dincidence = -y[9] + (((kappa * y[6] * y[1]) / y[7])/y[7])
    
    res <- c(dS,dcd4_1,dcd4_2,dcd4_3,dcd4_4,dInfec,dN,dprev,dincidence)
    list(res)
    
    
    })

}


################################################################################################################################
## Now we have written the function we will try to use the ode solver in deSolve to iterate through ############################
################################################################################################################################

out_epp_cd4_hiv <- ode(inits_hiv_cd4, times, epp_cd4, params_cd4_hiv, method="ode45")

out_epp_cd4_hiv<-data.frame(out_epp_cd4_hiv)

names(out_epp_cd4_hiv)<-c("time","susceptible","cd4 > 500","500 > cd4 > 350",
                          "350 > cd4 > 200", "200 > cd4", "infected", "n", "prevalence","incidence")
out_epp_cd4_hiv$prev_percent<-out_epp_cd4_hiv$prevalence * 100

out_epp_cd4_hiv

################################################################################################################################
## Now we have some data we can plot our prevalence curve and incidence curve ##################################################
################################################################################################################################

prev_plot<-ggplot(data = out_epp_cd4_hiv, aes(x=time,y=prev_percent))+geom_line(colour="blue")+
  labs(x="Time",y="Prevalence (%)")

total_infected_plot<-ggplot(data = out_epp_cd4_hiv,aes(x=time,y=incidence))+geom_line(colour="red")+
  labs(x="Time",y="Incidence rate")

total_pop_plot<-ggplot(data = out_epp_cd4_hiv,aes(x=time,y=n))+geom_line(colour="forest green")+
  labs(x="Time",y="Total Population size")

require(gridExtra)
require(grid)

grid::grid.draw(rbind(ggplotGrob(prev_plot), ggplotGrob(total_infected_plot), ggplotGrob(total_pop_plot),  size = "last"))

require(reshape2)
melt_incidence_prev<-out_epp_cd4_hiv[,-c(2:8,11)]
melted_incidence_prevalence<-melt(melt_incidence_prev,id="time")

incidence_prevalence_plot<-ggplot(data = melted_incidence_prevalence)+geom_line(aes(x=time,y=value,colour=variable),size=1.1)+
  labs(x="Time")
plot(incidence_prevalence_plot)
