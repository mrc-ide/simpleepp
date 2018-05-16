###################################################################################################################################
## Testing out how the new art diseased model work ################################################################################
###################################################################################################################################

require(ggplot2)
require(rstan)
require(reshape2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=T)

expose_stan_functions("hiv_project/simpleepp/stan_files/chunks/ART_DIAG_model.STAN")

rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

step_vector<-seq(1970.1,by = 0.1,length.out = 500)
kappa_params<-c(log(0.5),log(0.1),0.3,1995)

kappa<- exp(rlogistic(step_vector,kappa_params))                       ## This is our transmission parameter for the output

iota_val<-0.0001                                                       ## This is our initial proportion infected

mu <- 1/35                                                             ## Population level death rate
            
sigma<- c(1/3.16,1/2.13,1/3.2)                                         ## Progression along Cd4 stages among untreated

mu_i <- c(0.003, 0.008, 0.035, 0.27)                                   ## Mortality by stage, no ART

mu_d <- c(0.003,0.008,0.035,0.27)                                      ## Mortality be stage on diagnosed

mu_a <- c(0.002, 0.006, 0.006, 0.03)                                   ## Mortality by stage, on ART

omega <- 0.9                                                          ## Reduction in transmissability on art

theta <- 0.2                                                           ## Reduction when know diagnosed

dt <- 0.1

start<- 1970

diag_start<- 1982

art_start<-1996

sigmoid_curve_function<-function(x,lp){                               ## This is our function to produce the change in 
  (1 / (1 + exp(-lp[1]*x))) * lp[2]                                   ## diag rates and art uptake rates through time 
}

zero_time_art<-rep(0,((art_start-start)*10))                          ## How long initially art is not available       

zero_time_diag<-rep(0,((diag_start-start)*10))                        ## How long initially diagnoses not available

art_length<-length(kappa)+1 - length(zero_time_art)                   ## how long art is available for 

diag_time<-length(kappa)+1 - length(zero_time_diag)                   ## how long diag is available for

elll<-c(0.05,0.7)                                                     ## first number represents the curve on the uptake,
ellm<-c(0.05,0.7)                                                     ## second number represents the peak fraction on ART per time step
elln<-c(0.05,0.8)
ello<-c(0.1,0.9)

art_col_1<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),l = elll)) ## These are our columns for the Alpha matirx
art_col_2<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),ellm))     ## 1 is cd4>500, 4 is cd4<200
art_col_3<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),elln))
art_col_4<-c(zero_time_art, sigmoid_curve_function(seq(1,art_length),ello))

diag_params<-c(0.1,0.9)
diag_params_aids<-c(0.5,0.9)

diag_col_1<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))
diag_col_2<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))   ## Proportion to be diagnosed at each infected cd4 stage
diag_col_3<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params))
diag_col_4<-c(zero_time_diag, sigmoid_curve_function(seq(1,diag_time),diag_params_aids))

alpha<-cbind(art_col_1,art_col_2,art_col_3,art_col_4)
diag<-cbind(diag_col_1,diag_col_2,diag_col_3,diag_col_4)

art_prog<-c(1/10,1/10,1/10)

obem<-(60*24*365)/(60*10^6)*5

test_run<-simpleepp_art_diag(kappa = kappa,iota = iota,alpha = alpha,mu = mu,sigma = sigma, mu_i = mu_i, mu_d = mu_d,mu_a = mu_a,
                             omega = omega,theta = theta,dt = dt,start = start,diag_start = diag_start,art_start = art_start,
                             diag = diag,art_prog = art_prog,onem = obem)
sim_data_art<-data.frame(test_run)
sim_data_art$time<-seq(start,by = dt,length.out = nrow(alpha))
names(sim_data_art)[1:7]<-c("kappa","art","diag","N","ART_inc", "incidence","prevalence")
sim_data_art$prev_percent<-sim_data_art$prevalence * 100

sim_plot<-function(sim_df,col_ids_to_get_rid){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=lambda),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  melt_df<-sim_df[,-c(col_ids_to_get_rid)]
  melted_df_funky<-melt(melt_df,id="time")
  
  three_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot,incidence=incidence_plot,kappa=kappa_plot,whole=three_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_data_art,c(2:5,9))
plot(plotted_sim$whole)

sample_range<-1970:2015
sample_years<-46
sample_n<-1000


sample_function<-function(year_range,number_of_years_to_sample,people_t0_sample,simulated_df,prevalence_column_id){
  sample_years_hiv <- number_of_years_to_sample # number of days sampled throughout the epidemic
  sample_n <- people_t0_sample # number of host individuals sampled per day
  
  # Choose which days the samples were taken. 
  # Ideally this would be daily, but we all know that is difficult.
  sample_time_hiv = sort(sample(year_range, sample_years_hiv, replace=F))
  
  # Extract the "true" fraction of the population that is infected on each of the sampled days:
  sample_propinf_hiv = simulated_df[simulated_df$time %in% sample_time_hiv, prevalence_column_id]
  
  ## this just samples our prevalence, to get a probability that the sample we take is HIV infected then we need to divide
  ## by 100
  
  #sample_propinf_hiv<-sample_propinf_hiv/100
  
  # Generate binomially distributed data.
  # So, on each day we sample a given number of people (sample_n), and measure how many are infected.
  # We expect binomially distributed error in this estimate, hence the random number generation.
  sample_y_hiv_prev = rbinom(sample_years_hiv, sample_n, sample_propinf_hiv)
  sample_prev_hiv_percentage<-(sample_y_hiv_prev/sample_n)*100
  
  ## lets have a ggplot of the y (infected) and out sample of Y over time 
  sample_df_100<-data.frame(cbind(sample_time_hiv,sample_prev_hiv_percentage,sample_y_hiv_prev,sample_time_hiv,sample_propinf_hiv))
  return(sample_df_100)  
}



sample_df_1000_second<-sample_function(sample_range,sample_years,sample_n,
                                       simulated_df = sim_data_art,prevalence_column_id = 7)


plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = sim_data_art,sample_df = sample_df_1000_second)

diagnosed_cd4_sample<-function(year_range_to_sample_from,integer_number_of_years_we_sampled_from,simulated_data,
                              col_id_time_and_inc){
  data_output<-NULL
  
  for(i in year_range_to_sample_from[1]:year_range_to_sample_from[integer_number_of_years_we_sampled_from]){
    overall_year<- i
    tot_for_indep_year<-0
    for(j in 1:10){
      add_on<-(j-1)/10
      overall_year_dec<-overall_year + add_on
      art_inc<-simulated_data[simulated_data$time %in% overall_year_dec, col_id_time_and_inc$inc]
      tot_for_indep_year<-tot_for_indep_year + art_inc
    }
    
    year_data<-cbind(overall_year,tot_for_indep_year)
    data_output<-rbind.data.frame(data_output,year_data)
  }
  
  poisson_samples<-rpois(nrow(data_output),data_output$tot_for_indep_year)
  
  data_output$samples<-poisson_samples
  names(data_output)<-c("year","true","sampled")
  
  sample_plot<-ggplot(data = data_output)+geom_line(aes(x=year,y=true),colour="midnightblue")+
    geom_point(aes(x=year,y=sampled),colour="red")+labs(x="year",y="Number of ART starting with CD4 > 500",
                                                        title="Sampled numbers of new ART starters with CD4 > 500")
  
  return(list(sample_data=data_output,sample_plot=sample_plot))
  
}

a<- art_start:2015
b<-length(a)
col_ids<-list(time=8,inc=5)

poisson_sampled_data<-diagnosed_cd4_sample(a,b,simulated_data = sim_data_art,col_id_time_and_inc = col_ids)

poisson_sampled_data$sample_plot

















