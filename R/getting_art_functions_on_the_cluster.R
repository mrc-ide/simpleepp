##################################################################################################################################
## ART functions on the cluster ##################################################################################################
##################################################################################################################################

setwd("Y:")
options(didehpc.username = "jd2117",didehpc.home = "Y:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts"

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),
                            sources = "simpleepp/R/ART_model_simple_functions_cluster.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)


###############################################################################################################################
## NOw we've got the queue set up we need to generate the stan_data for using in the cluster runs #############################
###############################################################################################################################
require(rstan)
require(ggplot2)
require(reshape2)
require(ggpubr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=T)

expose_stan_functions("simpleepp/stan_files/chunks/ART_diag_prev_fitting.stan")

sim_art_epidemic<-function(kappa_params,params){

  rlogistic <- function(t, p) {
   p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
   }

step_vector<-seq(1970.1,by = 0.1,length.out = 500)
kappa_params<-c(log(kappa_params$peak),log(kappa_params$low),kappa_params$grad,kappa_params$decline)

kappa<- exp(rlogistic(step_vector,kappa_params))                       ## This is our transmission parameter for the output

iota_val<-params$iota                                                       ## This is our initial proportion infected

mu <- params$mu                                                             ## Population level death rate

sigma<- params$sigma                                         ## Progression along Cd4 stages among untreated

mu_i <- params$mu_i                                   ## Mortality by stage, no ART

mu_d <- params$mu_d                                      ## Mortality be stage on diagnosed

mu_a <- params$mu_a                                   ## Mortality by stage, on ART

omega <- params$omega                                                          ## Reduction in transmissability on art

theta <- params$theta                                                           ## Reduction when know diagnosed

dt <- params$dt

start<- params$start

diag_start<- params$diag_start

art_start<-params$art_start

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

art_prog<-params$art_prog

obem<-params$obem

test_run<-simpleepp_art_diag(kappa = kappa,iota = iota,alpha = alpha,mu = mu,sigma = sigma, mu_i = mu_i, mu_d = mu_d,mu_a = mu_a,
                             omega = omega,theta = theta,dt = dt,start = start,diag_start = diag_start,art_start = art_start,
                             diag = diag,art_prog = art_prog,onem = obem)
sim_data_art_early_72_ART<-data.frame(test_run)
sim_data_art_early_72_ART$time<-seq(start,by = dt,length.out = nrow(alpha))
names(sim_data_art_early_72_ART)[1:7]<-c("kappa","art","diag","N","ART_inc", "incidence","prevalence")
sim_data_art_early_72_ART$prev_percent<-sim_data_art_early_72_ART$prevalence * 100

params$kappa<-kappa_params

return(list(df=sim_data_art_early_72_ART,params_used=params))

}

kappa_params<-list(peak=0.5,low=0.1,grad=0.3,decline=1995)
iota<-0.0001
mu <- 1/35                                                             ## Population level death rate
sigma<- c(1/3.16,1/2.13,1/3.2)                                         ## Progression along Cd4 stages among untreated
mu_i <- c(0.003, 0.008, 0.035, 0.27)                                   ## Mortality by stage, no ART
mu_d <- c(0.003,0.008,0.035,0.27)                                      ## Mortality be stage on diagnosed
mu_a <- c(0.002, 0.006, 0.006, 0.03)                                   ## Mortality by stage, on ART
omega <- 0.8                                                          ## Reduction in transmissability on art
theta <- 0.2                                                           ## Reduction when know diagnosed
dt <- 0.1                                                              ## Time step difference
start<- 1970                                                           ## start of the epidemic
diag_start<- 1971                                                      ## start of diagnosis
art_start<-1972                                                        ## start of ART
art_prog<-c(1/10,1/10,1/10)                                            ## Progression rates through CD4 stages on ART
obem<-(60*24*365)/(60*10^6)*5                                          ## Birth rate in the population

params<-list(iota=iota,mu=mu,sigma=sigma,mu_i=mu_i,mu_d=mu_d,mu_a=mu_a,
             omega=omega,theta=theta,dt=dt,start=start,diag_start=diag_start,
             art_start=art_start,art_prog=art_prog,obem=obem)

early_art_epidemic<-sim_art_epidemic(kappa_params = kappa_params,params = params)


sim_plot<-function(sim_df,col_ids_to_get_rid){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=incidence),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  melt_df<-sim_df[,-c(col_ids_to_get_rid)]
  melted_df_funky<-melt(melt_df,id="time")
  
  three_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot,incidence=incidence_plot,kappa=kappa_plot,whole=three_variable_plot))
  
  
}

plotted_sim<-sim_plot(early_art_epidemic$df,c(2:5,9))
plot(plotted_sim$whole)
plot(plotted_sim$incidence)
plot(plotted_sim$prevalence)

###############################################################################################################################
## Now we will perform our sampling from the population #######################################################################
###############################################################################################################################


sample_range<-1970:2015
sample_years<-46
sample_n<-10000


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



sample_df_prev<-sample_function(sample_range,sample_years,sample_n,
                                simulated_df = early_art_epidemic$df,prevalence_column_id = 7)


plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = early_art_epidemic$df,sample_df = sample_df_prev)

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

year_seq<- art_start:2015
number_year_obs<-length(year_seq)
col_ids<-list(time=8,inc=5)

poisson_sampled_data<-diagnosed_cd4_sample(year_seq,number_year_obs,simulated_data = early_art_epidemic$df,
                                           col_id_time_and_inc = col_ids)

poisson_sampled_data$sample_plot

sneaky_sample<-(art_start-start):(2015-start)*10+1
sneaky_poiss<-rpois(length(sneaky_sample),sim_data_art[sneaky_sample,5])

###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################
penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-0:45*10+1

stan_data_discrete_prev_rw<-list(
  n_obs = sample_years,
  y = as.array(sample_df_prev$sample_y_hiv_prev),
  n_sample = sample_n,
  time_steps_euler = length(xout)+1,
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

RW_art_prev_72_ART<-obj$enqueue(prev_data_fitting_RW(simulated_data = early_art_epidemic$df,
                                                            stan_data = stan_data_discrete_prev_rw,
                                                     params_used_for_sim = early_art_epidemic$params_used))
RW_art_prev_72_ART$status()
id_rw_art_72_prev<-RW_art_prev_72_ART$id
save(id_rw_art_72_prev,
     file = 
       "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/RW_ART_72_prev_fitting_15_34_20_MAY")


##############################################################################################################################
## Now lets run our RW poissson fit ##########################################################################################
##############################################################################################################################

penalty_order<-1
xout<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=penalty_order)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-seq((art_start-start)*10+1,(year_seq[length(year_seq)]-start+1)*10)

stan_data_discrete_count_rw<-list(
  n_obs = number_year_obs,
  y = as.array(poisson_sampled_data$sample_data$sampled),
  time_steps_euler = length(xout)+1,
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

RW_ART_count_72_art<-obj$enqueue(count_data_fitting_RW(early_art_epidemic$df,
                                                       stan_data_discrete_count_rw,
                                                       early_art_epidemic$params_used))
RW_ART_count_72_art$status()
id_rw_art_72_count<-RW_ART_count_72_art$id
save(id_rw_art_72_count,
     file = "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/RW_count_ART_72_15_35_MAY_20")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/RW_ART_COUNT_FITTING")

RW_art_count_check<-obj$task_get(id_rw_art_count)
a<-obj$task_times()
a$task_id[69]
RW_art_count_check<-obj$task_get(a$task_id[69])
RW_count_results<-RW_art_count_check$result()
plot(RW_count_results$df_output$median)
plot(RW_count_results$incidence_df$median)
plot(RW_count_results$r_fit_df$median)

##############################################################################################################################
## NOw for the splines #######################################################################################################
##############################################################################################################################

xout<-seq(1970,2020,0.1)
splines_creator<-function(knot_number,penalty_order){
  
  nk <- knot_number # number of knots
  dk <- diff(range(xout))/(nk-3)
  knots <- xstart + -3:nk*dk
  spline_matrix<- splineDesign(knots, step_vector, ord=4)
  penalty_matrix <- diff(diag(nk), diff=penalty_order)
  
  return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
  
}

knot_number= 7
penalty_order= 2
splines_matrices<-splines_creator(knot_number,penalty_order)
rows_to_evaluate<-0:45*10+1

stan_data_discrete_prev_spline<-list(
  n_obs = sample_years,
  y = as.array(sample_df_prev$sample_y_hiv_prev),
  n_sample = sample_n,
  knot_number = knot_number,
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

spline_art_72_prev<-obj$enqueue(prev_data_fitting_spline(early_art_epidemic$df,stan_data_discrete_prev_spline,
                                                         early_art_epidemic$params_used))
spline_art_72_prev$status()
spline_art_prev_id_72<-spline_art_72_prev$id
save(spline_art_prev_id_72,file =
       "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/spline_art_72_prev_15_37_MAY_20")


### NOw for the count data spline


rows_to_evaluate<-seq((art_start-start)*10+1,(year_seq[length(year_seq)]-start+1)*10)


stan_data_discrete_count_spline<-list(
  n_obs = nrow(poisson_sampled_data$sample_data),
  y = as.array(poisson_sampled_data$sample_data$sampled),
  n_sample = sample_n,
  knot_number = knot_number,
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  alpha = alpha,
  mu_d = mu_d,
  mu_a = mu_a,
  omega = omega,
  theta = theta,
  start = start,
  diag_start = diag_start,
  art_start = art_start,
  diag = diag,
  art_prog = art_prog,
  onem = obem
  
)

spline_art_count_72<-obj$enqueue(count_data_fitting_spline(early_art_epidemic$df,stan_data_discrete_count_spline,
                                                           early_art_epidemic$params_used))
spline_art_count_72$status()
spline_art_72_count_ID<-spline_art_count_72$id
save(spline_art_72_count_ID, file = 
       "C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/spline_art_COUNT_72_15_38_MAY_20")

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/SPLINE_COUNT_FIT_ART")

count_spline_check<-obj$task_get(spline_count_id)
count_spline_check$result()

###############################################################################################################################
## Loading up the results for each of the runs, first the normal 1996, 1982 art ad diag class #################################
###############################################################################################################################

a<-load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/RW_ART_COUNT_FITTING_THIRD_RUN")
b<-load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/RW_ART_PREVALENCE_FITTING_THIRD_RUN")
c<-load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/SPLINE_COUNT_FIT_ART_third_RUN")
d<-load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/SPLINE_PREV_FIT_ART_THIRD_RUN")
a
b
c
d

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/RW_ART_COUNT_FITTING_THIRD_RUN")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/RW_ART_PREVALENCE_FITTING_THIRD_RUN")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/SPLINE_COUNT_FIT_ART_third_RUN")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/SPLINE_PREV_FIT_ART_THIRD_RUN")

RW_ART_COUNT_96_ART<-obj$task_get(id_rw_art_count)
RW_ART_PREV_96_ART<-obj$task_get(id_rw_art_prev)
SPLINE_ART_COUNT_96_ART<-obj$task_get(spline_count_id)
SPLINE_ART_PREV_96_ART<-obj$task_get(spline_art_prev_id)

rw_art_96_count_results<-RW_ART_COUNT_96_ART$result()
rw_art_96_prev_results<-RW_ART_PREV_96_ART$result()
spline_96_count_results<-SPLINE_ART_COUNT_96_ART$result()
spline_96_prev_results<-SPLINE_ART_PREV_96_ART$result()

########
# NOw for the earlier ART results
########

load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/RW_ART_COUNT_FITTING")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/RW_ART_PREVALENCE_FITTING")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/SPLINE_COUNT_FIT_ART")
load("C:/Users/josh/Dropbox/hiv_project/analysis_of_cluster_run_datasets/art_simpleepp/cluster_run_ids/early_72_ART_ids/SPLINE_PREV_FIT_ART")

RW_ART_COUNT_72_ART<-obj$task_get(id_rw_art_count_early_art_72)
RW_ART_PREV_72_ART<-obj$task_get(id_rw_art_prev_early_72)
SPLINE_ART_COUNT_72_ART<-obj$task_get(spline_count_id_early_art_72)
SPLINE_ART_PREV_72_ART<-obj$task_get(spline_art_prev_id_early_art_72)

rw_art_72_count_results<-RW_ART_COUNT_72_ART$result()
rw_art_72_prev_results<-RW_ART_PREV_72_ART$result()
spline_72_count_results<-SPLINE_ART_COUNT_72_ART$result()
spline_72_prev_results<-SPLINE_ART_PREV_72_ART$result()





