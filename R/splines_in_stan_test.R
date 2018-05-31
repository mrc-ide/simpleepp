###################################################################################################################################################
## Working with the splines being created in stan itself ##########################################################################################
###################################################################################################################################################

require(rstan)
require(splines)
require(ggplot2)
require(reshape2)
rstan_options(auto_write=T)
options(mc.cores =  parallel::detectCores())

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/splines_in_stan_simpleepp.stan")



run_simulated_model<-function(params,times){
  
  mu <- params$mu                               # Non HIV mortality / exit from population
  sigma <- params$sigma           # Progression from stages of infection
  mu_i <- params$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params$iota
  
  
  dt <- times$dt                                     # time step
  nsteps <- as.integer(times$years/dt)                   # number of steps
  xstart <- times$start                                # the start of the epidemic
  step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  
  
  rlogistic <- function(t, p) {
    p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4]))) 
  }
  
  bell_curve_func<-function(a,b,c,x){
    f_t<-a*exp(-((x-b)^2)/(2*c^2))
    
  }
  
  kappa_params<-c(log(params$kappa[1]),log(params$kappa[2]),params$kappa[3],params$kappa[4])
  kappa<-exp(rlogistic(step_vector,kappa_params))
  
  #kappa<-bell_curve_func(0.5,0.2,0.2,seq(0,0.998,0.002))
  
  mod<-simpleepp_no_art(kappa = kappa,iota,mu,sigma,mu_i,dt)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa))
  
  
}

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20)           # Progression from stages of infection
mu_i <- c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
iota<-0.0001


dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_output_changed_to_bell_curve<-run_simulated_model(params_sim,times_sim)

sim_plot<-function(sim_df){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=lambda),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  melt_df<-sim_df[,-c(4)]
  melted_df_funky<-melt(melt_df,id="time")
  
  three_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot,incidence=incidence_plot,kappa=kappa_plot,whole=three_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output_changed_to_bell_curve$sim_df)
plot(plotted_sim$whole)

###################################################################################################################################
## Now we have simulated through our model we can extract some samples form it ####################################################
###################################################################################################################################

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
                                       simulated_df = sim_model_output$sim_df,prevalence_column_id = 3)




plot_sample<-function(sample_df,simulated_df){
  a<-ggplot(data = simulated_df,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
    geom_point(data=sample_df, aes(x=sample_time_hiv,y=sample_prev_hiv_percentage),colour="red",size=1)
  
  return(plot(a))
}

plot_sample(simulated_df = sim_model_output$sim_df,sample_df = sample_df_1000_second)


ggplot(data = sample_df_1000_second,aes(x=sample_time_hiv,y=sample_prev_hiv_percentage))+geom_point(colour="red",size=1.5)

###################################################################################################################################
## Now we have our sample from the simulated data we can call the stan script to sample from this data ############################
###################################################################################################################################

splines_creator<-function(knot_number,penalty_order){
  
  nk <- knot_number # number of knots
  dk <- diff(range(xout))/(nk-3)
  knots <- xstart + -3:nk*dk
  spline_matrix<- splineDesign(knots, step_vector, ord=4)
  penalty_matrix <- diff(diag(nk), diff=penalty_order)
  
  return(list(spline_matrix=spline_matrix,penalty_matrix=penalty_matrix))
  
}

knot_number= 10
penalty_order= 1

splines_matrices<-splines_creator(knot_number,penalty_order)
rows_to_evaluate<-0:45*10+1


num_knots <- 8
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
X <- seq(from=1970, to=2020, by=0.1)
num_data <- length(X)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
a0 <- 0.2
a <- rnorm(num_basis, 0, 5)
B_true <- t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
Y_true <- as.vector(a0*X + a%*%B_true)
Y <- Y_true + rnorm(length(X), 0, 0.2)
#splines_matrices$penalty_matrix<-t(splines_matrices$penalty_matrix)


stan_data_splines<-list(
  n_obs = sample_years,
  n_sample = sample_n,
  y = as.array(sample_df_1000_second$sample_y_hiv_prev),
  time_steps_euler = length(xout),
  penalty_order = penalty_order,
  knot_number = num_knots,
  estimate_years = 5,
  time_steps_year = 51,
  X_design = splines_matrices$spline_matrix,
  D_penalty = splines_matrices$penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  num_data = num_data,
  knots = knots,
  spline_degree = spline_degree,
  Y = Y,
  X = X
)

spline_mod<-stan_model("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/splines_in_stan_simpleepp.stan")
test<-sampling(spline_mod,stan_data_splines,iter = 500,chains = 3,control = list(adapt_delta=0.95))
test_extract<-rstan::extract(test)
test_extract$fitted_output
prev_median<-(apply(test_extract$fitted_output[,,1],2,median))
prev_low<-(apply(test_extract$fitted_output[,,1],2,quantile,probs=c(0.025)))
prev_high<-(apply(test_extract$fitted_output[,,1],2,quantile,probs=c(0.975)))
spline_produced_estimate_df<-cbind.data.frame(prev_low,prev_median,prev_high)
names(spline_produced_estimate_df)<-c("low","median","high")
spline_produced_estimate_df$time<-seq(1970,2020.1,0.1)
plot(sim_model_output$sim_df$kappa)

lines(spline_produced_estimate_df$median)
lines(spline_produced_estimate_df$low,col="red")
lines(spline_produced_estimate_df$high,col="red")

ggplot(data = spline_produced_estimate_df)+geom_line(aes(x=time,y=median),colour="dodgerblue",size=1.05)+
  geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="dodgerblue",fill="dodgerblue",alpha=0.15,size=1.02)+
  geom_line(data = sim_model_output$sim_df, aes(x=time,y=kappa),colour="firebrick1",size=1.05)

knots
spline_degree
num_knots
num_data
Y
X

expose_stan_functions("C:/Users/josh/Dropbox/hiv_project/simpleepp/stan_files/chunks/splines_in_stan_simpleepp.stan")


xout<-seq(1970,2020,0.1)
intro<-rep(knots[1],spline_degree)
outro<-rep(knots[length(knots)],spline_degree)
ext_knots<-c(intro,knots,outro)

a_tot<-NULL

for(i in 1:knot_number){

  a<-build_b_spline(xout,ext_knots,i,(spline_degree+1))
a_tot<-cbind(a_tot,a)
  
  }
a_tot[501,10]=1
dim(a_tot)
matplot(a_tot,type = "l")
matplot(spline_matrix,type = "l")
