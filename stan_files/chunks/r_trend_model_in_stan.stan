
functions{

matrix simpleepp_no_art(vector kappa, real iota, real mu, vector sigma,   // So in this case we have @kappa as our R parameter, @iota as our initial fraction infected
		 vector mu_i, real dt){                                           // @mu is our population death rate, @sigma is our progression rate through CD4 stages
                                                                          // @mu_i is our death rate for each of the cd4 stages, the length of this also determines the number of stages
                                                                          // and finally @dt is our time step, in Euler this is one tenth of the time step 																		  
  vector[rows(kappa)+1] S;                                                 
  matrix[rows(mu_i), rows(kappa)+1] I;
  vector[rows(kappa)+1] rho;                                              // rho is our I/N so when multipled by kappa it becomes our force of infection 
  vector[rows(kappa)+1] lambda;                                           // This is our force of infection 
  int DS;
  DS = rows(mu_i);

  // initial values
  S[1] = 1000000 * (1 - iota);
  I[1, 1] = 1000000 * iota;
  for(m in 2:DS)
    I[m, 1] = 0;
    
  rho[1] = 1 - S[1] / (S[1] + sum(I[,1]));
  lambda[1] = 0;

  for(t in 1:rows(kappa)) {
    real artcov;
    real deaths;
    real infections;
    real It;
    

    It = sum(I[,t]);
    
    
    lambda[t+1] = kappa[t] * rho[t];

    deaths = mu * (S[t] + It) + sum(mu_i .* I[,t]);
    S[t+1] = S[t] + dt*(-lambda[t+1] * S[t] - mu * S[t] + deaths);                            // 

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS]) * I[DS, t]);
    
        
    rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]) );
  }
  
  
  return(append_col(append_row(kappa[1], kappa),
		    append_col(lambda, rho)));
}

}

data {
 
 int<lower = 1> n_obs;                                                  // The number of time points we observe
 int<lower = 1> n_sample;                                               // The number of people we observe from the population
 int<lower = 0, upper = n_sample> y[n_obs];                             // This is the number of poeple infected from our random draw 
 int time_steps_euler;                                                  // This is our number of points over which to evaulate the function 
 int time_steps_year;                                                   // This is our number of years we evaluate
 int sample_step_euler;                                                 // Length of our prevalence data for Euler trend  
 vector[sample_step_euler] xout;                                        // This is our sequence of time points
 int estimate_years;                                                    // This is the number of years we predict data for
 real mu;                                                               // This is our population level death rate
 vector[3] sigma;                                                       // This is our vector of transition rates
 vector[4] mu_i;                                                        // This is our vector of death rates
 real dt;                                                               // This is our time step
 real dt_2;                                                             // this is our second time step
 int rows_to_interpret[time_steps_year - estimate_years ];              // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years
 vector[sample_step_euler] rho;                                         // This is our expanded y vector
 
 int beta_points;
}
 
 
 
parameters{
 
 real<lower = 0, upper=1> iota;                                         // The proportion of the population initially infected 
 vector<lower = 0>[beta_points] beta;                                   // This is the knot point values
 real<lower = 0, upper = 30> time_param;                            // This is the years it takes to stabilise the epidemic
 real<lower = 0, upper = 5> r_init;                                     // This is the R0 of the epidemic 
 
 
 
}

 
transformed parameters{
  
  vector[sample_step_euler - 1] r;                                          // Initializes the vector to store our r values through time
  matrix[size(rows_to_interpret), 3] y_hat;                             // This creates the matrix in which to store our transformed y_hat values
  real time_func;                                                       // Creates the values for which our time difference can take	
  real a;																// Creates the variable to store the log difference in r between time points 
  
  r[1] = r_init;                                                        // Links our r vector with the r0 we are letting the system estimate.
  
  for(i in 2:(time_steps_euler - 1)){
   time_func = xout[i-1] - 1970 - time_param;
    if(time_func <= 0)
	
	time_func = 0;
	
	else 
	
	time_func = time_func;
    
    a = (beta[2] * (beta[1] - r[i-1])) + (beta[3] * rho[i-1]) + (beta[4] * (rho[i] - rho[i-1]) * time_func / rho[i-1]);
	
    r[i] = exp(a) * r[i-1];
	
	}
  
  y_hat = simpleepp_no_art( r , iota, mu, sigma, mu_i, dt_2)[rows_to_interpret, ];

 
 }

 
model {
 
 iota ~ normal(0, 0.5);                                                  // This models our initial population size 
 
 beta ~ normal(0, 0.1);                                                  // So this models our beta values as a normal distribution, these can be < 0
 
 time_param ~ normal(15, 5);
 
 r_init ~ normal(0,0.5);

 y ~ binomial(n_sample, y_hat[ , 3]);                                   // This fits it to the binomially sampled data 

} 

generated quantities{
   
  vector[sample_step_euler -1 ] r_fit; 
  matrix[sample_step_euler, 3] fitted_model;
  real a_fit;
  real time_func_fit;
  
  
  r_fit[1] = r_init;
  
  for(i in 2:(time_steps_euler - 1)){
   time_func_fit = xout[i-1] - 1970 - time_param;
    if(time_func <= 0)
	
	time_func_fit = 0;
	
	else 
	
	time_func_fit = time_func_fit;
    
    a_fit = (beta[2] * (beta[1] - r[i-1])) + (beta[3] * rho[i-1]) + (beta[4] * (rho[i] - rho[i-1]) * time_func_fit / rho[i-1]);
	
    r_fit[i] = exp(a) * r[i-1];
	
	}
  
  fitted_model = simpleepp_no_art( r_fit , iota, mu, sigma, mu_i, dt_2);

} 










