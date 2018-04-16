
  functions{

  //Coding up the simple four bin cd4 count model with intially a fixed R
  // I will input the progression and death rates, leaving it to estimate only eta and kappa 
  
  real[] hiv_SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  real dydt[9];
  
  dydt[1] = (params[1] * 5000) - (params[2] * y[6] * y[1] / y[7] ) - ( 0.02857 * y[1] ) ;
  
  dydt[2] = (params[2] * y[6] * y[1] / y[7] ) - (0.003 * y[2] ) - (1 / 3.16 * y[2]);
  
  dydt[3] = (1 / 3.16 * y[2]) - (0.008 * y[3]) - (1 / 2.13 * y[3]);
  
  dydt[4] = (1 / 2.13 * y[3]) - (0.035 * y[4]) - (1 / 3.2 * y[4]);
  
  dydt[5] = (1 / 3.2 * y[4]) - (0.27 * y[5]);
  
  dydt[6] = dydt[2] + dydt[3] + dydt[4] + dydt[5];
  
  dydt[7] = dydt[1] + dydt[6];
  
  dydt[8] = (dydt[6] + y[6]) / (dydt[7] + y[7]) - y[8];
  
  dydt[9] = ((params[2] * y[6] * y[1] / y[7] ) / y[7]) - y[9];
  
  return dydt;
  
  
  
  }
  }

data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_sample; // Number of hosts sampled at each time point.
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  
  int y[n_obs]; // The binomially distributed data
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[0];
  int x_i[0];
}


parameters {
  real<lower = 0, upper = 1> params[n_params];
  real<lower = 0> I0; // Initial fraction of hosts susceptible
  }
  
transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real y0[n_difeq]; // Initial conditions for the system 
  
y0[1] = 1000000 - I0;      // initial susceptible population
y0[2] = I0;              // initial number with >500 CD4
y0[3] = 0;               // initial number with 500>cd4>350
y0[4] = 0;              // initial number with cd4 count >200 <350      
y0[5] = 0;               // initial number with dc4 count <200
y0[6] = I0;             // This is the initial number of infected individuals
y0[7] = 1000000;       // This is the total population size  
y0[8] = I0 / 1000000;  // This is the starting prevalence in the system
y0[9] = I0 / 1000000;  //This is the starting incidence in the system

y_hat = integrate_ode_rk45(hiv_SI, y0, t0, ts, params, x_r, x_i);


}

model {
  params ~ uniform(0, 1); //constrained at the lower limit of all the paramaters and the upper limit of the et paramter
  I0 ~ normal(0, 100 ); //constrained to be positive.
  
  y ~ binomial(n_sample, y_hat[, 8]  ); //y_hat[,8] is the prevalence as a fraction of the total population
  
  }
generated quantities {
  // Generate predicted data over the whole time series:
real fake_I[n_fake, n_difeq];



fake_I = integrate_ode_rk45(hiv_SI, y0, t0, fake_ts, params, x_r, x_i);

}


