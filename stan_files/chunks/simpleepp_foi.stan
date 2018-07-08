/**
 * This is a simplfied version of the EPP model. It removes all of the demographic
 * details and assumes a constant population size and simplifies the HIV natural
 * history and representation of antiretorviral treatment. A simple model seems a
 * good starting place to capture the essence of the problem.
 * @param kappa A vector of transmission rate at each time step.
 * @param iota Proportion initally infected.
 * @param alpha A vector of the proportion ART coverage at end of each time step.
 * @param mu Natural (non-HIV) rate of mortality / exiting population.
 * @param sigma Rate of progression through each infection stage (vector length m-1).
 * @param mu_i HIV mortality rate by infection stage (vector length m).
 * @param mu_a HIV mortality rate for persons on treatment.
 * @param omega Reduction in transmission rate for persons on ART
 * @param dt Time step for Euler integration.
 * @return An array with columns transmission rate, incidence, and prevalence.
 */
functions{
matrix simpleepp_foi(vector f_t, real iota , real mu, vector sigma,
		 vector mu_i, real dt, int foi_flag,real kappa_init){

  vector[rows(f_t)+1] S;
  matrix[rows(mu_i), rows(f_t)+1] I;
  vector[rows(f_t)+1] rho;
  vector[rows(f_t)+1] kappa;
  vector[rows(f_t)+1] lambda;
  matrix[rows(f_t)+1, 2] return_mat;
  int DS;
  DS = rows(mu_i);

  // initial values

  if(foi_flag == 0) {
    kappa[1] = f_t[1];
    S[1] = 1000 * (1 - iota);
    I[1, 1] = 1000 * iota;
  } else if(foi_flag == 1) {
    kappa[1] = kappa_init;
    S[1] = 1000 - (f_t[1] * (1000 / kappa[1]));
    I[1, 1] = f_t[1] * (1000 / kappa[1]);
  }

  for(m in 2:DS)
    I[m, 1] = 0;

  rho[1] = 1 - S[1] / (S[1] + sum(I[,1]) );
  lambda[1] = kappa[1] * rho[1];

  for(t in 1:rows(f_t)) {
    real deaths;
    real infections;
    real It;


    It = sum(I[,t]);

    if(foi_flag == 0){
      kappa[t+1] = f_t[t];
      lambda[t+1] = kappa[t+1] * rho[t];
    } else if(foi_flag == 1) {
      lambda[t+1] = f_t[t];
      kappa[t+1] = lambda[t+1] / (rho[t]);
    }

    deaths = mu * (S[t] + It) + sum(mu_i .* I[,t]);
    S[t+1] = S[t] + dt*(-lambda[t+1] * S[t] - mu * S[t] + deaths);

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS]) * I[DS, t]);

    rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]));
	}


  return_mat = append_col(kappa, lambda);

  return(append_col(return_mat, rho));
}

}

data {

 int<lower = 1> n_obs;                                                  // The number of time points we observe
 int<lower = 1> n_sample;                                               // The number of people we observe from the population
 int<lower = 0, upper = n_sample> y[n_obs];                             // This is the number of poeple infected from our random draw
 int time_steps_euler;                                                  // This is our number of points over which to evaulate the function
 int time_steps_year;                                                   // This is our number of years we evaluate
 int penalty_order;                                                     // This is the degree of penalty we place on the
 int knot_number;                                                       // Number of knots, our length of beta
 int estimate_years;                                                    // This is the number of years we predict data for
 matrix[time_steps_euler - 1 , knot_number  ] X_design;                 // This is our spline design matirx that we are modelling kappa with.
 matrix[knot_number - penalty_order , knot_number ] D_penalty;          // This is our penalty matrix, can be first or second order depending on the R code
 real mu;                                                               // This is our population level death rate
 vector[3] sigma;                                                       // This is our vector of transition rates
 vector[4] mu_i;                                                        // This is our vector of death rates
 real dt;                                                               // This is our time step
 real dt_2;                                                             // this is our second time step
 int rows_to_interpret[time_steps_year - estimate_years ];              // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years
 int foi_flag;
 real kappa_init;

}



parameters{

 real<lower = 0, upper=1> iota;                                         // The proportion of the population initially infected
 vector[knot_number] beta;                                   // This is the knot point values
 real<lower = 0, upper =1> sigma_pen;                                             // This is the penalty to apply to the spline to make it smooth


}


transformed parameters{

  vector[rows(X_design)] kappa;
  matrix[size(rows_to_interpret), 3] y_hat;

  kappa = exp(X_design * beta);
  y_hat = simpleepp_foi( kappa , iota, mu, sigma, mu_i, dt_2, foi_flag, kappa_init)[rows_to_interpret, ];


 }


model {



 iota ~ normal(0, 0.5);                                                  // This models our initial population size

 target += normal_lpdf( D_penalty * beta | 0, sigma_pen);                // So this models our penalized spline with a slightly altered distribution command


  y ~ binomial(n_sample, y_hat[ , 3]);                                   // This fits it to the binomially sampled data

}

generated quantities{

vector[rows(X_design)] kappa_fit;
matrix[time_steps_euler , 3] fitted_output;


kappa_fit = exp(X_design * beta);
fitted_output = simpleepp_foi( kappa_fit , iota, mu, sigma, mu_i, dt_2, foi_flag, kappa_init);       // This models the data after our fitting procedure

}
