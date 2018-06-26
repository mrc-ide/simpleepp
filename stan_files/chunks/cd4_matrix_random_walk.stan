
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

 int penalty_order;                                                     // This is our first or second order penalty

 int time_steps_year;                                                   // This is our number of years we evaluate

 int estimate_period;                                                   // Number of years we predict data for

 matrix[time_steps_euler - 1  , time_steps_year  ] X_design;                // This is our spline design matirx that we are modelling kappa with.

 matrix[time_steps_year - penalty_order , time_steps_year ] D_penalty;  // This is our penalty matrix, can be first or second order depending on the R code

 real mu;                                                               // This is our population level death rate

 vector[3] sigma;                                                       // This is our vector of transition rates

 vector[4] mu_i;                                                        // This is our vector of death rates

 real dt;                                                               // This is our time step

 real dt_2;                                                             // this is our second time step for generating the output from the fitted beta parameters

 int rows_to_interpret[n_obs];                             // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years

  }



parameters{

 real<lower = 0> iota;                                         // The proportion of the population initially infected

 vector[cols(X_design)] beta;                                  // This is the knot point values

 real<lower = 0> sigma_pen;                                    // This is the penalty to apply to the spline to make it smooth


 }


transformed parameters{

 vector[size(rows_to_interpret)] y_hat;

 y_hat = simpleepp_no_art( exp(X_design * beta) , iota, mu, sigma, mu_i, dt_2)[rows_to_interpret, 3];
}


model {



 iota ~ normal(0, 0.25);                                                  // This models our initial population size

 target += normal_lpdf( D_penalty * beta | 0, sigma_pen);                // So this models our penalized spline with a slightly altered distribution command

 y ~ binomial(n_sample, y_hat);                                          // This fits it to the binomially sampled data


}

generated quantities{

matrix[time_steps_euler , 3] fitted_output;
vector[rows(X_design)] fitted_kappa;

fitted_kappa = exp(X_design * beta);
fitted_output = simpleepp_no_art(fitted_kappa, iota, mu, sigma, mu_i, dt_2);       // This models the data after our fitting procedure

}
