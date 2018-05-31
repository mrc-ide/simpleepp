
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
 vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
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
 matrix[time_steps_euler - 1 , knot_number + 2 ] X_design;              // This is our spline design matirx that we are modelling kappa with. 
 matrix[knot_number + 2 - penalty_order  , knot_number +2 ] D_penalty;                      // This is our penalty matrix, can be first or second order depending on the R code 
 real mu;                                                               // This is our population level death rate
 vector[3] sigma;                                                       // This is our vector of transition rates
 vector[4] mu_i;                                                        // This is our vector of death rates
 real dt;                                                               // This is our time step
 real dt_2;                                                             // this is our second time step
 int rows_to_interpret[time_steps_year - estimate_years ];              // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years
 vector[knot_number] knots;                                             // the sequence of knots
 int spline_degree;                                                     // the degree of spline (is equal to order - 1)
 real X[time_steps_euler];                                              // Our vector of time points over which to evaluate the B spline
}
 
 transformed data {
  int num_basis = knot_number + spline_degree - 1;                        // total number of B-splines
  matrix[time_steps_euler, num_basis] B;                                  // matrix of B-splines
  vector[spline_degree + knot_number] ext_knots_temp;
  vector[2*spline_degree + knot_number] ext_knots;                        // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[knot_number], spline_degree));
  for (ind in 1:num_basis)
    B[:,ind] = (build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[time_steps_euler, knot_number + spline_degree - 1] = 1;
}
 
parameters{
 
 real<lower = 0, upper=1> iota;                                         // The proportion of the population initially infected 
 real<lower = 0> sigma_pen;                                             // This is the penalty to apply to the spline to make it smooth
 vector[num_basis] beta;
  
}

 
transformed parameters{

  matrix[size(rows_to_interpret), 3] y_hat;
  vector[time_steps_euler] splino;
  

  splino = B * beta;
  
    
  y_hat = simpleepp_no_art( splino , iota, mu, sigma, mu_i, dt_2)[rows_to_interpret, ];

 
 }

 
model {
  
    
  target += normal_lpdf( D_penalty * beta  | 0, sigma_pen);
 
  iota ~ normal(0, 0.5);                                                    // This models our initial population size 
 

  y ~ binomial(n_sample, y_hat[ , 3]);                                     // This fits it to the binomially sampled data 

} 

generated quantities{

matrix[time_steps_euler +1 , 3] fitted_output;
vector[time_steps_euler] splino_fit;

splino_fit = to_vector(B * beta);
 

fitted_output = simpleepp_no_art(splino_fit, iota, mu, sigma, mu_i, dt_2);       // This models the data after our fitting procedure

} 












