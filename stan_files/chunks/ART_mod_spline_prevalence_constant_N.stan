functions{
matrix simpleepp_art_diag(vector kappa, real iota, matrix alpha, real mu, vector sigma,     // so kappa is our transmission parameter, iota our starting prop infected,
		 vector mu_i, vector mu_d, vector mu_a, real omega, real theta, real dt, int start, // alpha our progression from diagnosed to ART, mu our pop death rate, 
		 int diag_start, int art_start, matrix diag, vector art_prog, real onem){           // sigma is our progression through cd4 classes, mu_i death from untreated HIV,
                                                                                            // mu_d is death from diagnosed hiv, mu_a death from art hiv, 
  vector[rows(kappa)+1] S;                                                                  // omega our reduction in foi on art, theta our reduction when diagnosed,   
  matrix[rows(mu_i), rows(kappa)+1] I;                                                      // dt is our time step, start is 1970, diag start will be 1981, art_start 1996
  matrix[rows(mu_i), rows(kappa)+1] D;                                                      // diag will be our diagnosis rates CD4, and art_prog progression through ART CD4 
  matrix[rows(mu_i), rows(kappa)+1] A;
  vector[rows(kappa)+1] ART_first;
  vector[rows(kappa)+1] Diag_first;
  matrix[rows(kappa)+1, 2] first_mat;
  matrix[rows(kappa)+1, 2] second_mat;
  matrix[rows(kappa)+1, 4] return_mat;
  vector[rows(kappa)+1] rho;
  vector[rows(kappa)+1] lambda;
  int DS;
  int delay_to_diag_start;
  int delay_to_ART_start;
  vector[rows(kappa)+1] N;
  matrix[rows(kappa)+1, 3] triplay_mat;
  vector[rows(kappa)+1] art_cd4_inc;
	
  DS = rows(mu_i);
  
  
  delay_to_diag_start = (diag_start - start) * 10;
  delay_to_ART_start = (art_start - start) * 10;

  // initial values
  S[1] = 1000000 * (1 - iota);
  I[1, 1] = 1000000 * iota;
  for(m in 2:DS)
    I[m, 1] = 0;
  for(i in 1:delay_to_diag_start){
	for(m in 1:DS)
    D[m, i] = 0;
	}
  for(i in 1:delay_to_ART_start){
    for(m in 1:DS)
    A[m, i] = 0;
	}
  
  rho[1] = 1 - S[1] / (S[1] + sum(I[,1]) + sum(A[,1]) + sum(D[,1]));
  lambda[1] = 0;

  for(t in 1:rows(kappa)) {
    real artcov;
    real deaths;
	real diagnosed;
    real infections;
    real It;
    real At;
	real Dt;
	matrix[rows(kappa)+1, DS] art_rates;
	matrix[rows(kappa)+1, DS] diag_rates;

	

    It = sum(I[,t]);
	Dt = sum(D[,t]);
    At = sum(A[,t]);

    artcov = At / (It + At + Dt);
	diagnosed = Dt / (It + At + Dt);
    lambda[t+1] = kappa[t] * rho[t] * (1 - omega * artcov) * (1 - theta * diagnosed);
    deaths = mu * (S[t] + It + At + Dt) + sum(mu_i .* I[,t]) + sum(mu_d .* D[, t]) + sum(mu_a .* A[, t]);
    
	
	// Now we are modelling the change in the different compartments 
	
	S[t+1] = S[t] + dt*( onem * S[t] -lambda[t+1] * S[t] - mu * S[t] + deaths);

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1] + diag[t, 1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m] + diag[t, m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS] + diag[t, DS]) * I[DS, t]);
	
	if (t >= delay_to_diag_start - 1){
	D[1,t+1] = D[1, t] + dt * (diag[t, 1] * I[1, t] - (mu + mu_d[1] + alpha[t, 1] + sigma[1]) * D[1, t]);
	for(m in 2:(DS-1))
	   D[m, t+1] = D[m, t] + dt * (diag[t, m] * I[m, t] + (sigma[m-1] * D[m-1, t]) - (mu + mu_d[m] + alpha[t, m] + sigma[m]) * D[m, t]);
	D[DS, t+1] = D[DS, t] + dt * (diag[t, DS] * I[DS, t] + (sigma[DS-1] * D[DS-1, t]) - (mu + mu_d[DS] + alpha[t, DS]) * D[DS, t]);
	}
    
	if( t >= delay_to_ART_start){
	A[1, t+1] = A[1, t] + dt * (alpha[t, 1] * D[1, t] - (mu + mu_a[1] + art_prog[1]) * A[1, t]);
	for(m in 2:(DS-1))
	  A[m, t+1] = A[m, t] + dt * (alpha[t, m] * D[m, t] + (art_prog[m-1] * A[m-1, t]) - (mu + mu_a[m] + art_prog[m]) * A[m, t]);
	A[DS, t+1] = A[DS, t] + dt * (alpha[t, DS] * D[DS, t] + (art_prog[DS-1] * A[DS-1, t]) - (mu + mu_a[DS]) * A[DS, t]);
	}
	
	art_cd4_inc[t] = dt * (alpha[t, 1] * D[1, t]);
    rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]) + sum(A[ ,t+1]) + sum(D[,t+1]));
  }
  
  ART_first = to_vector(A[1,]);
  Diag_first = to_vector(D[1,]);
  first_mat = append_col(ART_first,Diag_first);
  second_mat = append_col(lambda,rho);
  
  for(i in 1:rows(kappa)+1)
  N[i] = S[i] + sum(I[,i]) + sum(A[,i]) + sum(D[,i]);
  
  art_cd4_inc[501] = 0;
  
  triplay_mat = append_col(first_mat, N);
  return_mat = append_col(triplay_mat, art_cd4_inc);
  
  return(append_col(append_row(kappa[1], kappa),
		    append_col(return_mat,second_mat)));
}
}

data {
 
 int<lower = 1> n_obs;                                                  // The number of time points we observe
 int<lower = 0> y[n_obs];                                               // This is the number of poeple infected from our random draw
 int<lower = 1> n_sample;                                               // This is our number of people sampled from the population. 
 int time_steps_euler;                                                  // This is our number of points over which to evaulate the function 
 int penalty_order;                                                     // This is our first or second order penalty 
 int knot_number;                                                       // This is the number of knots in our spline
 int time_steps_year;                                                   // This is our number of years we evaluate  
 matrix[time_steps_euler - 1 , knot_number  ] X_design;                 // This is our spline design matirx that we are modelling kappa with. 
 matrix[knot_number - penalty_order , knot_number ] D_penalty;          // This is our penalty matrix, can be first or second order depending on the R code 
 real mu;                                                               // This is our population level death rate
 vector[3] sigma;                                                       // This is our vector of transition rates
 vector[4] mu_i;                                                        // This is our vector of death rates
 real dt_2;                                                             // this is our second time step for generating the output from the fitted beta parameters
 int rows_to_interpret[n_obs];                                          // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years
 matrix[time_steps_euler , rows(mu_i)] alpha;                           // This is our matrix of ART uptake rates from the diagnosed class
 vector[4] mu_d;                                                        // This is our vecotr of death rates among diagnosed individuals
 vector[4] mu_a;                                                        // This is our vector of death rates among the ART class
 real omega;															// Our reduction in transmissability from those on ART
 real theta;															// Our reduction in transmissability when diagnosed
 int start;																// The start year of our epidemic 
 int diag_start;                                                        // The start year of possible diagnosis 
 int art_start;                                                         // The start year of ART use
 matrix[time_steps_euler, rows(mu_i)] diag;								// Our matrix of rates moving from infected to diagnosed
 vector[3] art_prog;                                                    // Progression through the CD$ stages when on ART
 real onem;                                                             // Birth rate into susceptible population.
  }
  
 
 
 
parameters{
 
 real<lower = 0, upper=1> iota;                                         // The proportion of the population initially infected 
 
 vector<lower = 0, upper = 1>[knot_number] beta;                        // This is the knot point values
 
 real<lower = 0, upper = 1> sigma_pen;                                  // This is the penalty to apply to the spline to make it smooth
 
 
 }

 
transformed parameters{
 matrix[size(rows_to_interpret), 7] y_hat;                                         
 
 y_hat = simpleepp_art_diag( X_design * beta , iota, alpha, mu, sigma, mu_i, mu_d, mu_a, omega, theta, dt_2, start, diag_start, art_start, diag, art_prog, onem)[rows_to_interpret, ];
 
 }

 
model {
 
 
 
 iota ~ normal(0, 0.5);                                                  // This models our initial population size 
 
 target += normal_lpdf( D_penalty * beta | 0, sigma_pen);                // So this models our penalized spline with a slightly altered distribution command
 
  for( i in 1:n_obs){
  y[i] ~ binomial(n_sample ,y_hat[i, 7]);                                // This fits it to the binomially sampled data 
  }

} 

generated quantities{

matrix[time_steps_euler, 7] fitted_output;                               // This models the data after our fitting procedure


fitted_output = simpleepp_art_diag( X_design * beta , iota, alpha, mu, sigma, mu_i, mu_d, mu_a, omega, theta, dt_2, start, diag_start, art_start, diag, art_prog, onem);   
} 


