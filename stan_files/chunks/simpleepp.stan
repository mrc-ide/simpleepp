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

matrix simpleepp(vector f_t, real iota, vector alpha, real mu, vector sigma,
		 vector mu_i, vector mu_a, real omega, real dt, int foi_flag){

  vector[rows(kappa)+1] S;
  matrix[rows(mu_i), rows(kappa)+1] I;
  matrix[rows(mu_i), rows(kappa)+1] A;
  vector[rows(kappa)+1] rho;
  vector[rows(kappa)+1] kappa;
  vector[rows(kappa)+1] lambda;
  int DS;
  DS = rows(mu_i);

  // initial values

  if(foi_flag == 0) {
    kappa[1] = f_t[1];
    S[1] = 1000 * (1 - iota);
    I[1, 1] = 1000 * iota;
  } else if(foi_flag == 1) {
    kappa[1] = 0;
    S[1] = 1000;
    I[1, 1] = 0;
  }
    
  for(m in 2:DS)
    I[m, 1] = 0;
  for(m in 1:DS)
    A[m, 1] = 0;
  
  rho[1] = 1 - S[1] / (S[1] + sum(I[,1]) + sum(A[,1]));
  lambda[1] = 0;

  for(t in 1:rows(kappa)) {
    real artcov;
    real deaths;
    real infections;
    real It;
    real At;

    It = sum(I[,t]);
    At = sum(A[,t]);

    artcov = At / (It + At);
    if(foi_flag == 0){
      kappa[t+1] = f_t[t];
      lambda[t+1] = kappa[t+1] * rho[t] * (1 - omega * artcov);
    } else if(foi_flag == 1) {
      lambda[t+1] = f_t[t];
      kappa[t+1] = lambda[t+1] / (rho[t] * (1 - omega * artcov));
    }

    deaths = mu * (S[t] + It + At) + sum(mu_i .* I[,t]);
    S[t+1] = S[t] + dt*(-lambda[t+1] * S[t] - mu * S[t] + deaths);

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS]) * I[DS, t]);
    
    for(m in 1:DS)
      A[,t+1] = A[,t] + dt * mu_a .* A[,t];
    
    rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]) + sum(A[ ,t+1]));
  }
  
  return(append_col(append_row(kappa[1], kappa),
		    append_col(lambda, rho)));
}



/**
 * This is a simplfied version of the four parameter 'EPP classic' model. The population size
 * is assumed to be constant, disease progression is simplified to progression through four
 * stages, and ART is not modelled.
 * @param r Initial epidemic growth rate.
 * @param f0 Proportion of the population initially at risk.
 * @param iota Proportion initially infected.
 * @param phi Behavioural response, determines proportion who enter 'at risk' group.
 * @param mu Natural (non-HIV) rate of mortality / exiting population.
 * @param sigma Rate of progression through each infection stage (vector length m-1).
 * @param mu_i HIV mortality rate by infection stage (vector length m).
 * @param dt Time step for Euler integration.
 * @return An array with columns transmission rate, incidence, and prevalence.
 */

matrix simpleepp_classic(real r, real f0, real iota, real phi, real mu, vector sigma,
			 vector mu_i, real dt, int nsteps){

  vector[nsteps+1] X;  // not-at-risk population
  vector[nsteps+1] Z;  // at-risk population
  matrix[rows(mu_i), nsteps+1] I;
  vector[nsteps+1] rho;
  vector[nsteps+1] lambda;
  vector[nsteps] kappa;
  int DS;
  DS = rows(mu_i);

  // initial values
  X[1] = 1000 * (1 - f0);
  Z[1] = 1000 * (f0 - iota);
  I[1, 1] = 1000 * iota;
  for(m in 2:DS)
    I[m, 1] = 0;
  
  rho[1] = 1 - (X[1] + Z[1]) / (X[1] + Z[1] + sum(I[,1]));
  lambda[1] = 0;

  for(t in 1:nsteps) {
    real deaths;
    real infections;
    real It;
    real St;
    real Nt;
    real ft;

    It = sum(I[,t]);
    St = X[t] + Z[t];
    Nt = St + It;

    ft = exp(phi * X[t] / Nt - (1-f0)) / (exp(phi * X[t] / Nt - (1-f0)) - 1 + 1/f0);
    lambda[t+1] = r * rho[t] * Z[t]/St;
    kappa[t] = lambda[t+1] / rho[t];

    deaths = mu * Nt + sum(mu_i .* I[,t]);
    X[t+1] = X[t] + dt*(-mu * X[t] + (1 - ft) * deaths);
    Z[t+1] = Z[t] + dt*(-lambda[t+1] * St - mu*Z[t] + ft * deaths);

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * St - (mu + mu_i[1] + sigma[1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS]) * I[DS, t]);

    rho[t+1] = 1 - (X[t+1] + Z[t+1]) / (X[t+1] + Z[t+1] + sum(I[,t+1]));
  }
  
  return(append_col(append_row(kappa[1], kappa),
		    append_col(lambda, rho)));
}
