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
matrix simpleepp_monitor_N(vector kappa, real iota,  real mu, vector sigma,
		 vector mu_i, real omega, real dt, int foi_flag){

  vector[rows(kappa)+1] S;
  matrix[rows(mu_i), rows(kappa)+1] I;
  matrix[rows(mu_i), rows(kappa)+1] A;
  vector[rows(kappa)+1] rho;
  vector[rows(kappa)+1] lambda;
  int DS;
  vector[rows(kappa)+1] N;
  matrix[rows(kappa)+1, 2] first_mat;
  matrix[rows(kappa)+1, 3] return_mat;  

  
  DS = rows(mu_i);
  
  // initial values

   S[1] = 1000000 * (1 - iota);
   I[1, 1] = 1000000 * iota;
    
  for(m in 2:DS)
    I[m, 1] = 0;
  
  rho[1] = 1 - S[1] / (S[1] + sum(I[,1]) );
  lambda[1] = 0;

  for(t in 1:rows(kappa)) {
    real artcov;
    real deaths;
    real infections;
    real It;
    real At;

    It = sum(I[,t]);
    
    lambda[t+1] = kappa[t] * rho[t];
	deaths = mu * (S[t] + It) + sum(mu_i .* I[,t]);
    
	S[t+1] = S[t] + dt*(-lambda[t+1] * S[t] - mu * S[t] + deaths);

    I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1]) * I[1, t]);
    for(m in 2:(DS-1))
      I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]) * I[m, t]);
    I[DS, t+1] = I[DS, t] + dt*(sigma[DS-1] * I[DS-1, t] - (mu + mu_i[DS]) * I[DS, t]);
    
    rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]));
  }
  
  for(i in 1:rows(kappa)+1)
  N[i] = S[i] + sum(I[,i]);
  
  first_mat = append_col(lambda, rho);
  return_mat = append_col(first_mat, N);
  
  
  return(append_col(append_row(kappa[1], kappa), return_mat));
}

}
