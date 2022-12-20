// ####################################################################### #
// PROJECT: [PhD; INTRINSIC FITNESS STAN MODELS] 
// CONTENTS: 
// 	- STAN model code for execution of joint individual fitness (NDD) and response-effect (RIM) models
// AUTHOR: 
//	- Malyon Bimler
//		- Conceptualisation and Coding
//		- https://orcid.org/0000-0003-0059-2360
// 	- Erik Kusch
//		- Adaptation from original Bimler code, restructuring and commenting
//		- erik.kusch@bio.au.dk
//		- https://orcid.org/0000-0002-4984-7646
// ####################################################################### #
  
// ====================================================================================================================

// DATA BLOCK - where we place all of the data that is fed to the model and declare the dimensions of each data object

data {
  int<lower=1> S;          // total number of focal species (s in the manuscript) 
  int<lower=1> N;          // total number of observations (rows in model matrix) (n in the manuscript)
  int<lower=0> T;          // total number of interaction partners (columns in model matrix) (t in the manuscript)
  int<lower=0> I;          // total number of identifiable interactions 
  int<lower=0> P;          // total number of unique plots ++++
  
  int<lower=0> species_ID[N];   // index matching observations to focal species (d in the manuscript)
  int<lower=0> plot_ID[N];   // index matching observations to plot identity ++++
  real<lower=0> perform[N];      // response variable - CHANGE FROM MEE CODE
    
  int<lower=0> icol[I];  // indices matching pairwise inferrable to location in interaction matrix
  int<lower=0> irow[I];
  
  matrix[N,T] X;         // neighbour abundances (the model matrix)

} 

parameters {
  
  vector[S] gamma_i;    // species-specific intercept 
  
  vector<lower=0>[P-1] epsilon_p;  // random effect for plot up to P-1 ++++

  real<lower=0> sigma_p;   // random effect precision ++++
    
  vector<lower=0>[S] sigma; // species-specific scale parameter for a lognormal distribution
    
  vector[I] beta_ij;     // vector of interactions which are inferrable by the NDDM 
    
  real<lower=0> weight;    // weighting value controlling the average strength of interactions (must be positive)
  unit_vector[S] response; // species-specific effect parameter (r in the manuscript but without weight)
  unit_vector[T] effect;   // species-specific effect parameter (e in the manuscript)

} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
 
  vector[N] mu;              // the RIM linear predictor for perform 
  vector[N] mu2;             // the joint model linear predictor for perform 
  
  matrix[S, T] ri_betaij;   // interaction matrix for response - impact estimates (re in the manuscript)
  matrix[S, T] ndd_betaij;  // interaction matrix for joint model interaction estimates (B in the manuscript)
  real epsilon_p1;        // random effect for last plot ++++
  vector[P] epsilon;      // random effect for plot ++++
  
  epsilon_p1 = -sum(epsilon_p);   // random effect for plot[P] = - sum(random effects[1:(P-1)] ++++
  epsilon = append_row(epsilon_p, epsilon_p1); // ++++S
      
  // get RIM interactions
  ri_betaij = weight * response*effect';  // the apostrophe transposes the effect vector
   
  // Estimate response-impact interactions
  for(n in 1:N) {
       mu[n] = gamma_i[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]);  // CHANGE FROM MEE CODE 
  }
  
  // NDDM estimates identifiable interactions, and uses RIM estimates when non-identifiable:
  ndd_betaij = ri_betaij; // initialise nddm interaction matrix to rim estimates
  // match identifiable interactions parameters to the correct position in the interaction matrix
  for(i in 1:I) {
    ndd_betaij[irow[i], icol[i]] = beta_ij[i];
  }

  // estimate identifiable interactions
  for(n in 1:N) {
        mu2[n] = gamma_i[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]);  // CHANGE FROM MEE CODE 
   }
   
} 

model {

  // priors
  gamma_i ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  sigma ~ exponential(1);  // prior for scale CHANGE FROM MEE CODE
  beta_ij ~ normal(0,1);    // prior for interactions inferred by the NDDM
  weight ~ normal(0, 10); // constrained by parameter definition to be positive
  // no prior needed for response or effect as we can use the default prior for the unit_vector
  epsilon_p ~ normal(0, sigma_p); // ++++
  sigma_p ~ normal(0, 10); // ++++
       
  // maximising the likelihood for the joint model 
  for(n in 1:N) {
  
    // maximise the likelihood of the RIM for all observations
    perform[n] ~ lognormal(mu[n] + epsilon[plot_ID[n]], sigma[species_ID[n]]); // CHANGE FROM MEE CODE ++++
    // maximise the likelihood of the NDDM over identifiable interactions
    target += lognormal_lpdf(perform[n] | mu2[n] + epsilon[plot_ID[n]], sigma[species_ID[n]]); // CHANGE FROM MEE CODE ++++
  }
  
  
}

generated quantities {

  vector[N] log_lik_rim;	// log-likelihood for the RIM
  vector[N] log_lik_nddm;	// log-likelihood for the joint model
  
  for(n in 1:N) {	
    log_lik_rim[n] = lognormal_lpdf(perform[n] | mu[n]  + epsilon[plot_ID[n]], sigma[species_ID[n]]); // CHANGE FROM MEE CODE
	
    log_lik_nddm[n] = lognormal_lpdf(perform[n] | mu2[n]  + epsilon[plot_ID[n]], sigma[species_ID[n]]); // CHANGE FROM MEE CODE
  }
  
}
