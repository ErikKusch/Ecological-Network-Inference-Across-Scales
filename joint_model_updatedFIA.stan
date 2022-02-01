// ####################################################################### #
// PROJECT: [PhD; INTRINSIC FITNESS STAN MODELS] 
// CONTENTS: 
// 	- STAN model code for execution of joint individual fitness (IF) and response-effect (RE) models
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

	// DATA DIMENSIONS ----------------------------------------------------------------
	int<lower=1> S;          								// number of species (elements) 
	int<lower=1> N;          								// number of observations (rows in model matrix)
	int<lower=0> K;          								// number of neighbours (columns in model matrix)
	int<lower=0> I;          								// number of realised interactions 
  
	int<lower=0> species_ID[N];	   							// index matching species to observations

	// INDICES MATCHING SPECIES TO INTERACTIONS
	int<lower=0> istart[S];       								// for indexing of icol and irow, records where a new actor-species is recorded	
	int<lower=0> iend[S];									// for indexing of icol and irow, records where the last of an actor-species is recorded 
				// istart = 1; iend = 21 indicates that the first 20 elements of icol and irow target species 1 (irow = 1)
	int<lower=0> icol[I];									// for indexing of observed interactions, records column of interaction partner in model matrix
	int<lower=0> irow[I];									// for indexing of observed interactions, records row of interaction partner in model matrix
				// icol = 1; irow = 2 indicates that there is an interaction between species 1 and 2


	// DATA ---------------------------------------------------------------------------
	real<lower=0> fitness[N]; 	     							// response variable 
	matrix[N,K] X;         									// neighbour abundances (the model matrix)

} 

// ====================================================================================================================
// PARAMETER BLOCK - where we declare model parameters

parameters {
  
	// IF MODEL PARAMETERS -------------------------------------------------------------
	vector[S] beta_i0;    									// species-specific intercept 
	vector[I] beta_ij;  									// vector of interactions which have been realised
	vector<lower=0>[S] sigma; 								// species-specific dispersion deviation parameter

	// RE MODEL PARAMETERS ------------------------------------------------------------
	vector<lower=0>[1] response1; 								// competitive response parameter; >= 0 to avoid bimodality in response and effect  
	vector[S-1] responseSm1;								// other species-specific response parameters
	unit_vector[K] effect; 									// competitive effect parameter; each species has the same effect on all other species, can be facilitative (-) or competitive (+)
	// real<lower=0> sigma;									// response-effect deviation, scale parameter for the logistic distribution (used to estimate re's)
} 

// ====================================================================================================================
// TRANSFORMED PARAMETERS - unknown parameters which can be computed given the known parameters in the parameters block; saved in the output

transformed parameters {

	// IF MODEL PARAMETERS ------------------------------------------------------------  
  	vector[N] mu2;             // the NDDM linear predictor for perform (here seed production)
        matrix[S, K] ndd_betaij;  // interaction matrix for the NDDM							// the linear predictor for fitness

	// RE MODEL PARAMETERS ------------------------------------------------------------
	vector[N] mu;              // the RIM linear predictor for perform (here seed production)
	vector[S] response;        								// combined vector of species-specific responses
  	matrix[S, K] ri_betaij;   // interaction matrix for the RIM

	// stitch together the response values
	  response = append_row(response1, responseSm1);
 	 // get RIM interactions
 	 ri_betaij = response*effect';
   
	  // RIM estimates all interactions
	  for(n in 1:N) {
	       mu[n] = beta_i0[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]);  
	  }
	  
	  // NDDM estimates inferrable interactions, and uses RIM estimates when non-inferrable:
  
	  ndd_betaij = ri_betaij; // initialise nddm interaction matrix to rim estimates
	  // match inferrable interactions parameters to the correct position in the interaction matrix
	  for(s in 1:S) {
	    for(i in istart[s]:iend[s]) {
	      ndd_betaij[irow[i], icol[i]] = beta_ij[i];
	   }
	  }
	  // estimate inferrable interactions
	  for(n in 1:N) {
	        mu2[n] = beta_i0[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]);  
	   }											// end of interaction-loop

} 

// ====================================================================================================================
// MODEL BLOCK - priors and likelihood evaluations

model {

	// IF PRIORS ----------------------------------------------------------------------
	beta_i0 ~ cauchy(0,10);   								// prior for the intercept following Gelman 2008
	//disp_dev ~ cauchy(0, 1);  								// safer to place prior on disp_dev than on phi
  	beta_ij ~ normal(0,1);    // prior for interactions inferred by the NDDM

	// RE PRIORS ----------------------------------------------------------------------  
	response1 ~ normal(0, 1);   								// defining a lower limit in parameters block truncates the normal distribution
	responseSm1 ~ normal(0,1);								// defining a lower limit in parameters block truncates the normal distribution
	sigma ~ cauchy(0, 1);									// response-effect deviation
	// no prior needed for effect as we can use the default prior for the unit_vector

	// RESPONSE EFFECT MODEL LIKELIHOOD ---------------------------------------------
	for(n in 1:N) {										// observation-loop; loop over all observations
		fitness[n] ~ lognormal(mu[n], sigma[species_ID[n]]);					// lognormal outcome distribution, individual mu for each observation and species-specific dispersion parameter
	}											// end of observation loop

	  // add likelihood for the NDDM over inferrable interactions only.
	  for (n in 1:N) {
	    target += lognormal_lpdf(fitness[n] | mu2[n], sigma[species_ID[n]]);
	  }											// end of interaction loop
  
}
