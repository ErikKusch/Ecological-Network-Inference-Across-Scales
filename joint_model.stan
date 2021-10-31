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
	vector<lower=0>[S] disp_dev; 								// species-specific dispersion deviation parameter

	// RE MODEL PARAMETERS ------------------------------------------------------------
	vector<lower=0>[1] response1; 								// competitive response parameter; >= 0 to avoid bimodality in response and effect  
	vector[S-1] responseSm1;								// other species-specific response parameters
	unit_vector[K] effect; 									// competitive effect parameter; each species has the same effect on all other species, can be facilitative (-) or competitive (+)
	real<lower=0> sigma;									// response-effect deviation, scale parameter for the logistic distribution (used to estimate re's)
} 

// ====================================================================================================================
// TRANSFORMED PARAMETERS - unknown parameters which can be computed given the known parameters in the parameters block; saved in the output

transformed parameters {

	// IF MODEL PARAMETERS ------------------------------------------------------------  
	matrix[S, K] inter_mat;									// the community interaction matrix
	vector[N] mu;              								// the linear predictor for fitness

	// RE MODEL PARAMETERS ------------------------------------------------------------
	vector[S] response;        								// combined vector of species-specific responses
	vector[I] re;              								// interactions as calculated by the re model

	// OBSERVED INTERACTION PARAMETERS AND COMMUNITY MATRIX ---------------------------	// match beta_ij to the correct position of observed interactions in model matrix
	inter_mat = rep_matrix(0, S, K); 							// fill the community interaction matrix with 0 (instead of NA)
	for(s in 1:S) {										// species-loop, loop over all species
		for(i in istart[s]:iend[s]) {								// column-loop, loop over all interactions observed for this species
			inter_mat[irow[i], icol[i]] = beta_ij[i];						// update interaction matrix contents with interaction parameter
		}											// end of column-loop
	}											// end of species-loop

	// INDIVIDUAL FITNESS MODEL -------------------------------------------------------
	for(n in 1:N) {										// observation-loop, loop over all observations; this is where additional model variables should be placed
		mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], inter_mat[species_ID[n], ]));  
									// exp() link because outcome distribution is lognormal
									// matrices can be subset for rows by writing "matrix[row]"
									// dot_product() is a STAN function that runs the additive model of each sequential pair of the vectors contained in the function

	}											// end of observation-loop
	response = append_row(response1, responseSm1);						// stitch together the response values

	// RESPONSE EFFECT MODEL INTERACTION PARAMETER VECTOR -----------------------------
	for (i in 1:I) {									// interaction-loop, loop over all observed interactions
		re[i] = response[irow[i]]*effect[icol[i]];						// response-effect outcome is the product of the response and effect of the interaction partners
	}											// end of interaction-loop

} 

// ====================================================================================================================
// MODEL BLOCK - priors and likelihood evaluations

model {

	// IF PRIORS ----------------------------------------------------------------------
	beta_i0 ~ cauchy(0,10);   								// prior for the intercept following Gelman 2008
	disp_dev ~ cauchy(0, 1);  								// safer to place prior on disp_dev than on phi

	// RE PRIORS ----------------------------------------------------------------------  
	response1 ~ normal(0, 1);   								// defining a lower limit in parameters block truncates the normal distribution
	responseSm1 ~ normal(0,1);								// defining a lower limit in parameters block truncates the normal distribution
	sigma ~ cauchy(0, 1);									// response-effect deviation
	// no prior needed for effect as we can use the default prior for the unit_vector

	// INDIVIDUAL FITNESS (IF) LIKELIHOOD ---------------------------------------------
	for(n in 1:N) {										// observation-loop; loop over all observations
		fitness[n] ~ lognormal(mu[n], disp_dev[species_ID[n]]);					// lognormal outcome distribution, individual mu for each observation and species-specific dispersion parameter
	}											// end of observation loop

	// RESPONSE-EFFECT-INTERACTION (RE) LIKELIHOOD ------------------------------------
	for (i in 1:I) {									// interaction-loop, loop over all observed interactions
		target += logistic_lpdf(re[i] | beta_ij[i], sigma);					// evaluation of RE computation given IF interaction effect
	}											// end of interaction loop
  
}
