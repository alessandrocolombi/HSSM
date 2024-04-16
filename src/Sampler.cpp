#include "Sampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>

// Consrtructor
Sampler::Sampler( const GS_data & __gs_data, unsigned int seed ):gs_data(__gs_data), random_engine(seed) {

	//Initialize Full Conditional Objects
    		//auto gamma_ptr  = std::make_shared<FC_gamma>("gamma", h1, h2, pow, gs_data.d, adapt_var0, a1, b1, FixedGamma);
    		//auto beta_ptr   = std::make_shared<FC_beta>();
    		//auto lambda_ptr = std::make_shared<FC_Lambda>("lambda", a2, b2, FixedLambda);

    //Full Conditional vector that we will loop
    		//std::vector< std::shared_ptr<FullConditional> > fc{gamma_ptr, beta_ptr,lambda_ptr};
    		//std::swap(FullConditionals, fc);

    // Initialize return structure 
    // burn_in + n_iter*thin should be just n_iter, that is the number of values I want as output
    	//out.S.reserve(burn_in + n_iter*thin);
    	//out.mu.reserve(burn_in + n_iter*thin);
    	//out.sigma.reserve(burn_in + n_iter*thin);
}
