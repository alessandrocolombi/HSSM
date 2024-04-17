#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#define STRICT_R_HEADERS

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>      // For progress bar
#include <progress_bar.hpp>

//#include "GSL_wrappers.h"
//#include <gsl/gsl_sf.h>
//#include "ComponentPrior_factory.h"
#include "GSL_wrappers.h"
#include "FullConditional.h"
#include "GS_data_return_obj.h"
#include "FC_gamma.h"
#include "FC_Lambda.h"
#include "FC_beta.h"
#include "ComponentPrior_factory.h"
#include "species.h"


class Sampler {
public:

    void sample();

    // Constructor
    Sampler( const GS_data & __gs_data, unsigned int seed );

    // Data structure for the output
    return_obj out;

private:
    std::vector<std::shared_ptr<FullConditional> > FullConditionals; // vector of shared pointer to FC class
    sample::GSL_RNG random_engine; // GSL random engine to sample from random distribution
    GS_data gs_data; // data structure to store values that are updated during Gibbs Sampler
    //std::string model;
    //void store_params_values();
    //void store_tau();
    //void store_w_jk();
    //void GS_Step();
    //bool Partition_fixed;
};

//' Sampler
//'
//' @export
// [[Rcpp::export]]
Rcpp::List MCMC_Sampler_c(const Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& N__jk,
                          const std::vector<double>& X__j,
                          const std::vector<unsigned int>& n__j,
                          unsigned int r,
                          unsigned int niter, unsigned int nburn, unsigned int thin,
                          const Rcpp::List& option);

#endif
