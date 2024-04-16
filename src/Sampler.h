#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#define STRICT_R_HEADERS

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>      // For progress bar
#include <progress_bar.hpp>

//#include "GSL_wrappers.h"
//#include <gsl/gsl_sf.h>
//#include "ComponentPrior_factory.h"
#include "species.h"

struct GS_data
{
    // sampler settings
    unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;
    std::vector<double> X_j;
    unsigned int d;
    
    // Lambda hyperparam
    double Lambda0;
    double V_Lambda;
    // gamma_j hyperparam
    double sigma2;
    double a_gamma;
    double b_gamma;
    // gamma_j options
    std::vector<double> adapt_var_gamma_j;
    // beta hypeparam
    double sigma2_beta;
    // beta options
    double adapt_var_beta;
};

struct return_obj
{
    std::vector<Rcpp::NumericVector> gammas;
    Rcpp::NumericVector betas;
    Rcpp::NumericVector Lambdas;
    Rcpp::NumericVector logV; 
};

class Sampler {
public:

    void sample();

    // Constructor 
    Sampler( const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & N_jk,
             const GS_data & gs_data );
                
    // Data structure for the output
    return_obj out;
    
    // Get n_j, for output C in a better way
    //std::vector<unsigned int> get_nj() const {return gs_data.n_j;}
    
private:
    //std::vector<std::shared_ptr<FullConditional> > FullConditionals; // vector of shared pointer to FC class
    sample::GSL_RNG random_engine; // GSL random engine to sample from random distribution
    GS_data gs_data; // data structure to store values that are updated during Gibbs Sampler
    //std::string model;
    //void store_params_values();
    //void store_tau();
    //void store_w_jk();
    //void GS_Step();
    //bool Partition_fixed;
};


#endif //GDFMM_GIBBSSAMPLER_H
