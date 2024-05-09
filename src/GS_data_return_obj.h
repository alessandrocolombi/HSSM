#ifndef __HSSM_GSDATA_H__
#define __HSSM_GSDATA_H__

//#include "species.h"
#include "include_headers.h"
#include "recurrent_traits.h"

struct GS_data
{
    // data
    Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> N_jk;
    std::vector<double> X_j;
    std::vector<unsigned int> n_j;
    unsigned int d;
    unsigned int Kobs;

    // sampler settings
    unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;
    unsigned int iterations;

    // current values
    std::vector<double> gamma_j;
    double beta;
    double Lambda;
    double logV;

    // Lambda hyperparam
    double Lambda0;
    double V_Lambda;
    // Lambda options
    bool FixedLambda;
    double adapt_var_Lambda;

    // gamma_j hyperparam
    double sigma2;
    double a_gamma;
    std::vector<double> b_gamma;
    // gamma_j options
    std::vector<double> adapt_var_gamma_j;
    bool FixedGammas;

    // beta hypeparam
    double sigma2_beta;
    // beta options
    double adapt_var_beta;
    bool use_covariates;
    bool FixedBeta;

    // auxiliary options
    unsigned int M_max;
    bool dependentPrior;
};

struct return_obj
{
    std::vector<Rcpp::NumericVector> gammas;
    Rcpp::NumericVector betas;
    Rcpp::NumericVector Lambdas;
    Rcpp::NumericVector logV;
};

#endif
