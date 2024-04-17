#ifndef ___FC_BETA_H___
#define ___FC_BETA_H___

#include "FullConditional.h"
//#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_beta: public FullConditional {
protected:
    double log_full_beta( const double beta,
                          const std::vector<double>& gamma_j,
                          const std::vector<double>& x_j,
                          double sigma2, double sigma2_beta);
public:

    FC_beta(std::string _na, bool _keepfixed);

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

    ~FC_beta() {};



};


#endif
