#ifndef ___FC_GAMMAS_H___
#define ___FC_GAMMAS_H___

#include "FullConditional.h"
//#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_gamma: public FullConditional {
protected:
    /* METHODS */
    //double log_full_gamma(double x, double Lambda, unsigned int k, unsigned int M_star, const GDFMM_Traits::MatUnsCol & n_jk);
    //double sumlgamma(double x, const GDFMM_Traits::MatUnsCol& n_jk);
    //double l_dgamma(double gamma, double a, double b);

public:

    FC_gamma(std::string _na, bool _keepfixed);

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

    ~FC_gamma() {};
    


};


#endif 
