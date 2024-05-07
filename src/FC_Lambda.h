#ifndef ___FC_LAMBDA_H___
#define ___FC_LAMBDA_H___

#include "FullConditional.h"
#include "species.h"
//#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_Lambda: public FullConditional {
protected:
    double log_full_lambda(double x, double logV, double a_Lambda, double b_Lambda);
public:

    FC_Lambda(std::string _na, bool _keepfixed);

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

    ~FC_Lambda() {};



};


#endif
