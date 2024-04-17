#ifndef ___FC_BETA_H___
#define ___FC_BETA_H___

#include "FullConditional.h"
//#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_beta: public FullConditional {
protected:

public:

    FC_beta(std::string _na, bool _keepfixed);

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

    ~FC_beta() {};



};


#endif
