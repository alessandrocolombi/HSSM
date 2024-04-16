#ifndef ___HSSM_FULLCONDITIONAL_H___
#define ___HSSM_FULLCONDITIONAL_H___

#include <Rcpp.h>
//#include "species.h"
#include "GSL_wrappers.h"
#include "GS_data_return_obj.h"

class FullConditional{
public :
    FullConditional(const std::string& _name, bool _keepfixed) : name(_name), keep_fixed(_keepfixed){};
    FullConditional() = default;
    std::string name;
    bool keep_fixed;
    virtual void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) = 0;
    void print() const;
    bool binary_decision(double p1, const sample::GSL_RNG& engine) const;

    virtual ~FullConditional() = default;
};
#endif
