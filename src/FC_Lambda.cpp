#include "FC_Lambda.h"

FC_Lambda::FC_Lambda(std::string _na, bool _keepfixed) : FullConditional(_na,_keepfixed) {}

void FC_Lambda::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
  const std::vector<double>& gamma_j = gs_data.gamma_j;
  const std::vector<double>& b_gamma = gs_data.b_gamma;

  // Update hyperparameters
  double a_Lambda = (gs_data.Lambda0 * gs_data.Lambda0 )/(gs_data.V_Lambda);
  double b_Lambda = (gs_data.Lambda0)/(gs_data.V_Lambda);
  a_Lambda += (double)gs_data.d * gs_data.a_gamma;
  b_Lambda += std::inner_product( gamma_j.cbegin(),gamma_j.cend(),b_gamma.cbegin(), 0.0, std::plus<>(),
      			                      [](const double& gammaj, const double& b_gammaj){return ( gammaj * b_gammaj   );}
  	    					              );
    // Random sampler is created
    sample::rgamma Gamma;

    gs_data.Lambda = Gamma(gs_engine, a_Lambda, 1.0 /(b_Lambda) );
}
