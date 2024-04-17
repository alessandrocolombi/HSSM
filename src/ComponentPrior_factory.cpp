#include "ComponentPrior_factory.h"

// Creates the actual object using the string passed in the main
std::unique_ptr< ComponentPrior >
Select_ComponentPrior(std::string const & namePr, ComponentPrior_Parameters const & qM_params){
  // Select Prior
  if(namePr == "Poisson")
    return Create_qM<Mpriors::Poi1>(namePr, qM_params.Lambda);
  else if(namePr == "NegativeBinomial")
    return Create_qM<Mpriors::NegBin1>(namePr, qM_params.p, qM_params.n_succ);
  else
    throw std::runtime_error("Error, only possible priors for the number of components are Poi1 and NegBin1");
}
