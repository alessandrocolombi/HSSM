#include "FC_beta.h"

FC_beta::FC_beta(std::string _na, bool _keepfixed) : FullConditional(_na,_keepfixed) {};

double FC_beta::log_full_beta( const double beta,
                               const std::vector<double>& gamma_j,
                               const std::vector<double>& x_j,
                               double sigma2, double sigma2_beta)
{
  double temp = std::inner_product( gamma_j.cbegin(),gamma_j.cend(),x_j.cbegin(), 0.0, std::plus<>(),
      			                        [&beta](const double& gammaj, const double& xj){return ( gammaj * std::exp(-xj*beta) );}
  	    					                );
  double res = -0.5 * ( 2.0/sigma2 * temp + (beta*beta)/(sigma2_beta) );
  return res;
}

void FC_beta::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
  // Sampling functions
  sample::rnorm rnorm;
  sample::runif runif;
  // Basic quantities
  double ww_g = std::pow( (double)(gs_data.iterations + 1), -0.7);
  std::vector<double>& gamma_j = gs_data.gamma_j;
  double beta_old   = gs_data.beta;
  // Draw proposed value
  double beta_prime =  rnorm(gs_engine, beta_old, std::sqrt(gs_data.adapt_var_beta));
  // Accept or reject it
  double log_acc = log_full_beta(beta_prime, gamma_j, gs_data.X_j, gs_data.sigma2, gs_data.sigma2_beta) -
                   log_full_beta(beta_old,   gamma_j, gs_data.X_j, gs_data.sigma2, gs_data.sigma2_beta);
  double acc = std::exp(log_acc);
  double u_acc = runif(gs_engine);
  if( u_acc  < std::min(1.0, acc) ){
    gs_data.beta = beta_prime;
    for(unsigned int j = 0; j < gs_data.d; j++)
     gs_data.b_gamma[j] = gs_data.a_gamma*std::exp(-gs_data.X_j[j]*gs_data.beta);
  }


  gs_data.adapt_var_beta *=  std::exp(  ww_g * ( std::min(1.0, acc) - 0.44 )  );

}
