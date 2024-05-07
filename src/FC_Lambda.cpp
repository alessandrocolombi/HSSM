#include "FC_Lambda.h"

FC_Lambda::FC_Lambda(std::string _na, bool _keepfixed) : FullConditional(_na,_keepfixed) {}

double FC_Lambda::log_full_lambda(double x, double logV, double a_Lambda, double b_Lambda)
{
  double res{logV};
  res += a_Lambda*std::log(x) - x*b_Lambda;
  return res;
}

void FC_Lambda::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
  const std::vector<double>& gamma_j = gs_data.gamma_j;
  const std::vector<double>& b_gamma = gs_data.b_gamma;
  // Random sampler is created
  sample::rgamma Gamma;
  // Sampling functions
  sample::rnorm rnorm;
  sample::runif runif;
  double ww_g = std::pow( (double)(gs_data.iterations + 1), -0.7);
  // Update hyperparameters
  double a_Lambda = (gs_data.Lambda0 * gs_data.Lambda0 )/(gs_data.V_Lambda);
  double b_Lambda = (gs_data.Lambda0)/(gs_data.V_Lambda);
  a_Lambda += (double)gs_data.d * gs_data.a_gamma;
  b_Lambda += std::inner_product( gamma_j.cbegin(),gamma_j.cend(),b_gamma.cbegin(), 0.0, std::plus<>(),
      			                      [](const double& gammaj, const double& b_gammaj){return ( gammaj * b_gammaj   );}
  	    					              );

  // MCMC corretto
  double Lambda_old = gs_data.Lambda;
  // Draw proposed value
  double log_Lambda_prime =  rnorm(gs_engine, std::log(Lambda_old), std::sqrt(gs_data.adapt_var_Lambda));
  double Lambda_prime     =  std::exp(log_Lambda_prime);
  // Compute logVprime
  Rcpp::List prior_param =  Rcpp::List::create( Rcpp::Named("lambda") = Lambda_prime);
  auto qM_ptr = Wrapper_ComponentPrior("Poisson", prior_param);
  ComponentPrior& qM(*qM_ptr);
  double logVprime = compute_log_Vprior(gs_data.Kobs, gs_data.n_j, gs_data.gamma_j, qM, gs_data.M_max);
  // Accept or reject it
  double log_acc = log_full_lambda( Lambda_prime, logVprime,    a_Lambda, b_Lambda  ) -
                   log_full_lambda( Lambda_old,   gs_data.logV, a_Lambda, b_Lambda  ) +
                   log_Lambda_prime - std::log(Lambda_old);
  double acc = std::exp(log_acc);
  double u_acc = runif(gs_engine);
  if( u_acc  < std::min(1.0, acc) ){
    gs_data.logV = logVprime;
    gs_data.Lambda = Lambda_prime;
  }
  //Rcpp::Rcout<<"gs_data.adapt_var_Lambda: Prima"<<std::endl;
  //Rcpp::Rcout<<gs_data.adapt_var_Lambda<<std::endl;
  gs_data.adapt_var_Lambda *=  std::exp(  ww_g * ( std::min(1.0, acc) - 0.44 )  );
  //Rcpp::Rcout<<"Dopo:"<<std::endl;
  //Rcpp::Rcout<<gs_data.adapt_var_Lambda<<std::endl;
  if(std::isnan(gs_data.Lambda)){
    throw std::runtime_error("Error in FC_Lambda::update, there is a nan");
  } 



    //gs_data.Lambda = Gamma(gs_engine, a_Lambda, 1.0 /(b_Lambda) );
}
