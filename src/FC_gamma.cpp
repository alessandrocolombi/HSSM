#include "FC_gamma.h"

FC_gamma::FC_gamma(std::string _na, bool _keepfixed) : FullConditional(_na,_keepfixed) {};


double FC_gamma::log_full_gamma_j(double gammaj, double Lambda, double a_gamma, double b_gammaj, double logV, const VecUnsCol& n_jk )
{
  double temp{0.0};
  for(unsigned int k = 0; k < n_jk.size(); k++){
    temp += log_raising_factorial(n_jk(k), gammaj);
  }
  double res = (a_gamma - 1.0)*std::log(gammaj) - Lambda*b_gammaj*gammaj + logV + temp;
  return res;
}

void FC_gamma::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine)
{
  // Sampling functions
  sample::rnorm rnorm;
  sample::runif runif;
  // Basic quantities
  double ww_g = std::pow( (double)(gs_data.iterations + 1), -0.7);

  for(unsigned int j = 0; j < gs_data.d; j++){
    double gammaj_old = gs_data.gamma_j[j];
    // Draw proposed value
    double log_gammaj_prime =  rnorm(gs_engine, std::log(gammaj_old), std::sqrt(gs_data.adapt_var_gamma_j[j]));
    double gammaj_prime     = std::exp(log_gammaj_prime);
    // Compute logVprime
    Rcpp::List prior_param = 	Rcpp::List::create( Rcpp::Named("lambda") = gs_data.Lambda);
    auto qM_ptr = Wrapper_ComponentPrior("Poisson", prior_param);
    ComponentPrior& qM(*qM_ptr);
    std::vector<double> gamma_prime_j = gs_data.gamma_j;
    gamma_prime_j[j] = gammaj_prime;
    double logVprime = compute_log_Vprior(gs_data.Kobs, gs_data.n_j, gamma_prime_j, qM, gs_data.M_max);
    // Accept or reject it
    double log_acc = log_full_gamma_j( gammaj_prime, gs_data.Lambda, gs_data.a_gamma, gs_data.b_gamma[j], logVprime, gs_data.N_jk.row(j) ) -
                     log_full_gamma_j( gammaj_old,   gs_data.Lambda, gs_data.a_gamma, gs_data.b_gamma[j], gs_data.logV, gs_data.N_jk.row(j) ) +
                     log_gammaj_prime - std::log(gammaj_old);
    double acc = std::exp(log_acc);
    double u_acc = runif(gs_engine);
    if( u_acc  < std::min(1.0, acc) ){
      gs_data.logV = logVprime;
      gs_data.gamma_j[j] = gammaj_prime;
    }
    gs_data.adapt_var_gamma_j[j] *=  std::exp(  ww_g * ( std::min(1.0, acc) - 0.44 )  );
  }

}
