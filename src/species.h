#ifndef __SPECIES_HPP__
#define __SPECIES_HPP__

#define STRICT_R_HEADERS

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>      // For progress bar
#include <progress_bar.hpp>

#include "GSL_wrappers.h"
#include <gsl/gsl_sf.h>
#include "ComponentPrior_factory.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------

// These functions computes log(sum_i(a_i)) using a stable formula for log values. Let us denote a* to the the maximum value of vector a which is attained when i = i*.
// ---> log(sum_i(a_i)) = log(a*) + log[ 1 + sum_{i not i*}(exp{log(a_i) - log(a*)}) ]
// See that only the logarithm of the elements of a are needed. Hence, it is likely that one has already computed them in log scale. If so, set is_log = T
// In this version of the function, the maximum and its position are passed as parameters. No check is performed to be sure that such information were true or not.
// It is assumed that the value in val_max is on the same scale of the values of a, i.e it is in log scale if is_log is set to true.
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max, const unsigned int& idx_max);

// As before but gets an iterator poining to the maximum value
double log_stable_sum(const std::vector<double>& a, const bool is_log, const HSSM_Traits::vector_d_citerator& it_max);

// In this specialized version of the function, the position of the max value is not provided. Hence, one additional operation is done. 
// Since exp( log(a*) - log(a*)) = 1, the previous formula becomes 
// ---> log(sum_i(a_i)) = log(a*) + log[ sum_{i}(exp{log(a_i) - log(a*)}) ]
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max);

// In this version of the formula, the maximum value is computed
double log_stable_sum(const std::vector<double>& a, const bool is_log);

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Factorials and Pochammer
//------------------------------------------------------------------------------------------------------------------------------------------------------


//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double raising_factorial(const unsigned int& n, const double& a);

//' log Raising Factorial
//'
//' \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double log_raising_factorial(const unsigned int& n, const double& a);

// NOTE: this is only for testing the best implementation. This function is not used.
// Raising Factorial
//
// \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
// \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
// The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
// \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double raising_factorial_poch(const unsigned int& n, const double& a);

// NOTE: this is only for testing the best implementation. This function is not used.
// log Raising Factorial
//
// \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
// See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details.
double log_raising_factorial_poch(const unsigned int& n, const double& a);


//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double my_falling_factorial(const unsigned int& n, const double& a);

//' log Falling Factorial
//'
//' \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double my_log_falling_factorial(const unsigned int& n, const double& a);

// Falling Factorial
//
// \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
// The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
// \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double my_falling_factorial_old(const unsigned int& n, const double& a);

// log Falling Factorial
//
// \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
// See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details.
double my_log_falling_factorial_old(const unsigned int& n, const double& a);

//' Pochhammer Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol,
//' \mjsdeqn{(a)^{x} = \frac{\Gamma(a+x)}{\Gamma(a)}}.
//' Where \code{x} is a real number. When x is an integer, such a function coincides with the rising factorial defined in \code{\link{raising_factorial}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double compute_Pochhammer(const unsigned int& x, const double& a);

//' Pochhammer log Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol in log form. See \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double compute_log_Pochhammer(const unsigned int& x, const double& a);


//' Zeta Riemann function in log scale
//'
//' @export
// [[Rcpp::export]]
double log_zeta_Riemann(double s);
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	C numbers
//------------------------------------------------------------------------------------------------------------------------------------------------------


// NOTE: this is only for testing the best implementation. This function is not used.
// Build matrix of logC numbers
//
// This is the recursive function called by \code{\link{my_logC}} that builds the matrix containing all the log(|C(n,k)|) numbers.
// It gets as input the element (n,k) to build, the scale s and location r (here defined as positive and non-negative numbers) and the
// matrix res to be build. The matrix is constructed diagonal by diagonal, starting from the bottom.
// Important remark, note that log(|C(n,0)|) = log(raising factorial (n,r)).
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res);

// NOTE: this is only for testing the best implementation. This function is not used.
// My logC
//
//  This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides.
//  It returns an (n+1 x n+1) matrix containing all the log(|C(nn,k)|) numbers, for nn = 0,...,n+1 and k = 0,...,nn.
//  scale and location must be negative and non-positive, respectively.
//  As a consequence, it is memory expensive.
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
my_logC(const unsigned int& n, const double& scale, const double& location);

// NOTE: this is only for testing the best implementation. This function is not used.
// Compute log of absolute values of non Central C number
//
// \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
// It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n.
// scale and location must be negative and non-positive, respectively.
// It uses eigen objects, apparetly it is slower than using Rcpp vectors.
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location);

// NOTE: this is only for testing the best implementation. This function is not used.
// My logC2 - central
//
// This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides for central numbers.
// This is the specialization of \code{\link{my_logC2}} for central C numbers.
// scale must be negative.
Eigen::VectorXd my_logC2_central(const unsigned int& n, const double& scale);

// This is the implementation used in the code
//' compute_logC - Compute log of absolute values of non Central C number
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n.
//' This implementation uses Rcpp vectors.
//' @param scale must be strictly negative.
//' @param locatio must be non positive. Set to 0 for central C numbers.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location);


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Approximations for V
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' All terms of V
//'
//' This function returns the log of all terms in V, from m=0 up to M_max
//' @export
// [[Rcpp::export]]
std::vector<double> log_Vprior_long(const unsigned int& k, const std::vector<unsigned int>& n_i, 
									const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, 
									unsigned int M_max );

// log_Vprior_apprx2 pesca valori vicino alla moda di q_M ma esplode quando gamma
// diminuisce molto o n_j aumenta molto. Per questi casi difficili, log_Vprior_apprx3 funziona meglio
//
// log_Vprior_apprx1 è totalemente inutilizzabile


//' V approximation - Zeta Riemann
//'
//' Non usare questa funzione. Non ci sono funzioni che implementino il log della zeta-Riemann. Quindi ho problemi numerici e non va
//' @export
// [[Rcpp::export]]
int log_Vprior_apprx1(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
									  const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, 
									  unsigned int M_max );

//' V approximation - Under threshold
//'
//' This approximation has no theoretical guarantees
//' @export
// [[Rcpp::export]]
int log_Vprior_apprx2(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, 
						unsigned int M_max );

// Overloaded function
int log_Vprior_apprx2(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max );


//' V approximation - Compute Upper tail
//'
//' Non posso calcolare direttamente il log di P(M>N+K). Usala solo se gamma è molto piccolo o n_j è molto grande.
//' @export
// [[Rcpp::export]]
int log_Vprior_apprx3(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, 
						unsigned int M_max );


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A priori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------


// Da Testare (ma tanto è inutile, calcolo sempre solo il log)
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, 
						const ComponentPrior& qM, unsigned int M_max = 100 );

// Works for d>2
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );

//Direct formula per d=1 or d=2
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);

//Direct formula per d=1 or d=2 but the vector with C numbers is passed as input
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, 
									const std::vector<double>& gamma,
									const Rcpp::NumericVector& absC1, const Rcpp::NumericVector& absC2);

//Recursive formula for d>2
double compute_Kprior_unnormalized_recursive(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);


//Direct formula per d=1 or d=2
double compute_SK_prior_unnormalized(const unsigned int& k, const unsigned int& s, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);

//Direct formula per d=1 or d=2: C numbers are passed as an input
double compute_SK_prior_unnormalized(const unsigned int& k, const unsigned int& s, 
									 const std::vector<unsigned int>& n_i, const std::vector<double>& gamma,
									 const Rcpp::NumericVector& absC1, const Rcpp::NumericVector& absC2 );

//Recursive formula for d>2
double compute_SK_prior_unnormalized_recursive(const unsigned int& k, const unsigned int& s, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A posteriori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

Rcpp::NumericVector log_qM_post(const unsigned int& m, 
								const ComponentPrior& qM,
								const unsigned int& k, const std::vector<unsigned int>& n_j, 
								const std::vector<double>& gamma_j, double log_V, unsigned int M_max );

//' log qM - Posterior evaluation
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector log_qM_post(const unsigned int& m, 
								const Rcpp::String& prior, const Rcpp::List& prior_param,
								const unsigned int& k, const std::vector<unsigned int>& n_j, 
								const std::vector<double>& gamma_j, double log_V, unsigned int M_max );


std::vector<double> D_log_qM_post(	const ComponentPrior& qM,
									const unsigned int& k, const std::vector<unsigned int>& n_j, 
									const std::vector<double>& gamma_j, double log_V, unsigned int M_max );

// NON USARLA, è BUGGATA
std::vector<double> build_log_qM_post(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, 
									  const ComponentPrior& qM, unsigned int M_max = 100 );

// r are the distinct species in the new sample of size m_i
// k are the distinct species in the new sample of size n_i
double compute_log_Vpost(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );

double compute_log_Vpost_naive(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 	   const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100);

// r are the distinct species in the new sample of size m_i
// k are the distinct species in the new sample of size n_i
// Only for d=1 or d=2
double compute_Kpost_unnormalized(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 		  const std::vector<double>& gamma);

// r are the distinct species in the new sample of size m_i
// k are the distinct species in the new sample of size n_i
// This is for d>2
double compute_Kpost_unnormalized_recursive(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 		 		    const std::vector<double>& gamma);

// r = distinct in new sample of size m_i
// t = shared in new sample of size m_i
// k = distinct in observed sample of size n_i
// This is only for d=2 
double compute_SK_post_unnormalized(const unsigned int& r, const unsigned int& t, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 		    const std::vector<double>& gamma);


// r = distinct in new sample
// sigma = shared in new sample
// k = distinct in observed sample
//Recursive formula for d>2
double compute_SK_post_unnormalized_recursive(const unsigned int& r, const unsigned int& sigma, const unsigned int& k, const std::vector<unsigned int>& m_i, 
											  const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// The same code was repeated multiple times in different functions. Hence it has been collected in a single function.
// Wrap the call for the ComponentPrior from Rcpp objects to c++ object
std::unique_ptr< ComponentPrior > Wrapper_ComponentPrior(const Rcpp::String& prior, const Rcpp::List& prior_param);

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing probabilities
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' 
// [[Rcpp::export]] 
double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
							  const Rcpp::List& prior_param, unsigned int M_max  );

//' 
// [[Rcpp::export]] 
double p_shared_prior_c(const unsigned int& s, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
					 	const Rcpp::List& prior_param, unsigned int M_max  );


//' 
// [[Rcpp::export]] 
double p_distinct_posterior_c(const unsigned int& r, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, 
						      const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max );

//' 
// [[Rcpp::export]] 
double p_shared_posterior_c(const unsigned int& t, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, 
						    const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max );

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing expected values
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' 
// [[Rcpp::export]] 
Rcpp::List Expected_prior_c(const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& type, const Rcpp::String& prior, 
					    const Rcpp::List& prior_param, unsigned int M_max, double tol  );


//' 
// [[Rcpp::export]] 
Rcpp::List Expected_posterior_c(const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, 
						    const Rcpp::String& type, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, double tol);


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper bounds
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' 
// [[Rcpp::export]]
Rcpp::NumericVector Sums_logC(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j);

//' 
// [[Rcpp::export]] 
Rcpp::NumericVector UpperBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, 
								  const Rcpp::List& prior_param, unsigned int M_max  );

//' 
// [[Rcpp::export]] 
Rcpp::NumericVector LowerBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, 
								  const Rcpp::List& prior_param, unsigned int M_max  );


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Compute whole distributions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Compute the whole distribution for the prior number of distinct components
//' 
// [[Rcpp::export]] 
Rcpp::NumericVector D_distinct_prior_c( const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, 
										const Rcpp::String& prior, const Rcpp::List& prior_param, 
										unsigned int M_max, const int& Kstart, std::vector<double>& logV_vec   );

// Compute the whole joint distribution of the prior number of distinct and shared components
//' 
// [[Rcpp::export]] 
Rcpp::NumericMatrix D_joint_prior_c( const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, 
									 const Rcpp::String& prior, const Rcpp::List& prior_param, 
									 unsigned int M_max, const int& Kstart, std::vector<double>& logV_vec   );

// Compute MCMC for Prior
//' 
Rcpp::List  Distinct_Prior_MCMC( unsigned int Niter,
				                 const std::vector<unsigned int>& n_j,
				                 const std::vector<double>& gamma_j,
				                 const Rcpp::String& prior, const Rcpp::List& prior_param,
				                 unsigned int M_max = 1000,
				                 unsigned int seed = 0
				               );

// Compute MCMC for Prior
//' 
// [[Rcpp::export]] 
Rcpp::List Distinct_Prior_MCMC_c( Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& N,
								  Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> N_k,
								  const std::vector<unsigned int>& rho0,
								  const unsigned int& K0,
								  unsigned int Niter,
				                  const std::vector<unsigned int>& n_j,
				                  const std::vector<double>& gamma_j,
				                  const Rcpp::String& prior, const Rcpp::List& prior_param,
				                  unsigned int M_max,
				                  unsigned int seed
				                );

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Large n regime
//------------------------------------------------------------------------------------------------------------------------------------------------------

double log_factorial_mom(const int& r, const std::vector<double>& v_log_qM_post);

// ---- K post 
double p_distinct_post_largen(	const unsigned int& r, const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const std::vector<double>& v_log_qM_post );

double p_distinct_post_largen(	const unsigned int& r, const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
						      	double log_V, unsigned int M_max );

Rcpp::NumericVector D_distinct_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const std::vector<double>& v_log_qM_post );

Rcpp::NumericVector D_distinct_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
									      	double log_V, unsigned int M_max );

// ---- (K,S) post 
double p_jointKS_post_largen(	const unsigned int& r, const unsigned int& t,
								const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const std::vector<double>& v_log_qM_post );

double p_jointKS_post_largen(	const unsigned int& r, const unsigned int& t,
								const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
						      	double log_V, unsigned int M_max );

Rcpp::NumericMatrix D_jointKS_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const std::vector<double>& v_log_qM_post );

Rcpp::NumericMatrix D_jointKS_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
									      	double log_V, unsigned int M_max );

// ---- S post 
double p_shared_post_largen(	const unsigned int& t, const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const std::vector<double>& v_log_qM_post );

double p_shared_post_largen(	const unsigned int& t, const unsigned int& k, 
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						      	const std::vector<double>& gamma_j, 
						      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
						      	double log_V, unsigned int M_max );

Rcpp::NumericVector D_shared_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const std::vector<double>& v_log_qM_post );

Rcpp::NumericVector D_shared_post_largen(	const unsigned int& k, 
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
									      	const std::vector<double>& gamma_j, 
									      	const Rcpp::String& prior, const Rcpp::List& prior_param, 
									      	double log_V, unsigned int M_max );

// ---- K post - moments
double Kpost_mom1_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						   	const std::vector<double>& gamma_j, 
						   	const std::vector<double>& v_log_qM_post );

double Kpost_mom1_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						    const std::vector<double>& gamma_j, 
						    const Rcpp::String& prior, const Rcpp::List& prior_param, 
						    double log_V, unsigned int M_max );

double Kpost_mom2_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						   	const std::vector<double>& gamma_j, 
						   	const std::vector<double>& v_log_qM_post );

double Kpost_mom2_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						    const std::vector<double>& gamma_j, 
						    const Rcpp::String& prior, const Rcpp::List& prior_param, 
						    double log_V, unsigned int M_max );

double Kpost_var_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						   	const std::vector<double>& gamma_j, 
						   	const std::vector<double>& v_log_qM_post );

double Kpost_var_largen(	const unsigned int& k, 
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j, 
						    const std::vector<double>& gamma_j, 
						    const Rcpp::String& prior, const Rcpp::List& prior_param, 
						    double log_V, unsigned int M_max );

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Tests
//------------------------------------------------------------------------------------------------------------------------------------------------------


// This function computes prod_{i=1}^{n}( f(a_i*b_i) )
// Questo è solo un test, sarebbe carinio farla che prende in input anche la funzione f() generica da applicare
//double combined_product(const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const unsigned int& Mstar, const unsigned int& k){
//
	//return (
				//std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 1.0, std::multiplies<>(),
	    					   		//[&Mstar, &k](const double& nj, const double& gamma_j){return (  (double)nj + gamma_j*(double)(Mstar + k)  );} 
	    					   	  //)
			//);
//}
//
//// This function computes sum_{i=1}^{n}( f(ni_i*gamma_i) )
//// Questo è solo un test, sarebbe carinio farla che prende in input anche la funzione f() generica da applicare
//double combined_sum(const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const unsigned int& Mstar, const unsigned int& k){
//
	//return (
				//std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
	    					   		//[&Mstar, &k](const double& nj, const double& gamma_j){return (  (double)nj + gamma_j*(double)(Mstar + k)  );} 
	    					   	  //)
			//);
//}



//' Test ComponentPrior
//' @export
// [[Rcpp::export]]
void Test_Prior();


//' Test prod_sum
//' @export
// [[Rcpp::export]]
void Test_prod_sum();


#endif
