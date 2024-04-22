// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MCMC_Sampler_c
Rcpp::List MCMC_Sampler_c(const Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& N__jk, const std::vector<unsigned int>& n__j, unsigned int r, unsigned int niter, unsigned int nburn, unsigned int thin, const Rcpp::List& option);
RcppExport SEXP _HSSM_MCMC_Sampler_c(SEXP N__jkSEXP, SEXP n__jSEXP, SEXP rSEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP thinSEXP, SEXP optionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& >::type N__jk(N__jkSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n__j(n__jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type r(rSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type option(optionSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC_Sampler_c(N__jk, n__j, r, niter, nburn, thin, option));
    return rcpp_result_gen;
END_RCPP
}
// raising_factorial
double raising_factorial(const unsigned int& n, const double& a);
RcppExport SEXP _HSSM_raising_factorial(SEXP nSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(raising_factorial(n, a));
    return rcpp_result_gen;
END_RCPP
}
// log_raising_factorial
double log_raising_factorial(const unsigned int& n, const double& a);
RcppExport SEXP _HSSM_log_raising_factorial(SEXP nSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(log_raising_factorial(n, a));
    return rcpp_result_gen;
END_RCPP
}
// my_falling_factorial
double my_falling_factorial(const unsigned int& n, const double& a);
RcppExport SEXP _HSSM_my_falling_factorial(SEXP nSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(my_falling_factorial(n, a));
    return rcpp_result_gen;
END_RCPP
}
// my_log_falling_factorial
double my_log_falling_factorial(const unsigned int& n, const double& a);
RcppExport SEXP _HSSM_my_log_falling_factorial(SEXP nSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(my_log_falling_factorial(n, a));
    return rcpp_result_gen;
END_RCPP
}
// compute_Pochhammer
double compute_Pochhammer(const unsigned int& x, const double& a);
RcppExport SEXP _HSSM_compute_Pochhammer(SEXP xSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Pochhammer(x, a));
    return rcpp_result_gen;
END_RCPP
}
// compute_log_Pochhammer
double compute_log_Pochhammer(const unsigned int& x, const double& a);
RcppExport SEXP _HSSM_compute_log_Pochhammer(SEXP xSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log_Pochhammer(x, a));
    return rcpp_result_gen;
END_RCPP
}
// log_zeta_Riemann
double log_zeta_Riemann(double s);
RcppExport SEXP _HSSM_log_zeta_Riemann(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(log_zeta_Riemann(s));
    return rcpp_result_gen;
END_RCPP
}
// compute_logC
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location);
RcppExport SEXP _HSSM_compute_logC(SEXP nSEXP, SEXP scaleSEXP, SEXP locationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const double& >::type location(locationSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_logC(n, scale, location));
    return rcpp_result_gen;
END_RCPP
}
// log_Vprior_long
std::vector<double> log_Vprior_long(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_log_Vprior_long(SEXP kSEXP, SEXP n_iSEXP, SEXP gammaSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Vprior_long(k, n_i, gamma, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// log_Vprior_apprx1
int log_Vprior_apprx1(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol, const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_log_Vprior_apprx1(SEXP kSEXP, SEXP n_iSEXP, SEXP tolSEXP, SEXP gammaSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Vprior_apprx1(k, n_i, tol, gamma, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// log_Vprior_apprx2
int log_Vprior_apprx2(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol, const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_log_Vprior_apprx2(SEXP kSEXP, SEXP n_iSEXP, SEXP tolSEXP, SEXP gammaSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Vprior_apprx2(k, n_i, tol, gamma, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// log_Vprior_apprx3
int log_Vprior_apprx3(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol, const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_log_Vprior_apprx3(SEXP kSEXP, SEXP n_iSEXP, SEXP tolSEXP, SEXP gammaSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Vprior_apprx3(k, n_i, tol, gamma, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// log_qM_post
Rcpp::NumericVector log_qM_post(const unsigned int& m, const Rcpp::String& prior, const Rcpp::List& prior_param, const unsigned int& k, const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, double log_V, unsigned int M_max);
RcppExport SEXP _HSSM_log_qM_post(SEXP mSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP kSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP log_VSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< double >::type log_V(log_VSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(log_qM_post(m, prior, prior_param, k, n_j, gamma_j, log_V, M_max));
    return rcpp_result_gen;
END_RCPP
}
// p_distinct_prior_c
double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_p_distinct_prior_c(SEXP kSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(p_distinct_prior_c(k, n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// p_shared_prior_c
double p_shared_prior_c(const unsigned int& s, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_p_shared_prior_c(SEXP sSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(p_shared_prior_c(s, n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// p_joint_prior_c
double p_joint_prior_c(const unsigned int& k, const unsigned int& s, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_p_joint_prior_c(SEXP kSEXP, SEXP sSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(p_joint_prior_c(k, s, n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// p_distinct_posterior_c
double p_distinct_posterior_c(const unsigned int& r, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_p_distinct_posterior_c(SEXP rSEXP, SEXP kSEXP, SEXP m_jSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type m_j(m_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(p_distinct_posterior_c(r, k, m_j, n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// p_shared_posterior_c
double p_shared_posterior_c(const unsigned int& t, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_p_shared_posterior_c(SEXP tSEXP, SEXP kSEXP, SEXP m_jSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type m_j(m_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(p_shared_posterior_c(t, k, m_j, n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// Expected_prior_c
Rcpp::List Expected_prior_c(const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& type, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, double tol);
RcppExport SEXP _HSSM_Expected_prior_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP typeSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Expected_prior_c(n_j, gamma_j, type, prior, prior_param, M_max, tol));
    return rcpp_result_gen;
END_RCPP
}
// Expected_posterior_c
Rcpp::List Expected_posterior_c(const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& type, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, double tol);
RcppExport SEXP _HSSM_Expected_posterior_c(SEXP kSEXP, SEXP m_jSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP typeSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type m_j(m_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Expected_posterior_c(k, m_j, n_j, gamma_j, type, prior, prior_param, M_max, tol));
    return rcpp_result_gen;
END_RCPP
}
// Sums_logC
Rcpp::NumericVector Sums_logC(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j);
RcppExport SEXP _HSSM_Sums_logC(SEXP n_jSEXP, SEXP gamma_jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    rcpp_result_gen = Rcpp::wrap(Sums_logC(n_j, gamma_j));
    return rcpp_result_gen;
END_RCPP
}
// UpperBounds_c
Rcpp::NumericVector UpperBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_UpperBounds_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(UpperBounds_c(n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// LowerBounds_c
Rcpp::NumericVector LowerBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max);
RcppExport SEXP _HSSM_LowerBounds_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(LowerBounds_c(n_j, gamma_j, prior, prior_param, M_max));
    return rcpp_result_gen;
END_RCPP
}
// D_distinct_prior_c
Rcpp::NumericVector D_distinct_prior_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const int& Kstart, std::vector<double>& logV_vec);
RcppExport SEXP _HSSM_D_distinct_prior_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KstartSEXP, SEXP logV_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kstart(KstartSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type logV_vec(logV_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(D_distinct_prior_c(n_j, gamma_j, prior, prior_param, M_max, Kstart, logV_vec));
    return rcpp_result_gen;
END_RCPP
}
// D_distinct_prior_interval_c
Rcpp::NumericVector D_distinct_prior_interval_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const int& Kmin, const int& Kmax, std::vector<double>& logV_vec, bool print);
RcppExport SEXP _HSSM_D_distinct_prior_interval_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KminSEXP, SEXP KmaxSEXP, SEXP logV_vecSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type logV_vec(logV_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(D_distinct_prior_interval_c(n_j, gamma_j, prior, prior_param, M_max, Kmin, Kmax, logV_vec, print));
    return rcpp_result_gen;
END_RCPP
}
// D_joint_prior_c
Rcpp::NumericMatrix D_joint_prior_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const int& Kstart, std::vector<double>& logV_vec);
RcppExport SEXP _HSSM_D_joint_prior_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KstartSEXP, SEXP logV_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kstart(KstartSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type logV_vec(logV_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(D_joint_prior_c(n_j, gamma_j, prior, prior_param, M_max, Kstart, logV_vec));
    return rcpp_result_gen;
END_RCPP
}
// D_joint_prior_square_c
Rcpp::NumericMatrix D_joint_prior_square_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const int& Kmin, const int& Kmax, const int& Smin, const int& Smax, std::vector<double>& logV_vec, bool print);
RcppExport SEXP _HSSM_D_joint_prior_square_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KminSEXP, SEXP KmaxSEXP, SEXP SminSEXP, SEXP SmaxSEXP, SEXP logV_vecSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Smin(SminSEXP);
    Rcpp::traits::input_parameter< const int& >::type Smax(SmaxSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type logV_vec(logV_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(D_joint_prior_square_c(n_j, gamma_j, prior, prior_param, M_max, Kmin, Kmax, Smin, Smax, logV_vec, print));
    return rcpp_result_gen;
END_RCPP
}
// Distinct_Prior_MCMC_c
Rcpp::List Distinct_Prior_MCMC_c(Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& N, Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> N_k, const std::vector<unsigned int>& rho0, const unsigned int& K0, unsigned int Niter, const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, unsigned int seed);
RcppExport SEXP _HSSM_Distinct_Prior_MCMC_c(SEXP NSEXP, SEXP N_kSEXP, SEXP rho0SEXP, SEXP K0SEXP, SEXP NiterSEXP, SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> >::type N_k(N_kSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type rho0(rho0SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(Distinct_Prior_MCMC_c(N, N_k, rho0, K0, Niter, n_j, gamma_j, prior, prior_param, M_max, seed));
    return rcpp_result_gen;
END_RCPP
}
// D_joint_post_square_c
Rcpp::NumericMatrix D_joint_post_square_c(const std::vector<unsigned int>& m_j, const std::vector<double>& gamma_j, const std::vector<unsigned int>& n_j, const unsigned int& k, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const int& Kmin, const int& Kmax, const int& Smin, const int& Smax, std::vector<double>& logVpost_vec, bool print);
RcppExport SEXP _HSSM_D_joint_post_square_c(SEXP m_jSEXP, SEXP gamma_jSEXP, SEXP n_jSEXP, SEXP kSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KminSEXP, SEXP KmaxSEXP, SEXP SminSEXP, SEXP SmaxSEXP, SEXP logVpost_vecSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type m_j(m_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< const int& >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< const int& >::type Smin(SminSEXP);
    Rcpp::traits::input_parameter< const int& >::type Smax(SmaxSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type logVpost_vec(logVpost_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(D_joint_post_square_c(m_j, gamma_j, n_j, k, prior, prior_param, M_max, Kmin, Kmax, Smin, Smax, logVpost_vec, print));
    return rcpp_result_gen;
END_RCPP
}
// Test_Prior
void Test_Prior();
RcppExport SEXP _HSSM_Test_Prior() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Test_Prior();
    return R_NilValue;
END_RCPP
}
// Test_prod_sum
void Test_prod_sum();
RcppExport SEXP _HSSM_Test_prod_sum() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Test_prod_sum();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HSSM_MCMC_Sampler_c", (DL_FUNC) &_HSSM_MCMC_Sampler_c, 7},
    {"_HSSM_raising_factorial", (DL_FUNC) &_HSSM_raising_factorial, 2},
    {"_HSSM_log_raising_factorial", (DL_FUNC) &_HSSM_log_raising_factorial, 2},
    {"_HSSM_my_falling_factorial", (DL_FUNC) &_HSSM_my_falling_factorial, 2},
    {"_HSSM_my_log_falling_factorial", (DL_FUNC) &_HSSM_my_log_falling_factorial, 2},
    {"_HSSM_compute_Pochhammer", (DL_FUNC) &_HSSM_compute_Pochhammer, 2},
    {"_HSSM_compute_log_Pochhammer", (DL_FUNC) &_HSSM_compute_log_Pochhammer, 2},
    {"_HSSM_log_zeta_Riemann", (DL_FUNC) &_HSSM_log_zeta_Riemann, 1},
    {"_HSSM_compute_logC", (DL_FUNC) &_HSSM_compute_logC, 3},
    {"_HSSM_log_Vprior_long", (DL_FUNC) &_HSSM_log_Vprior_long, 6},
    {"_HSSM_log_Vprior_apprx1", (DL_FUNC) &_HSSM_log_Vprior_apprx1, 7},
    {"_HSSM_log_Vprior_apprx2", (DL_FUNC) &_HSSM_log_Vprior_apprx2, 7},
    {"_HSSM_log_Vprior_apprx3", (DL_FUNC) &_HSSM_log_Vprior_apprx3, 7},
    {"_HSSM_log_qM_post", (DL_FUNC) &_HSSM_log_qM_post, 8},
    {"_HSSM_p_distinct_prior_c", (DL_FUNC) &_HSSM_p_distinct_prior_c, 6},
    {"_HSSM_p_shared_prior_c", (DL_FUNC) &_HSSM_p_shared_prior_c, 6},
    {"_HSSM_p_joint_prior_c", (DL_FUNC) &_HSSM_p_joint_prior_c, 7},
    {"_HSSM_p_distinct_posterior_c", (DL_FUNC) &_HSSM_p_distinct_posterior_c, 8},
    {"_HSSM_p_shared_posterior_c", (DL_FUNC) &_HSSM_p_shared_posterior_c, 8},
    {"_HSSM_Expected_prior_c", (DL_FUNC) &_HSSM_Expected_prior_c, 7},
    {"_HSSM_Expected_posterior_c", (DL_FUNC) &_HSSM_Expected_posterior_c, 9},
    {"_HSSM_Sums_logC", (DL_FUNC) &_HSSM_Sums_logC, 2},
    {"_HSSM_UpperBounds_c", (DL_FUNC) &_HSSM_UpperBounds_c, 5},
    {"_HSSM_LowerBounds_c", (DL_FUNC) &_HSSM_LowerBounds_c, 5},
    {"_HSSM_D_distinct_prior_c", (DL_FUNC) &_HSSM_D_distinct_prior_c, 7},
    {"_HSSM_D_distinct_prior_interval_c", (DL_FUNC) &_HSSM_D_distinct_prior_interval_c, 9},
    {"_HSSM_D_joint_prior_c", (DL_FUNC) &_HSSM_D_joint_prior_c, 7},
    {"_HSSM_D_joint_prior_square_c", (DL_FUNC) &_HSSM_D_joint_prior_square_c, 11},
    {"_HSSM_Distinct_Prior_MCMC_c", (DL_FUNC) &_HSSM_Distinct_Prior_MCMC_c, 11},
    {"_HSSM_D_joint_post_square_c", (DL_FUNC) &_HSSM_D_joint_post_square_c, 13},
    {"_HSSM_Test_Prior", (DL_FUNC) &_HSSM_Test_Prior, 0},
    {"_HSSM_Test_prod_sum", (DL_FUNC) &_HSSM_Test_prod_sum, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_HSSM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
