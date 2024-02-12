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
Rcpp::NumericVector D_distinct_prior_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, const unsigned int& Kexact);
RcppExport SEXP _HSSM_D_distinct_prior_c(SEXP n_jSEXP, SEXP gamma_jSEXP, SEXP priorSEXP, SEXP prior_paramSEXP, SEXP M_maxSEXP, SEXP KexactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior_param(prior_paramSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M_max(M_maxSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Kexact(KexactSEXP);
    rcpp_result_gen = Rcpp::wrap(D_distinct_prior_c(n_j, gamma_j, prior, prior_param, M_max, Kexact));
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
    {"_HSSM_raising_factorial", (DL_FUNC) &_HSSM_raising_factorial, 2},
    {"_HSSM_log_raising_factorial", (DL_FUNC) &_HSSM_log_raising_factorial, 2},
    {"_HSSM_my_falling_factorial", (DL_FUNC) &_HSSM_my_falling_factorial, 2},
    {"_HSSM_my_log_falling_factorial", (DL_FUNC) &_HSSM_my_log_falling_factorial, 2},
    {"_HSSM_compute_Pochhammer", (DL_FUNC) &_HSSM_compute_Pochhammer, 2},
    {"_HSSM_compute_log_Pochhammer", (DL_FUNC) &_HSSM_compute_log_Pochhammer, 2},
    {"_HSSM_compute_logC", (DL_FUNC) &_HSSM_compute_logC, 3},
    {"_HSSM_p_distinct_prior_c", (DL_FUNC) &_HSSM_p_distinct_prior_c, 6},
    {"_HSSM_p_shared_prior_c", (DL_FUNC) &_HSSM_p_shared_prior_c, 6},
    {"_HSSM_p_distinct_posterior_c", (DL_FUNC) &_HSSM_p_distinct_posterior_c, 8},
    {"_HSSM_p_shared_posterior_c", (DL_FUNC) &_HSSM_p_shared_posterior_c, 8},
    {"_HSSM_Expected_prior_c", (DL_FUNC) &_HSSM_Expected_prior_c, 7},
    {"_HSSM_Expected_posterior_c", (DL_FUNC) &_HSSM_Expected_posterior_c, 9},
    {"_HSSM_Sums_logC", (DL_FUNC) &_HSSM_Sums_logC, 2},
    {"_HSSM_UpperBounds_c", (DL_FUNC) &_HSSM_UpperBounds_c, 5},
    {"_HSSM_LowerBounds_c", (DL_FUNC) &_HSSM_LowerBounds_c, 5},
    {"_HSSM_D_distinct_prior_c", (DL_FUNC) &_HSSM_D_distinct_prior_c, 6},
    {"_HSSM_Test_Prior", (DL_FUNC) &_HSSM_Test_Prior, 0},
    {"_HSSM_Test_prod_sum", (DL_FUNC) &_HSSM_Test_prod_sum, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_HSSM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
