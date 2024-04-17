#' Set options
#'
#' @export
set_options_sampler = function( X_j = NULL,
                                gamma0 = c(1,1), beta0 = 0,Lambda0 = 0,
                                UpdateLambda = TRUE, UpdateGamma = TRUE, UpdateBeta = TRUE,
                                use_covariates = TRUE,
                                L0 = 100, V_Lambda = 100,
                                sigma2_beta = 100, adapt_var_gamma_j = c(1,1),
                                adapt_var_beta = 1,
                                M_max = 1000,
                                gamma_guess = 1 )
{
  d = length(gamma0)
  sigma2 = (L0*L0)/(d*V_Lambda)
  a_gamma = 1/sigma2

  if(use_covariates){
    if(is.null(X_j))
      stop("Error: if use_covariates is TRUE, X_j can not be NULL")
    if(length(X_j)!=d)
      stop("Error: if use_covariates is TRUE, X_j must have length equal to d")
  }else{
    X_j = rep(1,d)
    UpdateBeta = FALSE
    beta0 = log(gamma_guess)+log(L0)
  }
  b_gamma = a_gamma*exp(-X_j*beta0)
  res = list( "X_j" = X_j,
              "gamma0" = gamma0, "beta0" = beta0,"Lambda0" = Lambda0,
              "UpdateLambda" = UpdateLambda, "UpdateGamma" = UpdateGamma, "UpdateBeta" = UpdateBeta,
              "use_covariates" = use_covariates,
              "L0" = L0, "V_Lambda" = V_Lambda,
              "sigma2" = sigma2,
              "sigma2_beta" = sigma2_beta, "adapt_var_gamma_j" = adapt_var_gamma_j,
              "adapt_var_beta" = adapt_var_beta, "a_gamma" = a_gamma, "b_gamma" = b_gamma,
              "M_max" = M_max)
  return(res)
}
