#' Set options
#'
#' @export
set_options_sampler = function( X_j = NULL,
                                gamma0 = c(1,1), beta0 = 0,Lambda0 = 10,
                                UpdateLambda = TRUE, UpdateGamma = TRUE, UpdateBeta = TRUE,
                                use_covariates = TRUE,
                                L0 = 100, V_Lambda = 100,
                                sigma2_beta = 100, adapt_var_gamma_j = c(1,1),
                                adapt_var_Lambda = 1,
                                adapt_var_beta = 1,
                                dependentPrior = TRUE,
                                M_max = 100, seed = 12345,
                                gamma_guess = 1 )
{
  d = length(gamma0)
  sigma2 = (d*V_Lambda)/(L0*L0)
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
              "UpdateLambda" = UpdateLambda, "UpdateGamma" = UpdateGamma,
              "UpdateBeta" = UpdateBeta,
              "use_covariates" = use_covariates,
              "L0" = L0, "V_Lambda" = V_Lambda,
              "sigma2" = sigma2,
              "sigma2_beta" = sigma2_beta, "adapt_var_gamma_j" = adapt_var_gamma_j,
              "adapt_var_beta" = adapt_var_beta, "adapt_var_Lambda" = adapt_var_Lambda,
              "a_gamma" = a_gamma, "b_gamma" = b_gamma,
              "M_max" = M_max, "seed" = seed,"dependentPrior" = dependentPrior)
  return(res)
}

#' Set options
#'
#' @export
set_options_MCMC = function( gamma0 = c(1,1), Lambda0 = 10,
                             UpdateLambda = TRUE, UpdateGamma = TRUE,
                             L0 = 100, V_Lambda = 100,
                             a_gamma = 1, b_gamma = c(1,1),
                             adapt_var_gamma_j = c(1,1),
                             adapt_var_Lambda = 1,
                             dependentPrior = FALSE,
                             M_max = 100, seed = 12345,
                             gamma_guess = 1 )
{
  d = length(gamma0)

  X_j = rep(1,d)
  UpdateBeta = FALSE
  beta0 = log(gamma_guess)+log(L0)



  res = list( "X_j" = X_j,
              "gamma0" = gamma0, "beta0" = 0,"Lambda0" = Lambda0,
              "UpdateLambda" = UpdateLambda, "UpdateGamma" = UpdateGamma,
              "UpdateBeta" = FALSE,
              "use_covariates" = FALSE,
              "L0" = L0, "V_Lambda" = V_Lambda,
              "sigma2" = 1,
              "sigma2_beta" = 1, "adapt_var_gamma_j" = adapt_var_gamma_j,
              "adapt_var_beta" = 1, "adapt_var_Lambda" = adapt_var_Lambda,
              "a_gamma" = a_gamma, "b_gamma" = b_gamma,
              "M_max" = M_max, "seed" = seed, "dependentPrior" = dependentPrior)
  return(res)
}

#' Run MCMC Sampler
#'
#' @export
MCMC_Sampler = function(counts_jk,niter,nburn = 0,thin = 1,option)
{
  d = nrow(counts_jk)
  r = ncol(counts_jk)
  n =  sum(counts_jk)
  n_j = c(apply(counts_jk,1,sum))
  n_K = c(apply(counts_jk,2,sum))
  if( any(n_K== 0) )
    stop("Error: there are empty columns in counts_jk")
  if(length(n_j)!=length(option$gamma0))
    stop("Error, the length of n_j (group sizes) and gamma has to be equal")
  if(option$use_covariates){
    if(mean(option$X_j) > 1e-10 )
      stop("Error: if use_covariates is TRUE, X_j must have zero mean")
  }

  return( MCMC_Sampler_c(counts_jk,n_j,r,niter,nburn,thin,option) )
}
