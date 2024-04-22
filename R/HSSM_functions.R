#' p_distinct_prior
#'
#' This function computes the a priori probability that the number of distinct species is equal to \code{k}.
#' @param k integer, the number of distinct species whose probability has to be computed.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_distinct_prior = function(k,n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter is expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters are expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(k<0)
    stop("Error, the number of distinct species k can not be negative")
  if(k==0)
    return (0)
  if(k > sum(n_j))
    return (0)

  # Compute non trivial cases
  return (  p_distinct_prior_c(k,n_j,gamma,prior,prior_params,Max_iter)  )
}


#' p_shared_prior
#'
#' This function computes the a priori probability that the number of shared species is equal to \code{s}.
#' @param s integer, the number of shared species whose probability has to be computed.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_shared_prior = function(s,n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(s<0)
    stop("Error, the number of shared species s can not be negative")
  if(s > min(n_j))
    return (0)



  # Compute non trivial cases
  return (  p_shared_prior_c(s,n_j,gamma,prior,prior_params,Max_iter)  )
}

#' p_joint_prior
#'
#' This function computes the joint prior probability that the number of distinct and shared species is equal to \code{(k,s)}.
#' @param k integer, the number of distinct species whose probability has to be computed.
#' @param s integer, the number of shared species whose probability has to be computed.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_joint_prior = function(k,s,n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(s<0)
    stop("Error, the number of shared species s can not be negative")
  if(s > min(n_j))
    return (0)
  if(s > k)
    return (0)

  # Compute non trivial cases
  return (  p_joint_prior_c(k,s,n_j,gamma,prior,prior_params,Max_iter)  )
}


#' p_distinct_posterior
#'
#' This function computes the probability a posteriori that the number of distinct species is equal to \code{r}.
#' @param r integer, the number of distinct species whose probability has to be computed.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j a positive scalar or a vector of positive integers that defines the size of the groups that have been alread observed. It must be of the same size of \code{m_j}.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_distinct_posterior = function(r, k, m_j, n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if(length(n_j)!=length(m_j))
    stop("The length of m_j must be equal to the length of n_j")
  if( any(m_j<0) || any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check Special cases
  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(sum(n_j) == 0)
    cat("p_distinct_posterior has been called but vector n_i of previous observations is made of all zeros. Call the p_distinct_prior function instead. \n")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")
  if(r > sum(m_j) )
    return (0)
  if( r == 0 && sum(m_j) == 0 ){
    return (1) # just for coherence
  }

  # Compute non trivial cases
  return (  p_distinct_posterior_c(r,k,m_j,n_j,gamma,prior,prior_params,Max_iter)  )
}

#' p_shared_posterior
#'
#' This function computes the posterior probability that the number of shared species is equal to \code{t}.
#' @param t integer, the number of shared species whose probability has to be computed.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j a positive scalar or a vector of positive integers that defines the size of the groups that have been alread observed. It must be of the same size of \code{m_j}.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_shared_posterior = function(t, k, m_j, n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if(length(n_j)!=length(m_j))
    stop("The length of m_j must be equal to the length of n_j")
  if( any(m_j<0) || any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check Special cases
  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")
  if(sum(n_j) == 0)
    cat("p_distinct_posterior has been called but vector n_i of previous observations is made of all zeros. Call the p_distinct_prior function instead. \n")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(t<0)
    stop("Error, the number of shared species t can not be negative")
  if(t > min(m_j))
    return (0)

  # Compute non trivial cases
  return (p_shared_posterior_c(t, k, m_j, n_j, gamma, prior, prior_params, Max_iter ) )
}



#' Expected_prior
#'
#' This function computes the expected value and the variance of the a priori probability of distinct or shared species.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param type string, select \code{distinct} to compute the expected value of the number of distinct species or select \code{shared} to compute the expected value of the number of shared species.
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
Expected_prior = function(n_j, gamma, type, prior = "Poisson", ..., Max_iter = 100, tol = 10^-12){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(type != "distinct" & type != "shared")
    stop("type can only be distinct or shared")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(sum(n_j) == 0)
    stop("All values in n_j are 0.")

  if(type == "shared" & length(n_j) < 2 )
    stop("At least two groups are needed to compute the expected values of shared components")


  # Compute non trivial cases
  return( Expected_prior_c(n_j, gamma, type, prior, prior_params, Max_iter, tol  )  )
}


#' Expected_posterior
#'
#' This function computes the expected value and the variance of the a posteriori probability of distinct or shared species.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param type string, select \code{distinct} to compute the expected value of the number of distinct species or select \code{shared} to compute the expected value of the number of shared species.
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
Expected_posterior = function(k, m_j, n_j, gamma, type, prior = "Poisson", ..., Max_iter = 100, tol = 10^-12){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(type != "distinct" & type != "shared")
    stop("type can only be distinct or shared")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(sum(n_j) == 0){
    stop("All values in n_j are 0. Call the Expected_prior function instead")
  }
  if(sum(m_j) == 0)
    stop("All values in m_j are 0.")

  if(type == "shared" & length(m_j) < 2 )
    stop("At least two groups are needed to compute the expected values of shared components")

  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")

  # Compute non trivial cases
  return( Expected_posterior_c(k, m_j, n_j, gamma, type, prior, prior_params, Max_iter, tol)  )
}



#' Upper Bounds
#' @export
UpperBounds = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  UpperBounds_c(n_j,gamma,prior,prior_params,Max_iter)  )
}


#' Lower Bounds
#' @export
LowerBounds = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  LowerBounds_c(n_j,gamma,prior,prior_params,Max_iter)  )
}

#' Compute prior distribution
#' @export
D_distinct_prior = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100, Kstart = 1, logV_vec = NULL){
  l = list(...)
  L = length(l)
  n = sum(n_j)

  if(is.null(logV_vec)){
    logV_vec = rep(-Inf, n+1)
  }

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kstart<1)
    stop("The starting value for K must be >= 1")
  if(length(logV_vec)!=(n+1))
    stop("Length of logV_vec must be equal to n+1. logV_vec=NULL is a valid option.")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_distinct_prior_c(n_j,gamma,prior,prior_params,Max_iter,Kstart,logV_vec)  )
}

#' Compute prior distribution
#' @export
D_distinct_prior_interval = function( n_j, gamma, prior = "Poisson", ..., Max_iter = 100, 
                                      Kmin, Kmax, logV_vec = NULL, print = TRUE)
{
  l = list(...)
  L = length(l)
  n = sum(n_j)

  if(is.null(logV_vec)){
    logV_vec = rep(-Inf, n+1)
  }

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kmin<1)
    stop("The starting value for K must be >= 1")
  if(Kmax>n)
    stop("The maximum value for K must be <= n")
  if(length(logV_vec)!=(n+1))
    stop("Length of logV_vec must be equal to n+1. logV_vec=NULL is a valid option.")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_distinct_prior_interval_c(n_j,gamma,prior,prior_params,Max_iter,Kmin,Kmax,logV_vec,print)  )
}


#' Compute prior distribution
#' @export
D_joint_prior = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100, Kstart = 1, logV_vec = NULL){
  l = list(...)
  L = length(l)
  n = sum(n_j)

  if(is.null(logV_vec)){
    logV_vec = rep(-Inf, n+1)
  }

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kstart<1)
    stop("The starting value for K must be >= 1")
  if(length(logV_vec)!=(n+1))
    stop("Length of logV_vec must be equal to n+1. logV_vec=NULL is a valid option.")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_joint_prior_c(n_j,gamma,prior,prior_params,Max_iter,Kstart,logV_vec)  )
}

#' Compute prior distribution
#' @export
D_joint_prior_square = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100, 
                                Kmin, Kmax, Smin, Smax, 
                                logV_vec = NULL, print = TRUE)
{
  l = list(...)
  L = length(l)
  n = sum(n_j)

  if(is.null(logV_vec)){
    logV_vec = rep(-Inf, n+1)
  }

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kmin<1)
    stop("The starting value for K must be >= 1")
  if(Kmax>n)
    stop("The maximum value for K must be <= n")
  if(Smin<0)
    stop("The starting value for S must be >= 0")
  if(Smax>n)
    stop("The maximum value for S must be <= n")
  if(length(logV_vec)!=(n+1))
    stop("Length of logV_vec must be equal to n+1. logV_vec=NULL is a valid option.")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_joint_prior_square_c(n_j,gamma,prior,prior_params,Max_iter,Kmin,Kmax,Smin,Smax,logV_vec,print)  )
}

#' Compute posterior distribution
#' @export
D_joint_post_square = function(m_j, n_j, k, gamma, prior = "Poisson", ..., Max_iter = 100, 
                                Kmin, Kmax, Smin, Smax, 
                                logVpost_vec = NULL, print = TRUE)
{
  l = list(...)
  L = length(l)
  n = sum(n_j)
  m = sum(m_j)

  if(is.null(logVpost_vec)){
    logVpost_vec = rep(-Inf, m+1)
  }

  #checks
  if(length(m_j)!=length(gamma))
    stop("The length of m_j must be equal to the length of gamma")
  if(length(n_j)!=length(gamma))
      stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) || any(m_j<0))
    stop("The elements of n_j,m_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kmin<0)
    stop("The starting value for K must be >= 0")
  if(Kmax>n)
    stop("The maximum value for K must be <= n")
  if(Smin<0)
    stop("The starting value for S must be >= 0")
  if(Smax>n)
    stop("The maximum value for S must be <= n")
  if(length(logVpost_vec)!=(m+1))
    stop("Length of logVpost_vec must be equal to m+1. logVpost_vec=NULL is a valid option.")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( ! all( names(l) %in% names(prior_params) ) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_joint_post_square_c( m_j,gamma,n_j,k,prior,prior_params,Max_iter,Kmin,Kmax,Smin,Smax,logVpost_vec,print) )
}

#' arrange_partition
#'
#' This function takes as input a partition and fix it according to the notation used to define partitions in the sampler.
#' @param partition [vector] the input partition.
#' @return [vector] containing the partition following the requirement of the sampler.
#' @export
arrange_partition = function(partition){

  num = sort(unique(partition)) #get passed indices
  idx = seq(0, length(num)-1, by=1) #create sequence with correct indices

  for(h in 1:length(num)){
    if(num[h]!=idx[h])
      partition[partition == num[h]] = idx[h]
  }

  # A possible way to do this operation using apply is
  #apply(as.array(partition), MARGIN = 1, FUN = function(x){
  #h = which(num == x)
  #if(num[h]!=idx[h]) return (idx[h])
  #else return (x)
  #})
  # but it looks more expensive to me
  return (partition)
}


#' Distinct_Prior_MCMC
#'
#' @export
Distinct_Prior_MCMC = function( Niter, n_j, gamma_j, part_init = NULL,
                                prior, prior_param, M_max = 1000, seed = 0)
{
  d = length(n_j)
  n = sum(n_j)

  if(is.null(part_init))
    part_init = rep(0,n)
  if(Niter <= 0)
    stop("Niter must be positive")
  if(length(n_j)!=length(gamma_j))
    stop("length of n_j differs from length gamma_j")
  if(length(part_init) != n)
    stop("Number of points in part_init is not equal to n")

  part_init = arrange_partition(part_init)
  idx_temp  = cumsum(c(0,n_j))
  idx_start = idx_temp[1:d]+1
  idx_end   = idx_temp[2:(d+1)]

  # handle initial partition: compute N_k and N_jk
  global_clustering = table(part_init)
  K0 = length(global_clustering)
  cluster_labels = names(global_clustering)
  initial_values = data.frame( matrix(0,nrow = d, ncol = K0) )
  names(initial_values) = cluster_labels

  for(j in 1:d){
    temp = table(part_init[idx_start[j]:idx_end[j]])
    insert = unlist(lapply(as.list(cluster_labels), FUN = function(el){el %in% names(temp)}))
    initial_values[j,insert] = c(temp)
  }

  N_k = c(apply( initial_values, 2, sum))

  return( Distinct_Prior_MCMC_c( as.matrix(initial_values), N_k, part_init, K0,
                                 Niter, n_j, gamma_j,
                                 prior, prior_param,
                                 M_max, seed )
        )

}




#' Negative Binomial parametrization
#'
#' @export
NegBin_parametrization = function(mean,variance)
{
  if(variance <= 0)
    stop("Error in NegBin_parametrization: the variance must be strictly positive")

  res = list("r"=-1,"p"=-1)
  res$p = (mean)/(variance)
  res$r = (mean*mean)/(variance - mean)
  return(res)
}

#' Train-Test dataset
#'
#' @export
Ants_train_test = function(counts_long_all,d = 2,keep = 0.5,seed = 1234)
{
  # initialize long list with repeated labels
  species_long_all = lapply(1:d, function(x){c()})
  names(species_long_all) = levels(counts_long_all$site)
  for(i in 1:nrow(counts_long_all)){
    
    if(counts_long_all[i,1] > 0){
      counts = counts_long_all[i,1] # get number of repetitions
      species = counts_long_all[i,2] # get species name
      site = counts_long_all[i,3] # get area
      
      # repeat "species" for "counts" times and concatenate with past values in the same area
      species_long_all[[site]] = c(species_long_all[[site]],
                                   rep(as.character(species),counts))
    }
    
  }
  
  # define sizes of training and testind dataset
  n_j_train = ceiling(keep * n_j_all)
  m_j_test  = n_j_all - n_j_train
  
  
  # sample indexes of training sets in both populations
  set.seed(seed) # set the seed for reproducibility
  idx_train1 = sort(sample( seq_along(species_long_all[[1]]),
                            size = n_j_train[1], replace = FALSE ))
  idx_test1  = setdiff(seq_along(species_long_all[[1]]), idx_train1)
  idx_train2 = sort(sample( seq_along(species_long_all[[2]]),
                            size = n_j_train[2], replace = FALSE ))
  idx_test2  = setdiff(seq_along(species_long_all[[2]]), idx_train1)
  
  
  # initialize return objects
  data_training = counts_all;data_training[,-1] = 0
  data_test     = counts_all;data_test[,-1] = 0
  
  # compute frequecies in training and fill the return object
  x = species_long_all[[1]][idx_train1];Tx = table(x);names_species1 = names(Tx)
  y = species_long_all[[2]][idx_train1];Ty = table(y);names_species2 = names(Ty)
  for(i in 1:length(Tx)){
    data_training[1,names_species1[i]] = Tx[i]
  }
  for(i in 1:length(Ty)){
    data_training[2,names_species2[i]] = Ty[i]
  }
  
  # compute frequecies in test and fill the return object  
  x = species_long_all[[1]][idx_test1];Tx = table(x);names_species1 = names(Tx)
  y = species_long_all[[2]][idx_test2];Ty = table(y);names_species2 = names(Ty)
  for(i in 1:length(Tx)){
    data_test[1,names_species1[i]] = Tx[i]
  }
  for(i in 1:length(Ty)){
    data_test[2,names_species2[i]] = Ty[i]
  }
  
  # eliminate empty columns
  data_training = data_training[,c(TRUE,apply(data_training[,-1], 2, sum)>0)]
  data_test     = data_test[,c(TRUE,apply(data_test[,-1], 2, sum)>0)]
  
  res = list("data_training" = data_training,
             "data_test" = data_test)
  return(res)
}

