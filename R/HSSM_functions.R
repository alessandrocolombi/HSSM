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
D_distinct_post_interval = function( m_j, n_j, k, gamma, prior = "Poisson", ..., Max_iter = 100,
                                     Kmin, Kmax,
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
  if(Kmax>m)
    stop("The maximum value for K must be <= m")
  if(length(logVpost_vec) > (m+1))
    stop("Length of logVpost_vec must be equal to m+1 (or Kmax-Kmin+1). logVpost_vec=NULL is a valid option.")

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
  return (  D_distinct_post_interval_c(m_j,gamma,n_j,k,prior,prior_params,Max_iter,Kmin,Kmax,logVpost_vec,print)  )
}

#' Compute posterior distribution
#' @export
D_joint_post_square = function( m_j, n_j, k, gamma, prior = "Poisson", ..., Max_iter = 100,
                                Kmin, Kmax, Smin, Smax,
                                approxC = TRUE,
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
  if(Kmax>m)
    stop("The maximum value for K must be <= m")
  if(Smin<0)
    stop("The starting value for S must be >= 0")
  if(Smax>m)
    stop("The maximum value for S must be <= m")
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
  return (  D_joint_post_square_c( m_j,gamma,n_j,k,prior,prior_params,Max_iter,Kmin,Kmax,Smin,Smax,approxC,logVpost_vec,print) )
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



#' p_shared_post_largen
#'
#' @export
p_shared_post_largen = function(t, k, m_j, n_j, gamma, prior = "Poisson", ..., log_Vprior, Max_iter = 100){
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
  return (p_shared_post_largen_c(t, k, m_j, n_j, gamma, prior, prior_params, log_Vprior, Max_iter ) )
}


#' Compute posterior distribution large n
#' @export
D_joint_post_largen = function( m_j, n_j, k, gamma, prior = "Poisson", ..., Max_iter = 100,
                                Kmin, Kmax, Smin, Smax,
                                log_Vprior, print = TRUE)
{
  l = list(...)
  L = length(l)
  n = sum(n_j)
  m = sum(m_j)

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
  if(Kmax>m)
    stop("The maximum value for K must be <= m")
  if(Smin<0)
    stop("The starting value for S must be >= 0")
  if(Smax>m)
    stop("The maximum value for S must be <= m")

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
  return (  D_jointKS_post_largen_c( k,m_j,n_j,gamma,prior,prior_params,Kmin,Kmax,Smin,Smax,log_Vprior,Max_iter ) )
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
Ants_train_test = function( counts_long_all,n_j_all,counts_all,
                            d = 2,keep = 0.5,seed = 1234 )
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

### credible intervals
Find_Credible_Int = function(pmf, q = c(0.025,0.975))
{
  res = c(0,0)
  cdf = cumsum(pmf)
  stop_lower <- stop_upper <- stop <- FALSE
  i = 1
  while(!stop){
    if(cdf[i] > q[1] & stop_lower == FALSE){
      res[1] = i
      stop_lower = TRUE
    }
    if(cdf[i] > q[2] & stop_upper == FALSE){
      res[2] = i
      stop_upper = TRUE
    }
    if(stop_lower & stop_upper)
      stop = TRUE

    i = i+1
    if(i > length(cdf))
      stop = TRUE
  }
  return(res)
}




#' Moment estimator - BO
#'
#' @export
BO_MomEst = function(n_j,
                     K12,S12,K1,K2,
                     Kmin = 1,Kmax = 100,Smin = 0,Smax = 100,
                     gamma_LB , gamma_UB,
                     lambda_LB, lambda_UB,
                     BO_iterations = 50,
                     Max_iter = 100,
                     perr = 2, pesi= c(0,1/3,1/3,1/3), normalize = FALSE )
{
  library(mlrMBO)
  obj.fun <- makeSingleObjectiveFunction(
    # give a name to the objective function
    name = "Mom_Estimator",
    fn = function(x) {
      gamma_j = c(x[1],x[2])
      lambda = x[3]
      prior = "Poisson"

      ## Global quantities
      Kmax12 = min(Kmax,sum(n_j))
      Smax12 = min(Smax,sum(n_j))
      SK_joint = D_joint_prior_square(n_j = n_j, gamma = gamma_j,
                                      prior = prior, lambda = lambda,
                                      Max_iter = Max_iter,
                                      Kmin = Kmin, Kmax = Kmax12,
                                      Smin = Smin, Smax = Smax12,
                                      logV_vec = NULL, print = FALSE)
      marginals = vector("list",2)
      names(marginals) = c("K","S")
      marginals$S = apply(SK_joint,1,sum)
      marginals$K = apply(SK_joint,2,sum)

      ExpK = mean( sample(Kmin:Kmax12, size = 10000,replace = TRUE, prob = marginals$K[(Kmin:Kmax12)]) )
      ExpS = mean( sample(Smin:Kmax12, size = 10000,replace = TRUE, prob = marginals$S[(Smin:Kmax12)+1]) )

      check = sapply(marginals,sum)
      if( check[1]<0.99 || check[2]<0.99 ){
        cat("\n gamma_j = ",gamma_j,"lambda = ",lambda," \n")
        cat("NON sommano a 1")
        ExpK = 1000*K12
        ExpS = 1000*S12
      }

      ## Local quantities
      Kmax1 = min(Kmax,n_j[1])
      K_prior = D_distinct_prior_interval(n_j = n_j[1], gamma = gamma_j[1],
                                          prior = prior, lambda = lambda,
                                          Max_iter = 100,
                                          Kmin = Kmin, Kmax = Kmax1,
                                          logV_vec = NULL, print = FALSE)

      xK = sample(Kmin:Kmax1, size = 10000,replace = TRUE, prob = K_prior[(Kmin:Kmax1)])
      qK = quantile(xK, probs = c(0.025,0.975))
      ExpK1 = mean( xK )


      Kmax2 = min(Kmax,n_j[2])
      K_prior = D_distinct_prior_interval(n_j = n_j[2], gamma = gamma_j[2],
                                          prior = prior, lambda = lambda,
                                          Max_iter = 100,
                                          Kmin = Kmin, Kmax = Kmax2,
                                          logV_vec = NULL, print = FALSE)

      xK = sample(Kmin:Kmax2, size = 10000,replace = TRUE, prob = K_prior[(Kmin:Kmax2)])
      qK = quantile(xK, probs = c(0.025,0.975))
      ExpK2 = mean( xK )

      v_err = c( (ExpK-K12),
                 (ExpS-S12),
                 (ExpK1-K1),
                 (ExpK2-K2) )
      if(normalize){
        v_err[1] = v_err[1]/K12
        v_err[3] = v_err[3]/K1
        v_err[4] = v_err[4]/K2
        if(S12 > 0)
          v_err[2] = v_err[2]/S12
      }


      err = Inf
      if(perr == 1)
        err = sum( abs(pesi*v_err) )
      if(perr == 2)
        err = sum( (pesi*v_err)^2   )
      if(perr == Inf)
        err = max( abs(pesi*v_err)  )

      # err = ((ExpK-K12)^2 + (ExpS-S12)^2)/2
      # err = abs(ExpK-K12)/K12 + abs(ExpS-S12)/S12

      cat("\n ExpK = ",ExpK,"; ExpS = ",ExpS,"\n")
      cat("\n ExpK1 = ",ExpK1,"; ExpK2 = ",ExpK2,"\n")

      return( err )
    },

    # define if the objective function has to be minimized or maximized (i.e., accuracy must be maximized)
    minimize = T,

    # define the search space
    # nome - lower and upper bound
    par.set = makeParamSet(
      makeNumericParam("gamma1", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("gamma2", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("lambda", lower = lambda_LB,  upper = lambda_UB)
    )
  )


  # STEP 2: generation of the initial design
  des = generateDesign( n=5, getParamSet(obj.fun), fun=lhs::randomLHS )
  des$y = apply( des, 1, obj.fun )

  # STEP 4: Sequential process and acquisition
  control = makeMBOControl()
  control = setMBOControlTermination( control, iters=BO_iterations )

  # Run optimization
  res <- suppressWarnings(mbo(obj.fun, design=des, show.info=T , control = control))
  # print(res)
  best.seen <- getOptPathY(res$opt.path)
  result = list("gamma1" = res$x$gamma1,
                "gamma2" = res$x$gamma2,
                "lambda" = res$x$lambda,
                "error"  = res$y,
                "best.seen" = best.seen)
  return(result)

  # # retrieving and plotting results
  # best.seen <- getOptPathY(res$opt.path)
  # # best.seen <- c( max(best.seen[1:5]), best.seen[6:15] )
  # plot( cummin(best.seen), type="o", lwd=3, col="blue",
  #       ylim=c(0,10),
  #       ylab="best seen", xlab="trials")
  # lines( best.seen, type="o", lty=2, col="green", lwd=3 )
  # legend( "topright", legend=c("best seen","path"), col=c("blue","green"), lty=1:2, lwd=3, pch=1 )
}

#' Moment estimator - BO
#'
#' @export
BO_MomEst_multiple = function(n_j_list,v_K12,v_S12,
                              Kmin = 1,Kmax = 100,
                              Smin = 0,Smax = 100,
                              gamma_LB , gamma_UB,
                              lambda_LB, lambda_UB,
                              BO_iterations = 100,
                              Max_iter = 100,
                              perr = 2, c1 = 0.5, c2 = 0.5)
{
  L = length(n_j_list)
  if(length(v_K12)!=L)
    stop("length of v_K12 is not the same as n_j_list")
  if(length(v_S12)!=L)
    stop("length of v_S12 is not the same as n_j_list")
  if(perr != 2 & perr != 1 & perr != Inf)
    stop("perr can only be 1,2,Inf")
  if((c1+c2) != 1)
    stop("c1 and c2 must sum up to 1")

  library(mlrMBO)
  obj.fun <- makeSingleObjectiveFunction(
    # give a name to the objective function
    name = "Mom_Estimator",
    fn = function(x) {
      gamma_j = c(x[1],x[2])
      lambda = x[3]
      prior = "Poisson"
      v_errK <- v_errS <- v_ExpK <- v_ExpS <- rep(-1,L)

      for(it in 1:L){
        n_j_it = n_j_list[[it]]
        K12_it = v_K12[it]
        S12_it = v_S12[it]
        Kmin_it = min(Kmin,sum(n_j_it))
        Smin_it = min(Smin,sum(n_j_it))
        Kmax_it = min(Kmax,sum(n_j_it))
        Smax_it = min(Smax,sum(n_j_it))
        SK_joint = D_joint_prior_square(n_j = n_j_it, gamma = gamma_j,
                                        prior = prior, lambda = lambda,
                                        Max_iter = Max_iter,
                                        Kmin = Kmin_it, Kmax = Kmax_it,
                                        Smin = Smin_it, Smax = Smax_it,
                                        logV_vec = NULL, print = FALSE)

        marginals = vector("list",2)
        names(marginals) = c("K","S")
        marginals$S = apply(SK_joint,1,sum)
        marginals$K = apply(SK_joint,2,sum)

        check = sapply(marginals,sum)
        if( check[1]<0.99 || check[2]<0.99 ){
          cat("\n gamma_j = ",gamma_j,"lambda = ",lambda," \n")
          cat("NON sommano a 1")
          v_ExpK[it] = 1000*K12_it
          v_ExpS[it] = 1000*S12_it
        }else{
          v_ExpK[it] = mean( sample(Kmin_it:Kmax_it, size = 10000,replace = TRUE, prob = marginals$K[(Kmin_it:Kmax_it)]) )
          v_ExpS[it] = mean( sample(Smin_it:Smax_it, size = 10000,replace = TRUE, prob = marginals$S[(Smin_it:Smax_it)+1]) )
        }
        v_errK[it] = (v_ExpK[it]-K12_it)
        v_errS[it] = (v_ExpS[it]-S12_it)

        if(normalize){
          v_errK[it] = v_errK[it]/K12_it
          if(S12_it > 0)
            v_errS[it] = v_errS[it]/S12_it
        }
        # v_err[it] = abs(v_ExpK[it]-K12_it)/K12_it + abs(v_ExpS[it]-S12_it)/S12_it
      }
      cat("\n v_ExpK = ")
      print(v_ExpK)
      cat("\n v_ExpS = ")
      print(v_ExpS)

      err = Inf
      if(p == 1)
        err = sum( c1*abs(v_errK) + c2*abs(v_errS) )
      if(p == 2)
        err = sum( c1*(v_errK)^2 + c2*(v_errS)^2   )
      if(p == Inf)
        err = max( c1*abs(v_errK) + c2*abs(v_errS)   )
      return( err )
    },

    # define if the objective function has to be minimized or maximized (i.e., accuracy must be maximized)
    minimize = T,

    # define the search space
    # nome - lower and upper bound
    par.set = makeParamSet(
      makeNumericParam("gamma1", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("gamma2", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("lambda", lower = lambda_LB,  upper = lambda_UB)
    )
  )


  # STEP 2: generation of the initial design
  des = generateDesign( n=5, getParamSet(obj.fun), fun=lhs::randomLHS )
  des$y = apply( des, 1, obj.fun )

  # STEP 4: Sequential process and acquisition
  control = makeMBOControl()
  control = setMBOControlTermination( control, iters=BO_iterations )

  # Run optimization
  res <- suppressWarnings(mbo(obj.fun, design=des, show.info=T , control = control))
  # print(res)
  best.seen <- getOptPathY(res$opt.path)
  result = list("gamma1" = res$x$gamma1,
                "gamma2" = res$x$gamma2,
                "lambda" = res$x$lambda,
                "error"  = res$y,
                "best.seen" = best.seen)
  return(result)
  # # retrieving and plotting results

  # plot( cummin(best.seen), type="o", lwd=3, col="blue",
  #       ylim=c(0,10),
  #       ylab="best seen", xlab="trials")
  # lines( best.seen, type="o", lty=2, col="green", lwd=3 )
  # legend( "topright", legend=c("best seen","path"), col=c("blue","green"), lty=1:2, lwd=3, pch=1 )
}


#' Moment estimator - BO
#'
#' @export
BO_MomEst_NB = function(n_j,
                         K12,S12,K1,K2,
                         Kmin = 1,Kmax = 100,Smin = 0,Smax = 100,
                         gamma_LB , gamma_UB,
                         mu0_LB, mu0_UB, r_LB, r_UB,
                         BO_iterations = 50,
                         Max_iter = 100,
                         perr = 2, pesi= c(0.25,0.25,0.25,0.25) )
{
  library(mlrMBO)
  obj.fun <- makeSingleObjectiveFunction(
    # give a name to the objective function
    name = "Mom_Estimator",
    fn = function(x) {
      gamma_j = c(x[1],x[2])
      mu0 = x[3]
      r   = x[4]
      p = r/(mu0+r)
      prior = "NegativeBinomial"

      ## Global quantities
      SK_joint = D_joint_prior_square(n_j = n_j, gamma = gamma_j,
                                      prior = prior, p = p, r = r,
                                      Max_iter = Max_iter,
                                      Kmin = Kmin, Kmax = Kmax,
                                      Smin = Smin, Smax = Smax,
                                      logV_vec = NULL, print = FALSE)

      marginals = vector("list",2)
      names(marginals) = c("K","S")
      marginals$S = apply(SK_joint,1,sum)
      marginals$K = apply(SK_joint,2,sum)

      ExpK = mean( sample(Kmin:Kmax, size = 10000,replace = TRUE, prob = marginals$K[(Kmin:Kmax)]) )
      ExpS = mean( sample(Smin:Smax, size = 10000,replace = TRUE, prob = marginals$S[(Smin:Smax)+1]) )

      check = sapply(marginals,sum)
      if( check[1]<0.99 || check[2]<0.99 ){
        cat("\n gamma_j = ",gamma_j,"; mu0 = ",mu0,"; r = ",r," \n")
        cat("NON sommano a 1")
        ExpK = 1000*K12
        ExpS = 1000*S12
      }

      ## Local quantities
      K_prior = D_distinct_prior_interval(n_j = n_j[1], gamma = gamma_j[1],
                                          prior = prior, p = p, r = r,
                                          Max_iter = 100,
                                          Kmin = Kmin, Kmax = Kmax,
                                          logV_vec = NULL, print = FALSE)


      xK = sample(Kmin:Kmax, size = 10000,replace = TRUE, prob = K_prior[(Kmin:Kmax)])
      qK = quantile(xK, probs = c(0.025,0.975))
      ExpK1 = mean( xK )
      K_prior = D_distinct_prior_interval(n_j = n_j[2], gamma = gamma_j[2],
                                          prior = prior, p = p, r = r,
                                          Max_iter = 100,
                                          Kmin = Kmin, Kmax = Kmax,
                                          logV_vec = NULL, print = FALSE)


      xK = sample(Kmin:Kmax, size = 10000,replace = TRUE, prob = K_prior[(Kmin:Kmax)])
      qK = quantile(xK, probs = c(0.025,0.975))
      ExpK2 = mean( xK )

      v_err = c( (ExpK-K12),
                 (ExpS-S12),
                 (ExpK1-K1),
                 (ExpK2-K2) )
      if(normalize){
        v_err[1] = v_err[1]/K12
        v_err[3] = v_err[3]/K1
        v_err[4] = v_err[4]/K2
        if(S12 > 0)
          v_err[2] = v_err[2]/S12
      }

      err = Inf
      if(perr == 1)
        err = sum( abs(pesi*v_err) )
      if(perr == 2)
        err = sum( (pesi*v_err)^2   )
      if(perr == Inf)
        err = max( abs(pesi*v_err)  )

      # err = ((ExpK-K12)^2 + (ExpS-S12)^2)/2
      # err = abs(ExpK-K12)/K12 + abs(ExpS-S12)/S12

      cat("\n ExpK = ",ExpK,"; ExpS = ",ExpS,"\n")
      cat("\n ExpK1 = ",ExpK1,"; ExpK2 = ",ExpK2,"\n")
      return( err )
    },

    # define if the objective function has to be minimized or maximized (i.e., accuracy must be maximized)
    minimize = T,

    # define the search space
    # nome - lower and upper bound
    par.set = makeParamSet(
      makeNumericParam("gamma1", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("gamma2", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("mu0",    lower = mu0_LB,   upper = mu0_UB),
      makeNumericParam("r",      lower = r_LB,     upper = r_UB)
    )
  )


  # STEP 2: generation of the initial design
  des = generateDesign( n=5, getParamSet(obj.fun), fun=lhs::randomLHS )
  des$y = apply( des, 1, obj.fun )

  # STEP 4: Sequential process and acquisition
  control = makeMBOControl()
  control = setMBOControlTermination( control, iters=BO_iterations )

  # Run optimization
  res <- suppressWarnings(mbo(obj.fun, design=des, show.info=T , control = control))
  # print(res)
  best.seen <- getOptPathY(res$opt.path)
  result = list("gamma1" = res$x$gamma1,
                "gamma2" = res$x$gamma2,
                "mu0" = res$x$mu0,
                "r" = res$x$r,
                "error"  = res$y,
                "best.seen" = best.seen)
  return(result)
}

#' Moment estimator - BO
#'
#' @export
BO_MomEst_multiple_NB = function(n_j_list,v_K12,v_S12,
                              Kmin = 1,Kmax = 100,
                              Smin = 0,Smax = 100,
                              gamma_LB , gamma_UB,
                              mu0_LB, mu0_UB, r_LB, r_UB,
                              BO_iterations = 100,
                              Max_iter = 100,
                              perr = 2, c1 = 0.5, c2 = 0.5)
{
  L = length(n_j_list)
  if(length(v_K12)!=L)
    stop("length of v_K12 is not the same as n_j_list")
  if(length(v_S12)!=L)
    stop("length of v_S12 is not the same as n_j_list")
  if(perr != 2 & perr != 1 & perr != Inf)
    stop("p can only be 1,2,Inf")
  if((c1+c2) != 1)
    stop("c1 and c2 must sum up to 1")

  library(mlrMBO)
  obj.fun <- makeSingleObjectiveFunction(
    # give a name to the objective function
    name = "Mom_Estimator",
    fn = function(x) {
      gamma_j = c(x[1],x[2])
      mu0 = x[3]
      r   = x[4]
      p = r/(mu0+r)
      prior = "NegativeBinomial"
      v_errK <- v_errS <- v_ExpK <- v_ExpS <- rep(-1,L)

      for(it in 1:L){
        n_j_it = n_j_list[[it]]
        K12_it = v_K12[it]
        S12_it = v_S12[it]
        Kmin_it = min(Kmin,sum(n_j_it))
        Smin_it = min(Smin,sum(n_j_it))
        Kmax_it = min(Kmax,sum(n_j_it))
        Smax_it = min(Smax,sum(n_j_it))
        SK_joint = D_joint_prior_square(n_j = n_j_it, gamma = gamma_j,
                                        prior = prior, p = p, r = r,
                                        Max_iter = Max_iter,
                                        Kmin = Kmin_it, Kmax = Kmax_it,
                                        Smin = Smin_it, Smax = Smax_it,
                                        logV_vec = NULL, print = FALSE)

        marginals = vector("list",2)
        names(marginals) = c("K","S")
        marginals$S = apply(SK_joint,1,sum)
        marginals$K = apply(SK_joint,2,sum)

        check = sapply(marginals,sum)
        if( check[1]<0.99 || check[2]<0.99 ){
          cat("\n gamma_j = ",gamma_j,"; mu0 = ",mu0,"; r = ",r," \n")
          cat("NON sommano a 1")
          v_ExpK[it] = 1000*K12_it
          v_ExpS[it] = 1000*S12_it
        }else{
          v_ExpK[it] = mean( sample(Kmin_it:Kmax_it, size = 10000,replace = TRUE, prob = marginals$K[(Kmin_it:Kmax_it)]) )
          v_ExpS[it] = mean( sample(Smin_it:Smax_it, size = 10000,replace = TRUE, prob = marginals$S[(Smin_it:Smax_it)+1]) )
        }
        v_errK[it] = (v_ExpK[it]-K12_it)
        v_errS[it] = (v_ExpS[it]-S12_it)
        # v_err[it] = abs(v_ExpK[it]-K12_it)/K12_it + abs(v_ExpS[it]-S12_it)/S12_it
        if(normalize){
          v_errK[it] = v_errK[it]/K12_it
          if(S12_it > 0)
            v_errS[it] = v_errS[it]/S12_it
        }
      }
      cat("\n v_ExpK = ")
      print(v_ExpK)
      cat("\n v_ExpS = ")
      print(v_ExpS)

      err = Inf
      if(perr == 1)
        err = sum( c1*abs(v_errK) + c2*abs(v_errS) )
      if(perr == 2)
        err = sum( c1*(v_errK)^2 + c2*(v_errS)^2   )
      if(perr == Inf)
        err = max( c1*abs(v_errK) + c2*abs(v_errS)   )
      return( err )
    },

    # define if the objective function has to be minimized or maximized (i.e., accuracy must be maximized)
    minimize = T,

    # define the search space
    # nome - lower and upper bound
    par.set = makeParamSet(
      makeNumericParam("gamma1", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("gamma2", lower = gamma_LB, upper = gamma_UB),
      makeNumericParam("mu0",    lower = mu0_LB,   upper = mu0_UB),
      makeNumericParam("r",      lower = r_LB,     upper = r_UB)
    )
  )


  # STEP 2: generation of the initial design
  des = generateDesign( n=5, getParamSet(obj.fun), fun=lhs::randomLHS )
  des$y = apply( des, 1, obj.fun )

  # STEP 4: Sequential process and acquisition
  control = makeMBOControl()
  control = setMBOControlTermination( control, iters=BO_iterations )

  # Run optimization
  res <- suppressWarnings(mbo(obj.fun, design=des, show.info=T , control = control))
  # print(res)
  best.seen <- getOptPathY(res$opt.path)
  result = list("gamma1" = res$x$gamma1,
                "gamma2" = res$x$gamma2,
                "mu0" = res$x$mu0,
                "r" = res$x$r,
                "error"  = res$y,
                "best.seen" = best.seen)
  return(result)
  # # retrieving and plotting results

  # plot( cummin(best.seen), type="o", lwd=3, col="blue",
  #       ylim=c(0,10),
  #       ylab="best seen", xlab="trials")
  # lines( best.seen, type="o", lty=2, col="green", lwd=3 )
  # legend( "topright", legend=c("best seen","path"), col=c("blue","green"), lty=1:2, lwd=3, pch=1 )
}


#' Rarefaction curves
#'
#' @export
Rarefaction_curves = function(data, Nsort = 2, seed0 = 220424 ){

  add = function(x){Reduce("+", x)} # used to sum vectors/matrices stored in elements of a list
  data = data[apply(data,1,sum)>0,]
  d = ncol(data)
  r = nrow(data)
  n_j  = apply(data,2,sum)
  n    = sum(n_j)


  species_long_all = lapply(1:d, function(x){c()})
  names(species_long_all) = c("A1","A2")

  for(i in 1:r){
    for(j in 1:d){
      if(data[i,j] > 0){
        counts = data[i,j] # get number of repetitions
        species = i # get species name
        site = names(species_long_all)[j] # get area

        # repeat "species" for "counts" times and concatenate with past values in the same area
        species_long_all[[site]] = c(species_long_all[[site]],
                                     rep(as.character(species),counts))
      }
    }
  }

  set.seed(seed0)
  seeds = sample(1:99999,size = Nsort,replace = TRUE)
  K12obs_n <- K1obs_n <- K2obs_n <- S12obs_n <- matrix(0,n,Nsort)#lapply(1:Nsort,function(x){c()})

  for(it in 1:Nsort){
    cat("\n it = ",it,"\n")
    set.seed(seeds[it])
    new_idx1 = sample(1:n_j[1],size = n_j[1], replace = F)
    new_idx2 = sample(1:n_j[2],size = n_j[2], replace = F)
    new_idx12 = sample(1:n,size = n, replace = F)
    X1_reordered = species_long_all[[1]][new_idx1]
    X2_reordered = species_long_all[[2]][new_idx2]
    temp = tibble(species = c(X1_reordered,X2_reordered),
                  site    = c(rep(names(species_long_all)[1],n_j[1]),
                              rep(names(species_long_all)[2],n_j[2]))
    )
    temp_reordered = temp[new_idx12,]
    species_ordered_n = tibble(species = as.character(1:r),
                               counts_A1 = rep(0,r),
                               counts_A2 = rep(0,r) )

    pb = txtProgressBar(min = 0, max = n, initial = 0)
    for(i in 1:n){
      temp = temp_reordered[i,]
      j = -1
      if(temp$site == names(species_long_all)[1]){
        j = 1
      }else if(temp$site == names(species_long_all)[2]){
        j = 2
      }else{
        stop("ERRORE")
      }

      species_ordered_n[which(species_ordered_n$species == temp$species),j+1] = species_ordered_n[which(species_ordered_n$species == temp$species),j+1] + 1

      K12obs_n[i,it]  = nrow(species_ordered_n[apply(species_ordered_n[,-1],1,sum)>0,])
      K1obs_n[i,it]   = length(which(species_ordered_n[,2]>0))
      K2obs_n[i,it]   = length(which(species_ordered_n[,3]>0))
      #K12obs_n[[it]] = c(K12obs_n[[it]],nrow(species_ordered_n[apply(species_ordered_n[,-1],1,sum)>0,]))
      #K1obs_n[[it]]  = c(K1obs_n[[it]],length(which(species_ordered_n[,2]>0)))
      #K2obs_n[[it]]  = c(K2obs_n[[it]],length(which(species_ordered_n[,3]>0)))
      setTxtProgressBar(pb,i)
    }
    close(pb)
    S12obs_n[,it]  = K1obs_n[,it] + K2obs_n[,it] - K12obs_n[,it]
    # S12obs_n[[it]] = K1obs_n[[it]] + K2obs_n[[it]] - K12obs_n[[it]]
  }


  K12obs_summary = apply( K12obs_n,1,quantile, probs = c(0.025,0.5,0.975) )
  S12obs_summary = apply( S12obs_n,1,quantile, probs = c(0.025,0.5,0.975) )
  K12obs_mean    = apply( K12obs_n,1,mean )
  S12obs_mean    = apply( S12obs_n,1,mean )
  # K12obs_mean = add(K12obs_n)/Nsort
  # K1obs_mean = add(K1obs_n)/Nsort
  # K2obs_mean = add(K2obs_n)/Nsort
  # S12obs_mean = add(S12obs_n)/Nsort

  ## plot

  ## return
  res = list("K12obs_mean"  = K12obs_mean,
             "S12obs_mean" = S12obs_mean,
             "K12obs_summary" = K12obs_summary,
             "S12obs_summary" = S12obs_summary )
  return(res)
}




#' Rarefaction curves (d=1)
#'
#' @export
Rarefaction_curve_d1 = function(data, Nsort = 50, seed0 = 220424 ){

  add = function(x){Reduce("+", x)} # used to sum vectors/matrices stored in elements of a list
  data[which(data>0)]
  d = 1
  r = length(data)
  n = sum(data)


  species_long_all = c()

  for(i in 1:r){
    counts = data[i] # get number of repetitions
    species = i # get species name

    # repeat "species" for "counts" times and concatenate with past values in the same area
    species_long_all = c(species_long_all,
                         rep(as.character(species),counts))
  }

  set.seed(seed0)
  seeds = sample(1:99999,size = Nsort,replace = TRUE)
  Kobs_n = matrix(0,n,Nsort) #lapply(1:Nsort,function(x){c()})

  for(it in 1:Nsort){
    cat("\n it = ",it,"\n")
    set.seed(seeds[it])
    new_idx = sample(1:n,size = n, replace = F)
    X_reordered = species_long_all[new_idx]
    species_ordered_n = tibble(species = as.character(1:r),
                               counts = rep(0,r) )

    pb = txtProgressBar(min = 0, max = n, initial = 0)
    for(i in 1:n){
      x = X_reordered[i]
      species_ordered_n[which(species_ordered_n$species == x),2] = species_ordered_n[which(species_ordered_n$species == x),2] + 1

      Kobs_n[i,it]  = length(which(species_ordered_n[,2]>0))
      # Kobs_n[[it]]  = c(Kobs_n[[it]],length(which(species_ordered_n[,2]>0)))
      setTxtProgressBar(pb,i)
    }
    close(pb)
  }

  ## return
  Kobs_summary = apply( Kobs_n,1,quantile, probs = c(0.025,0.5,0.975) )
  Kobs_mean    = apply( Kobs_n,1,mean )
  res = list( "Kobs_mean"  = Kobs_mean,
              "Kobs_summary" = Kobs_summary )
  # Kobs_mean = add(Kobs_n)/Nsort
  return(res)
}









#' Training and Testing
#'
#' @export
Train_Test = function(data, keep = 0.5 , seed = 220424 ){

  data = data[apply(data,1,sum)>0,]
  d = ncol(data)
  r = nrow(data)
  n_j  = apply(data,2,sum)
  n    = sum(n_j)


  species_long_all = lapply(1:d, function(x){c()})
  names(species_long_all) = c("A1","A2")

  for(i in 1:r){
    for(j in 1:d){
      if(data[i,j] > 0){
        counts = data[i,j] # get number of repetitions
        species = i # get species name
        site = names(species_long_all)[j] # get area

        # repeat "species" for "counts" times and concatenate with past values in the same area
        species_long_all[[site]] = c(species_long_all[[site]],
                                     rep(as.character(species),counts))
      }
    }
  }

  set.seed(seed)
  # K12obs_n <- K1obs_n <- K2obs_n <- S12obs_n <- matrix(0,n,Nsort)#lapply(1:Nsort,function(x){c()})


  new_idx1 = sample(1:n_j[1],size = n_j[1], replace = F)
  new_idx2 = sample(1:n_j[2],size = n_j[2], replace = F)
  new_idx12 = sample(1:n,size = n, replace = F)
  X1_reordered = species_long_all[[1]][new_idx1]
  X2_reordered = species_long_all[[2]][new_idx2]
  temp = tibble(species = c(X1_reordered,X2_reordered),
                site    = c(rep(names(species_long_all)[1],n_j[1]),
                            rep(names(species_long_all)[2],n_j[2]))
  )
  temp_reordered = temp[new_idx12,]

  n_train = floor(keep * n)
  res = list("data_training" = temp_reordered[1:n_train,],
             "data_test" = temp_reordered[(n_train+1):(n),])

  # Sintetizzo il train set
  species = unique(res$data_training$species)
  res_tr = tibble("site" = c("A1","A2"))
  for(s in species){

    nA1 = nrow(res$data_training[which(res$data_training$species == s & res$data_training$site == "A1"),])
    nA2 = nrow(res$data_training[which(res$data_training$species == s & res$data_training$site == "A2"),])
    res_tr = res_tr %>% cbind( tibble(col = c(nA1,nA2)) )
    names(res_tr)[ncol(res_tr)] = s

  }
  # Sintetizzo il test set
  species = unique(res$data_test$species)
  res_te = tibble("site" = c("A1","A2"))
  for(s in species){

    nA1 = nrow(res$data_test[which(res$data_test$species == s & res$data_test$site == "A1"),])
    nA2 = nrow(res$data_test[which(res$data_test$species == s & res$data_test$site == "A2"),])
    res_te = res_te %>% cbind( tibble(col = c(nA1,nA2)) )
    names(res_te)[ncol(res_te)] = s
  }

  res = list("data_training" = res_tr,
             "data_test" = res_te)
  return(res)
}
