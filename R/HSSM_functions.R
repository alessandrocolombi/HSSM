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
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

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
D_distinct_prior = function(n_j, gamma, prior = "Poisson", ..., Max_iter = 100, Kexact = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(Kexact<=1)
    stop("The number of atoms whose probability is computed exactly must be strictly larger than 1");

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
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Compute non trivial cases
  return (  D_distinct_prior_c(n_j,gamma,prior,prior_params,Max_iter,Kexact)  )
}
