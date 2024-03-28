// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>

#include "species.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------

// These functions computes log(sum_i(a_i)) using a stable formula for log values. Let us denote a* to the the maximum value of vector a which is attained when i = i*.
// ---> log(sum_i(a_i)) = log(a*) + log[ 1 + sum_{i not i*}(exp{log(a_i) - log(a*)}) ]
// See that only the logarithm of the elements of a are needed. Hence, it is likely that one has already computed them in log scale. If so, set is_log = T
// In this versione of the function, the maximum and its position are passed as parameters. No check is performed to be sure that such information were true or not.
// It is assumed that the value in val_max is on the same scale of the values of a, i.e it is in log scale if is_log is set to true.
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max, const unsigned int& idx_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	if(is_log==TRUE){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
				std::log(1 +
					    std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
					    std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
				        )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;
		// Do not checks if values of a are strictly positive
		return ( std::log(val_max) +
				 std::log(1 +
					      std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )  +
					      std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )
			             )
			   );
	}
}

// As before but gets an iterator poining to the maximum value
double log_stable_sum(const std::vector<double>& a, const bool is_log, const HSSM_Traits::vector_d_citerator& it_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;
	double val_max{*it_max};
	if(is_log){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (    val_max +
					std::log(1 +
						    std::accumulate(   a.cbegin(), it_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
						    std::accumulate(   it_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
				            )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		// Do not checks if values of a are strictly positive
		return ( std::log(val_max) +
				 std::log(1 +
					      std::accumulate(   a.cbegin(), it_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )  +
					      std::accumulate(   it_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )
			             )
			   );
	}
}

// In this specialized version of the function, the position of the max value is not provided. Hence, one additional operation is done.
// Since exp( log(a*) - log(a*)) = 1, the previous formula becomes
// ---> log(sum_i(a_i)) = log(a*) + log[ sum_{i}(exp{log(a_i) - log(a*)}) ]
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	// Do not checks if it is really the max value
	if(is_log){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
					std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )   )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		return ( std::log(val_max) +
				 std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   ) )
			   );
	}

}

// In this version of the formula, the maximum value is computed
double log_stable_sum(const std::vector<double>& a, const bool is_log){
	if(a.size() == 0)
		return 0.0;

	// Computes maximum value
	HSSM_Traits::vector_d_citerator it_max{std::max_element(a.cbegin(), a.cend())};
	//double val_max{*it_max};
	// Calls the specialized version
	return log_stable_sum(a,is_log,it_max);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Factorials and Pochammer
//------------------------------------------------------------------------------------------------------------------------------------------------------

// NOTE: this is only for testing the best implementation. This function is not used.
double raising_factorial_poch(const unsigned int& n, const double& a)
{
	return( gsl_sf_poch(a, (double)n ) );
}

// NOTE: this is only for testing the best implementation. This function is not used.
double log_raising_factorial_poch(const unsigned int& n, const double& a) //secondo me troppo generale, può essere semplificata
{
	if(a<=0)
		throw std::runtime_error("Error in log_raising_factorial, can not compute the raising factorial of a negative number in log scale");

	return( gsl_sf_lnpoch(a, (double)n ) );
}

// This is used in the code
double log_raising_factorial(const unsigned int& n, const double& a)
{

	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	else{

		double val_max{std::log(a+n-1)};
		double res{1.0};
		if (n==1)
			return val_max;
		for(std::size_t i = 0; i <= n-2; ++i){
			res += std::log(a + (double)i) / val_max;
		}
		return val_max*res;

	}
}

// This is not used in the code but it is useful to be exported
double raising_factorial(const unsigned int& n, const double& a)
{
	if(n==0)
		return 1.0;
	if(n==1)
		return a;
	if(a<=0){
		double res{1.0};
		for(std::size_t i = 0; i <= n-1; ++i){
			res *= ( a + (double)i ) ;
		}
		return res;
	}
	else{
		return std::exp(log_raising_factorial(n,a));
	}

}


double my_falling_factorial_old(const unsigned int& n, const double& a)
{
	if(n%2 == 0) //n is even
		return( gsl_sf_poch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_poch(-a, (double)n ) );
}

double my_log_falling_factorial_old(const unsigned int& n, const double& a) //questo non va per a negativi!
{
	if(n%2 == 0) //n is even
		return( gsl_sf_lnpoch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_lnpoch(-a, (double)n ) );
}

// This is used in the code
double my_log_falling_factorial(const unsigned int& n, const double& a)
{
	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_falling_factorial, can not compute the falling factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	if(a-n+1<=0)
		throw std::runtime_error("Error in my_log_falling_factorial, can not compute the falling factorial (a)_n in log scale if a <= n-1");
	else{
		double val_max{std::log(a)};
		double res{1.0};
		if (n==1)
			return val_max;
		for(std::size_t i = 1; i <= n-1; ++i){
			res += std::log(a - (double)i) / val_max;
		}
		return val_max*res;
	}
}

// This is not used in the code but it is useful to be exported
double my_falling_factorial(const unsigned int& n, const double& a)
{
	if(n==0)
		return 1.0;
	if(n==1)
		return a;
	if(a<=0){
		double res{1.0};
		for(std::size_t i = 0; i <= n-1; ++i){
			res *= ( a + (double)i ) ;
		}
		return res;
	}
	else{
		return std::exp(my_log_falling_factorial(n,a));
	}
}

// This is not used in the code but it is useful to be exported
double compute_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_poch(a,x ) );
}

// This is not used in the code but it is useful to be exported
double compute_log_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_lnpoch(a,x ) );
}


double log_zeta_Riemann(double s){
	return( std::log(std::riemann_zetal(s)) );
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	C numbers
//------------------------------------------------------------------------------------------------------------------------------------------------------

// NOTE: this is only for testing the best implementation. This function is not used.
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res)
{

	//Rcpp::Rcout<<"Costruisco n = "<<n<<", k = "<<k<<std::endl;
	if(k>n)
		throw std::runtime_error("Can not call build_logC_matrix if k > n");

	// Boundary conditions
	if(n == 0 && k == 0)
		res(n,k) = 0;
	else if(k == 0)
		res(n,k) = std::log(raising_factorial(n,r));
	else if(n == k){
		build_logC_matrix(n-1, k-1, s, r, res);
		res(n,k) = std::log(s) + res(n-1,k-1);
	}
	else{
		double coef(s*k + r + n - 1);
		build_logC_matrix(n-1,k-1,s,r,res); // the matrix is always constructed diagonal by diagonal
		res(n,k) = std::log(coef) + res(n-1,k) + std::log( 1 + s/coef * std::exp( res(n-1,k-1) - res(n-1,k) ) );
	}

}

// NOTE: this is only for testing the best implementation. This function is not used.
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
my_logC(const unsigned int& n, const double& scale, const double& location)
{

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");

	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	MatCol res(MatCol::Constant(n+1,n+1,0.0)); //nxn matrix, n is in standard notation, i.e it starts from 1
	for(int k = n; k>=0; k--){
		//Rcpp::Rcout<<"Nel ciclo for, n = "<<n<<" , k = "<<k<<std::endl;
		build_logC_matrix(n,k,s,r,res);
		//Rcpp::Rcout<<"res finito k = "<<std::endl<<res<<std::endl;
	}
	return(res);
}

// NOTE: this is only for testing the best implementation. This function is not used.
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location)
{

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");

	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative


	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = std::log(raising_factorial(1,r));
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update(0) = std::log(raising_factorial(nn,r));
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	}
	return (LogC_old);
}

// NOTE: this is only for testing the best implementation. This function is not used.
Eigen::VectorXd my_logC2_central(const unsigned int& n, const double& scale)
{

	if(!( scale<0 ) )
		throw std::runtime_error("Error in my_logC2. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative.");

	const double& s = -scale; //s is strictly positive
	const double inf = std::numeric_limits<double>::infinity();

	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = -inf;
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	LogC_update(0) = -inf;
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	}
	return (LogC_old);
}


// This is used in the code
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");

	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative


	Rcpp::NumericVector LogC_old(n+1, 0.0);

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old[0] = log_raising_factorial(1,r);
	LogC_old[1] = std::log(s);

	//Rcpp::NumericVector LogC_update(LogC_old);
	Rcpp::NumericVector LogC_update(n+1, 0.0);
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update[0] = log_raising_factorial(nn,r);
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update[k] = std::log(coef) + LogC_old[k] + std::log( 1 + s/coef*std::exp( LogC_old[k-1] - LogC_old[k] ) );
		}
		LogC_update[nn] = nn*std::log(s); //update last element

		//std::copy(LogC_update.begin(),LogC_update.end(),LogC_old.begin());
		std::swap(LogC_old, LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	}
	return (LogC_old);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Approximations for V
//------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> log_Vprior_long(const unsigned int& k, const std::vector<unsigned int>& n_i,
									const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param,
									unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Initialize vector of results
	std::vector<double> log_vect_res(M_max+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
		// Formula implementation

		log_vect_res[Mstar] = log_raising_factorial(k,Mstar+1 ) +
							  qM.log_eval_prob(Mstar + k) -
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   			);

		// Check if it is the new maximum
        if(log_vect_res[Mstar]>val_max){
        	idx_max = Mstar;
        	val_max = log_vect_res[Mstar];
        }

	}
	// Formula to compute the log of all the sums in a stable way
	return log_vect_res;
}

int log_Vprior_apprx1(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
									  const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param,
									  unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Component prior preliminary operations
	// auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	// ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Constant term
	int n = std::accumulate(n_i.cbegin(),n_i.cend(),0.0);
	if(n-k < 2)
		throw std::runtime_error("Error in log_Vprior_apprx1: this approximation only works if K < n-1 ");
	double log_constant = 0.5*( std::log( std::exp(2.0)*k ) - std::log( 2.0*M_PI ) );
	log_constant -= std::inner_product( 	n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   	[](const unsigned int& nj, const double& gamma_j){
				       					   		return (nj * std::log(gamma_j) );
				       					   	}
			       					  );

	Rcpp::Rcout<<"log_constant:"<<std::endl<<log_constant<<std::endl;
	Rcpp::Rcout<<"std::log(tol)-log_constant:"<<std::endl<<std::log(tol)-log_constant<<std::endl;

	double log_z = std::log(std::riemann_zetal((double)(n-k)));

	unsigned int N{1};
	std::vector<double> resto;
	resto.reserve(M_max);
	double log_res{0.0};
	while(N<=M_max){
		Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
		resto.push_back(-std::log((double)N)*(double)(n-k));
		double log_resto = log_stable_sum(resto,TRUE);

		Rcpp::Rcout<<"log_z:"<<std::endl<<log_z<<std::endl;
		Rcpp::Rcout<<"log_resto:"<<std::endl<<log_resto<<std::endl;
		Rcpp::Rcout<<"diff:"<<std::endl<<std::exp(log_z)-std::exp(log_resto)<<std::endl;

		log_res = log_z + std::log(1-std::exp(log_resto-log_z) );
		if(log_res < std::log(tol)-log_constant){
			break;
		}
		N++;
	}

	// Formula to compute the log of all the sums in a stable way
	return N;
}

int log_Vprior_apprx2(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Constant term
	int n = std::accumulate(n_i.cbegin(),n_i.cend(),0.0);

	double log_constant = 0.5*( std::log( std::exp(2.0)*k ) - std::log( 2.0*M_PI ) );
	log_constant -= std::inner_product( 	n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   	[](const unsigned int& nj, const double& gamma_j){
				       					   		return (nj * std::log(gamma_j) );
				       					   	}
			       					  );

	unsigned int N{qM.get_mode()};
	while(N<=M_max){

		double log_res = -(double)(n-k)*log(N+k) + qM.log_eval_prob(N+k);
		if(log_res < std::log(tol)-log_constant){
			break;
		}
		N++;
	}

	// Formula to compute the log of all the sums in a stable way
	return N;
}

int log_Vprior_apprx2(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param,
						unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;
	return log_Vprior_apprx2(k,n_i,tol,gamma,qM,M_max );
}


int log_Vprior_apprx3(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Constant term
	int n = std::accumulate(n_i.cbegin(),n_i.cend(),0.0);
	double log_constant = 0.5*( std::log( std::exp(2.0)*k ) - std::log( 2.0*M_PI ) );
	log_constant -= std::inner_product( 	n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   	[](const unsigned int& nj, const double& gamma_j){
				       					   		return (nj * std::log(gamma_j) );
				       					   	}
			       					  );

	//Rcpp::Rcout<<"log_constant:"<<std::endl<<log_constant<<std::endl;
	//Rcpp::Rcout<<"std::log(tol)-log_constant:"<<std::endl<<std::log(tol)-log_constant<<std::endl;

	unsigned int N{1};
	while(N<=M_max){

		//Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
		double log_res = -(double)(n-k)*log(N+k) + std::log(qM.eval_upper_tail(N+k+1));
		if(log_res < std::log(tol)-log_constant){
			break;
		}
		N++;
	}

	// Formula to compute the log of all the sums in a stable way
	return N;
}


int log_Vprior_apprx3(const unsigned int& k, const std::vector<unsigned int>& n_i, const double& tol,
						const std::vector<double>& gamma, const Rcpp::String& prior, const Rcpp::List& prior_param,
						unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;
	return log_Vprior_apprx3(k,n_i,tol,gamma,qM,M_max) ;

}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A priori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------


// Da Testare (ma tanto è inutile, calcolo sempre solo il log)
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma,
					  const ComponentPrior& qM, unsigned int M_max )
{

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Vprior, the length of n_i (group sizes) and gamma has to be equal");
	double res(0.0);
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		res += raising_factorial(k, (double)(Mstar+1) ) * qM.eval_prob(Mstar + k) *
		       std::inner_product( n_i.cbegin(), n_i.cend(),gamma.cbegin(), 1.0, std::multiplies<>(),
		       					   [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return 1/raising_factorial( nj, gamma_j*((double)Mstar + (double)k));} );
		        // nj is an integer, this is just a raising factorial, not a pochammer
	}
	return res;
}


// Works for d>2
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma,
	                      const ComponentPrior& qM, unsigned int M_max )
{

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");

	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	// Compute the number of usefull iterations to be done
	int Max_M = 0; //log_Vprior_apprx3(k,n_i,1e-20,gamma,qM,M_max);
	Max_M = M_max;


	//Rcpp::Rcout<<"M_max:"<<std::endl<<M_max<<std::endl;
	//Rcpp::Rcout<<"Max_M:"<<std::endl<<Max_M<<std::endl;
	//Rcpp::Rcout<<"qM.get_mode():"<<std::endl<<qM.get_mode()<<std::endl;


	// Initialize vector of results
	std::vector<double> log_vect_res(Max_M+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= Max_M; ++Mstar){
		//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
		// Formula implementation

		log_vect_res[Mstar] = log_raising_factorial(k,Mstar+1 ) +
							  qM.log_eval_prob(Mstar + k) -
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   			);

		// Check if it is the new maximum
        if(log_vect_res[Mstar]>val_max){
        	idx_max = Mstar;
        	val_max = log_vect_res[Mstar];
        }

	}
	// Formula to compute the log of all the sums in a stable way
	return log_stable_sum(log_vect_res, TRUE, val_max, idx_max);
}


//questa è sola per 1 o 2 gruppi
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma)
{

	//Rcpp::Rcout<<"INPUT --> compute_Kprior_unnormalized: d = "<<n_i.size()<<"; K = "<<k<<std::endl;
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;

	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) must be equal to 1 or 2");

	if(k == 0){
		//Need to handle degenerate cases. The probability of k=0 is 1 when n_j = 0. It is needed to have coherence when some elements of n_i are equal to 0.
		if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			return 0.0;
		}
		else
			return -inf;
	}
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;



	if(n_i.size()==1){ // one group only
		Rcpp::NumericVector absC = compute_logC(n_i[0], -gamma[0], 0.0); //absC[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		return absC[k];
	}
	else{ // two groups case

		// Compute all C numbers required
		Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2

		const unsigned int start1 = std::max( 0, (int)k - (int)n_i[0] ); //max(0, k-n1)
		const unsigned int start2 = std::max( 0, (int)k - (int)n_i[1] ); //max(0, k-n2)
		const unsigned int end1   = std::min( (int)k, (int)n_i[1] );     //min(k,n2)
					//Rcpp::Rcout<<"start1 = "<<start1<<std::endl;
					//Rcpp::Rcout<<"start2 = "<<start2<<std::endl;
					//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;

					//std::vector<double> log_a(k+1, -inf);    // old version
		std::vector<double> log_a(end1-start1+1, -inf);    // This vector contains all the quantities that depend only on r1
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop
		unsigned int outer_indx{0};
		for(std::size_t r1=start1; r1 <= end1; ++r1){

			// Compute a_r1 using its definition
							//log_a[outer_indx] = gsl_sf_lnchoose(k,r1) - my_log_falling_factorial(r1,(double)k) + absC1[k-r1];  //<---- old version
			log_a[outer_indx] =  absC1[k-r1];

			// Prepare for computing the second term

			// Initialize vector of results
					//std::vector<double> log_vect_res(k-r1+1, -inf ); // old version
			std::vector<double> log_vect_res(k-r1-start2+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);

			// Inner loop on r2
			unsigned int inner_indx{0};
			for(std::size_t r2=start2; r2<= k-r1; ++r2){

				// Compute b_r2*c_r1r2
							//log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r1,r2) + std::lgamma(k-r2+1) +  absC2[k-r2];   //<---- old version
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r2,r1) + my_log_falling_factorial(k-r1-r2,(double)(k-r1)) +  absC2[k-r2];

				// Check if it is the new maximum of log_vect_res
	        	if(log_vect_res[inner_indx]>val_max2){
	        		idx_max2 = inner_indx;
	        		val_max2 = log_vect_res[inner_indx];
	        	}
	        			//Rcpp::Rcout<<"Computing for r1 = "<<r1<<" and r2 = "<<r2<<std::endl;
				inner_indx++;
			}


			// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
			log_a[outer_indx] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);

			// Check if it is the new maximum of log_a
	       	if(log_a[outer_indx]>val_max1){
	       		idx_max1 = outer_indx;
	       		val_max1 = log_a[outer_indx];
	       	}
	       	outer_indx++;
		}

		// Complete the sum over all elements in log_a
		//Rcpp::Rcout<<"OUTPUT --> compute_Kprior_unnormalized: K = "<<k<<"; d = "<<n_i.size()<<"; Prob = "<<log_stable_sum(log_a, TRUE, val_max1, idx_max1)<<std::endl;
		return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	}
}

//questa è sola per 1 o 2 gruppi e in più vengono passati i numeri C già calcolati
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i,
									const std::vector<double>& gamma,
									const Rcpp::NumericVector& absC1, const Rcpp::NumericVector& absC2)
{

	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) must be equal to 1 or 2");

	if(k == 0){
		//Need to handle degenerate cases. The probability of k=0 is 1 when n_j = 0. It is needed to have coherence when some elements of n_i are equal to 0.
		if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			return 0.0;
		}
		else
			return -inf;
	}
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;


	if(n_i.size()==1){ // one group only
		if(absC1.size() != (n_i[0]+1))
			throw std::runtime_error("Error in compute_Kprior_unnormalized: the length of absC1 is not compatible with n_1");
		return absC1[k];
	}
	else{ // two groups case

		// Check length of C numbers
		if(absC1.size() != (n_i[0]+1))
			throw std::runtime_error("Error in compute_Kprior_unnormalized: the length of absC1 is not compatible with n_1");
		if(absC2.size() != (n_i[1]+1))
			throw std::runtime_error("Error in compute_Kprior_unnormalized: the length of absC2 is not compatible with n_2");

		const unsigned int start1 = std::max( 0, (int)k - (int)n_i[0] ); //max(0, k-n1)
		const unsigned int start2 = std::max( 0, (int)k - (int)n_i[1] ); //max(0, k-n2)
		const unsigned int end1   = std::min( (int)k, (int)n_i[1] );     //min(k,n2)

		std::vector<double> log_a(end1-start1+1, -inf);    // This vector contains all the quantities that depend only on r1

		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop
		unsigned int outer_indx{0};
		for(std::size_t r1=start1; r1 <= end1; ++r1){

			// Compute a_r1 using its definition
			log_a[outer_indx] =  absC1[k-r1];

			// Prepare for computing the second term
			// Initialize vector of results
			std::vector<double> log_vect_res(k-r1-start2+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);

			// Inner loop on r2
			unsigned int inner_indx{0};
			for(std::size_t r2=start2; r2<= k-r1; ++r2){

				// Compute b_r2*c_r1r2
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r2,r1) + my_log_falling_factorial(k-r1-r2,(double)(k-r1)) +  absC2[k-r2];

				// Check if it is the new maximum of log_vect_res
	        	if(log_vect_res[inner_indx]>val_max2){
	        		idx_max2 = inner_indx;
	        		val_max2 = log_vect_res[inner_indx];
	        	}
				inner_indx++;
			}


			// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
			log_a[outer_indx] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);

			// Check if it is the new maximum of log_a
	       	if(log_a[outer_indx]>val_max1){
	       		idx_max1 = outer_indx;
	       		val_max1 = log_a[outer_indx];
	       	}
	       	outer_indx++;
		}

		// Complete the sum over all elements in log_a
		return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	}
}

//questa va bene per qualunque d, se d<=2 chiama compute_Kprior_unnormalized
double compute_Kprior_unnormalized_recursive(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma)
{

	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");
	if( n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior_unnormalized_recursive, the length of n_i (group sizes) must be greater than 0");
	if( n_i.size() <= 2)
		return compute_Kprior_unnormalized(k,n_i,gamma);
	if(k == 0)
		return -inf;
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		return -inf;

	//Rcpp::Rcout<<"INPUT --> compute_Kprior_unnormalized_recursive: d = "<<n_i.size()<<"; K = "<<k<<std::endl;
	//Rcpp::Rcout<<"Stampo n_i: ";
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;
	// if here, n_i.size()>2
	std::vector<double> log_a(k+1, -inf);
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);

	// Start for loop
	for(std::size_t k1=0; k1 <= k; ++k1){

				//Rcpp::Rcout<<"----> Dentro a k1 = "<<k1<<std::endl;
		// Compute recursive formula
		log_a[k1] = compute_Kprior_unnormalized_recursive(k1, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} );
				//Rcpp::Rcout<<"log_a["<<k1<<"] = "<<"NNP(K^(2) = "<<k1<<") = "<<log_a[k1]<<std::endl;
		// Prepare for computing the second term

		// Initialize vector of results
		std::vector<double> log_vect_res(k1+1, -inf );

		// Initialize quantities to find the maximum of log_vect_res
		unsigned int idx_max2(0);
		double val_max2(log_vect_res[idx_max2]);

		// Inner loop on r2
		unsigned int inner_indx{0};
		for(std::size_t k2=k-k1; k2<= k; ++k2){
					//Rcpp::Rcout<<"***** Dentro a k2 = "<<k2<<std::endl;

			// Compute last element
					//Rcpp::Rcout<<"n_i[n_i.size()-1] = "<<n_i[n_i.size()-1]<<std::endl;
					//Rcpp::Rcout<<"gamma[gamma.size()-1] = "<<gamma[gamma.size()-1]<<std::endl;
					//Rcpp::Rcout<<"NNP(K^(1) = "<<k2<<") = "<< compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  ) <<std::endl;
					//Rcpp::Rcout<<"coef = "<<std::exp( gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) )<<std::endl;
					//Rcpp::Rcout<<"coef giusto = "<< ( gsl_sf_fact(k1)*gsl_sf_fact(k2)  )/(gsl_sf_fact(k-k1)*gsl_sf_fact(k-k2)*gsl_sf_fact(k1+k2-k)  ) <<std::endl;
			log_vect_res[inner_indx] = gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) + compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );

		 	// Check if it is the new maximum of log_vect_res
	       	if(log_vect_res[inner_indx]>val_max2){
	       		idx_max2 = inner_indx;
	       		val_max2 = log_vect_res[inner_indx];
	       	}
	   		inner_indx++;
		}
		// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
					//Rcpp::Rcout<<"Vettore interno"<<std::endl;
					//for(auto __v : log_vect_res)
						//Rcpp::Rcout<<__v<<", ";
					//Rcpp::Rcout<<std::endl;
					//Rcpp::Rcout<<"Adding effect: log_stable_sum(log_vect_res) = "<<log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2)<<std::endl;
		log_a[k1] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
					//Rcpp::Rcout<<"Updated log_a["<<k1<<"] = "<<log_a[k1]<<std::endl;
		// Check if it is the new maximum of log_a
	   	if(log_a[k1]>val_max1){
	    	idx_max1 = k1;
	      	val_max1 = log_a[k1];
	    }
	}

			//Rcpp::Rcout<<"Vettore finale"<<std::endl;
			//for(auto __v : log_a)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
	// Complete the sum over all elements in log_a
	//Rcpp::Rcout<<"OUTPUT --> compute_Kprior_unnormalized_recursive: K = "<<k<<"; d = "<<n_i.size()<<"; Prob = "<<log_stable_sum(log_a, TRUE, val_max1, idx_max1)<<std::endl;
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);
}

//Direct formula per d=2
double compute_SK_prior_unnormalized(const unsigned int& k, const unsigned int& s,
									 const std::vector<unsigned int>& n_i, const std::vector<double>& gamma)
{

			//Rcpp::Rcout<<"Dentro compute_SK_prior_unnormalized"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) must be equal to 1 or 2");

	// d=1
	if(n_i.size() == 1)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) can not be zero");
	// k=0 case
	if(k==0){
		if( s==0 && *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			throw std::runtime_error("Error in compute_SK_prior_unnormalized, k and s can not be zero even if all n_j are zero");
			return 0.0;
		}
		else
			return -inf;
	}
	// k<s case, if here k>0
	if(k<s)
		return -inf;
	// k>n1+n2 case
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;
	// s > min(n1,n2) case
	if( s > *std::min_element(n_i.cbegin(),n_i.cend()) )
		return -inf;

	//If here, need to compute the formula
	// Compute all C numbers required
	Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2

	//Compute and check the range of r1
	const unsigned int start1 = std::max( 0, (int)k - (int)n_i[1] ); //max(0,k-n2)
	const unsigned int end1   = std::min( (int)(k-s), (int)n_i[0] - (int)s ); //min(k-s,n1-s)
	const int range_size = (int)end1-(int)start1+1;
				//Rcpp::Rcout<<"start1 = "<<start1<<std::endl;
				//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;
				//Rcpp::Rcout<<"Range size = "<<range_size<<std::endl;

	if(range_size <= 0) //All Cnumbers would be 0
		return -inf;

	std::vector<double> log_a(range_size, -inf);    // This vector contains all the quantities that depend only on r1
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);
	// Start for loop
	unsigned int outer_indx{0};
	for(std::size_t r1=start1; r1 <= end1; ++r1){
		// Compute a_r1 using its definition
		log_a[outer_indx] = gsl_sf_lnchoose(k-r1,s) + my_log_falling_factorial(s,(double)(s+r1)) + absC1[s+r1] + absC2[k-r1] ;

		// Check if it is the new maximum of log_a
	       if(log_a[outer_indx]>val_max1){
	       	idx_max1 = outer_indx;
	       	val_max1 = log_a[outer_indx];
	       }
	       outer_indx++;
	}
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	//return -1.0;
}

//Direct formula per d=2
double compute_SK_prior_unnormalized(const unsigned int& k, const unsigned int& s,
									 const std::vector<unsigned int>& n_i, const std::vector<double>& gamma,
									 const Rcpp::NumericVector& absC1, const Rcpp::NumericVector& absC2 )
{

			//Rcpp::Rcout<<"Dentro compute_SK_prior_unnormalized"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) must be equal to 1 or 2");

	// d=1
	if(n_i.size() == 1)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) can not be zero");
	// k=0 case
	if(k==0){
		if( s==0 && *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			throw std::runtime_error("Error in compute_SK_prior_unnormalized, k and s can not be zero even if all n_j are zero");
			return 0.0;
		}
		else
			return -inf;
	}
	// k<s case, if here k>0
	if(k<s)
		return -inf;
	// k>n1+n2 case
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;
	// s > min(n1,n2) case
	if( s > *std::min_element(n_i.cbegin(),n_i.cend()) )
		return -inf;

	//If here, need to compute the formula
	// C numbers are passed as input...

	//Compute and check the range of r1
	const unsigned int start1 = std::max( 0, (int)k - (int)n_i[1] ); //max(0,k-n2)
	const unsigned int end1   = std::min( (int)(k-s), (int)n_i[0] - (int)s ); //min(k-s,n1-s)
	const int range_size = (int)end1-(int)start1+1;
				//Rcpp::Rcout<<"start1 = "<<start1<<std::endl;
				//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;
				//Rcpp::Rcout<<"Range size = "<<range_size<<std::endl;

	if(range_size <= 0) //All Cnumbers would be 0
		return -inf;

	std::vector<double> log_a(range_size, -inf);    // This vector contains all the quantities that depend only on r1
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);
	// Start for loop
	unsigned int outer_indx{0};
	for(std::size_t r1=start1; r1 <= end1; ++r1){
		// Compute a_r1 using its definition
		log_a[outer_indx] = gsl_sf_lnchoose(k-r1,s) + my_log_falling_factorial(s,(double)(s+r1)) + absC1[s+r1] + absC2[k-r1] ;

		// Check if it is the new maximum of log_a
	       if(log_a[outer_indx]>val_max1){
	       	idx_max1 = outer_indx;
	       	val_max1 = log_a[outer_indx];
	       }
	       outer_indx++;
	}
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	//return -1.0;
}

//Recursive formula for d>2
double compute_SK_prior_unnormalized_recursive(const unsigned int& k, const unsigned int& s, const std::vector<unsigned int>& n_i,
											   const std::vector<double>& gamma)
{

			//Rcpp::Rcout<<"Dentro compute_SK_prior_unnormalized_recursive"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) must be equal to 1 or 2");
	// d=1
	if(n_i.size() == 1){
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) can not be zero");
		//Se finisco qua, una possibilità è chiamare compute_Kprior_unnormalized
	}
	// d=2 case
	if(n_i.size() == 2)
		return compute_SK_prior_unnormalized(k,s,n_i,gamma);
	// k=0 case
	if(k==0){
		if( s==0 && *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){
			Rcpp::Rcout<<"Sono in compute_SK_prior_unnormalized_recursive e sono finito nel caso P(K=0,S=0) con tutti i termini di n_i nulli"<<std::endl;
			return 0.0;
		}
		else
			return -inf;
	}
	// k<s case, if here k>0
	if(k<s)
		return -inf;
	// k>n1+...+nd case
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;
	// s > min(n1,...,nd) case
	if( s > *std::min_element(n_i.cbegin(),n_i.cend()) )
		return -inf;


	// if here, n_i.size()>2 in a non trivial case

	// Initialize quantities for sum on s1
	std::vector<double> log_s1(k-s+1, -inf);
	// Initialize quantities to find the maximum of log_s1
	unsigned int idx_max0(0);
	double val_max0(log_s1[idx_max0]);

	// Start loop on s1
	for(std::size_t s1 = 0; s1<=k-s; ++s1){
				//Rcpp::Rcout<<"^^^^ Dentro s1 = "<<s1<<std::endl;
		// Compute coeff that depends only on s1
		log_s1[s1] = my_log_falling_factorial(s,(double)(s+s1));

		// Initialize quantities for inner sums
		const unsigned int k_s_s1 = k-s-s1;
		// Compute the log of the inner sum on h and t
		std::vector<double> log_a(k_s_s1+1, -inf);
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop on h
		for(std::size_t h=0; h <= k_s_s1; ++h){

					//Rcpp::Rcout<<"----> Dentro a h = "<<h<<std::endl;
			// Compute recursive formula
			log_a[h] = compute_SK_prior_unnormalized_recursive(k-h, s+s1, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} );
					//Rcpp::Rcout<<"log_a["<<h<<"] = "<<"NNP(S^(2) = "<<s+s1<<", K^(2) = "<<k-h<<") = "<<log_a[h]<<std::endl;


			// Prepare for computing the second term
			const unsigned int k_s_s1_h = k-s-s1-h;
			// Initialize vector of results
			std::vector<double> log_vect_res(k_s_s1_h+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);

			// Inner loop on t
			unsigned int inner_indx{0};
			for(std::size_t t=0; t<= k_s_s1_h; ++t){
						//Rcpp::Rcout<<"***** Dentro a t = "<<t<<std::endl;

				// Compute last element
							//Rcpp::Rcout<<"********************"<<std::endl;
							//Rcpp::Rcout<<"Coef = "<<std::exp(my_log_falling_factorial(s,(double)(s+s1)) + gsl_sf_lnchoose(k_s_s1_h,t) + gsl_sf_lnchoose(s+h+t,s) + my_log_falling_factorial(t,(double)(h+t)))<<std::endl;
							//Rcpp::Rcout<<"Coef giusto = "<< ( gsl_sf_fact(s+s1)*gsl_sf_fact(k-s-s1-h)*gsl_sf_fact(s+h+t) )/(gsl_sf_fact(s)*gsl_sf_fact(s1)*gsl_sf_fact(h)*gsl_sf_fact(t)*gsl_sf_fact(k-s-s1-h-t)  ) <<std::endl;
							//Rcpp::Rcout<<"--------------------"<<std::endl;
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k_s_s1_h,t) + gsl_sf_lnchoose(s+h+t,s) + my_log_falling_factorial(t,(double)(h+t)) +
										   compute_Kprior_unnormalized( s+h+t, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );

			 	// Check if it is the new maximum of log_vect_res
		       	if(log_vect_res[inner_indx]>val_max2){
		       		idx_max2 = inner_indx;
		       		val_max2 = log_vect_res[inner_indx];
		       	}
		   		inner_indx++;
			}
			// Update outer vector
			log_a[h] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
						//Rcpp::Rcout<<"Updated log_a["<<h<<"] = "<<log_a[h]<<std::endl;
			// Check if it is the new maximum of log_a
		   	if(log_a[h]>val_max1){
		    	idx_max1 = h;
		      	val_max1 = log_a[h];
		    }
		}

		// End of loop on h. Update the sum on s1
		log_s1[s1] += log_stable_sum(log_a, TRUE, val_max1, idx_max1);

		// Check if it is the new maximum of log_s1
	   	if(log_s1[s1]>val_max0){
	    	idx_max0 = s1;
	      	val_max0 = log_s1[s1];
	    }

	}


	// Complete the sum over all elements in log_s1
	return log_stable_sum(log_s1, TRUE, val_max0, idx_max0);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A posteriori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------


Rcpp::NumericVector log_qM_post(const unsigned int& m,
								const ComponentPrior& qM,
								const unsigned int& k, const std::vector<unsigned int>& n_j,
								const std::vector<double>& gamma_j, double log_V, unsigned int M_max )
{
	// Basic quantities
	const unsigned int n{std::accumulate(n_j.cbegin(),n_j.cend(), 0)};
	const unsigned int d{gamma_j.size()};
	const double inf{std::numeric_limits<double>::infinity()};
	// Checks
	if(m < 0)
		throw std::runtime_error("Error in log_qM_post: m must be equal or greater than 0 ");
	if(k > n)
		throw std::runtime_error("Error in log_qM_post. It is not possible that k is higher than n.");
	if(n_j.size() != d)
		throw std::runtime_error("Error in log_qM_post. n_j does not have the same size of gamma_j.");

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );


	double logqMpost =	log_raising_factorial(k,(double)m+1.0) +
				  	  	qM.log_eval_prob(m+k) -
				      	std::inner_product( n_j.cbegin(),n_j.cend(),gamma_j.cbegin(), 0.0, std::plus<>(),
			       					[&m, &k](const unsigned int& nj, const double& gammaj){
			       						return log_raising_factorial( nj, gammaj*(m + k) );
			       					}
			       		);
	logqMpost -= log_V;
	Rcpp::NumericVector res = {logqMpost,log_V};
	return res;
}

Rcpp::NumericVector log_qM_post(const unsigned int& m,
								const Rcpp::String& prior, const Rcpp::List& prior_param,
								const unsigned int& k, const std::vector<unsigned int>& n_j,
								const std::vector<double>& gamma_j, double log_V, unsigned int M_max )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	return ( log_qM_post( m, qM, k, n_j, gamma_j, log_V, M_max ) );
}


std::vector<double> D_log_qM_post(	const ComponentPrior& qM,
									const unsigned int& k, const std::vector<unsigned int>& n_j,
									const std::vector<double>& gamma_j, double log_V, unsigned int M_max )
{
	double inf = std::numeric_limits<double>::infinity();

	unsigned int it_max{10000}; // set maximum number of evaluations for q_M post
	std::vector<double> log_res; // initialize return vector. elements are in log scale
	log_res.reserve(it_max);
	double max{-inf}; // initialize max value in log_res
	unsigned int idx_max{0}; // initialize position of max element in log_res
	bool repeat = TRUE;		 // flag to continue the while loop
	unsigned int mstar{0};   // initialize mstar
	unsigned int counter{0}; // initialize position of max element in log_res

	while( repeat ){
		// compute log_qM_post
		Rcpp::NumericVector temp = 	log_qM_post(mstar, qM, k, n_j, gamma_j, log_V, M_max );

		if(log_V == -inf){
			// if here, log_V was not provided and must be computed
			log_V = temp[1];
		}

		log_res.push_back(temp[0]); // store result

		// check for max
		if(log_res[mstar] > max){
			max = log_res[mstar];
			idx_max = mstar;
		}

		// compute cumulative probability. Note: this is not in log scale
		double cum_P = std::exp(log_stable_sum(log_res,TRUE,max,idx_max));

		// check if a sufficient amount of mass has been computed
		if( 1.0-cum_P < 1e-10)
			repeat = FALSE;
		// check if maximum number of iterations has been reached
		if(counter >= it_max){
			repeat = FALSE;
			Rcpp::Rcout<<"WARNING: v_log_qM_post reached "<<it_max<<" evaluations but it does not sum up to 1."<<std::endl;
		}

		counter++; // update counter
	}

	return log_res;
}

// NON USARLA, è BUGGATA
std::vector<double> build_log_qM_post(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma,
									  const ComponentPrior& qM, unsigned int M_max )
{
	// checks !!
	if(n_i.size() == 0)
		throw std::runtime_error("Error in build_log_qM_post, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in build_log_qM_post, the length of n_i (group sizes) and gamma has to be equal");


	if(k == 0){
		// The behaviour of such cases is indefined beacuse i am not sure if the need to be handled for coherence or not
		if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			throw std::runtime_error("Error in build_log_qM_post. The behaviuor of this function when k=0 and n_i is exactly zero is indefined.");
		}
		else
			throw std::runtime_error("Error in build_log_qM_post. It does not make any sense to evaluate this function when k=0 and n_i is not exactly zero. The behaviuor is indefined");
	}

	// Compute all terms from 0 to M_max+k
	const unsigned int end{M_max+k};
	// Initialize vector of results
	std::vector<double> log_vect_res(end+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= end; ++Mstar){

		// Formula implementation
		log_vect_res[Mstar] = log_raising_factorial(k,Mstar+1 ) +
							  qM.log_eval_prob(Mstar + k) -
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   			);

		// Check if it is the new maximum
	    if(log_vect_res[Mstar]>val_max){
	       	idx_max = Mstar;
	       	val_max = log_vect_res[Mstar];
	    }

	}

	// Formula to compute the log of all the sums in a stable way
	double log_norm_const{log_stable_sum(log_vect_res, TRUE, val_max, idx_max)};
			//Rcpp::Rcout<<"Stampo log_vect_res: ";
			//for(auto __v : log_vect_res)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
		//
			Rcpp::Rcout<<"log_norm_const:"<<std::endl<<log_norm_const<<std::endl;

	// Normalize all terms
	std::transform(log_vect_res.begin(), log_vect_res.end(), log_vect_res.begin(), [&log_norm_const](double& x){return x - log_norm_const;}); // MA NON HA MODIFICATO NIENTE!!

			Rcpp::Rcout<<"POST NORMALIZZAZIONE"<<std::endl;
				Rcpp::Rcout<<"Stampo log_vect_res: ";
			for(auto __v : log_vect_res)
				Rcpp::Rcout<<__v<<", ";
			Rcpp::Rcout<<std::endl;
	//return
	return log_vect_res;
}

double compute_log_Vpost(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i,
						 const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) and gamma has to be equal");

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){
		Rcpp::Rcout<<"The compute_log_Vpost has been called but vector n_i of previous observations is made of all zeros. Call the compute_log_Vprior function instead."<<std::endl;
		return compute_log_Vprior(r, m_i, gamma, qM, M_max );
	}
	if(k == 0){
		throw std::runtime_error("Error in compute_log_Vpost. If here, n_i is not exactly zero but k is set to 0. This makes no sense, such behaviuor is indefined");
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vpost. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	//if(r == 0) --> r may be 0 in posterior calculations
	const unsigned int start = r;
	const unsigned int end   = M_max;

	// Initialize vector of results and for log_Vprior
	std::vector<double> log_Vprior_vec(end+1, -std::numeric_limits<double>::infinity() );
	std::vector<double> log_vect_res(end-start+1, -std::numeric_limits<double>::infinity() );

	// Initialize quantities to find the maximum
	unsigned int idx_max_V{0};
	double val_max_V(log_Vprior_vec[idx_max_V]);

	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Can not use inner_product with 3 inputs. Need to compute in advance the vector of the arguments of log_raising_factorial
	std::vector<double> fact_argument(n_i.size(), 0.0);
	// Start the loop, let us compute all the elements
	unsigned int indx{0};
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){

		//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
		double coef_qM = log_raising_factorial(k,Mstar+1 ) +
					  	 qM.log_eval_prob(Mstar + k) -
						 std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   		 [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   		);
		// Assign value to compute log_Vprior
		log_Vprior_vec[Mstar] = coef_qM;
		// Check if it is maximum for log_Vprior
		if(log_Vprior_vec[Mstar]>val_max_V){
			idx_max_V = Mstar;
			val_max_V = log_Vprior_vec[Mstar];
		}

		// Formula implementation
		if(Mstar >= r){
			// Fill the elements of fact_argument
			std::transform(n_i.cbegin(),n_i.cend(),gamma.cbegin(),fact_argument.begin(),
						   [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return ( (double)nj + (double)(Mstar+k)*gamma_j );} );
			log_vect_res[indx] = my_log_falling_factorial(r, (double)Mstar ) +
								 coef_qM -
								 std::inner_product( m_i.cbegin(),m_i.cend(), fact_argument.cbegin(), 0.0, std::plus<>(),
								 				      [](const unsigned int& mj, const double& xj){return log_raising_factorial( mj, xj );}
								 				   );

			// Check if it is the new maximum
			if(log_vect_res[indx]>val_max){
			    idx_max = indx;
			    val_max = log_vect_res[indx];
			}

			indx++;

		}
	}

			//Rcpp::Rcout<<"Stampo log_vect_res: ";
			//for(auto __v : log_vect_res)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
			//Rcpp::Rcout<<"log_stable_sum(log_vect_res, TRUE, val_max, idx_max):"<<std::endl<<log_stable_sum(log_vect_res, TRUE, val_max, idx_max)<<std::endl;
			//Rcpp::Rcout<<"log_Vprior:"<<std::endl<<log_stable_sum(log_Vprior_vec, TRUE, val_max_V, idx_max_V)<<std::endl;

	// Formula to compute the log of all the sums in a stable way
	return (log_stable_sum(log_vect_res, TRUE, val_max, idx_max) - log_stable_sum(log_Vprior_vec, TRUE, val_max_V, idx_max_V) );
}

double compute_log_Vpost_naive(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i,
						 	   const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{
	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) and gamma has to be equal");

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){
		Rcpp::Rcout<<"The compute_log_Vpost has been called but vector n_i of previous observations is made of all zeros. Call the compute_log_Vprior function instead."<<std::endl;
		return compute_log_Vprior(r, m_i, gamma, qM, M_max );
	}
	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vpost. If here, n_i is not exactly zero but k is set to 0. This makes no sense, such behaviuor is indefined");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		throw std::runtime_error("Error in compute_log_Vpost. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	//if(r == 0) --> r may be 0 in posterior calculations

	double log_Vprior{compute_log_Vprior(k, n_i, gamma, qM, 1000)};
		// Initialize vector of results
	const unsigned int start = r;
	const unsigned int end   = M_max;

		std::vector<double> log_vect_res(end-start+1, -std::numeric_limits<double>::infinity() );
		// Initialize quantities to find the maximum
		unsigned int idx_max{0};
		double val_max(log_vect_res[idx_max]);

		// Start the loop, let us compute all the elements
		std::vector<double> fact_argument(n_i.size(), 0.0);
		unsigned int indx{0};
		for(std::size_t Mstar=start; Mstar <= end; ++Mstar){
			//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;

			// Formula implementation
			// Fill the elements of fact_argument
			std::transform(n_i.cbegin(),n_i.cend(),gamma.cbegin(),fact_argument.begin(),
						   [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return ( (double)nj + (double)(Mstar+k)*gamma_j );} );
			//Rcpp::Rcout<<"qM.log_eval_prob("<<Mstar + k<<") - log_Vprior = "<<std::endl<<qM.log_eval_prob(Mstar + k) - log_Vprior <<std::endl;
			log_vect_res[indx] =  gsl_sf_lnfact(Mstar+k) - gsl_sf_lnfact(Mstar-r) +
								  qM.log_eval_prob(Mstar + k) -
								  std::inner_product( m_i.cbegin(),m_i.cend(),fact_argument.cbegin(), 0.0, std::plus<>(),
								  				      [](const unsigned int& mj, const double& xj){return log_raising_factorial( mj, xj );}
								  				    ) -
								  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
				       					   			);
			// Check if it is the new maximum
	        if(log_vect_res[indx]>val_max){
	        	idx_max = indx;
	        	val_max = log_vect_res[indx];
	        }
	        indx++;
		}
		// Formula to compute the log of all the sums in a stable way
				//Rcpp::Rcout<<"Stampo log_vect_res: ";
				//for(auto __v : log_vect_res)
					//Rcpp::Rcout<<__v<<", ";
				//Rcpp::Rcout<<std::endl;
				//Rcpp::Rcout<<"log_stable_sum(log_vect_res, TRUE, val_max, idx_max):"<<std::endl<<log_stable_sum(log_vect_res, TRUE, val_max, idx_max)<<std::endl;
				//Rcpp::Rcout<<"log_Vprior:"<<std::endl<<log_Vprior<<std::endl;
		return ( log_stable_sum(log_vect_res, TRUE, val_max, idx_max) - log_Vprior );
}

//questa è sola per 1 o 2 gruppi
double compute_Kpost_unnormalized(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i,
						 		  const std::vector<double>& gamma)
{

	//Rcpp::Rcout<<"Dentro a compute_Kpost_unnormalized: ";
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;

	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (group sizes) and gamma has to be equal");


	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //check it is a prior case
		Rcpp::Rcout<<"The compute_Kpost_unnormalized has been called but vector n_i of previous observations is made of all zeros. Call the compute_Kprior_unnormalized function instead"<<std::endl;
		return compute_Kprior_unnormalized(r, m_i, gamma );
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) // check that the known data are coherent
		throw std::runtime_error("Error in compute_Kpost_unnormalized. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	if( r > std::accumulate(m_i.cbegin(), m_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;

	if( *std::max_element(m_i.cbegin(),m_i.cend()) == 0 ){
		    if( r == 0 )
		    	return 0.0;
		    else
		    	return -inf; //This case is handle just for completeness and to sake of clarity. However, the case m1=m2=0 and r>0 is handled in the previous if (r > m1+m2)
	}

	// if here, it is not possible that n1=n2=0.
	if( (n_i.size()==1) || ( n_i.size()==2 && m_i[1] == 0 ) ){ // one group only or m2=0
		Rcpp::NumericVector absC = compute_logC(  m_i[0], -gamma[0], - ( (double)k*gamma[0] + (double)n_i[0] )  ); //absC[i] = |C(m1,i,-gamma1,-(k*gamma1 + n1))| for i = 0,...,m1
		return absC[r];
	}
	else if( n_i.size()==2 && m_i[0] == 0 ){ // m1=0
		Rcpp::NumericVector absC = compute_logC(  m_i[1], -gamma[1], - ( (double)k*gamma[1] + (double)n_i[1] )  ); //absC[i] = |C(m2,i,-gamma2,-(k*gamma2 + n2))| for i = 0,...,m2
		return absC[r];
	}
	else{ // two groups case.

		// Compute all C numbers required
		Rcpp::NumericVector absC1 = compute_logC(  m_i[0], -gamma[0], - ( (double)k*gamma[0] + (double)n_i[0] )  ); //absC[i] = |C(m1,i,-gamma1,-(k*gamma1 + n1))| for i = 0,...,m1
		Rcpp::NumericVector absC2 = compute_logC(  m_i[1], -gamma[1], - ( (double)k*gamma[1] + (double)n_i[1] )  ); //absC[i] = |C(m2,i,-gamma2,-(k*gamma2 + n2))| for i = 0,...,m2

		const unsigned int start1 = std::max( 0, (int)r - (int)m_i[0] ); //max(0, r-m1)
		const unsigned int start2 = std::max( 0, (int)r - (int)m_i[1] ); //max(0, r-m2)
		const unsigned int end1   = std::min( (int)r, (int)m_i[1] );     //min(r,m2)
					//Rcpp::Rcout<<"start1 = "<<start1<<std::endl;
					//Rcpp::Rcout<<"start2 = "<<start2<<std::endl;
					//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;

		std::vector<double> log_a(end1-start1+1, -inf);    // This vector contains all the quantities that depend only on r1
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop
		unsigned int outer_indx{0};
		for(std::size_t r1=start1; r1 <= end1; ++r1){

			// Compute a_r1 using its definition
			log_a[outer_indx] =  absC1[r-r1];

			// Prepare for computing the second term
			const unsigned int end2{r-r1};     								//r-r1, this is just for coherence of notation, no operation is needed
			// Initialize vector of results
			std::vector<double> log_vect_res(end2-start2+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);

			// Inner loop on r2
			unsigned int inner_indx{0};
			for(std::size_t r2=start2; r2<= end2; ++r2){

				// Compute
				log_vect_res[inner_indx] = gsl_sf_lnchoose(r-r2,r1) + my_log_falling_factorial(r-r1-r2,(double)(r-r1)) +  absC2[r-r2];

				// Check if it is the new maximum of log_vect_res
	        	if(log_vect_res[inner_indx]>val_max2){
	        		idx_max2 = inner_indx;
	        		val_max2 = log_vect_res[inner_indx];
	        	}
	        			//Rcpp::Rcout<<"Computing for r1 = "<<r1<<" and r2 = "<<r2<<std::endl;
				inner_indx++;
			}


			// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
			log_a[outer_indx] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);

			// Check if it is the new maximum of log_a
	       	if(log_a[outer_indx]>val_max1){
	       		idx_max1 = outer_indx;
	       		val_max1 = log_a[outer_indx];
	       	}
	       	outer_indx++;
		}

		// Complete the sum over all elements in log_a
		return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	}
}

// Per d>2
double compute_Kpost_unnormalized_recursive(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i,
						 		 		    const std::vector<double>& gamma)
{
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kpost_unnormalized_recursive, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized_recursive, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");


	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //check it is a prior case
		Rcpp::Rcout<<"The compute_Kpost_unnormalized_recursive has been called but vector n_i of previous observations is made of all zeros. Call the compute_Kprior_unnormalized_recursive function instead"<<std::endl;
		return compute_Kprior_unnormalized_recursive(r, m_i, gamma );
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) // check that the known data are coherent
		throw std::runtime_error("Error in compute_Kpost_unnormalized_recursive. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	if( r > std::accumulate(m_i.cbegin(), m_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;

	if( *std::max_element(m_i.cbegin(),m_i.cend()) == 0 ){
		    if( r == 0 )
		    	return 0.0;
		    else
		    	return -inf; //This case is handle just for completeness and to sake of clarity. However, the case n1=n2=0 and r>0 is handled in the previous if (r > n1+n2)
	}

	// if d=1 or d=2
	if( n_i.size() <= 2)
		return compute_Kpost_unnormalized(r, k, m_i, n_i, gamma);

	// if here, n_i.size()>2
	std::vector<double> log_a(r+1, -inf);
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);

	// Start for loop
	for(std::size_t k1=0; k1 <= r; ++k1){
				//Rcpp::Rcout<<"----> Dentro a k1 = "<<k1<<std::endl;

		// Compute recursive formula
		log_a[k1] = compute_Kpost_unnormalized_recursive(k1, k, {m_i.cbegin(), m_i.cend()-1}, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} );


		// Prepare for computing the second term

		// Initialize vector of results
		std::vector<double> log_vect_res(k1+1, -inf );

		// Initialize quantities to find the maximum of log_vect_res
		unsigned int idx_max2(0);
		double val_max2(log_vect_res[idx_max2]);

		// Inner loop on r2
		unsigned int inner_indx{0};
		for(std::size_t k2=r-k1; k2<= r; ++k2){
					//Rcpp::Rcout<<"***** Dentro a k2 = "<<k2<<std::endl;

			// Compute coefficient and non-normalized prob
			log_vect_res[inner_indx] = gsl_sf_lnchoose(k2,r-k1) + my_log_falling_factorial(k1+k2-r, k1) +
									   compute_Kpost_unnormalized( k2, k, {m_i[m_i.size()-1]}, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );

		 	// Check if it is the new maximum of log_vect_res
	       	if(log_vect_res[inner_indx]>val_max2){
	       		idx_max2 = inner_indx;
	       		val_max2 = log_vect_res[inner_indx];
	       	}
	   		inner_indx++;
		}

		// Update outer vector
		log_a[k1] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);

		// Check if it is the new maximum of log_a
	   	if(log_a[k1]>val_max1){
	    	idx_max1 = k1;
	      	val_max1 = log_a[k1];
	    }
	}

	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);
}


//Direct formula per d=2
// r = distinct in new sample
// t = shared in new sample
// k = distinct in observed sample
double compute_SK_post_unnormalized(const unsigned int& r, const unsigned int& t, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i,
						 		    const std::vector<double>& gamma)
{

			//Rcpp::Rcout<<"Dentro compute_SK_post_unnormalized"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_post_unnormalized, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_SK_post_unnormalized, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_post_unnormalized, the length of n_i (group sizes) and gamma has to be equal");
	// d=1
	if(n_i.size() == 1)
		throw std::runtime_error("Error in compute_SK_post_unnormalized, two populations are needed to computed the probability of shared species");

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //check it is a prior case
		Rcpp::Rcout<<"The compute_SK_post_unnormalized has been called but vector n_i of previous observations is made of all zeros. Call the compute_SK_prior_unnormalized function instead"<<std::endl;
		return compute_SK_prior_unnormalized(r, t, m_i, gamma);
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) // check that the observed data are coherent
		throw std::runtime_error("Error in compute_SK_post_unnormalized. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	if( r > std::accumulate(m_i.cbegin(), m_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;

	if( *std::max_element(m_i.cbegin(),m_i.cend()) == 0 ){ // m1=m2=0 case
		    if( r == 0 && t == 0){
		    	throw std::runtime_error("Error in compute_SK_post_unnormalized, r and t can not be zero even if all m_j are zero. The behaviuor is indefined");
		    	return 0.0;
		    }
		    else
		    	return -inf; //This case is handle just for completeness and to sake of clarity. However, the case m1=m2=0 and r>0 is handled in the previous if (r > m1+m2)
	}

	// r<t case, more shared than distinct. This implies that the only possibility when r=0 is t=0
	if(r<t)
		return -inf;
	// t > min(m1,m2) case, more shared than observation in one group
	if( t > *std::min_element(m_i.cbegin(),m_i.cend()) )
		return -inf;


	//If here, need to compute the formula
	// Compute all C numbers required
	Rcpp::NumericVector absC1 = compute_logC(  m_i[0], -gamma[0], - ( (double)k*gamma[0] + (double)n_i[0] )  ); //absC[i] = |C(m1,i,-gamma1,-(k*gamma1 + n1))| for i = 0,...,m1
	Rcpp::NumericVector absC2 = compute_logC(  m_i[1], -gamma[1], - ( (double)k*gamma[1] + (double)n_i[1] )  ); //absC[i] = |C(m2,i,-gamma2,-(k*gamma2 + n2))| for i = 0,...,m2

	//Compute and check the range of r1
	const unsigned int start1 = std::max( 0, (int)r - (int)m_i[1] ); //max(0,r-m2)
	const unsigned int end1   = std::min( (int)(r-t), (int)m_i[0] - (int)t ); //min(r-t,m1-t)
	const int range_size = (int)end1-(int)start1+1;
				//Rcpp::Rcout<<"start1 = "<<start1<<std::endl;
				//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;
				//Rcpp::Rcout<<"Range size = "<<range_size<<std::endl;

	if(range_size <= 0) //All Cnumbers would be 0
		return -inf;

	std::vector<double> log_a(range_size, -inf);    // This vector contains all the quantities that depend only on r1
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);
	// Start for loop
	unsigned int outer_indx{0};
	for(std::size_t r1=start1; r1 <= end1; ++r1){
		// Compute a_r1 using its definition
		log_a[outer_indx] = gsl_sf_lnchoose(r-r1,t) + my_log_falling_factorial(t,(double)(t+r1)) + absC1[t+r1] + absC2[r-r1] ;

		// Check if it is the new maximum of log_a
	       if(log_a[outer_indx]>val_max1){
	       	idx_max1 = outer_indx;
	       	val_max1 = log_a[outer_indx];
	       }
	       outer_indx++;
	}
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	//return -1.0;
}

//Recursive formula for d>2
// r = distinct in new sample
// sigma = shared in new sample
// k = distinct in observed sample
double compute_SK_post_unnormalized_recursive(const unsigned int& r, const unsigned int& sigma, const unsigned int& k, const std::vector<unsigned int>& m_i,
											  const std::vector<unsigned int>& n_i, const std::vector<double>& gamma)
{

			//Rcpp::Rcout<<"Dentro compute_SK_post_unnormalized_recursive"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");

	// d=1
	if(n_i.size() == 1){
		throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive, two populations are needed to computed the probability of shared species");
		//Se finisco qua, una possibilità è chiamare compute_Kprior_unnormalized
	}
	// d=2 case
	if(n_i.size() == 2)
		return compute_SK_post_unnormalized(r, sigma, k, m_i, n_i, gamma);

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //check it is a prior case
		Rcpp::Rcout<<"The compute_SK_post_unnormalized_recursive has been called but vector n_i of previous observations is made of all zeros. Call the compute_SK_prior_unnormalized_recursive function instead"<<std::endl;
		return compute_SK_prior_unnormalized_recursive(r, sigma, m_i, gamma);
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) // check that the observed data are coherent
		throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	if( r > std::accumulate(m_i.cbegin(), m_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;

	if( *std::max_element(m_i.cbegin(),m_i.cend()) == 0 ){ // m1=m2=0 case
		    if( r == 0 && sigma == 0){
		    	throw std::runtime_error("Error in compute_SK_post_unnormalized_recursive, r and sigma can not be zero even if all m_j are zero. The behaviuor is indefined");
		    	return 0.0;
		    }
		    else
		    	return -inf; //This case is handle just for completeness and to sake of clarity. However, the case m1=m2=0 and r>0 is handled in the previous if (r > m1+m2)
	}

	// r<sigma case, more shared than distinct. This implies that the only possibility when r=0 is sigma=0
	if(r<sigma)
		return -inf;
	// sigma > min(m1,...,md) case, more shared than observation in one group
	if( sigma > *std::min_element(m_i.cbegin(),m_i.cend()) )
		return -inf;


	// if here, n_i.size()>2 in a non trivial case

	// Initialize quantities for sum on s1
	std::vector<double> log_s1(r-sigma+1, -inf);
	// Initialize quantities to find the maximum of log_s1
	unsigned int idx_max0(0);
	double val_max0(log_s1[idx_max0]);

	// Start loop on s1
	for(std::size_t s1 = 0; s1<=r-sigma; ++s1){
				//Rcpp::Rcout<<"^^^^ Dentro s1 = "<<s1<<std::endl;
		// Compute coeff that depends only on s1
		log_s1[s1] = my_log_falling_factorial(sigma,(double)(sigma+s1));

		// Initialize quantities for inner sums
		const unsigned int r_sigma_s1 = r-sigma-s1;
		// Compute the log of the inner sum on h and t
		std::vector<double> log_a(r_sigma_s1+1, -inf);
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop on h
		for(std::size_t h=0; h <= r_sigma_s1; ++h){

					//Rcpp::Rcout<<"----> Dentro a h = "<<h<<std::endl;
			// Compute recursive formula
			log_a[h] = compute_SK_post_unnormalized_recursive(r-h, sigma+s1, k,{m_i.cbegin(), m_i.cend()-1}, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} );
					//Rcpp::Rcout<<"log_a["<<h<<"] = "<<"NNP(S^(2) = "<<s+s1<<", K^(2) = "<<k-h<<") = "<<log_a[h]<<std::endl;


			// Prepare for computing the second term
			const unsigned int r_sigma_s1_h = r-sigma-s1-h;
			// Initialize vector of results
			std::vector<double> log_vect_res(r_sigma_s1_h+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);

			// Inner loop on t
			unsigned int inner_indx{0};
			for(std::size_t t=0; t<= r_sigma_s1_h; ++t){
						//Rcpp::Rcout<<"***** Dentro a t = "<<t<<std::endl;

				// Compute last element
							//Rcpp::Rcout<<"********************"<<std::endl;
							//Rcpp::Rcout<<"Coef = "<<std::exp(my_log_falling_factorial(s,(double)(s+s1)) + gsl_sf_lnchoose(r_sigma_s1_h,t) + gsl_sf_lnchoose(s+h+t,s) + my_log_falling_factorial(t,(double)(h+t)))<<std::endl;
							//Rcpp::Rcout<<"Coef giusto = "<< ( gsl_sf_fact(s+s1)*gsl_sf_fact(k-s-s1-h)*gsl_sf_fact(s+h+t) )/(gsl_sf_fact(s)*gsl_sf_fact(s1)*gsl_sf_fact(h)*gsl_sf_fact(t)*gsl_sf_fact(k-s-s1-h-t)  ) <<std::endl;
							//Rcpp::Rcout<<"--------------------"<<std::endl;

				log_vect_res[inner_indx] = gsl_sf_lnchoose(r_sigma_s1_h,t) + gsl_sf_lnchoose(sigma+h+t,sigma) + my_log_falling_factorial(t,(double)(h+t)) +
										   compute_Kpost_unnormalized( sigma+h+t, k, {m_i[m_i.size()-1]}, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );
			 	// Check if it is the new maximum of log_vect_res
		       	if(log_vect_res[inner_indx]>val_max2){
		       		idx_max2 = inner_indx;
		       		val_max2 = log_vect_res[inner_indx];
		       	}
		   		inner_indx++;
			}
			// Update outer vector
			log_a[h] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
						//Rcpp::Rcout<<"Updated log_a["<<h<<"] = "<<log_a[h]<<std::endl;
			// Check if it is the new maximum of log_a
		   	if(log_a[h]>val_max1){
		    	idx_max1 = h;
		      	val_max1 = log_a[h];
		    }
		}

		// End of loop on h. Update the sum on s1
		log_s1[s1] += log_stable_sum(log_a, TRUE, val_max1, idx_max1);

		// Check if it is the new maximum of log_s1
	   	if(log_s1[s1]>val_max0){
	    	idx_max0 = s1;
	      	val_max0 = log_s1[s1];
	    }

	}


	// Complete the sum over all elements in log_s1
	return log_stable_sum(log_s1, TRUE, val_max0, idx_max0);
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Prior selection
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Wrap the call for the ComponentPrior from Rcpp objects to c++ object
std::unique_ptr< ComponentPrior > Wrapper_ComponentPrior(const Rcpp::String& prior, const Rcpp::List& prior_param)
{
	// Component prior preliminary operations
	ComponentPrior_Parameters qM_params;
	if(prior == "Poisson"){
		qM_params.Lambda = prior_param["lambda"];
	}
	else if(prior == "NegativeBinomial"){
		qM_params.p = prior_param["p"];
		qM_params.n_succ = prior_param["r"];
	}
	else{
		throw std::runtime_error("Error in Wrapper_ComponentPrior, not implemented prior requested by R function");
	}

	//Rcpp::Rcout<<"Print ComponentPrior_Parameters: qM_params.Lambda = "<<qM_params.Lambda<<"; qM_params.p = "<<qM_params.p<<"; qM_params.n_succ = "<<qM_params.n_succ<<std::endl;
	return Select_ComponentPrior(prior, qM_params);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing probabilities
//------------------------------------------------------------------------------------------------------------------------------------------------------


double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior,
						  const Rcpp::List& prior_param, unsigned int M_max  )
{

	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Compute normalization constant
			//Rcpp::Rcout<<"Calcolo log_V:"<<std::endl;
	double log_V{ compute_log_Vprior(k, n_i, gamma, qM, M_max) };
			//Rcpp::Rcout<<"log_V = "<<log_V<<std::endl;

	// Compute unnormalized probability
			//Rcpp::Rcout<<"Calcolo log_K:"<<std::endl;
	double log_K{compute_Kprior_unnormalized_recursive(k, n_i, gamma)};
			//Rcpp::Rcout<<"log_K = "<<log_K<<std::endl;

	//return
	return std::exp(log_V + log_K);
}


double p_shared_prior_c(const unsigned int& s, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior,
					 	const Rcpp::List& prior_param, unsigned int M_max  )
{

	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Need to compute the probability for each possible K and marginalize it out
	double res{0.0};
	const unsigned int Kmax{ std::accumulate(n_i.cbegin(), n_i.cend(), 0) };

	const unsigned int start{ std::max(1,(int)s) }; //max(1,s)
	//start can not be 0! Why?
	// Because the problem is that when k=0, it does not make any sense to compute compute_log_Vprior(0), one gets NaN. Actually it rises an error at the moment.
	// However, it is possible that s=0 and k should range form 0 to n1+n2
	// Of course, P(K=0,S=0) is zero (at least, it is zero when n1+n2>0). This term is computed separately and should be 0 (or it can be omitted)
	if(s==0){ //can be omitted in this function
		res += std::exp(compute_SK_prior_unnormalized_recursive(0, s, n_i, gamma));
	}
	for(unsigned int k=start; k<=Kmax; k++){
				//Rcpp::Rcout<<"++++ k = "<<k<<std::endl;
		// Compute normalization constant
				//Rcpp::Rcout<<"Calcolo log_V:"<<std::endl;
		double log_V{ compute_log_Vprior(k, n_i, gamma, qM, M_max) };
				//Rcpp::Rcout<<"log_V = "<<log_V<<std::endl;

		// Compute unnormalized probability
				//Rcpp::Rcout<<"Calcolo log_SK:"<<std::endl;
		double log_SK{compute_SK_prior_unnormalized_recursive(k, s, n_i, gamma)};
				//Rcpp::Rcout<<"log_SK = "<<log_SK<<std::endl;

		// Compute P(K=k,S=s) and cumulate
		res += std::exp(log_V + log_SK);
				//Rcpp::Rcout<<"P(K="<<k<<", S="<<s<<") = "<<std::exp(log_V + log_SK)<<std::endl;
				//Rcpp::Rcout<<"res:"<<std::endl<<res<<std::endl;
	}

	//return
	return res;
}

double p_distinct_posterior_c(const unsigned int& r, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j,
						      const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> m_i   = Rcpp::as< std::vector<unsigned int> >(m_j);
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Compute normalization constant
			//Rcpp::Rcout<<"Calcolo log_V posteriori:"<<std::endl;
	double log_Vpost{ compute_log_Vpost(r, k, m_i, n_i, gamma, qM, M_max ) };
			//Rcpp::Rcout<<"log_Vpost = "<<log_Vpost<<std::endl;

			// That was just a check
			//Rcpp::Rcout<<"Calcolo NAIVE log_V posteriori:"<<std::endl;
			//double log_Vpost_naive{ compute_log_Vpost_naive(r, k, m_i, n_i, gamma, qM, M_max ) };
			//Rcpp::Rcout<<"log_Vpost_NAIVE = "<<log_Vpost_naive<<std::endl;

	// Compute unnormalized probability
			//Rcpp::Rcout<<"Calcolo log_K posteriori:"<<std::endl;
	double log_Kpost{compute_Kpost_unnormalized_recursive(r, k, m_i, n_i, gamma)};
			//Rcpp::Rcout<<"log_Kpost = "<<log_Kpost<<std::endl;

	//return
	//return std::exp(log_Vpost + log_Kpost);
	return std::exp(log_Vpost + log_Kpost);
}

double p_shared_posterior_c(const unsigned int& t, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j,
						    const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> m_i   = Rcpp::as< std::vector<unsigned int> >(m_j);
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Need to compute the probability for each possible r and marginalize it out
	double res{0.0};
	const unsigned int r_max{ std::accumulate(m_i.cbegin(), m_i.cend(), 0) };

	// r ranges from t to m1+...+md
	// when t=0, it is possible that r=0.
	for(unsigned int r=t; r<=r_max; ++r){
				//Rcpp::Rcout<<"++++ r = "<<r<<std::endl;

		// Compute normalization constant
				//Rcpp::Rcout<<"Calcolo log_Vpost:"<<std::endl;
		double log_Vpost{ compute_log_Vpost(r, k, m_i, n_i, gamma, qM, M_max ) };
				//Rcpp::Rcout<<"log_Vpost = "<<log_Vpost<<std::endl;

		// Compute unnormalized probability
				//Rcpp::Rcout<<"Calcolo log_SK_post:"<<std::endl;
		double log_SK_post{compute_SK_post_unnormalized_recursive(r, t, k, m_i, n_i, gamma)};
				//Rcpp::Rcout<<"log_SK_post = "<<log_SK_post<<std::endl;

		// Compute P(K=r,S=t|X) and cumulate
		res += std::exp(log_Vpost + log_SK_post);
				//Rcpp::Rcout<<"P(K="<<r<<", S="<<r<<" | X ) = "<<std::exp(log_Vpost + log_SK_post)<<std::endl;
				//Rcpp::Rcout<<"res:"<<std::endl<<res<<std::endl;
	}

	//return
	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing expected values
//------------------------------------------------------------------------------------------------------------------------------------------------------


Rcpp::List Expected_prior_c(const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& type, const Rcpp::String& prior,
					    const Rcpp::List& prior_param, unsigned int M_max, double tol  )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	double cumulative_prob{0.0};
	double first_moment{0.0};
	double second_moment{0.0};
	double variance{0.0};
	const double stop_criteria{1.0-tol};
	// distinct case
	if(type == "distinct"){
		const unsigned int Kmax = *std::max_element(n_i.cbegin(),n_i.cend()) ;
		unsigned int k{1};

		while(k<=Kmax && cumulative_prob < stop_criteria){
			const double prob = std::exp(compute_log_Vprior(k, n_i, gamma, qM, M_max) + compute_Kprior_unnormalized_recursive(k, n_i, gamma));
			cumulative_prob += prob;
			first_moment += k*prob;
			second_moment += k*k*prob;
			k++;
			//Check for User Interruption
			try{
			    Rcpp::checkUserInterrupt();
			}
			catch(Rcpp::internal::InterruptedException e){
			    //Print error and return
			    throw std::runtime_error("Execution stopped by the user");
			}
		}
	}
	else if(type == "shared"){ // shared case
		const unsigned int Smax = *std::min_element(n_i.cbegin(),n_i.cend()) ;
		unsigned int s{0}; // s=0 is useless from the point of view of computing the expected value but it is usefull to compute the cumulative_prob

		while(s<=Smax && cumulative_prob < stop_criteria){
			const double prob = p_shared_prior_c(s, n_j, gamma_j, prior, prior_param, M_max  ); // not the most efficient call
			cumulative_prob += prob;
			first_moment += s*prob;
			second_moment += s*s*prob;
			s++;
		}
	}
	else{
		throw std::runtime_error("Error in Expected_prior_c. type can only be equal to distinct or shared ");
		return -1.0;
	}
	variance = second_moment - (first_moment*first_moment);
	if(variance < 0)
		throw std::runtime_error("Error in Expected_prior_c, the variance can not be negative ");
	return Rcpp::List::create( Rcpp::Named("Mean") = first_moment, Rcpp::Named("Variance") = variance);
}


Rcpp::List Expected_posterior_c(const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j,
						    const Rcpp::String& type, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max, double tol)
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> m_i   = Rcpp::as< std::vector<unsigned int> >(m_j);
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	double cumulative_prob{0.0};
	double first_moment{0.0};
	double second_moment{0.0};
	double variance{0.0};
	const double stop_criteria{1.0-tol};

	// distinct case
	if(type == "distinct"){
		const unsigned int r_max = *std::max_element(m_i.cbegin(),m_i.cend()) ;
		unsigned int r{0}; // it is possible that r=0. This case is useless to compute the expected value but it is useful to compute cumulative_prob

		while(r<=r_max && cumulative_prob < stop_criteria){
			const double prob = std::exp(compute_log_Vpost(r, k, m_i, n_i, gamma, qM, M_max ) + compute_Kpost_unnormalized_recursive(r, k, m_i, n_i, gamma));
			cumulative_prob += prob;
			//Rcpp::Rcout<<"r:"<<std::endl<<r<<std::endl;
			//Rcpp::Rcout<<"prob:"<<std::endl<<prob<<std::endl;
			//Rcpp::Rcout<<"cumulative_prob:"<<std::endl<<cumulative_prob<<std::endl;
			first_moment += r*prob;
			second_moment += r*r*prob;
			r++;
		}
	}
	else if(type == "shared"){ // shared case
		const unsigned int sigma_max = *std::min_element(m_i.cbegin(),m_i.cend()) ;
		unsigned int sigma{0}; // sigma=0 is useless from the point of view of computing the expected value but it is usefull to compute the cumulative_prob

		while(sigma<=sigma_max && cumulative_prob < stop_criteria){

			const double prob = p_shared_posterior_c(sigma, k, m_j, n_j, gamma_j, prior, prior_param, M_max); // not the most efficient call
			cumulative_prob += prob;
			Rcpp::Rcout<<"sigma:"<<std::endl<<sigma<<std::endl;
			Rcpp::Rcout<<"prob:"<<std::endl<<prob<<std::endl;
			Rcpp::Rcout<<"cumulative_prob:"<<std::endl<<cumulative_prob<<std::endl;
			first_moment += sigma*prob;
			second_moment += sigma*sigma*prob;
			sigma++;
		}
	}
	else{
		throw std::runtime_error("Error in Expected_posterior_c. type can only be equal to distinct or shared ");
		return -1.0;
	}
	variance = second_moment - (first_moment*first_moment);
	if(variance < 0)
		throw std::runtime_error("Error in Expected_prior_c, the variance can not be negative ");
	return Rcpp::List::create( Rcpp::Named("Mean") = first_moment, Rcpp::Named("Variance") = variance);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper bounds
//------------------------------------------------------------------------------------------------------------------------------------------------------

Rcpp::NumericVector Sums_logC(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j)
{
	// checks
	if(n_j.size() != 2)
		throw std::runtime_error("Error in LowerBounds: only d=2 case is implemented");
	if(gamma_j.size() != 2)
		throw std::runtime_error("Error in LowerBounds: only d=2 case is implemented");

	const unsigned int n{n_j[0]+n_j[1]};
	unsigned int K{1};
	std::vector< std::vector<double> > a_old;
	std::vector<double> b;
	Rcpp::NumericVector res(n);
	double inf = std::numeric_limits<double>::infinity();
	// Compute all C numbers required
	Rcpp::NumericVector absC1 = compute_logC(n_j[0], -gamma_j[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_j[1], -gamma_j[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2
	std::vector<double> absC1v(n+1,-inf);
	std::copy(absC1.begin(),absC1.end(),absC1v.begin());

	//1) First case: K=1
	std::vector< std::vector<double> > a(K+1);
	a[0].resize(1);
	a[0][0] = absC2[1];

	a[1].resize(2);
	a[1][0] = absC2[1];
	a[1][1] = absC2[0];

	// compute log double sum
	b.clear();
	b.resize(K+1, -inf);
	std::transform(a.cbegin(),a.cend(),absC1v.cbegin(),b.begin(),[](const std::vector<double>& a_i, const double& lC_1i){ return(log_stable_sum(a_i,TRUE)+lC_1i); });
	res[K-1] = log_stable_sum(b,TRUE);

	// get ready for the next K
	K++;
	a_old = a;

	//2) All other cases
	while(K<=n){
		a.resize(K+1);
		for(unsigned int k = 0; k <= K; ++k){
			a[k].clear();
			a[k].resize(1+k, -inf); // resize correct number of terms and set all elements to -inf
			if(K <= n_j[1]) // fill the second element (if possible)
				a[k][0] = absC2[K];
			if(k > 0) // copy from old vector
				std::copy(a_old[k-1].begin(),a_old[k-1].end(), a[k].begin()+1);

			//Rcpp::Rcout<<"Stampo a["<<k<<"]: ";
				//for(auto __v : a[k])
					//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;

			//Check for User Interruption
			try{
			    Rcpp::checkUserInterrupt();
			}
			catch(Rcpp::internal::InterruptedException e){
			    //Print error and return
			    throw std::runtime_error("Execution stopped by the user");
			}

		}

		// compute log double sum
		b.clear();
		b.resize(K+1, -inf);
		std::transform(a.cbegin(),a.cend(),absC1v.cbegin(),b.begin(),[](const std::vector<double>& a_i, const double& lC_1i){ return(log_stable_sum(a_i,TRUE)+lC_1i); });
		res[K-1] = log_stable_sum(b,TRUE);
				//Rcpp::Rcout<<"res["<<K-1<<"] = "<<res[K-1]<<std::endl;
				//throw std::runtime_error("FERMO IO");
		// get ready for the next K
		K++;
		a_old = a;
	}

	return(res);
}

Rcpp::NumericVector UpperBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior,
					              const Rcpp::List& prior_param, unsigned int M_max  )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;
	const unsigned int n = n_j[0]+n_j[1];
	Rcpp::NumericVector res(n);

	// Compute V
	double log_V{0.0};
	double log_UBsums{0.0};
	for(unsigned int k=1;k<=n;++k){

		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max);
		// Compute upper bound for V
		//log_UBV = 0.5*std::log(k/(2.0*M_PI)) + 2.0 - n_j[0]*std::log(gamma_j[0]) - n_j[1]*std::log(gamma_j[1]);
		//if(k < n-1)
			//log_UBV += std::log(gsl_sf_zetam1_int(n-k));//gsl_sf_hzeta( (double)n - (double)k, (double)k );
		//else
			//log_UBV += std::log(qM.eval_upper_tail(k));

		log_UBsums = -gsl_sf_lnfact(k) + std::log(gamma_j[0]) + std::log(gamma_j[1]) + 2.0*std::log((double)k) +
						  ((double)n_j[0]-1.0)*std::log(gamma_j[0]*(double)k + n_j[0]) + ((double)n_j[1]-1.0)*std::log(gamma_j[1]*(double)k + n_j[1]);

		res[k-1] = log_V + log_UBsums;
	}

	//return
	return res;
}

Rcpp::NumericVector LowerBounds_c(const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j, const Rcpp::String& prior,
					              const Rcpp::List& prior_param, unsigned int M_max  )
{
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;
	const unsigned int n = n_j[0]+n_j[1];

	// Compute V
	Rcpp::NumericVector V_vect(n);
	for(unsigned int k=1;k<=n;++k){
		V_vect[k-1] = compute_log_Vprior(k, n_j, gamma_j, qM, M_max);
	}

	// Compute Sums C numbers
	Rcpp::NumericVector SumslogC = Sums_logC(n_j, gamma_j);
	Rcpp::NumericVector res = SumslogC + V_vect;

	//return
	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Compute whole distributions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Compute the whole distribution for the prior number of distinct components
Rcpp::NumericVector D_distinct_prior_c( const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j,
										const Rcpp::String& prior, const Rcpp::List& prior_param,
										unsigned int M_max,
										const int& Kstart,
										std::vector<double>& logV_vec   )
{
	double inf = std::numeric_limits<double>::infinity();
	// checks
	if(n_j.size() != 2)
		throw std::runtime_error("Error in D_distinct_prior_c: current implementation is only for d=2 groups ");
	const unsigned int n = std::accumulate(n_j.cbegin(),n_j.cend(), 0.0);

	// Initialize return quantities
	Rcpp::NumericVector res(n+1, 0.0);
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Compute all C numbers required
	Rcpp::Rcout<<"Compute C numbers ... ";
	Rcpp::NumericVector absC1 = compute_logC(n_j[0], -gamma_j[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_j[1], -gamma_j[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2
	Rcpp::Rcout<<" done! "<<std::endl;

	// Define grid search for k
	const unsigned int max_diff{std::max(Kstart-1,(int)n-Kstart)};
	std::vector<unsigned int> Ksearch;
	Ksearch.reserve(n);
	Ksearch.push_back( (unsigned int)Kstart );
	for(std::size_t i = 1; i <= max_diff; i++){
		if( (Kstart-i) > 0)
			Ksearch.push_back( (unsigned int)(Kstart-i) );
		if( (Kstart+i) <= n)
			Ksearch.push_back( (unsigned int)(Kstart+i) );
	}

	// Cycle for each k in Ksearch
	double log_V{0.0};
	double log_K{0.0};
	double cumulated{0.0};
	Progress progress_bar(n, TRUE); // Initialize progress bar
	for(unsigned int it=0; it<Ksearch.size(); ++it){
		unsigned int k = Ksearch[it];

		// compute V number if never computed before
		if(logV_vec[k] == -inf)
		    logV_vec[k] = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

		// get V number
		log_V = logV_vec[k];

		log_K = compute_Kprior_unnormalized(k, n_j, gamma_j, absC1, absC2);
		res[k] = std::exp(log_V + log_K);
		cumulated += res[k];
		//Rcpp::Rcout<<"P(K <= "<<k<<") = "<<cumulated<<std::endl;
		//Check for User Interruption
		try{
		    Rcpp::checkUserInterrupt();
		}
		catch(Rcpp::internal::InterruptedException e){
		    //Print error and return
		    throw std::runtime_error("Execution stopped by the user");
		}
		progress_bar.increment(); //update progress bar
		//Rcpp::Rcout<<"k = "<<k<<std::endl;

		if( 1.0 - cumulated < 1e-10 )
			break;
	}

	return (res);
}

Rcpp::List Distinct_Prior_MCMC( unsigned int Niter,
				                 const std::vector<unsigned int>& n_j,
				                 const std::vector<double>& gamma_j,
				                 const Rcpp::String& prior, const Rcpp::List& prior_param,
				                 unsigned int M_max,
				                 unsigned int seed
				              )
{
    // Standard definitions
    double inf = std::numeric_limits<double>::infinity();
    const unsigned int d{n_j.size()};
    unsigned int n = std::accumulate(n_j.cbegin(), n_j.cend(), 0);
    auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
    ComponentPrior& qM(*qM_ptr);

    // Set empty initial partition
    unsigned int K = 1;
    MatUnsCol N{MatUnsCol::Constant(d,K,0)};
    std::vector< std::vector<unsigned int> > z_ji(d);
    for (std::size_t j = 0; j < d; j++){
        std::vector<unsigned int> zeros(n_j[j],0);
        z_ji[j] = zeros;
        N(j,0) = n_j[j];
    }

    if(gamma_j.size() != d)
        throw std::runtime_error("Error in Distinct_Prior_MCMC: length of gamma_j must be equal to d, but it is not compatible with length of n_j");

    VecUnsCol N_k = N.colwise().sum(); // vector of length K, N_k[m] is the number of data assigned to cluster m

    // Define return objects
    std::vector<double> logV_vec(n+1,-inf); // Vector of logV, i.e., logV_vec[k] = log( V(n,k) )
    std::vector<unsigned int> K_MCMC(Niter+1,0); // Vector of MCMC values for K
    K_MCMC[0] = K;

    // Declare auxiliary quantities
    double log_probs_max;
    sample::sample_index sample_index;
    sample::GSL_RNG engine(seed); // initialize random engine with default random seed
    unsigned int new_z_ji;
    double logVnum{-inf};
    double logVden{-inf};

    Progress progress_bar(Niter, TRUE); // Initialize progress bar
    for(std::size_t it = 0; it < Niter; it++){
        //Rcpp::Rcout<<"###################################"<<std::endl;
        // Chinese Restaurant Franchise process allocation
        for(unsigned int j=0; j<d; j++)
        {
            for(unsigned int i=0; i<n_j[j]; i++)
            {

            			//Rcpp::Rcout<<"("<<j<<", "<<i<<")  ||  K = "<<K<<std::endl;
                // Shortcut to get cluster membership of observation (j,i)
                unsigned int C_ji = z_ji[j][i];

                // remove obs ji from its cluster. In this step, both the local and the global counts must be updated as well as the sums in that cluster
                N_k(C_ji)--; // decrease the global counts
                N(j,C_ji)--; // decrease the local counts

	                //Rcpp::Rcout<<"N_k = "<<N_k<<std::endl;
	                //Rcpp::Rcout<<"N_k(C_ji) = "<<N_k(C_ji)<<std::endl;

                // if the cluster becomes empty, then it must be removed.
                // This is achived by replacing the now empty cluster with the last one.
                if(N_k(C_ji) == 0){
                			//Rcpp::Rcout<<"------ Rimuovo un cluster -------"<<std::endl;
                    // last cluster to replace the now empty cluster.
                    N_k(C_ji) = N_k(K-1);        // update the global counts
                    N_k.conservativeResize(K-1);             // eliminate the last cluster, now empty
                    N.col(C_ji) = N.col(K-1);    // update the local counts
                    N.conservativeResize(d,K-1); // eliminate the last cluster, now empty
                    // change all labels according to new labeling. Data (there is only one and it is in position ji) with label C_ji is ruled out by setting its label equal to K-1
                    // while current data with label K-1 get label C_ji

                    for(unsigned int jj=0; jj<d; jj++){
                        for(unsigned int ii=0; ii<n_j[jj]; ii++){
                            if(z_ji[jj][ii] == K-1)
                                z_ji[jj][ii] = C_ji; // switch label for data currently in cluster K-1
                        }
                    }
                    z_ji[j][i] = K-1;

                    K--; // decrese the current number of clusters
                }


                VecRow log_probs_vec = VecRow::Constant(K+1, 0.0); // define vector of weights, currently in log-scale. lenght must be K+1

                // loop over current clusters
                if(K==0)
                    throw std::runtime_error("Error in Distinct_Prior_MCMC: K is 0, this should be impossible");

                for(std::size_t l = 0; l < K; l++){
                    log_probs_vec[l] =  std::log( (double)N(j,l) + gamma_j[j] ); // set prior term
                    		//Rcpp::Rcout<<"log_probs_vec["<<l<<"] = "<<log_probs_vec[l]<<std::endl;
                }

                // we are done looping over the already occupied clusters. Next, calculate the log probability of a new table.
                if(logV_vec[K+1] == -inf)
                    logV_vec[K+1] = compute_log_Vprior(K+1, n_j, gamma_j, qM, M_max );
                if(logV_vec[K] == -inf)
                    logV_vec[K] = compute_log_Vprior(K, n_j, gamma_j, qM, M_max );

                logVnum = logV_vec[K+1];
                logVden = logV_vec[K];

                		//Rcpp::Rcout<<"logVnum = "<<logVnum<<std::endl;
                		//Rcpp::Rcout<<"logVden = "<<logVden<<std::endl;
                log_probs_vec[K] =  std::log(gamma_j[j]) + logVnum - logVden;   // prior term
						//Rcpp::Rcout<<"log_probs_vec["<<K<<"] = "<<log_probs_vec[K]<<std::endl;

                // stable calculation of non-normalized weights in non-log scale
                log_probs_max = log_probs_vec.maxCoeff(); // get max value

                for(std::size_t m = 0; m < log_probs_vec.size(); m++){
                    log_probs_vec(m) = std::exp(log_probs_vec(m) - log_probs_max);
                     //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(log_probs_vec(m)))
                        throw std::runtime_error("Error in Distinct_Prior_MCMC, get a nan in probs_vec ");
                }

                		//Rcpp::Rcout<<"Probs: "<<log_probs_vec/log_probs_vec.sum()<<std::endl;
                // Draw a sample of which cluster customer ji should belong to
                new_z_ji = sample_index(engine, log_probs_vec); //values in log_probs_vec are no longer in log-scale here
                		//Rcpp::Rcout<<"new_z_ji = "<<new_z_ji<<std::endl;
                // set a new cluster, if necessary
                if( new_z_ji == K){
                			//Rcpp::Rcout<<"+++++++ Aggiungo un cluster +++++++"<<std::endl;
                			//throw std::runtime_error("FERMO IO - Aggiungo ");
                    N_k.conservativeResize(K+1);             // allocate space in global counts for the new cluster
                    N_k(K) = 0;                              // values are assigned later
                    N.conservativeResize(d,K+1);             // allocate space in local counts for the new cluster
                    N.col(K) = VecUnsCol::Constant(d,0);     // values are assigned later. Note that here we are setting empty tables in all the other restaurants
                    K++;
                }

                //Assign cluster membership and update counts
                z_ji[j][i] = new_z_ji;   // set new label
                N_k(new_z_ji)++;         // update global counts
                N(j,new_z_ji)++;         // update local counts

                if(N.sum() != n){
                	Rcpp::Rcout<<"N.sum():"<<std::endl<<N.sum()<<std::endl;
                	throw std::runtime_error("Error in Distinct_Prior_MCMC: N_k is not summing to n ");
                }
                if(N_k.sum() != n){
                	Rcpp::Rcout<<"N_k.sum():"<<std::endl<<N_k.sum()<<std::endl;
                	throw std::runtime_error("Error in Distinct_Prior_MCMC: N_k is not summing to n ");
                }
            }


        }

        progress_bar.increment(); //update progress bar

        if(K==0)
            throw std::runtime_error("Error in Distinct_Prior_MCMC: K is 0, this should be impossible");

        // Save MCMC
        K_MCMC[it+1] = K;

        //Check for User Interruption
        try{
            Rcpp::checkUserInterrupt();
        }
        catch(Rcpp::internal::InterruptedException e){
            //Print error and return
            throw std::runtime_error("Execution stopped by the user");
        }
    }

    return Rcpp::List::create( Rcpp::Named("K_MCMC") = K_MCMC, Rcpp::Named("logV_vec") = logV_vec);
}

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
				                )
{

    // Standard definitions
    double inf = std::numeric_limits<double>::infinity();
    const unsigned int d{n_j.size()};
    unsigned int n = std::accumulate(n_j.cbegin(), n_j.cend(), 0);
    auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
    ComponentPrior& qM(*qM_ptr);

    // Set empty initial partition
    unsigned int K{K0};
    std::vector< std::vector<unsigned int> > z_ji(d);
    std::vector<unsigned int> idx_temp{n_j};
    idx_temp.insert(idx_temp.begin(), 0);
    unsigned int cumsum = 0;
    for (std::size_t j = 0; j < d; j++){
    	unsigned int start = cumsum;
    	unsigned int end   = start + n_j[j];
    	cumsum += n_j[j];
    	z_ji[j].insert(z_ji[j].begin(), rho0.cbegin() + start, rho0.cbegin() + end);
    }
    if(gamma_j.size() != d)
        throw std::runtime_error("Error in Distinct_Prior_MCMC: length of gamma_j must be equal to d, but it is not compatible with length of n_j");

    //VecUnsCol N_k = N.colwise().sum(); // vector of length K, N_k[m] is the number of data assigned to cluster m

    // Define return objects
    std::vector<double> logV_vec(n+1,-inf); // Vector of logV, i.e., logV_vec[k] = log( V(n,k) )
    std::vector<unsigned int> K_MCMC(Niter+1,0); // Vector of MCMC values for K
    K_MCMC[0] = K;

    // Declare auxiliary quantities
    double log_probs_max;
    sample::sample_index sample_index;
    sample::GSL_RNG engine(seed); // initialize random engine with default random seed
    unsigned int new_z_ji;
    double logVnum{-inf};
    double logVden{-inf};

    Progress progress_bar(Niter, TRUE); // Initialize progress bar
    for(std::size_t it = 0; it < Niter; it++){
        //Rcpp::Rcout<<"###################################"<<std::endl;
        // Chinese Restaurant Franchise process allocation
        for(unsigned int j=0; j<d; j++)
        {
            for(unsigned int i=0; i<n_j[j]; i++)
            {

            			//Rcpp::Rcout<<"("<<j<<", "<<i<<")  ||  K = "<<K<<std::endl;
                // Shortcut to get cluster membership of observation (j,i)
                unsigned int C_ji = z_ji[j][i];

                // remove obs ji from its cluster. In this step, both the local and the global counts must be updated as well as the sums in that cluster
                N_k(C_ji)--; // decrease the global counts
                N(j,C_ji)--; // decrease the local counts

	                //Rcpp::Rcout<<"N_k = "<<N_k<<std::endl;
	                //Rcpp::Rcout<<"N_k(C_ji) = "<<N_k(C_ji)<<std::endl;

                // if the cluster becomes empty, then it must be removed.
                // This is achived by replacing the now empty cluster with the last one.
                if(N_k(C_ji) == 0){
                			//Rcpp::Rcout<<"------ Rimuovo un cluster -------"<<std::endl;
                    // last cluster to replace the now empty cluster.
                    N_k(C_ji) = N_k(K-1);        // update the global counts
                    N_k.conservativeResize(K-1);             // eliminate the last cluster, now empty
                    N.col(C_ji) = N.col(K-1);    // update the local counts
                    N.conservativeResize(d,K-1); // eliminate the last cluster, now empty
                    // change all labels according to new labeling. Data (there is only one and it is in position ji) with label C_ji is ruled out by setting its label equal to K-1
                    // while current data with label K-1 get label C_ji

                    for(unsigned int jj=0; jj<d; jj++){
                        for(unsigned int ii=0; ii<n_j[jj]; ii++){
                            if(z_ji[jj][ii] == K-1)
                                z_ji[jj][ii] = C_ji; // switch label for data currently in cluster K-1
                        }
                    }
                    z_ji[j][i] = K-1;

                    K--; // decrese the current number of clusters
                }


                VecRow log_probs_vec = VecRow::Constant(K+1, 0.0); // define vector of weights, currently in log-scale. lenght must be K+1

                // loop over current clusters
                if(K==0)
                    throw std::runtime_error("Error in Distinct_Prior_MCMC: K is 0, this should be impossible");

                for(std::size_t l = 0; l < K; l++){
                    log_probs_vec[l] =  std::log( (double)N(j,l) + gamma_j[j] ); // set prior term
                    		//Rcpp::Rcout<<"log_probs_vec["<<l<<"] = "<<log_probs_vec[l]<<std::endl;
                }

                // we are done looping over the already occupied clusters. Next, calculate the log probability of a new table.
                if(logV_vec[K+1] == -inf)
                    logV_vec[K+1] = compute_log_Vprior(K+1, n_j, gamma_j, qM, M_max );
                if(logV_vec[K] == -inf)
                    logV_vec[K] = compute_log_Vprior(K, n_j, gamma_j, qM, M_max );

                logVnum = logV_vec[K+1];
                logVden = logV_vec[K];

                		//Rcpp::Rcout<<"logVnum = "<<logVnum<<std::endl;
                		//Rcpp::Rcout<<"logVden = "<<logVden<<std::endl;
                log_probs_vec[K] =  std::log(gamma_j[j]) + logVnum - logVden;   // prior term
						//Rcpp::Rcout<<"log_probs_vec["<<K<<"] = "<<log_probs_vec[K]<<std::endl;

                // stable calculation of non-normalized weights in non-log scale
                log_probs_max = log_probs_vec.maxCoeff(); // get max value

                for(std::size_t m = 0; m < log_probs_vec.size(); m++){
                    log_probs_vec(m) = std::exp(log_probs_vec(m) - log_probs_max);
                     //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(log_probs_vec(m)))
                        throw std::runtime_error("Error in Distinct_Prior_MCMC, get a nan in probs_vec ");
                }

                		//Rcpp::Rcout<<"Probs: "<<log_probs_vec/log_probs_vec.sum()<<std::endl;
                // Draw a sample of which cluster customer ji should belong to
                new_z_ji = sample_index(engine, log_probs_vec); //values in log_probs_vec are no longer in log-scale here
                		//Rcpp::Rcout<<"new_z_ji = "<<new_z_ji<<std::endl;
                // set a new cluster, if necessary
                if( new_z_ji == K){
                			//Rcpp::Rcout<<"+++++++ Aggiungo un cluster +++++++"<<std::endl;
                			//throw std::runtime_error("FERMO IO - Aggiungo ");
                    N_k.conservativeResize(K+1);             // allocate space in global counts for the new cluster
                    N_k(K) = 0;                              // values are assigned later
                    N.conservativeResize(d,K+1);             // allocate space in local counts for the new cluster
                    N.col(K) = VecUnsCol::Constant(d,0);     // values are assigned later. Note that here we are setting empty tables in all the other restaurants
                    K++;
                }

                //Assign cluster membership and update counts
                z_ji[j][i] = new_z_ji;   // set new label
                N_k(new_z_ji)++;         // update global counts
                N(j,new_z_ji)++;         // update local counts

                if(N.sum() != n){
                	Rcpp::Rcout<<"N.sum():"<<std::endl<<N.sum()<<std::endl;
                	throw std::runtime_error("Error in Distinct_Prior_MCMC: N_k is not summing to n ");
                }
                if(N_k.sum() != n){
                	Rcpp::Rcout<<"N_k.sum():"<<std::endl<<N_k.sum()<<std::endl;
                	throw std::runtime_error("Error in Distinct_Prior_MCMC: N_k is not summing to n ");
                }
            }


        }

        progress_bar.increment(); //update progress bar

        if(K==0)
            throw std::runtime_error("Error in Distinct_Prior_MCMC: K is 0, this should be impossible");

        // Save MCMC
        K_MCMC[it+1] = K;

        //Check for User Interruption
        try{
            Rcpp::checkUserInterrupt();
        }
        catch(Rcpp::internal::InterruptedException e){
            //Print error and return
            throw std::runtime_error("Execution stopped by the user");
        }
    }

    return Rcpp::List::create( Rcpp::Named("K_MCMC") = K_MCMC, Rcpp::Named("logV_vec") = logV_vec);
}


Rcpp::NumericMatrix D_joint_prior_c( const std::vector<unsigned int>& n_j, const std::vector<double>& gamma_j,
									 const Rcpp::String& prior, const Rcpp::List& prior_param,
									 unsigned int M_max, const int& Kstart, std::vector<double>& logV_vec   )
{
	double inf = std::numeric_limits<double>::infinity();
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;
	const unsigned int n = std::accumulate(n_j.cbegin(),n_j.cend(), 0.0);
	// Compute all C numbers required
	Rcpp::Rcout<<"Compute C numbers ... ";
	Rcpp::NumericVector absC1 = compute_logC(n_j[0], -gamma_j[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_j[1], -gamma_j[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2
	Rcpp::Rcout<<" done! "<<std::endl;
	// Convert Rcpp vector
	Rcpp::NumericMatrix res(n+1,n+1);

	//VecCol marginal_K{VecCol::Constant(n+1)};
	//VecCol marginal_S{VecCol::Constant(n+1)};

	// Define grid search for k
	const unsigned int max_diff{std::max(Kstart-1,(int)n-Kstart)};
	std::vector<unsigned int> Ksearch;
	Ksearch.reserve(n);
	Ksearch.push_back( (unsigned int)Kstart );
	for(std::size_t i = 1; i <= max_diff; i++){
		if( (Kstart-i) > 0)
			Ksearch.push_back( (unsigned int)(Kstart-i) );
		if( (Kstart+i) <= n)
			Ksearch.push_back( (unsigned int)(Kstart+i) );
	}
	double log_V{0.0};
	double log_SK{0.0};
	double joint_cumulated{0.0};
	unsigned int nelem_max = n*(n+3)/2;
	Progress progress_bar(nelem_max, TRUE); // Initialize progress bar
	for(unsigned int it=0; it<Ksearch.size(); ++it){
		unsigned int k = Ksearch[it];
		// compute V number if never computed before
		if(logV_vec[k] == -inf)
		    logV_vec[k] = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );
		// get V number
		log_V = logV_vec[k];
		// Loop over all possible s values
		for(unsigned int s = 0; s <= k; s++){
					//Rcpp::Rcout<<"(k = "<<k<<", s = "<<s<<")"<<std::endl;
			// Compute unnormalized probability
					//Rcpp::Rcout<<"Calcolo log_SK:"<<std::endl;
			log_SK = compute_SK_prior_unnormalized(k, s, n_j, gamma_j, absC1, absC2);
					//Rcpp::Rcout<<"log_SK = "<<log_SK<<std::endl;

			res(s,k) = std::exp(log_V + log_SK); // save joint probability
			joint_cumulated += res(s,k);         // update joint cumulated distribution
					//Rcpp::Rcout<<"joint_cumulated:"<<std::endl<<joint_cumulated<<std::endl;
		}
		//Check for User Interruption
		try{
		    Rcpp::checkUserInterrupt();
		}
		catch(Rcpp::internal::InterruptedException e){
		    //Print error and return
		    throw std::runtime_error("Execution stopped by the user");
		}
		progress_bar.increment(); //update progress bar
		if( 1.0 - joint_cumulated < 1e-10 )
			break;
	}
	return res;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Large n regime
//------------------------------------------------------------------------------------------------------------------------------------------------------

double log_factorial_mom(const int& r, const std::vector<double>& v_log_qM_post)
{
	if(r == 0)
		return 0.0; // if here, Expected value is 1 and its log is 0

	double inf = std::numeric_limits<double>::infinity();
	unsigned int Natoms = v_log_qM_post.size(); // compute number of atoms in qM_post whose sum is enough to reach 1
	if(Natoms == 0)
		throw std::runtime_error("Error in log_factorial_mom: v_log_qM_post is empty");

	if(r >= Natoms)
		return -inf; // if here, Expected value is 0 and its log is -inf

	std::vector<double> log_res(Natoms-r,-inf); // initialize vector used to compute the log stable sum
	int counter{0}; // initialize counter
	double max{-inf}; // initialize max value in log_res
	unsigned int idx_max{0}; // initialize position of max element in log_res

	// start loop over all atomes
	for(std::size_t mstar = r; mstar < Natoms; mstar++){
		log_res[counter] = v_log_qM_post[mstar] + my_log_falling_factorial(r,mstar);

		// check for max
		if(log_res[counter] > max){
			max = log_res[counter];
			idx_max = counter;
		}
	}

	// compute log stable sum and return
	return log_stable_sum(log_res,TRUE,max,idx_max);
}

// ---- K post
double p_distinct_post_largen(	const unsigned int& r, const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const std::vector<double>& v_log_qM_post )
{
	double res{0.0};

	const unsigned int d = n_j.size();

	double log_fact_mom_r  = log_factorial_mom(r,v_log_qM_post);
	double log_fact_mom_r1 = log_factorial_mom(r+1,v_log_qM_post);
	//std::vector<double> A_j(d,0.0);
	std::vector<double> C_j(d,0.0);
	double C{0.0};
	for(std::size_t j = 0; j < d; j++){
		C_j[j] = gamma_j[j] * (double)m_j[j] / ((double)n_j[j]);
		C += C_j[j];
	}

	res = std::exp( log_fact_mom_r  + std::log(1.0 + 5.0/3.0 * (double)r * C ) ) -
		  std::exp( log_fact_mom_r1 + std::log(C) );
	if(res < 0){
		throw std::runtime_error("Error in p_distinct_post_largen: Probability is less than 0. ");
	}
	res *= std::exp( (double)r*std::log(3.0) - gsl_sf_lnfact((double)r) );
	return res;
}

double p_distinct_post_largen(	const unsigned int& r, const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const Rcpp::String& prior, const Rcpp::List& prior_param,
						      	double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return p_distinct_post_largen(r,k,m_j,n_j,gamma_j,v_log_qM_post);
}

Rcpp::NumericVector D_distinct_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const std::vector<double>& v_log_qM_post )
{
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);

	Rcpp::NumericVector res(m+1,0.0);
	double Pcum{0.0};

	for(std::size_t r=0; r <= m; r++){
		res[r] = p_distinct_post_largen(r, k, m_j, n_j, gamma_j, v_log_qM_post );
		Pcum += res[r];
		if(1.0 - Pcum < 1e-10)
			break;
	}
	return res;
}

Rcpp::NumericVector D_distinct_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const Rcpp::String& prior, const Rcpp::List& prior_param,
									      	double log_V, unsigned int M_max )
{

	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return D_distinct_post_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}


// ---- (K,S) post
double p_jointKS_post_largen(	const unsigned int& r, const unsigned int& t,
								const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const std::vector<double>& v_log_qM_post )
{
	if(t > r)
		return 0.0; // it is not possible to have more shared than distinct species

	double res{0.0};
	const unsigned int d = n_j.size();

	double log_fact_mom_r  = log_factorial_mom(r,v_log_qM_post);
	double log_fact_mom_r1 = log_factorial_mom(r+1,v_log_qM_post);
	//std::vector<double> A_j(d,0.0);
	std::vector<double> C_j(d,0.0);
	double C{0.0};
	for(std::size_t j = 0; j < d; j++){
		C_j[j] = gamma_j[j] * (double)m_j[j] / ((double)n_j[j]);
		C += C_j[j];
	}

	res = std::exp( -std::log(2) + log_fact_mom_r  + std::log(2.0 + (3.0 * (double)r +(double)t)*C ) ) -
		  std::exp( log_fact_mom_r1 + std::log(C) );
	if(res < 0){
		throw std::runtime_error("Error in p_jointKS_post_largen: Probability is less than 0. ");
	}
	res *= std::exp( (double)(r-t)*std::log(2.0) - gsl_sf_lnfact((double)t) - gsl_sf_lnfact((double)(r-t)) );
	return res;
}

double p_jointKS_post_largen(	const unsigned int& r, const unsigned int& t,
								const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const Rcpp::String& prior, const Rcpp::List& prior_param,
						      	double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return p_jointKS_post_largen(r,t,k,m_j,n_j,gamma_j,v_log_qM_post);
}

Rcpp::NumericMatrix D_jointKS_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const std::vector<double>& v_log_qM_post )
{
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);

	Rcpp::NumericMatrix res(m+1,m+1);
	double jointPcum{0.0};
	bool exit = FALSE;
	for(std::size_t r=0; r <= m; r++){

		if(exit)
			break;

		for(std::size_t t=0; t <= r; t++){
			res(t,r) = p_jointKS_post_largen(r, t, k, m_j, n_j, gamma_j, v_log_qM_post );

			jointPcum += res(t,r);
			if(1.0 - jointPcum < 1e-10)
				exit = TRUE;

			if(exit)
				break;

		}


	}

	return res;
}

Rcpp::NumericMatrix D_jointKS_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const Rcpp::String& prior, const Rcpp::List& prior_param,
									      	double log_V, unsigned int M_max )
{

	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return D_jointKS_post_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}

// ---- S post
double p_shared_post_largen(	const unsigned int& t, const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const std::vector<double>& v_log_qM_post )
{
	const double inf{std::numeric_limits<double>::infinity()};
	const unsigned int d = n_j.size();
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);
	const unsigned int Natoms = v_log_qM_post.size(); // compute number of atoms in qM_post whose sum is enough to reach 1

	if(t > m)
		return 0.0; // it is not possible to have more shared species than observations
	if(Natoms == 0)
		throw std::runtime_error("Error in p_shared_post_largen: v_log_qM_post is empty");

	if(t >= Natoms)
		return 0.0; // if here, Expected value is 0 beacuse q_M has no mass after mstar = t

	std::vector<double> log_res(Natoms-t,-inf); // initialize vector used to compute the log stable sum
	int counter{0}; // initialize counter
	double max{-inf}; // initialize max value in log_res
	unsigned int idx_max{0}; // initialize position of max element in log_res

	// Define C_j for each j = 1,...,d
	std::vector<double> C_j(d,0.0);
	double C{0.0};
	for(std::size_t j = 0; j < d; j++){
		C_j[j] = gamma_j[j] * (double)m_j[j] / ((double)n_j[j]);
		C += C_j[j];
	}

	// start loop over all atomes
	double temp{0.0};
	for(std::size_t mstar = t; mstar < Natoms; mstar++){
		log_res[counter] = v_log_qM_post[mstar] + (double)(mstar - t)*std::log(3.0) + gsl_sf_lnchoose(mstar,t);
		temp = (1.0 + (double)t*C)*gsl_cdf_binomial_P(m-t, 2.0/3.0, mstar-t);

		if( (mstar > t) & (mstar >= m) ){
			temp -= ((double)(mstar-t) * C)/(3.0) * gsl_ran_binomial_pdf(m-t, 2.0/3.0, mstar-t-1);
			if(temp < 0)
				throw std::runtime_error("Error in p_shared_post_largen: Probability is less than 0.");
		}
		log_res[counter] += std::log(temp);

		// check for max
		if(log_res[counter] > max){
			max = log_res[counter];
			idx_max = counter;
		}
	}

	// compute log stable sum and return
	return std::exp(log_stable_sum(log_res,TRUE,max,idx_max));
}

double p_shared_post_largen(	const unsigned int& t, const unsigned int& k,
								const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						      	const std::vector<double>& gamma_j,
						      	const Rcpp::String& prior, const Rcpp::List& prior_param,
						      	double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return p_shared_post_largen(t,k,m_j,n_j,gamma_j,v_log_qM_post);
}

Rcpp::NumericVector D_shared_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const std::vector<double>& v_log_qM_post )
{
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);

	Rcpp::NumericVector res(m+1,0.0);
	double Pcum{0.0};

	for(std::size_t t=0; t <= m; t++){
		res[t] = p_shared_post_largen(t, k, m_j, n_j, gamma_j, v_log_qM_post );
		Pcum += res[t];
		if(1.0 - Pcum < 1e-10)
			break;
	}
	return res;
}

Rcpp::NumericVector D_shared_post_largen(	const unsigned int& k,
											const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
									      	const std::vector<double>& gamma_j,
									      	const Rcpp::String& prior, const Rcpp::List& prior_param,
									      	double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return D_shared_post_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}

// ---- K post - moments
double Kpost_mom1_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						   	const std::vector<double>& gamma_j,
						   	const std::vector<double>& v_log_qM_post )
{
	const double inf{std::numeric_limits<double>::infinity()};

	const unsigned int d = n_j.size();
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);
	const unsigned int Natoms = v_log_qM_post.size(); // compute number of atoms in qM_post whose sum is enough to reach 1
	if(Natoms == 0)
		throw std::runtime_error("Error in Kpost_mom1_largen: v_log_qM_post is empty");

	std::vector<double> log_res(Natoms-1,-inf); // initialize vector used to compute the log stable sum
	int counter{0};   // initialize counter
	double max{-inf}; // initialize max value in log_res
	unsigned int idx_max{0}; // initialize position of max element in log_res

	// Define C_j for each j = 1,...,d
	std::vector<double> C_j(d,0.0);
	double C{0.0};
	for(std::size_t j = 0; j < d; j++){
		C_j[j] = gamma_j[j] * (double)m_j[j] / ((double)n_j[j]);
		C += C_j[j];
	}

	// start loop over all atomes
	double temp{0.0};
	for(std::size_t mstar = 1; mstar < Natoms; mstar++){
		log_res[counter] = v_log_qM_post[mstar] + (double)(mstar - 1)*std::log(4.0) + std::log(mstar);
		temp = (3.0 + 5.0*C)*gsl_cdf_binomial_P(m-1, 3.0/4.0, mstar-1);

		if( mstar > 1 ){
			temp += ( 3.0 * (double)(mstar-1) * C )/(4.0) *
					( 5.0 * gsl_cdf_binomial_P(m-2, 3.0/4.0, mstar-2) - gsl_cdf_binomial_P(m-1, 3.0/4.0, mstar-2) );
		}
		log_res[counter] += std::log(temp);

		// check for max
		if(log_res[counter] > max){
			max = log_res[counter];
			idx_max = counter;
		}
	}

	// compute log stable sum and return
	return std::exp(log_stable_sum(log_res,TRUE,max,idx_max));
}

double Kpost_mom1_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						    const std::vector<double>& gamma_j,
						    const Rcpp::String& prior, const Rcpp::List& prior_param,
						    double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return Kpost_mom1_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}

double Kpost_mom2_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						   	const std::vector<double>& gamma_j,
						   	const std::vector<double>& v_log_qM_post )
{
	const double inf{std::numeric_limits<double>::infinity()};

	const unsigned int d = n_j.size();
	const unsigned int m = std::accumulate(m_j.cbegin(), m_j.cend(), 0);
	const unsigned int Natoms = v_log_qM_post.size(); // compute number of atoms in qM_post whose sum is enough to reach 1
	if(Natoms == 0)
		throw std::runtime_error("Error in Kpost_mom2_largen: v_log_qM_post is empty");

	std::vector<double> log_res(Natoms-1,-inf); // initialize vector used to compute the log stable sum
	int counter{0};   // initialize counter
	double max{-inf}; // initialize max value in log_res
	unsigned int idx_max{0}; // initialize position of max element in log_res

	// Define C_j for each j = 1,...,d
	std::vector<double> C_j(d,0.0);
	double C{0.0};
	for(std::size_t j = 0; j < d; j++){
		C_j[j] = gamma_j[j] * (double)m_j[j] / ((double)n_j[j]);
		C += C_j[j];
	}

	// start loop over all atomes
	double temp{0.0};
	for(std::size_t mstar = 1; mstar < Natoms; mstar++){
		log_res[counter] = v_log_qM_post[mstar] + (double)(mstar - 1)*std::log(4.0) + std::log(mstar);
		temp = (3.0 + 5.0*C)*gsl_cdf_binomial_P(m-1, 3.0/4.0, mstar-1);

		if( mstar > 1 ){
			temp += ( 9.0 * (double)(mstar-1) )/(4.0) *
					( (1.0 + 5.0 * C) * gsl_cdf_binomial_P(m-2, 3.0/4.0, mstar-2) -
					  (1.0 / 3.0 * C  * gsl_cdf_binomial_P(m-1, 3.0/4.0, mstar-2) )
					);
		}
		if( mstar > 2 ){
			temp += ( 9.0 * (double)(mstar-1) * (double)(mstar-2) * C )/(16.0) *
					( 5.0 * gsl_cdf_binomial_P(m-3, 3.0/4.0, mstar-3) - gsl_cdf_binomial_P(m-1, 3.0/4.0, mstar-3) );
		}
		log_res[counter] += std::log(temp);

		// check for max
		if(log_res[counter] > max){
			max = log_res[counter];
			idx_max = counter;
		}
	}

	// compute log stable sum and return
	return std::exp(log_stable_sum(log_res,TRUE,max,idx_max));
}

double Kpost_mom2_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						    const std::vector<double>& gamma_j,
						    const Rcpp::String& prior, const Rcpp::List& prior_param,
						    double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return Kpost_mom2_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}

double Kpost_var_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						   	const std::vector<double>& gamma_j,
						   	const std::vector<double>& v_log_qM_post )
{
	double Exp_val = Kpost_mom1_largen(	k, m_j, n_j, gamma_j, v_log_qM_post );
	double Mom2    = Kpost_mom2_largen(	k, m_j, n_j, gamma_j, v_log_qM_post );
	return Mom2 - (Exp_val*Exp_val);
}

double Kpost_var_largen(	const unsigned int& k,
							const std::vector<unsigned int>& m_j, const std::vector<unsigned int>& n_j,
						    const std::vector<double>& gamma_j,
						    const Rcpp::String& prior, const Rcpp::List& prior_param,
						    double log_V, unsigned int M_max )
{
	const double inf{std::numeric_limits<double>::infinity()};
	// Component prior preliminary operations
	auto qM_ptr = Wrapper_ComponentPrior(prior, prior_param);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// The correctness of log_V is not checked. If it is -inf, then log_V is computed
	if(log_V == -inf)
		log_V = compute_log_Vprior(k, n_j, gamma_j, qM, M_max );

	std::vector<double> v_log_qM_post = D_log_qM_post(qM, k, n_j, gamma_j, log_V, M_max );

	return Kpost_var_largen(k,m_j,n_j,gamma_j,v_log_qM_post);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Tests
//------------------------------------------------------------------------------------------------------------------------------------------------------

// NON FUNZIONA, non usare
Rcpp::List Test_Riemann(double s, int n)
{
	double log_z = std::log(gsl_sf_zeta(s));
	std::vector<double> uno_N(n);
	std::iota(uno_N.begin(),uno_N.end(),1);

	Rcpp::Rcout<<"Stampo uno_N: ";
	for(auto __v : uno_N)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	// Per qualche motivo, non funziona!!
	std::for_each(uno_N.begin(),uno_N.end(), [&s](double h){
		Rcpp::Rcout<<"1.0/std::pow( (double)h,s):"<<std::endl<<1.0/std::pow( (double)h,s)<<std::endl;
		return (1.0/std::pow( (double)h,s) ) ; } );

	Rcpp::Rcout<<"Stampo uno_N: ";
	for(auto __v : uno_N)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;
	double log_resto = log_stable_sum(uno_N,TRUE);
	double log_res = log_z + std::log(1-std::exp(log_resto-log_z) );
	return Rcpp::List::create( Rcpp::Named("log_z") = log_z,
							   Rcpp::Named("log_resto") = log_resto,
							   Rcpp::Named("log_res") = log_res
							   );
}

void Test_Prior(){
	Rcpp::String form = "Poisson";
	double lambda(5.0);
	ComponentPrior_Parameters param;
	param.Lambda = lambda;
	auto qM_ptr = Select_ComponentPrior(form, param);
	ComponentPrior& qM(*qM_ptr);
	std::string i_am = qM_ptr->showMe();
	std::string i_am2 = qM.showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	Rcpp::Rcout<<"I am - II - : "<<i_am2<<std::endl;
	form = "NegativeBinomial";
	param.p = 0.45;
	param.n_succ = 2;
	auto qM_ptr2 = Select_ComponentPrior(form, param);
	i_am = qM_ptr2->showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	// Devo testare se le densità sono giuste!
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"============== Densità ============="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"Poisson1(0)="<<qM_ptr->eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"Poisson1(1)="<<qM_ptr->eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"Poisson1(2)="<<qM_ptr->eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"Poisson1(3)="<<qM_ptr->eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"Poisson1(4)="<<qM_ptr->eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(0)="<<qM_ptr2->eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(1)="<<qM_ptr2->eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(2)="<<qM_ptr2->eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial(3)="<<qM_ptr2->eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial(4)="<<qM_ptr2->eval_prob(4)<<std::endl;
	// Devo testare se le densità sono giuste!
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"============ Log Densità ==========="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(0)="<<qM_ptr->log_eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(1)="<<qM_ptr->log_eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(2)="<<qM_ptr->log_eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(3)="<<qM_ptr->log_eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(4)="<<qM_ptr->log_eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(0)="<<qM_ptr2->log_eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(1)="<<qM_ptr2->log_eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(2)="<<qM_ptr2->log_eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(3)="<<qM_ptr2->log_eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(4)="<<qM_ptr2->log_eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"================ Moda =============="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"La moda Poisson1 quando lambda = "<<lambda<<" è pari a "<<qM_ptr->get_mode()<<std::endl;
	Rcpp::Rcout<<"La moda NegativeBinomial1 quando p = "<<lambda<<" è pari a "<<param.p<<" e n è "<<param.n_succ<<" è pari a "<<qM_ptr2->get_mode()<<std::endl;
}


void Test_prod_sum(){

	double inf = std::numeric_limits<double>::infinity();

	std::vector<unsigned int> n{1,1,1};
	std::vector<double> gamma{1.0,2.0,3.0};
	unsigned int M{2};
	unsigned int k{0};

	std::vector<double> log_a{1.0, 2.0, 1.0};
	std::vector<double> a{ std::exp(1.0), std::exp(2.0), std::exp(1.0)};
	unsigned int idx{1};
	double max_log{log_a[idx]};
	double max_nolog{a[idx]};

	double res1 = combined_product(n,  gamma,  M,  k);
	double res2 = combined_sum(n,  gamma,  M,  k);

	// Test log_stable_sum
	double res3 = log_stable_sum(log_a, TRUE, max_log, idx);
	double res4 = log_stable_sum(a, FALSE, max_nolog, idx);

	double res5 = log_stable_sum(log_a, TRUE);
	double res6 = log_stable_sum(a, FALSE);

	Rcpp::Rcout<<"res combined_product = "<<res1<<std::endl;
	Rcpp::Rcout<<"res combined_sum = "<<res2<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (log and max) = "<<res3<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (no log and max) = "<<res4<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (log) = "<<res5<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum ( no log) = "<<res6<<std::endl;

	Rcpp::Rcout<<"--------------------------------------------"<<std::endl;
	std::vector<double> a_zero{ 0.0, 0.0, 0.0 };
	std::vector<double> log_a_inf{ -inf, -inf, -inf };
	Rcpp::Rcout<<"res log_stable_sum( a_zero ) = "<<log_stable_sum(a_zero, FALSE, 0.0, 0)<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum( log_a_inf ) = "<<log_stable_sum(log_a_inf, TRUE, -inf, 0)<<std::endl;
}
