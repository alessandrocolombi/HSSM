#include "Sampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>

// Consrtructor
Sampler::Sampler( const GS_data & __gs_data, unsigned int seed ):gs_data(__gs_data), random_engine(seed) {

	//Initialize Full Conditional Objects
  auto gamma_ptr  = std::make_shared<FC_gamma>("gamma",gs_data.FixedGammas);
  auto beta_ptr   = std::make_shared<FC_beta>("beta",gs_data.FixedBeta);
  auto Lambda_ptr = std::make_shared<FC_Lambda>("Lambda", gs_data.FixedLambda);

  //Full Conditional vector that we will loop over
  std::vector< std::shared_ptr<FullConditional> > fc{gamma_ptr, beta_ptr,Lambda_ptr};
  std::swap(FullConditionals, fc);

  // Initialize return structure
  // burn_in + n_iter*thin should be just n_iter, that is the number of values I want as output
    	//out.S.reserve(burn_in + n_iter*thin);
    	//out.mu.reserve(burn_in + n_iter*thin);
    	//out.sigma.reserve(burn_in + n_iter*thin);
}

void Sampler::sample() {

    Progress progress_bar(gs_data.burn_in + gs_data.n_iter*gs_data.thin, TRUE); // Initialize progress bar
    for(unsigned int it = 0; it < gs_data.burn_in + gs_data.n_iter*gs_data.thin; it++){

        // If we are in the right iteration store needed values
        // If burn_in is set to zero, the first values to be saved are the initial values.
        if(it>=gs_data.burn_in && (it-gs_data.burn_in)%gs_data.thin == 0){

                    //Rcpp::Rcout<<"it = "<<it<<std::endl;
                    //Rcpp::Rcout<<"Stampo gs_data.gamma_j: ";        
                    //for(auto __v : gs_data.gamma_j)
                        //Rcpp::Rcout<<__v<<", ";
                    //Rcpp::Rcout<<std::endl;
                    //Rcpp::Rcout<<"gs_data.Lambda = "<<gs_data.Lambda<<std::endl;
                    //Rcpp::Rcout<<"gs_data.logV = "<<gs_data.logV<<std::endl;

            out.gammas.push_back( Rcpp::NumericVector (gs_data.gamma_j.begin(),gs_data.gamma_j.end()) );  //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
            out.betas.push_back(gs_data.beta);
            out.Lambdas.push_back(gs_data.Lambda);
            out.logV.push_back(gs_data.logV);
        }

        // Sample from all full conditional
        for(auto full_cond: FullConditionals){

             //Rcpp::Rcout<<"-------------------------------------------"<<std::endl;
            if(!full_cond->keep_fixed){
                full_cond->update(gs_data, random_engine);
                //Rcpp::Rcout<<" --> done! "<<std::endl;
            }
            // Update logV
            Rcpp::List prior_param = 	Rcpp::List::create( Rcpp::Named("lambda") = gs_data.Lambda);
            auto qM_ptr = Wrapper_ComponentPrior("Poisson", prior_param);
            ComponentPrior& qM(*qM_ptr);
            gs_data.logV = compute_log_Vprior(gs_data.Kobs, gs_data.n_j, gs_data.gamma_j, qM, gs_data.M_max);
            //Rcpp::Rcout<<"********************************************"<<std::endl;

            //Check for User Interruption
            try{
                Rcpp::checkUserInterrupt();
            }
            catch(Rcpp::internal::InterruptedException e){
                //Print error and return
                throw std::runtime_error("Execution stopped by the user");
            }
        }

        //updating number of iterations necessary for MH algorithm
        gs_data.iterations = it;

        progress_bar.increment(); //update progress bar
    }
}

Rcpp::List MCMC_Sampler_c(const Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& N__jk,
                          const std::vector<unsigned int>& n__j,
                          unsigned int r,
                          unsigned int niter, unsigned int nburn, unsigned int thin,
                          const Rcpp::List& option)
{
  GS_data gs_data;
  // use covariates (model selection)
  gs_data.use_covariates = Rcpp::as<bool>(option["use_covariates"]);
  // data
  gs_data.N_jk = N__jk;
  gs_data.X_j = Rcpp::as<std::vector<double>>(option["X_j"]) ;
  gs_data.n_j = n__j;
  gs_data.Kobs = r;
  gs_data.d = gs_data.n_j.size();
  // basic options
  gs_data.n_iter = niter;
  gs_data.burn_in = nburn;
  gs_data.thin = thin;
  gs_data.iterations = 0;
  // initial values
  gs_data.gamma_j = Rcpp::as<std::vector<double>>(option["gamma0"]) ;
  gs_data.beta = Rcpp::as<double>(option["beta0"]);
  gs_data.Lambda = Rcpp::as<double>(option["Lambda0"]);

  Rcpp::List prior_param = 	Rcpp::List::create( Rcpp::Named("lambda") = gs_data.Lambda);
  auto qM_ptr = Wrapper_ComponentPrior("Poisson", prior_param);
  ComponentPrior& qM(*qM_ptr);
  gs_data.M_max = Rcpp::as<double>(option["M_max"]);
  gs_data.logV = compute_log_Vprior(gs_data.Kobs, gs_data.n_j, gs_data.gamma_j, qM, gs_data.M_max);


  // Fixed quantities
  gs_data.FixedLambda = !Rcpp::as<bool>(option["UpdateLambda"]);
  gs_data.FixedGammas = !Rcpp::as<bool>(option["UpdateGamma"]);
  gs_data.FixedBeta   = !Rcpp::as<bool>(option["UpdateBeta"]);
  // hyperparam
  gs_data.Lambda0     = Rcpp::as<double>(option["L0"]);
  gs_data.V_Lambda    = Rcpp::as<double>(option["V_Lambda"]);
  gs_data.sigma2      = Rcpp::as<double>(option["sigma2"]);
  gs_data.a_gamma     = Rcpp::as<double>(option["a_gamma"]);
  gs_data.b_gamma     = Rcpp::as<std::vector<double>>(option["b_gamma"]);
  gs_data.sigma2_beta = Rcpp::as<double>(option["sigma2_beta"]);
  // MCMC options
  gs_data.adapt_var_gamma_j = Rcpp::as<std::vector<double>>(option["adapt_var_gamma_j"]);
  gs_data.adapt_var_beta = Rcpp::as<double>(option["adapt_var_beta"]);



  Sampler sampler(gs_data, Rcpp::as<unsigned int>(option["seed"]));
  sampler.sample();
  Rcpp::List res = Rcpp::List::create( Rcpp::Named("gammas") = sampler.out.gammas,
                                     Rcpp::Named("betas") = sampler.out.betas,
                                     Rcpp::Named("Lambdas") = sampler.out.Lambdas,
                                     Rcpp::Named("logV") = sampler.out.logV
                                    );
  return res;
}
