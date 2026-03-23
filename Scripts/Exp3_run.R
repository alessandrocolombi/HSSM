wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[3] # <--- modify here
wd = paste0(choose_wd,"ScriptSpecies_shared/RevBA/")
setwd(wd)

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))

# Colors ------------------------------------------------------------------
mycol = hcl.colors(n = 5, palette = "viridis")

# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source("../SimStudyRaf/SimStudyRaf_functions.R")


# Options global ----------------------------------------------------------

seed0 = 271296
set.seed(seed0)
num_cores = 34 # <---
Nrep <- BB <- 100 # <---

n2grid = seq(50,400,by = 50) # <---
Lngrid = length(n2grid)   
n1_prop = 1
m_max = 200   # <---
mgrid = c(1,seq(10,m_max, by = 10)) # <---
Lmgrid = length(mgrid)

# (D,Z,G) -----------------------------------------------------------------


experiments = list(
  "Dirichelt" = c(0.5,0.5),
  "Zipfs" = c(2,2),
  "Geom"  = c(0.9,0.9)
)

M = 60
shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
pooled_indep = FALSE

save_all = TRUE

nome_base1 = paste0("save/Prob_msteps_")
# Assume that there are no area-specific species
names_shared <- names_species <- sapply(1:M, function(x){paste0("sh",as.character(x))})

## Run ---------------------------------------------------------------------
run_DZG = TRUE

if(run_DZG){
  igrid = c(1:3)
  i = 1
  for(i in igrid){
    name_exp = names(experiments)[i]
    cat("\n Start: ",name_exp,"\n")
    param = experiments[[i]]
    p1 = param[1]; p2 = param[2]
    
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
        
      # Generate true probabilities
      prob_true_mat1 <- prob_true_mat2 <- matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
      prob_true_mat1 = t(apply(prob_true_mat1, 1, function(x) { 
          pj = sim_generic(name_exp,M,p1)
          if(shuffle)
            pj = pj[sample(1:M,M)]
          pj  } ))  # (BB x M) matrix
      prob_true_mat2 = t(apply(prob_true_mat2, 1, function(x) { 
          pj = sim_generic(name_exp,M,p2)
          if(shuffle)
            pj = pj[sample(1:M,M)]
          pj } ))  # (BB x M) matrix
        
      ## Parallel run (no prints allowed)
      seeds = sample(1:999999, size = Nrep)
      idx_rep = 1:Nrep
      cluster <- makeCluster(num_cores, type = "SOCK")
      doSNOW::registerDoSNOW(cluster)
      clusterExport(cluster, list("ParEst_indicies", "ParEst_MLE", "ParEst_MLE_d1",
                                  "llik_VecFDP", "llik_FDP",
                                  "GiniSimpson_est","Mor_est",
                                  "compute_logC","PrSh0_c","PrDistinct0_c","PrDistinct0_d1_c"),
                    envir = environment())
      inner_result = parLapply( cl = cluster, x = idx_rep,
                                fun = Pr_msteps_fit_run,
                                n2grid = n2grid, names_species = names_species,
                                mgrid = mgrid,
                                probs_1 = prob_true_mat1, probs_2 = prob_true_mat2, seeds = seeds,
                                par_est_ind = par_est_ind, par_est_mle = par_est_mle,
                                pooled_indep = pooled_indep)
      stopCluster(cluster)
      
      ## REMARK: Here, Pr(S == 0|X),Pr(K == 0|X)... are saved
      save_name = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat") 
        
      if(save_all)
        save(inner_result,file = save_name)
    
  }
}

# Additive ----------------------------------------------------------------
set.seed(seed0)

## Dirichlet parameters
name_exp = "Dirichelt"
p1 = 0.5; p2 = 0.5
params = c(p1,p2)
masses = c(0,0.5,1)
params_grid = expand.grid(params,params,masses)
params_grid = params_grid[c(4,8,12),] # <--- select here the experiments to run

## Options
Q  = 80

shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
save_all = TRUE
pooled_indep = FALSE

nome_base1 = paste0("save/Prob_msteps_Additive_")


## Run ---------------------------------------------------------------------
run_add = TRUE
if(run_add){
  idx = 1
  Nrun = nrow(params_grid)
  for(idx in 1:Nrun ){
    ## Read parameters
    p1 = params_grid[idx,1]; p2 = params_grid[idx,2]; mass = params_grid[idx,3]
    cat("\n ",idx,"/",nrow(params_grid)," - p1 = ",p1," vs p2 = ",p2," || ","mass = ",mass,"\n")
    trim_p1 = get_first3digits(p1,4)
    trim_p2 = get_first3digits(p2,4)
    trim_mass = get_first3digits(mass,4)
    
    ## Set sizes (M,M1,M2)
    M = floor(mass*Q); M1 <- M2 <- floor( Q*(1-mass)/2 )
    if( (M+M1+M2) != Q )
      M = M + (Q-M-M1-M2)
    
    Q1 = M+M1; Q2 = M+M2
    
    # Species names
    names_shared <- names_distinct1 <- names_distinct2 <- c()
    if(M>0){
      names_shared <- sapply(1:M, function(x){paste0("sh",as.character(x))})
    }
    if(M1 > 0){
      names_distinct1 <- sapply(1:M1, function(x){paste0("d1_",as.character(x))})
    }
    if(M2 > 0){
      names_distinct2 <- sapply(1:M2, function(x){paste0("d2_",as.character(x))})
    }
    names_species <- c(names_shared,names_distinct1,names_distinct2)
    
    
    # Generate probabilities
    p_common <- matrix(0,nrow = BB, ncol = M)  # (BB x M) matrix
    p_only_1 <- matrix(0,nrow = BB, ncol = M1) # (BB x M1) matrix
    p_only_2 <- matrix(0,nrow = BB, ncol = M2) # (BB x M2) matrix
    p_common = t(apply(p_common, 1, function(x) { 
      pj = sim_generic(name_exp,M,p1)
      if(shuffle)
        pj = pj[sample(1:M,M)]
      pj = pj * mass
      pj  } ))
    p_only_1 = t(apply(p_only_1, 1, function(x) { 
      pj = sim_generic(name_exp,M1,p2)
      if(shuffle)
        pj = pj[sample(1:M1,M1)]
      pj = pj * (1-mass)
      pj  } ))
    p_only_2 = t(apply(p_only_2, 1, function(x) { 
      pj = sim_generic(name_exp,M2,p2)
      if(shuffle)
        pj = pj[sample(1:M2,M2)]
      pj = pj * (1-mass)
      pj  } ))
    
    
    prob_true_mat1 <- prob_true_mat2 <- matrix(0,nrow = BB, ncol = Q) # (BB x Q) matrix
    if(M > 0){
      prob_true_mat1[,1:M] <- prob_true_mat2[,1:M] <- p_common
    }
    if(M < Q){
      prob_true_mat1[,(M+1):(M+M1)] <- p_only_1
      prob_true_mat2[,(Q1+1):(Q)] <- p_only_2
    }
    
    ## Parallel run (no prints allowed)
    seeds = sample(1:999999, size = Nrep)
    idx_rep = 1:Nrep
    cluster <- makeCluster(num_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)
    clusterExport(cluster, list("ParEst_indicies", "ParEst_MLE", "ParEst_MLE_d1",
                                "llik_VecFDP", "llik_FDP",
                                "GiniSimpson_est","Mor_est",
                                "compute_logC","PrSh0_c","PrDistinct0_c","PrDistinct0_d1_c"),
                  envir = environment())
    inner_result = parLapply( cl = cluster, x = idx_rep,
                              fun = Pr_msteps_fit_run,
                              n2grid = n2grid, names_species = names_species,
                              mgrid = mgrid,
                              probs_1 = prob_true_mat1, probs_2 = prob_true_mat2, 
                              seeds = seeds,
                              par_est_ind = par_est_ind, par_est_mle = par_est_mle,
                              pooled_indep = pooled_indep)
    stopCluster(cluster)
    
    ## REMARK: Here, Pr(S == 0|X),Pr(K == 0|X)... are saved
    save_name = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat") 
    
    if(save_all)
      save(inner_result,file = save_name)
  }
  
}
