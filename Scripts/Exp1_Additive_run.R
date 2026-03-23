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


# Setting -----------------------------------------------------------------
seed0 = 271296
set.seed(seed0)
num_cores = 34

## Dirichlet parameters
name_exp = "Dirichelt"
p1 = 0.1; p2 = 0.5
params = c(p1,p2)
masses = c(0,0.5,1)
params_grid = expand.grid(params,params,masses)


## Options
Q  = 80
Nrep <- BB <- 100
m_j <- c(200,200)
nmax <- 800
ngrid <- seq(100,nmax,by = 100)
Lgrid = length(ngrid)
n1_prop = 1
shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
save_all = TRUE
pooled_indep = TRUE
run_chao = FALSE

nome_base1 = paste0("save/Pred_Additive_")

# Run ---------------------------------------------------------------------

idx = 2
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
    
  # Run
  results = lapply( X = ngrid, FUN = Prediction_fit,
                    n1_prop = n1_prop,
                    probs_1 = prob_true_mat1, probs_2 = prob_true_mat2,
                    m_j = m_j,
                    seed0 = seed0, num_cores = num_cores, 
                    par_est_ind = par_est_ind, par_est_mle = par_est_mle,
                    pooled_indep = pooled_indep, run_chao = run_chao)
    
  save_name = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat") 
    
  if(save_all)
    save(results,file = save_name)
    
  }
  


























