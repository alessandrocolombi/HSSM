wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"

wd = paste0(wd_g100,"ScriptSpecies_shared/RevBA")
setwd(wd)
# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source(paste0(wd_g100,"ScriptSpecies_shared/SimStudyRaf/SimStudyRaf_functions.R"))

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))

# Colors ------------------------------------------------------------------
mycol = hcl.colors(n = 5, palette = "viridis")

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
n2grid = seq(100,500,by = 50)
n1_prop = 1/5
Lgrid = length(n2grid)
shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
save_all = TRUE
pooled_indep = TRUE

nome_base1 = paste0("save/Pr1step_Add_n1small_")

# Run ---------------------------------------------------------------------

idx = 5
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
  
  ## Run
  Pr1step_ngrid = Pr1step_Nrep_fit(prob_true_mat1, prob_true_mat2, n2grid, 
                                   n1_prop = n1_prop,
                                   Nrep, num_cores, seed0, names_species,
                                   pooled_indep = pooled_indep,
                                   par_est_ind, par_est_mle)
  
  save_name = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat") 
  
  if(save_all)
    save(Pr1step_ngrid,file = save_name)
  
}



