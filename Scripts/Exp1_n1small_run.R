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

experiments = list(
  "Dirichelt" = c(0.1,0.5),
  "Zipfs" = c(1.3,2),
  "Geom"  = c(0.8,0.85,0.9)
)

M = 60
Nrep <- BB <- 100
m_j <- c(200,200)
# n2max <- 1000
# n2grid <- seq(100,n2max,by = 100)
n2grid = c(1000,2000)
n1_prop <- 1/10
Lgrid = length(n2grid)
shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
save_all = TRUE
pooled_indep = TRUE

nome_base1 = paste0("save/Pred_n1verysmall_")

# Run ---------------------------------------------------------------------

igrid = c(1:3)
i = 1
for(i in igrid){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
  # params_grid = expand.grid(params,params)
  params_grid <- rbind(
    cbind(params, params),
    t(combn(params, 2))
  )
  
  
  idx = 3
  for(idx in 1:nrow(params_grid) ){
    p1 = params_grid[idx,1]; p2 = params_grid[idx,2]
    cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2,"\n")
    
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
    
    # Run
    results = lapply( X = n2grid, FUN = Prediction_fit,
                      n1_prop = n1_prop,
                      probs_1 = prob_true_mat1, probs_2 = prob_true_mat2,
                      m_j = m_j,
                      seed0 = seed0, num_cores = num_cores, 
                      pooled_indep = pooled_indep,
                      par_est_ind = par_est_ind, par_est_mle = par_est_mle)
    
    save_name = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat") 
    
    if(save_all)
      save(results,file = save_name)
    
  }
  
}


























