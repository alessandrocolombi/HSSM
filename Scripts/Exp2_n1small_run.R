wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"

wd = paste0(wd_g100,"ScriptSpecies_shared/RevBA")
setwd(wd)
# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source(paste0(wd_g100,"ScriptSpecies_shared/SimStudyRaf/SimStudyRaf_functions.R"))

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
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

experiments = list(
  "Dirichelt" = c(0.1,0.5),
  "Zipfs" = c(1.3,2),
  "Geom"  = c(0.8,0.85,0.9)
)

M = 60
d = 2
Nrep <- BB <- 100
n1_prop = 1/5
shuffle = TRUE
par_est_ind = TRUE
par_est_mle = TRUE
pooled_indep = TRUE
save_all = TRUE

nome_base1 = paste0("save/Pr1step_n1small_")
n2grid = seq(100,500,by = 50)

# Assume that there are no area-specific species
names_shared <- names_species <- sapply(1:M, function(x){paste0("sh",as.character(x))})

# Run ---------------------------------------------------------------------

igrid = c(1:3)
i = 1
for(i in igrid){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
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
    
    ## Run
    Pr1step_ngrid = Pr1step_Nrep_fit(prob_true_mat1, prob_true_mat2, n2grid,
                                     n1_prop = n1_prop,
                                     Nrep, num_cores, seed0, names_species,
                                     pooled_indep = pooled_indep,
                                     par_est_ind, par_est_mle)
    
    name_save = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat") 
    if(save_all)
      save(Pr1step_ngrid,file = name_save)
  }
  
}


