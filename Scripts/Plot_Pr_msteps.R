wd_pc = "C:/Users/colom/"
wd = paste0(wd_pc,"ScriptSpecies_shared/RevBA")
setwd(wd)

# Libraries ---------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))


# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source("../SimStudyRaf/SimStudyRaf_functions.R")


# Plot options ------------------------------------------------------------

save_img = FALSE
width = 8; height = 6
cex.lab = 2
cex.axis = 2
cex = 2

# Options global ----------------------------------------------------------

seed0 = 271296
set.seed(seed0)
Nrep <- BB <- 100 # <---

n2grid = seq(50,400,by = 50) # <---
Lngrid = length(n2grid)   
n1_prop = 1
m_max = 200   # <---
mgrid = c(1,seq(10,m_max, by = 10)) # <---
Lmgrid = length(mgrid)

mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
# nomi_exp = 2*n2grid
# pos  = seq(1,length(nomi_exp),by = 1)
xpos = n2grid
xlabs = n2grid*2
idxs = c(1,4,7,10,14)
mgrid_pos = round(seq(min(mgrid),max(mgrid), length.out = 5))

img_folder = "img/Pr_m_steps/"#"img/Paper/"



# (D,Z,G) -----------------------------------------------------------------

experiments = list(
  "Dirichelt" = c(0.5,0.5),
  "Zipfs" = c(2,2),
  "Geom"  = c(0.9,0.9)
)

M = 60
nome_base1 = paste0("save/Prob_msteps_")
igrid = c(1:3)
i = 1

stat_names = c("Strue","ShInd","ShMLE","Ktrue","KInd","KMLE","Kpooled","K1true","K1Ind","K1MLE",
               "K1indep","K2true","K2Ind","K2MLE","K2indep")
qnt_names = c(rep("(S)",3),rep("(K)",4),rep("(K1)",4),rep("(K2)",4))
select_idx = list(c(1:3),c(4:6),c(8:10),c(12:14))
for(i in igrid){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  param = experiments[[i]]
  p1 = param[1]; p2 = param[2]
  
  trim_p1 = get_first3digits(p1,4)
  trim_p2 = get_first3digits(p2,4)
  
  filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat") 
  load(filename)
  
  hh = 1
  for(hh in 1:4){
    
    ## Read results
    res_all = vector("list", length(select_idx[[hh]]))
    names(res_all) = stat_names[select_idx[[hh]]]
    for(jj in seq_along(select_idx[[hh]])){
      temp <- simplify2array(lapply(inner_result, function(x) 1-x[[select_idx[[hh]][jj]]]))
      res_all[[jj]] <- apply(temp, c(1, 2), quantile, probs = c(0.5))
      
      img_name = paste0(img_folder,"Pr_m_steps_",name_exp,"_",stat_names[select_idx[[hh]][jj]],"_",trim_p1,"_",trim_p2,".pdf")
      ## Plot
      if(save_img)
        pdf(img_name, width = width, height = height)
      par(mfrow = c(1,1), mar = c(2.5,2.5,1,3), mgp=c(1.5,0.5,0), cex = 2)
      image( n2grid, mgrid,
             t(res_all[[jj]]),
             zlim = c(0,1),
             col = mycol,
             xlab = "n",
             ylab = paste0("Additional m - ",qnt_names[select_idx[[hh]][jj]]),
             main = " ",
             axes = FALSE )
      axis(1, at = xpos, labels = xlabs, cex.axis = 0.85 )
      axis(2, at = mgrid_pos, labels = mgrid_pos)
      box()
      fields::image.plot(
        n2grid, mgrid,t(res_all[[jj]]),
        zlim = c(0,1),
        col = mycol,
        legend.only = TRUE,
        horizontal = FALSE,
        legend.width = 1.2,            # controls legend thickness
        legend.shrink = 0.8,           # smaller legend
        legend.mar = 5.75,                # margin from image
        legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
      )
      if(save_img)
        dev.off()
    }


  }
  
}

# Additive -----------------------------------------------------------------

name_exp = "Dirichelt"
p1 = 0.5; p2 = 0.5
params = c(p1,p2)
masses = c(0,0.5,1)
params_grid = expand.grid(params,params,masses)
params_grid = params_grid[c(4,8,12),] # <--- select here the experiments to run

trim_p1 = get_first3digits(p1,4)
trim_p2 = get_first3digits(p2,4)

nome_base1 = paste0("save/Prob_msteps_Additive_")
igrid = c(1:3)
i = 1

stat_names = c("Strue","ShInd","ShMLE","Ktrue","KInd","KMLE","Kpooled","K1true","K1Ind","K1MLE",
               "K1indep","K2true","K2Ind","K2MLE","K2indep")
select_idx = list(c(1:3),c(4:6),c(8:10),c(12:14))
for(i in seq_along(masses)){
  mass = masses[i]
  trim_mass = get_first3digits(mass,4)
  filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat") 
  load(filename)
  
  hh = 1
  for(hh in 1:4){
    
    ## Read results
    res_all = vector("list", length(select_idx[[hh]]))
    names(res_all) = stat_names[select_idx[[hh]]]
    for(jj in seq_along(select_idx[[hh]])){
      temp <- simplify2array(lapply(inner_result, function(x) 1-x[[select_idx[[hh]][jj]]]))
      res_all[[jj]] <- apply(temp, c(1, 2), quantile, probs = c(0.5))
      
      img_name = paste0(img_folder,"Add_",name_exp,"_",stat_names[select_idx[[hh]][jj]],"_",trim_p1,"_",trim_p2,"_",trim_mass,".pdf")
      ## Plot
      if(save_img)
        pdf(img_name, width = width, height = height)
      par(mfrow = c(1,1), mar = c(2.5,2.5,1,3), mgp=c(1.5,0.5,0), cex = 2)
      image( n2grid, mgrid,
             t(res_all[[jj]]),
             zlim = c(0,1),
             col = mycol,
             xlab = "n",
             ylab = paste0("Additional m - ",qnt_names[select_idx[[hh]][jj]]),
             main = " ",
             axes = FALSE )
      axis(1, at = xpos, labels = xlabs, cex.axis = 0.85 )
      axis(2, at = mgrid_pos, labels = mgrid_pos)
      box()
      fields::image.plot(
        n2grid, mgrid,t(res_all[[jj]]),
        zlim = c(0,1),
        col = mycol,
        legend.only = TRUE,
        horizontal = FALSE,
        legend.width = 1.2,            # controls legend thickness
        legend.shrink = 0.8,           # smaller legend
        legend.mar = 5.75,                # margin from image
        legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
      )
      if(save_img)
        dev.off()
    }
    
    
  }
  
}







