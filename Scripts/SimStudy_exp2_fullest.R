# Init --------------------------------------------------------------------
# setwd("/home/lucia.paci/Lucia/Ale/ScriptSpecies_shared/SimStudyRaf")
# setwd("~/Documents/Projects/ScriptSpecies_shared/SimStudyRaf")
setwd("C:/Users/colom/ScriptSpecies_shared/SimStudyRaf")
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
mycol_pal1 = hcl.colors(n=6,palette = "Zissou1")
mycol_pal2 = hcl.colors(n=6,palette = "Temps")

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}

source("SimStudyRaf_functions.R")

# parametri comuni --------------------------------------------------------
seed0 = 271296
set.seed(seed0)


num_cores = 7 # parallel::detectCores() - 2
Nrep = 140
BOiter = 1
keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)
par_est_cor = TRUE; par_est_mom = FALSE

ylim_plot_simmodel = c(0,150)

# Yue6060 - n_j = c(400,400) ----------------------------------------------
ylim_plot_yue6060 = c(0,75)

## 8080  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue0808")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue0808")
nome_base_plot = "Yue6060 - 080v080 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.8
alpha2 = 0.8
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_0808_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

# pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
# dev.off()



## 8585  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue085085")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue085085")
nome_base_plot = "Yue6060 - 085v085 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha = 0.85
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_085085_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
dev.off()


## 9090  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue0909")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue0909")
nome_base_plot = "Yue6060 - 09v09 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha = 0.90
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })

final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_0909_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
dev.off()


## 8085  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue08085")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue08085")
nome_base_plot = "Yue6060 - 08v085 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.8
alpha2 = 0.85
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_08085_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

# pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
# dev.off()


## 8090  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue0809")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue0809")
nome_base_plot = "Yue6060 - 08v09 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.8
alpha2 = 0.9
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_0809_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

# pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
# dev.off()


## 8590  ---------------------------------------

macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/",macchina,"_res_fullestcor_Yue08509")
nome_base2 = paste0("save/",macchina,"_tabres_fullestcor_Yue08509")
nome_base_plot = "Yue6060 - 085v09 (n = "

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)

experiment = 400
cat("\n ... Sim from model, n_j = ",experiment," ... \n")


experiment_name = as.character(experiment)
n_j = c(experiment,experiment)

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.85
alpha2 = 0.90
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)

names_species = 1:M

sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

data_list = lapply(1:2, function(j){sample(names_species,size = n_j[j], prob = w_j[[j]], replace = TRUE)})
data = matrix(0,nrow = M, ncol = 2 )

for(m in 1:M){
  data[m,1] = length(which(data_list[[1]] == m))
  data[m,2] = length(which(data_list[[2]] == m))
}


data = data[apply(data, 1, function(x){x[1]+x[2]>0}),]

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,"_",experiment,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,"_",experiment,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = ylim_plot_yue6060
nome_fig_yue = paste0("img/",macchina,"_fullestcor_yue6060_08509_",experiment,".pdf")
titolo_plot = paste0(nome_base_plot,experiment,")")
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

pdf(nome_fig_yue,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = " ", ylab = " ", xaxt = "n",
     main = titolo_plot)
mtext(text = nomi_exp, side = 1, line = 0.3, at = pos, las = 1, cex = 0.8)
grid(lty = 1,lwd = 1, col = "gray90" )
abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  boxplot(results[[i]]$Stot_noi,
          at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol_pal2[1], pch = 16)
  # boxplot(results[[i]]$Stot_Chao,
  #         at = (2*i-1)+0.25, add = T, col = tail(mycol_pal2,1), pch = 16)
}
dev.off()
