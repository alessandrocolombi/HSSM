# Init --------------------------------------------------------------------
setwd("C:/Users/colom/ScriptSpecies_shared/DatiTrieste")

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(vegan)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("TrainTest_functions.R")
source("../SimStudyRaf/SimStudyRaf_functions.R")

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}


# colors ------------------------------------------------------------------


mycol_pal1 = hcl.colors(n=6,palette = "Zissou1")
mycol_pal2 = hcl.colors(n=6,palette = "Temps")



# load data ---------------------------------------------------------------


counts_all = read.table(file = "save/Ants_BB_OL.txt") # read data counts
# tolgo le due specie troppo abbondanti
counts_all = counts_all[,-c(10,22)]

d   = nrow(counts_all) # number of groups (must be two)
n_j_all = apply(counts_all[,-1],1,sum) # number of observations in each group
n_all   = sum(n_j_all) # total number of observations

K12_all = ncol(counts_all)-1 # global number of distinct species
K1_all  = ncol(counts_all[1,c(TRUE,apply(counts_all[1,-1], 2, sum)>0)])-1 # local number of distinct species (group 1)
K2_all  = ncol(counts_all[2,c(TRUE,apply(counts_all[2,-1], 2, sum)>0)])-1 # local number of distinct species (group 2)
S12_all = K1_all+K2_all-K12_all  # number of distinct shared species

data_all = as.matrix(counts_all[,-1])  # build data matrix for MCMC sampling
data = matrix(0,nrow = K12_all, ncol = 2)
colnames(data) = c("BB", "OL")
names_species = colnames(data_all)
row.names(data) = names_species
data[,1] = data_all[1,]; data[,2] = data_all[2,]
data_norm = apply(data,2,function(x){x/sum(x)})
w_j = lapply(1:2, function(j){data_norm[,j]})



# Parametri run -----------------------------------------------------------

seed0 = 271296
set.seed(seed0)


num_cores = 7 # parallel::detectCores() - 2
Nrep = 140
BOiter = 1
keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)
par_est_cor = TRUE; par_est_mom = FALSE



# Run ---------------------------------------------------------------------


macchina =  "PC" #"Catt" # "VMB"
nome_base1 = paste0("save/TrainTest_analysis/",macchina,"_Prediction_fullest_res")
nome_base2 = paste0("save/TrainTest_analysis/",macchina,"_Prediction_fullest_tabres")
nome_base_plot = "BB vs OL"

keeps = seq(0.1,0.9,by = 0.05)
Lkeeps = length(keeps)


sp1p2 = sum( w_j[[1]] * w_j[[2]] )
ssp1  = sum( w_j[[1]] * w_j[[1]] )
ssp2  = sum( w_j[[2]] * w_j[[2]] )

K12 = nrow(data) # vero numero di specie distinte
K1  = length(which(data[,1]>0)) # vero numero di distinte nella prima popolazione
K2  = length(which(data[,2]>0)) # vero numero di distinte nella seconda popolazione
S12 = K1 + K2 - K12 # vero numero di condivise tra le due popolazioni

results = lapply( X = keeps, FUN = TrainTest_fullest_param,
                  data = data, n_j = n_j_all,
                  seed0 = seed0, Nrep = Nrep, num_cores = num_cores, 
                  BOiter = BOiter, par_est_cor = par_est_cor, par_est_mom = par_est_mom,
                  sp1p2 = sp1p2, ssp1 = ssp1, ssp2 = ssp2)

name_save_1 = paste0(nome_base1,".Rdat")
save(results, file = name_save_1)


MC_Chao    = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_Chao)^2 ) })
MC_noi     = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi)^2 ) })
MC_noi_est = sapply(results, FUN = function(res_keep){ mean( (res_keep$Stot_true - res_keep$Stot_noi_est)^2 ) })
MC_noi_K   = sapply(results, FUN = function(res_keep){ mean( (res_keep$Ktot_true - res_keep$Ktot_noi)^2 ) })


final_resul = tibble("keep" = keeps, "MC_Chao" = MC_Chao, 
                     "MC_noi" = MC_noi,"MC_noi_est" = MC_noi_est,
                     "MC_noi_K" = MC_noi_K)


name_save_2 = paste0(nome_base2,".Rdat")
save(final_resul,file = name_save_2)


## plot
ylim_plot = c(0,20)
nome_fig = paste0("img/TrainTest_analysis/",macchina,"_Prediction_fullest.pdf")
titolo_plot = paste0(nome_base_plot)
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")


nomi_exp = as.character(final_resul$keep)
pos  = seq(1,2*length(nomi_exp),by = 2)

pdf(nome_fig,width = 8,height = 6)
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


