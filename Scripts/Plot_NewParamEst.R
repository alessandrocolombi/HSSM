setwd("C:/Users/colom/ScriptSpecies_shared/DatiTrieste")

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
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

name_folder_img = "G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/Script_HSSM/DatiTrieste/img/BB_OL/NewParamEst"

# colors ------------------------------------------------------------------


mycol_pal1 = hcl.colors(n=50,palette = "Zissou1")
mycol_pal2 = hcl.colors(n=6,palette = "Temps")

# load data ---------------------------------------------------------------

counts_all = read.table(file = "save/Ants_BB_OL.txt") # read data counts
# tolgo le due specie troppo abbondanti
# counts_all = counts_all[,-c(10,22)]

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

# Recode labels
old_labels = as.factor(1:K12_all)
new_labels1 = as.factor(sort(w_j[[1]], decreasing = TRUE, index.return = TRUE)$ix)
new_labels2 = as.factor(sort(w_j[[2]], decreasing = TRUE, index.return = TRUE)$ix)
names_species_ez = as.character(1:K12_all)
names_species_ez1 = as.character(new_labels1)
names_species_ez2 = as.character(new_labels2)
shared = apply(data, 1, function(x){x[1]>0 && x[2]>0})
names(shared) = names_species_ez
col_shared = rep("black", K12_all)
col_shared[shared] = "darkorange"

col_shared1 = col_shared[as.numeric(names_species_ez1)]
col_shared2 = col_shared[as.numeric(names_species_ez2)]
# recode_map <- setNames(new_labels, old_labels)
# VI_sara$cl <- recode(VI_sara$cl, !!!recode_map)

# Descriptive analysis ----------------------------------------------------

# All data
x1 = sort(data[,1], decreasing = TRUE)[1:K1_all]
names(x1) = names_species_ez1[1:K1_all]
x2 = sort(data[,2], decreasing = TRUE)[1:K2_all]
names(x2) = names_species_ez2[1:K2_all]
par(mfrow = c(2,1),mar = c(4,4,2,1), bty = "l")
barplot( x1, xlab = " ", ylab = "# of individuals", 
         cex.names = 0.7, col = col_shared1[1:K1_all] )
barplot( x2, xlab = " ", ylab = "# of individuals", 
         cex.names = 0.7, col = col_shared2[1:K2_all] )

# Remove most abundant species
x1 = sort(data[,1], decreasing = TRUE)[2:K1_all]
names(x1) = names_species_ez1[2:K1_all]
x2 = sort(data[,2], decreasing = TRUE)[2:K2_all]
names(x2) = names_species_ez2[2:K2_all]
par(mfrow = c(2,1),mar = c(4,4,2,1), bty = "l")
barplot( x1, xlab = " ", ylab = "# of individuals", 
         cex.names = 0.7, col = col_shared1[2:K1_all] )
barplot( x2, xlab = " ", ylab = "# of individuals", 
         cex.names = 0.7, col = col_shared2[2:K2_all] )

data_norm = apply(data[-c(9,21),],2,function(x){x/sum(x)})
w_j = lapply(1:2, function(j){data_norm[,j]})

# Species names transformation
ordered_names_1 = names(sort(w_j[[1]], decreasing = TRUE))
ordered_names_2 = names(sort(w_j[[2]], decreasing = TRUE))

my_trim_names = function(name) {
  # Split the name into two parts based on the underscore
  parts <- strsplit(name, "_")[[1]]
  # Extract the first 4 letters of each part
  part1 <- substr(parts[1], 1, 6)
  part2 <- substr(parts[2], 1, 4)
  # Combine the parts with a space
  paste(part1, part2, sep = " ")
}

trimmed_names_1 <- sapply(ordered_names_1, my_trim_names)
trimmed_names_2 <- sapply(ordered_names_2, my_trim_names)
shared_unordered = shared[-c(9,21)] 
names(shared_unordered) = names_species[-c(9,21)]
ordering_index_1 = sort(w_j[[1]], decreasing = TRUE, index.retur = TRUE)$ix
ordering_index_2 = sort(w_j[[2]], decreasing = TRUE, index.retur = TRUE)$ix
shared_ordered_1 = shared_unordered[ordering_index_1]
shared_ordered_2 = shared_unordered[ordering_index_2]
col_shared_1 = rep("black", K12_all)
col_shared_1[shared_ordered_1] = "red"
col_shared_2 = rep("black", K12_all)
col_shared_2[shared_ordered_2] = "red"



par(mfrow = c(1,2), mar = c(4.5,4,2,1), bty = "l")
plot(x = 1:K12_all, y = sort(w_j[[1]], decreasing = TRUE), 
     ylim = c(0,0.33), xaxt = "n", main = "Area 1 - Bosco Bovedo",
     type = "l", lty = 2, xlab = "", ylab = "Empirical species proportion")
points(x = 1:K12_all, y = sort(w_j[[1]], decreasing = TRUE), 
       type = "p", pch = 16, col = col_shared_1)
axis(side = 1, at = 1:K12_all, labels = FALSE)  # Suppress default labels
text(x = 1:K12_all, y = par("usr")[3] - 0.02,  # Adjust y-position to place labels
     labels = trimmed_names_1, col = col_shared_1, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)

plot(x = 1:K12_all, y = sort(w_j[[2]], decreasing = TRUE), 
     ylim = c(0,0.33), xaxt = "n", main = "Area 2 - Orto Lapidario",
     type = "l", lty = 2, xlab = "", ylab = "Empirical species proportion")
points(x = 1:K12_all, y = sort(w_j[[2]], decreasing = TRUE), 
       type = "p", pch = 16, col = col_shared_2)
axis(side = 1, at = 1:K12_all, labels = FALSE)  # Suppress default labels
text(x = 1:K12_all, y = par("usr")[3] - 0.02,  # Adjust y-position to place labels
     labels = trimmed_names_2, col = col_shared_2, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)

# Rarefaction -------------------------------------------------------------

# Nsort = 50
# seed  = 123
# Rar_curves  = Rarefaction_curves(data,  Nsort = Nsort, seed0 = seed         )
# Rar_curves1 = Rarefaction_curve_d1( data[,1],  Nsort = Nsort, seed0 = seed  )
# Rar_curves2 = Rarefaction_curve_d1( data[,2],  Nsort = Nsort, seed0 = seed  )
# 
# save(Rar_curves, file = paste0(name_folder_img,"Rar_curves.Rdat") )
# save(Rar_curves1, file = paste0(name_folder_img,"Rar_curves1.Rdat") )
# save(Rar_curves2, file = paste0(name_folder_img,"Rar_curves2.Rdat") )

load(paste0(name_folder_img,"Rar_curves.Rdat"))
load(paste0(name_folder_img,"Rar_curves1.Rdat"))
load(paste0(name_folder_img,"Rar_curves2.Rdat"))

# plot distinct
grid_curve = 1:length(Rar_curves$K12obs_summary[2,]) 
rar_curve  = Rar_curves$K12obs_summary[2,]
rar_curveLB  = Rar_curves$K12obs_summary[1,]
rar_curveUB  = Rar_curves$K12obs_summary[3,]
ylab_name = "# global distinct species"
col_curve = "black"
par(mfrow = c(1,1),mar = c(4,4,2,1), bty = "l")
plot(x = grid_curve, y = rar_curve, 
     type = "l", lwd = 1, xlab = "# obs.", ylab = ylab_name)
polygon( c(grid_curve, rev(grid_curve)),
         c(rar_curveLB, rev(rar_curveUB)),
         col = ACutils::t_col(col_curve, percent = 75),
         border = NA) # plot in-sample bands

# plot shared
grid_curve = 1:length(Rar_curves$S12obs_summary[2,]) 
rar_curve  = Rar_curves$S12obs_summary[2,]
rar_curveLB  = Rar_curves$S12obs_summary[1,]
rar_curveUB  = Rar_curves$S12obs_summary[3,]
ylab_name = "# shared species"
col_curve = "black"
par(mfrow = c(1,1),mar = c(4,4,2,1), bty = "l")
plot(x = grid_curve, y = rar_curve, 
     type = "l", lwd = 1, xlab = "# obs.", ylab = ylab_name)
polygon( c(grid_curve, rev(grid_curve)),
         c(rar_curveLB, rev(rar_curveUB)),
         col = ACutils::t_col(col_curve, percent = 75),
         border = NA) # plot in-sample bands

# plot local distinct 1
grid_curve = 1:length(Rar_curves1$Kobs_summary[2,]) 
rar_curve  = Rar_curves1$Kobs_summary[2,]
rar_curveLB  = Rar_curves1$Kobs_summary[1,]
rar_curveUB  = Rar_curves1$Kobs_summary[3,]
ylab_name = "# local distinct species - 1"
col_curve = "black"
par(mfrow = c(1,1),mar = c(4,4,2,1), bty = "l")
plot(x = grid_curve, y = rar_curve, 
     type = "l", lwd = 1, xlab = "# obs.", ylab = ylab_name)
polygon( c(grid_curve, rev(grid_curve)),
         c(rar_curveLB, rev(rar_curveUB)),
         col = ACutils::t_col(col_curve, percent = 75),
         border = NA) # plot in-sample bands

# plot local distinct 2
grid_curve = 1:length(Rar_curves2$Kobs_summary[2,]) 
rar_curve  = Rar_curves2$Kobs_summary[2,]
rar_curveLB  = Rar_curves2$Kobs_summary[1,]
rar_curveUB  = Rar_curves2$Kobs_summary[3,]
ylab_name = "# local distinct species - 2"
col_curve = "black"
par(mfrow = c(1,1),mar = c(4,4,2,1), bty = "l")
plot(x = grid_curve, y = rar_curve, 
     type = "l", lwd = 1, xlab = "# obs.", ylab = ylab_name)
polygon( c(grid_curve, rev(grid_curve)),
         c(rar_curveLB, rev(rar_curveUB)),
         col = ACutils::t_col(col_curve, percent = 75),
         border = NA) # plot in-sample bands


# Plot results  --------------------------------------------
## Rar TrainTest --------------------------------------------

name_folder_save = "C:/Users/colom/ScriptSpecies_shared/DatiTrieste/save/TrainTest_analysis/fullest/"
name_exp_base = "TrainTest_analysis_fullest"
keeps = c(0.2,0.25,0.5,0.75,0.8)

keep = keeps[2]

for(keep in keeps){
  name_exp = paste0(name_folder_save,name_exp_base,remove_points(keep),".Rdat")
  load(file = name_exp)
  
  run[[1]]$S_tr_mean
  run[[1]]$ExpS_tr
  Ntrain = length(run); Lcut1 = 100; Lcutprev = 100
  
  n_j_tr = run[[1]]$n_j_tr
  n_j_ts = run[[1]]$n_j_ts
  n_tr = sum(n_j_tr)
  n_ts = sum(n_j_ts)
  Lout_tr = length(run[[1]]$ExpK_tr)
  Lout1_tr = length(run[[1]]$ExpK1_tr)
  Lout2_tr = length(run[[1]]$ExpK2_tr)
  Lout_ts = length(run[[1]]$ExpKnew_ts)
  Lout1_ts = length(run[[1]]$ExpK1new_ts)
  Lout2_ts = length(run[[1]]$ExpK2new_ts)
  
  
  Kobs_rep  = matrix(0,nrow = Ntrain, ncol = n_tr)
  Sobs_rep  = matrix(0,nrow = Ntrain, ncol = n_tr)
  Kts_rep  = matrix(0,nrow = Ntrain, ncol = n_ts)
  Sts_rep  = matrix(0,nrow = Ntrain, ncol = n_ts)
  
  K1obs_rep = matrix(0,nrow = Ntrain, ncol = n_j_tr[1]+n_j_ts[1])
  K2obs_rep = matrix(0,nrow = Ntrain, ncol = n_j_tr[2]+n_j_ts[2])
  
  ExpKobs_rep  = matrix(0,nrow = Ntrain, ncol = Lout_tr)
  ExpSobs_rep  = matrix(0,nrow = Ntrain, ncol = Lout_tr)
  ExpK1obs_rep = matrix(0,nrow = Ntrain, ncol = Lout1_tr)
  ExpK2obs_rep = matrix(0,nrow = Ntrain, ncol = Lout2_tr)
  
  ExpKts_rep  = matrix(0,nrow = Ntrain, ncol = Lout_ts)
  ExpSts_rep  = matrix(0,nrow = Ntrain, ncol = Lout_ts)
  ExpK1ts_rep = matrix(0,nrow = Ntrain, ncol = Lout1_ts)
  ExpK2ts_rep = matrix(0,nrow = Ntrain, ncol = Lout2_ts)
  
  for(it in 1:Ntrain){
    
    Kobs_rep[it,] = run[[it]]$K_tr_mean
    Sobs_rep[it,] = run[[it]]$S_tr_mean
    Kts_rep[it,] = run[[it]]$K_ts_mean
    Sts_rep[it,] = run[[it]]$S_ts_mean
    
    K1obs_rep[it, 1:run[[it]]$n_j_tr[1] ] = run[[it]]$K1_tr_mean
    K1obs_rep[it, (run[[it]]$n_j_tr[1]+1):(n_j_tr[1]+n_j_ts[1]) ] = run[[it]]$K1_tr + run[[it]]$K1_ts_mean
    K2obs_rep[it, 1:run[[it]]$n_j_tr[2] ] = run[[it]]$K2_tr_mean
    K2obs_rep[it, (run[[it]]$n_j_tr[2]+1):(n_j_tr[2]+n_j_ts[2]) ] = run[[it]]$K2_tr + run[[it]]$K2_ts_mean
    
    ExpKobs_rep[it,] = run[[it]]$ExpK_tr
    ExpSobs_rep[it,] = run[[it]]$ExpS_tr
    ExpK1obs_rep[it,] = run[[it]]$ExpK1_tr
    ExpK2obs_rep[it,] = run[[it]]$ExpK2_tr
    
    ExpKts_rep[it,] = run[[it]]$ExpKnew_ts
    ExpSts_rep[it,] = run[[it]]$ExpSnew_ts
    ExpK1ts_rep[it,] = run[[it]]$ExpK1new_ts
    ExpK2ts_rep[it,] = run[[it]]$ExpK2new_ts
  }
  
  Kobs_summary = apply( Kobs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  Kobs_mean    = apply( Kobs_rep,2,mean )
  Sobs_summary = apply( Sobs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  Sobs_mean    = apply( Sobs_rep,2,mean )
  Kts_summary = apply( Kts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  Kts_mean    = apply( Kts_rep,2,mean )
  Sts_summary = apply( Sts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  Sts_mean    = apply( Sts_rep,2,mean )
  
  K1obs_summary = apply( K1obs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  K1obs_mean    = apply( K1obs_rep,2,mean )
  K2obs_summary = apply( K2obs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  K2obs_mean    = apply( K2obs_rep,2,mean )
  
  ExpKobs_summary = apply( ExpKobs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpKobs_mean    = apply( ExpKobs_rep,2,mean )
  ExpSobs_summary = apply( ExpSobs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpSobs_mean    = apply( ExpSobs_rep,2,mean )
  ExpK1obs_summary = apply( ExpK1obs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpK1obs_mean    = apply( ExpK1obs_rep,2,mean )
  ExpK2obs_summary = apply( ExpK2obs_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpK2obs_mean    = apply( ExpK2obs_rep,2,mean )
  
  ExpKts_summary = apply( ExpKts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpKts_mean    = apply( ExpKts_rep,2,mean )
  ExpSts_summary = apply( ExpSts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpSts_mean    = apply( ExpSts_rep,2,mean )
  ExpK1ts_summary = apply( ExpK1ts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpK1ts_mean    = apply( ExpK1ts_rep,2,mean )
  ExpK2ts_summary = apply( ExpK2ts_rep,2,quantile, probs = c(0.025,0.5,0.975) )
  ExpK2ts_mean    = apply( ExpK2ts_rep,2,mean )
  
  
  
  
  Ktr_mean = mean(sapply(run,function(x){x$K_tr}))
  Str_mean = mean(sapply(run,function(x){x$S_tr}))
  K1tr_mean = mean(sapply(run,function(x){x$K1_tr}))
  K2tr_mean = mean(sapply(run,function(x){x$K2_tr}))
  n_1_tr_mean  = floor(mean(sapply(run,function(x){x$n_j_tr[1]})))
  n_2_tr_mean  = floor(mean(sapply(run,function(x){x$n_j_tr[2]})))
  
  temp = tibble(site = c(rep("A1",n_1_tr_mean),rep("A2",n_2_tr_mean)))
  set.seed(123)
  new_idx = sample(1:(n_1_tr_mean+n_2_tr_mean),size = (n_1_tr_mean+n_2_tr_mean), replace = F)
  temp_reordered = temp[new_idx,]
  res = vector("list",length = 3)
  res[[1]] = rep(0,(n_1_tr_mean+n_2_tr_mean)+1)
  res[[2]] = rep(0,(n_1_tr_mean+n_2_tr_mean)+1)
  for(it in 2:((n_1_tr_mean+n_2_tr_mean)+1)){
    # cat("\n it = ",it,"\n")
    x = temp_reordered[it-1,]
    j = -1
    if(x$site == "A1"){
      j = 1
    }else if(x$site == "A2"){
      j = 2
    }else{
      stop("ERRORE")
    }
    
    res[[1]][it] = res[[1]][it-1]
    res[[2]][it] = res[[2]][it-1]
    res[[j]][it] = res[[j]][it] + 1
    
  }
  n_j_list_n = vector("list",(n_1_tr_mean+n_2_tr_mean))
  for(i in 1:(n_1_tr_mean+n_2_tr_mean)){
    n_j_list_n[[i]] = c(res[[1]][i+1],res[[2]][i+1])
  }
  
  cutpoints_train = floor( seq(10,(n_1_tr_mean+n_2_tr_mean), length.out = Lcut1 ) )
  Lgrid_train = length(cutpoints_train)
  
  cutpoints1_train = floor( seq(10,n_1_tr_mean, length.out = Lcut1 ) )
  Lgrid1_train = length(cutpoints1_train)
  
  cutpoints2_train = floor( seq(10,n_2_tr_mean, length.out = Lcut1 ) )
  Lgrid2_train = length(cutpoints2_train)
  
  
  ## Condivise
  perc = remove_points(keep)
  main_sh = paste0("Shared, train = ",perc)
  par(mfrow = c(1,1),mar = c(4,4,2,1), bty = "l")
  plot(x = 0, y = 0, type = "n",
       main = main_sh, xlab = "#obs.", ylab = " ",
       ylim = c(0,16),
       xlim = c(0,max(n_tr+(1:n_ts))+1),
       pch = 1) # init plot
  points(x = 1:n_tr, y = Sobs_mean, type = "l", lwd = 3) # plot mean obs
  
  abline(v = n_tr, lty = 2, col = "black") # plot train line
  polygon( c(1:n_tr, rev(1:n_tr)),
           c(Sobs_summary[1,], rev(Sobs_summary[3,])),
           col = ACutils::t_col("black", percent = 75),
           border = NA) # plot in-sample bands
  points(x = n_tr+(1:n_ts), y = Str_mean + Sts_mean, type = "l", lwd = 3) # plot mean test
  points(x = cutpoints_train, y =  ExpSobs_mean, type = "l", lwd = 3, col = "darkred") # plot in-sample mean
  points(x = n_tr+run[[1]]$m, Str_mean + ExpSts_mean,
         type = "l", lwd = 3, col = "darkgreen" ) # plot out-of-sample mean
  polygon( c(n_tr+run[[1]]$m, rev(n_tr+run[[1]]$m)),
           c(Str_mean + ExpSts_summary[1,], rev(Str_mean + ExpSts_summary[3,])),
           col = ACutils::t_col("darkgreen", percent = 75),
           border = NA) # plot out-of-sample bands
  polygon( c(n_tr+(1:n_ts), rev(n_tr+(1:n_ts))),
           c(Str_mean + Sts_summary[1,], rev(Str_mean + Sts_summary[3,])),
           col = ACutils::t_col("black", percent = 75),
           border = NA) # plot in-sample bands
  
  
}



## Prediction Test set -----------------------------------------------------
mycol = c("darkgreen","darkorange")

name_folder_save = "C:/Users/colom/ScriptSpecies_shared/DatiTrieste/save/TrainTest_analysis/"
name_exp_base = "PC_Prediction_fullest_res"
name_exp_base2 = "PC_Prediction_fullest_tabres"
load(file = paste0(name_folder_save,name_exp_base,".Rdat"))
load(file = paste0(name_folder_save,name_exp_base2,".Rdat"))

save_plots = TRUE

## plot
keeps = seq(0.1,0.9,by = 0.05)
ylim_plot = c(0,20)
titolo_plot = " "
S12_true = results[[1]]$Stot_true
nomi = c("Noi","Chao")

name_folder_img = "G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/Script_HSSM/DatiTrieste/img/BB_OL/NewParamEst/img/"
nome_fig = paste0(name_folder_img,"PredictionTestSet",".pdf")

# nomi_exp = as.character(keeps) # in decimale
nomi_exp = as.character(keeps * 100) # in percentuale
pos  = seq(1,2*length(nomi_exp),by = 2)

if(save_plots)
  pdf(nome_fig,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,2.5,1,1), bty = "l", mgp = c(2.5, 0.95, 0) )
plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
     xlab = "% training set", ylab = " ", 
     xaxt = "n", yaxt = "n",
     main = titolo_plot,
     cex.main = 2, cex.lab = 2, cex.axis = 1.75)
axis(side = 2, las = 1, cex.axis = 1.75)
axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 2)
grid(lty = 1,lwd = 1, col = "gray90" )
# abline(h = S12_true, lty = 2, col = "red")
for(i in 1:length(results)){
  # boxplot(results[[i]]$Stot_noi,
  #         at = (2*i-1)-0.25, add = T, col = mycol[3], pch = 1)
  boxplot(results[[i]]$Stot_noi_est,
          at = (2*i-1)+0.25, add = T, col = mycol[1], 
          pch = 1, yaxt = "n")
  boxplot(results[[i]]$Stot_Chao,
          at = (2*i-1)-0.25, add = T, col = tail(mycol,1), 
          pch = 1, yaxt = "n")
}
legend("topright", c("Proposed","Chao2000"), 
       col = mycol, pch = 16,  cex = 2)
dev.off()


## PrSh1step ---------------------------------------------------------------

name_folder_save = "C:/Users/colom/ScriptSpecies_shared/DatiTrieste/save/"
name_exp_base = "PrSh1step_all"
load(paste0(name_folder_save,name_exp_base,".Rdat"))
pchs  = c(17,rep(16,5))
mycol = hcl.colors(n=4,palette = "Zissou1")
mycol = c("black", mycol[1:3], "darkred", "darkgreen")

name_folder_img = "G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/Script_HSSM/DatiTrieste/img/BB_OL/NewParamEst/img/"
nome_fig = paste0(name_folder_img,"PrSh1step",".pdf")
save_plots = TRUE

ngrid = c(seq(50,600,by = 50))
shift = c(0,5,10,10,-10,-5)
titolo_plot = paste0("")

if(save_plots)
  pdf(nome_fig,width = 8,height = 6)
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l", mgp = c(2.5, 0.95, 0) )
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.1), 
     xlab = "n", ylab = "", # ylab = "P(Snew > 0)",
     xaxt = "n", yaxt = "n",
     cex.main = 2, cex.lab = 2, cex.axis = 1.8)
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, las = 1, cex.axis = 1.75)
axis(side = 1, at = ngrid, labels = as.character(2 * ngrid), las = 1, cex.axis = 1.75)
# for(kk in 1:length(PrSh_ngrid[[1]])){
for(kk in c(2,4,6)){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = pchs[kk], col = mycol[kk], lty = 2)
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 1 )
}
legend("topright", c("Yue","Chao2000","Proposed"),
       col = mycol[c(2,4,6)], pch = 16, cex = 2)
# legend("topright", c("True","Yue1","Yue2","Chao","Use prop.","Full est."),
#        col = mycol[1:6], pch = 16)
dev.off()









