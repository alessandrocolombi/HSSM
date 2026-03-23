wd = "C:/Users/colom/ScriptSpecies_shared/RevBA"
setwd(wd)

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))

# Image options -----------------------------------------------------------

save_img = TRUE
width = 56; height = 8
cex.lab = 2
cex.axis = 2
cex = 2

# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")

# a) (M=60) d = 2 -----------------------------------------------

experiments = list(
  "Zipfs" = c(2,2),
  "Geom"  = c(0.9,0.9),
  "Dirichelt" = c(0.5,0.5)
)

M = 60
n = 1000
ngrid = ceiling(seq(10,n, length.out = 100))
BB = 50

shuffle = TRUE

statistics = c("Shared","Distinct","Dist. 1","Dist. 2","Morosita")
ylabs = c("S","K","K1","K2")
save_plot_all = FALSE

i = 1
igrid = c(1:3)

for(i in igrid){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
  
  idx = 1
  for(idx in 1:1 ){
    # trimmed_param = get_first3digits(p,4)
    p1 = params[1]; p2 = params[2]
    cat("\n p1 = ",p1," vs p2 = ",p2,"\n")
    prob_true_mat1 <- prob_true_mat2 <- matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat1 = t(apply(prob_true_mat1, 1, function(x) { 
      pj = sim_generic(name_exp,M,p1)
      if(shuffle)
        pj = pj[sample(1:M,M)]
      pj  } ))
    prob_true_mat2 = t(apply(prob_true_mat2, 1, function(x) { 
      pj = sim_generic(name_exp,M,p2)
      if(shuffle)
        pj = pj[sample(1:M,M)]
      pj } ))
    
    mat <- lapply(ngrid, function(nn) {
      
      counts_list1 <- apply(prob_true_mat1, 1, function(pj) {
        rmultinom(1, size = nn, prob = pj)
      }) # (M x BB) matrix
      counts_list2 <- apply(prob_true_mat2, 1, function(pj) {
        rmultinom(1, size = nn, prob = pj)
      }) # (M x BB) matrix
      
      counts_stacked = rbind(counts_list1,counts_list2)
      mor = apply(counts_stacked, 2, function(x){
        data_mor = matrix(0,nrow = M, ncol = 2)
        data_mor[,1] = x[1:M]; data_mor[,2] = x[(M+1):(2*M)]
        Mor_est(data_mor)
      })
      
      counts_list = counts_list1 + counts_list2
      K = colSums(counts_list > 0)
      K1 = colSums(counts_list1 > 0)
      K2 = colSums(counts_list2 > 0)
      S = K1 + K2 - K
      cbind(S,K,K1,K2,mor)  # (BB x 5)
    })
    res <- simplify2array(mat) # (BB x 5 x ngrid)
    AccCrv_all = apply(res, c(2,3), quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)
    
    # Plot Accumulation curves
    ylim_max = max(AccCrv_all)
    if(save_img)
      pdf(paste0("img/AccumulationCurves/AccCrv_",name_exp,"_1row.pdf"), width = width, height = height)
    par(mfrow = c(1,4), bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,0.25,1), cex = 2, las = 1)
    for(hh in 1:4){
      main_plt = paste0(name_exp," - ",statistics[hh]," - ", p1,"vs",p2)
      xx = AccCrv_all[,hh,]
      plot(x = 0, y = 0, type = "n",
           main = "", xlab = "n", ylab = ylabs[hh],
           ylim = c(0,ylim_max+1),
           xlim = c(0,n+1),
           pch = 1) # init plot
      polygon( c(ngrid, rev(ngrid)),
               c(xx[1,], rev(xx[3,])),
               col = "grey75",
               border = NA) # plot in-sample bands
      points(x = ngrid, y = xx[2,], type = "l", lwd = 3) # plot mean obs
    }
    if(save_img)
      dev.off()
    
    
  }
}

# b) Additive - d = 2 ----------------------------------------------------

# M is the common part
# Mj is the jth specific part
# Qj = M + Mj is the total number of species in the jth group
# Q is the total

Q = 80
n = 1000
ngrid = ceiling(seq(10,n, length.out = 100))

experiments = list(
  "Dirichelt" = c(0.5,0.5)
)
Nexp = length(experiments)
masses = c(0,0.5,1)

BB = 50

shuffle = TRUE

statistics = c("Shared","Distinct","Dist. 1","Dist. 2","Morosita")
ylabs = c("S","K","K1","K2")

i = 1
igrid = c(1)

for(i in igrid){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
  
  idx = 1
  for(idx in seq_along(masses) ){
    # trimmed_param = get_first3digits(p,4)
    p1 = params[1]; p2 = params[2]; mass = masses[idx]
    cat("\n p1 = ",p1," vs p2 = ",p2," || ","mass = ",mass,"\n")
    M = floor(mass*Q); M1 <- M2 <- floor( Q*(1-mass)/2 )
    if( (M+M1+M2) != Q )
      M = M + (Q-M-M1-M2)
    
    Q1 = M+M1; Q2 = M+M2
    
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
    
    # Compute mass of common part
    mass_common = rep(0,BB)
    if(M > 0)
      mass_common = rowSums(prob_true_mat1[,1:M])
    
    mat <- lapply(ngrid, function(nn) {
      
      counts_list1 <- apply(prob_true_mat1, 1, function(pj) {
        rmultinom(1, size = nn, prob = pj)
      }) # (Q x BB) matrix
      counts_list2 <- apply(prob_true_mat2, 1, function(pj) {
        rmultinom(1, size = nn, prob = pj)
      }) # (Q x BB) matrix
      
      counts_stacked = rbind(counts_list1,counts_list2)
      mor = apply(counts_stacked, 2, function(x){
        data_mor = matrix(0,nrow = Q, ncol = 2)
        data_mor[,1] = x[1:Q]; data_mor[,2] = x[(Q+1):(2*Q)]
        Mor_est(data_mor)
      })
      
      counts_list = counts_list1 + counts_list2
      K = colSums(counts_list > 0)
      K1 = colSums(counts_list1 > 0)
      K2 = colSums(counts_list2 > 0)
      S = K1 + K2 - K
      cbind(S,K,K1, K2, mor)  # (BB x 5)
    })
    res <- simplify2array(mat) # (BB x 5 x ngrid)
    AccCrv_all = apply(res, c(2,3), quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)
    
    # Plot Accumulation curves
    ylim_max = max(AccCrv_all)
    if(save_img)
      pdf(paste0("img/AccumulationCurves/AccCrv_Add_",name_exp,"_",idx,"_1row.pdf"), width = width, height = height)
    par(mfrow = c(1,4), bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,0.25,1), cex = 2, las = 1)
    for(hh in 1:4){
      xx = AccCrv_all[,hh,]
      plot(x = 0, y = 0, type = "n",
           main = "", xlab = "n", ylab = ylabs[hh],
           ylim = c(0,ylim_max+1),
           xlim = c(0,n+1),
           pch = 1) # init plot
      polygon( c(ngrid, rev(ngrid)),
               c(xx[1,], rev(xx[3,])),
               col = "grey75",
               border = NA) # plot in-sample bands
      points(x = ngrid, y = xx[2,], type = "l", lwd = 3) # plot mean obs
    }
    if(save_img)
      dev.off()
  }
  
}











