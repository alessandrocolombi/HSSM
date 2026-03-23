wd = "C:/Users/colom/ScriptSpecies_shared/RevBA"
setwd(wd)


# Libraries --------------------------------------------------------------------
library(HSSM)
# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}

# Colors and legends ------------------------------------------------------------------
mycol_base = c("darkgreen","darkred","deeppink","aquamarine")
mycol_list = list(rep(c(mycol_base[2],mycol_base[1]),6), 
                  rep(c(mycol_base[2],mycol_base[1],mycol_base[3]),6),
                  rep(c(mycol_base[2],mycol_base[1],mycol_base[4]),6),
                  rep(c(mycol_base[2],mycol_base[1],mycol_base[4]),6)
)
legends = list(c("Bayes I","Bayes II"),
               c("Bayes I","Bayes II","Pooled"),
               c("Bayes I","Bayes II","Indep."),
               c("Bayes I","Bayes II","Indep."))

pos3 = seq(2,21,by = 7); pos6 = seq(2,21,by = 3.5)
# (D,Z,G) -----------------------------------------------------------------

## a) Basic options
save_plots = FALSE
name_fld_img_base = "img/Paper/Exp1/"
name_folder_save = "save/"
nome_base1 = paste0("save/Pred_")



experiments = list(
  "Dirichelt" = c(0.1,0.5),
  "Zipfs" = c(1.3,2),
  "Geom"  = c(0.8,0.85,0.9)
)

M = 60
Nrep <- BB <- 100
nmax <- 800
ngrid <- seq(100,nmax,by = 100)
Lgrid = length(ngrid)


igrid = c(1:3)
i = 1
for(i in igrid){ # for each experiment (D,G,Z)
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
  params_grid <- rbind( cbind(params, params), t(combn(params, 2)) )
  choose_n = 4
  for(choose_n in 1:Lgrid){ # for each sample size 
    ## Define main obj to be plot
    ExpPred = vector("list",4)
    names(ExpPred) = c("Sh","K","K1","K2")
    ExpPred = lapply(ExpPred, function(x) matrix(ncol = 0, nrow = Nrep))
    Params_BayesI <- Params_Bayes2 <- rep(-1,4)
    
    idx = 3
    for(idx in 1:nrow(params_grid) ){ # for each setting
      p1 = params_grid[idx,1]; p2 = params_grid[idx,2]
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2,"\n")
    
    trim_p1 = get_first3digits(p1,4)
    trim_p2 = get_first3digits(p2,4)
    filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat")
    load(file = filename)
    
    # Compute RMSE
    ## --> Compute SHARED
    RMSE_Sh_Ind = sapply(results, function(x) {
      x_true  = na.omit(x$Stot_true)
      x_Ind   = na.omit(x$Stot_est)
      rmse = sqrt((x_true - x_Ind)^2)
      rmse
    } )
    RMSE_Sh_MLE = sapply(results, function(x) {
      x_true  = na.omit(x$Stot_true)
      x_MLE   = na.omit(x$Stot_est_mle)
      rmse = sqrt((x_true - x_MLE)^2)
      rmse
    } )
    
    ## --> Compute DISTINCT
    RMSE_K_Ind = sapply(results, function(x) {
      x_true  = na.omit(x$Ktot_true)
      x_Ind   = na.omit(x$Ktot_est)
      rmse = sqrt((x_true - x_Ind)^2)
      rmse
    } )
    RMSE_K_MLE = sapply(results, function(x) {
      x_true  = na.omit(x$Ktot_true)
      x_MLE   = na.omit(x$Ktot_est_mle)
      rmse = sqrt((x_true - x_MLE)^2)
      rmse
    } )
    RMSE_K_pooled = sapply(results, function(x) {
      x_true  = na.omit(x$Ktot_true)
      x_pooled   = na.omit(x$Kpooled_est)
      rmse = sqrt((x_true - x_pooled)^2)
      rmse
    } )
    
    ## --> Compute LOCAL DISTINCT - Area 1
    RMSE_K1_Ind = sapply(results, function(x) {
      x_true  = na.omit(x$K1tot_true)
      x_Ind   = na.omit(x$K1tot_est)
      rmse = sqrt((x_true - x_Ind)^2)
      rmse
    } )
    RMSE_K1_MLE = sapply(results, function(x) {
      x_true  = na.omit(x$K1tot_true)
      x_MLE   = na.omit(x$K1tot_est_mle)
      rmse = sqrt((x_true - x_MLE)^2)
      rmse
    } )
    RMSE_K1_indep1 = sapply(results, function(x) {
      x_true  = na.omit(x$K1tot_true)
      x_indep1 = na.omit(x$K1indep_est)
      rmse = sqrt((x_true - x_indep1)^2)
      rmse
    } )
    
    ## --> Compute LOCAL DISTINCT - Area 2
    RMSE_K2_Ind = sapply(results, function(x) {
      x_true  = na.omit(x$K2tot_true)
      x_Ind   = na.omit(x$K2tot_est)
      rmse = sqrt((x_true - x_Ind)^2)
      rmse
    } )
    RMSE_K2_MLE = sapply(results, function(x) {
      x_true  = na.omit(x$K2tot_true)
      x_MLE   = na.omit(x$K2tot_est_mle)
      rmse = sqrt((x_true - x_MLE)^2)
      rmse
    } )
    RMSE_K2_indep2 = sapply(results, function(x) {
      x_true  = na.omit(x$K2tot_true)
      x_indep2 = na.omit(x$K2indep_est)
      rmse = sqrt((x_true - x_indep2)^2)
      rmse
    } )
    
    # Update global object
    ExpPred$Sh = cbind(ExpPred$Sh,RMSE_Sh_Ind[,choose_n],RMSE_Sh_MLE[,choose_n])
    ExpPred$K  = cbind(ExpPred$K,RMSE_K_Ind[,choose_n],RMSE_K_MLE[,choose_n],RMSE_K_pooled[,choose_n])
    ExpPred$K1 = cbind(ExpPred$K1,RMSE_K1_Ind[,choose_n],RMSE_K1_MLE[,choose_n],RMSE_K1_indep1[,choose_n])
    ExpPred$K2 = cbind(ExpPred$K2,RMSE_K2_Ind[,choose_n],RMSE_K2_MLE[,choose_n],RMSE_K2_indep2[,choose_n])
    
    
    # Save parameters
    Params_BayesI = c(median(results[[choose_n]]$gamma1_mle),
                      median(results[[choose_n]]$gamma2_mle),
                      median(results[[choose_n]]$lambda_mle),
                      median(results[[choose_n]]$Ktot_true))
    Params_BayesII = c(median(results[[choose_n]]$gamma1_ind),
                       median(results[[choose_n]]$gamma2_ind),
                       median(results[[choose_n]]$lambda_ind),
                       median(results[[choose_n]]$Ktot_true))
    
    ## Compute and plot qM_star
    Params_Bayes = list(Params_BayesI,Params_BayesII)
    qM_list = vector("list",2)
    for(idx_bayes in 1:2){
      gamma_vec = c(Params_Bayes[[idx_bayes]][1], Params_Bayes[[idx_bayes]][2])
      Lambda = Params_Bayes[[idx_bayes]][3]
      Kn = Params_Bayes[[idx_bayes]][4]
      qM_vec = c()
      for(mm in 0:1000){
        temp = HSSM::log_qM_post(m = mm, prior = "Poisson", prior_param = list("lambda"=Lambda),
                                 k = Kn, n_j = c(ngrid[choose_n],ngrid[choose_n]), 
                                 gamma_j = gamma_vec, log_V = -Inf, M_max = 200)[1]
        qM_vec = c(qM_vec, exp(temp))
        if(sum(qM_vec) >= (0.99))
          break
      }
      qM_list[[idx_bayes]] = qM_vec
    }
    LqMpost = max(sapply(qM_list,length))
    ymax = max(sapply(qM_list,max));ymin = 0
    ylim_plot = c(ymin,ymax)
    xgrid = 0:(LqMpost-1)
    xpos <- xlabs <- round(seq(xgrid[1],tail(xgrid,1),length.out = 5))
    
    nome_fig = paste0(name_fld_img_base,name_exp,"/n",ngrid[choose_n],"/","Mstar_",idx,".pdf")
    titolo_plot = ""
    
    if(save_plots)
      pdf(nome_fig,width = 14,height = 6)
    par(mfrow = c(1,1), mar = c(3,4,2,1), bty = "l", mgp = c(3, 0.95, 0), cex = 2, las = 1 )
    plot(0,0,type = "n",xlim = range(xgrid), ylim = ylim_plot,
         xlab = " ", ylab = paste0("Mstar"), 
         xaxt = "n", yaxt = "n",
         main = titolo_plot,
         cex.main = 1, cex.lab = 1, cex.axis = 1)
    axis(side = 2, las = 1, cex.axis = 1)
    axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = 1)
    grid(lty = 1,lwd = 1, col = "gray90" )
    points(x = xgrid[1:length(qM_list[[1]])], 
           y = qM_list[[1]], col = "darkred", 
           type = "b", lwd = 3, pch = 16)
    points(x = xgrid[1:length(qM_list[[2]])], 
           y = qM_list[[2]], 
           col = "darkgreen", type = "b", lwd = 3, pch = 16)
    legend("topright",c("Bayes I","Bayes II"), col = c("darkred","darkgreen"),
           pch = 16, cex = 1.25,  border = NA, bty = "n" )
    if(save_plots)
      dev.off()
  }
      
    ## Set plotting positions
    pos_ii = pos3
    if(name_exp == "Geom")
      pos_ii = pos6
    positions = list(pos_ii,pos_ii,pos_ii,pos_ii)
    shf1 = c(-0.25,0.25); shf2 = c(-0.5,0,0.5)
    shifts = list(shf1,shf2,shf2,shf2)
      
    ## Run plot
    hh = 2
    ylabs = c("S","K","K1","K2")
    for(hh in 1:4){ # for each quantity
      quantity = names(ExpPred)[hh]
      ymax = max(ExpPred[[hh]])
      ymin = 0
      ylim_plot = c(ymin,ymax)
        
      nome_fig = paste0(name_fld_img_base,name_exp,"/n",ngrid[choose_n],"/",quantity,".pdf")
      titolo_plot = ""
      nomi_exp = sapply(1:nrow(params_grid), function(idx) paste0(get_first3digits(name_exp,1),idx))
      pos = positions[[hh]]
      shift = shifts[[hh]]
      barplot_pos = c()
      for(l in 1:length(pos)){
        barplot_pos = c(barplot_pos, pos[l]+shift)
      }
        
      if(save_plots)
        pdf(nome_fig,width = 14,height = 6)
      par(mfrow = c(1,1), mar = c(3,3,2,1), bty = "l", mgp = c(2, 0.95, 0), cex = 2, las = 1 )
      plot(0,0,type = "n",xlim = c(0,21), ylim = ylim_plot,
           xlab = " ", ylab = paste0("RMSE (",ylabs[hh],")"), 
           xaxt = "n", yaxt = "n",
           main = titolo_plot,
           cex.main = 1, cex.lab = 1, cex.axis = 1)
      axis(side = 2, las = 1, cex.axis = 1)
      axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 1)
      grid(lty = 1,lwd = 1, col = "gray90" )
      for(i in 1:ncol(ExpPred[[hh]])){
        boxplot(ExpPred[[hh]][,i], at = barplot_pos[i], add = T, 
                col = mycol_list[[hh]][i], pch = 16, yaxt = "n", cex = 0.5)
      }
      legend("topleft",legends[[hh]], col = mycol_list[[hh]],
             pch = 16, cex = 1.25,  border = NA, bty = "n" )
      if(save_plots)
        dev.off()
    }
  }
  

  
  
  
}

# Additive -----------------------------------------------------------------

## a) Basic options
save_plots = FALSE
name_fld_img_base = "img/Paper/Exp1/Additive/"
name_folder_save = "save/Pred_Additive_"
nome_base1 = paste0("save/Pred_Additive_")


## Dirichlet parameters
name_exp = "Dirichelt"
p1 = 0.1; p2 = 0.5
params = c(p1,p2)
masses = c(0,0.5,1)
params_grid = expand.grid(params,params,masses)
idx_settings = lapply(1:4, function(j) seq(j,8+j,by = 4))
Lsettings = length(idx_settings)

Nrep <- BB <- 100
nmax <- 800
ngrid <- seq(100,nmax,by = 100)
Lgrid = length(ngrid)


i = 1
for(i in 1:Lsettings){ # for each experiment (A_i)
  idx_rep = i
  idxs = idx_settings[[i]]
  params_exp = params_grid[idxs,]
  
  choose_n = 1
  for(choose_n in 1:Lgrid){ # for each sample size 
    ## Define main obj to be plot
    ExpPred = vector("list",4)
    names(ExpPred) = c("Sh","K","K1","K2")
    ExpPred = lapply(ExpPred, function(x) matrix(ncol = 0, nrow = Nrep))
    Params_BayesI <- Params_Bayes2 <- rep(-1,4)
    
    idx = 1
    for(idx in 1:nrow(params_exp) ){ # for each setting
      p1 = params_exp[idx,1]; p2 = params_exp[idx,2]; mass = params_exp[idx,3]
      
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2," || mass = ",mass,"\n")
      
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
      trim_mass = get_first3digits(mass,4)
      
      filename = paste0(nome_base1,"Dirichelt","_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat")
      load(file = filename)
      
      # Compute RMSE
      ## --> Compute SHARED
      RMSE_Sh_Ind = sapply(results, function(x) {
        x_true  = na.omit(x$Stot_true)
        x_Ind   = na.omit(x$Stot_est)
        rmse = sqrt((x_true - x_Ind)^2)
        rmse
      } )
      RMSE_Sh_MLE = sapply(results, function(x) {
        x_true  = na.omit(x$Stot_true)
        x_MLE   = na.omit(x$Stot_est_mle)
        rmse = sqrt((x_true - x_MLE)^2)
        rmse
      } )
      
      ## --> Compute DISTINCT
      RMSE_K_Ind = sapply(results, function(x) {
        x_true  = na.omit(x$Ktot_true)
        x_Ind   = na.omit(x$Ktot_est)
        rmse = sqrt((x_true - x_Ind)^2)
        rmse
      } )
      RMSE_K_MLE = sapply(results, function(x) {
        x_true  = na.omit(x$Ktot_true)
        x_MLE   = na.omit(x$Ktot_est_mle)
        rmse = sqrt((x_true - x_MLE)^2)
        rmse
      } )
      RMSE_K_pooled = sapply(results, function(x) {
        x_true  = na.omit(x$Ktot_true)
        x_pooled   = na.omit(x$Kpooled_est)
        rmse = sqrt((x_true - x_pooled)^2)
        rmse
      } )
      
      ## --> Compute LOCAL DISTINCT - Area 1
      RMSE_K1_Ind = sapply(results, function(x) {
        x_true  = na.omit(x$K1tot_true)
        x_Ind   = na.omit(x$K1tot_est)
        rmse = sqrt((x_true - x_Ind)^2)
        rmse
      } )
      RMSE_K1_MLE = sapply(results, function(x) {
        x_true  = na.omit(x$K1tot_true)
        x_MLE   = na.omit(x$K1tot_est_mle)
        rmse = sqrt((x_true - x_MLE)^2)
        rmse
      } )
      RMSE_K1_indep1 = sapply(results, function(x) {
        x_true  = na.omit(x$K1tot_true)
        x_indep1 = na.omit(x$K1indep_est)
        rmse = sqrt((x_true - x_indep1)^2)
        rmse
      } )
      
      ## --> Compute LOCAL DISTINCT - Area 2
      RMSE_K2_Ind = sapply(results, function(x) {
        x_true  = na.omit(x$K2tot_true)
        x_Ind   = na.omit(x$K2tot_est)
        rmse = sqrt((x_true - x_Ind)^2)
        rmse
      } )
      RMSE_K2_MLE = sapply(results, function(x) {
        x_true  = na.omit(x$K2tot_true)
        x_MLE   = na.omit(x$K2tot_est_mle)
        rmse = sqrt((x_true - x_MLE)^2)
        rmse
      } )
      RMSE_K2_indep2 = sapply(results, function(x) {
        x_true  = na.omit(x$K2tot_true)
        x_indep2 = na.omit(x$K2indep_est)
        rmse = sqrt((x_true - x_indep2)^2)
        rmse
      } )
      
      # Update global object
      ExpPred$Sh = cbind(ExpPred$Sh,RMSE_Sh_Ind[,choose_n],RMSE_Sh_MLE[,choose_n])
      ExpPred$K  = cbind(ExpPred$K,RMSE_K_Ind[,choose_n],RMSE_K_MLE[,choose_n],RMSE_K_pooled[,choose_n])
      ExpPred$K1 = cbind(ExpPred$K1,RMSE_K1_Ind[,choose_n],RMSE_K1_MLE[,choose_n],RMSE_K1_indep1[,choose_n])
      ExpPred$K2 = cbind(ExpPred$K2,RMSE_K2_Ind[,choose_n],RMSE_K2_MLE[,choose_n],RMSE_K2_indep2[,choose_n])
      
      # Save parameters
      Params_BayesI = c(median(results[[choose_n]]$gamma1_mle),
                        median(results[[choose_n]]$gamma2_mle),
                        median(results[[choose_n]]$lambda_mle),
                        median(results[[choose_n]]$Ktot_true))
      Params_BayesII = c(median(results[[choose_n]]$gamma1_ind),
                         median(results[[choose_n]]$gamma2_ind),
                         median(results[[choose_n]]$lambda_ind),
                         median(results[[choose_n]]$Ktot_true))
      ## Compute and plot qM_star
      Params_Bayes = list(Params_BayesI,Params_BayesII)
      qM_list = vector("list",2)
      for(idx_bayes in 1:2){
        gamma_vec = c(Params_Bayes[[idx_bayes]][1], Params_Bayes[[idx_bayes]][2])
        Lambda = Params_Bayes[[idx_bayes]][3]
        Kn = Params_Bayes[[idx_bayes]][4]
        qM_vec = c()
        for(mm in 0:1000){
          temp = HSSM::log_qM_post(m = mm, prior = "Poisson", prior_param = list("lambda"=Lambda),
                                   k = Kn, n_j = c(ngrid[choose_n],ngrid[choose_n]), 
                                   gamma_j = gamma_vec, log_V = -Inf, M_max = 200)[1]
          qM_vec = c(qM_vec, exp(temp))
          if(sum(qM_vec) >= (0.99))
            break
        }
        qM_list[[idx_bayes]] = qM_vec
      }
      LqMpost = max(sapply(qM_list,length))
      ymax = max(sapply(qM_list,max));ymin = 0
      ylim_plot = c(ymin,ymax)
      xgrid = 0:(LqMpost-1)
      xpos <- xlabs <- round(seq(xgrid[1],tail(xgrid,1),length.out = 5))
      
      nome_fig = paste0(name_fld_img_base,"A",idx_rep,"/n",ngrid[choose_n],"/","Mstar_",idx,".pdf")
      titolo_plot = ""
      
      if(save_plots)
        pdf(nome_fig,width = 14,height = 6)
      par(mfrow = c(1,1), mar = c(3,4,2,1), bty = "l", mgp = c(3, 0.95, 0), cex = 2, las = 1 )
      plot(0,0,type = "n",xlim = range(xgrid), ylim = ylim_plot,
           xlab = " ", ylab = paste0("Mstar"), 
           xaxt = "n", yaxt = "n",
           main = titolo_plot,
           cex.main = 1, cex.lab = 1, cex.axis = 1)
      axis(side = 2, las = 1, cex.axis = 1)
      axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = 1)
      grid(lty = 1,lwd = 1, col = "gray90" )
      points(x = xgrid[1:length(qM_list[[1]])], 
             y = qM_list[[1]], col = "darkred", 
             type = "b", lwd = 3, pch = 16)
      points(x = xgrid[1:length(qM_list[[2]])], 
             y = qM_list[[2]], 
             col = "darkgreen", type = "b", lwd = 3, pch = 16)
      legend("topright",c("Bayes I","Bayes II"), col = c("darkred","darkgreen"),
             pch = 16, cex = 1.25,  border = NA, bty = "n" )
      if(save_plots)
        dev.off()
    }
    
    ## Set plotting positions
    pos_ii = pos3
    positions = list(pos_ii,pos_ii,pos_ii,pos_ii)
    shf1 = c(-0.25,0.25); shf2 = c(-0.5,0,0.5)
    shifts = list(shf1,shf2,shf2,shf2)
    
    ## Run plot
    hh = 1
    ylabs = c("S","K","K1","K2")
    for(hh in 1:4){ # for each quantity
      quantity = names(ExpPred)[hh]
      ymax = max(ExpPred[[hh]])
      ymin = 0
      ylim_plot = c(ymin,ymax)
      
      nome_fig = paste0(name_fld_img_base,"A",idx_rep,"/n",ngrid[choose_n],"/",quantity,".pdf")
      titolo_plot = ""
      nomi_exp = sapply(1:3, function(idx) paste0("A(",idx_rep,",",idx,")"))
      pos = positions[[hh]]
      shift = shifts[[hh]]
      barplot_pos = c()
      for(l in 1:length(pos)){
        barplot_pos = c(barplot_pos, pos[l]+shift)
      }
      
      if(save_plots)
        pdf(nome_fig,width = 14,height = 6)
      par(mfrow = c(1,1), mar = c(3,3,2,1), bty = "l", mgp = c(2, 0.95, 0), cex = 2, las = 1 )
      plot(0,0,type = "n",xlim = c(0,21), ylim = ylim_plot,
           xlab = " ", ylab = paste0("RMSE (",ylabs[hh],")"), 
           xaxt = "n", yaxt = "n",
           main = titolo_plot,
           cex.main = 1, cex.lab = 1, cex.axis = 1)
      axis(side = 2, las = 1, cex.axis = 1)
      axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 1)
      grid(lty = 1,lwd = 1, col = "gray90" )
      for(i in 1:ncol(ExpPred[[hh]])){
        boxplot(ExpPred[[hh]][,i], at = barplot_pos[i], add = T, 
                col = mycol_list[[hh]][i], pch = 16, yaxt = "n", cex = 0.5)
      }
      legend("topleft",legends[[hh]], col = mycol_list[[hh]],
             pch = 16, cex = 1.25,  border = NA, bty = "n" )
      if(save_plots)
        dev.off()
    }
  }
  

}

