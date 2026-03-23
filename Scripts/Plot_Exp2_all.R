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
mycol_base = c("black","darkred","darkgreen","#E8A419","#3B99B1","deeppink","aquamarine")
mycol_list = list(rep(c(mycol_base[1:5]),6), 
                  rep(c(mycol_base[c(1:3,6)]),6),
                  rep(c(mycol_base[c(1:3,7)]),6),
                  rep(c(mycol_base[c(1:3,7)]),6))

legends = list(c("True","Bayes I","Bayes II","Chao","Yue"),
               c("True","Bayes I","Bayes II","Pooled"),
               c("True","Bayes I","Bayes II","Indep."),
               c("True","Bayes I","Bayes II","Indep."))

legend_pos = c("topright","topright","topright","topright")

pch_list = list( rep(c(17,16,16,16,16),6),
                 rep(c(17,16,16,16),6),
                 rep(c(17,16,16,16),6),
                 rep(c(17,16,16,16),6))

pos3 = seq(2,21,by = 7); pos6 = seq(2,21,by = 3.5)
# (D,Z,G) -----------------------------------------------------------------

## a) Basic options
save_plots = TRUE
name_fld_img_base = "img/Paper/Exp2/"
name_folder_save = "save/"
nome_base1 = paste0("save/Pr1step_New_")



experiments = list(
  "Dirichelt" = c(0.1,0.5),
  "Zipfs" = c(1.3,2),
  "Geom"  = c(0.8,0.85,0.9)
)

M = 60
Nrep <- BB <- 100
nmax <- 800
ngrid = seq(50,400,by = 50)
Lgrid = length(ngrid)


igrid = c(1:3)
i = 1
for(i in igrid){ # for each experiment (D,G,Z)
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  params = experiments[[i]]
  params_grid <- rbind( cbind(params, params), t(combn(params, 2)) )
  choose_n = 3
  for(choose_n in 1:Lgrid){ # for each sample size 
    ## Define main obj to be plot
    ExpPr1 = vector("list",4)
    names(ExpPr1) = c("Sh","K","K1","K2")
    ExpPr1 = lapply(ExpPr1, function(x) matrix(ncol = 0, nrow = 3))
    
    idx = 3
    for(idx in 1:nrow(params_grid) ){ # for each setting
      p1 = params_grid[idx,1]; p2 = params_grid[idx,2]
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2,"\n")
      
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
      filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat")
      load(file = filename)
      
      # Update global object
      ExpPr1$Sh = cbind(ExpPr1$Sh,Pr1step_ngrid[[1]]$Strue[,choose_n], Pr1step_ngrid[[1]]$ShMLE[,choose_n],Pr1step_ngrid[[1]]$ShInd[,choose_n], Pr1step_ngrid[[1]]$Chao[,choose_n],Pr1step_ngrid[[1]]$Yue1[,choose_n])
      ExpPr1$K  = cbind(ExpPr1$K, Pr1step_ngrid[[1]]$Ktrue[,choose_n], Pr1step_ngrid[[1]]$KMLE[,choose_n], Pr1step_ngrid[[1]]$KInd[,choose_n],  Pr1step_ngrid[[1]]$Kpooled[,choose_n])
      ExpPr1$K1 = cbind(ExpPr1$K1,Pr1step_ngrid[[1]]$K1true[,choose_n],Pr1step_ngrid[[1]]$K1MLE[,choose_n],Pr1step_ngrid[[1]]$K1Ind[,choose_n], Pr1step_ngrid[[1]]$K1indep[,choose_n])
      ExpPr1$K2 = cbind(ExpPr1$K2,Pr1step_ngrid[[1]]$K2true[,choose_n],Pr1step_ngrid[[1]]$K2MLE[,choose_n],Pr1step_ngrid[[1]]$K2Ind[,choose_n], Pr1step_ngrid[[1]]$K2indep[,choose_n])
      
    }
    
    ## Set plotting positions
    pos_ii = pos3
    if(name_exp == "Geom")
      pos_ii = pos6
    positions = list(pos_ii,pos_ii,pos_ii,pos_ii)
    shf1 = c(-0.8,-0.4,0,0.4,0.8); shf2 = c(-1,-0.5,0,0.5)
    shifts = list(shf1,shf2,shf2,shf2)
    
    ## Run plot
    ylabs = c("S","K","K1","K2")
    hh = 1
    for(hh in 1:4){
      quantity = names(ExpPr1)[hh]
      ymax = max(ExpPr1[[hh]])
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
      par(mfrow = c(1,1), mar = c(2,4,2,1), bty = "l", mgp = c(3, 0.95, 0), cex = 2, las = 1 )
      plot(0,0,type = "n",xlim = c(0,23.5), ylim = ylim_plot,
           xlab = " ", ylab = paste0("Prob. next new (",ylabs[hh],")"), 
           xaxt = "n", yaxt = "n",
           main = titolo_plot,
           cex.main = 1, cex.lab = 1)
      axis(side = 2, las = 1)
      axis(side = 1, at = pos, labels = nomi_exp, las = 1)
      grid(lty = 1,lwd = 1, col = "gray90" )
      for(i in 1:ncol(ExpPr1[[hh]])){
        points(x = barplot_pos[i], y = ExpPr1[[hh]][2,i], 
               pch = pch_list[[hh]][i], col = mycol_list[[hh]][i])
        segments( x0 = barplot_pos[i], x1 = barplot_pos[i],
                  y0 = ExpPr1[[hh]][1,i], y1 = ExpPr1[[hh]][3,i],
                  col = mycol_list[[hh]][i], lty = 1, lwd = 7 )
      }
      legend(legend_pos[hh],legends[[hh]], col = mycol_list[[hh]],
             pch = 16, cex = 1, border = NA, bty = "n")
      if(save_plots)
        dev.off()
    }
  
    
    
  } # End choose n
} # End name-experiment

# Additive -----------------------------------------------------------------

## a) Basic options
save_plots = TRUE
name_fld_img_base = "img/Paper/Exp2/Additive/"
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
ngrid = seq(50,400,by = 50)
Lgrid = length(ngrid)


i = 1
for(i in 1:Lsettings){ # for each experiment (A_i)
  idx_rep = i
  idxs = idx_settings[[i]]
  params_exp = params_grid[idxs,]
  
  choose_n = 1
  for(choose_n in 1:Lgrid){ # for each sample size 
    ## Define main obj to be plot
    ExpPr1 = vector("list",4)
    names(ExpPr1) = c("Sh","K","K1","K2")
    ExpPr1 = lapply(ExpPr1, function(x) matrix(ncol = 0, nrow = 3))
    
    idx = 1
    for(idx in 1:nrow(params_exp) ){ # for each setting
      p1 = params_exp[idx,1]; p2 = params_exp[idx,2]; mass = params_exp[idx,3]
      
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2," || mass = ",mass,"\n")
      
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
      trim_mass = get_first3digits(mass,4)
      
      filename = paste0(nome_base1,"Dirichelt","_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat")
      load(file = filename)
      
      
      # Update global object
      ExpPr1$Sh = cbind(ExpPr1$Sh,Pr1step_ngrid[[1]]$Strue[,choose_n], Pr1step_ngrid[[1]]$ShMLE[,choose_n],Pr1step_ngrid[[1]]$ShInd[,choose_n], Pr1step_ngrid[[1]]$Chao[,choose_n],Pr1step_ngrid[[1]]$Yue1[,choose_n])
      ExpPr1$K  = cbind(ExpPr1$K, Pr1step_ngrid[[1]]$Ktrue[,choose_n], Pr1step_ngrid[[1]]$KMLE[,choose_n], Pr1step_ngrid[[1]]$KInd[,choose_n],  Pr1step_ngrid[[1]]$Kpooled[,choose_n])
      ExpPr1$K1 = cbind(ExpPr1$K1,Pr1step_ngrid[[1]]$K1true[,choose_n],Pr1step_ngrid[[1]]$K1MLE[,choose_n],Pr1step_ngrid[[1]]$K1Ind[,choose_n], Pr1step_ngrid[[1]]$K1indep[,choose_n])
      ExpPr1$K2 = cbind(ExpPr1$K2,Pr1step_ngrid[[1]]$K2true[,choose_n],Pr1step_ngrid[[1]]$K2MLE[,choose_n],Pr1step_ngrid[[1]]$K2Ind[,choose_n], Pr1step_ngrid[[1]]$K2indep[,choose_n])
      
    }
    
    ## Set plotting positions
    pos_ii = pos3
    positions = list(pos_ii,pos_ii,pos_ii,pos_ii)
    shf1 = c(-0.8,-0.4,0,0.4,0.8); shf2 = c(-1,-0.5,0,0.5)
    shifts = list(shf1,shf2,shf2,shf2)
    
    ## Run plot
    ylabs = c("S","K","K1","K2")
    hh = 1
    for(hh in 1:4){
      quantity = names(ExpPr1)[hh]
      ymax = max(ExpPr1[[hh]])
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
      par(mfrow = c(1,1), mar = c(2,4,2,1), bty = "l", mgp = c(3, 0.95, 0), cex = 2, las = 1 )
      plot(0,0,type = "n",xlim = c(0,23.5), ylim = ylim_plot,
           xlab = " ", ylab = paste0("Prob. next new (",ylabs[hh],")"), 
           xaxt = "n", yaxt = "n",
           main = titolo_plot,
           cex.main = 1, cex.lab = 1)
      axis(side = 2, las = 1)
      axis(side = 1, at = pos, labels = nomi_exp, las = 1)
      grid(lty = 1,lwd = 1, col = "gray90" )
      for(i in 1:ncol(ExpPr1[[hh]])){
        points(x = barplot_pos[i], y = ExpPr1[[hh]][2,i], 
               pch = pch_list[[hh]][i], col = mycol_list[[hh]][i])
        segments( x0 = barplot_pos[i], x1 = barplot_pos[i],
                  y0 = ExpPr1[[hh]][1,i], y1 = ExpPr1[[hh]][3,i],
                  col = mycol_list[[hh]][i], lty = 1, lwd = 7 )
      }
      legend(legend_pos[hh],legends[[hh]], col = mycol_list[[hh]],
             pch = 16, cex = 1, border = NA, bty = "n")
      if(save_plots)
        dev.off()
    }
    
    
    
  }
}

