wd = "C:/Users/colom/ScriptSpecies_shared/RevBA"
setwd(wd)


# Libraries --------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source("../SimStudyRaf/SimStudyRaf_functions.R")

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}

# Colors and legends ------------------------------------------------------


# Image options -----------------------------------------------------------

save_img = FALSE
width = 8; height = 6
cex.lab = 2
cex.axis = 2
cex = 2

# Figura Paper (D3-G6-Z3-A14 - A24 - A34) -----------------------------------------------------------------

## a) Basic options
choose_n = 4 # 4 is for n = 400
Nrep <- BB <- 100

name_fld_img = "img/Paper/Morosita"
name_folder_save = "save/"
nome_base1 = paste0("save/Pred_")

## b) D,Z,G: generate main obj
Mor_all = matrix(0,nrow = Nrep, ncol = 6)
names_experiments = c("Dirichelt","Geom","Zipfs")
params = list(c(0.5,0.5),c(0.9,0.9),c(2,2))
ii = 1
for(ii in 1:3){
  # Read
  p1 = params[[ii]][1]; p2 = params[[ii]][2]; trim_p1 = get_first3digits(p1,4); trim_p2 = get_first3digits(p2,4)
  filename = paste0(nome_base1,names_experiments[ii],"_",trim_p1,"_",trim_p2,".Rdat")
  load(file = filename)
  # Get Morisita values
  Mor_temp = sapply(results, function(x) x$Morosita )
  # Update global object
  Mor_all[,ii] = Mor_temp[,4]
}

## c) Additive: generate main obj
name_folder_save = "save/Pred_Additive_"
nome_base1 = paste0("save/Pred_Additive_")
masses = c(0,0.5,1)
p1 = 0.5; p2 = 0.5
trim_p1 = get_first3digits(p1,4); trim_p2 = get_first3digits(p2,4)

ii = 1
for(ii in 1:3){
  # Read
  mass = masses[ii]
  trim_mass = get_first3digits(mass,4)
  filename = paste0(nome_base1,"Dirichelt","_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat")
  load(file = filename)
  # Get Morisita values
  Mor_temp = sapply(results, function(x) x$Morosita )
  # Update global object
  Mor_all[,ii+3] = Mor_temp[,4]
}

## Set plotting positions
pos1 = seq(2,21,by = 3.5)
positions = pos1
## Run plot
ylabs = c("Morosita")
ymax = max(Mor_all)
ymin = 0
ylim_plot = c(ymin,ymax)
  
nome_fig = paste0(name_fld_img,".pdf")
titolo_plot = ""
  
nomi_exp = c(paste0("D",3),paste0("G",6),paste0("Z",3),
             paste0("A",14),paste0("A",24),paste0("A",34))
pos = positions
barplot_pos = c()
for(l in 1:length(pos)){
  barplot_pos = c(barplot_pos, pos[l])
}
  
if(save_img)
  pdf(nome_fig, width = width, height = height)
par(mfrow = c(1,1), mar = c(2,3,2,1), bty = "l", mgp = c(2, 0.95, 0), cex = 2, las = 1 )
plot(0,0,type = "n",xlim = c(0,21), ylim = ylim_plot,
     xlab = " ", ylab = ylabs, 
     xaxt = "n", yaxt = "n",
     main = titolo_plot)
axis(side = 2, las = 1, cex.axis = 1)
axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 1)
grid(lty = 1,lwd = 1, col = "gray90" )
for(i in 1:ncol(Mor_all)){
  boxplot(Mor_all[,i], at = barplot_pos[i], add = T, 
          col = "grey15", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()

# Figure Extra (all) -----------------------------------------------------------------

## Colors and legends ------------------------------------------------------------------
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
## (D,Z,G) -----------------------------------------------------------------

## a) Basic options
save_img = TRUE
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
  
  Mor_all = matrix(0,nrow = Nrep, ncol = nrow(params_grid))
  choose_n = 4
  for(choose_n in 1:Lgrid){ # for each sample size 
    idx = 3
    for(idx in 1:nrow(params_grid) ){ # for each setting
      p1 = params_grid[idx,1]; p2 = params_grid[idx,2]
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2,"\n")
        
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
      filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat")
      load(file = filename)
      # Get Morisita values
      Mor_temp = sapply(results, function(x) x$Morosita )
      Mor_all[,idx] = Mor_temp[,choose_n]
    }
    
    ## Set plotting positions
    pos_ii = pos3
    if(name_exp == "Geom")
      pos_ii = pos6
    ## Run plot
    ylabs = c("Morosita")
    ymax = max(Mor_all)
    ymin = 0
    ylim_plot = c(ymin,ymax)
    
    nome_fig = paste0(name_fld_img_base,name_exp,"/n",ngrid[choose_n],"/","Morosita",".pdf")
    titolo_plot = ""
    pos = pos_ii
    nomi_exp = sapply(1:nrow(params_grid), function(idx) paste0(get_first3digits(name_exp,1),idx))
    barplot_pos = c()
    for(l in 1:length(pos)){
      barplot_pos = c(barplot_pos, pos[l])
    }
    
    if(save_img)
      pdf(nome_fig, width = width, height = height)
    par(mfrow = c(1,1), mar = c(2,3,2,1), bty = "l", mgp = c(2, 0.95, 0), cex = 2, las = 1 )
    plot(0,0,type = "n",xlim = c(0,21), ylim = ylim_plot,
         xlab = " ", ylab = ylabs, 
         xaxt = "n", yaxt = "n",
         main = titolo_plot)
    axis(side = 2, las = 1, cex.axis = 1)
    axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 1)
    grid(lty = 1,lwd = 1, col = "gray90" )
    for(i in 1:ncol(Mor_all)){
      boxplot(Mor_all[,i], at = barplot_pos[i], add = T, 
              col = "grey15", pch = 16, yaxt = "n", cex = 0.5)
    }
    if(save_img)
      dev.off()
  }
}

# Additive -----------------------------------------------------------------

## a) Basic options
save_img = TRUE
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
  Mor_all = matrix(0,nrow = Nrep, ncol = nrow(params_exp))
  choose_n = 1
  for(choose_n in 1:Lgrid){ # for each sample size 
    
    idx = 1
    for(idx in 1:nrow(params_exp) ){ # for each setting
      p1 = params_exp[idx,1]; p2 = params_exp[idx,2]; mass = params_exp[idx,3]
      
      cat("\n",idx,"/",nrow(params_grid)," || p1 = ",p1," vs p2 = ",p2," || mass = ",mass,"\n")
      
      trim_p1 = get_first3digits(p1,4)
      trim_p2 = get_first3digits(p2,4)
      trim_mass = get_first3digits(mass,4)
      
      filename = paste0(nome_base1,"Dirichelt","_",trim_p1,"_",trim_p2,"_",trim_mass,".Rdat")
      load(file = filename)
      # Get Morisita values
      Mor_temp = sapply(results, function(x) x$Morosita )
      Mor_all[,idx] = Mor_temp[,choose_n]
    }
    
    ## Set plotting positions
    pos_ii = pos3
    ## Run plot
    ylabs = c("Morosita")
    ymax = max(Mor_all)
    ymin = 0
    ylim_plot = c(ymin,ymax)
    
    nome_fig = paste0(name_fld_img_base,"A",idx_rep,"/n",ngrid[choose_n],"/","Morosita",".pdf")
    titolo_plot = ""
    pos = pos_ii
    nomi_exp = sapply(1:3, function(idx) paste0("A(",idx_rep,",",idx,")"))
    barplot_pos = c()
    for(l in 1:length(pos)){
      barplot_pos = c(barplot_pos, pos[l])
    }
    
    if(save_img)
      pdf(nome_fig, width = width, height = height)
    par(mfrow = c(1,1), mar = c(2,3,2,1), bty = "l", mgp = c(2, 0.95, 0), cex = 2, las = 1 )
    plot(0,0,type = "n",xlim = c(0,21), ylim = ylim_plot,
         xlab = " ", ylab = ylabs, 
         xaxt = "n", yaxt = "n",
         main = titolo_plot)
    axis(side = 2, las = 1, cex.axis = 1)
    axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 1)
    grid(lty = 1,lwd = 1, col = "gray90" )
    for(i in 1:ncol(Mor_all)){
      boxplot(Mor_all[,i], at = barplot_pos[i], add = T, 
              col = "grey15", pch = 16, yaxt = "n", cex = 0.5)
    }
    if(save_img)
      dev.off()
    
  }
  
  
}


