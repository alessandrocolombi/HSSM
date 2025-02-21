# Init --------------------------------------------------------------------
# setwd("/home/lucia.paci/Lucia/Ale/ScriptSpecies_shared/SimStudyRaf")
# setwd("~/Documents/Projects/ScriptSpecies_shared/SimStudyRaf")
setwd("C:/Users/colom/ScriptSpecies_shared/SimStudyRaf")
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
# suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
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

mycol = hcl.colors(n=4,palette = "Zissou1")
mycol = c("black", mycol[1:3], "darkred","darkgreen")

name_fld_img = "G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/Script_HSSM/SimStudyRaf/img/NewParamEst_fullest/PrSh1step_fullest/"
name_folder_save = "C:/Users/colom/ScriptSpecies_shared/SimStudyRaf/save/"
ngrid = c(25,seq(50,400,by = 50))
shift = c(0,4,8,8,-8,-4)
pchs  = c(17,rep(16,5))
# Yue 6060 ----------------------------------------------------------------

my_ylim = lapply(1:6, function(x){c(0,0.15)} )
my_ylim[5:6] = lapply(1:2, function(x){c(0,0.25)} )
my_xlim = c(45,405)
alphas = c("08","085","09")
experiments = expand.grid(alphas,alphas)
experiments = experiments[-c(2,3,6),]
colnames(experiments) = c("alpha1","alpha2"); row.names(experiments) = as.character(1:nrow(experiments))

alphas_val = c(0.8,0.85,0.9)
experiments2 = expand.grid(alphas_val,alphas_val)
experiments2 = experiments2[-c(2,3,6),]

scenario_names = vector("list",nrow(experiments))
scenario_names[[1]] = "Scenario I - "
scenario_names[[2]] = "Scenario IV - "
scenario_names[[3]] = "Scenario II - "
scenario_names[[4]] = "Scenario V - "
scenario_names[[5]] = "Scenario VI - "
scenario_names[[6]] = "Scenario III - "

save_plots = TRUE

# pdf("img/Exp1_all.pdf")
for(ii in 1:nrow(experiments)){
  
  name_exp_base  = "PrSh1_Yue6060_fullest_"
  name_base_saveimg = "Yue6060_fullest_"
  name_exp = paste0(name_folder_save,name_exp_base,experiments[ii,1],experiments[ii,2],".Rdat")
  load(file = name_exp)
  
  
  ## plot
  ylim_plot = my_ylim[[ii]]
  nome_fig = paste0(name_fld_img,name_base_saveimg,experiments[ii,1],experiments[ii,2],".pdf")
  # titolo_plot = paste0("Yue6060 - ",experiments[ii,1],"v",experiments[ii,2])
  titolo_plot <- substitute(
    val0 ~ alpha[1] == val1 ~ "," ~ alpha[2] == val2,
    list(val1 = experiments2[ii,1], val2 = experiments2[ii,2], val0 = scenario_names[[ii]])
  )
  if(save_plots)
    pdf(nome_fig,width = 8,height = 6)
  par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l", mgp = c(2.5, 0.95, 0) )
  plot(0,0,type = "n",
       main = titolo_plot,
       ylim = ylim_plot, xlim = my_xlim,
       xaxt = "n", yaxt = "n",
       xlab = "n", ylab = "", #ylab = "P(Snew > 0)",
       cex.main = 2, cex.lab = 2, cex.axis = 2)
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, las = 1, cex.axis = 1.75)
  axis(side = 1, at = ngrid[2:9], 
       labels = as.character(2 * ngrid[2:9]), las = 1, 
       cex.axis = 1.75 )
  for(kk in c(1:2,4,6)){
    segments( x0 = ngrid[2:9]+shift[kk], x1 = ngrid[2:9]+shift[kk],
              y0 = PrSh_ngrid[[1]][[kk]][1,2:9], y1 = PrSh_ngrid[[1]][[kk]][3,2:9],
              col = mycol[kk], lty = 1 )
    points(x = ngrid[2:9]+shift[kk], 
           y = PrSh_ngrid[[1]][[kk]][2,2:9], type = "b", 
           pch = pchs[kk], col = mycol[kk], lty = 2)
  }
  legend("topright", c("True","Yue","ChaoSh","Proposed"),
         pch = pchs[c(1:2,4,6)],
         col = mycol[c(1:2,4,6)], cex = 2)
  if(save_plots)
    dev.off()
  
  
}
# dev.off()
