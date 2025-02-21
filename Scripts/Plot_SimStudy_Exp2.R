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

mycol = c("darkgreen","darkorange")
remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}

source("SimStudyRaf_functions.R")

name_fld_img = "G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/Script_HSSM/SimStudyRaf/img/NewParamEst_fullest/TrainTest_fullest/"
name_folder_save = "C:/Users/colom/ScriptSpecies_shared/SimStudyRaf/save/"
keeps = seq(0.1,0.9,by = 0.05)
# Yue 6060 ----------------------------------------------------------------

my_ylim = lapply(1:6, function(x){c(0,75)} )
my_ylim[1] = lapply(1:1, function(x){c(0,30)} )
my_ylim[2:4] = lapply(1:1, function(x){c(0,40)} )
my_ylim[5] = lapply(1:1, function(x){c(0,60)} )
alphas = c("08","085","09")#,"095")
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

# pdf("img/Exp2_all.pdf")
for(ii in 1:nrow(experiments)){
  
  name_exp_base  = "PC_res_fullestcor_Yue"
  name_exp_base2 = "PC_tabres_fullestcor_Yue"
  name_base_saveimg = "Yue6060_fullest_"
  name_exp = paste0(name_folder_save,name_exp_base,experiments[ii,1],experiments[ii,2],"_400.Rdat")
  name2_exp = paste0(name_folder_save,name_exp_base2,experiments[ii,1],experiments[ii,2],"_400.Rdat")
  load(file = name_exp);#load(file = name2_exp)
  
  
  
  
  ## plot
  ylim_plot = my_ylim[[ii]]
  nome_fig = paste0(name_fld_img,name_base_saveimg,experiments[ii,1],experiments[ii,2],"_400.pdf")
  # titolo_plot = paste0("Yue6060 - ",experiments[ii,1],"v",experiments[ii,2])
  titolo_plot <- substitute(
    val0 ~ alpha[1] == val1 ~ "," ~ alpha[2] == val2,
    list(val1 = experiments2[ii,1], val2 = experiments2[ii,2], val0 = scenario_names[[ii]])
  )
  S12_true = results[[1]]$Stot_true
  
  
  nomi_exp = as.character(keeps) # in decimale
  nomi_exp = as.character(keeps * 100) # in percentuale
  pos  = seq(1,2*length(nomi_exp),by = 2)
  
  if(save_plots)
    pdf(nome_fig,width = 8,height = 6)
  par(mfrow = c(1,1), mar = c(4,3,2,1), bty = "l", mgp = c(2.5, 0.95, 0) )
  plot(0,0,type = "n",xlim = c(0.5,2*length(nomi_exp)+0.5), ylim = ylim_plot,
       xlab = "% training set", ylab = " ", 
       xaxt = "n", yaxt = "n",
       main = titolo_plot,
       cex.main = 2, cex.lab = 2, cex.axis = 2)
  axis(side = 2, las = 1, cex.axis = 2)
  axis(side = 1, at = pos, labels = nomi_exp, las = 1, cex.axis = 2)
  grid(lty = 1,lwd = 1, col = "gray90" )
  abline(h = S12_true, lty = 4, col = "red")
  for(i in 1:length(results)){
    # boxplot(results[[i]]$Stot_noi,
    #         at = (2*i-1)-0.25, add = T, col = mycol_pal2[3], pch = 16)
    boxplot(results[[i]]$Stot_noi_est,
            at = (2*i-1)+0.25, add = T, col = mycol[1], pch = 1, 
            yaxt = "n")
    boxplot(results[[i]]$Stot_Chao,
            at = (2*i-1)-0.25, add = T, col = tail(mycol,1), pch = 1,
            yaxt = "n")
  }
  legend("topright",c("Chao2000","Proposed"), col = mycol[2:1], 
         pch = 16, cex = 1.5)
  if(save_plots)
    dev.off()
  
}
# dev.off()

# Sim model ----------------------------------------------------------------


my_ylim = c(0,150)

name_exp_base  = "PC_res_fullestcor_simmodel_400.Rdat"
name_exp_base2 = "PC_tabres_fullestcor_simmodel_400.Rdat"
name_base_saveimg = "simmodel_fullest_"
name_exp = paste0(name_folder_save,name_exp_base)
name2_exp = paste0(name_folder_save,name_exp_base2)
load(file = name_exp);load(file = name2_exp)




## plot
ylim_plot = my_ylim
nome_fig = paste0(name_fld_img,name_base_saveimg,"400.pdf")
titolo_plot = paste0("Sim. model")
S12_true = results[[1]]$Stot_true


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
  boxplot(results[[i]]$Stot_Chao,
          at = (2*i-1)+0.75, add = T, col = tail(mycol_pal2,1), pch = 16)
}
dev.off()








