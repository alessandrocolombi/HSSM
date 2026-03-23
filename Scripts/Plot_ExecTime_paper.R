wd = "C:/Users/colom/ScriptSpecies_shared/RevBA"
setwd(wd)


# Libraries --------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(iNEXT)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))


# Plot options ------------------------------------------------------------

mycol = c("darkred","darkgreen","darkred","darkgreen")
save_img = FALSE
width = 12; height = 6
cex.lab = 2
cex.axis = 2
cex = 2
cex.lgn = 2

# Custom functions --------------------------------------------------------
source("R/Rfunctions_aux.R")
source("C:/Users/colom/ScriptSpecies_shared/SimStudyRaf/SimStudyRaf_functions.R")

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}


# Paper - ExecTime - Dir -----------------------------------------------------------------

name_fld_img = "img/Paper/"
name_folder_save = "save/"
nmax = 800
M = 60
ngrid <- seq(100,nmax,by = 100)
Lgrid = length(ngrid)

experiments = list(
  "Dirichelt" = c(0.5,0.5)
)

nome_base1 = paste0("save/Pred_")

i = 1
name_exp = names(experiments)[i]
idx = 1
p1 = experiments[[i]][1]; p2 = experiments[[i]][2]
trim_p1 = get_first3digits(p1,4)
trim_p2 = get_first3digits(p2,4)
    
## Read
filename = paste0(nome_base1,name_exp,"_",trim_p1,"_",trim_p2,".Rdat") 
load(file = filename)

res_fit_time = lapply(results, function(x) matrix(c(x$time_fit_mle,x$time_fit_ind), nrow = 100, ncol = 2, byrow = FALSE) )
res_est_time = lapply(results, function(x) matrix(c(x$time_est_mle,x$time_est_ind), nrow = 100, ncol = 2, byrow = FALSE) )

## axis
xlim_plot = c(0,20)
xlabs = as.character(ngrid)
xpos = seq(1,19,length.out = length(ngrid))
ymax = 2.6;ymin = 0
ylim_plot = c(ymin,ymax)
ylabs = round(seq(ymin,ymax,length.out = 5),1)
titolo_plot = " "
    

## Boxplot positions
positions = xpos
shift = c(-0.25,0.25)
barplot_pos = vector("list", length = length(ngrid))
for(l in seq_along(positions)){
  barplot_pos[[l]] = positions[l]+shift
}

nome_fig = paste0(name_fld_img,name_exp,"_Time_est",".pdf")
if(save_img)
  pdf(nome_fig,width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,5,1,1), mgp=c(3.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,  yaxt = "n", xaxt = "n",
     xlab = "", ylab = "Time (sec.)",
     xlim = xlim_plot, ylim = ylim_plot,
     main = titolo_plot,
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, labels = ylabs, las = 1, cex.axis = cex.axis )
axis(1, at = xpos, labels = xlabs, cex.axis = cex.axis);mtext("n", side = 1, line = 2.5, cex = cex.axis)
for(h in 1:length(res_est_time)){
  for(i in 1:ncol(res_est_time[[h]])){
    boxplot(res_est_time[[h]][,i], 
            at = barplot_pos[[h]][i], col = mycol[i], 
            add = T, pch = 16, yaxt = "n", cex = 0.5)
  }
}
legend("topleft",c("Bayes I","Bayes II"), col = mycol,
       pch = 16, cex = cex.lgn,  border = NA, bty = "n" )
if(save_img)
  dev.off()


ymax = 0.4;ymin = 0
ylim_plot = c(ymin,ymax)
ylabs = round(seq(ymin,ymax,length.out = 5),1)
nome_fig = paste0(name_fld_img,name_exp,"_Time_fit",".pdf")
if(save_img)
  pdf(nome_fig,width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,5,1,1), mgp=c(3.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,  yaxt = "n", xaxt = "n",
     xlab = "", ylab = "Time (sec.)",
     xlim = xlim_plot, ylim = c(0,0.25),
     main = titolo_plot,
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, labels = ylabs, las = 1, cex.axis = cex.axis )
axis(1, at = xpos, labels = xlabs, cex.axis = cex.axis);mtext("n", side = 1, line = 2.5, cex = cex.axis)
for(h in 1:length(res_fit_time)){
  for(i in 1:ncol(res_fit_time[[h]])){
    boxplot(res_fit_time[[h]][,i], 
            at = barplot_pos[[h]][i], col = mycol[i], 
            add = T, pch = 16, yaxt = "n", cex = 0.5)
  }
}
legend("topleft",c("Bayes I","Bayes II"), col = mycol,
       pch = 16, cex = cex.lgn,  border = NA, bty = "n" )
if(save_img)
  dev.off()
    

