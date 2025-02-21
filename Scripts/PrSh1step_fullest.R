
setwd("C:/Users/colom/ScriptSpecies_shared/DatiTrieste")
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("TrainTest_functions.R")
source("../SimStudyRaf/SimStudyRaf_functions.R")

mycol_pal1 = hcl.colors(n=6,palette = "Zissou1")
mycol_pal2 = hcl.colors(n=6,palette = "Temps")
mycol = hcl.colors(n=4,palette = "Zissou1")
mycol = c("black", mycol[1:3], "darkred", "darkgreen")



# Load all data ---------------------------------------------------------------

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
data
data_norm = apply(data,2,function(x){x/sum(x)})
w_j = lapply(1:2, function(j){data_norm[,j]})
par(mfrow = c(1,2), mar = c(2,4,2,1), bty = "l")
plot(x = 1:K12_all, y = sort(w_j[[1]], decreasing = TRUE), type = "p", pch = 16)
plot(x = 1:K12_all, y = sort(w_j[[2]], decreasing = TRUE), type = "p", pch = 16)


# PrSh_1step --------------------------------------------------------------

## Curves (step by step est.) ------------------------------------------------------------------

# ngrid = c(25,seq(50,600,by = 50))
ngrid = c(seq(50,600,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species,
                                  gamma_j_est = NULL,
                                  lambda_est = NULL, estimate_param = TRUE)
# save(PrSh_ngrid, file = "save/PrSh1step_all.Rdat")

shift = c(0,5,10,10,-5,-10)
titolo_plot = paste0("")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.1), 
     ylab = "P(Snew > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}
legend("topright", c("True","Yue1","Yue2","Chao","Use prop.","Full est."),
       col = mycol[1:6], pch = 16)
