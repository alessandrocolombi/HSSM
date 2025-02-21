# Init --------------------------------------------------------------------
# setwd("/home/lucia.paci/Lucia/Ale/ScriptSpecies_shared/SimStudyRaf")
# setwd("~/Documents/Projects/ScriptSpecies_shared/SimStudyRaf")
setwd("C:/Users/colom/ScriptSpecies_shared/SimStudyRaf")
# setwd("C:/Users/AC883875/Documents/ScriptSpecies_shared/SimStudyRaf")


suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(HSSM)))
# suppressWarnings(suppressPackageStartupMessages(library(SpadeR)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

# Functions ---------------------------------------------------------------

remove_points <- function(x) {
  str_replace_all(x, "\\.", "")
}
source("SimStudyRaf_functions.R")


# colors ------------------------------------------------------------------

mycol_pal1 = hcl.colors(n=6,palette = "Zissou1")
mycol_pal2 = hcl.colors(n=6,palette = "Temps")
mycol = hcl.colors(n=4,palette = "Zissou1")
mycol = c("black", mycol[1:3], "darkred","darkgreen")


# Yue 6060 - 8080 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.8
alpha2 = 0.8
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)

shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.10), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}



# Yue 6060 - 8585 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.85
alpha2 = 0.85
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)


shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.25), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}


# Yue 6060 - 9090 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.9
alpha2 = 0.9
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)


shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.25), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}


# Yue 6060 - 8085 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.80
alpha2 = 0.85
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)


shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.25), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}


# Yue 6060 - 8090 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.80
alpha2 = 0.90
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)


shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.25), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}


# Yue 6060 - 8590 - SbS ---------------------------------------------------------

seed0 = 271296
set.seed(seed0)

## Setup
S1 <- S2 <- 60
S12 <- 60
M = S12 + (S1-S12) + (S2-S12)
names_shared = sapply(1:S12, function(x){paste0("sh",as.character(x))})
names_distinct = c()
if(M>S12)
  names_distinct = sapply(1:(M-S12), function(x){paste0("di",as.character(x))})
names_species = c(names_shared,names_distinct)
counts = matrix(0,nrow = 2,ncol = M)
d = 2

# Geometric function
alpha1 = 0.85
alpha2 = 0.90
w_temp_list = vector("list",length = 2)
w_temp_list[[1]] = sapply(1:S1,function(i){alpha1^i})
w_temp_list[[2]] = sapply(1:S1,function(i){alpha2^i})
w_temp_list = lapply(w_temp_list, function(w){
  res = w[sample(1:S12,S12)]
  if(M>S12)
    res = c(res,w[(S12+1):(S1)])
  return(res)
})

w_1 = c(w_temp_list[[1]],rep(0,M-S1))
w_2 = c(w_temp_list[[2]][1:S12],rep(0,S2-S12))
if(M>S12)
  w_2 = c(w_2,w_temp_list[[2]][(S12+1):S2])
w_j = vector("list",2)
w_j[[1]] = w_1/sum(w_1); w_j[[2]] = w_2/sum(w_2)


## Curves (step by step est.) ------------------------------------------------------------------

ngrid = c(25,seq(50,400,by = 50))
Nrep  = 140
num_cores = 7
PrSh_ngrid = ProbSh_1step_fullest(w_j, ngrid, Nrep, num_cores, seed0, names_species)

name_save = paste0("save/PrSh1_Yue6060_fullest_",remove_points(alpha1),remove_points(alpha2),".Rdat")
save(PrSh_ngrid,file = name_save)


shift = c(0,1,3,3,-1,-2)
titolo_plot = paste0("Yue6060 - ",remove_points(alpha1),"v",remove_points(alpha2), "- full est.")
par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
plot(0,0,type = "n",xlim = c( min(ngrid)+min(shift), max(ngrid)+max(shift) ),
     main = titolo_plot, ylim = c(0,0.25), 
     ylab = "P(S_new > 0)", xlab = "n")
for(kk in 1:length(PrSh_ngrid[[1]])){
  points(x = ngrid+shift[kk], 
         y = PrSh_ngrid[[1]][[kk]][2,], type = "b", 
         pch = 16, col = mycol[kk])
  segments( x0 = ngrid+shift[kk], x1 = ngrid+shift[kk],
            y0 = PrSh_ngrid[[1]][[kk]][1,], y1 = PrSh_ngrid[[1]][[kk]][3,],
            col = mycol[kk], lty = 2 )
}

