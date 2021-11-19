
#-----------Prortin DNA Map--------------------------------#


my.path<-"K:/TDP43_r2/"
setwd(my.path)

library(tidyverse)
library(Matrix)
library(dplyr)
library(magrittr)

source("MD_utility.R")
source("ML_utility.R")

system<-dir()[grep(".pdb", dir())]


files<-list()
for(kk in 1:3){
  
  print(paste("start", kk, "at", Sys.time(), sep=" "))
  
  data1<-read.table(system[4], sep="\t")
  data2<-read.table(system[kk], sep="\t")
  
  data<-rbind(data1, data2)
  
  xx<-mdml.matrix(data, type="ML")
  write.csv(xx,  paste0("distance_data_",system[kk],".csv", collapse = ""), row.names = F)
  
  files[[kk]]<-xx
  
}

#_________________________________________________________________#
#___________Global Feature_________________________________#

#----------Importing DATA------------------------#


my.files<- dir()[grep("distance", dir())][1:3] #reading the all interaction feature distance feature dataset

#----------Importing DATA------------------------#
wt<-read.csv(my.files[1], header = T)[1:2001,]
D169G<-read.csv(my.files[1], header = T)[2002:4002,]
Double_mut<-read.csv(my.files[2], header = T)[2002:4002,]
I168A<-read.csv(my.files[3], header = T)[2002:4002,]
data_mat<-rbind(wt, D169G, Double_mut, I168A)
data_mat<-filter(data_mat)

result<-ML(data_mat, 0.8, systems = 4, snaps = 2001)

#data_mat[1:100,1:100]

subDir<-"feature_imp"
dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
setwd(file.path(getwd(), subDir))
write.csv(result$system_feature[[3]], "Lasso_Global_imp_interaction_feature.csv", row.names = T)
write.csv(as.data.frame(result$system_feature[[2]][,1:2]), "Xgb_Global_imp_interaction_feature.csv", row.names = T)
write.csv(result$system_feature[[1]], "RF_Global_imp_interaction_feature.csv", row.names = T)
write.csv(result$system_feature[[4]], "Enet_Global_imp_interaction_feature.csv", row.names = T)
write.csv(result$system_feature[[5]], "PCA_Global_imp_interaction_feature.csv", row.names = T)




#-------------PLOT-----------------------------#
roc_plot(result$roc)

subDir<-"evaluation"
dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
setwd(file.path(getwd(), subDir))
#--------------------Model Evaluation--------------------------#
write.csv(conf.mat(result$tab_list), "model_met.csv", row.names = T)
write.csv(result$time, "model_time.csv", row.names = T)



#________________________Pairwise Features_____________________________#

setwd(my.path)
my.files<- dir()[grep("distance", dir())][1:3] #reading the all interaction feature distance feature dataset

mutation.system<-dir()[grep(".pdb", dir())]
mutation.system<-mutation.system[-grep("distance", mutation.system)]
mutation.system<-gsub(".pdb", "", mutation.system)
mutation.system<-mutation.system[-grep("WT_", mutation.system)]

for(i in 1:length(mutation.system)){

  setwd(my.path)
  data_mat<-read.csv(my.files[i], header = T)
  data_mat<-filter(data_mat)
  
  subDir<-"feature_imp"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
  setwd(file.path(getwd(), subDir))
  
  result<-ML(data_mat, 0.8, systems = 2, snaps = 2001, models = "ALL")
  #data_mat[1:100,1:100]
  
  write.csv(result$system_feature[[3]], paste0("Lasso_imp_int",system[i],".csv", collapse = ""), row.names = T)
  write.csv(as.data.frame(result$system_feature[[2]][,1:2]), paste0("XGB_imp_int",system[i],".csv", collapse = ""), row.names = T)
  write.csv(result$system_feature[[1]], paste0("RF_imp_int",system[i],".csv", collapse = ""), row.names = T)
  write.csv(result$system_feature[[4]], paste0("ENet_imp_int",system[i],".csv", collapse = ""), row.names = T)
  write.csv(result$system_feature[[5]], paste0("PCA_imp_int",system[i],".csv", collapse = ""), row.names = T)
  #write.csv(conf.mat(result$tab_list), paste0(system[i],"_model_met.csv", collapse = ""), row.names = T)
  roc_plot(result$roc)
  
}
