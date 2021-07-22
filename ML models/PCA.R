
setwd("C:/protein/TDP43/main/inter_distance")
set.seed(123)

library(tidyverse)
library(caret)
library(gplots)

system<-dir()[grep(".pdb", dir())]
#files<- dir()[grep("distance", dir())][1:3]
files<- dir()[grep("inverse", dir())][1:3] #!!!!!!! check inverse



zscore<-function(x){
  (x-mean(x))/sqrt(var(x))
}


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}


system_feature_pca<-list()
pca_intresting<-list()
tab_list<-list()

pca_imp_int<-list()
pca_imp_residue_list<-list()

for(k in 1:3){
  
  #--------------------------------Feature Importance---------------------------------------------------#
  dis_data<-read.csv(files[k], header = T)
  pca_data<-dis_data[sample(1:nrow(dis_data)),]
  eig<-eigen(cov(pca_data))
  pca_imp_int[[k]]<-eig$vectors[,1:10] %*% eig$values[1:10]
  row.names(pca_imp_int[[k]])<-names(pca_data)
  
  pca_imp_int[[k]]<-abs(pca_imp_int[[k]])
  pca_imp_int[[k]]<-min_max(pca_imp_int[[k]])*100
  
  residue_number<-as.numeric(unlist(regmatches( names(pca_data), 
                                                gregexpr("[[:digit:]]+",  names(pca_data)))))
  pca_imp<-cbind(pca_imp_int[[k]], matrix(residue_number, ncol=2, byrow = T))
  
  pca_imp_residue<-matrix(1:length(sort(unique(residue_number))), nrow=length(unique(residue_number)),2)
  
  for(i in 1:nrow(pca_imp_residue)){
    pca_imp_residue[i,2]<-sum(pca_imp[which(pca_imp==pca_imp_residue[i,1], arr.ind = T)[,1],1])
  }
  pca_imp_residue[,2]<-min_max(pca_imp_residue[,2])
  pca_imp_residue[,1]<- pca_imp_residue[,1]+14
  pca_imp_residue_list[[k]]<-pca_imp_residue
}



setwd("C:/protein/TDP43/main/importance_inverse") #!!!! check for inverse
for(x in 1:3){

  #write.csv(pca_imp_residue_list[[x]], 
            #paste0("pca_important_feature_",system[x],".csv", collapse = ""), 
            #row.names = F)
  
  write.csv(pca_imp_int[[x]], 
            paste0("pca_important_interaction_feature_inverse",system[x],".csv", collapse = ""), 
            row.names = T)
  
  
}










