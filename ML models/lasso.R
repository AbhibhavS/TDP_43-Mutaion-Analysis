

setwd("C:/protein/TDP43/main/inter_distance")
set.seed(123)

library(tidyverse)
library(caret)


#files<- dir()[grep("distance", dir())][1:3]
files<- dir()[grep("inverse", dir())][1:3] #!!!!!! check if inverse or not

zscore<-function(x){
  (x-mean(x))/sqrt(var(x))
}


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

system_feature<-list()
filtered_inter_distance_list<-list()
filtered_inverse_inter_distance_list<-list()
lasso_intresting<-list()
lasso_imp_residue_list<-list()

tab_list<-list()

for(k in 1:3){
  
  #--------------------------------Feature Importance---------------------------------------------------#
  
  data_mat<-read.csv(files[k], header = T)
  
  #mydata<-as.data.frame(apply(data_mat, 2, zscore))
  mydata<-as.data.frame(data_mat)
  
  mydata$label<-as.factor(c(rep(0,2001), rep(1,2001)))
  mydata_rand<-mydata[sample(1:nrow(mydata)),]
  
  part<-sort(sample(c(0,1), size = nrow(mydata), prob = c(0.8,0.2), replace = T))
  train<-mydata_rand[which(part==0),]
  test<-mydata_rand[which(part==1),]
  
  
  custom <- trainControl(method = "repeatedcv",
                         number = 10,
                         repeats = 5,
                         verboseIter = T)
  
  lasso <- train(label~.,
                 train, 
                 method = "glmnet",
                 tuneGrid = expand.grid(alpha=1,
                                        lambda= seq(10^(-10), 1, length = 5)),
                 trControl = custom)
  
  
  pred <- predict(lasso, test)
  tab<-table(Predicted= pred, Actual_en= test$label)
  tab_list[[k]]<-tab
  
  imp<-varImp(lasso, scale=T)
  imp<-as.matrix(imp$importance)
  system_feature[[k]]<-sort(imp[which(imp>0),], decreasing = T)
  
  residue_number<-regmatches( names(system_feature[[k]]), 
                              gregexpr("[[:digit:]]+",  names(system_feature[[k]])))
  lasso_intresting[[k]]<- sort(unique(as.numeric(unlist(residue_number))+14))
  
  aa<-cbind(system_feature[[k]], matrix(as.numeric(unlist(regmatches( names(system_feature[[k]]), 
                                                                      gregexpr("[[:digit:]]+",  
                                                                               names(system_feature[[k]]))))), 
                                        ncol=2, byrow = T))
  
  residue_imp<-c()
  for(i in 1:length(lasso_intresting[[k]])){
    residue_imp[i]<-sum(aa[which(aa==lasso_intresting[[k]][i]-14, arr.ind = T)[,1],1])
  }
  
  lasso_imp_residue_list[[k]]<-cbind(lasso_intresting[[k]], min_max(residue_imp))

}

setwd("C:/protein/TDP43/main/importance_inverse") # !!!!!! check where is it saving
for(x in 1:3){
  
  #write.csv(lasso_imp_residue_list[[x]], 
            #paste0("lasso_important_feature_",system[x],".csv", collapse = ""), 
            #row.names = F)
  
  write.csv(system_feature[[x]], 
            paste0("lasso_important_interaction_feature_inverse_",files[x],".csv", collapse = ""), 
            row.names = T)
  
}

