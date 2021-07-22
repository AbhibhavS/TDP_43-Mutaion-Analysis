
setwd("C:/protein/TDP43/main/inter_distance")
set.seed(123)

library(randomForest)
library(caret)


#files<- dir()[grep("distance", dir())][1:3]
files<- dir()[grep("inverse", dir())][1:3] #!!!!!!! check inverse

zscore<-function(x){
  (x-mean(x))/sqrt(var(x))
}


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

system_feature<-list()
filtered_inter_distance_list<-list()
filtered_inverse_inter_distance_list<-list()
RF_intresting<-list()
RF_imp_residue_list<-list()

tab_list<-list()

for(k in 1:3){
  
  #--------------------------------Feature Importance---------------------------------------------------#
  
  data_mat<-read.csv(files[k], header = T)
  
  mydata<-as.data.frame(apply(data_mat, 2, zscore))
  
  mydata$label<-as.factor(c(rep(0,2001), rep(1,2001)))
  mydata_rand<-mydata[sample(1:nrow(mydata)),]
  
  part<-sort(sample(c(0,1), size = nrow(mydata), prob = c(0.8,0.2), replace = T))
  train<-mydata_rand[which(part==0),]
  test<-mydata_rand[which(part==1),]
  
  model <- randomForest(label~.,
                        train, 
                        ntree = 300,
                        mtry = 8,
                        importance = TRUE,
                        proximity = TRUE,
                        do.trace=T)
  
  # Prediction & Confusion Matrix - train data
  
  pred <- predict(model, test)
  tab<-table(Predicted= pred, Actual_en= test$label)
  tab_list[[k]]<-tab
  
  imp<-as.data.frame(importance(model, type=2, scale = F))
  
  sorted_imp<-as.data.frame(imp[order(imp$MeanDecreaseGini, decreasing = T),])
  row.names(sorted_imp)<-row.names(imp)[order(imp$MeanDecreaseGini, decreasing = T)]
  names(sorted_imp)<-"importance"
  system_feature[[k]]<-min_max(sorted_imp)*100
  
  
}


setwd("C:/protein/TDP43/main/importance_inverse") #!!! check for the inverse
for(x in 1:3){
  
  #write.csv(RF_imp_residue_list[[x]], 
  #paste0("RF_important_feature_",system[x],".csv", collapse = ""), 
  #row.names = F)
  
  write.csv(system_feature[[x]], 
            paste0("RF_important_interaction_feature_inverse",files[x],".csv", collapse = ""), 
            row.names = T)
  
}

