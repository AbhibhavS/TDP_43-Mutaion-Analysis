
setwd("C:/protein/TDP43/main/inter_distance")
set.seed(123)

library(xgboost)
library(Matrix)
library(dplyr)
library(magrittr)

#files<- dir()[grep("distance", dir())][1:3]
files<- dir()[grep("inverse", dir())][1:3] #!!!!!!! check inverse

zscore<-function(x){
  (x-mean(x))/sqrt(var(x))
}

imp<-list()

for(i in 1:3){
  
  data<-read.csv(files[i], header = T)
  data<-apply(data,2, zscore)
  data<-as.data.frame(data)
  data$label<-c(rep(0,2001),rep(1,2001))
  data<-data[sample(1:nrow(data)),]
  
  
  #------------------Data Partition----------------#
  
  part<-sample(c(0,1), size=nrow(data), prob = c(0.8,0.2), replace = T)
  train<-data[part==0,]
  test<-data[part==1,]
  
  #-----------------Model Building----------------#
  
  trainm<-sparse.model.matrix(label ~ .-1, data = train)
  
  head(trainm)
  train_label <- train[,"label"]
  train_matrix <- xgb.DMatrix(data = as.matrix(trainm), label = train_label)
  
  testm <- sparse.model.matrix(label~.-1, data = test)
  test_label <- test[,"label"]
  test_matrix <- xgb.DMatrix(data = as.matrix(testm), label = test_label)
  
  
  #----------------Parameter------------------#
  nc <- length(unique(train_label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = nc)
  watchlist <- list(train = train_matrix, test = test_matrix)
  
  # ----------------eXtreme Gradient Boosting Model-----------------#
  
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = 1000,
                         watchlist = watchlist,
                         eta = 0.01,
                         max.depth = 3,
                         gamma = 0,
                         subsample = 1,
                         colsample_bytree = 1,
                         missing = NA,
                         seed = 333,
                         verbose = F)
  
  
  p <- predict(bst_model, newdata = test_matrix)
  pred <- matrix(p, nrow = nc, ncol = length(p)/nc) %>%
    t() %>%
    data.frame() %>%
    mutate(label = test_label, max_prob = max.col(., "last")-1)
  tab<-table(Prediction = pred$max_prob, Actual = pred$label)
  
  
  #----------------Feature Importance--------------------#
  
  # Feature importance
  imp[[i]] <- xgb.importance(colnames(train_matrix), model = bst_model)
  print(imp[[i]])
  #xgb.plot.importance(imp)
  
  #write.csv(imp, paste0("XGB_importance_",files[i],".csv", collapse = ""))
  
  print(tab)
  
}


setwd("C:/protein/TDP43/main/importance_inverse") #!!!!!! check for inverse
for(x in 1:3){
  
  #write.csv(lasso_imp_residue_list[[x]], 
  #paste0("lasso_important_feature_",system[x],".csv", collapse = ""), 
  #row.names = F)
  
  write.csv(imp[[x]][,1:2], 
            paste0("xgb_important_interaction_feature_inverse_",files[x],".csv", collapse = ""), 
            row.names = F)
  
}






