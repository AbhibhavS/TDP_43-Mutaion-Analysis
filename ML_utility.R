
zscore<-function(x){ 
  (x-mean(x))/sqrt(var(x))
} #standardization


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
} #Normalization


filter<-function(x){
  
  less_than_10<-unique(which(x <= 10, arr.ind = T)[,2])
  greater_than_10<-unique(which(x > 10, arr.ind = T)[,2])
  filtered <-intersect(less_than_10, greater_than_10)
  filtered_data_mat<-as.data.frame(x[,filtered])
  return(filtered_data_mat)
  
} #filtering the low influencing pairs


ml_split<-function(x, prop, z, systems, snaps){
  
  
  if(systems*snaps>nrow(x)){stop("Dataset has low number of snaps than provided")}
  x<-x[1:(systems*snaps),]
  part<-sort(sample(c(0,1), size = nrow(x), prob = c(prop,(1-prop)), replace = T)) #train, test proportion
  
  if(z==T){
    #----------inverse normalized (for RF and XGB)-----#
    mydata<-(as.data.frame(apply(x, 2, zscore)))^-1
    mydata$label<-as.factor(sort(c(rep(c(0:(systems-1)),snaps))))
    mydata_rand<-mydata[sample(1:nrow(mydata)),]
    train<-mydata_rand[which(part==0),]
    test<-mydata_rand[which(part==1),]
    a<-list(train, test)
    names(a)<-c("train", "test")
    return(a)
  }
  else{
    #----------inverse(for LASSO and ElasticNEt)-----#
    mydata<-(as.data.frame(x)^-1)
    mydata$label<-as.factor(sort(c(rep(c(0:(systems-1)),snaps))))
    mydata_rand<-mydata[sample(1:nrow(mydata)),]
    train<-mydata_rand[which(part==0),]
    test<-mydata_rand[which(part==1),]
    a<-list(train, test)
    names(a)<-c("train", "test")
    return(a)
  }
  
}


ML<-function(x, tr_size, systems, snaps, 
             models, ntree, mtry, xgb_params, lasso_method, enet_method){
  
  library(randomForest)
  library(caret)
  library(xgboost)
  library(Matrix)
  library(dplyr)
  library(magrittr)
  library(pROC)
  library(crayon)
  
  if(missing(models)){message(crayon::green$bold("models is missing: setting models to ALL as default")); models<-"ALL"}
  if(sum(models%in%c("ALL","RF","XGB","LASSO","ENet","PCA"))==0){stop("Invalid Model: try--> ALL, RF, XGB, LASSO, ENet, PCA")}
  
  
  ml_start_time <- Sys.time()
  #------------Intializing----------------------
  system_feature<-list() #to save important features
  tab_list<-list() #to save the confusion matix
  my_roc<-list() #to save the probabilities
  prob<-list()
  
  message(crayon::green$bold("lets crunch some data now..."))
  
  if(sum(models%in%c("RF","XGB","ALL","PCA"))>=1){
    
    b<-ml_split(x, prop = tr_size, z=T, systems = systems, snaps = snaps)
    train<-b$train
    test<-b$test
    if("PCA" %in% models | "ALL" %in% models){pca_data<-rbind(train, test)}
  }
  
  ite<-0
  name<-c()
  
  #__________________Random Forest_________________________#
  
  rf_start_time <- Sys.time()
  if("RF" %in% models | "ALL" %in% models ){
    message(crayon::green$bold("RF: Begins"))
    ite<-ite+1
    name<-c(name, "RF")
    
    if(missing(mtry)){
      message(crayon::red$bold("mtry is missing: setting mtry to 8 as default"))
      mtry<-8}
    
    if(missing(ntree)){
      message(crayon::red$bold("ntree is missing: setting ntree to 200 as default"))
      ntree<-200}
    
    RF <- randomForest(label~.,
                       train, 
                       ntree = ntree,
                       mtry = mtry,
                       importance = TRUE,
                       proximity = TRUE,
                       do.trace=T)
    
    # Prediction & Confusion Matrix - train data
    pred <- predict(RF, test) # predicted label
    pred1<-predict(RF, test, type="prob") #predicted probability
    tab<-table(Predicted= pred, Actual_en= test$label) #confusion matrix
    tab_list[[ite]]<-tab
    prob[[ite]]<-pred1
    my_roc[[ite]]<-roc(ifelse(test$label==0, 0,1),pred1[,1], percent = T)
    
    imp<-as.data.frame(importance(RF, type=2, scale = F)) #importance
    
    sorted_imp<-as.data.frame(imp[order(imp$MeanDecreaseGini, decreasing = T),])
    row.names(sorted_imp)<-row.names(imp)[order(imp$MeanDecreaseGini, decreasing = T)]
    names(sorted_imp)<-"importance"
    system_feature[[ite]]<-min_max(sorted_imp)*100
    message(crayon::yellow$bold("RF: Done"))
  }
  rf_end_time <- Sys.time()
  
  #_____________________________XGBOOST____________________________#
  
  xgb_start_time <- Sys.time()
  if("XGB" %in% models | "ALL" %in% models ){
    message(crayon::green$bold("XGBoost: Begins"))
    ite<-ite+1
    name<-c(name, "XGB")
    xgtrain<-train
    xgtrain$label<-as.numeric(train$label)-1
    
    xgtest<-test
    xgtest$label<-as.numeric(test$label)-1
    
    trainm<-sparse.model.matrix(label ~ .-1, data = xgtrain)
    train_label <- xgtrain[,"label"]
    train_matrix <- xgb.DMatrix(data = as.matrix(trainm), label = train_label)
    
    testm <- sparse.model.matrix(label~.-1, data = test)
    test_label <- xgtest[,"label"]
    test_matrix <- xgb.DMatrix(data = as.matrix(testm), label = test_label)
    
    
    #----------------Parameter------------------#
    
    if(missing(xgb_params)){
      message(crayon::red$bold("xgb_params is missing: setting:\n objective = multi:softprobe \n eval_metric = mlogloss\n nrounds = 1000, watchlist = watchlist\n eta = 0.01\n max.depth = 3 \n gamma = 0 \n subsample = 1 "))
      }
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
                           verbose = T)
    
    p <- predict(bst_model, newdata = test_matrix) #estimation of label
    pred <- matrix(p, nrow = nc, ncol = length(p)/nc) %>%
      t() %>%
      data.frame() %>%
      mutate(label = test_label, max_prob = max.col(., "last")-1) #prediction table
    tab<-table(Prediction = pred$max_prob, Actual = pred$label) #confusion matrix
    tab_list[[ite]]<-tab
    prob[[ite]]<-pred
    my_roc[[ite]]<-roc(ifelse(test$label==0, 0,1),pred[,1], percent = T)
    
    imp <- xgb.importance(colnames(train_matrix), model = bst_model)
    system_feature[[ite]]<-imp
    message(crayon::yellow$bold("XGBoost: Done"))
  }
  xgb_end_time <- Sys.time()
  
  #-------------------------------------------------------------------------
  
  if(sum(models%in%c("LASSO","ENet","ALL"))>=1){
    
    b<-ml_split(x, prop = tr_size, z=F, systems = systems, snaps = snaps)
    train<-b$train
    test<-b$test
  }
  
  #__________________________LASSO_______________________________#
  lasso_start_time <- Sys.time()
  if("LASSO" %in% models | "ALL" %in% models ){
    
    ite<-ite+1
    name<-c(name, "LASSO")
    message(crayon::green$bold("LASSO: Begins"))
    
    if(missing(lasso_method)){
      message(crayon::red$bold("lasso_method is missing: setting lasso_method to repeatedcv with n=10 and repeats =5 as default"))
      lasso_method<-"repeatedcv"}
    
    custom <- trainControl(method = lasso_method,
                           number = 10,
                           repeats = 5,
                           verboseIter = T)
    
    lasso <- train(label~.,
                   train, 
                   method = "glmnet",
                   tuneGrid = expand.grid(alpha=1,
                                          lambda= seq(10^(-10), 1, length = 5)),
                   trControl = custom)
    
    pred1 <- predict(lasso, test, type="prob") #predicted probability
    pred <- predict(lasso, test) # predicted label
    prob[[ite]]<-pred1 #confusion matrix
    my_roc[[ite]]<-roc(ifelse(test$label==0, 0,1),pred1[,1], percent = T)
    tab<-table(Predicted= pred, Actual_en= test$label)
    tab_list[[ite]]<-tab
    
    imp<-varImp(lasso, scale=T) #importance
    imp<-rowMeans(as.matrix(imp$importance))
    imp<-as.data.frame(sort(imp[which(imp>0)], decreasing = T))
    imp<-min_max(imp)*100
    names(imp)<-"importance"
    system_feature[[ite]]<-imp
    message(crayon::yellow$bold("LASSO: Done"))
  }
  lasso_end_time <- Sys.time()
  #_________________________Elastic Net_________________________
  
  enet_start_time <- Sys.time()
  if("ENet" %in% models | "ALL" %in% models ){
    
    ite<-ite+1
    name<-c(name, "ENet")
    message(crayon::green$bold("ENet: Begins"))
    
    if(missing(enet_method)){
      message(crayon::red$bold("enet_method is missing: setting enet_method to cv with n=3 as default"))
      enet_method<-"cv"}
    
    custom <- trainControl(method = enet_method,
                           number = 3,
                           verboseIter = T)
    
    ENet <- train(label~.,
                  train, 
                  method = "glmnet",
                  tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                         lambda= seq(0.0001, 1, length = 5)),
                  trControl = custom)
    
    
    
    pred1 <- predict(ENet, test, type="prob") #predicted probability
    pred <- predict(ENet, test)  # predicted label
    prob[[ite]]<-pred1 #confusion matrix
    my_roc[[ite]]<-roc(ifelse(test$label==0, 0,1),pred1[,1], percent = T)
    tab<-table(Predicted= pred, Actual_en= test$label)
    tab_list[[ite]]<-tab
    
    imp<-varImp(ENet, scale=T) #importance
    imp<-rowMeans(as.matrix(imp$importance))
    imp<-as.data.frame(sort(imp[which(imp>0)], decreasing = T))
    imp<-min_max(imp)*100
    names(imp)<-"importance"
    system_feature[[ite]]<-imp
    message(crayon::yellow$bold("ENet: Done"))
  }
  enet_end_time <- Sys.time()
  
  #_____________________________________PCA_______________________________________
  
  pca_start_time <- Sys.time()
  if("PCA" %in% models | "ALL" %in% models ){
    
    message(crayon::yellow$bold("PCA: Begins"))
    ite<-ite+1
    name<-c(name, "PCA")
    eig<-eigen(cov(pca_data[-ncol(pca_data)]))
    system_feature[[ite]]<-eig$vectors[,1:10] %*% eig$values[1:10]
    row.names(system_feature[[ite]])<-names(pca_data)[-length(pca_data)]
    system_feature[[ite]]<-abs(system_feature[[ite]])
    system_feature[[ite]]<-min_max(system_feature[[ite]])*100
    message(crayon::green$yellow("PCA: Done"))
    
  }
  pca_end_time <- Sys.time()
  
  ml_end_time <- Sys.time()
  
  run.time<-c(ml_end_time-ml_start_time,
              rf_end_time-rf_start_time,
              xgb_end_time-xgb_start_time,
              lasso_end_time-lasso_start_time,
              enet_end_time-enet_start_time,
              pca_end_time-pca_start_time)
  names(run.time)<-c("total","RF","XGB","LASSO","ENet","PCA")
  
  names(system_feature)<-name
  
  if(length(my_roc)>=1){
    names(my_roc)<-name[-length(name)]
    names(tab_list)<-name[-length(name)]}
  
  xx<-list(system_feature, tab_list, my_roc, run.time)
  names(xx)<-c("system_feature", "tab_list", "roc", "time")
  message(crayon::italic(crayon::white$bold("\n ....... \n All done \n.......\n bye bye")))
  return(xx)
  
}


#---roc plot function---

roc_plot<-function(list.roc){
  library(pROC)
  model_name<-c("Random Forest", "XGBoost", "LASSO", "Elastic Net")
  
  par(mfrow=c(2,2))
  for(i in 1:length(model_name)){
    
    plot(list.roc[[i]], 
         grid=T,
         xlab="False Positive Rate",
         ylab="True Positive Rate",
         main= model_name[i],
         #xaxt="n", yaxt="n",
         cex.lab=1.3,
         font.lab=9,
         print.auc.y=40,
         print.auc.x=24, 
         print.auc=T,
         legacy.axes=T,
         percent = T,
         col="red",
         lwd=3)
    
  }
}


#evaluation
metric<-function(tab){
  TP<-tab[2,2]; FN<-tab[1,2]; TN<-tab[1,1]; FP<-tab[2,1]
  ACC = (sum(diag(tab)))/sum(tab) #Accuracy
  Sen =TP/(TP+FN) #Sensitivity
  Spe =TN/(TN+FP) #Specificity
  Pre = TP/(TP+FP) #precision
  #MCC = (TP*TN - FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP))
  return(c(ACC, Sen, Spe, Pre))
} #evaluation  metric function


#Confusion Matrix
conf.mat<-function(list.tab){
  
  model_metric<-list()
  
  for(k in 1:4){
    tab1<-c(list.tab[[k]][1,1], 
            sum(list.tab[[k]][,1])-list.tab[[k]][1,1], 
            sum(list.tab[[k]][,2])-list.tab[[k]][2,2],
            list.tab[[k]][2,2])
    tab1<-matrix(tab1,2,2,byrow=F)
    
    tab2<-c(list.tab[[k]][1,1], 
            sum(list.tab[[k]][,1])-list.tab[[k]][1,1], 
            sum(list.tab[[k]][,3])-list.tab[[k]][3,3],
            list.tab[[k]][3,3])
    tab2<-matrix(tab2,2,2,byrow=F)
    
    tab3<-c(list.tab[[k]][1,1], 
            sum(list.tab[[k]][,1])-list.tab[[k]][1,1], 
            sum(list.tab[[k]][,4])-list.tab[[k]][4,4],
            list.tab[[k]][4,4])
    tab3<-matrix(tab3,2,2,byrow=F)
    
    model_metric[[k]]<-rbind(metric(tab1),metric(tab2),metric(tab3))
  }
  
  mod_met<- do.call(rbind, model_metric) #final evaluation table
  return(mod_met)
  #write.csv(mod_met, paste0(system[i],"model_metric", collapse = ""), row.names = F)
  
}


c(1,2,3)%in%c(1)

