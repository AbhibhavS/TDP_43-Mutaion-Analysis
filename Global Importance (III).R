


#setwd("C:/protein/TDP43/main/")
set.seed(123)

#--------Calling the required packages----------------
library(randomForest)
library(caret)
library(xgboost)
library(Matrix)
library(dplyr)
library(magrittr)
library(pROC)

protein2<-c("TSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNS") #TDP43 
protein2<-strsplit(protein2, "")[[1]]

DNA<-c("G","T","T","G","A","G","C","G","T","T") #DNA
sheet<-c(4:7,28:35, 42:50, 64:66, 69:74) #Sheet region in TDP43
helix<-c(14:23, 52:61) #Helix region in TDP43

#-------------Functions-----------
zscore<-function(x){ 
  (x-mean(x))/sqrt(var(x))
} #standardization


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
} #Normalization


#------------Intializing----------------------
system_feature<-list() #to save important features
#filtered_inter_distance_list<-list()
#filtered_inverse_inter_distance_list<-list()
#RF_intresting<-list()
#RF_imp_residue_list<-list()
tab_list<-list() #to save the confusion matix
prob<-list() #to save the probabilities



#--------------------------------DATA PREPARATION---------------------------------------------------#

#files<- dir()[grep(".csv", dir())][1:3] #reading the all interaction feature distance feature dataset

#----------Importing DATA------------------------#
wt_A169G<-read.csv(files[1], header = T)
Double_mut<-read.csv(files[2], header = T)[2002:4002,]
I168A<-read.csv(files[3], header = T)[2002:4002,]
data_mat<-rbind(wt_A169G, Double_mut, I168A)

#------------------Filtering---------------------#
less_than_10<-unique(which(data_mat <= 10, arr.ind = T)[,2])
greater_than_10<-unique(which(data_mat > 10, arr.ind = T)[,2])
filtered <-intersect(less_than_10, greater_than_10)
filtered_data_mat<-data_mat[,filtered]


#-----------------ML Splitting-----------------------#
part<-sort(sample(c(0,1), size = nrow(data_mat), prob = c(0.8,0.2), replace = T))

#----------inverse normalized (for RF and XGB)-----#
mydata<-(as.data.frame(apply(filtered_data_mat, 2, zscore)))^-1
mydata$label<-as.factor(sort(c(rep(c(0:3),2001))))
mydata_rand<-mydata[sample(1:nrow(mydata)),]
train<-mydata_rand[which(part==0),]
test<-mydata_rand[which(part==1),]

#----------inverse(for LASSO and ElasticNEt)-----#
mydata1<-(as.data.frame(filtered_data_mat)^-1)
mydata1$label<-as.factor(sort(c(rep(c(0:3),2001))))
mydata_rand1<-mydata1[sample(1:nrow(mydata1)),]
train1<-mydata_rand1[which(part==0),]
test1<-mydata_rand1[which(part==1),]

#-------------------------------------Modelling----------------------------#


#__________________Random Forest_________________________#

RF <- randomForest(label~.,
                   train, 
                   ntree = 200,
                   mtry = 8,
                   importance = TRUE,
                   proximity = TRUE,
                   do.trace=T)

# Prediction & Confusion Matrix - train data

pred <- predict(RF, test) # predicted label
pred1<-predict(RF, test, type="prob") #predicted probability
tab<-table(Predicted= pred, Actual_en= test$label) #confusion matrix
tab_list[[1]]<-tab
prob[[1]]<-pred1

imp<-as.data.frame(importance(RF, type=2, scale = F)) #importance

sorted_imp<-as.data.frame(imp[order(imp$MeanDecreaseGini, decreasing = T),])
row.names(sorted_imp)<-row.names(imp)[order(imp$MeanDecreaseGini, decreasing = T)]
names(sorted_imp)<-"importance"
system_feature[[1]]<-min_max(sorted_imp)*100



#_________________________Elastic Net_________________________

custom <- trainControl(method = "cv",
                       number = 3,
                       verboseIter = T)

ENet <- train(label~.,
              train1, 
              method = "glmnet",
              tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                     lambda= seq(0.0001, 1, length = 5)),
              trControl = custom)



pred1 <- predict(ENet, test1, type="prob") #predicted probability
pred <- predict(ENet, test1)  # predicted label
prob[[2]]<-pred1 #confusion matrix
tab<-table(Predicted= pred, Actual_en= test1$label)
tab_list[[2]]<-tab

imp<-varImp(ENet, scale=T) #importance
imp<-rowMeans(as.matrix(imp$importance))
imp<-as.data.frame(sort(imp[which(imp>0)], decreasing = T))
imp<-min_max(imp)*100
names(imp)<-"importance"
system_feature[[2]]<-imp


#__________________________LASSO_______________________________#


custom <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 5,
                       verboseIter = T)

lasso <- train(label~.,
               train1, 
               method = "glmnet",
               tuneGrid = expand.grid(alpha=1,
                                      lambda= seq(10^(-10), 1, length = 5)),
               trControl = custom)


pred1 <- predict(lasso, test1, type="prob") #predicted probability
pred <- predict(lasso, test1) # predicted label
prob[[3]]<-pred1 #confusion matrix
tab<-table(Predicted= pred, Actual_en= test1$label)
tab_list[[3]]<-tab

imp<-varImp(lasso, scale=T) #importance
imp<-rowMeans(as.matrix(imp$importance))
imp<-as.data.frame(sort(imp[which(imp>0)], decreasing = T))
imp<-min_max(imp)*100
names(imp)<-"importance"
system_feature[[3]]<-imp


#_____________________________XGBOOST____________________________#

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
tab_list[[4]]<-tab
prob[[4]]<-pred

imp <- xgb.importance(colnames(train_matrix), model = bst_model)
system_feature[[4]]<-imp


#----------------------------ROC AUC ----------------------------------------------

roc1<-roc(ifelse(test$label==0, 0,1),prob[[1]][,1], percent = T)
roc2<-roc(ifelse(test1$label==0, 0,1),prob[[2]][,1], percent = T)
roc3<-roc(ifelse(test1$label==0, 0,1),prob[[3]][,1], percent = T)
roc4<-roc(ifelse(test$label==0, 0,1),prob[[4]][,1], percent = T)

my_roc<-list(roc1, roc2, roc3, roc4)
model_name<-c("Random Forest", "Elastic Net", "LASSO", "XGBoost")


#-------------PLOT-----------------------------#
par(mfrow=c(2,2))
for(i in 1:4){
  
  plot(my_roc[[i]], 
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


#--------------------Model Evaluation--------------------------#

metric<-function(tab){
  TP<-tab[2,2]; FN<-tab[1,2]; TN<-tab[1,1]; FP<-tab[2,1]
  ACC = (sum(diag(tab)))/sum(tab) #Accuracy
  Sen =TP/(TP+FN) #Sensitivity
  Spe =TN/(TN+FP) #Specificity
  Pre = TP/(TP+FP) #precision
  #MCC = (TP*TN - FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP))
  return(c(ACC, Sen, Spe, Pre))
} #evaluation  metric function

names(tab_list)<-model_name
model_metric<-list()


for(k in 1:4){
tab1<-c(tab_list[[k]][1,1], 
        sum(tab_list[[k]][,1])-tab_list[[k]][1,1], 
        sum(tab_list[[k]][,2])-tab_list[[k]][2,2],
        tab_list[[k]][2,2])
tab1<-matrix(tab1,2,2,byrow=F)

tab2<-c(tab_list[[k]][1,1], 
        sum(tab_list[[k]][,1])-tab_list[[k]][1,1], 
        sum(tab_list[[k]][,3])-tab_list[[k]][3,3],
        tab_list[[k]][3,3])
tab2<-matrix(tab2,2,2,byrow=F)

tab3<-c(tab_list[[k]][1,1], 
        sum(tab_list[[k]][,1])-tab_list[[k]][1,1], 
        sum(tab_list[[k]][,4])-tab_list[[k]][4,4],
        tab_list[[k]][4,4])
tab3<-matrix(tab3,2,2,byrow=F)
        
model_metric[[k]]<-rbind(metric(tab1),metric(tab2),metric(tab3))
}

mod_met<- do.call(rbind, model_metric) #final evaluation table

#write.csv(mod_met, "model_metric.csv", row.names = F)

