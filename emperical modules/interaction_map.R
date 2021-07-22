
setwd("D:/Abhibhav/TDP43")
library(tidyverse)

system<-dir()[grep(".pdb", dir())]

system_imp<-list()

for(k in seq_along(system)){
  
  data<-read.table(system[k], sep="\t")
  
  start<- grep("MODEL", data[,1])
  end<-grep("ENDMDL", data[,1])
  
  snap<-list()
  
  for(i in seq_along(start)){
    
    current_snap <- str_squish(data[(start[i]+1):(end[i]-2),1])
    suppressWarnings(
    current_snap <- do.call(rbind, strsplit(current_snap, " ")))
    snap[[i]] <- as.data.frame(current_snap[which(current_snap[,3]=="CA"),])
    
  }
  
  euclid<- function(list_item){
    
    coord<-list_item[,7:9]
    distance_matrix<-matrix(NA, nrow(coord), nrow(coord))
    
    for(i in 1:nrow(coord)){
      for(j in 1:nrow(coord)){
        
        distance_matrix[i,j]<- sqrt(sum((as.numeric(coord[i,])-as.numeric(coord[j,]))^2))
        
      }
    }
    distance <- distance_matrix[upper.tri(distance_matrix)]
    
    return(distance)
  } #euclid distance function
  min_max<- function(x){
    (x-min(x))/(max(x)-min(x))
  }
  
  
  #----------------------Distance--------------------#
  distance_matrix <-lapply(snap, euclid)
  coord <- snap[[1]]
  interaction<-matrix(NA, nrow(coord), nrow(coord))
  
  for(i in 1:nrow(coord)){
    for(j in 1:nrow(coord)){
      
      interaction[i,j]<- paste0(coord$V4[i],"( atom-",as.character(i),"-)",
                                "-", coord$V4[j],"(atom-",as.character(j),"-)", 
                                collapse = "")
      
    }
  }
  interaction <- interaction[upper.tri(interaction)]
  
  
  inter_distance<-do.call(rbind, distance_matrix)
  colnames(inter_distance)<-interaction
  
  #--------------------------------Filtering-------------------------------------------------------#
  
  less_than_10<-unique(which(inter_distance <= 10, arr.ind = T)[,2])
  greater_than_10<-unique(which(inter_distance > 10, arr.ind = T)[,2])
  filtered <-intersect(less_than_10, greater_than_10)
  filtered_inter_distance<-inter_distance[,filtered]
  
  #--------------------------------Feature Importance---------------------------------------------------#
  
  eig<-eigen(cov(filtered_inter_distance))
  feature_imp<-eig$vectors[,1:10] %*% eig$values[1:10]
  feature_imp<- min_max(feature_imp)
  interaction_feature_df<-as.data.frame(cbind(interaction[filtered], feature_imp))
  
  interaction_imp_map<- matrix(0, nrow(coord), nrow(coord))
  
  for(x in 1:nrow(interaction_feature_df)){
    
    insta<-suppressWarnings(as.numeric(strsplit(interaction_feature_df$V1[x], 
                                                "-")[[1]])[!is.na(as.numeric(strsplit(interaction_feature_df$V1[x], "-")[[1]]))])
    interaction_imp_map[insta[1], insta[2]]<- as.numeric(interaction_feature_df$V2[x])
    
    interaction_imp_map[insta[2], insta[1]]<-interaction_imp_map[insta[1], insta[2]]
    
  }
  
  my_palette <-colorRampPalette(c("white", "black"))(n = 30)
  heatmap.2(interaction_imp_map, trace="none", 
            dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.6, cexCol = 0.7, 
            key.xlab = "importance",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            labRow = paste(c(coord$V4), 1:78),
            labCol = paste(c(coord$V4), 1:78),
            main = system[k]
            
  )
}















