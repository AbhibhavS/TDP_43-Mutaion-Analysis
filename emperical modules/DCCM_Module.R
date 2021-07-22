#---------------DCCM

setwd("K:/Abhibhav/TDP43")
library(tidyverse)
library(gplots)

system<-dir()[grep(".pdb", dir())]


my_palette <-colorRampPalette(c("red","white","blue"))(n = 400)

label<-list()

for(k in seq_along(system)){
  
  data<-read.table(system[k], sep="\t")
  data.protein<-data[grep(" CA ", data[,1]),] %>% str_squish %>% strsplit(" ")
  data.protein<-do.call(rbind, data.protein)
  label[[k]]<-data.protein[1:length(unique(data.protein[,6])), 4]
  
  start<- which(as.numeric(data.protein[,6])==1)
  end<- which(as.numeric(data.protein[,6])==length(unique(data.protein[,6])))
  
  snap<-list()
  
  for(i in seq_along(start)){
    
    current_snap <- data.protein[(start[i]):(end[i]),]
    snap[[i]] <- as.data.frame(current_snap)
    
  }
  
  snap_coord<- lapply(snap, function(x){cbind(as.numeric(x[,7]),
                                              as.numeric(x[,8]),
                                              as.numeric(x[,9]))})
  
  avg_r<- Reduce("+", snap_coord)/length(snap_coord)
  list_delta_r<-lapply(snap_coord, function(x){ x - avg_r})
  sqrt_l2_r<-sqrt(Reduce("+", lapply(list_delta_r,
                                     function(x){rowSums(x^2)}))/length(list_delta_r))
  
  protein_length<-nrow(snap_coord[[1]])
  
  DCCM<-matrix(0, protein_length, protein_length)
  
  for(i in 1:protein_length){
    for(j in 1:protein_length){
      
      numerator<- Reduce("+",lapply(list_delta_r, 
                                    function(x){ sum(x[i,]*x[j,])}))/length(list_delta_r)
      
      denominator<- sqrt_l2_r[i]*sqrt_l2_r[j]
      DCCM[i,j]<- numerator/denominator
      
    }
  }
  
  map<-DCCM
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Dynamic Cross Correlation",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = system[k],
            labCol = label[[k]],
            labRow = label[[k]])
  
  
}














