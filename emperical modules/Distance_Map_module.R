
#----------------------- Distance Map------------------------------------#

setwd("D:/Abhibhav/TDP43")
library(tidyverse)

system<-dir()[grep(".pdb", dir())]

euclid<- function(list_item){
  
  coord<-list_item[,7:9]
  distance_matrix<-matrix(NA, nrow(coord), nrow(coord))
  
  for(i in 1:nrow(coord)){
    for(j in 1:nrow(coord)){
      
      distance_matrix[i,j]<- sqrt(sum((as.numeric(coord[i,])-as.numeric(coord[j,]))^2))
      
    }
  }
  
  return(distance_matrix)
} #euclid distance function


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

distance_matrix<-list()
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
  
  #snap[[1]]
  #----------------------Distance--------------------#
  distance_matrix[[k]] <-lapply(snap, euclid)
  
}


distance.matrix<-list()

for(i in 1:4){ distance.matrix[[i]] <- Reduce("+", distance_matrix[[i]])/length(distance_matrix[[i]]) }
for(i in 1:4){ write.table(distance.matrix[[i]], paste0("distance_", system[i], ".txt", collapse = "")) }
names(distance.matrix)<- system

#------------------------------------------------------------------------#

library(gplots)

my_palette <-colorRampPalette(c("black","red","orange",
                                "blue","skyblue","white"))(n = 400)


for(i in 1:4){
  
  map<-distance.matrix[[i]]
  
  
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Avg. Distance (Angstrom)",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = system[i],
            labCol = label[[i]],
            labRow = label[[i]])
  
  
}


#--------------------------------------------------------------------


for(i in 1:3){
  map<-distance.matrix[[i]]- distance.matrix[[4]]
  
  my_palette <-colorRampPalette(c("red",
                                  "white","darkblue"))(n = 100)
  
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Avg. Distance Diff. (Angstrom)",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = paste(system[i],"-", system[4]),
            labCol = paste(label[[i]],"-",label[[4]]),
            labRow = paste(label[[i]],"-",label[[4]]))
  
}




















