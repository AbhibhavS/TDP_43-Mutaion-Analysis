
#----------------------- contact Map------------------------------------#

setwd("D:/Abhibhav/TDP43")
library(tidyverse)

system<-dir()[grep(".pdb", dir())]

contact.calculator<- function(list_item){
  
  coord<-list_item[,7:9]
  contact_matrix<-matrix(0, nrow(coord), nrow(coord))
  
  for(i in 1:nrow(coord)){
    for(j in 1:nrow(coord)){
      
      contact_matrix[i,j]<- ifelse(sqrt(sum((as.numeric(coord[i,])-as.numeric(coord[j,]))^2))<=10,
                                   contact_matrix[i,j]+1, contact_matrix[i,j]+0)
      
    }
  }
  
  return(contact_matrix)
} #contact.calculator contact function


min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

contact_matrix<-list()
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
  #----------------------contact--------------------#
  contact_matrix[[k]] <-lapply(snap, contact.calculator)
  
}


contact.matrix<-list()

for(i in 1:4){ contact.matrix[[i]] <- Reduce("+", contact_matrix[[i]])/length(contact_matrix[[i]]) }
for(i in 1:4){ write.table(contact.matrix[[i]], paste0("contact_", system[i], ".txt", collapse = "")) }
names(contact.matrix)<- system


#-------------------------------------------------------------------#



#------------------------------------------------------------------------#

library(gplots)

my_palette <-colorRampPalette(c("White",
                                "black"))(n = 50)


for(k in 1:4){
  
  current_map<-contact.matrix[[k]]
  neig<-4
  
  for(i in 1:ncol(current_map)){
    
    if(i-neig < 0){
      current_map[1:(i+neig),i]<-0
      next
    }
    else if(i+neig > ncol(current_map)){
      current_map[(i-neig):ncol(current_map),i]<-0
      next
    }
    else{
      current_map[(i-neig):(i+neig),i]<-0
      next
    }
    
  }
  
  heatmap.2(current_map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "contact intensity",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = system[k],
            labCol = label[[k]],
            labRow = label[[k]])
  
  
}


#--------------------------------------------------------------------























