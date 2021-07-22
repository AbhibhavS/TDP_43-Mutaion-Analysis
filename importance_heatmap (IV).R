
library(gplots)

setwd("C:/protein/TDP43/main/importance_inverse") #!!!check for inverse

files<-dir()[grep("csv", dir())]
system_mat<-matrix(files, ncol=4, byrow=T)
systems<-c("Global", "D169G", "D169G-I168A", "I168A")
mut<-matrix(c(67,66,67,67,67,66,66,66), ncol=4, byrow = F)

sheet<-c(4:7,28:35, 42:50, 64:66, 69:74)
helix<-c(14:23, 52:61)

protein1<-read.table("C:/protein/TDP43/main/protein_residue.txt", header = T)
protein2<-c("TSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNS")
protein2<-strsplit(protein2, "")[[1]]
DNA<-c("G","T","T","G","A","G","C","G","T","T")

colo<-c("orange", "blue", "yellow")
models<-c("Enet", "RF", "Average")

min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}


colo<-c("red", "blue", "yellow")

min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

system_map_DNA_prot<-list()
system_map_prot_prot<-list()

for(kk in seq_along(systems)){
  
  
  Enet<-read.csv(system_mat[,kk][1], row.names = 1)
  RF<-read.csv(system_mat[,kk][4], row.names = 1)
  #pca<-read.csv(system_mat[,kk][3], row.names = 1)
  
  model<-list(Enet, RF)#, pca)
  my_map_DNA_prot<-list()
  my_map_prot_prot<-list()
  
  for(k in seq_along(model)){
    
    residue_number<-as.numeric(unlist(regmatches( rownames(model[[k]]), 
                                                  gregexpr("[[:digit:]]+",  
                                                           rownames(model[[k]])))))
    
    imp_inter<-cbind(model[[k]], 
                     matrix(residue_number, 
                            ncol=2,
                            byrow = T))
    
    map<-matrix(0, 88, 88)
    for(i in 1:nrow(imp_inter)){
      
      map[imp_inter[i,2],imp_inter[i,3]] <- imp_inter[i,1]
    }
    
    map<-map+t(map)
    
    my_map_DNA_prot[[k]]<-map[79:88,1:78]
    my_map_prot_prot[[k]]<-map[1:78,1:78]
    
  }
  
  my_map_DNA_prot[[length(model)+1]]<-Reduce("+",my_map_DNA_prot)/length(model)
  my_map_prot_prot[[length(model)+1]]<-Reduce("+",my_map_prot_prot)/length(model)
  
  names(my_map_DNA_prot)<-models
  names(my_map_prot_prot)<-models
  
  system_map_DNA_prot[[kk]]<-my_map_DNA_prot
  system_map_prot_prot[[kk]]<-my_map_prot_prot
  
}

names(system_map_DNA_prot)<-systems
names(system_map_prot_prot)<-systems

#-------------------------------------------------------------------------

length(system_map_DNA_prot)

my_palette <-colorRampPalette(c("lightskyblue1", 
                                "navyblue"))(n = 100)
#my_palette <-colorspace::sequential_hcl

for(i in seq_along(system_map_DNA_prot)){
  for(j in seq_along(lengths(system_map_DNA_prot[[1]]))){
    
    interaction_DNA_prot<-system_map_DNA_prot[[i]][[j]]
    interaction_prot_prot<-system_map_prot_prot[[i]][[j]]
    
    heatmap.2( interaction_DNA_prot, trace="none", 
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
              ylab = "DNA_Residue",
              labRow = DNA,
              labCol = protein2,
              main = paste0(systems[i],"-",models[j],"-DNA-Prot")
              
    )
    
    heatmap.2( interaction_prot_prot, trace="none", 
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
               labRow = protein2,
               labCol = protein2,
               main = paste0(systems[i],"-",models[j],"-prot-Prot")
               
    )
    
  }
}





