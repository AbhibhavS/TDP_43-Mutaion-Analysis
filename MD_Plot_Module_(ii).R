
#-----------Prortin DNA Map--------------------------------#

setwd("K:/TDP43_r1")
set.seed(123)
source("MD_utility.R")
library(tidyverse)
library(gplots)

system<-dir()[grep(".pdb", dir())]
system<-system[-grep("distance", system)]

distance.matrix<-list()
contact.matrix<-list()


for(kk in 1:4){
  
  print(paste("start", kk, "at", Sys.time(), sep=" "))
  
  
  data<-read.table(system[kk], sep="\t")
  
  xx<-mdml.matrix(data, type = "MD")
  
  distance_matrix <- xx$distance_matrix
  contact_matrix <- xx$contact_matrix
  
  distance.matrix[[kk]] <- Reduce("+", distance_matrix)/length(distance_matrix)
  contact.matrix[[kk]] <- Reduce("+", contact_matrix)/length(contact_matrix)
  
}
names(distance.matrix)<- system
names(contact.matrix)<- system

#---------------------PLOTS---------------------------------#
#------------------------------------------------------------#
#------------------------------------------------------------#



#---------------Protein-Protein Average Inter-residual Distance-----------------------------------#

protein2<-c("TSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNS")
protein_wt<-strsplit(protein2, "")[[1]]
DNA_seq<-c("G","T","T","G","A","G","C","G","T","T")





for(i in 1:4){
  
  my_palette <-colorRampPalette(c("white","skyblue","blue"
                                  ,"orange","red","black"))(n = 400)
  map<-distance.matrix[[i]][1:78,1:78]
  
  
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Average Inter-residual Distance",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = paste(system[i]),
            labCol = protein_wt,
            labRow = protein_wt)
  
  
  #-----------------Protein-DNA Average Inter-residual Distance--------------------
  
  my_palette <-colorRampPalette(c("white","skyblue","blue"
                                  ,"orange","red","black"))(n = 400)
  
  map<-distance.matrix[[i]][79:88,1:78]
  
  
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Protein-DNA Avg residual Distance",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "DNA_Residue",
            main = system[i],
            labCol = protein_wt,
            labRow = DNA_seq)
  
  
  #------------------Contact Protein-Protein--------------#
  
  my_palette <-colorRampPalette(c("pink",
                                  "red", "black"))(n = 50)
  
  current_map<-contact.matrix[[i]][1:78,1:78]
  neig<-4
  
  for(j in 1:ncol(current_map)){
    
    if(j-neig < 0){
      current_map[1:(j+neig),j]<-0
      next
    }
    else if(j+neig > ncol(current_map)){
      current_map[(j-neig):ncol(current_map),j]<-0
      next
    }
    else{
      current_map[(j-neig):(j+neig),j]<-0
      next
    }
    
  }
  
  heatmap.2(current_map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Protein residues contact intensity",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = system[i],
            labCol = protein_wt,
            labRow = protein_wt)
  
  
  #------------------Contact Protein-DNA--------------#
  
  my_palette <-colorRampPalette(c("pink",
                                  "red", "black"))(n = 50)
  
  current_map<-contact.matrix[[i]][79:88,1:78]
  
  heatmap.2(current_map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.5, cexCol = 0.5, 
            key.xlab = "Protein-DNA residues contact intensity",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            ylab = "DNA_Residue",
            xlab = "Protein_Residue",
            main = system[i],
            labCol = protein_wt,
            labRow = DNA_seq)
  
  
}


#-------------------------Difference Distance Map---------------------

my_palette <-colorRampPalette(c("red",
                                "white","darkblue"))(n = 100)
for(i in 1:3){
  map<-distance.matrix[[i]][1:78,1:78]- distance.matrix[[4]][1:78,1:78]
  
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
            main = system[i],
            labCol = protein_wt,
            labRow = protein_wt)
}




















